"""Request-local diagnostics returned only to the requesting chat, never saved on disk.

A step's events are retained in session-scoped memory for recovery, including
completed steps whose response may have been lost. Reads are non-destructive
snapshots; entries expire after ten minutes without new events.
"""
from contextvars import ContextVar
from datetime import datetime, timezone
from functools import wraps
from inspect import signature
import re
import threading
import time
from time import monotonic

from pydantic import BaseModel
from app.services.assistant_errors import public_assistant_error

_active = ContextVar("assistant_trace", default=None)
_private = re.compile(r"api.?key|authorization|password|secret|credential|signature", re.I)

_PARTIAL: dict[tuple[str, str], dict] = {}
_PARTIAL_LOCK = threading.Lock()
_PARTIAL_TTL_SECONDS = 600


def _session_key(progress_id):
    from app.services.session_context import get_session_id
    return (get_session_id(), progress_id)


def _expire_partials(now: float) -> None:
    for key in [key for key, entry in _PARTIAL.items() if now - entry["updated"] > _PARTIAL_TTL_SECONDS]:
        _PARTIAL.pop(key, None)


def take_partial(progress_id: str) -> list:
    """A non-destructive snapshot of a running or completed step."""
    with _PARTIAL_LOCK:
        _expire_partials(time.monotonic())
        entry = _PARTIAL.get(_session_key(progress_id))
        return list(entry["events"]) if entry else []


def record_provider_reasoning(reply):
    """Record what the provider reply says about itself: token usage when the
    reply reports it (OpenAI/local `usage`, Anthropic `usage`, Gemini
    `usageMetadata`, Ollama eval counts) and reasoning text only when the
    provider returned some. Replies without either (model listings, capability
    probes) add nothing."""
    trace = _active.get()
    if trace is None or not isinstance(reply, dict):
        return
    usage = reply.get("usage") if isinstance(reply.get("usage"), dict) else reply.get("usageMetadata")
    if not isinstance(usage, dict):
        counts = {key: reply[key] for key in ("prompt_eval_count", "eval_count") if isinstance(reply.get(key), int)}
        usage = counts or None
    if usage:
        trace.record("Provider usage", {"model": reply.get("model"), "usage": usage})
    found = []

    def visit(value):
        if isinstance(value, dict):
            for key, item in value.items():
                if key in {"reasoning_content", "reasoning", "thinking"} and isinstance(item, str) and item:
                    found.append(item)
                elif key == "text" and value.get("thought") is True and isinstance(item, str):
                    found.append(item)
                elif key == "text" and value.get("type") in {"reasoning_text", "summary_text"} and isinstance(item, str):
                    found.append(item)
                elif isinstance(item, (dict, list)):
                    visit(item)
        elif isinstance(value, list):
            for item in value:
                visit(item)

    visit(reply)
    if found:
        trace.record("Provider-supplied reasoning", {"text": found})


def clean(value, secrets=()):
    if isinstance(value, BaseModel):
        value = value.model_dump(mode="json")
    if isinstance(value, dict):
        return {str(key): "[excluded]" if _private.search(str(key)) or key in
                {"image", "plotImage", "ollamaUrl", "localServerUrl"} else clean(item, secrets)
                for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(item, secrets) for item in value]
    if isinstance(value, str):
        for secret in secrets:
            if secret:
                value = value.replace(secret, "[redacted]")
        return re.sub(r"data:image/[^\s\"']+", "[image excluded]", value)
    if value is None or isinstance(value, (bool, int, float)):
        return value
    return f"[{type(value).__name__}]"


class Trace:
    def __init__(self, *secrets, trace_id: str | None = None):
        self.events: list[dict[str, object]] = []
        self.secrets = secrets
        self.trace_id = trace_id

    def record(self, kind, data):
        self.events.append({"time": datetime.now(timezone.utc).isoformat(),
                            "kind": kind, "data": clean(data, self.secrets)})
        if self.trace_id:
            now = time.monotonic()
            with _PARTIAL_LOCK:
                _expire_partials(now)
                entry = _PARTIAL.setdefault(_session_key(self.trace_id), {"events": []})
                entry["events"] = self.events
                entry["updated"] = now

    def __enter__(self):
        self.token = _active.set(self)
        return self

    def __exit__(self, *args):
        _active.reset(self.token)

    def finish(self) -> None:
        """Retain the snapshot: returning a response does not confirm delivery."""
        if self.trace_id:
            with _PARTIAL_LOCK:
                entry = _PARTIAL.get(_session_key(self.trace_id))
                if entry is not None:
                    entry["updated"] = time.monotonic()


def traced(function):
    """Record invoked assistant functions and outcomes, including failed calls."""
    parameters = signature(function)

    @wraps(function)
    def wrapper(*args, **kwargs):
        trace = _active.get()
        if trace is None:
            return function(*args, **kwargs)
        call_id = len(trace.events) + 1
        arguments = dict(parameters.bind(*args, **kwargs).arguments)
        # The step request is recorded once by the browser. Never duplicate credentials.
        for key in ("payload", "request", "schema", "check_cancelled"):
            arguments.pop(key, None)
        trace.record("function started", {"id": call_id, "function": function.__name__, "arguments": arguments})
        start = monotonic()
        try:
            result = function(*args, **kwargs)
        except Exception as exc:
            trace.record("function failed", {"id": call_id, "function": function.__name__,
                         "error": public_assistant_error(exc), "seconds": monotonic() - start})
            raise
        trace.record("function returned", {"id": call_id, "function": function.__name__,
                     "result": result, "seconds": monotonic() - start})
        return result
    return wrapper
