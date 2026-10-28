"""Diagnostic traces must preserve outcomes without mixing users or retaining credentials."""
import json
from concurrent.futures import ThreadPoolExecutor

import pytest

from app.services.assistant_trace import Trace, traced
from app.services.assistant_errors import DEFAULT_ERROR, public_assistant_error


@traced
def operation(value, fail=False):
    if fail:
        raise ValueError(value)
    return {"text": value, "api_key": "hidden", "thinking": "private", "image": "pixels"}


def test_trace_captures_order_and_redacts_nested_results_and_errors():
    with Trace("test-credential") as trace:
        operation("hello test-credential")
        with pytest.raises(ValueError):
            operation("error test-credential", fail=True)
    assert [e["kind"] for e in trace.events] == [
        "function started", "function returned", "function started", "function failed"]
    encoded = json.dumps(trace.events)
    assert "hello [redacted]" in encoded
    assert "error [redacted]" in encoded
    assert trace.events[-1]["data"]["error"] == DEFAULT_ERROR
    assert all(secret not in encoded for secret in ["test-credential", "hidden", "pixels"])
    assert "private" in encoded  # provider-returned thinking is retained
    operation("outside")
    assert len(trace.events) == 4


def test_parallel_traces_are_isolated():
    def run(value):
        with Trace() as trace:
            operation(value)
        return json.dumps(trace.events)
    with ThreadPoolExecutor(max_workers=2) as pool:
        first, second = list(pool.map(run, ["first-user", "second-user"]))
    assert "second-user" not in first
    assert "first-user" not in second


def test_provider_reasoning_including_streams_is_retained_separately():
    import io
    from app.services.provider_http import read_json_response
    from app.services.local_ai_server import _read_stream
    from app.services.assistant_trace import record_provider_reasoning

    with Trace() as trace:
        read_json_response(io.BytesIO(json.dumps({"content": [
            {"type": "thinking", "thinking": "provider thought", "signature": "not-exported"}
        ]}).encode()), "invalid")
        stream = b''.join(b'data: ' + json.dumps({"choices": [{"delta": {
            "reasoning_content": piece}}]}).encode() + b'\n' for piece in ["first ", "second"])
        reply = _read_stream(io.BytesIO(stream), None)
        record_provider_reasoning(reply)
        record_provider_reasoning({"candidates": [{"content": {"parts": [{"thought": True, "text": "summary"}]}}]})
        record_provider_reasoning({"choices": [{"message": {"content": "answer only"}}]})
    assert [event["data"]["text"] for event in trace.events] == [["provider thought"], ["first second"], ["summary"]]
    assert all(event["kind"] == "Provider-supplied reasoning" for event in trace.events), "a reply without reasoning adds no event"
    assert "not-exported" not in json.dumps(trace.events)


def test_provider_usage_is_recorded_in_every_shape(monkeypatch):
    from app.services.assistant_trace import record_provider_reasoning
    with Trace() as trace:
        record_provider_reasoning({"model": "gpt", "usage": {"prompt_tokens": 10, "completion_tokens": 3}})
        record_provider_reasoning({"usageMetadata": {"promptTokenCount": 7, "candidatesTokenCount": 2}})
        record_provider_reasoning({"model": "llama", "prompt_eval_count": 5, "eval_count": 1, "response": "x"})
        record_provider_reasoning({"models": [{"name": "listing"}]})
    assert [event["kind"] for event in trace.events] == ["Provider usage"] * 3
    assert trace.events[0]["data"] == {"model": "gpt", "usage": {"prompt_tokens": 10, "completion_tokens": 3}}
    assert trace.events[2]["data"]["usage"] == {"prompt_eval_count": 5, "eval_count": 1}


def test_running_and_completed_snapshots_survive_reads(monkeypatch):
    from app.services import assistant_trace
    with Trace("secret-value", trace_id="abandoned-step-1") as trace:
        operation("one secret-value")
        operation("two")
    partial = assistant_trace.take_partial("abandoned-step-1")
    assert [event["kind"] for event in partial] == ["function started", "function returned"] * 2
    assert "secret-value" not in json.dumps(partial) and "[redacted]" in json.dumps(partial)
    assert assistant_trace.take_partial("abandoned-step-1") == partial
    trace.record("later event", {})
    assert len(assistant_trace.take_partial("abandoned-step-1")) == len(partial) + 1
    with Trace(trace_id="finished-step-1") as trace:
        operation("three")
    trace.finish()
    assert len(assistant_trace.take_partial("finished-step-1")) == 2
    assert assistant_trace.take_partial("never-seen") == []


def test_route_returns_trace_on_success_and_failure(monkeypatch):
    from fastapi.testclient import TestClient
    from app.main import app
    from app.api import routes_chat
    from app.schemas.assistant import AssistantResponse

    monkeypatch.setattr(routes_chat, "ai_mode_enabled", lambda: True)

    def step(payload, cancelled):
        operation(payload.message, fail=payload.message == "fail")
        return AssistantResponse(model="test", revision=1, answer="done")

    monkeypatch.setattr(routes_chat, "assistant_step", step)
    with TestClient(app) as client:
        for message, expected in [("hello", 200), ("fail", 400)]:
            response = client.post("/api/chat/assistant", json={"message": message, "traceId": "completed-" + message,
                "view": {"module": "qc.coverage", "workflow": "protein", "controls": [], "destinations": {}}})
            assert response.status_code == expected, response.text
            assert len(response.json()["trace"]) == 2
            recovered = client.get("/api/chat/assistant-trace/completed-" + message)
            assert recovered.json()["trace"] == response.json()["trace"]
        # Recovery snapshots can be fetched repeatedly.
        monkeypatch.setattr(routes_chat, "assistant_step", lambda payload, cancelled: (_ for _ in ()).throw(RuntimeError("stopped")))
        from app.services.assistant_trace import Trace as _Trace
        with _Trace(trace_id="abcdefgh-1234") as trace:
            operation("partial")
        assert client.get("/api/chat/assistant-trace/abcdefgh-1234").json()["trace"] == trace.events
        assert client.get("/api/chat/assistant-trace/abcdefgh-1234").json()["trace"] == trace.events
        assert client.get("/api/chat/assistant-trace/bad id").status_code == 400


def test_recovery_is_session_scoped_and_expires(monkeypatch):
    from app.services import assistant_trace, session_context
    now = [1000.0]
    monkeypatch.setattr(assistant_trace.time, "monotonic", lambda: now[0])
    monkeypatch.setenv("COPYLOT_SESSION_ISOLATION", "1")
    token = session_context.set_session_id("first-session")
    try:
        with Trace(trace_id="shared-trace-id") as trace:
            trace.record("event", "first user")
        trace.finish()
        other = session_context.set_session_id("second-session")
        try:
            assert assistant_trace.take_partial("shared-trace-id") == []
        finally:
            session_context.reset_session_id(other)
        assert len(assistant_trace.take_partial("shared-trace-id")) == 1
        now[0] += 601
        assert assistant_trace.take_partial("shared-trace-id") == []
    finally:
        session_context.reset_session_id(token)


@pytest.mark.parametrize("exception_type", [ValueError, RuntimeError])
def test_internal_exception_details_never_reach_response_or_recovered_trace(monkeypatch, exception_type):
    from fastapi.testclient import TestClient
    from app.main import app
    from app.api import routes_chat

    internal = "parser failed at /srv/private/patient.tsv with credential=unregistered-secret"

    @traced
    def fail_inside_library():
        raise exception_type(internal)

    monkeypatch.setattr(routes_chat, "ai_mode_enabled", lambda: True)
    monkeypatch.setattr(routes_chat, "assistant_step", lambda *args: fail_inside_library())
    with TestClient(app) as client:
        response = client.post("/api/chat/assistant", json={
            "message": "Explain this plot", "traceId": "private-error-123",
            "view": {"module": "qc.coverage"},
        })
        assert response.status_code == (400 if exception_type is ValueError else 500)
        assert response.json()["trace"][-1]["data"]["error"] == DEFAULT_ERROR
        recovered = client.get("/api/chat/assistant-trace/private-error-123")
        assert recovered.json()["trace"] == response.json()["trace"]
        for body in (response.text, recovered.text):
            assert "patient.tsv" not in body
            assert "unregistered-secret" not in body
            assert "parser failed" not in body


@pytest.mark.parametrize("message", [
    "The analysis changed while the assistant was working. Please retry.",
    "The assistant request was cancelled.",
    "OpenAI request failed (HTTP 401).",
])
def test_approved_errors_remain_useful_but_cannot_carry_extra_details(message):
    assert public_assistant_error(ValueError(message)) == message
    assert public_assistant_error(ValueError(message + " /srv/private/patient.tsv")) == DEFAULT_ERROR
    assert public_assistant_error(ValueError(message, "private detail")) == DEFAULT_ERROR
    assert public_assistant_error(RuntimeError(message)) == DEFAULT_ERROR


def test_provider_timeouts_remain_actionable_without_copying_exception_text():
    from app.services.local_ai_server import _timeout_message
    assert "COPYLOT_LOCAL_AI_TIMEOUT" in public_assistant_error(ValueError(_timeout_message()))
    assert "retained" in public_assistant_error(ValueError(
        "Anthropic did not answer within 300 seconds; this step was stopped. Earlier completed actions are retained."
    ))
