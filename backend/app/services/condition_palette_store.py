from __future__ import annotations

import re

from app.schemas.annotation import AnnotationKind
from app.services.runtime_cache import invalidate_runtime_cache
from app.services.session_context import get_session_id

from app.services.copylot_core import COPYLOT_COLORS, UKJ_COLORS, MATPLOTLIB_TAB10

HEX_COLOR_RE = re.compile(r"^#[0-9a-fA-F]{6}$")


def _empty_palette_state() -> dict[AnnotationKind, dict[str, str]]:
    return {
        "protein": {},
        "phospho": {},
        "phosprot": {},
    }


_PALETTE_SESSIONS: dict[str, dict[AnnotationKind, dict[str, str]]] = {}


def _session_palettes() -> dict[AnnotationKind, dict[str, str]]:
    session_id = get_session_id()
    if session_id not in _PALETTE_SESSIONS:
        _PALETTE_SESSIONS[session_id] = _empty_palette_state()
    return _PALETTE_SESSIONS[session_id]


def _normalize_conditions(conditions: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for condition in conditions:
        normalized = str(condition).strip()
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        ordered.append(normalized)
    return ordered


def _generated_color_hex(index: int) -> str:
    return COPYLOT_COLORS[index % len(COPYLOT_COLORS)]


def _normalize_hex(value: str) -> str | None:
    text = str(value).strip()
    if not text:
        return None
    if not text.startswith("#"):
        text = f"#{text}"
    if not HEX_COLOR_RE.match(text):
        return None
    return text.upper()


def get_condition_palette(kind: AnnotationKind) -> dict[str, str]:
    return dict(_session_palettes().get(kind, {}))


def get_category_palette(kind: AnnotationKind) -> list[str]:
    """The selected preset in full, or the ordered custom condition colors."""
    overrides = get_condition_palette(kind)
    if not overrides:
        return list(COPYLOT_COLORS)
    from app.services.annotation_store import get_annotation

    annotation = get_annotation(kind)
    conditions = list(overrides)
    if annotation is not None and "condition" in annotation.metadata:
        conditions = _normalize_conditions(annotation.metadata["condition"].astype(str).tolist())
    colors = list(build_condition_color_map(kind, conditions).values())
    for preset in (COPYLOT_COLORS, UKJ_COLORS, MATPLOTLIB_TAB10):
        if colors and all(color == preset[i % len(preset)].upper() for i, color in enumerate(colors)):
            return list(preset)
    return colors or list(COPYLOT_COLORS)


def fill_category_palette(function, settings):
    """Resolve the GUI choice before rendering or capturing notebook settings."""
    if hasattr(function, "_category_palette_kind") and settings.get("plot_colors") is None:
        kind = settings.get("kind") or function._category_palette_kind
        settings["plot_colors"] = get_category_palette(kind)
    return settings


def apply_category_palettes(namespace, names, kind="protein"):
    """Bind session colors at renderer entry; exported plot bodies stay pure."""
    import inspect
    from functools import wraps

    def decorate(function):
        signature = inspect.signature(function)

        @wraps(function)
        def wrapped(*args, **kwargs):
            bound = signature.bind(*args, **kwargs)
            bound.apply_defaults()
            fill_category_palette(wrapped, bound.arguments)
            return function(*bound.args, **bound.kwargs)

        wrapped._category_palette_kind = kind
        return wrapped

    for name in names:
        namespace[name] = decorate(namespace[name])


def set_condition_palette(kind: AnnotationKind, palette: dict[str, str]) -> dict[str, str]:
    normalized: dict[str, str] = {}
    for condition, color in palette.items():
        condition_name = str(condition).strip()
        if not condition_name:
            continue
        hex_color = _normalize_hex(color)
        if hex_color is None:
            continue
        normalized[condition_name] = hex_color
    _session_palettes()[kind] = normalized
    invalidate_runtime_cache(f"palette:{kind}:updated")
    return dict(normalized)


def build_condition_color_map(kind: AnnotationKind, conditions: list[str]) -> dict[str, str]:
    ordered = _normalize_conditions(conditions)
    overrides = _session_palettes().get(kind, {})
    result: dict[str, str] = {}
    for index, condition in enumerate(ordered):
        override = _normalize_hex(overrides.get(condition, ""))
        result[condition] = override if override is not None else _generated_color_hex(index)
    return result


def reset_condition_palettes(*, all_sessions: bool = False) -> None:
    if all_sessions:
        _PALETTE_SESSIONS.clear()
    else:
        _PALETTE_SESSIONS.pop(get_session_id(), None)
    invalidate_runtime_cache("palette:reset-all", all_sessions=all_sessions)


def snapshot_condition_palette_session_state():
    """Detach-safe reference used before the session entry is atomically replaced."""
    return _PALETTE_SESSIONS.get(get_session_id())


def restore_condition_palette_session_state(snapshot) -> None:
    session_id = get_session_id()
    if snapshot is None:
        _PALETTE_SESSIONS.pop(session_id, None)
    else:
        _PALETTE_SESSIONS[session_id] = snapshot
    invalidate_runtime_cache("palette:transaction-restored")
