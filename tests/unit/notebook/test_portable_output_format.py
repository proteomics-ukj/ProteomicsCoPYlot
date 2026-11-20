from contextvars import ContextVar
from types import ModuleType

import pytest

from app.services.notebook_source import portable_output_format


def renderer():
    module = ModuleType("test_renderer")
    module.plot_output_format = ContextVar("test_output_format", default="png")
    return module


def test_restores_previous_format_after_rendering():
    module = renderer()
    module.plot_output_format.set("svg")
    with portable_output_format(module, "pdf"):
        assert module.plot_output_format.get() == "pdf"
    assert module.plot_output_format.get() == "svg"


def test_restores_format_and_propagates_rendering_error():
    module = renderer()
    error = RuntimeError("render failed")
    with pytest.raises(RuntimeError) as caught:
        with portable_output_format(module, "svg"):
            assert module.plot_output_format.get() == "svg"
            raise error
    assert caught.value is error
    assert module.plot_output_format.get() == "png"


def test_nested_failure_restores_outer_format_then_original():
    module = renderer()
    other = renderer()
    with portable_output_format(module, "svg"):
        with pytest.raises(ValueError, match="inner render failed"):
            with portable_output_format(module, "pdf"):
                assert module.plot_output_format.get() == "pdf"
                assert other.plot_output_format.get() == "png"
                raise ValueError("inner render failed")
        assert module.plot_output_format.get() == "svg"
        with portable_output_format(module, "png"):
            assert module.plot_output_format.get() == "png"
        assert module.plot_output_format.get() == "svg"
    assert module.plot_output_format.get() == "png"
