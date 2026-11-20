"""The backend follows the code-shape rules in tools/check_backend_rules.py.

The rules keep exception text out of HTTP responses, make every 500 leave a log
entry, and keep filesystem paths out of the messages a 400 shows the user. They
are about how code is written, so they are checked on the syntax tree, the way
frontend/eslint.config.js checks the frontend's. Each rule is also run on a
snippet that breaks it, so a checker that silently stopped matching fails here.
"""
from __future__ import annotations

import importlib.util
import sys

import pytest

from tests.repo_paths import ROOT


def _checker():
    spec = importlib.util.spec_from_file_location("check_backend_rules", ROOT / "tools" / "check_backend_rules.py")
    module = importlib.util.module_from_spec(spec)
    # Registered first: the checker's dataclass looks its module up by name.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_the_backend_follows_every_rule() -> None:
    problems = _checker().violations(ROOT)
    assert not problems, "\n".join(str(problem) for problem in problems)


ROUTE_HEADER = "import logging\nfrom fastapi import HTTPException\nLOGGER = logging.getLogger(__name__)\n"


@pytest.mark.parametrize("filename", ["routes_example.py", "plots_qc.py", "plot_responses.py", "nested/helper.py"])
def test_all_api_modules_are_scanned(tmp_path, filename):
    path = tmp_path / "backend" / "app" / "api" / filename
    path.parent.mkdir(parents=True)
    path.write_text(
        ROUTE_HEADER + "def route():\n    try:\n        pass\n"
        "    except Exception:\n        raise HTTPException(status_code=500, detail='Failed.')\n",
        encoding="utf-8",
    )
    problems = _checker().violations(tmp_path)
    assert [(problem.path, problem.rule) for problem in problems] == [
        (str(path.relative_to(tmp_path)), "unlogged-500"),
    ]


@pytest.mark.parametrize("body,rule", [
    # The exception's own text in the response, however it is spelled.
    ("    except Exception as e:\n        LOGGER.exception('x')\n"
     "        raise HTTPException(status_code=500, detail=f'Failed: {e}') from e\n", "exception-in-detail"),
    ("    except ValueError as exc:\n"
     "        raise HTTPException(status_code=400, detail=f'Bad input ({str(exc)})') from exc\n", "exception-in-detail"),
    # A 500 with nothing in the log.
    ("    except Exception as e:\n"
     "        raise HTTPException(status_code=500, detail='Failed.') from e\n", "unlogged-500"),
    # Logged only after the response was raised is not logged.
    ("    except Exception as e:\n        if e:\n"
     "            raise HTTPException(status_code=500, detail='Failed.') from e\n        LOGGER.exception('x')\n",
     "unlogged-500"),
])
def test_each_route_rule_catches_what_it_names(body: str, rule: str) -> None:
    source = ROUTE_HEADER + "def route():\n    try:\n        pass\n" + body
    assert rule in {problem.rule for problem in _checker().check_route_module(source, "routes_x.py")}


@pytest.mark.parametrize("body", [
    "    except Exception as e:\n        LOGGER.exception('x')\n"
    "        raise HTTPException(status_code=500, detail='Failed to store the file.') from e\n",
    "    except ValueError as e:\n        raise HTTPException(status_code=400, detail=str(e)) from e\n",
])
def test_the_route_rules_allow_the_house_pattern(body: str) -> None:
    source = ROUTE_HEADER + "def route():\n    try:\n        pass\n" + body
    assert _checker().check_route_module(source, "routes_x.py") == []


@pytest.mark.parametrize("message,flagged", [
    ("f'cannot read {path}'", True),
    ("f'cannot read {str(full_path)}'", True),
    ("f'searched {searched}'", True),
    ("f'missing in {candidate.parent}'", True),
    ("f'The file {path.name!r} is no longer available.'", False),
    ("'No dataset is loaded.'", False),
])
def test_the_path_rule_tells_a_location_from_a_file_name(message: str, flagged: bool) -> None:
    source = f"def load(path, full_path, searched, candidate):\n    raise ValueError({message})\n"
    problems = _checker().check_service_module(source, "services/x.py")
    assert bool(problems) is flagged, problems


@pytest.mark.parametrize("signature", [
    "def route(*, logger: logging.Logger)",
    "def route(logger: logging.Logger, /)",
    "async def route(logger: logging.Logger)",
])
def test_the_route_rules_allow_a_typed_injected_logger(signature: str) -> None:
    source = (
        "import logging\nfrom fastapi import HTTPException\n"
        f"{signature}:\n    try:\n        pass\n"
        "    except Exception as exc:\n        logger.exception('Failed')\n"
        "        raise HTTPException(status_code=500, detail='Failed.') from exc\n"
    )
    assert _checker().check_route_module(source, "helper.py") == []


@pytest.mark.parametrize("parameter,body", [
    ("logger", "        logger.exception('Failed')\n"),
    ("logger: str", "        logger.exception('Failed')\n"),
    ("logger: logging.Logger", "        pass\n"),
    ("logger: logging.Logger", "        other.exception('Failed')\n"),
    ("logger: logging.Logger", "        def later():\n            logger.exception('Failed')\n"),
    ("LOGGER: logging.Logger", "        def later():\n            LOGGER.exception('Failed')\n"),
])
def test_injected_logger_must_be_typed_and_called_in_the_handler(parameter: str, body: str) -> None:
    source = (
        "import logging\nfrom fastapi import HTTPException\n"
        f"def route({parameter}):\n    try:\n        pass\n    except Exception:\n"
        + body + "        raise HTTPException(status_code=500, detail='Failed.')\n"
    )
    assert "unlogged-500" in {p.rule for p in _checker().check_route_module(source, "helper.py")}


def test_injected_logger_must_precede_the_response() -> None:
    source = (
        "import logging\nfrom fastapi import HTTPException\n"
        "def route(logger: logging.Logger):\n    try:\n        pass\n    except Exception:\n"
        "        raise HTTPException(status_code=500, detail='Failed.')\n"
        "        logger.exception('Too late')\n"
    )
    assert "unlogged-500" in {p.rule for p in _checker().check_route_module(source, "helper.py")}


def test_a_typed_logger_in_another_function_does_not_authorize_logging() -> None:
    source = (
        "import logging\nfrom fastapi import HTTPException\n"
        "def unrelated(logger: logging.Logger):\n    pass\n"
        "def route(logger):\n    try:\n        pass\n    except Exception:\n"
        "        logger.exception('Failed')\n"
        "        raise HTTPException(status_code=500, detail='Failed.')\n"
    )
    assert "unlogged-500" in {p.rule for p in _checker().check_route_module(source, "helper.py")}


def test_a_typed_logger_does_not_allow_exception_text_in_the_response() -> None:
    source = (
        "import logging\nfrom fastapi import HTTPException\n"
        "def route(logger: logging.Logger):\n    try:\n        pass\n    except Exception as exc:\n"
        "        logger.exception('Failed')\n"
        "        raise HTTPException(status_code=500, detail=f'Failed: {exc}') from exc\n"
    )
    assert "exception-in-detail" in {p.rule for p in _checker().check_route_module(source, "helper.py")}
