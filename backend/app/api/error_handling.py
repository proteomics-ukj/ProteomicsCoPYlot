"""Common API error translation, retaining each caller's logger."""
from collections.abc import Iterator
from contextlib import contextmanager
import logging

from fastapi import HTTPException


@contextmanager
def api_errors(message: str, *, logger: logging.Logger) -> Iterator[None]:
    try:
        yield
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e
    except Exception as e:
        logger.exception(message)
        raise HTTPException(status_code=500, detail=f"{message}.") from e
