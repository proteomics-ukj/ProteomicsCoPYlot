"""Bounded, session-scoped preparation reuse for large peptide tables.

Values never escape without a deep copy. The key includes the session revision
and source file identity. Prepared inputs exclude cosmetic settings; rendered
image entries include their exact settings and output format.
"""
from collections import OrderedDict
from copy import deepcopy
from threading import RLock
from typing import Callable, TypeVar, cast
import sys

import pandas as pd

from app.services.runtime_cache import runtime_cache_version
from app.services.session_context import get_session_id

T = TypeVar('T')
_LOCK = RLock()
_CACHE: OrderedDict[tuple, tuple[object, int]] = OrderedDict()
_MAX_BYTES = 256 * 1024 * 1024
_MAX_ENTRIES = 16


def _size(value: object) -> int:
    if isinstance(value, pd.DataFrame):
        return int(value.memory_usage(deep=True).sum())
    if isinstance(value, pd.Series):
        return int(value.memory_usage(deep=True))
    if isinstance(value, dict):
        return sys.getsizeof(value) + sum(_size(k) + _size(v) for k, v in value.items())
    if isinstance(value, (tuple, list)):
        return sys.getsizeof(value) + sum(_size(v) for v in value)
    return sys.getsizeof(value)


def cached_peptide_preparation(name: str, options: tuple, factory: Callable[[], T]) -> T:
    from app.services.peptide_tools import _peptide_dataset_path
    from app.services.dataset_store import get_current_dataset

    stored = get_current_dataset('peptide')
    if stored is None:
        return factory()
    path = _peptide_dataset_path(stored)
    stat = path.stat()
    session = get_session_id()
    revision = runtime_cache_version()
    source = (str(path.resolve()), stat.st_mtime_ns, stat.st_size, stat.st_ctime_ns)
    key = (session, revision, source, name, options)
    with _LOCK:
        # Release obsolete preparations for this session, including external edits.
        for old in list(_CACHE):
            if old[0] == session and old[1:3] != key[1:3]:
                del _CACHE[old]
        if key in _CACHE:
            value, _ = _CACHE[key]
            _CACHE.move_to_end(key)
            return cast(T, deepcopy(value))
    value = factory()
    size = _size(value)
    current_stat = path.stat()
    unchanged = (current_stat.st_mtime_ns, current_stat.st_size, current_stat.st_ctime_ns) == source[1:]
    if size <= _MAX_BYTES and runtime_cache_version() == revision and unchanged:
        with _LOCK:
            _CACHE[key] = (value, size)
            while len(_CACHE) > _MAX_ENTRIES or sum(entry[1] for entry in _CACHE.values()) > _MAX_BYTES:
                _CACHE.popitem(last=False)
        return deepcopy(value)
    return value
