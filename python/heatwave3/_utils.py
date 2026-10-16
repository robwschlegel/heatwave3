"""Internal helpers: output paths, input resolution, dates and argument checks."""

from __future__ import annotations

import os
import warnings
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd

from . import _core

# Julian Day Number of 1970-01-01.
JD_UNIX_EPOCH = 2440588


# Output paths: a single `name` stem maps to the four product files.

def clim_path(name) -> str:
    return f"{os.fspath(name)}_clim.nc"


def events_path(name) -> str:
    return f"{os.fspath(name)}_events.nc"


def daily_path(name) -> str:
    return f"{os.fspath(name)}_events_daily.nc"


def proto_path(name) -> str:
    return f"{os.fspath(name)}_protoevents.nc"


def check_out_dir(path: str) -> None:
    out_dir = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(out_dir):
        raise FileNotFoundError(f"Output directory does not exist: {os.path.dirname(path)}")


# Inputs

def is_multi(file_in) -> bool:
    """True for a directory or a sequence of paths (as opposed to one file)."""
    return not isinstance(file_in, (str, os.PathLike))


def list_nc(directory, what: str = "the specified directory") -> list[str]:
    files = sorted(str(p) for p in Path(directory).iterdir()
                   if p.suffix in (".nc", ".nc4") and p.is_file())
    if not files:
        raise FileNotFoundError(f"No .nc or .nc4 files found in {what}.")
    return files


def resolve_input(file_in, missing_msg: str = "Input file does not exist: ",
                  dir_what: str = "the specified directory") -> tuple[list[str], bool]:
    """Resolve a file, directory or sequence of files to (files, multi_file)."""
    if not is_multi(file_in):
        path = os.fspath(file_in)
        if os.path.isdir(path):
            return list_nc(path, dir_what), True
        if not os.path.exists(path):
            raise FileNotFoundError(f"{missing_msg}{path}")
        return [path], False
    files = [os.fspath(f) for f in file_in]
    if len(files) == 1:
        return resolve_input(files[0], missing_msg, dir_what)
    return files, True


# Dates

def to_datetime64(dates) -> np.ndarray:
    """Convert date-like input(s) to datetime64[D]; raise ValueError if invalid."""
    try:
        with warnings.catch_warnings():
            # Invalid input raises below; don't also leak pandas' parse warnings.
            warnings.simplefilter("ignore")
            d = pd.to_datetime(pd.Series(np.atleast_1d(np.asarray(dates, dtype=object))))
    except (ValueError, TypeError) as e:
        raise ValueError(str(e)) from None
    if d.isna().any():
        raise ValueError("NaT")
    return d.to_numpy().astype("datetime64[D]")


def date_strings(dates) -> list[str]:
    return [str(d) for d in to_datetime64(dates)]


def to_jd(dates) -> list[int]:
    days = to_datetime64(dates).astype(np.int64)
    return [int(x) + JD_UNIX_EPOCH for x in days]


def jd_to_datetime(jd) -> np.ndarray:
    """Julian Day Numbers to datetime64[ns] (pandas' native resolution)."""
    days = np.asarray(jd, dtype=np.int64) - JD_UNIX_EPOCH
    return days.astype("datetime64[D]").astype("datetime64[ns]")


# Arguments

def pair(x, what: str) -> list[float] | None:
    """Validate an optional (min, max) pair."""
    if x is None:
        return None
    vals = [float(v) for v in np.atleast_1d(x)]
    if len(vals) != 2:
        raise ValueError(f"{what} must be a numeric vector of length 2.")
    return vals


def int_or_zero(x) -> int:
    """R's `if (is.numeric(x)) as.integer(x) else 0L`: False/None mean 0."""
    if x is None or isinstance(x, (bool, np.bool_)):
        return 0
    return int(x)


def nthreads(n_threads) -> int:
    """n_threads argument -> C++: None means the session default (0 tells the
    core to use getHW3threads())."""
    return 0 if n_threads is None else int(n_threads)


def compression_level(compression) -> int:
    """compression argument -> C++: one whole number from 0 to 9."""
    if (isinstance(compression, bool) or not isinstance(compression, (int, float))
            or compression != int(compression) or not 0 <= compression <= 9):
        raise ValueError("compression must be a whole number from 0 to 9.")
    return int(compression)


def choose(value: str, choices: tuple[str, ...], what: str) -> str:
    if value not in choices:
        raise ValueError(f"'{what}' must be one of: {', '.join(choices)}.")
    return value


def nullable_int(values, na_values=(0,)) -> pd.arrays.IntegerArray:
    """Integer column with the given sentinel values as missing (R NA_integer_)."""
    v = np.asarray(values, dtype=np.int64)
    mask = np.isin(v, na_values)
    return pd.arrays.IntegerArray(v.astype(np.int32), mask)


@contextmanager
def quiet_core(quiet: bool):
    """Silence the C++ core's messages (not warnings) on this thread while active."""
    old = _core.hw3_set_quiet(bool(quiet))
    try:
        yield
    finally:
        _core.hw3_set_quiet(old)
