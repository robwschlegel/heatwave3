"""Gridded climatology (port of R/ts2clm3.R)."""

from __future__ import annotations

from . import _core
from ._utils import (check_out_dir, clim_path, compression_level, date_strings, int_or_zero, nthreads, pair,
                     quiet_core, resolve_input)
from .export import console_summary


def ts2clm3(file_in, name, climatologyPeriod, lon_range=None, lat_range=None,
            time_range=None, depth=None, depth_range=None, var_name=None,
            maxPadLength=False, windowHalfWidth=5, pctile=90, smoothPercentile=True,
            smoothPercentileWidth=31, var=False, detrend=False, roundClm=4,
            n_threads=None, compression=1, skip_bad_files=False, quiet=False):
    """Compute a gridded climatology and write it to ``<name>_clim.nc``.

    The gridded equivalent of ``heatwaveR::ts2clm()``: for every pixel, the
    seasonal climatology (``seas``) and percentile threshold (``thresh``) for
    each day-of-year (1 to 366).

    Parameters
    ----------
    file_in : str, path-like, or list of these
        A NetCDF file of daily SST, a directory of daily ``.nc``/``.nc4``
        files, or a list of files.
    name : str or path-like
        Output path stem; the climatology is written to ``<name>_clim.nc``.
        The directory must exist.
    climatologyPeriod : (start, end) dates
        Base period, e.g. ``("1982-01-01", "2011-12-31")``.
    lon_range, lat_range : (min, max), optional
        Spatial subset.
    time_range : (start, end) dates, optional
        Temporal subset (single-file input only).
    depth : int, optional
        Legacy single depth-level index (0-based); the depth axis is dropped.
    depth_range : (min, max), optional
        Keep the depth levels in this range (metres) as a real dimension.
        Mutually exclusive with ``depth``.
    var_name : str, optional
        SST variable name; auto-detected when ``None``.
    maxPadLength : int or False
        Maximum run of missing days to interpolate; ``False`` means 0.
    windowHalfWidth : int
        Half-width (days) of the window pooled around each day-of-year.
    pctile : float
        Threshold percentile (90 for heatwaves, 10 for cold-spells).
    smoothPercentile : bool
        Smooth seas and thresh with a circular moving average?
    smoothPercentileWidth : int
        Width (days) of that moving average.
    var : bool
        Also compute the climatological variance.
    detrend : bool
        Remove a linear trend before computing the climatology.
    roundClm : int or False
        Decimal places for rounding; ``False`` means no rounding.
    n_threads : int, optional
        Threads for the per-pixel computation. ``None`` (the default) uses
        the session default: half the available cores, unless changed with
        :func:`setHW3threads` or the ``HEATWAVE3_NUM_THREADS`` environment
        variable.
    compression : int
        zlib compression level of the NetCDF output, 0 to 9: 1 (the default)
        compresses fastest, 9 gives the smallest files (a few percent
        smaller), and 0 writes uncompressed files, which are several times
        larger but much quicker to write. Gridded outputs are compressed on
        ``n_threads`` threads when heatwave3 was built with HDF5.
    skip_bad_files : bool
        For multi-file input, skip unreadable or mismatched files.
    quiet : bool
        Suppress all console output: the progress messages from the
        computation and the post-computation summary. Warnings, such as
        skipped files, are still shown.

    Returns
    -------
    str
        Path of the written climatology file.
    """
    if climatologyPeriod is None or len(climatologyPeriod) != 2:
        raise ValueError("climatologyPeriod must be a character vector of length 2.")
    if depth is not None and depth_range is not None:
        raise ValueError("Specify only one of 'depth' or 'depth_range', not both.")
    depth_range = pair(depth_range, "depth_range")

    file_out = clim_path(name)
    check_out_dir(file_out)

    files, multi = resolve_input(file_in)
    common = dict(
        file_out=file_out,
        var_name="" if var_name is None else var_name,
        climatologyPeriod=date_strings(climatologyPeriod),
        lon_range=pair(lon_range, "lon_range"),
        lat_range=pair(lat_range, "lat_range"),
        depth=-1 if depth is None else int(depth),
        depth_range=depth_range,
        maxPadLength=int_or_zero(maxPadLength),
        windowHalfWidth=int(windowHalfWidth),
        pctile=float(pctile),
        smoothPercentile=bool(smoothPercentile),
        smoothPercentileWidth=int(smoothPercentileWidth),
        compute_var=bool(var),
        roundClm=int_or_zero(roundClm),
        n_threads=nthreads(n_threads),
        compression=compression_level(compression),
        detrend=bool(detrend),
    )
    with quiet_core(quiet):
        if multi:
            _core.hw3_compute_clim_multi(files, skip_bad_files=bool(skip_bad_files), **common)
        else:
            time_range = None if time_range is None else date_strings(time_range)
            _core.hw3_compute_clim(files[0], time_range=time_range, **common)

    console_summary(file_out, quiet=quiet)
    return file_out
