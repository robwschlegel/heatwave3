"""Event detection (ports of R/detect_event3.R, R/detect3.R, R/exceedance3.R,
R/detect_event_array3.R)."""

from __future__ import annotations

import os
import tempfile

import numpy as np

from . import _core
from ._utils import (check_out_dir, choose, clim_path, compression_level, daily_path, date_strings,
                     events_path, is_multi, list_nc, nthreads, pair, proto_path, quiet_core,
                     resolve_input)
from .climatology import ts2clm3
from .export import console_summary


def detect_event3(file_in, name, clim_file=None, var_name=None, minDuration=5,
                  minDuration2=None, joinAcrossGaps=True, maxGap=2, maxGap2=None,
                  threshClim2=None, threshClim2_var_name=None, coldSpells=False,
                  category=False, hemisphere="south", roundRes=4, daily="none",
                  slim=False, protoEvent=False, n_threads=None, compression=1,
                  skip_bad_files=False, quiet=False):
    """Detect marine heatwaves (or cold-spells) in gridded data.

    The gridded equivalent of ``heatwaveR::detect_event()``, using a
    climatology from :func:`ts2clm3`.

    Parameters
    ----------
    file_in : str, path-like, or list of these
        SST NetCDF file, directory of daily files, or list of files.
    name : str or path-like
        Output path stem. Events go to ``<name>_events.nc``; per-day products
        to ``<name>_events_daily.nc`` and ``<name>_protoevents.nc``.
    clim_file : str, optional
        Climatology from :func:`ts2clm3`; defaults to ``<name>_clim.nc``.
    var_name : str, optional
        SST variable name; auto-detected when ``None``.
    minDuration : int
        Minimum event duration (days).
    minDuration2 : int, optional
        Minimum duration for the ``threshClim2`` pass; defaults to ``minDuration``.
    joinAcrossGaps : bool
        Join events separated by short gaps?
    maxGap : int
        Maximum gap (days) to join across.
    maxGap2 : int, optional
        Maximum gap for the ``threshClim2`` pass; defaults to ``maxGap``.
    threshClim2 : str, path-like, or list of these, optional
        Gridded logical criterion for a secondary pass (non-zero = True),
        aligned with ``file_in``.
    threshClim2_var_name : str, optional
        Variable name in ``threshClim2``; auto-detected when ``None``.
    coldSpells : bool
        Detect cold-spells instead of heatwaves.
    category : bool
        Compute Hobday et al. (2018) categories inline.
    hemisphere : {"south", "north"}
        Season naming (used when ``category=True``).
    roundRes : int
        Decimal places for event metrics.
    daily : {"none", "also", "only"}
        Also (or only) write the per-day series: temp, seas, thresh,
        intensity, event, event_no, and the daily Hobday (2018) category on
        event days above the threshold (missing on all other days).
    slim : bool
        With ``daily`` = ``"also"`` or ``"only"``, write only ``intensity`` and
        ``category`` to the per-day file (no temp, seas, thresh, event or
        event_no): about half the size, and quicker to write. The
        dropped variables can be rebuilt from the SST and climatology files.
    protoEvent : bool
        Write the per-day proto-event series instead of the event table.
        Mutually exclusive with ``daily``.
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
    dict
        Paths of the products written, keyed ``"events"``, ``"daily"`` and/or
        ``"protoevents"`` (the events file is skipped when there are no events).
        Use :func:`hw3_export` to read them.
    """
    minDuration2 = minDuration if minDuration2 is None else minDuration2
    maxGap2 = maxGap if maxGap2 is None else maxGap2
    daily = choose(daily, ("none", "also", "only"), "daily")
    if slim and daily == "none":
        raise ValueError("'slim' applies to the per-day file: set daily = \"also\" or \"only\".")
    if protoEvent and daily != "none":
        raise ValueError("Set only one of 'protoEvent' or 'daily' (\"also\"/\"only\"); "
                         "they are mutually exclusive.")

    if clim_file is None:
        clim_file = clim_path(name)
    clim_file = os.fspath(clim_file)
    if not os.path.exists(clim_file):
        raise FileNotFoundError(
            f"Climatology file does not exist: {clim_file}. It defaults to "
            "<name>_clim.nc; pass 'clim_file' to override.")
    check_out_dir(events_path(name))

    # A non-empty path tells the C++ layer to write that product; "" skips it.
    ev_file = "" if (protoEvent or daily == "only") else events_path(name)
    d_file = daily_path(name) if daily in ("also", "only") else ""
    p_file = proto_path(name) if protoEvent else ""

    south = choose(hemisphere, ("south", "north"), "hemisphere") == "south"
    vn = "" if var_name is None else var_name
    t2_vn = "" if threshClim2_var_name is None else threshClim2_var_name

    multi = is_multi(file_in) or os.path.isdir(file_in)
    if multi:
        files, _ = resolve_input(file_in)
    thresh2 = None
    if threshClim2 is not None:
        thresh2 = _resolve_nc_input(threshClim2, "threshClim2")
        if not multi and len(thresh2) > 1:
            raise ValueError("threshClim2 may only contain multiple files when file_in is "
                             "also a directory or vector of files.")

    args = dict(clim_file=clim_file, events_file=ev_file, var_name=vn,
                minDuration=int(minDuration), minDuration2=int(minDuration2),
                joinAcrossGaps=bool(joinAcrossGaps), maxGap=int(maxGap),
                maxGap2=int(maxGap2), coldSpells=bool(coldSpells), roundRes=int(roundRes),
                n_threads=nthreads(n_threads), compression=compression_level(compression),
                category=bool(category), southHemisphere=south,
                threshClim2_var_name=t2_vn, daily_file=d_file, proto_file=p_file,
                slim=bool(slim))
    if not multi and not os.path.exists(os.fspath(file_in)):
        raise FileNotFoundError(f"SST file does not exist: {os.fspath(file_in)}")
    with quiet_core(quiet):
        if multi:
            _core.hw3_detect_events_multi(files, threshClim2_files=thresh2,
                                          skip_bad_files=bool(skip_bad_files), **args)
        else:
            _core.hw3_detect_events(os.fspath(file_in),
                                    threshClim2_file="" if thresh2 is None else thresh2[0],
                                    **args)

    # Report the products actually written (no events file when no events).
    written = {k: v for k, v in (("events", ev_file), ("daily", d_file),
                                 ("protoevents", p_file)) if v and os.path.exists(v)}
    for f in written.values():
        console_summary(f, quiet=quiet)
    return written


def detect3(file_in, name, climatologyPeriod, var_name=None, lon_range=None,
            lat_range=None, time_range=None, depth=None, depth_range=None,
            maxPadLength=False, windowHalfWidth=5, pctile=90, smoothPercentile=True,
            smoothPercentileWidth=31, detrend=False, minDuration=5, minDuration2=None,
            joinAcrossGaps=True, maxGap=2, maxGap2=None, threshClim2=None,
            threshClim2_var_name=None, coldSpells=False, category=False,
            hemisphere="south", roundClm=4, roundRes=4, daily="none", slim=False,
            protoEvent=False, n_threads=None, compression=1, skip_bad_files=False,
            quiet=False):
    """Climatology and event detection in one call (the primary entry point).

    Runs :func:`ts2clm3` then :func:`detect_event3` with the same ``name``
    stem; see those functions for the parameters.

    Returns
    -------
    dict
        Paths written: ``"clim"`` plus those returned by :func:`detect_event3`.
    """
    daily = choose(daily, ("none", "also", "only"), "daily")
    if slim and daily == "none":
        raise ValueError("'slim' applies to the per-day file: set daily = \"also\" or \"only\".")
    clim_file = ts2clm3(file_in, name, climatologyPeriod, lon_range=lon_range,
                        lat_range=lat_range, time_range=time_range, depth=depth,
                        depth_range=depth_range, var_name=var_name,
                        maxPadLength=maxPadLength, windowHalfWidth=windowHalfWidth,
                        pctile=pctile, smoothPercentile=smoothPercentile,
                        smoothPercentileWidth=smoothPercentileWidth, detrend=detrend,
                        roundClm=roundClm, n_threads=n_threads, compression=compression,
                        skip_bad_files=skip_bad_files, quiet=quiet)
    written = detect_event3(file_in, name, clim_file=clim_file, var_name=var_name,
                            minDuration=minDuration, minDuration2=minDuration2,
                            joinAcrossGaps=joinAcrossGaps, maxGap=maxGap, maxGap2=maxGap2,
                            threshClim2=threshClim2,
                            threshClim2_var_name=threshClim2_var_name,
                            coldSpells=coldSpells, category=category,
                            hemisphere=hemisphere, roundRes=roundRes, daily=daily,
                            slim=slim, protoEvent=protoEvent, n_threads=n_threads,
                            compression=compression,
                            skip_bad_files=skip_bad_files, quiet=quiet)
    return {"clim": clim_file, **written}


def exceedance3(file_in, file_out, threshold, var_name=None, below=False, minDuration=5,
                joinAcrossGaps=True, maxGap=2, lon_range=None, lat_range=None,
                time_range=None, depth=None, roundRes=4, n_threads=None, compression=1,
                quiet=False):
    """Detect exceedances of a static threshold in gridded data.

    Parameters
    ----------
    file_in : str or path-like
        SST NetCDF file.
    file_out : str or path-like
        Path of the events NetCDF to write.
    threshold : float
        Static threshold (same units as the data).
    below : bool
        Detect runs below the threshold instead of above.
    minDuration, joinAcrossGaps, maxGap, roundRes, n_threads, compression
        As in :func:`detect_event3`.
    var_name, lon_range, lat_range, time_range, depth
        As in :func:`ts2clm3`.
    quiet : bool
        Suppress the progress messages from the computation. Warnings are
        still shown.

    Returns
    -------
    str
        ``file_out``.
    """
    if threshold is None:
        raise ValueError("A threshold value must be provided.")
    level = compression_level(compression)
    file_in, file_out = os.fspath(file_in), os.fspath(file_out)
    vn = "" if var_name is None else var_name
    fd, clim_file = tempfile.mkstemp(suffix=".nc")
    os.close(fd)
    try:
        with quiet_core(quiet):
            # A constant-value climatology (seas = thresh = threshold) reuses
            # the event detection machinery.
            _core.hw3_write_const_clim(
                file_in, clim_file, vn, float(threshold),
                lon_range=pair(lon_range, "lon_range"), lat_range=pair(lat_range, "lat_range"),
                time_range=None if time_range is None else date_strings(time_range),
                depth=-1 if depth is None else int(depth))
            _core.hw3_detect_events(
                file_in, clim_file, file_out, vn, minDuration=int(minDuration),
                minDuration2=int(minDuration), joinAcrossGaps=bool(joinAcrossGaps),
                maxGap=int(maxGap), maxGap2=int(maxGap), coldSpells=bool(below),
                roundRes=int(roundRes), n_threads=nthreads(n_threads),
                compression=level, category=False,
                southHemisphere=True)
    finally:
        os.remove(clim_file)
    return file_out


ARRAY_OUTPUTS = ("event", "event_no", "intensity", "category")


def detect_event_array3(temp, seas, thresh, minDuration=5, joinAcrossGaps=True, maxGap=2,
                        coldSpells=False, roundRes=4,
                        outputs=("event", "intensity", "category"), n_threads=None):
    """Detect marine heatwaves (or cold-spells) day by day in arrays.

    The in-memory counterpart of :func:`detect_event3` with ``daily="only"``:
    the same detection, but on arrays instead of files, and with the
    climatology already matched to each day. It returns only per-day results,
    because event metrics such as duration or cumulative intensity need the
    whole event, which a short series (a forecast window, for example) may
    cut off.

    Parameters
    ----------
    temp : array_like
        Daily temperature, time on the **last** axis, any number of leading
        axes (pixels). The days must be consecutive, with no gaps. Missing
        days are NaN.
    seas, thresh : array_like
        Seasonal climatology and threshold for each day of ``temp``, same
        shape. Match them to the dates beforehand, for example from a
        :func:`ts2clm3` climatology by day of year.
    minDuration : int
        Minimum event duration (days).
    joinAcrossGaps : bool
        Join events separated by short gaps?
    maxGap : int
        Maximum gap (days) to join across.
    coldSpells : bool
        Detect cold-spells instead of heatwaves.
    roundRes : int
        Decimal places for ``intensity``.
    outputs : str or sequence of str
        Which results to return, in this order: any of ``"event"``,
        ``"event_no"``, ``"intensity"`` and ``"category"``.
    n_threads : int, optional
        Threads for the per-pixel computation; ``None`` uses the session
        default (see :func:`setHW3threads`). Use 1 when the calls themselves
        run in parallel, as with dask.

    Returns
    -------
    numpy.ndarray or tuple of numpy.ndarray
        One array per name in ``outputs`` (a single array when only one is
        requested), each with the shape of ``temp``:

        - ``event`` (bool): the day belongs to an event.
        - ``event_no`` (int32): the event's number within its pixel, counted
          from the start of this call (0 outside events). It is not an
          identifier across calls: the same event in two overlapping series
          can have different numbers.
        - ``intensity`` (float32): ``temp - seas`` on every day with data,
          rounded to ``roundRes``; NaN where either is missing.
        - ``category`` (int8): the daily Hobday et al. (2018) category, 1 (I
          Moderate) to 4 (IV Extreme), on event days beyond the threshold;
          0 on all other days.

        A pixel with no temperature or no ``seas`` at all gets ``False``,
        0, NaN and 0. Results are identical to the per-day file of
        :func:`detect_event3` for the same days and climatology.

    Notes
    -----
    Whether a day belongs to an event depends on the days around it, so
    events are cut off at both ends of the series: a run of exceedances that
    started before the first day is shorter here, and one that continues
    past the last day may not reach ``minDuration``. Prepend earlier days
    when that matters. ``intensity`` does not depend on the other days, and
    ``category`` only through event membership.

    With xarray, time as the core dimension:

    >>> xarray.apply_ufunc(detect_event_array3, temp, seas, thresh,
    ...                    input_core_dims=[["time"]] * 3,
    ...                    output_core_dims=[["time"]] * 3,
    ...                    kwargs={"n_threads": 1},
    ...                    dask="parallelized",
    ...                    output_dtypes=[bool, np.float32, np.int8])  # doctest: +SKIP
    """
    if isinstance(outputs, str):
        outputs = (outputs,)
    outputs = tuple(outputs)
    if not outputs or len(set(outputs)) != len(outputs) or not set(outputs) <= set(ARRAY_OUTPUTS):
        raise ValueError(f"'outputs' must name one or more of {', '.join(ARRAY_OUTPUTS)}, "
                         f"each at most once; got {outputs!r}.")
    result = _core.hw3_detect_events_array(
        np.asarray(temp), np.asarray(seas), np.asarray(thresh),
        minDuration=int(minDuration), joinAcrossGaps=bool(joinAcrossGaps), maxGap=int(maxGap),
        coldSpells=bool(coldSpells), roundRes=int(roundRes), n_threads=nthreads(n_threads),
        **{name: name in outputs for name in ARRAY_OUTPUTS})
    arrays = tuple(result[name] for name in outputs)
    return arrays[0] if len(arrays) == 1 else arrays


def _resolve_nc_input(path, arg_name: str) -> list[str]:
    if not is_multi(path) and os.path.isdir(path):
        return list_nc(path, f"{arg_name} directory")
    files = [os.fspath(path)] if not is_multi(path) else [os.fspath(p) for p in path]
    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"{arg_name} file does not exist: {f}")
    return files
