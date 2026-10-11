"""Categories and summaries (ports of R/category3.R, R/category_daily3.R,
R/block_average3.R)."""

from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd

from . import _core
from ._frames import CAT_LABELS
from ._utils import choose, jd_to_datetime, nullable_int, resolve_input, to_jd
from .events import _resolve_nc_input

_SEASONS = np.array(["", "Summer", "Fall", "Winter", "Spring"], dtype=object)


def category3(event_file, clim_file=None, hemisphere="south", roundVal=4, S=None, sst_file=None,
              var_name=None):
    """Hobday et al. (2018) categories for each event.

    Categories follow ``heatwaveR::category()``. For each event day the
    category limits are the threshold and ``thresh + k * diff`` (k = 1, 2, 3),
    with ``diff = thresh - seas`` and each limit rounded to ``roundVal``
    decimals. ``p_moderate`` to ``p_extreme`` are the percentages of the
    event's duration (days, rounded to whole numbers) spent in each category,
    and the event's category is the highest one with a non-zero percentage.
    Cold-spells are recognised, as in heatwaveR, by their negative
    intensities. The season is that of the event's peak month.

    Reads the categories computed inline by ``detect_event3(category=True)``;
    for event files without them, computes them from ``clim_file`` and
    ``sst_file`` (the day counts need each event day's temperature).

    Parameters
    ----------
    event_file : str or path-like
        Events NetCDF from :func:`detect_event3`.
    clim_file : str or path-like, optional
        Climatology NetCDF; required when the event file has no categories.
    hemisphere : {"south", "north"}
        Season naming convention.
    roundVal : int
        Decimal places for the rounded outputs.
    S : bool, optional
        Deprecated; use ``hemisphere``.
    sst_file : str, path-like, or list of these, optional
        The SST used for detection (a file, a directory of daily files or a
        list of files); required when the event file has no categories.
    var_name : str, optional
        SST variable name; auto-detected when ``None``.

    Returns
    -------
    pandas.DataFrame
        One row per event: event_no, lon, lat, peak_date, category,
        intensity_max, duration, p_moderate, p_strong, p_severe, p_extreme
        (whole percentages), season.
    """
    if S is not None:
        warnings.warn("S is deprecated; use hemisphere instead.", DeprecationWarning,
                      stacklevel=2)
        hemisphere = "south" if S is True else "north"
    hemisphere = choose(hemisphere, ("south", "north"), "hemisphere")

    if event_file is None or not os.path.exists(event_file):
        raise FileNotFoundError(f"Event file does not exist: {event_file}")
    event_file = os.fspath(event_file)
    if clim_file is None:
        ev = _core.hw3_read_event_nc(event_file)
        if len(ev["category"]) == 0 or not (ev["category"] > 0).any():
            raise ValueError("clim_file and sst_file are required because the event file "
                             "does not contain pre-computed categories.")
        clim_file = ""
    elif not os.path.exists(clim_file):
        raise FileNotFoundError(f"Climatology file does not exist: {clim_file}")

    sst_files = None if sst_file is None else _resolve_nc_input(sst_file, "sst_file")
    r = _core.hw3_category(event_file, os.fspath(clim_file), hemisphere=hemisphere,
                           roundVal=int(roundVal), sst_files=sst_files,
                           var_name="" if var_name is None else var_name)
    return pd.DataFrame({
        "event_no": r["event_no"],
        "lon": r["lon"],
        "lat": r["lat"],
        "peak_date": jd_to_datetime(r["peak_jd"]),
        "category": CAT_LABELS[r["category"]],
        "intensity_max": r["intensity_max"],
        "duration": r["duration"],
        "p_moderate": r["p_moderate"],
        "p_strong": r["p_strong"],
        "p_severe": r["p_severe"],
        "p_extreme": r["p_extreme"],
        "season": _SEASONS[r["season"]],
    })


def category_daily3(sst_file, clim_file, event_file, time_range, var_name=None,
                    coldSpells=False, ice_thresh=-1.7, roundRes=2, skip_bad_files=False):
    """Per-pixel daily heatwave/cold-spell series for a date window.

    Built from SST, a climatology and an events file without re-running
    detection, so events that began before the window are still recognised.

    Parameters
    ----------
    sst_file : str, path-like, or list of these
        SST NetCDF file, directory of daily files, or list of files.
    clim_file, event_file : str or path-like
        Outputs of :func:`ts2clm3` and :func:`detect_event3`.
    time_range : (start, end) dates
        The window, e.g. ``("2026-05-10", "2026-05-20")``.
    var_name : str, optional
        SST variable name; auto-detected when ``None``.
    coldSpells : bool
        Cold-spell categories (also assigns the ice category).
    ice_thresh : float
        Threshold below which a cold-spell day is assigned category 5 (ice).
    roundRes : int
        Decimal places for ``intensity``.
    skip_bad_files : bool
        For multi-file input, skip unreadable or mismatched files.

    Returns
    -------
    pandas.DataFrame
        The full daily grid: lon, lat, t, temp, seas, thresh, intensity,
        event, event_no, category (and depth when depth-resolved).
        ``category`` is set on event-member exceedance days and missing
        elsewhere.
    """
    if not os.path.exists(clim_file):
        raise FileNotFoundError(f"Climatology file does not exist: {clim_file}")
    if not os.path.exists(event_file):
        raise FileNotFoundError(f"Events file does not exist: {event_file}")
    if time_range is None or len(time_range) != 2:
        raise ValueError("time_range must be a vector of two dates.")
    try:
        t_jd = to_jd(time_range)
    except ValueError:
        raise ValueError('time_range must be dates, e.g. ("2026-05-10", '
                         '"2026-05-20").') from None
    files, multi = resolve_input(sst_file, missing_msg="File does not exist: ",
                                 dir_what="the SST directory")
    vn = "" if var_name is None else var_name
    common = dict(t_jd_range=t_jd, var_name=vn, coldSpells=bool(coldSpells),
                  ice_thresh=float(ice_thresh), roundRes=int(roundRes))
    if multi:
        res = _core.hw3_category_daily_multi(files, os.fspath(clim_file),
                                             os.fspath(event_file),
                                             skip_bad_files=bool(skip_bad_files), **common)
    else:
        res = _core.hw3_category_daily(files[0], os.fspath(clim_file),
                                       os.fspath(event_file), **common)
    df = pd.DataFrame({
        "lon": res["lon"],
        "lat": res["lat"],
        "t": jd_to_datetime(res["jd"]),
        "temp": res["temp"],
        "seas": res["seas"],
        "thresh": res["thresh"],
        "intensity": res["intensity"],
        "event": res["event"] != 0,
        "event_no": nullable_int(res["event_no"], (_core.NA_INT,)),
        "category": nullable_int(res["category"], (_core.NA_INT,)),
    })
    if "depth" in res:
        df["depth"] = res["depth"]
    return df


def block_average3(event_file):
    """Yearly event-metric aggregation per pixel.

    Parameters
    ----------
    event_file : str or path-like
        Events NetCDF from :func:`detect_event3`.

    Returns
    -------
    pandas.DataFrame
        One row per (pixel, year of event start): lon, lat, year, count,
        duration_mean, duration_max, intensity_mean, intensity_max_mean,
        intensity_max_max, intensity_cumulative_mean, total_days, total_icum.
    """
    return pd.DataFrame(_core.hw3_block_average(os.fspath(event_file)))
