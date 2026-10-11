"""Long data frames from the C++ readers (ports of the R builders in export.R).

Gridded products use the pixel-major layout of the C++ core:
pixel = (ilon * nlat + ilat) * ndepth + idepth, then the third axis
(day-of-year or time) varies fastest within each pixel.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ._utils import jd_to_datetime, nullable_int

CAT_LABELS = np.array([None, "I Moderate", "II Strong", "III Severe", "IV Extreme"], dtype=object)
SEASON_LABELS = np.array([None, "Summer", "Fall", "Winter", "Spring"], dtype=object)

EVENT_METRICS = [
    "duration",
    "intensity_mean", "intensity_max", "intensity_var", "intensity_cumulative",
    "intensity_mean_relThresh", "intensity_max_relThresh",
    "intensity_var_relThresh", "intensity_cumulative_relThresh",
    "intensity_mean_abs", "intensity_max_abs", "intensity_var_abs",
    "intensity_cumulative_abs",
    "rate_onset", "rate_decline",
]


def _pixel_axes(nlon, nlat, ndepth, pixels):
    """(ilon, ilat, idepth) 0-based for 0-based pixel indices."""
    pixels = np.asarray(pixels)
    idepth = pixels % ndepth
    lonlat = pixels // ndepth
    return lonlat // nlat, lonlat % nlat, idepth


def clim_df(cd: dict, start: int = 0, stop: int | None = None) -> pd.DataFrame:
    """Climatology rows for pixels [start, stop): lon, lat, doy, seas, thresh, (depth)."""
    ndepth = max(1, cd["ndepth"])
    npix = cd["nlon"] * cd["nlat"] * ndepth
    stop = npix if stop is None else stop
    ndoy = cd["ndoy"]
    has_depth = ndepth > 1 and len(cd["depth"]) > 0
    pixels = np.arange(start, stop)
    ilon, ilat, idepth = _pixel_axes(cd["nlon"], cd["nlat"], ndepth, pixels)
    sl = slice(start * ndoy, stop * ndoy)
    df = pd.DataFrame({
        "lon": np.repeat(cd["lon"][ilon], ndoy),
        "lat": np.repeat(cd["lat"][ilat], ndoy),
        "doy": np.tile(np.arange(1, ndoy + 1, dtype=np.int32), len(pixels)),
        "seas": cd["seas"][sl],
        "thresh": cd["thresh"][sl],
    })
    if has_depth:
        df["depth"] = np.repeat(cd["depth"][idepth], ndoy)
    return df


def daily_df(dd: dict, start: int = 0, stop: int | None = None) -> pd.DataFrame:
    """Per-day rows (daily or protoevents) for pixels [start, stop)."""
    ndepth = max(1, dd["ndepth"])
    npix = dd["nlon"] * dd["nlat"] * ndepth
    stop = npix if stop is None else stop
    ntime = dd["ntime"]
    has_depth = ndepth > 1 and len(dd["depth"]) > 0
    pixels = np.arange(start, stop)
    ilon, ilat, idepth = _pixel_axes(dd["nlon"], dd["nlat"], ndepth, pixels)
    sl = slice(start * ntime, stop * ntime)
    has = lambda k: len(dd[k]) > 0
    df = pd.DataFrame({
        "lon": np.repeat(dd["lon"][ilon], ntime),
        "lat": np.repeat(dd["lat"][ilat], ntime),
        "t": np.tile(jd_to_datetime(dd["time_jd"]), len(pixels)),
        "temp": dd["temp"][sl],
        "seas": dd["seas"][sl],
        "thresh": dd["thresh"][sl],
    })
    if has("threshCriterion"):
        df["threshCriterion"] = dd["threshCriterion"][sl] != 0
    if has("durationCriterion"):
        df["durationCriterion"] = dd["durationCriterion"][sl] != 0
    if has("intensity"):
        df["intensity"] = dd["intensity"][sl]
    if has("event"):
        df["event"] = dd["event"][sl] != 0
    if has("event_no"):
        df["event_no"] = nullable_int(dd["event_no"][sl])
    if has("category"):
        df["category"] = nullable_int(dd["category"][sl])
    if has_depth:
        df["depth"] = np.repeat(dd["depth"][idepth], ntime)
    return df


def event_df(ev: dict, idx=None) -> pd.DataFrame:
    """Event rows (all, or the 0-based indices `idx`), as .event_df_from_indices()."""
    idx = np.arange(ev["nevents"]) if idx is None else np.asarray(idx, dtype=np.int64)
    ref = ev["ref_date_jd"]
    df = pd.DataFrame({
        "lon": ev["lon"][idx],
        "lat": ev["lat"][idx],
        "pixel_index": ev["pixel_index"][idx],
        "event_no": ev["event_no"][idx],
        "date_start": jd_to_datetime(ref + ev["date_start"][idx]),
        "date_peak": jd_to_datetime(ref + ev["date_peak"][idx]),
        "date_end": jd_to_datetime(ref + ev["date_end"][idx]),
    })
    for m in EVENT_METRICS:
        df[m] = ev[m][idx]
    if len(ev["depth"]) > 0:
        df["depth"] = ev["depth"][idx]
    if len(ev["category"]) > 0 and (ev["category"] > 0).any():
        df["category"] = CAT_LABELS[np.clip(ev["category"][idx], 0, 4)]
        for p in ("p_moderate", "p_strong", "p_severe", "p_extreme"):
            df[p] = ev[p][idx]
        df["season"] = SEASON_LABELS[np.clip(ev["season"][idx], 0, 4)]
    return df


def subset_df(sub: dict) -> pd.DataFrame:
    """Long frame from a hyperslab read: lon, lat, t|doy, vars..., (depth)."""
    nlon, nlat, n3 = sub["nlon"], sub["nlat"], sub["n3"]
    has_depth = bool(sub["has_depth"])
    ndepth = sub["ndepth"] if has_depth else 1
    vnames = list(sub["data"].keys())
    third_is_t = sub["third_name"] == "t"
    third = "t" if third_is_t else "doy"

    if nlon == 0 or nlat == 0 or (has_depth and ndepth == 0) or n3 == 0:
        cols = {"lon": np.array([], float), "lat": np.array([], float)}
        if has_depth:
            cols["depth"] = np.array([], float)
        cols[third] = (np.array([], "datetime64[ns]") if third_is_t
                       else np.array([], np.int32))
        for v in vnames:
            cols[v] = np.array([], float)
        return _finalize_subset(pd.DataFrame(cols), vnames)

    npix = nlon * nlat * ndepth
    df = pd.DataFrame({
        "lon": np.repeat(sub["lon"], nlat * ndepth * n3),
        "lat": np.tile(np.repeat(sub["lat"], ndepth * n3), nlon),
    })
    if third_is_t:
        df["t"] = np.tile(jd_to_datetime(sub["third"]), npix)
    else:
        df["doy"] = np.tile(sub["third"].astype(np.int32), npix)
    for v in vnames:
        df[v] = sub["data"][v]
    if has_depth:
        df["depth"] = np.tile(np.repeat(sub["depth"], n3), nlon * nlat)
    return _finalize_subset(df, vnames)


def _finalize_subset(df: pd.DataFrame, vnames: list[str]) -> pd.DataFrame:
    for v in ("threshCriterion", "durationCriterion", "event"):
        if v in vnames:
            df[v] = df[v].to_numpy() != 0
    for v in ("event_no", "category"):
        if v in vnames:
            df[v] = nullable_int(df[v].to_numpy())
    return df
