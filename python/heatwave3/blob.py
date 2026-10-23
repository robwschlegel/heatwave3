"""Spatially connected heatwave blobs (port of R/detect_blob3.R)."""

from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd

from . import _core
from ._utils import jd_to_datetime, list_nc

_RETURNS = ("event", "daily", "voxel")


def detect_blob3(sst_file, clim_file, var_name=None, wrapDateline=False, minVoxels=1,
                 topN=None, rankBy="cumI_km2_day", coldSpells=False,
                 return_=("event", "daily"), skip_bad_files=False, connectivity=6):
    """Detect spatially connected marine heatwave (or cold-spell) blobs.

    Builds the exceedance mask (``temp > thresh``, or ``temp < thresh`` for
    cold-spells) from the SST and the climatology, then labels connected
    voxels across longitude, latitude, (optionally) depth and time with a
    union-find algorithm. Neighbours differ by at most one grid step along
    each axis, and ``connectivity`` sets which of them in the lon/lat/time
    cube count. With the default (6, faces only) a voxel has no horizontal
    or space-time diagonal neighbours, so a blob continues from one day to
    the next only through pixels that exceed on both days. When
    ``clim_file`` is depth-resolved (from ``ts2clm3(depth_range=...)``), a
    voxel also joins the same pixel on the same day one level above and
    below (never diagonally in depth), and volume and depth metrics are
    added.

    Parameters
    ----------
    sst_file : str, path-like, or list of these
        SST NetCDF file (as passed to :func:`ts2clm3`), directory of daily
        files, or list of files.
    clim_file : str or path-like
        Climatology from :func:`ts2clm3`; its grid sets the region analysed.
    var_name : str, optional
        SST variable name; auto-detected when ``None``.
    wrapDateline : bool
        Let blobs connect across the antimeridian (longitude only).
    minVoxels : int
        Minimum number of voxels for a blob to be kept.
    topN : int, optional
        Keep only the top ``topN`` blobs ranked by ``rankBy``.
    rankBy : str
        Ranking metric: ``"cumI_km2_day"`` (default), ``"peakArea_km2"``,
        ``"duration"``, ``"meanArea_km2"``, and for depth-resolved input
        ``"totalVolume_km3"``, ``"peakVolume_km3"``, ``"meanVolume_km3"``.
    coldSpells : bool
        Detect cold-spell blobs (needs a low-percentile climatology, e.g.
        ``pctile=10``).
    return_ : sequence of {"event", "daily", "voxel"}
        Components to return. Named ``return`` in R; ``return`` is a Python
        keyword.
    skip_bad_files : bool
        For multi-file SST input, skip unreadable or mismatched files.
    connectivity : {6, 10, 18, 26}
        Which neighbours in the 3 x 3 x 3 lon/lat/time cube join a blob:
        6 (default) faces only; 10 adds the horizontal diagonals on the same
        day; 18 adds all edges, including space-time diagonals (a pixel joins
        its neighbour on the next day even if the two never exceed on the
        same day); 26 adds the corners. With depth the cube applies on each
        level and depth levels join only straight up and down, so 8, 12, 20
        or 28 neighbours in all. Higher settings only merge blobs. Adjacency
        is by grid index, so a space-time diagonal stands for a displacement
        of one grid cell per day.

    Returns
    -------
    dict of pandas.DataFrame
        The requested components: ``"event"`` (one row per blob),
        ``"daily"`` (one row per blob and day) and ``"voxel"`` (one row per
        voxel, the full footprint).

    Notes
    -----
    Blobs with exactly equal values of ``rankBy`` (common for small blobs)
    are ranked in blob order. On CPUs whose compilers fuse multiply-adds by
    default (e.g. ARM, including Apple Silicon) the accumulated metrics can
    differ from x86-64 in the last bit, so tied blobs can be ordered and
    ranked differently, and ``topN`` can select a different blob at a tie.
    Values agree to about 1e-15 relative.
    """
    if sst_file is None or clim_file is None:
        raise ValueError("sst_file and clim_file must both be provided.")
    clim_file = os.fspath(clim_file)
    if not os.path.exists(clim_file):
        raise FileNotFoundError(f"Climatology file does not exist: {clim_file}")

    ret = [return_] if isinstance(return_, str) else list(return_)
    bad = [r for r in ret if r not in _RETURNS]
    if bad:
        raise ValueError(f"return contains unsupported value: {bad[0]}")
    ret = list(dict.fromkeys(ret))

    if connectivity not in (6, 10, 18, 26) or isinstance(connectivity, bool):
        raise ValueError("connectivity must be one of 6, 10, 18 or 26.")

    vn = "" if var_name is None else var_name
    cd = _core.hw3_read_clim_nc(clim_file)
    if len(cd["lon"]) < 2 or len(cd["lat"]) < 2:
        raise ValueError("detect_blob3 requires at least two longitude and two latitude "
                         "cells to compute cell area.")
    ndepth = max(1, cd["ndepth"])
    has_depth = ndepth > 1 and len(cd["depth"]) > 0
    depth_range = [float(cd["depth"].min()), float(cd["depth"].max())] if has_depth else None

    valid_rank = ["cumI_km2_day", "peakArea_km2", "duration", "meanArea_km2"]
    if has_depth:
        valid_rank += ["totalVolume_km3", "peakVolume_km3", "meanVolume_km3"]
    if rankBy not in valid_rank:
        raise ValueError(f"rankBy must be one of: {', '.join(valid_rank)}")

    sst_files, multi = _resolve_sst(sst_file)
    extent = dict(lon_range=[float(cd["lon"].min()), float(cd["lon"].max())],
                  lat_range=[float(cd["lat"].min()), float(cd["lat"].max())],
                  depth_range=depth_range)
    if multi:
        gd = _core.hw3_read_sst_multi(sst_files, vn, skip_bad_files=bool(skip_bad_files), **extent)
    else:
        gd = _core.hw3_read_sst(sst_files[0], vn, **extent)

    nlon, nlat, ntime = gd["nlon"], gd["nlat"], gd["ntime"]
    if nlon != cd["nlon"] or nlat != cd["nlat"] or gd["ndepth"] != cd["ndepth"]:
        raise ValueError(f"Grid mismatch: SST is {nlon} x {nlat} x {gd['ndepth']} depth but "
                         f"climatology is {cd['nlon']} x {cd['nlat']} x {cd['ndepth']} depth")

    # Day-of-year (1..366, heatwaveR convention) for each time step
    doy = np.array([_core.hw3_jd_to_doy(int(jd)) for jd in gd["time_days"]], dtype=np.int64)

    # Exceedance mask [lon, lat, depth, time]. The SST and climatology are
    # pixel-major with pixel = (ilon * nlat + ilat) * ndepth + idepth, so a
    # C-order reshape gives [lon, lat, depth, time|doy]. Comparisons with NaN
    # are False, matching R's skip of NA temp or thresh.
    sst = gd["sst"].reshape(nlon, nlat, ndepth, ntime)
    seas = cd["seas"].reshape(nlon, nlat, ndepth, 366)
    thresh = cd["thresh"].reshape(nlon, nlat, ndepth, 366)[..., doy - 1]
    exceed = sst < thresh if coldSpells else sst > thresh
    mask = np.asfortranarray(exceed.astype(np.int32))

    # Connected components with the lon/lat/time neighbours set by connectivity
    # (6: faces only); depth levels join only straight up and down
    labels, n_components = _core.label_components(mask, nlon, nlat, ndepth, ntime,
                                                  bool(wrapDateline), int(connectivity))
    if n_components == 0:
        return {r: pd.DataFrame() for r in _RETURNS if r in ret}

    # Labelled voxels in R's column-major order, as 1-based indices
    lab_flat = labels.ravel(order="F")
    w0 = np.flatnonzero(lab_flat > 0)
    lab_vec = lab_flat[w0].astype(np.int32)
    ix = (w0 % nlon + 1).astype(np.int32)
    iy = ((w0 // nlon) % nlat + 1).astype(np.int32)
    iz = ((w0 // (nlon * nlat)) % ndepth + 1).astype(np.int32)
    it = (w0 // (nlon * nlat * ndepth) + 1).astype(np.int32)

    # Anomaly (delta), sign-flipped for cold-spells
    delta = sst[ix - 1, iy - 1, iz - 1, it - 1] - seas[ix - 1, iy - 1, iz - 1, doy[it - 1] - 1]
    if coldSpells:
        delta = -delta

    # Cell area by latitude (km^2); math.cos is the C library cos, as in R
    dlon = abs(cd["lon"][1] - cd["lon"][0])
    cell_height_km = 111.32 * abs(cd["lat"][1] - cd["lat"][0])
    cell_area_lat = np.array([111.32 * math.cos(la * math.pi / 180) * dlon * cell_height_km
                              for la in cd["lat"]])

    # Vertical layer thickness (km) from depth cell centres: interior
    # boundaries are midpoints, top/bottom mirror the nearest interior one.
    if has_depth:
        d = cd["depth"]
        mid = (d[:-1] + d[1:]) / 2
        bounds = np.concatenate([[2 * d[0] - mid[0]], mid, [2 * d[-1] - mid[-1]]])
        layer_thickness_km = np.diff(bounds) / 1000
    else:
        layer_thickness_km = np.array([0.0])

    daily = pd.DataFrame(_core.blob_daily_summary(
        lab_vec, ix, iy, it, delta, cell_area_lat, cd["lon"], cd["lat"], iz,
        cd["depth"] if has_depth else np.array([], dtype=float), layer_thickness_km,
        has_depth))
    dates = jd_to_datetime(gd["time_days"])
    daily["date"] = dates[daily["t_idx"].to_numpy() - 1]

    event = _event_summary(daily, has_depth)

    voxel = None
    if "voxel" in ret:
        voxel = pd.DataFrame({
            "blob": lab_vec,
            "lon": cd["lon"][ix - 1],
            "lat": cd["lat"][iy - 1],
            "t_idx": it,
            "date": dates[it - 1],
            "delta": delta,
        })
        if has_depth:
            voxel["depth"] = cd["depth"][iz - 1]

    # Filter by minVoxels
    if minVoxels > 1:
        blobs, counts = np.unique(lab_vec, return_counts=True)
        keep = blobs[counts >= minVoxels]
        event, daily = event[event["blob"].isin(keep)], daily[daily["blob"].isin(keep)]
        if voxel is not None:
            voxel = voxel[voxel["blob"].isin(keep)]

    # Rank (stable, as R's order()) and keep the top N
    if len(event):
        event = event.sort_values(rankBy, ascending=False, kind="stable")
        event["rank"] = np.arange(1, len(event) + 1, dtype=np.int32)
        if topN is not None and topN < len(event):
            event = event.iloc[:int(topN)]
            keep = event["blob"].to_numpy()
            daily = daily[daily["blob"].isin(keep)]
            if voxel is not None:
                voxel = voxel[voxel["blob"].isin(keep)]

    out = {"event": event, "daily": daily, "voxel": voxel}
    return {r: out[r].reset_index(drop=True) for r in _RETURNS if r in ret}


def _event_summary(daily: pd.DataFrame, has_depth: bool) -> pd.DataFrame:
    """Per-blob event metrics from the per-(blob, day) table (daily is sorted
    by blob then day, so each blob's rows are in date order)."""
    g = daily.groupby("blob", sort=True)
    first_max = lambda col: daily.loc[g[col].idxmax().to_numpy()]
    peak = first_max("area_km2")
    w = daily["area_km2"]
    ev = pd.DataFrame({
        "blob": g.size().index.to_numpy(dtype=np.int32),
        "duration": g.size().to_numpy(dtype=np.int32),
        "date_start": g["date"].min().to_numpy(),
        "date_end": g["date"].max().to_numpy(),
        "date_peak": peak["date"].to_numpy(),
        "cumI_km2_day": (daily["mean_delta"] * w).groupby(daily["blob"]).sum().to_numpy(),
        "peakArea_km2": g["area_km2"].max().to_numpy(),
        "meanArea_km2": g["area_km2"].mean().to_numpy(),
        "centroid_lon": ((daily["centroid_lon"] * w).groupby(daily["blob"]).sum()
                         / g["area_km2"].sum()).to_numpy(),
        "centroid_lat": ((daily["centroid_lat"] * w).groupby(daily["blob"]).sum()
                         / g["area_km2"].sum()).to_numpy(),
    })
    if has_depth:
        ev["depth_min_m"] = g["depth_min_m"].min().to_numpy()
        ev["depth_max_m"] = g["depth_max_m"].max().to_numpy()
        ev["totalVolume_km3"] = g["volume_km3"].sum().to_numpy()
        ev["peakVolume_km3"] = g["volume_km3"].max().to_numpy()
        ev["meanVolume_km3"] = g["volume_km3"].mean().to_numpy()
        ev["depthOfPeakIntensity_m"] = first_max("max_delta")["depth_at_max_delta_m"].to_numpy()
    return ev


def _resolve_sst(sst_file) -> tuple[list[str], bool]:
    """As R's .resolve_nc_input(sst_file, "sst_file") plus the multi-file test."""
    if isinstance(sst_file, (str, os.PathLike)):
        path = os.fspath(sst_file)
        if os.path.isdir(path):
            return list_nc(path, "sst_file directory"), True
        files = [path]
    else:
        files = [os.fspath(f) for f in sst_file]
    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"sst_file file does not exist: {f}")
    return files, len(files) > 1
