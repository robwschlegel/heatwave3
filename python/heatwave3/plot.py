"""Plotting with matplotlib (ports of R/geom_flame3.R, R/geom_lolli3.R,
R/event_line3.R and R/plot_spatial3.R).

R's geom_flame3() and geom_lolli3() are ggplot2 layers; here they are
functions that draw onto a matplotlib Axes. event_line3() and plot_metric3()
keep the R arguments, add a Python-only ``ax``, and return the Axes.
"""

from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd

from . import _core
from ._utils import jd_to_datetime, list_nc, to_datetime64
from ._zarr import clim_input, local_dir, nc_input

# R colours used by the R package
_SALMON = "#FA8072"       # "salmon"
_STEELBLUE3 = "#4F94CD"   # "steelblue3"
_GREY35 = "#595959"       # "grey35"
_GREY40 = "#666666"       # "grey40"
_DARKGREEN = "#006400"    # "darkgreen"


def _plt():
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError("Plotting needs matplotlib: pip install matplotlib "
                          "(or pip install heatwave3[plot])") from None
    return plt


def _column(v, data, like=None):
    """A column of `data` by name, or the values themselves; scalars are
    repeated to the length of `like` (as ggplot2 recycles a constant)."""
    out = np.asarray(data[v]) if (data is not None and isinstance(v, str)) else np.asarray(v)
    return np.broadcast_to(out, np.shape(like)) if (like is not None and out.ndim == 0) else out


def flame_runs(y, y2, n=0, n_gap=0, reverse=False):
    """Index ranges [start, end] (inclusive) of the flame polygons for one series.

    The run logic of R's GeomFlame3$draw_group(): runs where y > y2 (y < y2
    when reverse) of at least n + 1 points are kept, and gaps of at most
    n_gap points between two kept runs are kept too (drawn as their own
    polygon). Comparisons with missing values count as not exceeding.
    """
    y, y2 = np.asarray(y, dtype=float), np.asarray(y2, dtype=float)
    exceed = y < y2 if reverse else y > y2
    if exceed.size == 0:
        return []
    change = np.flatnonzero(exceed[1:] != exceed[:-1]) + 1
    starts = np.concatenate([[0], change])
    ends = np.concatenate([change, [exceed.size]]) - 1
    values = exceed[starts]
    lengths = ends - starts + 1
    keep = values & (lengths >= max(1, n + 1))
    if n_gap > 0 and keep.sum() > 1:
        for i in range(1, len(keep) - 1):
            if (not keep[i] and not values[i] and lengths[i] <= n_gap
                    and keep[i - 1] and keep[i + 1]):
                keep[i] = True
    return [(int(s), int(e)) for s, e, k in zip(starts, ends, keep) if k]


def geom_flame3(x, y, y2, ax=None, n=0, n_gap=0, reverse=False, data=None, group=None,
                **kwargs):
    """Fill between a series and its threshold where it exceeds it ("flames").

    The matplotlib counterpart of the R ``geom_flame3()`` layer. Each run of
    exceedance becomes one polygon: along ``y`` over the run and back along
    ``y2``.

    Parameters
    ----------
    x, y, y2 : array-like, or column names when ``data`` is given
        x values (e.g. dates), the series (e.g. temperature) and the threshold.
    ax : matplotlib Axes, optional
        Axes to draw on; defaults to the current Axes.
    n : int
        Keep only runs of more than ``n`` consecutive exceeding points.
    n_gap : int
        Also fill gaps of at most ``n_gap`` points between two kept runs.
    reverse : bool
        Fill where ``y < y2`` (cold-spells) instead of ``y > y2``.
    data : pandas.DataFrame, optional
        Source of the columns named by ``x``, ``y``, ``y2`` and ``group``.
    group : array-like or column name, optional
        Series identifier, e.g. depth: each group is treated as its own
        series (rows keep their order within a group). R folds a mapped
        ``depth`` aesthetic into the grouping in the same way.
    **kwargs
        Passed to :class:`matplotlib.collections.PolyCollection`; the default
        is ``facecolor="salmon"`` with no edge.

    Returns
    -------
    matplotlib.collections.PolyCollection or None
        The flames, or ``None`` when nothing exceeds.
    """
    plt = _plt()
    from matplotlib.collections import PolyCollection

    ax = plt.gca() if ax is None else ax
    xv = _column(x, data)
    yv, y2v = _column(y, data, like=xv), _column(y2, data, like=xv)
    groups = (np.zeros(len(xv), dtype=int) if group is None
              else pd.factorize(pd.Series(_column(group, data)))[0])

    ax.xaxis.update_units(xv)
    xn = np.asarray(ax.xaxis.convert_units(xv), dtype=float)
    verts = []
    for g in range(groups.max() + 1 if len(groups) else 0):
        idx = np.flatnonzero(groups == g)
        for s, e in flame_runs(yv[idx], y2v[idx], n, n_gap, reverse):
            seg = idx[s:e + 1]
            upper, lower = (y2v[seg], yv[seg]) if reverse else (yv[seg], y2v[seg])
            verts.append(np.column_stack([np.concatenate([xn[seg], xn[seg][::-1]]),
                                          np.concatenate([upper, lower[::-1]])]))
    if not verts:
        return None
    kw = {"facecolor": _SALMON, "edgecolor": "none"} | kwargs
    coll = PolyCollection(verts, **kw)
    ax.add_collection(coll)
    ax.autoscale_view()
    return coll


def geom_lolli3(x, y, ax=None, data=None, color=_GREY35, size=1.0, **kwargs):
    """Lollipops (a stem from 0 to ``y`` and a point) for event metrics.

    The matplotlib counterpart of the R ``geom_lolli3()`` layer.

    Parameters
    ----------
    x, y : array-like, or column names when ``data`` is given
        Positions (e.g. peak dates) and values (e.g. maximum intensity).
    ax : matplotlib Axes, optional
        Axes to draw on; defaults to the current Axes.
    data : pandas.DataFrame, optional
        Source of the columns named by ``x`` and ``y``.
    color : colour
        Stem and point colour (R default ``"grey35"``).
    size : float
        Point diameter in millimetres, as in R (default 1).
    **kwargs
        Passed to :meth:`matplotlib.axes.Axes.scatter` for the points.

    Returns
    -------
    (LineCollection, PathCollection)
        The stems and the points.
    """
    plt = _plt()
    ax = plt.gca() if ax is None else ax
    xv = _column(x, data)
    yv = _column(y, data, like=xv)
    stems = ax.vlines(xv, 0, yv, colors=color, linewidth=0.5 * 72.27 / 25.4)
    pts_diam = size * 72 / 25.4  # mm -> points
    points = ax.scatter(xv, yv, s=pts_diam ** 2, color=color, zorder=3, **kwargs)
    return stems, points


def event_line3(sst_file, clim_file, lon, lat, depth=None, var_name=None, start_date=None,
                end_date=None, spread=150, metric="intensity_cumulative", event_file=None,
                coldSpells=False, ax=None):
    """Temperature, climatology and threshold at one pixel, with event flames.

    Parameters
    ----------
    sst_file : str, path-like, or directory
        SST NetCDF file, a directory of daily files, or a Zarr store (see
        :func:`ts2clm3`).
    clim_file : str or path-like
        Climatology from :func:`ts2clm3`, NetCDF or Zarr (see
        :func:`detect_event3`).
    lon, lat : float
        The pixel (the nearest grid cell is used).
    depth : float, optional
        Depth (metres) for a depth-resolved climatology; required then, and
        not allowed otherwise. The nearest level is used.
    var_name : str, optional
        SST variable name; auto-detected when ``None``.
    start_date, end_date : dates, optional
        Plot window (both must be given).
    spread : int
        Without a window: days either side of the peak of the pixel's
        strongest event in ``event_file``.
    metric : str
        Event metric used to choose that event.
    event_file : str, optional
        Events from :func:`detect_event3`, used to centre the window.
    coldSpells : bool
        Fill below the threshold, in blue, for cold-spells.
    ax : matplotlib Axes, optional
        Axes to draw on (Python only); a new figure is created otherwise.

    Returns
    -------
    matplotlib.axes.Axes
    """
    plt = _plt()
    eps = 0.01
    tmp_dir = local_dir(clim_file, event_file, sst_file)
    with clim_input(os.fspath(clim_file), tmp_dir) as f:
        cd = _core.hw3_read_clim_nc(f)
    lon_idx = int(np.argmin(np.abs(cd["lon"] - lon)))
    lat_idx = int(np.argmin(np.abs(cd["lat"] - lat)))
    ndepth = max(1, cd["ndepth"])
    if ndepth > 1 and depth is None:
        raise ValueError(f"clim_file is depth-resolved (ndepth = {ndepth}); specify 'depth' "
                         "(metres) to pick a level.")
    if ndepth <= 1 and depth is not None:
        raise ValueError("clim_file has no depth dimension; omit 'depth'.")

    # Resolve to the exact depth stored in clim_file so the SST read matches it
    idepth, resolved_depth, depth_range = 0, None, None
    if depth is not None:
        idepth = int(np.argmin(np.abs(cd["depth"] - depth)))
        resolved_depth = float(cd["depth"][idepth])
        depth_range = [resolved_depth, resolved_depth]
    px = (lon_idx * cd["nlat"] + lat_idx) * ndepth + idepth

    window = dict(lon_range=[lon - eps, lon + eps], lat_range=[lat - eps, lat + eps],
                  depth_range=depth_range)
    # A point outside the SST grid stops the reader before any data is read
    try:
        with nc_input(sst_file, tmp_dir, var_name, **window) as (src, vn, _):
            vn = "" if vn is None else vn
            if isinstance(src, (str, os.PathLike)) and os.path.isdir(src):
                files = list_nc(src, os.fspath(src))
                gd = _core.hw3_read_sst_multi(files, vn, **window)
            else:
                gd = _core.hw3_read_sst(os.fspath(src), vn, **window)
    except RuntimeError as e:
        if "No coordinates found" not in str(e):
            raise
        gd = {"nlon": 0, "nlat": 0}
    if gd["nlon"] * gd["nlat"] == 0:
        raise ValueError(f"No data found at lon = {lon}, lat = {lat}")

    temp = gd["sst"][:gd["ntime"]].copy()
    doy = np.array([_core.hw3_jd_to_doy(int(j)) for j in gd["time_days"]])
    seas = cd["seas"][px * 366 + doy - 1].copy()
    thresh = cd["thresh"][px * 366 + doy - 1].copy()
    if np.nanmedian(temp) > 200:  # kelvin
        temp, seas, thresh = temp - 273.15, seas - 273.15, thresh - 273.15

    df = pd.DataFrame({"t": jd_to_datetime(gd["time_days"]), "temp": temp, "seas": seas,
                       "thresh": thresh}).dropna()

    if start_date is not None and end_date is not None:
        s, e = to_datetime64([start_date, end_date]).astype("datetime64[ns]")
        df = df[(df["t"] >= s) & (df["t"] <= e)]
    elif event_file is not None:
        ev = _core.hw3_read_event_nc(os.fspath(event_file))
        peaks = jd_to_datetime(ev["ref_date_jd"] + ev["date_peak"])
        mask = (np.abs(ev["lon"] - lon) < eps) & (np.abs(ev["lat"] - lat) < eps)
        if resolved_depth is not None and len(ev["depth"]) > 0:
            mask &= np.abs(ev["depth"] - resolved_depth) < 1e-3
        if mask.any():
            # Cold-spell intensities are negative: the strongest is the minimum
            pick = np.argmin if coldSpells else np.argmax
            best = np.flatnonzero(mask)[pick(ev[metric][mask])]
            centre, span = peaks[best], np.timedelta64(int(spread), "D")
            df = df[(df["t"] >= centre - span) & (df["t"] <= centre + span)]
    if len(df) == 0:
        raise ValueError("No data in the specified date range.")

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4))
    geom_flame3("t", "temp", "thresh", ax=ax, data=df, reverse=bool(coldSpells),
                facecolor=_STEELBLUE3 if coldSpells else _SALMON)
    ax.plot(df["t"], df["temp"], color="black", linewidth=0.85, label="temp")
    ax.plot(df["t"], df["seas"], color=_GREY40, linestyle="--", linewidth=0.85, label="seas")
    ax.plot(df["t"], df["thresh"], color=_DARKGREEN, linewidth=0.85, label="thresh")
    import matplotlib.dates as mdates
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
    ax.set_ylabel("Temperature (°C)")
    ax.set_title(f"lon = {lon:.3f}, lat = {lat:.3f}" if resolved_depth is None
                 else f"lon = {lon:.3f}, lat = {lat:.3f}, depth = {resolved_depth:.1f} m")
    ax.grid(True, color="0.92")
    ax.set_axisbelow(True)
    return ax


def plot_metric3(event_file, metric="intensity_max", summary="mean", coastline=True,
                 ax=None, **kwargs):
    """Map of an event metric aggregated per pixel.

    Parameters
    ----------
    event_file : str or path-like
        Events from :func:`detect_event3`.
    metric : str
        Event metric, e.g. ``"intensity_max"`` (default), ``"duration"``,
        ``"intensity_cumulative"``, ``"rate_onset"``, or any ``relThresh``/
        ``abs`` variant.
    summary : {"mean", "max", "min", "sum", "count"}
        Aggregation across each pixel's events (computed in C++).
    coastline : bool
        Draw the 1:50m Natural Earth coastline (the one R draws), with
        ``cartopy``, which ``pip install "heatwave3[plot]"`` installs. Skipped
        when cartopy is missing, as R skips it without ``rnaturalearth``,
        ``rnaturalearthdata`` and ``sf``.
    ax : matplotlib Axes, optional
        Axes to draw on (Python only); a new figure is created otherwise.
    **kwargs
        Passed to :meth:`matplotlib.axes.Axes.pcolormesh` (e.g. ``cmap``);
        R passes ``...`` to ``scale_fill_viridis_c()``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    plt = _plt()
    agg = pd.DataFrame(_core.hw3_read_metric_summary(os.fspath(event_file), metric, summary))
    lons, lats = np.sort(agg["lon"].unique()), np.sort(agg["lat"].unique())
    grid = (agg.pivot_table(index="lat", columns="lon", values="value", aggfunc="first")
            .reindex(index=lats, columns=lons).to_numpy())

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))
    mesh = ax.pcolormesh(lons, lats, grid, shading="nearest", **({"cmap": "viridis"} | kwargs))
    ax.figure.colorbar(mesh, ax=ax, label=f"{metric} ({summary})")
    ax.set(xlabel="Longitude", ylabel="Latitude")
    if len(lats):  # R's coord_quickmap(): 1 degree of latitude vs longitude at mid-latitude
        ax.set_aspect(1 / math.cos(math.radians((lats.min() + lats.max()) / 2)))
    if coastline:
        _add_coastline(ax)
    return ax


def _add_coastline(ax) -> None:
    try:
        import cartopy.feature as cfeature
        from cartopy.mpl.geoaxes import GeoAxes
    except ImportError:
        return
    if isinstance(ax, GeoAxes):
        ax.add_feature(cfeature.COASTLINE, edgecolor="0.19", linewidth=0.85)
    else:
        from cartopy.io import shapereader
        from shapely.geometry import box
        xlim, ylim = ax.get_xlim(), ax.get_ylim()  # keep the data extent
        view = box(xlim[0], ylim[0], xlim[1], ylim[1])
        # Only the coastline inside the map, so saved vector files do not carry
        # every coastline on Earth
        for geom in shapereader.Reader(
                shapereader.natural_earth("50m", "physical", "coastline")).geometries():
            part = geom.intersection(view)
            for line in getattr(part, "geoms", [part]):
                if line.geom_type == "LineString" and not line.is_empty:
                    x, y = line.xy
                    ax.plot(x, y, color="0.19", linewidth=0.85)
        ax.set(xlim=xlim, ylim=ylim)
