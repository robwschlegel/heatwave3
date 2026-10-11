"""Plotting: each Python plot draws exactly the data the R plot draws.

R references (from make_golden.R): flame polygons from R's own
GeomFlame3$draw_group(), event_line3()'s plotted data and title, and
plot_metric3()'s raster data.
"""

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import heatwave3 as hw3  # noqa: E402
import test_parity_r as tp  # noqa: E402

PLOT_GOLDEN = {"flame_input", "flame_default", "flame_n3", "flame_gap", "flame_reverse",
               "line_event", "line_event_title", "line_window", "line_window_title",
               "line_cold", "line_cold_title", "line_depth", "line_depth_title",
               "metric_default", "metric_duration_max"}

pytestmark = pytest.mark.skipif(not (tp.SST.exists() and tp.GOLDEN.exists()),
                                reason="test data or golden files not available")


def golden_df(name):
    n, cols = tp.load_golden(name)
    return pd.DataFrame({c: (v if t == "character" else np.asarray(v)) for c, (t, v) in cols.items()})


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    wd = tmp_path_factory.mktemp("plot")
    s, k, g = (str(wd / x) for x in "skg")
    hw3.detect3(tp.SST, s, tp.CP, quiet=True)
    hw3.detect3(tp.SST, k, tp.CP, pctile=10, coldSpells=True, quiet=True)
    hw3.detect3(tp.GLO, g, ("1993-01-01", "2019-12-31"), depth_range=(0, 10), quiet=True)
    return {"s": s, "k": k, "g": g}


@pytest.mark.parametrize("name, kwargs", [
    ("flame_default", {}), ("flame_n3", {"n": 3}), ("flame_gap", {"n": 1, "n_gap": 2}),
    ("flame_reverse", {"reverse": True})])
def test_flame_polygons_match_r(name, kwargs):
    d = golden_df("flame_input")
    _, ax = plt.subplots()
    coll = hw3.geom_flame3(d["x"], d["y"], d["y2"], ax=ax, **kwargs)
    ref = golden_df(name)
    want = [ref.loc[ref["id"] == i, ["x", "y"]].to_numpy() for i in sorted(ref["id"].unique())]
    got = [p.vertices[:len(w)] for p, w in zip(coll.get_paths(), want)]
    assert len(coll.get_paths()) == len(want) > 0
    for gv, wv in zip(got, want):
        assert np.array_equal(gv, wv)


def test_flame_with_data_columns_and_groups():
    d = pd.DataFrame({"x": np.tile(np.arange(10.0), 2), "y": np.tile([0, 1, 1, 0] * 2 + [1, 1], 2),
                      "thr": 0.5, "depth": np.repeat([1.0, 5.0], 10)})
    _, ax = plt.subplots()
    # Grouped by depth, each depth is its own series: 3 runs per depth
    grouped = hw3.geom_flame3("x", "y", "thr", ax=ax, data=d, group="depth")
    assert len(grouped.get_paths()) == 6
    assert hw3.geom_flame3("x", "y", 9, ax=ax, data=d) is None  # nothing exceeds


def test_flame_accepts_dates():
    t = pd.date_range("2020-01-01", periods=6)
    _, ax = plt.subplots()
    coll = hw3.geom_flame3(t, [0, 2, 2, 0, 2, 0], 1, ax=ax)
    assert len(coll.get_paths()) == 2


def test_lolli_stems_and_points():
    x, y = np.array([1.0, 2.0, 3.0]), np.array([0.5, 2.0, 1.0])
    _, ax = plt.subplots()
    stems, points = hw3.geom_lolli3(x, y, ax=ax)
    segs = stems.get_segments()
    assert [tuple(s[:, 1]) for s in segs] == [(0, 0.5), (0, 2.0), (0, 1.0)]
    assert np.array_equal(points.get_offsets(), np.column_stack([x, y]))


LINES = {
    "line_event": lambda r: hw3.event_line3(tp.SST, r["s"] + "_clim.nc", 26.525, -34.125,
                                            event_file=r["s"] + "_events.nc"),
    "line_window": lambda r: hw3.event_line3(tp.SST, r["s"] + "_clim.nc", 26.575, -34.025,
                                             start_date="2015-01-01", end_date="2016-06-30"),
    "line_cold": lambda r: hw3.event_line3(tp.SST, r["k"] + "_clim.nc", 26.525, -34.075,
                                           event_file=r["k"] + "_events.nc", coldSpells=True,
                                           metric="intensity_max", spread=60),
    "line_depth": lambda r: hw3.event_line3(tp.GLO, r["g"] + "_clim.nc", 25.16667, -34.91667,
                                            depth=5, start_date="2010-01-01",
                                            end_date="2010-12-31"),
}


@pytest.mark.parametrize("name", sorted(LINES))
def test_event_line3_draws_the_r_data(runs, name):
    ax = LINES[name](runs)
    ref = golden_df(name)
    days = np.asarray(ax.lines[0].get_xdata(), dtype="datetime64[D]").astype(np.int64)
    assert np.array_equal(days, ref["t"].to_numpy())
    for line, col in zip(ax.lines, ["temp", "seas", "thresh"]):
        assert line.get_label() == col
        assert np.array_equal(np.asarray(line.get_ydata(), float), ref[col].to_numpy())
    assert ax.get_title() == golden_df(f"{name}_title")["title"][0]
    flames = [c for c in ax.collections if isinstance(c, matplotlib.collections.PolyCollection)]
    if flames:
        want = "#4f94cd" if name == "line_cold" else "#fa8072"
        assert matplotlib.colors.to_hex(flames[0].get_facecolor()[0]) == want


def test_event_line3_errors(runs):
    with pytest.raises(RuntimeError, match="No coordinates found"):
        hw3.event_line3(tp.SST, runs["s"] + "_clim.nc", 999, 999)
    with pytest.raises(ValueError, match="No data in the specified date range"):
        hw3.event_line3(tp.SST, runs["s"] + "_clim.nc", 26.525, -34.125,
                        start_date="1900-01-01", end_date="1900-01-02")
    with pytest.raises(ValueError, match="specify 'depth'"):
        hw3.event_line3(tp.GLO, runs["g"] + "_clim.nc", 25.16667, -34.91667)
    with pytest.raises(ValueError, match="omit 'depth'"):
        hw3.event_line3(tp.SST, runs["s"] + "_clim.nc", 26.525, -34.125, depth=5)


@pytest.mark.parametrize("name, product, kwargs, label", [
    ("metric_default", "s", {}, "intensity_max (mean)"),
    ("metric_duration_max", "g", {"metric": "duration", "summary": "max", "cmap": "magma"},
     "duration (max)")])
def test_plot_metric3_maps_the_r_data(runs, name, product, kwargs, label):
    ax = hw3.plot_metric3(runs[product] + "_events.nc", coastline=False, **kwargs)
    ref = golden_df(name)
    lons, lats = np.sort(ref["lon"].unique()), np.sort(ref["lat"].unique())
    want = (ref.pivot_table(index="lat", columns="lon", values="value", aggfunc="first")
            .reindex(index=lats, columns=lons).to_numpy())
    got = np.ma.filled(ax.collections[0].get_array().astype(float), np.nan).reshape(want.shape)
    assert np.array_equal(got, want, equal_nan=True)
    assert ax.figure.axes[-1].get_ylabel() == label


def test_plot_metric3_coastline_is_optional(runs):
    pytest.importorskip("cartopy", reason="cartopy installed: coastlines would be drawn")
    hw3.plot_metric3(runs["s"] + "_events.nc", coastline=True)


def test_plot_metric3_without_cartopy_skips_coastline(runs):
    try:
        import cartopy  # noqa: F401
        pytest.skip("cartopy is installed")
    except ImportError:
        pass
    ax = hw3.plot_metric3(runs["s"] + "_events.nc", coastline=True)
    assert len(ax.lines) == 0
