"""Parity with the R package: run the same pipeline as make_golden.R through
the Python API and compare every output table with the R reference exactly.

The golden files are produced by python/tests/make_golden.R; regenerate them
after any change that alters outputs.
"""

import gzip
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import heatwave3 as hw3

ROOT = Path(__file__).resolve().parents[2]
EXT = ROOT / "inst" / "extdata"
GOLDEN = Path(__file__).resolve().parent / "data" / "golden"
SST, GLO, GLO_DIR = EXT / "sst_test.nc", EXT / "glorys_depth_test.nc", EXT / "glorys_depth_test_annual"
CP = ("1982-01-01", "2011-12-31")
INT_MIN = np.iinfo(np.int32).min

# detect_blob3() output is compared to rtol 1e-14 rather than exactly, for two
# reasons (decision 2026-10-11): event columns aggregated with sum()/mean() use
# 80-bit long double in R but not in pandas (1-2 ulp), and C++ blob metrics
# differ in the last bit between CPU types (<= ~4e-16 relative), because
# compilers fuse multiply-adds (FMA) by default on ARM but not on x86-64. The
# golden files may come from either CPU type, and exactly tied blobs can be
# ranked in a different order, so on every platform blobs are matched by id
# and compared without their rank. Everything else must match exactly.
BLOB_KEYS = {"event": ["blob"], "daily": ["blob", "t_idx"], "voxel": ["blob", "t_idx", "lon", "lat"]}

pytestmark = pytest.mark.skipif(not (SST.exists() and GOLDEN.exists()),
                                reason="test data or golden files not available")


def load_golden(name):
    meta = json.loads((GOLDEN / f"{name}.json").read_text())
    n = meta["nrow"]
    raw = gzip.decompress((GOLDEN / f"{name}.bin.gz").read_bytes())
    pos, cols = 0, {}
    for c in meta["columns"]:
        if c["type"] == "character":
            cols[c["name"]] = ("character", [None if v is None else v for v in (c["values"] or [])])
            continue
        size = 8 if c["type"] == "double" else 4
        dtype = "<f8" if size == 8 else "<i4"
        cols[c["name"]] = (c["type"], np.frombuffer(raw, dtype, n, pos))
        pos += size * n
    assert pos == len(raw)
    return n, cols


def as_r(series, rtype):
    """A pandas column in the golden file's representation."""
    if rtype == "double":
        assert series.dtype == np.float64, series.dtype
        return series.to_numpy()
    if rtype == "Date":
        assert np.issubdtype(series.dtype, np.datetime64), series.dtype
        return series.to_numpy().astype("datetime64[D]").astype(np.int64)
    if rtype == "logical":
        assert series.dtype == bool, series.dtype
        return series.to_numpy().astype(np.int64)
    if rtype == "integer":
        assert str(series.dtype) in ("int32", "Int32"), series.dtype
        return series.astype("Int64").fillna(INT_MIN).to_numpy(dtype=np.int64)
    if rtype == "character":
        return [None if pd.isna(v) else v for v in series]
    raise AssertionError(rtype)


def assert_matches_r(df, name):
    n, cols = load_golden(name)
    assert list(df.columns) == list(cols), f"{name}: columns differ"
    assert len(df) == n, f"{name}: {len(df)} rows vs {n} in R"
    for col, (rtype, ref) in cols.items():
        got = as_r(df[col], rtype)
        if rtype == "character":
            assert got == ref, f"{name}.{col} differs"
        elif rtype == "double":
            assert np.array_equal(got, ref, equal_nan=True), f"{name}.{col} differs"
        else:
            assert np.array_equal(got, ref.astype(np.int64)), f"{name}.{col} differs"


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    wd = tmp_path_factory.mktemp("parity")
    s, p, k, g, m = (wd / x for x in "spkgm")
    hw3.ts2clm3(SST, s, CP, quiet=True)
    hw3.detect_event3(SST, s, category=True, daily="also", quiet=True)
    shutil.copy(f"{s}_clim.nc", f"{p}_clim.nc")
    hw3.detect_event3(SST, p, protoEvent=True, quiet=True)
    hw3.detect_event3(SST, p, quiet=True)
    hw3.detect3(SST, k, CP, pctile=10, coldSpells=True, category=True, quiet=True)
    hw3.exceedance3(SST, wd / "exceed.nc", threshold=296)
    hw3.ts2clm3(GLO, g, ("1993-01-01", "2019-12-31"), depth_range=(0, 10), quiet=True)
    hw3.detect_event3(GLO, g, category=True, daily="also", quiet=True)
    hw3.detect3(GLO_DIR, m, ("1993-01-01", "2019-12-31"), depth_range=(0, 10), quiet=True)
    return {x: str(wd / x) for x in "spkgm"} | {"exceed": str(wd / "exceed.nc")}


CASES = {
    "clim": lambda r: hw3.hw3_export(r["s"] + "_clim.nc"),
    "clim_sub": lambda r: hw3.hw3_export(r["s"] + "_clim.nc", vars="thresh",
                                         lon_range=(26.5, 26.55), n=500),
    "events": lambda r: hw3.hw3_export(r["s"] + "_events.nc"),
    "events_sub": lambda r: hw3.hw3_export(r["s"] + "_events.nc",
                                           vars=["duration", "intensity_max"],
                                           lat_range=(-34.2, -34.05),
                                           time_range=("2010-01-01", "2012-12-31")),
    "daily": lambda r: hw3.hw3_export(r["s"] + "_events_daily.nc"),
    "daily_sub": lambda r: hw3.hw3_export(r["s"] + "_events_daily.nc",
                                          vars=["temp", "event", "category"],
                                          time_range=("2015-01-01", "2015-12-31")),
    "category3": lambda r: hw3.category3(r["s"] + "_events.nc"),
    "block_average3": lambda r: hw3.block_average3(r["s"] + "_events.nc"),
    "category_daily3": lambda r: hw3.category_daily3(SST, r["s"] + "_clim.nc",
                                                     r["s"] + "_events.nc",
                                                     ("2019-06-01", "2019-12-31")),
    "proto_sub": lambda r: hw3.hw3_export(r["p"] + "_protoevents.nc",
                                          time_range=("2016-01-01", "2016-12-31")),
    "category3_fallback": lambda r: hw3.category3(r["p"] + "_events.nc", r["p"] + "_clim.nc",
                                                  hemisphere="north", sst_file=SST),
    "cold_events": lambda r: hw3.hw3_export(r["k"] + "_events.nc"),
    "category_daily3_cold": lambda r: hw3.category_daily3(SST, r["k"] + "_clim.nc",
                                                          r["k"] + "_events.nc",
                                                          ("2005-01-01", "2005-12-31"),
                                                          coldSpells=True),
    "exceedance3": lambda r: hw3.hw3_export(r["exceed"]),
    "depth_clim": lambda r: hw3.hw3_export(r["g"] + "_clim.nc"),
    "depth_events": lambda r: hw3.hw3_export(r["g"] + "_events.nc"),
    "depth_daily_sub": lambda r: hw3.hw3_export(r["g"] + "_events_daily.nc", depth_range=(0, 3),
                                                time_range=("2010-01-01", "2010-03-01")),
    "depth_events_sub": lambda r: hw3.hw3_export(r["g"] + "_events.nc", depth_range=(0, 3)),
    "multi_events": lambda r: hw3.hw3_export(r["m"] + "_events.nc"),
}

BLOBS = {
    "blob": lambda r: hw3.detect_blob3(SST, r["s"] + "_clim.nc",
                                       return_=("event", "daily", "voxel")),
    "blob_filtered": lambda r: hw3.detect_blob3(SST, r["s"] + "_clim.nc", wrapDateline=True,
                                                minVoxels=20, topN=5, rankBy="peakArea_km2"),
    "blob_cold": lambda r: hw3.detect_blob3(SST, r["k"] + "_clim.nc", coldSpells=True),
    "blob_depth": lambda r: hw3.detect_blob3(GLO, r["g"] + "_clim.nc", minVoxels=10,
                                             rankBy="totalVolume_km3",
                                             return_=("event", "daily", "voxel")),
    "blob_multi": lambda r: hw3.detect_blob3(GLO_DIR, r["g"] + "_clim.nc", topN=3),
}
BLOB_CASES = [(b, part) for b, parts in [("blob", "event daily voxel"),
                                         ("blob_filtered", "event daily"),
                                         ("blob_cold", "event daily"),
                                         ("blob_depth", "event daily voxel"),
                                         ("blob_multi", "event daily")]
              for part in parts.split()]


def test_every_golden_file_has_a_case():
    expected = set(CASES) | {f"{b}_{part}" for b, part in BLOB_CASES}
    from test_plot import PLOT_GOLDEN  # plot references are tested in test_plot.py
    assert {p.name.removesuffix(".json") for p in GOLDEN.glob("*.json")} == expected | PLOT_GOLDEN


@pytest.mark.parametrize("name", sorted(CASES))
def test_matches_r(runs, name):
    assert_matches_r(CASES[name](runs), name)


@pytest.fixture(scope="module")
def blob_runs(runs):
    return {name: f(runs) for name, f in BLOBS.items()}


def assert_blob_matches_r_unordered(df, name, keys):
    """Blob output against R, matched on `keys`, ignoring order and rank."""
    n, cols = load_golden(name)
    assert list(df.columns) == list(cols) and len(df) == n, name
    ref = pd.DataFrame({c: (v if t == "character" else np.asarray(v)) for c, (t, v) in cols.items()})
    got = pd.DataFrame({c: as_r(df[c], cols[c][0]) for c in df.columns})
    keys = keys + (["depth"] if "depth" in ref and "voxel" in name else [])
    ref, got = (x.sort_values(keys, kind="stable").reset_index(drop=True) for x in (ref, got))
    for c, (rtype, _) in cols.items():
        if c == "rank":
            continue
        if rtype == "double":
            np.testing.assert_allclose(got[c], ref[c], rtol=1e-14, atol=0, err_msg=f"{name}.{c}")
        else:
            assert np.array_equal(got[c].to_numpy(), ref[c].to_numpy()), f"{name}.{c} differs"


@pytest.mark.parametrize("blob, part", BLOB_CASES)
def test_blob_matches_r(blob_runs, blob, part):
    res = blob_runs[blob]
    assert set(res) == {p for b, p in BLOB_CASES if b == blob}
    assert_blob_matches_r_unordered(res[part], f"{blob}_{part}", BLOB_KEYS[part])
