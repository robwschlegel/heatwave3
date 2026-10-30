"""Python-specific behaviour of the public API: argument handling, errors,
export formats and return values. Numerical parity with R is covered by
test_parity_r.py."""

import datetime as dt
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import heatwave3 as hw3

EXT = Path(__file__).resolve().parents[2] / "inst" / "extdata"
SST = EXT / "sst_test.nc"
CP = ("1982-01-01", "2011-12-31")

pytestmark = pytest.mark.skipif(not SST.exists(), reason="sst_test.nc not available")


@pytest.fixture(scope="module")
def stem(tmp_path_factory):
    s = tmp_path_factory.mktemp("api") / "s"
    hw3.detect3(SST, s, CP, category=True, daily="also", quiet=True)
    return str(s)


def test_detect3_returns_written_paths(stem):
    out = hw3.detect3(SST, stem + "x", CP, quiet=True)
    assert list(out) == ["clim", "events"]
    assert all(Path(p).exists() for p in out.values())


def test_accepts_path_objects_and_date_objects(tmp_path):
    a = hw3.ts2clm3(SST, tmp_path / "a", CP, quiet=True)
    b = hw3.ts2clm3(str(SST), str(tmp_path / "b"), (dt.date(1982, 1, 1), np.datetime64("2011-12-31")),
                    quiet=True)
    pd.testing.assert_frame_equal(hw3.hw3_export(a), hw3.hw3_export(b))


def test_false_means_zero_like_r(tmp_path):
    # R: if (is.numeric(x)) as.integer(x) else 0L. In Python bool is an int,
    # so True must not silently become 1.
    ref = hw3.hw3_export(hw3.ts2clm3(SST, tmp_path / "zero", CP, maxPadLength=0, quiet=True))
    for flag in (False, True):
        got = hw3.hw3_export(hw3.ts2clm3(SST, tmp_path / f"f{flag}", CP, maxPadLength=flag,
                                         quiet=True))
        pd.testing.assert_frame_equal(got, ref)


def test_roundClm_false_disables_rounding(tmp_path):
    raw = hw3.hw3_export(hw3.ts2clm3(SST, tmp_path / "raw", CP, roundClm=False, quiet=True))
    rounded = hw3.hw3_export(hw3.ts2clm3(SST, tmp_path / "rnd", CP, quiet=True))
    assert not np.array_equal(raw["seas"], rounded["seas"])
    np.testing.assert_allclose(raw["seas"], rounded["seas"], atol=5e-5)


def test_directory_and_file_list_inputs_agree(tmp_path):
    glo_dir = EXT / "glorys_depth_test_annual"
    files = sorted(glo_dir.glob("*.nc"))
    a = hw3.ts2clm3(glo_dir, tmp_path / "dir", ("1993-01-01", "2019-12-31"), quiet=True)
    b = hw3.ts2clm3(files, tmp_path / "list", ("1993-01-01", "2019-12-31"), quiet=True)
    pd.testing.assert_frame_equal(hw3.hw3_export(a), hw3.hw3_export(b))


def test_no_events_skips_events_file(stem, tmp_path):
    out = hw3.detect_event3(SST, tmp_path / "none", clim_file=stem + "_clim.nc",
                            minDuration=100_000, quiet=True)
    assert out == {}


def test_column_dtypes(stem):
    daily = hw3.hw3_export(stem + "_events_daily.nc", n=10)
    assert str(daily["event_no"].dtype) == "Int32"
    assert str(daily["category"].dtype) == "Int32"
    assert daily["event"].dtype == bool
    assert np.issubdtype(daily["t"].dtype, np.datetime64)
    cat = hw3.category3(stem + "_events.nc")
    assert np.issubdtype(cat["peak_date"].dtype, np.datetime64)
    assert cat["event_no"].dtype == np.int32


@pytest.mark.parametrize("call, exc, msg", [
    (lambda s: hw3.ts2clm3(SST, s, ("1982-01-01",)), ValueError, "climatologyPeriod"),
    (lambda s: hw3.ts2clm3(SST, s, CP, depth=0, depth_range=(0, 5)), ValueError, "only one of"),
    (lambda s: hw3.ts2clm3(SST, "/no/such/dir/x", CP), FileNotFoundError, "Output directory"),
    (lambda s: hw3.ts2clm3("/no/such.nc", s, CP), FileNotFoundError, "Input file does not exist"),
    (lambda s: hw3.detect_event3(SST, s, daily="sometimes"), ValueError, "'daily' must be"),
    (lambda s: hw3.detect_event3(SST, s, daily="also", protoEvent=True), ValueError, "mutually"),
    (lambda s: hw3.detect_event3(SST, s, slim=True), ValueError, "'slim' applies"),
    (lambda s: hw3.detect3(SST, s + "slim", CP, slim=True), ValueError, "'slim' applies"),
    (lambda s: hw3.detect_event3(SST, s + "nope"), FileNotFoundError, "Climatology file"),
    (lambda s: hw3.detect_event3(SST, s, hemisphere="east"), ValueError, "'hemisphere'"),
    (lambda s: hw3.category3(s + "_events.nc", hemisphere="up"), ValueError, "'hemisphere'"),
    (lambda s: hw3.category_daily3(SST, s + "_clim.nc", s + "_events.nc", ("2019-01-01",)),
     ValueError, "two dates"),
    (lambda s: hw3.category_daily3(SST, s + "_clim.nc", s + "_events.nc", ("x", "y")),
     ValueError, "must be dates"),
    (lambda s: hw3.hw3_export(s + "_clim.nc", vars="nope"), RuntimeError, "not in this"),
    (lambda s: hw3.hw3_export(s + "_clim.nc", file_out=s + ".rds"), ValueError, "R format"),
    (lambda s: hw3.hw3_export(s + "_clim.nc", file_out=s + ".xlsx"), ValueError, ".parquet or .zarr"),
])
def test_errors(stem, call, exc, msg):
    with pytest.raises(exc, match=msg):
        call(stem)


def test_clim_time_range_warns(stem):
    with pytest.warns(UserWarning, match="ignored for a climatology"):
        hw3.hw3_export(stem + "_clim.nc", time_range=("2000-01-01", "2000-12-31"), n=5)


def test_category3_S_is_deprecated(stem):
    with pytest.warns(DeprecationWarning):
        a = hw3.category3(stem + "_events.nc", S=False)
    pd.testing.assert_frame_equal(a, hw3.category3(stem + "_events.nc", hemisphere="north"))


@pytest.mark.parametrize("product", ["_clim.nc", "_events.nc", "_events_daily.nc"])
def test_chunked_csv_matches_whole_table(stem, tmp_path, product):
    whole = hw3.hw3_export(stem + product)
    out = hw3.hw3_export(stem + product, file_out=tmp_path / "x.csv", chunk_size=1000)
    back = pd.read_csv(out)
    assert list(back.columns) == list(whole.columns)
    assert len(back) == len(whole)
    num = [c for c in whole.columns if whole[c].dtype == np.float64]
    np.testing.assert_allclose(back[num].to_numpy(), whole[num].to_numpy(), rtol=1e-12)


@pytest.mark.parametrize("product", ["_clim.nc", "_events.nc", "_events_daily.nc"])
def test_chunked_parquet_round_trips_exactly(stem, tmp_path, product):
    pytest.importorskip("pyarrow")
    whole = hw3.hw3_export(stem + product)
    out = hw3.hw3_export(stem + product, file_out=tmp_path / "x.parquet", chunk_size=1000)
    back = pd.read_parquet(out)
    for c in whole.columns:
        if pd.api.types.is_datetime64_any_dtype(whole[c]):
            back[c] = back[c].astype(whole[c].dtype)
    pd.testing.assert_frame_equal(back, whole)


def test_subset_export_writes_file(stem, tmp_path):
    out = hw3.hw3_export(stem + "_events.nc", file_out=tmp_path / "e.csv", n=7)
    assert len(pd.read_csv(out)) == 7


def test_threads():
    old = hw3.getHW3threads()
    try:
        assert hw3.setHW3threads(3) == 3
        assert hw3.getHW3threads() == 3
    finally:
        hw3.setHW3threads(old)


def test_quiet_false_prints_summary(stem, capsys):
    hw3.ts2clm3(SST, stem + "q", CP)
    assert "Climatology written to" in capsys.readouterr().out


def test_quiet_silences_core_messages(tmp_path, capfd):
    capfd.readouterr()
    hw3.detect3(SST, tmp_path / "q", CP, daily="also", quiet=True)
    out, err = capfd.readouterr()
    assert out == "" and err == ""
    hw3.detect_event3(SST, tmp_path / "q", quiet=False)
    out, _ = capfd.readouterr()
    assert "Detecting events" in out and "Events written to" in out


def test_quiet_still_reports_warnings(tmp_path, capfd):
    files = [SST, tmp_path / "missing.nc", SST]
    capfd.readouterr()
    hw3.ts2clm3(files, tmp_path / "w", CP, skip_bad_files=True, quiet=True)
    out, err = capfd.readouterr()
    assert out == ""
    assert "failed to read" in err and "missing.nc" in err


def test_quiet_is_restored_after_an_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        hw3.detect_event3(SST, tmp_path / "x", clim_file=tmp_path / "none.nc", quiet=True)
    with pytest.raises(RuntimeError):
        hw3.ts2clm3(SST, tmp_path / "y", CP, var_name="nope", quiet=True)
    assert hw3._core.hw3_set_quiet(False) is False


def test_daily_categories_only_on_event_days(stem):
    d = hw3.hw3_export(stem + "_events_daily.nc")
    assert d.loc[~d["event"], "category"].isna().all()
    assert d.loc[d["event"], "category"].notna().any()
    win = ("2015-01-01", "2016-12-31")
    cd = hw3.category_daily3(SST, stem + "_clim.nc", stem + "_events.nc", win)
    dw = d[(d["t"] >= win[0]) & (d["t"] <= win[1])]
    key = ["lon", "lat", "t"]
    a, b = dw.sort_values(key).reset_index(drop=True), cd.sort_values(key).reset_index(drop=True)
    pd.testing.assert_series_equal(a["category"], b["category"])
    pd.testing.assert_series_equal(a["event"], b["event"])


def test_exceedance3_quiet(tmp_path, capfd):
    capfd.readouterr()
    out = hw3.exceedance3(SST, tmp_path / "e.nc", threshold=296, quiet=True)
    assert capfd.readouterr().out == ""
    assert Path(out).exists()
    hw3.exceedance3(SST, tmp_path / "f.nc", threshold=296)
    assert "Detecting events" in capfd.readouterr().out


def test_detect_blob3_return_handling(stem):
    clim = stem + "_clim.nc"
    assert list(hw3.detect_blob3(SST, clim, return_="daily")) == ["daily"]
    # Components come back in a fixed order, duplicates dropped (as in R)
    res = hw3.detect_blob3(SST, clim, return_=["voxel", "event", "voxel"])
    assert list(res) == ["event", "voxel"]
    ev = res["event"]
    assert list(ev["rank"]) == list(range(1, len(ev) + 1))
    assert ev["cumI_km2_day"].is_monotonic_decreasing
    assert np.issubdtype(ev["date_peak"].dtype, np.datetime64)


@pytest.mark.parametrize("kwargs, exc, msg", [
    (dict(return_=("event", "blobs")), ValueError, "unsupported value: blobs"),
    (dict(rankBy="totalVolume_km3"), ValueError, "rankBy must be one of"),
    (dict(connectivity=8), ValueError, "connectivity must be one of"),
    (dict(sst_file="/no/such.nc"), FileNotFoundError, "sst_file file does not exist"),
])
def test_detect_blob3_errors(stem, kwargs, exc, msg):
    args = dict(sst_file=SST, clim_file=stem + "_clim.nc") | kwargs
    with pytest.raises(exc, match=msg):
        hw3.detect_blob3(**args)


def test_detect_blob3_without_blobs_returns_empty_frames(tmp_path):
    clim = str(tmp_path / "never_clim.nc")
    hw3._core.hw3_write_const_clim(str(SST), clim, "", 1000.0)  # nothing exceeds 1000 K
    res = hw3.detect_blob3(SST, clim, return_=("event", "voxel"))
    assert list(res) == ["event", "voxel"]
    assert all(df.empty for df in res.values())


def test_category3_needs_sst_without_inline_categories(tmp_path):
    hw3.detect3(SST, tmp_path / "nc", CP, quiet=True)
    ev, clim = str(tmp_path / "nc_events.nc"), str(tmp_path / "nc_clim.nc")
    with pytest.raises(ValueError, match="clim_file and sst_file are required"):
        hw3.category3(ev)
    with pytest.raises(RuntimeError, match="sst_file are required"):
        hw3.category3(ev, clim)
    cats = hw3.category3(ev, clim, sst_file=SST)
    assert cats["category"].notna().all()
    # whole percentages of the duration; gap days below the threshold are not
    # counted, so (as in heatwaveR) they may add up to less than 100
    p = cats[["p_moderate", "p_strong", "p_severe", "p_extreme"]]
    assert ((p % 1 == 0) & (p >= 0)).all().all() and (p.sum(axis=1) <= 101).all()


def test_compression_must_be_0_to_9(tmp_path):
    for bad in (-1, 10, 1.5, "1", True, None):
        with pytest.raises(ValueError, match="compression must be a whole number from 0 to 9"):
            hw3.ts2clm3(SST, tmp_path / "x", CP, compression=bad, quiet=True)


def _products(stem):
    out = {k: hw3.hw3_export(f"{stem}{x}") for k, x in
           (("clim", "_clim.nc"), ("events", "_events.nc"), ("daily", "_events_daily.nc"))}
    out["events"] = out["events"].sort_values(["lon", "lat", "event_no"], ignore_index=True)
    return out


def test_compression_levels_give_same_values(tmp_path):
    stems = {lv: tmp_path / f"l{lv}" for lv in (0, 1, 9)}
    for lv, stem in stems.items():
        hw3.detect3(SST, stem, CP, category=True, daily="also", n_threads=1,
                    compression=lv, quiet=True)
    ref = _products(stems[1])
    for lv in (0, 9):
        for k, df in _products(stems[lv]).items():
            pd.testing.assert_frame_equal(df, ref[k])
    size = {lv: os.path.getsize(f"{s}_events_daily.nc") for lv, s in stems.items()}
    assert size[0] > 2 * size[1] and size[9] <= size[1]


def test_parallel_writer_matches_serial_writer(tmp_path, monkeypatch):
    def run(serial, stem, lv):
        monkeypatch.setenv("HEATWAVE3_SERIAL_WRITE", "1" if serial else "0")
        hw3.detect3(SST, stem, CP, category=True, daily="also", n_threads=4,
                    compression=lv, quiet=True)
        hw3.detect_event3(SST, f"{stem}p", clim_file=f"{stem}_clim.nc", protoEvent=True,
                          n_threads=4, compression=lv, quiet=True)
        out = _products(stem)
        out["proto"] = hw3.hw3_export(f"{stem}p_protoevents.nc")
        return out
    for lv in (1, 6):
        par = run(False, tmp_path / f"par{lv}", lv)
        ser = run(True, tmp_path / f"ser{lv}", lv)
        for k in par:
            pd.testing.assert_frame_equal(par[k], ser[k])


def test_slim_daily_has_only_intensity_and_category(stem, tmp_path):
    full = hw3.hw3_export(stem + "_events_daily.nc")
    for daily in ("also", "only"):
        out = hw3.detect_event3(SST, tmp_path / daily, clim_file=stem + "_clim.nc",
                                category=True, daily=daily, slim=True, quiet=True)
        slim = hw3.hw3_export(out["daily"])
        assert list(slim.columns) == ["lon", "lat", "t", "intensity", "category"]
        pd.testing.assert_frame_equal(slim, full[list(slim.columns)])
    assert "events" not in out  # daily="only"
    assert os.path.getsize(tmp_path / "also_events_daily.nc") < 0.7 * os.path.getsize(
        stem + "_events_daily.nc")


def test_slim_leaves_the_events_table_unchanged(stem, tmp_path):
    hw3.detect_event3(SST, tmp_path / "s", clim_file=stem + "_clim.nc", category=True,
                      daily="also", slim=True, quiet=True)
    srt = ["lon", "lat", "event_no"]
    pd.testing.assert_frame_equal(
        hw3.hw3_export(tmp_path / "s_events.nc").sort_values(srt, ignore_index=True),
        hw3.hw3_export(stem + "_events.nc").sort_values(srt, ignore_index=True))


def test_climatology_long_names_are_exact(tmp_path):
    netCDF4 = pytest.importorskip("netCDF4")
    clim = hw3.ts2clm3(SST, tmp_path / "c", CP, var=True, quiet=True)
    with netCDF4.Dataset(clim) as nc:
        # exact text: nothing read past the end of the string literal
        assert nc["seas"].getncattr("long_name") == "seasonal climatology"
        assert nc["thresh"].getncattr("long_name") == "threshold climatology"
        assert nc["var"].getncattr("long_name") == "variance climatology"
