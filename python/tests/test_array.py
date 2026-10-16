"""detect_event_array3(): per-day detection on arrays, checked against the
per-day product of detect_event3() on the same days and climatology."""

from pathlib import Path

import numpy as np
import pytest

import heatwave3 as hw3
from heatwave3 import _core

EXT = Path(__file__).resolve().parents[2] / "inst" / "extdata"
SST = EXT / "sst_test.nc"
CP = ("1982-01-01", "2011-12-31")
ALL = ("event", "event_no", "intensity", "category")

pytestmark = pytest.mark.skipif(not SST.exists(), reason="sst_test.nc not available")


def per_day_climatology(clim_file, time_jd):
    """seas and thresh of a ts2clm3 file matched to each day: (pixels, days)."""
    clim = _core.hw3_read_clim_nc(str(clim_file))
    npix = clim["nlon"] * clim["nlat"]
    doy = np.array([_core.hw3_jd_to_doy(int(jd)) for jd in time_jd]) - 1
    return (clim["seas"].reshape(npix, 366)[:, doy], clim["thresh"].reshape(npix, 366)[:, doy])


def read_daily(daily_file):
    daily = _core.hw3_read_daily_nc(str(daily_file))
    shape = (daily["nlon"] * daily["nlat"], daily["ntime"])
    return {name: daily[name].reshape(shape) for name in ALL}


def assert_same_as_daily(arrays, daily):
    event, event_no, intensity, category = arrays
    np.testing.assert_array_equal(event, daily["event"].astype(bool))
    np.testing.assert_array_equal(event_no, daily["event_no"])
    # the per-day file stores intensity as float, as the array output
    np.testing.assert_array_equal(intensity, daily["intensity"].astype(np.float32))
    np.testing.assert_array_equal(category, daily["category"])


@pytest.fixture(scope="module", params=[False, True], ids=["heatwaves", "cold-spells"])
def case(request, tmp_path_factory):
    cold = request.param
    stem = tmp_path_factory.mktemp("array") / ("cold" if cold else "heat")
    clim = hw3.ts2clm3(SST, stem, CP, pctile=10 if cold else 90, quiet=True)
    written = hw3.detect_event3(SST, stem, clim_file=clim, coldSpells=cold, category=True,
                                daily="only", quiet=True)
    sst = _core.hw3_read_sst(str(SST))
    temp = sst["sst"].reshape(sst["nlon"] * sst["nlat"], sst["ntime"])
    seas, thresh = per_day_climatology(clim, sst["time_days"])
    return {"cold": cold, "clim": clim, "temp": temp, "seas": seas, "thresh": thresh,
            "time_days": sst["time_days"], "daily": read_daily(written["daily"])}


def test_full_series_matches_detect_event3(case):
    arrays = hw3.detect_event_array3(case["temp"], case["seas"], case["thresh"],
                                     coldSpells=case["cold"], outputs=ALL)
    assert arrays[0].any() and arrays[3].any()  # the case has events and categories
    assert_same_as_daily(arrays, case["daily"])


def write_window(path, start, stop, missing_days=()):
    """Copy days [start, stop) of sst_test.nc byte for byte; missing_days
    (relative to start) are set to the variable's _FillValue."""
    netCDF4 = pytest.importorskip("netCDF4")
    with netCDF4.Dataset(SST) as src, netCDF4.Dataset(path, "w") as dst:
        src.set_auto_maskandscale(False)
        for name, dim in src.dimensions.items():
            dst.createDimension(name, stop - start if name == "time" else len(dim))
        dst.setncatts({k: src.getncattr(k) for k in src.ncattrs()})
        for name, var in src.variables.items():
            attrs = {k: var.getncattr(k) for k in var.ncattrs()}
            out = dst.createVariable(name, var.dtype, var.dimensions,
                                     fill_value=attrs.pop("_FillValue", None))
            out.set_auto_maskandscale(False)
            out.setncatts(attrs)
            data = var[start:stop] if "time" in var.dimensions else var[:]
            if name == "analysed_sst":
                data = data.copy()
                data[list(missing_days)] = var.getncattr("_FillValue")
            out[:] = data


def longest_event(event):
    """(start, length) of the longest run of event days in any pixel."""
    best = (0, 0)
    for series in event:
        edges = np.flatnonzero(np.diff(np.concatenate([[0], series.astype(np.int8), [0]])))
        for start, stop in zip(edges[::2], edges[1::2]):
            if stop - start > best[1]:
                best = (int(start), int(stop - start))
    return best


@pytest.mark.parametrize("length", [17, 60])
def test_windows_match_detect_event3_on_the_window(case, tmp_path, length):
    start, n = longest_event(case["daily"]["event"])
    assert n >= 10
    # (first day, missing days within the window, whether it must hold an event)
    windows = [
        (max(start - 3, 0), (), True),                   # the event starts inside
        (start + 2, (), True),                           # cut off at its start
        (start + n - 3 - length, (), True),              # cut off at its end
        (start + 2, (length - 3, length - 2), True),     # cut off, with missing days
        (0, (), False),                                  # the first days of the series
    ]
    for k, (first, missing, has_event) in enumerate(windows):
        stop = first + length
        window = tmp_path / f"w{k}.nc"
        write_window(window, first, stop, missing)
        written = hw3.detect_event3(window, tmp_path / f"w{k}", clim_file=case["clim"],
                                    coldSpells=case["cold"], category=True, daily="only",
                                    quiet=True)
        temp = case["temp"][:, first:stop].copy()
        temp[:, list(missing)] = np.nan
        arrays = hw3.detect_event_array3(temp, case["seas"][:, first:stop],
                                         case["thresh"][:, first:stop],
                                         coldSpells=case["cold"], outputs=ALL)
        assert_same_as_daily(arrays, read_daily(written["daily"]))
        if has_event:
            assert arrays[0].any(), f"window {k} has no event"


def test_any_leading_shape_and_input_types(case):
    temp, seas, thresh = case["temp"], case["seas"], case["thresh"]
    flat = hw3.detect_event_array3(temp, seas, thresh, outputs=ALL)
    grid = hw3.detect_event_array3(temp.reshape(2, 3, -1), seas.reshape(2, 3, -1),
                                   thresh.reshape(2, 3, -1), outputs=ALL)
    one = hw3.detect_event_array3(temp[4], seas[4], thresh[4], outputs=ALL)
    strided = hw3.detect_event_array3(temp[::2], seas[::2], thresh[::2], outputs=ALL)
    for a, g, o, s in zip(flat, grid, one, strided):
        np.testing.assert_array_equal(g.reshape(a.shape), a)
        np.testing.assert_array_equal(o, a[4])
        np.testing.assert_array_equal(s, a[::2])
    # float32 input is compared in double, like the same values given as float64
    t32 = temp.astype(np.float32)
    for a, b in zip(hw3.detect_event_array3(t32, seas, thresh, outputs=ALL),
                    hw3.detect_event_array3(t32.astype(np.float64), seas, thresh, outputs=ALL)):
        np.testing.assert_array_equal(a, b)


def test_outputs_selection_order_and_dtypes(case):
    temp, seas, thresh = case["temp"], case["seas"], case["thresh"]
    full = dict(zip(ALL, hw3.detect_event_array3(temp, seas, thresh, outputs=ALL)))
    assert {k: v.dtype for k, v in full.items()} == {
        "event": np.bool_, "event_no": np.int32, "intensity": np.float32, "category": np.int8}
    default = hw3.detect_event_array3(temp, seas, thresh)
    assert len(default) == 3
    for got, name in zip(default, ("event", "intensity", "category")):
        np.testing.assert_array_equal(got, full[name])
    single = hw3.detect_event_array3(temp, seas, thresh, outputs="intensity")
    np.testing.assert_array_equal(single, full["intensity"])
    swapped = hw3.detect_event_array3(temp, seas, thresh, outputs=["category", "event"])
    np.testing.assert_array_equal(swapped[0], full["category"])
    np.testing.assert_array_equal(swapped[1], full["event"])


def test_pixels_without_data_or_climatology_get_the_defaults():
    rng = np.random.default_rng(0)
    temp = rng.uniform(19.0, 20.4, (3, 30))  # below the threshold, except where set
    seas, thresh = np.full_like(temp, 20.0), np.full_like(temp, 20.5)
    temp[0, :] = np.nan
    seas[1, :] = np.nan
    temp[2, 5:15] = 25.0
    event, event_no, intensity, category = hw3.detect_event_array3(temp, seas, thresh, outputs=ALL)
    for px in (0, 1):
        assert not event[px].any() and not event_no[px].any() and not category[px].any()
        assert np.isnan(intensity[px]).all()
    assert event[2, 5:15].all() and (event_no[2, 5:15] == 1).all()
    np.testing.assert_array_equal(intensity[2], np.round(temp[2] - 20.0, 4).astype(np.float32))


def test_parameters_reach_the_detection():
    temp = np.full((1, 20), 20.0)
    temp[0, 2:6] = 25.0     # 4 days
    temp[0, 8:14] = 25.0    # 6 days, after a 2-day gap
    seas, thresh = np.full_like(temp, 20.0), np.full_like(temp, 21.0)
    event = hw3.detect_event_array3(temp, seas, thresh, outputs="event")[0]
    assert event[8:14].all() and not event[2:6].any()     # the short run alone is no event
    event = hw3.detect_event_array3(temp, seas, thresh, minDuration=4, outputs="event")[0]
    assert event[2:14].all()                               # both qualify, gap joined
    event = hw3.detect_event_array3(temp, seas, thresh, minDuration=4, maxGap=1,
                                    outputs="event")[0]
    assert event[2:6].all() and not event[6:8].any() and event[8:14].all()
    cold = hw3.detect_event_array3(-temp, -seas, -thresh, coldSpells=True, outputs="event")[0]
    np.testing.assert_array_equal(cold, hw3.detect_event_array3(temp, seas, thresh,
                                                                outputs="event")[0])


@pytest.mark.parametrize("bad, message", [
    ({"outputs": ()}, "'outputs' must name"),
    ({"outputs": ("event", "event")}, "'outputs' must name"),
    ({"outputs": ("events",)}, "'outputs' must name"),
    ({"minDuration": 0}, "minDuration must be at least 1"),
    ({"maxGap": -1}, "maxGap must not be negative"),
])
def test_invalid_arguments(bad, message):
    x = np.zeros((2, 10))
    with pytest.raises((ValueError, RuntimeError), match=message):
        hw3.detect_event_array3(x, x, x, **bad)


def test_shapes_must_match():
    x = np.zeros((2, 10))
    with pytest.raises(ValueError, match="same shape"):
        hw3.detect_event_array3(x, x[:, :9], x)
    with pytest.raises(ValueError, match="same shape"):
        hw3.detect_event_array3(x, x, x.ravel())


def test_threads_and_concurrent_calls_agree(case):
    from concurrent.futures import ThreadPoolExecutor
    temp, seas, thresh = case["temp"], case["seas"], case["thresh"]
    ref = hw3.detect_event_array3(temp, seas, thresh, outputs=ALL, n_threads=1)
    multi = hw3.detect_event_array3(temp, seas, thresh, outputs=ALL, n_threads=4)
    with ThreadPoolExecutor(4) as pool:
        parallel = list(pool.map(
            lambda px: hw3.detect_event_array3(temp[px], seas[px], thresh[px], outputs=ALL,
                                               n_threads=1), range(temp.shape[0])))
    for k in range(len(ALL)):
        np.testing.assert_array_equal(multi[k], ref[k])
        np.testing.assert_array_equal(np.stack([p[k] for p in parallel]), ref[k])


def test_xarray_apply_ufunc_with_dask_matches_one_call(case):
    xr = pytest.importorskip("xarray")
    pytest.importorskip("dask")
    temp, seas, thresh = case["temp"], case["seas"], case["thresh"]
    names = ("event", "intensity", "category")
    ref = hw3.detect_event_array3(temp, seas, thresh)

    def as_da(x):
        return xr.DataArray(x, dims=("pixel", "time")).chunk({"pixel": 2, "time": -1})

    out = xr.apply_ufunc(hw3.detect_event_array3, as_da(temp), as_da(seas), as_da(thresh),
                         input_core_dims=[["time"]] * 3, output_core_dims=[["time"]] * 3,
                         kwargs={"n_threads": 1}, dask="parallelized",
                         output_dtypes=[bool, np.float32, np.int8])
    for got, expected, name in zip(out, ref, names):
        assert got.chunks[0] == (2, 2, 2), name
        np.testing.assert_array_equal(got.values, expected, err_msg=name)
