"""Zarr input (Python only): every SST-reading function gives the same result
from a Zarr store as from a NetCDF holding the same data."""

import functools
import http.server
import os
import threading
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

xr = pytest.importorskip("xarray")
pytest.importorskip("zarr")
pytest.importorskip("netCDF4")
pytest.importorskip("cftime")

import heatwave3 as hw3  # noqa: E402
from heatwave3 import _zarr  # noqa: E402

EXT = Path(__file__).resolve().parents[2] / "inst" / "extdata"
SST = EXT / "sst_test.nc"
DEPTH = EXT / "glorys_depth_test.nc"
CP = ("1982-01-01", "2011-12-31")
CP_DEPTH = ("1993-01-01", "2022-12-31")

pytestmark = pytest.mark.skipif(not SST.exists() or not DEPTH.exists(),
                                reason="bundled test data not available")


def raw(path):
    """A NetCDF as stored: packed integers, raw time, no Zarr encodings."""
    ds = xr.open_dataset(path, mask_and_scale=False, decode_times=False)
    for v in ds.variables:
        ds[v].encoding = {}
    return ds


def no_temp_files(d):
    return not [f for f in os.listdir(d) if f.startswith(".hw3_zarr_")]


@pytest.fixture(scope="module")
def stores(tmp_path_factory):
    d = tmp_path_factory.mktemp("zarr")
    s = {}
    ds = raw(SST)
    for fmt in (2, 3):
        s[f"v{fmt}"] = str(d / f"sst_v{fmt}.zarr")
        ds.to_zarr(s[f"v{fmt}"], zarr_format=fmt, consolidated=fmt == 2, mode="w")
    # float32 data, in both formats
    f32 = xr.open_dataset(SST).astype({"analysed_sst": "float32"})
    s["f32_nc"] = str(d / "sst_f32.nc")
    f32.to_netcdf(s["f32_nc"])
    f32 = f32.copy()
    for v in f32.variables:
        f32[v].encoding = {}
    s["f32"] = str(d / "sst_f32.zarr")
    f32.to_zarr(s["f32"], consolidated=False, mode="w")
    # Copernicus ARCO layout: negative elevations, deepest first
    g = raw(DEPTH)
    g = g.isel(depth=slice(None, None, -1)).rename(depth="elevation")
    g = g.assign_coords(elevation=("elevation", -g["elevation"].values, g["elevation"].attrs))
    s["arco"] = str(d / "glorys_arco.zarr")
    g.to_zarr(s["arco"], zarr_format=2, mode="w")
    return s


def same_product(a, b):
    pd.testing.assert_frame_equal(hw3.hw3_export(a), hw3.hw3_export(b))


@pytest.mark.parametrize("store", ["v2", "v3"])
def test_detect3_equals_netcdf(stores, tmp_path, store):
    ref = hw3.detect3(SST, tmp_path / "nc", CP, category=True, daily="also", quiet=True)
    got = hw3.detect3(stores[store], tmp_path / "z", CP, category=True, daily="also",
                      quiet=True)
    assert list(got) == list(ref)
    for k in ref:
        same_product(ref[k], got[k])
    assert no_temp_files(tmp_path)


def test_float32_store_equals_float32_netcdf(stores, tmp_path):
    ref = hw3.detect3(stores["f32_nc"], tmp_path / "nc", CP, quiet=True)
    got = hw3.detect3(stores["f32"], tmp_path / "z", CP, quiet=True)
    for k in ref:
        same_product(ref[k], got[k])


def test_subset_equals_netcdf(stores, tmp_path):
    kw = dict(lon_range=(26.5, 26.53), lat_range=(-34.1, -34.0),
              time_range=("1982-01-01", "2015-12-31"))
    same_product(hw3.ts2clm3(SST, tmp_path / "nc", CP, quiet=True, **kw),
                 hw3.ts2clm3(stores["v2"], tmp_path / "z", CP, quiet=True, **kw))


def test_detect_event3_reads_the_climatology_extent(stores, tmp_path):
    clim = hw3.ts2clm3(SST, tmp_path / "c", CP, lon_range=(26.5, 26.53), quiet=True)
    a = hw3.detect_event3(SST, tmp_path / "nc", clim_file=clim, protoEvent=True, quiet=True)
    b = hw3.detect_event3(stores["v3"], tmp_path / "z", clim_file=clim, protoEvent=True,
                          quiet=True)
    same_product(a["protoevents"], b["protoevents"])


def test_threshclim2_from_zarr(stores, tmp_path):
    clim = hw3.ts2clm3(SST, tmp_path / "c", CP, quiet=True)
    a = hw3.detect_event3(SST, tmp_path / "nc", clim_file=clim, threshClim2=SST,
                          quiet=True)
    b = hw3.detect_event3(SST, tmp_path / "z", clim_file=clim, threshClim2=stores["v2"],
                          quiet=True)
    same_product(a["events"], b["events"])


def test_exceedance3_equals_netcdf(stores, tmp_path):
    a = hw3.exceedance3(SST, tmp_path / "nc.nc", 295.0, quiet=True)
    b = hw3.exceedance3(stores["v2"], tmp_path / "z.nc", 295.0, quiet=True)
    same_product(a, b)


def test_category_functions_equal_netcdf(stores, tmp_path):
    out = hw3.detect3(SST, tmp_path / "s", CP, quiet=True)
    pd.testing.assert_frame_equal(
        hw3.category3(out["events"], out["clim"], sst_file=SST),
        hw3.category3(out["events"], out["clim"], sst_file=stores["v2"]))
    tr = ("2010-01-01", "2010-03-31")
    pd.testing.assert_frame_equal(
        hw3.category_daily3(SST, out["clim"], out["events"], tr),
        hw3.category_daily3(stores["v3"], out["clim"], out["events"], tr))
    assert no_temp_files(tmp_path)


def test_detect_blob3_equals_netcdf(stores, tmp_path):
    clim = hw3.ts2clm3(SST, tmp_path / "c", CP, quiet=True)
    a = hw3.detect_blob3(SST, clim, return_=["event", "daily"])
    b = hw3.detect_blob3(stores["v2"], clim, return_=["event", "daily"])
    for k in a:
        pd.testing.assert_frame_equal(a[k], b[k])


def test_event_line3_draws_the_same_data(stores, tmp_path):
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    out = hw3.detect3(SST, tmp_path / "s", CP, quiet=True)
    kw = dict(lon=26.525, lat=-34.125, event_file=out["events"])
    a = hw3.event_line3(SST, out["clim"], **kw)
    b = hw3.event_line3(stores["v2"], out["clim"], **kw)
    for la, lb in zip(a.lines, b.lines):
        np.testing.assert_array_equal(la.get_xydata(), lb.get_xydata())
    with pytest.raises(ValueError, match="No data found"):
        hw3.event_line3(stores["v2"], out["clim"], lon=0, lat=0)


def test_depth_levels_count_from_the_surface(stores, tmp_path):
    # The ARCO store holds the same levels negated and deepest first
    for depth in (0, 3):
        same_product(hw3.ts2clm3(DEPTH, tmp_path / f"nc{depth}", CP_DEPTH, depth=depth,
                                 quiet=True),
                     hw3.ts2clm3(stores["arco"], tmp_path / f"z{depth}", CP_DEPTH,
                                 depth=depth, quiet=True))
    ref = hw3.detect3(DEPTH, tmp_path / "ncr", CP_DEPTH, depth_range=(0, 20), quiet=True)
    got = hw3.detect3(stores["arco"], tmp_path / "zr", CP_DEPTH, depth_range=(0, 20),
                      quiet=True)
    for k in ref:
        same_product(ref[k], got[k])
    with pytest.raises(ValueError, match="out of range"):
        hw3.ts2clm3(stores["arco"], tmp_path / "bad", CP_DEPTH, depth=99, quiet=True)


def test_source_file_names_the_store(stores, tmp_path):
    import netCDF4
    out = hw3.detect3(stores["v2"], tmp_path / "z", CP, quiet=True)
    for f in out.values():
        with netCDF4.Dataset(f) as nc:
            assert nc.source_file == stores["v2"]


def test_temp_file_removed_on_error(stores, tmp_path, monkeypatch):
    with pytest.raises(RuntimeError, match="not found"):
        hw3.ts2clm3(stores["v2"], tmp_path / "z", CP, var_name="nope", quiet=True)

    def fail(*a, **k):
        raise RuntimeError("boom")

    monkeypatch.setattr(hw3._core, "hw3_compute_clim", fail)
    with pytest.raises(RuntimeError, match="boom"):
        hw3.ts2clm3(stores["v2"], tmp_path / "z", CP, quiet=True)
    assert no_temp_files(tmp_path)


def test_missing_dependency_message(monkeypatch, tmp_path):
    import builtins
    real = builtins.__import__

    def fake(name, *a, **k):
        if name == "zarr":
            raise ImportError("no zarr")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake)
    with pytest.raises(ImportError, match=r"heatwave3\[zarr\]"):
        hw3.ts2clm3(tmp_path / "x.zarr", tmp_path / "z", CP, quiet=True)


def test_is_zarr(tmp_path, stores):
    assert _zarr.is_zarr("https://example.org/data/geoChunked.zarr")
    assert _zarr.is_zarr(stores["v3"] + "/")
    assert not _zarr.is_zarr(SST)
    assert not _zarr.is_zarr(EXT)  # a directory of NetCDF files
    assert not _zarr.is_zarr([str(SST)])


# Climatologies stored as Zarr

def clim_to_zarr(clim_nc, store, transpose=None):
    ds = xr.open_dataset(clim_nc)
    if transpose:
        ds = ds.transpose(*transpose)
    for v in ds.variables:
        ds[v].encoding = {}
    ds.to_zarr(store, consolidated=False, mode="w")
    return str(store)


@pytest.fixture(scope="module")
def clims(tmp_path_factory):
    d = tmp_path_factory.mktemp("clim")
    c = {"nc": hw3.ts2clm3(SST, d / "s", CP, quiet=True)}
    c["zarr"] = clim_to_zarr(c["nc"], d / "s_clim.zarr")
    c["zarr_t"] = clim_to_zarr(c["nc"], d / "t_clim.zarr", ("doy", "lat", "lon"))
    c["ev"] = hw3.detect_event3(SST, d / "s", clim_file=c["nc"], quiet=True)["events"]
    return c


@pytest.mark.parametrize("which", ["zarr", "zarr_t"])
def test_detect_event3_with_zarr_climatology(clims, stores, tmp_path, which):
    kw = dict(category=True, daily="also", quiet=True)
    ref = hw3.detect_event3(SST, tmp_path / "nc", clim_file=clims["nc"], **kw)
    for sst in (SST, stores["v2"]):
        got = hw3.detect_event3(sst, tmp_path / "z", clim_file=clims[which], **kw)
        for k in ref:
            same_product(ref[k], got[k])
    assert no_temp_files(tmp_path)


def test_other_functions_with_zarr_climatology(clims, stores, tmp_path):
    ev, nc, z = clims["ev"], clims["nc"], clims["zarr"]
    pd.testing.assert_frame_equal(hw3.category3(ev, nc, sst_file=SST),
                                  hw3.category3(ev, z, sst_file=stores["v2"]))
    tr = ("2010-01-01", "2010-03-31")
    pd.testing.assert_frame_equal(hw3.category_daily3(SST, nc, ev, tr),
                                  hw3.category_daily3(stores["v3"], z, ev, tr))
    a = hw3.detect_blob3(SST, nc, return_=["event", "daily"])
    b = hw3.detect_blob3(stores["v2"], z, return_=["event", "daily"])
    for k in a:
        pd.testing.assert_frame_equal(a[k], b[k])
    pytest.importorskip("matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    kw = dict(lon=26.525, lat=-34.125, event_file=ev)
    for la, lb in zip(hw3.event_line3(SST, nc, **kw).lines,
                      hw3.event_line3(stores["v2"], z, **kw).lines):
        np.testing.assert_array_equal(la.get_xydata(), lb.get_xydata())
    assert no_temp_files(os.path.dirname(ev))


def test_depth_climatology_as_zarr(stores, tmp_path):
    # Depth-resolved (depth dimension) and single-level (depth_value attribute)
    for name, kw in (("range", dict(depth_range=(0, 20))), ("level", dict(depth=2))):
        clim = hw3.ts2clm3(DEPTH, tmp_path / name, CP_DEPTH, quiet=True, **kw)
        zclim = clim_to_zarr(clim, tmp_path / f"{name}_clim.zarr")
        a = hw3.detect_event3(DEPTH, tmp_path / f"{name}nc", clim_file=clim, quiet=True)
        b = hw3.detect_event3(stores["arco"], tmp_path / f"{name}z", clim_file=zclim,
                              quiet=True)
        same_product(a["events"], b["events"])


def test_region_of_a_zarr_climatology(clims, stores, tmp_path):
    # A climatology for the whole grid serving a regional run equals one made
    # for the region alone, from NetCDF or Zarr inputs
    reg = dict(lon_range=(26.5, 26.53), lat_range=(-34.1, -34.0))
    clim_r = hw3.ts2clm3(SST, tmp_path / "r", CP, quiet=True, **reg)
    ref = hw3.detect_event3(SST, tmp_path / "r", clim_file=clim_r, category=True,
                            n_threads=1, quiet=True)
    for sst, clim in ((SST, clims["nc"]), (stores["v2"], clims["zarr"]),
                      (stores["v3"], clims["nc"])):
        got = hw3.detect_event3(sst, tmp_path / "g", clim_file=clim, category=True,
                                n_threads=1, quiet=True, **reg)
        same_product(ref["events"], got["events"])
    tr = ("2000-01-01", "2009-12-31")
    a = hw3.detect_event3(SST, tmp_path / "t", clim_file=clims["nc"], time_range=tr,
                          n_threads=1, quiet=True)
    b = hw3.detect_event3(stores["v2"], tmp_path / "tz", clim_file=clims["zarr"],
                          time_range=tr, n_threads=1, quiet=True)
    same_product(a["events"], b["events"])
    ev = hw3.hw3_export(a["events"])
    assert ev["date_start"].min() >= np.datetime64(tr[0])
    assert ev["date_end"].max() <= np.datetime64(tr[1])
    with pytest.raises(ValueError, match="single-file input only"):
        hw3.detect_event3([SST, SST], tmp_path / "m", clim_file=clims["nc"], time_range=tr,
                          quiet=True)


def test_transposed_netcdf_climatology_is_refused(clims, tmp_path):
    ds = xr.open_dataset(clims["nc"]).transpose("doy", "lat", "lon").load()
    ds.to_netcdf(tmp_path / "t_clim.nc")
    with pytest.raises(RuntimeError, match=r"stored as \(doy, lat, lon\)"):
        hw3.detect_event3(SST, tmp_path / "x", clim_file=tmp_path / "t_clim.nc", quiet=True)


def test_not_a_climatology(stores, tmp_path):
    with pytest.raises(ValueError, match="not a heatwave3 climatology"):
        hw3.detect_event3(SST, tmp_path / "x", clim_file=stores["v2"], quiet=True)
    assert no_temp_files(tmp_path)


# Remote stores: a local HTTP server that, like an S3 bucket that forbids
# listing, answers 403 for any key that does not exist.

class _Handler(http.server.SimpleHTTPRequestHandler):
    fail = None  # path answered with 500

    def log_message(self, *a):
        pass

    def send_head(self):
        if self.fail and self.path.endswith(self.fail):
            self.send_error(500)
            return None
        if not os.path.exists(self.translate_path(self.path)):
            self.send_error(403)
            return None
        return super().send_head()


@pytest.fixture
def server(tmp_path):
    # One pixel is all fill (land), and its chunk is therefore not stored
    ds = raw(SST)
    sst = ds["analysed_sst"].values.copy()
    sst[:, 0, 0] = ds["analysed_sst"].attrs["_FillValue"]
    ds["analysed_sst"].values = sst
    ds["analysed_sst"].encoding = {"chunks": (14276, 1, 1)}
    ds.to_zarr(tmp_path / "land.zarr", zarr_format=2, mode="w")
    ds.to_netcdf(tmp_path / "land.nc")
    assert not (tmp_path / "land.zarr" / "analysed_sst" / "0.0.0").exists()
    handler = functools.partial(_Handler, directory=str(tmp_path))
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    yield tmp_path, f"http://127.0.0.1:{httpd.server_address[1]}/land.zarr"
    httpd.shutdown()
    _Handler.fail = None


def test_remote_store_403_is_a_missing_chunk(server, tmp_path):
    d, url = server
    ref = hw3.detect3(d / "land.nc", tmp_path / "nc", CP, quiet=True)
    got = hw3.detect3(url, tmp_path / "z", CP, quiet=True)
    for k in ref:
        same_product(ref[k], got[k])


def test_remote_climatology(server, tmp_path):
    d, url = server
    clim = hw3.ts2clm3(d / "land.nc", d / "land", CP, quiet=True)
    clim_to_zarr(clim, d / "land_clim.zarr")
    base = url.rsplit("/", 1)[0]
    a = hw3.detect_event3(d / "land.nc", tmp_path / "nc", clim_file=clim, quiet=True)
    b = hw3.detect_event3(url, tmp_path / "z", clim_file=f"{base}/land_clim.zarr",
                          quiet=True)
    same_product(a["events"], b["events"])
    assert no_temp_files(tmp_path)


def test_remote_store_server_error_raises(server, tmp_path):
    _, url = server
    _Handler.fail = "analysed_sst/0.0.1"
    with pytest.raises(Exception, match="500"):
        hw3.ts2clm3(url, tmp_path / "z", CP, quiet=True)
    assert no_temp_files(tmp_path)


GLORYS = ("https://s3.waw3-1.cloudferro.com/mdl-arco-geo-025/arco/GLOBAL_MULTIYEAR_PHY_001_030/"
          "cmems_mod_glo_phy_my_0.083deg_P1D-m_202311/geoChunked.zarr")


@pytest.mark.skipif(os.environ.get("HW3_NETWORK_TESTS") != "1",
                    reason="set HW3_NETWORK_TESTS=1 to read the public GLORYS store")
def test_glorys_arco_store(tmp_path):
    out = hw3.detect3(GLORYS, tmp_path / "g", CP_DEPTH, var_name="thetao",
                      lon_range=(-5.0, -4.9), lat_range=(45.0, 45.1), category=True,
                      hemisphere="north", quiet=True)
    import netCDF4
    with netCDF4.Dataset(out["clim"]) as nc:
        assert abs(nc.depth_value - 0.494) < 1e-3  # level 0 is the surface
    assert len(hw3.hw3_export(out["events"])) > 0


# ---- Zarr output ----

def same_as_netcdf(nc, store):
    """The store decodes as the NetCDF does, and holds the same raw values,
    types, attributes and chunks."""
    import zarr
    a, b = xr.open_dataset(nc), xr.open_zarr(store)
    assert a.identical(b.load())
    ra, rb = xr.open_dataset(nc, decode_cf=False), xr.open_zarr(store, decode_cf=False)
    assert ra.identical(rb.load())
    import netCDF4
    g = zarr.open_group(store, mode="r")
    assert g.metadata.zarr_format == 3
    with netCDF4.Dataset(nc) as d:
        for name, v in d.variables.items():
            ch = v.chunking()
            assert g[name].chunks == (v.shape if ch == "contiguous" else tuple(ch))
            assert g[name].dtype == v.dtype


@pytest.fixture(scope="module")
def products(tmp_path_factory):
    d = tmp_path_factory.mktemp("products")
    out = hw3.detect3(SST, d / "t", CP, category=True, daily="also", quiet=True)
    out["protoevents"] = hw3.detect_event3(SST, d / "p", clim_file=out["clim"],
                                           protoEvent=True, quiet=True)["protoevents"]
    out["depth"] = hw3.detect3(DEPTH, d / "dz", CP_DEPTH, depth_range=(0, 20),
                               daily="also", quiet=True)["daily"]
    return out


@pytest.mark.parametrize("product", ["clim", "events", "daily", "protoevents", "depth"])
def test_hw3_export_to_zarr(products, tmp_path, product):
    store = str(tmp_path / f"{product}.zarr")
    assert hw3.hw3_export(products[product], file_out=store) == store
    same_as_netcdf(products[product], store)
    # Replaces an existing store
    hw3.hw3_export(products[product], file_out=store)
    same_as_netcdf(products[product], store)


def test_blocks_cover_every_chunk(tmp_path):
    import netCDF4
    nc = str(tmp_path / "b.nc")
    with netCDF4.Dataset(nc, "w") as d:
        for n, size in (("lon", 7), ("lat", 5), ("time", 11)):
            d.createDimension(n, size)
        v = d.createVariable("x", "f4", ("lon", "lat", "time"), chunksizes=(3, 2, 4))
        v[:] = np.arange(7 * 5 * 11, dtype="f4").reshape(7, 5, 11)
        v.units = "kelvin"
    store = _zarr.write_zarr(nc, str(tmp_path / "b.zarr"))
    same_as_netcdf(nc, store)


def test_zarr_climatology_output_used_as_clim_file(products, tmp_path):
    store = hw3.hw3_export(products["clim"], file_out=tmp_path / "c.zarr")
    a = hw3.detect_event3(SST, tmp_path / "nc", clim_file=products["clim"], quiet=True)
    b = hw3.detect_event3(SST, tmp_path / "z", clim_file=store, quiet=True)
    same_product(a["events"], b["events"])


def test_zarr_output_to_url(products):
    import zarr
    url = "memory://hw3/events.zarr"
    hw3.hw3_export(products["events"], file_out=url)
    a = xr.open_dataset(products["events"])
    b = xr.open_zarr(zarr.storage.FsspecStore.from_url(url, read_only=True))
    assert a.identical(b.load())


def test_zarr_out_copies_every_product(stores, tmp_path):
    zdir = tmp_path / "zarr"
    zdir.mkdir()
    out = hw3.detect3(stores["v3"], tmp_path / "x", CP, category=True, daily="also",
                      quiet=True, zarr_out=zdir)
    assert sorted(out) == ["clim", "daily", "events"]
    for f in out.values():
        assert f.endswith(".nc") and os.path.exists(f)
        same_as_netcdf(f, str(zdir / (Path(f).stem + ".zarr")))
    # source_file names the input store in the copies too
    assert xr.open_zarr(zdir / "x_clim.zarr").attrs["source_file"] == stores["v3"]
    assert no_temp_files(tmp_path)

    clim = hw3.ts2clm3(SST, tmp_path / "y", CP, quiet=True, zarr_out=zdir)
    same_as_netcdf(clim, str(zdir / "y_clim.zarr"))
    p = hw3.detect_event3(SST, tmp_path / "y", protoEvent=True, quiet=True, zarr_out=zdir)
    same_as_netcdf(p["protoevents"], str(zdir / "y_protoevents.zarr"))


def test_zarr_output_errors(products, tmp_path):
    with pytest.raises(ValueError, match="whole product"):
        hw3.hw3_export(products["daily"], file_out=tmp_path / "d.zarr", n=5)
    with pytest.raises(FileNotFoundError, match="Output directory"):
        hw3.hw3_export(products["clim"], file_out=tmp_path / "no" / "c.zarr")
    # zarr_out is checked before anything is computed
    with pytest.raises(FileNotFoundError, match="zarr_out"):
        hw3.ts2clm3(SST, tmp_path / "e", CP, quiet=True, zarr_out=tmp_path / "missing")
    assert not os.path.exists(tmp_path / "e_clim.nc")
