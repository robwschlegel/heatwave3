"""Zarr input (Python only; decision 9 in CLAUDE.md).

heatwave3 reads its input in C++ through libnetcdf. For a Zarr store, local or
remote, the part of the store a call needs is copied to a temporary NetCDF next
to the outputs, in the store's own storage type and packing (packed integers
keep their scale_factor/add_offset/_FillValue, float32 stays float32), and the
C++ reads that file as usual. The copy is the requested region padded by one
cell on each side, so the C++ makes its normal selection from it and the
results are those of a NetCDF holding the same data.

The vertical axis is ordered from the surface down and stored as positive
depths in metres, whatever the store's convention (Copernicus Marine ARCO
stores hold negative elevations, deepest first), so depth level 0 is always
the surface.

A climatology may also be a Zarr store holding a heatwave3 climatology
(``seas``/``thresh`` on lon, lat, doy and optionally depth). It is copied
(the requested region plus one cell, or whole) to a temporary NetCDF in the
dimension order the C++ reads.
"""

from __future__ import annotations

import datetime as dt
import os
import uuid
import warnings
from contextlib import contextmanager

import numpy as np

URL_SCHEMES = ("http://", "https://", "s3://", "gs://", "gcs://", "az://", "abfs://")

# Coordinate tolerance of the C++ reader (find_subset_indices() in netcdf_io.cpp)
EPS = 1e-4

# Variable auto-detection, as detect_sst_variable() in netcdf_io.cpp
SST_NAMES = ("analysed_sst", "sst", "SST", "sea_surface_temperature", "temperature",
             "temp", "thetao", "tos", "t_an", "sst_mean")
SST_STANDARD_NAMES = ("sea_surface_temperature", "sea_surface_foundation_temperature",
                      "sea_surface_skin_temperature", "sea_water_temperature")


def is_zarr(path) -> bool:
    """True for a path or URL ending in ``.zarr``, or a local Zarr directory."""
    if not isinstance(path, (str, os.PathLike)):
        return False
    p = os.fspath(path).rstrip("/")
    if p.endswith(".zarr"):
        return True
    return os.path.isdir(p) and any(os.path.exists(os.path.join(p, m))
                                    for m in (".zgroup", ".zmetadata", "zarr.json"))


def _xarray():
    try:
        import cftime  # noqa: F401
        import netCDF4  # noqa: F401
        import xarray as xr
        import zarr  # noqa: F401
    except ImportError as e:
        raise ImportError("Reading a Zarr store needs the optional dependencies: "
                          "pip install 'heatwave3[zarr]'") from e
    return xr


def _is_403(e) -> bool:
    # aiohttp (http/https) carries the status; s3fs raises PermissionError
    return getattr(e, "status", None) == 403 or isinstance(e, PermissionError)


def _store(path: str):
    if not path.startswith(URL_SCHEMES):
        return path
    from zarr.storage import FsspecStore

    class MissingAs403(FsspecStore):
        # A bucket that forbids listing answers 403 for a key that does not
        # exist, and chunks holding only fill values (e.g. all land) are often
        # not stored. Treat 403 alone as a missing chunk (fill value); any
        # other error still raises, so a server fault never becomes NaN.
        async def get(self, key, prototype, byte_range=None):
            try:
                return await super().get(key, prototype, byte_range)
            except Exception as e:
                if _is_403(e):
                    return None
                raise

    return MissingAs403.from_url(path, read_only=True)


def _role(name: str, attrs) -> str | None:
    # As classify_coord() in read_sst_netcdf(), plus "elevation"
    ax, sn, un = (str(attrs.get(k, "")) for k in ("axis", "standard_name", "units"))
    if (ax == "X" or sn == "longitude" or un in ("degrees_east", "degree_east") or
            name in ("lon", "longitude", "x", "xt_ocean", "nav_lon", "geolon", "ni", "xh")):
        return "lon"
    if (ax == "Y" or sn == "latitude" or un in ("degrees_north", "degree_north") or
            name in ("lat", "latitude", "y", "yt_ocean", "nav_lat", "geolat", "nj", "yh")):
        return "lat"
    if (ax == "T" or sn == "time" or "since" in un or
            name in ("time", "t", "time_counter")):
        return "time"
    if (ax == "Z" or sn in ("depth", "altitude", "height") or
            name in ("depth", "level", "lev", "z", "zlev", "elevation")):
        return "depth"
    return None


def _axes(ds, var: str, path: str) -> dict[str, str]:
    """Map role (lon, lat, time, depth) to the variable's dimension name."""
    roles = {}
    for d in ds[var].dims:
        cands = [d] if d in ds.variables else []
        cands += [n for n, v in ds.variables.items() if v.dims == (d,) and n != d]
        for c in cands:
            r = _role(c, ds[c].attrs)
            if r and r not in roles:
                roles[r] = d
                break
    if not {"lon", "lat", "time"} <= roles.keys():
        raise RuntimeError(f"Cannot identify lon/lat/time dimensions in {path}. Ensure the "
                           "store has coordinate variables with CF-compliant axis, "
                           "standard_name, or units attributes.")
    return roles


def _find_var(ds, path: str) -> str:
    for name, v in ds.variables.items():
        sn = str(v.attrs.get("standard_name", ""))
        if name in SST_NAMES or any(s in sn for s in SST_STANDARD_NAMES):
            return name
    raise RuntimeError(f"Could not auto-detect SST variable in {path}")


def _window(coord, rng) -> slice:
    """Index slice covering ``rng`` (with the C++ tolerance) plus one cell on
    each side; the whole axis without a range."""
    n = len(coord)
    if rng is None:
        return slice(0, n)
    inside = np.flatnonzero((coord >= rng[0] - EPS) & (coord <= rng[1] + EPS))
    if inside.size == 0:
        raise RuntimeError("No coordinates found within specified range")
    return slice(max(int(inside[0]) - 1, 0), min(int(inside[-1]) + 2, n))


def _time_window(tvar, time_range) -> slice:
    """Time steps covering ``time_range`` (dates) plus one day and one step on
    each side, compared in the store's own time units."""
    n = tvar.shape[0]
    if time_range is None:
        return slice(0, n)
    import cftime
    units = tvar.attrs["units"]
    cal = str(tvar.attrs.get("calendar", "standard"))
    day = dt.timedelta(days=1)
    lo, hi = (cftime.datetime(*map(int, str(d)[:10].split("-")), calendar=cal)
              for d in time_range)
    bounds = cftime.date2num([lo - day, hi + day], units, calendar=cal)
    return _window(np.asarray(tvar.values, dtype=float), (bounds[0], bounds[1]))


def _levels(values, depth, depth_range):
    """Stored indices of the levels to copy (surface first) and their depths
    (positive metres). ``depth`` is a level index counted from the surface."""
    z = np.asarray(values, dtype=float)
    z = -z if np.all(z <= 0) else z  # elevation (negative up) -> depth
    order = np.argsort(np.abs(z), kind="stable")
    if depth is not None:
        if not 0 <= depth < len(order):
            raise ValueError(f"depth index {depth} is out of range: the store has "
                             f"{len(order)} levels (0 = surface).")
        keep = order[depth:depth + 1]
    elif depth_range is not None:
        keep = order[_window(z[order], depth_range)]
    else:
        keep = order[:1]  # the C++ default: level 0
    return keep, z[keep]


def _open(path: str):
    """The store as stored: lazy, packed values and raw time."""
    xr = _xarray()
    with warnings.catch_warnings():
        # Consolidated metadata is used when present; its absence (usual for
        # Zarr v3) only costs a few extra metadata reads
        warnings.filterwarnings("ignore", "Failed to open Zarr store with consolidated")
        return xr, xr.open_zarr(_store(path), consolidated=None, chunks=None,
                                mask_and_scale=False, decode_times=False)


def _to_netcdf(sub, tmp_dir, fills):
    """Write ``sub`` to a temporary NetCDF in ``tmp_dir``, storage types as
    they are, ``fills`` the _FillValue of each data variable."""
    for name in sub.variables:
        sub[name].attrs.pop("_FillValue", None)
        sub[name].encoding = {}
    enc = {name: {"_FillValue": None} for name in sub.coords}
    enc.update({n: {"_FillValue": f, "dtype": sub[n].dtype} for n, f in fills.items()})
    tmp = os.path.join(tmp_dir, f".hw3_zarr_{os.getpid()}_{uuid.uuid4().hex[:8]}.nc")
    try:
        sub.to_netcdf(tmp, engine="netcdf4", format="NETCDF4", encoding=enc)
    except BaseException:
        _remove(tmp)
        raise
    return tmp


def _fill(ds, name):
    return ds[name].attrs.get("_FillValue", ds[name].encoding.get("_FillValue"))


def _write_subset(path, tmp_dir, var_name, lon_range, lat_range, time_range, depth,
                  depth_range, quiet):
    xr, ds = _open(path)
    var = var_name or _find_var(ds, path)
    if var not in ds.data_vars:
        raise RuntimeError(f"Variable '{var}' not found in {path}")
    ax = _axes(ds, var, path)
    sel = {ax["lon"]: _window(np.asarray(ds[ax["lon"]].values, dtype=float), lon_range),
           ax["lat"]: _window(np.asarray(ds[ax["lat"]].values, dtype=float), lat_range),
           ax["time"]: _time_window(ds[ax["time"]], time_range)}
    zdim = ax.get("depth")
    if zdim is not None:
        keep, zvals = _levels(ds[zdim].values, depth, depth_range)
        sel[zdim] = keep
    v = ds[var].isel(sel)
    # 1-D variables along the copied dimensions carry the coordinates
    coords = {n: ds[n].isel({ds[n].dims[0]: sel[ds[n].dims[0]]}).variable
              for n in ds.variables if ds[n].ndim == 1 and ds[n].dims[0] in sel}
    sub = xr.Dataset({var: v.variable}, coords=coords)
    if zdim is not None:
        sub = sub.drop_vars(zdim)
        if zdim != "depth":
            sub = sub.rename({zdim: "depth"})
        sub = sub.assign_coords(depth=("depth", zvals, {
            "units": "m", "positive": "down", "axis": "Z", "standard_name": "depth"}))
    sub.attrs = dict(ds.attrs)
    if not quiet:
        print(f"Reading {var} from Zarr store {path}...")
        if zdim is not None:
            print(f"Depth levels from the surface down: {', '.join(f'{z:g}' for z in zvals)} m")
    # Storage type and packing exactly as stored
    return _to_netcdf(sub, tmp_dir, {var: _fill(ds, var)}), var


# Climatology files: the layout write_clim_netcdf() writes and
# read_clim_netcdf() reads as one block, so the order matters
CLIM_DIMS = ("lon", "lat", "depth", "doy")


def _write_clim(path, tmp_dir, lon_range=None, lat_range=None):
    _, ds = _open(path)
    bad = f"{path} is not a heatwave3 climatology"
    missing = [v for v in ("seas", "thresh") if v not in ds.data_vars]
    if missing:
        raise ValueError(f"{bad}: no {' or '.join(missing)} variable.")
    dims = ds["seas"].dims
    order = tuple(d for d in CLIM_DIMS if d in dims)
    if (set(dims) != set(order) or not {"lon", "lat", "doy"} <= set(dims) or
            ds.sizes["doy"] != 366 or ds["thresh"].dims != dims):
        raise ValueError(f"{bad}: seas and thresh must have dimensions lon, lat, doy (366) "
                         f"and optionally depth; seas has {dims}.")
    names = [v for v in ("seas", "thresh", "var") if v in ds.data_vars]
    sub = ds[names].drop_vars([c for c in ds.coords if c not in CLIM_DIMS])
    # The region plus one cell on each side; the C++ then selects it exactly
    sub = sub.isel(lon=_window(np.asarray(ds["lon"].values, dtype=float), lon_range),
                   lat=_window(np.asarray(ds["lat"].values, dtype=float), lat_range))
    sub = sub.transpose(*order)
    sub.attrs = dict(ds.attrs)
    return _to_netcdf(sub, tmp_dir, {n: _fill(ds, n) for n in names})


def _remove(path: str) -> None:
    try:
        os.remove(path)
    except FileNotFoundError:
        pass


@contextmanager
def nc_input(file_in, tmp_dir, var_name=None, lon_range=None, lat_range=None,
             time_range=None, depth=None, depth_range=None, quiet=True):
    """Yield ``(file, var_name, depth)`` for the C++ reader.

    NetCDF input passes through unchanged. A Zarr store is copied (the region
    plus one cell on each side, in its own packing) to a temporary NetCDF in
    ``tmp_dir``, removed on exit; a ``depth`` level index then becomes 0, the
    one level copied. Ranges are (min, max) pairs, ``time_range`` dates.
    """
    if not is_zarr(file_in):
        yield file_in, var_name, depth
        return
    tmp, var = _write_subset(os.fspath(file_in), tmp_dir, var_name, lon_range, lat_range,
                             time_range, depth, depth_range, quiet)
    try:
        yield tmp, var, (None if depth is None else 0)
    finally:
        _remove(tmp)


@contextmanager
def clim_input(clim_file, tmp_dir, lon_range=None, lat_range=None):
    """Yield the climatology path for the C++ reader: a NetCDF unchanged, a
    Zarr store (a heatwave3 climatology) as a temporary NetCDF in
    ``tmp_dir`` (only the lon/lat ranges given, plus one cell on each side),
    removed on exit."""
    if not is_zarr(clim_file):
        yield clim_file
        return
    tmp = _write_clim(os.fspath(clim_file), tmp_dir, lon_range, lat_range)
    try:
        yield tmp
    finally:
        _remove(tmp)


def local_dir(*paths) -> str:
    """Directory for temporary files: that of the first local path given,
    else the working directory."""
    for p in paths:
        if p and isinstance(p, (str, os.PathLike)) and not os.fspath(p).startswith(URL_SCHEMES):
            return os.path.dirname(os.path.abspath(os.fspath(p)))
    return os.getcwd()


def check_clim(clim_file) -> None:
    if not is_zarr(clim_file) and not os.path.exists(clim_file):
        raise FileNotFoundError(f"Climatology file does not exist: {clim_file}")


def set_source(files, source) -> None:
    """Record the Zarr store, not the temporary copy, as ``source_file`` of
    the written outputs."""
    if not is_zarr(source):
        return
    import netCDF4
    for f in files:
        if f and os.path.exists(f):
            with netCDF4.Dataset(f, "a") as nc:
                nc.setncattr("source_file", os.fspath(source))


def clim_extent(clim_file) -> dict:
    """``lon_range``, ``lat_range`` and ``depth_range`` (``None`` without depth)
    of a climatology file: the region an SST read is matched to."""
    _xarray()
    import netCDF4
    with netCDF4.Dataset(os.fspath(clim_file)) as nc:
        lon = np.asarray(nc["lon"][:], dtype=float)
        lat = np.asarray(nc["lat"][:], dtype=float)
        if "depth" in nc.variables:
            z = np.asarray(nc["depth"][:], dtype=float)
        elif "depth_value" in nc.ncattrs():
            z = np.array([float(nc.getncattr("depth_value"))])
        else:
            z = None
    return {"lon_range": [lon.min(), lon.max()], "lat_range": [lat.min(), lat.max()],
            "depth_range": None if z is None else [z.min(), z.max()]}
