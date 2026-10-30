"""Read or export a heatwave3 NetCDF product (port of R/export.R)."""

from __future__ import annotations

import os
import warnings

import numpy as np
import pandas as pd

from . import _core
from ._frames import clim_df, daily_df, event_df, subset_df
from ._utils import choose, pair, to_datetime64, to_jd
from ._zarr import is_url, is_zarr, write_zarr

# Row count above which whole-table reads warn and inline previews are skipped.
PREVIEW_LIMIT = 5_000_000

_TYPE_TO_PRODUCT = {"clim": "climatology", "event": "events",
                    "daily": "daily", "protoevents": "protoevents"}
_PRODUCT_TO_CHUNKTYPE = {"climatology": "clim", "events": "event",
                         "daily": "daily", "protoevents": "daily"}


def hw3_export(file, file_out=None, vars=None, lon_range=None, lat_range=None,
               time_range=None, depth_range=None, n=None, type=None,
               chunk_size=1_000_000):
    """Read or export a heatwave3 NetCDF product.

    Reads any heatwave3 NetCDF output (climatology, events, daily series or
    proto-events) with the C++ reader, and either returns it as a long
    ``pandas.DataFrame`` or writes it to a flat companion file. The product
    type is detected automatically.

    Parameters
    ----------
    file : str or path-like
        A heatwave3 NetCDF output file.
    file_out : str or path-like, optional
        Path for a flat companion file ending in ``.csv`` or ``.parquet``
        (the latter needs ``pyarrow``), or a Zarr store ending in ``.zarr``:
        a local path or a URL such as ``s3://bucket/x_clim.zarr`` (needs
        ``pip install 'heatwave3[zarr]'``; ``s3://`` also needs ``s3fs``).
        A Zarr store holds the whole product as stored in the NetCDF (same
        variables, types, attributes and chunks, Zarr v3), replacing any
        store there; it cannot be combined with the subset arguments. When
        given, the product is written there and the path is returned. When
        ``None`` (the default), the product is returned as a DataFrame.
    vars : list of str, optional
        Data variables to read, e.g. ``["seas", "thresh"]``. ``None`` reads all.
    lon_range, lat_range : (min, max), optional
        Only overlapping grid cells are read from disk (a hyperslab read). For
        events these filter by each event's coordinate.
    time_range : (start, end) dates, optional
        Subsets daily and proto-event products; for events, keeps events
        overlapping the range. Ignored (with a warning) for a climatology.
    depth_range : (min, max), optional
        Depth window in metres, for depth-resolved products.
    n : int, optional
        Return only the first ``n`` rows; for gridded products this also caps
        how much is read from disk.
    type : {"clim", "event", "daily", "protoevents"}, optional
        Override the auto-detected product.
    chunk_size : int
        Rows per write chunk for a whole (unsubset) ``.csv``/``.parquet`` export.

    Returns
    -------
    pandas.DataFrame or str
        The product (or requested subset), or ``file_out`` when writing.
    """
    file = os.fspath(file)
    if not os.path.exists(file):
        raise FileNotFoundError(f"File does not exist: {file}")

    if type is None:
        product = _core.hw3_file_meta(file)["product"]
    else:
        product = _TYPE_TO_PRODUCT[choose(type, tuple(_TYPE_TO_PRODUCT), "type")]
    ctype = _product_to_chunktype(product)
    subset_requested = any(x is not None for x in
                           (vars, lon_range, lat_range, time_range, depth_range, n))

    if file_out is None:
        if subset_requested:
            return _read_subset_df(file, product, vars, lon_range, lat_range,
                                   time_range, depth_range, n)
        nrows = _core.hw3_file_meta(file)["nrows"]
        if nrows > PREVIEW_LIMIT:
            warnings.warn(
                f"Reading {nrows:,} rows into memory; pass vars/lon_range/lat_range/"
                "time_range/depth_range/n to read a subset, or file_out to export instead.",
                stacklevel=2)
        return read_product_df(file, product)

    file_out = os.fspath(file_out)
    if is_zarr(file_out):
        if subset_requested:
            raise ValueError("A Zarr file_out holds the whole product: remove "
                             "vars/lon_range/lat_range/time_range/depth_range/n.")
        if not is_url(file_out):
            _validate_out_dir(file_out)
        return write_zarr(file, file_out)
    fmt = _validate_save_file(file_out)

    if subset_requested:
        _write_df(_read_subset_df(file, product, vars, lon_range, lat_range,
                                  time_range, depth_range, n), file_out, fmt)
        return file_out

    chunk_size = int(chunk_size)
    if chunk_size < 1:
        raise ValueError("chunk_size must be a positive integer.")
    if os.path.exists(file_out):
        os.remove(file_out)
    if fmt == "csv":
        first = True
        for chunk in _iterate_chunks(file, ctype, chunk_size):
            chunk.to_csv(file_out, mode="w" if first else "a", header=first, index=False)
            first = False
    else:
        pq = _require_pyarrow()
        import pyarrow as pa
        writer = None
        try:
            for chunk in _iterate_chunks(file, ctype, chunk_size):
                tab = pa.Table.from_pandas(chunk, preserve_index=False)
                if writer is None:
                    writer = pq.ParquetWriter(file_out, tab.schema)
                writer.write_table(tab)
        finally:
            if writer is not None:
                writer.close()
    return file_out


# ---- Readers ----

def read_product_df(file: str, product: str) -> pd.DataFrame:
    if product == "climatology":
        return clim_df(_core.hw3_read_clim_nc(file))
    if product == "events":
        return event_df(_core.hw3_read_event_nc(file))
    if product in ("daily", "protoevents"):
        return daily_df(_core.hw3_read_daily_nc(file))
    raise ValueError(f"Unrecognised product: {product}")


def _product_to_chunktype(product: str) -> str:
    try:
        return _PRODUCT_TO_CHUNKTYPE[product]
    except KeyError:
        raise ValueError(f"Unrecognised product: {product}") from None


def _read_subset_df(file, product, vars=None, lon_range=None, lat_range=None,
                    time_range=None, depth_range=None, n=None) -> pd.DataFrame:
    vars = _as_list(vars)
    if product == "events":
        return _read_event_subset_df(file, vars, lon_range, lat_range, time_range,
                                     depth_range, n)
    t_jd = None
    if time_range is not None:
        if product == "climatology":
            warnings.warn("time_range is ignored for a climatology (no time dimension).",
                          stacklevel=3)
        else:
            try:
                t_jd = to_jd(time_range)
            except ValueError:
                raise ValueError('time_range must be dates, e.g. ("2010-01-01", '
                                 '"2010-12-31").') from None
    sub = _core.hw3_read_subset(file, lon_range=pair(lon_range, "lon_range"),
                                lat_range=pair(lat_range, "lat_range"), t_jd_range=t_jd,
                                depth_range=pair(depth_range, "depth_range"), vars=vars,
                                max_rows=-1 if n is None else int(n))
    df = subset_df(sub)
    if n is not None:
        df = df.head(int(n))
    return df


def _read_event_subset_df(file, vars=None, lon_range=None, lat_range=None,
                          time_range=None, depth_range=None, n=None) -> pd.DataFrame:
    df = event_df(_core.hw3_read_event_nc(file))
    keep = np.ones(len(df), dtype=bool)
    if lon_range is not None:
        lo, hi = min(lon_range), max(lon_range)
        keep &= (df["lon"] >= lo).to_numpy() & (df["lon"] <= hi).to_numpy()
    if lat_range is not None:
        lo, hi = min(lat_range), max(lat_range)
        keep &= (df["lat"] >= lo).to_numpy() & (df["lat"] <= hi).to_numpy()
    if time_range is not None:
        try:
            d = to_datetime64(time_range).astype("datetime64[ns]")
        except ValueError:
            raise ValueError("time_range must be dates.") from None
        keep &= ((df["date_end"] >= d.min()).to_numpy()
                 & (df["date_start"] <= d.max()).to_numpy())
    if depth_range is not None:
        if "depth" not in df:
            warnings.warn("depth_range is ignored for events with no depth column "
                          "(not from a depth-resolved run).", stacklevel=4)
        else:
            lo, hi = min(depth_range), max(depth_range)
            keep &= (df["depth"] >= lo).to_numpy() & (df["depth"] <= hi).to_numpy()
    df = df[keep]
    if vars is not None:
        miss = [v for v in vars if v not in df.columns]
        if miss:
            raise ValueError(f"Variable(s) not in events: {', '.join(miss)}. "
                             f"Available: {', '.join(df.columns)}")
        cols = ["lon", "lat", "event_no"] + [v for v in vars
                                             if v not in ("lon", "lat", "event_no")]
        df = df[list(dict.fromkeys(cols))]
    if n is not None:
        df = df.head(int(n))
    return df.reset_index(drop=True)


def _iterate_chunks(file: str, ctype: str, chunk_size: int):
    if ctype == "clim":
        cd = _core.hw3_read_clim_nc(file)
        npix = cd["nlon"] * cd["nlat"] * max(1, cd["ndepth"])
        step = max(1, chunk_size // cd["ndoy"])
        for start in range(0, npix, step):
            yield clim_df(cd, start, min(npix, start + step))
    elif ctype == "daily":
        dd = _core.hw3_read_daily_nc(file)
        npix = dd["nlon"] * dd["nlat"] * max(1, dd["ndepth"])
        step = max(1, chunk_size // max(1, dd["ntime"]))
        for start in range(0, npix, step):
            yield daily_df(dd, start, min(npix, start + step))
    else:
        ev = _core.hw3_read_event_nc(file)
        if ev["nevents"] == 0:
            yield event_df(ev, [])
            return
        for start in range(0, ev["nevents"], chunk_size):
            yield event_df(ev, np.arange(start, min(ev["nevents"], start + chunk_size)))


# ---- Writers ----

def _validate_out_dir(file_out: str) -> None:
    out_dir = os.path.dirname(os.path.abspath(file_out))
    if not os.path.isdir(out_dir):
        raise FileNotFoundError(f"Output directory does not exist: {os.path.dirname(file_out)}")


def _validate_save_file(file_out: str) -> str:
    _validate_out_dir(file_out)
    ext = os.path.splitext(file_out)[1].lower().lstrip(".")
    if ext == "rds":
        raise ValueError(".rds is an R format; use .csv or .parquet from Python.")
    if ext not in ("csv", "parquet"):
        raise ValueError("file_out must end in .csv, .parquet or .zarr.")
    return ext


def _write_df(df: pd.DataFrame, file_out: str, fmt: str) -> None:
    if fmt == "csv":
        df.to_csv(file_out, index=False)
    else:
        _require_pyarrow()
        df.to_parquet(file_out, index=False)


def _require_pyarrow():
    try:
        import pyarrow.parquet as pq
    except ImportError:
        raise ImportError("The 'pyarrow' package is required for Parquet export. "
                          "Install with: pip install pyarrow") from None
    return pq


def _as_list(x):
    if x is None:
        return None
    if isinstance(x, str):
        return [x]
    return list(x)


# ---- Post-compute console summary ----

def console_summary(file: str, quiet: bool = False, n: int = 5) -> None:
    """Head, tail and light statistics of a written product (R's .hw3_console_summary)."""
    if quiet or not os.path.exists(file):
        return
    meta = _core.hw3_file_meta(file)
    product = meta["product"]
    label = {"climatology": "Climatology", "events": "Events",
             "daily": "Daily series", "protoevents": "Proto-events"}.get(product, product)
    bar = "-" * 66
    print(f"\n{bar}\n{label} written to: {file}")
    line = f"Rows (long format): {meta['nrows']:,}"
    if product != "events":
        line += f"   grid: {meta['nlon']} lon x {meta['nlat']} lat"
    print(line)
    if meta["nrows"] > PREVIEW_LIMIT:
        print(f"(too large to preview inline; use hw3_export() to read or export)\n{bar}")
        return
    df = read_product_df(file, product)
    print("\nHead:")
    print(df.head(n).to_string())
    print("\nTail:")
    print(df.tail(n).to_string())
    print("\nSummary:")
    _print_stats(df, product)
    print(f'\nExamine with  hw3_export("{file}", n=20)\n'
          f'or export with hw3_export("{file}", file_out="out.csv")  (.csv/.parquet)\n{bar}')


def _print_stats(df: pd.DataFrame, product: str) -> None:
    def rng(x):
        x = np.asarray(x, dtype=float)
        x = x[np.isfinite(x)]
        return (np.nan, np.nan) if x.size == 0 else (x.min(), x.max())

    def fmt(v):
        return " to ".join(f"{x:.4g}" for x in v)

    def cats(x):
        x = pd.Series(x).dropna()
        if len(x):
            tb = x.value_counts().sort_index()
            print("  category: " + "  ".join(f"{k}={v}" for k, v in tb.items()))

    if product == "climatology":
        vp = df.loc[np.isfinite(df["seas"]), ["lon", "lat"]].drop_duplicates()
        print(f"  ocean pixels (valid climatology): {len(vp)}")
        print(f"  seas:   {fmt(rng(df['seas']))}")
        print(f"  thresh: {fmt(rng(df['thresh']))}")
    elif product == "events":
        print(f"  events: {len(df)}   pixels with events: "
              f"{len(df[['lon', 'lat']].drop_duplicates())}")
        if len(df):
            print(f"  dates:  {df['date_start'].min():%Y-%m-%d} to "
                  f"{df['date_end'].max():%Y-%m-%d}")
            print(f"  duration (days): {fmt(rng(df['duration']))}")
            print(f"  intensity_max:   {fmt(rng(df['intensity_max']))}")
            if "category" in df:
                cats(df["category"])
    else:
        if "t" in df:
            print(f"  dates: {df['t'].min():%Y-%m-%d} to {df['t'].max():%Y-%m-%d}")
        if "event" in df:
            print(f"  event-days: {int(df['event'].sum())}")
        if "threshCriterion" in df:
            print(f"  threshCriterion days: {int(df['threshCriterion'].sum())}")
        if "category" in df:
            cats(df["category"])
