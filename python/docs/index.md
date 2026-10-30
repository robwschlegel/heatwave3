# heatwave3

**heatwave3** detects marine heatwaves (MHWs) and marine cold-spells (MCSs)
directly in gridded NetCDF data, following the definitions of Hobday et al.
(2016, 2018). The climatology, event detection, categorisation and blob
detection run in multi-threaded C++ with direct NetCDF input and output, so
whole grids are processed without converting them into per-pixel time series.

The Python package is the companion of the
[heatwave3 R package](https://robwschlegel.github.io/heatwave3/). Both call
the same C++ core, so for the same input they produce the same results. The
function names, arguments and output columns are also the same, which means
code moves between the two languages almost line for line.

```{note}
The Python and R packages share one version number: heatwave3 1.5.0 has the
same functions and outputs in both languages, and the test suite checks every
output against the R package. The one difference is Zarr input, which only
Python reads (see "Differences" in [R to Python](r-to-python.md)).
```

## Installation

```bash
pip install heatwave3              # computation and NetCDF/CSV export
pip install "heatwave3[plot]"      # + plotting functions (matplotlib, cartopy)
pip install "heatwave3[parquet]"   # + Parquet export (pyarrow)
```

Pre-built wheels are available for Python 3.11 to 3.14 on Linux (x86-64 and
ARM64) and macOS (Apple Silicon and Intel, macOS 11 or later). They include
their own copy of the NetCDF and HDF5 libraries, so nothing else needs to be
installed. On other platforms, pip builds the package from source, which needs
a C++17 compiler, CMake and the NetCDF-C library (for example
`apt install libnetcdf-dev` or `brew install netcdf`).

## A first example

```python
import heatwave3 as hw3

# Climatology, then events with Hobday et al. (2018) categories. The outputs
# are NetCDF files named from the stem: results/run1_clim.nc, ..._events.nc
out = hw3.detect3("sst.nc", "results/run1", ("1982-01-01", "2011-12-31"),
                  category=True, n_threads=4)

events = hw3.hw3_export(out["events"])      # a pandas DataFrame
hw3.hw3_export(out["events"], file_out="events.csv")
```

The [Getting started](getting-started.md) guide walks through the full
workflow. If you know the R package, the [R to Python](r-to-python.md) page
lists every function and the few differences.

```{toctree}
:maxdepth: 2
:hidden:

getting-started
r-to-python
r-and-python-notebook
benchmark
api
```
