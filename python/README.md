# heatwave3

Detect marine heatwaves (MHWs) and marine cold-spells (MCSs) directly in
gridded NetCDF data, following the Hobday et al. (2016, 2018) definitions.

`heatwave3` is the Python interface to the same C++ core as the
[heatwave3 R package](https://robwschlegel.github.io/heatwave3/): climatology,
event detection, categorisation and blob detection run in multi-threaded C++
with direct NetCDF I/O. Function names, arguments and outputs match the R
package, so code translates line for line.

```python
import heatwave3 as hw3

# Climatology + events (+ categories, + daily series) in one call
out = hw3.detect3("sst.nc", "results/run1", ("1982-01-01", "2011-12-31"),
                  category=True, daily="also", n_threads=4)

events = hw3.hw3_export(out["events"])                       # pandas DataFrame
daily = hw3.hw3_export(out["daily"], time_range=("2023-01-01", "2023-12-31"))
hw3.hw3_export(out["events"], file_out="events.parquet")     # or .csv
```

Also available: `ts2clm3()`, `detect_event3()`, `category3()`,
`category_daily3()`, `block_average3()`, `exceedance3()`, `getHW3threads()`
and `setHW3threads()`.

This package is in early development (alpha).
