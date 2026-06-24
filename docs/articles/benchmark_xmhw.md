# Performance benchmark: heatwave3 vs xmhw

## Overview

This vignette benchmarks `heatwave3` against the python module `xmhw` on
the OSTIA gridded SST dataset for the French coastline of the
Mediterranean Sea (104 lon × 60 lat, 6,240 ocean pixels, 15,706 daily
time steps, 1982–2024).

We compare three configurations:

1.  **heatwaveR**, serial per-pixel R + C++ (the standard approach)
2.  **heatwave3 (1 thread)**, pure C++ with single-threaded execution
3.  **xmhw (1 thread)**, pure python without dask
4.  **heatwave3 (12 threads)**, pure C++ with multi-threading
5.  **xmhw (12 threads)**, pure python with dask enabled

``` r

# Necessary R packages
library(reticulate)
library(ncdf4)
library(heatwaveR)
library(heatwave3)
```

## Test data

``` r

# Choose bounding box
bbox <- c(
  2.5,   # °E  (just west of Cerbère)
  41.0,  # °N  (south of Corsica, open sea)
  7.7,   # °E  (Menton / Ligurian coast)
  44.0   # °N  (northern fringe, pre-Alps shelf)
)

# Select daterange
time_range <- c(
  "1982-01-01",
  "2024-12-31"
)

# Simply download the data proscribed here via the UI
# https://data.marine.copernicus.eu/product/SST_GLO_SST_L4_REP_OBSERVATIONS_010_011/download?dataset=METOFFICE-GLO-SST-L4-REP-OBS-SST_202003

# The NetCDF file for further use
sst_file <- "/home/calanus/data/OSTIA/METOFFICE-GLO-SST-L4-REP-OBS-SST_French_Med.nc"
# 104 lon × 60 lat at 0.05° resolution
# 2.5–7.7°E, 41–44°N (French Mediterranean coastline)
# 15,706 daily time steps (1982-01-01 to 2024-12-31)
# 6,240 ocean pixels (remainder is land)
```

## Python environment

We will control the comparison across R and python by using an existing
python v3.12 environemnt created for this tutorial with the following
modules installed: `netCDF4`, `xarray`, `dask`, and `xmhw`.

``` r

# Activate the conda environment and load functionality into R environment
use_condaenv("RiOMar", required = TRUE)
nc4 <- import("netCDF4")
xmhw <- import("xmhw.xmhw")
xr <- import("xarray")
dask <- import("dask")
builtins <- import_builtins()

# Pull out the two key functions
x_threshold <- xmhw$threshold
x_detect <- xmhw$detect
py_slice <- builtins$slice
```

## heatwaveR: serial baseline

`heatwaveR` processes one pixel at a time. We time a 20-pixel sample
(data pre-loaded to isolate computation from I/O) and extrapolate to the
full grid.

``` r

nc <- nc_open(sst_file)
time_raw <- ncvar_get(nc, "time")
dates <- as.Date("1970-01-01") + time_raw / 86400

# Pre-load one longitude column (20 lat pixels)
sst_col <- ncvar_get(nc, "analysed_sst",
                     start = c(1, 1, 1), 
                     count = c(1, 20, -1))
nc_close(nc)

ocean_idx <- which(!is.na(sst_col[, 1]))
n_test <- 20L

t_hwr <- system.time({
  for (j in ocean_idx[1:n_test]) {
    clm <- ts2clm(data.frame(t = dates, temp = as.numeric(sst_col[j, ])),
                  climatologyPeriod = c("1982-01-01", "2011-12-31"))
    ev <- detect_event(clm)
  }
})

per_pixel <- t_hwr[3] / n_test
n_ocean <- 6240  # pre-counted
estimated <- per_pixel * n_ocean

cat("Per pixel:      ", round(per_pixel, 3), "sec\n")
cat("Estimated total:", round(estimated), "sec (",
    round(estimated / 60, 1), "min)\n")
```

## xmhw: single-thread

``` r

# Open sst file with xarray and extract the same subset of data as the heatwaveR example above
x_ds <- xr$open_dataset(sst_file)
x_sst <- x_ds[["analysed_sst"]]$isel(
  longitude = 0L, # first longitude only (1 pixel)
  latitude  = py_slice(0L, 20L) # first 20 latitudes
)

# Check output
# print(x_sst)
# print(x_sst$shape) 

# Determine climatology data
x_clim_time <- system.time(
  x_clim <- x_threshold(x_sst)
)

# Detect events
x_event_time <- system.time(
  x_event <- x_detect(x_sst, x_clim['thresh'], x_clim['seas'])
)

# Time estimates
x_clim_est <- x_clim_time[3]/20*6240
x_event_est <- x_event_time[3]/20*6240

cat("Climatology (estimated):", round(x_clim_est), "sec (",
    round((x_clim_est) / 60, 1), "min)\n")
cat("Detection (estimated):  ", round(x_event_est), "sec (",
    round((x_event_est) / 60, 1), "min)\n")
cat("Total (estimated):", round(x_clim_est + x_event_est), "sec (",
    round((x_clim_est + x_event_est) / 60, 1), "min)\n")
```

## heatwave3: single-threaded

``` r

stem <- file.path(tempdir(), "bench")

t_clim <- system.time(
  ts2clm3(sst_file, name = stem,
          climatologyPeriod = c("1982-01-01", "2011-12-31"),
          var_name = "analysed_sst", n_threads = 1L, quiet = TRUE)
)

t_event <- system.time(
  detect_event3(sst_file, name = stem,
                var_name = "analysed_sst", n_threads = 1L, quiet = TRUE)
)

cat("Climatology:", round(t_clim[3], 1), "sec\n")
cat("Detection:  ", round(t_event[3], 1), "sec\n")
cat("Total:      ", round(t_clim[3] + t_event[3], 1), "sec\n")
```

## heatwave3: 12 threads

``` r

t_clim_12 <- system.time(
  ts2clm3(sst_file, name = stem,
          climatologyPeriod = c("1982-01-01", "2011-12-31"),
          var_name = "analysed_sst", n_threads = 12L, quiet = TRUE)
)

t_event_12 <- system.time(
  detect_event3(sst_file, name = stem,
                var_name = "analysed_sst", n_threads = 12L, quiet = TRUE)
)

cat("Climatology:", round(t_clim_12[3], 1), "sec\n")
cat("Detection:  ", round(t_event_12[3], 1), "sec\n")
cat("Total:      ", round(t_clim_12[3] + t_event_12[3], 1), "sec\n")
```

## xmhw: 12 threads

Note that the `xmhw` will also work with `dask`, which itself can have
parellelism established. This is however outside of the scope of the
capacity of an Rmarkdown file and `reticulate` to be able to replicate
reasonably. The interested researcher is advised to read the tutorial
available on the `xmhw` website:
<https://xmhw.readthedocs.io/en/latest/dask.html>

## Results

Benchmarked on Ubuntu 24.04, 16 cores, R 4.6.0.

| Method | Climatology | Detection | Total |  |
|----|----|----|----|----|
| heatwaveR (serial) | n/a | n/a | ca. 1,262 sec (21 min) | 1× |
| xmhw (serial) | ca. 2,115 sec | ca. 300 sec | ca. 2,415 sec (40.2 min) | 0.5× |
| heatwave3 (1 thread) | 17.5 sec | 5 sec | **22.5 sec** | **56×** |
| heatwave3 (12 threads) | 3.7 sec | 3.7 sec | **7.4 sec** | **170×** |

Note then that the base `heatwaveR` script is already almost twice as
fast as `xmhw` if this analysis has not been optimised for `dask`. Also
note that as part of the operation of `heatwave3` the output of the
files are already saved to local disk in an orderly fashion. The saving
of the output of `xmhw` would add additional total run-time.
