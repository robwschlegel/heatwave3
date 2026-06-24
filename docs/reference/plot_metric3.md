# Spatial map of event metrics

Creates a map showing the spatial distribution of a chosen event metric
from the event NetCDF output. The per-pixel aggregation is performed in
C++ for efficiency; the R layer handles only the ggplot2 rendering.

## Usage

``` r
plot_metric3(
  event_file,
  metric = "intensity_max",
  summary = "mean",
  coastline = TRUE,
  ...
)
```

## Arguments

- event_file:

  Path to the event NetCDF file from
  [`detect_event3`](https://robwschlegel.github.io/heatwave3/index.html/reference/detect_event3.md).

- metric:

  Character. The event metric to map. Options include `"intensity_max"`
  (default), `"intensity_mean"`, `"intensity_cumulative"`, `"duration"`,
  `"rate_onset"`, `"rate_decline"`, and all `relThresh`/`abs` variants.

- summary:

  Character. How to aggregate across events per pixel. One of `"mean"`
  (default), `"max"`, `"min"`, `"sum"`, or `"count"`.

- coastline:

  Logical. Add a coastline layer? Requires the `rnaturalearth` and `sf`
  packages. Default `TRUE`.

- ...:

  Additional arguments passed to
  [`ggplot2::scale_fill_viridis_c`](https://ggplot2.tidyverse.org/reference/scale_viridis.html).

## Value

A ggplot object. The underlying data is accessible via
[`ggplot2::layer_data()`](https://ggplot2.tidyverse.org/reference/ggplot_build.html)
or by calling `hw3_read_metric_summary()` directly.

## Examples

``` r
# \donttest{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
clim_file <- tempfile(fileext = ".nc")
event_file <- tempfile(fileext = ".nc")

detect3(sst_file, clim_file, event_file,
        climatologyPeriod = c("1982-01-01", "2011-12-31"))
#> Reading SST data from /tmp/RtmptnP3kw/temp_libpath10f40785ac7f0/heatwave3/extdata/sst_test.nc...
#> Error: NetCDF error in inq varid /tmp/RtmpbqjzHu/file162601d0df178.nc: NetCDF: Variable not found

plot_metric3(event_file, metric = "intensity_max", summary = "mean")
#> Error: NetCDF error in open /tmp/RtmpbqjzHu/file162601d0df178.nc: No such file or directory
# }
```
