<div id="main" class="col-md-9" role="main">

# Spatial map of event metrics

<div class="ref-description section level2">

Creates a map showing the spatial distribution of a chosen event metric
from the event NetCDF output. The per-pixel aggregation is performed in
C++ for efficiency; the R layer handles only the ggplot2 rendering.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
plot_metric3(
  event_file,
  metric = "intensity_max",
  summary = "mean",
  coastline = TRUE,
  ...
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   event\_file:

    Path to the event NetCDF file from `detect_event3`.

-   metric:

    Character. The event metric to map. Options include
    `"intensity_max"` (default), `"intensity_mean"`,
    `"intensity_cumulative"`, `"duration"`, `"rate_onset"`,
    `"rate_decline"`, and all `relThresh`/`abs` variants.

-   summary:

    Character. How to aggregate across events per pixel. One of `"mean"`
    (default), `"max"`, `"min"`, `"sum"`, or `"count"`.

-   coastline:

    Logical. Add a coastline layer? Requires the `rnaturalearth` and
    `sf` packages. Default `TRUE`.

-   ...:

    Additional arguments passed to `ggplot2::scale_fill_viridis_c`.

</div>

<div class="section level2">

## Value

A ggplot object. The underlying data is accessible via
`ggplot2::layer_data()` or by calling `hw3_read_metric_summary()`
directly.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# \donttest{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
clim_file <- tempfile(fileext = ".nc")
event_file <- tempfile(fileext = ".nc")

detect3(sst_file, clim_file, event_file,
        climatologyPeriod = c("1982-01-01", "2011-12-31"))
#> Reading SST data from /private/var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T/RtmpEGZ5mo/temp_libpath719a3f2010e4/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Computing climatology with 1 thread(s)...
#> 
  1/6 pixels (16%)
  2/6 pixels (33%)
  3/6 pixels (50%)
  4/6 pixels (66%)
  5/6 pixels (83%)
  6/6 pixels (100%)
#> Writing climatology to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a45e8ea76.nc...
#> Done.
#> Reading climatology from /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a45e8ea76.nc...
#> Reading SST data from /private/var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T/RtmpEGZ5mo/temp_libpath719a3f2010e4/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Detecting events with 1 thread(s)...
#> 
  1/6 pixels (16%)
  2/6 pixels (33%)
  3/6 pixels (50%)
  4/6 pixels (66%)
  5/6 pixels (83%)
  6/6 pixels (100%)
#> Found 610 events across 6 pixels
#> Writing events to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a3e790cd5.nc...
#> Done.

plot_metric3(event_file, metric = "intensity_max", summary = "mean")

# }
```

</div>

</div>

</div>
