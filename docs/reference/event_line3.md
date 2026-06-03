<div id="main" class="col-md-9" role="main">

# Plot a per-pixel MHW time series from NetCDF output

<div class="ref-description section level2">

Extracts a single pixel's time series from the SST and climatology
NetCDF files and produces an event\_line-style plot with flame polygons.
The SST and climatology data are read via the C++ backend; no external
NetCDF packages are required.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
event_line3(
  sst_file,
  clim_file,
  lon,
  lat,
  var_name = NULL,
  start_date = NULL,
  end_date = NULL,
  spread = 150,
  metric = "intensity_cumulative",
  event_file = NULL,
  coldSpells = FALSE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   sst\_file:

    Path to the SST NetCDF file (or directory of daily files).

-   clim\_file:

    Path to the climatology NetCDF from `ts2clm3`.

-   lon:

    Longitude of the pixel to plot.

-   lat:

    Latitude of the pixel to plot.

-   var\_name:

    SST variable name. If `NULL`, auto-detected.

-   start\_date, end\_date:

    Optional date range for the plot window (character, for example
    `"2010-01-01"`). If both are `NULL`, the plot is centred on the most
    intense event (see `spread`).

-   spread:

    Number of days before and after the peak event to show. Default
    `150`. Only used when `start_date`/`end_date` are not set.

-   metric:

    Event metric to use for selecting the peak event when `event_file`
    is supplied. Default `"intensity_cumulative"`.

-   event\_file:

    Optional path to event NetCDF for centring the window on the most
    intense event.

-   coldSpells:

    Logical. Render cold-spell (blue) or heatwave (red) flames? Default
    `FALSE`.

</div>

<div class="section level2">

## Value

A ggplot object that can be further customised with the standard `+`
operator.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
# \donttest{
sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
clim_file <- tempfile(fileext = ".nc")

ts2clm3(sst_file, clim_file,
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
#> Writing climatology to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//RtmpEGZ5mo/file719a6f51e8ac.nc...
#> Done.

event_line3(sst_file, clim_file,
            lon = 26.525, lat = -34.125,
            start_date = "2010-01-01", end_date = "2012-12-31")

# }
```

</div>

</div>

</div>
