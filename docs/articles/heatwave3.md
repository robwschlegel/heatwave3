<div id="main" class="col-md-9" role="main">

# Getting started with heatwave3

<div class="section level2">

## Introduction

**heatwave3** detects marine heatwaves (MHWs) and cold-spells directly
on gridded NetCDF data using the Hobday et al. (2016, 2018) definition.
All core algorithms are implemented in C++ with OpenMP parallelism,
making it orders of magnitude faster than per-pixel processing with
heatwaveR.

The typical workflow has three steps:

1.  `ts2clm3()`, to compute seasonal and threshold climatologies
2.  `detect_event3()`, to detect events from the climatologies
3.  Analyse/visualise with `category3()`, `block_average3()`,
    `detect_blob3()`, or the plotting functions

</div>

<div class="section level2">

## Data input

heatwave3 accepts three kinds of SST input:

-   **A single multi-timestep NetCDF file** (such as a merged time
    series)
-   **A vector of daily NetCDF file paths**
-   **A directory** containing daily `.nc` or `.nc4` files

NetCDF variable and dimension names are auto-detected using CF
conventions (`axis`, `standard_name`, `units` attributes) and common
naming patterns. This means the package works with GHRSST, OISST, OSTIA,
ERA5, CMIP6, and other datasets without needing to specify dimension
names.

<div id="cb1" class="sourceCode">

``` r
library(heatwave3)

# Single multi-timestep file
sst_file <- "/path/to/sst_merged.nc"
clim_file <- tempfile(fileext = ".nc")
ts2clm3(file_in = sst_file,
        file_out = clim_file,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        lon_range = c(15, 35), lat_range = c(-38, -28),
        n_threads = 4)

# Directory of daily files
ts2clm3(file_in = "/path/to/daily_sst/",
        file_out = clim_file,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        n_threads = 4)

# Explicit vector of file paths
daily_files <- list.files("/path/to/daily/", pattern = "[.]nc4$",
                          full.names = TRUE)
ts2clm3(file_in = daily_files,
        file_out = clim_file,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        n_threads = 4)
```

</div>

</div>

<div class="section level2">

## Basic pipeline

Here we show the full pipeline using the bundled OSTIA test dataset (3×2
pixels off the Agulhas coast, 1982–2021, 40 years of daily SST).

<div id="cb2" class="sourceCode">

``` r
library(heatwave3)

sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
clim_file <- tempfile(fileext = ".nc")
event_file <- tempfile(fileext = ".nc")

# All-in-one: climatology + event detection
detect3(file_in = sst_file,
        file_out_clim = clim_file,
        file_out_event = event_file,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        n_threads = 2)
#> Reading SST data from /private/var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T/RtmpEGZ5mo/temp_libpath719a3f2010e4/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Computing climatology with 2 thread(s)...
#>   1/6 pixels (16%)  2/6 pixels (33%)  3/6 pixels (50%)  4/6 pixels (66%)  5/6 pixels (83%)  6/6 pixels (100%)
#> Writing climatology to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//Rtmp5p68X5/file758b29b8123e.nc...
#> Done.
#> Reading climatology from /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//Rtmp5p68X5/file758b29b8123e.nc...
#> Reading SST data from /private/var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T/RtmpEGZ5mo/temp_libpath719a3f2010e4/heatwave3/extdata/sst_test.nc...
#> Grid: 2 lon x 3 lat x 14276 time = 6 pixels
#> Detecting events with 2 thread(s)...
#>   1/6 pixels (16%)  2/6 pixels (33%)  3/6 pixels (50%)  4/6 pixels (66%)  5/6 pixels (83%)  6/6 pixels (100%)
#> Found 610 events across 6 pixels
#> Writing events to /var/folders/3w/nmplbnm109b9903rx8z9q0kc0000gn/T//Rtmp5p68X5/file758b73b5b9d7.nc...
#> Done.
```

</div>

</div>

<div class="section level2">

## Per-pixel time series plot

Extract a single pixel and produce a heatwaveR-style event line plot
with flame polygons showing events above the threshold.

<div id="cb3" class="sourceCode">

``` r
event_line3(sst_file = sst_file,
            clim_file = clim_file,
            lon = 26.525, lat = -34.125,
            start_date = "2010-01-01",
            end_date = "2012-12-31")
```

</div>

![](heatwave3_files/figure-html/event-line-1.png)

</div>

<div class="section level2">

## Spatial maps

Visualise the spatial distribution of any event metric across the grid.

<div id="cb4" class="sourceCode">

``` r
plot_metric3(event_file, metric = "intensity_max", summary = "mean")
```

</div>

![](heatwave3_files/figure-html/spatial-map-1.png)

</div>

<div class="section level2">

## Event categories and yearly aggregation

<div id="cb5" class="sourceCode">

``` r
cats <- category3(event_file, clim_file)
head(cats)
#>   event_no    lon     lat  peak_date category intensity_max duration p_moderate
#> 1        1 26.525 -34.125 1982-11-14     <NA>        3.1531       19          0
#> 2        2 26.525 -34.125 1983-04-20     <NA>        3.6150        5          0
#> 3        3 26.525 -34.125 1983-05-30     <NA>        2.0416        6          0
#> 4        4 26.525 -34.125 1983-06-25     <NA>        2.4626        7          0
#> 5        5 26.525 -34.125 1983-07-12     <NA>        1.8050        6          0
#> 6        6 26.525 -34.125 1983-07-21     <NA>        1.9993        7          0
#>   p_strong p_severe p_extreme season
#> 1        0        0         0       
#> 2        0        0         0       
#> 3        0        0         0       
#> 4        0        0         0       
#> 5        0        0         0       
#> 6        0        0         0
table(cats$category)
#> < table of extent 0 >
```

</div>

<div id="cb6" class="sourceCode">

``` r
ba <- block_average3(event_file)
head(ba)
#>      lon     lat year count duration_mean duration_max intensity_mean
#> 1 26.525 -34.125 1982     1      19.00000           19       2.209900
#> 2 26.525 -34.125 1983     5       6.20000            7       2.130060
#> 3 26.525 -34.125 1984     2       7.00000            9       2.016500
#> 4 26.525 -34.125 1985     2       5.00000            5       1.632900
#> 5 26.525 -34.125 1986     5       8.40000           14       2.530060
#> 6 26.525 -34.125 1987     3      11.66667           15       2.399233
#>   intensity_max_mean intensity_max_max intensity_cumulative_mean total_days
#> 1           3.153100            3.1531                  41.98740         19
#> 2           2.384700            3.6150                  12.80478         31
#> 3           2.448550            2.7094                  13.67920         14
#> 4           1.723400            1.9355                   8.16455         10
#> 5           3.153800            3.8774                  21.01420         42
#> 6           2.840733            3.1123                  28.25203         35
#>   total_icum
#> 1    41.9874
#> 2    64.0239
#> 3    27.3584
#> 4    16.3291
#> 5   105.0710
#> 6    84.7561
```

</div>

</div>

<div class="section level2">

## Detrended climatology

By default, `ts2clm3()` follows Hobday et al. (2016) with a fixed
baseline. To remove the long-term warming trend before computing the
climatology (following Jacox et al. 2020), set `detrend = TRUE`:

<div id="cb7" class="sourceCode">

``` r
ts2clm3(file_in = sst_file,
        file_out = tempfile(fileext = ".nc"),
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        detrend = TRUE)
```

</div>

</div>

<div class="section level2">

## Output formats

The primary output is always NetCDF. Additional companion files can be
written alongside in CSV, RDA, or Parquet format:

<div id="cb8" class="sourceCode">

``` r
# Write events as Parquet alongside the NetCDF
detect_event3(file_in = sst_file,
              clim_file = clim_file,
              file_out = event_file,
              save_format = "parquet")

# Or export after the fact
hw3_export(event_file, format = "csv", type = "event")
hw3_export(clim_file, format = "rda", type = "clim")
```

</div>

</div>

<div class="section level2">

## Cold spells

Detect marine cold-spells by setting `pctile = 10` for the climatology
and `coldSpells = TRUE` for detection.

<div id="cb9" class="sourceCode">

``` r
clim_cold <- tempfile(fileext = ".nc")
event_cold <- tempfile(fileext = ".nc")

ts2clm3(file_in = sst_file,
        file_out = clim_cold,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        pctile = 10, n_threads = 2)

detect_event3(file_in = sst_file,
              clim_file = clim_cold,
              file_out = event_cold,
              coldSpells = TRUE, n_threads = 2)
```

</div>

</div>

<div class="section level2">

## Spatial blob detection

The examples below reproduce the spatial blob analysis from Schlegel et
al., using the OSTIA dataset covering the South African coast (15–35°E,
28–38°S). They show heatwave blob footprints, daily area evolution,
centroid trajectories, persistence maps, and the cold-spell equivalents.

All six figure types correspond to the example outputs in
`heatwaveR/dev_doc/figs/`.

<div class="section level3">

### Setup: climatology and blob detection

<div id="cb10" class="sourceCode">

``` r
library(heatwave3)
library(ggplot2)

sst_file <- "/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc"
clim_file <- tempfile(fileext = ".nc")
xlim <- c(15, 35); ylim <- c(-38, -28)

ts2clm3(file_in = sst_file,
        file_out = clim_file,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        lon_range = xlim, lat_range = ylim,
        var_name = "analysed_sst",
        n_threads = 4)

# Detect spatial blobs, including voxel data for footprint maps
blobs <- detect_blob3(file_in = sst_file,
                      clim_file = clim_file,
                      var_name = "analysed_sst",
                      minVoxels = 200L,
                      topN = 10L,
                      rankBy = "cumI_km2_day",
                      return = c("event", "daily", "voxel"))

# Coastline for map overlays
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
land_clip <- sf::st_crop(land, xmin = xlim[1], xmax = xlim[2],
                         ymin = ylim[1], ymax = ylim[2])
```

</div>

</div>

<div class="section level3">

### Figure 1: Top 6 blob peak-day footprints

<div id="cb11" class="sourceCode">

``` r
top6 <- blobs$event[1:6, ]
top6$blob_lbl <- paste0("#", top6$rank, " | ", top6$date_peak)
vox6 <- blobs$voxel[blobs$voxel$blob %in% top6$blob, ]
vox6$blob_lbl <- top6$blob_lbl[match(vox6$blob, top6$blob)]

# Filter to peak date for each blob
vox6_peak <- do.call(rbind, lapply(split(vox6, vox6$blob), function(sub) {
  peak <- top6$date_peak[top6$blob == sub$blob[1]]
  sub[sub$date == peak, ]
}))

ggplot() +
  geom_sf(data = land_clip, fill = "grey92", colour = "grey60", linewidth = 0.2) +
  geom_raster(data = vox6_peak, aes(lon, lat, fill = delta)) +
  scale_fill_viridis_c(name = expression(Delta ~ "(K)"),
                       option = "inferno", direction = 1) +
  facet_wrap(~ blob_lbl, ncol = 3) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "OSTIA spatial heatwave blobs: top 6 by |cumI|, peak-day footprints",
       subtitle = "Region 15-35 E, 38-28 S",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))
```

</div>

</div>

<div class="section level3">

### Figure 2: Daily area of top 6 blobs

<div id="cb12" class="sourceCode">

``` r
daily6 <- blobs$daily[blobs$daily$blob %in% top6$blob, ]
daily6$blob_lbl <- factor(top6$blob_lbl[match(daily6$blob, top6$blob)],
                          levels = top6$blob_lbl)

ggplot(daily6, aes(date, area_km2 / 1000, colour = blob_lbl)) +
  geom_line(linewidth = 0.8) +
  scale_colour_viridis_d(name = "rank", option = "turbo") +
  labs(title = "Daily area of top 6 spatial heatwave blobs",
       x = NULL, y = expression(Area ~ (10^3 ~ km^2))) +
  theme_minimal(base_size = 11)
```

</div>

</div>

<div class="section level3">

### Figure 3: Centroid trajectories of top 3 blobs

<div id="cb13" class="sourceCode">

``` r
top3 <- blobs$event[1:3, ]
top3$blob_lbl <- paste0("#", top3$rank, " (", top3$date_peak, ")")
daily3 <- blobs$daily[blobs$daily$blob %in% top3$blob, ]
daily3$blob_lbl <- factor(top3$blob_lbl[match(daily3$blob, top3$blob)],
                          levels = top3$blob_lbl)

for (b in unique(daily3$blob)) {
  idx <- daily3$blob == b
  daily3$day_from_start[idx] <- as.numeric(
    daily3$date[idx] - min(daily3$date[idx])
  )
}

ggplot() +
  geom_sf(data = land_clip, fill = "grey92", colour = "grey60", linewidth = 0.2) +
  geom_path(data = daily3, aes(centroid_lon, centroid_lat,
                                colour = day_from_start), linewidth = 0.6) +
  geom_point(data = daily3, aes(centroid_lon, centroid_lat,
                                 colour = day_from_start,
                                 size = area_km2 / 1000), alpha = 0.7) +
  scale_colour_viridis_c(name = "day from start") +
  scale_size_continuous(name = expression(Area ~ (10^3 ~ km^2)),
                        range = c(0.5, 5)) +
  facet_wrap(~ blob_lbl, ncol = 3) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "Centroid trajectories of top 3 blobs",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 10)
```

</div>

</div>

<div class="section level3">

### Figure 4: Persistence, days under heatwave per pixel

<div id="cb14" class="sourceCode">

``` r
vox3 <- blobs$voxel[blobs$voxel$blob %in% top3$blob, ]
vox3$blob_lbl <- factor(top3$blob_lbl[match(vox3$blob, top3$blob)],
                        levels = top3$blob_lbl)

# Count days per (blob, lon, lat)
vox3_persist <- aggregate(date ~ blob_lbl + lon + lat, data = vox3,
                          FUN = length)
names(vox3_persist)[4] <- "n_days"

ggplot() +
  geom_sf(data = land_clip, fill = "grey92", colour = "grey60", linewidth = 0.2) +
  geom_raster(data = vox3_persist, aes(lon, lat, fill = n_days)) +
  scale_fill_viridis_c(name = "days under event", option = "plasma") +
  facet_wrap(~ blob_lbl, ncol = 3) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "Footprint persistence: days under heatwave per pixel, top 3 blobs",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 10)
```

</div>

</div>

<div class="section level3">

### Figure 5: Cold-spell blob footprints

<div id="cb15" class="sourceCode">

``` r
clim_cold <- tempfile(fileext = ".nc")
ts2clm3(file_in = sst_file,
        file_out = clim_cold,
        climatologyPeriod = c("1982-01-01", "2011-12-31"),
        lon_range = xlim, lat_range = ylim,
        var_name = "analysed_sst",
        pctile = 10, n_threads = 4)

mcs_blobs <- detect_blob3(file_in = sst_file,
                          clim_file = clim_cold,
                          var_name = "analysed_sst",
                          coldSpells = TRUE,
                          minVoxels = 200L,
                          topN = 10L,
                          rankBy = "cumI_km2_day",
                          return = c("event", "daily", "voxel"))
```

</div>

<div id="cb16" class="sourceCode">

``` r
mcs6 <- mcs_blobs$event[1:6, ]
mcs6$blob_lbl <- paste0("#", mcs6$rank, " | ", mcs6$date_peak)
mvox6 <- mcs_blobs$voxel[mcs_blobs$voxel$blob %in% mcs6$blob, ]
mvox6$blob_lbl <- mcs6$blob_lbl[match(mvox6$blob, mcs6$blob)]

mvox6_peak <- do.call(rbind, lapply(split(mvox6, mvox6$blob), function(sub) {
  peak <- mcs6$date_peak[mcs6$blob == sub$blob[1]]
  sub[sub$date == peak, ]
}))
mvox6_peak$delta <- -mvox6_peak$delta  # show negative anomalies

ggplot() +
  geom_sf(data = land_clip, fill = "grey92", colour = "grey60", linewidth = 0.2) +
  geom_raster(data = mvox6_peak, aes(lon, lat, fill = delta)) +
  scale_fill_viridis_c(name = expression(Delta ~ "(K)"),
                       option = "mako", direction = 1) +
  facet_wrap(~ blob_lbl, ncol = 3) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "OSTIA cold-spell blobs: top 6 by |cumI|, peak-day footprints",
       subtitle = "Region 15-35 E, 38-28 S",
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold"))
```

</div>

</div>

<div class="section level3">

### Figure 6: Cold-spell daily area

<div id="cb17" class="sourceCode">

``` r
mdaily6 <- mcs_blobs$daily[mcs_blobs$daily$blob %in% mcs6$blob, ]
mdaily6$blob_lbl <- factor(mcs6$blob_lbl[match(mdaily6$blob, mcs6$blob)],
                           levels = mcs6$blob_lbl)

ggplot(mdaily6, aes(date, area_km2 / 1000, colour = blob_lbl)) +
  geom_line(linewidth = 0.8) +
  scale_colour_viridis_d(name = "rank", option = "mako") +
  labs(title = "Daily area of top 6 spatial cold-spell blobs",
       x = NULL, y = expression(Area ~ (10^3 ~ km^2))) +
  theme_minimal(base_size = 11)
```

</div>

</div>

</div>

<div class="section level2">

## Performance

heatwave3 runs much faster than per-pixel processing with heatwaveR. On
a 20×20 grid (400 pixels, 14,276 daily time steps each):

| Method                         | Time       | Speedup     |
|--------------------------------|------------|-------------|
| heatwaveR (serial, 400 pixels) | ca. 38 sec | 1×          |
| heatwave3 (4 threads)          | 0.53 sec   | **ca. 71×** |

The speedup comes from:

1.  **C++ implementation** of all core algorithms (climatology + event
    detection)
2.  **OpenMP parallelism** across pixels
3.  **Direct NetCDF I/O** via libnetcdf (no R overhead)

</div>

</div>
