# Tests for plot_metric3(): spatial map of aggregated event metrics
# (hw3_read_metric_summary() does the C++ aggregation; this file just
# renders it with ggplot2).

test_that("plot_metric3 returns a ggplot object for the default metric", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  detect3(sst_file, name = stem,
          climatologyPeriod = c("1982-01-01", "2011-12-31"), quiet = TRUE)

  p <- plot_metric3(paste0(stem, "_events.nc"), coastline = FALSE)
  expect_s3_class(p, "ggplot")
  expect_equal(p$scales$scales[[1]]$name, "intensity_max (mean)")
})

test_that("plot_metric3 accepts alternate metric/summary combinations and extra scale args", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  detect3(sst_file, name = stem,
          climatologyPeriod = c("1982-01-01", "2011-12-31"), quiet = TRUE)
  events_file <- paste0(stem, "_events.nc")

  p_dur <- plot_metric3(events_file, metric = "duration", summary = "max",
                        coastline = FALSE, option = "magma")
  expect_s3_class(p_dur, "ggplot")
  expect_equal(p_dur$scales$scales[[1]]$name, "duration (max)")

  p_count <- plot_metric3(events_file, metric = "intensity_mean", summary = "count",
                          coastline = FALSE)
  expect_s3_class(p_count, "ggplot")
})

test_that("plot_metric3 adds a coastline layer when coastline = TRUE and dependencies are present", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  skip_if_not_installed("rnaturalearth")
  skip_if_not_installed("sf")

  stem <- tempfile()
  detect3(sst_file, name = stem,
          climatologyPeriod = c("1982-01-01", "2011-12-31"), quiet = TRUE)

  p <- plot_metric3(paste0(stem, "_events.nc"), coastline = TRUE)
  expect_s3_class(p, "ggplot")
  is_sf_layer <- vapply(p$layers, function(l) inherits(l$geom, "GeomSf"), logical(1))
  expect_true(any(is_sf_layer))
})
