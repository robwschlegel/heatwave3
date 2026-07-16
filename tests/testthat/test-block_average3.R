# Tests for block_average3(): yearly per-pixel aggregation of event metrics.

test_that("block_average3 aggregates event metrics by pixel and year", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  detect3(sst_file, name = stem,
          climatologyPeriod = c("1982-01-01", "2011-12-31"), quiet = TRUE)

  ev <- hw3_export(paste0(stem, "_events.nc"))
  ba <- block_average3(paste0(stem, "_events.nc"))

  expect_s3_class(ba, "data.frame")
  expect_named(ba, c("lon", "lat", "year", "count", "duration_mean",
                     "duration_max", "intensity_mean", "intensity_max_mean",
                     "intensity_max_max", "intensity_cumulative_mean",
                     "total_days", "total_icum"))
  expect_true(nrow(ba) > 0)
  # one row per (pixel, year) with at least one event
  expect_equal(sum(ba$count), nrow(ev))
  expect_true(all(ba$duration_max >= ba$duration_mean))
})
