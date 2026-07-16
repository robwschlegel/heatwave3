# Tests for exceedance3(): static-threshold event detection that bypasses
# ts2clm3 by writing a constant climatology (hw3_write_const_clim) and running
# the standard C++ detector on it.

test_that("exceedance3 detects periods above a fixed threshold", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  event_file <- tempfile(fileext = ".nc")
  out <- exceedance3(sst_file, event_file, threshold = 295, var_name = "analysed_sst")

  expect_equal(out, event_file)
  expect_true(file.exists(event_file))
  expect_equal(hw3_file_meta(event_file)$product, "events")

  ev <- hw3_export(event_file)
  expect_true(nrow(ev) > 0)
  expect_true(all(ev$intensity_max > 0))
})

test_that("exceedance3 below = TRUE detects periods under the threshold", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  event_file <- tempfile(fileext = ".nc")
  exceedance3(sst_file, event_file, threshold = 290, below = TRUE)

  ev <- hw3_export(event_file)
  expect_true(nrow(ev) > 0)
  # heatwaveR convention: cold-spell intensities are reported positive
  expect_true(all(ev$intensity_mean > 0, na.rm = TRUE))
})

test_that("exceedance3 requires a threshold", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  expect_error(exceedance3(sst_file, tempfile(fileext = ".nc")),
               "threshold value must be provided")
})

test_that("exceedance3 respects lon/lat/time subsetting and minDuration", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  event_file <- tempfile(fileext = ".nc")
  exceedance3(sst_file, event_file, threshold = 295,
              lon_range = c(26.5, 26.6), lat_range = c(-34.15, -34.05),
              minDuration = 3L, maxGap = 1L, joinAcrossGaps = FALSE,
              roundRes = 2L)

  ev <- hw3_export(event_file)
  expect_true(all(ev$duration >= 3))
})
