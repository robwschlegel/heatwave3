# Gaps in event_line3() not covered by test-depth_range_daily.R: the
# directory/multi-file sst_file branch, event_file-based window centring
# (only reached when start_date/end_date are both absent -- the existing
# depth test passes all three, so the `if` takes the start_date/end_date
# branch and never reaches `else if (!is.null(event_file))`), and the two
# "no data" error paths.

glorys_file <- system.file("extdata/glorys_depth_test.nc", package = "heatwave3")
glorys_dir  <- system.file("extdata/glorys_depth_test_annual", package = "heatwave3")

test_that("event_line3 centres the window on the most intense event when only event_file is given", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  detect3(sst_file, name = stem, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)

  p <- event_line3(sst_file, paste0(stem, "_clim.nc"), lon = 26.525, lat = -34.125,
                   event_file = paste0(stem, "_events.nc"))
  expect_s3_class(p, "ggplot")

  # A pixel with no matching events: event_file is ignored, full series plotted.
  p_full <- event_line3(sst_file, paste0(stem, "_clim.nc"), lon = 26.525, lat = -34.125)
  expect_s3_class(p_full, "ggplot")
})

test_that("event_line3 accepts a directory of daily files (hw3_read_sst_multi)", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")
  skip_if(glorys_dir == "", "glorys_depth_test_annual not available")

  stem <- tempfile()
  ts2clm3(glorys_file, name = stem, climatologyPeriod = c("1993-01-01", "2019-12-31"),
          depth = 0, quiet = TRUE)

  p <- event_line3(glorys_dir, paste0(stem, "_clim.nc"),
                   lon = 25.16667, lat = -34.91667,
                   start_date = "2000-01-01", end_date = "2001-12-31")
  expect_s3_class(p, "ggplot")
})

test_that("event_line3 errors clearly when no data is found", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)
  clim_file <- paste0(stem, "_clim.nc")

  # hw3_read_sst() itself throws for an out-of-range lon/lat before
  # event_line3's own "No data found" check is ever reached (that check
  # appears to be unreachable given current hw3_read_sst() behaviour).
  expect_error(
    event_line3(sst_file, clim_file, lon = 999, lat = 999),
    "No coordinates found"
  )
  expect_error(
    event_line3(sst_file, clim_file, lon = 26.525, lat = -34.125,
               start_date = "1900-01-01", end_date = "1900-01-02"),
    "No data in the specified date range"
  )
})
