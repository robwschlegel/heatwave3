# Input checks that stop before any computation: missing files, empty
# directories and malformed arguments.

sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
cp <- c("1982-01-01", "2011-12-31")
empty_dir <- function() { d <- tempfile(); dir.create(d); d }

test_that("ts2clm3 checks its inputs", {
  skip_if(sst_file == "", "sst_test.nc not available")
  expect_error(ts2clm3(sst_file, tempfile(), cp, depth_range = 5, quiet = TRUE),
               "depth_range must be a numeric vector of length 2")
  expect_error(ts2clm3(sst_file, file.path(tempfile(), "x"), cp, quiet = TRUE),
               "Output directory does not exist")
  expect_error(ts2clm3(empty_dir(), tempfile(), cp, quiet = TRUE),
               "No .nc or .nc4 files found")
})

test_that("detect_event3 checks its inputs", {
  skip_if(sst_file == "", "sst_test.nc not available")
  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = cp, quiet = TRUE)
  clim <- paste0(stem, "_clim.nc")
  expect_error(detect_event3(sst_file, file.path(tempfile(), "x"), clim_file = clim,
                             quiet = TRUE),
               "Output directory does not exist")
  expect_error(detect_event3(empty_dir(), tempfile(), clim_file = clim, quiet = TRUE),
               "No .nc or .nc4 files found")
  expect_error(detect_event3(sst_file, tempfile(), clim_file = clim,
                             threshClim2 = c(clim, clim), quiet = TRUE),
               "threshClim2 may only contain multiple files")
  expect_error(detect_event3(tempfile(fileext = ".nc"), tempfile(), clim_file = clim,
                             quiet = TRUE),
               "SST file does not exist")
})

test_that("category_daily3 checks its inputs", {
  skip_if(sst_file == "", "sst_test.nc not available")
  stem <- tempfile()
  detect3(sst_file, name = stem, climatologyPeriod = cp, quiet = TRUE)
  clim <- paste0(stem, "_clim.nc"); ev <- paste0(stem, "_events.nc")
  window <- c("2010-01-01", "2010-01-31")
  expect_error(category_daily3(sst_file, tempfile(fileext = ".nc"), ev, window),
               "Climatology file does not exist")
  expect_error(category_daily3(sst_file, clim, tempfile(fileext = ".nc"), window),
               "Events file does not exist")
  expect_error(category_daily3(sst_file, clim, ev, c("2010-01-01", "not a date")),
               "time_range must be dates")
  expect_error(category_daily3(empty_dir(), clim, ev, window),
               "No .nc or .nc4 files found in the SST directory")
})

test_that("event_line3 checks its inputs", {
  skip_if(sst_file == "", "sst_test.nc not available")
  skip_if_not_installed("ggplot2")
  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = cp, quiet = TRUE)
  clim <- paste0(stem, "_clim.nc")
  expect_error(event_line3(empty_dir(), clim, lon = 26.525, lat = -34.125),
               "No .nc/.nc4 files in")
  expect_error(event_line3(sst_file, clim, lon = 0, lat = 0),
               "No coordinates found within specified range")
})
