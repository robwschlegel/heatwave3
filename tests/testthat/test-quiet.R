# quiet = TRUE silences all console output from ts2clm3()/detect_event3()/
# detect3(), including the C++ progress messages, but not warnings.

test_that("quiet = TRUE silences the C++ messages and the summary", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  stem <- tempfile()
  cp <- c("1982-01-01", "2011-12-31")

  expect_identical(capture.output(ts2clm3(sst_file, name = stem, climatologyPeriod = cp,
                                          quiet = TRUE)), character(0))
  expect_identical(capture.output(detect_event3(sst_file, name = stem, daily = "also",
                                                quiet = TRUE)), character(0))
  expect_identical(capture.output(detect3(sst_file, name = tempfile(), climatologyPeriod = cp,
                                          quiet = TRUE)), character(0))

  # The setting is restored afterwards, so the next call is verbose again
  expect_false(hw3_set_quiet(FALSE))
  out <- capture.output(detect_event3(sst_file, name = stem, quiet = FALSE))
  expect_true(any(grepl("Detecting events", out)))
  expect_true(any(grepl("Events written to", out)))
})

test_that("quiet = TRUE still reports warnings about skipped files", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  files <- c(sst_file, file.path(tempdir(), "missing_hw3_quiet.nc"), sst_file)
  msgs <- capture.output(
    out <- capture.output(ts2clm3(files, name = tempfile(),
                                  climatologyPeriod = c("1982-01-01", "2011-12-31"),
                                  skip_bad_files = TRUE, quiet = TRUE)),
    type = "message")
  expect_identical(out, character(0))
  expect_true(any(grepl("failed to read .*missing_hw3_quiet.nc", msgs)))
})

test_that("quiet is restored even when the call fails", {
  expect_error(detect_event3(system.file("extdata/sst_test.nc", package = "heatwave3"),
                             name = tempfile(), clim_file = tempfile(fileext = ".nc"),
                             quiet = TRUE))
  expect_false(hw3_set_quiet(FALSE))
})
