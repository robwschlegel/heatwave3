# Regression test for CF time parsing: a daily product stamped at noon
# (e.g. OSTIA/GHRSST, 12:00:00) must be dated to the calendar day that contains
# it, not rounded forward a day.

test_that("noon-stamped daily files land on the correct calendar day", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  skip_if_not_installed("ncdf4")

  cp <- c("1982-01-01", "2011-12-31")

  # As-shipped fixture is stamped at midnight.
  stem0 <- tempfile()
  detect3(sst_file, name = stem0, climatologyPeriod = cp, quiet = TRUE)
  ev0 <- hw3_export(paste0(stem0, "_events.nc"))

  # Shift every timestamp by +12 h (the fixture's time units are seconds).
  noon <- tempfile(fileext = ".nc")
  file.copy(sst_file, noon)
  nc <- ncdf4::nc_open(noon, write = TRUE)
  ncdf4::ncvar_put(nc, "time", ncdf4::ncvar_get(nc, "time") + 43200)
  ncdf4::nc_close(nc)

  stem1 <- tempfile()
  detect3(noon, name = stem1, climatologyPeriod = cp, quiet = TRUE)
  ev1 <- hw3_export(paste0(stem1, "_events.nc"))

  # Identical detection and, crucially, identical dates: the half-day stamp
  # must not move events to the next day.
  expect_equal(nrow(ev0), nrow(ev1))
  expect_equal(ev0$date_start, ev1$date_start)
  expect_equal(ev0$date_peak, ev1$date_peak)
  expect_equal(ev0$date_end, ev1$date_end)
})
