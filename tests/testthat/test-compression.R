# compression: zlib level of the outputs, and the parallel (HDF5 direct chunk)
# writer against netCDF's own serial writer.

cp <- c("1982-01-01", "2011-12-31")
sst <- system.file("extdata/sst_test.nc", package = "heatwave3")

products <- function(stem) {
  lapply(c(clim = "_clim.nc", events = "_events.nc", daily = "_events_daily.nc"),
         function(x) hw3_export(paste0(stem, x)))
}

sorted_events <- function(ev) ev[order(ev$lon, ev$lat, ev$event_no), ]

test_that("compression must be a whole number from 0 to 9", {
  for (bad in list(-1, 10, 1.5, "1", NA, c(1, 2))) {
    expect_error(ts2clm3(sst, tempfile(), cp, compression = bad, quiet = TRUE),
                 "compression must be a whole number from 0 to 9")
  }
  expect_error(exceedance3(sst, tempfile(fileext = ".nc"), threshold = 296,
                           compression = 10, quiet = TRUE),
               "compression must be a whole number from 0 to 9")
})

test_that("every compression level gives the same values", {
  stems <- c(l0 = tempfile(), l1 = tempfile(), l9 = tempfile())
  for (lv in c(0L, 1L, 9L)) {
    detect3(sst, stems[[paste0("l", lv)]], cp, category = TRUE, daily = "also",
            n_threads = 1L, compression = lv, quiet = TRUE)
  }
  ref <- products(stems[["l1"]])
  for (s in stems[c("l0", "l9")]) expect_identical(products(s), ref)
  daily <- paste0(stems, "_events_daily.nc")
  expect_gt(file.size(daily[1]), 2 * file.size(daily[2]))
  expect_lte(file.size(daily[3]), file.size(daily[2]))
})

test_that("the parallel writer gives the same files as netCDF's serial writer", {
  run <- function(serial, ...) {
    stem <- tempfile()
    withr_env <- Sys.getenv("HEATWAVE3_SERIAL_WRITE", unset = NA)
    Sys.setenv(HEATWAVE3_SERIAL_WRITE = if (serial) "1" else "0")
    on.exit(if (is.na(withr_env)) Sys.unsetenv("HEATWAVE3_SERIAL_WRITE")
            else Sys.setenv(HEATWAVE3_SERIAL_WRITE = withr_env))
    detect3(sst, stem, cp, category = TRUE, daily = "also", n_threads = 2L,
            quiet = TRUE, ...)
    detect_event3(sst, paste0(stem, "p"), clim_file = paste0(stem, "_clim.nc"),
                  protoEvent = TRUE, n_threads = 2L, quiet = TRUE, ...)
    p <- products(stem)
    p$events <- sorted_events(p$events)
    p$proto <- hw3_export(paste0(stem, "p_protoevents.nc"))
    p
  }
  for (lv in c(1L, 6L)) {
    par <- run(FALSE, compression = lv)
    ser <- run(TRUE, compression = lv)
    rownames(par$events) <- rownames(ser$events) <- NULL
    expect_identical(par, ser)
  }
})
