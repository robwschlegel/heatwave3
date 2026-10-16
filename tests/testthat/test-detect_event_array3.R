# detect_event_array3(): per-day detection on arrays, checked against the
# per-day product of detect_event3() on the same days and climatology.

all_outputs <- c("event", "event_no", "intensity", "category")

# The bundled series as [time, pixel] matrices, with the climatology matched
# to each day, and detect_event3()'s per-day product for it.
array_case <- function(coldSpells) {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  stem <- tempfile()
  clim <- ts2clm3(sst_file, name = stem, pctile = if (coldSpells) 10 else 90,
                  climatologyPeriod = c("1982-01-01", "2011-12-31"), quiet = TRUE)
  detect_event3(sst_file, stem, clim_file = clim, coldSpells = coldSpells,
                category = TRUE, daily = "only", quiet = TRUE)
  g <- hw3_read_sst(sst_file, "analysed_sst")
  cl <- hw3_read_clim_nc(clim)
  doy <- vapply(g$time_days, hw3_jd_to_doy, integer(1))
  list(sst_file = sst_file, clim = clim, time_days = g$time_days, lon = g$lon, lat = g$lat,
       temp = matrix(g$sst, nrow = g$ntime),
       seas = matrix(cl$seas, nrow = 366)[doy, ],
       thresh = matrix(cl$thresh, nrow = 366)[doy, ],
       daily = read_daily(paste0(stem, "_events_daily.nc")))
}

read_daily <- function(daily_file) {
  d <- hw3_read_daily_nc(daily_file)
  lapply(d[all_outputs], function(x) matrix(x, nrow = d$ntime))
}

expect_same_as_daily <- function(res, daily) {
  expect_identical(res$event, daily$event == 1L)
  expect_identical(res$event_no, daily$event_no)
  # the internal reader returns the core's NaN for missing days; the array
  # function returns R's NA, as heatwave3's user-facing functions do
  daily$intensity[is.nan(daily$intensity)] <- NA_real_
  expect_identical(res$intensity, daily$intensity)
  expect_identical(res$category, daily$category)
}

for (cold in c(FALSE, TRUE)) {
  test_that(paste("detect_event_array3 matches detect_event3 on the full series,",
                  if (cold) "cold-spells" else "heatwaves"), {
    case <- array_case(cold)
    res <- detect_event_array3(case$temp, case$seas, case$thresh, coldSpells = cold,
                               outputs = all_outputs)
    expect_true(any(res$event) && any(res$category > 0))
    expect_same_as_daily(res, case$daily)
  })

  test_that(paste("detect_event_array3 matches detect_event3 on a window,",
                  if (cold) "cold-spells" else "heatwaves"), {
    skip_if_not_installed("ncdf4")
    case <- array_case(cold)
    # A 17-day window that starts 2 days into the longest event (cut off at
    # the start), with missing days near its end
    runs <- do.call(rbind, lapply(seq_len(ncol(case$temp)), function(px) {
      r <- rle(case$daily$event[, px] == 1L)
      data.frame(start = cumsum(r$lengths) - r$lengths + 1, length = r$lengths,
                 event = r$values)
    }))
    longest <- runs[runs$event, ][which.max(runs$length[runs$event]), ]
    expect_gt(longest$length, 9)
    days <- seq(longest$start + 2, length.out = 17)
    temp <- case$temp[days, , drop = FALSE]
    temp[c(15, 16), ] <- NA

    window <- tempfile(fileext = ".nc")
    nc_dims <- list(ncdf4::ncdim_def("longitude", "degrees_east", case$lon),
                    ncdf4::ncdim_def("latitude", "degrees_north", case$lat),
                    ncdf4::ncdim_def("time", "days since 1970-01-01",
                                     case$time_days[days] - 2440588L))
    var <- ncdf4::ncvar_def("analysed_sst", "kelvin", nc_dims, missval = NA,
                            prec = "double")
    nc <- ncdf4::nc_create(window, var)
    # [time, pixel] -> [lon, lat, time]; pixel = lon * nlat + lat, as the core
    ncdf4::ncvar_put(nc, var, aperm(array(temp, c(17, 3, 2)), c(3, 2, 1)))
    ncdf4::nc_close(nc)

    stem <- tempfile()
    detect_event3(window, stem, clim_file = case$clim, coldSpells = cold,
                  category = TRUE, daily = "only", quiet = TRUE)
    res <- detect_event_array3(temp, case$seas[days, ], case$thresh[days, ],
                               coldSpells = cold, outputs = all_outputs)
    expect_true(any(res$event))
    expect_same_as_daily(res, read_daily(paste0(stem, "_events_daily.nc")))
  })
}

test_that("detect_event_array3 accepts vectors, matrices and arrays alike", {
  case <- array_case(FALSE)
  flat <- detect_event_array3(case$temp, case$seas, case$thresh, outputs = all_outputs)
  grid <- detect_event_array3(array(case$temp, c(nrow(case$temp), 3, 2)),
                              array(case$seas, c(nrow(case$temp), 3, 2)),
                              array(case$thresh, c(nrow(case$temp), 3, 2)),
                              outputs = all_outputs)
  one <- detect_event_array3(case$temp[, 5], case$seas[, 5], case$thresh[, 5],
                             outputs = all_outputs)
  for (out in all_outputs) {
    expect_identical(dim(grid[[out]]), c(nrow(case$temp), 3L, 2L))
    expect_identical(as.vector(grid[[out]]), as.vector(flat[[out]]))
    expect_null(dim(one[[out]]))
    expect_identical(one[[out]], flat[[out]][, 5])
  }
})

test_that("detect_event_array3 returns the requested outputs, in order, with dimnames", {
  temp <- matrix(20, 30, 2, dimnames = list(NULL, c("a", "b")))
  temp[6:15, 1] <- 25
  seas <- matrix(20, 30, 2)
  thresh <- matrix(21, 30, 2)
  res <- detect_event_array3(temp, seas, thresh)
  expect_named(res, c("event", "intensity", "category"))
  expect_type(res$event, "logical")
  expect_type(res$intensity, "double")
  expect_type(res$category, "integer")
  expect_identical(dimnames(res$event), dimnames(temp))
  expect_identical(which(res$event[, "a"]), 6:15)
  expect_identical(unique(res$category[6:15, 1]), 4L)  # 25 is beyond 4x (21 + 3 x 1)
  single <- detect_event_array3(temp, seas, thresh, outputs = "event_no")
  expect_true(is.matrix(single))
  expect_identical(single[6:15, 1], rep(1L, 10))
  expect_named(detect_event_array3(temp, seas, thresh, outputs = c("category", "event")),
               c("category", "event"))
})

test_that("detect_event_array3 gives empty pixels the defaults", {
  temp <- matrix(20, 30, 3)
  temp[, 1] <- NA
  temp[6:15, 3] <- 25
  seas <- matrix(20, 30, 3)
  seas[, 2] <- NA
  thresh <- matrix(21, 30, 3)
  res <- detect_event_array3(temp, seas, thresh, outputs = all_outputs)
  for (px in 1:2) {
    expect_false(any(res$event[, px]))
    expect_true(all(res$event_no[, px] == 0L) && all(res$category[, px] == 0L))
    expect_true(all(is.na(res$intensity[, px])))
  }
  expect_true(all(res$event[6:15, 3]))
})

test_that("detect_event_array3 passes the detection parameters on", {
  temp <- rep(20, 20)
  temp[3:6] <- 25    # 4 days
  temp[9:14] <- 25   # 6 days, after a 2-day gap
  seas <- rep(20, 20)
  thresh <- rep(21, 20)
  ev <- function(...) detect_event_array3(temp, seas, thresh, outputs = "event", ...)
  expect_identical(which(ev()), 9:14)
  expect_identical(which(ev(minDuration = 4)), 3:14)
  expect_identical(which(ev(minDuration = 4, maxGap = 1)), c(3:6, 9:14))
  expect_identical(
    detect_event_array3(-temp, -seas, -thresh, coldSpells = TRUE, outputs = "event"), ev())
})

test_that("detect_event_array3 rejects invalid input", {
  x <- matrix(0, 10, 2)
  expect_error(detect_event_array3(x, x, x, outputs = character(0)), "'outputs' must name")
  expect_error(detect_event_array3(x, x, x, outputs = c("event", "event")), "'outputs' must name")
  expect_error(detect_event_array3(x, x, x, outputs = "events"), "'outputs' must name")
  expect_error(detect_event_array3(x, x[1:9, ], x), "same dimensions")
  expect_error(detect_event_array3(x, x, as.vector(x)), "same dimensions")
  expect_error(detect_event_array3(x, x, x, minDuration = 0), "minDuration must be at least 1")
  expect_error(detect_event_array3(x, x, x, maxGap = -1), "maxGap must not be negative")
})

test_that("detect_event_array3 gives the same results on any number of threads", {
  case <- array_case(FALSE)
  one <- detect_event_array3(case$temp, case$seas, case$thresh, outputs = all_outputs,
                             n_threads = 1)
  four <- detect_event_array3(case$temp, case$seas, case$thresh, outputs = all_outputs,
                              n_threads = 4)
  expect_identical(four, one)
})
