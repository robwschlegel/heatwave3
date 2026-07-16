# hw3_pixel_clim() computes a single-pixel climatology directly from plain
# numeric vectors (no NetCDF I/O). It backs no R wrapper -- ts2clm3() always
# goes through the gridded NetCDF path -- so it needs a direct test to be
# exercised at all.

.jd <- function(date) as.integer(as.Date(date)) + 2440588L

synth_pixel <- function(start = "1990-01-01", end = "2009-12-31", seed = 1) {
  set.seed(seed)
  dates <- seq(as.Date(start), as.Date(end), by = "day")
  yday_frac <- as.integer(format(dates, "%j")) / 365
  temp <- 290 + 5 * sin(2 * pi * yday_frac) + stats::rnorm(length(dates), 0, 0.5)
  list(temp = temp, time_jd = .jd(dates), start = start, end = end)
}

test_that("hw3_pixel_clim computes a plausible seasonal climatology", {
  px <- synth_pixel()
  res <- hw3_pixel_clim(px$temp, px$time_jd, px$start, px$end)

  expect_named(res, c("seas", "thresh", "valid"))
  expect_true(res$valid)
  expect_length(res$seas, 366)
  expect_length(res$thresh, 366)
  expect_false(anyNA(res$seas))
  expect_false(anyNA(res$thresh))
  # thresh (90th pctile) must sit above seas (the mean) everywhere
  expect_true(all(res$thresh > res$seas))
  # seasonal cycle recovered within tolerance of the synthetic amplitude
  expect_equal(max(res$seas) - min(res$seas), 10, tolerance = 1)
})

test_that("hw3_pixel_clim respects detrend, smoothPercentile, and maxPadLength arguments", {
  px <- synth_pixel()

  res_raw <- hw3_pixel_clim(px$temp, px$time_jd, px$start, px$end,
                            smoothPercentile = FALSE)
  expect_true(res_raw$valid)
  expect_false(anyNA(res_raw$thresh))

  res_detrend <- hw3_pixel_clim(px$temp, px$time_jd, px$start, px$end,
                                detrend = TRUE)
  expect_true(res_detrend$valid)

  res_pad <- hw3_pixel_clim(px$temp, px$time_jd, px$start, px$end,
                            maxPadLength = 3L)
  expect_true(res_pad$valid)

  res_pctile10 <- hw3_pixel_clim(px$temp, px$time_jd, px$start, px$end, pctile = 10)
  expect_true(all(res_pctile10$thresh < res_pctile10$seas))
})

test_that("hw3_pixel_clim reports valid = FALSE for insufficient data", {
  dates <- seq(as.Date("2020-01-01"), as.Date("2020-01-10"), by = "day")
  temp <- rep(290, length(dates))
  res <- hw3_pixel_clim(temp, .jd(dates), "2020-01-01", "2020-01-10")
  expect_false(res$valid)
})
