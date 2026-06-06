# Tests for detect_event3 — name-stem outputs and the daily / protoEvent
# permutations. A climatology is built first under the same stem.

cp <- c("1982-01-01", "2011-12-31")

make_clim <- function() {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = cp, quiet = TRUE)
  list(sst = sst_file, stem = stem)
}

test_that("detect_event3 cold-spell intensities are negative (heatwaveR convention)", {
  fx <- make_clim()
  stem <- tempfile()
  ts2clm3(fx$sst, name = stem, climatologyPeriod = cp, pctile = 10, quiet = TRUE)
  detect_event3(fx$sst, name = stem, coldSpells = TRUE, quiet = TRUE)
  ev <- hw3_export(paste0(stem, "_events.nc"))
  # heatwaveR reports cold-spell intensities with the sign of the anomaly
  # retained, i.e. negative; ranking/severity use magnitude.
  expect_true(all(ev$intensity_mean < 0, na.rm = TRUE))
  expect_true(all(ev$intensity_max < 0, na.rm = TRUE))
  expect_true(all(ev$intensity_cumulative < 0, na.rm = TRUE))
  expect_true(all(ev$intensity_var >= 0, na.rm = TRUE))        # variance positive
  expect_true(all(ev$intensity_mean_abs > 0, na.rm = TRUE))    # absolute temps positive
})

test_that("detect_event3 cold-spell intensities match heatwaveR like-for-like", {
  skip_if_not_installed("heatwaveR")
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = cp, pctile = 10, quiet = TRUE)
  detect_event3(sst_file, name = stem, coldSpells = TRUE, quiet = TRUE)
  ev <- hw3_export(paste0(stem, "_events.nc"))

  # heatwaveR on the first pixel of the same grid
  gd <- hw3_read_sst(sst_file, "")
  nt <- gd$ntime
  df <- data.frame(t = as.Date(gd$time_days - 2440588L, origin = "1970-01-01"),
                   temp = gd$sst[1:nt])
  tsR <- heatwaveR::ts2clm(df, climatologyPeriod = cp, pctile = 10)
  R <- heatwaveR::detect_event(tsR, coldSpells = TRUE)$event

  H <- ev[ev$lon == gd$lon[1] & ev$lat == gd$lat[1], ]
  H <- H[order(H$date_start), ]
  n <- min(nrow(R), nrow(H))
  expect_gt(n, 0)
  for (col in c("intensity_mean", "intensity_max", "intensity_cumulative",
                "intensity_mean_relThresh", "intensity_max_relThresh",
                "intensity_cumulative_relThresh", "intensity_var",
                "rate_onset", "rate_decline")) {
    expect_equal(H[[col]][seq_len(n)], R[[col]][seq_len(n)], tolerance = 1e-3,
                 info = col)
  }
})

test_that("detect_event3 heatwave intensities and rates match heatwaveR", {
  skip_if_not_installed("heatwaveR")
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = cp, quiet = TRUE)
  detect_event3(sst_file, name = stem, quiet = TRUE)
  ev <- hw3_export(paste0(stem, "_events.nc"))

  gd <- hw3_read_sst(sst_file, "")
  nt <- gd$ntime
  df <- data.frame(t = as.Date(gd$time_days - 2440588L, origin = "1970-01-01"),
                   temp = gd$sst[1:nt])
  tsR <- heatwaveR::ts2clm(df, climatologyPeriod = cp)
  R <- heatwaveR::detect_event(tsR)$event

  H <- ev[ev$lon == gd$lon[1] & ev$lat == gd$lat[1], ]
  H <- H[order(H$date_start), ]
  n <- min(nrow(R), nrow(H))
  expect_gt(n, 0)
  for (col in c("intensity_mean", "intensity_max", "intensity_cumulative",
                "intensity_var", "rate_onset", "rate_decline")) {
    expect_equal(H[[col]][seq_len(n)], R[[col]][seq_len(n)], tolerance = 1e-3,
                 info = col)
  }
})

test_that("daily = 'none' (default) writes only the events file", {
  fx <- make_clim()
  out <- detect_event3(fx$sst, name = fx$stem, quiet = TRUE)

  expect_named(out, "events")
  expect_true(file.exists(paste0(fx$stem, "_events.nc")))
  expect_false(file.exists(paste0(fx$stem, "_events_daily.nc")))
  expect_equal(hw3_file_meta(out[["events"]])$product, "events")
})

test_that("daily = 'also' writes events + a self-consistent per-day file", {
  fx <- make_clim()
  out <- detect_event3(fx$sst, name = fx$stem, category = TRUE,
                       daily = "also", quiet = TRUE)

  expect_named(out, c("events", "daily"))
  expect_true(all(file.exists(out)))

  ev <- hw3_export(out[["events"]])
  dd <- hw3_export(out[["daily"]])
  expect_named(dd, c("lon", "lat", "t", "temp", "seas", "thresh",
                     "intensity", "event", "event_no", "category"))
  # daily category is NA off-threshold and 1-4 on-threshold
  expect_true(all(is.na(dd$category) | dd$category %in% 1:4))
  exceed <- !is.na(dd$temp) & !is.na(dd$thresh) & dd$temp > dd$thresh
  expect_equal(!is.na(dd$category), exceed)
  # event-days equal the total event duration
  expect_equal(sum(!is.na(dd$event_no)), sum(ev$duration))
})

test_that("daily = 'only' writes the per-day file and no events table", {
  fx <- make_clim()
  out <- detect_event3(fx$sst, name = fx$stem, daily = "only", quiet = TRUE)

  expect_named(out, "daily")
  expect_true(file.exists(paste0(fx$stem, "_events_daily.nc")))
  expect_false(file.exists(paste0(fx$stem, "_events.nc")))
  expect_equal(hw3_file_meta(out[["daily"]])$product, "daily")
})

test_that("protoEvent writes the heatwaveR-style proto series only", {
  fx <- make_clim()
  out <- detect_event3(fx$sst, name = fx$stem, protoEvent = TRUE, quiet = TRUE)

  expect_named(out, "protoevents")
  expect_true(file.exists(paste0(fx$stem, "_protoevents.nc")))
  expect_false(file.exists(paste0(fx$stem, "_events.nc")))
  expect_equal(hw3_file_meta(out[["protoevents"]])$product, "protoevents")

  proto <- hw3_export(out[["protoevents"]])
  expect_named(proto, c("lon", "lat", "t", "temp", "seas", "thresh",
                        "threshCriterion", "durationCriterion",
                        "event", "event_no"))
  expect_type(proto$event, "logical")
  # duration-criterion days are a subset of event days
  expect_true(all(proto$event[proto$durationCriterion]))
})

test_that("clim_file defaults to <name>_clim.nc but can be overridden", {
  fx <- make_clim()

  # default: picks up fx$stem_clim.nc automatically
  expect_no_error(detect_event3(fx$sst, name = tempfile(),
                                clim_file = paste0(fx$stem, "_clim.nc"),
                                quiet = TRUE))

  # missing default climatology -> informative error
  expect_error(detect_event3(fx$sst, name = tempfile(), quiet = TRUE),
               "Climatology file does not exist")
})

test_that("detect_event3 argument validation works", {
  fx <- make_clim()

  expect_error(detect_event3(fx$sst, name = fx$stem,
                             protoEvent = TRUE, daily = "also", quiet = TRUE),
               "mutually exclusive")
  expect_error(detect_event3(name = fx$stem), "file_in and name")
  expect_error(detect_event3(fx$sst, name = fx$stem, daily = "sometimes",
                             quiet = TRUE),
               "should be one of")
})

test_that("hw3_export round-trips every product to CSV and previews", {
  fx <- make_clim()
  detect_event3(fx$sst, name = fx$stem, category = TRUE, daily = "also",
                quiet = TRUE)
  detect_event3(fx$sst, name = fx$stem, protoEvent = TRUE, quiet = TRUE)

  files <- c(
    clim = paste0(fx$stem, "_clim.nc"),
    events = paste0(fx$stem, "_events.nc"),
    daily = paste0(fx$stem, "_events_daily.nc"),
    protoevents = paste0(fx$stem, "_protoevents.nc")
  )
  for (f in files) {
    df_whole <- hw3_export(f)
    expect_s3_class(df_whole, "data.frame")
    expect_equal(nrow(hw3_export(f, n = 3)), min(3, nrow(df_whole)))
    csv <- tempfile(fileext = ".csv")
    hw3_export(f, file_out = csv)
    expect_equal(nrow(utils::read.csv(csv)), nrow(df_whole))
  }
})
