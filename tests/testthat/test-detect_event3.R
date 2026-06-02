# Tests for detect_event3 — compare against heatwaveR pixel by pixel

test_that("detect_event3 matches heatwaveR::detect_event for a single pixel", {
  skip_if_not(file.exists("/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc"),
              "SWIO test data not available")

  sst_file <- "/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc"
  clim_file <- tempfile(fileext = ".nc")
  event_file <- tempfile(fileext = ".nc")
  on.exit(unlink(c(clim_file, event_file)))

  # Run heatwave3
  ts2clm3(file_in = sst_file, file_out = clim_file,
          climatologyPeriod = c("1982-01-01", "2011-12-31"),
          lon_range = c(25.025, 25.075),
          lat_range = c(-33.975, -33.925),
          var_name = "analysed_sst")

  detect_event3(file_in = sst_file, clim_file = clim_file,
                file_out = event_file, var_name = "analysed_sst")

  expect_true(file.exists(event_file))

  # Read hw3 events
  nc <- ncdf4::nc_open(event_file)
  hw3_dur <- ncdf4::ncvar_get(nc, "duration")
  hw3_im <- ncdf4::ncvar_get(nc, "intensity_mean")
  hw3_ix <- ncdf4::ncvar_get(nc, "intensity_max")
  hw3_ic <- ncdf4::ncvar_get(nc, "intensity_cumulative")
  hw3_lon <- ncdf4::ncvar_get(nc, "lon")
  hw3_lat <- ncdf4::ncvar_get(nc, "lat")
  ncdf4::nc_close(nc)

  # Filter to first valid pixel
  nc_cl <- ncdf4::nc_open(clim_file)
  cl_lon <- ncdf4::ncvar_get(nc_cl, "lon")
  cl_lat <- ncdf4::ncvar_get(nc_cl, "lat")
  ncdf4::nc_close(nc_cl)

  px_lon <- cl_lon[1]; px_lat <- cl_lat[1]
  px1 <- which(abs(hw3_lon - px_lon) < 0.01 & abs(hw3_lat - px_lat) < 0.01)

  # Run heatwaveR on the same pixel
  nc2 <- ncdf4::nc_open(sst_file)
  lon <- ncdf4::ncvar_get(nc2, "longitude")
  lat <- ncdf4::ncvar_get(nc2, "latitude")
  time_raw <- ncdf4::ncvar_get(nc2, "time")
  li <- which.min(abs(lon - px_lon))
  lj <- which.min(abs(lat - px_lat))
  sst <- ncdf4::ncvar_get(nc2, "analysed_sst", start = c(li, lj, 1), count = c(1, 1, -1))
  ncdf4::nc_close(nc2)

  dates <- as.Date("1981-01-01") + time_raw / 86400
  clim_r <- heatwaveR::ts2clm(data.frame(t = dates, temp = as.numeric(sst)),
                               climatologyPeriod = c("1982-01-01", "2011-12-31"))
  ev_r <- heatwaveR::detect_event(clim_r)$event

  # Compare
  expect_equal(length(px1), nrow(ev_r))
  expect_equal(as.integer(hw3_dur[px1]), ev_r$duration)
  expect_equal(max(abs(hw3_im[px1] - ev_r$intensity_mean)), 0, tolerance = 2e-04)
  expect_equal(max(abs(hw3_ix[px1] - ev_r$intensity_max)), 0, tolerance = 2e-04)
})

test_that("detect_event3 input validation works", {
  expect_error(detect_event3(), "file_in")
})

test_that("detect_event3 writes a custom companion file", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  out_dir <- tempfile()
  dir.create(out_dir)
  clim_file <- file.path(out_dir, "clim.nc")
  event_file <- file.path(out_dir, "events.nc")
  csv_file <- file.path(out_dir, "custom_events.csv")

  ts2clm3(
    file_in = sst_file,
    file_out = clim_file,
    climatologyPeriod = c("1982-01-01", "2011-12-31")
  )

  event_df <- detect_event3(
    file_in = sst_file,
    clim_file = clim_file,
    file_out = event_file,
    save_file = csv_file,
    return_df = TRUE
  )

  expect_s3_class(event_df, "data.frame")
  expect_true(file.exists(csv_file))
  expect_gt(nrow(event_df), 0)

  csv_df <- utils::read.csv(csv_file, stringsAsFactors = FALSE)
  date_cols <- c("date_start", "date_peak", "date_end")
  csv_df[date_cols] <- lapply(csv_df[date_cols], as.Date)
  expect_equal(csv_df, event_df, ignore_attr = TRUE)

  nc <- ncdf4::nc_open(event_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_equal(ncdf4::ncatt_get(nc, "intensity_mean", "units")$value, "kelvin")
  expect_equal(ncdf4::ncatt_get(nc, "intensity_cumulative", "units")$value, "kelvin days")
  expect_equal(ncdf4::ncatt_get(nc, "rate_onset", "units")$value, "kelvin/day")
})

test_that("detect_event3 threshClim2 matches heatwaveR secondary event pass", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  skip_if_not_installed("heatwaveR")
  skip_if_not_installed("ncdf4")

  out_dir <- tempfile()
  dir.create(out_dir)
  clim_file <- file.path(out_dir, "clim.nc")
  event_file <- file.path(out_dir, "events.nc")
  mask_file <- file.path(out_dir, "thresh2.nc")

  ts2clm3(
    file_in = sst_file,
    file_out = clim_file,
    climatologyPeriod = c("1982-01-01", "2011-12-31"),
    lon_range = c(26.525, 26.525),
    lat_range = c(-34.125, -34.125)
  )

  nc <- ncdf4::nc_open(sst_file)
  lon <- ncdf4::ncvar_get(nc, "longitude")
  lat <- ncdf4::ncvar_get(nc, "latitude")
  time_raw <- ncdf4::ncvar_get(nc, "time")
  time_units <- ncdf4::ncatt_get(nc, "time", "units")$value
  time_calendar <- ncdf4::ncatt_get(nc, "time", "calendar")$value
  sst <- as.numeric(ncdf4::ncvar_get(
    nc, "analysed_sst",
    start = c(1, 1, 1),
    count = c(1, 1, -1)
  ))
  ncdf4::nc_close(nc)

  dates <- as.Date("1981-01-01") + time_raw / 86400
  mask_vec <- dates %in% as.Date(c(
    "1982-11-06", "1982-11-07",
    "1982-11-09", "1982-11-10", "1982-11-11"
  ))

  lon_dim <- ncdf4::ncdim_def("longitude", "degrees_east", lon)
  lat_dim <- ncdf4::ncdim_def("latitude", "degrees_north", lat)
  time_dim <- ncdf4::ncdim_def("time", time_units, time_raw)
  mask_var <- ncdf4::ncvar_def(
    "thresh2", "1", list(lon_dim, lat_dim, time_dim),
    missval = NA_real_, prec = "double"
  )
  nc_mask <- ncdf4::nc_create(mask_file, mask_var)
  if (!is.null(time_calendar) && !is.na(time_calendar)) {
    ncdf4::ncatt_put(nc_mask, "time", "calendar", time_calendar)
  }
  mask_grid <- array(0, dim = c(length(lon), length(lat), length(time_raw)))
  mask_grid[1, 1, ] <- as.numeric(mask_vec)
  ncdf4::ncvar_put(nc_mask, "thresh2", mask_grid)
  ncdf4::nc_close(nc_mask)

  ev_hw3 <- detect_event3(
    file_in = sst_file,
    clim_file = clim_file,
    file_out = event_file,
    threshClim2 = mask_file,
    threshClim2_var_name = "thresh2",
    minDuration2 = 2L,
    maxGap2 = 1L,
    return_df = TRUE
  )

  clim_r <- heatwaveR::ts2clm(
    data.frame(t = dates, temp = sst),
    climatologyPeriod = c("1982-01-01", "2011-12-31")
  )
  ev_r <- heatwaveR::detect_event(
    clim_r,
    threshClim2 = mask_vec,
    minDuration2 = 2L,
    maxGap2 = 1L
  )$event
  ev_r <- ev_r[!is.na(ev_r$event_no), , drop = FALSE]

  expect_equal(nrow(ev_hw3), nrow(ev_r))
  expect_equal(ev_hw3$date_start, ev_r$date_start)
  expect_equal(ev_hw3$date_peak, ev_r$date_peak)
  expect_equal(ev_hw3$date_end, ev_r$date_end)
  expect_equal(ev_hw3$duration, ev_r$duration)
  expect_equal(ev_hw3$intensity_max, ev_r$intensity_max, tolerance = 2e-04)
  expect_equal(
    ev_hw3$intensity_cumulative,
    ev_r$intensity_cumulative,
    tolerance = 2e-04
  )
})

test_that("Cold spell detection works", {
  skip_if_not(file.exists("/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc"),
              "SWIO test data not available")

  sst_file <- "/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc"
  clim_file <- tempfile(fileext = ".nc")
  event_file <- tempfile(fileext = ".nc")
  on.exit(unlink(c(clim_file, event_file)))

  ts2clm3(file_in = sst_file, file_out = clim_file,
          climatologyPeriod = c("1982-01-01", "2011-12-31"),
          lon_range = c(25.025, 25.075),
          lat_range = c(-33.975, -33.925),
          var_name = "analysed_sst",
          pctile = 10)

  detect_event3(file_in = sst_file, clim_file = clim_file,
                file_out = event_file, var_name = "analysed_sst",
                coldSpells = TRUE)

  nc <- ncdf4::nc_open(event_file)
  im <- ncdf4::ncvar_get(nc, "intensity_mean")
  ncdf4::nc_close(nc)

  # Cold spell intensities should be positive (negated)
  expect_true(all(im > 0, na.rm = TRUE))
})
