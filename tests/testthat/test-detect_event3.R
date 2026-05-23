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
