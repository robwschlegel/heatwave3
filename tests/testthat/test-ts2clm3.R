# Tests for ts2clm3 — compare against heatwaveR pixel by pixel

test_that("ts2clm3 matches heatwaveR::ts2clm for a single pixel", {
  skip_if_not(file.exists("/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc"),
              "SWIO test data not available")

  sst_file <- "/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc"
  clim_file <- tempfile(fileext = ".nc")
  on.exit(unlink(clim_file))

  ts2clm3(file_in = sst_file, file_out = clim_file,
          climatologyPeriod = c("1982-01-01", "2011-12-31"),
          lon_range = c(25.025, 25.075),
          lat_range = c(-33.975, -33.925),
          var_name = "analysed_sst")

  expect_true(file.exists(clim_file))

  # Read hw3 output
  nc <- ncdf4::nc_open(clim_file)
  seas_hw3 <- ncdf4::ncvar_get(nc, "seas")
  thresh_hw3 <- ncdf4::ncvar_get(nc, "thresh")
  hw3_lon <- ncdf4::ncvar_get(nc, "lon")
  hw3_lat <- ncdf4::ncvar_get(nc, "lat")
  ncdf4::nc_close(nc)

  # Read pixel for heatwaveR
  nc2 <- ncdf4::nc_open(sst_file)
  lon <- ncdf4::ncvar_get(nc2, "longitude")
  lat <- ncdf4::ncvar_get(nc2, "latitude")
  time_raw <- ncdf4::ncvar_get(nc2, "time")
  li <- which.min(abs(lon - hw3_lon[1]))
  lj <- which.min(abs(lat - hw3_lat[1]))
  sst <- ncdf4::ncvar_get(nc2, "analysed_sst", start = c(li, lj, 1), count = c(1, 1, -1))
  ncdf4::nc_close(nc2)

  dates <- as.Date("1981-01-01") + time_raw / 86400

  # Run heatwaveR
  clim_r <- heatwaveR::ts2clm(data.frame(t = dates, temp = as.numeric(sst)),
                               climatologyPeriod = c("1982-01-01", "2011-12-31"),
                               clmOnly = TRUE)

  # Compare — NetCDF dims are [doy, lat, lon] in R
  hw3_seas_px <- seas_hw3[, 1, 1]
  hw3_thresh_px <- thresh_hw3[, 1, 1]

  expect_equal(length(hw3_seas_px), 366)
  expect_equal(max(abs(hw3_seas_px - clim_r$seas)), 1e-04, tolerance = 1e-04)
  expect_equal(max(abs(hw3_thresh_px - clim_r$thresh)), 0, tolerance = 1e-04)
})

test_that("ts2clm3 input validation works", {
  expect_error(ts2clm3(), "file_in")
  expect_error(ts2clm3("nonexistent.nc", "out.nc", c("1982-01-01", "2011-12-31")),
               "does not exist")
  expect_error(ts2clm3("in.nc", "out.nc", "1982-01-01"),
               "length 2")
  expect_error(ts2clm3("in.nc", "out.nc", c("1982-01-01", "2011-12-31"),
                       save_file = "custom.txt"),
               "save_file")
})

test_that("ts2clm3 can return a data frame and write a custom companion file", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  out_dir <- tempfile()
  dir.create(out_dir)
  clim_file <- file.path(out_dir, "clim.nc")
  csv_file <- file.path(out_dir, "custom_climatology.csv")

  clim_df <- ts2clm3(
    file_in = sst_file,
    file_out = clim_file,
    climatologyPeriod = c("1982-01-01", "2011-12-31"),
    save_file = csv_file,
    return_df = TRUE
  )

  expect_s3_class(clim_df, "data.frame")
  expect_named(clim_df, c("lon", "lat", "doy", "seas", "thresh"))
  expect_true(file.exists(csv_file))
  expect_equal(nrow(clim_df), 2 * 3 * 366)

  csv_df <- utils::read.csv(csv_file)
  expect_equal(csv_df, clim_df, ignore_attr = TRUE)
})

test_that("hw3_export writes chunked climatology parquet", {
  skip_if_not_installed("arrow")

  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  out_dir <- tempfile()
  dir.create(out_dir)
  clim_file <- file.path(out_dir, "clim.nc")
  parquet_file <- file.path(out_dir, "clim.parquet")

  clim_df <- ts2clm3(
    file_in = sst_file,
    file_out = clim_file,
    climatologyPeriod = c("1982-01-01", "2011-12-31"),
    return_df = TRUE
  )

  hw3_export(clim_file, parquet_file, type = "clim", chunk_size = 100L)
  pq_df <- as.data.frame(arrow::read_parquet(parquet_file))

  expect_equal(pq_df, clim_df, ignore_attr = TRUE)
})

test_that("multi-file readers fail fast unless skipping is explicit", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  missing_file <- file.path(tempdir(), "missing_hw3_file.nc")
  clim_file <- tempfile(fileext = ".nc")

  expect_error(
    ts2clm3(
      file_in = c(sst_file, missing_file),
      file_out = clim_file,
      climatologyPeriod = c("1982-01-01", "2011-12-31")
    ),
    "No such file|open"
  )

  expect_no_error(
    ts2clm3(
      file_in = c(missing_file, sst_file),
      file_out = clim_file,
      climatologyPeriod = c("1982-01-01", "2011-12-31"),
      skip_bad_files = TRUE
    )
  )
})

test_that("ts2clm3 preserves input temperature units in climatology output", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  clim_file <- tempfile(fileext = ".nc")
  ts2clm3(
    file_in = sst_file,
    file_out = clim_file,
    climatologyPeriod = c("1982-01-01", "2011-12-31")
  )

  nc <- ncdf4::nc_open(clim_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  expect_equal(ncdf4::ncatt_get(nc, "seas", "units")$value, "kelvin")
  expect_equal(ncdf4::ncatt_get(nc, "thresh", "units")$value, "kelvin")
})

test_that("ts2clm3 rejects unsupported CF calendars explicitly", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  sst_copy <- tempfile(fileext = ".nc")
  file.copy(sst_file, sst_copy)

  nc <- ncdf4::nc_open(sst_copy, write = TRUE)
  ncdf4::ncatt_put(nc, "time", "calendar", "noleap")
  ncdf4::nc_close(nc)

  expect_error(
    ts2clm3(
      file_in = sst_copy,
      file_out = tempfile(fileext = ".nc"),
      climatologyPeriod = c("1982-01-01", "2011-12-31")
    ),
    "Unsupported CF calendar"
  )
})
