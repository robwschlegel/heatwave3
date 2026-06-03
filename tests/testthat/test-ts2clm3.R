# Tests for ts2clm3 — name-stem output, console summary, and heatwaveR parity

test_that("ts2clm3 matches heatwaveR::ts2clm for a single pixel", {
  skip_if_not(file.exists("/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc"),
              "SWIO test data not available")

  sst_file <- "/Volumes/OceanData/OSTIA_East_Coast_MHW/SWIO_Jan1982-Dec2021.nc"
  stem <- tempfile()
  clim_file <- paste0(stem, "_clim.nc")
  on.exit(unlink(clim_file))

  ts2clm3(file_in = sst_file, name = stem,
          climatologyPeriod = c("1982-01-01", "2011-12-31"),
          lon_range = c(25.025, 25.075),
          lat_range = c(-33.975, -33.925),
          var_name = "analysed_sst", quiet = TRUE)

  expect_true(file.exists(clim_file))

  nc <- ncdf4::nc_open(clim_file)
  seas_hw3 <- ncdf4::ncvar_get(nc, "seas")
  thresh_hw3 <- ncdf4::ncvar_get(nc, "thresh")
  hw3_lon <- ncdf4::ncvar_get(nc, "lon")
  hw3_lat <- ncdf4::ncvar_get(nc, "lat")
  ncdf4::nc_close(nc)

  nc2 <- ncdf4::nc_open(sst_file)
  lon <- ncdf4::ncvar_get(nc2, "longitude")
  lat <- ncdf4::ncvar_get(nc2, "latitude")
  time_raw <- ncdf4::ncvar_get(nc2, "time")
  li <- which.min(abs(lon - hw3_lon[1]))
  lj <- which.min(abs(lat - hw3_lat[1]))
  sst <- ncdf4::ncvar_get(nc2, "analysed_sst", start = c(li, lj, 1), count = c(1, 1, -1))
  ncdf4::nc_close(nc2)

  dates <- as.Date("1981-01-01") + time_raw / 86400
  clim_r <- heatwaveR::ts2clm(data.frame(t = dates, temp = as.numeric(sst)),
                               climatologyPeriod = c("1982-01-01", "2011-12-31"),
                               clmOnly = TRUE)

  hw3_seas_px <- seas_hw3[, 1, 1]
  hw3_thresh_px <- thresh_hw3[, 1, 1]

  expect_equal(length(hw3_seas_px), 366)
  expect_equal(max(abs(hw3_seas_px - clim_r$seas)), 1e-04, tolerance = 1e-04)
  expect_equal(max(abs(hw3_thresh_px - clim_r$thresh)), 0, tolerance = 1e-04)
})

test_that("ts2clm3 input validation works", {
  expect_error(ts2clm3(), "file_in and name")
  expect_error(ts2clm3("nonexistent.nc", tempfile(), c("1982-01-01", "2011-12-31")),
               "does not exist")
  expect_error(ts2clm3("in.nc", tempfile(), "1982-01-01"),
               "length 2")
})

test_that("ts2clm3 writes <name>_clim.nc and is readable via hw3_export", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  out <- ts2clm3(sst_file, name = stem,
                 climatologyPeriod = c("1982-01-01", "2011-12-31"), quiet = TRUE)

  expect_equal(out, paste0(stem, "_clim.nc"))
  expect_true(file.exists(out))

  clim_df <- hw3_export(out)
  expect_s3_class(clim_df, "data.frame")
  expect_named(clim_df, c("lon", "lat", "doy", "seas", "thresh"))
  expect_equal(nrow(clim_df), 2 * 3 * 366)

  # auto-detected product
  expect_equal(hw3_file_meta(out)$product, "climatology")

  # CSV export round-trips
  csv_file <- tempfile(fileext = ".csv")
  hw3_export(out, file_out = csv_file)
  expect_equal(utils::read.csv(csv_file), clim_df, ignore_attr = TRUE)

  # portion read
  expect_equal(nrow(hw3_export(out, n = 10)), 10)
})

test_that("ts2clm3 prints a console summary unless quiet", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  out <- utils::capture.output(
    ts2clm3(sst_file, name = stem, climatologyPeriod = c("1982-01-01", "2011-12-31"))
  )
  joined <- paste(out, collapse = "\n")
  expect_match(joined, "Climatology written to")
  expect_match(joined, "hw3_export")

  # quiet suppresses the summary block
  stem2 <- tempfile()
  out2 <- utils::capture.output(
    ts2clm3(sst_file, name = stem2, climatologyPeriod = c("1982-01-01", "2011-12-31"),
            quiet = TRUE)
  )
  expect_false(any(grepl("Climatology written to", out2)))
})

test_that("hw3_export writes chunked climatology parquet", {
  skip_if_not_installed("arrow")

  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  ts2clm3(sst_file, name = stem,
          climatologyPeriod = c("1982-01-01", "2011-12-31"), quiet = TRUE)
  clim_file <- paste0(stem, "_clim.nc")

  parquet_file <- tempfile(fileext = ".parquet")
  hw3_export(clim_file, parquet_file, chunk_size = 100L)
  pq_df <- as.data.frame(arrow::read_parquet(parquet_file))

  expect_equal(pq_df, hw3_export(clim_file), ignore_attr = TRUE)
})

test_that("multi-file readers fail fast unless skipping is explicit", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  missing_file <- file.path(tempdir(), "missing_hw3_file.nc")

  expect_error(
    ts2clm3(file_in = c(sst_file, missing_file), name = tempfile(),
            climatologyPeriod = c("1982-01-01", "2011-12-31"), quiet = TRUE),
    "No such file|open"
  )

  expect_no_error(
    ts2clm3(file_in = c(missing_file, sst_file), name = tempfile(),
            climatologyPeriod = c("1982-01-01", "2011-12-31"),
            skip_bad_files = TRUE, quiet = TRUE)
  )
})

test_that("ts2clm3 preserves input temperature units in climatology output", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  ts2clm3(sst_file, name = stem,
          climatologyPeriod = c("1982-01-01", "2011-12-31"), quiet = TRUE)

  nc <- ncdf4::nc_open(paste0(stem, "_clim.nc"))
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
    ts2clm3(file_in = sst_copy, name = tempfile(),
            climatologyPeriod = c("1982-01-01", "2011-12-31"), quiet = TRUE),
    "Unsupported CF calendar"
  )
})
