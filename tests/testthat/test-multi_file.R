# Multi-file ingestion paths (hw3_detect_events_multi, hw3_category_daily_multi,
# hw3_read_sst_multi in heatwave3_init.cpp) are only reached when file_in/
# sst_file has length > 1 or is a directory. No other test file exercises that
# branch for detect_event3()/category_daily3()/detect_blob3() -- ts2clm3()'s
# multi-file path is already covered in test-depth_range.R. Uses the annual
# split of the small GLORYS fixture at inst/extdata/glorys_depth_test_annual/
# (one file per year, 1993-2024) alongside the single merged file
# inst/extdata/glorys_depth_test.nc.

glorys_file <- system.file("extdata/glorys_depth_test.nc", package = "heatwave3")
glorys_dir  <- system.file("extdata/glorys_depth_test_annual", package = "heatwave3")
cp <- c("1993-01-01", "2019-12-31")

test_that("detect_event3 multi-file directory ingestion matches the single-file result", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")
  skip_if(glorys_dir == "", "glorys_depth_test_annual not available")

  stem <- tempfile()
  ts2clm3(glorys_file, name = stem, climatologyPeriod = cp, depth = 0, quiet = TRUE)
  clim_file <- paste0(stem, "_clim.nc")

  stem_single <- tempfile()
  detect_event3(glorys_file, name = stem_single, clim_file = clim_file, quiet = TRUE)

  stem_multi <- tempfile()
  detect_event3(glorys_dir, name = stem_multi, clim_file = clim_file, quiet = TRUE)

  ev_single <- hw3_export(paste0(stem_single, "_events.nc"))
  ev_multi  <- hw3_export(paste0(stem_multi, "_events.nc"))

  ord_s <- order(ev_single$lon, ev_single$lat, ev_single$date_start)
  ord_m <- order(ev_multi$lon, ev_multi$lat, ev_multi$date_start)
  ev_single <- ev_single[ord_s, ]; rownames(ev_single) <- NULL
  ev_multi  <- ev_multi[ord_m, ];  rownames(ev_multi)  <- NULL

  expect_equal(ev_single, ev_multi)
})

test_that("detect_event3 accepts an explicit vector of files, and skip_bad_files controls failure", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")
  skip_if(glorys_dir == "", "glorys_depth_test_annual not available")

  stem <- tempfile()
  ts2clm3(glorys_file, name = stem, climatologyPeriod = cp, depth = 0, quiet = TRUE)
  clim_file <- paste0(stem, "_clim.nc")

  files <- sort(list.files(glorys_dir, pattern = "\\.nc$", full.names = TRUE))
  bad_file <- file.path(tempdir(), "bad_glorys.nc")
  file.create(bad_file)
  on.exit(unlink(bad_file))

  expect_error(
    detect_event3(c(files, bad_file), name = tempfile(), clim_file = clim_file, quiet = TRUE),
    "NetCDF|open|Invalid"
  )

  expect_no_error(
    detect_event3(c(bad_file, files), name = tempfile(), clim_file = clim_file,
                  skip_bad_files = TRUE, quiet = TRUE)
  )
})

# category_daily3's multi-file reader applies the same time window to every
# file passed in (it is designed for one-time-step-per-day archives where the
# caller pre-filters to the files covering time_range -- see its @param
# sst_file docs). A file with no overlap throws unless skip_bad_files = TRUE.
# A window spanning the 1999/2000 file boundary exercises the genuine
# cross-file merge; skip_bad_files is exercised separately below.
test_that("category_daily3 multi-file SST source merges across a file boundary and matches the single-file result", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")
  skip_if(glorys_dir == "", "glorys_depth_test_annual not available")

  stem <- tempfile()
  ts2clm3(glorys_file, name = stem, climatologyPeriod = cp, depth = 0, quiet = TRUE)
  detect_event3(glorys_file, name = stem, quiet = TRUE)
  clim_file  <- paste0(stem, "_clim.nc")
  event_file <- paste0(stem, "_events.nc")

  win <- c("1999-12-20", "2000-01-10")
  res_single <- category_daily3(glorys_file, clim_file, event_file, time_range = win)

  files <- sort(list.files(glorys_dir, pattern = "\\.nc$", full.names = TRUE))
  boundary_files <- files[grepl("1999\\.nc$|2000\\.nc$", files)]
  res_multi <- category_daily3(boundary_files, clim_file, event_file, time_range = win)

  ord_s <- order(res_single$lon, res_single$lat, res_single$t)
  ord_m <- order(res_multi$lon, res_multi$lat, res_multi$t)
  res_single <- res_single[ord_s, ]; rownames(res_single) <- NULL
  res_multi  <- res_multi[ord_m, ];  rownames(res_multi)  <- NULL

  expect_equal(res_single, res_multi)
})

test_that("category_daily3 skip_bad_files skips both unreadable and out-of-window files", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")
  skip_if(glorys_dir == "", "glorys_depth_test_annual not available")

  stem <- tempfile()
  ts2clm3(glorys_file, name = stem, climatologyPeriod = cp, depth = 0, quiet = TRUE)
  detect_event3(glorys_file, name = stem, quiet = TRUE)
  clim_file  <- paste0(stem, "_clim.nc")
  event_file <- paste0(stem, "_events.nc")
  win <- c("2000-06-01", "2000-08-01")

  # Passing the whole directory fails by default: most files have no time
  # step inside win, and skip_bad_files defaults to FALSE.
  expect_error(
    category_daily3(glorys_dir, clim_file, event_file, time_range = win),
    "No time steps found"
  )

  # With skip_bad_files = TRUE, the 31 non-overlapping files are silently
  # skipped and the result matches the single-file read for the same window.
  res_single <- category_daily3(glorys_file, clim_file, event_file, time_range = win)
  res_dir <- category_daily3(glorys_dir, clim_file, event_file, time_range = win,
                             skip_bad_files = TRUE)
  ord_s <- order(res_single$lon, res_single$lat, res_single$t)
  ord_d <- order(res_dir$lon, res_dir$lat, res_dir$t)
  expect_equal(res_single[ord_s, ], res_dir[ord_d, ], ignore_attr = TRUE)

  # A genuinely unreadable file is skipped the same way.
  bad_file <- file.path(tempdir(), "bad_glorys2.nc")
  file.create(bad_file)
  on.exit(unlink(bad_file))
  files <- sort(list.files(glorys_dir, pattern = "\\.nc$", full.names = TRUE))
  boundary_files <- files[grepl("2000\\.nc$", files)]

  expect_error(
    category_daily3(c(bad_file, boundary_files), clim_file, event_file, time_range = win),
    "NetCDF|Unknown file format"
  )
  expect_no_error(
    category_daily3(c(bad_file, boundary_files), clim_file, event_file, time_range = win,
                    skip_bad_files = TRUE)
  )
})

test_that("detect_blob3 accepts a directory of daily files (hw3_read_sst_multi)", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")
  skip_if(glorys_dir == "", "glorys_depth_test_annual not available")

  stem <- tempfile()
  ts2clm3(glorys_file, name = stem, climatologyPeriod = cp, depth_range = c(0, 10),
          quiet = TRUE)
  clim_file <- paste0(stem, "_clim.nc")

  single <- detect_blob3(glorys_file, clim_file, return = "event")
  multi  <- detect_blob3(glorys_dir, clim_file, return = "event")

  expect_equal(nrow(single$event), nrow(multi$event))
  expect_equal(sort(single$event$cumI_km2_day), sort(multi$event$cumI_km2_day),
               tolerance = 1e-6)
})
