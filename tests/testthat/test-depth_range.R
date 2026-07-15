# Tests for depth_range (native depth-resolved climatology/events) and the
# legacy single-index depth squeeze, using the small 2x2 pixel x 5 depth-level
# GLORYS fixture (inst/extdata/glorys_depth_test.nc and its annual split in
# inst/extdata/glorys_depth_test_annual/). See dev/4d-depth-range-progress.md.

glorys_file <- system.file("extdata/glorys_depth_test.nc", package = "heatwave3")
glorys_dir  <- system.file("extdata/glorys_depth_test_annual", package = "heatwave3")
clim_period <- c("1993-01-01", "2019-12-31")

test_that("ts2clm3 depth_range reads all 5 levels with no NA and a real depth dimension", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")

  stem <- tempfile()
  ts2clm3(glorys_file, name = stem, climatologyPeriod = clim_period,
          depth_range = c(0, 10), quiet = TRUE)

  clim_df <- hw3_export(paste0(stem, "_clim.nc"))
  expect_true("depth" %in% names(clim_df))
  expect_equal(nrow(clim_df), 2 * 2 * 5 * 366)
  expect_equal(sort(unique(clim_df$depth)),
               c(0.494025, 1.541375, 2.645669, 3.819495, 5.078224),
               tolerance = 1e-5)
  expect_false(anyNA(clim_df$seas))
  expect_false(anyNA(clim_df$thresh))
})

test_that("depth and depth_range are mutually exclusive", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")

  expect_error(
    ts2clm3(glorys_file, name = tempfile(), climatologyPeriod = clim_period,
            depth = 2, depth_range = c(0, 10), quiet = TRUE),
    "only one of 'depth' or 'depth_range'"
  )
})

test_that("detect3 depth_range produces events with a matching depth column", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")

  stem <- tempfile()
  detect3(glorys_file, name = stem, climatologyPeriod = clim_period,
          depth_range = c(0, 10), category = TRUE, quiet = TRUE)

  ev <- hw3_export(paste0(stem, "_events.nc"))
  expect_true(nrow(ev) > 0)
  expect_true("depth" %in% names(ev))
  expect_equal(sort(unique(ev$depth)),
               c(0.494025, 1.541375, 2.645669, 3.819495, 5.078224),
               tolerance = 1e-5)
  expect_true("category" %in% names(ev))
})


test_that("legacy single-index depth squeeze matches its own level, not the surface", {
  # Regression test for the pre-existing bug documented in
  # dev/4d-depth-range-progress.md: detect_event3() used to always compare
  # against the surface (depth index 0) SST regardless of which depth index
  # ts2clm3() was squeezed to. Fixed by recording the resolved depth value on
  # the climatology file (as a global attribute for the 1-level squeeze case)
  # so detect_event3() can match its SST read to the same level.
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")

  stem3 <- tempfile()
  ts2clm3(glorys_file, name = stem3, climatologyPeriod = clim_period,
          depth = 3, quiet = TRUE)
  detect_event3(glorys_file, name = stem3, quiet = TRUE)
  ev3 <- hw3_export(paste0(stem3, "_events.nc"))

  stem0 <- tempfile()
  ts2clm3(glorys_file, name = stem0, climatologyPeriod = clim_period,
          depth = 0, quiet = TRUE)
  detect_event3(glorys_file, name = stem0, quiet = TRUE)
  ev0 <- hw3_export(paste0(stem0, "_events.nc"))

  # depth index 3 (3.819495 m) must differ from the surface (depth index 0) --
  # proves detect_event3() is not silently reading surface SST for k != 0.
  expect_false(isTRUE(all.equal(ev3, ev0)))

  # ...and must bit-match the known-correct depth_range single-level path.
  stem_dr <- tempfile()
  ts2clm3(glorys_file, name = stem_dr, climatologyPeriod = clim_period,
          depth_range = c(3.819495, 3.819495), quiet = TRUE)
  detect_event3(glorys_file, name = stem_dr, quiet = TRUE)
  ev_dr <- hw3_export(paste0(stem_dr, "_events.nc"))

  ord3  <- order(ev3$lon, ev3$lat, ev3$date_start)
  orddr <- order(ev_dr$lon, ev_dr$lat, ev_dr$date_start)
  common_cols <- setdiff(intersect(names(ev3), names(ev_dr)), "depth")

  ev3_sorted  <- ev3[ord3, common_cols]; rownames(ev3_sorted) <- NULL
  evdr_sorted <- ev_dr[orddr, common_cols]; rownames(evdr_sorted) <- NULL

  expect_equal(ev3_sorted, evdr_sorted)

  # Legacy squeeze keeps the pre-existing 3D schema: no depth column.
  expect_false("depth" %in% names(ev3))
})

test_that("multi-file directory ingestion with depth_range matches the single-file result", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")
  skip_if(glorys_dir == "", "glorys_depth_test_annual not available")

  stem_single <- tempfile()
  ts2clm3(glorys_file, name = stem_single, climatologyPeriod = clim_period,
          depth_range = c(0, 10), quiet = TRUE)

  stem_multi <- tempfile()
  ts2clm3(glorys_dir, name = stem_multi, climatologyPeriod = clim_period,
          depth_range = c(0, 10), quiet = TRUE)

  cl_single <- hw3_export(paste0(stem_single, "_clim.nc"))
  cl_multi  <- hw3_export(paste0(stem_multi, "_clim.nc"))

  ord_single <- order(cl_single$lon, cl_single$lat, cl_single$depth, cl_single$doy)
  ord_multi  <- order(cl_multi$lon, cl_multi$lat, cl_multi$depth, cl_multi$doy)
  cl_single <- cl_single[ord_single, ]; rownames(cl_single) <- NULL
  cl_multi  <- cl_multi[ord_multi, ]; rownames(cl_multi) <- NULL

  expect_equal(cl_single, cl_multi)
})
