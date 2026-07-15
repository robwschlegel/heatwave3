# Tests for detect_blob3()'s connected-component labelling, including the
# 4D (lon/lat/depth/time) depth-connectivity extension. No prior coverage
# existed for blob detection before this file; see
# dev/4d-water-column-mhw-proposal.qmd Phase 2.

# ---- label_components_cpp(): hand-built masks with known component counts ----

test_that("label_components_cpp reproduces 3D behaviour when nz = 1", {
  nx <- 3; ny <- 3; nz <- 1; nt <- 2
  mask <- array(0L, dim = c(nx, ny, nz, nt))
  mask[1, 1, 1, 1] <- 1L
  mask[3, 3, 1, 2] <- 1L
  labels <- label_components_cpp(as.integer(mask), nx, ny, nz, nt, FALSE)
  expect_equal(attr(labels, "n_components"), 2L)

  mask2 <- array(0L, dim = c(nx, ny, nz, nt))
  mask2[1, 1, 1, 1] <- 1L
  mask2[1, 1, 1, 2] <- 1L  # same lon/lat, adjacent time -> connects
  labels2 <- label_components_cpp(as.integer(mask2), nx, ny, nz, nt, FALSE)
  expect_equal(attr(labels2, "n_components"), 1L)
})

test_that("label_components_cpp connects voxels across adjacent depth indices", {
  mask <- array(0L, dim = c(2, 2, 3, 2))
  mask[1, 1, 1, 1] <- 1L
  mask[1, 1, 2, 1] <- 1L  # same lon/lat/time, depth index 1 and 2 (adjacent)
  labels <- label_components_cpp(as.integer(mask), 2, 2, 3, 2, FALSE)
  expect_equal(attr(labels, "n_components"), 1L)
})

test_that("label_components_cpp does not connect voxels across a depth gap", {
  mask <- array(0L, dim = c(2, 2, 3, 2))
  mask[1, 1, 1, 1] <- 1L
  mask[1, 1, 3, 1] <- 1L  # depth index 1 and 3 -- gap at index 2, not adjacent
  labels <- label_components_cpp(as.integer(mask), 2, 2, 3, 2, FALSE)
  expect_equal(attr(labels, "n_components"), 2L)
})

test_that("label_components_cpp's dateline wrap only applies to longitude, not depth", {
  mask <- array(0L, dim = c(3, 2, 3, 1))
  mask[1, 1, 1, 1] <- 1L  # lon index 1 (first), depth 1
  mask[3, 1, 1, 1] <- 1L  # lon index 3 (last), depth 1 -- wraps to the first voxel
  mask[1, 1, 3, 1] <- 1L  # lon index 1, depth index 3 (last) -- no depth wrap
  labels <- label_components_cpp(as.integer(mask), 3, 2, 3, 1, TRUE)
  expect_equal(attr(labels, "n_components"), 2L)
})

test_that("label_components_cpp connects a single blob chained across all 4 axes", {
  mask <- array(0L, dim = c(2, 2, 2, 2))
  mask[1, 1, 1, 1] <- 1L
  mask[2, 1, 1, 1] <- 1L  # +lon
  mask[2, 2, 1, 1] <- 1L  # +lat
  mask[2, 2, 2, 1] <- 1L  # +depth
  mask[2, 2, 2, 2] <- 1L  # +time
  labels <- label_components_cpp(as.integer(mask), 2, 2, 2, 2, FALSE)
  expect_equal(attr(labels, "n_components"), 1L)
})

# ---- blob_daily_summary_cpp(): hand-computed depth/volume fields ----

test_that("blob_daily_summary_cpp computes volume/depth fields correctly", {
  cell_area_lat <- c(100, 200)
  lons <- c(10, 20)
  lats <- c(-30, -20)
  depths <- c(1, 5, 10)
  layer_thickness_km <- c(0.002, 0.004, 0.006)

  labels <- c(1L, 1L)
  ix <- c(1L, 1L); iy <- c(1L, 1L); it <- c(1L, 1L); iz <- c(1L, 2L)
  delta <- c(1.0, 3.0)

  res <- blob_daily_summary_cpp(labels, ix, iy, it, delta, cell_area_lat, lons, lats, 1L,
                                iz, depths, layer_thickness_km, TRUE)

  expect_equal(res$volume_km3, 100 * 0.002 + 100 * 0.004)
  expect_equal(res$depth_min_m, 1)
  expect_equal(res$depth_max_m, 5)
  expect_equal(res$depth_at_max_delta_m, 5)  # delta = 3.0 (the max) was at depth index 2
  expect_equal(res$max_delta, 3)
  expect_equal(res$area_km2, 200)  # unaffected by depth logic
})

test_that("blob_daily_summary_cpp omits depth columns when has_depth is FALSE", {
  cell_area_lat <- c(100, 200)
  labels <- c(1L, 1L)
  ix <- c(1L, 1L); iy <- c(1L, 1L); it <- c(1L, 1L)
  delta <- c(1.0, 3.0)

  res <- blob_daily_summary_cpp(labels, ix, iy, it, delta, cell_area_lat,
                                c(10, 20), c(-30, -20), 1L,
                                integer(0), numeric(0), numeric(0), FALSE)
  expect_false("volume_km3" %in% names(res))
  expect_equal(res$area_km2, 200)
})

# ---- detect_blob3(): end-to-end, both 3D and depth-resolved ----

test_that("detect_blob3 is unaffected by depth support for ordinary 3D input", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)
  blobs <- detect_blob3(sst_file, paste0(stem, "_clim.nc"),
                        return = c("event", "daily", "voxel"))

  expect_true(nrow(blobs$event) > 0)
  expect_false(any(grepl("depth|[Vv]olume", names(blobs$event))))
  expect_false(any(grepl("depth|[Vv]olume", names(blobs$daily))))
  expect_false("depth" %in% names(blobs$voxel))
  expect_false(anyNA(blobs$event$cumI_km2_day))
})

test_that("detect_blob3 connects blobs through depth for depth-resolved input", {
  glorys_file <- system.file("extdata/glorys_depth_test.nc", package = "heatwave3")
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")

  stem <- tempfile()
  ts2clm3(glorys_file, name = stem, climatologyPeriod = c("1993-01-01", "2019-12-31"),
          depth_range = c(0, 10), quiet = TRUE)
  blobs <- detect_blob3(glorys_file, paste0(stem, "_clim.nc"),
                        return = c("event", "daily", "voxel"))

  expect_true(nrow(blobs$event) > 0)
  expect_true(all(c("depth_min_m", "depth_max_m", "totalVolume_km3", "peakVolume_km3",
                    "meanVolume_km3", "depthOfPeakIntensity_m") %in% names(blobs$event)))
  expect_true(all(c("volume_km3", "depth_min_m", "depth_max_m", "depth_at_max_delta_m")
                  %in% names(blobs$daily)))
  expect_true("depth" %in% names(blobs$voxel))
  expect_false(anyNA(blobs$event$depth_min_m))
  expect_false(anyNA(blobs$event$totalVolume_km3))
  expect_true(all(blobs$event$depth_min_m <= blobs$event$depth_max_m))
  expect_true(all(blobs$event$totalVolume_km3 > 0))

  # The whole point of 4D connectivity: most blobs should span more than one
  # depth level (voxels at adjacent depths get unioned into the same blob).
  expect_true(mean(blobs$event$depth_min_m != blobs$event$depth_max_m) > 0.5)
})

test_that("detect_blob3 rankBy accepts volume metrics only for depth-resolved input", {
  glorys_file <- system.file("extdata/glorys_depth_test.nc", package = "heatwave3")
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")

  stem <- tempfile()
  ts2clm3(glorys_file, name = stem, climatologyPeriod = c("1993-01-01", "2019-12-31"),
          depth_range = c(0, 10), quiet = TRUE)

  blobs <- detect_blob3(glorys_file, paste0(stem, "_clim.nc"),
                        topN = 5, rankBy = "totalVolume_km3", return = "event")
  expect_equal(nrow(blobs$event), 5)
  expect_false(is.unsorted(rev(blobs$event$totalVolume_km3)))

  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  stem3d <- tempfile()
  ts2clm3(sst_file, name = stem3d, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)
  expect_error(
    detect_blob3(sst_file, paste0(stem3d, "_clim.nc"), rankBy = "totalVolume_km3"),
    "rankBy must be one of"
  )
})

test_that("detect_blob3 minVoxels filters depth-resolved blobs", {
  glorys_file <- system.file("extdata/glorys_depth_test.nc", package = "heatwave3")
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")

  stem <- tempfile()
  ts2clm3(glorys_file, name = stem, climatologyPeriod = c("1993-01-01", "2019-12-31"),
          depth_range = c(0, 10), quiet = TRUE)

  all_blobs <- detect_blob3(glorys_file, paste0(stem, "_clim.nc"), minVoxels = 1, return = "event")
  filtered <- detect_blob3(glorys_file, paste0(stem, "_clim.nc"), minVoxels = 50, return = "event")
  expect_true(nrow(filtered$event) <= nrow(all_blobs$event))
})
