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

# ---- connectivity: which neighbours join a blob ----

# Number of components of a two-voxel mask at each connectivity setting
n_pair <- function(dims, a, b, wrap = FALSE) {
  mask <- array(0L, dim = dims)
  mask[matrix(a, 1)] <- 1L
  mask[matrix(b, 1)] <- 1L
  vapply(c(6L, 10L, 18L, 26L), function(conn)
    attr(label_components_cpp(as.integer(mask), dims[1], dims[2], dims[3], dims[4],
                              wrap, conn), "n_components"), integer(1))
}

test_that("connectivity sets which diagonal neighbours join (6/10/18/26)", {
  d3 <- c(3, 3, 1, 3)
  # Horizontal diagonal, same day: joined from 10 up
  expect_equal(n_pair(d3, c(1, 1, 1, 1), c(2, 2, 1, 1)), c(2L, 1L, 1L, 1L))
  # Space-time diagonal (lon + time): joined from 18 up
  expect_equal(n_pair(d3, c(1, 1, 1, 1), c(2, 1, 1, 2)), c(2L, 2L, 1L, 1L))
  # Corner (lon + lat + time): joined only at 26
  expect_equal(n_pair(d3, c(1, 1, 1, 1), c(2, 2, 1, 2)), c(2L, 2L, 2L, 1L))
  # Two steps apart: never joined
  expect_equal(n_pair(d3, c(1, 1, 1, 1), c(3, 1, 1, 1)), c(2L, 2L, 2L, 2L))

  d4 <- c(3, 3, 3, 3)
  # The lon/lat/time cube applies on each depth level, named as in 3D
  expect_equal(n_pair(d4, c(1, 1, 2, 1), c(2, 2, 2, 1)), c(2L, 1L, 1L, 1L))
  expect_equal(n_pair(d4, c(1, 1, 2, 1), c(2, 2, 2, 2)), c(2L, 2L, 2L, 1L))
  # Same pixel and day, adjacent level: always joined
  expect_equal(n_pair(d4, c(1, 1, 1, 1), c(1, 1, 2, 1)), c(1L, 1L, 1L, 1L))
  # Depth diagonals (with lon, or with time, or all four axes): never joined
  expect_equal(n_pair(d4, c(1, 1, 1, 1), c(2, 1, 2, 1)), c(2L, 2L, 2L, 2L))
  expect_equal(n_pair(d4, c(1, 1, 1, 1), c(1, 1, 2, 2)), c(2L, 2L, 2L, 2L))
  expect_equal(n_pair(d4, c(1, 1, 1, 1), c(2, 2, 2, 2)), c(2L, 2L, 2L, 2L))
})

test_that("dateline wrap joins diagonal neighbours across the boundary", {
  d <- c(4, 3, 1, 2)
  # First and last longitude columns, one row apart, same day
  expect_equal(n_pair(d, c(1, 1, 1, 1), c(4, 2, 1, 1), wrap = TRUE), c(2L, 1L, 1L, 1L))
  expect_equal(n_pair(d, c(1, 1, 1, 1), c(4, 2, 1, 1), wrap = FALSE), c(2L, 2L, 2L, 2L))
  # Across the boundary and one day later
  expect_equal(n_pair(d, c(4, 1, 1, 1), c(1, 1, 1, 2), wrap = TRUE), c(2L, 2L, 1L, 1L))
})

# Brute-force labelling: flood fill from each unlabelled voxel in raster
# order, so labels are numbered as label_components_cpp() numbers them.
label_bruteforce <- function(mask, conn, wrap) {
  d <- dim(mask)
  # lon/lat/time cube on one level, plus straight up and down in depth
  steps <- as.matrix(expand.grid(-1:1, -1:1, 0, -1:1))
  axes <- rowSums(steps != 0)
  keep <- axes == 1 | (conn == 10 & axes == 2 & steps[, 4] == 0) |
    (conn == 18 & axes == 2) | (conn == 26 & axes > 0)
  steps <- steps[keep, , drop = FALSE]
  if (d[3] > 1) steps <- rbind(steps, c(0, 0, -1, 0), c(0, 0, 1, 0))
  lab <- array(0L, d)
  k <- 0L
  for (v in which(mask == 1L)) {
    if (lab[v] > 0L) next
    k <- k + 1L
    lab[v] <- k
    queue <- v
    while (length(queue)) {
      p <- arrayInd(queue[1], d); queue <- queue[-1]
      nb <- sweep(steps, 2, p, "+")
      if (wrap) nb[, 1] <- (nb[, 1] - 1) %% d[1] + 1
      ok <- apply(nb, 1, function(x) all(x >= 1 & x <= d))
      nb <- nb[ok, , drop = FALSE]
      if (!nrow(nb)) next
      idx <- 1 + colSums((t(nb) - 1) * cumprod(c(1, d[-4])))
      idx <- idx[mask[idx] == 1L & lab[idx] == 0L]
      lab[idx] <- k
      queue <- c(queue, idx)
    }
  }
  lab
}

test_that("label_components_cpp matches a brute-force labelling at every connectivity", {
  set.seed(42)
  for (d in list(c(6, 5, 1, 7), c(4, 4, 3, 5))) {
    for (p in c(0.15, 0.35)) {
      mask <- array(as.integer(runif(prod(d)) < p), d)
      for (conn in c(6L, 10L, 18L, 26L)) for (wrap in c(FALSE, TRUE)) {
        got <- label_components_cpp(as.integer(mask), d[1], d[2], d[3], d[4], wrap, conn)
        ref <- label_bruteforce(mask, conn, wrap)
        expect_identical(as.vector(got), as.vector(ref),
                         info = sprintf("dims %s, p = %g, conn = %d, wrap = %s",
                                        paste(d, collapse = "x"), p, conn, wrap))
        expect_equal(attr(got, "n_components"), max(ref))
      }
    }
  }
})

test_that("label_components_cpp rejects other connectivity values", {
  expect_error(label_components_cpp(rep(0L, 8), 2, 2, 1, 2, FALSE, 8L),
               "connectivity must be 6, 10, 18 or 26")
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

test_that("detect_blob3 connectivity only merges blobs, and 6 is the default", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")

  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)
  clim <- paste0(stem, "_clim.nc")
  run <- function(...) detect_blob3(sst_file, clim, return = "voxel", ...)$voxel
  key <- function(v) paste(v$lon, v$lat, v$t_idx)

  expect_identical(run(), run(connectivity = 6))
  v6 <- run(connectivity = 6)
  for (conn in c(10, 18, 26)) {
    v <- run(connectivity = conn)
    # Same exceedance voxels; each 6-connected blob lies inside one blob
    expect_setequal(key(v), key(v6))
    merged <- v$blob[match(key(v6), key(v))]
    expect_true(all(tapply(merged, v6$blob, function(x) length(unique(x))) == 1))
    expect_lte(length(unique(v$blob)), length(unique(v6$blob)))
  }

  expect_error(detect_blob3(sst_file, clim, connectivity = 8),
               "connectivity must be one of 6, 10, 18 or 26")
})
