# detect_blob3(): v1 spatial-events labeller. The mask is the Hobday
# duration-filtered per-pixel event flag from detect_event3(daily = "also").
# Tests cover the cell-area integral, mask-from-event reshape, 6-connectivity
# labelling, the C++ metric reducer on a hand-computable case, sign retention /
# magnitude ranking, and an end-to-end run on the bundled fixture.

test_that("cell area: exact spherical-cap integrates to 4 pi R^2", {
  lat <- seq(-90 + 0.125, 90 - 0.125, by = 0.25)   # 720 rows tiling the sphere
  nlon <- 1440
  R <- 6371
  target <- 4 * pi * R^2

  a_exact <- heatwave3:::.hw3_cell_area_lat(lat, 0.25, 0.25, "exact")
  expect_equal(nlon * sum(a_exact), target, tolerance = 1e-9)  # relative

  a_coslat <- heatwave3:::.hw3_cell_area_lat(lat, 0.25, 0.25, "coslat")
  expect_lt(abs(nlon * sum(a_coslat) - target) / target, 1e-3) # approx, <0.1%
})

test_that("mask comes from the event flag, not exceedance", {
  # pixel-major event[px*ntime+t], px = i*nlat + j ; nlon=2, nlat=1, ntime=2
  event <- c(1L, 0L,   # px0 (i0): event on t0 only
             0L, 1L)   # px1 (i1): event on t1 only
  mask <- heatwave3:::hw3_blob_mask_from_event(event, 2L, 1L, 2L)
  # column-major v = i + nlon*j + nlon*nlat*t
  expect_equal(mask, c(1L, 0L, 0L, 1L))  # (i0,t0) and (i1,t1) only
})

test_that("label_components_3d_cpp matches a hand-built oracle (6-connectivity)", {
  nx <- 3; ny <- 3; nt <- 1
  mask <- integer(nx * ny * nt)
  idx <- function(i, j, k) i + nx * j + nx * ny * k + 1L  # 0-based i,j,k
  mask[idx(0, 0, 0)] <- 1L
  mask[idx(1, 0, 0)] <- 1L      # adjacent to (0,0) -> same blob
  mask[idx(2, 1, 0)] <- 1L      # corner-touch only -> separate blob
  lab <- heatwave3:::label_components_3d_cpp(mask, nx, ny, nt, FALSE)
  expect_equal(attr(lab, "n_components"), 2L)
  expect_equal(lab[idx(0, 0, 0)], lab[idx(1, 0, 0)])
  expect_false(lab[idx(0, 0, 0)] == lab[idx(2, 1, 0)])
})

test_that("hw3_blob_reduce computes metrics correctly on a synthetic blob", {
  nlon <- 2L; nlat <- 2L; ntime <- 2L
  ncell <- nlon * nlat
  lon <- c(0, 10); lat <- c(0, 20)
  cell_area_lat <- c(10, 20)                       # constant per latitude row
  # per-(pixel, time) arrays: temp[px*ntime+t], seas = 0, thresh = 1.
  temp <- numeric(ncell * ntime)
  seas <- numeric(ncell * ntime)
  thresh <- rep(1, ncell * ntime)
  setT <- function(i, j, t, val) temp[(i * nlat + j) * ntime + t + 1L] <<- val
  setT(0, 0, 0, 2)   # delta 2  (j=0, area 10)
  setT(1, 0, 0, 3)   # delta 3
  setT(0, 0, 1, 5)   # delta 5

  labels <- integer(ncell * ntime)
  v <- function(i, j, t) i + nlon * j + nlon * nlat * t + 1L
  labels[v(0, 0, 0)] <- 1L; labels[v(1, 0, 0)] <- 1L; labels[v(0, 0, 1)] <- 1L

  red <- heatwave3:::hw3_blob_reduce(labels, 1L, temp, seas, thresh,
                                     lon, lat, cell_area_lat,
                                     nlon, nlat, ntime, TRUE)
  e <- red$event
  expect_equal(e$n_voxels, 3L)
  expect_equal(e$duration_days, 2L)
  expect_equal(e$cumI_km2_d, 100)              # 2*10 + 3*10 + 5*10
  expect_equal(e$volume_km2_d, 30)
  expect_equal(e$total_area_km2, 20)           # two unique cells, area 10 each
  expect_equal(e$peak_area_km2, 20)
  expect_equal(e$mean_area_km2, 15)            # 30 / 2 active days
  expect_equal(e$mean_intensity, 100 / 30)
  expect_equal(e$max_intensity, 5)             # signed extreme
  expect_equal(e$peak_severity, 5)             # diff = 1, sev = delta
  expect_equal(e$frac_strong, 0.5)             # peak day: delta 2 -> Strong
  expect_equal(e$frac_severe, 0.5)             #           delta 3 -> Severe
  expect_equal(e$frac_moderate, 0)
  expect_equal(e$frac_extreme, 0)
  expect_equal(e$centroid_lon, 5, tolerance = 1e-6)   # great-circle of (0,0),(10,0)
  expect_equal(e$centroid_lat, 0, tolerance = 1e-6)
  expect_equal(length(red$voxel$event_no), 3L)
})

test_that("cold-spell metrics retain sign; ranking uses magnitude", {
  nlon <- 2L; nlat <- 1L; ntime <- 1L
  lon <- c(0, 10); lat <- 0; cell_area_lat <- 10
  seas <- numeric(nlon * nlat * ntime)
  thresh <- rep(-1, nlon * nlat * ntime)           # cold threshold below seas
  temp <- numeric(nlon * nlat * ntime)
  temp[(0 * nlat + 0) * ntime + 1L] <- -2          # blob 1: delta -2
  temp[(1 * nlat + 0) * ntime + 1L] <- -5          # blob 2: delta -5
  labels <- c(1L, 2L)
  red <- heatwave3:::hw3_blob_reduce(labels, 2L, temp, seas, thresh,
                                     lon, lat, cell_area_lat, nlon, nlat, ntime, FALSE)
  e <- red$event
  expect_true(all(e$cumI_km2_d < 0))               # sign retained
  expect_equal(e$max_intensity[e$event_no == 2], -5)
  expect_true(all(e$peak_severity > 0))            # severity positive both polarities
  ord <- order(-abs(e$cumI_km2_d))
  expect_equal(e$event_no[ord][1], 2L)             # |cumI| of blob 2 > blob 1
})

test_that("detect_blob3 runs end-to-end from a daily product (v1 vocabulary)", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)
  detect_event3(sst_file, name = stem, daily = "also", quiet = TRUE)
  daily_file <- paste0(stem, "_events_daily.nc")

  res <- detect_blob3(daily_file, cellAreaMethod = "exact",
                      return = c("event", "daily"))
  expect_true(all(c("event", "daily") %in% names(res)))
  cols <- c("event_no", "date_start", "date_end", "date_peak", "duration_days",
            "n_voxels", "peak_area_km2", "mean_area_km2", "total_area_km2",
            "volume_km2_d", "cumI_km2_d", "mean_intensity", "max_intensity",
            "peak_severity", "frac_moderate", "frac_strong", "frac_severe",
            "frac_extreme", "centroid_lon", "centroid_lat",
            "bbox_lon_min", "bbox_lon_max", "bbox_lat_min", "bbox_lat_max", "rank")
  expect_named(res$event, cols)
  if (nrow(res$event) > 0) {
    expect_s3_class(res$event$date_peak, "Date")
    expect_true(all(res$event$peak_area_km2 > 0))
    expect_true(all(res$event$volume_km2_d >= res$event$peak_area_km2))
    expect_equal(res$event$rank, seq_len(nrow(res$event)))
  }

  res2 <- detect_blob3(daily_file, rankBy = "volume_km2_day", return = "event")
  expect_true("event" %in% names(res2))
})

test_that("detect_blob3 validates inputs", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  stem <- tempfile()
  ts2clm3(sst_file, name = stem, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)
  detect_event3(sst_file, name = stem, daily = "also", quiet = TRUE)
  daily_file <- paste0(stem, "_events_daily.nc")

  expect_error(detect_blob3(daily_file, connectivity = 26L), "connectivity = 6")
  expect_error(detect_blob3(daily_file, rankBy = "nonsense"), "rankBy must be")
  expect_error(detect_blob3("nope.nc"), "does not exist")
  # an events file (not a daily product) lacks the gridded event variable
  expect_error(detect_blob3(paste0(stem, "_events.nc")))
})
