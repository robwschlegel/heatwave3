# Tests for the second wave of depth_range support: depth-resolved
# daily/protoEvent output, hw3_export()'s streaming subset, category_daily3(),
# event_line3(), and geom_flame3()'s depth-aware grouping. Builds on the
# fixture/tests in test-depth_range.R -- see dev/4d-depth-range-progress.md.

glorys_file <- system.file("extdata/glorys_depth_test.nc", package = "heatwave3")
clim_period <- c("1993-01-01", "2019-12-31")
depth_vals <- c(0.494025, 1.541375, 2.645669, 3.819495, 5.078224)

test_that("detect3 depth_range + daily='also' writes a depth-resolved daily product", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")

  stem <- tempfile()
  detect3(glorys_file, name = stem, climatologyPeriod = clim_period,
          depth_range = c(0, 10), daily = "also", quiet = TRUE)

  df <- hw3_export(paste0(stem, "_events_daily.nc"))
  expect_true("depth" %in% names(df))
  expect_equal(nrow(df), 2L * 2L * 5L * 11688L)
  expect_equal(sort(unique(df$depth)), depth_vals, tolerance = 1e-5)
  expect_false(anyNA(df$temp))
})

test_that("detect3 depth_range + protoEvent writes a depth-resolved protoevents product", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")

  stem <- tempfile()
  detect3(glorys_file, name = stem, climatologyPeriod = clim_period,
          depth_range = c(0, 10), protoEvent = TRUE, quiet = TRUE)

  df <- hw3_export(paste0(stem, "_protoevents.nc"))
  expect_true("depth" %in% names(df))
  expect_true("threshCriterion" %in% names(df))
  expect_equal(sort(unique(df$depth)), depth_vals, tolerance = 1e-5)
})

test_that("hw3_export streaming subset is depth-aware for climatology and daily", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")

  stem <- tempfile()
  detect3(glorys_file, name = stem, climatologyPeriod = clim_period,
          depth_range = c(0, 10), daily = "also", quiet = TRUE)

  clim_whole  <- hw3_export(paste0(stem, "_clim.nc"))
  clim_sub    <- hw3_export(paste0(stem, "_clim.nc"), lon_range = c(25.1, 25.2),
                            depth_range = c(0, 3))
  clim_whole_sub <- clim_whole[clim_whole$lon >= 25.1 & clim_whole$lon <= 25.2 &
                                clim_whole$depth <= 3, ]
  o1 <- order(clim_sub$lon, clim_sub$lat, clim_sub$depth, clim_sub$doy)
  o2 <- order(clim_whole_sub$lon, clim_whole_sub$lat, clim_whole_sub$depth, clim_whole_sub$doy)
  a <- clim_sub[o1, ]; rownames(a) <- NULL
  b <- clim_whole_sub[o2, ]; rownames(b) <- NULL
  expect_equal(a, b)

  daily_sub <- hw3_export(paste0(stem, "_events_daily.nc"),
                          time_range = c("1993-01-01", "1993-01-10"),
                          depth_range = c(0, 3))
  expect_equal(nrow(daily_sub), 2L * 2L * 3L * 10L)
  expect_equal(sort(unique(daily_sub$depth)), depth_vals[1:3], tolerance = 1e-5)

  # n= preview still works and stays depth-aware
  prev <- hw3_export(paste0(stem, "_clim.nc"), n = 5)
  expect_equal(nrow(prev), 5L)
  expect_true("depth" %in% names(prev))
})

test_that("category_daily3 auto-surfaces a depth column for depth-resolved inputs", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")

  stem <- tempfile()
  detect3(glorys_file, name = stem, climatologyPeriod = clim_period,
          depth_range = c(0, 10), quiet = TRUE)

  res <- category_daily3(glorys_file, paste0(stem, "_clim.nc"), paste0(stem, "_events.nc"),
                         time_range = c("2020-01-01", "2020-01-10"))
  expect_true("depth" %in% names(res))
  expect_equal(nrow(res), 2L * 2L * 5L * 10L)
  expect_equal(sort(unique(res$depth)), depth_vals, tolerance = 1e-5)

  # Ordinary 3D clim/event pair: unaffected, no depth column.
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  stem3d <- tempfile()
  ts2clm3(sst_file, name = stem3d, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)
  detect_event3(sst_file, name = stem3d, quiet = TRUE)
  res3d <- category_daily3(sst_file, paste0(stem3d, "_clim.nc"), paste0(stem3d, "_events.nc"),
                           time_range = c("2010-01-01", "2010-01-05"))
  expect_false("depth" %in% names(res3d))
})

test_that("event_line3 supports a single depth level and errors clearly otherwise", {
  skip_if(glorys_file == "", "glorys_depth_test.nc not available")

  stem <- tempfile()
  detect3(glorys_file, name = stem, climatologyPeriod = clim_period,
          depth_range = c(0, 10), quiet = TRUE)

  p <- event_line3(glorys_file, paste0(stem, "_clim.nc"),
                   lon = 25.16667, lat = -34.91667, depth = 3,
                   event_file = paste0(stem, "_events.nc"),
                   start_date = "2015-01-01", end_date = "2016-12-31")
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "depth = ")

  expect_error(
    event_line3(glorys_file, paste0(stem, "_clim.nc"), lon = 25.16667, lat = -34.91667),
    "depth-resolved"
  )

  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  stem3d <- tempfile()
  ts2clm3(sst_file, name = stem3d, climatologyPeriod = c("1982-01-01", "2011-12-31"),
          quiet = TRUE)
  expect_error(
    event_line3(sst_file, paste0(stem3d, "_clim.nc"), lon = 26.525, lat = -34.125, depth = 10),
    "no depth dimension"
  )
  # 3D regression: still works with no depth argument.
  p3d <- event_line3(sst_file, paste0(stem3d, "_clim.nc"), lon = 26.525, lat = -34.125,
                     start_date = "2010-01-01", end_date = "2012-12-31")
  expect_s3_class(p3d, "ggplot")
})

test_that("geom_flame3 folds a mapped depth aesthetic into row grouping", {
  df <- data.frame(
    x = rep(1:10, 2),
    y = c(1, 1, 3, 3, 1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1, 1, 3, 3, 1, 1),
    y2 = 2,
    depth = rep(c(10, 20), each = 10)
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, y2 = y2, depth = depth)) +
    geom_flame3()
  built <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(sort(unique(built$group)), c(1, 2))

  # No depth column: grouping is untouched (still the ggplot2 default).
  df2 <- df[df$depth == 10, ]
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = x, y = y, y2 = y2)) + geom_flame3()
  built2 <- ggplot2::ggplot_build(p2)$data[[1]]
  expect_equal(unique(built2$group), -1)

  # An existing group aesthetic combines with depth rather than being replaced.
  df3 <- df
  df3$region <- rep(c("A", "B"), times = 10)
  p3 <- ggplot2::ggplot(df3, ggplot2::aes(x = x, y = y, y2 = y2, depth = depth, group = region)) +
    geom_flame3()
  built3 <- ggplot2::ggplot_build(p3)$data[[1]]
  expect_equal(sort(unique(built3$group)), 1:4)
})

test_that("geom_lolli3 renders unaffected with a depth column present", {
  df <- data.frame(
    x = as.Date("2020-01-01") + rep(seq(0, 90, by = 30), 2),
    y = c(1.2, 2.1, 0.8, 1.5, 0.9, 1.7, 2.3, 1.1),
    depth = rep(c(10, 20), each = 4)
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = factor(depth))) +
    geom_lolli3()
  built <- ggplot2::ggplot_build(p)$data[[1]]
  expect_equal(nrow(built), nrow(df))
})
