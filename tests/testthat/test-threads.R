# Tests for thread-count management (getHW3threads/setHW3threads, wrapping
# hw3_get_threads()/hw3_set_threads() in heatwave3_init.cpp) and for forcing
# the multi-thread worker-pool branch of hw3::parallel_for() (hw3_omp.h).
# No other test calls these directly, and every other test uses the default
# n_threads = 1, so neither the pool branch nor hw3::hw3_default_threads()'s
# env-var branch is otherwise reached.

test_that("getHW3threads/setHW3threads round-trip and validate", {
  old <- getHW3threads()
  on.exit(setHW3threads(old), add = TRUE)

  expect_type(old, "integer")
  expect_gte(old, 1L)

  expect_equal(setHW3threads(3L), 3L)
  expect_equal(getHW3threads(), 3L)

  # 0 resets to the auto-detected default (50% of cores)
  expect_gte(setHW3threads(0L), 1L)
})

test_that("hw3_default_threads() honours R_HEATWAVE3_NUM_THREADS when set", {
  old <- getHW3threads()
  old_env <- Sys.getenv("R_HEATWAVE3_NUM_THREADS", unset = NA)
  on.exit({
    setHW3threads(old)
    if (is.na(old_env)) Sys.unsetenv("R_HEATWAVE3_NUM_THREADS")
    else Sys.setenv(R_HEATWAVE3_NUM_THREADS = old_env)
  }, add = TRUE)

  Sys.setenv(R_HEATWAVE3_NUM_THREADS = "3")
  expect_equal(setHW3threads(0L), 3L)

  Sys.unsetenv("R_HEATWAVE3_NUM_THREADS")
  expect_gte(setHW3threads(0L), 1L)
})

test_that("n_threads > 1 exercises the parallel_for worker-pool branch and matches n_threads = 1", {
  sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
  skip_if(sst_file == "", "sst_test.nc not available")
  cp <- c("1982-01-01", "2011-12-31")

  stem1 <- tempfile()
  ts2clm3(sst_file, name = stem1, climatologyPeriod = cp, n_threads = 1L, quiet = TRUE)
  stem2 <- tempfile()
  ts2clm3(sst_file, name = stem2, climatologyPeriod = cp, n_threads = 2L, quiet = TRUE)

  expect_equal(hw3_export(paste0(stem1, "_clim.nc")), hw3_export(paste0(stem2, "_clim.nc")))

  detect_event3(sst_file, name = stem1, clim_file = paste0(stem1, "_clim.nc"),
                n_threads = 1L, quiet = TRUE)
  detect_event3(sst_file, name = stem2, clim_file = paste0(stem2, "_clim.nc"),
                n_threads = 2L, quiet = TRUE)

  # parallel_for dispatches pixels dynamically for load balancing, so worker
  # threads finish (and so append to the merged event list) in a different
  # order than the single-threaded pixel-ascending order -- same events, only
  # the row order differs. Sort before comparing, as elsewhere in this suite.
  ev1 <- hw3_export(paste0(stem1, "_events.nc"))
  ev2 <- hw3_export(paste0(stem2, "_events.nc"))
  ord1 <- order(ev1$lon, ev1$lat, ev1$date_start)
  ord2 <- order(ev2$lon, ev2$lat, ev2$date_start)
  ev1 <- ev1[ord1, ]; rownames(ev1) <- NULL
  ev2 <- ev2[ord2, ]; rownames(ev2) <- NULL
  expect_equal(ev1, ev2)
})

test_that(".onLoad initializes fork safety", {
  # .onLoad() runs automatically once, at namespace-load time -- before covr's
  # line tracing is attached to that code -- so it's otherwise never credited
  # even though it genuinely executes on every package load. Call it directly
  # as an ordinary function; both arguments are ignored by the body.
  expect_no_error(heatwave3:::.onLoad("x", "heatwave3"))
})

test_that("hw3_version reports package metadata", {
  v <- hw3_version()
  expect_type(v, "list")
  expect_equal(v$package, "heatwave3")
  expect_match(v$version, "^[0-9]+\\.[0-9]+\\.[0-9]+$")
  expect_equal(v$cpp_standard, "C++17")
})
