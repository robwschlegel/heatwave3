# Regenerate the R reference outputs ("golden" files) that the Python tests
# compare against. Run from the repository root after any change to the C++
# core or the R/Python front-ends that alters outputs:
#
#   Rscript python/tests/make_golden.R
#
# Each data.frame is saved exactly: <name>.json holds the column names, types
# and any character columns; <name>.bin.gz holds the numeric, integer, logical
# and Date columns as little-endian binary, in column order. Integer and
# logical NA are INT_MIN; Dates are days since 1970-01-01.

suppressMessages(devtools::load_all(".", quiet = TRUE))
out <- file.path("python", "tests", "data", "golden")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(out, full.names = TRUE))

save_golden <- function(df, name) {
  stopifnot(is.data.frame(df))
  cols <- lapply(names(df), function(nm) {
    x <- df[[nm]]
    type <- if (inherits(x, "Date")) "Date" else if (is.character(x)) "character" else
      if (is.logical(x)) "logical" else if (is.integer(x)) "integer" else
      if (is.double(x)) "double" else stop("unsupported column type: ", nm)
    list(name = nm, type = type,
         values = if (type == "character") as.list(x) else NULL)
  })
  con <- gzfile(file.path(out, paste0(name, ".bin.gz")), "wb")
  for (nm in names(df)) {
    x <- df[[nm]]
    if (inherits(x, "Date")) writeBin(as.integer(x), con, size = 4, endian = "little")
    else if (is.logical(x) || is.integer(x)) writeBin(as.integer(x), con, size = 4, endian = "little")
    else if (is.double(x)) writeBin(x, con, size = 8, endian = "little")
  }
  close(con)
  jsonlite::write_json(list(nrow = nrow(df), columns = cols),
                       file.path(out, paste0(name, ".json")),
                       auto_unbox = TRUE, null = "null", na = "null", digits = NA)
}

ext <- file.path("inst", "extdata")
sst <- file.path(ext, "sst_test.nc")
glo <- file.path(ext, "glorys_depth_test.nc")
glo_dir <- file.path(ext, "glorys_depth_test_annual")
wd <- tempfile("golden"); dir.create(wd)
q <- function(expr) invisible(capture.output(expr))
cp <- c("1982-01-01", "2011-12-31")

# Single-file heatwaves (single-threaded, so event order is deterministic)
s <- file.path(wd, "s")
q(ts2clm3(sst, name = s, climatologyPeriod = cp, quiet = TRUE))
q(detect_event3(sst, name = s, category = TRUE, daily = "also", quiet = TRUE))
save_golden(hw3_export(paste0(s, "_clim.nc")), "clim")
save_golden(hw3_export(paste0(s, "_clim.nc"), vars = "thresh",
                       lon_range = c(26.5, 26.55), n = 500), "clim_sub")
save_golden(hw3_export(paste0(s, "_events.nc")), "events")
save_golden(hw3_export(paste0(s, "_events.nc"), vars = c("duration", "intensity_max"),
                       lat_range = c(-34.2, -34.05),
                       time_range = c("2010-01-01", "2012-12-31")), "events_sub")
save_golden(hw3_export(paste0(s, "_events_daily.nc")), "daily")
save_golden(hw3_export(paste0(s, "_events_daily.nc"), vars = c("temp", "event", "category"),
                       time_range = c("2015-01-01", "2015-12-31")), "daily_sub")
save_golden(category3(paste0(s, "_events.nc")), "category3")
save_golden(block_average3(paste0(s, "_events.nc")), "block_average3")
save_golden(category_daily3(sst, paste0(s, "_clim.nc"), paste0(s, "_events.nc"),
                            time_range = c("2019-06-01", "2019-12-31")), "category_daily3")

# Proto-events and the category fallback (no inline categories, north)
p <- file.path(wd, "p")
invisible(file.copy(paste0(s, "_clim.nc"), paste0(p, "_clim.nc")))
q(detect_event3(sst, name = p, protoEvent = TRUE, quiet = TRUE))
save_golden(hw3_export(paste0(p, "_protoevents.nc"), time_range = c("2016-01-01", "2016-12-31")),
            "proto_sub")
q(detect_event3(sst, name = p, quiet = TRUE))
save_golden(category3(paste0(p, "_events.nc"), paste0(p, "_clim.nc"), hemisphere = "north"),
            "category3_fallback")

# Cold-spells via detect3()
k <- file.path(wd, "k")
q(detect3(sst, name = k, climatologyPeriod = cp, pctile = 10, coldSpells = TRUE,
          category = TRUE, quiet = TRUE))
save_golden(hw3_export(paste0(k, "_events.nc")), "cold_events")
save_golden(category_daily3(sst, paste0(k, "_clim.nc"), paste0(k, "_events.nc"),
                            time_range = c("2005-01-01", "2005-12-31"), coldSpells = TRUE),
            "category_daily3_cold")

# Static-threshold exceedance
ex <- file.path(wd, "exceed.nc")
q(exceedance3(sst, ex, threshold = 296))
save_golden(hw3_export(ex), "exceedance3")

# Depth-resolved (single file) and multi-file directory input
g <- file.path(wd, "g")
q(ts2clm3(glo, name = g, climatologyPeriod = c("1993-01-01", "2019-12-31"),
          depth_range = c(0, 10), quiet = TRUE))
q(detect_event3(glo, name = g, category = TRUE, daily = "also", quiet = TRUE))
save_golden(hw3_export(paste0(g, "_clim.nc")), "depth_clim")
save_golden(hw3_export(paste0(g, "_events.nc")), "depth_events")
save_golden(hw3_export(paste0(g, "_events_daily.nc"), depth_range = c(0, 3),
                       time_range = c("2010-01-01", "2010-03-01")), "depth_daily_sub")
save_golden(hw3_export(paste0(g, "_events.nc"), depth_range = c(0, 3)), "depth_events_sub")
m <- file.path(wd, "m")
q(detect3(glo_dir, name = m, climatologyPeriod = c("1993-01-01", "2019-12-31"),
          depth_range = c(0, 10), quiet = TRUE))
save_golden(hw3_export(paste0(m, "_events.nc")), "multi_events")

# Blobs: 3D (all components), wrapped + filtered + ranked, cold-spells,
# depth-resolved (4D) and multi-file input
save_blob <- function(res, name) for (k in names(res)) save_golden(res[[k]], paste0(name, "_", k))
save_blob(detect_blob3(sst, paste0(s, "_clim.nc"), return = c("event", "daily", "voxel")), "blob")
save_blob(detect_blob3(sst, paste0(s, "_clim.nc"), wrapDateline = TRUE, minVoxels = 20L,
                       topN = 5L, rankBy = "peakArea_km2"), "blob_filtered")
save_blob(detect_blob3(sst, paste0(k, "_clim.nc"), coldSpells = TRUE), "blob_cold")
save_blob(detect_blob3(glo, paste0(g, "_clim.nc"), minVoxels = 10L, rankBy = "totalVolume_km3",
                       return = c("event", "daily", "voxel")), "blob_depth")
q(save_blob(detect_blob3(glo_dir, paste0(g, "_clim.nc"), topN = 3L), "blob_multi"))

# Plotting: the data each R plot draws. Flame polygons come straight from
# GeomFlame3$draw_group() with an identity coordinate transform.
flame_poly <- function(d, ...) {
  d <- cbind(d, colour = NA, fill = "salmon", linewidth = 0.5, linetype = 1, alpha = NA, group = 1)
  g <- GeomFlame3$draw_group(d, NULL, list(transform = function(data, panel_params) data), ...)
  if (inherits(g, "null")) return(data.frame(x = numeric(0), y = numeric(0), id = integer(0)))
  data.frame(x = as.numeric(g$x), y = as.numeric(g$y), id = as.integer(g$id))
}
set.seed(1)
fx <- seq(0, 30, length.out = 400)
flame_in <- data.frame(x = fx, y = sin(fx) + rnorm(400, sd = 0.3), y2 = 0.4 + 0.1 * cos(fx / 3))
flame_in$y[c(57, 200:203)] <- NA
save_golden(flame_in, "flame_input")
save_golden(flame_poly(flame_in), "flame_default")
save_golden(flame_poly(flame_in, n = 3), "flame_n3")
save_golden(flame_poly(flame_in, n = 1, n_gap = 2), "flame_gap")
save_golden(flame_poly(flame_in, reverse = TRUE), "flame_reverse")

save_line <- function(p, name) {
  save_golden(p$data, name)
  save_golden(data.frame(title = p$labels$title), paste0(name, "_title"))
}
q(save_line(event_line3(sst, paste0(s, "_clim.nc"), lon = 26.525, lat = -34.125,
                        event_file = paste0(s, "_events.nc")), "line_event"))
q(save_line(event_line3(sst, paste0(s, "_clim.nc"), lon = 26.575, lat = -34.025,
                        start_date = "2015-01-01", end_date = "2016-06-30"), "line_window"))
q(save_line(event_line3(sst, paste0(k, "_clim.nc"), lon = 26.525, lat = -34.075,
                        event_file = paste0(k, "_events.nc"), coldSpells = TRUE,
                        metric = "intensity_max", spread = 60), "line_cold"))
q(save_line(event_line3(glo, paste0(g, "_clim.nc"), lon = 25.16667, lat = -34.91667, depth = 5,
                        start_date = "2010-01-01", end_date = "2010-12-31"), "line_depth"))
save_golden(plot_metric3(paste0(s, "_events.nc"), coastline = FALSE)$data, "metric_default")
save_golden(plot_metric3(paste0(g, "_events.nc"), metric = "duration", summary = "max",
                         coastline = FALSE)$data, "metric_duration_max")

cat("Golden files written to", out, "\n")

# Signatures of the exported R functions, for the API-parity test
sigs <- lapply(setNames(nm = getNamespaceExports("heatwave3")), function(f) {
  fn <- get(f, envir = asNamespace("heatwave3"))
  if (!is.function(fn)) return(NULL)
  fm <- formals(fn)
  lapply(fm, function(x) if (identical(x, quote(expr = ))) "<required>" else paste(deparse(x), collapse = " "))
})
jsonlite::write_json(Filter(Negate(is.null), sigs), file.path(out, "..", "r_signatures.json"),
                     auto_unbox = TRUE, pretty = TRUE)
cat("R signatures written\n")
