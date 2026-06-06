#' Daily marine heatwave / cold-spell categories for a date window
#'
#' Builds the per-pixel, per-day series for a window of dates directly from an
#' SST file, a heatwave3 climatology, and a heatwave3 events file, using the C++
#' readers throughout. It returns the full daily grid for the window with the
#' same columns as the \code{\link{detect_event3}} \code{daily = "also"} product
#' (\code{lon}, \code{lat}, \code{t}, \code{temp}, \code{seas}, \code{thresh},
#' \code{intensity}, \code{event}, \code{event_no}, \code{category}).
#'
#' Crucially, event membership is read from the events file rather than
#' re-detected on the window, so events that began before \code{time_range} are
#' still recognised. Re-running \code{detect_event3()} on a short recent window
#' instead would silently drop ongoing events whose in-window portion is shorter
#' than the minimum duration. Only the requested time window of the SST file is
#' read (a hyperslab), which makes it cheap to call repeatedly for the most
#' recent days, as a daily operational pipeline does.
#'
#' The Marine Heatwave Tracker's \code{load_sub_cat_clim()} subset (event-member
#' days with a category) is recovered with
#' \code{subset(x, !is.na(category))[c("t", "lon", "lat", "event_no",
#' "intensity", "category")]}.
#'
#' @param sst_file SST source: a single multi-time NetCDF (for example a
#'   per-longitude OISST time series), \strong{or} a directory / character vector
#'   of daily files (one time step each, e.g. daily-global OSTIA/GHRSST). Only
#'   files covering \code{time_range} should be passed for a directory/vector.
#'   The variable is auto-detected unless \code{var_name} is given.
#' @param clim_file Path to the climatology NetCDF from \code{\link{ts2clm3}}.
#' @param event_file Path to the events NetCDF from \code{\link{detect_event3}}.
#' @param time_range Character vector of length 2, \code{c(start, end)} dates,
#'   for example \code{c("2026-05-10", "2026-05-20")}.
#' @param lon_range,lat_range Optional numeric vectors of length 2,
#'   \code{c(min, max)}, restricting the output to a spatial sub-region. Like
#'   \code{time_range}, these subset the SST \strong{read} (a hyperslab),
#'   intersected with the climatology extent, so only the requested window is
#'   pulled from disk and the climatology and events are aligned to it. Default
#'   \code{NULL} (the full climatology extent).
#' @param var_name Name of the SST variable. If \code{NULL}, auto-detected.
#' @param coldSpells Logical. Categorise marine cold-spells instead of heatwaves?
#'   Default \code{FALSE}.
#' @param ice_thresh Numeric. For cold-spells only, days whose threshold is below
#'   this value (in the SST's units) are assigned category 5 ("ice"), following
#'   the Marine Heatwave Tracker convention. Default \code{-1.7}. Ignored for
#'   heatwaves.
#' @param roundRes Number of decimal places for the returned \code{intensity}.
#'   Default \code{2}, mirroring the Marine Heatwave Tracker's
#'   \code{load_sub_cat_clim()} (\code{round(temp - seas, 2)}). Pass
#'   \code{roundRes = 4} to match \code{\link{detect_event3}}'s
#'   \code{daily = "also"} product exactly.
#' @param skip_bad_files Logical. When \code{sst_file} is a directory or vector,
#'   skip unreadable files or files with mismatched grids instead of failing.
#'   Default \code{FALSE}.
#'
#' @return A \code{data.frame} with one row per pixel per day in the window:
#'   \code{lon}, \code{lat}, \code{t} (Date), \code{temp}, \code{seas},
#'   \code{thresh}, \code{intensity} (\code{temp - seas}), \code{event} (logical;
#'   is the day inside a detected event?), \code{event_no} (integer, \code{NA}
#'   off-event), and \code{category} (integer on event-member exceedance days;
#'   1 = I Moderate ... 4 = IV Extreme, 5 = ice for cold-spells; \code{NA}
#'   otherwise). Columns match the \code{detect_event3} \code{daily = "also"}
#'   product.
#'
#' @details
#' The SST, climatology, and events files must share the same spatial grid (as
#' they do when produced from one heatwave3 run for a given region or longitude).
#' Event membership and \code{event_no} are taken from \code{event_file}, so the
#' events must already cover the requested window.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' category_daily3(
#'   sst_file   = "oisst-avhrr-v02r01.ts.0001.nc",
#'   clim_file  = "MHW_0001_1982-2011_clim.nc",
#'   event_file = "MHW_0001_1982-2011_events.nc",
#'   time_range = c("2026-05-10", "2026-05-20")
#' )
#' }
category_daily3 <- function(sst_file, clim_file, event_file, time_range,
                            lon_range = NULL, lat_range = NULL,
                            var_name = NULL, coldSpells = FALSE,
                            ice_thresh = -1.7, roundRes = 2L,
                            skip_bad_files = FALSE) {
  if (!file.exists(clim_file)) {
    stop("Climatology file does not exist: ", clim_file, call. = FALSE)
  }
  if (!file.exists(event_file)) {
    stop("Events file does not exist: ", event_file, call. = FALSE)
  }
  if (length(time_range) != 2) {
    stop("time_range must be a vector of two dates.", call. = FALSE)
  }
  d <- as.Date(time_range)
  if (anyNA(d)) {
    stop("time_range must be dates, e.g. c(\"2026-05-10\", \"2026-05-20\").",
         call. = FALSE)
  }
  t_jd <- as.integer(d) + 2440588L

  if (!is.null(lon_range) && length(lon_range) != 2) {
    stop("lon_range must be NULL or a vector of two longitudes.", call. = FALSE)
  }
  if (!is.null(lat_range) && length(lat_range) != 2) {
    stop("lat_range must be NULL or a vector of two latitudes.", call. = FALSE)
  }
  lon_r <- if (is.null(lon_range)) numeric(0) else as.double(lon_range)
  lat_r <- if (is.null(lat_range)) numeric(0) else as.double(lat_range)
  vn <- if (is.null(var_name)) "" else var_name

  # SST source: single multi-time file, directory, or vector of daily files
  multi <- FALSE
  if (length(sst_file) == 1 && dir.exists(sst_file)) {
    sst_file <- sort(list.files(sst_file, pattern = "\\.(nc|nc4)$",
                                full.names = TRUE))
    if (length(sst_file) == 0) {
      stop("No .nc or .nc4 files found in the SST directory.", call. = FALSE)
    }
    multi <- TRUE
  } else if (length(sst_file) > 1) {
    multi <- TRUE
  } else if (!file.exists(sst_file)) {
    stop("File does not exist: ", sst_file, call. = FALSE)
  }

  res <- if (multi) {
    hw3_category_daily_multi(
      sst_file, clim_file, event_file, t_jd_range = t_jd,
      lon_range = lon_r, lat_range = lat_r, var_name = vn,
      coldSpells = coldSpells, ice_thresh = as.double(ice_thresh),
      roundRes = as.integer(roundRes), skip_bad_files = skip_bad_files
    )
  } else {
    hw3_category_daily(
      sst_file, clim_file, event_file, t_jd_range = t_jd,
      lon_range = lon_r, lat_range = lat_r, var_name = vn,
      coldSpells = coldSpells, ice_thresh = as.double(ice_thresh),
      roundRes = as.integer(roundRes)
    )
  }

  data.frame(
    lon = res$lon,
    lat = res$lat,
    t = as.Date(res$jd - 2440588L, origin = "1970-01-01"),
    temp = res$temp,
    seas = res$seas,
    thresh = res$thresh,
    intensity = res$intensity,
    event = res$event,
    event_no = res$event_no,
    category = res$category,
    stringsAsFactors = FALSE
  )
}
