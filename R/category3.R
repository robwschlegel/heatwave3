#' Convenience wrapper for marine heatwave event categories
#'
#' \code{category3()} is a convenience wrapper for adding or reading Hobday
#' et al. (2018) categories (I Moderate, II Strong, III Severe, IV Extreme)
#' for an event NetCDF file. New analyses should usually set
#' \code{category = TRUE} in \code{\link{detect_event3}} or
#' \code{\link{detect3}}, which computes the same categories during event
#' detection and avoids a second pass over the files.
#'
#' Categories follow \code{heatwaveR::category()}. For each event day the
#' category limits are the threshold and \code{thresh + k * diff} (k = 1, 2, 3),
#' with \code{diff = thresh - seas} and each limit rounded to \code{roundVal}
#' decimals. \code{p_moderate} to \code{p_extreme} are the percentages of the
#' event's duration (days, rounded to whole numbers) spent in each category,
#' and the event's category is the highest one with a non-zero percentage.
#' Works for both marine heatwaves and cold-spells (cold-spells are recognised,
#' as in heatwaveR, by their negative intensities). The season is that of the
#' event's peak month.
#'
#' If the event file already contains categories, \code{category3()} reads
#' those values directly, and \code{clim_file} and \code{sst_file} may be
#' omitted. Otherwise both are required, because the percentages are counts of
#' days, which need each event day's temperature.
#'
#' @param event_file Path to the event NetCDF file from \code{\link{detect_event3}}.
#' @param clim_file Path to the climatology NetCDF from \code{\link{ts2clm3}}.
#'   Only required when the event file does not contain pre-computed categories.
#' @param hemisphere Character. Either \code{"south"} (default, austral:
#'   DJF = Summer) or \code{"north"} (boreal: DJF = Winter).
#' @param roundVal Decimal places for rounding. Default \code{4}.
#' @param S Deprecated. Use \code{hemisphere} instead.
#' @param sst_file The SST used for detection (a file, a directory of daily
#'   files or a vector of files). Only required when the event file does not
#'   contain pre-computed categories.
#' @param var_name Name of the SST variable in \code{sst_file}; auto-detected
#'   when \code{NULL}.
#'
#' @return A data.frame with columns: event_no, lon, lat, peak_date, category,
#'   intensity_max, duration, p_moderate, p_strong, p_severe, p_extreme (whole
#'   percentages), season.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
#' stem <- file.path(tempdir(), "demo")
#'
#' detect3(sst_file, name = stem,
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"))
#'
#' cats <- category3(paste0(stem, "_events.nc"), paste0(stem, "_clim.nc"))
#' table(cats$category)
#' }
category3 <- function(event_file, clim_file = NULL,
                      hemisphere = "south",
                      roundVal = 4L,
                      S = NULL,
                      sst_file = NULL,
                      var_name = NULL) {

  if (!is.null(S)) {
    warning("S is deprecated; use hemisphere instead.", call. = FALSE)
    hemisphere <- if (isTRUE(S)) "south" else "north"
  }
  hemisphere <- match.arg(hemisphere, c("south", "north"))

  if (missing(event_file) || !file.exists(event_file)) {
    stop("Event file does not exist: ", event_file, call. = FALSE)
  }

  if (is.null(clim_file)) {
    ev <- hw3_read_event_nc(event_file)
    if (length(ev$category) == 0 || !any(ev$category > 0)) {
      stop(
        "clim_file and sst_file are required because the event file does not ",
        "contain pre-computed categories.",
        call. = FALSE
      )
    }
    clim_file <- ""
  } else if (!file.exists(clim_file)) {
    stop("Climatology file does not exist: ", clim_file, call. = FALSE)
  }

  sst_files <- if (is.null(sst_file)) NULL else .resolve_nc_input(sst_file, "sst_file")
  hw3_category(event_file, clim_file, hemisphere = hemisphere,
               roundVal = as.integer(roundVal), sst_files = sst_files,
               var_name = if (is.null(var_name)) "" else var_name)
}
