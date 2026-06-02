#' Convenience wrapper for marine heatwave event categories
#'
#' \code{category3()} is a convenience wrapper for adding or reading Hobday
#' et al. (2018) categories (I Moderate, II Strong, III Severe, IV Extreme)
#' for an event NetCDF file. New analyses should usually set
#' \code{category = TRUE} in \code{\link{detect_event3}} or
#' \code{\link{detect3}}, which computes the same categories during event
#' detection and avoids a second pass over the files.
#'
#' If the event file already contains categories, \code{category3()} reads
#' those values directly and \code{clim_file} may be omitted. If the event file
#' does not contain categories, \code{clim_file} is required so the category
#' width can be computed from the seasonal and threshold climatologies.
#'
#' Works for both marine heatwaves and cold-spells: the category width is
#' \code{|seas - thresh|}, which is always positive regardless of whether the
#' threshold is above (heatwave) or below (cold-spell) the seasonal mean.
#'
#' @param event_file Path to the event NetCDF file from \code{\link{detect_event3}}.
#' @param clim_file Path to the climatology NetCDF from \code{\link{ts2clm3}}.
#'   Only required when the event file does not contain pre-computed categories.
#' @param hemisphere Character. Either \code{"south"} (default, austral:
#'   DJF = Summer) or \code{"north"} (boreal: DJF = Winter).
#' @param roundVal Decimal places for rounding. Default \code{4}.
#' @param S Deprecated. Use \code{hemisphere} instead.
#'
#' @return A data.frame with columns: event_no, lon, lat, peak_date, category,
#'   intensity_max, duration, p_moderate, p_strong, p_severe, p_extreme, season.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
#' clim_file <- tempfile(fileext = ".nc")
#' event_file <- tempfile(fileext = ".nc")
#'
#' detect3(sst_file, clim_file, event_file,
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"))
#'
#' cats <- category3(event_file, clim_file)
#' table(cats$category)
#' }
category3 <- function(event_file, clim_file = NULL,
                      hemisphere = "south",
                      roundVal = 4L,
                      S = NULL) {

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
        "clim_file is required because the event file does not contain ",
        "pre-computed categories.",
        call. = FALSE
      )
    }
    clim_file <- ""
  } else if (!file.exists(clim_file)) {
    stop("Climatology file does not exist: ", clim_file, call. = FALSE)
  }

  hw3_category(event_file, clim_file, hemisphere = hemisphere,
               roundVal = roundVal)
}
