#' Categorise marine heatwave events
#'
#' Assigns Hobday et al. (2018) categories (I Moderate, II Strong, III Severe,
#' IV Extreme) to events. If the event file was produced with
#' \code{category = TRUE} in \code{\link{detect_event3}}, the pre-computed
#' categories are read directly. Otherwise they are computed from the
#' climatology file.
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
    hemisphere <- if (isTRUE(S)) "south" else "north"
  }
  hemisphere <- match.arg(hemisphere, c("south", "north"))

  if (is.null(clim_file)) {
    clim_file <- ""
  }

  hw3_category(event_file, clim_file, hemisphere = hemisphere,
               roundVal = roundVal)
}
