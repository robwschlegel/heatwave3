#' All-in-one marine heatwave detection for gridded data
#'
#' The primary entry point for \code{heatwave3}. Runs
#' \code{\link{ts2clm3}} followed by \code{\link{detect_event3}} in a
#' single call, with optional inline category computation.
#'
#' @param file_in Path to input NetCDF file (or directory/vector of daily files).
#' @param file_out_clim Path for the output climatology NetCDF file.
#' @param file_out_event Path for the output event NetCDF file.
#' @param climatologyPeriod Character vector of length 2 with start and end
#'   dates of the baseline period, e.g. \code{c("1991-01-01", "2020-12-31")}.
#' @param var_name Name of the SST variable. If \code{NULL}, auto-detected.
#' @param lon_range Optional \code{c(min, max)} longitude range.
#' @param lat_range Optional \code{c(min, max)} latitude range.
#' @param time_range Optional \code{c("start", "end")} date range.
#' @param depth Optional depth/level index for 4D data.
#' @param maxPadLength Max consecutive NAs to interpolate. Default \code{FALSE}.
#' @param windowHalfWidth Half-width of climatology window. Default \code{5}.
#' @param pctile Percentile for threshold. Default \code{90} (heatwaves);
#'   use \code{10} for cold-spells.
#' @param smoothPercentile Apply rolling mean smoothing? Default \code{TRUE}.
#' @param smoothPercentileWidth Smoothing window width. Default \code{31}.
#' @param minDuration Minimum event duration in days. Default \code{5}.
#' @param joinAcrossGaps Join events across short gaps? Default \code{TRUE}.
#' @param maxGap Maximum gap length to join. Default \code{2}.
#' @param coldSpells Detect cold-spells? Default \code{FALSE}.
#' @param category Logical. Compute Hobday et al. (2018) severity categories
#'   inline? Default \code{FALSE}.
#' @param hemisphere Character. \code{"south"} (default) or \code{"north"}.
#' @param roundClm Decimal places for climatology rounding. Default \code{4}.
#' @param roundRes Decimal places for event metric rounding. Default \code{4}.
#' @param return_df Logical. Return the event table as a \code{data.frame}?
#'   Default \code{FALSE}.
#' @param n_threads Number of OpenMP threads. Default \code{1}.
#'
#' @return If \code{return_df = FALSE} (the default), invisibly returns a list
#'   with \code{clim_file} and \code{event_file} paths. If
#'   \code{return_df = TRUE}, returns a \code{data.frame} of event metrics.
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
#' }
detect3 <- function(file_in, file_out_clim, file_out_event,
                    climatologyPeriod,
                    var_name = NULL,
                    lon_range = NULL, lat_range = NULL,
                    time_range = NULL, depth = NULL,
                    maxPadLength = FALSE,
                    windowHalfWidth = 5L,
                    pctile = 90,
                    smoothPercentile = TRUE,
                    smoothPercentileWidth = 31L,
                    minDuration = 5L,
                    joinAcrossGaps = TRUE,
                    maxGap = 2L,
                    coldSpells = FALSE,
                    category = FALSE,
                    hemisphere = "south",
                    roundClm = 4L,
                    roundRes = 4L,
                    return_df = FALSE,
                    n_threads = 1L) {

  ts2clm3(file_in = file_in, file_out = file_out_clim,
          climatologyPeriod = climatologyPeriod,
          lon_range = lon_range, lat_range = lat_range,
          time_range = time_range, depth = depth,
          var_name = var_name,
          maxPadLength = maxPadLength,
          windowHalfWidth = windowHalfWidth,
          pctile = pctile,
          smoothPercentile = smoothPercentile,
          smoothPercentileWidth = smoothPercentileWidth,
          roundClm = roundClm,
          n_threads = n_threads)

  result <- detect_event3(
    file_in = file_in, clim_file = file_out_clim,
    file_out = file_out_event,
    var_name = var_name,
    minDuration = minDuration,
    joinAcrossGaps = joinAcrossGaps,
    maxGap = maxGap,
    coldSpells = coldSpells,
    category = category,
    hemisphere = hemisphere,
    roundRes = roundRes,
    return_df = return_df,
    n_threads = n_threads
  )

  if (return_df) return(result)

  invisible(list(clim_file = file_out_clim, event_file = file_out_event))
}
