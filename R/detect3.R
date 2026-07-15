#' All-in-one marine heatwave detection for gridded data
#'
#' The primary entry point for \code{heatwave3}. Runs
#' \code{\link{ts2clm3}} followed by \code{\link{detect_event3}} in a
#' single call, with optional inline category computation.
#'
#' @param file_in Path to input NetCDF file (or directory/vector of daily files).
#' @param name Output base name (a path stem) shared by both stages. Produces
#'   \code{paste0(name, "_clim.nc")} and \code{paste0(name, "_events.nc")} (plus
#'   the per-day products; see \code{daily} and \code{protoEvent}). The
#'   directory must exist.
#' @param climatologyPeriod Character vector of length 2 with start and end
#'   dates of the baseline period, for example \code{c("1991-01-01", "2020-12-31")}.
#' @param var_name Name of the SST variable. If \code{NULL}, auto-detected.
#' @param lon_range Optional \code{c(min, max)} longitude range.
#' @param lat_range Optional \code{c(min, max)} latitude range.
#' @param time_range Optional \code{c("start", "end")} date range.
#' @param depth Optional depth/level index for 4D data.
#' @param depth_range Optional numeric vector of length 2:
#'   \code{c(min_depth, max_depth)}. Runs the climatology and event detection
#'   natively across this contiguous band of depth levels instead of a
#'   single squeezed level -- see \code{\link{ts2clm3}}'s \code{depth_range}
#'   for details. Mutually exclusive with \code{depth}. Compatible with
#'   \code{daily}/\code{protoEvent}: the per-day product gains a \code{depth}
#'   dimension the same way the climatology and event files do.
#' @param maxPadLength Max consecutive NAs to interpolate. Default \code{FALSE}.
#' @param windowHalfWidth Half-width of climatology window. Default \code{5}.
#' @param pctile Percentile for threshold. Default \code{90} (heatwaves);
#'   use \code{10} for cold-spells.
#' @param smoothPercentile Apply rolling mean smoothing? Default \code{TRUE}.
#' @param smoothPercentileWidth Smoothing window width. Default \code{31}.
#' @param detrend Logical. Remove a linear trend before climatology calculation?
#'   Default \code{FALSE}.
#' @param minDuration Minimum event duration in days. Default \code{5}.
#' @param minDuration2 Secondary minimum duration. See
#'   \code{\link{detect_event3}}.
#' @param joinAcrossGaps Join events across short gaps? Default \code{TRUE}.
#' @param maxGap Maximum gap length to join. Default \code{2}.
#' @param maxGap2 Secondary maximum gap length. See
#'   \code{\link{detect_event3}}.
#' @param threshClim2 Optional gridded NetCDF logical criterion for secondary
#'   event detection. See \code{\link{detect_event3}}.
#' @param threshClim2_var_name Name of the \code{threshClim2} variable. If
#'   \code{NULL}, it is auto-detected.
#' @param coldSpells Detect cold-spells? Default \code{FALSE}.
#' @param category Logical. Compute Hobday et al. (2018) severity categories
#'   inline? Default \code{FALSE}.
#' @param hemisphere Character. \code{"south"} (default) or \code{"north"}.
#' @param roundClm Decimal places for climatology rounding. Default \code{4}.
#' @param roundRes Decimal places for event metric rounding. Default \code{4}.
#' @param daily Per-day output control: \code{"none"} (default), \code{"also"},
#'   or \code{"only"}. See \code{\link{detect_event3}}.
#' @param protoEvent Logical. Write the per-day proto-event series
#'   (\code{paste0(name, "_protoevents.nc")}) instead of the event table. See
#'   \code{\link{detect_event3}}. Default \code{FALSE}.
#' @param n_threads Number of OpenMP threads. Default \code{1}.
#' @param skip_bad_files Logical. For multi-file inputs, skip unreadable files
#'   or files with mismatched grids instead of failing. Default \code{FALSE}.
#' @param quiet Logical. Suppress the post-computation console summaries?
#'   Default \code{FALSE}.
#'
#' @return Invisibly returns a named character vector of the NetCDF files
#'   written (\code{clim} plus a subset of \code{events}, \code{daily},
#'   \code{protoevents}). Use \code{\link{hw3_export}} to read or export any of
#'   them.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
#'
#' detect3(sst_file, name = file.path(tempdir(), "demo"),
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"))
#' }
detect3 <- function(file_in, name,
                    climatologyPeriod,
                    var_name = NULL,
                    lon_range = NULL, lat_range = NULL,
                    time_range = NULL, depth = NULL,
                    depth_range = NULL,
                    maxPadLength = FALSE,
                    windowHalfWidth = 5L,
                    pctile = 90,
                    smoothPercentile = TRUE,
                    smoothPercentileWidth = 31L,
                    detrend = FALSE,
                    minDuration = 5L,
                    minDuration2 = minDuration,
                    joinAcrossGaps = TRUE,
                    maxGap = 2L,
                    maxGap2 = maxGap,
                    threshClim2 = NULL,
                    threshClim2_var_name = NULL,
                    coldSpells = FALSE,
                    category = FALSE,
                    hemisphere = "south",
                    roundClm = 4L,
                    roundRes = 4L,
                    daily = c("none", "also", "only"),
                    protoEvent = FALSE,
                    n_threads = 1L,
                    skip_bad_files = FALSE,
                    quiet = FALSE) {

  daily <- match.arg(daily)

  clim_file <- ts2clm3(file_in = file_in, name = name,
          climatologyPeriod = climatologyPeriod,
          lon_range = lon_range, lat_range = lat_range,
          time_range = time_range, depth = depth,
          depth_range = depth_range,
          var_name = var_name,
          maxPadLength = maxPadLength,
          windowHalfWidth = windowHalfWidth,
          pctile = pctile,
          smoothPercentile = smoothPercentile,
          smoothPercentileWidth = smoothPercentileWidth,
          detrend = detrend,
          roundClm = roundClm,
          n_threads = n_threads,
          skip_bad_files = skip_bad_files,
          quiet = quiet)

  written <- detect_event3(
    file_in = file_in, name = name, clim_file = clim_file,
    var_name = var_name,
    minDuration = minDuration,
    minDuration2 = minDuration2,
    joinAcrossGaps = joinAcrossGaps,
    maxGap = maxGap,
    maxGap2 = maxGap2,
    threshClim2 = threshClim2,
    threshClim2_var_name = threshClim2_var_name,
    coldSpells = coldSpells,
    category = category,
    hemisphere = hemisphere,
    roundRes = roundRes,
    daily = daily,
    protoEvent = protoEvent,
    n_threads = n_threads,
    skip_bad_files = skip_bad_files,
    quiet = quiet
  )

  invisible(c(clim = clim_file, written))
}
