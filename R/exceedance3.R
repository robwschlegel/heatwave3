#' Static threshold exceedance detection for gridded data
#'
#' Detects periods where SST exceeds (or falls below) a fixed threshold
#' for each pixel, without requiring a climatology baseline.
#'
#' @param file_in Path to input NetCDF file.
#' @param file_out Path for output event NetCDF file.
#' @param threshold The static threshold value.
#' @param var_name SST variable name. If \code{NULL}, auto-detected.
#' @param below Logical. If \code{TRUE}, detect exceedances below the threshold.
#'   Default \code{FALSE}.
#' @param minDuration Minimum event duration in days. Default \code{5}.
#' @param joinAcrossGaps Join events across short gaps? Default \code{TRUE}.
#' @param maxGap Maximum gap to join. Default \code{2}.
#' @param lon_range Optional \code{c(min, max)} longitude range.
#' @param lat_range Optional \code{c(min, max)} latitude range.
#' @param time_range Optional \code{c("start", "end")} date range.
#' @param depth Optional depth index.
#' @param roundRes Decimal places for rounding. Default \code{4}.
#' @param n_threads Number of threads for parallel computation. Default \code{1}.
#'
#' @return Invisibly returns the output file path.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
#' event_file <- tempfile(fileext = ".nc")
#'
#' # Detect periods where SST exceeds 295 K
#' exceedance3(sst_file, event_file, threshold = 295)
#' }
exceedance3 <- function(file_in, file_out, threshold,
                        var_name = NULL,
                        below = FALSE,
                        minDuration = 5L,
                        joinAcrossGaps = TRUE,
                        maxGap = 2L,
                        lon_range = NULL, lat_range = NULL,
                        time_range = NULL, depth = NULL,
                        roundRes = 4L,
                        n_threads = 1L) {

  if (missing(threshold))
    stop("A threshold value must be provided.", call. = FALSE)

  # Create a dummy climatology where seas = threshold and thresh = threshold
  # This reuses the event detection machinery
  clim_file <- tempfile(fileext = ".nc")
  on.exit(unlink(clim_file), add = TRUE)

  # Create a constant climatology file using ts2clm3 with a dummy period,
  # then overwrite seas/thresh with the static threshold value.
  # For simplicity, run ts2clm3 first then modify the output.
  vn <- if (is.null(var_name)) "" else var_name

  # Write a constant-value climatology NetCDF via the C++ helper
  hw3_write_const_clim(
    file_in = file_in,
    clim_file = clim_file,
    var_name = vn,
    threshold = as.double(threshold),
    lon_range = lon_range,
    lat_range = lat_range,
    time_range = time_range,
    depth = if (is.null(depth)) -1L else as.integer(depth)
  )

  # Call the C++ detector directly so exceedance3 keeps its exact-path file_out
  # API (the detect_event3 wrapper derives paths from a name stem instead).
  hw3_detect_events(
    file_in = file_in, clim_file = clim_file, events_file = file_out,
    var_name = vn,
    minDuration = as.integer(minDuration),
    minDuration2 = as.integer(minDuration),
    joinAcrossGaps = joinAcrossGaps,
    maxGap = as.integer(maxGap),
    maxGap2 = as.integer(maxGap),
    coldSpells = below,
    roundRes = as.integer(roundRes),
    n_threads = as.integer(n_threads),
    category = FALSE,
    southHemisphere = TRUE,
    threshClim2_file = "",
    threshClim2_var_name = "",
    daily_file = "",
    proto_file = ""
  )

  invisible(file_out)
}
