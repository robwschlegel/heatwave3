#' Detect marine heatwave events in gridded data
#'
#' Detects marine heatwave (or cold-spell) events for each pixel in a gridded
#' dataset, using climatologies computed by \code{\link{ts2clm3}}. This is the
#' gridded equivalent of \code{heatwaveR::detect_event()}.
#'
#' @param file_in Path to the input NetCDF file containing SST data.
#' @param clim_file Path to the climatology NetCDF file produced by
#'   \code{\link{ts2clm3}}.
#' @param file_out Path for the output NetCDF file containing detected events.
#' @param var_name Name of the SST variable. If \code{NULL}, auto-detected.
#' @param minDuration Minimum duration (days) for an event. Default \code{5}.
#' @param minDuration2 Minimum duration for secondary threshold events.
#'   Default equals \code{minDuration}.
#' @param joinAcrossGaps Logical. Join events separated by short gaps?
#'   Default \code{TRUE}.
#' @param maxGap Maximum gap length (days) to join across. Default \code{2}.
#' @param maxGap2 Maximum gap for secondary threshold. Default equals \code{maxGap}.
#' @param coldSpells Logical. Detect cold-spells instead of heatwaves?
#'   Default \code{FALSE}.
#' @param roundRes Number of decimal places for rounding event metrics.
#'   Default \code{4}.
#' @param save_format Additional output format. Default \code{NULL} (NetCDF only).
#'   Options: \code{"csv"}, \code{"rda"}, \code{"parquet"}.
#' @param n_threads Number of OpenMP threads. Default \code{1}.
#'
#' @return Invisibly returns the path to the output file. Events are written as
#'   a CF ragged array NetCDF with 19 event metrics per event.
#'
#' @export
#'
#' @examples
#' \donttest{
#' sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
#' clim_file <- tempfile(fileext = ".nc")
#' event_file <- tempfile(fileext = ".nc")
#'
#' ts2clm3(sst_file, clim_file,
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"))
#'
#' detect_event3(sst_file, clim_file, event_file)
#'
#' # Read the result
#' nc <- ncdf4::nc_open(event_file)
#' cat("Events detected:", nc$dim$event$len, "\n")
#' ncdf4::nc_close(nc)
#' }
detect_event3 <- function(file_in, clim_file, file_out,
                          var_name = NULL,
                          minDuration = 5L,
                          minDuration2 = minDuration,
                          joinAcrossGaps = TRUE,
                          maxGap = 2L,
                          maxGap2 = maxGap,
                          coldSpells = FALSE,
                          roundRes = 4L,
                          save_format = NULL,
                          n_threads = 1L) {

  if (missing(file_in) || missing(clim_file) || missing(file_out))
    stop("file_in, clim_file, and file_out must all be provided.", call. = FALSE)
  if (!file.exists(clim_file))
    stop("Climatology file does not exist: ", clim_file, call. = FALSE)

  vn <- if (is.null(var_name)) "" else var_name

  # Determine input mode: single file, vector of files, or directory
  multi_file <- FALSE
  if (length(file_in) == 1 && dir.exists(file_in)) {
    file_in <- sort(list.files(file_in, pattern = "\\.(nc|nc4)$",
                               full.names = TRUE))
    if (length(file_in) == 0)
      stop("No .nc or .nc4 files found in the specified directory.", call. = FALSE)
    multi_file <- TRUE
  } else if (length(file_in) > 1) {
    multi_file <- TRUE
  }

  if (multi_file) {
    hw3_detect_events_multi(
      files = file_in,
      clim_file = clim_file,
      file_out = file_out,
      var_name = vn,
      minDuration = as.integer(minDuration),
      minDuration2 = as.integer(minDuration2),
      joinAcrossGaps = joinAcrossGaps,
      maxGap = as.integer(maxGap),
      maxGap2 = as.integer(maxGap2),
      coldSpells = coldSpells,
      roundRes = as.integer(roundRes),
      n_threads = as.integer(n_threads)
    )
  } else {
    if (!file.exists(file_in))
      stop("SST file does not exist: ", file_in, call. = FALSE)
    hw3_detect_events(
      file_in = file_in,
      clim_file = clim_file,
      file_out = file_out,
      var_name = vn,
      minDuration = as.integer(minDuration),
      minDuration2 = as.integer(minDuration2),
      joinAcrossGaps = joinAcrossGaps,
      maxGap = as.integer(maxGap),
      maxGap2 = as.integer(maxGap2),
      coldSpells = coldSpells,
      roundRes = as.integer(roundRes),
      n_threads = as.integer(n_threads)
    )
  }

  if (!is.null(save_format)) {
    hw3_export(file_out, format = save_format, type = "event")
  }

  invisible(file_out)
}
