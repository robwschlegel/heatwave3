#' Detect marine heatwave events in gridded data
#'
#' Detects marine heatwave (or cold-spell) events for each pixel in a gridded
#' dataset, using climatologies computed by \code{\link{ts2clm3}}. This is the
#' gridded equivalent of \code{heatwaveR::detect_event()}.
#'
#' @param file_in Path to the input NetCDF file containing SST data.
#'   May also be a character vector of file paths or a directory path
#'   containing daily \code{.nc}/\code{.nc4} files.
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
#' @param category Logical. Compute Hobday et al. (2018) severity categories
#'   (I Moderate through IV Extreme) inline during detection? Categories are
#'   written to the event NetCDF as additional variables. Default \code{FALSE}.
#' @param hemisphere Character. Season-naming convention: \code{"south"}
#'   (default, austral: DJF = Summer) or \code{"north"} (boreal:
#'   DJF = Winter). Only used when \code{category = TRUE}.
#' @param roundRes Number of decimal places for rounding event metrics.
#'   Default \code{4}.
#' @param return_df Logical. If \code{TRUE}, return the event table as a
#'   \code{data.frame} in addition to writing the NetCDF file. Default
#'   \code{FALSE} (returns the file path invisibly).
#' @param save_format Additional output format. Default \code{NULL} (NetCDF only).
#'   Options: \code{"csv"}, \code{"rda"}, \code{"parquet"}.
#' @param n_threads Number of OpenMP threads. Default \code{1}.
#'
#' @return If \code{return_df = FALSE} (the default), invisibly returns the
#'   path to the output file. If \code{return_df = TRUE}, returns a
#'   \code{data.frame} of event metrics (the NetCDF is still written).
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
#' }
detect_event3 <- function(file_in, clim_file, file_out,
                          var_name = NULL,
                          minDuration = 5L,
                          minDuration2 = minDuration,
                          joinAcrossGaps = TRUE,
                          maxGap = 2L,
                          maxGap2 = maxGap,
                          coldSpells = FALSE,
                          category = FALSE,
                          hemisphere = "south",
                          roundRes = 4L,
                          return_df = FALSE,
                          save_format = NULL,
                          n_threads = 1L) {

  if (missing(file_in) || missing(clim_file) || missing(file_out))
    stop("file_in, clim_file, and file_out must all be provided.", call. = FALSE)
  if (!file.exists(clim_file))
    stop("Climatology file does not exist: ", clim_file, call. = FALSE)

  hemisphere <- match.arg(hemisphere, c("south", "north"))
  south <- hemisphere == "south"
  vn <- if (is.null(var_name)) "" else var_name

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
      files = file_in, clim_file = clim_file, file_out = file_out,
      var_name = vn,
      minDuration = as.integer(minDuration),
      minDuration2 = as.integer(minDuration2),
      joinAcrossGaps = joinAcrossGaps,
      maxGap = as.integer(maxGap), maxGap2 = as.integer(maxGap2),
      coldSpells = coldSpells, roundRes = as.integer(roundRes),
      n_threads = as.integer(n_threads),
      category = category, southHemisphere = south
    )
  } else {
    if (!file.exists(file_in))
      stop("SST file does not exist: ", file_in, call. = FALSE)
    hw3_detect_events(
      file_in = file_in, clim_file = clim_file, file_out = file_out,
      var_name = vn,
      minDuration = as.integer(minDuration),
      minDuration2 = as.integer(minDuration2),
      joinAcrossGaps = joinAcrossGaps,
      maxGap = as.integer(maxGap), maxGap2 = as.integer(maxGap2),
      coldSpells = coldSpells, roundRes = as.integer(roundRes),
      n_threads = as.integer(n_threads),
      category = category, southHemisphere = south
    )
  }

  if (!is.null(save_format)) {
    hw3_export(file_out, format = save_format, type = "event")
  }

  if (return_df) {
    return(.read_event_df(file_out))
  }

  invisible(file_out)
}

# Internal helper: read event NetCDF into a tidy data.frame
.read_event_df <- function(event_file) {
  ev <- hw3_read_event_nc(event_file)
  ref_date <- as.Date("1970-01-01") + (ev$ref_date_jd - 2440588L)

  cat_labels <- c("I Moderate", "II Strong", "III Severe", "IV Extreme")
  sea_labels <- c("Summer", "Fall", "Winter", "Spring")

  df <- data.frame(
    event_no             = ev$event_no,
    lon                  = ev$lon,
    lat                  = ev$lat,
    date_start           = ref_date + ev$date_start,
    date_peak            = ref_date + ev$date_peak,
    date_end             = ref_date + ev$date_end,
    duration             = ev$duration,
    intensity_mean       = ev$intensity_mean,
    intensity_max        = ev$intensity_max,
    intensity_var        = ev$intensity_var,
    intensity_cumulative = ev$intensity_cumulative,
    rate_onset           = ev$rate_onset,
    rate_decline         = ev$rate_decline,
    stringsAsFactors     = FALSE
  )

  if (length(ev$category) > 0 && any(ev$category > 0)) {
    df$category <- factor(cat_labels[ev$category], levels = cat_labels)
    df$p_moderate <- ev$p_moderate
    df$p_strong   <- ev$p_strong
    df$p_severe   <- ev$p_severe
    df$p_extreme  <- ev$p_extreme
    df$season     <- factor(sea_labels[ev$season], levels = sea_labels)
  }

  df
}
