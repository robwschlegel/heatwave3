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
#' @param minDuration2 Minimum duration for events that also satisfy
#'   \code{threshClim2}. Used only when \code{threshClim2} is supplied.
#' @param joinAcrossGaps Logical. Join events separated by short gaps?
#'   Default \code{TRUE}.
#' @param maxGap Maximum gap length (days) to join across. Default \code{2}.
#' @param maxGap2 Maximum gap length for the secondary \code{threshClim2}
#'   criterion. Used only when \code{threshClim2} is supplied.
#' @param threshClim2 Optional gridded NetCDF logical criterion for the
#'   secondary event pass, equivalent to \code{heatwaveR::detect_event()}'s
#'   \code{threshClim2}. The file must align with \code{file_in}; non-zero and
#'   non-missing values are treated as \code{TRUE}. May be a single NetCDF file,
#'   a vector of files, or a directory of daily \code{.nc}/\code{.nc4} files.
#' @param threshClim2_var_name Name of the secondary criterion variable. If
#'   \code{NULL}, it is auto-detected.
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
#' @param save_file Optional path for an additional output file. The extension
#'   determines the format and must be one of \code{.csv}, \code{.rds}, or
#'   \code{.parquet}. If \code{NULL}, no companion file is written.
#' @param return_df Logical. If \code{TRUE}, return the event table as a
#'   \code{data.frame} in addition to writing the NetCDF file. Default
#'   \code{FALSE} (returns the file path invisibly).
#' @param n_threads Number of OpenMP threads. Default \code{1}.
#' @param skip_bad_files Logical. For multi-file inputs, skip unreadable files
#'   or files with mismatched grids instead of failing. Default \code{FALSE}.
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
detect_event3 <- function(
  file_in,
  clim_file,
  file_out,
  var_name = NULL,
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
  roundRes = 4L,
  return_df = FALSE,
  save_file = NULL,
  n_threads = 1L,
  skip_bad_files = FALSE
) {
  if (missing(file_in) || missing(clim_file) || missing(file_out)) {
    stop(
      "file_in, clim_file, and file_out must all be provided.",
      call. = FALSE
    )
  }
  if (!file.exists(clim_file)) {
    stop("Climatology file does not exist: ", clim_file, call. = FALSE)
  }
  if (!is.null(save_file)) {
    .validate_save_file(save_file)
  }

  hemisphere <- match.arg(hemisphere, c("south", "north"))
  south <- hemisphere == "south"
  vn <- if (is.null(var_name)) "" else var_name
  t2_vn <- if (is.null(threshClim2_var_name)) "" else threshClim2_var_name

  multi_file <- FALSE
  if (length(file_in) == 1 && dir.exists(file_in)) {
    file_in <- sort(list.files(
      file_in,
      pattern = "\\.(nc|nc4)$",
      full.names = TRUE
    ))
    if (length(file_in) == 0) {
      stop(
        "No .nc or .nc4 files found in the specified directory.",
        call. = FALSE
      )
    }
    multi_file <- TRUE
  } else if (length(file_in) > 1) {
    multi_file <- TRUE
  }

  thresh2_files <- NULL
  if (!is.null(threshClim2)) {
    thresh2_files <- .resolve_nc_input(threshClim2, "threshClim2")
    if (!multi_file && length(thresh2_files) > 1) {
      stop(
        "threshClim2 may only contain multiple files when file_in is also ",
        "a directory or vector of files.",
        call. = FALSE
      )
    }
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
      n_threads = as.integer(n_threads),
      category = category,
      southHemisphere = south,
      threshClim2_files = thresh2_files,
      threshClim2_var_name = t2_vn,
      skip_bad_files = skip_bad_files
    )
  } else {
    if (!file.exists(file_in)) {
      stop("SST file does not exist: ", file_in, call. = FALSE)
    }
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
      n_threads = as.integer(n_threads),
      category = category,
      southHemisphere = south,
      threshClim2_file = if (is.null(thresh2_files)) "" else thresh2_files[1],
      threshClim2_var_name = t2_vn
    )
  }

  if (!is.null(save_file)) {
    hw3_export(file_out, file_out = save_file, type = "event")
  }

  if (return_df) {
    return(.read_event_df(file_out))
  }

  invisible(file_out)
}

.resolve_nc_input <- function(path, arg_name) {
  if (length(path) == 1 && dir.exists(path)) {
    files <- sort(list.files(path, pattern = "\\.(nc|nc4)$", full.names = TRUE))
    if (length(files) == 0) {
      stop(
        "No .nc or .nc4 files found in ", arg_name, " directory.",
        call. = FALSE
      )
    }
    return(files)
  }

  missing_files <- path[!file.exists(path)]
  if (length(missing_files) > 0) {
    stop(
      arg_name, " file does not exist: ", missing_files[1],
      call. = FALSE
    )
  }

  path
}

# Internal helper: read event NetCDF into a tidy data.frame
.read_event_df <- function(event_file) {
  ev <- hw3_read_event_nc(event_file)
  .event_df_from_indices(ev, seq_len(ev$nevents))
}

.event_df_from_indices <- function(ev, idx) {
  ref_date <- as.Date("1970-01-01") + (ev$ref_date_jd - 2440588L)

  cat_labels <- c("I Moderate", "II Strong", "III Severe", "IV Extreme")
  sea_labels <- c("Summer", "Fall", "Winter", "Spring")

  df <- data.frame(
    lon = ev$lon[idx],
    lat = ev$lat[idx],
    pixel_index = ev$pixel_index[idx],
    event_no = ev$event_no[idx],
    date_start = ref_date + ev$date_start[idx],
    date_peak = ref_date + ev$date_peak[idx],
    date_end = ref_date + ev$date_end[idx],
    duration = ev$duration[idx],
    intensity_mean = ev$intensity_mean[idx],
    intensity_max = ev$intensity_max[idx],
    intensity_var = ev$intensity_var[idx],
    intensity_cumulative = ev$intensity_cumulative[idx],
    intensity_mean_relThresh = ev$intensity_mean_relThresh[idx],
    intensity_max_relThresh = ev$intensity_max_relThresh[idx],
    intensity_var_relThresh = ev$intensity_var_relThresh[idx],
    intensity_cumulative_relThresh = ev$intensity_cumulative_relThresh[idx],
    intensity_mean_abs = ev$intensity_mean_abs[idx],
    intensity_max_abs = ev$intensity_max_abs[idx],
    intensity_var_abs = ev$intensity_var_abs[idx],
    intensity_cumulative_abs = ev$intensity_cumulative_abs[idx],
    rate_onset = ev$rate_onset[idx],
    rate_decline = ev$rate_decline[idx],
    stringsAsFactors = FALSE
  )

  if (length(ev$category) > 0 && any(ev$category > 0)) {
    df$category <- cat_labels[ev$category[idx]]
    df$p_moderate <- ev$p_moderate[idx]
    df$p_strong <- ev$p_strong[idx]
    df$p_severe <- ev$p_severe[idx]
    df$p_extreme <- ev$p_extreme[idx]
    df$season <- sea_labels[ev$season[idx]]
  }

  df
}
