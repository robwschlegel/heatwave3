#' Compute climatology for gridded NetCDF data
#'
#' Calculates seasonal and threshold climatologies for each pixel in a gridded
#' NetCDF file, following the Hobday et al. (2016) methodology. This is the
#' gridded equivalent of \code{heatwaveR::ts2clm()}.
#'
#' @param file_in Path to a single multi-timestep NetCDF file, **or** a character
#'   vector of daily NetCDF file paths, **or** a directory path containing daily
#'   NetCDF files (matched by \code{.nc} or \code{.nc4} extension).
#' @param name Output base name (a path stem). The climatology is written to
#'   \code{paste0(name, "_clim.nc")}, so \code{name = "results/cape_coast"}
#'   produces \code{results/cape_coast_clim.nc}. The directory must exist.
#' @param climatologyPeriod A character vector of length 2 specifying the start and
#'   end dates of the baseline period, for example \code{c("1982-01-01", "2011-12-31")}.
#' @param lon_range Optional numeric vector of length 2: \code{c(min_lon, max_lon)}.
#'   If \code{NULL}, all longitudes are used.
#' @param lat_range Optional numeric vector of length 2: \code{c(min_lat, max_lat)}.
#'   If \code{NULL}, all latitudes are used.
#' @param time_range Optional character vector of length 2: \code{c("start", "end")}.
#'   If \code{NULL}, all time steps are read.
#' @param depth Optional integer depth/level index for 4D data. Default \code{NULL}
#'   (no depth subsetting).
#' @param var_name Name of the SST variable in the NetCDF file. If \code{NULL},
#'   the variable is auto-detected.
#' @param maxPadLength Maximum number of consecutive missing days to interpolate.
#'   Default \code{FALSE} (no interpolation). Set to an integer to enable.
#' @param windowHalfWidth Half-width of the sliding window for climatology
#'   calculation. Default \code{5} (11-day window).
#' @param pctile Percentile for the threshold climatology. Default \code{90}.
#' @param smoothPercentile Logical. Apply rolling mean smoothing to the climatology?
#'   Default \code{TRUE}.
#' @param smoothPercentileWidth Width of the rolling mean window for smoothing.
#'   Default \code{31}.
#' @param var Logical. Compute variance climatology? Default \code{FALSE}.
#' @param detrend Logical. Remove a linear trend from each pixel's time series
#'   before computing the climatology? Default \code{FALSE} (fixed-baseline,
#'   Hobday et al. 2016). Set to \code{TRUE} to apply the detrended-baseline
#'   approach (Jacox et al. 2020).
#' @param roundClm Number of decimal places for rounding. Default \code{4}.
#'   Set to \code{FALSE} to disable.
#' @param n_threads Number of OpenMP threads for parallel computation. Default \code{1}.
#' @param skip_bad_files Logical. For multi-file inputs, skip unreadable files
#'   or files with mismatched grids instead of failing. Default \code{FALSE}.
#' @param quiet Logical. Suppress the post-computation console summary (head,
#'   tail, and summary statistics of the climatology)? Default \code{FALSE}.
#'
#' @return Invisibly returns the path to the climatology NetCDF
#'   (\code{paste0(name, "_clim.nc")}). Use \code{\link{hw3_export}} to read it
#'   into a \code{data.frame} or export it to CSV/RDS/Parquet.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ts2clm3(file_in = "path/to/sst.nc",
#'         name = file.path(tempdir(), "cape_coast"),
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"),
#'         lon_range = c(25, 26), lat_range = c(-34, -33))
#' }
ts2clm3 <- function(file_in, name,
                    climatologyPeriod,
                    lon_range = NULL,
                    lat_range = NULL,
                    time_range = NULL,
                    depth = NULL,
                    var_name = NULL,
                    maxPadLength = FALSE,
                    windowHalfWidth = 5L,
                    pctile = 90,
                    smoothPercentile = TRUE,
                    smoothPercentileWidth = 31L,
                    var = FALSE,
                    detrend = FALSE,
                    roundClm = 4L,
                    n_threads = 1L,
                    skip_bad_files = FALSE,
                    quiet = FALSE) {

  if (missing(file_in) || missing(name))
    stop("Both file_in and name must be provided.", call. = FALSE)
  if (missing(climatologyPeriod) || length(climatologyPeriod) != 2)
    stop("climatologyPeriod must be a character vector of length 2.", call. = FALSE)

  file_out <- .hw3_clim_path(name)
  out_dir <- dirname(file_out)
  if (!dir.exists(out_dir))
    stop("Output directory does not exist: ", out_dir, call. = FALSE)

  pad <- if (is.numeric(maxPadLength)) as.integer(maxPadLength) else 0L
  rnd <- if (is.numeric(roundClm)) as.integer(roundClm) else 0L
  vn <- if (is.null(var_name)) "" else var_name
  dp <- if (is.null(depth)) -1L else as.integer(depth)

  # Determine input mode: single file, vector of files, or directory
  multi_file <- FALSE
  if (length(file_in) == 1 && dir.exists(file_in)) {
    # Directory: list all .nc and .nc4 files
    file_in <- sort(list.files(file_in, pattern = "\\.(nc|nc4)$",
                               full.names = TRUE))
    if (length(file_in) == 0)
      stop("No .nc or .nc4 files found in the specified directory.", call. = FALSE)
    multi_file <- TRUE
  } else if (length(file_in) > 1) {
    multi_file <- TRUE
  } else {
    if (!file.exists(file_in))
      stop("Input file does not exist: ", file_in, call. = FALSE)
  }

  if (multi_file) {
    hw3_compute_clim_multi(
      files = file_in,
      file_out = file_out,
      var_name = vn,
      climatologyPeriod = climatologyPeriod,
      lon_range = lon_range,
      lat_range = lat_range,
      depth = dp,
      maxPadLength = pad,
      windowHalfWidth = as.integer(windowHalfWidth),
      pctile = as.double(pctile),
      smoothPercentile = smoothPercentile,
      smoothPercentileWidth = as.integer(smoothPercentileWidth),
      compute_var = var,
      roundClm = rnd,
      n_threads = as.integer(n_threads),
      detrend = detrend,
      skip_bad_files = skip_bad_files
    )
  } else {
    hw3_compute_clim(
      file_in = file_in,
      file_out = file_out,
      var_name = vn,
      climatologyPeriod = climatologyPeriod,
      lon_range = lon_range,
      lat_range = lat_range,
      time_range = time_range,
      depth = dp,
      maxPadLength = pad,
      windowHalfWidth = as.integer(windowHalfWidth),
      pctile = as.double(pctile),
      smoothPercentile = smoothPercentile,
      smoothPercentileWidth = as.integer(smoothPercentileWidth),
      compute_var = var,
      roundClm = rnd,
      n_threads = as.integer(n_threads),
      detrend = detrend
    )
  }

  .hw3_console_summary(file_out, quiet = quiet)
  invisible(file_out)
}
