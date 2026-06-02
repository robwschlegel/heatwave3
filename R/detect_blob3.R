#' Detect spatially connected marine heatwave blobs
#'
#' Identifies spatially contiguous marine heatwave (or cold-spell) events
#' ("blobs") by performing 3D connected-component labeling on the exceedance
#' mask across longitude, latitude, and time.
#'
#' @section How it works:
#' This function does \strong{not} read a pre-computed event file. Instead, it
#' takes the raw SST data and the climatology, and internally constructs a 3D
#' boolean exceedance mask by comparing each pixel's daily SST against its
#' threshold climatology (\code{temp > thresh} for heatwaves, or
#' \code{temp < thresh} for cold-spells). Connected voxels in this 3D mask are
#' then labeled as coherent spatial events using a union-find algorithm.
#'
#' This is the same approach used by \code{heatwaveR::detect_blob3()}. The
#' reason an event file from \code{\link{detect_event3}} cannot be used
#' directly is that the event file stores per-pixel event \emph{summaries}
#' (start/end dates, intensity metrics), not the daily per-pixel boolean mask
#' needed for spatial connectivity analysis.
#'
#' @param sst_file Path to the SST NetCDF file (the same file passed to
#'   \code{\link{ts2clm3}}).
#' @param clim_file Path to the climatology NetCDF produced by
#'   \code{\link{ts2clm3}}. The spatial grid of the climatology determines
#'   the region analysed.
#' @param var_name Name of the SST variable in \code{sst_file}. If \code{NULL},
#'   auto-detected from CF attributes.
#' @param connectivity Integer. Voxel connectivity for labeling. Default
#'   \code{6} (face-adjacent in 3D: left/right, up/down, forward/backward in
#'   time).
#' @param wrapDateline Logical. Wrap the longitude axis so that blobs can
#'   connect across the antimeridian? Default \code{FALSE}.
#' @param minVoxels Minimum number of space-time voxels for a blob to be
#'   retained. Default \code{1}. Set higher (e.g. \code{200}) to filter out
#'   small, short-lived events.
#' @param topN Return only the top N blobs ranked by \code{rankBy}. Default
#'   \code{NULL} (return all).
#' @param rankBy Character. Metric to rank blobs by. Default
#'   \code{"cumI_km2_day"} (cumulative intensity-weighted area). Other options
#'   include \code{"peakArea_km2"}, \code{"duration"}, \code{"meanArea_km2"}.
#' @param coldSpells Logical. If \code{TRUE}, detect cold-spell blobs
#'   (\code{temp < thresh}). Requires a climatology computed with a low
#'   percentile (e.g. \code{pctile = 10} in \code{\link{ts2clm3}}).
#'   Default \code{FALSE}.
#' @param return Character vector specifying which components to include in the
#'   output. One or more of:
#'   \describe{
#'     \item{\code{"event"}}{One row per blob with summary metrics (duration,
#'       peak area, cumulative intensity, centroid, etc.)}
#'     \item{\code{"daily"}}{One row per (blob, date) with daily area, mean/max
#'       anomaly, centroid, and bounding box}
#'     \item{\code{"voxel"}}{One row per (blob, lon, lat, date) — the full 3D
#'       footprint. Required for spatial footprint maps and persistence
#'       analysis. Can be large.}
#'   }
#'   Default \code{c("event", "daily")}.
#' @param skip_bad_files Logical. For multi-file SST inputs, skip unreadable
#'   files or files with mismatched grids instead of failing. Default
#'   \code{FALSE}.
#'
#' @return A named list containing the requested components (\code{event},
#'   \code{daily}, and/or \code{voxel} data.frames).
#'
#' @seealso \code{\link{ts2clm3}} for computing the climatology,
#'   \code{\link{detect_event3}} for per-pixel event detection (different from
#'   spatial blob detection).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sst_file <- "/path/to/sst.nc"
#' clim_file <- tempfile(fileext = ".nc")
#'
#' # Step 1: compute climatology
#' ts2clm3(sst_file, clim_file,
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"),
#'         n_threads = 4)
#'
#' # Step 2: detect spatial blobs (heatwaves)
#' blobs <- detect_blob3(
#'   sst_file  = sst_file,
#'   clim_file = clim_file,
#'   minVoxels = 200,
#'   topN      = 10,
#'   return    = c("event", "daily", "voxel")
#' )
#'
#' head(blobs$event)   # summary per blob
#' head(blobs$daily)   # daily progression
#' head(blobs$voxel)   # full spatial footprint
#'
#' # Step 2b: detect cold-spell blobs (use pctile = 10 climatology)
#' clim_cold <- tempfile(fileext = ".nc")
#' ts2clm3(sst_file, clim_cold,
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"),
#'         pctile = 10, n_threads = 4)
#'
#' mcs <- detect_blob3(
#'   sst_file   = sst_file,
#'   clim_file  = clim_cold,
#'   coldSpells = TRUE,
#'   minVoxels  = 200,
#'   topN       = 6
#' )
#' }
#'
detect_blob3 <- function(sst_file, clim_file,
                         var_name = NULL,
                         connectivity = 6L,
                         wrapDateline = FALSE,
                         minVoxels = 1L,
                         topN = NULL,
                         rankBy = "cumI_km2_day",
                         coldSpells = FALSE,
                         return = c("event", "daily"),
                         skip_bad_files = FALSE) {

  if (missing(sst_file) || missing(clim_file)) {
    stop("sst_file and clim_file must both be provided.", call. = FALSE)
  }
  if (!file.exists(clim_file)) {
    stop("Climatology file does not exist: ", clim_file, call. = FALSE)
  }

  connectivity <- as.integer(connectivity)
  if (!identical(connectivity, 6L)) {
    stop("Only connectivity = 6 is currently supported.", call. = FALSE)
  }

  valid_return <- c("event", "daily", "voxel")
  bad_return <- setdiff(return, valid_return)
  if (length(bad_return) > 0) {
    stop("return contains unsupported value: ", bad_return[1], call. = FALSE)
  }
  return <- unique(return)

  valid_rank <- c("cumI_km2_day", "peakArea_km2", "duration", "meanArea_km2")
  if (!rankBy %in% valid_rank) {
    stop("rankBy must be one of: ", paste(valid_rank, collapse = ", "),
         call. = FALSE)
  }

  # Read SST and climatology
  vn <- if (is.null(var_name)) "" else var_name
  cd <- hw3_read_clim_nc(clim_file)
  if (length(cd$lon) < 2 || length(cd$lat) < 2) {
    stop(
      "detect_blob3 requires at least two longitude and two latitude cells ",
      "to compute cell area.",
      call. = FALSE
    )
  }

  sst_files <- .resolve_nc_input(sst_file, "sst_file")
  if (length(sst_files) > 1 || (length(sst_file) == 1 && dir.exists(sst_file))) {
    gd <- hw3_read_sst_multi(
      sst_files, vn,
      lon_range = c(min(cd$lon), max(cd$lon)),
      lat_range = c(min(cd$lat), max(cd$lat)),
      skip_bad_files = skip_bad_files
    )
  } else {
    gd <- hw3_read_sst(
      sst_files[1], vn,
      lon_range = c(min(cd$lon), max(cd$lon)),
      lat_range = c(min(cd$lat), max(cd$lat))
    )
  }

  nlon <- gd$nlon; nlat <- gd$nlat; ntime <- gd$ntime
  if (nlon != cd$nlon || nlat != cd$nlat) {
    stop(
      "Grid mismatch: SST is ", nlon, " x ", nlat,
      " but climatology is ", cd$nlon, " x ", cd$nlat,
      call. = FALSE
    )
  }
  sst <- gd$sst

  # Map DOY for each time step
  doy_arr <- vapply(gd$time_days, function(jd) {
    # Use the C++ function if available
    hw3_jd_to_doy(jd)
  }, integer(1))

  # Build 3D event mask: 1 = above threshold, 0 = not
  mask <- integer(nlon * nlat * ntime)
  for (px in seq_len(nlon * nlat)) {
    px_offset <- (px - 1L) * ntime
    ilon <- ((px - 1L) %/% nlat) + 1L
    ilat <- ((px - 1L) %% nlat) + 1L
    px_clim_offset <- (px - 1L) * 366

    for (t in seq_len(ntime)) {
      d <- doy_arr[t]
      s_val <- cd$seas[px_clim_offset + d]
      th_val <- cd$thresh[px_clim_offset + d]
      temp_val <- sst[px_offset + t]

      if (is.na(temp_val) || is.na(th_val)) next

      exceed <- if (coldSpells) temp_val < th_val else temp_val > th_val
      if (exceed) {
        # Column-major index for R 3D array [ilon, ilat, t]
        idx <- ilon + (ilat - 1L) * nlon + (t - 1L) * nlon * nlat
        mask[idx] <- 1L
      }
    }
  }

  dim(mask) <- c(nlon, nlat, ntime)

  # 3D connected-component labeling
  labels <- label_components_3d_cpp(as.integer(mask), nlon, nlat, ntime, wrapDateline)
  n_components <- attr(labels, "n_components")

  if (n_components == 0) {
    empty_result <- list()
    if ("event" %in% return) empty_result$event <- data.frame()
    if ("daily" %in% return) empty_result$daily <- data.frame()
    if ("voxel" %in% return) empty_result$voxel <- data.frame()
    return(empty_result)
  }

  # Extract voxel indices for labeled voxels
  which_labeled <- which(labels > 0)
  lab_vec <- labels[which_labeled]

  # Convert linear indices to (i, j, k) — 1-based
  ix <- ((which_labeled - 1L) %% nlon) + 1L
  iy <- (((which_labeled - 1L) %/% nlon) %% nlat) + 1L
  it <- ((which_labeled - 1L) %/% (nlon * nlat)) + 1L

  # Compute anomaly (delta)
  delta <- numeric(length(which_labeled))
  for (v in seq_along(which_labeled)) {
    px <- (ix[v] - 1L) * nlat + iy[v] - 1L
    t_idx <- it[v]
    d <- doy_arr[t_idx]
    temp_val <- sst[px * ntime + t_idx]
    seas_val <- cd$seas[px * 366 + d]
    delta[v] <- temp_val - seas_val
    if (coldSpells) delta[v] <- -delta[v]
  }

  # Cell area by latitude (km^2)
  cell_width_km <- 111.32 * cos(cd$lat * pi / 180) * abs(diff(cd$lon[1:2]))
  cell_height_km <- 111.32 * abs(diff(cd$lat[1:2]))
  cell_area_lat <- cell_width_km * cell_height_km

  # Compute daily summaries
  daily <- blob_daily_summary_cpp(
    lab_vec, ix, iy, it, delta,
    cell_area_lat, cd$lon, cd$lat, ntime
  )
  daily_df <- as.data.frame(daily, stringsAsFactors = FALSE)

  # Convert time indices to dates
  daily_df$date <- as.Date("1970-01-01") + (gd$time_days[daily_df$t_idx] - 2440588L)

  # Aggregate to event-level summaries
  blobs <- unique(daily_df$blob)
  event_list <- lapply(blobs, function(b) {
    sub <- daily_df[daily_df$blob == b, , drop = FALSE]
    data.frame(
      blob = b,
      duration = nrow(sub),
      date_start = min(sub$date),
      date_end = max(sub$date),
      date_peak = sub$date[which.max(sub$area_km2)],
      cumI_km2_day = sum(sub$mean_delta * sub$area_km2),
      peakArea_km2 = max(sub$area_km2),
      meanArea_km2 = mean(sub$area_km2),
      centroid_lon = stats::weighted.mean(sub$centroid_lon, sub$area_km2),
      centroid_lat = stats::weighted.mean(sub$centroid_lat, sub$area_km2),
      stringsAsFactors = FALSE
    )
  })
  event_df <- do.call(rbind, event_list)

  # Build voxel data if requested
  voxel_df <- NULL
  if ("voxel" %in% return) {
    voxel_dates <- as.Date("1970-01-01") + (gd$time_days[it] - 2440588L)
    voxel_df <- data.frame(
      blob = lab_vec,
      lon = cd$lon[ix],
      lat = cd$lat[iy],
      t_idx = it,
      date = voxel_dates,
      delta = delta,
      stringsAsFactors = FALSE
    )
  }

  # Filter by minVoxels
  if (minVoxels > 1) {
    voxel_counts <- table(lab_vec)
    keep <- as.integer(names(voxel_counts[voxel_counts >= minVoxels]))
    event_df <- event_df[event_df$blob %in% keep, , drop = FALSE]
    daily_df <- daily_df[daily_df$blob %in% keep, , drop = FALSE]
    if (!is.null(voxel_df))
      voxel_df <- voxel_df[voxel_df$blob %in% keep, , drop = FALSE]
  }

  # Rank and filter
  if (nrow(event_df) > 0) {
    event_df <- event_df[order(-event_df[[rankBy]]), , drop = FALSE]
    event_df$rank <- seq_len(nrow(event_df))
    if (!is.null(topN) && topN < nrow(event_df)) {
      keep_blobs <- event_df$blob[1:topN]
      event_df <- event_df[1:topN, , drop = FALSE]
      daily_df <- daily_df[daily_df$blob %in% keep_blobs, , drop = FALSE]
      if (!is.null(voxel_df))
        voxel_df <- voxel_df[voxel_df$blob %in% keep_blobs, , drop = FALSE]
    }
  }

  rownames(event_df) <- NULL
  rownames(daily_df) <- NULL

  result <- list()
  if ("event" %in% return) result$event <- event_df
  if ("daily" %in% return) result$daily <- daily_df
  if ("voxel" %in% return) { rownames(voxel_df) <- NULL; result$voxel <- voxel_df }

  result
}
