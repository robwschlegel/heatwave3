#' Detect spatially connected marine heatwave / cold-spell blobs
#'
#' Identifies spatially contiguous marine heatwave (or cold-spell) events
#' ("blobs") by 3-D connected-component labelling of the per-pixel Hobday event
#' mask across longitude, latitude, and time, and reduces each blob to the v1
#' event-table metric vocabulary. This is the heatwave3 implementation of the
#' v1 spatial-events labeller (Smit design note; spatially-connected MHW review,
#' Section 15.2.1).
#'
#' @section How it works:
#' The function reads the per-day product written by
#' \code{detect_event3(daily = "also")}, which already carries the
#' \strong{Hobday duration-filtered per-pixel \code{event} flag} together with
#' \code{temp}, \code{seas}, and \code{thresh}. Per-pixel detection is therefore
#' the user's responsibility (it happened in \code{\link{detect_event3}}); the
#' mask is the \code{event} flag, not a raw daily exceedance. Connected voxels
#' in the 3-D mask are labelled with a union-find algorithm (6-connectivity:
#' face-adjacency in \eqn{(x, y, t)}), and the labelled voxels are reduced to
#' per-event and per-day metrics in a single C++ pass. All speed-critical work
#' is in C++.
#'
#' Polarity (heatwave vs cold-spell) is inferred from the sign of the anomaly
#' \eqn{\Delta = temp - seas}; no \code{coldSpells} switch is needed. Signed
#' metrics retain their sign (cold spells are negative) and ranking uses
#' magnitude, so the same defaults behave correctly for both.
#'
#' @param daily_file Path to the per-day NetCDF product from
#'   \code{detect_event3(name = ..., daily = "also")} (the
#'   \code{<name>_events_daily.nc} file). It must contain the \code{event},
#'   \code{temp}, \code{seas}, and \code{thresh} variables.
#' @param connectivity Integer. Voxel connectivity for labelling. Only \code{6}
#'   (face-adjacency in 3-D) is supported in v1; 10, 18, and 26 are reserved for
#'   v2.
#' @param wrapDateline Logical. Join the first and last longitude columns so
#'   blobs can connect across the antimeridian. Default \code{FALSE}. Polar wrap
#'   is not handled.
#' @param cellAreaMethod Character. Cell-area model. \code{"coslat"} (default)
#'   uses the cosine-of-latitude approximation
#'   \eqn{A = R^2 \Delta\lambda \Delta\phi \cos\phi}; \code{"exact"} uses the
#'   spherical-cap area
#'   \eqn{A = R^2 \Delta\lambda [\sin(\phi + \Delta\phi/2) - \sin(\phi -
#'   \Delta\phi/2)]}. \eqn{R = 6371} km.
#' @param minVoxels Integer. Discard blobs with fewer than this many space-time
#'   voxels. Default \code{1} (no filter).
#' @param minArea Numeric. Discard blobs whose peak daily area (km\eqn{^2}) is
#'   below this. Default \code{0} (no filter). Published reference values for
#'   orientation: a 75th-percentile-of-pooled-areas cut (Scannell et al. 2024),
#'   \eqn{5^\circ \times 5^\circ} (Sun et al. 2023). Single-product spatial-event
#'   counts inherit the product's coastal biases; a non-default filter triggers a
#'   reminder.
#' @param minDuration Integer. Discard blobs shorter than this many days
#'   (inclusive span). Default \code{0} (no filter).
#' @param topN Return only the top N blobs ranked by \code{rankBy}. Default
#'   \code{NULL} (return all).
#' @param rankBy Character. Metric to rank blobs by, applied to its
#'   \strong{magnitude} so MHWs and MCSs rank consistently. One of
#'   \code{"cumI_km2_day"} (default; cumulative intensity-weighted area, column
#'   \code{cumI_km2_d}), \code{"peakArea_km2"} (\code{peak_area_km2}),
#'   \code{"volume_km2_day"} (\code{volume_km2_d}), or \code{"duration"}
#'   (\code{duration_days}).
#' @param return Character vector of components to include: \code{"event"} (one
#'   row per blob), \code{"daily"} (one row per blob-day), and/or \code{"voxel"}
#'   (one row per blob-cell-day; can be large). Default \code{c("event",
#'   "daily")}.
#'
#' @return A named list with the requested components. The \code{event}
#'   data.frame carries the v1 vocabulary: \code{event_no}, \code{date_start},
#'   \code{date_end}, \code{date_peak}, \code{duration_days}, \code{n_voxels},
#'   \code{peak_area_km2}, \code{mean_area_km2}, \code{total_area_km2}
#'   (union of unique cells), \code{volume_km2_d} (sum of daily areas),
#'   \code{cumI_km2_d} (area-weighted cumulative anomaly), \code{mean_intensity},
#'   \code{max_intensity} (signed extreme), \code{peak_severity} (max Hobday
#'   multiplier), \code{frac_moderate}/\code{frac_strong}/\code{frac_severe}/
#'   \code{frac_extreme} (peak-day area fractions by Hobday category),
#'   \code{centroid_lon}/\code{centroid_lat} (peak-day great-circle
#'   area-weighted centroid), \code{bbox_lon_min}/\code{bbox_lon_max}/
#'   \code{bbox_lat_min}/\code{bbox_lat_max}, and \code{rank}.
#'
#' @seealso \code{\link{detect_event3}} (run with \code{daily = "also"}) to
#'   produce the input, \code{\link{ts2clm3}} for the climatology.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' detect_event3(sst_file, name = "X", daily = "also")
#' blobs <- detect_blob3("X_events_daily.nc", topN = 10, cellAreaMethod = "exact",
#'                       return = c("event", "daily", "voxel"))
#' head(blobs$event)
#' }
detect_blob3 <- function(daily_file,
                         connectivity = 6L,
                         wrapDateline = FALSE,
                         cellAreaMethod = c("coslat", "exact"),
                         minVoxels = 1L,
                         minArea = 0,
                         minDuration = 0L,
                         topN = NULL,
                         rankBy = "cumI_km2_day",
                         return = c("event", "daily")) {

  if (missing(daily_file)) {
    stop("daily_file must be provided (a detect_event3(daily = \"also\") ",
         "product).", call. = FALSE)
  }
  if (!file.exists(daily_file)) {
    stop("daily_file does not exist: ", daily_file, call. = FALSE)
  }

  connectivity <- as.integer(connectivity)
  if (!identical(connectivity, 6L)) {
    stop("Only connectivity = 6 is supported in v1 (10/18/26 are reserved ",
         "for v2).", call. = FALSE)
  }

  cellAreaMethod <- match.arg(cellAreaMethod)

  valid_return <- c("event", "daily", "voxel")
  bad_return <- setdiff(return, valid_return)
  if (length(bad_return) > 0) {
    stop("return contains unsupported value: ", bad_return[1], call. = FALSE)
  }
  return <- unique(return)

  rank_map <- c(cumI_km2_day = "cumI_km2_d", peakArea_km2 = "peak_area_km2",
                volume_km2_day = "volume_km2_d", duration = "duration_days")
  if (!rankBy %in% names(rank_map)) {
    stop("rankBy must be one of: ", paste(names(rank_map), collapse = ", "),
         call. = FALSE)
  }
  rank_col <- rank_map[[rankBy]]

  # Read the per-day product (C++ reader).
  dd <- hw3_read_daily_nc(daily_file)
  if (length(dd$event) == 0 || length(dd$temp) == 0 ||
      length(dd$seas) == 0 || length(dd$thresh) == 0) {
    stop("daily_file must contain event, temp, seas, and thresh variables. ",
         "Produce it with detect_event3(..., daily = \"also\").", call. = FALSE)
  }
  nlon <- dd$nlon; nlat <- dd$nlat; ntime <- dd$ntime
  if (nlon < 2 || nlat < 2) {
    stop("detect_blob3 requires at least two longitude and two latitude cells ",
         "to compute cell area.", call. = FALSE)
  }

  # Hobday event flag -> 3-D mask (C++), then 6-connectivity labelling (C++).
  mask <- hw3_blob_mask_from_event(dd$event, nlon, nlat, ntime)
  labels <- label_components_3d_cpp(mask, nlon, nlat, ntime, wrapDateline)
  n_components <- attr(labels, "n_components")

  if (n_components == 0) {
    empty <- list()
    if ("event" %in% return) empty$event <- data.frame()
    if ("daily" %in% return) empty$daily <- data.frame()
    if ("voxel" %in% return) empty$voxel <- data.frame()
    return(empty)
  }

  # Per-latitude cell area (km^2) under the chosen method.
  dlon <- stats::median(abs(diff(dd$lon)))
  dlat <- stats::median(abs(diff(dd$lat)))
  cell_area_lat <- .hw3_cell_area_lat(dd$lat, dlon, dlat, cellAreaMethod)

  # Metric reduction (C++, one pass over labelled voxels).
  want_voxel <- "voxel" %in% return
  red <- hw3_blob_reduce(labels, n_components, dd$temp, dd$seas, dd$thresh,
                         dd$lon, dd$lat, cell_area_lat,
                         nlon, nlat, ntime, want_voxel)

  to_date <- function(t_idx) {
    as.Date(dd$time_jd[t_idx] - 2440588L, origin = "1970-01-01")
  }

  ev <- red$event
  event_df <- data.frame(
    event_no = ev$event_no,
    date_start = to_date(ev$t_start),
    date_end = to_date(ev$t_end),
    date_peak = to_date(ev$t_peak),
    duration_days = ev$duration_days,
    n_voxels = ev$n_voxels,
    peak_area_km2 = ev$peak_area_km2,
    mean_area_km2 = ev$mean_area_km2,
    total_area_km2 = ev$total_area_km2,
    volume_km2_d = ev$volume_km2_d,
    cumI_km2_d = ev$cumI_km2_d,
    mean_intensity = ev$mean_intensity,
    max_intensity = ev$max_intensity,
    peak_severity = ev$peak_severity,
    frac_moderate = ev$frac_moderate,
    frac_strong = ev$frac_strong,
    frac_severe = ev$frac_severe,
    frac_extreme = ev$frac_extreme,
    centroid_lon = ev$centroid_lon,
    centroid_lat = ev$centroid_lat,
    bbox_lon_min = ev$bbox_lon_min,
    bbox_lon_max = ev$bbox_lon_max,
    bbox_lat_min = ev$bbox_lat_min,
    bbox_lat_max = ev$bbox_lat_max,
    stringsAsFactors = FALSE
  )

  dl <- red$daily
  daily_df <- data.frame(
    event_no = dl$event_no,
    date = to_date(dl$t_idx),
    area_km2 = dl$area_km2,
    mean_intensity = dl$mean_intensity,
    max_intensity = dl$max_intensity,
    centroid_lon = dl$centroid_lon,
    centroid_lat = dl$centroid_lat,
    bbox_lon_min = dl$bbox_lon_min,
    bbox_lon_max = dl$bbox_lon_max,
    bbox_lat_min = dl$bbox_lat_min,
    bbox_lat_max = dl$bbox_lat_max,
    n_cells = dl$n_cells,
    stringsAsFactors = FALSE
  )

  voxel_df <- NULL
  if (want_voxel) {
    vx <- red$voxel
    voxel_df <- data.frame(
      event_no = vx$event_no,
      lon = dd$lon[vx$i + 1L],
      lat = dd$lat[vx$j + 1L],
      date = to_date(vx$t_idx + 1L),
      intensity = vx$delta,
      stringsAsFactors = FALSE
    )
  }

  # Polarity: infer from the sign of the anomaly; warn on mixed-sign input.
  sgn <- sign(event_df$mean_intensity)
  if (any(sgn > 0, na.rm = TRUE) && any(sgn < 0, na.rm = TRUE)) {
    warning("Input mixes positive and negative anomalies; detect_blob3 expects ",
            "a single polarity (a stitched MHW + MCS event flag?).",
            call. = FALSE)
  }

  # Filters: minVoxels, minArea (peak area), minDuration. Default = no filter.
  filtered <- FALSE
  keep <- rep(TRUE, nrow(event_df))
  if (minVoxels > 1) { keep <- keep & event_df$n_voxels >= minVoxels; filtered <- TRUE }
  if (minArea > 0) { keep <- keep & event_df$peak_area_km2 >= minArea; filtered <- TRUE }
  if (minDuration > 0) { keep <- keep & event_df$duration_days >= minDuration; filtered <- TRUE }
  if (filtered) {
    event_df <- event_df[keep, , drop = FALSE]
    warning("Spatial-event counts inherit the SST product's coastal and ",
            "resolution biases; min* filters change the count accordingly.",
            call. = FALSE)
  }

  # Rank by magnitude, then top-N.
  if (nrow(event_df) > 0) {
    event_df <- event_df[order(-abs(event_df[[rank_col]])), , drop = FALSE]
    event_df$rank <- seq_len(nrow(event_df))
    if (!is.null(topN) && topN < nrow(event_df)) {
      event_df <- event_df[seq_len(topN), , drop = FALSE]
    }
  }
  keep_ids <- event_df$event_no

  daily_df <- daily_df[daily_df$event_no %in% keep_ids, , drop = FALSE]
  if (!is.null(voxel_df)) {
    voxel_df <- voxel_df[voxel_df$event_no %in% keep_ids, , drop = FALSE]
  }

  rownames(event_df) <- NULL
  rownames(daily_df) <- NULL

  result <- list()
  if ("event" %in% return) result$event <- event_df
  if ("daily" %in% return) result$daily <- daily_df
  if ("voxel" %in% return) { rownames(voxel_df) <- NULL; result$voxel <- voxel_df }
  result
}

# Per-latitude cell area (km^2) for a regular lon-lat grid. dlon_deg, dlat_deg
# are the (positive) grid spacings in degrees; R = 6371 km.
.hw3_cell_area_lat <- function(lat, dlon_deg, dlat_deg,
                               method = c("coslat", "exact")) {
  method <- match.arg(method)
  R <- 6371
  d2r <- pi / 180
  dlon <- abs(dlon_deg) * d2r
  half <- abs(dlat_deg) * d2r / 2
  phi <- lat * d2r
  if (method == "exact") {
    R^2 * dlon * (sin(phi + half) - sin(phi - half))
  } else {
    R^2 * dlon * (2 * half) * cos(phi)
  }
}
