#' Categorise marine heatwave events
#'
#' Assigns Hobday et al. (2018) categories (Moderate, Strong, Severe, Extreme)
#' to events based on intensity relative to the threshold-to-seasonal difference.
#'
#' @param event_file Path to the event NetCDF file from \code{\link{detect_event3}}.
#' @param clim_file Path to the climatology NetCDF from \code{\link{ts2clm3}}.
#'   Required for computing category thresholds.
#' @param S Logical. Use the Southern Hemisphere season convention? Default \code{TRUE}.
#' @param name Character string to prepend to event names. Default \code{"Event"}.
#' @param roundVal Decimal places for rounding. Default \code{4}.
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
category3 <- function(event_file, clim_file,
                      S = TRUE, name = "Event",
                      roundVal = 4L) {

  if (!requireNamespace("ncdf4", quietly = TRUE))
    stop("ncdf4 package is required for category3()", call. = FALSE)

  nc_ev <- ncdf4::nc_open(event_file)
  on.exit(ncdf4::nc_close(nc_ev), add = TRUE)

  ev_lon <- ncdf4::ncvar_get(nc_ev, "lon")
  ev_lat <- ncdf4::ncvar_get(nc_ev, "lat")
  i_max <- ncdf4::ncvar_get(nc_ev, "intensity_max")
  dur <- ncdf4::ncvar_get(nc_ev, "duration")
  eno <- ncdf4::ncvar_get(nc_ev, "event_no")

  # Read date_peak and convert to actual dates
  dp <- ncdf4::ncvar_get(nc_ev, "date_peak")
  dp_units <- ncdf4::ncatt_get(nc_ev, "date_peak", "units")$value
  ref_date <- as.Date(sub("days since ", "", dp_units))
  peak_dates <- ref_date + dp

  # Read climatology to compute category thresholds
  nc_cl <- ncdf4::nc_open(clim_file)
  cl_lon <- ncdf4::ncvar_get(nc_cl, "lon")
  cl_lat <- ncdf4::ncvar_get(nc_cl, "lat")
  seas <- ncdf4::ncvar_get(nc_cl, "seas")
  thresh <- ncdf4::ncvar_get(nc_cl, "thresh")
  ncdf4::nc_close(nc_cl)

  nevents <- length(i_max)
  category <- character(nevents)
  p_mod <- p_str <- p_sev <- p_ext <- numeric(nevents)

  for (i in seq_len(nevents)) {
    # Find matching pixel in climatology
    lon_idx <- which.min(abs(cl_lon - ev_lon[i]))
    lat_idx <- which.min(abs(cl_lat - ev_lat[i]))

    # Get peak DOY
    peak_doy <- as.integer(format(peak_dates[i], "%j"))
    # Adjust for non-leap year
    yr <- as.integer(format(peak_dates[i], "%Y"))
    is_leap <- (yr %% 4 == 0 & yr %% 100 != 0) | (yr %% 400 == 0)
    if (!is_leap && peak_doy >= 60) peak_doy <- peak_doy + 1L

    # dims are [doy, lat, lon] in R due to NetCDF reversal
    s <- seas[peak_doy, lat_idx, lon_idx]
    th <- thresh[peak_doy, lat_idx, lon_idx]
    diff <- th - s

    if (is.na(diff) || diff <= 0) {
      category[i] <- NA_character_
      next
    }

    # Category boundaries: 1x, 2x, 3x, 4x the threshold-seas difference
    cats <- c(diff, 2 * diff, 3 * diff, 4 * diff)

    # Proportion of time in each category
    ix <- i_max[i]
    if (ix >= cats[4]) {
      category[i] <- "IV Extreme"
    } else if (ix >= cats[3]) {
      category[i] <- "III Severe"
    } else if (ix >= cats[2]) {
      category[i] <- "II Strong"
    } else {
      category[i] <- "I Moderate"
    }

    p_mod[i] <- round(min(1, max(0, ix / diff)), roundVal)
    p_str[i] <- round(min(1, max(0, (ix - diff) / diff)), roundVal)
    p_sev[i] <- round(min(1, max(0, (ix - 2 * diff) / diff)), roundVal)
    p_ext[i] <- round(min(1, max(0, (ix - 3 * diff) / diff)), roundVal)
  }

  # Season assignment
  month <- as.integer(format(peak_dates, "%m"))
  if (S) {
    season <- ifelse(month %in% c(12, 1, 2), "Summer",
              ifelse(month %in% 3:5, "Fall",
              ifelse(month %in% 6:8, "Winter", "Spring")))
  } else {
    season <- ifelse(month %in% c(12, 1, 2), "Winter",
              ifelse(month %in% 3:5, "Spring",
              ifelse(month %in% 6:8, "Summer", "Fall")))
  }

  data.frame(
    event_no = eno,
    lon = round(ev_lon, 4),
    lat = round(ev_lat, 4),
    peak_date = peak_dates,
    category = category,
    intensity_max = round(i_max, roundVal),
    duration = dur,
    p_moderate = p_mod,
    p_strong = p_str,
    p_severe = p_sev,
    p_extreme = p_ext,
    season = season,
    stringsAsFactors = FALSE
  )
}
