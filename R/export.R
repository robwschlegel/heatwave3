#' Export heatwave3 output to additional formats
#'
#' Reads a heatwave3 NetCDF output file (climatology or events) and writes it
#' in an additional format. The NetCDF file is always the primary output; this
#' function creates a companion file alongside it.
#'
#' @param nc_file Path to the heatwave3 NetCDF output file.
#' @param format Character. One of \code{"csv"}, \code{"rda"}, or \code{"parquet"}.
#' @param type Character. One of \code{"clim"} or \code{"event"}, indicating the
#'   type of output file.
#'
#' @return Invisibly returns the path to the exported file.
#'
#' @details
#' \itemize{
#'   \item \code{"csv"}: Writes a flat CSV file. Not recommended for large grids
#'     as it can produce very large files.
#'   \item \code{"rda"}: Writes an R data file (.rda) containing a list named
#'     \code{hw3_data}.
#'   \item \code{"parquet"}: Writes an Apache Parquet file. Requires the
#'     \pkg{arrow} package.
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
#' event_file <- tempfile(fileext = ".nc")
#' clim_file <- tempfile(fileext = ".nc")
#'
#' detect3(sst_file, clim_file, event_file,
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"))
#'
#' # Export events as CSV
#' csv_path <- hw3_export(event_file, format = "csv", type = "event")
#' head(read.csv(csv_path))
#'
#' # Export climatology as RDA
#' hw3_export(clim_file, format = "rda", type = "clim")
#' }
hw3_export <- function(nc_file, format = c("csv", "rda", "parquet"),
                       type = c("clim", "event")) {

  format <- match.arg(format)
  type <- match.arg(type)

  if (!requireNamespace("ncdf4", quietly = TRUE))
    stop("ncdf4 package is required for hw3_export()", call. = FALSE)

  nc <- ncdf4::nc_open(nc_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  if (type == "clim") {
    lon <- ncdf4::ncvar_get(nc, "lon")
    lat <- ncdf4::ncvar_get(nc, "lat")
    doy <- ncdf4::ncvar_get(nc, "doy")
    seas   <- ncdf4::ncvar_get(nc, "seas",   collapse_degen = FALSE)
    thresh <- ncdf4::ncvar_get(nc, "thresh", collapse_degen = FALSE)

    # Reshape from [doy, lat, lon] to long data.frame
    rows <- list()
    for (i in seq_along(lon)) {
      for (j in seq_along(lat)) {
        rows[[length(rows) + 1]] <- data.frame(
          lon = lon[i], lat = lat[j], doy = doy,
          seas = seas[, j, i], thresh = thresh[, j, i],
          stringsAsFactors = FALSE
        )
      }
    }
    hw3_data <- do.call(rbind, rows)

  } else {
    # Event file — all variables are 1D along the event dimension
    varnames <- names(nc$var)
    hw3_data <- as.data.frame(
      lapply(varnames, function(v) ncdf4::ncvar_get(nc, v)),
      stringsAsFactors = FALSE
    )
    names(hw3_data) <- varnames

    # Convert date fields to actual dates
    for (dvar in c("date_start", "date_peak", "date_end")) {
      if (dvar %in% names(hw3_data)) {
        du <- ncdf4::ncatt_get(nc, dvar, "units")
        if (du$hasatt) {
          ref <- as.Date(sub("days since ", "", du$value))
          hw3_data[[dvar]] <- ref + hw3_data[[dvar]]
        }
      }
    }

    # Convert integer category and season codes to Hobday et al. (2018)
    # labels (matches what category3() returns).
    cat_labels    <- c("I Moderate", "II Strong", "III Severe", "IV Extreme")
    season_labels <- c("Summer", "Fall", "Winter", "Spring")
    if ("category" %in% names(hw3_data)) {
      ci <- hw3_data$category
      hw3_data$category <- ifelse(ci >= 1 & ci <= 4, cat_labels[ci], NA_character_)
    }
    if ("season" %in% names(hw3_data)) {
      si <- hw3_data$season
      hw3_data$season <- ifelse(si >= 1 & si <= 4, season_labels[si], NA_character_)
    }
  }

  # Write in requested format
  base <- sub("\\.[^.]+$", "", nc_file)

  if (format == "csv") {
    out_file <- paste0(base, ".csv")
    if (nrow(hw3_data) > 1e6)
      message("Note: writing ", nrow(hw3_data), " rows to CSV. ",
              "Consider using 'parquet' for large datasets.")
    utils::write.csv(hw3_data, out_file, row.names = FALSE)
  } else if (format == "rda") {
    out_file <- paste0(base, ".rda")
    save(hw3_data, file = out_file)
  } else if (format == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE))
      stop("The 'arrow' package is required for Parquet export. ",
           "Install with: install.packages('arrow')", call. = FALSE)
    out_file <- paste0(base, ".parquet")
    arrow::write_parquet(hw3_data, out_file)
  }

  invisible(out_file)
}
