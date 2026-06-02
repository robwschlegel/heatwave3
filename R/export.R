#' Export heatwave3 output to additional formats
#'
#' Reads a heatwave3 NetCDF output file (climatology or events) and writes it
#' to a companion file. The output format is inferred from the extension of
#' \code{file_out}.
#'
#' @param nc_file Path to the heatwave3 NetCDF output file.
#' @param file_out Path for the exported companion file. The extension must be
#'   one of \code{.csv}, \code{.rds}, or \code{.parquet}.
#' @param type Character. One of \code{"clim"} or \code{"event"}, indicating the
#'   type of output file.
#'
#' @return Invisibly returns the path to the exported file.
#'
#' @details
#' \itemize{
#'   \item \code{"csv"}: Writes a flat CSV file. Not recommended for large grids
#'     as it can produce very large files.
#'   \item \code{"rds"}: Writes the data frame with \code{saveRDS()}.
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
#' csv_path <- hw3_export(event_file, file_out = tempfile(fileext = ".csv"),
#'                        type = "event")
#' head(read.csv(csv_path))
#'
#' # Export climatology as RDS
#' hw3_export(clim_file, file_out = tempfile(fileext = ".rds"), type = "clim")
#' }
hw3_export <- function(nc_file, file_out,
                       type = c("clim", "event")) {

  type <- match.arg(type)
  format <- .validate_save_file(file_out)

  if (type == "clim") {
    hw3_data <- .read_clim_df(nc_file)

  } else {
    hw3_data <- .read_event_df(nc_file)
  }

  if (format == "csv") {
    if (nrow(hw3_data) > 1e6)
      message("Note: writing ", nrow(hw3_data), " rows to CSV. ",
              "Consider using 'parquet' for large datasets.")
    utils::write.csv(hw3_data, file_out, row.names = FALSE)
  } else if (format == "rds") {
    saveRDS(hw3_data, file = file_out)
  } else if (format == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE))
      stop("The 'arrow' package is required for Parquet export. ",
           "Install with: install.packages('arrow')", call. = FALSE)
    arrow::write_parquet(hw3_data, file_out)
  }

  invisible(file_out)
}

# Internal helper: validate and infer export format from a path
.validate_save_file <- function(save_file) {
  if (missing(save_file) || is.null(save_file))
    stop("save_file must be provided.", call. = FALSE)

  out_dir <- dirname(save_file)
  if (!dir.exists(out_dir))
    stop("Output directory does not exist: ", out_dir, call. = FALSE)

  ext <- tolower(tools::file_ext(save_file))
  if (!ext %in% c("csv", "rds", "parquet")) {
    stop("save_file must end in .csv, .rds, or .parquet.", call. = FALSE)
  }

  ext
}

# Internal helper: read climatology NetCDF into a long data.frame
.read_clim_df <- function(clim_file) {
  cd <- hw3_read_clim_nc(clim_file)
  npixels <- cd$nlon * cd$nlat
  ndoy <- cd$ndoy

  rows <- vector("list", npixels)
  for (px in seq_len(npixels)) {
    ilon <- ((px - 1L) %/% cd$nlat) + 1L
    ilat <- ((px - 1L) %% cd$nlat) + 1L
    offset <- (px - 1L) * ndoy
    idx <- offset + seq_len(ndoy)

    rows[[px]] <- data.frame(
      lon = cd$lon[ilon],
      lat = cd$lat[ilat],
      doy = seq_len(ndoy),
      seas = cd$seas[idx],
      thresh = cd$thresh[idx],
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}
