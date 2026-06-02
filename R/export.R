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
#' @param chunk_size Number of rows per CSV or Parquet write chunk. Default
#'   \code{1000000}.
#'
#' @return Invisibly returns the path to the exported file.
#'
#' @details
#' \itemize{
#'   \item \code{"csv"}: Writes a flat CSV file in row chunks.
#'   \item \code{"rds"}: Writes the data frame with \code{saveRDS()}.
#'   \item \code{"parquet"}: Writes a single Apache Parquet file in row groups.
#'     Requires the \pkg{arrow} package.
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
                       type = c("clim", "event"),
                       chunk_size = 1000000L) {

  type <- match.arg(type)
  format <- .validate_save_file(file_out)
  chunk_size <- as.integer(chunk_size)
  if (is.na(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be a positive integer.", call. = FALSE)
  }

  if (format == "csv") {
    .write_hw3_csv(nc_file, file_out, type, chunk_size)
  } else if (format == "rds") {
    hw3_data <- if (type == "clim") .read_clim_df(nc_file) else .read_event_df(nc_file)
    saveRDS(hw3_data, file = file_out)
  } else if (format == "parquet") {
    .write_hw3_parquet(nc_file, file_out, type, chunk_size)
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
  .clim_df_from_pixels(cd, seq_len(cd$nlon * cd$nlat))
}

.clim_df_from_pixels <- function(cd, pixels) {
  npixels <- cd$nlon * cd$nlat
  ndoy <- cd$ndoy

  rows <- vector("list", length(pixels))
  for (i in seq_along(pixels)) {
    px <- pixels[i]
    if (px < 1L || px > npixels) {
      stop("Pixel index out of range.", call. = FALSE)
    }
    ilon <- ((px - 1L) %/% cd$nlat) + 1L
    ilat <- ((px - 1L) %% cd$nlat) + 1L
    offset <- (px - 1L) * ndoy
    idx <- offset + seq_len(ndoy)

    rows[[i]] <- data.frame(
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

.iterate_hw3_chunks <- function(nc_file, type, chunk_size, callback) {
  if (type == "clim") {
    cd <- hw3_read_clim_nc(nc_file)
    npixels <- cd$nlon * cd$nlat
    pixels_per_chunk <- max(1L, floor(chunk_size / cd$ndoy))
    for (start in seq.int(1L, npixels, by = pixels_per_chunk)) {
      end <- min(npixels, start + pixels_per_chunk - 1L)
      callback(.clim_df_from_pixels(cd, start:end))
    }
  } else {
    ev <- hw3_read_event_nc(nc_file)
    if (ev$nevents == 0L) {
      callback(.event_df_from_indices(ev, integer(0)))
      return(invisible(NULL))
    }
    for (start in seq.int(1L, ev$nevents, by = chunk_size)) {
      end <- min(ev$nevents, start + chunk_size - 1L)
      callback(.event_df_from_indices(ev, start:end))
    }
  }
  invisible(NULL)
}

.write_hw3_csv <- function(nc_file, file_out, type, chunk_size) {
  if (file.exists(file_out)) {
    unlink(file_out)
  }
  first <- TRUE
  .iterate_hw3_chunks(nc_file, type, chunk_size, function(chunk) {
    utils::write.table(
      chunk,
      file = file_out,
      sep = ",",
      row.names = FALSE,
      col.names = first,
      append = !first,
      qmethod = "double"
    )
    first <<- FALSE
  })
  invisible(file_out)
}

.write_hw3_parquet <- function(nc_file, file_out, type, chunk_size) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("The 'arrow' package is required for Parquet export. ",
         "Install with: install.packages('arrow')", call. = FALSE)
  }
  if (file.exists(file_out)) {
    unlink(file_out)
  }

  writer <- NULL
  sink <- NULL
  wrote_chunk <- FALSE
  close_writer <- function() {
    if (!is.null(writer)) writer$Close()
    if (!is.null(sink) && "close" %in% names(sink)) sink$close()
  }
  on.exit(close_writer(), add = TRUE)

  .iterate_hw3_chunks(nc_file, type, chunk_size, function(chunk) {
    tab <- arrow::as_arrow_table(chunk)
    if (is.null(writer)) {
      sink <<- arrow::FileOutputStream$create(file_out)
      writer <<- arrow::ParquetFileWriter$create(
        tab$schema,
        sink,
        properties = arrow::ParquetWriterProperties$create(names(chunk))
      )
    }
    writer$WriteTable(tab, chunk_size = max(1L, nrow(chunk)))
    wrote_chunk <<- TRUE
  })

  if (!wrote_chunk) {
    hw3_data <- if (type == "clim") .read_clim_df(nc_file) else .read_event_df(nc_file)
    arrow::write_parquet(hw3_data, file_out)
  }

  invisible(file_out)
}
