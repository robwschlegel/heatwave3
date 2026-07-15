#' Read or export a heatwave3 NetCDF product
#'
#' Reads any heatwave3 NetCDF output (climatology, events, daily series, or
#' proto-events) using the fast C++ reader, and either returns it as a long
#' \code{data.frame} for examination or writes it to a flat companion file.
#' The product type is detected automatically from the file.
#'
#' @param file Path to a heatwave3 NetCDF output file (any of the products
#'   written by \code{\link{ts2clm3}} or \code{\link{detect_event3}}).
#' @param file_out Optional path for a flat companion file. When supplied, the
#'   product is written there and the path is returned invisibly; the extension
#'   must be one of \code{.csv}, \code{.rds}, or \code{.parquet}. When
#'   \code{NULL} (the default), the product is returned as a \code{data.frame}.
#' @param vars Optional character vector of data variables to read. \code{NULL}
#'   (the default) reads all variables present. For example \code{c("seas",
#'   "thresh")} for a climatology, or \code{c("temp", "category")} for a daily
#'   product. An unknown name raises an error listing the available variables.
#' @param lon_range,lat_range Optional numeric \code{c(min, max)} windows. Only
#'   the overlapping grid cells are read from disk (a true hyperslab read, so
#'   the rest of the file is never loaded). For the events product these filter
#'   events by their per-event coordinate.
#' @param time_range Optional \code{c("start", "end")} dates (e.g.
#'   \code{c("2015-01-01", "2015-12-31")}). Subsets the time dimension of daily
#'   and proto-event products; for events it keeps events overlapping the range.
#'   Ignored for a climatology (which has a day-of-year axis, not time).
#' @param depth_range Optional numeric \code{c(min_depth, max_depth)} window
#'   (metres). Only meaningful for a depth-resolved climatology or daily/
#'   proto-event product (written with \code{ts2clm3(depth_range = ...)}); for
#'   events it filters by each event's \code{depth}. Ignored (with a warning)
#'   if the file has no depth dimension/column.
#' @param n Optional integer. Return only the first \code{n} rows. For gridded
#'   products this also caps how much is read from disk, so it is a cheap preview
#'   even for very large files. \code{NULL} (the default) returns the whole
#'   (possibly subset) table, with a warning if an unsubset whole read is large.
#' @param type Optional override of the auto-detected product, one of
#'   \code{"clim"}, \code{"event"}, \code{"daily"}, or \code{"protoevents"}.
#'   Normally left \code{NULL}.
#' @param chunk_size Number of rows per CSV or Parquet write chunk for a whole
#'   (unsubset) export. Default \code{1000000}.
#'
#' @return If \code{file_out} is supplied, invisibly returns its path. Otherwise
#'   returns a \code{data.frame} (the whole product, or the requested subset).
#'
#' @details
#' Any of \code{vars}, \code{lon_range}, \code{lat_range}, \code{time_range},
#' \code{depth_range}, or \code{n} triggers a streaming hyperslab read: for the
#' gridded products only the requested window and variables are read from the
#' NetCDF, so memory and I/O scale with the subset, not the file. A whole
#' (unsubset) export to \code{file_out} is written in row chunks.
#'
#' Writing formats (chosen by the \code{file_out} extension): \code{.csv},
#' \code{.rds}, or \code{.parquet} (the last requires the \pkg{arrow} package).
#'
#' @export
#'
#' @examples
#' \donttest{
#' sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
#' stem <- file.path(tempdir(), "demo")
#' detect3(sst_file, name = stem,
#'         climatologyPeriod = c("1982-01-01", "2011-12-31"), daily = "also")
#'
#' # First rows of the events product (cheap preview)
#' head(hw3_export(paste0(stem, "_events.nc"), n = 10))
#'
#' # Only seas/thresh for a lon-lat window of the climatology
#' hw3_export(paste0(stem, "_clim.nc"), vars = c("seas", "thresh"),
#'            lon_range = c(26.5, 26.6))
#'
#' # A date window and chosen variables of the daily product, to CSV
#' hw3_export(paste0(stem, "_events_daily.nc"),
#'            vars = c("temp", "category"), time_range = c("2010-01-01", "2010-12-31"),
#'            file_out = tempfile(fileext = ".csv"))
#' }
hw3_export <- function(file, file_out = NULL, vars = NULL,
                       lon_range = NULL, lat_range = NULL, time_range = NULL,
                       depth_range = NULL,
                       n = NULL, type = NULL, chunk_size = 1000000L) {

  if (!file.exists(file)) {
    stop("File does not exist: ", file, call. = FALSE)
  }

  product <- if (is.null(type)) hw3_file_meta(file)$product else .type_to_product(type)
  ctype <- .product_to_chunktype(product)
  subset_requested <- !is.null(vars) || !is.null(lon_range) ||
    !is.null(lat_range) || !is.null(time_range) || !is.null(depth_range) ||
    !is.null(n)

  # Return a data.frame
  if (is.null(file_out)) {
    if (subset_requested) {
      return(.read_subset_df(file, product, vars, lon_range, lat_range,
                             time_range, depth_range, n))
    }
    nrows <- hw3_file_meta(file)$nrows
    if (nrows > .HW3_PREVIEW_LIMIT) {
      warning("Reading ", format(nrows, big.mark = ","), " rows into memory; ",
              "pass vars/lon_range/lat_range/time_range/depth_range/n to read ",
              "a subset, or file_out to export instead.", call. = FALSE)
    }
    return(.read_product_df(file, product))
  }

  # Write a flat companion file
  format <- .validate_save_file(file_out)

  if (subset_requested) {
    # The subset is bounded by the user; read it and write in one pass.
    .write_df(.read_subset_df(file, product, vars, lon_range, lat_range,
                              time_range, depth_range, n),
              file_out, format)
    return(invisible(file_out))
  }

  chunk_size <- as.integer(chunk_size)
  if (is.na(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be a positive integer.", call. = FALSE)
  }
  if (format == "csv") {
    .write_hw3_csv(file, file_out, ctype, chunk_size)
  } else if (format == "rds") {
    saveRDS(.read_product_df(file, product), file = file_out)
  } else if (format == "parquet") {
    .write_hw3_parquet(file, file_out, ctype, chunk_size)
  }

  invisible(file_out)
}

# Write a single in-memory data.frame to a flat file by format.
.write_df <- function(df, file_out, format) {
  if (format == "csv") {
    utils::write.csv(df, file_out, row.names = FALSE)
  } else if (format == "rds") {
    saveRDS(df, file_out)
  } else if (format == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("The 'arrow' package is required for Parquet export.", call. = FALSE)
    }
    arrow::write_parquet(df, file_out)
  }
  invisible(file_out)
}

# Streaming subset read -> long data.frame. Gridded products use the C++
# hyperslab reader; events (sparse, small) are filtered after a whole read.
.read_subset_df <- function(file, product, vars = NULL, lon_range = NULL,
                            lat_range = NULL, time_range = NULL,
                            depth_range = NULL, n = NULL) {
  if (product == "events") {
    return(.read_event_subset_df(file, vars, lon_range, lat_range, time_range,
                                 depth_range, n))
  }
  t_jd <- NULL
  if (!is.null(time_range)) {
    if (product == "climatology") {
      warning("time_range is ignored for a climatology (no time dimension).",
              call. = FALSE)
    } else {
      d <- as.Date(time_range)
      if (anyNA(d)) stop("time_range must be dates, e.g. c(\"2010-01-01\", \"2010-12-31\").",
                         call. = FALSE)
      t_jd <- as.integer(d) + 2440588L
    }
  }
  sub <- hw3_read_subset(file, lon_range = lon_range, lat_range = lat_range,
                         t_jd_range = t_jd, depth_range = depth_range, vars = vars,
                         max_rows = if (is.null(n)) -1L else as.integer(n))
  df <- .subset_to_df(sub)
  if (!is.null(n)) df <- utils::head(df, as.integer(n))
  df
}

# Build a long data.frame from the C++ hyperslab list (gridded products).
# has_depth: an extra depth axis sits between lat and the third (doy/time)
# axis in the flat [lon][lat][depth][third] hyperslab, matching the on-disk
# layout (see read_subset_netcdf()); reduces to the plain [lon][lat][third]
# decode below when has_depth is FALSE.
.subset_to_df <- function(sub) {
  nlon <- sub$nlon; nlat <- sub$nlat; n3 <- sub$n3
  has_depth <- isTRUE(sub$has_depth)
  ndepth <- if (has_depth) sub$ndepth else 1L
  vnames <- names(sub$data)
  third_is_t <- identical(sub$third_name, "t")

  if (nlon == 0L || nlat == 0L || (has_depth && ndepth == 0L) || n3 == 0L) {
    df <- data.frame(lon = numeric(0), lat = numeric(0))
    if (has_depth) df$depth <- numeric(0)
    df[[if (third_is_t) "t" else "doy"]] <-
      if (third_is_t) as.Date(character(0)) else integer(0)
    for (v in vnames) df[[v]] <- numeric(0)
    return(.finalize_subset_cols(df, vnames))
  }

  np <- nlon * nlat * ndepth
  if (has_depth) {
    lon_col <- rep(sub$lon, each = nlat * ndepth * n3)
    lat_col <- rep(rep(sub$lat, each = ndepth * n3), times = nlon)
    depth_col <- rep(rep(sub$depth, each = n3), times = nlon * nlat)
  } else {
    lon_col <- rep(sub$lon, each = nlat * n3)
    lat_col <- rep(rep(sub$lat, each = n3), times = nlon)
  }
  if (third_is_t) {
    third_col <- rep(as.Date("1970-01-01") + (sub$third - 2440588L), times = np)
    df <- data.frame(lon = lon_col, lat = lat_col, t = third_col, stringsAsFactors = FALSE)
  } else {
    df <- data.frame(lon = lon_col, lat = lat_col, doy = rep(sub$third, times = np),
                     stringsAsFactors = FALSE)
  }
  for (v in vnames) df[[v]] <- sub$data[[v]]
  if (has_depth) df$depth <- depth_col
  .finalize_subset_cols(df, vnames)
}

.finalize_subset_cols <- function(df, vnames) {
  for (v in intersect(c("threshCriterion", "durationCriterion", "event"), vnames)) {
    df[[v]] <- df[[v]] != 0L
  }
  for (v in intersect(c("event_no", "category"), vnames)) {
    x <- df[[v]]; x[x == 0L] <- NA_integer_; df[[v]] <- x
  }
  df
}

# Events: filter a whole (small) read by coordinate/time/depth and select columns.
.read_event_subset_df <- function(file, vars = NULL, lon_range = NULL,
                                  lat_range = NULL, time_range = NULL,
                                  depth_range = NULL, n = NULL) {
  df <- .read_event_df(file)
  keep <- rep(TRUE, nrow(df))
  if (!is.null(lon_range)) keep <- keep & df$lon >= min(lon_range) & df$lon <= max(lon_range)
  if (!is.null(lat_range)) keep <- keep & df$lat >= min(lat_range) & df$lat <= max(lat_range)
  if (!is.null(time_range)) {
    d <- as.Date(time_range)
    if (anyNA(d)) stop("time_range must be dates.", call. = FALSE)
    keep <- keep & df$date_end >= min(d) & df$date_start <= max(d)
  }
  if (!is.null(depth_range)) {
    if (is.null(df$depth)) {
      warning("depth_range is ignored for events with no depth column ",
              "(not from a depth-resolved run).", call. = FALSE)
    } else {
      keep <- keep & df$depth >= min(depth_range) & df$depth <= max(depth_range)
    }
  }
  df <- df[keep, , drop = FALSE]
  if (!is.null(vars)) {
    miss <- setdiff(vars, names(df))
    if (length(miss)) {
      stop("Variable(s) not in events: ", paste(miss, collapse = ", "),
           ". Available: ", paste(names(df), collapse = ", "), call. = FALSE)
    }
    df <- df[, union(c("lon", "lat", "event_no"), vars), drop = FALSE]
  }
  if (!is.null(n)) df <- utils::head(df, as.integer(n))
  rownames(df) <- NULL
  df
}

# Row count above which whole-table reads emit a warning / inline previews are
# skipped. Generous enough for per-slice workflows.
.HW3_PREVIEW_LIMIT <- 5000000L

# Output-path helpers: a single 'name' stem maps to the four product files.
.hw3_clim_path   <- function(name) paste0(name, "_clim.nc")
.hw3_events_path <- function(name) paste0(name, "_events.nc")
.hw3_daily_path  <- function(name) paste0(name, "_events_daily.nc")
.hw3_proto_path  <- function(name) paste0(name, "_protoevents.nc")

# Map the stored product label / a 'type' override to the chunk-writer type.
.product_to_chunktype <- function(product) {
  switch(product,
    climatology = "clim",
    events = "event",
    daily = "daily",
    protoevents = "daily",
    stop("Unrecognised product: ", product, call. = FALSE))
}

.type_to_product <- function(type) {
  type <- match.arg(type, c("clim", "event", "daily", "protoevents"))
  switch(type, clim = "climatology", event = "events",
         daily = "daily", protoevents = "protoevents")
}

.read_product_df <- function(file, product) {
  switch(product,
    climatology = .read_clim_df(file),
    events = .read_event_df(file),
    daily = .read_daily_df(file),
    protoevents = .read_daily_df(file),
    stop("Unrecognised product: ", product, call. = FALSE))
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
  ndepth <- max(1L, cd$ndepth)
  .clim_df_from_pixels(cd, seq_len(cd$nlon * cd$nlat * ndepth))
}

# pixel = (ilon * nlat + ilat) * ndepth + idepth (0-based), matching the C++
# pixel-major layout. Reduces to the plain ilon*nlat+ilat decode when
# ndepth == 1 (ordinary 3D climatology).
.clim_df_from_pixels <- function(cd, pixels) {
  ndepth <- max(1L, cd$ndepth)
  npixels <- cd$nlon * cd$nlat * ndepth
  ndoy <- cd$ndoy
  has_depth <- ndepth > 1L && length(cd$depth) > 0L

  rows <- vector("list", length(pixels))
  for (i in seq_along(pixels)) {
    px <- pixels[i]
    if (px < 1L || px > npixels) {
      stop("Pixel index out of range.", call. = FALSE)
    }
    p0 <- px - 1L
    idepth <- p0 %% ndepth
    lonlat0 <- p0 %/% ndepth
    ilon <- (lonlat0 %/% cd$nlat) + 1L
    ilat <- (lonlat0 %% cd$nlat) + 1L
    offset <- p0 * ndoy
    idx <- offset + seq_len(ndoy)

    df <- data.frame(
      lon = cd$lon[ilon],
      lat = cd$lat[ilat],
      doy = seq_len(ndoy),
      seas = cd$seas[idx],
      thresh = cd$thresh[idx],
      stringsAsFactors = FALSE
    )
    if (has_depth) df$depth <- cd$depth[idepth + 1L]
    rows[[i]] <- df
  }

  do.call(rbind, rows)
}

# Internal helper: read a per-day NetCDF (protoEvent / daily) into a long
# data.frame, one row per pixel per time step. Columns present depend on the
# file's variable set.
.read_daily_df <- function(daily_file) {
  dd <- hw3_read_daily_nc(daily_file)
  ndepth <- max(1L, dd$ndepth)
  .daily_df_from_pixels(dd, seq_len(dd$nlon * dd$nlat * ndepth))
}

# pixel = (ilon * nlat + ilat) * ndepth + idepth (0-based), matching the C++
# pixel-major layout -- see .clim_df_from_pixels(). Reduces to the plain
# ilon*nlat+ilat decode when ndepth == 1 (ordinary 3D daily/protoevents).
.daily_df_from_pixels <- function(dd, pixels) {
  ndepth <- max(1L, dd$ndepth)
  npixels <- dd$nlon * dd$nlat * ndepth
  ntime <- dd$ntime
  has_depth <- ndepth > 1L && length(dd$depth) > 0L
  date <- as.Date("1970-01-01") + (dd$time_jd - 2440588L)
  has <- function(x) length(dd[[x]]) > 0L

  rows <- vector("list", length(pixels))
  for (i in seq_along(pixels)) {
    px <- pixels[i]
    if (px < 1L || px > npixels) {
      stop("Pixel index out of range.", call. = FALSE)
    }
    p0 <- px - 1L
    idepth <- p0 %% ndepth
    lonlat0 <- p0 %/% ndepth
    ilon <- (lonlat0 %/% dd$nlat) + 1L
    ilat <- (lonlat0 %% dd$nlat) + 1L
    idx <- p0 * ntime + seq_len(ntime)

    df <- data.frame(
      lon = dd$lon[ilon],
      lat = dd$lat[ilat],
      t = date,
      temp = dd$temp[idx],
      seas = dd$seas[idx],
      thresh = dd$thresh[idx],
      stringsAsFactors = FALSE
    )
    if (has("threshCriterion")) df$threshCriterion <- dd$threshCriterion[idx] != 0L
    if (has("durationCriterion")) df$durationCriterion <- dd$durationCriterion[idx] != 0L
    if (has("intensity")) df$intensity <- dd$intensity[idx]
    if (has("event")) df$event <- dd$event[idx] != 0L
    if (has("event_no")) {
      eno <- dd$event_no[idx]
      eno[eno == 0L] <- NA_integer_
      df$event_no <- eno
    }
    if (has("category")) {
      ct <- dd$category[idx]
      ct[ct == 0L] <- NA_integer_
      df$category <- ct
    }
    if (has_depth) df$depth <- dd$depth[idepth + 1L]
    rows[[i]] <- df
  }

  do.call(rbind, rows)
}

.iterate_hw3_chunks <- function(nc_file, type, chunk_size, callback) {
  if (type == "clim") {
    cd <- hw3_read_clim_nc(nc_file)
    npixels <- cd$nlon * cd$nlat * max(1L, cd$ndepth)
    pixels_per_chunk <- max(1L, floor(chunk_size / cd$ndoy))
    for (start in seq.int(1L, npixels, by = pixels_per_chunk)) {
      end <- min(npixels, start + pixels_per_chunk - 1L)
      callback(.clim_df_from_pixels(cd, start:end))
    }
  } else if (type == "daily") {
    dd <- hw3_read_daily_nc(nc_file)
    npixels <- dd$nlon * dd$nlat * max(1L, dd$ndepth)
    pixels_per_chunk <- max(1L, floor(chunk_size / max(1L, dd$ntime)))
    for (start in seq.int(1L, npixels, by = pixels_per_chunk)) {
      end <- min(npixels, start + pixels_per_chunk - 1L)
      callback(.daily_df_from_pixels(dd, start:end))
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
    hw3_data <- switch(type,
      clim = .read_clim_df(nc_file),
      daily = .read_daily_df(nc_file),
      .read_event_df(nc_file))
    arrow::write_parquet(hw3_data, file_out)
  }

  invisible(file_out)
}

# Post-compute console summary: head + tail + light stats + a reminder to use
# hw3_export(). Reads via the C++ reader; skips the inline preview for very
# large products (shows metadata only).
.hw3_console_summary <- function(file, quiet = FALSE, n = 5L) {
  if (isTRUE(quiet) || !file.exists(file)) return(invisible(NULL))
  meta <- hw3_file_meta(file)
  product <- meta$product
  label <- switch(product,
    climatology = "Climatology", events = "Events",
    daily = "Daily series", protoevents = "Proto-events", product)
  bar <- strrep("-", 66)

  cat("\n", bar, "\n", label, " written to: ", file, "\n", sep = "")
  cat("Rows (long format): ", format(meta$nrows, big.mark = ","), sep = "")
  if (product != "events") {
    cat("   grid: ", meta$nlon, " lon x ", meta$nlat, " lat", sep = "")
  }
  cat("\n")

  if (meta$nrows > .HW3_PREVIEW_LIMIT) {
    cat("(too large to preview inline; use hw3_export() to read or export)\n",
        bar, "\n", sep = "")
    return(invisible(NULL))
  }

  df <- .read_product_df(file, product)
  cat("\nHead:\n"); print(utils::head(df, n))
  cat("\nTail:\n"); print(utils::tail(df, n))
  cat("\nSummary:\n"); .hw3_print_stats(df, product)
  cat("\nExamine with  hw3_export(\"", file, "\", n = 20)\n",
      "or export with hw3_export(\"", file,
      "\", file_out = \"out.csv\")  (.csv/.rds/.parquet)\n", bar, "\n", sep = "")
  invisible(NULL)
}

.hw3_print_stats <- function(df, product) {
  rng <- function(x) { x <- x[is.finite(x)]; if (!length(x)) c(NA, NA) else range(x) }
  fmt <- function(v) paste(formatC(v, digits = 4, format = "g"), collapse = " to ")
  cats <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) return(invisible(NULL))
    tb <- table(x)
    cat("  category: ", paste(names(tb), tb, sep = "=", collapse = "  "), "\n", sep = "")
  }
  if (product == "climatology") {
    vp <- unique(df[is.finite(df$seas), c("lon", "lat")])
    cat("  ocean pixels (valid climatology): ", nrow(vp), "\n", sep = "")
    cat("  seas:   ", fmt(rng(df$seas)), "\n", sep = "")
    cat("  thresh: ", fmt(rng(df$thresh)), "\n", sep = "")
  } else if (product == "events") {
    cat("  events: ", nrow(df), "   pixels with events: ",
        nrow(unique(df[, c("lon", "lat")])), "\n", sep = "")
    if (nrow(df) > 0) {
      cat("  dates:  ", format(min(df$date_start)), " to ",
          format(max(df$date_end)), "\n", sep = "")
      cat("  duration (days): ", fmt(rng(df$duration)), "\n", sep = "")
      cat("  intensity_max:   ", fmt(rng(df$intensity_max)), "\n", sep = "")
      if (!is.null(df$category)) cats(df$category)
    }
  } else { # daily / protoevents
    if (!is.null(df$t)) {
      cat("  dates: ", format(min(df$t)), " to ", format(max(df$t)), "\n", sep = "")
    }
    if (!is.null(df$event)) {
      cat("  event-days: ", sum(df$event, na.rm = TRUE), "\n", sep = "")
    }
    if (!is.null(df$threshCriterion)) {
      cat("  threshCriterion days: ", sum(df$threshCriterion, na.rm = TRUE), "\n", sep = "")
    }
    if (!is.null(df$category)) cats(df$category)
  }
  invisible(NULL)
}
