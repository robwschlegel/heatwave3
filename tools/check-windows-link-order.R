# Guard the static-library order required by the rwinlib/netcdf 4.9.0 bundle.
# This runs in Windows CI before R CMD check so an accidental reorder produces
# a direct diagnostic instead of the bundle's less obvious duplicate-symbol
# linker failure.

makevars_path <- file.path("src", "Makevars.win")
makevars <- paste(readLines(makevars_path, warn = FALSE), collapse = " ")

xml_position <- regexpr("-lxml2", makevars, fixed = TRUE)[[1L]]
crypto_position <- regexpr("-lcrypto", makevars, fixed = TRUE)[[1L]]

if (xml_position < 1L || crypto_position < 1L) {
  stop("src/Makevars.win must link both libxml2 and libcrypto.", call. = FALSE)
}

if (xml_position > crypto_position) {
  stop(
    paste(
      "src/Makevars.win must link -lxml2 before -lcrypto.",
      "The rwinlib bundle otherwise extracts two incompatible pathtools.o objects."
    ),
    call. = FALSE
  )
}

message("Windows NetCDF static-library order is valid: -lxml2 precedes -lcrypto.")
