# This script documents how inst/extdata/glorys_depth_test.nc and
# inst/extdata/glorys_depth_test_annual/*.nc were built from the full GLORYS12
# Port Elizabeth download (dev/data/, not tracked in git). Not run automatically;
# kept for provenance/reproducibility. Requires cdo and nccopy (netcdf-c) on PATH.

# Source ----------------------------------------------------------------

src <- "dev/data/cmems_mod_glo_phy_my_0.083deg_P1D-m_port_elizabeth_20260707T153618.nc"
# 61 lon x 25 lat x 50 depth x 11688 time (daily, 1993-01-01 to 2024-12-31).

# Pixel selection ---------------------------------------------------------

# A 2x2 open-ocean block was chosen by computing a land/sea mask via
# `cdo timmin` at depth level 5 (5.08 m) and again at level 30 (380 m), then
# searching for a contiguous 2x2 block (plus a 1-pixel buffer ring) with no
# missing values at either depth -- i.e. genuinely open ocean, not just a
# shallow shelf pixel. Grid indices 27:28 (lon) x 2:3 (lat), 1-based, were
# selected: 25.16667-25.25 degrees E, -34.91667 to -34.83333 degrees N.
# Verified with `cdo info` (Miss = 0) across the full time series at all of
# the first 5 depth levels (0.49-5.08 m).

# Single fixture: inst/extdata/glorys_depth_test.nc -----------------------

system(paste(
  "cdo -s -O -z zip9 selindexbox,27,28,2,3 -sellevidx,1/5 -selname,thetao",
  src, "inst/extdata/glorys_depth_test.nc"
))

# cdo writes one HDF5 chunk per timestep by default, which balloons file size
# for a long, spatially tiny series (per-chunk overhead >> data). Rechunk to a
# single chunk covering the whole series and re-deflate.
system(paste(
  "nccopy -u -d9 -s -c thetao:11688,5,2,2",
  "inst/extdata/glorys_depth_test.nc /tmp/glorys_onechunk.nc",
  "&& mv /tmp/glorys_onechunk.nc inst/extdata/glorys_depth_test.nc"
))

nc <- ncdf4::nc_open("inst/extdata/glorys_depth_test.nc", write = TRUE)
ncdf4::ncatt_put(nc, 0, "title",
  paste("GLORYS12V1 potential temperature subset for heatwave3 depth_range",
        "testing: 2x2 open-ocean pixels off Port Elizabeth (SW Indian Ocean),",
        "first 5 depth levels (0.49-5.08 m), full 1993-2024 daily record"))
ncdf4::ncatt_put(nc, 0, "source",
  paste("MERCATOR GLORYS12V1 reanalysis, subset of the Port Elizabeth GLORYS",
        "box used in dev/4d-depth-range-progress.md"))
ncdf4::ncatt_put(nc, 0, "comment",
  paste("Test fixture for heatwave3 depth_range coverage. All 4 pixels are",
        "ocean (no NA/land) at all 5 depth levels across the full time record."))
ncdf4::nc_close(nc)

# Annual split: inst/extdata/glorys_depth_test_annual/*.nc ----------------

# One file per calendar year, for exercising the multi-file (directory of
# daily/annual files) + depth_range ingestion path.
dir.create("inst/extdata/glorys_depth_test_annual", showWarnings = FALSE)
system(paste(
  "cdo -s -O -z zip9 splityear inst/extdata/glorys_depth_test.nc",
  "inst/extdata/glorys_depth_test_annual/glorys_depth_test_"
))

annual_files <- list.files("inst/extdata/glorys_depth_test_annual",
                            full.names = TRUE)
for (f in annual_files) {
  ntime <- as.integer(system(paste("cdo -s ntime", f), intern = TRUE))
  tmp <- tempfile(fileext = ".nc")
  system(paste0("nccopy -u -d9 -s -c thetao:", ntime, ",5,2,2 ", f, " ", tmp))
  file.copy(tmp, f, overwrite = TRUE)
  unlink(tmp)
}
