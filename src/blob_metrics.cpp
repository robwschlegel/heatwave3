// R wrapper for hw3::blob_daily_summary() (blob_core.cpp), the per-(blob, day)
// metric reduction for spatial events. Returns a long list of columns
// suitable for constructing a data.frame on the R side.
//
// Inputs are parallel vectors describing event voxels only:
//   labels[i]  : blob label (>= 1)
//   ix[i]      : x-axis index (1-based)
//   iy[i]      : y-axis index (1-based)
//   it[i]      : time index (1-based)
//   delta[i]   : anomaly value
// plus the lookup tables cell_area_lat (length ny), lons (length nx),
// lats (length ny).
//
// When has_depth is true, three more parallel/lookup inputs are used:
//   iz[i]                          : depth-axis index (1-based)
//   depths (length nz)             : depth coordinate values (metres)
//   layer_thickness_km (length nz) : vertical thickness of each depth level
// and the output list gains volume_km3/depth_min_m/depth_max_m/
// depth_at_max_delta_m alongside the existing area-only columns.

#include <Rcpp.h>
#include "blob_core.h"
#include "heatwave3_types.h"

using namespace Rcpp;

// [[Rcpp::export]]
List blob_daily_summary_cpp(IntegerVector labels,
                            IntegerVector ix,
                            IntegerVector iy,
                            IntegerVector it,
                            NumericVector delta,
                            NumericVector cell_area_lat,
                            NumericVector lons,
                            NumericVector lats,
                            int nt,
                            IntegerVector iz,
                            NumericVector depths,
                            NumericVector layer_thickness_km,
                            bool has_depth) {
  const R_xlen_t n = labels.size();
  if (ix.size() != n || iy.size() != n || it.size() != n || delta.size() != n) {
    stop("Input vectors must all have the same length.");
  }
  if (has_depth && iz.size() != n) {
    stop("iz must have the same length as labels when has_depth is TRUE.");
  }

  hw3::BlobDailySummary s = hw3::blob_daily_summary(
    labels.begin(), ix.begin(), iy.begin(), it.begin(), delta.begin(), n,
    cell_area_lat.begin(), lons.begin(), lats.begin(),
    iz.begin(), depths.begin(), layer_thickness_km.begin(), has_depth);

  (void)nt;  // currently unused; reserved for future bounds checking

  List out = List::create(
    _["blob"] = s.blob,
    _["t_idx"] = s.t_idx,
    _["area_km2"] = s.area_km2,
    _["mean_delta"] = s.mean_delta,
    _["max_delta"] = s.max_delta,
    _["centroid_lon"] = s.centroid_lon,
    _["centroid_lat"] = s.centroid_lat,
    _["lon_min"] = s.lon_min,
    _["lon_max"] = s.lon_max,
    _["lat_min"] = s.lat_min,
    _["lat_max"] = s.lat_max,
    _["n_cells"] = s.n_cells
  );
  if (has_depth) {
    // The core marks "no depth recorded" with NaN; R expects NA.
    NumericVector at_max = wrap(s.depth_at_max_delta_m);
    for (R_xlen_t i = 0; i < at_max.size(); ++i)
      if (hw3::is_na(at_max[i])) at_max[i] = NA_REAL;
    out["volume_km3"] = s.volume_km3;
    out["depth_min_m"] = s.depth_min_m;
    out["depth_max_m"] = s.depth_max_m;
    out["depth_at_max_delta_m"] = at_max;
  }
  return out;
}
