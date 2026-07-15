// Per-(blob, day) metric reduction for heatwaveR spatial events.
// Single pass over event voxels. Returns a long list of columns suitable for
// constructing a data.table on the R side.
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
// depth_at_max_delta_m alongside the existing area-only columns (computed
// identically to the has_depth = false case either way, so 3D output is
// unchanged bit-for-bit).

#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cstdint>

using namespace Rcpp;

namespace {
struct DailyAcc {
  double area;
  double sum_delta_area;
  double sum_lon_area;
  double sum_lat_area;
  double lon_min, lon_max;
  double lat_min, lat_max;
  double max_delta;
  int n_cells;
  // Depth-only fields (left at their initial values when !has_depth).
  double volume;
  double depth_min, depth_max;
  double depth_at_max_delta;
  DailyAcc()
    : area(0.0), sum_delta_area(0.0), sum_lon_area(0.0), sum_lat_area(0.0),
      lon_min(R_PosInf), lon_max(R_NegInf),
      lat_min(R_PosInf), lat_max(R_NegInf),
      max_delta(R_NegInf), n_cells(0),
      volume(0.0), depth_min(R_PosInf), depth_max(R_NegInf),
      depth_at_max_delta(NA_REAL) {}
};
} // namespace

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

  std::unordered_map<int64_t, DailyAcc> acc;
  acc.reserve(static_cast<size_t>(n / 4) + 16);

  for (R_xlen_t k = 0; k < n; ++k) {
    int lab = labels[k];
    int t   = it[k];
    int64_t key = (static_cast<int64_t>(lab) << 32) | static_cast<uint32_t>(t);
    DailyAcc& e = acc[key];
    double area = cell_area_lat[iy[k] - 1];
    double lon = lons[ix[k] - 1];
    double lat = lats[iy[k] - 1];
    double d   = delta[k];
    e.area += area;
    e.sum_delta_area += d * area;
    e.sum_lon_area += lon * area;
    e.sum_lat_area += lat * area;
    if (lon < e.lon_min) e.lon_min = lon;
    if (lon > e.lon_max) e.lon_max = lon;
    if (lat < e.lat_min) e.lat_min = lat;
    if (lat > e.lat_max) e.lat_max = lat;
    if (d > e.max_delta) {
      e.max_delta = d;
      if (has_depth) e.depth_at_max_delta = depths[iz[k] - 1];
    }
    ++e.n_cells;
    if (has_depth) {
      int dz = iz[k] - 1;
      double dep = depths[dz];
      double vol = area * layer_thickness_km[dz];
      e.volume += vol;
      if (dep < e.depth_min) e.depth_min = dep;
      if (dep > e.depth_max) e.depth_max = dep;
    }
  }

  const R_xlen_t m = acc.size();

  // Sort keys so output is ordered by (blob, t_idx).
  std::vector<int64_t> keys;
  keys.reserve(m);
  for (auto& kv : acc) keys.push_back(kv.first);
  std::sort(keys.begin(), keys.end());

  IntegerVector out_blob(m), out_t(m), out_n(m);
  NumericVector out_area(m), out_mean_delta(m), out_max_delta(m);
  NumericVector out_centroid_lon(m), out_centroid_lat(m);
  NumericVector out_lon_min(m), out_lon_max(m), out_lat_min(m), out_lat_max(m);
  NumericVector out_volume(has_depth ? m : 0);
  NumericVector out_depth_min(has_depth ? m : 0), out_depth_max(has_depth ? m : 0);
  NumericVector out_depth_at_max(has_depth ? m : 0);

  for (R_xlen_t i = 0; i < m; ++i) {
    const auto& e = acc.find(keys[i])->second;
    int lab = static_cast<int>(keys[i] >> 32);
    int t   = static_cast<int>(keys[i] & 0xFFFFFFFFLL);
    out_blob[i] = lab;
    out_t[i] = t;
    out_area[i] = e.area;
    out_mean_delta[i] = e.sum_delta_area / e.area;
    out_max_delta[i] = e.max_delta;
    out_centroid_lon[i] = e.sum_lon_area / e.area;
    out_centroid_lat[i] = e.sum_lat_area / e.area;
    out_lon_min[i] = e.lon_min;
    out_lon_max[i] = e.lon_max;
    out_lat_min[i] = e.lat_min;
    out_lat_max[i] = e.lat_max;
    out_n[i] = e.n_cells;
    if (has_depth) {
      out_volume[i] = e.volume;
      out_depth_min[i] = e.depth_min;
      out_depth_max[i] = e.depth_max;
      out_depth_at_max[i] = e.depth_at_max_delta;
    }
  }

  (void)nt;  // currently unused; reserved for future bounds checking

  List out = List::create(
    _["blob"] = out_blob,
    _["t_idx"] = out_t,
    _["area_km2"] = out_area,
    _["mean_delta"] = out_mean_delta,
    _["max_delta"] = out_max_delta,
    _["centroid_lon"] = out_centroid_lon,
    _["centroid_lat"] = out_centroid_lat,
    _["lon_min"] = out_lon_min,
    _["lon_max"] = out_lon_max,
    _["lat_min"] = out_lat_min,
    _["lat_max"] = out_lat_max,
    _["n_cells"] = out_n
  );
  if (has_depth) {
    out["volume_km3"] = out_volume;
    out["depth_min_m"] = out_depth_min;
    out["depth_max_m"] = out_depth_max;
    out["depth_at_max_delta_m"] = out_depth_at_max;
  }
  return out;
}
