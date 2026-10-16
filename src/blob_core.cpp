// Blob detection core, shared by the R and Python bindings.
//
// label_components(): 3D/4D connected-component labelling for heatwave-event
// masks. Two-pass scan with union-find. Face-adjacent connectivity: 6
// (lon/lat/time) when nz == 1, or 8 (lon/lat/depth/time) when nz > 1.
// Optional dateline wrap on the x-axis only (depth is never periodic).
//
// Mask layout matches R's column-major 4D array:
//   m[i, j, kd, kt] -> mask[i + nx*j + nx*ny*kd + nx*ny*nz*kt]
// with i in [0, nx), j in [0, ny), kd in [0, nz), kt in [0, nt).
// Reduces exactly to the plain 3D layout (m[i,j,kt] -> i + nx*j + nx*ny*kt)
// when nz == 1. Returned labels are dense, 1..K. Background voxels carry 0.
//
// blob_daily_summary(): per-(blob, day) metric reduction for spatial events.
// Single pass over event voxels. When has_depth is true, volume and depth
// columns are added alongside the area-only columns (computed identically to
// the has_depth = false case either way, so 3D output is unchanged).

#include "blob_core.h"
#include <unordered_map>
#include <algorithm>
#include <limits>

namespace hw3 {

static inline int uf_find(std::vector<int>& parent, int x) {
  while (parent[x] != x) {
    parent[x] = parent[parent[x]]; // path halving
    x = parent[x];
  }
  return x;
}

static inline void uf_union(std::vector<int>& parent,
                            std::vector<int>& rank_,
                            int a, int b) {
  int ra = uf_find(parent, a);
  int rb = uf_find(parent, b);
  if (ra == rb) return;
  if (rank_[ra] < rank_[rb]) {
    parent[ra] = rb;
  } else if (rank_[ra] > rank_[rb]) {
    parent[rb] = ra;
  } else {
    parent[rb] = ra;
    rank_[ra]++;
  }
}

int label_components(const int* mask, int* labels,
                     int nx, int ny, int nz, int nt,
                     bool wrap_dateline) {
  const int64_t n = static_cast<int64_t>(nx) * ny * nz * nt;

  // Union-find. Reserve to amortise growth.
  std::vector<int> parent;
  std::vector<int> rank_;
  parent.reserve(1024);
  rank_.reserve(1024);
  // sentinel for index 0 (unused; labels start at 1)
  parent.push_back(0);
  rank_.push_back(0);

  int next_label = 1;

  // Pass 1: raster scan in (i, j, kd, kt) order, column-major indexing.
  const int64_t stride_y = nx;
  const int64_t stride_z = static_cast<int64_t>(nx) * ny;
  const int64_t stride_t = static_cast<int64_t>(nx) * ny * nz;

  for (int kt = 0; kt < nt; ++kt) {
    const int64_t base_t = static_cast<int64_t>(kt) * stride_t;
    for (int kd = 0; kd < nz; ++kd) {
      const int64_t base_zt = base_t + static_cast<int64_t>(kd) * stride_z;
      for (int j = 0; j < ny; ++j) {
        const int64_t base_yzt = base_zt + static_cast<int64_t>(j) * stride_y;
        for (int i = 0; i < nx; ++i) {
          const int64_t v = base_yzt + i;
          if (mask[v] == 0) continue;

          int lx = (i > 0) ? labels[v - 1] : 0;
          int ly = (j > 0) ? labels[v - stride_y] : 0;
          int lz = (kd > 0) ? labels[v - stride_z] : 0;
          int lt = (kt > 0) ? labels[v - stride_t] : 0;

          int neigh[4] = {lx, ly, lz, lt};
          int chosen = 0;
          for (int s = 0; s < 4; ++s) {
            if (neigh[s] > 0) {
              if (chosen == 0 || neigh[s] < chosen) chosen = neigh[s];
            }
          }

          if (chosen == 0) {
            // start a new component
            parent.push_back(next_label);
            rank_.push_back(0);
            labels[v] = next_label;
            ++next_label;
          } else {
            labels[v] = chosen;
            for (int s = 0; s < 4; ++s) {
              if (neigh[s] > 0 && neigh[s] != chosen) {
                uf_union(parent, rank_, chosen, neigh[s]);
              }
            }
          }
        }
      }
    }
  }

  // Dateline wrap: union edges across the x boundary, then re-resolve.
  if (wrap_dateline && nx >= 2) {
    for (int kt = 0; kt < nt; ++kt) {
      const int64_t base_t = static_cast<int64_t>(kt) * stride_t;
      for (int kd = 0; kd < nz; ++kd) {
        const int64_t base_zt = base_t + static_cast<int64_t>(kd) * stride_z;
        for (int j = 0; j < ny; ++j) {
          const int64_t base_yzt = base_zt + static_cast<int64_t>(j) * stride_y;
          const int64_t v0 = base_yzt;
          const int64_t v1 = base_yzt + (nx - 1);
          if (labels[v0] > 0 && labels[v1] > 0) {
            uf_union(parent, rank_, labels[v0], labels[v1]);
          }
        }
      }
    }
  }

  // Pass 2: resolve labels via union-find and compactify to 1..K.
  // First, find the canonical root for each label; build a dense remap.
  std::vector<int> remap(parent.size(), 0);
  int next_dense = 1;
  for (int64_t v = 0; v < n; ++v) {
    int lab = labels[v];
    if (lab == 0) continue;
    int root = uf_find(parent, lab);
    int dense = remap[root];
    if (dense == 0) {
      dense = next_dense++;
      remap[root] = dense;
    }
    labels[v] = dense;
  }

  return next_dense - 1;
}

namespace {
constexpr double POS_INF = std::numeric_limits<double>::infinity();
constexpr double NEG_INF = -std::numeric_limits<double>::infinity();

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
      lon_min(POS_INF), lon_max(NEG_INF),
      lat_min(POS_INF), lat_max(NEG_INF),
      max_delta(NEG_INF), n_cells(0),
      volume(0.0), depth_min(POS_INF), depth_max(NEG_INF),
      depth_at_max_delta(std::numeric_limits<double>::quiet_NaN()) {}
};
} // namespace

BlobDailySummary blob_daily_summary(const int* labels, const int* ix,
                                    const int* iy, const int* it,
                                    const double* delta, int64_t n,
                                    const double* cell_area_lat,
                                    const double* lons, const double* lats,
                                    const int* iz, const double* depths,
                                    const double* layer_thickness_km,
                                    bool has_depth) {
  std::unordered_map<int64_t, DailyAcc> acc;
  acc.reserve(static_cast<size_t>(n / 4) + 16);

  for (int64_t k = 0; k < n; ++k) {
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

  const size_t m = acc.size();

  // Sort keys so output is ordered by (blob, t_idx).
  std::vector<int64_t> keys;
  keys.reserve(m);
  for (auto& kv : acc) keys.push_back(kv.first);
  std::sort(keys.begin(), keys.end());

  BlobDailySummary out;
  out.blob.resize(m); out.t_idx.resize(m); out.n_cells.resize(m);
  out.area_km2.resize(m); out.mean_delta.resize(m); out.max_delta.resize(m);
  out.centroid_lon.resize(m); out.centroid_lat.resize(m);
  out.lon_min.resize(m); out.lon_max.resize(m);
  out.lat_min.resize(m); out.lat_max.resize(m);
  if (has_depth) {
    out.volume_km3.resize(m);
    out.depth_min_m.resize(m); out.depth_max_m.resize(m);
    out.depth_at_max_delta_m.resize(m);
  }

  for (size_t i = 0; i < m; ++i) {
    const auto& e = acc.find(keys[i])->second;
    out.blob[i] = static_cast<int>(keys[i] >> 32);
    out.t_idx[i] = static_cast<int>(keys[i] & 0xFFFFFFFFLL);
    out.area_km2[i] = e.area;
    out.mean_delta[i] = e.sum_delta_area / e.area;
    out.max_delta[i] = e.max_delta;
    out.centroid_lon[i] = e.sum_lon_area / e.area;
    out.centroid_lat[i] = e.sum_lat_area / e.area;
    out.lon_min[i] = e.lon_min;
    out.lon_max[i] = e.lon_max;
    out.lat_min[i] = e.lat_min;
    out.lat_max[i] = e.lat_max;
    out.n_cells[i] = e.n_cells;
    if (has_depth) {
      out.volume_km3[i] = e.volume;
      out.depth_min_m[i] = e.depth_min;
      out.depth_max_m[i] = e.depth_max;
      out.depth_at_max_delta_m[i] = e.depth_at_max_delta;
    }
  }

  return out;
}

} // namespace hw3
