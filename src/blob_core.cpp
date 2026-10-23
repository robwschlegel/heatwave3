// Blob detection core, shared by the R and Python bindings.
//
// label_components(): 3D/4D connected-component labelling for heatwave-event
// masks. Two-pass scan with union-find. `connectivity` counts neighbours in
// the 3 x 3 x 3 lon/lat/time cube on one depth level:
//   6  : faces only (one axis changes)
//   10 : faces + lon/lat diagonals on the same day
//   18 : faces + edges (up to two axes change)
//   26 : faces + edges + corners
// With depth (nz > 1), a voxel also joins the same pixel on the same day one
// level above and below, whatever `connectivity` is (no depth diagonals), so
// 8, 12, 20 or 28 neighbours in all. Optional dateline wrap on the x-axis
// only, with the same neighbours across the boundary (depth is never
// periodic).
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
#include <stdexcept>

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

namespace {

struct Offset { int dx, dy, dz, dt; };

// Neighbour offsets allowed by `connectivity` (see the header comment).
std::vector<Offset> neighbour_offsets(int connectivity, bool has_z) {
  if (connectivity != 6 && connectivity != 10 &&
      connectivity != 18 && connectivity != 26) {
    throw std::runtime_error("connectivity must be 6, 10, 18 or 26.");
  }
  std::vector<Offset> out;
  for (int dt = -1; dt <= 1; ++dt)
    for (int dy = -1; dy <= 1; ++dy)
      for (int dx = -1; dx <= 1; ++dx) {
        int axes = (dx != 0) + (dy != 0) + (dt != 0);
        if (axes == 0) continue;
        bool keep = axes == 1 ||
          (connectivity == 10 && axes == 2 && dt == 0) ||
          (connectivity == 18 && axes == 2) ||
          connectivity == 26;
        if (keep) out.push_back({dx, dy, 0, dt});
      }
  if (has_z) {
    out.push_back({0, 0, -1, 0});
    out.push_back({0, 0, 1, 0});
  }
  return out;
}

} // namespace

int label_components(const int* mask, int* labels,
                     int nx, int ny, int nz, int nt,
                     bool wrap_dateline, int connectivity) {
  const int64_t n = static_cast<int64_t>(nx) * ny * nz * nt;
  const std::vector<Offset> all = neighbour_offsets(connectivity, nz > 1);

  // Pass 1 looks only at neighbours already visited by the raster scan,
  // namely those whose first non-zero step in (dt, dz, dy, dx) is -1.
  std::vector<Offset> back;
  for (const Offset& o : all) {
    int first = o.dt != 0 ? o.dt : o.dz != 0 ? o.dz : o.dy != 0 ? o.dy : o.dx;
    if (first < 0) back.push_back(o);
  }

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
  auto index = [&](int i, int j, int kd, int kt) {
    return kt * stride_t + kd * stride_z + j * stride_y + i;
  };
  auto inside = [&](int j, int kd, int kt) {
    return j >= 0 && j < ny && kd >= 0 && kd < nz && kt >= 0 && kt < nt;
  };
  std::vector<int64_t> back_step;
  for (const Offset& o : back) back_step.push_back(index(o.dx, o.dy, o.dz, o.dt));

  for (int kt = 0; kt < nt; ++kt) {
    for (int kd = 0; kd < nz; ++kd) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          const int64_t v = index(i, j, kd, kt);
          if (mask[v] == 0) continue;

          int chosen = 0;
          // Away from the domain edges every visited neighbour is in bounds
          // (time steps in `back` are never +1).
          const bool interior = i > 0 && i < nx - 1 && j > 0 && j < ny - 1 &&
            (nz == 1 || (kd > 0 && kd < nz - 1)) && kt > 0;
          for (size_t s = 0; s < back.size(); ++s) {
            const Offset& o = back[s];
            if (!interior) {
              int ii = i + o.dx;
              if (ii < 0 || ii >= nx || !inside(j + o.dy, kd + o.dz, kt + o.dt)) continue;
            }
            int l = labels[v + back_step[s]];
            if (l == 0) continue;
            if (chosen == 0) chosen = l;
            else if (l != chosen) uf_union(parent, rank_, chosen, l);
          }

          if (chosen == 0) {
            // start a new component
            parent.push_back(next_label);
            rank_.push_back(0);
            labels[v] = next_label;
            ++next_label;
          } else {
            labels[v] = chosen;
          }
        }
      }
    }
  }

  // Dateline wrap: union the x = 0 column with its neighbours in the
  // x = nx - 1 column (every pair across the boundary, both directions).
  if (wrap_dateline && nx >= 2) {
    for (int kt = 0; kt < nt; ++kt) {
      for (int kd = 0; kd < nz; ++kd) {
        for (int j = 0; j < ny; ++j) {
          const int l0 = labels[index(0, j, kd, kt)];
          if (l0 == 0) continue;
          for (const Offset& o : all) {
            if (o.dx != -1 || !inside(j + o.dy, kd + o.dz, kt + o.dt)) continue;
            int l1 = labels[index(nx - 1, j + o.dy, kd + o.dz, kt + o.dt)];
            if (l1 > 0) uf_union(parent, rank_, l0, l1);
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
