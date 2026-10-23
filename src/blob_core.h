#ifndef HW3_BLOB_CORE_H
#define HW3_BLOB_CORE_H

// Language-neutral blob routines behind label_components_cpp() and
// blob_daily_summary_cpp() (blob_label.cpp / blob_metrics.cpp hold the R
// wrappers). Input lengths are validated by the callers.

#include <cstdint>
#include <vector>

namespace hw3 {

// 3D/4D connected-component labelling. `mask` and `labels` both hold
// nx*ny*nz*nt values in column-major [i, j, kd, kt] order; `labels` must be
// zero-initialised. `connectivity` is 6, 10, 18 or 26, counted in the
// lon/lat/time cube; depth always joins only straight up and down (see
// blob_core.cpp). Other values throw. Writes dense labels 1..K and returns K.
int label_components(const int* mask, int* labels,
                     int nx, int ny, int nz, int nt,
                     bool wrap_dateline, int connectivity);

// Per-(blob, day) metrics, ordered by (blob, t_idx). The depth vectors are
// filled only when has_depth is true.
struct BlobDailySummary {
    std::vector<int> blob, t_idx, n_cells;
    std::vector<double> area_km2, mean_delta, max_delta;
    std::vector<double> centroid_lon, centroid_lat;
    std::vector<double> lon_min, lon_max, lat_min, lat_max;
    std::vector<double> volume_km3, depth_min_m, depth_max_m, depth_at_max_delta_m;
};

// Voxel inputs are parallel arrays of length n with 1-based ix/iy/it/iz.
// depth_at_max_delta_m is NaN where no depth was recorded.
BlobDailySummary blob_daily_summary(const int* labels, const int* ix,
                                    const int* iy, const int* it,
                                    const double* delta, int64_t n,
                                    const double* cell_area_lat,
                                    const double* lons, const double* lats,
                                    const int* iz, const double* depths,
                                    const double* layer_thickness_km,
                                    bool has_depth);

} // namespace hw3

#endif
