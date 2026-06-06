// Spatial-event ("blob") mask construction and metric reduction for
// heatwave3::detect_blob3().
//
// Two speed-critical steps live here, both single-pass over the grid:
//   hw3_blob_build_mask()  builds the 3-D daily-exceedance mask (column-major,
//                          ready for label_components_3d_cpp()).
//   hw3_blob_reduce()      reduces labelled voxels to the v1 event and daily
//                          tables: areas, volume, cumulative intensity, signed
//                          intensities, peak Hobday severity and category-area
//                          fractions, great-circle area-weighted peak-day
//                          centroid, and bounding box.
//
// Layout conventions (shared with detect_blob3.R). The inputs are the gridded
// daily product from detect_event3(daily = "also"); the mask is the Hobday
// duration-filtered per-pixel event flag, not a raw exceedance.
//   temp/seas/thresh/event   pixel-major   x[px*ntime + t],  px = i*nlat + j
//   mask / labels            column-major  v = i + nlon*j + nlon*nlat*t
//
// Anomaly Delta = temp - seas keeps its sign (negative for cold spells). The
// Hobday severity multiplier s = Delta / (thresh - seas) is >= 1 for every
// event voxel and positive for both polarities.

#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <climits>

using namespace Rcpp;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Reshape the pixel-major Hobday event flag into the column-major mask consumed
// by label_components_3d_cpp(). event[px*ntime + t] > 0 marks an event voxel.
// [[Rcpp::export]]
IntegerVector hw3_blob_mask_from_event(IntegerVector event,
                                       int nlon, int nlat, int ntime) {
  const R_xlen_t ncell = static_cast<R_xlen_t>(nlon) * nlat;
  IntegerVector mask(ncell * ntime);  // zero-initialised, column-major
  const int* ev = INTEGER(event);
  int* maskp = INTEGER(mask);
  for (int i = 0; i < nlon; ++i) {
    for (int j = 0; j < nlat; ++j) {
      const R_xlen_t px = static_cast<R_xlen_t>(i) * nlat + j;
      const int* ep = ev + px * ntime;
      for (int t = 0; t < ntime; ++t) {
        if (ep[t] != NA_INTEGER && ep[t] > 0) {
          const R_xlen_t v = static_cast<R_xlen_t>(i)
                           + static_cast<R_xlen_t>(nlon) * j
                           + static_cast<R_xlen_t>(nlon) * nlat * t;
          maskp[v] = 1;
        }
      }
    }
  }
  return mask;
}

namespace {
struct DailyAcc {
  double area, sum_delta_area;
  double gx, gy, gz;                 // great-circle centroid accumulators
  double lon_min, lon_max, lat_min, lat_max;
  double max_delta;                  // signed, by magnitude
  double cat_area[5];                // index 1..4 = Moderate..Extreme
  int n_cells;
  bool max_set;
  DailyAcc()
    : area(0), sum_delta_area(0), gx(0), gy(0), gz(0),
      lon_min(R_PosInf), lon_max(R_NegInf), lat_min(R_PosInf), lat_max(R_NegInf),
      max_delta(0), n_cells(0), max_set(false) {
    cat_area[0] = cat_area[1] = cat_area[2] = cat_area[3] = cat_area[4] = 0.0;
  }
};
} // namespace

// Reduce labelled voxels to the v1 event and daily tables. labels are the dense
// component labels (1..n_blobs) in column-major order; cell_area_lat[j] is the
// per-latitude cell area (km^2) precomputed on the R side per cellAreaMethod.
// [[Rcpp::export]]
List hw3_blob_reduce(IntegerVector labels, int n_blobs,
                     NumericVector temp, NumericVector seas, NumericVector thresh,
                     NumericVector lon, NumericVector lat,
                     NumericVector cell_area_lat,
                     int nlon, int nlat, int ntime, bool want_voxel) {
  const R_xlen_t ncell = static_cast<R_xlen_t>(nlon) * nlat;
  const R_xlen_t nvox_total = ncell * ntime;
  const int* lab = INTEGER(labels);
  const double* tempp = REAL(temp);
  const double* seasp = REAL(seas);
  const double* thrp = REAL(thresh);
  const double* lonp = REAL(lon);
  const double* latp = REAL(lat);
  const double* areap = REAL(cell_area_lat);
  const double DEG = M_PI / 180.0;
  const int K = n_blobs;

  // Per-blob accumulators (index 1..K).
  std::vector<double> b_cumI(K + 1, 0.0), b_vol(K + 1, 0.0), b_total(K + 1, 0.0);
  std::vector<double> b_maxdelta(K + 1, 0.0), b_maxsev(K + 1, R_NegInf);
  std::vector<int> b_nvox(K + 1, 0), b_tmin(K + 1, INT_MAX), b_tmax(K + 1, INT_MIN);
  std::vector<double> b_lonmin(K + 1, R_PosInf), b_lonmax(K + 1, R_NegInf),
                      b_latmin(K + 1, R_PosInf), b_latmax(K + 1, R_NegInf);
  std::vector<char> b_maxset(K + 1, 0);

  std::unordered_map<int64_t, DailyAcc> dmap;
  std::unordered_set<int64_t> seen_cells;     // (blob, cell) -> total_area dedupe

  std::vector<int> vx_blob, vx_i, vx_j, vx_t;
  std::vector<double> vx_delta;

  for (R_xlen_t v = 0; v < nvox_total; ++v) {
    const int L = lab[v];
    if (L == 0) continue;
    const int i = static_cast<int>(v % nlon);
    const R_xlen_t rem = v / nlon;
    const int j = static_cast<int>(rem % nlat);
    const int t = static_cast<int>(rem / nlat);
    const R_xlen_t px = static_cast<R_xlen_t>(i) * nlat + j;
    const R_xlen_t pt = px * ntime + t;
    const double tv = tempp[pt];
    const double s = seasp[pt];
    const double th = thrp[pt];
    if (ISNAN(tv) || ISNAN(s)) continue;
    const double delta = tv - s;
    const double diff = th - s;
    const double area = areap[j];
    const double lonv = lonp[i];
    const double latv = latp[j];
    const double sev = (!ISNAN(th) && diff != 0.0) ? delta / diff : NA_REAL;

    // Per-blob reduction.
    b_cumI[L] += delta * area;
    b_vol[L] += area;
    ++b_nvox[L];
    if (t < b_tmin[L]) b_tmin[L] = t;
    if (t > b_tmax[L]) b_tmax[L] = t;
    if (!b_maxset[L] || std::fabs(delta) > std::fabs(b_maxdelta[L])) {
      b_maxdelta[L] = delta; b_maxset[L] = 1;
    }
    if (!ISNAN(sev) && sev > b_maxsev[L]) b_maxsev[L] = sev;
    if (lonv < b_lonmin[L]) b_lonmin[L] = lonv;
    if (lonv > b_lonmax[L]) b_lonmax[L] = lonv;
    if (latv < b_latmin[L]) b_latmin[L] = latv;
    if (latv > b_latmax[L]) b_latmax[L] = latv;
    const int64_t ckey = static_cast<int64_t>(L) * ncell + (i + static_cast<int64_t>(nlon) * j);
    if (seen_cells.insert(ckey).second) b_total[L] += area;

    // Per-(blob, day) reduction.
    const int64_t dkey = static_cast<int64_t>(L) * ntime + t;
    DailyAcc& e = dmap[dkey];
    e.area += area;
    e.sum_delta_area += delta * area;
    const double phi = latv * DEG, lam = lonv * DEG, cphi = std::cos(phi);
    e.gx += area * cphi * std::cos(lam);
    e.gy += area * cphi * std::sin(lam);
    e.gz += area * std::sin(phi);
    if (lonv < e.lon_min) e.lon_min = lonv;
    if (lonv > e.lon_max) e.lon_max = lonv;
    if (latv < e.lat_min) e.lat_min = latv;
    if (latv > e.lat_max) e.lat_max = latv;
    if (!e.max_set || std::fabs(delta) > std::fabs(e.max_delta)) {
      e.max_delta = delta; e.max_set = true;
    }
    ++e.n_cells;
    if (!ISNAN(sev)) {
      int c = static_cast<int>(std::floor(sev));
      if (c < 1) c = 1; if (c > 4) c = 4;
      e.cat_area[c] += area;
    }

    if (want_voxel) {
      vx_blob.push_back(L); vx_i.push_back(i); vx_j.push_back(j);
      vx_t.push_back(t); vx_delta.push_back(delta);
    }
  }

  // Daily table, ordered by (blob, t); also derive per-blob peak day.
  std::vector<int64_t> dkeys;
  dkeys.reserve(dmap.size());
  for (auto& kv : dmap) dkeys.push_back(kv.first);
  std::sort(dkeys.begin(), dkeys.end());
  const R_xlen_t M = dkeys.size();

  IntegerVector d_blob(M), d_t(M), d_n(M);
  NumericVector d_area(M), d_mi(M), d_md(M), d_clon(M), d_clat(M),
                d_lonmin(M), d_lonmax(M), d_latmin(M), d_latmax(M);

  std::vector<int> b_ndays(K + 1, 0), pk_t(K + 1, -1);
  std::vector<double> pk_area(K + 1, R_NegInf), pk_clon(K + 1, NA_REAL),
                      pk_clat(K + 1, NA_REAL), pk_cat(4 * (K + 1), 0.0);

  for (R_xlen_t m = 0; m < M; ++m) {
    const DailyAcc& e = dmap[dkeys[m]];
    const int L = static_cast<int>(dkeys[m] / ntime);
    const int t = static_cast<int>(dkeys[m] % ntime);
    double clon = NA_REAL, clat = NA_REAL;
    const double horiz = std::sqrt(e.gx * e.gx + e.gy * e.gy);
    if (horiz > 0 || e.gz != 0) {
      clon = std::atan2(e.gy, e.gx) / DEG;
      clat = std::atan2(e.gz, horiz) / DEG;
    }
    d_blob[m] = L; d_t[m] = t + 1;
    d_area[m] = e.area;
    d_mi[m] = e.area > 0 ? e.sum_delta_area / e.area : NA_REAL;
    d_md[m] = e.max_set ? e.max_delta : NA_REAL;
    d_clon[m] = clon; d_clat[m] = clat;
    d_lonmin[m] = e.lon_min; d_lonmax[m] = e.lon_max;
    d_latmin[m] = e.lat_min; d_latmax[m] = e.lat_max;
    d_n[m] = e.n_cells;

    ++b_ndays[L];
    if (e.area > pk_area[L]) {
      pk_area[L] = e.area; pk_t[L] = t; pk_clon[L] = clon; pk_clat[L] = clat;
      for (int c = 0; c < 4; ++c) pk_cat[L * 4 + c] = e.cat_area[c + 1];
    }
  }

  // Event table.
  IntegerVector e_no(K), e_ts(K), e_te(K), e_tp(K), e_dur(K), e_nvox(K);
  NumericVector e_parea(K), e_marea(K), e_tarea(K), e_vol(K), e_cumI(K),
                e_mi(K), e_maxI(K), e_psev(K),
                e_fmod(K), e_fstr(K), e_fsev(K), e_fext(K),
                e_clon(K), e_clat(K), e_lonmin(K), e_lonmax(K), e_latmin(K), e_latmax(K);
  for (int L = 1; L <= K; ++L) {
    const int idx = L - 1;
    e_no[idx] = L;
    e_ts[idx] = b_tmin[L] + 1;
    e_te[idx] = b_tmax[L] + 1;
    e_tp[idx] = pk_t[L] + 1;
    e_dur[idx] = b_tmax[L] - b_tmin[L] + 1;
    e_nvox[idx] = b_nvox[L];
    e_parea[idx] = pk_area[L];
    e_marea[idx] = b_ndays[L] > 0 ? b_vol[L] / b_ndays[L] : NA_REAL;
    e_tarea[idx] = b_total[L];
    e_vol[idx] = b_vol[L];
    e_cumI[idx] = b_cumI[L];
    e_mi[idx] = b_vol[L] != 0.0 ? b_cumI[L] / b_vol[L] : NA_REAL;
    e_maxI[idx] = b_maxset[L] ? b_maxdelta[L] : NA_REAL;
    e_psev[idx] = (b_maxsev[L] == R_NegInf) ? NA_REAL : b_maxsev[L];
    const double pa = pk_area[L];
    e_fmod[idx] = pa > 0 ? pk_cat[L * 4 + 0] / pa : NA_REAL;
    e_fstr[idx] = pa > 0 ? pk_cat[L * 4 + 1] / pa : NA_REAL;
    e_fsev[idx] = pa > 0 ? pk_cat[L * 4 + 2] / pa : NA_REAL;
    e_fext[idx] = pa > 0 ? pk_cat[L * 4 + 3] / pa : NA_REAL;
    e_clon[idx] = pk_clon[L]; e_clat[idx] = pk_clat[L];
    e_lonmin[idx] = b_lonmin[L]; e_lonmax[idx] = b_lonmax[L];
    e_latmin[idx] = b_latmin[L]; e_latmax[idx] = b_latmax[L];
  }

  const char* enm[] = {"event_no", "t_start", "t_end", "t_peak", "duration_days",
    "n_voxels", "peak_area_km2", "mean_area_km2", "total_area_km2", "volume_km2_d",
    "cumI_km2_d", "mean_intensity", "max_intensity", "peak_severity",
    "frac_moderate", "frac_strong", "frac_severe", "frac_extreme",
    "centroid_lon", "centroid_lat", "bbox_lon_min", "bbox_lon_max",
    "bbox_lat_min", "bbox_lat_max"};
  List event(24);
  event[0] = e_no; event[1] = e_ts; event[2] = e_te; event[3] = e_tp;
  event[4] = e_dur; event[5] = e_nvox; event[6] = e_parea; event[7] = e_marea;
  event[8] = e_tarea; event[9] = e_vol; event[10] = e_cumI; event[11] = e_mi;
  event[12] = e_maxI; event[13] = e_psev; event[14] = e_fmod; event[15] = e_fstr;
  event[16] = e_fsev; event[17] = e_fext; event[18] = e_clon; event[19] = e_clat;
  event[20] = e_lonmin; event[21] = e_lonmax; event[22] = e_latmin; event[23] = e_latmax;
  event.attr("names") = CharacterVector(enm, enm + 24);

  const char* dnm[] = {"event_no", "t_idx", "area_km2", "mean_intensity",
    "max_intensity", "centroid_lon", "centroid_lat", "bbox_lon_min",
    "bbox_lon_max", "bbox_lat_min", "bbox_lat_max", "n_cells"};
  List daily(12);
  daily[0] = d_blob; daily[1] = d_t; daily[2] = d_area; daily[3] = d_mi;
  daily[4] = d_md; daily[5] = d_clon; daily[6] = d_clat; daily[7] = d_lonmin;
  daily[8] = d_lonmax; daily[9] = d_latmin; daily[10] = d_latmax; daily[11] = d_n;
  daily.attr("names") = CharacterVector(dnm, dnm + 12);

  if (want_voxel) {
    List voxel = List::create(
      _["event_no"] = wrap(vx_blob), _["i"] = wrap(vx_i), _["j"] = wrap(vx_j),
      _["t_idx"] = wrap(vx_t), _["delta"] = wrap(vx_delta));
    return List::create(_["event"] = event, _["daily"] = daily, _["voxel"] = voxel);
  }
  return List::create(_["event"] = event, _["daily"] = daily);
}
