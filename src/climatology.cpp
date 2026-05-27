#include "climatology.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cstring>
#include <atomic>
#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace hw3 {

// ---- Date utilities ----

int ymd_to_jd(int y, int m, int d) {
    int a = (14 - m) / 12;
    int yy = y + 4800 - a;
    int mm = m + 12 * a - 3;
    return d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045;
}

void jd_to_ymd(int jd, int& y, int& m, int& d) {
    int a = jd + 32044;
    int b = (4 * a + 3) / 146097;
    int c = a - (146097 * b) / 4;
    int dd = (4 * c + 3) / 1461;
    int e = c - (1461 * dd) / 4;
    int mm = (5 * e + 2) / 153;
    d = e - (153 * mm + 2) / 5 + 1;
    m = mm + 3 - 12 * (mm / 10);
    y = 100 * b + dd - 4800 + mm / 10;
}

static bool is_leap_year(int y) {
    return (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0);
}

bool is_leap_year_from_jd(int jd) {
    int y, m, d;
    jd_to_ymd(jd, y, m, d);
    return is_leap_year(y);
}

// Returns DOY in 1..366 range, matching heatwaveR convention:
// - Leap years: DOY 60 = Feb 29
// - Non-leap years: DOYs skip 60 (go from 59 to 61)
int jd_to_doy_366(int jd) {
    int y, m, d;
    jd_to_ymd(jd, y, m, d);

    // Standard DOY (1-based)
    int jan1_jd = ymd_to_jd(y, 1, 1);
    int doy = jd - jan1_jd + 1;

    // Non-leap years: shift DOY >= 60 up by 1 to skip DOY 60
    if (!is_leap_year(y) && doy >= 60) {
        doy += 1;
    }
    return doy;
}

// ---- NA interpolation (linear, up to maxPadLength consecutive NAs) ----

static void interpolate_na(double* vals, int n, int maxPadLength) {
    if (maxPadLength <= 0) return;

    int i = 0;
    while (i < n) {
        if (is_na(vals[i])) {
            int gap_start = i;
            while (i < n && is_na(vals[i])) ++i;
            int gap_len = i - gap_start;

            if (gap_len <= maxPadLength) {
                double v_before = (gap_start > 0) ? vals[gap_start - 1] : NA_DOUBLE;
                double v_after = (i < n) ? vals[i] : NA_DOUBLE;

                if (!is_na(v_before) && !is_na(v_after)) {
                    for (int j = 0; j < gap_len; ++j) {
                        double frac = static_cast<double>(j + 1) / (gap_len + 1);
                        vals[gap_start + j] = v_before + frac * (v_after - v_before);
                    }
                } else if (!is_na(v_before)) {
                    for (int j = 0; j < gap_len; ++j) vals[gap_start + j] = v_before;
                } else if (!is_na(v_after)) {
                    for (int j = 0; j < gap_len; ++j) vals[gap_start + j] = v_after;
                }
            }
        } else {
            ++i;
        }
    }
}

// ---- Type-7 quantile (R's default) ----

static double quantile_type7(const double* sorted, int n, double p) {
    if (n == 0) return NA_DOUBLE;
    if (n == 1) return sorted[0];

    double index = (n - 1) * p;
    int lo = static_cast<int>(std::floor(index));
    int hi = lo + 1;
    if (hi >= n) hi = n - 1;
    double hhi = index - lo;
    double hlo = 1.0 - hhi;
    return hlo * sorted[lo] + hhi * sorted[hi];
}

// ---- Rolling mean with circular padding ----

static void smooth_circular(double* vals, int n, int window) {
    if (window <= 1 || n == 0) return;
    int half = window / 2;

    // Create padded array: [last half values] [values] [first half values]
    int padded_len = n + 2 * half;
    std::vector<double> padded(padded_len);
    for (int i = 0; i < half; ++i) padded[i] = vals[n - half + i];
    for (int i = 0; i < n; ++i) padded[half + i] = vals[i];
    for (int i = 0; i < half; ++i) padded[half + n + i] = vals[i];

    // Compute rolling mean (centered, width = window)
    // The RcppRoll::roll_mean with n=31 produces output of length = input - n + 1
    // After padding by half on each side, the valid output aligns correctly
    std::vector<double> result(n);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < window; ++j) {
            sum += padded[i + j];
        }
        result[i] = sum / window;
    }
    std::memcpy(vals, result.data(), n * sizeof(double));
}

// ---- Main per-pixel climatology ----

PixelClimResult compute_pixel_climatology(
    const double* temp,
    const int* time_jd,
    int ntime,
    int clim_start_jd,
    int clim_end_jd,
    int windowHalfWidth,
    double pctile,
    bool smoothPercentile,
    int smoothPercentileWidth,
    bool compute_var,
    int maxPadLength,
    bool detrend) {

    PixelClimResult result;
    result.valid = false;
    for (int i = 0; i < 366; ++i) {
        result.seas[i] = NA_DOUBLE;
        result.thresh[i] = NA_DOUBLE;
        result.var[i] = NA_DOUBLE;
    }

    // Step 1: Build continuous daily series with DOY assignment
    // Find date range
    int first_jd = time_jd[0];
    int last_jd = time_jd[ntime - 1];
    int total_days = last_jd - first_jd + 1;

    // Create full daily array (fill gaps with NA)
    std::vector<double> daily(total_days, NA_DOUBLE);
    std::vector<int> daily_jd(total_days);
    std::vector<int> daily_doy(total_days);

    for (int i = 0; i < total_days; ++i) {
        daily_jd[i] = first_jd + i;
        daily_doy[i] = jd_to_doy_366(daily_jd[i]);
    }

    // Fill in observed values
    for (int i = 0; i < ntime; ++i) {
        int idx = time_jd[i] - first_jd;
        if (idx >= 0 && idx < total_days) {
            daily[idx] = temp[i];
        }
    }

    // Step 2: Interpolate NAs if requested
    if (maxPadLength > 0) {
        interpolate_na(daily.data(), total_days, maxPadLength);
    }

    // Step 2b: Linear detrend within climatology period (Jacox et al. 2020)
    if (detrend) {
        int cs_idx = std::max(0, clim_start_jd - first_jd);
        int ce_idx = std::min(total_days - 1, clim_end_jd - first_jd);

        // Least-squares fit: temp = a + b * t
        double sum_t = 0, sum_y = 0, sum_tt = 0, sum_ty = 0;
        int n_valid = 0;
        for (int i = cs_idx; i <= ce_idx; ++i) {
            if (is_na(daily[i])) continue;
            double t = static_cast<double>(i - cs_idx);
            sum_t += t;
            sum_y += daily[i];
            sum_tt += t * t;
            sum_ty += t * daily[i];
            ++n_valid;
        }

        if (n_valid > 2) {
            double mean_t = sum_t / n_valid;
            double mean_y = sum_y / n_valid;
            double slope = (sum_ty - n_valid * mean_t * mean_y) /
                           (sum_tt - n_valid * mean_t * mean_t);

            // Subtract trend from ALL days (not just clim period)
            // Keep intercept at the midpoint of the climatology period
            double mid_t = static_cast<double>(ce_idx - cs_idx) / 2.0;
            for (int i = 0; i < total_days; ++i) {
                if (!is_na(daily[i])) {
                    double t = static_cast<double>(i - cs_idx);
                    daily[i] -= slope * (t - mid_t);
                }
            }
        }
    }

    // Step 3: Extract climatology period
    int clim_start_idx = clim_start_jd - first_jd;
    int clim_end_idx = clim_end_jd - first_jd;
    if (clim_start_idx < 0) clim_start_idx = 0;
    if (clim_end_idx >= total_days) clim_end_idx = total_days - 1;

    // Determine years in climatology period
    int y1, m1, d1, y2, m2, d2;
    jd_to_ymd(clim_start_jd, y1, m1, d1);
    jd_to_ymd(clim_end_jd, y2, m2, d2);
    int n_years = y2 - y1 + 1;

    if (n_years < 3) return result;

    // Step 4: Spread into DOY × year matrix
    // Matrix: 366 rows (DOY 1-366) × n_years columns
    std::vector<double> spread(366 * n_years, NA_DOUBLE);

    for (int i = clim_start_idx; i <= clim_end_idx; ++i) {
        int y, m, d;
        jd_to_ymd(daily_jd[i], y, m, d);
        int col = y - y1;
        int row = daily_doy[i] - 1; // 0-based
        if (row >= 0 && row < 366 && col >= 0 && col < n_years) {
            spread[row * n_years + col] = daily[i];
        }
    }

    // Step 5: Fill Feb 29 (DOY 60, index 59) for non-leap years
    // Interpolate from DOY 59 and DOY 61
    for (int col = 0; col < n_years; ++col) {
        if (is_na(spread[59 * n_years + col])) {
            double v59 = spread[58 * n_years + col]; // DOY 59
            double v61 = spread[60 * n_years + col]; // DOY 61
            if (!is_na(v59) && !is_na(v61)) {
                spread[59 * n_years + col] = (v59 + v61) / 2.0;
            } else if (!is_na(v59)) {
                spread[59 * n_years + col] = v59;
            } else if (!is_na(v61)) {
                spread[59 * n_years + col] = v61;
            }
        }
    }

    // Step 6: Add circular padding (windowHalfWidth rows at each end)
    int padded_rows = 366 + 2 * windowHalfWidth;
    std::vector<double> padded(padded_rows * n_years);

    // Copy tail padding (last windowHalfWidth DOYs)
    for (int r = 0; r < windowHalfWidth; ++r) {
        int src_row = 366 - windowHalfWidth + r;
        std::memcpy(&padded[r * n_years],
                     &spread[src_row * n_years],
                     n_years * sizeof(double));
    }
    // Copy main data
    std::memcpy(&padded[windowHalfWidth * n_years],
                 spread.data(),
                 366 * n_years * sizeof(double));
    // Copy head padding (first windowHalfWidth DOYs)
    for (int r = 0; r < windowHalfWidth; ++r) {
        std::memcpy(&padded[(windowHalfWidth + 366 + r) * n_years],
                     &spread[r * n_years],
                     n_years * sizeof(double));
    }

    // Step 7: Compute climatology using sliding window
    double p = pctile / 100.0;
    int window_rows = 2 * windowHalfWidth + 1;
    int max_window_elements = window_rows * n_years;
    std::vector<double> window_vals(max_window_elements);

    for (int doy = 0; doy < 366; ++doy) {
        int center = doy + windowHalfWidth;

        // Extract window: rows [center - windowHalfWidth, center + windowHalfWidth]
        int n_finite = 0;
        for (int r = center - windowHalfWidth; r <= center + windowHalfWidth; ++r) {
            for (int c = 0; c < n_years; ++c) {
                double v = padded[r * n_years + c];
                if (!is_na(v)) {
                    window_vals[n_finite++] = v;
                }
            }
        }

        if (n_finite == 0) continue;

        // Sort for quantile and mean
        std::sort(window_vals.begin(), window_vals.begin() + n_finite);

        // Mean (seasonal climatology)
        double sum = 0.0;
        for (int i = 0; i < n_finite; ++i) sum += window_vals[i];
        result.seas[doy] = sum / n_finite;

        // Threshold (percentile)
        result.thresh[doy] = quantile_type7(window_vals.data(), n_finite, p);

        // Variance (optional)
        if (compute_var) {
            double mean = result.seas[doy];
            double var_sum = 0.0;
            for (int i = 0; i < n_finite; ++i) {
                double diff = window_vals[i] - mean;
                var_sum += diff * diff;
            }
            result.var[doy] = (n_finite > 1) ?
                std::sqrt(var_sum / (n_finite - 1)) : 0.0;
        }
    }

    // Step 8: Smooth with rolling mean (circular padding)
    if (smoothPercentile) {
        smooth_circular(result.seas, 366, smoothPercentileWidth);
        smooth_circular(result.thresh, 366, smoothPercentileWidth);
        if (compute_var) {
            smooth_circular(result.var, 366, smoothPercentileWidth);
        }
    }

    // Check if we got any valid values
    for (int i = 0; i < 366; ++i) {
        if (!is_na(result.seas[i])) {
            result.valid = true;
            break;
        }
    }

    return result;
}

// ---- Grid-level computation with OpenMP ----

void compute_climatology_grid(
    const double* sst,
    const int* time_jd,
    int npixels,
    int ntime,
    int clim_start_jd,
    int clim_end_jd,
    int windowHalfWidth,
    double pctile,
    bool smoothPercentile,
    int smoothPercentileWidth,
    bool compute_var,
    int maxPadLength,
    bool detrend,
    int n_threads,
    double* seas_out,
    double* thresh_out,
    double* var_out) {

    // Initialize output to NA
    size_t grid_size = static_cast<size_t>(npixels) * 366;
    for (size_t i = 0; i < grid_size; ++i) {
        seas_out[i] = NA_DOUBLE;
        thresh_out[i] = NA_DOUBLE;
        if (var_out) var_out[i] = NA_DOUBLE;
    }

#ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
#endif

    std::atomic<int> done_pixels{0};
    int report_interval = std::max(1, npixels / 20);

    #pragma omp parallel for schedule(dynamic)
    for (int px = 0; px < npixels; ++px) {
        const double* pixel_temp = sst + static_cast<size_t>(px) * ntime;

        // Skip all-NA pixels
        bool has_data = false;
        for (int t = 0; t < ntime; ++t) {
            if (!is_na(pixel_temp[t])) { has_data = true; break; }
        }
        if (!has_data) {
            int cur = ++done_pixels;
            if (cur % report_interval == 0 || cur == npixels) {
                #pragma omp critical
                {
                    Rcpp::Rcout << "\r  " << cur << "/" << npixels << " pixels ("
                                << (100 * cur / npixels) << "%)" << std::flush;
                }
            }
            continue;
        }

        PixelClimResult cr = compute_pixel_climatology(
            pixel_temp, time_jd, ntime,
            clim_start_jd, clim_end_jd,
            windowHalfWidth, pctile,
            smoothPercentile, smoothPercentileWidth,
            compute_var, maxPadLength, detrend
        );

        if (cr.valid) {
            size_t offset = static_cast<size_t>(px) * 366;
            std::memcpy(seas_out + offset, cr.seas, 366 * sizeof(double));
            std::memcpy(thresh_out + offset, cr.thresh, 366 * sizeof(double));
            if (var_out && compute_var) {
                std::memcpy(var_out + offset, cr.var, 366 * sizeof(double));
            }
        }

        int cur = ++done_pixels;
        if (cur % report_interval == 0 || cur == npixels) {
            #pragma omp critical
            {
                Rcpp::Rcout << "\r  " << cur << "/" << npixels << " pixels ("
                            << (100 * cur / npixels) << "%)" << std::flush;
            }
        }
    }

    Rcpp::Rcout << std::endl;
}

} // namespace hw3
