#include "event_detect.h"
#include "climatology.h"
#include "hw3_omp.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cstring>
#include <mutex>
#include <atomic>
#include <Rcpp.h>

namespace hw3 {

// ---- Proto-event detection (RLE + duration filter + gap joining) ----

struct ProtoEvent {
    int start;
    int end;
    int event_no;
};

static std::vector<ProtoEvent> proto_event(
    const std::vector<bool>& thresh_criterion,
    int n,
    int minDuration,
    bool joinAcrossGaps,
    int maxGap) {

    // Phase 1: Run-length encoding — find consecutive TRUE runs
    std::vector<ProtoEvent> runs;
    int i = 0;
    while (i < n) {
        if (thresh_criterion[i]) {
            int start = i;
            while (i < n && thresh_criterion[i]) ++i;
            int len = i - start;
            if (len >= minDuration) {
                ProtoEvent pe;
                pe.start = start;
                pe.end = i - 1;
                pe.event_no = 0;
                runs.push_back(pe);
            }
        } else {
            ++i;
        }
    }

    if (runs.empty()) return runs;

    // Phase 2: Join across gaps
    if (joinAcrossGaps && maxGap > 0 && runs.size() > 1) {
        std::vector<ProtoEvent> merged;
        merged.push_back(runs[0]);
        for (size_t r = 1; r < runs.size(); ++r) {
            int gap = runs[r].start - merged.back().end - 1;
            if (gap <= maxGap) {
                merged.back().end = runs[r].end;
            } else {
                merged.push_back(runs[r]);
            }
        }
        runs = merged;
    }

    // Phase 3: Re-check duration after merging
    std::vector<ProtoEvent> final_events;
    int eno = 1;
    for (auto& pe : runs) {
        int dur = pe.end - pe.start + 1;
        if (dur >= minDuration) {
            pe.event_no = eno++;
            final_events.push_back(pe);
        }
    }

    return final_events;
}

// ---- Compute event metrics ----

static double round_to(double val, int decimals) {
    if (is_na(val)) return val;
    double factor = std::pow(10.0, decimals);
    return std::round(val * factor) / factor;
}

PixelEvents detect_pixel_events(
    const double* temp,
    const double* seas,
    const double* thresh,
    const int* time_jd,
    const int* doy,
    int ntime,
    int minDuration,
    int minDuration2,
    bool joinAcrossGaps,
    int maxGap,
    int maxGap2,
    bool coldSpells,
    int roundRes,
    bool category,
    bool southHemisphere) {

    PixelEvents result;
    result.valid = false;

    if (ntime <= 0) return result;

    // Step 0: Expand the (possibly sparse) input series to a dense daily
    // array, filling missing dates with NA. This mirrors heatwaveR's
    // make_whole() and ensures that proto_event(), which walks indices,
    // does not glue runs together across calendar-day gaps in the input.
    int first_jd = time_jd[0];
    int last_jd = time_jd[ntime - 1];
    int dntime = last_jd - first_jd + 1;

    std::vector<double> dense_temp(dntime, NA_DOUBLE);
    std::vector<int>    dense_jd(dntime);
    std::vector<int>    dense_doy(dntime);
    for (int i = 0; i < dntime; ++i) {
        dense_jd[i] = first_jd + i;
        dense_doy[i] = jd_to_doy_366(dense_jd[i]);
    }
    for (int t = 0; t < ntime; ++t) {
        int idx = time_jd[t] - first_jd;
        if (idx >= 0 && idx < dntime) {
            dense_temp[idx] = temp[t];
        }
    }
    // From here on we work entirely on the dense arrays. The original
    // sparse `temp`, `time_jd`, `doy` and `ntime` are NOT used below;
    // `index_*` fields in EventResult therefore index the dense array.
    (void)doy;

    // Build per-day seas and thresh from 366-day climatology
    std::vector<double> daily_seas(dntime), daily_thresh(dntime);
    for (int t = 0; t < dntime; ++t) {
        int d = dense_doy[t] - 1; // 0-based DOY index
        if (d >= 0 && d < 366) {
            daily_seas[t] = seas[d];
            daily_thresh[t] = thresh[d];
        } else {
            daily_seas[t] = NA_DOUBLE;
            daily_thresh[t] = NA_DOUBLE;
        }
    }

    // Step 1: Threshold criterion
    // Using char instead of bool to allow .data() access
    std::vector<char> thresh_raw(dntime, 0);
    std::vector<double> temp_filled(dntime);

    for (int t = 0; t < dntime; ++t) {
        double tv = dense_temp[t];
        if (is_na(tv)) {
            // Missing values set to seas (breaks event continuity)
            temp_filled[t] = daily_seas[t];
            thresh_raw[t] = 0;
        } else if (is_na(daily_thresh[t])) {
            temp_filled[t] = tv;
            thresh_raw[t] = 0;
        } else {
            temp_filled[t] = tv;
            if (coldSpells) {
                thresh_raw[t] = (tv < daily_thresh[t]) ? 1 : 0;
            } else {
                thresh_raw[t] = (tv > daily_thresh[t]) ? 1 : 0;
            }
        }
    }

    // Convert to bool array for proto_event
    std::vector<bool> thresh_criterion(dntime);
    for (int t = 0; t < dntime; ++t) thresh_criterion[t] = (thresh_raw[t] != 0);

    // Step 2: Detect proto-events — use manual loop since vector<bool> lacks .data()
    auto events = proto_event(thresh_criterion, dntime,
                              minDuration, joinAcrossGaps, maxGap);

    if (events.empty()) {
        result.valid = true;
        return result;
    }

    // Step 3: Compute metrics for each event
    for (auto& pe : events) {
        EventResult er;
        er.event_no = pe.event_no;
        er.category = 0;
        er.p_moderate = er.p_strong = er.p_severe = er.p_extreme = 0.0;
        er.season = 0;
        er.index_start = pe.start;
        er.index_end = pe.end;
        er.jd_start = dense_jd[pe.start];
        er.jd_end   = dense_jd[pe.end];
        er.duration = pe.end - pe.start + 1;

        // Intensity relative to seas
        double sum_rel_seas = 0.0, max_rel_seas = -1e30;
        double sum_sq_rel_seas = 0.0;
        int peak_idx = pe.start;

        // Intensity relative to thresh
        double sum_rel_thresh = 0.0, max_rel_thresh = -1e30;
        double sum_sq_rel_thresh = 0.0;

        // Absolute intensity
        double sum_abs = 0.0, max_abs = -1e30;
        double sum_sq_abs = 0.0;

        int n_valid = 0;

        for (int t = pe.start; t <= pe.end; ++t) {
            double rel_seas = temp_filled[t] - daily_seas[t];
            double rel_thresh = temp_filled[t] - daily_thresh[t];
            double abs_val = temp_filled[t];

            if (coldSpells) {
                rel_seas = -rel_seas;
                rel_thresh = -rel_thresh;
            }

            sum_rel_seas += rel_seas;
            sum_sq_rel_seas += rel_seas * rel_seas;
            if (rel_seas > max_rel_seas) {
                max_rel_seas = rel_seas;
                peak_idx = t;
            }

            sum_rel_thresh += rel_thresh;
            sum_sq_rel_thresh += rel_thresh * rel_thresh;
            if (rel_thresh > max_rel_thresh) max_rel_thresh = rel_thresh;

            sum_abs += abs_val;
            sum_sq_abs += abs_val * abs_val;
            if (coldSpells) {
                if (abs_val < max_abs || max_abs < -1e29) max_abs = abs_val;
            } else {
                if (abs_val > max_abs) max_abs = abs_val;
            }

            ++n_valid;
        }

        er.index_peak = peak_idx;
        er.jd_peak    = dense_jd[peak_idx];

        if (n_valid > 0) {
            er.intensity_mean = sum_rel_seas / n_valid;
            er.intensity_max = max_rel_seas;
            er.intensity_cumulative = sum_rel_seas;

            double var_seas = (n_valid > 1) ?
                (sum_sq_rel_seas - sum_rel_seas * sum_rel_seas / n_valid) / (n_valid - 1) : 0.0;
            er.intensity_var = std::sqrt(std::max(0.0, var_seas));

            er.intensity_mean_relThresh = sum_rel_thresh / n_valid;
            er.intensity_max_relThresh = max_rel_thresh;
            er.intensity_cumulative_relThresh = sum_rel_thresh;
            double var_thresh = (n_valid > 1) ?
                (sum_sq_rel_thresh - sum_rel_thresh * sum_rel_thresh / n_valid) / (n_valid - 1) : 0.0;
            er.intensity_var_relThresh = std::sqrt(std::max(0.0, var_thresh));

            er.intensity_mean_abs = sum_abs / n_valid;
            er.intensity_max_abs = max_abs;
            er.intensity_cumulative_abs = sum_abs;
            double var_abs = (n_valid > 1) ?
                (sum_sq_abs - sum_abs * sum_abs / n_valid) / (n_valid - 1) : 0.0;
            er.intensity_var_abs = std::sqrt(std::max(0.0, var_abs));

            // Rate of onset and decline (half-day interpolation)
            double start_rel = temp_filled[pe.start] - daily_seas[pe.start];
            double end_rel = temp_filled[pe.end] - daily_seas[pe.end];
            if (coldSpells) {
                start_rel = -start_rel;
                end_rel = -end_rel;
            }

            int days_to_peak = peak_idx - pe.start;
            int days_from_peak = pe.end - peak_idx;
            er.rate_onset = (days_to_peak > 0) ?
                (max_rel_seas - start_rel) / (days_to_peak + 0.5) : 0.0;
            er.rate_decline = (days_from_peak > 0) ?
                (max_rel_seas - end_rel) / (days_from_peak + 0.5) : 0.0;
        }

        // Round
        if (roundRes > 0) {
            er.intensity_mean = round_to(er.intensity_mean, roundRes);
            er.intensity_max = round_to(er.intensity_max, roundRes);
            er.intensity_var = round_to(er.intensity_var, roundRes);
            er.intensity_cumulative = round_to(er.intensity_cumulative, roundRes);
            er.intensity_mean_relThresh = round_to(er.intensity_mean_relThresh, roundRes);
            er.intensity_max_relThresh = round_to(er.intensity_max_relThresh, roundRes);
            er.intensity_var_relThresh = round_to(er.intensity_var_relThresh, roundRes);
            er.intensity_cumulative_relThresh = round_to(er.intensity_cumulative_relThresh, roundRes);
            er.intensity_mean_abs = round_to(er.intensity_mean_abs, roundRes);
            er.intensity_max_abs = round_to(er.intensity_max_abs, roundRes);
            er.intensity_var_abs = round_to(er.intensity_var_abs, roundRes);
            er.intensity_cumulative_abs = round_to(er.intensity_cumulative_abs, roundRes);
            er.rate_onset = round_to(er.rate_onset, roundRes);
            er.rate_decline = round_to(er.rate_decline, roundRes);
        }

        // Category computation (inline, zero extra I/O)
        if (category) {
            int peak_doy_idx = dense_doy[er.index_peak] - 1;
            double s = seas[peak_doy_idx];
            double th = thresh[peak_doy_idx];
            double diff = std::abs(s - th);

            if (!is_na(s) && !is_na(th) && diff > 0.0) {
                double ix = er.intensity_max;
                if (ix >= 4.0 * diff) er.category = 4;
                else if (ix >= 3.0 * diff) er.category = 3;
                else if (ix >= 2.0 * diff) er.category = 2;
                else er.category = 1;

                double round_mult = std::pow(10.0, roundRes > 0 ? roundRes : 4);
                auto rnd = [&](double v) { return std::round(v * round_mult) / round_mult; };
                er.p_moderate = rnd(std::min(1.0, std::max(0.0, ix / diff)));
                er.p_strong   = rnd(std::min(1.0, std::max(0.0, (ix - diff) / diff)));
                er.p_severe   = rnd(std::min(1.0, std::max(0.0, (ix - 2.0 * diff) / diff)));
                er.p_extreme  = rnd(std::min(1.0, std::max(0.0, (ix - 3.0 * diff) / diff)));
            }

            // Season from peak date (use dense JD; the original `time_jd[]`
            // is the caller's sparse array and would index the wrong day
            // when the input series had gaps).
            int peak_jd = er.jd_peak;
            int y, m, d;
            jd_to_ymd(peak_jd, y, m, d);
            if (southHemisphere) {
                if (m == 12 || m <= 2) er.season = 1; // Summer
                else if (m <= 5) er.season = 2; // Fall
                else if (m <= 8) er.season = 3; // Winter
                else er.season = 4; // Spring
            } else {
                if (m == 12 || m <= 2) er.season = 3; // Winter
                else if (m <= 5) er.season = 4; // Spring
                else if (m <= 8) er.season = 1; // Summer
                else er.season = 2; // Fall
            }
        }

        result.events.push_back(er);
    }

    result.valid = true;
    return result;
}

// ---- Grid-level event detection with OpenMP ----

void detect_events_grid(
    const double* sst,
    const double* seas_clim,
    const double* thresh_clim,
    const int* time_jd,
    int npixels,
    int ntime,
    int minDuration,
    int minDuration2,
    bool joinAcrossGaps,
    int maxGap,
    int maxGap2,
    bool coldSpells,
    int roundRes,
    int n_threads,
    bool category,
    bool southHemisphere,
    std::vector<double>& event_lon,
    std::vector<double>& event_lat,
    std::vector<int>& pixel_index,
    std::vector<EventResult>& all_events,
    std::vector<int>& date_start,
    std::vector<int>& date_peak,
    std::vector<int>& date_end,
    const std::vector<double>& lon,
    const std::vector<double>& lat,
    int nlon, int nlat) {

    // Compute DOY for each time step
    std::vector<int> doy_arr(ntime);
    for (int t = 0; t < ntime; ++t) {
        doy_arr[t] = jd_to_doy_366(time_jd[t]);
    }

    // Per-pixel results (collected thread-safely)
    struct PixelResult {
        int px;
        PixelEvents pe;
    };

    int nt = (n_threads > 0) ? n_threads : hw3::get_threads(npixels);

    std::vector<std::vector<PixelResult>> thread_results;
    thread_results.resize(nt);

    std::atomic<int> done_pixels{0};
    int report_interval = std::max(1, npixels / 20);

    #pragma omp parallel for schedule(dynamic) num_threads(nt)
    for (int px = 0; px < npixels; ++px) {
        const double* pixel_temp = sst + static_cast<size_t>(px) * ntime;
        const double* pixel_seas = seas_clim + static_cast<size_t>(px) * 366;
        const double* pixel_thresh = thresh_clim + static_cast<size_t>(px) * 366;

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

        // Skip pixels with no valid climatology
        bool has_clim = false;
        for (int d = 0; d < 366; ++d) {
            if (!is_na(pixel_seas[d])) { has_clim = true; break; }
        }
        if (!has_clim) {
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

        PixelEvents pe = detect_pixel_events(
            pixel_temp, pixel_seas, pixel_thresh,
            time_jd, doy_arr.data(), ntime,
            minDuration, minDuration2,
            joinAcrossGaps, maxGap, maxGap2,
            coldSpells, roundRes,
            category, southHemisphere
        );

        if (pe.valid && !pe.events.empty()) {
            int tid = omp_get_thread_num();
            PixelResult pr;
            pr.px = px;
            pr.pe = std::move(pe);
            thread_results[tid].push_back(std::move(pr));
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

    // Merge thread results
    for (auto& tr : thread_results) {
        for (auto& pr : tr) {
            int px = pr.px;
            int ilon = px / nlat;
            int ilat = px % nlat;

            for (auto& ev : pr.pe.events) {
                event_lon.push_back(lon[ilon]);
                event_lat.push_back(lat[ilat]);
                pixel_index.push_back(px);
                date_start.push_back(ev.jd_start);
                date_peak.push_back(ev.jd_peak);
                date_end.push_back(ev.jd_end);
                all_events.push_back(ev);
            }
        }
    }
}

} // namespace hw3
