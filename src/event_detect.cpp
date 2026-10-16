#include "event_detect.h"
#include "hw3_round.h"
#include "climatology.h"
#include "hw3_omp.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <cstring>
#include <mutex>
#include <atomic>
#include "hw3_log.h"

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
    int maxGap,
    signed char* duration_mask = nullptr) {

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
                // Mark the duration criterion (runs >= minDuration, pre-join)
                if (duration_mask) {
                    for (int k = start; k <= i - 1; ++k) duration_mask[k] = 1;
                }
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

// R's round(), as heatwaveR (see hw3_round.h)
static double round_to(double val, int decimals) {
    return r_round(val, decimals);
}

// ---- Per-day Hobday (2018) category ----
// 0 = no exceedance (or undefined), 1 = I Moderate ... 4 = IV Extreme.
// Mirrors the inline event-category thresholds: exceedance measured in
// multiples of |thresh - seas| above (heatwave) or below (cold-spell) seas.
//
// The arithmetic lives in the one-day helpers below, shared by the
// day-of-year form (CategoryLevels, for 366-day climatologies) and the
// per-day form used by detect_dense_events(), so both give identical levels.

// One day's four levels: the threshold, then thresh + k x diff for k = 1..3,
// diff = round(thresh - seas, digits), each level rounded and built
// cumulatively, as heatwaveR. All NaN when seas or thresh is missing.
static void day_levels(double s, double th, int digits, double* lv) {
    if (is_na(s) || is_na(th)) {
        for (int k = 0; k < 4; ++k) lv[k] = NA_DOUBLE;
        return;
    }
    double diff = round_to(th - s, digits);
    lv[0] = th;
    for (int k = 1; k < 4; ++k) lv[k] = round_to(lv[k - 1] + diff, digits);
}

// The daily category of one day from its four levels (see daily_category())
static int day_category(double temp, const double* lv, bool coldSpells) {
    if (is_na(temp)) return 0;
    for (int k = 0; k < 4; ++k)
        if (is_na(lv[k])) return 0;
    for (int k = 3; k >= 0; --k)
        if (coldSpells ? (temp < lv[k]) : (temp > lv[k])) return k + 1;
    return 0;
}

// The event category over n days (see event_category()). levels_of_day(i, lv)
// writes day i's four levels to lv and returns false when the day has none.
template <typename LevelsOfDay>
static EventCategory event_category_from(const double* temp, int n, int duration,
                                         LevelsOfDay levels_of_day, bool coldSpells) {
    int count[4] = {0, 0, 0, 0};
    double lv[4];
    for (int i = 0; i < n; ++i) {
        double t = temp[i];
        if (is_na(t) || !levels_of_day(i, lv)) continue;
        for (int k = 0; k < 4; ++k) {
            if (!is_na(lv[k]) && (coldSpells ? (t <= lv[k]) : (t >= lv[k]))) ++count[k];
        }
    }
    EventCategory ec;
    if (count[0] == 0) {
        ec.p_moderate = ec.p_strong = ec.p_severe = ec.p_extreme = NA_DOUBLE;
        return ec;
    }
    // R: round((moderate - strong) / duration * 100); nearbyint rounds half
    // to even, as R does.
    auto pct = [&](int c) { return std::nearbyint(static_cast<double>(c) / duration * 100.0); };
    ec.p_moderate = pct(count[0] - count[1]);
    ec.p_strong   = pct(count[1] - count[2]);
    ec.p_severe   = pct(count[2] - count[3]);
    ec.p_extreme  = pct(count[3]);
    ec.category = ec.p_extreme > 0 ? 4 : ec.p_severe > 0 ? 3 : ec.p_strong > 0 ? 2 : 1;
    return ec;
}

void category_levels(const double* seas, const double* thresh, int digits,
                     CategoryLevels& out) {
    double lv[4];
    for (int d = 0; d < 366; ++d) {
        day_levels(seas[d], thresh[d], digits, lv);
        for (int k = 0; k < 4; ++k) out.level[k][d] = lv[k];
    }
}

// Day d's four levels from the day-of-year form; false when d is out of range
static bool doy_levels(const CategoryLevels& L, int d, double* lv) {
    if (d < 0 || d >= 366) return false;
    for (int k = 0; k < 4; ++k) lv[k] = L.level[k][d];
    return true;
}

int daily_category(double temp, int doy, const CategoryLevels& L, bool coldSpells) {
    double lv[4];
    if (is_na(temp) || !doy_levels(L, doy - 1, lv)) return 0;
    return day_category(temp, lv, coldSpells);
}

EventCategory event_category(const double* temp, const int* doy, int n, int duration,
                             const CategoryLevels& L, bool coldSpells) {
    return event_category_from(temp, n, duration,
                               [&](int i, double* lv) { return doy_levels(L, doy[i] - 1, lv); },
                               coldSpells);
}

// ---- Event detection on a dense daily series ----
//
// Everything after the dense daily expansion of detect_pixel_events():
// threshold exceedance, proto-events, per-day output and event metrics.
// Inputs are on a dense daily grid of n days: dense_temp (missing days NaN),
// dense_thresh2 (0/1, or nullptr), and each day's own seas and thresh.
// dense_jd (Julian days) is needed only for event dates and seasons; with
// dense_jd == nullptr, want_events must be false.
// Per-day output goes to the caller's grid of out_n days: out_temp[t] is that
// day's temperature and out_index[t] its dense index (nullptr: identity).
// With want_events == false the per-day output is filled and no event
// metrics are computed (the returned events are empty).
// levels_of_day(t, lv) writes dense day t's four category levels to lv and
// returns false when the day has none; it is called only when a category is
// wanted, and only on event days. detect_pixel_events() looks the levels up
// in a 366-day table built once per pixel, which is cheaper than rounding
// four levels for every day of a long series.
template <typename LevelsOfDay>
static PixelEvents detect_dense_events(
    const double* dense_temp,
    const char* dense_thresh2,
    const double* daily_seas,
    const double* daily_thresh,
    const int* dense_jd,
    int dntime,
    const double* out_temp,
    const int* out_index,
    int out_n,
    int minDuration,
    int minDuration2,
    bool joinAcrossGaps,
    int maxGap,
    int maxGap2,
    bool coldSpells,
    int roundRes,
    bool category,
    bool southHemisphere,
    const DailyBuffers* daily,
    bool want_events,
    LevelsOfDay levels_of_day) {

    PixelEvents result;
    result.valid = false;

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

    // Optional dense per-day masks for daily output (filled on the dense grid,
    // then mapped back to the caller's input grid below).
    const bool want_daily = (daily != nullptr);
    std::vector<signed char> dense_dur, dense_evt;
    std::vector<int> dense_eno;
    if (want_daily) {
        dense_dur.assign(dntime, 0);
        dense_evt.assign(dntime, 0);
        dense_eno.assign(dntime, 0);
    }

    // Step 2: Detect proto-events — use manual loop since vector<bool> lacks .data()
    // The duration_mask is captured from the FINAL proto_event pass (the
    // secondary one when threshClim2 is supplied, otherwise the primary).
    signed char* dur_mask = (want_daily && dense_thresh2 == nullptr) ? dense_dur.data() : nullptr;
    auto events = proto_event(thresh_criterion, dntime,
                              minDuration, joinAcrossGaps, maxGap, dur_mask);

    if (dense_thresh2 != nullptr && !events.empty()) {
        std::vector<bool> primary_event(dntime, false);
        for (const auto& pe : events) {
            for (int t = pe.start; t <= pe.end; ++t) {
                primary_event[t] = true;
            }
        }

        std::vector<bool> secondary_criterion(dntime, false);
        for (int t = 0; t < dntime; ++t) {
            secondary_criterion[t] = primary_event[t] && (dense_thresh2[t] != 0);
        }

        events = proto_event(secondary_criterion, dntime,
                             minDuration2, joinAcrossGaps, maxGap2,
                             want_daily ? dense_dur.data() : nullptr);
    }

    // Fill the per-day output buffers (when requested) before any early return,
    // so that valid pixels with zero events still receive a complete daily
    // series (seas/thresh/threshCriterion present, event flags all 0).
    if (want_daily) {
        for (const auto& pe : events) {
            for (int k = pe.start; k <= pe.end && k < dntime; ++k) {
                dense_evt[k] = 1;
                dense_eno[k] = pe.event_no;
            }
        }
        for (int t = 0; t < out_n; ++t) {
            int di = out_index ? out_index[t] : t;
            double s = NA_DOUBLE, th = NA_DOUBLE, intens = NA_DOUBLE;
            signed char tc = 0, dc = 0, ev = 0, cat = 0;
            int eno = 0;
            if (di >= 0 && di < dntime) {
                s = daily_seas[di];
                th = daily_thresh[di];
                tc = thresh_raw[di];
                dc = dense_dur[di];
                ev = dense_evt[di];
                eno = dense_eno[di];
                double tv = out_temp[t];
                if (!is_na(tv) && !is_na(s)) {
                    intens = tv - s;
                    if (roundRes > 0) intens = round_to(intens, roundRes);
                }
                // Hobday (2018) categories describe events: only event days
                // are categorised, by heatwaveR's rule (category_daily3() too).
                double lv[4];
                if (ev && daily->category && levels_of_day(di, lv))
                    cat = static_cast<signed char>(day_category(tv, lv, coldSpells));
            }
            if (daily->seas) daily->seas[t] = static_cast<float>(s);
            if (daily->thresh) daily->thresh[t] = static_cast<float>(th);
            if (daily->threshCriterion) daily->threshCriterion[t] = tc;
            if (daily->durationCriterion) daily->durationCriterion[t] = dc;
            if (daily->event) daily->event[t] = ev;
            if (daily->event_no) daily->event_no[t] = eno;
            if (daily->intensity) daily->intensity[t] = static_cast<float>(intens);
            if (daily->category) daily->category[t] = cat;
        }
    }

    if (events.empty() || !want_events) {
        result.valid = true;
        return result;
    }

    // Step 3: Compute metrics for each event
    for (auto& pe : events) {
        EventResult er{};
        er.event_no = pe.event_no;
        er.category = 0;
        er.p_moderate = er.p_strong = er.p_severe = er.p_extreme = 0.0;
        er.season = 0;
        er.index_start = pe.start;
        er.index_end = pe.end;
        er.jd_start = dense_jd[pe.start];
        er.jd_end   = dense_jd[pe.end];
        er.duration = pe.end - pe.start + 1;

        // Event days' intensities relative to seas and thresh, and absolute
        // (sign-flipped for cold-spells, as heatwaveR), summarised with R's
        // mean(), sum() and sd() so the rounded metrics match heatwaveR.
        int nd = pe.end - pe.start + 1;
        std::vector<double> rel_seas_v(nd), rel_thresh_v(nd), abs_v(nd);
        double max_rel_seas = -1e30, max_rel_thresh = -1e30, max_abs = -1e30;
        int peak_idx = pe.start;
        for (int t = pe.start; t <= pe.end; ++t) {
            double rel_seas = temp_filled[t] - daily_seas[t];
            double rel_thresh = temp_filled[t] - daily_thresh[t];
            double abs_val = temp_filled[t];
            if (coldSpells) {
                rel_seas = -rel_seas;
                rel_thresh = -rel_thresh;
            }
            rel_seas_v[t - pe.start] = rel_seas;
            rel_thresh_v[t - pe.start] = rel_thresh;
            abs_v[t - pe.start] = abs_val;
            if (rel_seas > max_rel_seas) {
                max_rel_seas = rel_seas;
                peak_idx = t;
            }
            if (rel_thresh > max_rel_thresh) max_rel_thresh = rel_thresh;
            if (coldSpells) {
                if (abs_val < max_abs || max_abs < -1e29) max_abs = abs_val;
            } else {
                if (abs_val > max_abs) max_abs = abs_val;
            }
        }
        int n_valid = nd;

        er.index_peak = peak_idx;
        er.jd_peak    = dense_jd[peak_idx];

        if (n_valid > 0) {
            er.intensity_mean = r_mean(rel_seas_v.data(), nd);
            er.intensity_max = max_rel_seas;
            er.intensity_cumulative = r_sum(rel_seas_v.data(), nd);
            er.intensity_var = r_sd(rel_seas_v.data(), nd);

            er.intensity_mean_relThresh = r_mean(rel_thresh_v.data(), nd);
            er.intensity_max_relThresh = max_rel_thresh;
            er.intensity_cumulative_relThresh = r_sum(rel_thresh_v.data(), nd);
            er.intensity_var_relThresh = r_sd(rel_thresh_v.data(), nd);

            er.intensity_mean_abs = r_mean(abs_v.data(), nd);
            er.intensity_max_abs = max_abs;
            er.intensity_cumulative_abs = r_sum(abs_v.data(), nd);
            er.intensity_var_abs = r_sd(abs_v.data(), nd);

            // Rates of onset and decline, as heatwaveR (and Hobday et al.
            // 2016): from the midpoint between the first day and the day
            // before it, and to the midpoint between the last day and the day
            // after it; NA when that day is outside the series or missing.
            // Computed on the sign-flipped series for cold-spells, as heatwaveR.
            const double sgn = coldSpells ? -1.0 : 1.0;
            auto rel = [&](int t) { return sgn * temp_filled[t] - sgn * daily_seas[t]; };
            auto mid = [&](int edge, int outside) {
                // 0.5 * (A + B - C), in heatwaveR's order of operations
                return 0.5 * (rel(edge) + sgn * dense_temp[outside] - sgn * daily_seas[outside]);
            };
            int days_to_peak = peak_idx - pe.start;
            int days_from_peak = pe.end - peak_idx;
            er.rate_onset = (pe.start > 0) ?
                (max_rel_seas - mid(pe.start, pe.start - 1)) / (days_to_peak + 0.5) : NA_DOUBLE;
            er.rate_decline = (pe.end < dntime - 1) ?
                (max_rel_seas - mid(pe.end, pe.end + 1)) / (days_from_peak + 0.5) : NA_DOUBLE;

            // Cold-spells: intensities and rates are reported negative, as
            // heatwaveR (the variances and the _abs metrics are unchanged).
            if (coldSpells) {
                er.intensity_mean = -er.intensity_mean;
                er.intensity_max = -er.intensity_max;
                er.intensity_cumulative = -er.intensity_cumulative;
                er.intensity_mean_relThresh = -er.intensity_mean_relThresh;
                er.intensity_max_relThresh = -er.intensity_max_relThresh;
                er.intensity_cumulative_relThresh = -er.intensity_cumulative_relThresh;
                er.rate_onset = -er.rate_onset;
                er.rate_decline = -er.rate_decline;
            }
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
            // heatwaveR::category(): days at or beyond each level, over the
            // event's days, against the original (unfilled) temperatures
            const int start = pe.start;
            EventCategory ec = event_category_from(
                dense_temp + start, pe.end - start + 1, er.duration,
                [&](int i, double* lv) { return levels_of_day(start + i, lv); },
                coldSpells);
            er.category = ec.category;
            er.p_moderate = ec.p_moderate;
            er.p_strong = ec.p_strong;
            er.p_severe = ec.p_severe;
            er.p_extreme = ec.p_extreme;

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

PixelEvents detect_pixel_events(
    const double* temp,
    const double* thresh2,
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
    bool southHemisphere,
    const DailyBuffers* daily) {

    if (ntime <= 0) {
        PixelEvents result;
        result.valid = false;
        return result;
    }

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
    std::vector<char> dense_thresh2;
    if (thresh2 != nullptr) {
        dense_thresh2.assign(dntime, 0);
        for (int t = 0; t < ntime; ++t) {
            int idx = time_jd[t] - first_jd;
            if (idx >= 0 && idx < dntime) {
                dense_thresh2[idx] = (!is_na(thresh2[t]) && thresh2[t] != 0.0) ? 1 : 0;
            }
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

    // Per-day output is written on the caller's (possibly sparse) grid
    std::vector<int> out_index;
    if (daily != nullptr) {
        out_index.resize(ntime);
        for (int t = 0; t < ntime; ++t) out_index[t] = time_jd[t] - first_jd;
    }

    // heatwaveR category levels for this pixel (per day of year), when any
    // category output is wanted
    const bool want_levels = category || (daily != nullptr && daily->category != nullptr);
    CategoryLevels levels;
    if (want_levels) category_levels(seas, thresh, 4, levels);

    return detect_dense_events(
        dense_temp.data(), thresh2 != nullptr ? dense_thresh2.data() : nullptr,
        daily_seas.data(), daily_thresh.data(), dense_jd.data(), dntime,
        temp, out_index.empty() ? nullptr : out_index.data(), ntime,
        minDuration, minDuration2, joinAcrossGaps, maxGap, maxGap2,
        coldSpells, roundRes, category, southHemisphere, daily,
        /*want_events=*/true,
        [&](int t, double* lv) { return doy_levels(levels, dense_doy[t] - 1, lv); });
}

// ---- Grid-level event detection with OpenMP ----

void detect_events_grid(
    const double* sst,
    const double* thresh2,
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
    int nlon, int nlat,
    DailyBuffers* daily_grid,
    int ndepth,
    const std::vector<double>& depth,
    std::vector<double>* event_depth) {

    // Compute DOY for each time step
    std::vector<int> doy_arr(ntime);
    for (int t = 0; t < ntime; ++t) {
        doy_arr[t] = jd_to_doy_366(time_jd[t]);
    }

    // Per-day output buffers arrive uninitialised. Each worker first sets its
    // pixel's slice to the defaults (NA for the real-valued fields, 0 for the
    // flag/integer fields), so skipped (land / no-climatology) pixels keep them
    // and valid pixels overwrite them in detect_pixel_events. Done per pixel on
    // the worker threads: a serial fill of the whole grid took ~0.3 s on 80M cells.
    const bool want_daily = (daily_grid != nullptr);
    const float NA_F = std::numeric_limits<float>::quiet_NaN();

    // Per-pixel results (collected thread-safely)
    struct PixelResult {
        int px;
        PixelEvents pe;
    };

    int nt = (n_threads > 0) ? n_threads : hw3::get_threads(npixels);

    // Each worker accumulates into its own results vector (indexed by worker id),
    // merged after the parallel region; no locking in the hot loop.
    std::vector<std::vector<PixelResult>> thread_results(nt);

    auto report = [&](int cur) {
        hw3::out() << "\r  " << cur << "/" << npixels << " pixels ("
                    << (100 * cur / npixels) << "%)" << std::flush;
    };

    hw3::parallel_for(npixels, nt, [&](int wid, int px) {
        const double* pixel_temp = sst + static_cast<size_t>(px) * ntime;
        const double* pixel_thresh2 = thresh2 ? thresh2 + static_cast<size_t>(px) * ntime : nullptr;
        const double* pixel_seas = seas_clim + static_cast<size_t>(px) * 366;
        const double* pixel_thresh = thresh_clim + static_cast<size_t>(px) * 366;

        if (want_daily) {
            size_t off = static_cast<size_t>(px) * ntime;
            auto fill = [&](auto* buf, auto v) { if (buf) std::fill(buf + off, buf + off + ntime, v); };
            fill(daily_grid->seas, NA_F);
            fill(daily_grid->thresh, NA_F);
            fill(daily_grid->intensity, NA_F);
            fill(daily_grid->threshCriterion, static_cast<signed char>(0));
            fill(daily_grid->durationCriterion, static_cast<signed char>(0));
            fill(daily_grid->event, static_cast<signed char>(0));
            fill(daily_grid->category, static_cast<signed char>(0));
            fill(daily_grid->event_no, 0);
        }

        // Skip all-NA pixels
        bool has_data = false;
        for (int t = 0; t < ntime; ++t) {
            if (!is_na(pixel_temp[t])) { has_data = true; break; }
        }
        if (!has_data) return;

        // Skip pixels with no valid climatology
        bool has_clim = false;
        for (int d = 0; d < 366; ++d) {
            if (!is_na(pixel_seas[d])) { has_clim = true; break; }
        }
        if (!has_clim) return;

        // Build this pixel's slice of the grid-wide daily buffers (if any).
        DailyBuffers pdaily;
        DailyBuffers* pdaily_ptr = nullptr;
        if (want_daily) {
            size_t off = static_cast<size_t>(px) * ntime;
            if (daily_grid->seas)      pdaily.seas = daily_grid->seas + off;
            if (daily_grid->thresh)    pdaily.thresh = daily_grid->thresh + off;
            if (daily_grid->intensity) pdaily.intensity = daily_grid->intensity + off;
            if (daily_grid->threshCriterion)   pdaily.threshCriterion = daily_grid->threshCriterion + off;
            if (daily_grid->durationCriterion) pdaily.durationCriterion = daily_grid->durationCriterion + off;
            if (daily_grid->event)     pdaily.event = daily_grid->event + off;
            if (daily_grid->category)  pdaily.category = daily_grid->category + off;
            if (daily_grid->event_no)  pdaily.event_no = daily_grid->event_no + off;
            pdaily_ptr = &pdaily;
        }

        PixelEvents pe = detect_pixel_events(
            pixel_temp, pixel_thresh2, pixel_seas, pixel_thresh,
            time_jd, doy_arr.data(), ntime,
            minDuration, minDuration2,
            joinAcrossGaps, maxGap, maxGap2,
            coldSpells, roundRes,
            category, southHemisphere,
            pdaily_ptr
        );

        if (pe.valid && !pe.events.empty()) {
            PixelResult pr;
            pr.px = px;
            pr.pe = std::move(pe);
            thread_results[wid].push_back(std::move(pr));
        }
    }, report);

    hw3::out() << std::endl;

    // Merge thread results
    bool want_event_depth = (ndepth > 1) && (event_depth != nullptr) && !depth.empty();
    for (auto& tr : thread_results) {
        for (auto& pr : tr) {
            int px = pr.px;
            int per_lonlat = nlat * ndepth;
            int ilon = px / per_lonlat;
            int rem = px % per_lonlat;
            int ilat = rem / ndepth;
            int idepth = rem % ndepth;

            for (auto& ev : pr.pe.events) {
                event_lon.push_back(lon[ilon]);
                event_lat.push_back(lat[ilat]);
                pixel_index.push_back(px);
                date_start.push_back(ev.jd_start);
                date_peak.push_back(ev.jd_peak);
                date_end.push_back(ev.jd_end);
                all_events.push_back(ev);
                if (want_event_depth) event_depth->push_back(depth[idepth]);
            }
        }
    }
}

} // namespace hw3
