#ifndef HEATWAVE3_EVENT_DETECT_H
#define HEATWAVE3_EVENT_DETECT_H

#include <vector>
#include "heatwave3_types.h"

namespace hw3 {

struct PixelEvents {
    std::vector<EventResult> events;
    bool valid;
};

// Optional per-pixel daily output. Each non-null pointer addresses a slice of
// length == ntime (the caller's input time grid) into a grid-wide buffer; the
// caller offsets the pointers per pixel before the call. detect_pixel_events
// fills only the non-null arrays. Used for protoEvent and per-day output.
struct DailyBuffers {
    double* seas = nullptr;              // daily seasonal climatology
    double* thresh = nullptr;            // daily threshold climatology
    signed char* threshCriterion = nullptr;   // temp exceeds threshold (0/1)
    signed char* durationCriterion = nullptr; // in a >= minDuration run (0/1)
    signed char* event = nullptr;        // in a final, gap-joined event (0/1)
    int* event_no = nullptr;             // event number (0 = none)
    double* intensity = nullptr;         // temp - seas (rounded to roundRes)
    signed char* category = nullptr;     // daily Hobday category (0 = none, 1-4)
};

// ---- Hobday et al. (2018) categories, as heatwaveR::category() ----
//
// Category levels for one pixel, per day of year (index 0-365): the threshold
// and thresh + k x diff for k = 1, 2, 3, where diff = round(thresh - seas,
// digits) and each level is rounded to `digits`, built cumulatively as in
// heatwaveR. NaN where the climatology is missing.
struct CategoryLevels {
    double level[4][366];  // [0] thresh, [1] 2x, [2] 3x, [3] 4x
};
void category_levels(const double* seas, const double* thresh, int digits,
                     CategoryLevels& out);

// Daily category, as heatwaveR::category(climatology = TRUE): 1-4 when temp
// is strictly beyond the threshold .. 4x level (below it for cold-spells),
// 0 otherwise, when temp is missing or when any level is missing. doy is
// 1-366.
int daily_category(double temp, int doy, const CategoryLevels& L, bool coldSpells);

// Event category, as heatwaveR::category(): over the event's n days, count
// the days at or beyond each level (missing temperatures do not count); p_*
// are integer percentages of `duration` (rounded half to even, like R's
// round()), and the category is IV if p_extreme > 0, else III if
// p_severe > 0, else II if p_strong > 0, else I. When no day reaches the
// threshold the category is 0 (none) and p_* are NaN, as heatwaveR's NA.
struct EventCategory {
    int category = 0;
    double p_moderate, p_strong, p_severe, p_extreme;
};
EventCategory event_category(const double* temp, const int* doy, int n, int duration,
                             const CategoryLevels& L, bool coldSpells);

// Detect events for a single pixel
PixelEvents detect_pixel_events(
    const double* temp,        // daily temperature
    const double* thresh2,     // optional secondary logical criterion
    const double* seas,        // seasonal climatology (366 DOY)
    const double* thresh,      // threshold climatology (366 DOY)
    const int* time_jd,        // Julian Day numbers
    const int* doy,            // DOY for each time step (1-366)
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
    const DailyBuffers* daily = nullptr  // optional per-day output (this pixel's slice)
);

// Detect events for entire grid with OpenMP
void detect_events_grid(
    const double* sst,         // [npixels * ntime] pixel-major
    const double* thresh2,     // optional [npixels * ntime] secondary criterion
    const double* seas_clim,   // [npixels * 366]
    const double* thresh_clim, // [npixels * 366]
    const int* time_jd,        // [ntime]
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
    // Outputs (populated by function)
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
    // Optional grid-wide per-day output buffers. When non-null, each member
    // points to an [npixels * ntime] array (pixel-major) that is filled for
    // every valid pixel. Members left null are skipped. Drives protoEvent and
    // per-day output.
    DailyBuffers* daily_grid = nullptr,
    // Depth axis: pixel = (ilon * nlat + ilat) * ndepth + idepth. ndepth == 1
    // (depth empty) for ordinary 3D data, reducing pixel to ilon*nlat+ilat.
    // When ndepth > 1, event_depth is populated (one value per event,
    // metres) in parallel with event_lon/event_lat; left empty otherwise.
    int ndepth = 1,
    const std::vector<double>& depth = {},
    std::vector<double>* event_depth = nullptr
);

} // namespace hw3

#endif
