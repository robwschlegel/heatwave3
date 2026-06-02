#ifndef HEATWAVE3_EVENT_DETECT_H
#define HEATWAVE3_EVENT_DETECT_H

#include <vector>
#include "heatwave3_types.h"

namespace hw3 {

struct PixelEvents {
    std::vector<EventResult> events;
    bool valid;
};

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
    bool southHemisphere
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
    int nlon, int nlat
);

} // namespace hw3

#endif
