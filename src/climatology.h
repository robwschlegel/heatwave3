#ifndef HEATWAVE3_CLIMATOLOGY_H
#define HEATWAVE3_CLIMATOLOGY_H

#include <vector>
#include "heatwave3_types.h"

namespace hw3 {

struct PixelClimResult {
    double seas[366];
    double thresh[366];
    double var[366];
    bool valid;
};

// Day-of-year for a Gregorian date (Julian Day Number)
int jd_to_doy_366(int jd);
bool is_leap_year_from_jd(int jd);
void jd_to_ymd(int jd, int& y, int& m, int& d);
int ymd_to_jd(int y, int m, int d);

// Compute climatology for a single pixel
PixelClimResult compute_pixel_climatology(
    const double* temp,        // temperature time series
    const int* time_jd,        // Julian Day numbers
    int ntime,                 // length of time series
    int clim_start_jd,         // start of climatology period (JD)
    int clim_end_jd,           // end of climatology period (JD)
    int windowHalfWidth,       // default 5
    double pctile,             // default 90
    bool smoothPercentile,     // default true
    int smoothPercentileWidth, // default 31
    bool compute_var,          // default false
    int maxPadLength,          // max consecutive NAs to interpolate, 0 = no interpolation
    bool detrend               // remove linear trend before climatology calc
);

// Compute climatology for entire grid with OpenMP
void compute_climatology_grid(
    const double* sst,         // [npixels * ntime] pixel-major
    const int* time_jd,        // [ntime]
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
    double* seas_out,          // [npixels * 366]
    double* thresh_out,        // [npixels * 366]
    double* var_out            // [npixels * 366] (may be nullptr)
);

} // namespace hw3

#endif
