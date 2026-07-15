#ifndef HEATWAVE3_TYPES_H
#define HEATWAVE3_TYPES_H

#include <vector>
#include <string>
#include <cmath>
#include <limits>

namespace hw3 {

constexpr double NA_DOUBLE = std::numeric_limits<double>::quiet_NaN();

inline bool is_na(double x) {
    return std::isnan(x);
}

struct SubsetSpec {
    double lon_min = -std::numeric_limits<double>::infinity();
    double lon_max =  std::numeric_limits<double>::infinity();
    double lat_min = -std::numeric_limits<double>::infinity();
    double lat_max =  std::numeric_limits<double>::infinity();
    double time_min = -std::numeric_limits<double>::infinity();
    double time_max =  std::numeric_limits<double>::infinity();
    // Legacy single-level squeeze: select one depth index and drop the depth
    // dimension entirely (existing behaviour, unchanged).
    int depth_index = -1;
    // New: keep a contiguous *range* of depth levels (inclusive, by depth
    // value in the file's own units, typically metres) as a real dimension
    // in the output. Mutually exclusive with depth_index; ignored unless
    // finite. Default (both infinite) preserves the legacy squeeze-to-index
    // behaviour above.
    double depth_min = -std::numeric_limits<double>::infinity();
    double depth_max =  std::numeric_limits<double>::infinity();
};

struct PixelClim {
    std::vector<double> seas;
    std::vector<double> thresh;
    std::vector<double> var;
};

struct EventResult {
    int event_no;
    int index_start;
    int index_peak;
    int index_end;
    int jd_start;
    int jd_peak;
    int jd_end;
    int duration;
    double intensity_mean;
    double intensity_max;
    double intensity_var;
    double intensity_cumulative;
    double intensity_mean_relThresh;
    double intensity_max_relThresh;
    double intensity_var_relThresh;
    double intensity_cumulative_relThresh;
    double intensity_mean_abs;
    double intensity_max_abs;
    double intensity_var_abs;
    double intensity_cumulative_abs;
    double rate_onset;
    double rate_decline;
    int category;       // 0=unset/NA, 1=I Moderate, 2=II Strong, 3=III Severe, 4=IV Extreme
    double p_moderate;
    double p_strong;
    double p_severe;
    double p_extreme;
    int season;         // 0=unset, 1=Summer, 2=Fall, 3=Winter, 4=Spring (hemisphere-dependent)
};

struct EventData {
    std::vector<double> lon;
    std::vector<double> lat;
    // Per-event depth (metres). Empty when the source data had no depth
    // range (ndepth == 1) -- matches how `category`/`season` are also only
    // populated when present in the file.
    std::vector<double> depth;
    std::vector<int> pixel_index;
    std::vector<int> event_no;
    std::vector<int> date_start;
    std::vector<int> date_peak;
    std::vector<int> date_end;
    std::vector<int> duration;
    std::vector<double> intensity_mean;
    std::vector<double> intensity_max;
    std::vector<double> intensity_var;
    std::vector<double> intensity_cumulative;
    std::vector<double> intensity_mean_relThresh;
    std::vector<double> intensity_max_relThresh;
    std::vector<double> intensity_var_relThresh;
    std::vector<double> intensity_cumulative_relThresh;
    std::vector<double> intensity_mean_abs;
    std::vector<double> intensity_max_abs;
    std::vector<double> intensity_var_abs;
    std::vector<double> intensity_cumulative_abs;
    std::vector<double> rate_onset;
    std::vector<double> rate_decline;
    std::vector<int> category;
    std::vector<double> p_moderate;
    std::vector<double> p_strong;
    std::vector<double> p_severe;
    std::vector<double> p_extreme;
    std::vector<int> season;
    int ref_date_jd;
    int nevents;
};

} // namespace hw3

#endif
