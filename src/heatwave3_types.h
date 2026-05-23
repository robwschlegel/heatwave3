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
    int depth_index = -1;
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
};

} // namespace hw3

#endif
