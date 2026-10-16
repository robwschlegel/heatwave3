#ifndef HW3_API_H
#define HW3_API_H

// Language-neutral entry points shared by the R bindings (heatwave3_init.cpp)
// and the Python bindings. Plain structs and std::vector in and out; no R or
// Python types. Errors are thrown as std::runtime_error, which both Rcpp and
// pybind11 translate into native errors. Missing values are NaN for doubles
// and NA_INT for integers.

#include <climits>
#include <string>
#include <vector>
#include "heatwave3_types.h"
#include "netcdf_io.h"
#include "climatology.h"

namespace hw3 {

// Integer missing value; the same bit pattern as R's NA_INTEGER.
constexpr int NA_INT = INT_MIN;

// Parse "YYYY-MM-DD" to a Julian Day Number.
int parse_date_to_jd(const std::string& s);

// Build a SubsetSpec. Each range is either empty (unset) or {min, max};
// time_range holds "YYYY-MM-DD" strings. `depth` (>= 0) is the legacy
// single-level index; it is mutually exclusive with depth_range.
SubsetSpec make_subset(const std::vector<double>& lon_range,
                       const std::vector<double>& lat_range,
                       const std::vector<std::string>& time_range,
                       int depth,
                       const std::vector<double>& depth_range);

// ---- SST reading ----

GridData read_sst(const std::string& file_in, std::string var_name,
                  const SubsetSpec& ss);

GridData read_sst_multi(const std::vector<std::string>& files,
                        std::string var_name, const SubsetSpec& ss,
                        bool skip_bad_files);

// ---- Climatology ----

struct ClimParams {
    int maxPadLength = 0;
    int windowHalfWidth = 5;
    double pctile = 90;
    bool smoothPercentile = true;
    int smoothPercentileWidth = 31;
    bool compute_var = false;
    int roundClm = 4;
    int n_threads = 1;
    bool detrend = false;
    int compression = 1;  // zlib level of the output (0 = none)
};

PixelClimResult pixel_clim(const std::vector<double>& temp,
                           const std::vector<int>& time_jd,
                           const std::string& clim_start,
                           const std::string& clim_end,
                           const ClimParams& p);

void compute_clim(const std::string& file_in, const std::string& file_out,
                  std::string var_name,
                  const std::string& clim_start, const std::string& clim_end,
                  const SubsetSpec& ss, const ClimParams& p);

void compute_clim_multi(const std::vector<std::string>& files,
                        const std::string& file_out, std::string var_name,
                        const std::string& clim_start, const std::string& clim_end,
                        const SubsetSpec& ss, const ClimParams& p,
                        bool skip_bad_files);

void write_const_clim(const std::string& file_in, const std::string& clim_file,
                      std::string var_name, double threshold,
                      const SubsetSpec& ss);

// ---- Event detection ----

struct DetectParams {
    int minDuration = 5;
    int minDuration2 = 5;
    bool joinAcrossGaps = true;
    int maxGap = 2;
    int maxGap2 = 2;
    bool coldSpells = false;
    int roundRes = 4;
    int n_threads = 1;
    bool category = false;
    bool southHemisphere = true;
    int compression = 1;  // zlib level of the outputs (0 = none)
    bool slim = false;    // daily file: only intensity + category
};

// Output is selected by which paths are non-empty (events_file, daily_file,
// proto_file); an empty threshClim2 path/vector skips the secondary pass.
void detect_events(const std::string& file_in, const std::string& clim_file,
                   const std::string& events_file, std::string var_name,
                   const DetectParams& p,
                   const std::string& threshClim2_file,
                   const std::string& threshClim2_var_name,
                   const std::string& daily_file, const std::string& proto_file);

void detect_events_multi(const std::vector<std::string>& files,
                         const std::string& clim_file,
                         const std::string& events_file, std::string var_name,
                         const DetectParams& p,
                         const std::vector<std::string>& threshClim2_files,
                         const std::string& threshClim2_var_name,
                         bool skip_bad_files,
                         const std::string& daily_file, const std::string& proto_file);

// ---- Per-day detection on arrays (detect_event_array3) ----

// Output buffers of detect_events_array(), each [npixels * ntime] and owned by
// the caller; a null member is not computed. event is 0/1 in one byte (so a
// numpy bool array can be passed); category is 0 (none) to 4.
struct ArrayDetectOutputs {
    signed char* event = nullptr;
    int* event_no = nullptr;
    float* intensity = nullptr;
    signed char* category = nullptr;
};

// temp, seas and thresh: [npixels * ntime] pixel-major, gap-free daily series,
// each day paired with its own seas and thresh. Uses p.minDuration,
// p.joinAcrossGaps, p.maxGap, p.coldSpells, p.roundRes and p.n_threads
// (0 = session default). No files and no console output.
void detect_events_array(const double* temp, const double* seas, const double* thresh,
                         int npixels, int ntime, const DetectParams& p,
                         const ArrayDetectOutputs& out);

// ---- Daily categories for a date window ----

// Full per-pixel per-day grid. event is 0/1; event_no and category are NA_INT
// where not applicable; depth is filled only when has_depth.
struct DailyCatResult {
    bool has_depth = false;
    std::vector<int> jd;
    std::vector<double> lon, lat, depth;
    std::vector<double> temp, seas, thresh, intensity;
    std::vector<int> event, event_no, category;
};

DailyCatResult category_daily(const std::string& sst_file,
                              const std::string& clim_file,
                              const std::string& event_file,
                              const std::vector<int>& t_jd_range,
                              std::string var_name, bool coldSpells,
                              double ice_thresh, int roundRes);

DailyCatResult category_daily_multi(const std::vector<std::string>& files,
                                    const std::string& clim_file,
                                    const std::string& event_file,
                                    const std::vector<int>& t_jd_range,
                                    std::string var_name, bool coldSpells,
                                    double ice_thresh, int roundRes,
                                    bool skip_bad_files);

// ---- Event-file summaries ----

// One row per event. category: 0 = NA, 1..4 = I Moderate..IV Extreme;
// p_*: integer percentages of the duration (NaN = NA);
// season: 0 = unset, 1..4 = Summer, Fall, Winter, Spring.
struct CategoryResult {
    std::vector<int> event_no, peak_jd, category, duration, season;
    std::vector<double> lon, lat, intensity_max;
    std::vector<double> p_moderate, p_strong, p_severe, p_extreme;
};

// Categories precomputed in the events file are returned as stored;
// otherwise they are computed as heatwaveR::category() from the climatology
// and the SST (sst_files: one file or several daily files).
CategoryResult category(const std::string& event_file,
                        const std::string& clim_file,
                        const std::vector<std::string>& sst_files,
                        std::string var_name, bool skip_bad_files,
                        const std::string& hemisphere, int roundVal);

// One row per (pixel, year), ordered by (pixel_index, year).
struct BlockAverageResult {
    std::vector<double> lon, lat;
    std::vector<int> year, count, duration_max, total_days;
    std::vector<double> duration_mean, intensity_mean, intensity_max_mean;
    std::vector<double> intensity_max_max, intensity_cumulative_mean, total_icum;
};

BlockAverageResult block_average(const std::string& event_file);

// One row per pixel, ordered by pixel_index.
struct MetricSummaryResult {
    std::vector<double> lon, lat, value;
};

MetricSummaryResult metric_summary(const std::string& event_file,
                                   const std::string& metric,
                                   const std::string& summary);

} // namespace hw3

#endif
