// R bindings. Converts between Rcpp types and the language-neutral API in
// hw3_api.h; all algorithm and orchestration logic lives there so the Python
// bindings share it.

#include <Rcpp.h>
#include "heatwave3_types.h"
#include "hw3_api.h"
#include "hw3_log.h"
#include "hw3_omp.h"
#include "netcdf_io.h"
#include "climatology.h"

// ---- Message sinks for the shared core (hw3_log.h) ----

namespace hw3 {
std::ostream& out() { return Rcpp::Rcout; }
std::ostream& err() { return Rcpp::Rcerr; }
}

// ---- Conversion helpers ----

static std::vector<double> opt_dbl(Rcpp::Nullable<Rcpp::NumericVector> x) {
    if (x.isNull()) return {};
    return Rcpp::as<std::vector<double>>(Rcpp::NumericVector(x));
}

static std::vector<std::string> opt_str(Rcpp::Nullable<Rcpp::CharacterVector> x) {
    if (x.isNull()) return {};
    return Rcpp::as<std::vector<std::string>>(Rcpp::CharacterVector(x));
}

// The core marks missing doubles with NaN; R's NA_real_ is a distinct NaN.
static Rcpp::NumericVector na_real(const std::vector<double>& v) {
    Rcpp::NumericVector out(v.size());
    for (size_t i = 0; i < v.size(); ++i)
        out[i] = hw3::is_na(v[i]) ? NA_REAL : v[i];
    return out;
}

static Rcpp::List grid_list(const hw3::GridData& gd) {
    return Rcpp::List::create(
        Rcpp::Named("lon") = gd.lon,
        Rcpp::Named("lat") = gd.lat,
        Rcpp::Named("nlon") = gd.nlon,
        Rcpp::Named("nlat") = gd.nlat,
        Rcpp::Named("ntime") = gd.ntime,
        Rcpp::Named("time_days") = gd.time_days,
        Rcpp::Named("ndepth") = gd.ndepth,
        Rcpp::Named("depth") = gd.depth,
        Rcpp::Named("sst") = gd.sst
    );
}

static hw3::ClimParams clim_params(int maxPadLength, int windowHalfWidth,
                                   double pctile, bool smoothPercentile,
                                   int smoothPercentileWidth, bool compute_var,
                                   int roundClm, int n_threads, bool detrend) {
    hw3::ClimParams p;
    p.maxPadLength = maxPadLength;
    p.windowHalfWidth = windowHalfWidth;
    p.pctile = pctile;
    p.smoothPercentile = smoothPercentile;
    p.smoothPercentileWidth = smoothPercentileWidth;
    p.compute_var = compute_var;
    p.roundClm = roundClm;
    p.n_threads = n_threads;
    p.detrend = detrend;
    return p;
}

static hw3::DetectParams detect_params(int minDuration, int minDuration2,
                                       bool joinAcrossGaps, int maxGap, int maxGap2,
                                       bool coldSpells, int roundRes, int n_threads,
                                       bool category, bool southHemisphere) {
    hw3::DetectParams p;
    p.minDuration = minDuration;
    p.minDuration2 = minDuration2;
    p.joinAcrossGaps = joinAcrossGaps;
    p.maxGap = maxGap;
    p.maxGap2 = maxGap2;
    p.coldSpells = coldSpells;
    p.roundRes = roundRes;
    p.n_threads = n_threads;
    p.category = category;
    p.southHemisphere = southHemisphere;
    return p;
}

// ---- Thread management (exported to R) ----

// [[Rcpp::export]]
int hw3_get_threads() {
    return hw3::get_threads();
}

// [[Rcpp::export]]
int hw3_set_threads(int n) {
    hw3::set_threads(n);
    return hw3::get_threads();
}

// [[Rcpp::export]]
void hw3_init_fork_safety() {
    hw3::hw3_init_fork_handler();
}

// [[Rcpp::export]]
Rcpp::List hw3_version() {
    return Rcpp::List::create(
        Rcpp::Named("package") = "heatwave3",
        Rcpp::Named("version") = "1.0.0",
        Rcpp::Named("cpp_standard") = "C++17"
    );
}

// [[Rcpp::export]]
Rcpp::List hw3_read_sst(std::string file_in,
                        std::string var_name,
                        Rcpp::Nullable<Rcpp::NumericVector> lon_range = R_NilValue,
                        Rcpp::Nullable<Rcpp::NumericVector> lat_range = R_NilValue,
                        Rcpp::Nullable<Rcpp::CharacterVector> time_range = R_NilValue,
                        int depth = -1,
                        Rcpp::Nullable<Rcpp::NumericVector> depth_range = R_NilValue) {
    hw3::SubsetSpec ss = hw3::make_subset(opt_dbl(lon_range), opt_dbl(lat_range),
                                          opt_str(time_range), depth, opt_dbl(depth_range));
    return grid_list(hw3::read_sst(file_in, var_name, ss));
}

// Read and merge multiple daily NetCDF files
// [[Rcpp::export]]
Rcpp::List hw3_read_sst_multi(Rcpp::CharacterVector files,
                              std::string var_name,
                              Rcpp::Nullable<Rcpp::NumericVector> lon_range = R_NilValue,
                              Rcpp::Nullable<Rcpp::NumericVector> lat_range = R_NilValue,
                              int depth = -1,
                              bool skip_bad_files = false,
                              Rcpp::Nullable<Rcpp::NumericVector> depth_range = R_NilValue) {
    hw3::SubsetSpec ss = hw3::make_subset(opt_dbl(lon_range), opt_dbl(lat_range),
                                          {}, depth, opt_dbl(depth_range));
    return grid_list(hw3::read_sst_multi(Rcpp::as<std::vector<std::string>>(files),
                                         var_name, ss, skip_bad_files));
}

// Debug: compute climatology for a single pixel from R vectors
// [[Rcpp::export]]
Rcpp::List hw3_pixel_clim(Rcpp::NumericVector temp,
                          Rcpp::IntegerVector time_jd,
                          std::string clim_start,
                          std::string clim_end,
                          int windowHalfWidth = 5,
                          double pctile = 90,
                          bool smoothPercentile = true,
                          int smoothPercentileWidth = 31,
                          int maxPadLength = 0,
                          bool detrend = false) {
    int n = temp.size();
    std::vector<double> t(n);
    std::vector<int> jd(n);
    for (int i = 0; i < n; ++i) {
        t[i] = Rcpp::NumericVector::is_na(temp[i]) ? hw3::NA_DOUBLE : temp[i];
        jd[i] = time_jd[i];
    }

    hw3::ClimParams p = clim_params(maxPadLength, windowHalfWidth, pctile,
                                    smoothPercentile, smoothPercentileWidth,
                                    false, 0, 1, detrend);
    hw3::PixelClimResult cr = hw3::pixel_clim(t, jd, clim_start, clim_end, p);

    Rcpp::NumericVector seas(366), thresh(366);
    for (int i = 0; i < 366; ++i) {
        seas[i] = hw3::is_na(cr.seas[i]) ? NA_REAL : cr.seas[i];
        thresh[i] = hw3::is_na(cr.thresh[i]) ? NA_REAL : cr.thresh[i];
    }

    return Rcpp::List::create(
        Rcpp::Named("seas") = seas,
        Rcpp::Named("thresh") = thresh,
        Rcpp::Named("valid") = cr.valid
    );
}

// Single-file climatology entry point
// [[Rcpp::export]]
void hw3_compute_clim(std::string file_in,
                      std::string file_out,
                      std::string var_name,
                      Rcpp::CharacterVector climatologyPeriod,
                      Rcpp::Nullable<Rcpp::NumericVector> lon_range = R_NilValue,
                      Rcpp::Nullable<Rcpp::NumericVector> lat_range = R_NilValue,
                      Rcpp::Nullable<Rcpp::CharacterVector> time_range = R_NilValue,
                      int depth = -1,
                      int maxPadLength = 0,
                      int windowHalfWidth = 5,
                      double pctile = 90,
                      bool smoothPercentile = true,
                      int smoothPercentileWidth = 31,
                      bool compute_var = false,
                      int roundClm = 4,
                      int n_threads = 1,
                      bool detrend = false,
                      Rcpp::Nullable<Rcpp::NumericVector> depth_range = R_NilValue) {
    hw3::SubsetSpec ss = hw3::make_subset(opt_dbl(lon_range), opt_dbl(lat_range),
                                          opt_str(time_range), depth, opt_dbl(depth_range));
    hw3::ClimParams p = clim_params(maxPadLength, windowHalfWidth, pctile,
                                    smoothPercentile, smoothPercentileWidth,
                                    compute_var, roundClm, n_threads, detrend);
    hw3::compute_clim(file_in, file_out, var_name,
                      Rcpp::as<std::string>(climatologyPeriod[0]),
                      Rcpp::as<std::string>(climatologyPeriod[1]), ss, p);
}

// Multi-file climatology entry point
// [[Rcpp::export]]
void hw3_compute_clim_multi(Rcpp::CharacterVector files,
                            std::string file_out,
                            std::string var_name,
                            Rcpp::CharacterVector climatologyPeriod,
                            Rcpp::Nullable<Rcpp::NumericVector> lon_range = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericVector> lat_range = R_NilValue,
                            int depth = -1,
                            int maxPadLength = 0,
                            int windowHalfWidth = 5,
                            double pctile = 90,
                            bool smoothPercentile = true,
                            int smoothPercentileWidth = 31,
                            bool compute_var = false,
                            int roundClm = 4,
                            int n_threads = 1,
                            bool detrend = false,
                            bool skip_bad_files = false,
                            Rcpp::Nullable<Rcpp::NumericVector> depth_range = R_NilValue) {
    hw3::SubsetSpec ss = hw3::make_subset(opt_dbl(lon_range), opt_dbl(lat_range),
                                          {}, depth, opt_dbl(depth_range));
    hw3::ClimParams p = clim_params(maxPadLength, windowHalfWidth, pctile,
                                    smoothPercentile, smoothPercentileWidth,
                                    compute_var, roundClm, n_threads, detrend);
    hw3::compute_clim_multi(Rcpp::as<std::vector<std::string>>(files), file_out, var_name,
                            Rcpp::as<std::string>(climatologyPeriod[0]),
                            Rcpp::as<std::string>(climatologyPeriod[1]),
                            ss, p, skip_bad_files);
}

// [[Rcpp::export]]
void hw3_write_const_clim(std::string file_in,
                          std::string clim_file,
                          std::string var_name,
                          double threshold,
                          Rcpp::Nullable<Rcpp::NumericVector> lon_range = R_NilValue,
                          Rcpp::Nullable<Rcpp::NumericVector> lat_range = R_NilValue,
                          Rcpp::Nullable<Rcpp::CharacterVector> time_range = R_NilValue,
                          int depth = -1) {
    hw3::SubsetSpec ss = hw3::make_subset(opt_dbl(lon_range), opt_dbl(lat_range),
                                          opt_str(time_range), depth, {});
    hw3::write_const_clim(file_in, clim_file, var_name, threshold, ss);
}

// [[Rcpp::export]]
int hw3_jd_to_doy(int jd) {
    return hw3::jd_to_doy_366(jd);
}

// [[Rcpp::export]]
Rcpp::List hw3_read_clim_nc(std::string clim_file) {
    hw3::ClimData cd = hw3::read_clim_netcdf(clim_file);
    return Rcpp::List::create(
        Rcpp::Named("lon") = cd.lon,
        Rcpp::Named("lat") = cd.lat,
        Rcpp::Named("nlon") = cd.nlon,
        Rcpp::Named("nlat") = cd.nlat,
        Rcpp::Named("ndoy") = cd.ndoy,
        Rcpp::Named("ndepth") = cd.ndepth,
        Rcpp::Named("depth") = cd.depth,
        Rcpp::Named("seas") = cd.seas,
        Rcpp::Named("thresh") = cd.thresh
    );
}

// [[Rcpp::export]]
Rcpp::List hw3_read_subset(std::string file,
                           Rcpp::Nullable<Rcpp::NumericVector> lon_range = R_NilValue,
                           Rcpp::Nullable<Rcpp::NumericVector> lat_range = R_NilValue,
                           Rcpp::Nullable<Rcpp::IntegerVector> t_jd_range = R_NilValue,
                           Rcpp::Nullable<Rcpp::NumericVector> depth_range = R_NilValue,
                           Rcpp::Nullable<Rcpp::CharacterVector> vars = R_NilValue,
                           int max_rows = -1) {
    std::vector<int> tr;
    if (t_jd_range.isNotNull()) tr = Rcpp::as<std::vector<int>>(Rcpp::IntegerVector(t_jd_range));
    hw3::SubsetResult s = hw3::read_subset_netcdf(file, opt_dbl(lon_range), opt_dbl(lat_range),
                                                  tr, opt_dbl(depth_range), opt_str(vars),
                                                  max_rows);
    Rcpp::List data;
    for (const auto& v : s.data) {
        if (v.is_double) data[v.name] = v.dbl;
        else             data[v.name] = v.intg;
    }
    return Rcpp::List::create(
        Rcpp::Named("product") = s.product,
        Rcpp::Named("lon") = s.lon,
        Rcpp::Named("lat") = s.lat,
        Rcpp::Named("has_depth") = s.has_depth,
        Rcpp::Named("depth") = s.depth,
        Rcpp::Named("third") = s.third,
        Rcpp::Named("third_name") = s.third_name,
        Rcpp::Named("nlon") = s.nlon,
        Rcpp::Named("nlat") = s.nlat,
        Rcpp::Named("ndepth") = s.ndepth,
        Rcpp::Named("n3") = s.n3,
        Rcpp::Named("data") = data
    );
}

static Rcpp::List daily_cat_list(const hw3::DailyCatResult& r) {
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("jd") = r.jd,
        Rcpp::Named("lon") = r.lon,
        Rcpp::Named("lat") = r.lat,
        Rcpp::Named("temp") = na_real(r.temp),
        Rcpp::Named("seas") = na_real(r.seas),
        Rcpp::Named("thresh") = na_real(r.thresh),
        Rcpp::Named("intensity") = na_real(r.intensity),
        Rcpp::Named("event") = Rcpp::LogicalVector(r.event.begin(), r.event.end()),
        Rcpp::Named("event_no") = r.event_no,
        Rcpp::Named("category") = r.category
    );
    if (r.has_depth) out["depth"] = r.depth;
    return out;
}

// Daily marine-heatwave / cold-spell categories for a date window, built from a
// single multi-time SST file + a heatwave3 climatology + a heatwave3 events
// file. Only the SST time-window is read (hyperslab); no detection is re-run.
// [[Rcpp::export]]
Rcpp::List hw3_category_daily(std::string sst_file,
                              std::string clim_file,
                              std::string event_file,
                              std::vector<int> t_jd_range,   // [jd0, jd1]
                              std::string var_name = "",
                              bool coldSpells = false,
                              double ice_thresh = -1.7,
                              int roundRes = 2) {
    return daily_cat_list(hw3::category_daily(sst_file, clim_file, event_file, t_jd_range,
                                              var_name, coldSpells, ice_thresh, roundRes));
}

// As above, but reading the SST window from a vector of daily files (one time
// step each), e.g. daily-global OSTIA/GHRSST files for the requested window.
// [[Rcpp::export]]
Rcpp::List hw3_category_daily_multi(Rcpp::CharacterVector files,
                                    std::string clim_file,
                                    std::string event_file,
                                    std::vector<int> t_jd_range,
                                    std::string var_name = "",
                                    bool coldSpells = false,
                                    double ice_thresh = -1.7,
                                    int roundRes = 2,
                                    bool skip_bad_files = false) {
    return daily_cat_list(hw3::category_daily_multi(
        Rcpp::as<std::vector<std::string>>(files), clim_file, event_file, t_jd_range,
        var_name, coldSpells, ice_thresh, roundRes, skip_bad_files));
}

// [[Rcpp::export]]
Rcpp::List hw3_file_meta(std::string file) {
    hw3::FileMeta m = hw3::read_file_meta(file);
    return Rcpp::List::create(
        Rcpp::Named("product") = m.product,
        Rcpp::Named("nlon") = m.nlon,
        Rcpp::Named("nlat") = m.nlat,
        Rcpp::Named("n3") = m.n3,
        Rcpp::Named("nrows") = static_cast<double>(m.nrows)
    );
}

// [[Rcpp::export]]
Rcpp::List hw3_read_daily_nc(std::string daily_file) {
    hw3::DailyData dd = hw3::read_daily_netcdf(daily_file);
    auto schar_to_int = [](const std::vector<signed char>& v) {
        Rcpp::IntegerVector out(v.size());
        for (size_t i = 0; i < v.size(); ++i) out[i] = v[i];
        return out;
    };
    return Rcpp::List::create(
        Rcpp::Named("lon") = dd.lon,
        Rcpp::Named("lat") = dd.lat,
        Rcpp::Named("time_jd") = dd.time_jd,
        Rcpp::Named("nlon") = dd.nlon,
        Rcpp::Named("nlat") = dd.nlat,
        Rcpp::Named("ntime") = dd.ntime,
        Rcpp::Named("ndepth") = dd.ndepth,
        Rcpp::Named("depth") = dd.depth,
        Rcpp::Named("temp") = dd.temp,
        Rcpp::Named("seas") = dd.seas,
        Rcpp::Named("thresh") = dd.thresh,
        Rcpp::Named("intensity") = dd.intensity,
        Rcpp::Named("threshCriterion") = schar_to_int(dd.threshCriterion),
        Rcpp::Named("durationCriterion") = schar_to_int(dd.durationCriterion),
        Rcpp::Named("event") = schar_to_int(dd.event),
        Rcpp::Named("category") = schar_to_int(dd.category),
        Rcpp::Named("event_no") = dd.event_no
    );
}

// Single-file event detection
// [[Rcpp::export]]
void hw3_detect_events(std::string file_in,
                       std::string clim_file,
                       std::string events_file,
                       std::string var_name,
                       int minDuration = 5,
                       int minDuration2 = 5,
                       bool joinAcrossGaps = true,
                       int maxGap = 2,
                       int maxGap2 = 2,
                       bool coldSpells = false,
                       int roundRes = 4,
                       int n_threads = 1,
                       bool category = false,
                       bool southHemisphere = true,
                       std::string threshClim2_file = "",
                       std::string threshClim2_var_name = "",
                       std::string daily_file = "",
                       std::string proto_file = "") {
    hw3::DetectParams p = detect_params(minDuration, minDuration2, joinAcrossGaps,
                                        maxGap, maxGap2, coldSpells, roundRes,
                                        n_threads, category, southHemisphere);
    hw3::detect_events(file_in, clim_file, events_file, var_name, p,
                       threshClim2_file, threshClim2_var_name, daily_file, proto_file);
}

// Multi-file event detection
// [[Rcpp::export]]
void hw3_detect_events_multi(Rcpp::CharacterVector files,
                             std::string clim_file,
                             std::string events_file,
                             std::string var_name,
                             int minDuration = 5,
                             int minDuration2 = 5,
                             bool joinAcrossGaps = true,
                             int maxGap = 2,
                             int maxGap2 = 2,
                             bool coldSpells = false,
                             int roundRes = 4,
                             int n_threads = 1,
                             bool category = false,
                             bool southHemisphere = true,
                             Rcpp::Nullable<Rcpp::CharacterVector> threshClim2_files = R_NilValue,
                             std::string threshClim2_var_name = "",
                             bool skip_bad_files = false,
                             std::string daily_file = "",
                             std::string proto_file = "") {
    hw3::DetectParams p = detect_params(minDuration, minDuration2, joinAcrossGaps,
                                        maxGap, maxGap2, coldSpells, roundRes,
                                        n_threads, category, southHemisphere);
    hw3::detect_events_multi(Rcpp::as<std::vector<std::string>>(files), clim_file,
                             events_file, var_name, p, opt_str(threshClim2_files),
                             threshClim2_var_name, skip_bad_files, daily_file, proto_file);
}

// [[Rcpp::export]]
Rcpp::List hw3_read_event_nc(std::string event_file) {
    hw3::EventData ed = hw3::read_event_netcdf(event_file);
    return Rcpp::List::create(
        Rcpp::Named("lon") = ed.lon,
        Rcpp::Named("lat") = ed.lat,
        Rcpp::Named("depth") = ed.depth,
        Rcpp::Named("pixel_index") = ed.pixel_index,
        Rcpp::Named("event_no") = ed.event_no,
        Rcpp::Named("date_start") = ed.date_start,
        Rcpp::Named("date_peak") = ed.date_peak,
        Rcpp::Named("date_end") = ed.date_end,
        Rcpp::Named("duration") = ed.duration,
        Rcpp::Named("intensity_mean") = ed.intensity_mean,
        Rcpp::Named("intensity_max") = ed.intensity_max,
        Rcpp::Named("intensity_var") = ed.intensity_var,
        Rcpp::Named("intensity_cumulative") = ed.intensity_cumulative,
        Rcpp::Named("intensity_mean_relThresh") = ed.intensity_mean_relThresh,
        Rcpp::Named("intensity_max_relThresh") = ed.intensity_max_relThresh,
        Rcpp::Named("intensity_var_relThresh") = ed.intensity_var_relThresh,
        Rcpp::Named("intensity_cumulative_relThresh") = ed.intensity_cumulative_relThresh,
        Rcpp::Named("intensity_mean_abs") = ed.intensity_mean_abs,
        Rcpp::Named("intensity_max_abs") = ed.intensity_max_abs,
        Rcpp::Named("intensity_var_abs") = ed.intensity_var_abs,
        Rcpp::Named("intensity_cumulative_abs") = ed.intensity_cumulative_abs,
        Rcpp::Named("rate_onset") = ed.rate_onset,
        Rcpp::Named("rate_decline") = ed.rate_decline,
        Rcpp::Named("category") = ed.category,
        Rcpp::Named("p_moderate") = ed.p_moderate,
        Rcpp::Named("p_strong") = ed.p_strong,
        Rcpp::Named("p_severe") = ed.p_severe,
        Rcpp::Named("p_extreme") = ed.p_extreme,
        Rcpp::Named("season") = ed.season,
        Rcpp::Named("ref_date_jd") = ed.ref_date_jd,
        Rcpp::Named("nevents") = ed.nevents
    );
}

// [[Rcpp::export]]
Rcpp::DataFrame hw3_category(std::string event_file,
                             std::string clim_file,
                             std::string hemisphere = "south",
                             int roundVal = 4) {
    hw3::CategoryResult r = hw3::category(event_file, clim_file, hemisphere, roundVal);

    const char* cat_labels[] = {"", "I Moderate", "II Strong", "III Severe", "IV Extreme"};
    const char* season_labels[] = {"", "Summer", "Fall", "Winter", "Spring"};

    int n = static_cast<int>(r.event_no.size());
    Rcpp::CharacterVector category(n, NA_STRING);
    Rcpp::CharacterVector season(n);
    Rcpp::DateVector peak_date(n);
    for (int i = 0; i < n; ++i) {
        // Peak date as R Date (days since 1970-01-01); JD of 1970-01-01 = 2440588
        peak_date[i] = Rcpp::Date(r.peak_jd[i] - 2440588);
        if (r.category[i] >= 1) category[i] = cat_labels[r.category[i]];
        season[i] = season_labels[r.season[i]];
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("event_no") = r.event_no,
        Rcpp::Named("lon") = r.lon,
        Rcpp::Named("lat") = r.lat,
        Rcpp::Named("peak_date") = peak_date,
        Rcpp::Named("category") = category,
        Rcpp::Named("intensity_max") = r.intensity_max,
        Rcpp::Named("duration") = r.duration,
        Rcpp::Named("p_moderate") = r.p_moderate,
        Rcpp::Named("p_strong") = r.p_strong,
        Rcpp::Named("p_severe") = r.p_severe,
        Rcpp::Named("p_extreme") = r.p_extreme,
        Rcpp::Named("season") = season,
        Rcpp::Named("stringsAsFactors") = false
    );
}

// [[Rcpp::export]]
Rcpp::DataFrame hw3_block_average(std::string event_file) {
    hw3::BlockAverageResult r = hw3::block_average(event_file);
    return Rcpp::DataFrame::create(
        Rcpp::Named("lon") = r.lon,
        Rcpp::Named("lat") = r.lat,
        Rcpp::Named("year") = r.year,
        Rcpp::Named("count") = r.count,
        Rcpp::Named("duration_mean") = r.duration_mean,
        Rcpp::Named("duration_max") = r.duration_max,
        Rcpp::Named("intensity_mean") = r.intensity_mean,
        Rcpp::Named("intensity_max_mean") = r.intensity_max_mean,
        Rcpp::Named("intensity_max_max") = r.intensity_max_max,
        Rcpp::Named("intensity_cumulative_mean") = r.intensity_cumulative_mean,
        Rcpp::Named("total_days") = r.total_days,
        Rcpp::Named("total_icum") = r.total_icum,
        Rcpp::Named("stringsAsFactors") = false
    );
}

// [[Rcpp::export]]
Rcpp::DataFrame hw3_read_metric_summary(std::string event_file,
                                        std::string metric = "intensity_max",
                                        std::string summary = "mean") {
    hw3::MetricSummaryResult r = hw3::metric_summary(event_file, metric, summary);
    return Rcpp::DataFrame::create(
        Rcpp::Named("lon") = r.lon,
        Rcpp::Named("lat") = r.lat,
        Rcpp::Named("value") = r.value,
        Rcpp::Named("stringsAsFactors") = false
    );
}
