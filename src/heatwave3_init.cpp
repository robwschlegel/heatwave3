#include <Rcpp.h>
#include <cstdlib>
#include <map>
#include "heatwave3_types.h"
#include "netcdf_io.h"
#include "climatology.h"
#include "event_detect.h"

// Parse "YYYY-MM-DD" to Julian Day Number
static int parse_date_to_jd(const std::string& s) {
    int y = std::atoi(s.substr(0, 4).c_str());
    int m = std::atoi(s.substr(5, 2).c_str());
    int d = std::atoi(s.substr(8, 2).c_str());
    return hw3::ymd_to_jd(y, m, d);
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
                        int depth = -1) {
    hw3::SubsetSpec ss;
    if (lon_range.isNotNull()) {
        Rcpp::NumericVector lr(lon_range);
        ss.lon_min = lr[0]; ss.lon_max = lr[1];
    }
    if (lat_range.isNotNull()) {
        Rcpp::NumericVector lr(lat_range);
        ss.lat_min = lr[0]; ss.lat_max = lr[1];
    }
    if (time_range.isNotNull()) {
        Rcpp::CharacterVector tr(time_range);
        ss.time_min = parse_date_to_jd(Rcpp::as<std::string>(tr[0]));
        ss.time_max = parse_date_to_jd(Rcpp::as<std::string>(tr[1]));
    }
    ss.depth_index = depth;

    if (var_name.empty()) {
        var_name = hw3::detect_sst_variable(file_in);
    }

    hw3::GridData gd = hw3::read_sst_netcdf(file_in, var_name, ss);

    return Rcpp::List::create(
        Rcpp::Named("lon") = gd.lon,
        Rcpp::Named("lat") = gd.lat,
        Rcpp::Named("nlon") = gd.nlon,
        Rcpp::Named("nlat") = gd.nlat,
        Rcpp::Named("ntime") = gd.ntime,
        Rcpp::Named("time_days") = gd.time_days,
        Rcpp::Named("sst") = gd.sst
    );
}

// Read and merge multiple daily NetCDF files
// [[Rcpp::export]]
Rcpp::List hw3_read_sst_multi(Rcpp::CharacterVector files,
                              std::string var_name,
                              Rcpp::Nullable<Rcpp::NumericVector> lon_range = R_NilValue,
                              Rcpp::Nullable<Rcpp::NumericVector> lat_range = R_NilValue,
                              int depth = -1) {
    hw3::SubsetSpec ss;
    if (lon_range.isNotNull()) {
        Rcpp::NumericVector lr(lon_range);
        ss.lon_min = lr[0]; ss.lon_max = lr[1];
    }
    if (lat_range.isNotNull()) {
        Rcpp::NumericVector lr(lat_range);
        ss.lat_min = lr[0]; ss.lat_max = lr[1];
    }
    ss.depth_index = depth;

    std::vector<std::string> file_vec;
    for (int i = 0; i < files.size(); ++i) {
        file_vec.push_back(Rcpp::as<std::string>(files[i]));
    }

    if (var_name.empty() && !file_vec.empty()) {
        var_name = hw3::detect_sst_variable(file_vec[0]);
    }

    Rcpp::Rcout << "Reading " << file_vec.size() << " files..." << std::endl;
    hw3::GridData gd = hw3::read_sst_multi_netcdf(file_vec, var_name, ss);
    Rcpp::Rcout << "Merged: " << gd.nlon << " lon x " << gd.nlat << " lat x "
                << gd.ntime << " time" << std::endl;

    return Rcpp::List::create(
        Rcpp::Named("lon") = gd.lon,
        Rcpp::Named("lat") = gd.lat,
        Rcpp::Named("nlon") = gd.nlon,
        Rcpp::Named("nlat") = gd.nlat,
        Rcpp::Named("ntime") = gd.ntime,
        Rcpp::Named("time_days") = gd.time_days,
        Rcpp::Named("sst") = gd.sst
    );
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

    int cs = parse_date_to_jd(clim_start);
    int ce = parse_date_to_jd(clim_end);

    hw3::PixelClimResult cr = hw3::compute_pixel_climatology(
        t.data(), jd.data(), n, cs, ce,
        windowHalfWidth, pctile, smoothPercentile, smoothPercentileWidth,
        false, maxPadLength, detrend
    );

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

// Shared helper: compute and write climatology from a pre-loaded GridData
static void compute_and_write_clim(hw3::GridData& gd,
                                   const std::string& file_out,
                                   const std::string& source_label,
                                   const std::string& cp0, const std::string& cp1,
                                   int maxPadLength, int windowHalfWidth,
                                   double pctile, bool smoothPercentile,
                                   int smoothPercentileWidth,
                                   bool compute_var, int roundClm, int n_threads,
                                   bool detrend = false) {
    int clim_start_jd = parse_date_to_jd(cp0);
    int clim_end_jd = parse_date_to_jd(cp1);
    int npixels = gd.nlon * gd.nlat;
    Rcpp::Rcout << "Grid: " << gd.nlon << " lon x " << gd.nlat << " lat x "
                << gd.ntime << " time = " << npixels << " pixels" << std::endl;

    // Allocate output
    size_t grid_size = static_cast<size_t>(npixels) * 366;
    std::vector<double> seas(grid_size);
    std::vector<double> thresh(grid_size);
    std::vector<double> var_clim;
    if (compute_var) var_clim.resize(grid_size);

    Rcpp::Rcout << "Computing climatology with " << n_threads << " thread(s)..." << std::endl;

    hw3::compute_climatology_grid(
        gd.sst.data(), gd.time_days.data(),
        npixels, gd.ntime,
        clim_start_jd, clim_end_jd,
        windowHalfWidth, pctile,
        smoothPercentile, smoothPercentileWidth,
        compute_var, maxPadLength, detrend, n_threads,
        seas.data(), thresh.data(),
        compute_var ? var_clim.data() : nullptr
    );

    // Round
    if (roundClm > 0) {
        double factor = std::pow(10.0, roundClm);
        for (size_t i = 0; i < grid_size; ++i) {
            if (!hw3::is_na(seas[i]))
                seas[i] = std::round(seas[i] * factor) / factor;
            if (!hw3::is_na(thresh[i]))
                thresh[i] = std::round(thresh[i] * factor) / factor;
            if (compute_var && !hw3::is_na(var_clim[i]))
                var_clim[i] = std::round(var_clim[i] * factor) / factor;
        }
    }

    Rcpp::Rcout << "Writing climatology to " << file_out << "..." << std::endl;

    hw3::write_clim_netcdf(
        file_out, gd.lon, gd.lat, gd.nlon, gd.nlat,
        seas, thresh, var_clim, compute_var,
        source_label, cp0, cp1,
        pctile, windowHalfWidth, smoothPercentileWidth
    );

    Rcpp::Rcout << "Done." << std::endl;
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
                      bool detrend = false) {

    hw3::SubsetSpec ss;
    if (lon_range.isNotNull()) {
        Rcpp::NumericVector lr(lon_range);
        ss.lon_min = lr[0]; ss.lon_max = lr[1];
    }
    if (lat_range.isNotNull()) {
        Rcpp::NumericVector lr(lat_range);
        ss.lat_min = lr[0]; ss.lat_max = lr[1];
    }
    if (time_range.isNotNull()) {
        Rcpp::CharacterVector tr(time_range);
        ss.time_min = parse_date_to_jd(Rcpp::as<std::string>(tr[0]));
        ss.time_max = parse_date_to_jd(Rcpp::as<std::string>(tr[1]));
    }
    ss.depth_index = depth;
    if (var_name.empty()) var_name = hw3::detect_sst_variable(file_in);

    Rcpp::Rcout << "Reading SST data from " << file_in << "..." << std::endl;
    hw3::GridData gd = hw3::read_sst_netcdf(file_in, var_name, ss);

    std::string cp0 = Rcpp::as<std::string>(climatologyPeriod[0]);
    std::string cp1 = Rcpp::as<std::string>(climatologyPeriod[1]);
    compute_and_write_clim(gd, file_out, file_in, cp0, cp1,
                           maxPadLength, windowHalfWidth, pctile,
                           smoothPercentile, smoothPercentileWidth,
                           compute_var, roundClm, n_threads, detrend);
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
                            bool detrend = false) {

    hw3::SubsetSpec ss;
    if (lon_range.isNotNull()) {
        Rcpp::NumericVector lr(lon_range);
        ss.lon_min = lr[0]; ss.lon_max = lr[1];
    }
    if (lat_range.isNotNull()) {
        Rcpp::NumericVector lr(lat_range);
        ss.lat_min = lr[0]; ss.lat_max = lr[1];
    }
    ss.depth_index = depth;

    std::vector<std::string> file_vec;
    for (int i = 0; i < files.size(); ++i)
        file_vec.push_back(Rcpp::as<std::string>(files[i]));

    if (var_name.empty() && !file_vec.empty())
        var_name = hw3::detect_sst_variable(file_vec[0]);

    Rcpp::Rcout << "Reading " << file_vec.size() << " daily files..." << std::endl;
    hw3::GridData gd = hw3::read_sst_multi_netcdf(file_vec, var_name, ss);

    std::string cp0 = Rcpp::as<std::string>(climatologyPeriod[0]);
    std::string cp1 = Rcpp::as<std::string>(climatologyPeriod[1]);
    std::string source_label = std::to_string(file_vec.size()) + " daily files";
    compute_and_write_clim(gd, file_out, source_label, cp0, cp1,
                           maxPadLength, windowHalfWidth, pctile,
                           smoothPercentile, smoothPercentileWidth,
                           compute_var, roundClm, n_threads, detrend);
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
    hw3::SubsetSpec ss;
    if (lon_range.isNotNull()) {
        Rcpp::NumericVector lr(lon_range);
        ss.lon_min = lr[0]; ss.lon_max = lr[1];
    }
    if (lat_range.isNotNull()) {
        Rcpp::NumericVector lr(lat_range);
        ss.lat_min = lr[0]; ss.lat_max = lr[1];
    }
    if (time_range.isNotNull()) {
        Rcpp::CharacterVector tr(time_range);
        ss.time_min = parse_date_to_jd(Rcpp::as<std::string>(tr[0]));
        ss.time_max = parse_date_to_jd(Rcpp::as<std::string>(tr[1]));
    }
    ss.depth_index = depth;

    if (var_name.empty()) var_name = hw3::detect_sst_variable(file_in);
    hw3::GridData gd = hw3::read_sst_netcdf(file_in, var_name, ss);

    int npixels = gd.nlon * gd.nlat;
    size_t grid_size = static_cast<size_t>(npixels) * 366;
    std::vector<double> seas(grid_size, threshold);
    std::vector<double> thresh(grid_size, threshold);

    hw3::write_clim_netcdf(clim_file, gd.lon, gd.lat, gd.nlon, gd.nlat,
                           seas, thresh, {}, false,
                           file_in, "static", "static",
                           90.0, 5, 31);
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
        Rcpp::Named("seas") = cd.seas,
        Rcpp::Named("thresh") = cd.thresh
    );
}

// Shared helper: detect events from pre-loaded GridData + ClimData
static void detect_and_write_events(hw3::GridData& gd, hw3::ClimData& cd,
                                    const std::string& file_out,
                                    const std::string& source_label,
                                    const std::string& clim_file,
                                    int minDuration, int minDuration2,
                                    bool joinAcrossGaps, int maxGap, int maxGap2,
                                    bool coldSpells, int roundRes, int n_threads,
                                    bool category, bool southHemisphere) {
    if (gd.nlon != cd.nlon || gd.nlat != cd.nlat) {
        Rcpp::stop("Grid mismatch: SST is %d x %d but climatology is %d x %d",
                   gd.nlon, gd.nlat, cd.nlon, cd.nlat);
    }

    int npixels = gd.nlon * gd.nlat;
    Rcpp::Rcout << "Grid: " << gd.nlon << " lon x " << gd.nlat << " lat x "
                << gd.ntime << " time = " << npixels << " pixels" << std::endl;

    std::vector<double> event_lon, event_lat;
    std::vector<int> pixel_idx;
    std::vector<hw3::EventResult> all_events;
    std::vector<int> ds, dp, de;

    Rcpp::Rcout << "Detecting events with " << n_threads << " thread(s)..." << std::endl;

    hw3::detect_events_grid(
        gd.sst.data(), cd.seas.data(), cd.thresh.data(),
        gd.time_days.data(), npixels, gd.ntime,
        minDuration, minDuration2,
        joinAcrossGaps, maxGap, maxGap2,
        coldSpells, roundRes, n_threads,
        category, southHemisphere,
        event_lon, event_lat, pixel_idx, all_events,
        ds, dp, de,
        gd.lon, gd.lat, gd.nlon, gd.nlat
    );

    Rcpp::Rcout << "Found " << all_events.size() << " events across "
                << npixels << " pixels" << std::endl;

    if (all_events.empty()) {
        Rcpp::Rcout << "No events detected. Skipping output." << std::endl;
        return;
    }

    int ref_jd = gd.time_days[0];
    std::vector<int> ds_rel(ds.size()), dp_rel(dp.size()), de_rel(de.size());
    for (size_t i = 0; i < ds.size(); ++i) {
        ds_rel[i] = ds[i] - ref_jd;
        dp_rel[i] = dp[i] - ref_jd;
        de_rel[i] = de[i] - ref_jd;
    }

    Rcpp::Rcout << "Writing events to " << file_out << "..." << std::endl;

    hw3::write_event_netcdf(
        file_out, event_lon, event_lat, pixel_idx, all_events,
        ds_rel, dp_rel, de_rel, ref_jd,
        source_label, clim_file,
        minDuration, maxGap, coldSpells,
        southHemisphere
    );

    Rcpp::Rcout << "Done." << std::endl;
}

// Single-file event detection
// [[Rcpp::export]]
void hw3_detect_events(std::string file_in,
                       std::string clim_file,
                       std::string file_out,
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
                       bool southHemisphere = true) {

    Rcpp::Rcout << "Reading climatology from " << clim_file << "..." << std::endl;
    hw3::ClimData cd = hw3::read_clim_netcdf(clim_file);

    hw3::SubsetSpec ss;
    ss.lon_min = *std::min_element(cd.lon.begin(), cd.lon.end());
    ss.lon_max = *std::max_element(cd.lon.begin(), cd.lon.end());
    ss.lat_min = *std::min_element(cd.lat.begin(), cd.lat.end());
    ss.lat_max = *std::max_element(cd.lat.begin(), cd.lat.end());

    if (var_name.empty()) var_name = hw3::detect_sst_variable(file_in);

    Rcpp::Rcout << "Reading SST data from " << file_in << "..." << std::endl;
    hw3::GridData gd = hw3::read_sst_netcdf(file_in, var_name, ss);

    detect_and_write_events(gd, cd, file_out, file_in, clim_file,
                            minDuration, minDuration2,
                            joinAcrossGaps, maxGap, maxGap2,
                            coldSpells, roundRes, n_threads,
                            category, southHemisphere);
}

// Multi-file event detection
// [[Rcpp::export]]
void hw3_detect_events_multi(Rcpp::CharacterVector files,
                             std::string clim_file,
                             std::string file_out,
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
                             bool southHemisphere = true) {

    Rcpp::Rcout << "Reading climatology from " << clim_file << "..." << std::endl;
    hw3::ClimData cd = hw3::read_clim_netcdf(clim_file);

    hw3::SubsetSpec ss;
    ss.lon_min = *std::min_element(cd.lon.begin(), cd.lon.end());
    ss.lon_max = *std::max_element(cd.lon.begin(), cd.lon.end());
    ss.lat_min = *std::min_element(cd.lat.begin(), cd.lat.end());
    ss.lat_max = *std::max_element(cd.lat.begin(), cd.lat.end());

    std::vector<std::string> file_vec;
    for (int i = 0; i < files.size(); ++i)
        file_vec.push_back(Rcpp::as<std::string>(files[i]));

    if (var_name.empty() && !file_vec.empty())
        var_name = hw3::detect_sst_variable(file_vec[0]);

    Rcpp::Rcout << "Reading " << file_vec.size() << " daily files..." << std::endl;
    hw3::GridData gd = hw3::read_sst_multi_netcdf(file_vec, var_name, ss);

    std::string source_label = std::to_string(file_vec.size()) + " daily files";
    detect_and_write_events(gd, cd, file_out, source_label, clim_file,
                            minDuration, minDuration2,
                            joinAcrossGaps, maxGap, maxGap2,
                            coldSpells, roundRes, n_threads,
                            category, southHemisphere);
}

// [[Rcpp::export]]
Rcpp::List hw3_read_event_nc(std::string event_file) {
    hw3::EventData ed = hw3::read_event_netcdf(event_file);
    return Rcpp::List::create(
        Rcpp::Named("lon") = ed.lon,
        Rcpp::Named("lat") = ed.lat,
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
    hw3::EventData ed = hw3::read_event_netcdf(event_file);
    bool southHemisphere = (hemisphere == "south");

    int n = ed.nevents;

    Rcpp::CharacterVector category(n, NA_STRING);
    Rcpp::NumericVector p_moderate(n, 0.0);
    Rcpp::NumericVector p_strong(n, 0.0);
    Rcpp::NumericVector p_severe(n, 0.0);
    Rcpp::NumericVector p_extreme(n, 0.0);
    Rcpp::CharacterVector season(n);
    Rcpp::DateVector peak_date(n);

    double round_mult = std::pow(10.0, roundVal);

    // Check if event file already has category data
    bool has_precomputed = !ed.category.empty();

    if (has_precomputed) {
        // Use pre-computed category data from event file
        const char* cat_labels[] = {"", "I Moderate", "II Strong", "III Severe", "IV Extreme"};
        const char* season_labels[] = {"", "Summer", "Fall", "Winter", "Spring"};

        for (int i = 0; i < n; ++i) {
            int peak_jd = ed.ref_date_jd + ed.date_peak[i];
            peak_date[i] = Rcpp::Date(peak_jd - 2440588);

            int cat = ed.category[i];
            if (cat >= 1 && cat <= 4) {
                category[i] = cat_labels[cat];
            }
            p_moderate[i] = ed.p_moderate[i];
            p_strong[i] = ed.p_strong[i];
            p_severe[i] = ed.p_severe[i];
            p_extreme[i] = ed.p_extreme[i];

            int sea = ed.season[i];
            if (sea >= 1 && sea <= 4) {
                season[i] = season_labels[sea];
            }
        }
    } else {
        if (clim_file.empty())
            Rcpp::stop("Event file has no pre-computed categories; clim_file is required.");
        hw3::ClimData cd = hw3::read_clim_netcdf(clim_file);

        for (int i = 0; i < n; ++i) {
            // Convert relative date to absolute JD
            int peak_jd = ed.ref_date_jd + ed.date_peak[i];

            // Peak date as R Date (days since 1970-01-01)
            // JD of 1970-01-01 = 2440588
            peak_date[i] = Rcpp::Date(peak_jd - 2440588);

            // DOY in 1..366 (heatwaveR convention: non-leap skips 60)
            int doy = hw3::jd_to_doy_366(peak_jd);

            // Find nearest pixel in climatology grid
            double ev_lon = ed.lon[i];
            double ev_lat = ed.lat[i];
            int lon_idx = 0, lat_idx = 0;
            double best_lon = std::abs(cd.lon[0] - ev_lon);
            for (int j = 1; j < cd.nlon; ++j) {
                double d = std::abs(cd.lon[j] - ev_lon);
                if (d < best_lon) { best_lon = d; lon_idx = j; }
            }
            double best_lat = std::abs(cd.lat[0] - ev_lat);
            for (int j = 1; j < cd.nlat; ++j) {
                double d = std::abs(cd.lat[j] - ev_lat);
                if (d < best_lat) { best_lat = d; lat_idx = j; }
            }

            // Climatology layout: [pixel * 366 + doy-1] where pixel = ilon * nlat + ilat
            int px = lon_idx * cd.nlat + lat_idx;
            double s = cd.seas[px * 366 + (doy - 1)];
            double th = cd.thresh[px * 366 + (doy - 1)];

            // Category width: |seas - thresh|, always positive
            double diff = std::abs(s - th);
            if (hw3::is_na(s) || hw3::is_na(th) || diff <= 0.0) continue;

            double ix = ed.intensity_max[i];

            if (ix >= 4.0 * diff) {
                category[i] = "IV Extreme";
            } else if (ix >= 3.0 * diff) {
                category[i] = "III Severe";
            } else if (ix >= 2.0 * diff) {
                category[i] = "II Strong";
            } else {
                category[i] = "I Moderate";
            }

            auto rnd = [&](double v) {
                return std::round(v * round_mult) / round_mult;
            };
            p_moderate[i] = rnd(std::min(1.0, std::max(0.0, ix / diff)));
            p_strong[i]   = rnd(std::min(1.0, std::max(0.0, (ix - diff) / diff)));
            p_severe[i]   = rnd(std::min(1.0, std::max(0.0, (ix - 2.0 * diff) / diff)));
            p_extreme[i]  = rnd(std::min(1.0, std::max(0.0, (ix - 3.0 * diff) / diff)));

            // Season from peak month
            int y, m, d;
            hw3::jd_to_ymd(peak_jd, y, m, d);
            if (southHemisphere) {
                if (m == 12 || m <= 2)      season[i] = "Summer";
                else if (m >= 3 && m <= 5)  season[i] = "Fall";
                else if (m >= 6 && m <= 8)  season[i] = "Winter";
                else                        season[i] = "Spring";
            } else {
                if (m == 12 || m <= 2)      season[i] = "Winter";
                else if (m >= 3 && m <= 5)  season[i] = "Spring";
                else if (m >= 6 && m <= 8)  season[i] = "Summer";
                else                        season[i] = "Fall";
            }
        }
    }

    // Round lon/lat for display
    Rcpp::NumericVector out_lon(n), out_lat(n), out_imax(n);
    for (int i = 0; i < n; ++i) {
        out_lon[i] = std::round(ed.lon[i] * round_mult) / round_mult;
        out_lat[i] = std::round(ed.lat[i] * round_mult) / round_mult;
        out_imax[i] = std::round(ed.intensity_max[i] * round_mult) / round_mult;
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("event_no") = ed.event_no,
        Rcpp::Named("lon") = out_lon,
        Rcpp::Named("lat") = out_lat,
        Rcpp::Named("peak_date") = peak_date,
        Rcpp::Named("category") = category,
        Rcpp::Named("intensity_max") = out_imax,
        Rcpp::Named("duration") = ed.duration,
        Rcpp::Named("p_moderate") = p_moderate,
        Rcpp::Named("p_strong") = p_strong,
        Rcpp::Named("p_severe") = p_severe,
        Rcpp::Named("p_extreme") = p_extreme,
        Rcpp::Named("season") = season,
        Rcpp::Named("stringsAsFactors") = false
    );
}

// [[Rcpp::export]]
Rcpp::DataFrame hw3_block_average(std::string event_file) {
    hw3::EventData ed = hw3::read_event_netcdf(event_file);

    int n = ed.nevents;
    if (n == 0) {
        return Rcpp::DataFrame::create(
            Rcpp::Named("lon") = Rcpp::NumericVector(0),
            Rcpp::Named("lat") = Rcpp::NumericVector(0),
            Rcpp::Named("year") = Rcpp::IntegerVector(0),
            Rcpp::Named("count") = Rcpp::IntegerVector(0),
            Rcpp::Named("duration_mean") = Rcpp::NumericVector(0),
            Rcpp::Named("duration_max") = Rcpp::IntegerVector(0),
            Rcpp::Named("intensity_mean") = Rcpp::NumericVector(0),
            Rcpp::Named("intensity_max_mean") = Rcpp::NumericVector(0),
            Rcpp::Named("intensity_max_max") = Rcpp::NumericVector(0),
            Rcpp::Named("intensity_cumulative_mean") = Rcpp::NumericVector(0),
            Rcpp::Named("total_days") = Rcpp::IntegerVector(0),
            Rcpp::Named("total_icum") = Rcpp::NumericVector(0),
            Rcpp::Named("stringsAsFactors") = false
        );
    }

    // Compute year from date_start for each event
    std::vector<int> year(n);
    for (int i = 0; i < n; ++i) {
        int jd = ed.ref_date_jd + ed.date_start[i];
        int y, m, d;
        hw3::jd_to_ymd(jd, y, m, d);
        year[i] = y;
    }

    // Group by (pixel_index, year) using a map
    struct BlockAccum {
        double lon, lat;
        int year;
        int count = 0;
        int duration_sum = 0;
        int duration_max = 0;
        double imean_sum = 0.0;
        double imax_sum = 0.0;
        double imax_max = -1e30;
        double icum_sum = 0.0;
    };
    std::map<std::pair<int,int>, BlockAccum> groups;

    for (int i = 0; i < n; ++i) {
        auto key = std::make_pair(ed.pixel_index[i], year[i]);
        auto& g = groups[key];
        if (g.count == 0) {
            g.lon = ed.lon[i];
            g.lat = ed.lat[i];
            g.year = year[i];
        }
        g.count++;
        g.duration_sum += ed.duration[i];
        if (ed.duration[i] > g.duration_max) g.duration_max = ed.duration[i];
        g.imean_sum += ed.intensity_mean[i];
        g.imax_sum += ed.intensity_max[i];
        if (ed.intensity_max[i] > g.imax_max) g.imax_max = ed.intensity_max[i];
        g.icum_sum += ed.intensity_cumulative[i];
    }

    int ng = static_cast<int>(groups.size());
    Rcpp::NumericVector o_lon(ng), o_lat(ng), o_dur_mean(ng);
    Rcpp::NumericVector o_imean(ng), o_imax_mean(ng), o_imax_max(ng);
    Rcpp::NumericVector o_icum_mean(ng), o_total_icum(ng);
    Rcpp::IntegerVector o_year(ng), o_count(ng), o_dur_max(ng), o_total_days(ng);

    int idx = 0;
    for (auto& [key, g] : groups) {
        o_lon[idx] = g.lon;
        o_lat[idx] = g.lat;
        o_year[idx] = g.year;
        o_count[idx] = g.count;
        o_dur_mean[idx] = static_cast<double>(g.duration_sum) / g.count;
        o_dur_max[idx] = g.duration_max;
        o_imean[idx] = g.imean_sum / g.count;
        o_imax_mean[idx] = g.imax_sum / g.count;
        o_imax_max[idx] = g.imax_max;
        o_icum_mean[idx] = g.icum_sum / g.count;
        o_total_days[idx] = g.duration_sum;
        o_total_icum[idx] = g.icum_sum;
        idx++;
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("lon") = o_lon,
        Rcpp::Named("lat") = o_lat,
        Rcpp::Named("year") = o_year,
        Rcpp::Named("count") = o_count,
        Rcpp::Named("duration_mean") = o_dur_mean,
        Rcpp::Named("duration_max") = o_dur_max,
        Rcpp::Named("intensity_mean") = o_imean,
        Rcpp::Named("intensity_max_mean") = o_imax_mean,
        Rcpp::Named("intensity_max_max") = o_imax_max,
        Rcpp::Named("intensity_cumulative_mean") = o_icum_mean,
        Rcpp::Named("total_days") = o_total_days,
        Rcpp::Named("total_icum") = o_total_icum,
        Rcpp::Named("stringsAsFactors") = false
    );
}

// [[Rcpp::export]]
Rcpp::DataFrame hw3_read_metric_summary(std::string event_file,
                                        std::string metric = "intensity_max",
                                        std::string summary = "mean") {
    hw3::EventData ed = hw3::read_event_netcdf(event_file);
    int n = ed.nevents;
    if (n == 0) {
        return Rcpp::DataFrame::create(
            Rcpp::Named("lon") = Rcpp::NumericVector(0),
            Rcpp::Named("lat") = Rcpp::NumericVector(0),
            Rcpp::Named("value") = Rcpp::NumericVector(0),
            Rcpp::Named("stringsAsFactors") = false
        );
    }

    // Select the metric vector
    const std::vector<double>* vals = nullptr;
    const std::vector<int>* ivals = nullptr;
    if (metric == "intensity_max") vals = &ed.intensity_max;
    else if (metric == "intensity_mean") vals = &ed.intensity_mean;
    else if (metric == "intensity_var") vals = &ed.intensity_var;
    else if (metric == "intensity_cumulative") vals = &ed.intensity_cumulative;
    else if (metric == "intensity_mean_relThresh") vals = &ed.intensity_mean_relThresh;
    else if (metric == "intensity_max_relThresh") vals = &ed.intensity_max_relThresh;
    else if (metric == "intensity_var_relThresh") vals = &ed.intensity_var_relThresh;
    else if (metric == "intensity_cumulative_relThresh") vals = &ed.intensity_cumulative_relThresh;
    else if (metric == "intensity_mean_abs") vals = &ed.intensity_mean_abs;
    else if (metric == "intensity_max_abs") vals = &ed.intensity_max_abs;
    else if (metric == "intensity_var_abs") vals = &ed.intensity_var_abs;
    else if (metric == "intensity_cumulative_abs") vals = &ed.intensity_cumulative_abs;
    else if (metric == "rate_onset") vals = &ed.rate_onset;
    else if (metric == "rate_decline") vals = &ed.rate_decline;
    else if (metric == "duration") ivals = &ed.duration;
    else Rcpp::stop("Unknown metric: %s", metric.c_str());

    // Group by pixel_index and aggregate
    struct Accum { double lon, lat, sum; double mx, mn; int count; };
    std::map<int, Accum> groups;
    for (int i = 0; i < n; ++i) {
        double v = vals ? (*vals)[i] : static_cast<double>((*ivals)[i]);
        auto& g = groups[ed.pixel_index[i]];
        if (g.count == 0) {
            g.lon = ed.lon[i]; g.lat = ed.lat[i];
            g.sum = v; g.mx = v; g.mn = v; g.count = 1;
        } else {
            g.sum += v;
            if (v > g.mx) g.mx = v;
            if (v < g.mn) g.mn = v;
            g.count++;
        }
    }

    int ng = static_cast<int>(groups.size());
    Rcpp::NumericVector o_lon(ng), o_lat(ng), o_val(ng);
    int idx = 0;
    for (auto& [key, g] : groups) {
        o_lon[idx] = g.lon;
        o_lat[idx] = g.lat;
        if (summary == "mean") o_val[idx] = g.sum / g.count;
        else if (summary == "max") o_val[idx] = g.mx;
        else if (summary == "min") o_val[idx] = g.mn;
        else if (summary == "sum") o_val[idx] = g.sum;
        else if (summary == "count") o_val[idx] = static_cast<double>(g.count);
        else Rcpp::stop("Unknown summary: %s", summary.c_str());
        idx++;
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("lon") = o_lon,
        Rcpp::Named("lat") = o_lat,
        Rcpp::Named("value") = o_val,
        Rcpp::Named("stringsAsFactors") = false
    );
}
