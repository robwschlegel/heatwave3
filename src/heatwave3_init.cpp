#include <Rcpp.h>
#include <cstdlib>
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
                                    bool coldSpells, int roundRes, int n_threads) {
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
        minDuration, maxGap, coldSpells
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
                       int n_threads = 1) {

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
                            coldSpells, roundRes, n_threads);
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
                             int n_threads = 1) {

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
                            coldSpells, roundRes, n_threads);
}
