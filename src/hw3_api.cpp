#include "hw3_api.h"
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <stdexcept>
#include "hw3_log.h"
#include "hw3_round.h"
#include "hw3_omp.h"
#include "event_detect.h"

namespace hw3 {

// printf-style message for std::runtime_error.
static std::string strf(const char* fmt, ...) {
    char buf[512];
    va_list args;
    va_start(args, fmt);
    std::vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    return std::string(buf);
}

// Writer options for one call; the zlib level must be 0 (none) to 9
static WriteOpts write_opts(int compression, int n_threads) {
    if (compression < 0 || compression > 9)
        throw std::runtime_error(strf("compression must be an integer from 0 to 9, not %d",
                                      compression));
    WriteOpts wo;
    wo.compression = compression;
    wo.n_threads = n_threads;
    return wo;
}

int parse_date_to_jd(const std::string& s) {
    int y = std::atoi(s.substr(0, 4).c_str());
    int m = std::atoi(s.substr(5, 2).c_str());
    int d = std::atoi(s.substr(8, 2).c_str());
    return ymd_to_jd(y, m, d);
}

SubsetSpec make_subset(const std::vector<double>& lon_range,
                       const std::vector<double>& lat_range,
                       const std::vector<std::string>& time_range,
                       int depth,
                       const std::vector<double>& depth_range) {
    SubsetSpec ss;
    if (!lon_range.empty()) {
        ss.lon_min = lon_range[0]; ss.lon_max = lon_range[1];
    }
    if (!lat_range.empty()) {
        ss.lat_min = lat_range[0]; ss.lat_max = lat_range[1];
    }
    if (!time_range.empty()) {
        ss.time_min = parse_date_to_jd(time_range[0]);
        ss.time_max = parse_date_to_jd(time_range[1]);
    }
    // A legacy single-level `depth` index and a `depth_range` are mutually
    // exclusive.
    if (depth >= 0 && !depth_range.empty()) {
        throw std::runtime_error("Specify only one of 'depth' or 'depth_range', not both.");
    }
    if (!depth_range.empty()) {
        ss.depth_min = depth_range[0]; ss.depth_max = depth_range[1];
    } else {
        ss.depth_index = depth;
    }
    return ss;
}

// ---- SST reading ----

GridData read_sst(const std::string& file_in, std::string var_name,
                  const SubsetSpec& ss) {
    if (var_name.empty()) {
        var_name = detect_sst_variable(file_in);
    }
    return read_sst_netcdf(file_in, var_name, ss);
}

GridData read_sst_multi(const std::vector<std::string>& files,
                        std::string var_name, const SubsetSpec& ss,
                        bool skip_bad_files) {
    if (var_name.empty() && !skip_bad_files && !files.empty()) {
        var_name = detect_sst_variable(files[0]);
    }

    out() << "Reading " << files.size() << " files..." << std::endl;
    GridData gd = read_sst_multi_netcdf(files, var_name, ss, skip_bad_files);
    out() << "Merged: " << gd.nlon << " lon x " << gd.nlat << " lat x "
          << gd.ntime << " time" << std::endl;
    return gd;
}

// ---- Climatology ----

PixelClimResult pixel_clim(const std::vector<double>& temp,
                           const std::vector<int>& time_jd,
                           const std::string& clim_start,
                           const std::string& clim_end,
                           const ClimParams& p) {
    int cs = parse_date_to_jd(clim_start);
    int ce = parse_date_to_jd(clim_end);

    return compute_pixel_climatology(
        temp.data(), time_jd.data(), static_cast<int>(temp.size()), cs, ce,
        p.windowHalfWidth, p.pctile, p.smoothPercentile, p.smoothPercentileWidth,
        p.compute_var, p.maxPadLength, p.detrend
    );
}

// Shared helper: compute and write climatology from a pre-loaded GridData
static void compute_and_write_clim(GridData& gd,
                                   const std::string& file_out,
                                   const std::string& source_label,
                                   const std::string& cp0, const std::string& cp1,
                                   const ClimParams& p) {
    WriteOpts wo = write_opts(p.compression, p.n_threads);
    int clim_start_jd = parse_date_to_jd(cp0);
    int clim_end_jd = parse_date_to_jd(cp1);
    int npixels = gd.nlon * gd.nlat * gd.ndepth;
    out() << "Grid: " << gd.nlon << " lon x " << gd.nlat << " lat";
    if (gd.ndepth > 1) out() << " x " << gd.ndepth << " depth";
    out() << " x " << gd.ntime << " time = " << npixels << " pixels" << std::endl;

    // Allocate output
    size_t grid_size = static_cast<size_t>(npixels) * 366;
    std::vector<double> seas(grid_size);
    std::vector<double> thresh(grid_size);
    std::vector<double> var_clim;
    if (p.compute_var) var_clim.resize(grid_size);

    out() << "Computing climatology with "
          << (p.n_threads > 0 ? p.n_threads : get_threads(npixels)) << " thread(s)..." << std::endl;

    compute_climatology_grid(
        gd.sst.data(), gd.time_days.data(),
        npixels, gd.ntime,
        clim_start_jd, clim_end_jd,
        p.windowHalfWidth, p.pctile,
        p.smoothPercentile, p.smoothPercentileWidth,
        p.compute_var, p.maxPadLength, p.detrend, p.n_threads,
        seas.data(), thresh.data(),
        p.compute_var ? var_clim.data() : nullptr
    );

    // Round, with R's round() as heatwaveR's ts2clm(roundClm = )
    if (p.roundClm > 0) {
        for (size_t i = 0; i < grid_size; ++i) {
            seas[i] = r_round(seas[i], p.roundClm);
            thresh[i] = r_round(thresh[i], p.roundClm);
            if (p.compute_var) var_clim[i] = r_round(var_clim[i], p.roundClm);
        }
    }

    out() << "Writing climatology to " << file_out << "..." << std::endl;

    write_clim_netcdf(
        file_out, gd.lon, gd.lat, gd.nlon, gd.nlat,
        seas, thresh, var_clim, p.compute_var,
        source_label, cp0, cp1,
        p.pctile, p.windowHalfWidth, p.smoothPercentileWidth,
        gd.temp_units, gd.ndepth, gd.depth, wo
    );

    out() << "Done." << std::endl;
}

void compute_clim(const std::string& file_in, const std::string& file_out,
                  std::string var_name,
                  const std::string& clim_start, const std::string& clim_end,
                  const SubsetSpec& ss, const ClimParams& p) {
    if (var_name.empty()) var_name = detect_sst_variable(file_in);

    out() << "Reading SST data from " << file_in << "..." << std::endl;
    GridData gd = read_sst_netcdf(file_in, var_name, ss);

    compute_and_write_clim(gd, file_out, file_in, clim_start, clim_end, p);
}

void compute_clim_multi(const std::vector<std::string>& files,
                        const std::string& file_out, std::string var_name,
                        const std::string& clim_start, const std::string& clim_end,
                        const SubsetSpec& ss, const ClimParams& p,
                        bool skip_bad_files) {
    if (var_name.empty() && !skip_bad_files && !files.empty())
        var_name = detect_sst_variable(files[0]);

    out() << "Reading " << files.size() << " daily files..." << std::endl;
    GridData gd = read_sst_multi_netcdf(files, var_name, ss, skip_bad_files);

    std::string source_label = std::to_string(files.size()) + " daily files";
    compute_and_write_clim(gd, file_out, source_label, clim_start, clim_end, p);
}

void write_const_clim(const std::string& file_in, const std::string& clim_file,
                      std::string var_name, double threshold,
                      const SubsetSpec& ss) {
    if (var_name.empty()) var_name = detect_sst_variable(file_in);
    GridData gd = read_sst_netcdf(file_in, var_name, ss);

    int npixels = gd.nlon * gd.nlat;
    size_t grid_size = static_cast<size_t>(npixels) * 366;
    std::vector<double> seas(grid_size, threshold);
    std::vector<double> thresh(grid_size, threshold);

    write_clim_netcdf(clim_file, gd.lon, gd.lat, gd.nlon, gd.nlat,
                      seas, thresh, {}, false,
                      file_in, "static", "static",
                      90.0, 5, 31, gd.temp_units);
}

// ---- Event detection ----

// SubsetSpec spanning the climatology grid extent, plus an optional [jd0,jd1]
// time window (Julian Days) for the SST read. Depth-resolved climatology
// (ts2clm3(depth_range=...)) or a legacy single-level squeeze (ts2clm3(depth=k),
// recorded as a one-element cd.depth): match the SST read to the same depth(s)
// so detection stays consistent with the climatology it's compared against,
// with no separate depth argument needed.
static SubsetSpec clim_extent_subset(const ClimData& cd,
                                     const std::vector<int>& t_jd_range = {}) {
    SubsetSpec ss;
    ss.lon_min = *std::min_element(cd.lon.begin(), cd.lon.end());
    ss.lon_max = *std::max_element(cd.lon.begin(), cd.lon.end());
    ss.lat_min = *std::min_element(cd.lat.begin(), cd.lat.end());
    ss.lat_max = *std::max_element(cd.lat.begin(), cd.lat.end());
    if (!cd.depth.empty()) {
        ss.depth_min = *std::min_element(cd.depth.begin(), cd.depth.end());
        ss.depth_max = *std::max_element(cd.depth.begin(), cd.depth.end());
    }
    if (t_jd_range.size() == 2) {
        ss.time_min = std::min(t_jd_range[0], t_jd_range[1]);
        ss.time_max = std::max(t_jd_range[0], t_jd_range[1]);
    }
    return ss;
}

// Shared helper: detect events from pre-loaded GridData + ClimData
// Output is selected by which paths are non-empty:
//   events_file : event-metrics table (ragged NetCDF)
//   daily_file  : per-day series with intensity + daily category
//   proto_file  : per-day proto series with threshCriterion + durationCriterion
// Any combination is valid; an empty string skips that product. Detection
// always runs (the daily/proto series need event membership).
static void detect_and_write_events(GridData& gd, ClimData& cd,
                                    const GridData* thresh2_gd,
                                    const std::string& events_file,
                                    const std::string& source_label,
                                    const std::string& clim_file,
                                    const DetectParams& p,
                                    const std::string& daily_file,
                                    const std::string& proto_file) {
    if (gd.nlon != cd.nlon || gd.nlat != cd.nlat || gd.ndepth != cd.ndepth) {
        throw std::runtime_error(strf(
            "Grid mismatch: SST is %d x %d x %d depth but climatology is %d x %d x %d depth",
            gd.nlon, gd.nlat, gd.ndepth, cd.nlon, cd.nlat, cd.ndepth));
    }
    if (thresh2_gd != nullptr) {
        if (gd.nlon != thresh2_gd->nlon || gd.nlat != thresh2_gd->nlat ||
            gd.ntime != thresh2_gd->ntime) {
            throw std::runtime_error(strf(
                "Grid mismatch: SST is %d x %d x %d but threshClim2 is %d x %d x %d",
                gd.nlon, gd.nlat, gd.ntime,
                thresh2_gd->nlon, thresh2_gd->nlat, thresh2_gd->ntime));
        }
        for (int t = 0; t < gd.ntime; ++t) {
            if (gd.time_days[t] != thresh2_gd->time_days[t]) {
                throw std::runtime_error(strf(
                    "Time mismatch between SST and threshClim2 at time index %d", t + 1));
            }
        }
    }

    WriteOpts wo = write_opts(p.compression, p.n_threads);
    int npixels = gd.nlon * gd.nlat * gd.ndepth;
    out() << "Grid: " << gd.nlon << " lon x " << gd.nlat << " lat";
    if (gd.ndepth > 1) out() << " x " << gd.ndepth << " depth";
    out() << " x " << gd.ntime << " time = " << npixels << " pixels" << std::endl;

    std::vector<double> event_lon, event_lat, event_depth;
    std::vector<int> pixel_idx;
    std::vector<EventResult> all_events;
    std::vector<int> ds, dp, de;

    // Optional per-day output buffers. The proto set carries the
    // threshold/duration flags; the daily set carries intensity and the daily
    // category. Both share temp/seas/thresh and the event flag/number (not in
    // a slim daily file, which has intensity + category only). temp is written
    // directly from gd.sst.
    //
    // Buffers are allocated uninitialised: detect_events_grid() sets every
    // pixel's slice on its worker thread (std::vector would first zero-fill
    // them on this thread, ~0.4 s for 80M cells).
    const bool want_proto = !proto_file.empty();
    const bool want_daily_series = !daily_file.empty();
    const bool want_daily = want_proto || want_daily_series;
    const bool full_daily = want_daily_series && !p.slim;
    size_t gs = static_cast<size_t>(npixels) * gd.ntime;
    std::unique_ptr<float[]> b_seas, b_thresh, b_int;
    std::unique_ptr<signed char[]> b_tc, b_dc, b_ev, b_cat;
    std::unique_ptr<int[]> b_eno;
    DailyBuffers daily_grid;
    DailyBuffers* daily_grid_ptr = nullptr;
    if (want_daily) {
        if (want_proto || full_daily) {
            b_seas.reset(new float[gs]); b_thresh.reset(new float[gs]);
            b_ev.reset(new signed char[gs]); b_eno.reset(new int[gs]);
            daily_grid.seas = b_seas.get();
            daily_grid.thresh = b_thresh.get();
            daily_grid.event = b_ev.get();
            daily_grid.event_no = b_eno.get();
        }
        if (want_proto) {
            b_tc.reset(new signed char[gs]); b_dc.reset(new signed char[gs]);
            daily_grid.threshCriterion = b_tc.get();
            daily_grid.durationCriterion = b_dc.get();
        }
        if (want_daily_series) {
            b_int.reset(new float[gs]); b_cat.reset(new signed char[gs]);
            daily_grid.intensity = b_int.get();
            daily_grid.category = b_cat.get();
        }
        daily_grid_ptr = &daily_grid;
    }

    out() << "Detecting events with "
          << (p.n_threads > 0 ? p.n_threads : get_threads(npixels)) << " thread(s)..." << std::endl;

    detect_events_grid(
        gd.sst.data(),
        thresh2_gd ? thresh2_gd->sst.data() : nullptr,
        cd.seas.data(), cd.thresh.data(),
        gd.time_days.data(), npixels, gd.ntime,
        p.minDuration, p.minDuration2,
        p.joinAcrossGaps, p.maxGap, p.maxGap2,
        p.coldSpells, p.roundRes, p.n_threads,
        p.category, p.southHemisphere,
        event_lon, event_lat, pixel_idx, all_events,
        ds, dp, de,
        gd.lon, gd.lat, gd.nlon, gd.nlat,
        daily_grid_ptr,
        gd.ndepth, gd.depth, &event_depth
    );

    out() << "Found " << all_events.size() << " events across "
          << npixels << " pixels" << std::endl;

    // Per-day proto-event series (threshCriterion / durationCriterion set).
    // Written even when no events were detected.
    if (want_proto) {
        out() << "Writing per-day proto-event series to " << proto_file << "..." << std::endl;
        write_daily_netcdf(
            proto_file, gd.lon, gd.lat, gd.nlon, gd.nlat, gd.ntime, gd.time_days,
            gd.sst.data(), b_seas.get(), b_thresh.get(),
            b_tc.get(), b_dc.get(), b_ev.get(), b_eno.get(),
            nullptr, nullptr,
            "protoevents", source_label, clim_file, p.minDuration, p.maxGap,
            p.coldSpells, p.southHemisphere, gd.temp_units,
            gd.ndepth, gd.depth, wo);
        out() << "Done." << std::endl;
    }

    // Per-day series with intensity + daily category (and, unless slim, temp,
    // seas, thresh and the event flag/number).
    if (want_daily_series) {
        out() << "Writing per-day series" << (p.slim ? " (slim)" : "") << " to "
              << daily_file << "..." << std::endl;
        write_daily_netcdf(
            daily_file, gd.lon, gd.lat, gd.nlon, gd.nlat, gd.ntime, gd.time_days,
            p.slim ? nullptr : gd.sst.data(),
            p.slim ? nullptr : b_seas.get(), p.slim ? nullptr : b_thresh.get(),
            nullptr, nullptr,
            p.slim ? nullptr : b_ev.get(), p.slim ? nullptr : b_eno.get(),
            b_int.get(), b_cat.get(),
            "daily", source_label, clim_file, p.minDuration, p.maxGap,
            p.coldSpells, p.southHemisphere, gd.temp_units,
            gd.ndepth, gd.depth, wo);
        out() << "Done." << std::endl;
    }

    // Event-metrics table.
    if (!events_file.empty()) {
        if (all_events.empty()) {
            out() << "No events detected. Skipping event output." << std::endl;
        } else {
            int ref_jd = gd.time_days[0];
            std::vector<int> ds_rel(ds.size()), dp_rel(dp.size()), de_rel(de.size());
            for (size_t i = 0; i < ds.size(); ++i) {
                ds_rel[i] = ds[i] - ref_jd;
                dp_rel[i] = dp[i] - ref_jd;
                de_rel[i] = de[i] - ref_jd;
            }
            out() << "Writing events to " << events_file << "..." << std::endl;
            write_event_netcdf(
                events_file, event_lon, event_lat, pixel_idx, all_events,
                ds_rel, dp_rel, de_rel, ref_jd,
                source_label, clim_file,
                p.minDuration, p.maxGap, p.coldSpells,
                p.southHemisphere, gd.temp_units, event_depth, wo
            );
            out() << "Done." << std::endl;
        }
    }
}

void detect_events(const std::string& file_in, const std::string& clim_file,
                   const std::string& events_file, std::string var_name,
                   const DetectParams& p,
                   const std::string& threshClim2_file,
                   const std::string& threshClim2_var_name,
                   const std::string& daily_file, const std::string& proto_file) {
    out() << "Reading climatology from " << clim_file << "..." << std::endl;
    ClimData cd = read_clim_netcdf(clim_file);
    SubsetSpec ss = clim_extent_subset(cd);

    if (var_name.empty()) var_name = detect_sst_variable(file_in);

    out() << "Reading SST data from " << file_in << "..." << std::endl;
    GridData gd = read_sst_netcdf(file_in, var_name, ss);

    std::unique_ptr<GridData> thresh2_gd;
    if (!threshClim2_file.empty()) {
        std::string t2_var = threshClim2_var_name;
        if (t2_var.empty()) t2_var = detect_sst_variable(threshClim2_file);
        out() << "Reading secondary threshold from " << threshClim2_file << "..." << std::endl;
        thresh2_gd.reset(new GridData(read_sst_netcdf(threshClim2_file, t2_var, ss)));
    }

    detect_and_write_events(gd, cd, thresh2_gd.get(), events_file, file_in, clim_file,
                            p, daily_file, proto_file);
}

void detect_events_multi(const std::vector<std::string>& files,
                         const std::string& clim_file,
                         const std::string& events_file, std::string var_name,
                         const DetectParams& p,
                         const std::vector<std::string>& threshClim2_files,
                         const std::string& threshClim2_var_name,
                         bool skip_bad_files,
                         const std::string& daily_file, const std::string& proto_file) {
    out() << "Reading climatology from " << clim_file << "..." << std::endl;
    ClimData cd = read_clim_netcdf(clim_file);
    SubsetSpec ss = clim_extent_subset(cd);

    if (var_name.empty() && !skip_bad_files && !files.empty())
        var_name = detect_sst_variable(files[0]);

    out() << "Reading " << files.size() << " daily files..." << std::endl;
    GridData gd = read_sst_multi_netcdf(files, var_name, ss, skip_bad_files);

    std::unique_ptr<GridData> thresh2_gd;
    if (!threshClim2_files.empty()) {
        std::string t2_var = threshClim2_var_name;
        if (t2_var.empty() && !skip_bad_files)
            t2_var = detect_sst_variable(threshClim2_files[0]);
        out() << "Reading " << threshClim2_files.size()
              << " secondary threshold files..." << std::endl;
        thresh2_gd.reset(new GridData(read_sst_multi_netcdf(threshClim2_files, t2_var, ss,
                                                            skip_bad_files)));
    }

    std::string source_label = std::to_string(files.size()) + " daily files";
    detect_and_write_events(gd, cd, thresh2_gd.get(), events_file, source_label, clim_file,
                            p, daily_file, proto_file);
}

// ---- Per-day detection on arrays ----

void detect_events_array(const double* temp, const double* seas, const double* thresh,
                         int npixels, int ntime, const DetectParams& p,
                         const ArrayDetectOutputs& out) {
    if (p.minDuration < 1) throw std::runtime_error("minDuration must be at least 1.");
    if (p.maxGap < 0) throw std::runtime_error("maxGap must not be negative.");
    DailyBuffers buffers;
    buffers.event = out.event;
    buffers.event_no = out.event_no;
    buffers.intensity = out.intensity;
    buffers.category = out.category;
    detect_events_array_grid(temp, seas, thresh, npixels, ntime,
                             p.minDuration, p.joinAcrossGaps, p.maxGap,
                             p.coldSpells, p.roundRes, p.n_threads, buffers);
}

// ---- Daily categories for a date window ----

// Core of the daily-category computation: given an already-read SST window
// (GridData), climatology, and events, return the full per-pixel per-day grid
// for the window. Event membership (event / event_no) is read from the events
// file -- no detection is re-run, so events that began before the window are
// still recognised. Mirrors the columns of the detect_event3 daily ("also")
// product: temp, seas, thresh, intensity for every day; event/event_no from the
// event ranges; and the Hobday (2018) category on event-member exceedance days
// (NA otherwise). Cold-spells additionally get category 5 (ice) where the
// threshold falls below ice_thresh.
static DailyCatResult daily_cat_core(const GridData& gd, const ClimData& cd,
                                     const EventData& ed,
                                     bool coldSpells, double ice_thresh, int roundRes) {
    if (gd.nlon != cd.nlon || gd.nlat != cd.nlat || gd.ndepth != cd.ndepth) {
        throw std::runtime_error(strf(
            "Grid mismatch: SST is %d x %d x %d depth but climatology is %d x %d x %d depth",
            gd.nlon, gd.nlat, gd.ndepth, cd.nlon, cd.nlat, cd.ndepth));
    }
    bool has_depth = (gd.ndepth > 1) && !gd.depth.empty();
    int npixels = gd.nlon * gd.nlat * gd.ndepth;

    // Per-pixel event ranges in absolute Julian Days.
    struct Ev { int ds; int de; int no; };
    std::vector<std::vector<Ev>> by_pixel(npixels);
    for (int i = 0; i < ed.nevents; ++i) {
        int px = ed.pixel_index[i];
        if (px < 0 || px >= npixels) continue;
        by_pixel[px].push_back({ ed.ref_date_jd + ed.date_start[i],
                                 ed.ref_date_jd + ed.date_end[i],
                                 ed.event_no[i] });
    }

    auto rnd = [&](double v) { return (roundRes > 0) ? r_round(v, roundRes) : v; };

    // Full daily grid: every pixel x every day in the window.
    size_t nrow = static_cast<size_t>(npixels) * gd.ntime;
    DailyCatResult res;
    res.has_depth = has_depth;
    res.jd.resize(nrow); res.event_no.resize(nrow); res.category.resize(nrow);
    res.lon.resize(nrow); res.lat.resize(nrow); res.temp.resize(nrow);
    res.seas.resize(nrow); res.thresh.resize(nrow); res.intensity.resize(nrow);
    res.event.resize(nrow);
    res.depth.resize(has_depth ? nrow : 0);

    size_t r = 0;
    for (int px = 0; px < npixels; ++px) {
        const auto& evs = by_pixel[px];
        int idepth = px % gd.ndepth;
        int lonlat0 = px / gd.ndepth;
        int ilon = lonlat0 / gd.nlat, ilat = lonlat0 % gd.nlat;
        double plon = gd.lon[ilon], plat = gd.lat[ilat];
        double pdepth = has_depth ? gd.depth[idepth] : NA_DOUBLE;
        const double* temp = gd.sst.data() + static_cast<size_t>(px) * gd.ntime;
        const double* seas = cd.seas.data() + static_cast<size_t>(px) * 366;
        const double* thr  = cd.thresh.data() + static_cast<size_t>(px) * 366;
        CategoryLevels levels;
        category_levels(seas, thr, 4, levels);

        for (int t = 0; t < gd.ntime; ++t, ++r) {
            int jd = gd.time_days[t];
            res.jd[r]  = jd;
            res.lon[r] = plon;
            res.lat[r] = plat;
            if (has_depth) res.depth[r] = pdepth;

            double tv = temp[t];
            int doy = jd_to_doy_366(jd) - 1;
            bool clim_ok = (doy >= 0 && doy < 366);
            double s = clim_ok ? seas[doy] : NA_DOUBLE;
            double th = clim_ok ? thr[doy] : NA_DOUBLE;
            bool tv_ok = !is_na(tv);
            bool s_ok  = clim_ok && !is_na(s);
            bool th_ok = clim_ok && !is_na(th);

            res.temp[r]      = tv_ok ? tv : NA_DOUBLE;
            res.seas[r]      = s_ok  ? s  : NA_DOUBLE;
            res.thresh[r]    = th_ok ? th : NA_DOUBLE;
            res.intensity[r] = (tv_ok && s_ok) ? rnd(tv - s) : NA_DOUBLE;

            // Event membership from the full-record events file.
            int eno = 0;
            for (const auto& e : evs) {
                if (jd >= e.ds && jd <= e.de) { eno = e.no; break; }
            }
            res.event[r]    = (eno != 0);
            res.event_no[r] = (eno != 0) ? eno : NA_INT;

            // Category on event days, as heatwaveR::category(climatology =
            // TRUE); for cold-spells, every event day whose threshold is below
            // ice_thresh is ice (5), as heatwaveR's MCSice.
            int cat = NA_INT;
            if (eno != 0) {
                int c = daily_category(tv, doy + 1, levels, coldSpells);
                if (c > 0) cat = c;
                if (coldSpells && th_ok && th < ice_thresh) cat = 5;
            }
            res.category[r] = cat;
        }
    }
    return res;
}

// Daily marine-heatwave / cold-spell categories for a date window, built from a
// single multi-time SST file + a heatwave3 climatology + a heatwave3 events
// file. Only the SST time-window is read (hyperslab); no detection is re-run.
DailyCatResult category_daily(const std::string& sst_file,
                              const std::string& clim_file,
                              const std::string& event_file,
                              const std::vector<int>& t_jd_range,
                              std::string var_name, bool coldSpells,
                              double ice_thresh, int roundRes) {
    ClimData cd = read_clim_netcdf(clim_file);
    SubsetSpec ss = clim_extent_subset(cd, t_jd_range);
    if (var_name.empty()) var_name = detect_sst_variable(sst_file);
    GridData gd = read_sst_netcdf(sst_file, var_name, ss);
    EventData ed = read_event_netcdf(event_file);
    return daily_cat_core(gd, cd, ed, coldSpells, ice_thresh, roundRes);
}

// As above, but reading the SST window from a vector of daily files (one time
// step each), e.g. daily-global OSTIA/GHRSST files for the requested window.
DailyCatResult category_daily_multi(const std::vector<std::string>& files,
                                    const std::string& clim_file,
                                    const std::string& event_file,
                                    const std::vector<int>& t_jd_range,
                                    std::string var_name, bool coldSpells,
                                    double ice_thresh, int roundRes,
                                    bool skip_bad_files) {
    ClimData cd = read_clim_netcdf(clim_file);
    SubsetSpec ss = clim_extent_subset(cd, t_jd_range);
    if (var_name.empty() && !skip_bad_files && !files.empty())
        var_name = detect_sst_variable(files[0]);
    GridData gd = read_sst_multi_netcdf(files, var_name, ss, skip_bad_files);
    EventData ed = read_event_netcdf(event_file);
    return daily_cat_core(gd, cd, ed, coldSpells, ice_thresh, roundRes);
}

// ---- Event-file summaries ----

CategoryResult category(const std::string& event_file,
                        const std::string& clim_file,
                        const std::vector<std::string>& sst_files,
                        std::string var_name, bool skip_bad_files,
                        const std::string& hemisphere, int roundVal) {
    EventData ed = read_event_netcdf(event_file);
    bool southHemisphere = (hemisphere == "south");

    int n = ed.nevents;

    CategoryResult res;
    res.category.assign(n, 0);
    res.p_moderate.assign(n, NA_DOUBLE);
    res.p_strong.assign(n, NA_DOUBLE);
    res.p_severe.assign(n, NA_DOUBLE);
    res.p_extreme.assign(n, NA_DOUBLE);
    res.season.assign(n, 0);
    res.peak_jd.assign(n, 0);

    for (int i = 0; i < n; ++i) res.peak_jd[i] = ed.ref_date_jd + ed.date_peak[i];

    // Check if event file already has category data. The "category" NetCDF
    // variable is always defined by write_event_netcdf() -- filled with 0
    // when detect_event3() was called with category = FALSE -- so emptiness
    // alone can't distinguish "no data" from "real all-zero placeholder".
    // Require at least one real (1-4) category code, mirroring the R-level
    // check in category3().
    bool has_precomputed = false;
    for (int c : ed.category) {
        if (c >= 1 && c <= 4) { has_precomputed = true; break; }
    }

    if (has_precomputed) {
        // Use pre-computed category data from event file
        for (int i = 0; i < n; ++i) {
            int cat = ed.category[i];
            if (cat >= 1 && cat <= 4) res.category[i] = cat;
            res.p_moderate[i] = ed.p_moderate[i];
            res.p_strong[i] = ed.p_strong[i];
            res.p_severe[i] = ed.p_severe[i];
            res.p_extreme[i] = ed.p_extreme[i];
            int sea = ed.season[i];
            if (sea >= 1 && sea <= 4) res.season[i] = sea;
        }
    } else {
        // Compute as heatwaveR::category(), which needs each event day's
        // temperature: read the SST over the climatology grid.
        if (clim_file.empty() || sst_files.empty())
            throw std::runtime_error("Event file has no pre-computed categories; clim_file and "
                                     "sst_file are required to compute them.");
        ClimData cd = read_clim_netcdf(clim_file);
        SubsetSpec ss = clim_extent_subset(cd);
        GridData gd;
        if (sst_files.size() == 1) {
            if (var_name.empty()) var_name = detect_sst_variable(sst_files[0]);
            gd = read_sst_netcdf(sst_files[0], var_name, ss);
        } else {
            if (var_name.empty() && !skip_bad_files) var_name = detect_sst_variable(sst_files[0]);
            gd = read_sst_multi_netcdf(sst_files, var_name, ss, skip_bad_files);
        }
        if (gd.nlon != cd.nlon || gd.nlat != cd.nlat || gd.ndepth != cd.ndepth) {
            throw std::runtime_error(strf(
                "Grid mismatch: SST is %d x %d x %d depth but climatology is %d x %d x %d depth",
                gd.nlon, gd.nlat, gd.ndepth, cd.nlon, cd.nlat, cd.ndepth));
        }
        const int npixels = gd.nlon * gd.nlat * gd.ndepth;

        // Cold-spells, as heatwaveR: any negative maximum intensity
        bool coldSpells = false;
        for (double v : ed.intensity_max) if (v < 0) { coldSpells = true; break; }

        // Events grouped by pixel, so each pixel's levels are computed once
        std::vector<int> order(n);
        for (int i = 0; i < n; ++i) order[i] = i;
        std::stable_sort(order.begin(), order.end(),
                         [&](int x, int y) { return ed.pixel_index[x] < ed.pixel_index[y]; });

        CategoryLevels levels;
        int levels_px = -1;
        std::vector<double> temp;
        std::vector<int> doy;
        for (int i : order) {
            int px = ed.pixel_index[i];
            if (px < 0 || px >= npixels) continue;
            if (px != levels_px) {
                category_levels(cd.seas.data() + static_cast<size_t>(px) * 366,
                                cd.thresh.data() + static_cast<size_t>(px) * 366,
                                roundVal, levels);
                levels_px = px;
            }
            // The event's days on a dense daily axis (missing days are NaN)
            int js = ed.ref_date_jd + ed.date_start[i];
            int je = ed.ref_date_jd + ed.date_end[i];
            temp.assign(je - js + 1, NA_DOUBLE);
            doy.resize(je - js + 1);
            for (int jd = js; jd <= je; ++jd) doy[jd - js] = jd_to_doy_366(jd);
            const double* series = gd.sst.data() + static_cast<size_t>(px) * gd.ntime;
            auto it = std::lower_bound(gd.time_days.begin(), gd.time_days.end(), js);
            for (; it != gd.time_days.end() && *it <= je; ++it)
                temp[*it - js] = series[it - gd.time_days.begin()];

            EventCategory ec = event_category(temp.data(), doy.data(), je - js + 1,
                                              ed.duration[i], levels, coldSpells);
            res.category[i] = ec.category;
            res.p_moderate[i] = ec.p_moderate;
            res.p_strong[i] = ec.p_strong;
            res.p_severe[i] = ec.p_severe;
            res.p_extreme[i] = ec.p_extreme;

            // Season from peak month (codes: 1 Summer, 2 Fall, 3 Winter, 4 Spring)
            int y, m, d;
            jd_to_ymd(res.peak_jd[i], y, m, d);
            if (southHemisphere) {
                if (m == 12 || m <= 2)      res.season[i] = 1;
                else if (m >= 3 && m <= 5)  res.season[i] = 2;
                else if (m >= 6 && m <= 8)  res.season[i] = 3;
                else                        res.season[i] = 4;
            } else {
                if (m == 12 || m <= 2)      res.season[i] = 3;
                else if (m >= 3 && m <= 5)  res.season[i] = 4;
                else if (m >= 6 && m <= 8)  res.season[i] = 1;
                else                        res.season[i] = 2;
            }
        }
    }

    // Round lon/lat for display; i_max as heatwaveR (round(intensity_max, roundVal))
    res.event_no = ed.event_no;
    res.duration = ed.duration;
    res.lon.resize(n); res.lat.resize(n); res.intensity_max.resize(n);
    for (int i = 0; i < n; ++i) {
        res.lon[i] = r_round(ed.lon[i], roundVal);
        res.lat[i] = r_round(ed.lat[i], roundVal);
        res.intensity_max[i] = r_round(ed.intensity_max[i], roundVal);
    }
    return res;
}

BlockAverageResult block_average(const std::string& event_file) {
    EventData ed = read_event_netcdf(event_file);
    BlockAverageResult res;

    int n = ed.nevents;
    if (n == 0) return res;

    // Compute year from date_start for each event
    std::vector<int> year(n);
    for (int i = 0; i < n; ++i) {
        int jd = ed.ref_date_jd + ed.date_start[i];
        int y, m, d;
        jd_to_ymd(jd, y, m, d);
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

    for (auto& [key, g] : groups) {
        res.lon.push_back(g.lon);
        res.lat.push_back(g.lat);
        res.year.push_back(g.year);
        res.count.push_back(g.count);
        res.duration_mean.push_back(static_cast<double>(g.duration_sum) / g.count);
        res.duration_max.push_back(g.duration_max);
        res.intensity_mean.push_back(g.imean_sum / g.count);
        res.intensity_max_mean.push_back(g.imax_sum / g.count);
        res.intensity_max_max.push_back(g.imax_max);
        res.intensity_cumulative_mean.push_back(g.icum_sum / g.count);
        res.total_days.push_back(g.duration_sum);
        res.total_icum.push_back(g.icum_sum);
    }
    return res;
}

MetricSummaryResult metric_summary(const std::string& event_file,
                                   const std::string& metric,
                                   const std::string& summary) {
    EventData ed = read_event_netcdf(event_file);
    MetricSummaryResult res;
    int n = ed.nevents;
    if (n == 0) return res;

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
    else throw std::runtime_error("Unknown metric: " + metric);

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

    for (auto& [key, g] : groups) {
        res.lon.push_back(g.lon);
        res.lat.push_back(g.lat);
        if (summary == "mean") res.value.push_back(g.sum / g.count);
        else if (summary == "max") res.value.push_back(g.mx);
        else if (summary == "min") res.value.push_back(g.mn);
        else if (summary == "sum") res.value.push_back(g.sum);
        else if (summary == "count") res.value.push_back(static_cast<double>(g.count));
        else throw std::runtime_error("Unknown summary: " + summary);
    }
    return res;
}

} // namespace hw3
