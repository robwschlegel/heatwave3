#ifndef HEATWAVE3_NETCDF_IO_H
#define HEATWAVE3_NETCDF_IO_H

#include <string>
#include <vector>
#include <Rcpp.h>
#include "heatwave3_types.h"

namespace hw3 {

struct DimInfo {
    std::string name;
    size_t len;
    std::vector<double> values;
};

struct GridData {
    std::vector<double> lon;
    std::vector<double> lat;
    std::vector<double> time_raw;
    std::vector<int> time_days;
    std::string time_units;
    std::string time_calendar;
    std::string temp_units;
    int nlon;
    int nlat;
    int ntime;
    std::vector<double> sst;  // [nlon * nlat * ntime], pixel-major order
};

struct ClimData {
    std::vector<double> lon;
    std::vector<double> lat;
    int nlon;
    int nlat;
    int ndoy;
    std::vector<double> seas;   // [nlon * nlat * 366]
    std::vector<double> thresh; // [nlon * nlat * 366]
    std::vector<double> var;    // [nlon * nlat * 366], may be empty
};

// Per-day, per-pixel output (protoEvent / daily). Variable vectors are filled
// only when present in the file; absent ones are left empty. All data vectors
// are [nlon * nlat * ntime] in pixel-major ([lon][lat][time]) order.
struct DailyData {
    std::vector<double> lon;
    std::vector<double> lat;
    std::vector<int> time_jd;          // Julian Day per time step
    int nlon = 0;
    int nlat = 0;
    int ntime = 0;
    std::vector<double> temp;
    std::vector<double> seas;
    std::vector<double> thresh;
    std::vector<double> intensity;
    std::vector<signed char> threshCriterion;
    std::vector<signed char> durationCriterion;
    std::vector<signed char> event;
    std::vector<signed char> category;
    std::vector<int> event_no;
};

// Lightweight metadata for any heatwave3 product file (reads attributes and
// dimension lengths only, never the data arrays).
struct FileMeta {
    std::string product;   // "climatology" | "events" | "daily" | "protoevents"
    int nlon = 0;
    int nlat = 0;
    int n3 = 0;            // ndoy, ntime, or nevents depending on product
    long long nrows = 0;   // long-table row count
};

void nc_check(int status, const std::string& context);

FileMeta read_file_meta(const std::string& file);

GridData read_sst_netcdf(const std::string& file_in,
                         const std::string& var_name,
                         const SubsetSpec& subset);

std::string detect_sst_variable(const std::string& file_in);

void parse_cf_time(const std::string& units,
                   const std::string& calendar,
                   const std::vector<double>& time_raw,
                   std::vector<int>& julian_days);

void write_clim_netcdf(const std::string& file_out,
                       const std::vector<double>& lon,
                       const std::vector<double>& lat,
                       int nlon, int nlat,
                       const std::vector<double>& seas,
                       const std::vector<double>& thresh,
                       const std::vector<double>& var,
                       bool has_var,
                       const std::string& source_file,
                       const std::string& clim_period_start,
                       const std::string& clim_period_end,
                       double pctile,
                       int windowHalfWidth,
                       int smoothPercentileWidth,
                       const std::string& temp_units);

void write_event_netcdf(const std::string& file_out,
                        const std::vector<double>& event_lon,
                        const std::vector<double>& event_lat,
                        const std::vector<int>& pixel_index,
                        const std::vector<EventResult>& events,
                        const std::vector<int>& date_start,
                        const std::vector<int>& date_peak,
                        const std::vector<int>& date_end,
                        int ref_date_jd,
                        const std::string& source_file,
                        const std::string& clim_file,
                        int minDuration,
                        int maxGap,
                        bool coldSpells,
                        bool southHemisphere,
                        const std::string& temp_units);

ClimData read_clim_netcdf(const std::string& clim_file);

EventData read_event_netcdf(const std::string& event_file);

// Write a per-day, per-pixel NetCDF ([lon, lat, time]). temp and the daily
// buffers are [nlon*nlat*ntime] pixel-major; any null buffer pointer omits
// that variable (this selects the protoEvent vs per-day variable sets).
void write_daily_netcdf(const std::string& file_out,
                        const std::vector<double>& lon,
                        const std::vector<double>& lat,
                        int nlon, int nlat, int ntime,
                        const std::vector<int>& time_days,
                        const double* temp,
                        const double* seas,
                        const double* thresh,
                        const signed char* threshCriterion,
                        const signed char* durationCriterion,
                        const signed char* event,
                        const int* event_no,
                        const double* intensity,
                        const signed char* category,
                        const std::string& product,
                        const std::string& source_file,
                        const std::string& clim_file,
                        int minDuration,
                        int maxGap,
                        bool coldSpells,
                        bool southHemisphere,
                        const std::string& temp_units);

DailyData read_daily_netcdf(const std::string& daily_file);

// Streaming hyperslab read of a gridded product (climatology / daily /
// protoevents). Reads only the requested lon/lat/time index window and only
// the requested data variables. Empty range vectors mean "full extent"; an
// empty `vars` means "all data variables present". `max_rows` (>= 0) caps the
// number of long-table rows by limiting how many longitude columns are read
// (<0 = no cap). Returns an Rcpp::List consumed by the R builder. Throws if the
// file is an events product (those are filtered in R).
Rcpp::List read_subset_netcdf(const std::string& file,
                              const std::vector<double>& lon_range,
                              const std::vector<double>& lat_range,
                              const std::vector<int>& t_jd_range,
                              const std::vector<std::string>& vars,
                              int max_rows);

// Read and merge multiple daily NetCDF files into a single GridData
GridData read_sst_multi_netcdf(const std::vector<std::string>& files,
                               const std::string& var_name,
                               const SubsetSpec& subset,
                               bool skip_bad_files = false);

} // namespace hw3

#endif
