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

void nc_check(int status, const std::string& context);

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
                       int smoothPercentileWidth);

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
                        bool coldSpells);

ClimData read_clim_netcdf(const std::string& clim_file);

// Read and merge multiple daily NetCDF files into a single GridData
GridData read_sst_multi_netcdf(const std::vector<std::string>& files,
                               const std::string& var_name,
                               const SubsetSpec& subset);

} // namespace hw3

#endif
