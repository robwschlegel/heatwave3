#include "netcdf_io.h"
#include <netcdf.h>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <regex>

namespace hw3 {

void nc_check(int status, const std::string& context) {
    if (status != NC_NOERR) {
        throw std::runtime_error("NetCDF error in " + context + ": " +
                                 nc_strerror(status));
    }
}

// Read a text attribute that may be NC_CHAR or NC_STRING
static bool nc_get_att_str(int ncid, int varid, const char* name, std::string& out) {
    nc_type atype;
    size_t alen;
    if (nc_inq_att(ncid, varid, name, &atype, &alen) != NC_NOERR) return false;

    if (atype == NC_STRING) {
        char* sp = nullptr;
        if (nc_get_att_string(ncid, varid, name, &sp) != NC_NOERR) return false;
        out = sp;
        nc_free_string(1, &sp);
    } else {
        std::vector<char> buf(alen + 1, 0);
        if (nc_get_att_text(ncid, varid, name, buf.data()) != NC_NOERR) return false;
        out = std::string(buf.data(), alen);
    }
    // Trim trailing nulls
    while (!out.empty() && out.back() == '\0') out.pop_back();
    return true;
}

// ---- CF time parsing ----

static int date_to_jd(int y, int m, int d) {
    // Julian Day Number from Gregorian calendar date
    int a = (14 - m) / 12;
    int yy = y + 4800 - a;
    int mm = m + 12 * a - 3;
    return d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045;
}

static void jd_to_date(int jd, int& y, int& m, int& d) {
    int a = jd + 32044;
    int b = (4 * a + 3) / 146097;
    int c = a - (146097 * b) / 4;
    int dd = (4 * c + 3) / 1461;
    int e = c - (1461 * dd) / 4;
    int mm = (5 * e + 2) / 153;
    d = e - (153 * mm + 2) / 5 + 1;
    m = mm + 3 - 12 * (mm / 10);
    y = 100 * b + dd - 4800 + mm / 10;
}

void parse_cf_time(const std::string& units,
                   const std::string& calendar,
                   const std::vector<double>& time_raw,
                   std::vector<int>& julian_days) {
    // Parse "seconds since YYYY-MM-DD" or "days since YYYY-MM-DD" etc.
    std::regex re(R"((seconds|minutes|hours|days)\s+since\s+(\d{4})-(\d{1,2})-(\d{1,2}))");
    std::smatch match;
    if (!std::regex_search(units, match, re)) {
        throw std::runtime_error("Cannot parse CF time units: " + units);
    }

    std::string unit_type = match[1].str();
    int ref_y = std::stoi(match[2].str());
    int ref_m = std::stoi(match[3].str());
    int ref_d = std::stoi(match[4].str());
    int ref_jd = date_to_jd(ref_y, ref_m, ref_d);

    double scale = 1.0;
    if (unit_type == "seconds") scale = 1.0 / 86400.0;
    else if (unit_type == "minutes") scale = 1.0 / 1440.0;
    else if (unit_type == "hours") scale = 1.0 / 24.0;
    else if (unit_type == "days") scale = 1.0;

    julian_days.resize(time_raw.size());
    for (size_t i = 0; i < time_raw.size(); ++i) {
        julian_days[i] = ref_jd + static_cast<int>(std::round(time_raw[i] * scale));
    }
}

// ---- Variable detection ----

std::string detect_sst_variable(const std::string& file_in) {
    int ncid;
    nc_check(nc_open(file_in.c_str(), NC_NOWRITE, &ncid), "open " + file_in);

    int nvars;
    nc_check(nc_inq_nvars(ncid, &nvars), "inq nvars");

    // Priority list of standard names / variable names
    const std::vector<std::string> sst_names = {
        "analysed_sst", "sst", "SST", "sea_surface_temperature",
        "temperature", "temp", "thetao", "tos", "t_an", "sst_mean"
    };
    const std::vector<std::string> sst_standard_names = {
        "sea_surface_temperature",
        "sea_surface_foundation_temperature",
        "sea_surface_skin_temperature",
        "sea_water_temperature"
    };

    std::string found_name;

    for (int v = 0; v < nvars; ++v) {
        char vname[NC_MAX_NAME + 1];
        nc_check(nc_inq_varname(ncid, v, vname), "inq varname");

        // Check standard_name attribute
        std::string sn;
        if (nc_get_att_str(ncid, v, "standard_name", sn)) {
            for (auto& s : sst_standard_names) {
                if (sn.find(s) != std::string::npos) {
                    found_name = vname;
                    break;
                }
            }
        }
        if (!found_name.empty()) break;

        // Check variable name
        std::string vn(vname);
        for (auto& s : sst_names) {
            if (vn == s) {
                found_name = vn;
                break;
            }
        }
        if (!found_name.empty()) break;
    }

    nc_close(ncid);

    if (found_name.empty()) {
        throw std::runtime_error("Could not auto-detect SST variable in " + file_in);
    }
    return found_name;
}

// ---- Subset index helpers ----

static void find_subset_indices(const std::vector<double>& coord,
                                double min_val, double max_val,
                                size_t& start, size_t& count) {
    if (std::isinf(min_val) && std::isinf(max_val)) {
        start = 0;
        count = coord.size();
        return;
    }
    // Use tolerance for float-precision coordinates
    double eps = 1e-4;
    size_t first = coord.size(), last = 0;

    for (size_t i = 0; i < coord.size(); ++i) {
        if (coord[i] >= (min_val - eps) && coord[i] <= (max_val + eps)) {
            if (i < first) first = i;
            if (i > last) last = i;
        }
    }

    if (first >= coord.size()) {
        throw std::runtime_error("No coordinates found within specified range");
    }
    start = first;
    count = last - first + 1;
}

// ---- Main read function ----

GridData read_sst_netcdf(const std::string& file_in,
                         const std::string& var_name,
                         const SubsetSpec& subset) {
    GridData gd;
    int ncid;
    nc_check(nc_open(file_in.c_str(), NC_NOWRITE, &ncid), "open " + file_in);

    // Get variable id
    int varid;
    nc_check(nc_inq_varid(ncid, var_name.c_str(), &varid), "inq varid " + var_name);

    // Get variable dimensions
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    nc_check(nc_inq_var(ncid, varid, nullptr, nullptr, &ndims, dimids, nullptr),
             "inq var dims");

    // ---- Identify lon, lat, time dimensions ----
    // Two-pass strategy:
    //   Pass 1: match dimension/coordinate variable names against known patterns
    //   Pass 2: fallback via CF axis, standard_name, and units attributes
    int lon_dimid = -1, lat_dimid = -1, time_dimid = -1, depth_dimid = -1;
    int lon_varid = -1, lat_varid = -1, time_varid = -1;

    // Helper: classify a role from name, axis, standard_name, units
    auto classify_coord = [&](int did, int cvid, const std::string& dn) {
        std::string ax, sn, un;
        nc_get_att_str(ncid, cvid, "axis", ax);
        nc_get_att_str(ncid, cvid, "standard_name", sn);
        nc_get_att_str(ncid, cvid, "units", un);

        // Longitude
        if (ax == "X" || sn == "longitude" ||
            un == "degrees_east" || un == "degree_east" ||
            dn == "lon" || dn == "longitude" || dn == "x" ||
            dn == "xt_ocean" || dn == "nav_lon" || dn == "geolon" ||
            dn == "ni" || dn == "xh") {
            if (lon_dimid < 0) { lon_dimid = did; lon_varid = cvid; }
        }
        // Latitude
        else if (ax == "Y" || sn == "latitude" ||
                 un == "degrees_north" || un == "degree_north" ||
                 dn == "lat" || dn == "latitude" || dn == "y" ||
                 dn == "yt_ocean" || dn == "nav_lat" || dn == "geolat" ||
                 dn == "nj" || dn == "yh") {
            if (lat_dimid < 0) { lat_dimid = did; lat_varid = cvid; }
        }
        // Time
        else if (ax == "T" || sn == "time" ||
                 un.find("since") != std::string::npos ||
                 dn == "time" || dn == "t" || dn == "time_counter") {
            if (time_dimid < 0) { time_dimid = did; time_varid = cvid; }
        }
        // Depth
        else if (ax == "Z" || sn == "depth" || sn == "altitude" ||
                 dn == "depth" || dn == "level" || dn == "lev" ||
                 dn == "z" || dn == "zlev") {
            if (depth_dimid < 0) depth_dimid = did;
        }
    };

    // Pass 1: look for coordinate variables whose name matches the dimension name
    for (int d = 0; d < ndims; ++d) {
        char dname[NC_MAX_NAME + 1];
        size_t dlen;
        nc_check(nc_inq_dim(ncid, dimids[d], dname, &dlen), "inq dim");
        std::string dn(dname);

        int coord_varid;
        if (nc_inq_varid(ncid, dname, &coord_varid) == NC_NOERR) {
            classify_coord(dimids[d], coord_varid, dn);
        }
    }

    // Pass 2: for any unidentified dimensions, scan ALL 1-D variables
    // that use that dimension for CF attributes
    if (lon_varid < 0 || lat_varid < 0 || time_varid < 0) {
        int total_vars;
        nc_inq_nvars(ncid, &total_vars);

        for (int d = 0; d < ndims; ++d) {
            int did = dimids[d];
            // Skip already-identified dims
            if (did == lon_dimid || did == lat_dimid ||
                did == time_dimid || did == depth_dimid) continue;

            char dname[NC_MAX_NAME + 1];
            size_t dlen;
            nc_inq_dim(ncid, did, dname, &dlen);
            std::string dn(dname);

            for (int v = 0; v < total_vars; ++v) {
                int vndims;
                int vdimids[NC_MAX_VAR_DIMS];
                nc_inq_var(ncid, v, nullptr, nullptr, &vndims, vdimids, nullptr);
                if (vndims == 1 && vdimids[0] == did) {
                    classify_coord(did, v, dn);
                    if (did == lon_dimid || did == lat_dimid ||
                        did == time_dimid || did == depth_dimid) break;
                }
            }
        }
    }

    if (lon_varid < 0 || lat_varid < 0 || time_varid < 0) {
        nc_close(ncid);
        throw std::runtime_error("Cannot identify lon/lat/time dimensions in " + file_in +
                                 ". Ensure the file has coordinate variables with CF-compliant "
                                 "axis, standard_name, or units attributes.");
    }

    // Read full coordinate arrays
    size_t lon_len, lat_len, time_len;
    nc_check(nc_inq_dimlen(ncid, lon_dimid, &lon_len), "inq lon len");
    nc_check(nc_inq_dimlen(ncid, lat_dimid, &lat_len), "inq lat len");
    nc_check(nc_inq_dimlen(ncid, time_dimid, &time_len), "inq time len");

    std::vector<double> full_lon(lon_len), full_lat(lat_len), full_time(time_len);
    nc_check(nc_get_var_double(ncid, lon_varid, full_lon.data()), "get lon");
    nc_check(nc_get_var_double(ncid, lat_varid, full_lat.data()), "get lat");
    nc_check(nc_get_var_double(ncid, time_varid, full_time.data()), "get time");

    // Get time units
    if (!nc_get_att_str(ncid, time_varid, "units", gd.time_units)) {
        nc_close(ncid);
        throw std::runtime_error("Cannot read time units attribute from " + file_in);
    }

    if (!nc_get_att_str(ncid, time_varid, "calendar", gd.time_calendar)) {
        gd.time_calendar = "proleptic_gregorian";
    }

    // Apply spatial subsetting
    size_t lon_start, lon_count, lat_start, lat_count;
    find_subset_indices(full_lon, subset.lon_min, subset.lon_max, lon_start, lon_count);
    find_subset_indices(full_lat, subset.lat_min, subset.lat_max, lat_start, lat_count);

    // Apply time subsetting via Julian days
    std::vector<int> full_jd;
    parse_cf_time(gd.time_units, gd.time_calendar, full_time, full_jd);

    size_t time_start = 0, time_count = time_len;
    if (!std::isinf(subset.time_min) || !std::isinf(subset.time_max)) {
        int tmin_jd = std::isinf(subset.time_min) ? full_jd.front() :
                      static_cast<int>(subset.time_min);
        int tmax_jd = std::isinf(subset.time_max) ? full_jd.back() :
                      static_cast<int>(subset.time_max);
        size_t first = time_len, last = 0;
        for (size_t i = 0; i < time_len; ++i) {
            if (full_jd[i] >= tmin_jd && full_jd[i] <= tmax_jd) {
                if (i < first) first = i;
                if (i > last) last = i;
            }
        }
        if (first >= time_len) {
            nc_close(ncid);
            throw std::runtime_error("No time steps found within specified range");
        }
        time_start = first;
        time_count = last - first + 1;
    }

    gd.nlon = static_cast<int>(lon_count);
    gd.nlat = static_cast<int>(lat_count);
    gd.ntime = static_cast<int>(time_count);

    gd.lon.assign(full_lon.begin() + lon_start,
                  full_lon.begin() + lon_start + lon_count);
    gd.lat.assign(full_lat.begin() + lat_start,
                  full_lat.begin() + lat_start + lat_count);
    gd.time_raw.assign(full_time.begin() + time_start,
                       full_time.begin() + time_start + time_count);
    gd.time_days.assign(full_jd.begin() + time_start,
                        full_jd.begin() + time_start + time_count);

    // Get scale_factor, add_offset, _FillValue
    double scale_factor = 1.0, add_offset = 0.0;
    nc_get_att_double(ncid, varid, "scale_factor", &scale_factor);
    nc_get_att_double(ncid, varid, "add_offset", &add_offset);

    nc_type var_type;
    nc_check(nc_inq_vartype(ncid, varid, &var_type), "inq vartype");

    short fill_short = NC_FILL_SHORT;
    float fill_float = NC_FILL_FLOAT;
    double fill_double = NC_FILL_DOUBLE;
    nc_get_att_short(ncid, varid, "_FillValue", &fill_short);
    nc_get_att_float(ncid, varid, "_FillValue", &fill_float);
    nc_get_att_double(ncid, varid, "_FillValue", &fill_double);

    // Build hyperslab start/count arrays
    // Dimension ordering from the variable's dimids
    std::vector<size_t> start(ndims, 0), count(ndims, 1);
    for (int d = 0; d < ndims; ++d) {
        if (dimids[d] == lon_dimid)  { start[d] = lon_start;  count[d] = lon_count; }
        if (dimids[d] == lat_dimid)  { start[d] = lat_start;  count[d] = lat_count; }
        if (dimids[d] == time_dimid) { start[d] = time_start; count[d] = time_count; }
        if (dimids[d] == depth_dimid) {
            start[d] = (subset.depth_index >= 0) ?
                       static_cast<size_t>(subset.depth_index) : 0;
            count[d] = 1;
        }
    }

    size_t total = lon_count * lat_count * time_count;
    gd.sst.resize(total);

    if (var_type == NC_SHORT || var_type == NC_BYTE) {
        std::vector<short> raw(total);
        nc_check(nc_get_vara_short(ncid, varid, start.data(), count.data(), raw.data()),
                 "get sst data");
        for (size_t i = 0; i < total; ++i) {
            if (raw[i] == fill_short) {
                gd.sst[i] = NA_DOUBLE;
            } else {
                gd.sst[i] = static_cast<double>(raw[i]) * scale_factor + add_offset;
            }
        }
    } else if (var_type == NC_FLOAT) {
        std::vector<float> raw(total);
        nc_check(nc_get_vara_float(ncid, varid, start.data(), count.data(), raw.data()),
                 "get sst data");
        for (size_t i = 0; i < total; ++i) {
            if (raw[i] == fill_float || std::isnan(raw[i])) {
                gd.sst[i] = NA_DOUBLE;
            } else {
                gd.sst[i] = static_cast<double>(raw[i]) * scale_factor + add_offset;
            }
        }
    } else {
        nc_check(nc_get_vara_double(ncid, varid, start.data(), count.data(), gd.sst.data()),
                 "get sst data");
        for (size_t i = 0; i < total; ++i) {
            if (gd.sst[i] == fill_double || std::isnan(gd.sst[i])) {
                gd.sst[i] = NA_DOUBLE;
            }
        }
    }

    // Reorder to pixel-major: sst[pixel][time] where pixel = lon_idx * nlat + lat_idx
    // Input is in file dimension order (typically [lon, lat, time])
    // We want contiguous time series per pixel for efficient per-pixel processing
    std::vector<double> reordered(total);
    // Determine dimension order
    (void)depth_dimid; // may be unused if no depth dim
    // If depth dim was squeezed, adjust
    int eff_idx = 0;
    size_t eff_dims[3] = {1, 1, 1};
    int dim_map[3] = {-1, -1, -1}; // map eff position to role: 0=lon,1=lat,2=time
    for (int d = 0; d < ndims; ++d) {
        if (dimids[d] == depth_dimid && count[d] == 1) continue;
        eff_dims[eff_idx] = count[d];
        if (dimids[d] == lon_dimid)  dim_map[eff_idx] = 0;
        if (dimids[d] == lat_dimid)  dim_map[eff_idx] = 1;
        if (dimids[d] == time_dimid) dim_map[eff_idx] = 2;
        eff_idx++;
    }

    // Reorder: target layout is [pixel][time] where pixel = ilon * nlat + ilat
    size_t stride2 = 1;
    size_t stride1 = eff_dims[2];
    size_t stride0 = eff_dims[1] * eff_dims[2];

    for (size_t i0 = 0; i0 < eff_dims[0]; ++i0) {
        for (size_t i1 = 0; i1 < eff_dims[1]; ++i1) {
            for (size_t i2 = 0; i2 < eff_dims[2]; ++i2) {
                size_t src = i0 * stride0 + i1 * stride1 + i2 * stride2;
                // Map to lon/lat/time indices
                size_t idx[3] = {i0, i1, i2};
                size_t ilon = 0, ilat = 0, itime = 0;
                for (int e = 0; e < eff_idx; ++e) {
                    if (dim_map[e] == 0) ilon  = idx[e];
                    if (dim_map[e] == 1) ilat  = idx[e];
                    if (dim_map[e] == 2) itime = idx[e];
                }
                size_t pixel = ilon * static_cast<size_t>(gd.nlat) + ilat;
                size_t dst = pixel * static_cast<size_t>(gd.ntime) + itime;
                reordered[dst] = gd.sst[src];
            }
        }
    }
    gd.sst = std::move(reordered);

    nc_close(ncid);
    return gd;
}

// ---- Write climatology NetCDF ----

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
                       int smoothPercentileWidth) {
    int ncid;
    nc_check(nc_create(file_out.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid),
             "create " + file_out);

    // Define dimensions
    int lon_dimid, lat_dimid, doy_dimid;
    nc_check(nc_def_dim(ncid, "lon", nlon, &lon_dimid), "def dim lon");
    nc_check(nc_def_dim(ncid, "lat", nlat, &lat_dimid), "def dim lat");
    nc_check(nc_def_dim(ncid, "doy", 366, &doy_dimid), "def dim doy");

    // Define coordinate variables
    int lon_varid, lat_varid, doy_varid;
    nc_check(nc_def_var(ncid, "lon", NC_DOUBLE, 1, &lon_dimid, &lon_varid), "def var lon");
    nc_check(nc_def_var(ncid, "lat", NC_DOUBLE, 1, &lat_dimid, &lat_varid), "def var lat");
    nc_check(nc_def_var(ncid, "doy", NC_INT, 1, &doy_dimid, &doy_varid), "def var doy");

    nc_put_att_text(ncid, lon_varid, "units", 12, "degrees_east");
    nc_put_att_text(ncid, lon_varid, "axis", 1, "X");
    nc_put_att_text(ncid, lat_varid, "units", 13, "degrees_north");
    nc_put_att_text(ncid, lat_varid, "axis", 1, "Y");
    nc_put_att_text(ncid, doy_varid, "long_name", 11, "day_of_year");

    // Define data variables: [lon, lat, doy]
    int dims3[3] = {lon_dimid, lat_dimid, doy_dimid};
    int seas_varid, thresh_varid, var_varid;
    nc_check(nc_def_var(ncid, "seas", NC_DOUBLE, 3, dims3, &seas_varid), "def var seas");
    nc_check(nc_def_var(ncid, "thresh", NC_DOUBLE, 3, dims3, &thresh_varid), "def var thresh");

    nc_put_att_text(ncid, seas_varid, "long_name", 21, "seasonal climatology");
    nc_put_att_text(ncid, seas_varid, "units", 6, "degC");
    nc_put_att_text(ncid, thresh_varid, "long_name", 30, "threshold climatology");
    nc_put_att_text(ncid, thresh_varid, "units", 6, "degC");

    // Enable compression
    nc_def_var_deflate(ncid, seas_varid, 1, 1, 4);
    nc_def_var_deflate(ncid, thresh_varid, 1, 1, 4);

    if (has_var) {
        nc_check(nc_def_var(ncid, "var", NC_DOUBLE, 3, dims3, &var_varid), "def var var");
        nc_put_att_text(ncid, var_varid, "long_name", 21, "variance climatology");
        nc_def_var_deflate(ncid, var_varid, 1, 1, 4);
    }

    // Global attributes
    nc_put_att_text(ncid, NC_GLOBAL, "source_file",
                    source_file.size(), source_file.c_str());
    nc_put_att_text(ncid, NC_GLOBAL, "climatologyPeriod_start",
                    clim_period_start.size(), clim_period_start.c_str());
    nc_put_att_text(ncid, NC_GLOBAL, "climatologyPeriod_end",
                    clim_period_end.size(), clim_period_end.c_str());
    nc_put_att_double(ncid, NC_GLOBAL, "pctile", NC_DOUBLE, 1, &pctile);
    nc_put_att_int(ncid, NC_GLOBAL, "windowHalfWidth", NC_INT, 1, &windowHalfWidth);
    nc_put_att_int(ncid, NC_GLOBAL, "smoothPercentileWidth", NC_INT, 1, &smoothPercentileWidth);
    const char* conv = "CF-1.8";
    nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, conv);
    const char* created_by = "heatwave3";
    nc_put_att_text(ncid, NC_GLOBAL, "created_by", 9, created_by);

    nc_check(nc_enddef(ncid), "enddef");

    // Write coordinates
    nc_check(nc_put_var_double(ncid, lon_varid, lon.data()), "put lon");
    nc_check(nc_put_var_double(ncid, lat_varid, lat.data()), "put lat");
    std::vector<int> doy_vals(366);
    for (int i = 0; i < 366; ++i) doy_vals[i] = i + 1;
    nc_check(nc_put_var_int(ncid, doy_varid, doy_vals.data()), "put doy");

    // Write data — input layout is [pixel][doy] where pixel = ilon*nlat + ilat
    // NetCDF layout is [lon][lat][doy]
    // These are the same order, so write directly
    nc_check(nc_put_var_double(ncid, seas_varid, seas.data()), "put seas");
    nc_check(nc_put_var_double(ncid, thresh_varid, thresh.data()), "put thresh");
    if (has_var) {
        nc_check(nc_put_var_double(ncid, var_varid, var.data()), "put var");
    }

    nc_check(nc_close(ncid), "close " + file_out);
}

// ---- Read climatology NetCDF ----

ClimData read_clim_netcdf(const std::string& clim_file) {
    ClimData cd;
    int ncid;
    nc_check(nc_open(clim_file.c_str(), NC_NOWRITE, &ncid), "open " + clim_file);

    int lon_varid, lat_varid, seas_varid, thresh_varid;
    nc_check(nc_inq_varid(ncid, "lon", &lon_varid), "inq lon");
    nc_check(nc_inq_varid(ncid, "lat", &lat_varid), "inq lat");
    nc_check(nc_inq_varid(ncid, "seas", &seas_varid), "inq seas");
    nc_check(nc_inq_varid(ncid, "thresh", &thresh_varid), "inq thresh");

    int lon_dimid, lat_dimid, doy_dimid;
    nc_check(nc_inq_dimid(ncid, "lon", &lon_dimid), "inq dim lon");
    nc_check(nc_inq_dimid(ncid, "lat", &lat_dimid), "inq dim lat");
    nc_check(nc_inq_dimid(ncid, "doy", &doy_dimid), "inq dim doy");

    size_t nlon, nlat, ndoy;
    nc_check(nc_inq_dimlen(ncid, lon_dimid, &nlon), "inq nlon");
    nc_check(nc_inq_dimlen(ncid, lat_dimid, &nlat), "inq nlat");
    nc_check(nc_inq_dimlen(ncid, doy_dimid, &ndoy), "inq ndoy");

    cd.nlon = static_cast<int>(nlon);
    cd.nlat = static_cast<int>(nlat);
    cd.ndoy = static_cast<int>(ndoy);

    cd.lon.resize(nlon);
    cd.lat.resize(nlat);
    nc_check(nc_get_var_double(ncid, lon_varid, cd.lon.data()), "get lon");
    nc_check(nc_get_var_double(ncid, lat_varid, cd.lat.data()), "get lat");

    size_t total = nlon * nlat * ndoy;
    cd.seas.resize(total);
    cd.thresh.resize(total);
    nc_check(nc_get_var_double(ncid, seas_varid, cd.seas.data()), "get seas");
    nc_check(nc_get_var_double(ncid, thresh_varid, cd.thresh.data()), "get thresh");

    // Check for var
    int var_varid;
    if (nc_inq_varid(ncid, "var", &var_varid) == NC_NOERR) {
        cd.var.resize(total);
        nc_check(nc_get_var_double(ncid, var_varid, cd.var.data()), "get var");
    }

    nc_close(ncid);
    return cd;
}

// ---- Write event NetCDF (ragged array) ----

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
                        bool coldSpells) {
    int ncid;
    nc_check(nc_create(file_out.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid),
             "create " + file_out);

    size_t nevents = events.size();

    // Define event dimension
    int event_dimid;
    nc_check(nc_def_dim(ncid, "event", nevents, &event_dimid), "def dim event");

    // Helper to define and register a double variable on event dim
    auto def_dvar = [&](const char* name, const char* long_name, const char* units) -> int {
        int vid;
        nc_check(nc_def_var(ncid, name, NC_DOUBLE, 1, &event_dimid, &vid), name);
        nc_put_att_text(ncid, vid, "long_name", strlen(long_name), long_name);
        if (units) nc_put_att_text(ncid, vid, "units", strlen(units), units);
        nc_def_var_deflate(ncid, vid, 1, 1, 4);
        return vid;
    };
    auto def_ivar = [&](const char* name, const char* long_name, const char* units) -> int {
        int vid;
        nc_check(nc_def_var(ncid, name, NC_INT, 1, &event_dimid, &vid), name);
        nc_put_att_text(ncid, vid, "long_name", strlen(long_name), long_name);
        if (units) nc_put_att_text(ncid, vid, "units", strlen(units), units);
        nc_def_var_deflate(ncid, vid, 1, 1, 4);
        return vid;
    };

    int v_lon = def_dvar("lon", "longitude", "degrees_east");
    int v_lat = def_dvar("lat", "latitude", "degrees_north");
    int v_pixel = def_ivar("pixel_index", "pixel index (ilon * nlat + ilat)", nullptr);
    int v_eno = def_ivar("event_no", "event number within pixel", nullptr);
    int v_ds = def_ivar("date_start", "start date", "days since reference");
    int v_dp = def_ivar("date_peak", "peak date", "days since reference");
    int v_de = def_ivar("date_end", "end date", "days since reference");
    int v_dur = def_ivar("duration", "event duration", "days");
    int v_im = def_dvar("intensity_mean", "mean intensity", "degC");
    int v_ix = def_dvar("intensity_max", "maximum intensity", "degC");
    int v_iv = def_dvar("intensity_var", "intensity variability", "degC");
    int v_ic = def_dvar("intensity_cumulative", "cumulative intensity", "degC days");
    int v_imrt = def_dvar("intensity_mean_relThresh", "mean intensity rel. threshold", "degC");
    int v_ixrt = def_dvar("intensity_max_relThresh", "max intensity rel. threshold", "degC");
    int v_ivrt = def_dvar("intensity_var_relThresh", "intensity var rel. threshold", "degC");
    int v_icrt = def_dvar("intensity_cumulative_relThresh", "cumulative intensity rel. threshold", "degC days");
    int v_ima = def_dvar("intensity_mean_abs", "mean absolute intensity", "degC");
    int v_ixa = def_dvar("intensity_max_abs", "max absolute intensity", "degC");
    int v_iva = def_dvar("intensity_var_abs", "absolute intensity variability", "degC");
    int v_ica = def_dvar("intensity_cumulative_abs", "cumulative absolute intensity", "degC days");
    int v_ro = def_dvar("rate_onset", "onset rate", "degC/day");
    int v_rd = def_dvar("rate_decline", "decline rate", "degC/day");

    // Reference date attribute
    int ref_y, ref_m, ref_d;
    jd_to_date(ref_date_jd, ref_y, ref_m, ref_d);
    char ref_str[32];
    snprintf(ref_str, sizeof(ref_str), "%04d-%02d-%02d", ref_y, ref_m, ref_d);
    std::string ref_units = std::string("days since ") + ref_str;
    nc_put_att_text(ncid, v_ds, "units", ref_units.size(), ref_units.c_str());
    nc_put_att_text(ncid, v_dp, "units", ref_units.size(), ref_units.c_str());
    nc_put_att_text(ncid, v_de, "units", ref_units.size(), ref_units.c_str());

    // Global attributes
    nc_put_att_text(ncid, NC_GLOBAL, "source_file",
                    source_file.size(), source_file.c_str());
    nc_put_att_text(ncid, NC_GLOBAL, "climatology_file",
                    clim_file.size(), clim_file.c_str());
    nc_put_att_int(ncid, NC_GLOBAL, "minDuration", NC_INT, 1, &minDuration);
    nc_put_att_int(ncid, NC_GLOBAL, "maxGap", NC_INT, 1, &maxGap);
    int cs = coldSpells ? 1 : 0;
    nc_put_att_int(ncid, NC_GLOBAL, "coldSpells", NC_INT, 1, &cs);
    const char* conv = "CF-1.8";
    nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, conv);
    const char* created_by = "heatwave3";
    nc_put_att_text(ncid, NC_GLOBAL, "created_by", 9, created_by);

    nc_check(nc_enddef(ncid), "enddef");

    // Write data
    nc_check(nc_put_var_double(ncid, v_lon, event_lon.data()), "put lon");
    nc_check(nc_put_var_double(ncid, v_lat, event_lat.data()), "put lat");
    nc_check(nc_put_var_int(ncid, v_pixel, pixel_index.data()), "put pixel");

    // Extract fields into contiguous arrays
    std::vector<int> eno(nevents), dur(nevents);
    std::vector<double> im(nevents), ix(nevents), iv(nevents), ic(nevents);
    std::vector<double> imrt(nevents), ixrt(nevents), ivrt(nevents), icrt(nevents);
    std::vector<double> ima(nevents), ixa(nevents), iva(nevents), ica(nevents);
    std::vector<double> ro(nevents), rd(nevents);

    for (size_t i = 0; i < nevents; ++i) {
        const auto& e = events[i];
        eno[i] = e.event_no;
        dur[i] = e.duration;
        im[i] = e.intensity_mean;
        ix[i] = e.intensity_max;
        iv[i] = e.intensity_var;
        ic[i] = e.intensity_cumulative;
        imrt[i] = e.intensity_mean_relThresh;
        ixrt[i] = e.intensity_max_relThresh;
        ivrt[i] = e.intensity_var_relThresh;
        icrt[i] = e.intensity_cumulative_relThresh;
        ima[i] = e.intensity_mean_abs;
        ixa[i] = e.intensity_max_abs;
        iva[i] = e.intensity_var_abs;
        ica[i] = e.intensity_cumulative_abs;
        ro[i] = e.rate_onset;
        rd[i] = e.rate_decline;
    }

    nc_check(nc_put_var_int(ncid, v_eno, eno.data()), "put event_no");
    nc_check(nc_put_var_int(ncid, v_ds, date_start.data()), "put date_start");
    nc_check(nc_put_var_int(ncid, v_dp, date_peak.data()), "put date_peak");
    nc_check(nc_put_var_int(ncid, v_de, date_end.data()), "put date_end");
    nc_check(nc_put_var_int(ncid, v_dur, dur.data()), "put duration");
    nc_check(nc_put_var_double(ncid, v_im, im.data()), "put intensity_mean");
    nc_check(nc_put_var_double(ncid, v_ix, ix.data()), "put intensity_max");
    nc_check(nc_put_var_double(ncid, v_iv, iv.data()), "put intensity_var");
    nc_check(nc_put_var_double(ncid, v_ic, ic.data()), "put intensity_cumulative");
    nc_check(nc_put_var_double(ncid, v_imrt, imrt.data()), "put imrt");
    nc_check(nc_put_var_double(ncid, v_ixrt, ixrt.data()), "put ixrt");
    nc_check(nc_put_var_double(ncid, v_ivrt, ivrt.data()), "put ivrt");
    nc_check(nc_put_var_double(ncid, v_icrt, icrt.data()), "put icrt");
    nc_check(nc_put_var_double(ncid, v_ima, ima.data()), "put ima");
    nc_check(nc_put_var_double(ncid, v_ixa, ixa.data()), "put ixa");
    nc_check(nc_put_var_double(ncid, v_iva, iva.data()), "put iva");
    nc_check(nc_put_var_double(ncid, v_ica, ica.data()), "put ica");
    nc_check(nc_put_var_double(ncid, v_ro, ro.data()), "put rate_onset");
    nc_check(nc_put_var_double(ncid, v_rd, rd.data()), "put rate_decline");

    nc_check(nc_close(ncid), "close " + file_out);
}

// ---- Read and merge multiple daily NetCDF files ----

GridData read_sst_multi_netcdf(const std::vector<std::string>& files,
                               const std::string& var_name,
                               const SubsetSpec& subset) {
    if (files.empty()) {
        throw std::runtime_error("No files provided to read_sst_multi_netcdf");
    }

    // Read first file to get spatial grid and metadata
    GridData first = read_sst_netcdf(files[0], var_name, subset);
    int nlon = first.nlon;
    int nlat = first.nlat;
    int npixels = nlon * nlat;

    // Collect all time steps and data
    struct TimeSlice {
        int jd;
        std::vector<double> data; // [npixels]
    };
    std::vector<TimeSlice> slices;
    slices.reserve(files.size());

    // Add slices from first file
    for (int t = 0; t < first.ntime; ++t) {
        TimeSlice ts;
        ts.jd = first.time_days[t];
        ts.data.resize(npixels);
        for (int px = 0; px < npixels; ++px) {
            ts.data[px] = first.sst[static_cast<size_t>(px) * first.ntime + t];
        }
        slices.push_back(std::move(ts));
    }

    // Read remaining files
    for (size_t f = 1; f < files.size(); ++f) {
        try {
            GridData gd = read_sst_netcdf(files[f], var_name, subset);

            // Verify spatial grid matches
            if (gd.nlon != nlon || gd.nlat != nlat) {
                Rcpp::Rcerr << "Warning: grid mismatch in " << files[f]
                            << " (" << gd.nlon << "x" << gd.nlat
                            << " vs " << nlon << "x" << nlat << "), skipping"
                            << std::endl;
                continue;
            }

            for (int t = 0; t < gd.ntime; ++t) {
                TimeSlice ts;
                ts.jd = gd.time_days[t];
                ts.data.resize(npixels);
                for (int px = 0; px < npixels; ++px) {
                    ts.data[px] = gd.sst[static_cast<size_t>(px) * gd.ntime + t];
                }
                slices.push_back(std::move(ts));
            }
        } catch (const std::exception& e) {
            Rcpp::Rcerr << "Warning: failed to read " << files[f]
                        << ": " << e.what() << ", skipping" << std::endl;
        }
    }

    // Sort by Julian Day
    std::sort(slices.begin(), slices.end(),
              [](const TimeSlice& a, const TimeSlice& b) { return a.jd < b.jd; });

    // Remove duplicate dates (keep first occurrence)
    {
        auto it = std::unique(slices.begin(), slices.end(),
                              [](const TimeSlice& a, const TimeSlice& b) {
                                  return a.jd == b.jd;
                              });
        slices.erase(it, slices.end());
    }

    // Assemble into GridData
    GridData result;
    result.lon = first.lon;
    result.lat = first.lat;
    result.nlon = nlon;
    result.nlat = nlat;
    result.ntime = static_cast<int>(slices.size());
    result.time_units = first.time_units;
    result.time_calendar = first.time_calendar;

    result.time_days.resize(result.ntime);
    for (int t = 0; t < result.ntime; ++t) {
        result.time_days[t] = slices[t].jd;
    }

    // Repack to pixel-major order: sst[pixel * ntime + t]
    size_t total = static_cast<size_t>(npixels) * result.ntime;
    result.sst.resize(total);
    for (int t = 0; t < result.ntime; ++t) {
        for (int px = 0; px < npixels; ++px) {
            result.sst[static_cast<size_t>(px) * result.ntime + t] = slices[t].data[px];
        }
    }

    return result;
}

} // namespace hw3
