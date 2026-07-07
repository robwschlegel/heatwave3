#include "netcdf_io.h"
#include <netcdf.h>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <regex>
#include <cctype>

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
    std::string cal = calendar;
    std::transform(cal.begin(), cal.end(), cal.begin(), ::tolower);
    if (cal.empty()) cal = "proleptic_gregorian";
    if (cal != "gregorian" && cal != "standard" && cal != "proleptic_gregorian") {
        throw std::runtime_error(
            "Unsupported CF calendar '" + calendar + "'. heatwave3 currently "
            "supports only gregorian, standard, and proleptic_gregorian calendars."
        );
    }

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
        // Map each timestamp to the calendar day that contains it (days since
        // the reference midnight), i.e. floor, matching heatwaveR's as.Date()
        // truncation. This keeps noon-stamped daily products (OSTIA/GHRSST,
        // stamped at 12:00:00) on the correct day rather than rounding them
        // forward. The small epsilon absorbs floating-point error at midnight.
        julian_days[i] = ref_jd + static_cast<int>(std::floor(time_raw[i] * scale + 1e-6));
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
    int lon_varid = -1, lat_varid = -1, time_varid = -1, depth_varid = -1;

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
            if (depth_dimid < 0) { depth_dimid = did; depth_varid = cvid; }
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

    // Depth coordinate (only present for 4D data). Read eagerly whenever the
    // dimension exists so both the legacy single-index squeeze and the new
    // depth_min/depth_max range mode below can use it.
    std::vector<double> full_depth;
    if (depth_dimid >= 0 && depth_varid >= 0) {
        size_t depth_len;
        nc_check(nc_inq_dimlen(ncid, depth_dimid, &depth_len), "inq depth len");
        full_depth.resize(depth_len);
        nc_check(nc_get_var_double(ncid, depth_varid, full_depth.data()), "get depth");
    }
    bool depth_range_requested = !std::isinf(subset.depth_min) || !std::isinf(subset.depth_max);
    if (depth_range_requested && depth_dimid >= 0 && full_depth.empty()) {
        nc_close(ncid);
        throw std::runtime_error("depth_range was requested but " + file_in +
                                 " has a depth dimension with no readable coordinate "
                                 "variable to select a range from.");
    }

    // Get time units
    if (!nc_get_att_str(ncid, time_varid, "units", gd.time_units)) {
        nc_close(ncid);
        throw std::runtime_error("Cannot read time units attribute from " + file_in);
    }

    if (!nc_get_att_str(ncid, time_varid, "calendar", gd.time_calendar)) {
        gd.time_calendar = "proleptic_gregorian";
    }
    if (!nc_get_att_str(ncid, varid, "units", gd.temp_units)) {
        gd.temp_units = "unknown";
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

    // Resolve the depth hyperslab window (only meaningful if depth_dimid >= 0):
    //   - depth_range_requested: contiguous index range covering
    //     [depth_min, depth_max] in the file's own depth coordinate values.
    //   - else legacy squeeze: subset.depth_index (or 0 by default), count 1.
    size_t depth_start = 0, depth_count = 1;
    if (depth_dimid >= 0) {
        if (depth_range_requested) {
            find_subset_indices(full_depth, subset.depth_min, subset.depth_max,
                                depth_start, depth_count);
        } else {
            depth_start = (subset.depth_index >= 0) ?
                          static_cast<size_t>(subset.depth_index) : 0;
            depth_count = 1;
        }
    }

    // Build hyperslab start/count arrays
    // Dimension ordering from the variable's dimids
    std::vector<size_t> start(ndims, 0), count(ndims, 1);
    for (int d = 0; d < ndims; ++d) {
        if (dimids[d] == lon_dimid)  { start[d] = lon_start;  count[d] = lon_count; }
        if (dimids[d] == lat_dimid)  { start[d] = lat_start;  count[d] = lat_count; }
        if (dimids[d] == time_dimid) { start[d] = time_start; count[d] = time_count; }
        if (dimids[d] == depth_dimid) {
            start[d] = depth_start;
            count[d] = depth_count;
        }
    }

    size_t total = lon_count * lat_count * time_count * depth_count;
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

    // Set the depth axis on the output before reordering (needed for the
    // pixel formula below).
    gd.ndepth = (depth_dimid >= 0) ? static_cast<int>(depth_count) : 1;
    if (depth_dimid >= 0 && !full_depth.empty()) {
        gd.depth.assign(full_depth.begin() + depth_start,
                        full_depth.begin() + depth_start + depth_count);
    } else {
        gd.depth.clear();
    }

    // Reorder to pixel-major: sst[pixel][time] where
    // pixel = (ilon * nlat + ilat) * ndepth + idepth (idepth always 0 when
    // ndepth == 1, reducing to the plain 3D layout).
    // Input is in file dimension order (e.g. [time, depth, lat, lon]); walk
    // it generically over up to 4 axes (lon/lat/time always present, depth
    // only when depth_dimid >= 0) so this works whether the file is 3D or 4D
    // and regardless of the file's own dimension ordering.
    std::vector<double> reordered(total);
    int eff_idx = 0;
    size_t eff_dims[4] = {1, 1, 1, 1};
    int dim_map[4] = {-1, -1, -1, -1}; // 0=lon, 1=lat, 2=time, 3=depth
    for (int d = 0; d < ndims; ++d) {
        eff_dims[eff_idx] = count[d];
        if (dimids[d] == lon_dimid)   dim_map[eff_idx] = 0;
        if (dimids[d] == lat_dimid)   dim_map[eff_idx] = 1;
        if (dimids[d] == time_dimid)  dim_map[eff_idx] = 2;
        if (dimids[d] == depth_dimid) dim_map[eff_idx] = 3;
        eff_idx++;
    }

    size_t total_eff = 1;
    for (int e = 0; e < eff_idx; ++e) total_eff *= eff_dims[e];

    for (size_t lin = 0; lin < total_eff; ++lin) {
        // Decompose the flat (row-major, last axis fastest) source index
        // into per-axis indices.
        size_t rem = lin;
        size_t idx[4] = {0, 0, 0, 0};
        for (int e = eff_idx - 1; e >= 0; --e) {
            idx[e] = rem % eff_dims[e];
            rem /= eff_dims[e];
        }

        size_t ilon = 0, ilat = 0, itime = 0, idepth = 0;
        for (int e = 0; e < eff_idx; ++e) {
            if (dim_map[e] == 0) ilon   = idx[e];
            if (dim_map[e] == 1) ilat   = idx[e];
            if (dim_map[e] == 2) itime  = idx[e];
            if (dim_map[e] == 3) idepth = idx[e];
        }
        size_t pixel = (ilon * static_cast<size_t>(gd.nlat) + ilat) *
                       static_cast<size_t>(gd.ndepth) + idepth;
        size_t dst = pixel * static_cast<size_t>(gd.ntime) + itime;
        reordered[dst] = gd.sst[lin];
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
                       int smoothPercentileWidth,
                       const std::string& temp_units,
                       int ndepth,
                       const std::vector<double>& depth) {
    int ncid;
    nc_check(nc_create(file_out.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid),
             "create " + file_out);

    bool has_depth = (ndepth > 1) && !depth.empty();

    // Define dimensions
    int lon_dimid, lat_dimid, doy_dimid, depth_dimid = -1;
    nc_check(nc_def_dim(ncid, "lon", nlon, &lon_dimid), "def dim lon");
    nc_check(nc_def_dim(ncid, "lat", nlat, &lat_dimid), "def dim lat");
    nc_check(nc_def_dim(ncid, "doy", 366, &doy_dimid), "def dim doy");
    if (has_depth) {
        nc_check(nc_def_dim(ncid, "depth", ndepth, &depth_dimid), "def dim depth");
    }

    // Define coordinate variables
    int lon_varid, lat_varid, doy_varid, depth_varid = -1;
    nc_check(nc_def_var(ncid, "lon", NC_DOUBLE, 1, &lon_dimid, &lon_varid), "def var lon");
    nc_check(nc_def_var(ncid, "lat", NC_DOUBLE, 1, &lat_dimid, &lat_varid), "def var lat");
    nc_check(nc_def_var(ncid, "doy", NC_INT, 1, &doy_dimid, &doy_varid), "def var doy");
    if (has_depth) {
        nc_check(nc_def_var(ncid, "depth", NC_DOUBLE, 1, &depth_dimid, &depth_varid), "def var depth");
        nc_put_att_text(ncid, depth_varid, "units", 1, "m");
        nc_put_att_text(ncid, depth_varid, "positive", 4, "down");
        nc_put_att_text(ncid, depth_varid, "axis", 1, "Z");
        nc_put_att_text(ncid, depth_varid, "standard_name", 5, "depth");
    }

    nc_put_att_text(ncid, lon_varid, "units", 12, "degrees_east");
    nc_put_att_text(ncid, lon_varid, "axis", 1, "X");
    nc_put_att_text(ncid, lat_varid, "units", 13, "degrees_north");
    nc_put_att_text(ncid, lat_varid, "axis", 1, "Y");
    nc_put_att_text(ncid, doy_varid, "long_name", 11, "day_of_year");

    // Define data variables: [lon, lat, doy] (3D) or [lon, lat, depth, doy] (4D)
    int dims3[3] = {lon_dimid, lat_dimid, doy_dimid};
    int dims4[4] = {lon_dimid, lat_dimid, depth_dimid, doy_dimid};
    int ndims_data = has_depth ? 4 : 3;
    int* dims_data = has_depth ? dims4 : dims3;
    int seas_varid, thresh_varid, var_varid;
    nc_check(nc_def_var(ncid, "seas", NC_DOUBLE, ndims_data, dims_data, &seas_varid), "def var seas");
    nc_check(nc_def_var(ncid, "thresh", NC_DOUBLE, ndims_data, dims_data, &thresh_varid), "def var thresh");

    nc_put_att_text(ncid, seas_varid, "long_name", 21, "seasonal climatology");
    nc_put_att_text(ncid, seas_varid, "units", temp_units.size(), temp_units.c_str());
    nc_put_att_text(ncid, thresh_varid, "long_name", 30, "threshold climatology");
    nc_put_att_text(ncid, thresh_varid, "units", temp_units.size(), temp_units.c_str());

    // Enable compression
    nc_def_var_deflate(ncid, seas_varid, 1, 1, 4);
    nc_def_var_deflate(ncid, thresh_varid, 1, 1, 4);

    if (has_var) {
        nc_check(nc_def_var(ncid, "var", NC_DOUBLE, ndims_data, dims_data, &var_varid), "def var var");
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
    nc_put_att_text(ncid, NC_GLOBAL, "hw3_product", 11, "climatology");

    nc_check(nc_enddef(ncid), "enddef");

    // Write coordinates
    nc_check(nc_put_var_double(ncid, lon_varid, lon.data()), "put lon");
    nc_check(nc_put_var_double(ncid, lat_varid, lat.data()), "put lat");
    if (has_depth) {
        nc_check(nc_put_var_double(ncid, depth_varid, depth.data()), "put depth");
    }
    std::vector<int> doy_vals(366);
    for (int i = 0; i < 366; ++i) doy_vals[i] = i + 1;
    nc_check(nc_put_var_int(ncid, doy_varid, doy_vals.data()), "put doy");

    // Write data — input layout is [pixel][doy] where
    // pixel = (ilon*nlat + ilat)*ndepth + idepth, matching NetCDF layout
    // [lon][lat][doy] (3D) or [lon][lat][depth][doy] (4D). Same order in
    // both cases, so write directly with no reordering.
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

    // Optional depth dimension (depth-resolved climatology written with
    // ndepth > 1). Absent in ordinary 3D climatology files.
    int depth_dimid;
    if (nc_inq_dimid(ncid, "depth", &depth_dimid) == NC_NOERR) {
        size_t ndepth;
        nc_check(nc_inq_dimlen(ncid, depth_dimid, &ndepth), "inq ndepth");
        cd.ndepth = static_cast<int>(ndepth);
        int depth_varid;
        nc_check(nc_inq_varid(ncid, "depth", &depth_varid), "inq depth");
        cd.depth.resize(ndepth);
        nc_check(nc_get_var_double(ncid, depth_varid, cd.depth.data()), "get depth");
    }

    size_t total = nlon * nlat * static_cast<size_t>(cd.ndepth) * ndoy;
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

// ---- Read event NetCDF (ragged array) ----

EventData read_event_netcdf(const std::string& event_file) {
    EventData ed;
    int ncid;
    nc_check(nc_open(event_file.c_str(), NC_NOWRITE, &ncid), "open " + event_file);

    int event_dimid;
    nc_check(nc_inq_dimid(ncid, "event", &event_dimid), "inq dim event");
    size_t nevents;
    nc_check(nc_inq_dimlen(ncid, event_dimid, &nevents), "inq nevents");
    ed.nevents = static_cast<int>(nevents);

    auto read_dvec = [&](const char* name, std::vector<double>& v) {
        int vid;
        nc_check(nc_inq_varid(ncid, name, &vid), std::string("inq ") + name);
        v.resize(nevents);
        nc_check(nc_get_var_double(ncid, vid, v.data()), std::string("get ") + name);
    };
    auto read_ivec = [&](const char* name, std::vector<int>& v) {
        int vid;
        nc_check(nc_inq_varid(ncid, name, &vid), std::string("inq ") + name);
        v.resize(nevents);
        nc_check(nc_get_var_int(ncid, vid, v.data()), std::string("get ") + name);
    };

    read_dvec("lon", ed.lon);
    read_dvec("lat", ed.lat);
    read_ivec("pixel_index", ed.pixel_index);
    read_ivec("event_no", ed.event_no);
    read_ivec("date_start", ed.date_start);
    read_ivec("date_peak", ed.date_peak);
    read_ivec("date_end", ed.date_end);
    read_ivec("duration", ed.duration);
    read_dvec("intensity_mean", ed.intensity_mean);
    read_dvec("intensity_max", ed.intensity_max);
    read_dvec("intensity_var", ed.intensity_var);
    read_dvec("intensity_cumulative", ed.intensity_cumulative);
    read_dvec("intensity_mean_relThresh", ed.intensity_mean_relThresh);
    read_dvec("intensity_max_relThresh", ed.intensity_max_relThresh);
    read_dvec("intensity_var_relThresh", ed.intensity_var_relThresh);
    read_dvec("intensity_cumulative_relThresh", ed.intensity_cumulative_relThresh);
    read_dvec("intensity_mean_abs", ed.intensity_mean_abs);
    read_dvec("intensity_max_abs", ed.intensity_max_abs);
    read_dvec("intensity_var_abs", ed.intensity_var_abs);
    read_dvec("intensity_cumulative_abs", ed.intensity_cumulative_abs);
    read_dvec("rate_onset", ed.rate_onset);
    read_dvec("rate_decline", ed.rate_decline);

    // Optional category fields (may not exist in older files)
    auto try_read_ivec = [&](const char* name, std::vector<int>& v) {
        int vid;
        if (nc_inq_varid(ncid, name, &vid) == NC_NOERR) {
            v.resize(nevents);
            nc_get_var_int(ncid, vid, v.data());
        }
    };
    auto try_read_dvec = [&](const char* name, std::vector<double>& v) {
        int vid;
        if (nc_inq_varid(ncid, name, &vid) == NC_NOERR) {
            v.resize(nevents);
            nc_get_var_double(ncid, vid, v.data());
        }
    };

    try_read_ivec("category", ed.category);
    try_read_dvec("p_moderate", ed.p_moderate);
    try_read_dvec("p_strong", ed.p_strong);
    try_read_dvec("p_severe", ed.p_severe);
    try_read_dvec("p_extreme", ed.p_extreme);
    try_read_ivec("season", ed.season);
    try_read_dvec("depth", ed.depth);

    // Parse reference date from date_start units attribute
    int ds_vid;
    nc_check(nc_inq_varid(ncid, "date_start", &ds_vid), "inq date_start");
    std::string units_str;
    nc_get_att_str(ncid, ds_vid, "units", units_str);
    // Parse "days since YYYY-MM-DD"
    std::regex re(R"(days\s+since\s+(\d{4})-(\d{1,2})-(\d{1,2}))");
    std::smatch match;
    if (std::regex_search(units_str, match, re)) {
        int ry = std::stoi(match[1].str());
        int rm = std::stoi(match[2].str());
        int rd = std::stoi(match[3].str());
        ed.ref_date_jd = date_to_jd(ry, rm, rd);
    } else {
        ed.ref_date_jd = 0;
    }

    nc_close(ncid);
    return ed;
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
                        bool coldSpells,
                        bool southHemisphere,
                        const std::string& temp_units,
                        const std::vector<double>& event_depth) {
    int ncid;
    nc_check(nc_create(file_out.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid),
             "create " + file_out);

    size_t nevents = events.size();
    bool has_depth = !event_depth.empty();

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
    std::string cumulative_units = temp_units + " days";
    std::string rate_units = temp_units + "/day";
    int v_im = def_dvar("intensity_mean", "mean intensity", temp_units.c_str());
    int v_ix = def_dvar("intensity_max", "maximum intensity", temp_units.c_str());
    int v_iv = def_dvar("intensity_var", "intensity variability", temp_units.c_str());
    int v_ic = def_dvar("intensity_cumulative", "cumulative intensity", cumulative_units.c_str());
    int v_imrt = def_dvar("intensity_mean_relThresh", "mean intensity rel. threshold", temp_units.c_str());
    int v_ixrt = def_dvar("intensity_max_relThresh", "max intensity rel. threshold", temp_units.c_str());
    int v_ivrt = def_dvar("intensity_var_relThresh", "intensity var rel. threshold", temp_units.c_str());
    int v_icrt = def_dvar("intensity_cumulative_relThresh", "cumulative intensity rel. threshold", cumulative_units.c_str());
    int v_ima = def_dvar("intensity_mean_abs", "mean absolute intensity", temp_units.c_str());
    int v_ixa = def_dvar("intensity_max_abs", "max absolute intensity", temp_units.c_str());
    int v_iva = def_dvar("intensity_var_abs", "absolute intensity variability", temp_units.c_str());
    int v_ica = def_dvar("intensity_cumulative_abs", "cumulative absolute intensity", cumulative_units.c_str());
    int v_ro = def_dvar("rate_onset", "onset rate", rate_units.c_str());
    int v_rd = def_dvar("rate_decline", "decline rate", rate_units.c_str());
    int v_cat = def_ivar("category", "event category", nullptr);
    int v_pm = def_dvar("p_moderate", "proportion moderate", nullptr);
    int v_ps = def_dvar("p_strong", "proportion strong", nullptr);
    int v_pv = def_dvar("p_severe", "proportion severe", nullptr);
    int v_pe = def_dvar("p_extreme", "proportion extreme", nullptr);
    int v_sea = def_ivar("season", "peak season", nullptr);
    int v_depth = has_depth ? def_dvar("depth", "depth", "m") : -1;

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
    const char* hemi = southHemisphere ? "south" : "north";
    nc_put_att_text(ncid, NC_GLOBAL, "hemisphere", strlen(hemi), hemi);
    nc_put_att_text(ncid, NC_GLOBAL, "hw3_product", 6, "events");

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
    std::vector<int> cat(nevents), sea(nevents);
    std::vector<double> pm(nevents), ps(nevents), pv2(nevents), pex(nevents);

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
        cat[i] = e.category;
        pm[i] = e.p_moderate;
        ps[i] = e.p_strong;
        pv2[i] = e.p_severe;
        pex[i] = e.p_extreme;
        sea[i] = e.season;
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
    nc_check(nc_put_var_int(ncid, v_cat, cat.data()), "put category");
    nc_check(nc_put_var_double(ncid, v_pm, pm.data()), "put p_moderate");
    nc_check(nc_put_var_double(ncid, v_ps, ps.data()), "put p_strong");
    nc_check(nc_put_var_double(ncid, v_pv, pv2.data()), "put p_severe");
    nc_check(nc_put_var_double(ncid, v_pe, pex.data()), "put p_extreme");
    nc_check(nc_put_var_int(ncid, v_sea, sea.data()), "put season");
    if (has_depth) {
        nc_check(nc_put_var_double(ncid, v_depth, event_depth.data()), "put depth");
    }

    nc_check(nc_close(ncid), "close " + file_out);
}

// ---- Write per-day NetCDF ([lon, lat, time]) ----

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
                        const std::string& temp_units) {
    int ncid;
    nc_check(nc_create(file_out.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid),
             "create " + file_out);

    int lon_dimid, lat_dimid, time_dimid;
    nc_check(nc_def_dim(ncid, "lon", nlon, &lon_dimid), "def dim lon");
    nc_check(nc_def_dim(ncid, "lat", nlat, &lat_dimid), "def dim lat");
    nc_check(nc_def_dim(ncid, "time", ntime, &time_dimid), "def dim time");

    int lon_varid, lat_varid, time_varid;
    nc_check(nc_def_var(ncid, "lon", NC_DOUBLE, 1, &lon_dimid, &lon_varid), "def var lon");
    nc_check(nc_def_var(ncid, "lat", NC_DOUBLE, 1, &lat_dimid, &lat_varid), "def var lat");
    nc_check(nc_def_var(ncid, "time", NC_DOUBLE, 1, &time_dimid, &time_varid), "def var time");
    nc_put_att_text(ncid, lon_varid, "units", 12, "degrees_east");
    nc_put_att_text(ncid, lon_varid, "axis", 1, "X");
    nc_put_att_text(ncid, lat_varid, "units", 13, "degrees_north");
    nc_put_att_text(ncid, lat_varid, "axis", 1, "Y");
    // Time axis as "days since 1970-01-01" (jd 2440588), derived from the
    // Julian-day vector so it works identically for single- and multi-file
    // inputs (whose raw time units may differ).
    const char* time_units = "days since 1970-01-01";
    const char* time_calendar = "standard";
    nc_put_att_text(ncid, time_varid, "units", strlen(time_units), time_units);
    nc_put_att_text(ncid, time_varid, "calendar", strlen(time_calendar), time_calendar);
    nc_put_att_text(ncid, time_varid, "axis", 1, "T");

    int dims3[3] = {lon_dimid, lat_dimid, time_dimid};

    auto def_fvar = [&](const char* name, const char* long_name, const char* units) -> int {
        int vid;
        nc_check(nc_def_var(ncid, name, NC_FLOAT, 3, dims3, &vid), name);
        nc_put_att_text(ncid, vid, "long_name", strlen(long_name), long_name);
        if (units) nc_put_att_text(ncid, vid, "units", strlen(units), units);
        nc_def_var_deflate(ncid, vid, 1, 1, 4);
        return vid;
    };
    auto def_bvar = [&](const char* name, const char* long_name) -> int {
        int vid;
        nc_check(nc_def_var(ncid, name, NC_BYTE, 3, dims3, &vid), name);
        nc_put_att_text(ncid, vid, "long_name", strlen(long_name), long_name);
        nc_def_var_deflate(ncid, vid, 1, 1, 4);
        return vid;
    };
    auto def_ivar = [&](const char* name, const char* long_name) -> int {
        int vid;
        nc_check(nc_def_var(ncid, name, NC_INT, 3, dims3, &vid), name);
        nc_put_att_text(ncid, vid, "long_name", strlen(long_name), long_name);
        nc_def_var_deflate(ncid, vid, 1, 1, 4);
        return vid;
    };

    int v_temp = -1, v_seas = -1, v_thresh = -1, v_tc = -1, v_dc = -1;
    int v_ev = -1, v_eno = -1, v_int = -1, v_cat = -1;
    if (temp)   v_temp   = def_fvar("temp", "sea water temperature", temp_units.c_str());
    if (seas)   v_seas   = def_fvar("seas", "seasonal climatology", temp_units.c_str());
    if (thresh) v_thresh = def_fvar("thresh", "threshold climatology", temp_units.c_str());
    if (threshCriterion)
        v_tc = def_bvar("threshCriterion", "temperature exceeds threshold (1=yes)");
    if (durationCriterion)
        v_dc = def_bvar("durationCriterion", "within a run of at least minDuration days (1=yes)");
    if (event)  v_ev = def_bvar("event", "day within a detected event (1=yes)");
    if (event_no) v_eno = def_ivar("event_no", "event number within pixel (0=none)");
    if (intensity)
        v_int = def_fvar("intensity", "temperature anomaly (temp - seas)", temp_units.c_str());
    if (category)
        v_cat = def_bvar("category", "daily Hobday category (0=none, 1=I ... 4=IV)");

    nc_put_att_text(ncid, NC_GLOBAL, "source_file", source_file.size(), source_file.c_str());
    nc_put_att_text(ncid, NC_GLOBAL, "climatology_file", clim_file.size(), clim_file.c_str());
    nc_put_att_int(ncid, NC_GLOBAL, "minDuration", NC_INT, 1, &minDuration);
    nc_put_att_int(ncid, NC_GLOBAL, "maxGap", NC_INT, 1, &maxGap);
    int cs = coldSpells ? 1 : 0;
    nc_put_att_int(ncid, NC_GLOBAL, "coldSpells", NC_INT, 1, &cs);
    const char* conv = "CF-1.8";
    nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 6, conv);
    const char* created_by = "heatwave3";
    nc_put_att_text(ncid, NC_GLOBAL, "created_by", 9, created_by);
    const char* hemi = southHemisphere ? "south" : "north";
    nc_put_att_text(ncid, NC_GLOBAL, "hemisphere", strlen(hemi), hemi);
    nc_put_att_text(ncid, NC_GLOBAL, "hw3_product", product.size(), product.c_str());

    nc_check(nc_enddef(ncid), "enddef");

    nc_check(nc_put_var_double(ncid, lon_varid, lon.data()), "put lon");
    nc_check(nc_put_var_double(ncid, lat_varid, lat.data()), "put lat");
    std::vector<double> time_vals(ntime);
    for (int t = 0; t < ntime; ++t)
        time_vals[t] = static_cast<double>(time_days[t] - 2440588);  // jd -> days since 1970
    nc_check(nc_put_var_double(ncid, time_varid, time_vals.data()), "put time");

    // Data buffers are pixel-major [pixel][time] = [lon][lat][time], which is
    // exactly the declared dimension order, so write directly. Real-valued
    // variables are NC_FLOAT; nc_put_var_double converts (NaN stays NaN -> NA).
    if (temp)   nc_check(nc_put_var_double(ncid, v_temp, temp), "put temp");
    if (seas)   nc_check(nc_put_var_double(ncid, v_seas, seas), "put seas");
    if (thresh) nc_check(nc_put_var_double(ncid, v_thresh, thresh), "put thresh");
    if (threshCriterion)
        nc_check(nc_put_var_schar(ncid, v_tc, threshCriterion), "put threshCriterion");
    if (durationCriterion)
        nc_check(nc_put_var_schar(ncid, v_dc, durationCriterion), "put durationCriterion");
    if (event)  nc_check(nc_put_var_schar(ncid, v_ev, event), "put event");
    if (event_no) nc_check(nc_put_var_int(ncid, v_eno, event_no), "put event_no");
    if (intensity) nc_check(nc_put_var_double(ncid, v_int, intensity), "put intensity");
    if (category) nc_check(nc_put_var_schar(ncid, v_cat, category), "put category");

    nc_check(nc_close(ncid), "close " + file_out);
}

// ---- Streaming hyperslab subset of a gridded product ----

Rcpp::List read_subset_netcdf(const std::string& file,
                              const std::vector<double>& lon_range,
                              const std::vector<double>& lat_range,
                              const std::vector<int>& t_jd_range,
                              const std::vector<std::string>& vars,
                              int max_rows) {
    int ncid;
    nc_check(nc_open(file.c_str(), NC_NOWRITE, &ncid), "open " + file);

    std::string product;
    nc_get_att_str(ncid, NC_GLOBAL, "hw3_product", product);

    auto has_dim = [&](const char* name) {
        int d; return nc_inq_dimid(ncid, name, &d) == NC_NOERR;
    };
    if (product.empty()) {
        if (has_dim("event")) product = "events";
        else if (has_dim("doy")) product = "climatology";
        else if (has_dim("time")) product = "daily";
    }
    if (product == "events") {
        nc_close(ncid);
        Rcpp::stop("read_subset_netcdf: events are subset in R, not via hyperslab.");
    }
    bool is_clim = (product == "climatology");
    if (is_clim && has_dim("depth")) {
        nc_close(ncid);
        Rcpp::stop("hw3_export()'s streaming subset (vars/lon_range/lat_range/"
                   "time_range/n) is not yet supported for a depth-resolved "
                   "(depth_range) climatology. Call hw3_export(file) with no "
                   "subset arguments for a full read instead.");
    }
    const char* third_dim = is_clim ? "doy" : "time";

    auto get_dimlen = [&](const char* name) -> size_t {
        int d; nc_check(nc_inq_dimid(ncid, name, &d), std::string("dim ") + name);
        size_t n; nc_check(nc_inq_dimlen(ncid, d, &n), std::string("dimlen ") + name);
        return n;
    };
    size_t nlon = get_dimlen("lon");
    size_t nlat = get_dimlen("lat");
    size_t n3   = get_dimlen(third_dim);

    std::vector<double> lon(nlon), lat(nlat);
    int vid_lon, vid_lat;
    nc_check(nc_inq_varid(ncid, "lon", &vid_lon), "inq lon");
    nc_check(nc_inq_varid(ncid, "lat", &vid_lat), "inq lat");
    nc_check(nc_get_var_double(ncid, vid_lon, lon.data()), "get lon");
    nc_check(nc_get_var_double(ncid, vid_lat, lat.data()), "get lat");

    // Third coordinate values: DOY 1..366 for climatology, else Julian Days.
    std::vector<int> third_vals(n3);
    if (is_clim) {
        for (size_t k = 0; k < n3; ++k) third_vals[k] = static_cast<int>(k) + 1;
    } else {
        int vid_t; nc_check(nc_inq_varid(ncid, "time", &vid_t), "inq time");
        std::vector<double> traw(n3);
        nc_check(nc_get_var_double(ncid, vid_t, traw.data()), "get time");
        std::string tunits, tcal;
        nc_get_att_str(ncid, vid_t, "units", tunits);
        nc_get_att_str(ncid, vid_t, "calendar", tcal);
        if (tcal.empty()) tcal = "standard";
        parse_cf_time(tunits, tcal, traw, third_vals);
    }

    // Contiguous index window from a monotonic coordinate and [min,max].
    auto coord_window = [](const std::vector<double>& c,
                           const std::vector<double>& rng,
                           int& i0, int& i1) {
        int n = static_cast<int>(c.size());
        if (rng.size() != 2) { i0 = 0; i1 = n - 1; return; }
        double lo = std::min(rng[0], rng[1]), hi = std::max(rng[0], rng[1]);
        i0 = n; i1 = -1;
        for (int i = 0; i < n; ++i) {
            if (c[i] >= lo && c[i] <= hi) { if (i < i0) i0 = i; if (i > i1) i1 = i; }
        }
    };

    int li0, li1, lj0, lj1;
    coord_window(lon, lon_range, li0, li1);
    coord_window(lat, lat_range, lj0, lj1);

    int t0 = 0, t1 = static_cast<int>(n3) - 1;
    if (!is_clim && t_jd_range.size() == 2) {
        int lo = std::min(t_jd_range[0], t_jd_range[1]);
        int hi = std::max(t_jd_range[0], t_jd_range[1]);
        t0 = static_cast<int>(n3); t1 = -1;
        for (int k = 0; k < static_cast<int>(n3); ++k) {
            if (third_vals[k] >= lo && third_vals[k] <= hi) {
                if (k < t0) t0 = k; if (k > t1) t1 = k;
            }
        }
    }

    // Candidate data variables per product (read those that are present).
    std::vector<std::string> avail = is_clim
        ? std::vector<std::string>{"seas", "thresh", "var"}
        : std::vector<std::string>{"temp", "seas", "thresh", "intensity",
                                    "threshCriterion", "durationCriterion",
                                    "event", "event_no", "category"};
    std::vector<std::string> present;
    for (const auto& nm : avail) {
        int v; if (nc_inq_varid(ncid, nm.c_str(), &v) == NC_NOERR) present.push_back(nm);
    }
    std::vector<std::string> wanted = vars.empty() ? present : vars;
    for (const auto& w : wanted) {
        if (std::find(present.begin(), present.end(), w) == present.end()) {
            std::string avs;
            for (size_t i = 0; i < present.size(); ++i)
                avs += (i ? ", " : "") + present[i];
            nc_close(ncid);
            Rcpp::stop("Variable '%s' is not in this %s file. Available: %s",
                       w, product, avs);
        }
    }

    // Empty window (no overlap) -> return a zero-row result.
    bool empty = (li1 < li0) || (lj1 < lj0) || (t1 < t0);
    int nlon_sub = empty ? 0 : (li1 - li0 + 1);
    int nlat_sub = empty ? 0 : (lj1 - lj0 + 1);
    int n3_sub   = empty ? 0 : (t1 - t0 + 1);

    // Row cap: limit how many longitude columns we read.
    int nlon_read = nlon_sub;
    if (!empty && max_rows >= 0) {
        long long per = static_cast<long long>(nlat_sub) * n3_sub; // rows per lon column
        if (per > 0) {
            long long lon_need = (static_cast<long long>(max_rows) + per - 1) / per;
            if (lon_need < nlon_read) nlon_read = std::max(1, static_cast<int>(lon_need));
        }
    }

    Rcpp::List data;
    if (!empty) {
        size_t start[3] = {static_cast<size_t>(li0), static_cast<size_t>(lj0),
                           static_cast<size_t>(t0)};
        size_t count[3] = {static_cast<size_t>(nlon_read),
                           static_cast<size_t>(nlat_sub),
                           static_cast<size_t>(n3_sub)};
        size_t slab = static_cast<size_t>(nlon_read) * nlat_sub * n3_sub;
        for (const auto& w : wanted) {
            int vid; nc_check(nc_inq_varid(ncid, w.c_str(), &vid), "inq " + w);
            nc_type vt; nc_check(nc_inq_vartype(ncid, vid, &vt), "type " + w);
            if (vt == NC_FLOAT || vt == NC_DOUBLE) {
                std::vector<double> buf(slab);
                nc_check(nc_get_vara_double(ncid, vid, start, count, buf.data()), "get " + w);
                data[w] = buf;
            } else if (vt == NC_INT) {
                std::vector<int> buf(slab);
                nc_check(nc_get_vara_int(ncid, vid, start, count, buf.data()), "get " + w);
                data[w] = buf;
            } else { // NC_BYTE flags
                std::vector<signed char> b(slab);
                nc_check(nc_get_vara_schar(ncid, vid, start, count, b.data()), "get " + w);
                data[w] = std::vector<int>(b.begin(), b.end());
            }
        }
    }
    nc_close(ncid);

    std::vector<double> lon_s, lat_s;
    std::vector<int> third_s;
    for (int i = 0; i < nlon_read; ++i) lon_s.push_back(lon[li0 + i]);
    for (int j = 0; j < nlat_sub; ++j)  lat_s.push_back(lat[lj0 + j]);
    for (int k = 0; k < n3_sub; ++k)    third_s.push_back(third_vals[t0 + k]);

    return Rcpp::List::create(
        Rcpp::Named("product") = product,
        Rcpp::Named("lon") = lon_s,
        Rcpp::Named("lat") = lat_s,
        Rcpp::Named("third") = third_s,
        Rcpp::Named("third_name") = std::string(is_clim ? "doy" : "t"),
        Rcpp::Named("nlon") = nlon_read,
        Rcpp::Named("nlat") = nlat_sub,
        Rcpp::Named("n3") = n3_sub,
        Rcpp::Named("data") = data
    );
}

// ---- Lightweight product metadata (attrs + dim lengths only) ----

FileMeta read_file_meta(const std::string& file) {
    FileMeta m;
    int ncid;
    nc_check(nc_open(file.c_str(), NC_NOWRITE, &ncid), "open " + file);

    nc_get_att_str(ncid, NC_GLOBAL, "hw3_product", m.product);

    auto dimlen = [&](const char* name, size_t& out) -> bool {
        int did;
        if (nc_inq_dimid(ncid, name, &did) != NC_NOERR) return false;
        return nc_inq_dimlen(ncid, did, &out) == NC_NOERR;
    };

    size_t nlon = 0, nlat = 0, ndoy = 0, ntime = 0, nevent = 0, ndepth = 1;
    bool has_event = dimlen("event", nevent);
    bool has_doy = dimlen("doy", ndoy);
    bool has_time = dimlen("time", ntime);
    dimlen("lon", nlon);
    dimlen("lat", nlat);
    // Depth-resolved climatology (ts2clm3(depth_range=...)): fold into nrows.
    size_t depth_dimlen;
    if (dimlen("depth", depth_dimlen)) ndepth = depth_dimlen;

    // Infer product for files written before the hw3_product attribute existed.
    if (m.product.empty()) {
        if (has_event) m.product = "events";
        else if (has_doy) m.product = "climatology";
        else if (has_time) m.product = "daily";
    }

    if (m.product == "events") {
        m.n3 = static_cast<int>(nevent);
        m.nrows = static_cast<long long>(nevent);
    } else if (m.product == "climatology") {
        m.nlon = static_cast<int>(nlon);
        m.nlat = static_cast<int>(nlat);
        m.n3 = static_cast<int>(ndoy);
        m.nrows = static_cast<long long>(nlon) * nlat * ndepth * ndoy;
    } else { // daily / protoevents
        m.nlon = static_cast<int>(nlon);
        m.nlat = static_cast<int>(nlat);
        m.n3 = static_cast<int>(ntime);
        m.nrows = static_cast<long long>(nlon) * nlat * ntime;
    }

    nc_close(ncid);
    return m;
}

// ---- Read per-day NetCDF ----

DailyData read_daily_netcdf(const std::string& daily_file) {
    DailyData dd;
    int ncid;
    nc_check(nc_open(daily_file.c_str(), NC_NOWRITE, &ncid), "open " + daily_file);

    int lon_dimid, lat_dimid, time_dimid;
    nc_check(nc_inq_dimid(ncid, "lon", &lon_dimid), "inq dim lon");
    nc_check(nc_inq_dimid(ncid, "lat", &lat_dimid), "inq dim lat");
    nc_check(nc_inq_dimid(ncid, "time", &time_dimid), "inq dim time");

    size_t nlon, nlat, ntime;
    nc_check(nc_inq_dimlen(ncid, lon_dimid, &nlon), "inq nlon");
    nc_check(nc_inq_dimlen(ncid, lat_dimid, &nlat), "inq nlat");
    nc_check(nc_inq_dimlen(ncid, time_dimid, &ntime), "inq ntime");
    dd.nlon = static_cast<int>(nlon);
    dd.nlat = static_cast<int>(nlat);
    dd.ntime = static_cast<int>(ntime);

    int lon_varid, lat_varid, time_varid;
    nc_check(nc_inq_varid(ncid, "lon", &lon_varid), "inq lon");
    nc_check(nc_inq_varid(ncid, "lat", &lat_varid), "inq lat");
    nc_check(nc_inq_varid(ncid, "time", &time_varid), "inq time");

    dd.lon.resize(nlon);
    dd.lat.resize(nlat);
    nc_check(nc_get_var_double(ncid, lon_varid, dd.lon.data()), "get lon");
    nc_check(nc_get_var_double(ncid, lat_varid, dd.lat.data()), "get lat");

    std::vector<double> time_raw(ntime);
    nc_check(nc_get_var_double(ncid, time_varid, time_raw.data()), "get time");
    std::string tunits, tcal;
    nc_get_att_str(ncid, time_varid, "units", tunits);
    nc_get_att_str(ncid, time_varid, "calendar", tcal);
    if (tcal.empty()) tcal = "standard";
    parse_cf_time(tunits, tcal, time_raw, dd.time_jd);

    size_t total = nlon * nlat * ntime;
    auto rd_f = [&](const char* name, std::vector<double>& v) {
        int vid;
        if (nc_inq_varid(ncid, name, &vid) == NC_NOERR) {
            v.resize(total);
            nc_check(nc_get_var_double(ncid, vid, v.data()), std::string("get ") + name);
        }
    };
    auto rd_b = [&](const char* name, std::vector<signed char>& v) {
        int vid;
        if (nc_inq_varid(ncid, name, &vid) == NC_NOERR) {
            v.resize(total);
            nc_check(nc_get_var_schar(ncid, vid, v.data()), std::string("get ") + name);
        }
    };
    auto rd_i = [&](const char* name, std::vector<int>& v) {
        int vid;
        if (nc_inq_varid(ncid, name, &vid) == NC_NOERR) {
            v.resize(total);
            nc_check(nc_get_var_int(ncid, vid, v.data()), std::string("get ") + name);
        }
    };

    rd_f("temp", dd.temp);
    rd_f("seas", dd.seas);
    rd_f("thresh", dd.thresh);
    rd_f("intensity", dd.intensity);
    rd_b("threshCriterion", dd.threshCriterion);
    rd_b("durationCriterion", dd.durationCriterion);
    rd_b("event", dd.event);
    rd_b("category", dd.category);
    rd_i("event_no", dd.event_no);

    nc_close(ncid);
    return dd;
}

// ---- Read and merge multiple daily NetCDF files ----

GridData read_sst_multi_netcdf(const std::vector<std::string>& files,
                               const std::string& var_name,
                               const SubsetSpec& subset,
                               bool skip_bad_files) {
    if (files.empty()) {
        throw std::runtime_error("No files provided to read_sst_multi_netcdf");
    }

    // Read first usable file to get spatial grid and metadata.
    std::string read_var_name = var_name;
    size_t first_file_index = 0;
    GridData first;
    bool have_first = false;
    for (; first_file_index < files.size(); ++first_file_index) {
        try {
            std::string file_var_name = read_var_name;
            if (file_var_name.empty()) {
                file_var_name = detect_sst_variable(files[first_file_index]);
            }
            first = read_sst_netcdf(files[first_file_index], file_var_name, subset);
            if (read_var_name.empty()) {
                read_var_name = file_var_name;
            }
            have_first = true;
            break;
        } catch (const std::exception& e) {
            if (!skip_bad_files) {
                throw;
            }
            Rcpp::Rcerr << "Warning: failed to read " << files[first_file_index]
                        << ": " << e.what() << ", skipping" << std::endl;
        }
    }
    if (!have_first) {
        throw std::runtime_error("No readable files provided to read_sst_multi_netcdf");
    }

    int nlon = first.nlon;
    int nlat = first.nlat;
    int ndepth = first.ndepth;
    int npixels = nlon * nlat * ndepth;

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
    for (size_t f = first_file_index + 1; f < files.size(); ++f) {
        try {
            GridData gd = read_sst_netcdf(files[f], read_var_name, subset);

            // Verify spatial (and depth) grid matches
            if (gd.nlon != nlon || gd.nlat != nlat || gd.ndepth != ndepth) {
                std::string msg = "Grid mismatch in " + files[f] +
                    " (" + std::to_string(gd.nlon) + "x" + std::to_string(gd.nlat) +
                    "x" + std::to_string(gd.ndepth) +
                    " vs " + std::to_string(nlon) + "x" + std::to_string(nlat) +
                    "x" + std::to_string(ndepth) + ")";
                if (skip_bad_files) {
                    Rcpp::Rcerr << "Warning: " << msg << ", skipping" << std::endl;
                    continue;
                }
                throw std::runtime_error(msg);
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
            if (skip_bad_files) {
                Rcpp::Rcerr << "Warning: failed to read " << files[f]
                            << ": " << e.what() << ", skipping" << std::endl;
                continue;
            }
            throw;
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
    result.ndepth = ndepth;
    result.depth = first.depth;
    result.ntime = static_cast<int>(slices.size());
    result.time_units = first.time_units;
    result.time_calendar = first.time_calendar;
    result.temp_units = first.temp_units;

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
