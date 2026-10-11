// Python bindings (pybind11). Converts between Python/numpy types and the
// language-neutral API in src/hw3_api.h, mirroring the R exports in
// src/heatwave3_init.cpp name for name. Missing values: NaN for floats,
// INT_MIN (hw3::NA_INT) for integers.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "hw3_api.h"
#include "hw3_log.h"
#include "hw3_omp.h"
#include "blob_core.h"
#include "climatology.h"
#include "netcdf_io.h"

namespace py = pybind11;

#ifndef HW3_VERSION
#define HW3_VERSION "0.0.0"
#endif

// ---- Message sinks for the shared core (hw3_log.h) ----
//
// Text is buffered and forwarded to sys.stdout / sys.stderr on each flush
// (std::endl, std::flush), reacquiring the GIL, so messages appear in Jupyter
// as well as terminals. The core only writes from the calling thread.

namespace {
class PyStreamBuf : public std::stringbuf {
public:
    explicit PyStreamBuf(const char* name) : name_(name) {}
protected:
    int sync() override {
        std::string s = str();
        str("");
        if (s.empty()) return 0;
        try {
            py::gil_scoped_acquire gil;
            py::object stream = py::module_::import("sys").attr(name_);
            if (!stream.is_none()) {
                stream.attr("write")(s);
                stream.attr("flush")();
            }
        } catch (...) {
            // Never let a logging failure abort a computation.
        }
        return 0;
    }
private:
    const char* name_;
};
} // namespace

namespace hw3 {
std::ostream& out() {
    static PyStreamBuf buf("stdout");
    static std::ostream os(&buf);
    return os;
}
std::ostream& err() {
    static PyStreamBuf buf("stderr");
    static std::ostream os(&buf);
    return os;
}
} // namespace hw3

// ---- Conversion helpers ----

namespace {

using OptDbl = std::optional<std::vector<double>>;
using OptStr = std::optional<std::vector<std::string>>;

// Hand a std::vector to numpy without copying; the array owns the vector.
template <class T>
py::array_t<T> np(std::vector<T>&& v) {
    auto* p = new std::vector<T>(std::move(v));
    py::capsule owner(p, [](void* q) { delete static_cast<std::vector<T>*>(q); });
    return py::array_t<T>(static_cast<py::ssize_t>(p->size()), p->data(), owner);
}

template <class T>
std::vector<T> val(const std::optional<std::vector<T>>& x) {
    return x ? *x : std::vector<T>{};
}

py::dict grid_dict(hw3::GridData&& gd) {
    py::dict d;
    d["lon"] = np(std::move(gd.lon));
    d["lat"] = np(std::move(gd.lat));
    d["nlon"] = gd.nlon;
    d["nlat"] = gd.nlat;
    d["ntime"] = gd.ntime;
    d["time_days"] = np(std::move(gd.time_days));
    d["ndepth"] = gd.ndepth;
    d["depth"] = np(std::move(gd.depth));
    d["sst"] = np(std::move(gd.sst));
    return d;
}

hw3::ClimParams clim_params(int maxPadLength, int windowHalfWidth, double pctile,
                            bool smoothPercentile, int smoothPercentileWidth,
                            bool compute_var, int roundClm, int n_threads,
                            bool detrend) {
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

hw3::DetectParams detect_params(int minDuration, int minDuration2, bool joinAcrossGaps,
                                int maxGap, int maxGap2, bool coldSpells, int roundRes,
                                int n_threads, bool category, bool southHemisphere) {
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

py::dict daily_cat_dict(hw3::DailyCatResult&& r) {
    py::dict d;
    d["jd"] = np(std::move(r.jd));
    d["lon"] = np(std::move(r.lon));
    d["lat"] = np(std::move(r.lat));
    d["temp"] = np(std::move(r.temp));
    d["seas"] = np(std::move(r.seas));
    d["thresh"] = np(std::move(r.thresh));
    d["intensity"] = np(std::move(r.intensity));
    d["event"] = np(std::move(r.event));
    d["event_no"] = np(std::move(r.event_no));
    d["category"] = np(std::move(r.category));
    if (r.has_depth) d["depth"] = np(std::move(r.depth));
    return d;
}

} // namespace

// ---- Module ----

PYBIND11_MODULE(_core, m) {
    m.doc() = "heatwave3 C++ core (internal; use the heatwave3 package API).";
    m.attr("NA_INT") = hw3::NA_INT;

    // Thread management and metadata

    m.def("hw3_get_threads", []() { return hw3::get_threads(); });
    m.def("hw3_set_threads", [](int n) { hw3::set_threads(n); return hw3::get_threads(); },
          py::arg("n"));
    m.def("hw3_version", []() {
        py::dict d;
        d["package"] = "heatwave3";
        d["version"] = HW3_VERSION;
        d["cpp_standard"] = "C++17";
        return d;
    });
    m.def("hw3_jd_to_doy", [](int jd) { return hw3::jd_to_doy_366(jd); }, py::arg("jd"));

    // SST reading

    m.def("hw3_read_sst",
          [](const std::string& file_in, const std::string& var_name,
             const OptDbl& lon_range, const OptDbl& lat_range, const OptStr& time_range,
             int depth, const OptDbl& depth_range) {
              hw3::GridData gd;
              {
                  py::gil_scoped_release nogil;
                  hw3::SubsetSpec ss = hw3::make_subset(val(lon_range), val(lat_range),
                                                        val(time_range), depth, val(depth_range));
                  gd = hw3::read_sst(file_in, var_name, ss);
              }
              return grid_dict(std::move(gd));
          },
          py::arg("file_in"), py::arg("var_name") = "", py::arg("lon_range") = py::none(),
          py::arg("lat_range") = py::none(), py::arg("time_range") = py::none(),
          py::arg("depth") = -1, py::arg("depth_range") = py::none());

    m.def("hw3_read_sst_multi",
          [](const std::vector<std::string>& files, const std::string& var_name,
             const OptDbl& lon_range, const OptDbl& lat_range, int depth,
             bool skip_bad_files, const OptDbl& depth_range) {
              hw3::GridData gd;
              {
                  py::gil_scoped_release nogil;
                  hw3::SubsetSpec ss = hw3::make_subset(val(lon_range), val(lat_range), {},
                                                        depth, val(depth_range));
                  gd = hw3::read_sst_multi(files, var_name, ss, skip_bad_files);
              }
              return grid_dict(std::move(gd));
          },
          py::arg("files"), py::arg("var_name") = "", py::arg("lon_range") = py::none(),
          py::arg("lat_range") = py::none(), py::arg("depth") = -1,
          py::arg("skip_bad_files") = false, py::arg("depth_range") = py::none());

    // Climatology

    m.def("hw3_pixel_clim",
          [](py::array_t<double, py::array::c_style | py::array::forcecast> temp,
             py::array_t<int, py::array::c_style | py::array::forcecast> time_jd,
             const std::string& clim_start, const std::string& clim_end,
             int windowHalfWidth, double pctile, bool smoothPercentile,
             int smoothPercentileWidth, int maxPadLength, bool detrend) {
              if (temp.size() != time_jd.size())
                  throw std::invalid_argument("temp and time_jd must have the same length.");
              std::vector<double> t(temp.data(), temp.data() + temp.size());
              std::vector<int> jd(time_jd.data(), time_jd.data() + time_jd.size());
              hw3::ClimParams p = clim_params(maxPadLength, windowHalfWidth, pctile,
                                              smoothPercentile, smoothPercentileWidth,
                                              false, 0, 1, detrend);
              hw3::PixelClimResult cr;
              {
                  py::gil_scoped_release nogil;
                  cr = hw3::pixel_clim(t, jd, clim_start, clim_end, p);
              }
              py::dict d;
              d["seas"] = np(std::vector<double>(cr.seas, cr.seas + 366));
              d["thresh"] = np(std::vector<double>(cr.thresh, cr.thresh + 366));
              d["valid"] = cr.valid;
              return d;
          },
          py::arg("temp"), py::arg("time_jd"), py::arg("clim_start"), py::arg("clim_end"),
          py::arg("windowHalfWidth") = 5, py::arg("pctile") = 90.0,
          py::arg("smoothPercentile") = true, py::arg("smoothPercentileWidth") = 31,
          py::arg("maxPadLength") = 0, py::arg("detrend") = false);

    m.def("hw3_compute_clim",
          [](const std::string& file_in, const std::string& file_out,
             const std::string& var_name, const std::vector<std::string>& climatologyPeriod,
             const OptDbl& lon_range, const OptDbl& lat_range, const OptStr& time_range,
             int depth, int maxPadLength, int windowHalfWidth, double pctile,
             bool smoothPercentile, int smoothPercentileWidth, bool compute_var,
             int roundClm, int n_threads, bool detrend, const OptDbl& depth_range) {
              if (climatologyPeriod.size() != 2)
                  throw std::invalid_argument("climatologyPeriod must have two dates.");
              py::gil_scoped_release nogil;
              hw3::SubsetSpec ss = hw3::make_subset(val(lon_range), val(lat_range),
                                                    val(time_range), depth, val(depth_range));
              hw3::compute_clim(file_in, file_out, var_name,
                                climatologyPeriod[0], climatologyPeriod[1], ss,
                                clim_params(maxPadLength, windowHalfWidth, pctile,
                                            smoothPercentile, smoothPercentileWidth,
                                            compute_var, roundClm, n_threads, detrend));
          },
          py::arg("file_in"), py::arg("file_out"), py::arg("var_name"),
          py::arg("climatologyPeriod"), py::arg("lon_range") = py::none(),
          py::arg("lat_range") = py::none(), py::arg("time_range") = py::none(),
          py::arg("depth") = -1, py::arg("maxPadLength") = 0, py::arg("windowHalfWidth") = 5,
          py::arg("pctile") = 90.0, py::arg("smoothPercentile") = true,
          py::arg("smoothPercentileWidth") = 31, py::arg("compute_var") = false,
          py::arg("roundClm") = 4, py::arg("n_threads") = 1, py::arg("detrend") = false,
          py::arg("depth_range") = py::none());

    m.def("hw3_compute_clim_multi",
          [](const std::vector<std::string>& files, const std::string& file_out,
             const std::string& var_name, const std::vector<std::string>& climatologyPeriod,
             const OptDbl& lon_range, const OptDbl& lat_range, int depth, int maxPadLength,
             int windowHalfWidth, double pctile, bool smoothPercentile,
             int smoothPercentileWidth, bool compute_var, int roundClm, int n_threads,
             bool detrend, bool skip_bad_files, const OptDbl& depth_range) {
              if (climatologyPeriod.size() != 2)
                  throw std::invalid_argument("climatologyPeriod must have two dates.");
              py::gil_scoped_release nogil;
              hw3::SubsetSpec ss = hw3::make_subset(val(lon_range), val(lat_range), {},
                                                    depth, val(depth_range));
              hw3::compute_clim_multi(files, file_out, var_name,
                                      climatologyPeriod[0], climatologyPeriod[1], ss,
                                      clim_params(maxPadLength, windowHalfWidth, pctile,
                                                  smoothPercentile, smoothPercentileWidth,
                                                  compute_var, roundClm, n_threads, detrend),
                                      skip_bad_files);
          },
          py::arg("files"), py::arg("file_out"), py::arg("var_name"),
          py::arg("climatologyPeriod"), py::arg("lon_range") = py::none(),
          py::arg("lat_range") = py::none(), py::arg("depth") = -1,
          py::arg("maxPadLength") = 0, py::arg("windowHalfWidth") = 5,
          py::arg("pctile") = 90.0, py::arg("smoothPercentile") = true,
          py::arg("smoothPercentileWidth") = 31, py::arg("compute_var") = false,
          py::arg("roundClm") = 4, py::arg("n_threads") = 1, py::arg("detrend") = false,
          py::arg("skip_bad_files") = false, py::arg("depth_range") = py::none());

    m.def("hw3_write_const_clim",
          [](const std::string& file_in, const std::string& clim_file,
             const std::string& var_name, double threshold, const OptDbl& lon_range,
             const OptDbl& lat_range, const OptStr& time_range, int depth) {
              py::gil_scoped_release nogil;
              hw3::SubsetSpec ss = hw3::make_subset(val(lon_range), val(lat_range),
                                                    val(time_range), depth, {});
              hw3::write_const_clim(file_in, clim_file, var_name, threshold, ss);
          },
          py::arg("file_in"), py::arg("clim_file"), py::arg("var_name"), py::arg("threshold"),
          py::arg("lon_range") = py::none(), py::arg("lat_range") = py::none(),
          py::arg("time_range") = py::none(), py::arg("depth") = -1);

    // Event detection

    m.def("hw3_detect_events",
          [](const std::string& file_in, const std::string& clim_file,
             const std::string& events_file, const std::string& var_name,
             int minDuration, int minDuration2, bool joinAcrossGaps, int maxGap, int maxGap2,
             bool coldSpells, int roundRes, int n_threads, bool category,
             bool southHemisphere, const std::string& threshClim2_file,
             const std::string& threshClim2_var_name, const std::string& daily_file,
             const std::string& proto_file) {
              py::gil_scoped_release nogil;
              hw3::detect_events(file_in, clim_file, events_file, var_name,
                                 detect_params(minDuration, minDuration2, joinAcrossGaps,
                                               maxGap, maxGap2, coldSpells, roundRes,
                                               n_threads, category, southHemisphere),
                                 threshClim2_file, threshClim2_var_name, daily_file,
                                 proto_file);
          },
          py::arg("file_in"), py::arg("clim_file"), py::arg("events_file"),
          py::arg("var_name"), py::arg("minDuration") = 5, py::arg("minDuration2") = 5,
          py::arg("joinAcrossGaps") = true, py::arg("maxGap") = 2, py::arg("maxGap2") = 2,
          py::arg("coldSpells") = false, py::arg("roundRes") = 4, py::arg("n_threads") = 1,
          py::arg("category") = false, py::arg("southHemisphere") = true,
          py::arg("threshClim2_file") = "", py::arg("threshClim2_var_name") = "",
          py::arg("daily_file") = "", py::arg("proto_file") = "");

    m.def("hw3_detect_events_multi",
          [](const std::vector<std::string>& files, const std::string& clim_file,
             const std::string& events_file, const std::string& var_name,
             int minDuration, int minDuration2, bool joinAcrossGaps, int maxGap, int maxGap2,
             bool coldSpells, int roundRes, int n_threads, bool category,
             bool southHemisphere, const OptStr& threshClim2_files,
             const std::string& threshClim2_var_name, bool skip_bad_files,
             const std::string& daily_file, const std::string& proto_file) {
              py::gil_scoped_release nogil;
              hw3::detect_events_multi(files, clim_file, events_file, var_name,
                                       detect_params(minDuration, minDuration2,
                                                     joinAcrossGaps, maxGap, maxGap2,
                                                     coldSpells, roundRes, n_threads,
                                                     category, southHemisphere),
                                       val(threshClim2_files), threshClim2_var_name,
                                       skip_bad_files, daily_file, proto_file);
          },
          py::arg("files"), py::arg("clim_file"), py::arg("events_file"),
          py::arg("var_name"), py::arg("minDuration") = 5, py::arg("minDuration2") = 5,
          py::arg("joinAcrossGaps") = true, py::arg("maxGap") = 2, py::arg("maxGap2") = 2,
          py::arg("coldSpells") = false, py::arg("roundRes") = 4, py::arg("n_threads") = 1,
          py::arg("category") = false, py::arg("southHemisphere") = true,
          py::arg("threshClim2_files") = py::none(), py::arg("threshClim2_var_name") = "",
          py::arg("skip_bad_files") = false, py::arg("daily_file") = "",
          py::arg("proto_file") = "");

    // Product readers

    m.def("hw3_file_meta", [](const std::string& file) {
        hw3::FileMeta fm;
        {
            py::gil_scoped_release nogil;
            fm = hw3::read_file_meta(file);
        }
        py::dict d;
        d["product"] = fm.product;
        d["nlon"] = fm.nlon;
        d["nlat"] = fm.nlat;
        d["n3"] = fm.n3;
        d["nrows"] = fm.nrows;
        return d;
    }, py::arg("file"));

    m.def("hw3_read_clim_nc", [](const std::string& clim_file) {
        hw3::ClimData cd;
        {
            py::gil_scoped_release nogil;
            cd = hw3::read_clim_netcdf(clim_file);
        }
        py::dict d;
        d["lon"] = np(std::move(cd.lon));
        d["lat"] = np(std::move(cd.lat));
        d["nlon"] = cd.nlon;
        d["nlat"] = cd.nlat;
        d["ndoy"] = cd.ndoy;
        d["ndepth"] = cd.ndepth;
        d["depth"] = np(std::move(cd.depth));
        d["seas"] = np(std::move(cd.seas));
        d["thresh"] = np(std::move(cd.thresh));
        return d;
    }, py::arg("clim_file"));

    m.def("hw3_read_event_nc", [](const std::string& event_file) {
        hw3::EventData ed;
        {
            py::gil_scoped_release nogil;
            ed = hw3::read_event_netcdf(event_file);
        }
        py::dict d;
        d["lon"] = np(std::move(ed.lon));
        d["lat"] = np(std::move(ed.lat));
        d["depth"] = np(std::move(ed.depth));
        d["pixel_index"] = np(std::move(ed.pixel_index));
        d["event_no"] = np(std::move(ed.event_no));
        d["date_start"] = np(std::move(ed.date_start));
        d["date_peak"] = np(std::move(ed.date_peak));
        d["date_end"] = np(std::move(ed.date_end));
        d["duration"] = np(std::move(ed.duration));
        d["intensity_mean"] = np(std::move(ed.intensity_mean));
        d["intensity_max"] = np(std::move(ed.intensity_max));
        d["intensity_var"] = np(std::move(ed.intensity_var));
        d["intensity_cumulative"] = np(std::move(ed.intensity_cumulative));
        d["intensity_mean_relThresh"] = np(std::move(ed.intensity_mean_relThresh));
        d["intensity_max_relThresh"] = np(std::move(ed.intensity_max_relThresh));
        d["intensity_var_relThresh"] = np(std::move(ed.intensity_var_relThresh));
        d["intensity_cumulative_relThresh"] = np(std::move(ed.intensity_cumulative_relThresh));
        d["intensity_mean_abs"] = np(std::move(ed.intensity_mean_abs));
        d["intensity_max_abs"] = np(std::move(ed.intensity_max_abs));
        d["intensity_var_abs"] = np(std::move(ed.intensity_var_abs));
        d["intensity_cumulative_abs"] = np(std::move(ed.intensity_cumulative_abs));
        d["rate_onset"] = np(std::move(ed.rate_onset));
        d["rate_decline"] = np(std::move(ed.rate_decline));
        d["category"] = np(std::move(ed.category));
        d["p_moderate"] = np(std::move(ed.p_moderate));
        d["p_strong"] = np(std::move(ed.p_strong));
        d["p_severe"] = np(std::move(ed.p_severe));
        d["p_extreme"] = np(std::move(ed.p_extreme));
        d["season"] = np(std::move(ed.season));
        d["ref_date_jd"] = ed.ref_date_jd;
        d["nevents"] = ed.nevents;
        return d;
    }, py::arg("event_file"));

    m.def("hw3_read_daily_nc", [](const std::string& daily_file) {
        hw3::DailyData dd;
        {
            py::gil_scoped_release nogil;
            dd = hw3::read_daily_netcdf(daily_file);
        }
        py::dict d;
        d["lon"] = np(std::move(dd.lon));
        d["lat"] = np(std::move(dd.lat));
        d["time_jd"] = np(std::move(dd.time_jd));
        d["nlon"] = dd.nlon;
        d["nlat"] = dd.nlat;
        d["ntime"] = dd.ntime;
        d["ndepth"] = dd.ndepth;
        d["depth"] = np(std::move(dd.depth));
        d["temp"] = np(std::move(dd.temp));
        d["seas"] = np(std::move(dd.seas));
        d["thresh"] = np(std::move(dd.thresh));
        d["intensity"] = np(std::move(dd.intensity));
        d["threshCriterion"] = np(std::move(dd.threshCriterion));
        d["durationCriterion"] = np(std::move(dd.durationCriterion));
        d["event"] = np(std::move(dd.event));
        d["category"] = np(std::move(dd.category));
        d["event_no"] = np(std::move(dd.event_no));
        return d;
    }, py::arg("daily_file"));

    m.def("hw3_read_subset",
          [](const std::string& file, const OptDbl& lon_range, const OptDbl& lat_range,
             const std::optional<std::vector<int>>& t_jd_range, const OptDbl& depth_range,
             const OptStr& vars, int max_rows) {
              hw3::SubsetResult s;
              {
                  py::gil_scoped_release nogil;
                  s = hw3::read_subset_netcdf(file, val(lon_range), val(lat_range),
                                              val(t_jd_range), val(depth_range), val(vars),
                                              max_rows);
              }
              py::dict data;
              for (auto& v : s.data) {
                  if (v.is_double) data[py::str(v.name)] = np(std::move(v.dbl));
                  else             data[py::str(v.name)] = np(std::move(v.intg));
              }
              py::dict d;
              d["product"] = s.product;
              d["lon"] = np(std::move(s.lon));
              d["lat"] = np(std::move(s.lat));
              d["has_depth"] = s.has_depth;
              d["depth"] = np(std::move(s.depth));
              d["third"] = np(std::move(s.third));
              d["third_name"] = s.third_name;
              d["nlon"] = s.nlon;
              d["nlat"] = s.nlat;
              d["ndepth"] = s.ndepth;
              d["n3"] = s.n3;
              d["data"] = data;
              return d;
          },
          py::arg("file"), py::arg("lon_range") = py::none(), py::arg("lat_range") = py::none(),
          py::arg("t_jd_range") = py::none(), py::arg("depth_range") = py::none(),
          py::arg("vars") = py::none(), py::arg("max_rows") = -1);

    // Daily categories for a date window

    m.def("hw3_category_daily",
          [](const std::string& sst_file, const std::string& clim_file,
             const std::string& event_file, const std::vector<int>& t_jd_range,
             const std::string& var_name, bool coldSpells, double ice_thresh, int roundRes) {
              hw3::DailyCatResult r;
              {
                  py::gil_scoped_release nogil;
                  r = hw3::category_daily(sst_file, clim_file, event_file, t_jd_range,
                                          var_name, coldSpells, ice_thresh, roundRes);
              }
              return daily_cat_dict(std::move(r));
          },
          py::arg("sst_file"), py::arg("clim_file"), py::arg("event_file"),
          py::arg("t_jd_range"), py::arg("var_name") = "", py::arg("coldSpells") = false,
          py::arg("ice_thresh") = -1.7, py::arg("roundRes") = 2);

    m.def("hw3_category_daily_multi",
          [](const std::vector<std::string>& files, const std::string& clim_file,
             const std::string& event_file, const std::vector<int>& t_jd_range,
             const std::string& var_name, bool coldSpells, double ice_thresh, int roundRes,
             bool skip_bad_files) {
              hw3::DailyCatResult r;
              {
                  py::gil_scoped_release nogil;
                  r = hw3::category_daily_multi(files, clim_file, event_file, t_jd_range,
                                                var_name, coldSpells, ice_thresh, roundRes,
                                                skip_bad_files);
              }
              return daily_cat_dict(std::move(r));
          },
          py::arg("files"), py::arg("clim_file"), py::arg("event_file"),
          py::arg("t_jd_range"), py::arg("var_name") = "", py::arg("coldSpells") = false,
          py::arg("ice_thresh") = -1.7, py::arg("roundRes") = 2,
          py::arg("skip_bad_files") = false);

    // Event-file summaries

    m.def("hw3_category",
          [](const std::string& event_file, const std::string& clim_file,
             const std::string& hemisphere, int roundVal) {
              hw3::CategoryResult r;
              {
                  py::gil_scoped_release nogil;
                  r = hw3::category(event_file, clim_file, hemisphere, roundVal);
              }
              py::dict d;
              d["event_no"] = np(std::move(r.event_no));
              d["lon"] = np(std::move(r.lon));
              d["lat"] = np(std::move(r.lat));
              d["peak_jd"] = np(std::move(r.peak_jd));
              d["category"] = np(std::move(r.category));   // 0 = NA, 1..4
              d["intensity_max"] = np(std::move(r.intensity_max));
              d["duration"] = np(std::move(r.duration));
              d["p_moderate"] = np(std::move(r.p_moderate));
              d["p_strong"] = np(std::move(r.p_strong));
              d["p_severe"] = np(std::move(r.p_severe));
              d["p_extreme"] = np(std::move(r.p_extreme));
              d["season"] = np(std::move(r.season));       // 0 = unset, 1..4
              return d;
          },
          py::arg("event_file"), py::arg("clim_file") = "", py::arg("hemisphere") = "south",
          py::arg("roundVal") = 4);

    m.def("hw3_block_average", [](const std::string& event_file) {
        hw3::BlockAverageResult r;
        {
            py::gil_scoped_release nogil;
            r = hw3::block_average(event_file);
        }
        py::dict d;
        d["lon"] = np(std::move(r.lon));
        d["lat"] = np(std::move(r.lat));
        d["year"] = np(std::move(r.year));
        d["count"] = np(std::move(r.count));
        d["duration_mean"] = np(std::move(r.duration_mean));
        d["duration_max"] = np(std::move(r.duration_max));
        d["intensity_mean"] = np(std::move(r.intensity_mean));
        d["intensity_max_mean"] = np(std::move(r.intensity_max_mean));
        d["intensity_max_max"] = np(std::move(r.intensity_max_max));
        d["intensity_cumulative_mean"] = np(std::move(r.intensity_cumulative_mean));
        d["total_days"] = np(std::move(r.total_days));
        d["total_icum"] = np(std::move(r.total_icum));
        return d;
    }, py::arg("event_file"));

    m.def("hw3_read_metric_summary",
          [](const std::string& event_file, const std::string& metric,
             const std::string& summary) {
              hw3::MetricSummaryResult r;
              {
                  py::gil_scoped_release nogil;
                  r = hw3::metric_summary(event_file, metric, summary);
              }
              py::dict d;
              d["lon"] = np(std::move(r.lon));
              d["lat"] = np(std::move(r.lat));
              d["value"] = np(std::move(r.value));
              return d;
          },
          py::arg("event_file"), py::arg("metric") = "intensity_max",
          py::arg("summary") = "mean");

    // Blob detection

    m.def("label_components",
          [](py::array_t<int, py::array::f_style | py::array::forcecast> mask,
             int nx, int ny, int nz, int nt, bool wrap_dateline) {
              const py::ssize_t n = static_cast<py::ssize_t>(nx) * ny * nz * nt;
              if (mask.size() != n)
                  throw std::invalid_argument("mask length does not match nx * ny * nz * nt");
              py::array_t<int, py::array::f_style> labels(
                  std::vector<py::ssize_t>{nx, ny, nz, nt});
              std::fill(labels.mutable_data(), labels.mutable_data() + n, 0);
              int k;
              {
                  py::gil_scoped_release nogil;
                  k = hw3::label_components(mask.data(), labels.mutable_data(),
                                            nx, ny, nz, nt, wrap_dateline);
              }
              return py::make_tuple(labels, k);
          },
          "Column-major [i, j, kd, kt] labelling; returns (labels, n_components).",
          py::arg("mask"), py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("nt"),
          py::arg("wrap_dateline") = false);

    m.def("blob_daily_summary",
          [](py::array_t<int, py::array::c_style | py::array::forcecast> labels,
             py::array_t<int, py::array::c_style | py::array::forcecast> ix,
             py::array_t<int, py::array::c_style | py::array::forcecast> iy,
             py::array_t<int, py::array::c_style | py::array::forcecast> it,
             py::array_t<double, py::array::c_style | py::array::forcecast> delta,
             py::array_t<double, py::array::c_style | py::array::forcecast> cell_area_lat,
             py::array_t<double, py::array::c_style | py::array::forcecast> lons,
             py::array_t<double, py::array::c_style | py::array::forcecast> lats,
             py::array_t<int, py::array::c_style | py::array::forcecast> iz,
             py::array_t<double, py::array::c_style | py::array::forcecast> depths,
             py::array_t<double, py::array::c_style | py::array::forcecast> layer_thickness_km,
             bool has_depth) {
              const py::ssize_t n = labels.size();
              if (ix.size() != n || iy.size() != n || it.size() != n || delta.size() != n)
                  throw std::invalid_argument("Input vectors must all have the same length.");
              if (has_depth && iz.size() != n)
                  throw std::invalid_argument(
                      "iz must have the same length as labels when has_depth is True.");
              hw3::BlobDailySummary s;
              {
                  py::gil_scoped_release nogil;
                  s = hw3::blob_daily_summary(labels.data(), ix.data(), iy.data(), it.data(),
                                              delta.data(), n, cell_area_lat.data(),
                                              lons.data(), lats.data(), iz.data(),
                                              depths.data(), layer_thickness_km.data(),
                                              has_depth);
              }
              py::dict d;
              d["blob"] = np(std::move(s.blob));
              d["t_idx"] = np(std::move(s.t_idx));
              d["area_km2"] = np(std::move(s.area_km2));
              d["mean_delta"] = np(std::move(s.mean_delta));
              d["max_delta"] = np(std::move(s.max_delta));
              d["centroid_lon"] = np(std::move(s.centroid_lon));
              d["centroid_lat"] = np(std::move(s.centroid_lat));
              d["lon_min"] = np(std::move(s.lon_min));
              d["lon_max"] = np(std::move(s.lon_max));
              d["lat_min"] = np(std::move(s.lat_min));
              d["lat_max"] = np(std::move(s.lat_max));
              d["n_cells"] = np(std::move(s.n_cells));
              if (has_depth) {
                  d["volume_km3"] = np(std::move(s.volume_km3));
                  d["depth_min_m"] = np(std::move(s.depth_min_m));
                  d["depth_max_m"] = np(std::move(s.depth_max_m));
                  d["depth_at_max_delta_m"] = np(std::move(s.depth_at_max_delta_m));
              }
              return d;
          },
          "1-based voxel indices, as in the R blob_daily_summary_cpp().",
          py::arg("labels"), py::arg("ix"), py::arg("iy"), py::arg("it"), py::arg("delta"),
          py::arg("cell_area_lat"), py::arg("lons"), py::arg("lats"), py::arg("iz"),
          py::arg("depths"), py::arg("layer_thickness_km"), py::arg("has_depth"));
}
