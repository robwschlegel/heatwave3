#include "hw3_h5chunk.h"

#ifndef HW3_HAVE_HDF5

namespace hw3 {
bool h5_direct_available() { return false; }
std::vector<bool> h5_write_chunks(const std::string&, const std::vector<ChunkSrc>& vars, int) {
    return std::vector<bool>(vars.size(), false);
}
} // namespace hw3

#else

#include "hw3_omp.h"
#include <hdf5.h>
#include <zlib.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

namespace hw3 {

bool h5_direct_available() {
    // HEATWAVE3_SERIAL_WRITE=1 selects netCDF's own (single-threaded) writes,
    // e.g. to check that both paths give the same files
    const char* serial = std::getenv("HEATWAVE3_SERIAL_WRITE");
    if (serial && *serial && std::strcmp(serial, "0") != 0) return false;
    // Headers and library must be the same HDF5 release: struct layouts and
    // identifier sizes differ between releases.
    unsigned maj = 0, min = 0, rel = 0;
    if (H5get_libversion(&maj, &min, &rel) < 0) return false;
    return maj == H5_VERS_MAJOR && min == H5_VERS_MINOR && rel == H5_VERS_RELEASE;
}

namespace {

// One variable, as stored in the file
struct Plan {
    hid_t did = -1;
    const ChunkSrc* src = nullptr;
    int rank = 0;
    hsize_t dims[H5S_MAX_RANK], cdims[H5S_MAX_RANK];
    size_t esz = 0;          // bytes per element in the file
    bool to_float = false;   // double -> float32
    bool shuffle = false;
    int level = 0;           // deflate level
    size_t nchunks = 1;
    size_t chunk_elems = 1;
};

// Checks the variable can be written here (chunked, native byte order, a type
// matching the source, and only shuffle and/or deflate, shuffle first)
bool make_plan(hid_t fid, const ChunkSrc& v, Plan& p) {
    p.src = &v;
    p.did = H5Dopen2(fid, v.name.c_str(), H5P_DEFAULT);
    if (p.did < 0) return false;
    bool ok = false;
    hid_t dcpl = H5Dget_create_plist(p.did);
    hid_t type = H5Dget_type(p.did);
    hid_t space = H5Dget_space(p.did);
    do {
        if (dcpl < 0 || type < 0 || space < 0) break;
        if (H5Pget_layout(dcpl) != H5D_CHUNKED) break;
        p.rank = H5Sget_simple_extent_ndims(space);
        if (p.rank < 1 || H5Pget_chunk(dcpl, p.rank, p.cdims) != p.rank) break;
        H5Sget_simple_extent_dims(space, p.dims, nullptr);

        H5T_class_t cls = H5Tget_class(type);
        p.esz = H5Tget_size(type);
        if (H5Tget_order(type) != H5Tget_order(H5T_NATIVE_INT)) break;
        if (v.kind == ChunkSrc::DOUBLE) {
            if (cls != H5T_FLOAT || (p.esz != 4 && p.esz != 8)) break;
            p.to_float = (p.esz == 4);
        } else {
            size_t want = (v.kind == ChunkSrc::SCHAR) ? 1 : 4;
            if (cls != H5T_INTEGER || p.esz != want || H5Tget_sign(type) != H5T_SGN_2) break;
        }

        int nf = H5Pget_nfilters(dcpl);
        bool good = nf >= 1;
        for (int i = 0; i < nf && good; ++i) {
            unsigned flags = 0, cd[8] = {0}, config = 0;
            size_t ncd = 8;
            H5Z_filter_t f = H5Pget_filter2(dcpl, i, &flags, &ncd, cd, 0, nullptr, &config);
            if (f == H5Z_FILTER_SHUFFLE && i == 0) p.shuffle = true;
            else if (f == H5Z_FILTER_DEFLATE && i == nf - 1 && ncd >= 1) p.level = static_cast<int>(cd[0]);
            else good = false;
        }
        if (!good || p.level < 1) break;

        for (int k = 0; k < p.rank; ++k) {
            if (p.cdims[k] == 0) { good = false; break; }
            p.nchunks *= (p.dims[k] + p.cdims[k] - 1) / p.cdims[k];
            p.chunk_elems *= p.cdims[k];
        }
        ok = good;
    } while (false);
    if (space >= 0) H5Sclose(space);
    if (type >= 0) H5Tclose(type);
    if (dcpl >= 0) H5Pclose(dcpl);
    if (!ok) { H5Dclose(p.did); p.did = -1; }
    return ok;
}

// Element offset of chunk `ci` (row-major over the chunk grid)
void chunk_origin(const Plan& p, size_t ci, hsize_t* off) {
    for (int k = p.rank - 1; k >= 0; --k) {
        size_t n = (p.dims[k] + p.cdims[k] - 1) / p.cdims[k];
        off[k] = (ci % n) * p.cdims[k];
        ci /= n;
    }
}

// Copies one chunk out of the source array, converted to the file type.
// Edge chunks are zero padded (HDF5 always stores full chunks).
void fill_chunk(const Plan& p, const hsize_t* off, unsigned char* out) {
    std::memset(out, 0, p.chunk_elems * p.esz);
    int r = p.rank;
    hsize_t ext[H5S_MAX_RANK], stride[H5S_MAX_RANK], cstride[H5S_MAX_RANK], idx[H5S_MAX_RANK] = {0};
    for (int k = 0; k < r; ++k) ext[k] = std::min(p.cdims[k], p.dims[k] - off[k]);
    stride[r - 1] = 1; cstride[r - 1] = 1;
    for (int k = r - 2; k >= 0; --k) {
        stride[k] = stride[k + 1] * p.dims[k + 1];
        cstride[k] = cstride[k + 1] * p.cdims[k + 1];
    }
    size_t run = ext[r - 1];
    for (;;) {
        size_t s = off[r - 1], d = 0;
        for (int k = 0; k < r - 1; ++k) { s += (off[k] + idx[k]) * stride[k]; d += idx[k] * cstride[k]; }
        switch (p.src->kind) {
        case ChunkSrc::DOUBLE: {
            const double* src = static_cast<const double*>(p.src->data) + s;
            if (p.to_float) {
                float* dst = reinterpret_cast<float*>(out) + d;
                for (size_t i = 0; i < run; ++i) dst[i] = static_cast<float>(src[i]);
            } else {
                std::memcpy(reinterpret_cast<double*>(out) + d, src, run * sizeof(double));
            }
            break;
        }
        case ChunkSrc::SCHAR:
            std::memcpy(out + d, static_cast<const signed char*>(p.src->data) + s, run);
            break;
        case ChunkSrc::INT:
            std::memcpy(reinterpret_cast<int*>(out) + d, static_cast<const int*>(p.src->data) + s,
                        run * sizeof(int));
            break;
        }
        int k = r - 2;
        while (k >= 0 && ++idx[k] == ext[k]) idx[k--] = 0;
        if (k < 0) break;
    }
}

// HDF5's shuffle filter: byte j of every element, for j = 0 .. esz - 1
void shuffle_bytes(const unsigned char* in, unsigned char* out, size_t n, size_t esz) {
    for (size_t j = 0; j < esz; ++j) {
        unsigned char* o = out + j * n;
        for (size_t i = 0; i < n; ++i) o[i] = in[i * esz + j];
    }
}

} // namespace

std::vector<bool> h5_write_chunks(const std::string& file,
                                  const std::vector<ChunkSrc>& vars,
                                  int n_threads) {
    std::vector<bool> written(vars.size(), false);
    if (vars.empty() || !h5_direct_available()) return written;

    // Silence HDF5's error stack printing while probing; restored below
    H5E_auto2_t old_func = nullptr;
    void* old_data = nullptr;
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_data);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

    hid_t fid = H5Fopen(file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (fid < 0) {
        H5Eset_auto2(H5E_DEFAULT, old_func, old_data);
        return written;
    }

    std::vector<Plan> plans;
    std::vector<size_t> which;  // index into vars
    for (size_t v = 0; v < vars.size(); ++v) {
        Plan p;
        if (make_plan(fid, vars[v], p)) { plans.push_back(p); which.push_back(v); }
    }

    // Tasks are (variable, chunk) pairs, compressed in batches so that only a
    // bounded number of compressed chunks is held before being stored
    std::vector<std::pair<int, size_t>> tasks;
    size_t max_bytes = 0;
    for (size_t i = 0; i < plans.size(); ++i) {
        for (size_t c = 0; c < plans[i].nchunks; ++c) tasks.emplace_back(static_cast<int>(i), c);
        max_bytes = std::max(max_bytes, plans[i].chunk_elems * plans[i].esz);
    }
    int nt = std::max(1, n_threads > 0 ? n_threads : get_threads());
    size_t batch = static_cast<size_t>(nt) * 4;
    std::vector<std::vector<unsigned char>> raw(nt, std::vector<unsigned char>(max_bytes));
    std::vector<std::vector<unsigned char>> shuf(nt, std::vector<unsigned char>(max_bytes));
    std::vector<std::vector<unsigned char>> comp(batch);
    std::vector<int> zerr(batch, Z_OK);

    std::string error;
    for (size_t b0 = 0; b0 < tasks.size() && error.empty(); b0 += batch) {
        size_t nb = std::min(batch, tasks.size() - b0);
        parallel_for(static_cast<int>(nb), nt, [&](int wid, int j) {
            const Plan& p = plans[tasks[b0 + j].first];
            hsize_t off[H5S_MAX_RANK];
            chunk_origin(p, tasks[b0 + j].second, off);
            size_t nbytes = p.chunk_elems * p.esz;
            fill_chunk(p, off, raw[wid].data());
            const unsigned char* in = raw[wid].data();
            if (p.shuffle && p.esz > 1) {
                shuffle_bytes(in, shuf[wid].data(), p.chunk_elems, p.esz);
                in = shuf[wid].data();
            }
            uLongf len = compressBound(static_cast<uLong>(nbytes));
            comp[j].resize(len);
            zerr[j] = compress2(comp[j].data(), &len, in, static_cast<uLong>(nbytes), p.level);
            comp[j].resize(len);
        }, [](int) {});

        for (size_t j = 0; j < nb && error.empty(); ++j) {
            const Plan& p = plans[tasks[b0 + j].first];
            hsize_t off[H5S_MAX_RANK];
            chunk_origin(p, tasks[b0 + j].second, off);
            if (zerr[j] != Z_OK) {
                error = "zlib compression failed for " + p.src->name;
            } else if (H5Dwrite_chunk(p.did, H5P_DEFAULT, 0, off, comp[j].size(),
                                      comp[j].data()) < 0) {
                error = "HDF5 chunk write failed for " + p.src->name;
            }
            std::vector<unsigned char>().swap(comp[j]);
        }
    }

    for (auto& p : plans) H5Dclose(p.did);
    herr_t closed = H5Fclose(fid);
    H5Eset_auto2(H5E_DEFAULT, old_func, old_data);
    if (!error.empty()) throw std::runtime_error(error + " in " + file);
    if (closed < 0) throw std::runtime_error("HDF5 could not close " + file);
    for (size_t v : which) written[v] = true;
    return written;
}

} // namespace hw3

#endif
