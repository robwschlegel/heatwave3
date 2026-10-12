#ifndef HW3_H5CHUNK_H
#define HW3_H5CHUNK_H

// Parallel compression of gridded NetCDF-4 variables.
//
// HDF5 (which stores NetCDF-4 files) compresses chunks one at a time on the
// calling thread, so writing compressed per-day products is limited to one
// core. Here the variables are defined and the file closed by netCDF as
// usual; the file is then reopened with HDF5, each chunk is converted,
// shuffled and deflated on worker threads exactly as HDF5's own filters would
// do it, and the compressed chunks are stored with H5Dwrite_chunk(). Values
// and file layout are identical to writing through netCDF.
//
// Compiled in only when the build found HDF5 >= 1.10.3 and zlib
// (HW3_HAVE_HDF5); otherwise, and for any variable whose filters, type or
// layout are not the plain shuffle + deflate this code reproduces, the caller
// writes with netCDF instead.

#include <string>
#include <vector>

namespace hw3 {

struct ChunkSrc {
    enum Kind { DOUBLE, SCHAR, INT };
    std::string name;   // variable name in the file
    const void* data;   // full array, row-major in the variable's dimension order
    Kind kind;
};

// True when this build can write chunks directly (HDF5 found at build time
// and the runtime library matches its headers), unless the environment
// variable HEATWAVE3_SERIAL_WRITE is set (to anything but 0).
bool h5_direct_available();

// Writes each variable of `vars` into the closed NetCDF-4 file `file`.
// Returns, per variable, whether it was written; variables not written must
// be written with netCDF. Throws on HDF5 errors after writing has begun.
std::vector<bool> h5_write_chunks(const std::string& file,
                                  const std::vector<ChunkSrc>& vars,
                                  int n_threads);

} // namespace hw3

#endif
