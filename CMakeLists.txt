# Python build of heatwave3 (driven by scikit-build-core via pyproject.toml).
# Compiles the language-neutral C++ core in src/ plus the pybind11 bindings in
# python/src/. The R package's Rcpp files in src/ are not compiled here.

cmake_minimum_required(VERSION 3.18...3.30)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(Threads REQUIRED)

pybind11_add_module(_core
  python/src/bindings.cpp
  src/netcdf_io.cpp
  src/climatology.cpp
  src/event_detect.cpp
  src/hw3_api.cpp
  src/blob_core.cpp
  src/hw3_h5chunk.cpp
)
target_include_directories(_core PRIVATE src)
target_compile_definitions(_core PRIVATE HW3_VERSION="${SKBUILD_PROJECT_VERSION_FULL}")
target_link_libraries(_core PRIVATE Threads::Threads)
# Keep the path to libnetcdf in the installed module's rpath, so the wheel
# repair tools (auditwheel, delocate) can find the library to bundle.
set_target_properties(_core PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

# libnetcdf: netCDFConfig.cmake first, then nc-config, then pkg-config
# (mirrors the R package's configure script).
find_package(netCDF CONFIG QUIET)
if(netCDF_FOUND)
  message(STATUS "netCDF ${netCDF_VERSION} via netCDFConfig.cmake")
  target_link_libraries(_core PRIVATE netCDF::netcdf)
else()
  find_program(NC_CONFIG nc-config)
  if(NC_CONFIG)
    message(STATUS "netCDF via ${NC_CONFIG}")
    execute_process(COMMAND ${NC_CONFIG} --cflags OUTPUT_VARIABLE NC_CFLAGS
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${NC_CONFIG} --libs OUTPUT_VARIABLE NC_LIBS
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    separate_arguments(NC_CFLAGS UNIX_COMMAND "${NC_CFLAGS}")
    separate_arguments(NC_LIBS UNIX_COMMAND "${NC_LIBS}")
    target_compile_options(_core PRIVATE ${NC_CFLAGS})
    target_link_libraries(_core PRIVATE ${NC_LIBS})
  else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET netcdf)
    message(STATUS "netCDF ${NETCDF_VERSION} via pkg-config")
    target_link_libraries(_core PRIVATE PkgConfig::NETCDF)
  endif()
endif()

# Optional: compress gridded outputs on all threads (src/hw3_h5chunk.cpp),
# which needs HDF5 >= 1.10.3 and zlib. Use the HDF5 that libnetcdf uses (the
# wheels point HDF5_ROOT at the prefix built by python/tools/build_netcdf.sh)
# and check that its headers match the library; otherwise outputs are written
# by netCDF on one thread. Disable with -DHW3_HDF5=OFF.
option(HW3_HDF5 "Compress gridded outputs on all threads (HDF5 >= 1.10.3 + zlib)" ON)
set(HW3_HDF5_USED OFF)
if(HW3_HDF5)
  find_package(ZLIB QUIET)
  find_package(HDF5 QUIET COMPONENTS C)
  if(HDF5_FOUND AND ZLIB_FOUND)
    # Only libhdf5 itself (FindHDF5 also lists its dependencies, e.g. curl)
    if(HDF5_C_LIBRARY_hdf5)
      set(HW3_H5_LIBS ${HDF5_C_LIBRARY_hdf5})
    else()
      set(HW3_H5_LIBS ${HDF5_C_LIBRARIES})
    endif()
    if(NOT HW3_H5_LIBS)
      set(HW3_H5_LIBS ${HDF5_LIBRARIES})
    endif()
    include(CheckCXXSourceRuns)
    set(CMAKE_REQUIRED_INCLUDES ${HDF5_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${HW3_H5_LIBS} ZLIB::ZLIB)
    set(CMAKE_REQUIRED_DEFINITIONS ${HDF5_DEFINITIONS})
    check_cxx_source_runs([=[
#include <hdf5.h>
#include <zlib.h>
int main() {
  unsigned a, b, c;
  auto wc = &H5Dwrite_chunk;
  auto cz = &compress2;
  if (!wc || !cz || H5get_libversion(&a, &b, &c) < 0) return 1;
  return (a == H5_VERS_MAJOR && b == H5_VERS_MINOR && c == H5_VERS_RELEASE) ? 0 : 2;
}]=] HW3_HDF5_WORKS)
    unset(CMAKE_REQUIRED_INCLUDES)
    unset(CMAKE_REQUIRED_LIBRARIES)
    unset(CMAKE_REQUIRED_DEFINITIONS)
    if(HW3_HDF5_WORKS)
      target_compile_definitions(_core PRIVATE HW3_HAVE_HDF5 ${HDF5_DEFINITIONS})
      target_include_directories(_core PRIVATE ${HDF5_INCLUDE_DIRS})
      target_link_libraries(_core PRIVATE ${HW3_H5_LIBS} ZLIB::ZLIB)
      set(HW3_HDF5_USED ON)
    endif()
  endif()
endif()
if(HW3_HDF5_USED)
  message(STATUS "HDF5 ${HDF5_VERSION} for parallel compression: ${HW3_H5_LIBS}")
else()
  message(STATUS "HDF5 for parallel compression: not used (gridded outputs are compressed on one thread)")
endif()

install(TARGETS _core LIBRARY DESTINATION heatwave3)
