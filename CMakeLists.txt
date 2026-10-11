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

install(TARGETS _core LIBRARY DESTINATION heatwave3)
