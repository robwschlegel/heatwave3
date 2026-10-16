#!/usr/bin/env bash
# Build pinned HDF5 and netcdf-c as shared libraries into PREFIX, for the
# Python wheels (cibuildwheel's before-all on Linux and macOS). auditwheel /
# delocate then bundle them into each wheel.
#
#   python/tools/build_netcdf.sh PREFIX
#
# Features are kept to what heatwave3 needs (local NetCDF-3/4 files, deflate
# compression via HDF5): DAP/byte-range (curl), NCZarr, plugins, SZIP and the
# optional bzip2/blosc/zstd filters are off, so nothing else is linked. Checksums match
# the conda-forge recipes for the same versions. Re-running with the same
# versions already installed in PREFIX is a no-op.

set -euo pipefail

PREFIX=${1:?usage: build_netcdf.sh PREFIX}
HDF5_VERSION=1.14.6
HDF5_SHA256=e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b
NETCDF_VERSION=4.9.3
NETCDF_SHA256=990f46d49525d6ab5dc4249f8684c6deeaf54de6fec63a187e9fb382cc0ffdff

STAMP="$PREFIX/.heatwave3-netcdf-$HDF5_VERSION-$NETCDF_VERSION"
if [[ -f "$STAMP" ]]; then
  echo "HDF5 $HDF5_VERSION and netcdf-c $NETCDF_VERSION already in $PREFIX"
  exit 0
fi

JOBS=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

fetch() {  # url sha256 -> extracted into $WORK
  local url=$1 sha=$2 out
  out="$WORK/$(basename "$url")"
  curl -fsSL --retry 3 -o "$out" "$url"
  if command -v sha256sum >/dev/null; then
    echo "$sha  $out" | sha256sum -c -
  else
    echo "$sha  $out" | shasum -a 256 -c -
  fi
  tar -xzf "$out" -C "$WORK"
}

# Record where the libraries live (rpath on Linux, absolute install names on
# macOS) so auditwheel / delocate can resolve libnetcdf -> libhdf5_hl ->
# libhdf5 without LD_LIBRARY_PATH / DYLD_LIBRARY_PATH.
COMMON=(
  -DCMAKE_BUILD_TYPE=Release
  -DCMAKE_INSTALL_PREFIX="$PREFIX"
  -DCMAKE_INSTALL_LIBDIR=lib
  -DCMAKE_PREFIX_PATH="$PREFIX"
  -DCMAKE_INSTALL_RPATH="$PREFIX/lib"
  -DCMAKE_INSTALL_NAME_DIR="$PREFIX/lib"
  -DBUILD_SHARED_LIBS=ON
)
# On macOS never pick up Homebrew libraries: they target the runner's own
# macOS version, which breaks the wheel's deployment target.
if [[ "$(uname)" == "Darwin" ]]; then
  COMMON+=(-DCMAKE_IGNORE_PREFIX_PATH="/opt/homebrew;/usr/local")
fi

fetch "https://github.com/HDFGroup/hdf5/releases/download/hdf5_${HDF5_VERSION}/hdf5-${HDF5_VERSION}.tar.gz" "$HDF5_SHA256"
cmake -S "$WORK/hdf5-$HDF5_VERSION" -B "$WORK/build-hdf5" "${COMMON[@]}" \
  -DBUILD_STATIC_LIBS=OFF -DBUILD_TESTING=OFF \
  -DHDF5_BUILD_TOOLS=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_BUILD_HL_LIB=ON \
  -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_FORTRAN=OFF -DHDF5_BUILD_JAVA=OFF \
  -DHDF5_BUILD_DOC=OFF -DHDF5_ENABLE_Z_LIB_SUPPORT=ON -DHDF5_ENABLE_SZIP_SUPPORT=OFF
cmake --build "$WORK/build-hdf5" --parallel "$JOBS"
cmake --install "$WORK/build-hdf5"

fetch "https://github.com/Unidata/netcdf-c/archive/refs/tags/v${NETCDF_VERSION}.tar.gz" "$NETCDF_SHA256"
cmake -S "$WORK/netcdf-c-$NETCDF_VERSION" -B "$WORK/build-netcdf" "${COMMON[@]}" \
  -DHDF5_ROOT="$PREFIX" \
  -DNETCDF_ENABLE_NETCDF_4=ON -DNETCDF_ENABLE_HDF5=ON \
  -DNETCDF_ENABLE_DAP=OFF -DNETCDF_ENABLE_BYTERANGE=OFF -DNETCDF_ENABLE_NCZARR=OFF \
  -DNETCDF_ENABLE_PLUGINS=OFF -DNETCDF_ENABLE_FILTER_TESTING=OFF \
  -DNETCDF_ENABLE_FILTER_SZIP=OFF -DNETCDF_ENABLE_FILTER_BZ2=OFF \
  -DNETCDF_ENABLE_FILTER_BLOSC=OFF -DNETCDF_ENABLE_FILTER_ZSTD=OFF \
  -DNETCDF_ENABLE_LIBXML2=OFF -DNETCDF_ENABLE_HDF4=OFF \
  -DNETCDF_ENABLE_TESTS=OFF -DNETCDF_BUILD_UTILITIES=OFF -DNETCDF_ENABLE_EXAMPLES=OFF
cmake --build "$WORK/build-netcdf" --parallel "$JOBS"
cmake --install "$WORK/build-netcdf"

touch "$STAMP"
echo "Installed HDF5 $HDF5_VERSION and netcdf-c $NETCDF_VERSION into $PREFIX"
