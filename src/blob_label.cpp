// R wrapper for hw3::label_components() (blob_core.cpp), the 3D/4D
// connected-component labelling of heatwave-event masks. The mask is R's
// column-major 4D array m[i, j, kd, kt]; see blob_core.cpp for the layout
// and connectivity rules.
//
// Returned labels are dense, 1..K. Background voxels carry 0.

#include <Rcpp.h>
#include "blob_core.h"

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector label_components_cpp(IntegerVector mask,
                                   int nx, int ny, int nz, int nt,
                                   bool wrap_dateline = false,
                                   int connectivity = 6) {
  const R_xlen_t n = static_cast<R_xlen_t>(nx) * ny * nz * nt;
  if (mask.size() != n) {
    stop("mask length does not match nx * ny * nz * nt");
  }

  IntegerVector labels(n);  // initialised to 0
  int k = hw3::label_components(mask.begin(), labels.begin(),
                                nx, ny, nz, nt, wrap_dateline, connectivity);

  labels.attr("dim") = IntegerVector::create(nx, ny, nz, nt);
  labels.attr("n_components") = k;
  return labels;
}
