#ifndef HW3_LOG_H
#define HW3_LOG_H

// Message sinks for the language-neutral core. Declared here, defined by each
// language binding: the R package returns Rcpp::Rcout / Rcpp::Rcerr
// (heatwave3_init.cpp); the Python package returns streams that forward to
// sys.stdout / sys.stderr. Only call these from the main thread (e.g. the
// parallel_for report callback), never from worker threads.

#include <ostream>

namespace hw3 {

std::ostream& out();
std::ostream& err();

} // namespace hw3

#endif
