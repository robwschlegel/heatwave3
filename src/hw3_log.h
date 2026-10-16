#ifndef HW3_LOG_H
#define HW3_LOG_H

// Message sinks for the language-neutral core. out_sink()/err_sink() are
// declared here and defined by each language binding: the R package returns
// Rcpp::Rcout / Rcpp::Rcerr (heatwave3_init.cpp); the Python package returns
// streams that forward to sys.stdout / sys.stderr. Only call out()/err() from
// the main thread (e.g. the parallel_for report callback), never from worker
// threads.

#include <ostream>

namespace hw3 {

std::ostream& out_sink();
std::ostream& err_sink();

// Per-thread quiet flag, set by the language wrappers for the duration of a
// call (quiet = TRUE). Messages are only written from the calling thread, so
// a thread-local flag silences exactly the call that set it.
inline bool& quiet_flag() {
    thread_local bool quiet = false;
    return quiet;
}

// Informational messages (progress, files read and written): discarded when
// quiet.
inline std::ostream& out() {
    if (quiet_flag()) {
        thread_local std::ostream discard(nullptr);
        return discard;
    }
    return out_sink();
}

// Warnings (e.g. skipped files): always reported.
inline std::ostream& err() {
    return err_sink();
}

} // namespace hw3

#endif
