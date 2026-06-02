#ifndef HW3_OMP_H
#define HW3_OMP_H

// Centralised OpenMP compatibility layer for heatwave3.
//
// Provides:
//   - Stub macros when _OPENMP is not defined (eliminates #ifdef guards)
//   - Private thread manager (never touches omp_set_num_threads global state)
//   - Fork safety via pthread_atfork
//
// Modelled on data.table's openmp-utils.c approach.

#ifdef _OPENMP
  #include <omp.h>
#else
  // Stubs so source code can call omp_* without #ifdef guards
  #define omp_get_max_threads()  1
  #define omp_get_thread_num()   0
  #define omp_get_num_procs()    1
#endif

#include <algorithm>
#include <cstdlib>
#include <cstring>

namespace hw3 {

// ---- Private thread manager ------------------------------------------------
//
// All parallel regions use:
//   #pragma omp parallel for schedule(static, 1) num_threads(hw3::get_threads(n))
//
// static (not dynamic) scheduling is deliberate: it emits only
// __kmpc_for_static_init, which every libomp provides, whereas
// schedule(dynamic) needs __kmpc_dispatch_deinit, absent from R's bundled
// libomp. The chunk size of 1 round-robins pixels for load balance.
//
// This avoids calling omp_set_num_threads(), which is a global side effect
// that changes thread counts for every OpenMP-using package in the R session.

inline int& hw3_threads_ref() {
    // Default: 0 means "auto" (resolved at first use)
    static int threads = 0;
    return threads;
}

// Resolve the default thread count (called once, lazily)
inline int hw3_default_threads() {
    // Check environment variable first
    const char* env = std::getenv("R_HEATWAVE3_NUM_THREADS");
    if (env && std::strlen(env) > 0) {
        int n = std::atoi(env);
        if (n > 0) return n;
    }
    // Default to 50% of available cores (polite, leaves room for other work)
    int procs = omp_get_num_procs();
    return std::max(1, procs / 2);
}

// Get the thread count for a parallel region with `n` iterations.
// Applies throttling: at least 100 iterations per thread to avoid overhead.
inline int get_threads(int n_iterations = 0) {
    int& t = hw3_threads_ref();
    if (t == 0) t = hw3_default_threads();

    int threads = t;

    // Throttle: don't spawn more threads than useful
    if (n_iterations > 0) {
        int max_useful = std::max(1, n_iterations / 100);
        threads = std::min(threads, max_useful);
    }
    return std::max(1, threads);
}

// Set the thread count (called from R via setHW3threads)
inline void set_threads(int n) {
    hw3_threads_ref() = (n <= 0) ? hw3_default_threads() : n;
}

// ---- Fork safety -----------------------------------------------------------
//
// OpenMP inside a forked child (e.g. parallel::mclapply) can deadlock.
// We register a pthread_atfork handler that drops to 1 thread before fork
// and restores afterward.

#if defined(_OPENMP) && !defined(_WIN32)
  #include <pthread.h>

  inline int& hw3_pre_fork_threads() {
      static int saved = 0;
      return saved;
  }

  inline void hw3_before_fork() {
      hw3_pre_fork_threads() = hw3_threads_ref();
      hw3_threads_ref() = 1;
  }

  inline void hw3_after_fork() {
      hw3_threads_ref() = hw3_pre_fork_threads();
  }

  inline void hw3_init_fork_handler() {
      pthread_atfork(hw3_before_fork, hw3_after_fork, nullptr);
  }
#else
  inline void hw3_init_fork_handler() {}
#endif

} // namespace hw3

#endif
