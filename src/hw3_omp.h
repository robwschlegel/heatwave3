#ifndef HW3_OMP_H
#define HW3_OMP_H

// heatwave3 parallelism backend.
//
// Uses C++ std::thread rather than OpenMP. On macOS a dlopen-ed libomp can fail
// to allocate its per-worker thread-local storage (the gtid slot) in a process
// that already carries a large native footprint -- terra/GDAL, or the
// Positron/ark kernel that embeds R -- and then crashes in __kmp_suspend the
// instant a parallel region starts. std::thread uses pthreads and pthread-key
// TLS, so it is immune to that exhaustion, works regardless of which thread R
// runs on, and needs no special OpenMP runtime. Each parallel region joins its
// workers before returning, so it is also fork-safe under parallel::mclapply.

#include <thread>
#include <vector>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <algorithm>
#include <cstdlib>
#include <cstring>

namespace hw3 {

// ---- Thread-count management (no global side effects) ----------------------

inline int& hw3_threads_ref() {
    static int threads = 0;  // 0 = "auto", resolved lazily
    return threads;
}

inline int hw3_default_threads() {
    for (const char* var : {"R_HEATWAVE3_NUM_THREADS", "HEATWAVE3_NUM_THREADS"}) {
        const char* env = std::getenv(var);
        if (env && std::strlen(env) > 0) {
            int n = std::atoi(env);
            if (n > 0) return n;
        }
    }
    unsigned hc = std::thread::hardware_concurrency();
    int procs = (hc > 0) ? static_cast<int>(hc) : 1;
    int n = std::max(1, procs / 2);  // 50% of cores, leaving room for other work
    // CRAN allows at most 2 cores: R CMD check --as-cran sets
    // _R_CHECK_LIMIT_CORES_ (read as parallel::detectCores() reads it), and
    // check machines may set OMP_THREAD_LIMIT.
    const char* chk = std::getenv("_R_CHECK_LIMIT_CORES_");
    if (chk && std::strlen(chk) > 0 && std::strcmp(chk, "false") != 0)
        n = std::min(n, 2);
    const char* lim = std::getenv("OMP_THREAD_LIMIT");
    if (lim && std::atoi(lim) > 0) n = std::min(n, std::atoi(lim));
    return n;
}

// Thread count for a region with `n_iterations`. Throttles to at least ~100
// iterations per thread to avoid spawning threads for trivial work.
inline int get_threads(int n_iterations = 0) {
    int& t = hw3_threads_ref();
    if (t == 0) t = hw3_default_threads();
    int threads = t;
    if (n_iterations > 0) {
        int max_useful = std::max(1, n_iterations / 100);
        threads = std::min(threads, max_useful);
    }
    return std::max(1, threads);
}

inline void set_threads(int n) {
    hw3_threads_ref() = (n <= 0) ? hw3_default_threads() : n;
}

// std::thread workers are joined within every parallel region, so no heatwave3
// threads are alive across a fork; no atfork handler is needed.
inline void hw3_init_fork_handler() {}

// ---- Parallel for ----------------------------------------------------------
//
// Runs body(worker_id, i) for each i in [0, n) across `nthreads` workers
// (worker_id in [0, nthreads)). Items are dispatched dynamically (atomic
// counter) for load balance across the land/ocean cost mix.
//
// IMPORTANT: body() runs on worker threads and must NOT touch the R API.
// report(done_count) is invoked only on the calling (main) thread -- before,
// during (periodically), and at the end of the run -- so it may use
// hw3::out() (hw3_log.h) for progress reporting.
template <class Body, class Report>
inline void parallel_for(int n, int nthreads, Body body, Report report) {
    if (n <= 0) return;
    nthreads = std::max(1, std::min(nthreads, n));

    if (nthreads == 1) {
        int interval = std::max(1, n / 100);
        for (int i = 0; i < n; ++i) {
            body(0, i);
            if ((i + 1) % interval == 0 || i + 1 == n) report(i + 1);
        }
        return;
    }

    // The main thread reports progress every 80 ms, but is woken as soon as
    // the last item finishes, so short regions do not wait out the interval.
    std::atomic<int> next{0};
    std::atomic<int> done{0};
    std::mutex mtx;
    std::condition_variable finished;
    auto worker = [&](int wid) {
        int i;
        while ((i = next.fetch_add(1, std::memory_order_relaxed)) < n) {
            body(wid, i);
            if (done.fetch_add(1, std::memory_order_acq_rel) + 1 == n) {
                std::lock_guard<std::mutex> lock(mtx);
                finished.notify_one();
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    for (int w = 0; w < nthreads; ++w) pool.emplace_back(worker, w);

    int last = 0;
    for (;;) {
        int d = done.load(std::memory_order_acquire);
        if (d != last) { report(d); last = d; }
        if (d >= n) break;
        std::unique_lock<std::mutex> lock(mtx);
        finished.wait_for(lock, std::chrono::milliseconds(80), [&] {
            return done.load(std::memory_order_acquire) >= n;
        });
    }
    for (auto& th : pool) th.join();
}

} // namespace hw3

#endif
