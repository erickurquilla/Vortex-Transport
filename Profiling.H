#ifndef PROFILING_H
#define PROFILING_H

// Lightweight, opt-in instrumentation for identifying hotspots before
// parallelizing the code.
//
// Usage, per function:
//   - PROFILE_DECLARE("name");  at file scope, registers the name so it shows
//     up in the summary table even if the function is never called (0 calls)
//   - PROFILE_SCOPE("name");    at the top of the function (or any scope),
//     accumulates the elapsed wall time and call count of that scope
//
// A summary table (calls, total time, % of wall time) plus the peak memory
// usage (VmPeak/VmHWM from /proc/self/status) is printed to stderr at exit.
//
// The macros expand to nothing unless ENABLE_TIMERS is defined, so the
// default build is unaffected. ENABLE_TIMERS is defined automatically by the
// profiling build: `make profile` (or `make PROFILE=1`).
//
// Note: timers are inclusive. A timed scope that calls other timed scopes
// (e.g. rk4 -> evolve_elem -> compute_numerical_flux) counts the same wall
// time in every enclosing timer, so the percentages do not sum to 100.
// For very hot micro-functions (millions of calls) the two clock reads per
// call add measurable overhead; trust perf for those, and the table for the
// coarser functions.

#ifdef ENABLE_TIMERS

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace profiling {

struct TimerData {
    double total_seconds;
    long long calls;
    TimerData() : total_seconds(0.0), calls(0) {}
};

class Registry {
public:
    static Registry& instance() {
        static Registry registry;
        return registry;
    }

    // returns the stats slot for `name`, creating it (0 calls) if new.
    // std::map guarantees the reference stays valid for the registry lifetime.
    TimerData& timer(const std::string& name) { return timers_[name]; }

    ~Registry() { report(); }

private:
    Registry() : program_start_(std::chrono::steady_clock::now()) {}
    Registry(const Registry&);            // non-copyable (C++11-compatible)
    Registry& operator=(const Registry&);

    // return e.g. "VmHWM:      33196 kB" from /proc/self/status
    static std::string proc_status_line(const std::string& key) {
        std::ifstream status_file("/proc/self/status");
        std::string line;
        while (std::getline(status_file, line)) {
            if (line.compare(0, key.size(), key) == 0) return line;
        }
        return key + "\tunavailable (no /proc/self/status)";
    }

    void report() const {
        double wall_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - program_start_).count();

        // sort timers by total time (largest first), then by name
        std::vector<std::pair<std::string, TimerData> > sorted(timers_.begin(), timers_.end());
        std::sort(sorted.begin(), sorted.end(), order_for_report);

        std::fprintf(stderr, "\n");
        std::fprintf(stderr, "==================== ENABLE_TIMERS profiling summary ====================\n");
        std::fprintf(stderr, "%-45s %12s %12s %8s\n", "function", "calls", "total (s)", "% wall");
        std::fprintf(stderr, "--------------------------------------------------------------------------\n");
        for (size_t i = 0; i < sorted.size(); ++i) {
            double percent = wall_seconds > 0.0 ? 100.0 * sorted[i].second.total_seconds / wall_seconds : 0.0;
            std::fprintf(stderr, "%-45s %12lld %12.3f %8.2f\n",
                         sorted[i].first.c_str(),
                         sorted[i].second.calls,
                         sorted[i].second.total_seconds,
                         percent);
        }
        std::fprintf(stderr, "--------------------------------------------------------------------------\n");
        std::fprintf(stderr, "wall time: %.3f s\n", wall_seconds);
        std::fprintf(stderr, "note: timers are inclusive, nested scopes overlap; %%s need not sum to 100\n");
        std::fprintf(stderr, "peak memory: %s | %s\n",
                     proc_status_line("VmHWM").c_str(),   // peak resident set size
                     proc_status_line("VmPeak").c_str()); // peak virtual memory
        std::fprintf(stderr, "==========================================================================\n");
    }

    static bool order_for_report(const std::pair<std::string, TimerData>& a,
                                 const std::pair<std::string, TimerData>& b) {
        if (a.second.total_seconds != b.second.total_seconds)
            return a.second.total_seconds > b.second.total_seconds;
        return a.first < b.first;
    }

    std::chrono::steady_clock::time_point program_start_;
    std::map<std::string, TimerData> timers_;
};

class ScopedTimer {
public:
    explicit ScopedTimer(TimerData& data)
        : data_(data), start_(std::chrono::steady_clock::now()) {}

    ~ScopedTimer() {
        data_.total_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_).count();
        data_.calls += 1;
    }

private:
    TimerData& data_;
    std::chrono::steady_clock::time_point start_;
};

} // namespace profiling

#define PROFILING_CONCAT_IMPL(a, b) a##b
#define PROFILING_CONCAT(a, b) PROFILING_CONCAT_IMPL(a, b)

// file-scope registration: makes `name` appear in the summary (with 0 calls)
// even if the function is never executed. Also forces the registry to be
// created during static initialization, before main() starts.
#define PROFILE_DECLARE(name) \
    static profiling::TimerData& PROFILING_CONCAT(profiling_declared_timer_, __LINE__) \
        __attribute__((unused)) = profiling::Registry::instance().timer(name)

// scope timing: the map lookup happens only once (function-local static),
// so each call only costs two clock reads and two additions.
#define PROFILE_SCOPE(name) \
    static profiling::TimerData& PROFILING_CONCAT(profiling_timer_data_, __LINE__) = \
        profiling::Registry::instance().timer(name); \
    profiling::ScopedTimer PROFILING_CONCAT(profiling_scoped_timer_, __LINE__)(PROFILING_CONCAT(profiling_timer_data_, __LINE__))

#else // ENABLE_TIMERS not defined: macros compile to nothing

#define PROFILE_SCOPE(name)
#define PROFILE_DECLARE(name)

#endif // ENABLE_TIMERS

#endif // PROFILING_H
