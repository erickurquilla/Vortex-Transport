#ifndef PROFILING_H
#define PROFILING_H

// Lightweight, opt-in instrumentation for identifying hotspots before
// parallelizing the code.
//
// Usage: place PROFILE_SCOPE("some name"); at the top of a function (or any
// scope). The elapsed wall time of that scope is accumulated per name and a
// summary table (total time, % of wall time, call count) plus the peak memory
// usage (VmPeak/VmHWM from /proc/self/status) is printed to stderr at exit.
//
// The macros expand to nothing unless ENABLE_TIMERS is defined, so the
// default build is unaffected. ENABLE_TIMERS is defined automatically by the
// profiling build: `make profile` (or `make PROFILE=1`).
//
// Note: timers are inclusive. A timed scope that calls other timed scopes
// (e.g. rk4 -> evolve_elem -> compute_numerical_flux) counts the same wall
// time in every enclosing timer, so the percentages do not sum to 100.

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

    void add(const char* name, double seconds) {
        TimerData& data = timers_[name];
        data.total_seconds += seconds;
        data.calls += 1;
    }

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

        // sort timers by total time, largest first
        std::vector<std::pair<std::string, TimerData> > sorted(timers_.begin(), timers_.end());
        std::sort(sorted.begin(), sorted.end(), by_total_time_descending);

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

    static bool by_total_time_descending(const std::pair<std::string, TimerData>& a,
                                         const std::pair<std::string, TimerData>& b) {
        return a.second.total_seconds > b.second.total_seconds;
    }

    std::chrono::steady_clock::time_point program_start_;
    std::map<std::string, TimerData> timers_;
};

class ScopedTimer {
public:
    explicit ScopedTimer(const char* name)
        : name_(name), start_(std::chrono::steady_clock::now()) {}

    ~ScopedTimer() {
        double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_).count();
        Registry::instance().add(name_, seconds);
    }

private:
    const char* name_;
    std::chrono::steady_clock::time_point start_;
};

} // namespace profiling

#define PROFILING_CONCAT_IMPL(a, b) a##b
#define PROFILING_CONCAT(a, b) PROFILING_CONCAT_IMPL(a, b)
#define PROFILE_SCOPE(name) \
    profiling::ScopedTimer PROFILING_CONCAT(profiling_scoped_timer_, __LINE__)(name)

#else // ENABLE_TIMERS not defined: macros compile to nothing

#define PROFILE_SCOPE(name)

#endif // ENABLE_TIMERS

#endif // PROFILING_H
