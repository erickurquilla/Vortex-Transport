#ifndef PROFILING_H
#define PROFILING_H

// Lightweight, opt-in instrumentation reporting SELF and INCLUSIVE time per
// function, designed to stay correct once the code is parallelized
// (OpenMP/MPI ranks: each thread accumulates into its own storage, merged
// only at exit, so there is no lock contention on the hot path).
//
// Usage, per function:
//   - PROFILE_DECLARE("name");  at file scope, registers the name so it shows
//     up in the summary tables even if the function is never called (0 calls)
//   - PROFILE_SCOPE("name");    at the top of the function (or any scope)
//
// How self time is computed: scoped timers form a per-thread call stack.
// When a timer is destroyed it knows its full elapsed time (inclusive) and
// how much of that was spent inside nested timed scopes (child time, credited
// to it by its children on their destruction). Self = inclusive - children.
// The timer then credits its own elapsed time to its parent's child counter.
//
// Interpretation:
//   - SELF %       sums to ~100% across all functions in a single-threaded,
//                  fully instrumented run ("where is CPU time actually
//                  spent"). With multiple threads the sum can exceed 100% of
//                  wall time (threads run concurrently).
//   - INCLUSIVE %  does NOT sum to 100%: nested calls overlap (rk4 contains
//                  evolve_elem, which contains the compute_* kernels), and
//                  recursive calls double-count. This is expected. It answers
//                  "which high-level code path is expensive".
//
// The macros expand to nothing unless ENABLE_TIMERS is defined, so the
// default build is unaffected. ENABLE_TIMERS is defined automatically by the
// profiling build: `make profile` (or `make PROFILE=1`).
//
// Per-call overhead when enabled: two clock reads, a handful of additions,
// and two thread-local pointer accesses. The name lookup happens only once
// per function per thread (cached in a static thread_local pointer). Still,
// for micro-functions called millions of times (e.g. numerical_flux) the
// overhead is measurable; trust perf for those.

#ifdef ENABLE_TIMERS

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace profiling {

struct ProfileEntry {
    long long calls;
    double self_seconds;       // elapsed time minus time spent in nested timed scopes
    double incl_seconds;       // full elapsed time (children included)
    double min_call_seconds;   // min/max inclusive time of a single call
    double max_call_seconds;

    ProfileEntry()
        : calls(0), self_seconds(0.0), incl_seconds(0.0),
          min_call_seconds(std::numeric_limits<double>::infinity()),
          max_call_seconds(0.0) {}

    void merge(const ProfileEntry& other) {
        calls += other.calls;
        self_seconds += other.self_seconds;
        incl_seconds += other.incl_seconds;
        if (other.min_call_seconds < min_call_seconds) min_call_seconds = other.min_call_seconds;
        if (other.max_call_seconds > max_call_seconds) max_call_seconds = other.max_call_seconds;
    }
};

class ThreadState;

// Global registry: holds declared names, data merged from finished threads,
// and pointers to still-live thread states. Its mutex is only taken at
// registration/merge/report time, never while timing.
class GlobalRegistry {
public:
    static GlobalRegistry& instance() {
        static GlobalRegistry registry;
        return registry;
    }

    void declare(const std::string& name) {
        std::lock_guard<std::mutex> lock(mutex_);
        merged_[name]; // creates a default (0 calls) entry if not present
    }

    void register_thread(ThreadState* state) {
        std::lock_guard<std::mutex> lock(mutex_);
        live_threads_.push_back(state);
    }

    void merge_and_unregister(ThreadState* state); // defined below ThreadState

    ~GlobalRegistry() { report(); }

private:
    GlobalRegistry() : program_start_(std::chrono::steady_clock::now()) {}
    GlobalRegistry(const GlobalRegistry&);            // non-copyable (C++11-compatible)
    GlobalRegistry& operator=(const GlobalRegistry&);

    // return e.g. "VmHWM:      33196 kB" from /proc/self/status
    static std::string proc_status_line(const std::string& key) {
        std::ifstream status_file("/proc/self/status");
        std::string line;
        while (std::getline(status_file, line)) {
            if (line.compare(0, key.size(), key) == 0) return line;
        }
        return key + "\tunavailable (no /proc/self/status)";
    }

    typedef std::pair<std::string, ProfileEntry> NamedEntry;

    static bool by_self_time(const NamedEntry& a, const NamedEntry& b) {
        if (a.second.self_seconds != b.second.self_seconds)
            return a.second.self_seconds > b.second.self_seconds;
        return a.first < b.first;
    }

    static bool by_incl_time(const NamedEntry& a, const NamedEntry& b) {
        if (a.second.incl_seconds != b.second.incl_seconds)
            return a.second.incl_seconds > b.second.incl_seconds;
        return a.first < b.first;
    }

    void print_table(const std::vector<NamedEntry>& rows, double wall_seconds,
                     bool sorted_by_self) const {
        if (sorted_by_self) {
            std::fprintf(stderr, "%-45s %10s %11s %8s %11s %8s %14s\n",
                         "function", "calls", "self (s)", "self %", "incl (s)", "incl %", "avg self (ms)");
        } else {
            std::fprintf(stderr, "%-45s %10s %11s %8s %11s %8s %13s %13s\n",
                         "function", "calls", "incl (s)", "incl %", "self (s)", "self %", "min/call (ms)", "max/call (ms)");
        }
        std::fprintf(stderr, "--------------------------------------------------------------------------------------------------------------\n");
        for (size_t i = 0; i < rows.size(); ++i) {
            const std::string& name = rows[i].first;
            const ProfileEntry& e = rows[i].second;
            double self_pct = wall_seconds > 0.0 ? 100.0 * e.self_seconds / wall_seconds : 0.0;
            double incl_pct = wall_seconds > 0.0 ? 100.0 * e.incl_seconds / wall_seconds : 0.0;
            if (sorted_by_self) {
                double avg_self_ms = e.calls > 0 ? 1000.0 * e.self_seconds / e.calls : 0.0;
                std::fprintf(stderr, "%-45s %10lld %11.3f %8.2f %11.3f %8.2f %14.6f\n",
                             name.c_str(), e.calls, e.self_seconds, self_pct,
                             e.incl_seconds, incl_pct, avg_self_ms);
            } else {
                double min_ms = e.calls > 0 ? 1000.0 * e.min_call_seconds : 0.0;
                double max_ms = e.calls > 0 ? 1000.0 * e.max_call_seconds : 0.0;
                std::fprintf(stderr, "%-45s %10lld %11.3f %8.2f %11.3f %8.2f %13.6f %13.6f\n",
                             name.c_str(), e.calls, e.incl_seconds, incl_pct,
                             e.self_seconds, self_pct, min_ms, max_ms);
            }
        }
    }

    void report(); // defined below ThreadState (needs its definition)

    std::chrono::steady_clock::time_point program_start_;
    std::mutex mutex_;
    std::map<std::string, ProfileEntry> merged_;   // declared names + data from finished threads
    std::vector<ThreadState*> live_threads_;
};

// Per-thread accumulator. Created lazily on each thread's first timed scope,
// merged into the global registry when the thread finishes (or at report
// time for threads still alive at program exit, e.g. OpenMP pool threads).
class ThreadState {
public:
    static ThreadState& instance() {
        static thread_local ThreadState state;
        return state;
    }

    // pointer stays valid for the thread's lifetime (unordered_map nodes are
    // stable under rehash), so PROFILE_SCOPE caches it per thread
    ProfileEntry* entry(const std::string& name) { return &entries_[name]; }

    const std::unordered_map<std::string, ProfileEntry>& entries() const { return entries_; }

    ThreadState() { GlobalRegistry::instance().register_thread(this); }
    ~ThreadState() { GlobalRegistry::instance().merge_and_unregister(this); }

private:
    std::unordered_map<std::string, ProfileEntry> entries_;
};

inline void GlobalRegistry::merge_and_unregister(ThreadState* state) {
    std::lock_guard<std::mutex> lock(mutex_);
    live_threads_.erase(std::remove(live_threads_.begin(), live_threads_.end(), state),
                        live_threads_.end());
    typedef std::unordered_map<std::string, ProfileEntry>::const_iterator Iter;
    for (Iter it = state->entries().begin(); it != state->entries().end(); ++it) {
        merged_[it->first].merge(it->second);
    }
}

inline void GlobalRegistry::report() {
    double wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - program_start_).count();

    // snapshot: declared + finished-thread data, plus live threads (at exit
    // those should be quiescent; only main typically remains)
    std::map<std::string, ProfileEntry> total;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        total = merged_;
        for (size_t t = 0; t < live_threads_.size(); ++t) {
            typedef std::unordered_map<std::string, ProfileEntry>::const_iterator Iter;
            const std::unordered_map<std::string, ProfileEntry>& e = live_threads_[t]->entries();
            for (Iter it = e.begin(); it != e.end(); ++it) total[it->first].merge(it->second);
        }
    }

    std::vector<NamedEntry> rows(total.begin(), total.end());

    std::fprintf(stderr, "\n");
    std::fprintf(stderr, "======================== ENABLE_TIMERS profiling summary ========================\n");
    std::fprintf(stderr, "wall time: %.3f s\n", wall_seconds);
    std::fprintf(stderr, "peak memory: %s | %s\n",
                 proc_status_line("VmHWM").c_str(),   // peak resident set size
                 proc_status_line("VmPeak").c_str()); // peak virtual memory
    std::fprintf(stderr, "\n");

    // Table 1: sorted by SELF time — "where is CPU time actually spent".
    // Self percentages sum to ~100% in a fully instrumented single-threaded
    // run; with multiple threads the sum can exceed 100% of wall time.
    std::fprintf(stderr, "--- by SELF time (exclusive; where CPU time is actually spent) ---\n");
    std::sort(rows.begin(), rows.end(), by_self_time);
    print_table(rows, wall_seconds, true);
    std::fprintf(stderr, "\n");

    // Table 2: sorted by INCLUSIVE time — "which high-level path is
    // expensive". Inclusive percentages do NOT sum to 100%: nested calls
    // overlap and recursion double-counts. This is expected.
    std::fprintf(stderr, "--- by INCLUSIVE time (children included; which code path is expensive) ---\n");
    std::sort(rows.begin(), rows.end(), by_incl_time);
    print_table(rows, wall_seconds, false);
    std::fprintf(stderr, "==================================================================================\n");
}

// RAII timer forming a per-thread call stack (via a thread-local "current
// timer" pointer). On destruction it updates its function's entry and
// credits its elapsed time to the parent so the parent can compute self time.
class ScopedTimer {
public:
    explicit ScopedTimer(ProfileEntry* entry)
        : entry_(entry),
          parent_(current()),
          child_seconds_(0.0),
          start_(std::chrono::steady_clock::now()) {
        current() = this;
    }

    ~ScopedTimer() {
        double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_).count();
        double self = elapsed - child_seconds_;
        if (self < 0.0) self = 0.0; // guard against clock jitter
        entry_->calls += 1;
        entry_->incl_seconds += elapsed;
        entry_->self_seconds += self;
        if (elapsed < entry_->min_call_seconds) entry_->min_call_seconds = elapsed;
        if (elapsed > entry_->max_call_seconds) entry_->max_call_seconds = elapsed;
        current() = parent_;
        if (parent_) parent_->child_seconds_ += elapsed;
    }

private:
    static ScopedTimer*& current() {
        static thread_local ScopedTimer* current_timer = 0;
        return current_timer;
    }

    ProfileEntry* entry_;
    ScopedTimer* parent_;
    double child_seconds_;
    std::chrono::steady_clock::time_point start_;
};

// helper for PROFILE_DECLARE: registers a name at static-initialization time
struct DeclareRegistrar {
    explicit DeclareRegistrar(const char* name) {
        GlobalRegistry::instance().declare(name);
    }
};

} // namespace profiling

#define PROFILING_CONCAT_IMPL(a, b) a##b
#define PROFILING_CONCAT(a, b) PROFILING_CONCAT_IMPL(a, b)

// file-scope registration: makes `name` appear in the summary (with 0 calls)
// even if the function is never executed. Also forces the global registry to
// be created during static initialization, before main() starts.
#define PROFILE_DECLARE(name) \
    static profiling::DeclareRegistrar PROFILING_CONCAT(profiling_declared_, __LINE__)(name)

// scope timing: the name lookup happens only once per function per thread
// (static thread_local), so each call costs two clock reads, a few adds,
// and two thread-local pointer accesses.
#define PROFILE_SCOPE(name) \
    static thread_local profiling::ProfileEntry* PROFILING_CONCAT(profiling_entry_, __LINE__) = \
        profiling::ThreadState::instance().entry(name); \
    profiling::ScopedTimer PROFILING_CONCAT(profiling_scoped_timer_, __LINE__)(PROFILING_CONCAT(profiling_entry_, __LINE__))

#else // ENABLE_TIMERS not defined: macros compile to nothing

#define PROFILE_SCOPE(name)
#define PROFILE_DECLARE(name)

#endif // ENABLE_TIMERS

#endif // PROFILING_H
