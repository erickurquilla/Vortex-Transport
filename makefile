CC = g++
CFLAGS = -std=c++11 -Wall
LDFLAGS = -I/usr/include/eigen3

# ------------------------------------------------------------------
# OpenMP (opt-in, disabled by default)
#
#   make ENABLE_OPENMP=1    OpenMP build (-fopenmp, defines VORTEX_USE_OPENMP)
#   make ENABLE_OPENMP=0    serial build (default)
#
# Run `make clean` when switching between the two modes, since both
# configurations produce object files with the same names.
# ------------------------------------------------------------------
ENABLE_OPENMP ?= FALSE
ifeq ($(ENABLE_OPENMP),TRUE)
CFLAGS  += -fopenmp -DVORTEX_USE_OPENMP
LDFLAGS += -fopenmp
endif

# ------------------------------------------------------------------
# Git commit hash baked into the binary at compile time (printed by main).
# Appends "-dirty" if the working tree has uncommitted changes.
# ------------------------------------------------------------------
GIT_COMMIT_HASH := $(shell git rev-parse --short=12 HEAD 2>/dev/null || echo unknown)$(shell git diff --quiet 2>/dev/null || echo -dirty)
GIT_FLAGS = -DGIT_COMMIT_HASH=\"$(GIT_COMMIT_HASH)\"

# ------------------------------------------------------------------
# Profiling build (opt-in, does not affect the default build)
#
#   make profile        (or: make PROFILE=1)
#
# builds a separate executable `Vortex-Transport-prof` from separate
# object files (*.prof.o), compiled with:
#   -O2                     keep realistic optimization
#   -g                      debug symbols so perf can resolve names
#   -fno-omit-frame-pointer accurate perf call stacks with `perf record -g`
#   -DENABLE_TIMERS         enable the Profiling.H timing/memory summary
#   -pthread                thread support for the (thread-safe) timers
# ------------------------------------------------------------------
PROFILE_FLAGS = -O2 -g -fno-omit-frame-pointer -DENABLE_TIMERS -pthread

SRCS = main.cpp Parameters.cpp Utilities.cpp Meshgeneration.cpp Element.cpp Lagrangebasis.cpp Quadraturerule.cpp Preevolve.cpp Evolve.cpp Numericalflux.cpp Timestepping.cpp
OBJS = $(SRCS:.cpp=.o)
PROF_OBJS = $(SRCS:.cpp=.prof.o)
EXEC = Vortex-Transport
PROF_EXEC = Vortex-Transport-prof

ifeq ($(PROFILE),1)
.DEFAULT_GOAL := $(PROF_EXEC)
endif

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

$(PROF_EXEC): $(PROF_OBJS)
	$(CC) $(CFLAGS) $(PROFILE_FLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

# objects also depend on the headers so edits to any .H trigger a rebuild
HDRS = $(wildcard *.H)

# force rebuild of main whenever the commit/dirty state may have changed
.git/HEAD .git/index:
	@true

%.o: %.cpp $(HDRS)
	$(CC) $(CFLAGS) $(LDFLAGS) -c $< -o $@

%.prof.o: %.cpp $(HDRS)
	$(CC) $(CFLAGS) $(PROFILE_FLAGS) $(LDFLAGS) -c $< -o $@

# only main needs the hash define; rebuild it when HEAD/index change
main.o: main.cpp $(HDRS) .git/HEAD .git/index
	$(CC) $(CFLAGS) $(GIT_FLAGS) $(LDFLAGS) -c main.cpp -o $@

main.prof.o: main.cpp $(HDRS) .git/HEAD .git/index
	$(CC) $(CFLAGS) $(PROFILE_FLAGS) $(GIT_FLAGS) $(LDFLAGS) -c main.cpp -o $@

.PHONY: profile clean
profile: $(PROF_EXEC)

clean:
	rm -f $(OBJS) $(PROF_OBJS) $(EXEC) $(PROF_EXEC)
	rm -rf grid
	rm -rf output
	rm -rf flamegraph.svg perf.data* perf_report.txt
	rm -rf plottools/*.pdf plottools/*.mp4
	rm -rf Vortex-Transport*
