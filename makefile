CC = g++
CFLAGS = -std=c++11 -Wall
LDFLAGS = -I/usr/include/eigen3

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

%.o: %.cpp $(HDRS)
	$(CC) $(CFLAGS) $(LDFLAGS) -c $< -o $@

%.prof.o: %.cpp $(HDRS)
	$(CC) $(CFLAGS) $(PROFILE_FLAGS) $(LDFLAGS) -c $< -o $@

.PHONY: profile clean
profile: $(PROF_EXEC)

clean:
	rm -f $(OBJS) $(PROF_OBJS) $(EXEC) $(PROF_EXEC)
	rm -rf grid
	rm -rf output
	rm -rf flamegraph.svg perf.data* perf_report.txt
	rm -rf plottools/*.pdf plottools/*.mp4
	rm -rf Vortex-Transport*
