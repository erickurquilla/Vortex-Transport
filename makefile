CC = g++
CFLAGS = -std=c++11 -Wall
LDFLAGS = -I/usr/include/eigen3

SRCS = main.cpp Parameters.cpp Utilities.cpp Meshgeneration.cpp Element.cpp Lagrangebasis.cpp Quadraturerule.cpp Preevolve.cpp Evolve.cpp Numericalflux.cpp Timestepping.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = Vortex-Transport

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LIBS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJS) $(EXEC)
	rm -r grid
	rm -r output