# If you do not have ROOT installed on your system, remove the following
# definition and replace it with anything needed to include and link to
# your own uniform random number generator.

INCLUDES = -I $(shell root-config --incdir)
LIBS = $(shell root-config --libs)
GCCFLAGS = -std=c++11 -g

all: ex

ex: ex.cc AdaptiveSampler.cc
	g++ $(GCCFLAGS) -o $@ $(INCLUDES) $^ $(LIBS)
