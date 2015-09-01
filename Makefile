CXXFLAGS = -g -O2 -Wall -std=c++0x -static-libgcc
LIB = -lm -lz -lpthread

GCC_VER := $(shell echo `$(CXX) -dumpversion | cut -f1-2 -d.`)

ifeq "4.5" "$(word 1, $(sort 4.5 $(GCC_VER)))"
	CXXFLAGS += -static-libstdc++
endif

all: samcmp

samcmp: samcmp.cpp
	$(CXX) $(CXXFLAGS) samcmp.cpp -o samcmp