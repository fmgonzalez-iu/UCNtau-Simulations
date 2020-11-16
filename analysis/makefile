OPT=-O3
CPP=mpic++
CPPFLAGS=$(OPT) -std=c++11
CC=mpicc
CFLAGS=$(OPT)
CPPFLAGS=$(OPT) -std=c++11 -fopenmp `root-config --cflags`
LFLAGS=-L/N/soft/cle7/gsl/2.6/lib
CPPLFLAGS=-fopenmp -lgsl -lgslcblas `root-config --libs`

MPICPP := $(shell command -v mpic++ 2> /dev/null)
MPICC := $(shell command -v mpicc 2> /dev/null)

ifndef MPICPP
	CPP=CC
endif
ifndef MPICC
	CC=cc
endif


all: analyze

debug: OPT=-g
debug: all


bin/xorshift.o: src/xorshift.c inc/xorshift.h
	$(CC) $(CFLAGS) -c -o bin/xorshift.o src/xorshift.c

bin/chisq_spectrum.o: src/chisq_spectrum.cpp inc/chisq_spectrum.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/chisq_spectrum.o src/chisq_spectrum.cpp

bin/file_reader.o: src/file_reader.cpp inc/file_reader.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/file_reader.o src/file_reader.cpp

bin/generate_weighting.o: src/generate_weighting.cpp inc/generate_weighting.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/generate_weighting.o src/generate_weighting.cpp

bin/quant_refl.o: src/quant_refl.cpp inc/quant_refl.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/quant_refl.o src/quant_refl.cpp

bin/stitcher.o: src/stitcher.cpp inc/stitcher.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/stitcher.o src/stitcher.cpp
	
bin/surface_multilayer.o: src/surface_multilayer.cpp inc/surface_multilayer.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/surface_multilayer.o src/surface_multilayer.cpp
bin/2det_plotting.o: src/2det_plotting.cpp inc/2det_plotting.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/2det_plotting.o src/2det_plotting.cpp

bin/analyze.o: analyze.cpp 
	$(CPP) $(CPPFLAGS) -c -o bin/analyze.o analyze.cpp

analyze: bin/analyze.o bin/xorshift.o bin/chisq_spectrum.o bin/file_reader.o bin/generate_weighting.o bin/quant_refl.o bin/stitcher.o bin/surface_multilayer.o bin/2det_plotting.o
	$(CPP) $(LFLAGS) -o analyze bin/analyze.o bin/xorshift.o bin/chisq_spectrum.o bin/file_reader.o bin/generate_weighting.o bin/quant_refl.o bin/stitcher.o bin/surface_multilayer.o bin/2det_plotting.o $(CPPLFLAGS)

clean:
	find ./bin/ -type f -name '*.o' -delete
	rm -rf analyze
