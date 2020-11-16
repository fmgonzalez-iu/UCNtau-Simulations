OPT=-O3
CPP=mpic++
CPPFLAGS=$(OPT) -std=c++11
CC=mpicc
CFLAGS=$(OPT)

MPICPP := $(shell command -v mpic++ 2> /dev/null)
MPICC := $(shell command -v mpicc 2> /dev/null)

ifndef MPICPP
	CPP=CC
endif
ifndef MPICC
	CC=cc
endif

all: sim

debug: OPT=-g
debug: all

bin/fields_nate.o: src/fields_nate.c inc/fields_nate.h
	$(CC) $(CFLAGS) -c -o bin/fields_nate.o src/fields_nate.c

bin/xorshift.o: src/xorshift.c inc/xorshift.h
	$(CC) $(CFLAGS) -c -o bin/xorshift.o src/xorshift.c

bin/quant_refl.o: src/quant_refl.cpp inc/quant_refl.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/quant_refl.o src/quant_refl.cpp

bin/symplectic.o: src/symplectic.cpp inc/symplectic.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/symplectic.o src/symplectic.cpp

bin/geometry.o: src/geometry.cpp inc/geometry.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/geometry.o src/geometry.cpp

bin/track_gen.o: src/track_gen.cpp inc/track_gen.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/track_gen.o src/track_gen.cpp

bin/trackUCN.o: src/trackUCN.cpp inc/trackUCN.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/trackUCN.o src/trackUCN.cpp

bin/lyap.o: src/lyap.cpp inc/lyap.hpp
	$(CPP) $(CPPFLAGS) -c -o bin/lyap.o src/lyap.cpp

bin/sim.o: sim.cpp
	$(CPP) $(CPPFLAGS) -c -o bin/sim.o sim.cpp

sim: bin/sim.o bin/trackUCN.o bin/track_gen.o bin/geometry.o bin/symplectic.o bin/quant_refl.o bin/xorshift.o bin/fields_nate.o bin/lyap.o
	$(CPP) $(LFLAGS) -o sim bin/sim.o bin/trackUCN.o bin/track_gen.o bin/geometry.o bin/symplectic.o bin/quant_refl.o bin/xorshift.o bin/fields_nate.o bin/lyap.o

bin/find_min.o: find_min.c
	$(CC) $(CFLAGS) -c find_min.c -o bin/find_min.o


find_min: bin/find_min.o bin/fields_nate.o
	$(CC) $(CFLAGS) bin/fields_nate.o bin/find_min.o -o find_min $(LFLAGS) -lgsl

clean:
	find ./bin/ -type f -name '*.o' -delete
	rm -rf sim find_min
