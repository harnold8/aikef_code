SRCS := $(wildcard src/*.cpp)
OBJS := $(SRCS:src/%.cpp=obj/%.o)

# TYPE is used to switch between compilers
# usage: make TYPE=home_icc
#        make TYPE=home_gcc
#        make TYPE=home_clang
#        etc.

#TYPE := home_gcc
TYPE := home_icc

LIB = -L lib_64/lib_silo lib_64/lib_silo/libsiloh5.a lib_64/lib_silo/libhdf5.a lib_64/lib_silo/libsz.a lib_64/lib_silo/libz.a lib_64/lib_gsl/libgsl.a -ldl -lm #-lstdc++

ifeq ($(TYPE),home_gcc)
	OMPI_CC := gcc
	OMPI_CXX := g++
	ARCH := -march=native -mtune=native -m64 -mfpmath=sse -msse4 -msse4.1 -msse4.2 
	OPTIMIZE := $(ARCH) -O2 -fno-builtin -fno-rtti -ffunction-sections -fno-signed-zeros -ftree-vectorize 
	STANDARD := # -ansi -pedantic -std=c++0x 
else ifeq ($(TYPE),home_gcc_O3)
	OMPI_CC := gcc
	OMPI_CXX := g++
	ARCH := -march=native -mtune=native -m64 -mfpmath=sse -msse4 -msse4.1 -msse4.2 
	OPTIMIZE := $(ARCH) -O3 -mvzeroupper -fno-builtin -fno-rtti -ffunction-sections -fno-signed-zeros -ftree-vectorize 
	STANDARD := -Wall -Wextra -ansi -pedantic -std=c++11 
else ifeq ($(TYPE),home_icc)
	OMPI_CC := icc
	OMPI_CXX := icpc
	ARCH := -march=native -mtune=native -m64 -mfpmath=sse -msse4.2 #-axsse4.2 
	OPTIMIZE := $(ARCH) -O2 -fno-builtin -fno-rtti -ffunction-sections -fno-signed-zeros -ftree-vectorize
#	STANDARD := -Wall -Wextra -ansi -pedantic -std=c++11 
else ifeq ($(TYPE),home_icc_O3)
	OMPI_CC := icc
	OMPI_CXX := icpc
	ARCH := -march=native -mtune=native -m64 -mfpmath=sse -msse4.2 -axsse4.2 
	OPTIMIZE := $(ARCH) -O3 -fno-builtin -fno-rtti -ffunction-sections
	STANDARD := -Wall -Wextra -ansi -pedantic -std=c++11 
else ifeq ($(TYPE),home_clang)
	OMPI_CC := clang
	OMPI_CXX := clang++
	ARCH := -march=native -mtune=native -m64 -msse4 -msse4.1 -msse4.2 
	OPTIMIZE := $(ARCH) -O2 -fno-builtin -fno-rtti -ffunction-sections -fno-signed-zeros 
	STANDARD := -Wall -Wextra -pedantic -std=c++11 
else ifeq ($(TYPE),home_clang_O3)
	OMPI_CC := clang
	OMPI_CXX := clang++
	ARCH := -march=native -mtune=native -m64 -msse4 -msse4.1 -msse4.2 
	OPTIMIZE := $(ARCH) -O3 -fno-builtin -fno-rtti -ffunction-sections -fno-signed-zeros 
	STANDARD := -Wall -Wextra -pedantic -std=c++11 
else
$(error Variable TYPE not recognised.)
endif

export OMPI_CC
export OMPI_CXX

COMPILER := mpicxx
INCLUDE  := -Isrc/vectorclass


bin/aikef_mpi: $(OBJS)
	$(COMPILER) $(INCLUDE) $(OBJS) -o bin/aikef_mpi $(OPTIMIZE) $(STANDARD) $(LIB)

obj/%.o: src/%.cpp
	$(COMPILER) $(INCLUDE) -c $< -o $@ $(OPTIMIZE) $(STANDARD)


clean:
	rm -f $(OBJS) bin/aikef_mpi

rmstate:
	rm -r bin/State
	rm -r bin/Last_State
	rm $(wildcard bin/process_*.log)

rmdata:
	rm -r  data/gnuplot	
	mkdir data/gnuplot

	rm -r  data/trajectories
	mkdir data/trajectories

	rm -r data/native
	mkdir data/native

	rm -r  data/silo
	mkdir data/silo

	rm -r  data/silo_3D
	mkdir data/silo_3D

	rm -r data/vtk
	mkdir data/vtk

	rm -r data/pTracer
	mkdir data/pTracer

	rm -r data/particle_tracks
	mkdir data/particle_tracks

	rm -r data/particle_detector
	mkdir data/particle_detector

	rm -r data/lineout
	mkdir data/lineout

	rm -rf bin/State
	rm -rf bin/Last_State

	rm -r bin/process_0.log
	rm -r bin/*.[o][0123456789][0123456789][01234567890]*

rmdir:
	rm -r bin data obj

info:
	@echo 
	@echo COMPILER = $(COMPILER)
	@echo OMPI_CC  = $(OMPI_CC)
	@echo OMPI_CXX = $(OMPI_CXX)
	@echo 
	@echo INCLUDE  = $(INCLUDE)
	@echo 
	@echo OPTIMIZE = $(OPTIMIZE)
	@echo 
	@echo STANDARD = $(STANDARD)
	@echo 
	@echo LIB  = $(LIB)
	@echo 
	@echo SRCS = $(SRCS)
	@echo 
	@echo OBJS = $(OBJS)
	@echo 
	@echo 

mkdir:
	-mkdir bin
	-mkdir data
	-mkdir obj
	-mkdir data/lineout
	-mkdir data/particle_detector
	-mkdir data/particle_tracks
	-mkdir data/uniform_output
	-mkdir data/silo
	-mkdir data/silo_3D
	-mkdir data/trajectories
