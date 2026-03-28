# HPSC 2026 Particle Simulator Makefile
FC = gfortran
FFLAGS = -O3 -fopenmp -Wall

# This is the 'all' target - what happens when you just type 'make'
all: simulator

# We list the physics module first so its .mod file is created before main.f90 needs it
simulator: src/physics_data.f90 src/main.f90
	$(FC) $(FFLAGS) src/physics_data.f90 src/main.f90 -o simulator

# This cleans up the directory
clean:
	rm -f simulator *.mod
