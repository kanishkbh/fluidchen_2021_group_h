#!/bin/bash
mpirun -np 1 ../../../build/fluidchen ./32_32/LidDrivenCavity.dat >> logs.dat 
mpirun -np 2  ../../../build/fluidchen ./64_32/LidDrivenCavity.dat >> logs.dat 
mpirun -np 3 ../../../build/fluidchen ./96_32/LidDrivenCavity.dat >> logs.dat 
mpirun -np 4 ../../../build/fluidchen ./128_32/LidDrivenCavity.dat >> logs.dat 
