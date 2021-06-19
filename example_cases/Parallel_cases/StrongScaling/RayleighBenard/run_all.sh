#!/bin/bash 
mpirun -np 1 ../../../../build/fluidchen ./1_1/RayleighBenard.dat >> logs.dat 
mpirun -np 2 ../../../../build/fluidchen ./2_1/RayleighBenard.dat >> logs.dat 
mpirun -np 3 ../../../../build/fluidchen ./3_1/RayleighBenard.dat >> logs.dat 
mpirun -np 4 ../../../../build/fluidchen ./4_1/RayleighBenard.dat >> logs.dat 
