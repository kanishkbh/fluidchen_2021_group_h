---
title: Documentation
tags: []
---

# CFD Lab - Final Project (Group H)
- Boris Martin 
- Kanishk Bhatia 
- Rahul Manavalan
## Iterative Krylov solvers and Preconditioners  

### Motivation 

The "vanilla" fluidchen solver uses Gauss Seidel method for solving the Pressure Poisson Equation. While the GS is respectable in its own right, it does not scale well if the domain is refined as this is after all a stationary iterative scheme. 

In this project, we seek faster means to solve the PPE, as this would go a long way in exploring finer mesh sizes for the problems that we may setup using our fluidchen solver. And the best iterative scheme there is for symmetric systems is the CG method. 

In addition to exploring the traditional CG, we also implement some preconditioners for the same to reduce the number of iterations that might take for the desired convergence. 

### Additions to "vanilla" fluidchen 

- CG solver - Conjugate Gradient 
- SD solver - Steepest Descent 
- Jacobi Preconditioned CG 
- Richardson CG 
- Symmetric Gauss  Seidel preconditioned CG 
- IDLU preconditioned CG 

All of the implementations have parallel versions, where applicable. By this we mean that, the notion of a parallel GS preconditioner or a IDLU preconditioner does not really make sense. This is something for the user's consideration. 

### Implementation Details 

The implementation exploits the runtime polymorphism that is furnished by the PressureSolver type and builds on it to include more solvers. 

Note that the CG solvers and the preconditioners themselves could be made more efficient by using runtime polymorphism again, but this was something that was not attempted owing to a lack of time. 

### Building 

Building the code for the project should be no different than the one for the "vanilla" fluidchen solver. Nonethless for the sake of completion we mention 
        
        mkdir build 
        cd build 
        cmake .. 
        make 

### Interface 

The solver can be chosen by simply specifying the solver on the config file using the **pressure_solver** keyword. The options that the user has at disposal are : 

- CG 
- SD 
- CG_Jacobi 
- CG_Richardson
- CG_GS (Serial)
- CG_IDLU (Serial)




        pressure_solver {solver_param}
        # Example : 
        pressure_solver CG_GS

### Execution 

The first few solvers offer themselves to parallelization, but the preconditioned GS and IDLU are outliers and can only be run serially. 

Serial (Withing build) 
    
    ./fluidchen ../example_cases/Parallel_Cases/1_4/LidDrivenCavity.dat
    
Parallel (Within build)

    mpirun -np 4 ./fluidchen ../example_cases/Parallel_Cases/1_4/LidDrivenCavity.dat

### Why is CG better ? (Disclaimar : Involves theory) 

As mentioned, the SOR solver is a stationary iterative method and the rate of convergence is quite long. This is particularly true for fine refinements as we had seen in the instances in WS-3. 

CG is based on a different philosophy, wherein in order to solve the LSE, we minimize a certain cost function 

$ f(x) = frac{x^{T}Ax}{2} $
