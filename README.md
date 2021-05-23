![](FluidchenLogo.png)

### What this (solver) is and is not.  
This is a solver for incompressible Navier Stokes equations that is not essentially dimensionless. (More on the equation to follow.)
It includes provision for doing coupled temperature simulation, in the domain of flow where advection is negligible. 
Any scenario besides the ones mentioned, is beyond the scope of this solver. 

## Building
After cloning the repo, navigate to the directory 'fluidchen-skeleton' and build the files as follows:
```shell
mkdir build
cd build
cmake ..
make
```
### Guidelines to generate Geometry file
Create a PGM file that models the bitwise representation of the domain for simulation.

    Use the following legend as reference : 
    0 - Fluid element 
    1 - Inflow element 
    2 - Outflow element
    3 - Fixed wall (with a predefined temperature or not)
    4 - Fixed wall (with a predefined temperature or not)
    5 - Adiabatic wall 
    6 - Miscellaneous wall 1  
    7 - Miscellaneous wall 2 
    8 - Moving wall 

### Guidelines to generate the parameter data file 
For a simulation that involves energy considerations (with negligible advection), inflow and outflows BCs : 

1) Define initial velocities UIN and VIN at the .dat file. (Inflow boundary conditions).
2) The outflow is dealt with using homogenous Neumann BCs.  
3) Define initial temperatures wall_temp_i. 
4) For other simulation parameters refer to one of the example_cases. 
5) It is essential to include the filename.ppm in the .dat file. (See example)

Other simulations that are subsets of the described case should be easily adaptable. (When in doubt look at the examples.) 

### Running
In order to run **Fluidchen**, the case file should be given as input parameter. Some default case files are located in the `example_cases` directory.  
After building, navigate back to the directory 'fluidchen-skeleton' and build the files as follows:

```shell
./build/fluidchen ./example_cases/LidDrivenCavity/LidDrivenCavity.dat
```

This will run the case file and create the output folder `./example_cases/LidDrivenCavity/LidDrivenCavity_Output` which holds the `.vtk` files of the solution. The output folder is created in the same location as your case file. Note that this may require write permissions in the given directory.

If input file does not contain a geometry file, fluidchen will run lid-driven cavity case with given parameters.

### Predefined batches of input
The second worksheet deals with several cases that include energy and non-energy simulations. To run all the files together use **./run_all.sh **
**Disclaimar:** This may take several minutes, hence use with caution. 

### Time stepping

Time stepping in the `.dat` files is given by `dt` and `tau`. **If `tau` is postive, adapatative time stepping (with `tau` as safety factor) is used and `dt` is discarded. Otherwise, fixed time stepping is used with step size `dt`.**
  
