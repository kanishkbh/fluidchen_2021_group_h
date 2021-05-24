![](FluidchenLogo.png)

### What this solver is and is not.  
This is a finite differences solver for incompressible Navier-Stokes equations (in 2D) with heat transfer included.
It includes provision for doing coupled temperature simulations with the Boussinesq approximation : density is an affine function of temperature for the buyancy force, but the fluid is still assumed incompressible for the other parts of the equation. Furthermore, the fluid is Newtonian (viscosity is constant and isotropic), heat diffusivity is constant too, and heat dissipation due to friction is neglected. 

## Building
After cloning the repo, navigate to the directory 'fluidchen-skeleton' and build the files as follows:
```shell
mkdir build
cd build
cmake ..
make
```
### Guidelines to generate Geometry file (.pgm)
Create a PGM file that models the bitwise representation of the domain for simulation.

    Use the following legend as reference : 
    0 - Fluid element 
    1 - Inflow element 
    2 - Outflow element
    3 - Fixed wall
    4 - Fixed wall
    5 - Fixed wall
    6 - Fixed wall
    7 - Fixed wall
    8 - Moving wall 

Fixed walls with IDs in range 3-7 allow to have up to 5 different kinds of wall, with different temperature boundary conditions (for instance, an adiabatic wall and up to 4 different fixed temperatures). These boundary conditions are defined in the data file.
Note that the grid must include the ghost cells, so the dimensions must be (imax + 2) * (jmax + 2) with imax and jmax the values used in the .dat file.

### Guidelines to generate the parameter data file (.dat)
For a simulation that involves energy considerations (with negligible advection), inflow and outflows BCs : 

1) Define initial velocities UIN and VIN at the .dat file. (Inflow boundary conditions).
2) The outflow is dealt with using homogenous Neumann BCs.  
3) Define boundary temperatures wall_temp_i. Demperature at -1 is interpreted as adiabatic wall. (for i in range 3-8) 
4) For other simulation parameters refer to one of the example_cases. 
5) It is essential to include the name of the .ppm file, in local path relative to the .dat location. If it is not provided, a Lid-Driven Cavity with dimensions and time stepping as described is generated instead.

Other simulations should be easy to design from the provided exemples.

### Running
In order to run **Fluidchen**, the case file should be given as input parameter. Some default case files are located in the `example_cases` directory.  
After building, navigate back to the directory 'fluidchen-skeleton' and build the files as follows:

```shell
./build/fluidchen ./example_cases/LidDrivenCavity/LidDrivenCavity.dat
```

This will run the case file and create the output folder `./example_cases/LidDrivenCavity/LidDrivenCavity_Output` which holds the `.vtk` files of the solution. The output folder is created in the same location as your case file. Note that this may require write permissions in the given directory.
A log file in the output folder is also created and updated in real time, with file name CASE_NAME_log.txt. It is recommended to read it in command line with, `tail -f ./example_cases/CASE_NAME/CASE_NAME_Output/CASE_NAME_log.txt`. At the moment, it includes time steps counter, simulation time, time step, number of SOR iterations and pressure residual. Changing its extension to .csv makes it possible to open it with spreadsheet software (Excel, ...), and it can be read with Matlab (no matter the extension) with the `readtable` function.

If input file does not contain a geometry file, fluidchen will run lid-driven cavity case with given parameters.

### Predefined batches of input
The second worksheet deals with several cases that include energy and non-energy simulations. To run all the files together use **./run_all.sh **
**Disclaimer:** This may take several minutes, hence use with caution. 

### Time stepping

Time stepping in the `.dat` files is given by `dt` and `tau`. **If `tau` is postive, adapatative time stepping (with `tau` as safety factor) is used and `dt` is discarded. Otherwise, fixed time stepping is used with step size `dt`.**
  
