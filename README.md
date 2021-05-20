![](FluidchenLogo.png)

## Building
After cloning the repo, navigate to the directory 'fluidchen-skeleton' and build the files as follows:
```shell
mkdir build
cd build
cmake ..
make
```

## Running

In order to run **Fluidchen**, the case file should be given as input parameter. Some default case files are located in the `example_cases` directory.  
After building, navigate back to the directory 'fluidchen-skeleton' and build the files as follows:

```shell
./build/fluidchen ./example_cases/LidDrivenCavity/LidDrivenCavity.dat
```

This will run the case file and create the output folder `./example_cases/LidDrivenCavity/LidDrivenCavity_Output` which holds the `.vtk` files of the solution. The output folder is created in the same location as your case file. Note that this may require write permissions in the given directory. A log file in the output folder is also created and updated in real time, with file name CASE_NAME_log.txt. It is recommended to read it in command line with, `tail -f ./example_cases/CASE_NAME/CASE_NAME_Output/CASE_NAME_log.txt`. At the moment, it includes time steps counter, simulation time, time step, number of SOR iterations and pressure residual.

If input file does not contain a geometry file, fluidchen will run lid-driven cavity case with given parameters.

### Predefined batches of input

For the first worksheet, we've studied effects of changing viscosity, time stepping and mesh size. There is a subfolder for each parametric study in the `example_cases`, which contains a bash script `run_all.sh` so that you can queue up all simulations an go grab a coffe until the results show up.

### Time stepping

Time stepping in the `.dat` files is given by `dt` and `tau`. **If `tau` is postive, adapatative time stepping (with `tau` as safety factor) is used and `dt` is discarded. Otherwise, fixed time stepping is used with step size `dt`.**

### Defining Geometry File 

As prescribed, we use a pgm file to define the geometry of the domain that we simulate. Use the following as a reference should you choose to define your geometry for the simulation domain.  

FLUID           - 0 
INFLOW          - 1 
OUTFLOW         - 2
FIXED WALL      - 3 
INSULATED WALL  - 4 
COLD WALL       - 5 
HEATED WALL     - 6 
FIXEDWALL OTHER - 7
MOVING WALL     - 8 