![](FluidchenLogo.png)

## Running

In order to run **Fluidchen**, the case file should be given as input parameter. Some default case files are located in the `example_cases` directory. If you installed **Fluidchen**, you can execute them from anywhere you want as
For Serial:

```shell
fluidchen /path/to/fluidchen/example_cases/LidDrivenCavity/LidDrivenCavity.dat
```

This will run the case file and create the output folder `/path/to/case/case_name_Output` which holds the `.vtk` files of the solution. The output folder is created in the same location as your case file. Note that this may require write permissions in the given directory. The folder will also contain a `case_name_log.txt` that records the history of time steps and SOR iterations required.

If input file does not contain a geometry file, fluidchen will run lid-driven cavity case with given parameters.

### Predefined batches of input

For the first worksheet, we've studied effects of changing viscosity, time stepping and mesh size. There is a subfolder for each parametric study in the `example_cases`, which contains a bash script `run_all.sh` so that you can queue up all simulations an go grab a coffe until the results show up.

### Time stepping

Time stepping in the `.dat` files is given by `dt` and `tau`. If `tau` is postive, adapatative time stepping (with `tau` as safety factor) is used and `dt` is discarded. Otherwise, fixed time stepping is used with step size `dt`.
