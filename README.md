# IMAU-ICE

Contact: Tijn Berends (c.j.berends@uu.nl), Jorjo Bernales (j.a.bernalesconcha@uu.nl)

This repository contains the source code of IMAU-ICE, which can be found in the `src/` directory.
A few template files for compiling the code can be found in `templates/`
A few example model configuration files can be found in `config-files/`

## Getting the model

First you'll need a local copy of the model in your machine. You have two main options:
1. Use Git: open your beloved terminal, go to wherever you want the model to reside, and type `git clone https://github.com/IMAU-paleo/IMAU-ICE.git`. Assuming you got some internet, you should shortly have a copy of the entire thing there.

2. Click the green "Code" button on GitHub, and click "Download ZIP". Then, decompress.

We strongly recommend to use Git though, as this will give you easy access to all other branches besides the `main` branch, which is always kinda behind all other developmental branches that include the latests additions, improvements, and fixes.

## Compiling the model

To use the model, you'll need to compile it into a file that can be executed by your machine, using a few software that you might or might not have on there already. Here's the list:

1. A Fortran compiler: the classic and free `gfortran` works well.
2. OpenMPI: the model is decently parallelised, so you'll need this one too (even if you plan on running it serial mode).
3. PETSc: used to lift the heaviest weight in the model: system of partial differential equations.
4. NetCDF: to read and write input data and model output files. You'll need both the basic and Fortran libraries.
5. LAPACK: to solve linear systems of equations.
6. GNU `make`: to deal with interdependencies between different model Fortran modules during compilation.

If you can install all of these libraries through some package manager (e.g. Homebrew) your life will be easier. That way their versions and configurations will be consistent. 99% of compilation errors are due to wrongly set paths to library files. Make sure your terminal actually knows where all those libraries and files are, and that those files are actually there. Don't underestimate the importance of that last line.

Once you have everything, you'll need to modify any of the `Makefile_include_*` files in `src/` (or create your own) so that the paths to each of the libraries declared there (and a couple of other details) match the paths on your machine. Once your custom `Makefile_include_` file is looking good, you'll need to modify the `src/Makefile`, so it points to your custom file.

Now you are ready to run the `compile_all.csh` script in the main folder (or in `templates/`; in this case copy it to the main folder), which basically goes into `src/`, and runs `make` there. `make` will search for the Makefile, which in turn checks your `Makefile_include_*` file, and then reads the paths to all the required libraries there (and a couple of other details). Then, `make` uses the declared Fortran compiler there (if you got `gfortran` and OpenMPI, you should have access to `mpifort`, so declare that one) and compiles the code.

## Running the model

If the compilation was successful, you should get no errors and a lovely executable file in your main model directory called `IMAU_ICE_program`. That's the one you run in order to use the model. In addition, you'll need the `mpiexec` command and any model configuration file (those found in `config-files`). As an example, you would run the model like this:

`mpiexec -n 2 IMAU_ICE_program config-files/config_test`

The `-n` option tells `mpiexec` how many cores (or processes) you want to use in parallel (here 2). `config-files/config_test` is the path to a configuration file you will use, which tells the model all the details about your simulation. In this case, `config_test` will run a short 1000-year schematic experiment. More "realistic" (e.g. Greenland and Antarctica) experiments will require a number of input data files containing observation-based topographic/bathymetric (initial conditions) and forcing data (boundary conditions and time-dependent series), plus other herbs. Examples are Bedmachine and ERA5 data. You can always contact us (Tijn and Jorjo) to get these: c.j.berends@uu.nl & j.a.bernalesconcha@uu.nl , or get in touch with IMAU to know our current whereabouts. This also applies to any other questions you might have; as always, the harder you tried and the more you googled before getting in touch, the happier we'll be and the faster your issue will be solved.

## Checking the output

Model output uses NetCDF files. Python, Matlab, and the rest all have a way to extract data from these files, so check your preferred visualisation software. A quick, simple, and somewhat grumpy way is to use the `ncview` tool.

## Payment methods

Just joking. Good luck have fun.

Last updated: 2023-08-14 by Jorjo Bernales.