#!/bin/bash

# This is the script you run from the terminal to start an IMAU-ICE simulation.
# Doesn't need any arguments, those are listed below.
# The number of cores and requested walltime are set in run_IMAU_ICE_Snellius.csh

# Load relevant modules
module purge

module load 2021                                   # load the 2021 software environment
module load foss/2021a                             # load the foss toolchain
module load netCDF/4.8.0-gompi-2021a
module load netCDF-Fortran/4.5.3-gompi-2021a       # load the netcdf modules built with the intel compilers and MPI
module load PETSc/3.15.1-foss-2021a

# Submit the run script with sbatch
sbatch run_IMAU_ICE_Snellius.csh config-files/config_ANT_relax_100yr_40km test
