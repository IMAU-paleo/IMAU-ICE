#!/bin/bash

# This is the script you run from the terminal to start an IMAU-ICE simulation.
# Doesn't need any arguments, those are listed below.
# The number of cores and requested walltime are set in run_IMAU_ICE_Cartesius.csh

# Clean up the scratch directory
rm -rf /scratch-local/berends/*

# Load relevant modules
module purge
module load 2020                                  # load the 2020 software environment
module load intel/2020a                           # load the intel toolchain (intel compilers, intel MPI, MKL, etc)
module load netCDF/4.7.4-iimpi-2020a              # load the netcdf modules built with the intel compilers and MPI
module load netCDF-Fortran/4.5.2-iimpi-2020a      # load the netcdf modules built with the intel compilers and MPI
module load PETSc/3.12.4-intel-2020a-Python-3.8.2
module load impi/2019.8.254-iccifort-2020.1.217

# Submit the run script with sbatch
sbatch run_IMAU_ICE_Cartesius.csh src_dev/config-files/config_test test
