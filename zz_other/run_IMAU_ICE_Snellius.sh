#!/bin/bash
#Set job requirements
#SBATCH --time=10:00
#SBATCH --partition=thin
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1

# Load modules for MPI and other parallel libraries
module load 2021
module load eb/4.5.2
eblocalinstall PETSc-3.15.1-foss-2021a.eb
module load foss/2021a
module load netCDF-Fortran/4.5.3-gompi-2021a
module load PETSc/3.15.1-foss-2021a

# Execute the program
srun IMAU_ICE_program config-files/config_test
