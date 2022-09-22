#!/bin/bash
#Set job requirements
#SBATCH --time=8:00:00
#SBATCH --partition=thin
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
# #SBATCH --mem=50G

# Load modules for MPI, NetCDF, PETSc, and other required libraries
module load 2021
module load eb/4.5.2
module load foss/2021a
module load netCDF-Fortran/4.5.3-gompi-2021a
module load PETSc/3.15.1-foss-2021a
module load imkl/2021.2.0-iompi-2021a
module load OpenMPI/4.1.1-GCC-10.3.0

# Execute the program
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE1_UKESM1-0-LL-ssp585_MARv3.12_Rlow
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE1_UKESM1-0-LL-ssp585_MARv3.12_Rmed
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE1_UKESM1-0-LL-ssp585_MARv3.12_Rhigh
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE2_UKESM1-0-LL-ssp585_MARv3.12_Rlow
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE2_UKESM1-0-LL-ssp585_MARv3.12_Rmed
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE2_UKESM1-0-LL-ssp585_MARv3.12_Rhigh
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE3_UKESM1-0-LL-ssp585_MARv3.12_Rlow
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE3_UKESM1-0-LL-ssp585_MARv3.12_Rmed
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE3_UKESM1-0-LL-ssp585_MARv3.12_Rhigh
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE4_UKESM1-0-LL-ssp585_MARv3.12_Rlow
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE4_UKESM1-0-LL-ssp585_MARv3.12_Rmed
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE4_UKESM1-0-LL-ssp585_MARv3.12_Rhigh
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE6_UKESM1-0-LL-ssp585_MARv3.12_Rlow
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE6_UKESM1-0-LL-ssp585_MARv3.12_Rmed
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE6_UKESM1-0-LL-ssp585_MARv3.12_Rhigh
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE7_UKESM1-0-LL-ssp585_MARv3.12_Rlow
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE7_UKESM1-0-LL-ssp585_MARv3.12_Rmed
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE7_UKESM1-0-LL-ssp585_MARv3.12_Rhigh
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE8_UKESM1-0-LL-ssp585_MARv3.12_Rlow
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE8_UKESM1-0-LL-ssp585_MARv3.12_Rmed
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE8_UKESM1-0-LL-ssp585_MARv3.12_Rhigh
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE5_UKESM1-0-LL-ssp585_MARv3.12_Rlow
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE5_UKESM1-0-LL-ssp585_MARv3.12_Rmed
srun IMAU_ICE_program PROTECT_projections/config-files/config_IMAUICE5_UKESM1-0-LL-ssp585_MARv3.12_Rhigh
