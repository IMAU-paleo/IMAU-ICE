module purge

module load 2021                                   # load the 2021 software environment
module load foss/2021a                             # load the foss toolchain
module load netCDF/4.8.0-gompi-2021a
module load netCDF-Fortran/4.5.3-gompi-2021a       # load the netcdf modules built with the intel compilers and MPI
module load PETSc/3.15.1-foss-2021a

cd src

make clean
make all

cd ..

rm -f IMAU_ICE_program

mv src/IMAU_ICE_program .
