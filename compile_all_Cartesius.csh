module purge

module load 2020                                   # load the 2020 software environment
module load intel/2020a                            # load the intel toolchain (intel compilers, intel MPI, MKL, etc)
module load netCDF/4.7.4-iimpi-2020a               # load the netcdf modules built with the intel compilers and MPI
module load netCDF-Fortran/4.5.2-iimpi-2020a       # load the netcdf modules built with the intel compilers and MPI
module load PETSc/3.12.4-intel-2020a-Python-3.8.2

cd src_dev

#make clean
make all

cd ..

rm -f IMAU_ICE_program

mv src_dev/IMAU_ICE_program .
