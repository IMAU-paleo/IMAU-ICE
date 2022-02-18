#!/bin/sh -f

module purge

module load 2021
module load eb/4.5.2
eblocalinstall PETSc-3.15.1-foss-2021a.eb
module load foss/2021a
module load netCDF-Fortran/4.5.3-gompi-2021a
module load PETSc/3.15.1-foss-2021a

cd src

make clean
make all

cd ..

rm -f IMAU_ICE_program

mv src/IMAU_ICE_program .
