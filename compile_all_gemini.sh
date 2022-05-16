module load mpi/openmpi-x86_64
module load petsc/3.16.3

cd src
mkdir module-files
mkdir object-files
#make clean
make all
cd ..

rm -f IMAU_ICE_program

mv src/IMAU_ICE_program .
