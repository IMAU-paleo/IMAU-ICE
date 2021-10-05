rm -f src/liblapack.so
ln -s /usr/lib64/liblapack.so.3.4.2 src/liblapack.so

module load mpi/openmpi-x86_64

cd src
make clean
make all
cd ..

rm -f IMAU_ICE_program

mv src/IMAU_ICE_program .
