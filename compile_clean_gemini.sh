rm -f src_v1.1.0/liblapack.so
ln -s /usr/lib64/liblapack.so.3.4.2 src_v1.1.0/liblapack.so

module load mpi/openmpi-x86_64

cd src_v1.1.0
make clean
make all
cd ..

rm -f IMAU_ICE_program

mv src_v1.1.0/IMAU_ICE_program .