rm -f src_dev/liblapack.so
ln -s /usr/lib64/liblapack.so.3.4.2 src_dev/liblapack.so

module load mpi/openmpi-x86_64

cd src_dev
#make clean
make all
cd ..

rm -f IMAU_ICE_program

mv src_dev/IMAU_ICE_program .
