#! /bin/csh -f

cd src_v1.1.1

#make clean
make all

cd ..

rm -f IMAU_ICE_program

mv src_v1.1.1/IMAU_ICE_program .
