#! /bin/csh -f

cd src

#make clean
make all

cd ..

rm -f IMAU_ICE_program

mv src/IMAU_ICE_program .
