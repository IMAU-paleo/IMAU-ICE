#! /bin/csh -f

cd src_dev

make clean
make all

cd ..

rm -f IMAU_ICE_program

mv src_dev/IMAU_ICE_program .
