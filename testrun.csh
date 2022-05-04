#! /bin/csh -f

./compile_all.csh

rm -rf results_20*

mpiexec -n 2 IMAU_ICE_program   config-files/config_test