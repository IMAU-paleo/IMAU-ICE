#! /bin/csh -f

./compile_all.csh

rm -rf results_20*
rm -rf results_test_GRL

mpiexec -n 2 IMAU_ICE_program   config-files/config_test