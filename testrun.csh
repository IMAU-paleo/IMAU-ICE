#! /bin/csh -f

./compile_all.csh

rm -rf results_2021*

mpiexec -n 2 IMAU_ICE_program   src/config-files/config_test2 EISMINT1/config-files/var_50km