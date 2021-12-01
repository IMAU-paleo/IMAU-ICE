#! /bin/csh -f

./compile_all.csh

rm -rf results_2021*

mpiexec -n 2 IMAU_ICE_program   EISMINT1/config-files/EISMINT1_spinup_fixed EISMINT1/config-files/var_50km