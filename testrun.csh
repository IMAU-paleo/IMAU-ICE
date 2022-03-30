#! /bin/csh -f

./compile_all.csh

rm -rf results_202*

mpiexec -n 2 IMAU_ICE_program   config-files/config_test    MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO