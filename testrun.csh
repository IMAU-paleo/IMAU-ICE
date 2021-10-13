#! /bin/csh -f

./compile_all.csh

rm -rf results_2021*

mpiexec -n 2 IMAU_ICE_program src/config-files/config_MISMIPplus_init_5km_slid3