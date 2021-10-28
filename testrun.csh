#! /bin/csh -f

./compile_all.csh

rm -rf results_2021*

mpiexec -n 2 IMAU_ICE_program MISOMIPplus/config-files/config_MISOMIPplus_IceOcean0_5km_slid1_Favier2019_lin