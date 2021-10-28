#! /bin/csh -f

./compile_all.csh

rm -rf results_2021*
rm -rf MISOMIPplus/MISOMIPplus_IceOcean0_2km_slid1_PICO

mpiexec -n 2 IMAU_ICE_program   MISOMIPplus/config-files/template_MISOMIPplus_IceOcean0   MISOMIPplus/config-files/var_2km_slid1_PICO