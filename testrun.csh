#! /bin/csh -f

./compile_all.csh

rm -rf results_202*

mpiexec -n 2 IMAU_ICE_program   config-files/config_ABUM_40km