#! /bin/csh -f

./compile_all.csh

rm -rf results_2021*

mpiexec -n 2 IMAU_ICE_program   config-files/config_ANT_thermo_240kyr_40km