#! /bin/csh -f

./compile_all.csh

rm -rf results_2021*
rm -rf MISMIP_mod_spinup_40km

mpiexec -n 2 IMAU_ICE_program   config-files/config_MISMIP_mod_spinup_40km