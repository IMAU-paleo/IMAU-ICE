#! /bin/csh -f

./compile_all.csh

rm -rf results_ANT_relax_100yr_40km
#rm -rf results_ANT_thermo_LGC_40km

mpiexec -n 2 IMAU_ICE_program src/config-files/config_ANT_relax_100yr_40km
#mpiexec -n 2 IMAU_ICE_program src/config-files/config_ANT_thermo_240kyr_40km