#! /bin/csh -f

./compile_all.csh

rm -rf results_GRL_relax_100yr_40km
rm -rf results_GRL_thermo_LGC_40km

mpiexec -n 2 IMAU_ICE_program src_dev/config-files/config_GRL_relax_100yr_40km
mpiexec -n 2 IMAU_ICE_program src_dev/config-files/config_GRL_thermo_240kyr_40km