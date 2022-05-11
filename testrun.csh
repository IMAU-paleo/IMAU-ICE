#! /bin/csh -f

./compile_all.csh

rm -rf results_20*

#mpiexec -n 2 IMAU_ICE_program   config-files/config_spinup01_relax_fixedshelf_100yr_40km
#mpiexec -n 2 IMAU_ICE_program   config-files/config_spinup02_thermo_fixedgeo_240kyr_40km
#mpiexec -n 2 IMAU_ICE_program   config-files/config_spinup03_relax_fixedshelf_100kyr_40km
rm -rf results_spinup04_free_schematicocean_10kyr_40km
mpiexec -n 2 IMAU_ICE_program   config-files/config_spinup04_free_schematicocean_10kyr_40km
