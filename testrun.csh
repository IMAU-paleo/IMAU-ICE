#! /bin/csh -f

./compile_all.csh

rm -rf results_20*

rm -rf results_invert_bed_fixedshelf_40km
mpiexec -n 2 IMAU_ICE_program   config-files/config_invert_bed_fixedshelf_40km

#rm -rf results_invert_BMB_fixedsheet_40km
#mpiexec -n 2 IMAU_ICE_program   config-files/config_invert_BMB_fixedsheet_40km