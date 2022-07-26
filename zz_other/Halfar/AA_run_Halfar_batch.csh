#! /bin/csh -f

cd ..
./compile_all.csh

#rm -rf results_2021*
#rm -rf Halfar_*

#mpiexec -n 2 IMAU_ICE_program   Halfar/config-files/config_template_Halfar  Halfar/config-files/config_var_50km
mpiexec -n 2 IMAU_ICE_program   Halfar/config-files/config_template_Halfar  Halfar/config-files/config_var_25km
mpiexec -n 2 IMAU_ICE_program   Halfar/config-files/config_template_Halfar  Halfar/config-files/config_var_16km
mpiexec -n 2 IMAU_ICE_program   Halfar/config-files/config_template_Halfar  Halfar/config-files/config_var_8km