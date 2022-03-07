#! /bin/csh -f

rm -rf SSA_icestream_*

cd ..
./compile_all.csh

#rm -rf results_2021*

mpiexec -n 2 IMAU_ICE_program   SSA_icestream/config-files/config_template_SSA_icestream  SSA_icestream/config-files/config_var_40km
mpiexec -n 2 IMAU_ICE_program   SSA_icestream/config-files/config_template_SSA_icestream  SSA_icestream/config-files/config_var_32km
mpiexec -n 2 IMAU_ICE_program   SSA_icestream/config-files/config_template_SSA_icestream  SSA_icestream/config-files/config_var_20km
mpiexec -n 2 IMAU_ICE_program   SSA_icestream/config-files/config_template_SSA_icestream  SSA_icestream/config-files/config_var_16km
mpiexec -n 2 IMAU_ICE_program   SSA_icestream/config-files/config_template_SSA_icestream  SSA_icestream/config-files/config_var_10km