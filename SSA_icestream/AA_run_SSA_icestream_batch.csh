#! /bin/csh -f

cd ..
./compile_all.csh

#rm -rf results_2021*
rm -rf SSA_icestream_*

mpiexec -n 2 IMAU_ICE_program   SSA_icestream/config-files/config_template_SSA_icestream  SSA_icestream/config-files/config_var_5km
mpiexec -n 2 IMAU_ICE_program   SSA_icestream/config-files/config_template_SSA_icestream  SSA_icestream/config-files/config_var_4km
mpiexec -n 2 IMAU_ICE_program   SSA_icestream/config-files/config_template_SSA_icestream  SSA_icestream/config-files/config_var_3km
#mpiexec -n 2 IMAU_ICE_program   SSA_icestream/config-files/config_template_SSA_icestream  SSA_icestream/config-files/config_var_2km
#mpiexec -n 2 IMAU_ICE_program   SSA_icestream/config-files/config_template_SSA_icestream  SSA_icestream/config-files/config_var_1km
#mpiexec -n 2 IMAU_ICE_program   SSA_icestream/config-files/config_template_SSA_icestream  SSA_icestream/config-files/config_var_500m
#mpiexec -n 2 IMAU_ICE_program   SSA_icestream/config-files/config_template_SSA_icestream  SSA_icestream/config-files/config_var_250m