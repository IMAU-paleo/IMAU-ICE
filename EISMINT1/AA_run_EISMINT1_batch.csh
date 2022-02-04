#! /bin/csh -f

cd ..
./compile_all.csh

#rm -rf results_2021*
#rm -rf EISMINT1_*

# spin-up simulations
mpiexec -n 2 IMAU_ICE_program   EISMINT1/config-files/config_template_EISMINT1_A  EISMINT1/config-files/config_var_50km
mpiexec -n 2 IMAU_ICE_program   EISMINT1/config-files/config_template_EISMINT1_D  EISMINT1/config-files/config_var_50km

# glacial cycle simulations
mpiexec -n 2 IMAU_ICE_program   EISMINT1/config-files/config_template_EISMINT1_B  EISMINT1/config-files/config_var_50km
mpiexec -n 2 IMAU_ICE_program   EISMINT1/config-files/config_template_EISMINT1_C  EISMINT1/config-files/config_var_50km
mpiexec -n 2 IMAU_ICE_program   EISMINT1/config-files/config_template_EISMINT1_E  EISMINT1/config-files/config_var_50km
mpiexec -n 2 IMAU_ICE_program   EISMINT1/config-files/config_template_EISMINT1_F  EISMINT1/config-files/config_var_50km