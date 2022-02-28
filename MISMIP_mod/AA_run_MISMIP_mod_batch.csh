#! /bin/csh -f

cd ..
./compile_all.csh

#mpiexec -n 2 IMAU_ICE_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA          MISMIP_mod/config-files/config_var_40km
mpiexec -n 2 IMAU_ICE_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA          MISMIP_mod/config-files/config_var_32km
mpiexec -n 2 IMAU_ICE_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA          MISMIP_mod/config-files/config_var_20km
mpiexec -n 2 IMAU_ICE_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA          MISMIP_mod/config-files/config_var_16km
mpiexec -n 2 IMAU_ICE_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA          MISMIP_mod/config-files/config_var_10km