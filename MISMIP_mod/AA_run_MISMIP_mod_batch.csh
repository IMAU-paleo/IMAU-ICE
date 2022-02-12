#! /bin/csh -f

cd ..
./compile_all.csh

# spin-up simulations
#mpiexec -n 2 IMAU_ICE_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA          MISMIP_mod/config-files/config_var_10km
#mpiexec -n 2 IMAU_ICE_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid        MISMIP_mod/config-files/config_var_10km
mpiexec -n 2 IMAU_ICE_program   MISMIP_mod/config-files/config_template_MISMIP_mod_DIVA_sans     MISMIP_mod/config-files/config_var_10km
#mpiexec -n 2 IMAU_ICE_program   MISMIP_mod/config-files/config_template_MISMIP_mod_hybrid_sans   MISMIP_mod/config-files/config_var_10km