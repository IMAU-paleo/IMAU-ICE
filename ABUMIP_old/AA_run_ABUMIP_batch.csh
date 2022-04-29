#! /bin/csh -f

cd ..
./compile_all.csh

#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_template_ABUMIP_spinup   ABUMIP/config-files/config_var_40km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_template_ABUC            ABUMIP/config-files/config_var_40km
mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_template_ABUM            ABUMIP/config-files/config_var_40km
mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_template_ABUK            ABUMIP/config-files/config_var_40km