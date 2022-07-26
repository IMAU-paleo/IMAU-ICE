#! /bin/csh -f

cd ..
./compile_all.csh

rm -rf ABUMIP/ABUMIP_*

#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUC_40km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUC_32km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUC_20km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUC_16km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUC_10km

#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUM_40km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUM_32km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUM_20km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUM_16km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUM_10km

mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUK_40km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUK_32km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUK_20km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUK_16km
#mpiexec -n 2 IMAU_ICE_program   ABUMIP/config-files/config_ABUMIP_ABUCK10km