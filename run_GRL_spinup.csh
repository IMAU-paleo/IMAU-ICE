#! /bin/csh -f

mpiexec -n 1 IMAU_ICE_program config-files/config_GRL_40km_spin01_calib_csm
mpiexec -n 1 IMAU_ICE_program config-files/config_GRL_40km_spin02_glacial_csm
mpiexec -n 1 IMAU_ICE_program config-files/config_GRL_40km_spin03_termina_csm
mpiexec -n 1 IMAU_ICE_program config-files/config_GRL_40km_spin04_holocene_csm

echo 'Done :)'
echo ''
