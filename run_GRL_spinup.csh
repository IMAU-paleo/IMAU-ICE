#! /bin/csh -f

mpiexec -n 1 IMAU_ICE_program config-files/config_GRL_40km_spin01_calib
mpiexec -n 1 IMAU_ICE_program config-files/config_GRL_40km_spin02_glacial
mpiexec -n 1 IMAU_ICE_program config-files/config_GRL_40km_spin03_termina
mpiexec -n 1 IMAU_ICE_program config-files/config_GRL_40km_spin04_holocene

echo 'Done :)'
echo ''
