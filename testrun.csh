#! /bin/csh -f

./compile_all.csh

rm -rf results_202*

mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_20km_perfect
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_10km_perfect