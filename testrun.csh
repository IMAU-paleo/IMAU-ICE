#! /bin/csh -f

./compile_all.csh

rm -rf results_202*
rm -rf BIVMIP/BIVMIP_B_inv_20km_perfect

#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_perfect_40km

#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_inv_40km_perfect
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_inv_20km_perfect
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_inv_10km_perfect

#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_inv_40km_perfect_cont
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_inv_20km_perfect_cont
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_inv_10km_perfect_cont