#! /bin/csh -f

cd ..
./compile_all.csh

#rm -rf results_2021*
#rm -rf EISMINT1_*

## analytical-solution runs
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_perfect_40km
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_perfect_20km
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_perfect_10km

# perfect inversion runs
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_40km_perfect
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_20km_perfect
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_10km_perfect

# perturbed inversion runs
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_10km_visc_hi
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_10km_visc_lo
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_10km_SMB_hi
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_10km_SMB_lo
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_10km_p_hi
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_10km_p_lo
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_10km_ut_hi
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_10km_ut_lo