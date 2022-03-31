#! /bin/csh -f

cd ..
./compile_all.csh

# ====================
# ===== BIVMIP A =====
# ====================

# Perfect model runs
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_perfect_40km
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_perfect_20km
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_perfect_10km

# Perfect inversions
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_40km_perfect
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_20km_perfect
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_10km_perfect

# Perturbed inversions