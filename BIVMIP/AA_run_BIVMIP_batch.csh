#! /bin/csh -f

cd ..
./compile_all.csh

# ====================
# ===== BIVMIP A =====
# ====================

# Perfect model runs
#rm -rf BIVMIP/BIVMIP_A_perfect*
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_perfect_40km
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_perfect_20km
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_perfect_10km

# Perfect inversions
rm -rf BIVMIP/BIVMIP_A_inv*
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_40km_perfect
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_20km_perfect
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_A_inv_10km_perfect

# ====================
# ===== BIVMIP B =====
# ====================

# Perfect model runs
#rm -rf BIVMIP/BIVMIP_B_perfect*
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_perfect_40km
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_perfect_20km
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_perfect_10km

# Perfect inversions
rm -rf BIVMIP/BIVMIP_B_inv*
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_inv_40km_perfect
mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_inv_20km_perfect
#mpiexec -n 2 IMAU_ICE_program   BIVMIP/config-files/config_BIVMIP_B_inv_10km_perfect

