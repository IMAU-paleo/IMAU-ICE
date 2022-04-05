#! /bin/csh -f

cd ..
./compile_all.csh

# ========================
# ===== Experiment I =====
# ========================

# Target runs
rm -rf Berends2022_basal_inversion/exp_I_target_*
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_target_40km

