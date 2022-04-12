#! /bin/csh -f

cd ..
./compile_all.csh

# ========================
# ===== Experiment I =====
# ========================

# Target runs
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_target_40km
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_target_20km
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_target_10km

# NOTE: Before starting the inversion runs, be sure to create the target
#       velocity files using the Matlab script "AA_create_target_velocity_file.m"!

# Unperturbed inversion runs
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_40km_unperturbed
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_20km_unperturbed
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_10km_unperturbed

# NOTE: Before starting the perturbed topography inversion runs, be sure to create the
#       perturbed geometries by running the Matlab script "AA_create_experiment_I_perturbed_topo.m"
#       (located in input_data/experiment_I)!

# Perturbed inversion runs
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_20km_visc_hi
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_20km_visc_lo
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_20km_SMB_hi
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_20km_SMB_lo
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_20km_topo_hi
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_20km_topo_lo
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_20km_ut_hi
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_20km_ut_lo
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_20km_p_hi
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_inv_20km_p_lo

# Retreat runs
rm -rf Berends2022_basal_inversion/exp_I_retreat_20km_*
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_retreat_20km_target
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_retreat_20km_unperturbed
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_retreat_20km_visc_hi
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_retreat_20km_visc_lo
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_retreat_20km_SMB_hi
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_retreat_20km_SMB_lo
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_retreat_20km_topo_hi
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_retreat_20km_topo_lo
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_retreat_20km_ut_hi
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_retreat_20km_ut_lo
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_retreat_20km_p_hi
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_I_retreat_20km_p_lo