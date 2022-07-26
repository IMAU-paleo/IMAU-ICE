#! /bin/csh -f

cd ..
./compile_all.csh

# =========================
# ===== Experiment II =====
# =========================

# Tune flow factor to achieve stable grounding-line position at x = 450 km
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_tuneA_5km
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_tuneA_2km

# NOTE: The "tuneA" runs use a little iterative procedure to tune the flow factor A to achieve
#       a steady-state grounding-line position at x = 450 km, following the MISMIP+ protocol.
#       The values of A are written to the screen; be sure to use the last value for all
#       subsequent runs!

# Target runs
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_target_5km
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_target_2km

# NOTE: Before starting the inversion runs, be sure to create the target
#       velocity files using the Matlab script "AA_create_target_velocity_file.m"!

# Unperturbed inversion runs
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_unperturbed
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_2km_unperturbed

# NOTE: Before starting the perturbed topography inversion runs, be sure to create the
#       perturbed geometries by running the Matlab script "AA_create_experiment_II_perturbed_topo.m"
#       (located in input_data/experiment_II)!

# Perturbed inversion runs
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_visc_hi
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_visc_lo
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_SMB_hi
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_SMB_lo
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_BMB_hi
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_BMB_lo
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_topo_hi
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_topo_lo
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_ut_hi
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_ut_lo
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_p_hi
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_p_lo
#mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_inv_5km_noneq

# Retreat runs
rm -rf Berends2022_basal_inversion/exp_II_retreat_5km_*
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_target
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_unperturbed
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_visc_hi
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_visc_lo
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_SMB_hi
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_SMB_lo
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_BMB_hi
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_BMB_lo
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_topo_hi
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_topo_lo
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_ut_hi
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_ut_lo
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_p_hi
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_p_lo
mpiexec -n 2 IMAU_ICE_program   Berends2022_basal_inversion/config-files/config_exp_II_retreat_5km_noneq