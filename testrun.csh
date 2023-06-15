#! /bin/csh -f

./compile_all.csh

#rm -rf spinup_GRL/phase5_historicalperiod/hybrid_ZoetIverson_PMIP3ens_40km

#mpiexec -n 2 IMAU_ICE_program spinup_GRL/phase5_historicalperiod/config-files/config_hybrid_ZoetIverson_PMIP3ens_40km

rm -rf PROTECT_projections/IMAUICE1_ctrl
rm -rf PROTECT_projections/IMAUICE2_ctrl
rm -rf PROTECT_projections/IMAUICE3_ctrl
rm -rf PROTECT_projections/IMAUICE4_ctrl
rm -rf PROTECT_projections/IMAUICE5_ctrl
rm -rf PROTECT_projections/IMAUICE6_ctrl
rm -rf PROTECT_projections/IMAUICE7_ctrl
rm -rf PROTECT_projections/IMAUICE8_ctrl

mpiexec -n 2 IMAU_ICE_program PROTECT_projections/config-files/config_control_IMAUICE1
mpiexec -n 2 IMAU_ICE_program PROTECT_projections/config-files/config_control_IMAUICE2
mpiexec -n 2 IMAU_ICE_program PROTECT_projections/config-files/config_control_IMAUICE3
mpiexec -n 2 IMAU_ICE_program PROTECT_projections/config-files/config_control_IMAUICE4
mpiexec -n 2 IMAU_ICE_program PROTECT_projections/config-files/config_control_IMAUICE5
mpiexec -n 2 IMAU_ICE_program PROTECT_projections/config-files/config_control_IMAUICE6
mpiexec -n 2 IMAU_ICE_program PROTECT_projections/config-files/config_control_IMAUICE7
mpiexec -n 2 IMAU_ICE_program PROTECT_projections/config-files/config_control_IMAUICE8