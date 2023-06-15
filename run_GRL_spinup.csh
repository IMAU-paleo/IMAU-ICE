#! /bin/csh -f

rm -rf results_20220728_001

mpiexec -n 2 IMAU_ICE_program spinup_GRL/phase1_calib/config-files/config_spinup_GRL_phase1_calib_hybrid_ZoetIverson_PMIP3ens_32km_invTijn
