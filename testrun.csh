#! /bin/csh -f

./compile_all.csh

rm -rf spinup_GRL_new/phase2_glacialinception/hybrid_*

mpiexec -n 2 IMAU_ICE_program spinup_GRL_new/phase2_glacialinception/config-files/config_hybrid_ZoetIverson_PMIP3ens_40km