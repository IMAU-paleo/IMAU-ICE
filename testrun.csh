#! /bin/csh -f

./compile_all.csh

rm -rf results_20*
rm -rf extrapolated_ocean_files/*

mpiexec -n 2 IMAU_ICE_program   config-files/config_ANT_relax_100yr_40km