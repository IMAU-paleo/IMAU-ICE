#! /bin/csh -f

./compile_all.csh

rm -rf results_2021*

mpiexec -n 2 IMAU_ICE_program src_dev/config-files/config_test
#mpiexec -n 2 IMAU_ICE_program src_dev/config-files/config-files-benchmark/config_MISMIP_mod_64km