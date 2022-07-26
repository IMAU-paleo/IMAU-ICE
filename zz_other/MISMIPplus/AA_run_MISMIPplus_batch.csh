#! /bin/csh -f

cd ..
./compile_all.csh

#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_tuneGL_5km_Weertman_DIVA
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_tuneGL_5km_Tsai2015_DIVA
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_tuneGL_5km_Schoof2005_DIVA
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_tuneGL_5km_Weertman_hybrid
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_tuneGL_5km_Tsai2015_hybrid
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_tuneGL_5km_Schoof2005_hybrid
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_tuneGL_2km_Weertman_DIVA
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_tuneGL_2km_Tsai2015_DIVA
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_tuneGL_2km_Schoof2005_DIVA
#
#
#
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice0_5km_Weertman_DIVA
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice0_5km_Tsai2015_DIVA
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice0_5km_Schoof2005_DIVA
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice0_5km_Weertman_hybrid
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice0_5km_Tsai2015_hybrid
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice0_5km_Schoof2005_hybrid
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice0_2km_Weertman_DIVA
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice0_2km_Tsai2015_DIVA
#mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice0_2km_Schoof2005_DIVA



mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Weertman_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Tsai2015_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Schoof2005_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Weertman_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Tsai2015_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Schoof2005_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_2km_Weertman_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_2km_Tsai2015_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_2km_Schoof2005_DIVA_FCMP

mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Weertman_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Tsai2015_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Schoof2005_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Weertman_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Tsai2015_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Schoof2005_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_2km_Weertman_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_2km_Tsai2015_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_2km_Schoof2005_DIVA_PMP

mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Weertman_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Tsai2015_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Schoof2005_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Weertman_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Tsai2015_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_5km_Schoof2005_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_2km_Weertman_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_2km_Tsai2015_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1r_2km_Schoof2005_DIVA_NMP



mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Weertman_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Tsai2015_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Schoof2005_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Weertman_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Tsai2015_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Schoof2005_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_2km_Weertman_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_2km_Tsai2015_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_2km_Schoof2005_DIVA_FCMP

mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Weertman_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Tsai2015_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Schoof2005_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Weertman_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Tsai2015_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Schoof2005_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_2km_Weertman_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_2km_Tsai2015_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_2km_Schoof2005_DIVA_PMP

mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Weertman_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Tsai2015_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Schoof2005_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Weertman_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Tsai2015_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_5km_Schoof2005_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_2km_Weertman_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_2km_Tsai2015_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1rr_2km_Schoof2005_DIVA_NMP



mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Weertman_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Tsai2015_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Schoof2005_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Weertman_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Tsai2015_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Schoof2005_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_2km_Weertman_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_2km_Tsai2015_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_2km_Schoof2005_DIVA_FCMP

mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Weertman_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Tsai2015_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Schoof2005_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Weertman_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Tsai2015_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Schoof2005_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_2km_Weertman_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_2km_Tsai2015_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_2km_Schoof2005_DIVA_PMP

mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Weertman_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Tsai2015_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Schoof2005_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Weertman_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Tsai2015_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_5km_Schoof2005_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_2km_Weertman_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_2km_Tsai2015_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice1ra_2km_Schoof2005_DIVA_NMP



mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Weertman_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Tsai2015_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Schoof2005_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Weertman_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Tsai2015_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Schoof2005_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_2km_Weertman_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_2km_Tsai2015_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_2km_Schoof2005_DIVA_FCMP

mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Weertman_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Tsai2015_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Schoof2005_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Weertman_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Tsai2015_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Schoof2005_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_2km_Weertman_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_2km_Tsai2015_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_2km_Schoof2005_DIVA_PMP

mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Weertman_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Tsai2015_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Schoof2005_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Weertman_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Tsai2015_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_5km_Schoof2005_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_2km_Weertman_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_2km_Tsai2015_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2r_2km_Schoof2005_DIVA_NMP



mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Weertman_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Tsai2015_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Schoof2005_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Weertman_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Tsai2015_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Schoof2005_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_2km_Weertman_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_2km_Tsai2015_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_2km_Schoof2005_DIVA_FCMP

mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Weertman_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Tsai2015_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Schoof2005_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Weertman_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Tsai2015_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Schoof2005_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_2km_Weertman_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_2km_Tsai2015_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_2km_Schoof2005_DIVA_PMP

mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Weertman_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Tsai2015_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Schoof2005_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Weertman_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Tsai2015_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_5km_Schoof2005_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_2km_Weertman_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_2km_Tsai2015_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2rr_2km_Schoof2005_DIVA_NMP



mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Weertman_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Tsai2015_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Schoof2005_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Weertman_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Tsai2015_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Schoof2005_hybrid_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_2km_Weertman_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_2km_Tsai2015_DIVA_FCMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_2km_Schoof2005_DIVA_FCMP

mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Weertman_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Tsai2015_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Schoof2005_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Weertman_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Tsai2015_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Schoof2005_hybrid_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_2km_Weertman_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_2km_Tsai2015_DIVA_PMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_2km_Schoof2005_DIVA_PMP

mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Weertman_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Tsai2015_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Schoof2005_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Weertman_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Tsai2015_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_5km_Schoof2005_hybrid_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_2km_Weertman_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_2km_Tsai2015_DIVA_NMP
mpiexec -n 2 IMAU_ICE_program   MISMIPplus/config-files/config_MISMIPplus_ice2ra_2km_Schoof2005_DIVA_NMP