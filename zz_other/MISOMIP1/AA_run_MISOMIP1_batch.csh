#! /bin/csh -f

cd ..
./compile_all.csh

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean0_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP



mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1rr_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP



mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean1ra_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP



mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2rr_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP



mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_5km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Weertman     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Tsai2015     MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP

mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_lin
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_quad
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_Mplus
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_plume
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICO
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/MISOMIP1_template_IceOcean2ra_2km_Schoof2005   MISOMIP1/config-files/MISOMIP1_var_FCMP_PICOP