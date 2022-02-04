#! /bin/csh -f

cd ..
./compile_all.csh

#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_Mplus_FCMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_Mplus_NMP
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_Mplus_PMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_PICOP_FCMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_PICOP_NMP
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_PICOP_PMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_PICO_FCMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_PICO_NMP
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_PICO_PMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_lin_FCMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_lin_NMP
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_lin_PMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_plume_FCMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_plume_NMP
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_plume_PMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_quad_FCMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_quad_NMP
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_2km_quad_PMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_Mplus_FCMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_Mplus_NMP
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_Mplus_PMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_PICOP_FCMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_PICOP_NMP
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_PICOP_PMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_PICO_FCMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_PICO_NMP
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_PICO_PMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_lin_FCMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_lin_NMP
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_lin_PMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_plume_FCMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_plume_NMP
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_plume_PMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_quad_FCMP
#mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_quad_NMP
mpiexec -n 2 IMAU_ICE_program   MISOMIP1/config-files/subgridmelt_template   MISOMIP1/config-files/subgridmelt_var_5km_quad_PMP

