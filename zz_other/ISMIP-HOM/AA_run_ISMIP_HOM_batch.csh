#! /bin/csh -f

cd ..
./compile_all.csh

#rm -rf results_2021*
rm -rf ISMIP_HOM_*

#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_160_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_080_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_040_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_020_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_010_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_005_DIVA

#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_160_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_080_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_040_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_020_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_010_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_005_hybrid

#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_160_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_080_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_040_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_020_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_010_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_005_DIVA

#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_160_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_080_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_040_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_020_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_010_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_005_hybrid

#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_160_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_080_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_040_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_020_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_010_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_005_DIVA

#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_160_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_080_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_040_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_020_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_010_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_005_hybrid

#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_160_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_080_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_040_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_020_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_010_DIVA
mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_005_DIVA

#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_160_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_080_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_040_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_020_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_010_hybrid
mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_005_hybrid

#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_E_0_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_E_0_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_E_1_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_E_1_hybrid

#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_F_0_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_F_0_hybrid
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_F_1_DIVA
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_F_1_hybrid

#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_160_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_080_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_040_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_020_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_010_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_005_DIVA_sans
#
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_160_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_080_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_040_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_020_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_010_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_A  ISMIP-HOM/config-files/config_var_005_hybrid_sans
#
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_160_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_080_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_040_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_020_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_010_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_005_DIVA_sans
#
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_160_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_080_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_040_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_020_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_010_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_B  ISMIP-HOM/config-files/config_var_005_hybrid_sans
#
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_160_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_080_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_040_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_020_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_010_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_005_DIVA_sans
#
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_160_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_080_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_040_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_020_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_010_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_C  ISMIP-HOM/config-files/config_var_005_hybrid_sans
#
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_160_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_080_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_040_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_020_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_010_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_005_DIVA_sans
#
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_160_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_080_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_040_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_020_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_010_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_D  ISMIP-HOM/config-files/config_var_005_hybrid_sans
#
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_F_0_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_F_0_hybrid_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_F_1_DIVA_sans
#mpiexec -n 2 IMAU_ICE_program   ISMIP-HOM/config-files/config_template_ISMIP_HOM_F_1_hybrid_sans