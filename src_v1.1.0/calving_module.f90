MODULE calving_module

  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parallel_module,                 ONLY: par, sync
  USE data_types_module,               ONLY: type_grid, type_ice_model

  IMPLICIT NONE
    
CONTAINS

  ! Calving is called from the calculate_ice_thickness_change subroutine in the ice_dynamics module.

  ! Simple threshold thickness calving
  SUBROUTINE threshold_thickness_calving( grid, ice)
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(INOUT) :: ice 
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%mask_shelf_Aa( j,i) == 1 .AND. ice%Hi_Aa( j,i) < C%calving_threshold_thickness) ice%Hi_Aa( j,i) = 0._dp
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE threshold_thickness_calving
  
END MODULE calving_module
