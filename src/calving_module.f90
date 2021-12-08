MODULE calving_module

  ! Contains all the routines for calving.

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_reference_geometry
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             is_floating

  IMPLICIT NONE
  
CONTAINS

! == The main routines that should be called from the main ice model/program
! ==========================================================================

  SUBROUTINE apply_calving_law( grid, ice)
    ! Apply the selected calving law

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Apply the selected calving law
    IF     (C%choice_calving_law == 'none') THEN
      ! No calving at all
    ELSEIF (C%choice_calving_law == 'threshold_thickness') THEN
      CALL threshold_thickness_calving( grid, ice)
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_calving_law "', TRIM(C%choice_calving_law), '" not implemented in calculate_calving_flux!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE apply_calving_law
  
! == Routines for different calving laws
! ===================================

  SUBROUTINE threshold_thickness_calving( grid, ice)
    ! A simple threshold thickness calving law

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    ! Apply calving law
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%mask_cf_a( j,i) == 1 .AND. ice%Hi_eff_cf_a( j,i) < C%calving_threshold_thickness) THEN
        ice%Hi_a( j,i) = 0._dp
      END IF
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE threshold_thickness_calving
  
END MODULE calving_module
