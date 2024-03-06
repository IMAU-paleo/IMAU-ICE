MODULE calving_module

  ! Contains all the routines for calving.

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_reference_geometry
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             is_floating, map_square_to_square_cons_2nd_order_2D, transpose_dp_2D

  USE netcdf_debug_module,             ONLY: save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, save_variable_as_netcdf_int_3D, &
                                             save_variable_as_netcdf_dp_1D,  save_variable_as_netcdf_dp_2D,  save_variable_as_netcdf_dp_3D

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

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_calving_law'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Apply the selected calving law
    IF     (C%choice_calving_law == 'none') THEN
      ! No calving at all
    ELSEIF (C%choice_calving_law == 'threshold_thickness') THEN
      CALL threshold_thickness_calving( grid, ice)
    ELSE
      CALL crash('unknown choice_calving_law"' // TRIM(C%choice_calving_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_calving_law

! == Routines for different calving laws
! ======================================

  SUBROUTINE threshold_thickness_calving( grid, ice)
    ! A simple threshold thickness calving law

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'threshold_thickness_calving'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Apply calving law
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%mask_cf_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 1 .AND. ice%Hi_eff_cf_a( j,i) < C%calving_threshold_thickness) THEN
        ice%Calving( j,i) = ice%Calving( j,i) - ice%Hi_a( j,i) ! Add the ice thickness to the calving rate
        ice%Hi_a( j,i)    = 0._dp          ! Remove ice
      END IF
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE threshold_thickness_calving

END MODULE calving_module
