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
  USE data_types_netcdf_module,        ONLY: type_netcdf_prescribed_retreat_mask, type_netcdf_prescribed_retreat_mask_refice
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, inquire_prescribed_retreat_mask_file, &
                                             read_prescribed_retreat_mask_file, inquire_prescribed_retreat_mask_refice_file, &
                                             read_prescribed_retreat_mask_refice_file
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             is_floating, map_square_to_square_cons_2nd_order_2D, transpose_dp_2D

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
    ELSEIF (C%choice_calving_law == 'threshold_distance') THEN
      CALL threshold_distance_calving( grid, ice)
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
      IF (ice%mask_cf_a( j,i) == 1 .AND. ice%Hi_eff_cf_a( j,i) < C%calving_threshold_thickness) THEN
        ice%Hi_a( j,i) = 0._dp
      END IF
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE threshold_thickness_calving

  SUBROUTINE threshold_distance_calving( grid, ice)
    ! A simple threshold distance-from-center calving law

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'threshold_distance_calving'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: radius

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Apply calving law
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      radius = sqrt( grid%x(i)*grid%x(i) + grid%y(j)*grid%y(j))

      IF (ice%mask_cf_a( j,i) == 1 .AND. &
          ice%mask_shelf_a( j,i) == 1 .AND. &
          radius > C%calving_threshold_distance * 1000._dp) THEN
        ice%Hi_a( j,i) = 0._dp
      END IF
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE threshold_distance_calving

! == Routines for applying a prescribed retreat mask
! ==================================================

  SUBROUTINE apply_prescribed_retreat_mask( grid, ice, time)
    ! Apply a prescribed retreat mask (e.g. as in ISMIP6-Greenland or PROTECT-Greenland)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_prescribed_retreat_mask'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: wt0, wt1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time < ice%ice_fraction_retreat_mask_t0 .OR. time > ice%ice_fraction_retreat_mask_t1) THEN

      ! Find and read the two global time frames
      CALL sync
      CALL update_prescribed_retreat_mask_timeframes( grid, ice, time)

    END IF ! IF (time < ice%ice_fraction_retreat_mask_t0 .OR. time > ice%ice_fraction_retreat_mask_t1) THEN

    ! Interpolate the two timeframes in time
    wt0 = MAX( 0._dp, MIN( 1._dp, (ice%ice_fraction_retreat_mask_t1 - time) / (ice%ice_fraction_retreat_mask_t1 - ice%ice_fraction_retreat_mask_t0) ))
    wt1 = 1._dp - wt0

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ice%ice_fraction_retreat_mask( j,i) = (wt0 * ice%ice_fraction_retreat_mask0( j,i)) + &
                                            (wt1 * ice%ice_fraction_retreat_mask1( j,i))

    END DO
    END DO
    CALL sync

    ! Apply the retreat mask following the ISMIP protocol
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (ice%ice_fraction_retreat_mask( j,i) < 0.999_dp) THEN
        ice%Hi_a( j,i) = MIN( ice%Hi_a( j,i), ice%retreat_mask_Hi_ref( j,i) * ice%ice_fraction_retreat_mask( j,i))
      END IF

      ! Make sure it also works at the margin right from the start
      IF (i > 1 .AND. i < grid%nx .AND. j > 1 .AND. j < grid%ny) THEN
        IF (ice%ice_fraction_retreat_mask( j,i) >= 0.999_dp) THEN
          IF (ice%ice_fraction_retreat_mask( j-1,i-1) < 0.999_dp .OR. &
              ice%ice_fraction_retreat_mask( j-1,i  ) < 0.999_dp .OR. &
              ice%ice_fraction_retreat_mask( j-1,i+1) < 0.999_dp .OR. &
              ice%ice_fraction_retreat_mask( j  ,i-1) < 0.999_dp .OR. &
              ice%ice_fraction_retreat_mask( j  ,i+1) < 0.999_dp .OR. &
              ice%ice_fraction_retreat_mask( j+1,i-1) < 0.999_dp .OR. &
              ice%ice_fraction_retreat_mask( j+1,i  ) < 0.999_dp .OR. &
              ice%ice_fraction_retreat_mask( j+1,i+1) < 0.999_dp) THEN
            ice%Hi_a( j,i) = MIN( ice%Hi_a( j,i), ice%retreat_mask_Hi_ref( j,i) * ice%ice_fraction_retreat_mask( j,i))
          END IF
        END IF
      END IF

    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_prescribed_retreat_mask

  SUBROUTINE update_prescribed_retreat_mask_timeframes( grid, ice, time)
    ! Update the two timeframes of the prescribed retreat mask

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_prescribed_retreat_mask_timeframes'
    INTEGER                                            :: i,j
    TYPE(type_netcdf_prescribed_retreat_mask)          :: netcdf
    TYPE(type_grid)                                    :: grid_raw
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  ice_fraction_retreat_mask_raw0
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  ice_fraction_retreat_mask_raw1
    INTEGER                                            :: wice_fraction_retreat_mask_raw0
    INTEGER                                            :: wice_fraction_retreat_mask_raw1


    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire into the NetCDF file
    CALL allocate_shared_int_0D( grid_raw%nx, grid_raw%wnx)
    CALL allocate_shared_int_0D( grid_raw%ny, grid_raw%wny)
    CALL inquire_prescribed_retreat_mask_file( netcdf, grid_raw%nx, grid_raw%ny)
    CALL sync

    ! Allocate memory, read data
    CALL allocate_shared_dp_1D( grid_raw%nx,              grid_raw%x                    , grid_raw%wx                    )
    CALL allocate_shared_dp_1D(              grid_raw%ny, grid_raw%y                    , grid_raw%wy                    )
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, ice_fraction_retreat_mask_raw0, wice_fraction_retreat_mask_raw0)
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, ice_fraction_retreat_mask_raw1, wice_fraction_retreat_mask_raw1)
    CALL read_prescribed_retreat_mask_file( netcdf, grid_raw%x, grid_raw%y, time, ice_fraction_retreat_mask_raw0, ice_fraction_retreat_mask_raw1, &
      ice%ice_fraction_retreat_mask_t0, ice%ice_fraction_retreat_mask_t1)
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( ice_fraction_retreat_mask_raw0, 'ice_fraction_retreat_mask_raw0')
    CALL check_for_NaN_dp_2D( ice_fraction_retreat_mask_raw1, 'ice_fraction_retreat_mask_raw1')

    ! Since we want data represented as [j,i] internally, transpose the data we just read.
    CALL transpose_dp_2D( ice_fraction_retreat_mask_raw0, wice_fraction_retreat_mask_raw0)
    CALL transpose_dp_2D( ice_fraction_retreat_mask_raw1, wice_fraction_retreat_mask_raw1)

    ! Map (transposed) raw data to the model grid
    CALL map_square_to_square_cons_2nd_order_2D( grid_raw%nx, grid_raw%ny, grid_raw%x, grid_raw%y, grid%nx, grid%ny, grid%x, grid%y, &
      ice_fraction_retreat_mask_raw0, ice%ice_fraction_retreat_mask0)
    CALL map_square_to_square_cons_2nd_order_2D( grid_raw%nx, grid_raw%ny, grid_raw%x, grid_raw%y, grid%nx, grid%ny, grid%x, grid%y, &
      ice_fraction_retreat_mask_raw1, ice%ice_fraction_retreat_mask1)

    ! Limit ice fractions to [0,1]
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%ice_fraction_retreat_mask0( j,i) = MAX( 0._dp, MIN( 1._dp, ice%ice_fraction_retreat_mask0( j,i)))
      ice%ice_fraction_retreat_mask1( j,i) = MAX( 0._dp, MIN( 1._dp, ice%ice_fraction_retreat_mask1( j,i)))
    END DO
    END DO
    CALL sync

    ! Deallocate raw data
    CALL deallocate_shared( grid_raw%wnx                   )
    CALL deallocate_shared( grid_raw%wny                   )
    CALL deallocate_shared( grid_raw%wx                    )
    CALL deallocate_shared( grid_raw%wy                    )
    CALL deallocate_shared( wice_fraction_retreat_mask_raw0)
    CALL deallocate_shared( wice_fraction_retreat_mask_raw1)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_prescribed_retreat_mask_timeframes

  SUBROUTINE initialise_retreat_mask_refice( grid, ice)
    ! Initialise the reference ice thickness for a prescribed retreat mask

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_retreat_mask_refice'
    TYPE(type_netcdf_prescribed_retreat_mask_refice)   :: netcdf
    TYPE(type_grid)                                    :: grid_raw
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  Hi_raw
    INTEGER                                            :: wHi_raw

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire if all the required fields are present in the specified NetCDF file,
    ! and determine the dimensions of the memory to be allocated.
    CALL allocate_shared_int_0D( grid_raw%nx, grid_raw%wnx)
    CALL allocate_shared_int_0D( grid_raw%ny, grid_raw%wny)
    IF (par%master) THEN
      netcdf%filename = C%prescribed_retreat_mask_refice_filename
      CALL inquire_prescribed_retreat_mask_refice_file( netcdf, grid_raw%nx, grid_raw%ny)
    END IF
    CALL sync

    ! Allocate memory for raw data
    CALL allocate_shared_dp_1D( grid_raw%nx,              grid_raw%x , grid_raw%wx )
    CALL allocate_shared_dp_1D(              grid_raw%ny, grid_raw%y , grid_raw%wy )
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, Hi_raw     , wHi_raw     )

    ! Read raw data
    CALL read_prescribed_retreat_mask_refice_file( netcdf, grid_raw%x, grid_raw%y, Hi_raw)
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( Hi_raw, 'Hi_raw')

    ! Since we want data represented as [j,i] internally, transpose the data we just read.
    CALL transpose_dp_2D( Hi_raw, wHi_raw)

    ! Map (transposed) raw data to the model grid
    CALL map_square_to_square_cons_2nd_order_2D( grid_raw%nx, grid_raw%ny, grid_raw%x, &
      grid_raw%y, grid%nx, grid%ny, grid%x, grid%y, Hi_raw, ice%retreat_mask_Hi_ref)

    ! Deallocate raw data
    CALL deallocate_shared( grid_raw%wnx)
    CALL deallocate_shared( grid_raw%wny)
    CALL deallocate_shared( grid_raw%wx )
    CALL deallocate_shared( grid_raw%wy )
    CALL deallocate_shared( wHi_raw     )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_retreat_mask_refice

END MODULE calving_module
