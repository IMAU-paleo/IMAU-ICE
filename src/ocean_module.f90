MODULE ocean_module

  ! Contains all the routines for calculating the climate forcing.

  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_climate_matrix, type_subclimate_global, &
                                             type_climate_model, type_subclimate_region
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             error_function, smooth_Gaussian_2D, smooth_Shepard_2D, &
                                             map_glob_to_grid_2D, map_glob_to_grid_3D, &
                                             map_square_to_square_cons_2nd_order_2D, map_square_to_square_cons_2nd_order_3D

  IMPLICIT NONE
    
CONTAINS

  
  ! Initialising the region-specific ocean data
  SUBROUTINE initialise_oceans_regional( grid, climate, matrix)
    ! Allocate shared memory for the ocean data of the regional subclimates
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    TYPE(type_climate_matrix),           INTENT(IN)    :: matrix
    
    ! Local variables
    REAL(dp), DIMENSION(:    ), POINTER                :: z_ocean_dummy  ! DENK DROM
    INTEGER,                    POINTER                :: nz_ocean_dummy ! DENK DROM
    INTEGER :: wz_ocean_dummy, wnz_ocean_dummy                           ! DENK DROM
    INTEGER                                            :: k              ! DENK DROM
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
        ! Nothing to be done here
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_oceans_regional!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Create dummy vertical dimension - DENK DROM
    ! This should be replaced by copying the vertical coordinate of the World Ocean Atlas data!
    ! GCM data should be mapped from the GCM vertical grid to the WOA vertical grid on the global lat/lon grid immediately after reading the data.
    CALL allocate_shared_int_0D( nz_ocean_dummy, wnz_ocean_dummy)
    IF (par%master) nz_ocean_dummy = 100
    CALL sync
    CALL allocate_shared_dp_1D( nz_ocean_dummy, z_ocean_dummy, wz_ocean_dummy)
    IF (par%master) THEN
      ! Create a "dummy" vertical coordinate running from z=1.25 to z=5500 (similar to the world ocean atlas)
      DO k = 1, nz_ocean_dummy
        z_ocean_dummy( k) = 1.25_dp + (5500._dp - 1.25_dp) * REAL(k-1,dp) / REAL(nz_ocean_dummy-1,dp)
      END DO
    END IF
    CALL sync
    
    ! Allocate memory for ocean data in all the regional subclimates
    CALL allocate_subclimate_regional_oceans( grid, climate%PD_obs,   z_ocean_dummy, nz_ocean_dummy)
    CALL allocate_subclimate_regional_oceans( grid, climate%GCM_PI,   z_ocean_dummy, nz_ocean_dummy)
    CALL allocate_subclimate_regional_oceans( grid, climate%GCM_cold, z_ocean_dummy, nz_ocean_dummy)
    CALL allocate_subclimate_regional_oceans( grid, climate%GCM_warm, z_ocean_dummy, nz_ocean_dummy)
    CALL allocate_subclimate_regional_oceans( grid, climate%applied,  z_ocean_dummy, nz_ocean_dummy)
    
    ! Map ocean data from the global subclimates to the regional subclimates
    ! DENK DROM - these routines don't do anything yet!
    CALL map_ocean_data_global_to_regional( grid, matrix%PD_obs,   climate%PD_obs  )
    CALL map_ocean_data_global_to_regional( grid, matrix%GCM_PI,   climate%GCM_PI  )
    CALL map_ocean_data_global_to_regional( grid, matrix%GCM_cold, climate%GCM_cold)
    CALL map_ocean_data_global_to_regional( grid, matrix%GCM_warm, climate%GCM_warm)
    
    ! Correct regional ocean data for GCM bias
    ! (must be done on regional grid, since the GCM grid and the World Ocean Atlas grid are generally not the same!)
    ! DENK DROM - these routines don't do anything yet!
    CALL correct_GCM_bias_ocean( grid, climate, climate%GCM_PI)
    CALL correct_GCM_bias_ocean( grid, climate, climate%GCM_warm)
    CALL correct_GCM_bias_ocean( grid, climate, climate%GCM_cold)
    
    ! Initialise dummy ocean data fields - DENK DROM
    CALL initialise_dummy_ocean( grid, climate%applied)
  
  END SUBROUTINE initialise_oceans_regional
  SUBROUTINE allocate_subclimate_regional_oceans( grid, climate, z_ocean, nz_ocean)
    ! Allocate shared memory for the ocean data of a regional subclimate models
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: z_ocean
    INTEGER,                             INTENT(IN)    :: nz_ocean
    
    ! Set up the vertical coordinate
    CALL allocate_shared_int_0D( climate%nz_ocean, climate%wnz_ocean)
    IF (par%master) climate%nz_ocean = nz_ocean
    CALL sync
    CALL allocate_shared_dp_1D( climate%nz_ocean, climate%z_ocean, climate%wz_ocean)
    IF (par%master) climate%z_ocean = z_ocean
    CALL sync
    
    ! Allocate shared memory
    CALL allocate_shared_int_3D( climate%nz_ocean, grid%ny, grid%nx, climate%mask_ocean,       climate%wmask_ocean      )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%T_ocean,          climate%wT_ocean         )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%S_ocean,          climate%wS_ocean         )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%T_ocean_corr,     climate%wT_ocean_corr    )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%S_ocean_corr,     climate%wS_ocean_corr    )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%T_ocean_corr_ext, climate%wT_ocean_corr_ext)
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%S_ocean_corr_ext, climate%wS_ocean_corr_ext)
  
  END SUBROUTINE allocate_subclimate_regional_oceans
  SUBROUTINE map_ocean_data_global_to_regional( grid, clim_glob, clim_reg)
    ! Map ocean data for a single subclimate from the global grid to the regional grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_subclimate_global),        INTENT(IN)    :: clim_glob
    TYPE(type_subclimate_region),        INTENT(INOUT) :: clim_reg
    
    ! Nothing done here right now!
  
  END SUBROUTINE map_ocean_data_global_to_regional
  SUBROUTINE correct_GCM_bias_ocean( grid, climate, subclimate)
    ! Correct regional ocean data for GCM bias
    ! (must be done on regional grid, since the GCM grid and the World Ocean Atlas grid are generally not the same!)
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_model),            INTENT(IN)    :: climate
    TYPE(type_subclimate_region),        INTENT(INOUT) :: subclimate
    
    ! Local variables:
    INTEGER                                            :: i,j,k
    
    ! Nothing done here right now!
    
    ! Will probably need think a bit about differences in ice geometry between snapshots and PD?
  
  END SUBROUTINE correct_GCM_bias_ocean
  SUBROUTINE initialise_dummy_ocean( grid, climate)
    ! Initialise dummy ocean data fields - DENK DROM
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    
    ! Local variables
    INTEGER                                            :: i,j,k
    
    ! A simple temperature profile based on Antarctic waters from the World Ocean Atlas
  
  END SUBROUTINE initialise_dummy_ocean

END MODULE ocean_module