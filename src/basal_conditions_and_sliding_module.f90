MODULE basal_conditions_and_sliding_module

  ! Contains all the routines for calculating the basal conditions underneath the ice.

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_reference_geometry, type_BIV_target_velocity
  USE netcdf_output_module,            ONLY: create_BIV_bed_roughness_file
  USE netcdf_input_module,             ONLY: read_field_from_file_2D, read_field_from_xy_file_2D
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             SSA_Schoof2006_analytical_solution, smooth_Gaussian_2D, &
                                             extrapolate_Gaussian_floodfill, transpose_dp_2D, &
                                             map_square_to_square_cons_2nd_order_2D, interp_bilin_2D, &
                                             deallocate_grid
  USE general_ice_model_data_module,   ONLY: surface_elevation

  USE netcdf_debug_module,             ONLY: save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, save_variable_as_netcdf_int_3D, &
                                             save_variable_as_netcdf_dp_1D,  save_variable_as_netcdf_dp_2D,  save_variable_as_netcdf_dp_3D

  IMPLICIT NONE

CONTAINS

  ! The main routine, to be called from the ice_velocity_module
  SUBROUTINE calc_basal_conditions( grid, ice)
    ! Determine the basal conditions underneath the ice

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_conditions'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Basal hydrology
    CALL calc_basal_hydrology( grid, ice)

    ! Bed roughness
    CALL calc_bed_roughness( grid, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_conditions
  SUBROUTINE initialise_basal_conditions( grid, ice, region_name)
    ! Allocation and initialisation

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_conditions'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Basal hydrology
    CALL initialise_basal_hydrology( grid, ice)

    ! Bed roughness
    CALL initialise_bed_roughness( grid, ice, region_name)

    ! Basal inversion
    IF (C%do_BIVgeo) THEN
      CALL initialise_basal_inversion( grid, ice, region_name)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_conditions

! == Basal hydrology
! ==================

  SUBROUTINE calc_basal_hydrology( grid, ice)
    ! Calculate the pore water pressure and effective basal pressure

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_hydrology'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate pore water pressure using the chosen basal hydrology model
    ! ====================================================================

    IF     (C%choice_basal_hydrology == 'dry') THEN
      ! Assume zero subglacial water pressure, i.e. effective pressure is equal to overburden pressure everywhere
      ice%pore_water_pressure_a( :,grid%i1:grid%i2) = 0._dp
    ELSEIF (C%choice_basal_hydrology == 'saturated') THEN
      ! Assume all marine till is saturated (i.e. pore water pressure is equal to water pressure at depth everywhere)
      CALL calc_pore_water_pressure_saturated( grid, ice)
    ELSEIF (C%choice_basal_hydrology == 'Martin2011') THEN
      ! The Martin et al. (2011) parameterisation of pore water pressure
      CALL calc_pore_water_pressure_Martin2011( grid, ice)
    ELSE
      CALL crash('unknown choice_basal_hydrology' // TRIM(C%choice_basal_hydrology))
    END IF

    ! Calculate overburden and effective pressure
    ! ===========================================

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%overburden_pressure_a( j,i) = ice_density * grav * ice%Hi_a( j,i)
      ice%Neff_a(                j,i) = MAX( 0._dp, ice%overburden_pressure_a( j,i) - ice%pore_water_pressure_a( j,i))
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_hydrology
  SUBROUTINE initialise_basal_hydrology( grid, ice)
    ! Allocation and initialisation

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_hydrology'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    IF     (C%choice_basal_hydrology == 'dry') THEN
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%pore_water_pressure_a, ice%wpore_water_pressure_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%overburden_pressure_a, ice%woverburden_pressure_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%Neff_a               , ice%wNeff_a               )
    ELSEIF (C%choice_basal_hydrology == 'saturated') THEN
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%pore_water_pressure_a, ice%wpore_water_pressure_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%overburden_pressure_a, ice%woverburden_pressure_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%Neff_a               , ice%wNeff_a               )
    ELSEIF (C%choice_basal_hydrology == 'Martin2011') THEN
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%pore_water_pressure_a, ice%wpore_water_pressure_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%overburden_pressure_a, ice%woverburden_pressure_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%Neff_a               , ice%wNeff_a               )
    ELSE
      CALL crash('unknown choice_basal_hydrology' // TRIM(C%choice_basal_hydrology))
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_hydrology

  SUBROUTINE calc_pore_water_pressure_saturated( grid, ice)
    ! Calculate the pore water pressure
    !
    ! Assume all till is saturated, i.e. pore water pressure = -rho_w * g * Hb

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pore_water_pressure_saturated'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%pore_water_pressure_a( j,i) = -seawater_density * grav * (ice%SL_a( j,i) - ice%Hb_a( j,i))
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pore_water_pressure_saturated
  SUBROUTINE calc_pore_water_pressure_Martin2011( grid, ice)
    ! Calculate the pore water pressure
    !
    ! Use the parameterisation from Martin et al. (2011)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pore_water_pressure_Martin2011'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: lambda_p

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Pore water pressure scaling factor (Martin et al., 2011, Eq. 12)
      lambda_p = MIN( 1._dp, MAX( 0._dp, 1._dp - (ice%Hb_a( j,i) - ice%SL_a( j,i) - C%Martin2011_hydro_Hb_min) / (C%Martin2011_hydro_Hb_max - C%Martin2011_hydro_Hb_min) ))

      ! Pore water pressure (Martin et al., 2011, Eq. 11)
      ice%pore_water_pressure_a( j,i) = 0.96_dp * ice_density * grav * ice%Hi_a( j,i) * lambda_p

    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pore_water_pressure_Martin2011

! == Bed roughness
! ================

  SUBROUTINE calc_bed_roughness( grid, ice)
    ! Calculate the bed roughness

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! In case of no sliding or "idealised" sliding (e.g. ISMIP-HOM experiments), no bed roughness is required
    IF (C%choice_sliding_law == 'no_sliding' .OR. &
        C%choice_sliding_law == 'idealised') THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    IF (C%choice_basal_roughness == 'uniform') THEN
      ! Apply a uniform bed roughness - no need to do anything, as this has already been set in the initialise routine

    ELSEIF (C%choice_basal_roughness == 'parameterised') THEN
      ! Apply the chosen parameterisation of bed roughness

      IF     (C%choice_param_basal_roughness == 'Martin2011') THEN
        ! The Martin et al. (2011) parameterisation of basal roughness (specifically the till friction angle and till yield stress)
        CALL calc_bed_roughness_Martin2011( grid, ice)
      ELSEIF (C%choice_param_basal_roughness == 'SSA_icestream') THEN
        ! The basal roughness parameterisation in the SSA_icestream idealised-geometry experiment
        CALL calc_bed_roughness_SSA_icestream( grid, ice)
      ELSEIF (C%choice_param_basal_roughness == 'MISMIP+') THEN
        ! The basal roughness parameterisation in the MISMIP+ idealised-geometry experiment
        CALL calc_bed_roughness_MISMIPplus( grid, ice)
      ELSE
        CALL crash('unknown choice_param_basal_roughness' // TRIM(C%choice_param_basal_roughness))
      END IF

    ELSEIF (C%choice_basal_roughness == 'prescribed') THEN
      ! Basal roughness has been initialised from an external file; no need to do anything

    ELSE
      CALL crash('unknown choice_basal_roughness' // TRIM(C%choice_basal_roughness))
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness
  SUBROUTINE initialise_bed_roughness( grid, ice, region_name)

    ! Allocation and initialisation
    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! In case of no sliding or "idealised" sliding (e.g. ISMIP-HOM experiments), no bed roughness is required
    IF (C%choice_sliding_law == 'no_sliding' .OR. &
        C%choice_sliding_law == 'idealised') THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Allocate shared memory
    IF     (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%beta_sq_a , ice%wbeta_sq_a )
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Regularised Coulomb-type sliding law
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%phi_fric_a, ice%wphi_fric_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%tauc_a    , ice%wtauc_a    )
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%alpha_sq_a, ice%walpha_sq_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%beta_sq_a , ice%wbeta_sq_a )
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%alpha_sq_a, ice%walpha_sq_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%beta_sq_a , ice%wbeta_sq_a )
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%phi_fric_a, ice%wphi_fric_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%tauc_a    , ice%wtauc_a    )
    ELSE
      CALL crash('unknown choice_sliding_law' // TRIM(C%choice_sliding_law))
    END IF

    ! Initialise values
    IF (C%choice_basal_roughness == 'uniform') THEN
      ! Apply a uniform bed roughness

      IF     (C%choice_sliding_law == 'Weertman') THEN
        ice%beta_sq_a(   :,grid%i1:grid%i2) = C%uniform_Weertman_beta_sq
      ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
              C%choice_sliding_law == 'Coulomb_regularised') THEN
        ice%phi_fric_a(  :,grid%i1:grid%i2) = C%uniform_Coulomb_phi_fric
      ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
        ice%alpha_sq_a(  :,grid%i1:grid%i2) = C%uniform_Tsai2015_alpha_sq
        ice%beta_sq_a(   :,grid%i1:grid%i2) = C%uniform_Tsai2015_beta_sq
      ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
        ice%alpha_sq_a(  :,grid%i1:grid%i2) = C%uniform_Schoof2005_alpha_sq
        ice%beta_sq_a(   :,grid%i1:grid%i2) = C%uniform_Schoof2005_beta_sq
      ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
        ice%phi_fric_a(  :,grid%i1:grid%i2) = C%uniform_Coulomb_phi_fric
      ELSE
        CALL crash('unknown choice_sliding_law' // TRIM(C%choice_sliding_law))
      END IF

    ELSEIF (C%choice_basal_roughness == 'parameterised') THEN
      ! Apply the chosen parameterisation of bed roughness

      CALL calc_bed_roughness( grid, ice)

    ELSEIF (C%choice_basal_roughness == 'prescribed') THEN
      ! Initialise bed roughness from an external file

      CALL initialise_bed_roughness_from_file( grid, ice, region_name)

    ELSE
      CALL crash('unknown choice_basal_roughness' // TRIM(C%choice_basal_roughness))
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness
  SUBROUTINE initialise_bed_roughness_from_file( grid, ice, region_name)
    ! Initialise bed roughness with data from an external NetCDF file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_sliding_law == 'no_sliding' .OR. &
            C%choice_sliding_law == 'idealised') THEN
      ! No sliding allowed / sliding laws for some idealised experiments
      CALL crash('not defined for choice_sliding_law' // TRIM(C%choice_sliding_law))
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      CALL initialise_bed_roughness_from_file_Weertman( grid, ice, region_name)
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Coulomb-type sliding law
      CALL initialise_bed_roughness_from_file_Coulomb( grid, ice, region_name)
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      CALL initialise_bed_roughness_from_file_Tsai2015( grid, ice, region_name)
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      CALL initialise_bed_roughness_from_file_Schoof2005( grid, ice, region_name)
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      CALL initialise_bed_roughness_from_file_ZoetIverson( grid, ice, region_name)
    ELSE
      CALL crash('unknown choice_sliding_law' // TRIM(C%choice_sliding_law))
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file
  SUBROUTINE initialise_bed_roughness_from_file_Weertman( grid, ice, region_name)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Weertman-type sliding law: bed roughness described by beta_sq

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file_Weertman'
    CHARACTER(LEN=256)                                 :: filename
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename
    filename = C%basal_roughness_filename

    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( filename), '...'
    CALL sync

    ! Read grid & bed roughness data from file
    CALL read_BIV_bed_roughness_file( filename, ice, grid, region_name)
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( ice%beta_sq_a,  'ice%beta_sq_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file_Weertman
  SUBROUTINE initialise_bed_roughness_from_file_Coulomb( grid, ice, region_name)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Coulomb-type sliding law: bed roughness described by phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file_Coulomb'
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: i,j
    REAL(dp)                                           :: r_smooth

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename
    filename = C%basal_roughness_filename

    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( filename), '...'

    ! Read grid & bed roughness data from file
    CALL read_BIV_bed_roughness_file( filename, ice, grid, region_name)

    ! Safety
    CALL check_for_NaN_dp_2D( ice%phi_fric_a, 'ice%phi_fric_a')

    ! Smooth input bed roughness (for restarts with a different resolution)
    IF (C%do_smooth_phi_restart) THEN
      ! Smooth with a 2-D Gaussian filter with a standard deviation of 1/2 grid cell
      r_smooth = grid%dx * C%r_smooth_phi_restart

      ! Apply smoothing to the bed roughness
      CALL smooth_Gaussian_2D( grid, ice%phi_fric_a, r_smooth)

      ! Limit bed roughness
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%phi_fric_a( j,i) = MAX( 0.1_dp, MIN( 30._dp, ice%phi_fric_a( j,i)))
      END DO
      END DO
      CALL sync
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file_Coulomb
  SUBROUTINE initialise_bed_roughness_from_file_Tsai2015( grid, ice, region_name)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Tsai 2015 sliding law: bed roughness described by alpha_sq & beta_sq

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file_Tsai2015'
    CHARACTER(LEN=256)                                 :: filename

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename
    filename = C%basal_roughness_filename

    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( filename), '...'
    CALL sync

    ! Read grid & bed roughness data from file
    CALL read_BIV_bed_roughness_file( filename, ice, grid, region_name)
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( ice%alpha_sq_a, 'ice%alpha_sq_a')
    CALL check_for_NaN_dp_2D( ice%beta_sq_a,  'ice%beta_sq_a' )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file_Tsai2015
  SUBROUTINE initialise_bed_roughness_from_file_Schoof2005( grid, ice, region_name)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Schoof 2005 sliding law: bed roughness described by alpha_sq & beta_sq

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file_Schoof2005'
    CHARACTER(LEN=256)                                 :: filename

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename
    filename = C%basal_roughness_filename

    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( filename), '...'
    CALL sync

    ! Read grid & bed roughness data from file
    CALL read_BIV_bed_roughness_file( filename, ice, grid, region_name)
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( ice%alpha_sq_a, 'ice%alpha_sq_a')
    CALL check_for_NaN_dp_2D( ice%beta_sq_a,  'ice%beta_sq_a' )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file_Schoof2005
  SUBROUTINE initialise_bed_roughness_from_file_ZoetIverson( grid, ice, region_name)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Zoet-Iverson sliding law: bed roughness described by phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file_ZoetIverson'
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: i,j
    REAL(dp)                                           :: r_smooth

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename
    filename = C%basal_roughness_filename

    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( filename), '...'
    CALL sync

    ! Read grid & bed roughness data from file
    CALL read_BIV_bed_roughness_file( filename, ice, grid, region_name)
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( ice%phi_fric_a, 'ice%phi_fric_a')

    ! Smooth input bed roughness (for restarts with a different resolution)
    IF (C%do_smooth_phi_restart) THEN
      ! Smooth with a 2-D Gaussian filter with a standard deviation of 1/2 grid cell
      r_smooth = grid%dx * C%r_smooth_phi_restart

      ! Apply smoothing to the bed roughness
      CALL smooth_Gaussian_2D( grid, ice%phi_fric_a, r_smooth)

      ! Limit bed roughness
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%phi_fric_a( j,i) = MAX( 0.1_dp, MIN( 30._dp, ice%phi_fric_a( j,i)))
      END DO
      END DO
      CALL sync
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file_ZoetIverson

  ! The Martin et al. (2011) till parameterisation
  SUBROUTINE calc_bed_roughness_Martin2011( grid, ice)
    ! Calculate the till friction angle phi_fric and till yield stress tauc,
    ! using the till model by Martin et al. (2011).
    !
    ! Only applicable when choice_sliding_law = "Coulomb" or "Coulomb_regularised"

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness_Martin2011'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: w_Hb

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. (C%choice_sliding_law == 'Coulomb' .OR. C%choice_sliding_law == 'Coulomb_regularised' .OR. C%choice_sliding_law == 'Zoet-Iverson')) THEN
      CALL crash('only applicable when choice_sliding_law = "Coulomb", "Coulomb_regularised", or "Zoet-Iverson"!')
    END IF

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Martin et al. (2011) Eq. 10
      w_Hb = MIN( 1._dp, MAX( 0._dp, (ice%Hb_a( j,i) - C%Martin2011till_phi_Hb_min) / (C%Martin2011till_phi_Hb_max - C%Martin2011till_phi_Hb_min) ))
      ice%phi_fric_a( j,i) = (1._dp - w_Hb) * C%Martin2011till_phi_min + w_Hb * C%Martin2011till_phi_max

    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( ice%phi_fric_a, 'ice%phi_fric_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness_Martin2011

  ! Idealised cases
  SUBROUTINE calc_bed_roughness_SSA_icestream( grid, ice)
    ! Determine the basal conditions underneath the ice
    !
    ! Idealised case: SSA_icestream (i.e. the Schoof 2006 analytical solution)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness_SSA_icestream'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: dummy1

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      CALL SSA_Schoof2006_analytical_solution( ABS(ice%dHs_dx_a( j,i)), ice%Hi_a( j,i), ice%A_flow_vav_a( j,i), grid%y(j), dummy1, ice%tauc_a( j,i))
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness_SSA_icestream
  SUBROUTINE calc_bed_roughness_MISMIPplus( grid, ice)
    ! Determine the basal conditions underneath the ice
    !
    ! Idealised case: MISMIP+ (see Asay-Davis et al., 2016)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness_MISMIPplus'
    INTEGER                                            :: i,j
    REAL(dp), PARAMETER                                :: MISMIPplus_alpha_sq = 0.5_dp   ! Coulomb-law friction coefficient [unitless];         see Asay-Davis et al., 2016
    REAL(dp), PARAMETER                                :: MISMIPplus_beta_sq  = 1.0E4_dp ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3]; idem dito

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_sliding_law == 'Weertman') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the first (Weertman) sliding law option

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%beta_sq_a( j,i) = MISMIPplus_beta_sq
      END DO
      END DO
      CALL sync

    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the second (Tsai et al., 2015) sliding law option

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%alpha_sq_a( j,i) = MISMIPplus_alpha_sq
        ice%beta_sq_a(  j,i) = MISMIPplus_beta_sq
      END DO
      END DO
      CALL sync

    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the third (Schoof, 2005) sliding law option

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%alpha_sq_a( j,i) = MISMIPplus_alpha_sq
        ice%beta_sq_a(  j,i) = MISMIPplus_beta_sq
      END DO
      END DO
      CALL sync

    ELSE
      CALL crash('only defined when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness_MISMIPplus

! == Sliding laws
! ===============

  SUBROUTINE calc_sliding_law( grid, ice, u_a, v_a, beta_a)
    ! Calculate the sliding term beta in the SSA/DIVA using the specified sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No sliding allowed (choice of beta is trivial)
      beta_a( :,grid%i1:grid%i2) = 0._dp
      CALL sync
    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! Sliding laws for some idealised experiments
      CALL calc_sliding_law_idealised(           grid, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      CALL calc_sliding_law_Weertman(            grid, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Coulomb') THEN
      ! Coulomb-type sliding law
      CALL calc_sliding_law_Coulomb(             grid, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Regularised Coulomb-type sliding law
      CALL calc_sliding_law_Coulomb_regularised( grid, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      CALL calc_sliding_law_Tsai2015(            grid, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      CALL calc_sliding_law_Schoof2005(          grid, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      CALL calc_sliding_law_ZoetIverson(         grid, ice, u_a, v_a, beta_a)
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM(C%choice_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law

  SUBROUTINE calc_sliding_law_Weertman( grid, ice, u_a, v_a, beta_a)
    ! Weertman-type ("power law") sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Weertman'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)

      ! Asay-Davis et al. (2016), Eq. 6
      beta_a( j,i) = ice%beta_sq_a( j,i) * uabs_a ** (1._dp / C%slid_Weertman_m - 1._dp)

    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Weertman
  SUBROUTINE calc_sliding_law_Coulomb( grid, ice, u_a, v_a, beta_a)
    ! Coulomb-type sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Coulomb'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine the till yield stress
    IF (C%choice_basal_roughness == 'parameterised' .AND. C%choice_param_basal_roughness == 'SSA_icestream') THEN
      ! In this case, tauc has already been calculated
    ELSE
      ! Calculate the till yield stress from the till friction angle and the effective pressure

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%tauc_a( j,i) = TAN((pi / 180._dp) * ice%phi_fric_a( j,i)) * ice%Neff_a( j,i)
      END DO
      END DO

    END IF

    ! Calculate the basal friction coefficient
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)

      beta_a( j,i) = ice%tauc_a( j,i) / uabs_a

    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Coulomb
  SUBROUTINE calc_sliding_law_Coulomb_regularised( grid, ice, u_a, v_a, beta_a)
    ! Regularised Coulomb-type sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Coulomb_regularised'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Calculate the till yield stress from the till friction angle and the effective pressure
      ice%tauc_a( j,i) = TAN((pi / 180._dp) * ice%phi_fric_a( j,i)) * ice%Neff_a( j,i)

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)

      beta_a( j,i) = ice%tauc_a( j,i) * uabs_a ** (C%slid_Coulomb_reg_q_plastic - 1._dp) / (C%slid_Coulomb_reg_u_threshold ** C%slid_Coulomb_reg_q_plastic)

    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Coulomb_regularised
  SUBROUTINE calc_sliding_law_Tsai2015( grid, ice, u_a, v_a, beta_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Davis et al., 2016)
    !
    ! Asay-Davis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    !
    ! Tsai et al.: Marine ice-sheet profiles and stability under Coulomb basal conditions,
    ! Journal of Glaciology 61, 205–215, doi:10.3189/2015JoG14J221, 2015.

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Tsai2015'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)

      ! Asay-Davis et al. (2016), Eq. 7
      beta_a( j,i) = MIN( ice%alpha_sq_a( j,i) * ice%Neff_a( j,i), ice%beta_sq_a( j,i) * uabs_a ** (1._dp / C%slid_Weertman_m)) * uabs_a**(-1._dp)

    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Tsai2015
  SUBROUTINE calc_sliding_law_Schoof2005( grid, ice, u_a, v_a, beta_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Davis et al., 2016)
    !
    ! Asay-Davis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    !
    ! Schoof: The effect of cavitation on glacier sliding, P. Roy. Soc. A-Math. Phy., 461, 609–627, doi:10.1098/rspa.2004.1350, 2005

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Schoof2005'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a, alpha_sq, beta_sq

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)

      ! Abbreviations for shorter code
      alpha_sq = ice%alpha_sq_a( j,i)
      beta_sq  = ice%beta_sq_a(  j,i)

      ! Asay-Davis et al. (2016), Eq. 11
      beta_a( j,i) = ((beta_sq * uabs_a**(1._dp / C%slid_Weertman_m) * alpha_sq * ice%Neff_a( j,i)) / &
        ((beta_sq**C%slid_Weertman_m * uabs_a + (alpha_sq * ice%Neff_a( j,i))**C%slid_Weertman_m)**(1._dp / C%slid_Weertman_m))) * uabs_a**(-1._dp)

    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Schoof2005
  SUBROUTINE calc_sliding_law_ZoetIverson( grid, ice, u_a, v_a, beta_a)
    ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_ZoetIverson'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)

      ! Zoet & Iverson (2020), Eq. (3) (divided by u to give beta = tau_b / u)
      beta_a( j,i) = ice%Neff_a( j,i) * TAN((pi / 180._dp) * ice%phi_fric_a( j,i)) * (uabs_a**(1._dp / C%slid_ZI_p - 1._dp)) * ((uabs_a + C%slid_ZI_ut)**(-1._dp / C%slid_ZI_p))

    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_ZoetIverson

  SUBROUTINE calc_sliding_law_idealised( grid, ice, u_a, v_a, beta_a)
    ! Sliding laws for some idealised experiments

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised'
    REAL(dp) :: dummy_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings...
    dummy_dp = u_a( 1,1)
    dummy_dp = v_a( 1,1)

    IF     (C%choice_idealised_sliding_law == 'ISMIP_HOM_C') THEN
      ! ISMIP-HOM experiment C

      CALL calc_sliding_law_idealised_ISMIP_HOM_C( grid, beta_a)

    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_D') THEN
      ! ISMIP-HOM experiment D

      CALL calc_sliding_law_idealised_ISMIP_HOM_D( grid, beta_a)

    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_E') THEN
      ! ISMIP-HOM experiment E

      CALL calc_sliding_law_idealised_ISMIP_HOM_E( grid, beta_a)

    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_F') THEN
      ! ISMIP-HOM experiment F

      CALL calc_sliding_law_idealised_ISMIP_HOM_F( grid, ice, beta_a)

    ELSE
      CALL crash('unknown choice_idealised_sliding_law "' // TRIM(C%choice_idealised_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C( grid, beta_a)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment C

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_C'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      beta_a( j,i) = 1000._dp + 1000._dp * SIN( 2._dp * pi * grid%x( i) / C%ISMIP_HOM_L) * SIN( 2._dp * pi * grid%y( j) / C%ISMIP_HOM_L)
    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D( grid, beta_a)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment D

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_D'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      beta_a( j,i) = 1000._dp + 1000._dp * SIN( 2._dp * pi * grid%x( i) / C%ISMIP_HOM_L)
    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_E( grid, beta_a)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment E: use the externally prescribed slip zone in Haut Glacier d'Arolla

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_E'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: dummy1, dummy2, dummy3
    INTEGER                                            :: ios,slides

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN
      OPEN(   UNIT = 1337, FILE=C%ISMIP_HOM_E_Arolla_filename, ACTION='READ')
      DO i = 1, 51
        READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy1, dummy2, dummy3, slides
        IF (slides == 1) THEN
          DO j = 1, grid%ny
            beta_a( j,i) = 0._dp
          END DO
        ELSE
          DO j = 1, grid%ny
            beta_a( j,i) = 1.0E30_dp
          END DO
        END IF
        IF (ios /= 0) THEN
          CALL crash('length of text file "' // TRIM(C%ISMIP_HOM_E_Arolla_filename) // '" should be 51 lines!')
        END IF
      END DO
      CLOSE( UNIT  = 1337)
    END IF ! IF (par%master) THEN
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_E
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F( grid, ice, beta_a)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment F

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_F'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      beta_a( j,i) = (ice%A_flow_vav_a( j,i) * 1000._dp)**(-1._dp)
    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F

! == Basal inversion routines
! ===========================

  SUBROUTINE basal_inversion_geo( grid, ice, refgeo_PD, dt)
    ! Invert for bed roughness based on modelled vs. target geometry

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    REAL(dp),                            INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'basal_inversion_geo'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_BIVgeo_method == 'PDC2012') THEN
      ! Update the bed roughness according to Pollard & DeConto (2012)

      CALL update_bed_roughness_PDC2012( grid, ice, refgeo_PD)

    ELSEIF (C%choice_BIVgeo_method == 'Lipscomb2021') THEN
      ! Update the bed roughness according to Lipscomb et al. (2021)

      CALL update_bed_roughness_Lipscomb2021( grid, ice, refgeo_PD, dt)

    ELSEIF (C%choice_BIVgeo_method == 'CISM+') THEN
      ! Update the bed roughness according to the geometry+velocity-based CISM+ approach

      CALL update_bed_roughness_CISMplus( grid, ice, refgeo_PD, dt)

    ELSEIF (C%choice_BIVgeo_method == 'Berends2022') THEN
      ! PIEP

      CALL update_bed_roughness_Berends2022( grid, ice, refgeo_PD, dt)

    ELSEIF (C%choice_BIVgeo_method == 'Bernales2017') THEN
      ! Update the bed roughness according to Bernales et al. (2017)

      CALL update_bed_roughness_Bernales2017( grid, ice, refgeo_PD)

    ELSEIF (C%choice_BIVgeo_method == 'Pien2023') THEN
      ! Update the bed roughness according to Pien and friends (2023)

      CALL update_bed_roughness_Pien2023( grid, ice, refgeo_PD, dt)

    ELSE
      CALL crash('unknown choice_BIVgeo_method "' // TRIM(C%choice_BIVgeo_method) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE basal_inversion_geo

  SUBROUTINE update_bed_roughness_PDC2012( grid, ice, refgeo_PD)
    ! Update the bed roughness according to Pollard & DeConto (2012)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_bed_roughness_PDC2012'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_sliding_law == 'Coulomb' .OR. &
        C%choice_sliding_law == 'Coulomb_regularised') THEN

      CALL update_bed_roughness_PDC2012_Coulomb( grid, ice, refgeo_PD)

    ELSE
      CALL crash('not implemented for sliding law "' // TRIM(C%choice_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_bed_roughness_PDC2012

  SUBROUTINE update_bed_roughness_PDC2012_Coulomb( grid, ice, refgeo_PD)
    ! Update the bed roughness according to Pollard & DeConto (2012)
    !
    ! For the case of a (regularised) Coulomb sliding law: update phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_bed_roughness_PDC2012_Coulomb'
    INTEGER                                            :: i,j
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  mask
    INTEGER                                            :: wmask
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dz
    INTEGER                                            :: wdz
    REAL(dp), PARAMETER                                :: sigma = 500._dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx, mask, wmask)
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx, dz  , wdz  )

    ! Calculate dz (Pollard % DeConto (2012), just after Eq. 5)
    DO i = grid%i1, grid%I2
    DO j = 1, grid%ny

      IF (ice%mask_sheet_a( j,i) == 1) THEN

        mask( j,i) = 2

        dz( j,i) = MAX( -1.5_dp, MIN( 1.5_dp, (ice%Hs_a( j,i) - refgeo_PD%Hs( j,i)) / C%BIVgeo_PDC2012_hinv ))

      ELSE
        mask( j,i) = 1
      END IF

    END DO
    END DO
    CALL sync

    ! Extrapolate new values outside the ice sheet
    CALL extrapolate_updated_bed_roughness( grid, mask, dz)

    ! Apply regularisation
    CALL smooth_Gaussian_2D( grid, dz, sigma)

    ! Update bed roughness
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ! Pollard % DeConto (2012), Eq. 5
      ice%phi_fric_a( j,i) = MAX( 0.001_dp, MIN( 30._dp, ice%phi_fric_a( j,i) * (10._dp**(-dz( j,i))) ))
    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wdz  )
    CALL deallocate_shared( wmask)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_bed_roughness_PDC2012_Coulomb

  SUBROUTINE update_bed_roughness_Lipscomb2021( grid, ice, refgeo_PD, dt)
    ! Update the bed roughness according to Lipscomb et al. (2021)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    REAL(dp),                            INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_bed_roughness_Lipscomb2021'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_sliding_law == 'Coulomb' .OR. &
        C%choice_sliding_law == 'Coulomb_regularised') THEN

      CALL update_bed_roughness_Lipscomb2021_Coulomb( grid, ice, refgeo_PD, dt)

    ELSE
      CALL crash('not implemented for sliding law "' // TRIM(C%choice_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_bed_roughness_Lipscomb2021

  SUBROUTINE update_bed_roughness_Lipscomb2021_Coulomb( grid, ice, refgeo_PD, dt)
    ! Update the bed roughness according to Lipscomb et al. (2021)
    !
    ! For the case of a (regularised) Coulomb sliding law: update phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    REAL(dp),                            INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_bed_roughness_Lipscomb2021_Coulomb'
    INTEGER                                            :: i,j
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  mask
    INTEGER                                            :: wmask
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dCdt
    INTEGER                                            :: wdCdt
    REAL(dp), PARAMETER                                :: sigma = 500._dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx, mask, wmask)
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx, dCdt, wdCdt)

    ! Calculate dCdt
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (ice%mask_sheet_a( j,i) == 1) THEN

        mask( j,i) = 2

        ! Lipscomb et al. (2021), Eq. 8 (slightly altered for Coulomb friction)
        dCdt( j,i) = -ice%phi_fric_a( j,i) / C%BIVgeo_Lipscomb2021_H0 * ( (ice%Hi_a( j,i) - refgeo_PD%Hi( j,i)) / C%BIVgeo_Lipscomb2021_tauc + 2._dp * ice%dHi_dt_a( j,i))

      ELSE
        mask( j,i) = 1
      END IF

    END DO
    END DO
    CALL sync

    ! Extrapolate new values outside the ice sheet
    CALL extrapolate_updated_bed_roughness( grid, mask, dCdt)

    ! Apply regularisation
    CALL smooth_Gaussian_2D( grid, dCdt, sigma)

    ! Update bed roughness
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%phi_fric_a( j,i) = MAX( 1E-6_dp, MIN( 30._dp, ice%phi_fric_a( j,i) + dCdt( j,i) * dt ))
    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wdCdt)
    CALL deallocate_shared( wmask)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_bed_roughness_Lipscomb2021_Coulomb

  SUBROUTINE update_bed_roughness_CISMplus( grid, ice, refgeo_PD, dt)
    ! Update the bed roughness according to the geometry+velocity-based CISM+ approach

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    REAL(dp),                            INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_bed_roughness_CISMplus'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_sliding_law == 'Coulomb' .OR. &
        C%choice_sliding_law == 'Coulomb_regularised' .OR. &
        C%choice_sliding_law == 'Zoet-Iverson') THEN

      CALL update_bed_roughness_CISMplus_Coulomb( grid, ice, refgeo_PD, dt)

    ELSE
      CALL crash('not implemented for sliding law "' // TRIM(C%choice_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_bed_roughness_CISMplus

  SUBROUTINE update_bed_roughness_CISMplus_Coulomb( grid, ice, refgeo_PD, dt)
    ! Update the bed roughness according to the geometry+velocity-based CISM+ approach
    !
    ! For the case of a (regularised) Coulomb / Zoet-Iverson sliding law: update phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    REAL(dp),                            INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_bed_roughness_CISMplus_Coulomb'
    INTEGER                                            :: i,j
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  mask
    INTEGER                                            :: wmask
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dCdt
    INTEGER                                            :: wdCdt
    REAL(dp), PARAMETER                                :: sigma = 500._dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx, mask, wmask)
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx, dCdt, wdCdt)

    ! Calculate dCdt
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (ice%mask_sheet_a( j,i) == 1) THEN

        mask( j,i) = 2

        ! CISM+ inversion equation
        dCdt( j,i) = -ice%phi_fric_a( j,i) / C%BIVgeo_CISMplus_tauc * &
          ( C%BIVgeo_CISMplus_wH * (ice%Hi_a(        j,i) - refgeo_PD%Hi(           j,i)) / ( C%BIVgeo_CISMplus_H0) + &
            C%BIVgeo_CISMplus_wu * (ice%uabs_surf_a( j,i) - ice%BIV_uabs_surf_target( j,i)) / (-C%BIVgeo_CISMplus_u0) )

      ELSE
        mask( j,i) = 1
      END IF

    END DO
    END DO
    CALL sync

    ! Extrapolate new values outside the ice sheet
    CALL extrapolate_updated_bed_roughness( grid, mask, dCdt)

    ! Apply regularisation
    CALL smooth_Gaussian_2D( grid, dCdt, sigma)

    ! Update bed roughness
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%phi_fric_a( j,i) = MAX( 1E-6_dp, MIN( 30._dp, ice%phi_fric_a( j,i) + dCdt( j,i) * dt ))
    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wdCdt)
    CALL deallocate_shared( wmask)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_bed_roughness_CISMplus_Coulomb

  SUBROUTINE update_bed_roughness_Berends2022( grid, ice, refgeo_PD, dt)
    ! The geometry+velocity-based inversion from Berends et al. (2022)
    ! Here, the rate of change of the bed roughness is calculated
    ! by integrating the difference between modelled vs. target geometry+velocity
    ! along the flowline (both upstream and downstream).

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    REAL(dp),                            INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_bed_roughness_Tijn'
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  mask
    INTEGER                                            :: wmask
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dphi_dt
    INTEGER                                            :: wdphi_dt
    REAL(dp), PARAMETER                                :: dx_trace_rel    = 0.25_dp       ! Trace step size relative to grid resolution
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: trace_up, trace_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: w_up, w_down
    INTEGER                                            :: i,j,n_up,n_down,k
    REAL(dp), DIMENSION(2)                             :: p,pt
    REAL(dp)                                           :: Hs_mod, Hs_target, u_mod, u_target
    REAL(dp)                                           :: I1, I2, I3, I_tot
    REAL(dp)                                           :: R
    REAL(dp)                                           :: h_delta, h_dfrac
    REAL(dp)                                           :: sigma

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx, mask   , wmask   )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx, dphi_dt, wdphi_dt)

    ! Allocate memory for the up- and downstream traces
    ALLOCATE( trace_up(        2*MAX(grid%nx,grid%ny), 2))
    ALLOCATE( trace_down(      2*MAX(grid%nx,grid%ny), 2))

    ! Allocate memory for the linear scaling functions
    ALLOCATE( w_up(   2*MAX(grid%nx,grid%ny)   ))
    ALLOCATE( w_down( 2*MAX(grid%nx,grid%ny)   ))

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Obviously can only be done where there's (grounded) ice
      IF (ice%mask_sheet_a( j,i) == 0 .OR. ice%mask_margin_a( j,i) == 1 .OR. ice%Hi_a( j,i) < 1._dp .OR. &
          refgeo_PD%Hi( j,i) < 1._dp .OR. ice%BIV_uabs_surf_target( j,i) == 0._dp .OR. ice%mask_gl_a( j,i) == 1) THEN
        mask( j,i) = 1
        CYCLE
      ELSE
        mask( j,i) = 2
      END IF

      ! The point p
      p = [grid%x( i), grid%y( j)]

      ! Trace the flowline upstream and downstream (Berends et al., 2022, Eqs. 2)
      CALL trace_flowline_upstream(   grid, ice, p, dx_trace_rel, trace_up  , n_up               )
      CALL trace_flowline_downstream( grid, ice, p, dx_trace_rel, trace_down, n_down, refgeo_PD)

      ! Calculate linear scaling functions (Berends et al., Eqs. 5)
      w_up = 0._dp
      DO k = 1, n_up
        w_up( k) = REAL( n_up + 1 - k, dp)
      END DO
      w_up = w_up / SUM( w_up)

      w_down = 0._dp
      DO k = 1, n_down
        w_down( k) = REAL( n_down + 1 - k, dp)
      END DO
      w_down = w_down / SUM( w_down)

      ! Calculate upstream line integrals
      I1 = 0._dp
      I3 = 0._dp
      DO k = 1, n_up

        pt = trace_up( k,:)
        u_mod     = interp_bilin_2D( ice%uabs_surf_a         , grid%x, grid%y, pt(1), pt(2))
        u_target  = interp_bilin_2D( ice%BIV_uabs_surf_target, grid%x, grid%y, pt(1), pt(2))
        Hs_mod    = interp_bilin_2D( ice%Hs_a                , grid%x, grid%y, pt(1), pt(2))
        Hs_target = interp_bilin_2D( refgeo_PD%Hs          , grid%x, grid%y, pt(1), pt(2))

        ! If no target velocity data is available, assume zero difference
        IF (u_target /= u_target) u_target = u_mod

        I1 = I1 - ( u_mod -  u_target) * w_up( k) / C%BIVgeo_Berends2022_u0    ! Berends et al., (2022), Eq. 4a
        I3 = I3 + (Hs_mod - Hs_target) * w_up( k) / C%BIVgeo_Berends2022_H0    ! Berends et al., (2022), Eq. 4c

      END DO

      ! Calculate downstream line integral
      I2 = 0._dp
      DO k = 1, n_down

        pt = trace_down( k,:)
        u_mod     = interp_bilin_2D( ice%uabs_surf_a         , grid%x, grid%y, pt(1), pt(2))
        u_target  = interp_bilin_2D( ice%BIV_uabs_surf_target, grid%x, grid%y, pt(1), pt(2))

        ! If no target velocity data is available, assume zero difference
        IF (u_target /= u_target) u_target = u_mod

        I2 = I2 - ( u_mod -  u_target) * w_down( k) / C%BIVgeo_Berends2022_u0  ! Berends et al., (2022), Eq. 4b

      END DO

      ! Scale weights with local ice thickness * velocity
      ! (thinner and/or ice experiences less basal friction, so the solution is less sensitive to changes in bed roughness there)

      ! Berends et al. (2022), Eq. 7
      R = MAX( 0._dp, MIN( 1._dp, ((ice%uabs_surf_a( j,i) * ice%Hi_a( j,i)) / (C%BIVgeo_Berends2022_u_scale * C%BIVgeo_Berends2022_Hi_scale)) ))
      ! Berends et al. (2022), Eq. 6
      I_tot = (I1 + I2 + I3) * R

      ! Ice thickness difference w.r.t. reference thickness
      h_delta = ice%Hi_a(j,i) - refgeo_PD%Hi(j,i)
      ! Ratio between this difference and the reference ice thickness
      h_dfrac = h_delta / MAX( refgeo_PD%Hi(j,i), 1._dp)

      ! If the difference/fraction is outside the specified tolerance
      IF (ABS( h_delta) >= C%BIVgeo_Bernales2017_tol_diff .OR. &
          ABS( h_dfrac) >= C%BIVgeo_Bernales2017_tol_frac) THEN

        ! Further adjust only where the previous value is not improving the result
        IF ( (h_delta >= 0._dp .AND. ice%dHi_dt_a( j,i) >= 0.0_dp) .OR. &
             (h_delta <= 0._dp .AND. ice%dHi_dt_a( j,i) <= 0.0_dp) ) THEN

          ! Calculate rate of change of bed roughness (Berends et al. (2022), Eq. 8)
          dphi_dt( j,i) = -ice%phi_fric_a( j,i) * I_tot / C%BIVgeo_Berends2022_tauc

        END IF
      END IF

    END DO
    END DO
    CALL sync

    ! Extrapolate new values outside the ice sheet
    CALL extrapolate_updated_bed_roughness( grid, mask, dphi_dt)

    ! Update bed roughness (Berends et al. (2022), Eq. 9)

    ! First regularisation step (function F1 in Berends et al. (2022), Eq. 9)
    sigma = grid%dx / 1.5_dp
    CALL smooth_Gaussian_2D( grid, dphi_dt, sigma)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%phi_fric_a( j,i) = MAX( C%BIVgeo_Berends2022_phimin, MIN( C%BIVgeo_Berends2022_phimax, ice%phi_fric_a( j,i) + dphi_dt( j,i) * dt ))
    END DO
    END DO
    CALL sync

    ! Second regularisation step (function F2 in Berends et al. (2022), Eq. 9)
    sigma = grid%dx / 4._dp
    CALL smooth_Gaussian_2D( grid, ice%phi_fric_a, sigma)

    ! Clean up after yourself
    DEALLOCATE( trace_up  )
    DEALLOCATE( trace_down)
    DEALLOCATE( w_up      )
    DEALLOCATE( w_down    )
    CALL deallocate_shared( wmask   )
    CALL deallocate_shared( wdphi_dt)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_bed_roughness_Berends2022

  SUBROUTINE trace_flowline_upstream( grid, ice, p, dx_trace_rel, T, n)
    ! Trace the flowline passing through point p upstream.
    ! Returns a list T of n points on the flowline spaced dx_trace_rel * grid%dx apart.
    !
    ! Stop the trace when it encounters the ice divide (defined as an ice velocity lower than 1 m/yr)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    REAL(dp),                            INTENT(IN)    :: dx_trace_rel
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: T
    INTEGER,                             INTENT(OUT)   :: n

    ! Local variables:
    INTEGER                                            :: ii,jj,it,nmax
    REAL(dp), DIMENSION(2)                             :: pt, u

    nmax = SIZE( T,1)

    ! Initialise
    T  = 0._dp
    n  = 0
    pt = p
    it = 0

    DO WHILE (n < nmax)

      ! Find current tracer location grid indices
      ii = MAX( 1, MIN( grid%nx-1, 1 + FLOOR(( pt(1) + grid%dx/2._dp - grid%xmin) / grid%dx)))
      jj = MAX( 1, MIN( grid%ny-1, 1 + FLOOR(( pt(2) + grid%dx/2._dp - grid%ymin) / grid%dx)))

      ! Interpolate surface velocity to the tracer location
      u( 1) = interp_bilin_2D( ice%u_surf_a, grid%x, grid%y, pt(1), pt(2))
      u( 2) = interp_bilin_2D( ice%v_surf_a, grid%x, grid%y, pt(1), pt(2))

      ! If we've reached the ice divide, end the trace
      IF (NORM2( u) < 1._dp) EXIT

      ! Add current position to the traces
      n = n + 1
      T( n,:) = pt

      ! Normalise velocity vector
      u = u / NORM2( u)

      ! Move the tracer upstream
      pt = pt - u * dx_trace_rel * grid%dx

      ! Safety
      it = it + 1
      IF (it > nmax) EXIT

    END DO ! DO WHILE (n_up < MAX(grid%nx, grid%ny))

    ! Safety
    IF (n == 0) THEN
      n = 1
      T( 1,:) = p
    END IF

  END SUBROUTINE trace_flowline_upstream

  SUBROUTINE trace_flowline_downstream( grid, ice, p, dx_trace_rel, T, n, refgeo_PD)
    ! Trace the flowline passing through point p downstream.
    ! Returns a list T of n points on the flowline spaced dx_trace_rel * grid%dx apart.
    !
    ! Stop the trace when it encounters the ice margin (either in the modelled or the target geometry)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    REAL(dp),                            INTENT(IN)    :: dx_trace_rel
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: T
    INTEGER,                             INTENT(OUT)   :: n
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD

    ! Local variables:
    INTEGER                                            :: ii,jj,it,nmax
    REAL(dp), DIMENSION(2)                             :: pt, u

    nmax = SIZE( T,1)

    ! Initialise
    T  = 0._dp
    n  = 0
    pt = p
    it = 0

    DO WHILE (n < nmax)

      ! Find current tracer location grid indices
      ii = MAX( 1, MIN( grid%nx-1, 1 + FLOOR(( pt(1) + grid%dx/2._dp - grid%xmin) / grid%dx)))
      jj = MAX( 1, MIN( grid%ny-1, 1 + FLOOR(( pt(2) + grid%dx/2._dp - grid%ymin) / grid%dx)))

      ! If we've reached the ice margin, end the trace
      IF (ice%mask_ice_a( jj,ii) == 0 .OR. ice%mask_margin_a( jj,ii) == 1 .OR. ice%Hi_a( jj,ii) < 1._dp .OR. &
          refgeo_PD%Hi( jj,ii) < 1._dp .OR. ice%BIV_uabs_surf_target( jj,ii) == 0._dp) EXIT

      ! Interpolate surface velocity to the tracer location
      u( 1) = interp_bilin_2D( ice%u_surf_a, grid%x, grid%y, pt(1), pt(2))
      u( 2) = interp_bilin_2D( ice%v_surf_a, grid%x, grid%y, pt(1), pt(2))

      ! Safety
      IF (u( 1) == 0._dp .AND. u( 2) == 0._dp) EXIT

      ! Add current position to the traces
      n = n + 1
      T( n,:) = pt

      ! Normalise velocity vector
      u = u / NORM2( u)

      ! Move the tracer downstream
      pt = pt + u * dx_trace_rel * grid%dx

      ! Safety
      it = it + 1
      IF (it > nmax) EXIT

    END DO ! DO WHILE (n_down < MAX(grid%nx, grid%ny))

    ! Safety
    IF (n == 0) THEN
      n = 1
      T( 1,:) = p
    END IF

  END SUBROUTINE trace_flowline_downstream

  SUBROUTINE update_bed_roughness_Bernales2017( grid, ice, refgeo_PD)
    ! Update the bed roughness according to Bernales et al. (2017)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_bed_roughness_Bernales2017'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_sliding_law == 'Coulomb' .OR. &
        C%choice_sliding_law == 'Coulomb_regularised' .OR. &
        C%choice_sliding_law == 'Zoet-Iverson') THEN

      CALL update_bed_roughness_Bernales2017_CoulombZoetIverson( grid, ice, refgeo_PD)

    ELSE
      CALL crash('not implemented for sliding law "' // TRIM(C%choice_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_bed_roughness_Bernales2017

  SUBROUTINE update_bed_roughness_Bernales2017_CoulombZoetIverson( grid, ice, refgeo)
    ! Update the bed roughness according to Bernales et al. (2017)
    !
    ! For the case of a (regularised) Coulomb or Zoet-Iverson sliding law: update phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_bed_roughness_Bernales2017_CoulombZoetIverson'
    INTEGER                                            :: i,j
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  mask
    INTEGER                                            :: wmask
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dz
    INTEGER                                            :: wdz
    REAL(dp), PARAMETER                                :: sigma = 500._dp
    REAL(dp)                                           :: h_scale, h_delta, h_dfrac

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx, mask, wmask)
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx, dz  , wdz  )

    ! Define the ice thickness factor for scaling of inversion
    h_scale = 1.0_dp/C%BIVgeo_Bernales2017_hinv

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Ice thickness difference w.r.t. reference thickness
      h_delta = ice%Hi_a(j,i) - refgeo%Hi(j,i)
      ! Ratio between this difference and the reference ice thickness
      h_dfrac = h_delta / MAX(refgeo%Hi(j,i), 1._dp)

      ! Invert only where the model has grounded ice
      IF (ice%mask_sheet_a(j,i) == 1) THEN

        ! Check if grounded point is marginal or interior ice.
        ! If marginal, override its value during extrapolation.
        IF (ice%mask_margin_a(j,i) == 1 .OR. ice%mask_gl_a(j,i) == 1) THEN
          ! Mark this point as ice margin/grounding-line
          mask(j,i) = 1
        ELSE
          ! Mark this vertex as grounded ice
          mask(j,i) = 2
        END IF

        ! If the difference/fraction is outside the specified tolerance
        IF (ABS(h_delta) >= C%BIVgeo_Bernales2017_tol_diff .OR. &
            ABS(h_dfrac) >= C%BIVgeo_Bernales2017_tol_frac) THEN

          ! Further adjust only where the previous value is not improving the result
          IF ( (h_delta > 0._dp .AND. ice%dHi_dt_a(j,i) >= 0.0_dp) .OR. &
               (h_delta < 0._dp .AND. ice%dHi_dt_a(j,i) <= 0.0_dp) ) THEN

            ! Scale the difference and restrict it to the [-1.5 1.5] range
            dz(j,i) = MAX(-1.5_dp, MIN(1.5_dp, h_delta * h_scale))

          END IF ! else the fit is already improving for some other reason, so leave it alone

        END IF ! else the difference is within the specified tolerance, so leave it alone

      ELSE

        ! This vertex is not grounded ice sheet, so mark it for later extrapolation
        mask(j,i) = 1

      END IF ! (ice%mask_sheet_a(j,i) == 1)

    END DO
    END DO
    CALL sync

    ! Extrapolate new values outside the ice sheet
    CALL extrapolate_updated_bed_roughness( grid, mask, dz)

    ! Apply regularisation
    CALL smooth_Gaussian_2D( grid, dz, sigma)

    ! Update bed roughness
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ! Pollard % DeConto (2012), Eq. 5
      ice%phi_fric_a( j,i) = MAX( 0.1_dp, MIN( 30._dp, ice%phi_fric_a( j,i) * (10._dp**(-dz( j,i))) ))
    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wdz  )
    CALL deallocate_shared( wmask)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_bed_roughness_Bernales2017_CoulombZoetIverson

  SUBROUTINE update_bed_roughness_Pien2023( grid, ice, refgeo, dt)
    ! Update bed roughness according to Pien and friends, following
    ! Lipscomb et al. (2021) and van dem Akker et al. (202X).

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),               INTENT(IN)    :: grid
    TYPE(type_ice_model),          INTENT(INOUT) :: ice
    TYPE(type_reference_geometry), INTENT(IN)    :: refgeo
    REAL(dp),                      INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'update_bed_roughness_Pien2023'
    INTEGER                                      :: i,j
    INTEGER, DIMENSION(:,:), POINTER             ::  mask
    INTEGER                                      :: wmask
    REAL(dp), DIMENSION(:,:), POINTER            ::  dC_dt
    INTEGER                                      :: wdC_dt
    REAL(dp), PARAMETER                          :: sigma = 500._dp
    REAL(dp)                                     :: h_delta, c_ratio, reg_1st, reg_2nd, reg_3rd

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx, dC_dt, wdC_dt)
    CALL allocate_shared_int_2D( grid%ny, grid%nx, mask, wmask  )

    ! Loop over all grid points within this process
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Skip non-grounded, non-interior grid points
      IF (ice%mask_sheet_a(j,i) == 0 .OR. &
          ice%mask_gl_a(j,i) == 0) THEN
        ! Mark for extrapolation
        mask( j,i) = 1
        ! Go to next grid point
        CYCLE
      END IF

      ! == Metrics
      ! ==========

      ! Ice thickness difference w.r.t. observations (positive if too thick)
      h_delta = ice%Hi_a(j,i) - refgeo%Hi(j,i)

      ! Ratio between modelled friction and a target relaxation friction
      c_ratio = ice%phi_fric_a( j,i) !/ ice%phi_fric_relax_a( j,i)

      ! == Regularisation
      ! =================

      ! First regularisation term
      reg_1st = h_delta / C%BIVgeo_Pien2023_H0 / C%BIVgeo_Pien2023_tau

      ! Second regularisation term
      reg_2nd = ice%dHi_dt_a( j,i) * 2._dp / C%BIVgeo_Pien2023_H0

      ! Third regularisation term
      reg_3rd = LOG(c_ratio) * C%BIVgeo_Pien2023_r / C%BIVgeo_Pien2023_tau

      ! == Adjustment
      ! =============

      ! Compute adjustment to friction field
      dC_dt( j,i) = -ice%phi_fric_a( j,i) * (reg_1st + reg_2nd - reg_3rd)

      ! Use this point as seed during extrapolation
      mask( j,i) = 2

    END DO
    END DO
    CALL sync

    ! Extrapolate new values outside the ice sheet
    CALL extrapolate_updated_bed_roughness( grid, mask, dC_dt)

    ! Apply regularisation
    CALL smooth_Gaussian_2D( grid, dC_dt, sigma)

    ! Apply adjustment
    IF (par%master) THEN
      ice%phi_fric_a = ice%phi_fric_a + dC_dt * dt
    END IF
    CALL sync

    ! Limit bed roughness
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%phi_fric_a( j,i) = MAX( ice%phi_fric_a( j,i), 1._dp)
      ice%phi_fric_a( j,i) = MIN( ice%phi_fric_a( j,i), 30._dp)
    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wdC_dt)
    CALL deallocate_shared( wmask )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_bed_roughness_Pien2023

  SUBROUTINE extrapolate_updated_bed_roughness( grid, mask, d)
    ! The different geometry-based basal inversion routine only yield values
    ! beneath grounded ice; extrapolate new values outward to cover the entire domain.
    !
    ! Note about the mask:
    !    2 = data provided
    !    1 = no data provided, fill allowed
    !    0 = no fill allowed
    ! (so basically this routine extrapolates data from the area
    !  where mask == 2 into the area where mask == 1)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extrapolate_updated_bed_roughness'
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  mask_filled
    INTEGER                                            :: wmask_filled
    REAL(dp)                                           :: sigma

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx, mask_filled, wmask_filled)

    ! Perform the extrapolation
    sigma = 40000._dp
    IF (par%master) CALL extrapolate_Gaussian_floodfill( grid, mask, d, sigma, mask_filled)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wmask_filled)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE extrapolate_updated_bed_roughness

  SUBROUTINE initialise_basal_inversion( grid, ice, region_name)
    ! Fill in the initial guess for the bed roughness

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_inversion'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If needed, initialise target velocity fields
    ! ============================================

    IF     (C%choice_BIVgeo_method == 'PDC2012' .OR. &
            C%choice_BIVgeo_method == 'Lipscomb2021' .OR. &
            C%choice_BIVgeo_method == 'Bernales2017' .OR. &
            C%choice_BIVgeo_method == 'Pien2023') THEN
      ! Not needed in these methods
    ELSEIF (C%choice_BIVgeo_method == 'CISM+' .OR. &
            C%choice_BIVgeo_method == 'Berends2022') THEN
      ! Needed in these methods

      CALL initialise_basal_inversion_target_velocity( grid, ice, region_name)

    ELSE
      CALL crash('unknown choice_BIVgeo_method "' // TRIM(C%choice_BIVgeo_method) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_inversion

  SUBROUTINE initialise_basal_inversion_target_velocity( grid, ice, region_name)
    ! Initialise the target velocity fields used in a velocity-based basal inversion routine

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_inversion_target_velocity'
    TYPE(type_BIV_target_velocity)                     :: BIV_target
    INTEGER                                            :: i,j,i1,i2

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL warning('Reading basal inversion target velocity field has not been thoroughly tested. Please check if the data is read correctly.')

    ! Determine filename
    IF     (C%choice_BIVgeo_method == 'CISM+' .OR. &
            C%choice_BIVgeo_method == 'Berends2022') THEN
      BIV_target%filename = C%BIVgeo_target_velocity_filename
    ELSE
      CALL crash('unknown choice_BIVgeo_method "' // TRIM(C%choice_BIVgeo_method) // '"!')
    END IF

    IF (par%master) WRITE(0,*) '  Initialising basal inversion target velocity from file ', TRIM( BIV_target%filename), '...'

    ! Read input field
    CALL read_BIV_target_velocity( BIV_target, region_name)
    CALL sync

    ! NOTE: do not check for NaNs in the target velocity field. The Rignot 2011 Antarctica velocity product
    !       has a lot of missing data points, indicated by NaN values. This is acceptable, the basal inversion
    !       routine can handle that.

    ! CALL check_for_NaN_dp_2D( BIV_target%u_surf, 'BIV_target%u_surf')
    ! CALL check_for_NaN_dp_2D( BIV_target%v_surf, 'BIV_target%v_surf')

    CALL allocate_shared_dp_2D( BIV_target%grid%ny, BIV_target%grid%nx, BIV_target%uabs_surf, BIV_target%wuabs_surf)

    ! Get absolute velocity
    CALL partition_list( BIV_target%grid%nx, par%i, par%n, i1, i2)
    DO i = i1, i2
    DO j = 1, BIV_target%grid%ny
      BIV_target%uabs_surf( i,j) = SQRT( BIV_target%u_surf( i,j)**2 + BIV_target%v_surf( i,j)**2)
    END DO
    END DO
    CALL sync

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%BIV_uabs_surf_target, ice%wBIV_uabs_surf_target)

    ! Map (transposed) raw data to the model grid
    CALL map_square_to_square_cons_2nd_order_2D( BIV_target%grid%nx, BIV_target%grid%ny, BIV_target%grid%x, BIV_target%grid%y, grid%nx, grid%ny, grid%x, grid%y, BIV_target%uabs_surf, ice%BIV_uabs_surf_target)


    ! Deallocate raw data
    CALL deallocate_grid(   BIV_target%grid      )
    CALL deallocate_shared( BIV_target%wu_surf   )
    CALL deallocate_shared( BIV_target%wv_surf   )
    CALL deallocate_shared( BIV_target%wuabs_surf)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_inversion_target_velocity

  SUBROUTINE write_inverted_bed_roughness_to_file( grid, ice)
    ! Create a new NetCDF file and write the inverted bed roughness to it

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_inverted_bed_roughness_to_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL create_BIV_bed_roughness_file( grid, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_inverted_bed_roughness_to_file

  SUBROUTINE read_BIV_bed_roughness_file( filename, ice, grid, region_name)

    IMPLICIT NONE

    ! Read the extrapolated ocean data netcdf file
    CHARACTER(LEN=256),                  INTENT(IN)    :: filename
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Input variables:

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_BIV_bed_roughness_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_sliding_law == 'no_sliding') THEN
      CALL crash('not defined for choice_sliding_law = "no_sliding"!')
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      CALL read_field_from_file_2D(   filename, 'beta_sq' , grid, ice%beta_sq_a,   region_name)
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      CALL read_field_from_file_2D(   filename, 'phi_fric' , grid, ice%phi_fric_a, region_name)
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      CALL read_field_from_file_2D(   filename, 'alpha_sq', grid, ice%alpha_sq_a,  region_name)
      CALL read_field_from_file_2D(   filename, 'beta_sq' , grid, ice%beta_sq_a,   region_name)
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      CALL read_field_from_file_2D(   filename, 'alpha_sq', grid, ice%alpha_sq_a,  region_name)
      CALL read_field_from_file_2D(   filename, 'beta_sq' , grid, ice%beta_sq_a,   region_name)
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      CALL read_field_from_file_2D(   filename, 'phi_fric' , grid, ice%phi_fric_a, region_name)
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_BIV_bed_roughness_file

  SUBROUTINE read_BIV_target_velocity(    BIV_target, region_name)
    ! Read the target velocity NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_BIV_target_velocity), INTENT(INOUT) :: BIV_target
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_BIV_target_velocity'
    INTEGER                                       :: i,j
    LOGICAL                                       :: is_Rignot2011
    REAL(dp), DIMENSION(:    ), ALLOCATABLE       :: y_temp
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: u_temp, v_temp
    REAL(dp)                                      :: NaN

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    CALL read_field_from_xy_file_2D(           BIV_target%filename, 'u_surf', region_name, BIV_target%grid, BIV_target%u_surf,BIV_target%wu_surf)
    CALL read_field_from_xy_file_2D(           BIV_target%filename, 'v_surf', region_name, BIV_target%grid, BIV_target%v_surf,BIV_target%wv_surf)

    ! Exception: for some reason, the Rignot 2011 data has the y-axis reversed...
    is_Rignot2011 = .FALSE.
    DO i = 1, 256-36
      IF (BIV_target%filename(i:i+33) == 'antarctica_ice_velocity_450m_v2.nc') THEN
        is_Rignot2011 = .TRUE.
      END IF
    END DO
    IF (is_Rignot2011) THEN

      ! Allocate temporary memory for storing the upside-down data
      ALLOCATE( y_temp( BIV_target%grid%ny))
      ALLOCATE( u_temp( BIV_target%grid%nx, BIV_target%grid%ny))
      ALLOCATE( v_temp( BIV_target%grid%nx, BIV_target%grid%ny))

      ! Copy the upside-down data to temporary memory
      y_temp = BIV_target%grid%y
      u_temp = BIV_target%u_surf
      v_temp = BIV_target%v_surf

      ! Flip the data
      DO j = 1, BIV_target%grid%ny
        BIV_target%grid%y(   j) = y_temp(   BIV_target%grid%ny + 1 - j)
        BIV_target%u_surf( :,j) = u_temp( :,BIV_target%grid%ny + 1 - j)
        BIV_target%v_surf( :,j) = v_temp( :,BIV_target%grid%ny + 1 - j)
      END DO

      ! Deallocate temporary memory
      DEALLOCATE( y_temp)
      DEALLOCATE( u_temp)
      DEALLOCATE( v_temp)

      ! Set missing values to NaN
      NaN = 0._dp
      NaN = 0._dp / NaN
      DO i = 1, BIV_target%grid%nx
      DO j = 1, BIV_target%grid%ny
        IF (BIV_target%u_surf( i,j) == 0._dp .AND. BIV_target%v_surf( i,j) == 0._dp) THEN
          BIV_target%u_surf( i,j) = NaN
          BIV_target%v_surf( i,j) = NaN
        END IF
      END DO
      END DO

    END IF ! IF (is_Rignot2011) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_BIV_target_velocity

END MODULE basal_conditions_and_sliding_module
