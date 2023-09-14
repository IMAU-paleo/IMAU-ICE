MODULE climate_module

  ! Contains all the routines for calculating the climate forcing.

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_model_region, type_grid, type_grid_lonlat, type_ice_model, type_SMB_model, &
                                             type_reference_geometry, type_climate_model, type_climate_snapshot, &
                                             type_climate_model_PD_obs, type_climate_model_direct, type_climate_model_ISMIP_style
  USE netcdf_basic_module,             ONLY: inquire_var, inquire_var_multiple_options, field_name_options_wind_WE, field_name_options_wind_SN, &
                                             field_name_options_mask_ice, field_name_options_mask_ocean, field_name_options_Hs
  USE netcdf_input_module,             ONLY: read_field_from_file_2D_monthly, read_field_from_file_2D
  USE forcing_module,                  ONLY: forcing, get_insolation_at_time, update_CO2_at_model_time, get_summer_insolation_at_time
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             error_function, smooth_Gaussian_2D, smooth_Gaussian_3D, smooth_Shepard_2D, &
                                             map_glob_to_grid_2D, map_glob_to_grid_3D, transpose_dp_2D, transpose_dp_3D, &
                                             map_square_to_square_cons_2nd_order_2D, map_square_to_square_cons_2nd_order_3D
  USE derivatives_and_grids_module,    ONLY: ddx_a_to_a_2D, ddy_a_to_a_2D
  USE SMB_module,                      ONLY: run_SMB_model

  USE netcdf_debug_module,             ONLY: save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, save_variable_as_netcdf_int_3D, &
                                             save_variable_as_netcdf_dp_1D,  save_variable_as_netcdf_dp_2D,  save_variable_as_netcdf_dp_3D

  IMPLICIT NONE

CONTAINS

! == The routines that should be called from the main model
! =========================================================

  SUBROUTINE run_climate_model( region, time)
    ! Run the climate model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_climate_model == 'none') THEN
      ! No need to do anything.

    ELSEIF (C%choice_climate_model == 'idealised') THEN
      ! Assign some idealised temperature/precipitation

      CALL run_climate_model_idealised( region%grid, region%ice, region%climate, time)

    ELSEIF (C%choice_climate_model == 'PD_obs') THEN
      ! Keep the climate fixed to present-day observed conditions

      CALL run_climate_model_PD_obs( region%grid, region%ice, region%climate, region%name)

    ELSEIF (C%choice_climate_model == 'direct') THEN
      ! Use a directly prescribed global climate

      CALL run_climate_model_direct( region%grid, region%ice, region%climate, time, region%name)

    ELSEIF (C%choice_climate_model == 'direct_SMB') THEN
      ! Use a directly prescribed global SMB + 2-m air temperature

      CALL crash('choice_climate_model "'//TRIM( C%choice_climate_model)//'" is broken right now, must be fixed!')
      ! CALL run_climate_model_direct_SMB( region%grid, region%ice, region%climate, time)

    ELSEIF (C%choice_climate_model == 'matrix') THEN
      ! Use the warm/cold climate matrix (Berends et al., 2018)

      CALL run_climate_model_matrix( region%grid, region%ice, region%SMB, region%climate, region%refgeo_PD, region%name, region%time)

    ELSE
      CALL crash('unknown choice_climate_model"' // TRIM(C%choice_climate_model) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model
  SUBROUTINE initialise_climate_model( region)
    ! Initialise the climate model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Always allocate memory for the applied climate
    CALL allocate_shared_dp_3D( 12, region%grid%ny, region%grid%nx, region%climate%Q_TOA,   region%climate%wQ_TOA  )
    CALL allocate_shared_dp_3D( 12, region%grid%ny, region%grid%nx, region%climate%T2m,     region%climate%wT2m    )
    CALL allocate_shared_dp_3D( 12, region%grid%ny, region%grid%nx, region%climate%Precip,  region%climate%wPrecip )
    CALL allocate_shared_dp_3D( 12, region%grid%ny, region%grid%nx, region%climate%Wind_LR, region%climate%wWind_LR)
    CALL allocate_shared_dp_3D( 12, region%grid%ny, region%grid%nx, region%climate%Wind_DU, region%climate%wWind_DU)

    IF (par%master) WRITE (0,*) '  Initialising climate model "', TRIM(C%choice_climate_model), '"...'

    IF     (C%choice_climate_model == 'none') THEN
      ! No need to initialise anything other than the applied climate

    ELSEIF (C%choice_climate_model == 'idealised') THEN
      ! No need to initialise anything other than the applied climate

    ELSEIF (C%choice_climate_model == 'PD_obs') THEN
      ! Keep the climate fixed to present-day observed conditions

      CALL initialise_climate_model_PD_obs( region%grid, region%climate, region%name)

    ELSEIF (C%choice_climate_model == 'direct') THEN
      ! Use a directly prescribed climate

      CALL initialise_climate_model_direct( region%grid, region%climate%direct)

    ELSEIF (C%choice_climate_model == 'direct_SMB') THEN
      ! Use a directly prescribed climate

      CALL crash('choice_climate_model "'//TRIM( C%choice_climate_model)//'" is broken right now, must be fixed!')
!      CALL initialise_climate_model_direct_SMB( climate%direct)

    ELSEIF (C%choice_climate_model == 'matrix') THEN
      ! Use the warm/cold climate matrix (Berends et al., 2018)

      CALL initialise_climate_matrix( region%grid, region%climate, region%name, region%mask_noice)

    ELSE
      CALL crash('unknown choice_climate_model"' // TRIM(C%choice_climate_model) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model

! == Idealised climates
! =====================

  SUBROUTINE run_climate_model_idealised( grid, ice, climate, time)
    ! Run the climate model
    !
    ! Assign some idealised temperature/precipitation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_idealised_climate == 'EISMINT1_A' .OR. &
            C%choice_idealised_climate == 'EISMINT1_B' .OR. &
            C%choice_idealised_climate == 'EISMINT1_C' .OR. &
            C%choice_idealised_climate == 'EISMINT1_D' .OR. &
            C%choice_idealised_climate == 'EISMINT1_E' .OR. &
            C%choice_idealised_climate == 'EISMINT1_F') THEN
      CALL run_climate_model_idealised_EISMINT1( grid, ice, climate, time)
    ELSE
      CALL crash('unknown choice_idealised_climate"' // TRIM(C%choice_idealised_climate) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_idealised
  SUBROUTINE run_climate_model_idealised_EISMINT1( grid, ice, climate, time)
    ! The parameterised climates of the EISMINT1 idealised-geometry experiments

    USe parameters_module,           ONLY: T0, pi

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_idealised_EISMINT1'
    REAL(dp), PARAMETER                                :: lambda = -0.010_dp
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: dT_lapse, d, dT

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set precipitation to zero - SMB is parameterised anyway...
    climate%Precip( :,:,grid%i1:grid%i2) = 0._dp

    ! Surface temperature for fixed or moving margin experiments
    IF     (C%choice_idealised_climate == 'EISMINT1_A' .OR. &
            C%choice_idealised_climate == 'EISMINT1_B' .OR. &
            C%choice_idealised_climate == 'EISMINT1_C') THEN
      ! Moving margin

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny

        dT_lapse = ice%Hs_a(j,i) * lambda

        DO m = 1, 12
          climate%T2m(m,j,i) = 270._dp + dT_lapse
        END DO

      END DO
      END DO

    ELSEIF (C%choice_idealised_climate == 'EISMINT1_D' .OR. &
            C%choice_idealised_climate == 'EISMINT1_E' .OR. &
            C%choice_idealised_climate == 'EISMINT1_F') THEN
      ! Fixed margin

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny

        d = MAX( ABS(grid%x(i)/1000._dp), ABS(grid%y(j)/1000._dp))

        DO m = 1, 12
          climate%T2m(m,j,i) = 239._dp + (8.0E-08_dp * d**3)
        END DO

      END DO
      END DO

    ELSE
      CALL crash('unknown choice_idealised_climate"' // TRIM(C%choice_idealised_climate) // '"!')
    END IF
    CALL sync

    ! Glacial cycles
    IF     (C%choice_idealised_climate == 'EISMINT1_B' .OR. &
            C%choice_idealised_climate == 'EISMINT1_E') THEN
      IF (time > 0._dp) THEN
        dT = 10._dp * SIN(2._dp * pi * time / 20000._dp)
        climate%T2m( :,:,grid%i1:grid%i2) = climate%T2m( :,:,grid%i1:grid%i2) + dT
      END IF
    ELSEIF (C%choice_idealised_climate == 'EISMINT1_C' .OR. &
            C%choice_idealised_climate == 'EISMINT1_F') THEN
      IF (time > 0._dp) THEN
        dT = 10._dp * SIN(2._dp * pi * time / 40000._dp)
        climate%T2m( :,:,grid%i1:grid%i2) = climate%T2m( :,:,grid%i1:grid%i2) + dT
      END IF
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_idealised_EISMINT1

! == Static present-day observed climate
! ======================================

  SUBROUTINE run_climate_model_PD_obs( grid, ice, climate, region_name)
    ! Run the climate model
    !
    ! Keep the climate fixed to present-day observed conditions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_PD_obs'
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: lambda, dT_lapse
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  T2m_corr,  Precip_corr
    INTEGER                                            :: wT2m_corr, wPrecip_corr

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, T2m_corr   , wT2m_corr   )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, Precip_corr, wPrecip_corr)

    ! Keep the climate fixed to present-day observed conditions
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      climate%T2m(     m,j,i) = climate%PD_obs%snapshot%T2m(     m,j,i)
      climate%Precip(  m,j,i) = climate%PD_obs%snapshot%Precip(  m,j,i)
      climate%Wind_LR( m,j,i) = climate%PD_obs%snapshot%Wind_LR( m,j,i)
      climate%Wind_DU( m,j,i) = climate%PD_obs%snapshot%Wind_DU( m,j,i)
    END DO
    END DO
    END DO
    CALL sync

    ! Define lapse rate as a negative number
    lambda = -1._dp * ABS( C%constant_lapserate)

    ! Calculate lapse-rate-corrected temperature
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      dT_lapse = (ice%Hs_a( j,i) - climate%PD_obs%snapshot%Hs( j,i)) * lambda
      T2m_corr( m,j,i) = climate%T2m( m,j,i) + dT_lapse
    END DO
    END DO
    END DO

    ! Downscale precipitation from the coarse-resolution reference
    ! orography to the fine-resolution ice-model orography
    IF (region_name == 'NAM' .OR. region_name == 'EAS' .OR. region_name == 'PAT') THEN
      ! Use the Roe&Lindzen precipitation model to do this; Berends et al., 2018, Eqs. A3-A7

      CALL adapt_precip_Roe( grid, climate%PD_obs%snapshot%Hs, climate%PD_obs%snapshot%T2m, &
        climate%PD_obs%snapshot%Wind_LR, climate%PD_obs%snapshot%Wind_DU, climate%PD_obs%snapshot%Precip, &
        ice%Hs_a, T2m_corr, climate%Wind_LR, climate%Wind_DU, Precip_corr)

    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use a simpler temperature-based correction; Berends et al., 2018, Eq. 14

      CALL adapt_precip_CC( grid, ice%Hs_a, climate%PD_obs%snapshot%Hs, climate%PD_obs%snapshot%T2m, &
                            climate%PD_obs%snapshot%Precip, Precip_corr, region_name)

    END IF

    ! Store corrected climate
    climate%T2m(    :,:,grid%i1:grid%i2) = T2m_corr(    :,:,grid%i1:grid%i2)
    climate%Precip( :,:,grid%i1:grid%i2) = Precip_corr( :,:,grid%i1:grid%i2)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wT2m_corr   )
    CALL deallocate_shared( wPrecip_corr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_PD_obs
  SUBROUTINE initialise_climate_model_PD_obs( grid, climate,region_name)
    ! Initialise the climate model
    !
    ! Keep the climate fixed to present-day observed conditions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_model_PD_obs'
    LOGICAL                                            :: found_winds

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    CALL allocate_climate_snapshot( grid, climate%PD_obs%snapshot, 'PD_obs')

    ! Read climate from file
    CALL read_climate_snapshot( C%filename_PD_obs_climate, grid, climate%PD_obs%snapshot, found_winds, region_name)

    ! Safety
    IF (.NOT. found_winds) CALL crash('couldnt find wind fields for PD climate in file ' // TRIM( C%filename_PD_obs_climate))

    ! Initialise insolation at present-day (needed for the IMAU-ITM SMB model)
    CALL get_insolation_at_time( grid, 0.0_dp, climate%PD_obs%snapshot%Q_TOA)

    ! Copy all values to the applied climate
    climate%Q_TOA(   :,:,grid%i1:grid%i2) = climate%PD_obs%snapshot%Q_TOA(   :,:,grid%i1:grid%i2)
    climate%T2m(     :,:,grid%i1:grid%i2) = climate%PD_obs%snapshot%T2m(     :,:,grid%i1:grid%i2)
    climate%Precip(  :,:,grid%i1:grid%i2) = climate%PD_obs%snapshot%Precip(  :,:,grid%i1:grid%i2)
    climate%Wind_LR( :,:,grid%i1:grid%i2) = climate%PD_obs%snapshot%Wind_LR( :,:,grid%i1:grid%i2)
    climate%Wind_DU( :,:,grid%i1:grid%i2) = climate%PD_obs%snapshot%Wind_DU( :,:,grid%i1:grid%i2)
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_PD_obs

! == Directly prescribed time-varying climate
! ===========================================

  SUBROUTINE run_climate_model_direct( grid, ice, climate, time, region_name)
    ! Run the climate model
    !
    ! Use a directly prescribed climate

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    REAL(dp),                            INTENT(IN)    :: time
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_direct'
    REAL(dp)                                           :: wt0, wt1
    INTEGER                                            :: i,j,m
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  Hs
    INTEGER                                            :: wHs
    REAL(dp)                                           :: lambda, dT_lapse
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  T2m_corr,  Precip_corr
    INTEGER                                            :: wT2m_corr, wPrecip_corr

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct') THEN
      CALL crash('choice_climate_model should be "direct"!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, Hs, wHs)

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF files

    IF (time < climate%direct%t0 .OR. time > climate%direct%t1) THEN
      ! Find and read the two global time frames

      CALL update_direct_climate_timeframes( grid, climate, time, region_name)

    END IF ! IF (time < climate%direct%t0 .OR. time > climate%direct%t1) THEN

    ! Interpolate the two timeframes in time
    wt0 = (climate%direct%t1 - time) / (climate%direct%t1 - climate%direct%t0)
    wt1 = 1._dp - wt0

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      Hs(                  j,i) = (wt0 * climate%direct%timeframe0%Hs(        j,i)) + (wt1 * climate%direct%timeframe1%Hs(        j,i))
      climate%T2m(       m,j,i) = (wt0 * climate%direct%timeframe0%T2m(     m,j,i)) + (wt1 * climate%direct%timeframe1%T2m(     m,j,i))
      climate%Precip(    m,j,i) = (wt0 * climate%direct%timeframe0%Precip(  m,j,i)) + (wt1 * climate%direct%timeframe1%Precip(  m,j,i))
      climate%Wind_LR(   m,j,i) = (wt0 * climate%direct%timeframe0%Wind_LR( m,j,i)) + (wt1 * climate%direct%timeframe1%Wind_LR( m,j,i))
      climate%Wind_DU(   m,j,i) = (wt0 * climate%direct%timeframe0%Wind_DU( m,j,i)) + (wt1 * climate%direct%timeframe1%Wind_DU( m,j,i))
    END DO
    END DO
    END DO
    CALL sync

    ! Apply geometry corrections
    IF (C%do_direct_climate_geo_corr) THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, T2m_corr   , wT2m_corr   )
      CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, Precip_corr, wPrecip_corr)

      ! Define lapse rate as a negative number
      lambda = -1._dp * ABS( C%constant_lapserate)

      ! Calculate lapse-rate-corrected temperature
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
      DO m = 1, 12
        dT_lapse = (ice%Hs_a( j,i) - Hs( j,i)) * lambda
        T2m_corr( m,j,i) = climate%T2m( m,j,i) + dT_lapse
      END DO
      END DO
      END DO

      ! Downscale precipitation from the coarse-resolution reference
      ! orography to the fine-resolution ice-model orography
      IF (region_name == 'NAM' .OR. region_name == 'EAS' .OR. region_name == 'PAT') THEN
        ! Use the Roe&Lindzen precipitation model to do this; Berends et al., 2018, Eqs. A3-A7

        CALL adapt_precip_Roe( grid, Hs      , climate%T2m, climate%Wind_LR, climate%Wind_DU, climate%Precip, &
                                     ice%Hs_a, T2m_corr   , climate%Wind_LR, climate%Wind_DU, Precip_corr)

      ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
        ! Use a simpler temperature-based correction; Berends et al., 2018, Eq. 14

        CALL adapt_precip_CC( grid, ice%Hs_a, Hs, climate%T2m, climate%Precip, Precip_corr, region_name)

      END IF

      ! Store corrected climate
      climate%T2m(    :,:,grid%i1:grid%i2) = T2m_corr(    :,:,grid%i1:grid%i2)
      climate%Precip( :,:,grid%i1:grid%i2) = Precip_corr( :,:,grid%i1:grid%i2)
      CALL sync

      ! Clean up after yourself
      CALL deallocate_shared( wT2m_corr   )
      CALL deallocate_shared( wPrecip_corr)

    END IF ! IF (C%do_direct_climate_geo_corr) THEN

    ! Update insolation
    CALL get_insolation_at_time( grid, time, climate%Q_TOA)

    ! Clean up after yourself
    CALL deallocate_shared( wHs)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_direct
  SUBROUTINE update_direct_climate_timeframes( grid, climate, time, region_name)
    ! Update the two timeframes of the direct climate forcing

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    REAL(dp),                            INTENT(IN)    :: time
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_direct_climate_timeframes'
    INTEGER                                            :: year0, year1
    CHARACTER(LEN=4)                                   :: year0str, year1str
    CHARACTER(LEN=256)                                 :: filename0, filename1
    CHARACTER(LEN=256)                                 :: direct_climate_foldername, direct_climate_basefilename
    LOGICAL                                            :: found_winds

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Assume we have a single NetCDF file for each calendar year.
    ! Determine which two we need to read

    year0 = FLOOR(        time)
    year1 = MAX( CEILING( time), year0+1)

    WRITE( year0str,'(I4)') year0
    WRITE( year1str,'(I4)') year1

    ! Update timestamps
    climate%direct%t0 = REAL( year0,dp)
    climate%direct%t1 = REAL( year1,dp)

    ! Get filenames

    IF     (region_name == 'NAM') THEN
      direct_climate_foldername   = C%direct_climate_foldername_NAM
      direct_climate_basefilename = C%direct_climate_basefilename_NAM
    ELSEIF (region_name == 'EAS') THEN
      direct_climate_foldername   = C%direct_climate_foldername_EAS
      direct_climate_basefilename = C%direct_climate_basefilename_EAS
    ELSEIF (region_name == 'GRL') THEN
      direct_climate_foldername   = C%direct_climate_foldername_GRL
      direct_climate_basefilename = C%direct_climate_basefilename_GRL
    ELSEIF (region_name == 'ANT') THEN
      direct_climate_foldername   = C%direct_climate_foldername_ANT
      direct_climate_basefilename = C%direct_climate_basefilename_ANT
    ELSE
      CALL crash('unknown model region "' // TRIM( region_name) // '"!')
    END IF

    filename0 = TRIM( direct_climate_foldername) // '/' // TRIM( direct_climate_basefilename) // year0str // '.nc'
    filename1 = TRIM( direct_climate_foldername) // '/' // TRIM( direct_climate_basefilename) // year1str // '.nc'

    ! Prevent weird screen message output when using time display
    IF (par%master .AND. C%do_time_display .AND. time > C%start_time_of_run) WRITE(0,*) ''

    ! Read timeframes from files
    CALL read_climate_snapshot( filename0, grid, climate%direct%timeframe0, found_winds, region_name)
    ! IF (.NOT. found_winds) CALL crash('couldnt find wind fields for direct prescribed climate in file ' // TRIM( filename0)) !CvC
    CALL read_climate_snapshot( filename1, grid, climate%direct%timeframe1, found_winds, region_name)
    ! IF (.NOT. found_winds) CALL crash('couldnt find wind fields for direct prescribed climate in file ' // TRIM( filename1))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_direct_climate_timeframes
  SUBROUTINE initialise_climate_model_direct( grid, direct)
    ! Initialise the climate model
    !
    ! Use a directly prescribed climate

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_model_direct),     INTENT(INOUT) :: direct

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_model_direct'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct') THEN
      CALL crash('choice_climate_model should be "direct"!')
    END IF

    ! Set timeframes timestamps to impossible values, so that the first call to run_climate_model_direct
    ! is guaranteed to first read two new timeframes from the NetCDF file
    direct%t0 = C%start_time_of_run - 100._dp
    direct%t1 = C%start_time_of_run - 90._dp

    ! Allocate shared memory
    CALL allocate_climate_snapshot( grid, direct%timeframe0, 'climate_direct_timeframe0')
    CALL allocate_climate_snapshot( grid, direct%timeframe1, 'climate_direct_timeframe1')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_direct

! == Climate matrix
! ===========================

  ! Climate matrix with warm + cold snapshots, forced with CO2 (from record or from inverse routine) from Berends et al., 2018
  ! Generalised for different timeframes, L.B. Stap (2021)
  SUBROUTINE run_climate_model_matrix( grid, ice, SMB, climate, refgeo_PD, region_name, time)
    ! Use CO2 (either prescribed or inversely modelled) to force the 2-snapshot (PI-LGM) climate matrix (Berends et al., 2018)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_matrix'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Update insolation at model time
    CALL get_insolation_at_time( grid, time, climate%Q_TOA)

    ! Update 65N insolation at model time (perhaps replace with caloric summer)
    IF (C%do_combine_CO2_and_insolation) THEN
      CALL get_summer_insolation_at_time( time)
    END IF
  
   ! Update CO2 at model time
    IF (C%choice_forcing_method == 'CO2_direct') THEN
      CALL update_CO2_at_model_time( time)
    END IF
    
    ! Use the (CO2 + absorbed insolation)-based interpolation scheme for temperature
    CALL run_climate_model_matrix_temperature( grid, ice, SMB, climate, region_name)

    ! Use the (CO2 + ice-sheet geometry)-based interpolation scheme for precipitation
    CALL run_climate_model_matrix_precipitation( grid, ice, climate, refgeo_PD, region_name)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_matrix
  SUBROUTINE run_climate_model_matrix_temperature( grid, ice, SMB, climate, region_name)
    ! The (CO2 + absorbed insolation)-based matrix interpolation for temperature, from Berends et al. (2018)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_matrix_temperature'
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: CO2, w_CO2, w_QTOA
    REAL(dp)                                           :: w_ins_mean, w_ins_amplitude
    REAL(dp), DIMENSION(:,:  ), POINTER                :: w_ins_smooth  
    INTEGER                                            :: ww_ins_smooth
    REAL(dp)                                           :: w_ins_av
    REAL(dp), DIMENSION(:,:,:), POINTER                :: T_ref_GCM
    REAL(dp), DIMENSION(:,:  ), POINTER                :: Hs_GCM, lambda_GCM
    INTEGER                                            :: wT_ref_GCM, wHs_GCM, wlambda_GCM
    REAL(dp)                                           :: w_from_w_ins_av
    REAL(dp), PARAMETER                                :: w_cutoff = 0.5_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]
    REAL(dp), PARAMETER                                :: P_offset = 0.008_dp       ! Normalisation term in precipitation anomaly to avoid divide-by-nearly-zero

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_ins_smooth, ww_ins_smooth  )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, T_ref_GCM,    wT_ref_GCM     )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, Hs_GCM,       wHs_GCM    )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, lambda_GCM,   wlambda_GCM)

    ! Find CO2 interpolation weight (use either prescribed or modelled CO2)
    ! =====================================================================

    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      CO2 = forcing%CO2_obs
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CO2 = forcing%CO2_mod
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      CO2 = 0._dp
      CALL crash('must only be called with the correct forcing method, check your code!')
    ELSE
      CO2 = 0._dp
      CALL crash('unknown choice_forcing_method"' // TRIM(C%choice_forcing_method) // '"!')
    END IF

    w_CO2 = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (CO2 - C%matrix_low_CO2_level) / (C%matrix_high_CO2_level - C%matrix_low_CO2_level) ))   ! Berends et al., 2018 - Eq. 1


   ! Combine CO2 and insolation to obtain external forcing
   ! ======================================================================

    IF (C%do_combine_CO2_and_insolation) THEN
      ! Combine CO2 and insolation in the climate matrix method
      w_ins_mean      = 0._dp
      w_ins_amplitude = 0._dp

      ! Select the corresponding w_ins_mean and amplitude
      IF (region_name == 'NAM') THEN
        w_ins_mean      = C%insolation_weigth_mean_NAM
        w_ins_amplitude = C%insolation_weigth_amplitude_NAM
      ELSEIF (region_name == 'EAS') THEN
        w_ins_mean      = C%insolation_weigth_mean_EAS
        w_ins_amplitude = C%insolation_weigth_amplitude_EAS
      ELSEIF (region_name == 'GRL') THEN
        w_ins_mean      = C%insolation_weigth_mean_GRL
        w_ins_amplitude = C%insolation_weigth_amplitude_GRL
      ELSEIF (region_name == 'ANT') THEN
        w_ins_mean      = C%insolation_weigth_mean_ANT
        w_ins_amplitude = C%insolation_weigth_amplitude_ANT
      ELSE
        CALL crash('region_name "'//TRIM(region_name)//'" not found!')
      END IF
      
      ! Calculate the insolation weight based on 65N summer insolation
      w_QTOA = 0._dp

      ! Select the correct insolation for the Hemisphere
      IF (region_name == 'NAM' .OR. region_name == 'EAS' .OR. region_name == 'GRL') THEN
        w_QTOA = (forcing%Q_TOA_JJA_65N - w_ins_mean) / w_ins_amplitude
      ELSEIF (region_name == 'ANT') THEN
        w_QTOA = (forcing%Q_TOA_DJF_80S - w_ins_mean) / w_ins_amplitude
      ELSE
        CALL crash('region_name "'//TRIM(region_name)//'" not found!')
      END IF
      CALL sync

      ! Combine CO2 and insolation
      climate%matrix%w_EXT = w_CO2 + w_QTOA
      
    ELSE ! C%do_combine_CO2_and_insolation
      ! Only use CO2 in the climate matrix method
      climate%matrix%w_EXT = w_CO2
      
    END IF ! C%do_combine_CO2_and_insolation
    CALL sync

    ! Make sure w_EXT is not too large or small
    climate%matrix%w_EXT = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, climate%matrix%w_EXT))

    ! Find the interpolation weights based on absorbed insolation
    ! ===========================================================

    ! Calculate modelled absorbed insolation
    climate%matrix%I_abs( :,grid%i1:grid%i2) = 0._dp

    ! Calculate I_abs (annual total absorbed insolation) based on monthly albedo and insolation
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
        DO m = 1, 12
          climate%matrix%I_abs( j,i) = climate%matrix%I_abs( j,i) + climate%Q_TOA( m,j,i) * (1._dp - SMB%Albedo( m,j,i))  ! Berends et al., 2018 - Eq. 2
        END DO
    END DO
    END DO
    CALL sync

    ! Calculate the domain wide average insolation
    w_ins_av      = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (SUM( climate%matrix%I_abs         )      - SUM( climate%matrix%GCM_cold%I_abs)     ) / &
                                                           (SUM( climate%matrix%GCM_warm%I_abs)      - SUM( climate%matrix%GCM_cold%I_abs)     ) ))

    ! Account for local changes in albedo and insolation
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (climate%matrix%GCM_warm%mask_ocean( j,i) == 1 .AND. climate%matrix%GCM_cold%mask_ocean( j,i) == 1 ) THEN
        ! Ocean in both cases.
        !  Absorbed insolation in the ocean only depends on insolation, which can cause huge temporal and spatial
        !  gradients in w_ins. Therefore, use the domain-wide change in absorbed insolation instead
        climate%matrix%w_ins_T( j,i) = w_ins_av
      ELSE
        ! Land in at least one case, local I_abs can be used
        climate%matrix%w_ins_T( j,i) = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (     climate%matrix%I_abs(          j,i) -      climate%matrix%GCM_cold%I_abs( j,i)) / &  ! Berends et al., 2018 - Eq. 3
                                                                              (     climate%matrix%GCM_warm%I_abs( j,i) -      climate%matrix%GCM_cold%I_abs( j,i)) ))
        
        IF (climate%matrix%GCM_warm%mask_ice( j,i) == 0 .AND. climate%matrix%GCM_cold%mask_ice( j,i) == 0) THEN
          ! Small albedo differences
          !  If albedo in the warm and cold snapshot are very close together, I_abs will mostly vary due to insolation. This
          !  is a problem. Insolation may be substantially higher or lower compared to BOTH snapshots. Therefore,
          !  in these regions (e.g., dry tundra, warm regions) the w_ins_T will almost exclusively be too high or too low,
          !  which can affect local temperatures. In the tundra regions this may make it very difficult to melt ice.
          !
          !  Therefore, instead, we let w_ins_T be partly determined by w_ins_av
          climate%matrix%w_ins_T( j,i) = 0.5_dp * climate%matrix%w_ins_T( j,i) + 0.5_dp * w_ins_av 
        
        END IF
      END IF
    END DO
    END DO
    CALL sync


    ! Smooth the weighting field
    w_ins_smooth( :,grid%i1:grid%i2) = climate%matrix%w_ins_T( :,grid%i1:grid%i2)
    CALL smooth_Gaussian_2D( grid, w_ins_smooth, 200000._dp)

    ! Combine unsmoothed, smoothed, and regional average weighting fields (Berends et al., 2018, Eq. 4)
    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      climate%matrix%w_ice_T( :,grid%i1:grid%i2) = (1._dp * climate%matrix%w_ins_T(        :,grid%i1:grid%i2) + &
                                   3._dp * w_ins_smooth( :,grid%i1:grid%i2) + &
                                   3._dp * w_ins_av) / 7._dp
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      climate%matrix%w_ice_T( :,grid%i1:grid%i2) = (1._dp * w_ins_smooth( :,grid%i1:grid%i2) + &
                                   6._dp * w_ins_av) / 7._dp
    END IF


    ! Combine interpolation weights from absorbed insolation and CO2 into the final weights fields
    ! Berends et al., 2018 - Eqs. 5, 9 with weights 0.5 for NAM & EAS, and 0.75 for ANT
    ! Generalised: "switch" between matrix method and glacial index method by altering C%climate_matrix_CO2vsice_<region>
    IF         (region_name == 'NAM') THEN
      climate%matrix%w_tot_T( :,grid%i1:grid%i2) = (C%climate_matrix_CO2vsice_NAM * climate%matrix%w_EXT) + ((1._dp -  C%climate_matrix_CO2vsice_NAM) * climate%matrix%w_ice_T( :,grid%i1:grid%i2))
    ELSEIF     (region_name == 'EAS') THEN
      climate%matrix%w_tot_T( :,grid%i1:grid%i2) = (C%climate_matrix_CO2vsice_EAS * climate%matrix%w_EXT) + ((1._dp -  C%climate_matrix_CO2vsice_EAS) * climate%matrix%w_ice_T( :,grid%i1:grid%i2))
    ELSEIF     (region_name == 'GRL') THEN
      climate%matrix%w_tot_T( :,grid%i1:grid%i2) = (C%climate_matrix_CO2vsice_GRL * climate%matrix%w_EXT) + ((1._dp -  C%climate_matrix_CO2vsice_GRL) * climate%matrix%w_ice_T( :,grid%i1:grid%i2))
    ELSEIF     (region_name == 'ANT') THEN
      climate%matrix%w_tot_T( :,grid%i1:grid%i2) = (C%climate_matrix_CO2vsice_ANT * climate%matrix%w_EXT) + ((1._dp -  C%climate_matrix_CO2vsice_ANT) * climate%matrix%w_ice_T( :,grid%i1:grid%i2))
    END IF

    ! Interpolate between the GCM snapshots
    ! =====================================

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Find matrix-interpolated orography, lapse rate, and temperature
      Hs_GCM(       j,i) = (climate%matrix%w_tot_T( j,i) * climate%matrix%GCM_warm%Hs(       j,i)) + ((1._dp - climate%matrix%w_tot_T( j,i)) * climate%matrix%GCM_cold%Hs(       j,i))  ! Berends et al., 2018 - Eq. 8
      lambda_GCM(   j,i) = (climate%matrix%w_tot_T( j,i) * climate%matrix%GCM_warm%lambda(   j,i)) + ((1._dp - climate%matrix%w_tot_T( j,i)) * climate%matrix%GCM_cold%lambda(   j,i))  ! Not listed in the article, shame on me!
      T_ref_GCM(  :,j,i) = (climate%matrix%w_tot_T( j,i) * climate%matrix%GCM_warm%T2m(    :,j,i)) + ((1._dp - climate%matrix%w_tot_T( j,i)) * climate%matrix%GCM_cold%T2m(    :,j,i))  ! Berends et al., 2018 - Eq. 6

      ! Adapt temperature to model orography using matrix-derived lapse-rate
      DO m = 1, 12
        climate%T2m( m,j,i) = T_ref_GCM( m,j,i) - lambda_GCM( j,i) * (ice%Hs_a( j,i) - Hs_GCM( j,i))  ! Berends et al., 2018 - Eq. 11
      END DO

    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( ww_ins_smooth)
    CALL deallocate_shared( wT_ref_GCM)
    CALL deallocate_shared( wHs_GCM)
    CALL deallocate_shared( wlambda_GCM)
    
    ! Safety checks
    CALL check_safety_temperature( climate%T2m)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_matrix_temperature
  SUBROUTINE run_climate_model_matrix_precipitation( grid, ice, climate, refgeo_PD, region_name)
    ! The (CO2 + ice geometry)-based matrix interpolation for precipitation, from Berends et al. (2018)
    ! For NAM and EAS, this is based on local ice geometry and uses the Roe&Lindzen precipitation model for downscaling.
    ! For GRL and ANT, this is based on total ice volume,  and uses the simple CC   precipitation model for downscaling.
    ! The rationale for this difference is that glacial-interglacial differences in ice geometry are much more
    ! dramatic in NAM and EAS than they are in GRL and ANT.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_matrix_precipitation'
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: CO2
    REAL(dp), DIMENSION(:,:,:), POINTER                :: T_ref_GCM, T_ref_ice, P_ref_GCM
    REAL(dp), DIMENSION(:,:  ), POINTER                :: Hs_GCM, lambda_GCM
    REAL(dp), DIMENSION(:,:  ), POINTER                :: dHs_snapshots, dHs_ice
    INTEGER                                            :: wT_ref_GCM, wT_ref_ice, wlambda_GCM, wP_ref_GCM, wdHs_ice, wdHs_snapshots, wHs_GCM   
    REAL(dp), PARAMETER                                :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, T_ref_GCM,      wT_ref_GCM      )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, T_ref_ice,      wT_ref_ice      )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, P_ref_GCM,      wP_ref_GCM      )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, Hs_GCM,         wHs_GCM         )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, lambda_GCM,     wlambda_GCM     )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, dHs_ice,        wdHs_ice        )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, dHs_snapshots,  wdHs_snapshots  )

    ! Calculate interpolation weights based on ice geometry
    ! =====================================================

    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      CO2 = forcing%CO2_obs
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CO2 = forcing%CO2_mod
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      CO2 = 0._dp
      CALL crash('must only be called with the correct forcing method, check your code!')
    ELSE
      CO2 = 0._dp
      CALL crash('unknown choice_forcing_method"' // TRIM(C%choice_forcing_method) // '"!')
    END IF

    ! First calculate the total ice volume term (second term in the equation)
    ! - The climate matrix method interpolates precipitation with respect to the change in 
    !   topography between the warm and cold snapshot. We only consider this change in
    !   regions with ice-sheets. Therefore, we ignore topography changes due to e.g., sea level.
    
    dHs_snapshots( :,grid%i1:grid%i2) = 0._dp
    dHs_ice(       :,grid%i1:grid%i2) = 0._dp

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Topography change between the cold and/or warm periods
      IF ((climate%matrix%GCM_warm%mask_ice( j,i) == 1) .OR. (climate%matrix%GCM_cold%mask_ice( j,i) == 1)) THEN
        ! Check if the GCM_cold topography is higher than GCM_warm
        IF ((climate%matrix%GCM_cold%Hs( j,i) + 10._dp) > climate%matrix%GCM_warm%Hs( j,i)) THEN
          ! Calculate the difference between the modelled Hs and snapshot Hs
          dHs_snapshots( j,i) = climate%matrix%GCM_cold%Hs( j,i) - climate%matrix%GCM_warm%Hs( j,i)
        END IF
      END IF

      ! Ice sheet model topography change 
      IF (ice%mask_ice_a( j,i) == 1) THEN
        dHs_ice(       j,i) = ice%Hs_a( j,i)  - refgeo_PD%Hs( j,i)
      END IF

    END DO
    END DO
    CALL sync

    ! Calculate w_tot (the domain-wide temperature difference)
    IF (par%master) THEN
      climate%matrix%w_tot_P  = MAX(-w_cutoff, MIN(1._dp + w_cutoff, (SUM( dHs_ice) / SUM(dHs_snapshots))))
    END IF
    CALL sync

    ! Combine total and local topography
    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      ! Combine total + local topography; Adapted from Berends et al., 2018, Eq. 12
      ! Then the local ice thickness term
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
      
        ! Determine the local effect  due to precipitation for all four possible ice/noice combinations

        IF (.NOT.((climate%matrix%GCM_cold%Hs( j,i) + 10._dp) > climate%matrix%GCM_warm%Hs( j,i))) THEN
          ! GCM_cold topography MUST be higher than GCM_warm, otherwise the interpolation
          ! cannot work. Also, dividing by 0 will cause trouble. Therefore we use only total ice volume.
          climate%matrix%w_cold_P( j,i) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, climate%matrix%w_tot_P ))
          climate%matrix%w_warm_P( j,i) = 1._dp - climate%matrix%w_cold_P( j,i)
            
        ELSEIF ((    climate%matrix%GCM_warm%mask_ice( j,i) == 1) .AND. (climate%matrix%GCM_cold%mask_ice( j,i) == 1)) THEN
          ! Ice in both GCM states.  Linear inter- / extrapolation
          climate%matrix%w_cold_P( j,i) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, &
                                    (dHs_ice( j,i) / dHs_snapshots( j,i)) * climate%matrix%w_tot_P))
          climate%matrix%w_warm_P( j,i)  = 1._dp - climate%matrix%w_cold_P( j,i)  
          
        ELSEIF ((climate%matrix%GCM_warm%mask_ice( j,i) == 0) .AND. (climate%matrix%GCM_cold%mask_ice( j,i) == 0)) THEN
          ! No ice in any GCM state. Use only total ice volume.
          climate%matrix%w_cold_P( j,i) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, climate%matrix%w_tot_P ))
          climate%matrix%w_warm_P( j,i) = 1._dp - climate%matrix%w_cold_P( j,i)
            
        ELSEIF ((climate%matrix%GCM_warm%mask_ice( j,i) == 0) .AND. (climate%matrix%GCM_cold%mask_ice( j,i) == 1)) THEN
          ! No ice in warm climate, ice in cold climate. Linear inter- / extrapolation.
          climate%matrix%w_cold_P( j,i) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, &
                                 (dHs_ice( j,i) / dHs_snapshots( j,i)) * climate%matrix%w_tot_P))

          climate%matrix%w_warm_P( j,i)  = 1._dp - climate%matrix%w_cold_P( j,i)
            
        ELSEIF ((  climate%matrix%GCM_warm%mask_ice( j,i) == 1) .AND. (climate%matrix%GCM_cold%mask_ice( j,i) == 0)) THEN
          ! Ice in warm climate, no ice in cold climate. Use total ice volume.
          climate%matrix%w_cold_P( j,i) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, climate%matrix%w_tot_P ))
          climate%matrix%w_warm_P( j,i) = 1._dp - climate%matrix%w_cold_P( j,i)
        ELSE
          CALL crash('This should not be possible!')  
        END IF

      END DO
      END DO
      CALL sync

      ! Combine the local and domain-wide topography influence
      climate%matrix%w_cold_P( :,grid%i1:grid%i2) = climate%matrix%w_cold_P( :,grid%i1:grid%i2) * climate%matrix%w_tot_P

      ! Smooth the weighting field
      CALL smooth_Gaussian_2D( grid, climate%matrix%w_cold_P, 200000._dp)

      climate%matrix%w_warm_P( :,grid%i1:grid%i2) = 1._dp - climate%matrix%w_cold_P( :,grid%i1:grid%i2)

    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use only total ice volume and CO2; Berends et al., 2018, Eq. 13

      climate%matrix%w_cold_P( :,grid%i1:grid%i2) = climate%matrix%w_tot_P
      climate%matrix%w_warm_P( :,grid%i1:grid%i2) = 1._dp - climate%matrix%w_cold_P( :,grid%i1:grid%i2)

    END IF

    ! Glacial index method
    IF (C%switch_glacial_index_precip) THEN ! If a glacial index is used for the precipitation forcing, it will only depend on CO2
      climate%matrix%w_tot_P     = 1._dp - (MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (CO2 - C%matrix_low_CO2_level) / (C%matrix_high_CO2_level - C%matrix_low_CO2_level) )) )
      climate%matrix%w_cold_P( :,grid%i1:grid%i2) = climate%matrix%w_tot_P
      climate%matrix%w_warm_P( :,grid%i1:grid%i2) = 1._dp - climate%matrix%w_cold_P( :,grid%i1:grid%i2)
    END IF

    !  Limit w_cold and w_warm to prevent too much extrapolation
    IF (par%master) THEN
      climate%matrix%w_cold_P = MAX(-w_cutoff, MIN(1._dp + w_cutoff, climate%matrix%w_cold_P ))
      climate%matrix%w_warm_P = MAX(-w_cutoff, MIN(1._dp + w_cutoff, climate%matrix%w_warm_P ))
    END IF
    CALL sync

    ! Interpolate the GCM snapshots
    ! =============================

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      T_ref_GCM(  :,j,i) =      (climate%matrix%w_warm_P( j,i) *      climate%matrix%GCM_warm%T2m(    :,j,i))  + (climate%matrix%w_cold_P( j,i) *     climate%matrix%GCM_cold%T2m(    :,j,i))   ! Berends et al., 2018 - Eq. 6
      P_ref_GCM(  :,j,i) = EXP( (climate%matrix%w_warm_P( j,i) *  LOG(climate%matrix%GCM_warm%Precip( :,j,i))) + (climate%matrix%w_cold_P( j,i) * LOG(climate%matrix%GCM_cold%Precip( :,j,i)))) ! Berends et al., 2018 - Eq. 7
      Hs_GCM(       j,i) =      (climate%matrix%w_warm_P( j,i) *      climate%matrix%GCM_warm%Hs(       j,i))  + (climate%matrix%w_cold_P( j,i) *     climate%matrix%GCM_cold%Hs(       j,i))   ! Berends et al., 2018 - Eq. 8
      lambda_GCM(   j,i) =      (climate%matrix%w_warm_P( j,i) *      climate%matrix%GCM_warm%lambda(   j,i))  + (climate%matrix%w_cold_P( j,i) *     climate%matrix%GCM_cold%lambda(   j,i))   ! Not listed in the article, shame on me!

      ! Adapt reference temperature to model orography using matrix-derived lapse-rate
      DO m = 1, 12
        T_ref_ice( m,j,i) = T_ref_GCM( m,j,i) - lambda_GCM( j,i) * (ice%Hs_a( j,i) - Hs_GCM( j,i)) 
      END DO

    END DO
    END DO
    CALL sync

    ! If there is also an interpolation between wind
    IF (C%do_climate_matrix_wind) THEN
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      DO m = 1, 12
        climate%Wind_LR( m,j,i) =  (climate%matrix%w_warm_P( j,i) *climate%matrix%GCM_warm%Wind_LR(      m,j,i))  + (climate%matrix%w_cold_P(   j,i) *     climate%matrix%GCM_cold%Wind_LR(      m,j,i))
        climate%Wind_DU( m,j,i) =  (climate%matrix%w_warm_P( j,i) *climate%matrix%GCM_warm%Wind_DU(      m,j,i))  + (climate%matrix%w_cold_P(   j,i) *     climate%matrix%GCM_cold%Wind_DU(      m,j,i))
      END DO
    END DO
    END DO
    ELSE
        climate%Wind_LR( :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%Wind_LR( :,:,grid%i1:grid%i2) 
        climate%Wind_DU( :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%Wind_DU( :,:,grid%i1:grid%i2) 
    END IF

    ! Downscale precipitation from the coarse-resolution reference
    ! GCM orography to the fine-resolution ice-model orography
    ! ========================================================

    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      ! Use the Roe&Lindzen precipitation model to do this; Berends et al., 2018, Eqs. A3-A7
      CALL adapt_precip_Roe( grid, Hs_GCM,   T_ref_GCM,  climate%Wind_LR, climate%Wind_DU, P_ref_GCM, &
                                   ice%Hs_a, T_ref_ice,  climate%Wind_LR, climate%Wind_DU, climate%Precip)

   ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use a simpler temperature-based correction; Berends et al., 2018, Eq. 14
      CALL adapt_precip_CC( grid, ice%Hs_a, Hs_GCM, T_ref_GCM, P_ref_GCM, climate%Precip, region_name)
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wT_ref_GCM)
    CALL deallocate_shared( wT_ref_ice)
    CALL deallocate_shared( wP_ref_GCM)
    CALL deallocate_shared( wHs_GCM)
    CALL deallocate_shared( wlambda_GCM)
    CALL deallocate_shared( wdHs_ice)
    CALL deallocate_shared( wdHs_snapshots)
    
    ! Safety checks
    CALL check_safety_precipitation( climate%Precip)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_matrix_precipitation
  SUBROUTINE initialise_climate_matrix( grid, climate, region_name, mask_noice)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_matrix'
    INTEGER                                            :: i,j,m
    LOGICAL                                            :: found_winds_PD_obs, found_winds_PI, found_winds_warm, found_winds_cold

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, climate%matrix%I_abs           , climate%matrix%wI_abs           )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%matrix%GCM_bias_T2m    , climate%matrix%wGCM_bias_T2m    )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%matrix%GCM_bias_Precip , climate%matrix%wGCM_bias_Precip )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, climate%matrix%GCM_bias_Hs     , climate%matrix%wGCM_bias_Hs     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%matrix%GCM_bias_Wind_LR, climate%matrix%wGCM_bias_Wind_LR)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%matrix%GCM_bias_Wind_DU, climate%matrix%wGCM_bias_Wind_DU)

    ! Allocated shared memory for climate matrix interpolation
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate%matrix%w_ins_T,  climate%matrix%ww_ins_T)  
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate%matrix%w_ice_T,  climate%matrix%ww_ice_T)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate%matrix%w_tot_T,  climate%matrix%ww_tot_T)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate%matrix%w_warm_P, climate%matrix%ww_warm_P)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate%matrix%w_cold_P, climate%matrix%ww_cold_P)
    CALL allocate_shared_dp_0D(                   climate%matrix%w_tot_P,  climate%matrix%ww_tot_P)   
    CALL allocate_shared_dp_0D(                   climate%matrix%w_EXT,    climate%matrix%ww_EXT)

    ! Allocate memory for the regional ERA40 climate and the final applied climate
    CALL allocate_climate_snapshot( grid, climate%matrix%PD_obs,   name = 'PD_obs'  )
    CALL allocate_climate_snapshot( grid, climate%matrix%GCM_PI,   name = 'GCM_PI'  )
    CALL allocate_climate_snapshot( grid, climate%matrix%GCM_warm, name = 'GCM_warm')
    CALL allocate_climate_snapshot( grid, climate%matrix%GCM_cold, name = 'GCM_cold')

    ! Read climate data from files
    CALL read_climate_snapshot( C%filename_PD_obs_climate       , grid, climate%matrix%PD_obs  , found_winds_PD_obs, region_name)
    CALL read_climate_snapshot( C%filename_climate_snapshot_PI  , grid, climate%matrix%GCM_PI  , found_winds_PI    , region_name)
    CALL read_climate_snapshot( C%filename_climate_snapshot_warm, grid, climate%matrix%GCM_warm, found_winds_warm  , region_name)
    CALL read_climate_snapshot( C%filename_climate_snapshot_cold, grid, climate%matrix%GCM_cold, found_winds_cold  , region_name)

    ! Get the orbit time
    climate%matrix%GCM_PI%orbit_time   = 0._dp 
    climate%matrix%GCM_warm%orbit_time = C%matrix_warm_orbit_time
    climate%matrix%GCM_cold%orbit_time = C%matrix_cold_orbit_time

    ! Safety
    IF (.NOT. found_winds_PD_obs) CALL crash('couldnt find wind fields for PD climate in file ' // TRIM( C%filename_PD_obs_climate))

    ! If no wind fields are provided in the GCM snapshots (as is usually the case), use the present-day observed winds instead
    IF (.NOT. found_winds_PI) THEN
      IF (C%do_climate_matrix_wind) CALL crash('no wind fields found for PI climate snapshot, but do_climate_matrix_wind is TRUE.')
      IF (par%master) CALL warning('no wind fields found for PI   climate snapshot; using PD observed winds instead')
      climate%matrix%GCM_PI%Wind_LR(   :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%wind_LR( :,:,grid%i1:grid%i2)
      climate%matrix%GCM_PI%Wind_DU(   :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%Wind_DU( :,:,grid%i1:grid%i2)
    END IF
    IF (.NOT. found_winds_warm) THEN
      IF (C%do_climate_matrix_wind) CALL crash('no wind fields found for warm climate snapshot, but do_climate_matrix_wind is TRUE.')
      IF (par%master) CALL warning('no wind fields found for warm climate snapshot; using PD observed winds instead')
      climate%matrix%GCM_warm%Wind_LR( :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%wind_LR( :,:,grid%i1:grid%i2)
      climate%matrix%GCM_warm%Wind_DU( :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%Wind_DU( :,:,grid%i1:grid%i2)
    END IF
    IF (.NOT. found_winds_cold) THEN
      IF (C%do_climate_matrix_wind) CALL crash('no wind fields found for cold climate snapshot, but do_climate_matrix_wind is TRUE.')
      IF (par%master) CALL warning('no wind fields found for cold climate snapshot; using PD observed winds instead')
      climate%matrix%GCM_cold%Wind_LR( :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%wind_LR( :,:,grid%i1:grid%i2)
      climate%matrix%GCM_cold%Wind_DU( :,:,grid%i1:grid%i2) = climate%matrix%PD_obs%Wind_DU( :,:,grid%i1:grid%i2)
    END IF
    CALL sync
    
    ! Determine the snapshot ocean and ice mask
    CALL initialise_snapshot_mask(climate%matrix%PD_obs,   grid, mask_noice, region_name, C%filename_snapshot_mask_PD_obs  )
    CALL initialise_snapshot_mask(climate%matrix%GCM_PI,   grid, mask_noice, region_name, C%filename_snapshot_mask_GCM_PI  )
    CALL initialise_snapshot_mask(climate%matrix%GCM_warm, grid, mask_noice, region_name, C%filename_snapshot_mask_GCM_warm  )
    CALL initialise_snapshot_mask(climate%matrix%GCM_cold, grid, mask_noice, region_name, C%filename_snapshot_mask_GCM_cold  )

    ! Calculate spatially variable lapse rate

    ! Use a uniform value for the warm snapshot [this assumes "warm" is actually identical to PI!]
    climate%matrix%GCM_warm%lambda( :,grid%i1:grid%i2) = C%constant_lapserate
    CALL sync

    IF     (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      CALL initialise_matrix_calc_spatially_variable_lapserate( grid, climate%matrix%GCM_PI, climate%matrix%GCM_cold)
    ELSEIF (region_name == 'GLR' .OR. region_name == 'ANT') THEN
      climate%matrix%GCM_cold%lambda( :,grid%i1:grid%i2) = C%constant_lapserate
      CALL sync
    END IF

    ! Calculate GCM bias
    CALL initialise_matrix_calc_GCM_bias( grid, climate%matrix%GCM_PI, climate%matrix%PD_obs, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip, climate%matrix%GCM_bias_Hs, &
      climate%matrix%GCM_bias_Wind_LR,  climate%matrix%GCM_bias_Wind_DU, region_name)

    ! Apply bias correction on warm snapshot
    IF (C%climate_matrix_biascorrect_warm) CALL initialise_matrix_apply_bias_correction(   grid, climate%matrix%GCM_warm, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip, climate%matrix%GCM_bias_Hs, &
      climate%matrix%GCM_bias_Wind_LR,  climate%matrix%GCM_bias_Wind_DU )
     
    ! Apply bias correction on cold snapshot
    IF (C%climate_matrix_biascorrect_cold) CALL initialise_matrix_apply_bias_correction(    grid, climate%matrix%GCM_cold, &
      climate%matrix%GCM_bias_T2m, climate%matrix%GCM_bias_Precip, climate%matrix%GCM_bias_Hs, &
      climate%matrix%GCM_bias_Wind_LR,  climate%matrix%GCM_bias_Wind_DU)  

    ! Get reference absorbed insolation for the GCM snapshots
    CALL initialise_matrix_calc_absorbed_insolation( grid, climate%matrix%GCM_warm, region_name, mask_noice)
    CALL initialise_matrix_calc_absorbed_insolation( grid, climate%matrix%GCM_cold, region_name, mask_noice)

    ! Initialise applied climate with present-day observations
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      climate%T2m(     m,j,i) = climate%matrix%PD_obs%T2m(     m,j,i)
      climate%Precip(  m,j,i) = climate%matrix%PD_obs%Precip(  m,j,i)
      climate%Wind_LR( m,j,i) = climate%matrix%PD_obs%Wind_LR( m,j,i)
      climate%Wind_DU( m,j,i) = climate%matrix%PD_obs%Wind_DU( m,j,i)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_matrix
  SUBROUTINE initialise_matrix_calc_GCM_bias( grid, GCM_PI, PD_obs, GCM_bias_T2m, GCM_bias_Precip, GCM_bias_Hs, GCM_bias_Wind_LR, GCM_bias_Wind_DU, region_name)
    ! Calculate the GCM bias in temperature and precipitation
    !
    ! Account for the fact that the GCM PI snapshot has a lower resolution, and therefore
    ! a different surface elevation than the PD observed climatology!

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_snapshot),         INTENT(IN)    :: GCM_PI, PD_obs
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: GCM_bias_T2m
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: GCM_bias_Precip
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: GCM_bias_Hs
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: GCM_bias_Wind_LR
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: GCM_bias_Wind_DU

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_matrix_calc_GCM_bias'
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: T2m_SL_GCM, T2m_SL_obs
    REAL(dp), DIMENSION(:,:,:),          POINTER       :: PD_obs_corr
    INTEGER                                            :: wPD_obs_corr
    
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, PD_obs_corr, wPD_obs_corr)
    
    ! Calculate temperature bias
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
    ! === Topography === 
    GCM_bias_Hs(      j,i) = GCM_PI%Hs( j,i) - PD_obs%Hs( j,i)

    DO m = 1, 12

      ! === Temperature ===
      ! Scale modelled and observed temperature to sea level using a constant lapse rate
      T2m_SL_GCM = GCM_PI%T2m( m,j,i) + GCM_PI%Hs( j,i) * C%constant_lapserate
      T2m_SL_obs = PD_obs%T2m( m,j,i) + PD_obs%Hs( j,i) * C%constant_lapserate

      ! Calculate bias
      GCM_bias_T2m(    m,j,i) = T2m_SL_GCM            - T2m_SL_obs

      ! === Precipitation ===
      ! Apply bias correction
      GCM_bias_Precip( m,j,i) = GCM_PI%Precip( m,j,i) / PD_obs%Precip( m,j,i)
      
      ! === Wind ===
      ! Calculate the bias from wind.
      GCM_bias_Wind_LR( m,j,i) = GCM_PI%Wind_LR( m,j,i) - PD_obs%Wind_LR( m,j,i)
      GCM_bias_Wind_DU( m,j,i) = GCM_PI%Wind_DU( m,j,i) - PD_obs%Wind_DU( m,j,i)

    END DO
    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself 
    CALL deallocate_shared(wPD_obs_corr)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_matrix_calc_GCM_bias

   SUBROUTINE initialise_matrix_apply_bias_correction( grid, snapshot, bias_T2m, bias_Precip, bias_Hs, bias_Wind_LR, bias_Wind_DU)
    ! Apply a bias correction to this (GCM) snapshot

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_climate_snapshot),          INTENT(INOUT) :: snapshot
    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)    :: bias_T2m
    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)    :: bias_Precip
    REAL(dp), DIMENSION(:,:  ),           INTENT(IN)    :: bias_Hs
    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)    :: bias_Wind_LR
    REAL(dp), DIMENSION(:,:,:),           INTENT(IN)    :: bias_Wind_DU

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_matrix_apply_bias_correction'
    INTEGER                                             :: i,j,m

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%climate_matrix_biascorrect_cold) CALL crash('climate_matrix_biascorrect_cold  = .TRUE.: The bias correction in its current state should only be applied to PI/PD')

    ! Apply bias correction
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    snapshot%Hs( j,i) = snapshot%Hs( j,i) - bias_Hs( j,i)
    DO m = 1, 12
      snapshot%T2m(     m,j,i) = snapshot%T2m(     m,j,i) - bias_T2m(     m,j,i)
      snapshot%Precip(  m,j,i) = snapshot%Precip(  m,j,i) / bias_Precip(  m,j,i)
      snapshot%Wind_LR( m,j,i) = snapshot%Wind_LR( m,j,i) - bias_Wind_LR( m,j,i)
      snapshot%Wind_DU( m,j,i) = snapshot%Wind_DU( m,j,i) - bias_Wind_DU( m,j,i)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_matrix_apply_bias_correction
  
  SUBROUTINE initialise_matrix_calc_spatially_variable_lapserate( grid, snapshot_PI, snapshot)
    ! Calculate the spatially variable lapse-rate (for non-PI GCM climates; see Berends et al., 2018)
    ! Only meaningful for climates where there is ice (LGM, M2_Medium, M2_Large),
    ! and only intended for North America and Eurasia

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_climate_snapshot),          INTENT(IN)    :: snapshot_PI
    TYPE(type_climate_snapshot),          INTENT(INOUT) :: snapshot

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_matrix_calc_spatially_variable_lapserate'
    INTEGER                                             :: i,j,m
    INTEGER,  DIMENSION(:,:  ), POINTER                 ::  mask_calc_lambda
    INTEGER                                             :: wmask_calc_lambda
    REAL(dp)                                            :: dT_mean_nonice
    INTEGER                                             :: n_nonice, n_ice
    REAL(dp)                                            :: lambda_mean_ice

    REAL(dp), PARAMETER                                 :: lambda_min = 0.002_dp
    REAL(dp), PARAMETER                                 :: lambda_max = 0.05_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx, mask_calc_lambda, wmask_calc_lambda)

    ! Determine where the variable lapse rate should be calculated
    ! (i.e. where has the surface elevation increased substantially)
    ! ==============================================================

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (snapshot%Hs( j,i) > snapshot_PI%Hs( j,i) + 100._dp) THEN
        mask_calc_lambda( j,i) = 1
      ELSE
        mask_calc_lambda( j,i) = 0
      END IF

    END DO
    END DO
    CALL sync

    ! Calculate the regional average temperature change outside of the ice sheet
    ! ==========================================================================

    dT_mean_nonice = 0._dp
    n_nonice       = 0
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      IF (mask_calc_lambda( j,i) == 0) THEN
        dT_mean_nonice = dT_mean_nonice + snapshot%T2m( m,j,i) - snapshot_PI%T2m( m,j,i)
        n_nonice = n_nonice + 1
      END IF
    END DO
    END DO
    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dT_mean_nonice, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_nonice,       1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)

    dT_mean_nonice = dT_mean_nonice / REAL(n_nonice,dp)

    ! Calculate the lapse rate over the ice itself
    ! ============================================

    lambda_mean_ice = 0._dp
    n_ice           = 0

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (mask_calc_lambda( j,i) == 1) THEN

        DO m = 1, 12
          ! Berends et al., 2018 - Eq. 10
          snapshot%lambda( j,i) = snapshot%lambda( j,i) + 1/12._dp * MAX(lambda_min, MIN(lambda_max, &
            -(snapshot%T2m( m,j,i) - (snapshot_PI%T2m( m,j,i) + dT_mean_nonice)) / (snapshot%Hs( j,i) - snapshot_PI%Hs( j,i))))
        END DO

        lambda_mean_ice = lambda_mean_ice + snapshot%lambda( j,i)
        n_ice = n_ice + 1

      END IF

    END DO
    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, lambda_mean_ice, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_ice,           1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)

    lambda_mean_ice = lambda_mean_ice / n_ice

    ! Apply mean lapse-rate over ice to the rest of the region
    ! ========================================================

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (mask_calc_lambda( j,i) == 0) snapshot%lambda( j,i) = lambda_mean_ice
    END DO
    END DO
    CALL sync

    ! Smooth the lapse rate field with a 160 km Gaussian filter
    CALL smooth_Gaussian_2D( grid, snapshot%lambda, 160000._dp)

    ! Normalise the entire region to a mean lapse rate of 8 K /km
    snapshot%lambda( :,grid%i1:grid%i2) = snapshot%lambda( :,grid%i1:grid%i2) * (C%constant_lapserate / lambda_mean_ice)

    ! Clean up after yourself
    CALl deallocate_shared( wmask_calc_lambda)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_matrix_calc_spatially_variable_lapserate
  SUBROUTINE initialise_matrix_calc_absorbed_insolation( grid, snapshot, region_name, mask_noice)
    ! Calculate the yearly absorbed insolation for this (regional) GCM snapshot, to be used in the matrix interpolation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_climate_snapshot),          INTENT(INOUT) :: snapshot
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name
    INTEGER,  DIMENSION(:,:  ),           INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_matrix_calc_absorbed_insolation'
    INTEGER                                             :: i,j,m
    TYPE(type_ice_model)                                :: ice_dummy
    TYPE(type_climate_model)                            :: climate_dummy
    TYPE(type_SMB_model)                                :: SMB_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get insolation at the desired time from the insolation NetCDF file
    ! ==================================================================

    CALL get_insolation_at_time( grid, snapshot%orbit_time, snapshot%Q_TOA)

    ! Create temporary "dummy" climate, ice & SMB data structures,
    ! so we can run the SMB model and determine the reference albedo field
    ! ====================================================================

    ! Climate
    ! =======

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_dummy%T2m,    climate_dummy%wT2m)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_dummy%Precip, climate_dummy%wPrecip)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_dummy%Q_TOA,  climate_dummy%wQ_TOA)

    ! Copy climate fields
    climate_dummy%T2m(    :,:,grid%i1:grid%i2) = snapshot%T2m(    :,:,grid%i1:grid%i2)
    climate_dummy%Precip( :,:,grid%i1:grid%i2) = snapshot%Precip( :,:,grid%i1:grid%i2)
    climate_dummy%Q_TOA(  :,:,grid%i1:grid%i2) = snapshot%Q_TOA(  :,:,grid%i1:grid%i2)

    ! Ice
    ! ===

    CALL allocate_shared_int_2D(    grid%ny, grid%nx, ice_dummy%mask_ocean_a   , ice_dummy%wmask_ocean_a   )
    CALL allocate_shared_int_2D(    grid%ny, grid%nx, ice_dummy%mask_ice_a     , ice_dummy%wmask_ice_a     )
    CALL allocate_shared_int_2D(    grid%ny, grid%nx, ice_dummy%mask_shelf_a   , ice_dummy%wmask_shelf_a   )

    ice_dummy%mask_ice_a(    :,grid%i1:grid%i2) = snapshot%mask_ice(      :,grid%i1:grid%i2)
    ice_dummy%mask_ocean_a(  :,grid%i1:grid%i2) = snapshot%mask_ocean(    :,grid%i1:grid%i2)
    ice_dummy%mask_shelf_a(  :,grid%i1:grid%i2) = snapshot%mask_shelf(    :,grid%i1:grid%i2)

    ! SMB
    ! ===

    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB_dummy%AlbedoSurf      , SMB_dummy%wAlbedoSurf      )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB_dummy%MeltPreviousYear, SMB_dummy%wMeltPreviousYear)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB_dummy%FirnDepth       , SMB_dummy%wFirnDepth       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB_dummy%Rainfall        , SMB_dummy%wRainfall        )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB_dummy%Snowfall        , SMB_dummy%wSnowfall        )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB_dummy%AddedFirn       , SMB_dummy%wAddedFirn       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB_dummy%Melt            , SMB_dummy%wMelt            )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB_dummy%Refreezing      , SMB_dummy%wRefreezing      )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB_dummy%Refreezing_year , SMB_dummy%wRefreezing_year )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB_dummy%Runoff          , SMB_dummy%wRunoff          )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB_dummy%Albedo          , SMB_dummy%wAlbedo          )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB_dummy%Albedo_year     , SMB_dummy%wAlbedo_year     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB_dummy%SMB             , SMB_dummy%wSMB             )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB_dummy%SMB_year        , SMB_dummy%wSMB_year        )

    CALL allocate_shared_dp_0D( SMB_dummy%C_abl_constant, SMB_dummy%wC_abl_constant)
    CALL allocate_shared_dp_0D( SMB_dummy%C_abl_Ts,       SMB_dummy%wC_abl_Ts      )
    CALL allocate_shared_dp_0D( SMB_dummy%C_abl_Q,        SMB_dummy%wC_abl_Q       )
    CALL allocate_shared_dp_0D( SMB_dummy%C_refr,         SMB_dummy%wC_refr        )

    IF (par%master) THEN
      IF     (region_name == 'NAM') THEN
        SMB_dummy%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_NAM
        SMB_dummy%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_NAM
        SMB_dummy%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_NAM
        SMB_dummy%C_refr         = C%SMB_IMAUITM_C_refr_NAM
      ELSEIF (region_name == 'EAS') THEN
        SMB_dummy%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_EAS
        SMB_dummy%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_EAS
        SMB_dummy%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_EAS
        SMB_dummy%C_refr         = C%SMB_IMAUITM_C_refr_EAS
      ELSEIF (region_name == 'GRL') THEN
        SMB_dummy%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_GRL
        SMB_dummy%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_GRL
        SMB_dummy%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_GRL
        SMB_dummy%C_refr         = C%SMB_IMAUITM_C_refr_GRL
      ELSEIF (region_name == 'ANT') THEN
        SMB_dummy%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_ANT
        SMB_dummy%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_ANT
        SMB_dummy%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_ANT
        SMB_dummy%C_refr         = C%SMB_IMAUITM_C_refr_ANT
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Run the SMB model for 30 years for this particular climate
    ! (experimentally determined to be long enough to converge)
    DO i = 1, 30
      CALL run_SMB_model( grid, ice_dummy, climate_dummy, 0._dp, SMB_dummy, mask_noice)
    END DO
    CALL sync
    
    ! Calculate yearly total absorbed insolation
    snapshot%I_abs( :,grid%i1:grid%i2) = 0._dp
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      snapshot%I_abs( j,i) = snapshot%I_abs( j,i) + snapshot%Q_TOA( m,j,i) * (1._dp - SMB_dummy%Albedo( m,j,i))
    END DO
    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( ice_dummy%wmask_ocean_a)
    CALL deallocate_shared( ice_dummy%wmask_ice_a)
    CALL deallocate_shared( ice_dummy%wmask_shelf_a)
    CALL deallocate_shared( climate_dummy%wT2m)
    CALL deallocate_shared( climate_dummy%wPrecip)
    CALL deallocate_shared( climate_dummy%wQ_TOA)
    CALL deallocate_shared( SMB_dummy%wAlbedoSurf)
    CALL deallocate_shared( SMB_dummy%wMeltPreviousYear)
    CALL deallocate_shared( SMB_dummy%wFirnDepth)
    CALL deallocate_shared( SMB_dummy%wRainfall)
    CALL deallocate_shared( SMB_dummy%wSnowfall)
    CALL deallocate_shared( SMB_dummy%wAddedFirn)
    CALL deallocate_shared( SMB_dummy%wMelt)
    CALL deallocate_shared( SMB_dummy%wRefreezing)
    CALL deallocate_shared( SMB_dummy%wRefreezing_year)
    CALL deallocate_shared( SMB_dummy%wRunoff)
    CALL deallocate_shared( SMB_dummy%wAlbedo)
    CALL deallocate_shared( SMB_dummy%wAlbedo_year)
    CALL deallocate_shared( SMB_dummy%wSMB)
    CALL deallocate_shared( SMB_dummy%wSMB_year)
    CALL deallocate_shared( SMB_dummy%wC_abl_constant)
    CALL deallocate_shared( SMB_dummy%wC_abl_Ts)
    CALL deallocate_shared( SMB_dummy%wC_abl_Q)
    CALL deallocate_shared( SMB_dummy%wC_refr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_matrix_calc_absorbed_insolation

! == Some generally useful tools
! ==============================

  ! Allocate memory for a single climate snapshot
  SUBROUTINE allocate_climate_snapshot( grid, snapshot, name)
    ! Allocate shared memory for a single climate snapshot

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_snapshot),         INTENT(INOUT) :: snapshot
    CHARACTER(LEN=*),                    INTENT(IN)    :: name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_climate_snapshot'

    ! Add routine to path
    CALL init_routine( routine_name)

    snapshot%name = name

    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, snapshot%Hs,             snapshot%wHs            )
    CALL allocate_shared_int_2D(    grid%ny, grid%nx, snapshot%mask_ice,       snapshot%wmask_ice      )
    CALL allocate_shared_int_2D(    grid%ny, grid%nx, snapshot%mask_ocean,     snapshot%wmask_ocean    )
    CALL allocate_shared_int_2D(    grid%ny, grid%nx, snapshot%mask_shelf,     snapshot%wmask_shelf    )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%T2m,            snapshot%wT2m           )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Precip,         snapshot%wPrecip        )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Wind_WE,        snapshot%wWind_WE       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Wind_SN,        snapshot%wWind_SN       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Wind_LR,        snapshot%wWind_LR       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Wind_DU,        snapshot%wWind_DU       )
    

    CALL allocate_shared_dp_0D(                       snapshot%CO2,            snapshot%wCO2           )
    CALL allocate_shared_dp_0D(                       snapshot%orbit_time,     snapshot%worbit_time    )
    CALL allocate_shared_dp_0D(                       snapshot%orbit_ecc,      snapshot%worbit_ecc     )
    CALL allocate_shared_dp_0D(                       snapshot%orbit_obl,      snapshot%worbit_obl     )
    CALL allocate_shared_dp_0D(                       snapshot%orbit_pre,      snapshot%worbit_pre     )
    CALL allocate_shared_dp_0D(                       snapshot%sealevel,       snapshot%wsealevel      )

    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, snapshot%lambda,         snapshot%wlambda        )

    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Q_TOA,          snapshot%wQ_TOA         )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, snapshot%Albedo,         snapshot%wAlbedo        )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, snapshot%I_abs,          snapshot%wI_abs         )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_climate_snapshot
  SUBROUTINE read_climate_snapshot( filename, grid, snapshot, found_winds, region_name)
    ! Read a climate snapshot from a NetCDF file. Works both for global lon/lat files and regional x/y files

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                 INTENT(IN)    :: filename
    TYPE(type_grid),                    INTENT(IN)    :: grid
    TYPE(type_climate_snapshot),        INTENT(INOUT) :: snapshot
    LOGICAL,                            INTENT(OUT)   :: found_winds
    CHARACTER(LEN=3),                   INTENT(IN)    :: region_name
    

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'read_climate_snapshot'
    INTEGER                                           :: found_wind_WE, found_wind_SN, found_wind_LR, found_wind_DU
    LOGICAL                                           :: found_winds_LRDU 
    
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if wind fields are included in this file; if not, return -1
    CALL inquire_var_multiple_options( filename, field_name_options_wind_WE, found_wind_WE)
    CALL inquire_var_multiple_options( filename, field_name_options_wind_SN, found_wind_SN)

    ! Winds may also be LR or DU, check for that as well
    CALL inquire_var_multiple_options( filename, 'Wind_LR', found_wind_LR)
    CALL inquire_var_multiple_options( filename, 'Wind_DU', found_wind_DU)

	! Check if South-North / East-West winds exist
    IF (found_wind_WE /= -1 .AND. found_wind_SN /= -1) THEN
         found_winds = .TRUE.
    ELSE
         found_winds = .FALSE.
    END IF

    ! Check instead if the file contains Left-Right / Up-Down winds
    IF (found_wind_DU /= -1 .AND. found_wind_LR /= -1) THEN
         found_winds_LRDU = .TRUE.
    ELSE
         found_winds_LRDU = .FALSE.
    END IF

    CALL sync

    ! Read the climate snapshot regardless of x/y or lon/lat grid
    CALL read_field_from_file_2D(         filename, field_name_options_Hs , grid, snapshot%Hs,     region_name)
    CALL read_field_from_file_2D_monthly( filename, 'T2m'                 , grid, snapshot%T2m,    region_name)
    CALL read_field_from_file_2D_monthly( filename, 'Precip'              , grid, snapshot%Precip, region_name)

    IF (found_winds) THEN
         CALL read_field_from_file_2D_monthly( filename, field_name_options_wind_WE, grid, snapshot%Wind_WE, region_name)
         CALL read_field_from_file_2D_monthly( filename, field_name_options_wind_SN, grid, snapshot%Wind_SN, region_name)

         ! Make sure to project SN and WE winds to DU and LR winds
         CALL rotate_wind_to_model_grid( grid, snapshot%wind_WE, snapshot%wind_SN, snapshot%wind_LR, snapshot%wind_DU)
    ELSEIF (found_winds_LRDU) THEN 
         ! LR and DU winds are already correctly rotated, so these field need to be loaded.
         CALL read_field_from_file_2D_monthly( filename, 'Wind_LR', grid, snapshot%wind_LR, region_name)
         CALL read_field_from_file_2D_monthly( filename, 'Wind_DU', grid, snapshot%wind_DU, region_name)
    ELSE
        ! There are no wind fields in the file
    END IF

    ! Half wind to see if that works
    snapshot%wind_LR( :,:,grid%i1:grid%i2) = snapshot%wind_LR( :,:,grid%i1:grid%i2) / 2._dp
    snapshot%wind_DU( :,:,grid%i1:grid%i2) = snapshot%wind_DU( :,:,grid%i1:grid%i2) / 2._dp
    
    ! Safety checks
    CALL check_safety_temperature(   snapshot%T2m   )
    CALL check_safety_precipitation( snapshot%Precip)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_climate_snapshot

  SUBROUTINE initialise_snapshot_mask(snapshot, grid, mask_noice, region_name, filename_snapshot_mask)
    ! Load or create the mask belonging to the climate snapshot
    ! These mask are currently used in two places in the climate matrix method:
    ! 1) To calculate a reference albedo /absorbed insolation. 
    ! 2) To determine the location of ice sheets, to interpolate the precipitation.
    !
    ! Two methods can be used (see reference_mask_method config):
    ! 	- estimate: Estimate land and ice masks based on topography and temperature
    !	  (Which may go wrong if the topography is bumpy in the ocean, which happens quite
    ! 	   often.)
    !   -  file: Determine masks based on an external file, which fully circumnavigates 
    !	   problems with topography. The file does not need to have the same resolution
    !      as the climate forcing.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_climate_snapshot),          INTENT(INOUT) :: snapshot
    INTEGER,  DIMENSION(:,:  ),           INTENT(IN)    :: mask_noice
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name
    CHARACTER(LEN=*),                     INTENT(IN)    :: filename_snapshot_mask
    

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_snapshot_mask'
    INTEGER                                             :: i,j
    REAL(dp), DIMENSION(:,: ), POINTER                  ::  mask_ice_raw,  mask_ocean_raw
    INTEGER                                             :: wmask_ice_raw, wmask_ocean_raw        
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Determine ice mask
    ! =============================
    
    ! Estimate the albedo from the climate snapshots based on the topography
    IF (C%reference_mask_method == 'estimate') THEN

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny

        ! If surface topography is the same as lowest surface topography, it is likely ocean
        IF (snapshot%Hs( j,i) == MINVAL(snapshot%Hs)) THEN
          snapshot%mask_ocean( j,i) = 1
        ELSE
          snapshot%mask_ocean( j,i) = 0
        END IF

        ! If LGM-PI topography difference is large and the temperature is low, it is likely an ice sheet
        IF ((ABS(snapshot%Hs( j,i)) > 100._dp) .AND. SUM(snapshot%T2m( :,j,i)) / 12._dp < 0._dp) THEN
          snapshot%mask_ice(   j,i) = 1
        ELSE
          snapshot%mask_ice(   j,i) = 0
        END IF

        ! mask_shelf is used in the SMB model only to find open ocean; since mask_ocean
        ! in this case already marks only open ocean, no need to look for shelves
        snapshot%mask_shelf( j,i) = 0

      END DO
      END DO
      CALL sync

    ! Use a file to obtain the albedo reference ice masks
    !   NOTE: Make sure the input file uses 1 for ice (ocean) and 0 for ice-free (land).
    ELSEIF (C%reference_mask_method == 'file') THEN

      ! Allocate raw mask fields
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, mask_ice_raw,   wmask_ice_raw   )
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, mask_ocean_raw, wmask_ocean_raw ) 
     
      ! Read the files
      CALL read_field_from_file_2D( filename_snapshot_mask, field_name_options_mask_ice,   grid, mask_ice_raw,       region_name)
      CALL read_field_from_file_2D( filename_snapshot_mask, field_name_options_mask_ocean, grid, mask_ocean_raw,     region_name)
     
      ! Often GCM mask will use 100 (all ice) and 0 (no ice) to show how much ice is in the domain. We need to convert that to 1 and 0

      ! - Ice -
      IF ((MAXVAL(mask_ice_raw) >= 2._dp .AND. MAXVAL(mask_ice_raw) < 110._dp)) THEN
        ! We can resonably assume the GCM people used percentages
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
             IF (mask_ice_raw( j,i) > 50._dp) THEN
               snapshot%mask_ice( j,i) = 1
             ELSE
               snapshot%mask_ice( j,i) = 0
             END IF
        END DO
        END DO
        CALL sync

     ELSEIF (MAXVAL(mask_ice_raw) < 2._dp) THEN
        ! We can resonably assume the GCM people used fractions
        ! Note that this may go wrong in case a domain contains no ice. But then
        ! it should matter little.
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
             IF (mask_ice_raw( j,i) > 0.5_dp) THEN
               snapshot%mask_ice( j,i) = 1
             ELSE
               snapshot%mask_ice( j,i) = 0
             END IF
        END DO
        END DO
        CALL sync
     ELSE
        CALL crash('The ice mask file may not contain percentages or fractions')
     END IF

      ! - Ocean - 
      IF ((MAXVAL(mask_ocean_raw) >= 2._dp .AND. MAXVAL(mask_ocean_raw) < 110._dp)) THEN
        ! We can resonably assume the GCM people used percentages
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
             IF (mask_ocean_raw( j,i) < 50._dp) THEN
               snapshot%mask_ocean( j,i) = 1
             ELSE
               snapshot%mask_ocean( j,i) = 0
             END IF
        END DO
        END DO
        CALL sync
     ELSEIF (MAXVAL(mask_ocean_raw) < 2._dp) THEN
        ! We can resonably assume the GCM people used fractions
        ! Note that this may go wrong in case a domain contains no ice. But then 
        ! it should matter little.
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
             IF (mask_ocean_raw( j,i) < 0.5_dp) THEN
               snapshot%mask_ocean( j,i) = 1
             ELSE
               snapshot%mask_ocean( j,i) = 0
             END IF
        END DO
        END DO
        CALL sync
      ELSE
        CALL crash('The ocean mask file may not contain percentages or fractions')
      END IF

      ! Snapshots and ice can be seen as the same thing for albedo and precipitation topography.
      snapshot%mask_shelf( :, grid%i1:grid%i2) = 0
      
      ! Clean up your mess
      CALL deallocate_shared(wmask_ice_raw)
      CALL deallocate_shared(wmask_ocean_raw)

    ELSE
      ! A method was given that does not exist; crash
      CALL crash('reference_mask_method "'//TRIM( C%reference_mask_method)//'" not found!')
    END IF
    
    ! Remove ice where there should be none
    ! =====================================
    
    ! mask_noice (no ice) cells should be ice and land free
    !		For example, Greenland in the Eurasian domain
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (mask_noice( j,i) == 1) THEN
        snapshot%mask_ice(     j,i) = 0
        snapshot%mask_ocean(   j,i) = 1
        snapshot%mask_shelf(   j,i) = 0
      END IF
    END DO
    END DO
    CALL sync

    ! Safety checks
    CALL check_for_NaN_int_2D( snapshot%mask_ice  )
    CALL check_for_NaN_int_2D( snapshot%mask_ocean)
    CALL check_for_NaN_int_2D( snapshot%mask_shelf)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_snapshot_mask
  
  ! Two different parameterised precipitation models:
  ! - a simply Clausius-Clapeyron-based method            (used for GRL and ANT)
  ! - the Roe & Lindzen temperature/orography-based model (used for NAM and EAS)
  SUBROUTINE adapt_precip_CC(  grid, Hs, Hs_GCM, T_ref_GCM, P_ref_GCM, Precip_GCM, region_name)

    USE parameters_module, ONLY: T0

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Hs              ! Model orography (m)
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Hs_GCM          ! Reference orography (m)           - total ice-weighted
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: T_ref_GCM       ! Reference temperature (K)         - total ice-weighted
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: P_ref_GCM       ! Reference precipitation (m/month) - total ice-weighted
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Output variables:
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: Precip_GCM      ! Climate matrix precipitation

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'adapt_precip_CC'
    INTEGER                                            :: i,j,m
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  T_inv,  T_inv_ref
    INTEGER                                            :: wT_inv, wT_inv_ref

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, T_inv,     wT_inv    )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, T_inv_ref, wT_inv_ref)

    ! Calculate inversion layer temperatures
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      T_inv_ref( m,j,i) = 88.9_dp + 0.67_dp *  T_ref_GCM( m,j,i)
      T_inv(     m,j,i) = 88.9_dp + 0.67_dp * (T_ref_GCM( m,j,i) - 0.008_dp * (Hs( j,i) - Hs_GCM( j,i)))
    END DO
    END DO
    END DO
    CALL sync

    IF     (region_name == 'GRL') THEN
      ! Method of Jouzel and Merlivat (1984), see equation (4.82) in Huybrechts (1992)

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
      DO m = 1, 12
        Precip_GCM( m,j,i) = P_ref_GCM( m,j,i) * 1.04**(T_inv( m,j,i) - T_inv_ref( m,j,i))
      END DO
      END DO
      END DO
      CALL sync

    ELSEIF (region_name == 'ANT') THEN
      ! As with Lorius/Jouzel method (also Huybrechts, 2002

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
      DO m = 1, 12
        Precip_GCM( m,j,i) = P_ref_GCM( m,j,i) * (T_inv_ref( m,j,i) / T_inv( m,j,i))**2 * EXP(22.47_dp * (T0 / T_inv_ref( m,j,i) - T0 / T_inv( m,j,i)))
      END DO
      END DO
      END DO
      CALL sync

    ELSE
      CALL crash('adapt_precip_CC should only be used for Greenland and Antarctica!')
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wT_inv)
    CALL deallocate_shared( wT_inv_ref)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE adapt_precip_CC
  SUBROUTINE adapt_precip_Roe( grid, Hs1, T2m1, Wind_LR1, Wind_DU1, Precip1, &
                                     Hs2, T2m2, Wind_LR2, Wind_DU2, Precip2)
    ! Adapt precipitation from reference state 1 to model state 2, using the Roe&Lindzen precipitation model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Hs1,      Hs2
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: T2m1,     T2m2
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: Wind_LR1, Wind_LR2
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: Wind_DU1, Wind_DU2
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: Precip1
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: Precip2

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'adapt_precip_Roe'
    INTEGER                                            :: i,j,m
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dHs_dx1,  dHs_dx2
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dHs_dy1,  dHs_dy2
    INTEGER                                            :: wdHs_dx1, wdHs_dx2
    INTEGER                                            :: wdHs_dy1, wdHs_dy2
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  Precip_RL1,  Precip_RL2,  dPrecip_RL
    INTEGER                                            :: wPrecip_RL1, wPrecip_RL2, wdPrecip_RL

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, dHs_dx1,     wdHs_dx1   )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, dHs_dx2,     wdHs_dx2   )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, dHs_dy1,     wdHs_dy1   )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, dHs_dy2,     wdHs_dy2   )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, Precip_RL1,  wPrecip_RL1)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, Precip_RL2,  wPrecip_RL2)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, dPrecip_RL,  wdPrecip_RL)

    ! Calculate surface slopes for both states
    CALL ddx_a_to_a_2D( grid, Hs1, dHs_dx1)
    CALL ddx_a_to_a_2D( grid, Hs2, dHs_dx2)
    CALL ddy_a_to_a_2D( grid, Hs1, dHs_dy1)
    CALL ddy_a_to_a_2D( grid, Hs2, dHs_dy2)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12

      ! Calculate precipitation with the Roe&Lindzen model for both states
      CALL precipitation_model_Roe( T2m1( m,j,i), dHs_dx1( j,i), dHs_dy1( j,i), Wind_LR1( m,j,i), Wind_DU1( m,j,i), Precip_RL1( m,j,i))
      CALL precipitation_model_Roe( T2m2( m,j,i), dHs_dx2( j,i), dHs_dy2( j,i), Wind_LR2( m,j,i), Wind_DU2( m,j,i), Precip_RL2( m,j,i))

      ! Calculate the ratio between those two precipitation rates
      dPrecip_RL( m,j,i) = MAX(0.01_dp, MIN( 2._dp, Precip_RL2( m,j,i) / Precip_RL1( m,j,i) ))

      ! Applied model precipitation = (matrix-interpolated GCM reference precipitation) * RL ratio
      Precip2( m,j,i) = Precip1( m,j,i) * dPrecip_RL( m,j,i)

    END DO
    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wdHs_dx1)
    CALL deallocate_shared( wdHs_dx2)
    CALL deallocate_shared( wdHs_dy1)
    CALL deallocate_shared( wdHs_dy2)
    CALL deallocate_shared( wPrecip_RL1)
    CALL deallocate_shared( wPrecip_RL2)
    CALL deallocate_shared( wdPrecip_RL)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE adapt_precip_Roe
  SUBROUTINE precipitation_model_Roe( T2m, dHs_dx, dHs_dy, Wind_LR, Wind_DU, Precip)
    ! Precipitation model of Roe (J. Glac, 2002), integration from Roe and Lindzen (J. Clim. 2001)

    USE parameters_module, ONLY: T0, pi, sec_per_year

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: T2m                  ! 2-m air temperature [K]
    REAL(dp),                            INTENT(IN)    :: dHs_dx               ! Surface slope in the x-direction [m/m]
    REAL(dp),                            INTENT(IN)    :: dHs_dy               ! Surface slope in the y-direction [m/m]
    REAL(dp),                            INTENT(IN)    :: Wind_LR              ! Wind speed    in the x-direction [m/s]
    REAL(dp),                            INTENT(IN)    :: Wind_DU              ! Wind speed    in the y-direction [m/s]
    REAL(dp),                            INTENT(OUT)   :: Precip               ! Modelled precipitation

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'precipitation_model_Roe'
    REAL(dp)                                           :: upwind_slope         ! Upwind slope
    REAL(dp)                                           :: E_sat                ! Saturation vapour pressure as function of temperature [Pa]
    REAL(dp)                                           :: x0                   ! Integration parameter x0 [m s-1]
    REAL(dp)                                           :: err_in,err_out

    REAL(dp), PARAMETER                                :: e_sat0  = 611.2_dp   ! Saturation vapour pressure at 273.15 K [Pa]
    REAL(dp), PARAMETER                                :: c_one   = 17.67_dp   ! Constant c1 []
    REAL(dp), PARAMETER                                :: c_two   = 243.5_dp   ! Constant c2 [Celcius]

    REAL(dp), PARAMETER                                :: a_par   = 2.5E-11_dp ! Constant a [m2 s  kg-1] (from Roe et al., J. Clim. 2001)
    REAL(dp), PARAMETER                                :: b_par   = 5.9E-09_dp ! Constant b [m  s2 kg-1] (from Roe et al., J. Clim. 2001)
    REAL(dp), PARAMETER                                :: alpha   = 100.0_dp   ! Constant alpha [s m-1]

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the upwind slope
    upwind_slope = MAX(0._dp, Wind_LR * dHs_dx + Wind_DU * dHs_dy)

    ! Calculate the saturation vapour pressure E_sat:
    E_sat = e_sat0 * EXP( c_one * (T2m - T0) / (c_two + T2m - T0) )

    ! Calculate integration parameter x0 = a/b + w (with w = wind times slope)
    x0 = a_par / b_par + upwind_slope

    ! Calculate the error function (2nd term on the r.h.s.)
    err_in = alpha * ABS(x0)
    CALL error_function(err_in,err_out)

    ! Calculate precipitation rate as in Appendix of Roe et al. (J. Clim, 2001)
    Precip = ( b_par * E_sat ) * ( x0 / 2._dp + x0**2 * err_out / (2._dp * ABS(x0)) + &
                                         EXP (-alpha**2 * x0**2) / (2._dp * SQRT(pi) * alpha) ) * sec_per_year

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE precipitation_model_Roe

  ! Rotate wind_WE, wind_SN to wind_LR, wind_DU
  SUBROUTINE rotate_wind_to_model_grid( grid, wind_WE, wind_SN, wind_LR, wind_DU)
    ! Code copied from ANICE.

    USE parameters_module, ONLY: pi

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: wind_WE
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: wind_SN
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: wind_LR
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: wind_DU

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'rotate_wind_to_model_grid'
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: longitude_start, Uwind_x, Uwind_y, Vwind_x, Vwind_y

    ! Add routine to path
    CALL init_routine( routine_name)

    ! First find the first longitude which defines the start of quadrant I:
    longitude_start = grid%lambda_M - 90._dp

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12

      ! calculate x and y from the zonal wind
      Uwind_x =   wind_WE( m,j,i) * SIN((pi/180._dp) * (grid%lon( j,i) - longitude_start))
      Uwind_y = - wind_WE( m,j,i) * COS((pi/180._dp) * (grid%lon( j,i) - longitude_start))

      ! calculate x and y from the meridional winds
      Vwind_x =   wind_SN( m,j,i) * COS((pi/180._dp) * (grid%lon( j,i) - longitude_start))
      Vwind_y =   wind_SN( m,j,i) * SIN((pi/180._dp) * (grid%lon( j,i) - longitude_start))

      ! Sum up wind components resulting in ice-grid windfield
      wind_LR( m,j,i) = Uwind_x + Vwind_x   ! winds left to right
      wind_DU( m,j,i) = Uwind_y + Vwind_y   ! winds bottom to top

    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE rotate_wind_to_model_grid

  ! Safety checks
  SUBROUTINE check_safety_temperature( T2m)
    ! Safety checks on a monthly temperature field

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: T2m

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_safety_temperature'
    INTEGER                                            :: n1,n2,n3,i1,i2,i,j,k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get field size
    n1 = SIZE( T2m,1)
    n2 = SIZE( T2m,2)
    n3 = SIZE( T2m,3)

    ! Parallelisation
    CALL partition_list( n3, par%i, par%n, i1, i2)

    ! Perform safety check
    DO i = i1, i2
    DO j = 1, n2
    DO k = 1, n1

      ! Temperature errors
      IF     (T2m( k,j,i) < 0._dp) THEN
        CALL crash('negative temperature detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (ISNAN( T2m( k,j,i))) THEN
        CALL crash('NaN temperature detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      END IF

      ! Temperature warnings
      IF     (T2m( k,j,i) < 150._dp) THEN
        CALL warning('excessively low temperature (< 150 K) detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (T2m( k,j,i) > 350._dp) THEN
        CALL warning('excessively high temperature (> 350 K) detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      END IF

    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_safety_temperature
  SUBROUTINE check_safety_precipitation( Precip)
    ! Safety checks on a monthly precipitation field

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: Precip

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_safety_precipitation'
    INTEGER                                            :: n1,n2,n3,i1,i2,i,j,k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get field size
    n1 = SIZE( Precip,1)
    n2 = SIZE( Precip,2)
    n3 = SIZE( Precip,3)

    ! Parallelisation
    CALL partition_list( n3, par%i, par%n, i1, i2)

    ! Perform safety check
    DO i = i1, i2
    DO j = 1, n2
    DO k = 1, n1

      ! Precipitation errors
      IF     (Precip( k,j,i) <= 0._dp) THEN
        CALL crash('zero/negative precipitation detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (ISNAN(Precip( k,j,i))) THEN
        CALL crash('NaN precipitation detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      END IF

      ! Precipitation warnings
      IF     (Precip( k,j,i) > 10._dp) THEN
        CALL warning('excessively high precipitation (> 10 m/month) detected at [{int_01},{int_02},{int_03}]', int_01 = k, int_02 = j, int_03 = i)
      END IF

    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_safety_precipitation

!! == Directly prescribed regional SMB
!! ===================================
!
!  SUBROUTINE run_climate_model_direct_SMB_regional( grid, climate_matrix, time)
!    ! Run the regional climate model
!    !
!    ! Use a directly prescribed regional SMB
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    TYPE(type_grid),                          INTENT(IN)    :: grid
!    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix
!    REAL(dp),                                 INTENT(IN)    :: time
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'run_climate_model_direct_SMB_regional'
!    REAL(dp)                                                :: wt0, wt1
!    INTEGER                                                 :: i,j
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Safety
!    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
!      CALL crash('choice_SMB_model should be "direct_regional"!')
!    END IF
!
!    ! Check if the requested time is enveloped by the two timeframes;
!    ! if not, read the two relevant timeframes from the NetCDF file
!    IF (time < climate_matrix%SMB_direct%t0 .OR. time > climate_matrix%SMB_direct%t1) THEN
!
!      ! Find and read the two global time frames
!      CALL sync
!      CALL update_direct_regional_SMB_timeframes_from_file( grid, climate_matrix%SMB_direct, time)
!
!    END IF ! IF (time >= climate_matrix%SMB_direct%t0 .AND. time <= climate_matrix%SMB_direct%t1) THEN
!
!    ! Interpolate the two timeframes in time
!    wt0 = (climate_matrix%SMB_direct%t1 - time) / (climate_matrix%SMB_direct%t1 - climate_matrix%SMB_direct%t0)
!    wt1 = 1._dp - wt0
!
!    DO i = grid%i1, grid%i2
!    DO j = 1, grid%ny
!      climate_matrix%applied%T2m( :,j,i) = (wt0 * climate_matrix%SMB_direct%T2m_year0( j,i)) + (wt1 * climate_matrix%SMB_direct%T2m_year1( j,i))
!    END DO
!    END DO
!    CALL sync
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE run_climate_model_direct_SMB_regional
!  SUBROUTINE update_direct_regional_SMB_timeframes_from_file( grid, clim_reg, time)
!    ! Read the NetCDF file containing the regional climate forcing data. Only read the time
!    ! frames enveloping the current coupling timestep to save on memory usage.
!
!    IMPLICIT NONE
!
!    TYPE(type_grid),                            INTENT(IN)    :: grid
!    TYPE(type_direct_SMB_forcing_regional),     INTENT(INOUT) :: clim_reg
!    REAL(dp),                                   INTENT(IN)    :: time
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_direct_regional_SMB_timeframes_from_file'
!    INTEGER                                            :: ti0, ti1
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Safety
!    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
!      CALL crash('choice_SMB_model should be "direct_regional"!')
!    END IF
!
!    ! Find time indices to be read
!    IF (par%master) THEN
!
!      IF     (time < clim_reg%time( 1)) THEN
!
!        CALL warning('using constant start-of-record SMB when extrapolating!')
!        ti0 = 1
!        ti1 = 1
!        clim_reg%t0 = clim_reg%time( ti0) - 1._dp
!        clim_reg%t1 = clim_reg%time( ti1)
!
!      ELSEIF (time <= clim_reg%time( clim_reg%nyears)) THEN
!
!        ti1 = 1
!        DO WHILE (clim_reg%time( ti1) < time)
!          ti1 = ti1 + 1
!        END DO
!        ti0 = ti1 - 1
!
!        IF (ti0 == 0) THEN
!          ti0 = 1
!          ti1 = 2
!        ELSEIF (ti1 == clim_reg%nyears) THEN
!          ti0 = clim_reg%nyears - 1
!          ti1 = ti0 + 1
!        END IF
!
!        clim_reg%t0 = clim_reg%time( ti0)
!        clim_reg%t1 = clim_reg%time( ti1)
!
!      ELSE ! IF     (time < clim_reg%time( 1)) THEN
!
!        CALL warning('using constant end-of-record SMB when extrapolating!')
!        ti0 = clim_reg%nyears
!        ti1 = clim_reg%nyears
!        clim_reg%t0 = clim_reg%time( ti0) - 1._dp
!        clim_reg%t1 = clim_reg%time( ti1)
!
!      END IF ! IF     (time < clim_reg%time( 1)) THEN
!
!    END IF ! IF (par%master) THEN
!
!    ! Read new regional climate fields from the NetCDF file
!    IF (par%master) CALL read_direct_regional_SMB_file_timeframes( clim_reg, ti0, ti1)
!    CALL sync
!
!    ! Map the newly read data to the model grid
!    CALL map_square_to_square_cons_2nd_order_2D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%T2m_year0_raw, clim_reg%T2m_year0)
!    CALL map_square_to_square_cons_2nd_order_2D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%T2m_year1_raw, clim_reg%T2m_year1)
!    CALL map_square_to_square_cons_2nd_order_2D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%SMB_year0_raw, clim_reg%SMB_year0)
!    CALL map_square_to_square_cons_2nd_order_2D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%SMB_year1_raw, clim_reg%SMB_year1)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE update_direct_regional_SMB_timeframes_from_file
!  SUBROUTINE initialise_climate_model_direct_SMB_regional( grid, climate_matrix, region_name)
!    ! Initialise the regional climate model
!    !
!    ! Use a directly prescribed regional SMB
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    TYPE(type_grid),                          INTENT(IN)    :: grid
!    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix
!    CHARACTER(LEN=3),                         INTENT(IN)    :: region_name
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'initialise_climate_model_direct_SMB_regional'
!    INTEGER                                                 :: nx, ny
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Safety
!    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
!      CALL crash('choice_SMB_model should be "direct_regional"!')
!    END IF
!
!    ! The times at which we have climate fields from input, between which we'll interpolate
!    ! to find the climate at model time (t0 <= model_time <= t1)
!
!    CALL allocate_shared_dp_0D( climate_matrix%SMB_direct%t0, climate_matrix%SMB_direct%wt0)
!    CALL allocate_shared_dp_0D( climate_matrix%SMB_direct%t1, climate_matrix%SMB_direct%wt1)
!
!    IF (par%master) THEN
!      ! Give impossible values to timeframes, so that the first call to run_climate_model_direct_climate_global
!      ! is guaranteed to first read two new timeframes from the NetCDF file
!      climate_matrix%SMB_direct%t0 = C%start_time_of_run - 100._dp
!      climate_matrix%SMB_direct%t1 = C%start_time_of_run - 90._dp
!    END IF ! IF (par%master) THEN
!    CALL sync
!
!    ! Inquire into the direct global cliamte forcing netcdf file
!    CALL allocate_shared_int_0D( climate_matrix%SMB_direct%nyears, climate_matrix%SMB_direct%wnyears)
!    CALL allocate_shared_int_0D( climate_matrix%SMB_direct%nx_raw, climate_matrix%SMB_direct%wnx_raw)
!    CALL allocate_shared_int_0D( climate_matrix%SMB_direct%ny_raw, climate_matrix%SMB_direct%wny_raw)
!
!    ! Determine name of file to read data from
!    IF     (region_name == 'NAM') THEN
!      climate_matrix%SMB_direct%netcdf%filename = C%filename_direct_regional_SMB_NAM
!    ELSEIF (region_name == 'EAS') THEN
!      climate_matrix%SMB_direct%netcdf%filename = C%filename_direct_regional_SMB_EAS
!    ELSEIF (region_name == 'GRL') THEN
!      climate_matrix%SMB_direct%netcdf%filename = C%filename_direct_regional_SMB_GRL
!    ELSEIF (region_name == 'ANT') THEN
!      climate_matrix%SMB_direct%netcdf%filename = C%filename_direct_regional_SMB_ANT
!    END IF
!
!    IF (par%master) WRITE(0,*) '  Initialising direct regional SMB forcing from ', TRIM( climate_matrix%SMB_direct%netcdf%filename), '...'
!
!    IF (par%master) CALL inquire_direct_regional_SMB_forcing_file( climate_matrix%SMB_direct)
!    CALL sync
!
!    ! Abbreviations for shorter code
!    nx = climate_matrix%SMB_direct%nx_raw
!    ny = climate_matrix%SMB_direct%ny_raw
!
!    ! Allocate shared memory
!    CALL allocate_shared_dp_1D( climate_matrix%SMB_direct%nyears, climate_matrix%SMB_direct%time, climate_matrix%SMB_direct%wtime  )
!    CALL allocate_shared_dp_1D(               nx, climate_matrix%SMB_direct%x_raw,         climate_matrix%SMB_direct%wx_raw        )
!    CALL allocate_shared_dp_1D(      ny,          climate_matrix%SMB_direct%y_raw,         climate_matrix%SMB_direct%wy_raw        )
!
!    CALL allocate_shared_dp_2D(      ny,      nx, climate_matrix%SMB_direct%T2m_year0_raw, climate_matrix%SMB_direct%wT2m_year0_raw)
!    CALL allocate_shared_dp_2D(      ny,      nx, climate_matrix%SMB_direct%T2m_year1_raw, climate_matrix%SMB_direct%wT2m_year1_raw)
!    CALL allocate_shared_dp_2D(      ny,      nx, climate_matrix%SMB_direct%SMB_year0_raw, climate_matrix%SMB_direct%wSMB_year0_raw)
!    CALL allocate_shared_dp_2D(      ny,      nx, climate_matrix%SMB_direct%SMB_year1_raw, climate_matrix%SMB_direct%wSMB_year1_raw)
!
!    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%SMB_direct%T2m_year0,     climate_matrix%SMB_direct%wT2m_year0    )
!    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%SMB_direct%T2m_year1,     climate_matrix%SMB_direct%wT2m_year1    )
!    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%SMB_direct%SMB_year0,     climate_matrix%SMB_direct%wSMB_year0    )
!    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%SMB_direct%SMB_year1,     climate_matrix%SMB_direct%wSMB_year1    )
!
!    ! Read time and grid data
!    IF (par%master) CALL read_direct_regional_SMB_file_time_xy( climate_matrix%SMB_direct)
!    CALL sync
!
!    ! Lastly, allocate memory for the "applied" snapshot
!    CALL allocate_climate_snapshot( grid, climate_matrix%applied, name = 'applied')
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE initialise_climate_model_direct_SMB_regional

END MODULE climate_module
