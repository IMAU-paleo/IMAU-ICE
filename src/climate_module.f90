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
  USE data_types_module,               ONLY: type_model_region, type_grid, type_ice_model, type_SMB_model, type_climate_snapshot_global, &
                                             type_climate_matrix_global, type_climate_snapshot_regional, type_climate_matrix_regional, &
                                             type_direct_climate_forcing_global,   type_direct_SMB_forcing_global, type_reference_geometry, &
                                             type_direct_climate_forcing_regional, type_direct_SMB_forcing_regional
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, inquire_PD_obs_global_climate_file, read_PD_obs_global_climate_file, &
                                             inquire_GCM_global_climate_file, read_GCM_global_climate_file, inquire_direct_global_climate_forcing_file, &
                                             read_direct_global_climate_file_time_latlon, read_direct_global_climate_file_timeframes, &
                                             inquire_direct_regional_climate_forcing_file, read_direct_regional_climate_file_time_xy, &
                                             read_direct_regional_climate_file_timeframes, inquire_direct_global_SMB_forcing_file, &
                                             read_direct_global_SMB_file_time_latlon, read_direct_global_SMB_file_timeframes, &
                                             inquire_direct_regional_SMB_forcing_file, read_direct_regional_SMB_file_time_xy, &
                                             read_direct_regional_SMB_file_timeframes, inquire_ISMIP_forcing_SMB_baseline_file, &
                                             read_ISMIP_forcing_SMB_baseline_file, inquire_ISMIP_forcing_ST_baseline_file, &
                                             read_ISMIP_forcing_ST_baseline_file, inquire_ISMIP_forcing_aSMB_file, &
                                             read_ISMIP_forcing_aSMB_file, inquire_ISMIP_forcing_dSMBdz_file, &
                                             read_ISMIP_forcing_dSMBdz_file, inquire_ISMIP_forcing_aST_file, &
                                             read_ISMIP_forcing_aST_file, inquire_ISMIP_forcing_dSTdz_file, &
                                             read_ISMIP_forcing_dSTdz_file
  USE forcing_module,                  ONLY: forcing, get_insolation_at_time, update_CO2_at_model_time
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             error_function, smooth_Gaussian_2D, smooth_Shepard_2D, &
                                             map_glob_to_grid_2D, map_glob_to_grid_3D, transpose_dp_2D, transpose_dp_3D, &
                                             map_square_to_square_cons_2nd_order_2D, map_square_to_square_cons_2nd_order_3D
  USE derivatives_and_grids_module,    ONLY: ddx_a_to_a_2D, ddy_a_to_a_2D
  USE SMB_module,                      ONLY: run_SMB_model

  IMPLICIT NONE

CONTAINS

! == The main routines that should be called from the main ice model/program
! ==========================================================================

  SUBROUTINE run_climate_model( region, climate_matrix_global, time)
    ! Run the regional climate model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_climate_matrix_global),    INTENT(INOUT) :: climate_matrix_global
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_climate_model == 'none') THEN
      ! Possibly a direct SMB is used? Otherwise no need to do anything.

      IF     (C%choice_SMB_model == 'direct_global') THEN
        ! Use a directly prescribed global SMB
        CALL run_climate_model_direct_SMB_global( region%grid, climate_matrix_global%SMB_direct, region%climate_matrix, time)
      ELSEIF (C%choice_SMB_model == 'direct_regional') THEN
        ! Use a directly prescribed regional SMB
        CALL run_climate_model_direct_SMB_regional( region%grid, region%climate_matrix, time)
      END IF

    ELSEIF (C%choice_climate_model == 'idealised') THEN
      ! Assign some idealised temperature/precipitation

      CALL run_climate_model_idealised( region%grid, region%ice, region%climate_matrix%applied, time)

    ELSEIF (C%choice_climate_model == 'PD_obs') THEN
      ! Keep the climate fixed to present-day observed conditions

      CALL run_climate_model_PD_obs( region%grid, region%climate_matrix)

    ELSEIF (C%choice_climate_model == 'PD_dTglob') THEN
      ! Use the present-day climate plus a global temperature offset (de Boer et al., 2013)

      CALL run_climate_model_dT_glob( region%grid, region%ice, region%climate_matrix, region%name)

    ELSEIF (C%choice_climate_model == 'matrix_warm_cold') THEN
      ! Use the warm/cold climate matrix (Berends et al., 2018)

      CALL run_climate_model_matrix_warm_cold( region%grid, region%ice, region%SMB, region%climate_matrix, region%name, time)

    ELSEIF (C%choice_climate_model == 'direct_global') THEN
      ! Use a directly prescribed global climate

      CALL run_climate_model_direct_climate_global( region%grid, climate_matrix_global%direct, region%climate_matrix, time)

    ELSEIF (C%choice_climate_model == 'direct_regional') THEN
      ! Use a directly prescribed regional climate

      CALL run_climate_model_direct_climate_regional( region%grid, region%climate_matrix, time)

    ELSEIF (C%choice_climate_model == 'ISMIP_style_forcing') THEN
      ! Use the ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing

      CALL run_climate_model_ISMIP_forcing( region%grid, region%climate_matrix, time, region%refgeo_PD, region%ice)

    ELSE
      CALL crash('unknown choice_climate_model"' // TRIM(C%choice_climate_model) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model
  SUBROUTINE initialise_climate_model_global( climate_matrix)
    ! Initialise the global climate model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_climate_matrix_global),    INTENT(INOUT) :: climate_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_model_global'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) ' Initialising global climate model "', TRIM(C%choice_climate_model), '"...'

    IF     (C%choice_climate_model == 'none') THEN
      ! Possibly a direct SMB is used? Otherwise no need to do anything.

      IF     (C%choice_SMB_model == 'direct_global') THEN
        ! Use a directly prescribed global SMB
        CALL initialise_climate_model_direct_SMB_global( climate_matrix%SMB_direct)
      ELSEIF (C%choice_SMB_model == 'direct_regional') THEN
        ! Use a directly prescribed regional SMB, no need to initialise any global stuff
      END IF

    ELSEIF (C%choice_climate_model == 'idealised') THEN
      ! No need to initialise any global climate stuff

    ELSEIF (C%choice_climate_model == 'PD_obs') THEN
      ! Keep the climate fixed to present-day observed conditions

      CALL initialise_climate_model_global_PD_obs( climate_matrix)

    ELSEIF (C%choice_climate_model == 'PD_dTglob') THEN
      ! Use the present-day climate plus a global temperature offset (de Boer et al., 2013)

      CALL initialise_climate_model_global_PD_obs( climate_matrix)

    ELSEIF (C%choice_climate_model == 'matrix_warm_cold') THEN
      ! Use the warm/cold climate matrix (Berends et al., 2018)

      CALL initialise_climate_matrix_global( climate_matrix)

    ELSEIF (C%choice_climate_model == 'direct_global') THEN
      ! Use a directly prescribed global climate

      CALL initialise_climate_model_direct_climate_global( climate_matrix%direct)

    ELSEIF (C%choice_climate_model == 'direct_regional') THEN
      ! No need to initialise any global climate stuff

    ELSEIF (C%choice_climate_model == 'ISMIP_style_forcing') THEN
      ! Use the ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing

    ELSE
      CALL crash('unknown choice_climate_model"' // TRIM(C%choice_climate_model) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_global
  SUBROUTINE initialise_climate_model_regional( region, climate_matrix_global)
    ! Initialise the regional climate model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_climate_matrix_global),    INTENT(IN)    :: climate_matrix_global

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_model_regional'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) '  Initialising regional climate model "', TRIM(C%choice_climate_model), '"...'

    IF     (C%choice_climate_model == 'none') THEN
      ! Possibly a direct SMB is used? Otherwise no need to do anything.

      IF     (C%choice_SMB_model == 'direct_global') THEN
        ! Use a directly prescribed global SMB
        CALL initialise_climate_model_direct_SMB_global_regional( region%grid, region%climate_matrix)
      ELSEIF (C%choice_SMB_model == 'direct_regional') THEN
        ! Use a directly prescribed regional SMB, no need to initialise any global stuff
        CALL initialise_climate_model_direct_SMB_regional( region%grid, region%climate_matrix, region%name)
      END IF

    ELSEIF (C%choice_climate_model == 'idealised') THEN
      ! Only need to allocate memory for the "applied" regional snapshot

      CALL allocate_climate_snapshot_regional( region%grid, region%climate_matrix%applied, name = 'applied')

    ELSEIF (C%choice_climate_model == 'PD_obs') THEN
      ! Keep the climate fixed to present-day observed conditions

      CALL initialise_climate_model_regional_PD_obs( region%grid, climate_matrix_global, region%climate_matrix)

    ELSEIF (C%choice_climate_model == 'PD_dTglob') THEN
      ! Use the present-day climate plus a global temperature offset (de Boer et al., 2013)

      CALL initialise_climate_model_regional_PD_obs( region%grid, climate_matrix_global, region%climate_matrix)

    ELSEIF (C%choice_climate_model == 'matrix_warm_cold') THEN
      ! Use the warm/cold climate matrix (Berends et al., 2018)

      CALL initialise_climate_matrix_regional( region, climate_matrix_global)

    ELSEIF (C%choice_climate_model == 'direct_global') THEN
      ! Use a directly prescribed global climate

      CALL initialise_climate_model_direct_climate_global_regional( region%grid, region%climate_matrix)

    ELSEIF (C%choice_climate_model == 'direct_regional') THEN
      ! Use a directly prescribed global climate

      CALL initialise_climate_model_direct_climate_regional( region%grid, region%climate_matrix, region%name)

    ELSEIF (C%choice_climate_model == 'ISMIP_style_forcing') THEN
      ! Use the ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing

      CALL initialise_climate_model_ISMIP_forcing( region%grid, region%climate_matrix)

    ELSE
      CALL crash('unknown choice_climate_model"' // TRIM(C%choice_climate_model) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_regional

! == Idealised climates
! =====================

  SUBROUTINE run_climate_model_idealised( grid, ice, climate, time)
    ! Run the regional climate model
    !
    ! Assign some idealised temperature/precipitation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_climate_snapshot_regional), INTENT(INOUT) :: climate
    REAL(dp),                             INTENT(IN)    :: time

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
    TYPE(type_climate_snapshot_regional), INTENT(INOUT) :: climate
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

  SUBROUTINE run_climate_model_PD_obs( grid, climate_matrix)
    ! Run the regional climate model
    !
    ! Keep the climate fixed to present-day observed conditions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),  INTENT(INOUT) :: climate_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_PD_obs'
    INTEGER                                            :: i,j,m

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      climate_matrix%applied%T2m(     m,j,i) = climate_matrix%PD_obs%T2m(     m,j,i)
      climate_matrix%applied%Precip(  m,j,i) = climate_matrix%PD_obs%Precip(  m,j,i)
      climate_matrix%applied%Hs(        j,i) = climate_matrix%PD_obs%Hs(        j,i)
      climate_matrix%applied%Wind_LR( m,j,i) = climate_matrix%PD_obs%Wind_LR( m,j,i)
      climate_matrix%applied%Wind_DU( m,j,i) = climate_matrix%PD_obs%Wind_DU( m,j,i)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_PD_obs
  SUBROUTINE initialise_climate_model_global_PD_obs( climate_matrix)
    ! Initialise the global climate model
    !
    ! Keep the climate fixed to present-day observed conditions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_climate_matrix_global),    INTENT(INOUT) :: climate_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_model_global_PD_obs'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise the present-day observed global climate (e.g. ERA-40)
    CALL initialise_climate_PD_obs_global( climate_matrix%PD_obs, name = 'PD_obs')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_global_PD_obs
  SUBROUTINE initialise_climate_model_regional_PD_obs( grid, climate_matrix_global, climate_matrix)
    ! Initialise the regional climate model
    !
    ! Keep the climate fixed to present-day observed conditions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_matrix_global),    INTENT(IN)    :: climate_matrix_global
    TYPE(type_climate_matrix_regional),  INTENT(INOUT) :: climate_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_model_regional_PD_obs'
    INTEGER                                            :: i,j,m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise data structures for the regional ERA40 climate and the final applied climate
    CALL allocate_climate_snapshot_regional( grid, climate_matrix%PD_obs,  'PD_obs' )
    CALL allocate_climate_snapshot_regional( grid, climate_matrix%applied, 'applied')

    ! Map these subclimates from global grid to model grid
    CALL map_snapshot_to_grid( grid,  climate_matrix_global%PD_obs, climate_matrix%PD_obs)

    ! Initialise applied climate with present-day observations
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      climate_matrix%applied%T2m(     m,j,i) = climate_matrix%PD_obs%T2m(     m,j,i)
      climate_matrix%applied%Precip(  m,j,i) = climate_matrix%PD_obs%Precip(  m,j,i)
      climate_matrix%applied%Hs(        j,i) = climate_matrix%PD_obs%Hs(        j,i)
      climate_matrix%applied%Wind_LR( m,j,i) = climate_matrix%PD_obs%Wind_LR( m,j,i)
      climate_matrix%applied%Wind_DU( m,j,i) = climate_matrix%PD_obs%Wind_DU( m,j,i)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_regional_PD_obs

! == Present-day observed climate plus a global temperature offset (de Boer et al., 2013)
! =======================================================================================

  SUBROUTINE run_climate_model_dT_glob( grid, ice, climate_matrix, region_name)
    ! Use the climate parameterisation from de Boer et al., 2013 (global temperature offset calculated with the inverse routine,
    ! plus a precipitation correction based on temperature + orography changes (NAM & EAS; Roe & Lindzen model), or only temperature (GRL & ANT).
    ! (for more details, see de Boer, B., van de Wal, R., Lourens, L. J., Bintanja, R., and Reerink, T. J.:
    ! A continuous simulation of global ice volume over the past 1 million years with 3-D ice-sheet models, Climate Dynamics 41, 1365-1384, 2013)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_climate_matrix_regional),  INTENT(INOUT) :: climate_matrix
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_dT_glob'
    INTEGER                                            :: i,j,m
    REAL(dp), DIMENSION(:,:  ), POINTER                :: dHs_dx_ref, dHs_dy_ref
    REAL(dp)                                           :: dT_lapse
    REAL(dp), DIMENSION(:,:,:), POINTER                :: Precip_RL_ref, Precip_RL_mod, dPrecip_RL
    INTEGER                                            :: wdHs_dx_ref, wdHs_dy_ref, wPrecip_RL_ref, wPrecip_RL_mod, wdPrecip_RL
    REAL(dp), PARAMETER                                :: P_offset = 0.008_dp       ! Normalisation term in precipitation anomaly to avoid divide-by-nearly-zero

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, dHs_dx_ref   , wdHs_dx_ref   )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, dHs_dy_ref   , wdHs_dy_ref   )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, Precip_RL_ref, wPrecip_RL_ref)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, Precip_RL_mod, wPrecip_RL_mod)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, dPrecip_RL   , wdPrecip_RL   )

    ! Get surface slopes for the PD_obs reference orography
    CALL ddx_a_to_a_2D( grid, climate_matrix%PD_obs%Hs, dHs_dx_ref)
    CALL ddy_a_to_a_2D( grid, climate_matrix%PD_obs%Hs, dHs_dy_ref)

    ! Temperature: constant lapse rate plus global offset
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      dT_lapse = (ice%Hs_a( j,i) - climate_matrix%PD_obs%Hs( j,i)) * C%constant_lapserate
      DO m = 1, 12
        climate_matrix%applied%T2m( m,j,i) = climate_matrix%PD_obs%T2m( m,j,i) + dT_lapse + forcing%dT_glob_inverse
      END DO

    END DO
    END DO
    CALL sync

    ! Precipitation:
    ! NAM & EAS: Roe&Lindzen model to account for changes in orography and temperature
    ! GRL & ANT: simple correction based on temperature alone

    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
      DO m = 1, 12

        CALL precipitation_model_Roe( climate_matrix%PD_obs%T2m(  m,j,i), dHs_dx_ref(   j,i), dHs_dy_ref(   j,i), climate_matrix%PD_obs%Wind_LR( m,j,i), climate_matrix%PD_obs%Wind_DU( m,j,i), Precip_RL_ref( m,j,i))
        CALL precipitation_model_Roe( climate_matrix%applied%T2m( m,j,i), ice%dHs_dx_a( j,i), ice%dHs_dy_a( j,i), climate_matrix%PD_obs%Wind_LR( m,j,i), climate_matrix%PD_obs%Wind_DU( m,j,i), Precip_RL_mod( m,j,i))
        dPrecip_RL( m,j,i) = MAX(0.01_dp, MIN( 2._dp, Precip_RL_mod( m,j,i) / Precip_RL_ref( m,j,i) ))

        climate_matrix%applied%Precip( m,j,i) = climate_matrix%PD_obs%Precip( m,j,i) * dPrecip_RL( m,j,i)

      END DO
      END DO
      END DO
      CALL sync

    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN

      CALL adapt_precip_CC(  grid, ice%Hs_a, climate_matrix%PD_obs%Hs, climate_matrix%PD_obs%T2m, climate_matrix%PD_obs%Precip, climate_matrix%applied%Precip, region_name)

    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wdHs_dx_ref)
    CALL deallocate_shared( wdHs_dy_ref)
    CALL deallocate_shared( wPrecip_RL_ref)
    CALL deallocate_shared( wPrecip_RL_mod)
    CALL deallocate_shared( wdPrecip_RL)

    ! Safety
    CALL check_for_NaN_dp_3D( climate_matrix%applied%T2m   , 'climate_matrix%applied%T2m'   )
    CALL check_for_NaN_dp_3D( climate_matrix%applied%Precip, 'climate_matrix%applied%Precip')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_dT_glob

! == Warm/cold climate matrix
! ===========================

  ! Climate matrix with PI + LGM snapshots, forced with CO2 (from record or from inverse routine) from Berends et al., 2018
  ! Generalised for different timeframes, L.B. Stap (2021)
  SUBROUTINE run_climate_model_matrix_warm_cold( grid, ice, SMB, climate_matrix, region_name, time)
    ! Use CO2 (either prescribed or inversely modelled) to force the 2-snapshot (PI-LGM) climate matrix (Berends et al., 2018)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_matrix_regional),  INTENT(INOUT) :: climate_matrix
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_matrix_warm_cold'
    INTEGER                                            :: i,j,m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Update forcing at model time
    CALL get_insolation_at_time( grid, time, climate_matrix%applied%Q_TOA)
    CALL update_CO2_at_model_time( time)

    ! Use the (CO2 + absorbed insolation)-based interpolation scheme for temperature
    CALL run_climate_model_matrix_warm_cold_temperature( grid, ice, SMB, climate_matrix, region_name)

    ! Use the (CO2 + ice-sheet geometry)-based interpolation scheme for precipitation
    CALL run_climate_model_matrix_warm_cold_precipitation( grid, ice, climate_matrix, region_name)

    ! Safety checks
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      ! Temperature errors
      IF     (climate_matrix%applied%T2m( m,j,i) < 0._dp) THEN
        CALL crash('negative temperature detected at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      ELSEIF (climate_matrix%applied%T2m( m,j,i) /= climate_matrix%applied%T2m( m,j,i)) THEN
        CALL crash('NaN temperature detected at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      END IF
      ! Precipitation errors
      IF     (climate_matrix%applied%Precip( m,j,i) <= 0._dp) THEN
        CALL crash('zero/negative precipitation detected at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      ELSEIF (climate_matrix%applied%Precip( m,j,i) /= climate_matrix%applied%Precip( m,j,i)) THEN
        CALL crash('NaN precipitation detected at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      END IF
      ! Temperature warnings
      IF     (climate_matrix%applied%T2m( m,j,i) < 150._dp) THEN
        CALL warning('excessively low temperature (< 150 K) detected at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      ELSEIF (climate_matrix%applied%T2m( m,j,i) > 350._dp) THEN
        CALL warning('excessively high temperature (> 350 K) detected at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      END IF
      ! Precipitation warnings
      IF     (climate_matrix%applied%Precip( m,j,i) > 10._dp) THEN
        CALL warning('excessively high precipitation (> 10 m/month) detected at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      END IF
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_matrix_warm_cold
  SUBROUTINE run_climate_model_matrix_warm_cold_temperature( grid, ice, SMB, climate_matrix, region_name)
    ! The (CO2 + absorbed insolation)-based matrix interpolation for temperature, from Berends et al. (2018)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_matrix_regional),  INTENT(INOUT) :: climate_matrix
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_matrix_warm_cold_temperature'
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: CO2, w_CO2
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  w_ins,  w_ins_smooth,  w_ice,  w_tot
    INTEGER                                            :: ww_ins, ww_ins_smooth, ww_ice, ww_tot
    REAL(dp)                                           :: w_ins_av
    REAL(dp), DIMENSION(:,:,:), POINTER                :: T_ref_GCM
    REAL(dp), DIMENSION(:,:  ), POINTER                :: Hs_GCM, lambda_GCM
    INTEGER                                            :: wT_ref_GCM, wHs_GCM, wlambda_GCM

    REAL(dp), PARAMETER                                :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]
    REAL(dp), PARAMETER                                :: P_offset = 0.008_dp       ! Normalisation term in precipitation anomaly to avoid divide-by-nearly-zero

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_ins,        ww_ins         )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_ins_smooth, ww_ins_smooth  )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_ice,        ww_ice         )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_tot,        ww_tot         )
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

    ! Find the interpolation weights based on absorbed insolation
    ! ===========================================================

    ! Calculate modelled absorbed insolation
    climate_matrix%applied%I_abs( :,grid%i1:grid%i2) = 0._dp
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      climate_matrix%applied%I_abs( j,i) = climate_matrix%applied%I_abs( j,i) + climate_matrix%applied%Q_TOA( m,j,i) * (1._dp - SMB%Albedo( m,j,i))  ! Berends et al., 2018 - Eq. 2
    END DO
    END DO
    END DO
    CALL sync

    ! Calculate weighting field
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      w_ins( j,i) = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (    climate_matrix%applied%I_abs(  j,i) -     climate_matrix%GCM_cold%I_abs( j,i)) / &  ! Berends et al., 2018 - Eq. 3
                                                           (    climate_matrix%GCM_warm%I_abs( j,i) -     climate_matrix%GCM_cold%I_abs( j,i)) ))
    END DO
    END DO
    CALL sync
    w_ins_av      = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (SUM(climate_matrix%applied%I_abs )      - SUM(climate_matrix%GCM_cold%I_abs)     ) / &
                                                           (SUM(climate_matrix%GCM_warm%I_abs)      - SUM(climate_matrix%GCM_cold%I_abs)     ) ))

    ! Smooth the weighting field
    w_ins_smooth( :,grid%i1:grid%i2) = w_ins( :,grid%i1:grid%i2)
    CALL smooth_Gaussian_2D( grid, w_ins_smooth, 200000._dp)

    ! Combine unsmoothed, smoothed, and regional average weighting fields (Berends et al., 2018, Eq. 4)
    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      w_ice( :,grid%i1:grid%i2) = (1._dp * w_ins(        :,grid%i1:grid%i2) + &
                                   3._dp * w_ins_smooth( :,grid%i1:grid%i2) + &
                                   3._dp * w_ins_av) / 7._dp
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      w_ice( :,grid%i1:grid%i2) = (1._dp * w_ins_smooth( :,grid%i1:grid%i2) + &
                                   6._dp * w_ins_av) / 7._dp
    END IF


    ! Combine interpolation weights from absorbed insolation and CO2 into the final weights fields
    ! Berends et al., 2018 - Eqs. 5, 9 with weights 0.5 for NAM & EAS, and 0.75 for ANT
    ! Generalised: "switch" between matrix method and glacial index method by altering C%climate_matrix_CO2vsice_<region>
    IF         (region_name == 'NAM') THEN
      w_tot( :,grid%i1:grid%i2) = (C%climate_matrix_CO2vsice_NAM * w_CO2) + ((1._dp - C%climate_matrix_CO2vsice_NAM) * w_ice( :,grid%i1:grid%i2))
    ELSEIF     (region_name == 'EAS') THEN
      w_tot( :,grid%i1:grid%i2) = (C%climate_matrix_CO2vsice_EAS * w_CO2) + ((1._dp - C%climate_matrix_CO2vsice_EAS) * w_ice( :,grid%i1:grid%i2))
    ELSEIF     (region_name == 'GRL') THEN
      w_tot( :,grid%i1:grid%i2) = (C%climate_matrix_CO2vsice_GRL * w_CO2) + ((1._dp - C%climate_matrix_CO2vsice_GRL) * w_ice( :,grid%i1:grid%i2))
    ELSEIF     (region_name == 'ANT') THEN
      w_tot( :,grid%i1:grid%i2) = (C%climate_matrix_CO2vsice_ANT * w_CO2) + ((1._dp - C%climate_matrix_CO2vsice_ANT) * w_ice( :,grid%i1:grid%i2))
    END IF

    ! Interpolate between the GCM snapshots
    ! =====================================

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Find matrix-interpolated orography, lapse rate, and temperature
      Hs_GCM(       j,i) = (w_tot( j,i) * climate_matrix%GCM_warm%Hs(         j,i)) + ((1._dp - w_tot( j,i)) * climate_matrix%GCM_cold%Hs(         j,i))  ! Berends et al., 2018 - Eq. 8
      lambda_GCM(   j,i) = (w_tot( j,i) * climate_matrix%GCM_warm%lambda(     j,i)) + ((1._dp - w_tot( j,i)) * climate_matrix%GCM_cold%lambda(     j,i))  ! Not listed in the article, shame on me!
      T_ref_GCM(  :,j,i) = (w_tot( j,i) * climate_matrix%GCM_warm%T2m_corr( :,j,i)) + ((1._dp - w_tot( j,i)) * climate_matrix%GCM_cold%T2m_corr( :,j,i))  ! Berends et al., 2018 - Eq. 6

      ! Adapt temperature to model orography using matrix-derived lapse-rate
      DO m = 1, 12
        climate_matrix%applied%T2m( m,j,i) = T_ref_GCM( m,j,i) - lambda_GCM( j,i) * (ice%Hs_a( j,i) - Hs_GCM( j,i))  ! Berends et al., 2018 - Eq. 11
      END DO

    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( ww_ins)
    CALL deallocate_shared( ww_ins_smooth)
    CALL deallocate_shared( ww_ice)
    CALL deallocate_shared( ww_tot)
    CALL deallocate_shared( wT_ref_GCM)
    CALL deallocate_shared( wHs_GCM)
    CALL deallocate_shared( wlambda_GCM)

    ! Safety
    CALL check_for_NaN_dp_3D( climate_matrix%applied%T2m   , 'climate_matrix%applied%T2m')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_matrix_warm_cold_temperature
  SUBROUTINE run_climate_model_matrix_warm_cold_precipitation( grid, ice, climate_matrix, region_name)
    ! The (CO2 + ice geometry)-based matrix interpolation for precipitation, from Berends et al. (2018)
    ! For NAM and EAS, this is based on local ice geometry and uses the Roe&Lindzen precipitation model for downscaling.
    ! For GRL and ANT, this is based on total ice volume,  and uses the simple CC   precipitation model for downscaling.
    ! The rationale for this difference is that glacial-interglacial differences in ice geometry are much more
    ! dramatic in NAM and EAS than they are in GRL and ANT.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_climate_matrix_regional),  INTENT(INOUT) :: climate_matrix
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_matrix_warm_cold_precipitation'
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  w_warm,  w_cold
    INTEGER                                            :: ww_warm, ww_cold
    REAL(dp)                                           :: w_tot
    REAL(dp), DIMENSION(:,:,:), POINTER                :: T_ref_GCM, P_ref_GCM
    REAL(dp), DIMENSION(:,:  ), POINTER                :: Hs_GCM
    INTEGER                                            :: wT_ref_GCM, wP_ref_GCM, wHs_GCM

    REAL(dp), PARAMETER                                :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_warm,         ww_warm        )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_cold,         ww_cold        )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, T_ref_GCM,      wT_ref_GCM     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, P_ref_GCM,      wP_ref_GCM     )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, Hs_GCM,         wHs_GCM        )

    ! Calculate interpolation weights based on ice geometry
    ! =====================================================

    ! First calculate the total ice volume term (second term in the equation)
    w_tot = MAX(-w_cutoff, MIN(1._dp + w_cutoff, (SUM(ice%Hs_a) - SUM(climate_matrix%GCM_warm%Hs)) / (SUM(climate_matrix%GCM_cold%Hs) - SUM(climate_matrix%GCM_warm%Hs)) ))

    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      ! Combine total + local ice thicness; Berends et al., 2018, Eq. 12

      ! Then the local ice thickness term
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny

        IF (climate_matrix%GCM_warm%Hs( j,i) < climate_matrix%GCM_PI%Hs( j,i) + 50._dp) THEN
          IF (climate_matrix%GCM_cold%Hs( j,i) < climate_matrix%GCM_PI%Hs( j,i) + 50._dp) THEN
            ! No ice in any GCM state. Use only total ice volume.
            w_cold( j,i) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, w_tot ))
            w_warm( j,i) = 1._dp - w_cold( j,i)
          ELSE
            ! No ice in warm climate, ice in cold climate. Linear inter- / extrapolation.
            w_cold( j,i) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, ((ice%Hs_a( j,i) - climate_matrix%GCM_PI%Hs( j,i)) / (climate_matrix%GCM_cold%Hs( j,i) - climate_matrix%GCM_PI%Hs( j,i))) * w_tot ))
            w_warm( j,i)  = 1._dp - w_cold( j,i)
          END IF
        ELSE
          ! Ice in both GCM states.  Linear inter- / extrapolation
          w_cold( j,i) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, ((ice%Hs_a( j,i) - climate_matrix%GCM_PI%Hs( j,i)) / (climate_matrix%GCM_cold%Hs( j,i) - climate_matrix%GCM_PI%Hs( j,i))) * w_tot ))
          w_warm( j,i)  = 1._dp - w_cold( j,i)
        END IF

      END DO
      END DO
      CALL sync

      w_cold( :,grid%i1:grid%i2) = w_cold( :,grid%i1:grid%i2) * w_tot

      ! Smooth the weighting field
      CALL smooth_Gaussian_2D( grid, w_cold, 200000._dp)

      w_warm( :,grid%i1:grid%i2) = 1._dp - w_cold( :,grid%i1:grid%i2)

    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use only total ice volume and CO2; Berends et al., 2018, Eq. 13

      w_cold( :,grid%i1:grid%i2) = w_tot
      w_warm( :,grid%i1:grid%i2) = 1._dp - w_cold( :,grid%i1:grid%i2)

    END IF

    IF (C%switch_glacial_index_precip) THEN ! If a glacial index is used for the precipitation forcing, it will only depend on CO2
      w_tot = 1._dp - (MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (forcing%CO2_obs - C%matrix_low_CO2_level) / (C%matrix_high_CO2_level - C%matrix_low_CO2_level) )) )
      w_cold( :,grid%i1:grid%i2) = w_tot
      w_warm( :,grid%i1:grid%i2) = 1._dp - w_cold( :,grid%i1:grid%i2)
    END IF

    ! Interpolate the GCM snapshots
    ! =============================

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      T_ref_GCM(  :,j,i) =      (w_warm( j,i) *      climate_matrix%GCM_warm%T2m_corr(    :,j,i))  + (w_cold( j,i) *     climate_matrix%GCM_cold%T2m_corr(    :,j,i))   ! Berends et al., 2018 - Eq. 6
      P_ref_GCM(  :,j,i) = EXP( (w_warm( j,i) *  LOG(climate_matrix%GCM_warm%Precip_corr( :,j,i))) + (w_cold( j,i) * LOG(climate_matrix%GCM_cold%Precip_corr( :,j,i)))) ! Berends et al., 2018 - Eq. 7
      Hs_GCM(       j,i) =      (w_warm( j,i) *      climate_matrix%GCM_warm%Hs(            j,i))  + (w_cold( j,i) *     climate_matrix%GCM_cold%Hs(            j,i))   ! Berends et al., 2018 - Eq. 8

    END DO
    END DO
    CALL sync

    ! Downscale precipitation from the coarse-resolution reference
    ! GCM orography to the fine-resolution ice-model orography
    ! ========================================================

    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      ! Use the Roe&Lindzen precipitation model to do this; Berends et al., 2018, Eqs. A3-A7
      CALL adapt_precip_Roe( grid, Hs_GCM,   T_ref_GCM,                  climate_matrix%PD_obs%Wind_LR, climate_matrix%PD_obs%Wind_DU, P_ref_GCM, &
                                   ice%Hs_a, climate_matrix%applied%T2m, climate_matrix%PD_obs%Wind_LR, climate_matrix%PD_obs%Wind_DU, climate_matrix%applied%Precip)
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use a simpler temperature-based correction; Berends et al., 2018, Eq. 14
      CALL adapt_precip_CC( grid, ice%Hs_a, Hs_GCM, T_ref_GCM, P_ref_GCM, climate_matrix%applied%Precip, region_name)
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( ww_warm)
    CALL deallocate_shared( ww_cold)
    CALL deallocate_shared( wT_ref_GCM)
    CALL deallocate_shared( wP_ref_GCM)
    CALL deallocate_shared( wHs_GCM)

    ! Safety
    CALL check_for_NaN_dp_3D( climate_matrix%applied%Precip, 'climate_matrix%applied%Precip')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_matrix_warm_cold_precipitation

  ! Initialising the climate matrix, containing all the global subclimates
  ! (PD observations and GCM snapshots) on their own lat-lon grids
  SUBROUTINE initialise_climate_matrix_global( climate_matrix_global)
    ! Allocate shared memory for the global climate matrix

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_climate_matrix_global), INTENT(INOUT) :: climate_matrix_global

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_matrix_global'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise the present-day observed global climate (e.g. ERA-40)
    CALL initialise_climate_PD_obs_global(   climate_matrix_global%PD_obs,   name = 'PD_obs')

    ! Initialise the (GCM-modelled) global climate snapshots
    CALL initialise_climate_snapshot_global( climate_matrix_global%GCM_PI,   name = 'GCM_PI',   nc_filename = C%filename_climate_snapshot_PI,   CO2 = 280._dp,                 orbit_time = 0._dp                   )
    CALL initialise_climate_snapshot_global( climate_matrix_global%GCM_warm, name = 'GCM_warm', nc_filename = C%filename_climate_snapshot_warm, CO2 = C%matrix_high_CO2_level, orbit_time = C%matrix_warm_orbit_time)
    CALL initialise_climate_snapshot_global( climate_matrix_global%GCM_cold, name = 'GCM_cold', nc_filename = C%filename_climate_snapshot_cold, CO2 = C%matrix_low_CO2_level,  orbit_time = C%matrix_cold_orbit_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_matrix_global
  SUBROUTINE initialise_climate_PD_obs_global( PD_obs, name)
    ! Allocate shared memory for the global PD observed climate data fields (stored in the climate matrix),
    ! read them from the specified NetCDF file (latter only done by master process).

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_climate_snapshot_global), INTENT(INOUT) :: PD_obs
    CHARACTER(LEN=*),                   INTENT(IN)    :: name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_PD_obs_global'

    ! Add routine to path
    CALL init_routine( routine_name)

    PD_obs%name = name
    PD_obs%netcdf%filename   = C%filename_PD_obs_climate

    ! General forcing info (not relevant for PD_obs, but needed so that the same mapping routines as for GCM snapshots can be used)
    CALL allocate_shared_dp_0D( PD_obs%CO2,        PD_obs%wCO2       )
    CALL allocate_shared_dp_0D( PD_obs%orbit_time, PD_obs%worbit_time)
    CALL allocate_shared_dp_0D( PD_obs%orbit_ecc,  PD_obs%worbit_ecc )
    CALL allocate_shared_dp_0D( PD_obs%orbit_obl,  PD_obs%worbit_obl )
    CALL allocate_shared_dp_0D( PD_obs%orbit_pre,  PD_obs%worbit_pre )

    ! Inquire if all required variables are present in the NetCDF file, and read the grid size.
    CALL allocate_shared_int_0D( PD_obs%nlon, PD_obs%wnlon)
    CALL allocate_shared_int_0D( PD_obs%nlat, PD_obs%wnlat)
    IF (par%master) CALL inquire_PD_obs_global_climate_file( PD_obs)
    CALL sync

    ! Allocate memory
    CALL allocate_shared_dp_1D( PD_obs%nlon,                  PD_obs%lon,         PD_obs%wlon        )
    CALL allocate_shared_dp_1D(              PD_obs%nlat,     PD_obs%lat,         PD_obs%wlat        )
    CALL allocate_shared_dp_2D( PD_obs%nlon, PD_obs%nlat,     PD_obs%Hs,          PD_obs%wHs         )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%T2m,         PD_obs%wT2m        )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%Precip,      PD_obs%wPrecip     )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%Wind_WE,     PD_obs%wWind_WE    )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%Wind_SN,     PD_obs%wWind_SN    )

    ! Read data from the NetCDF file
    IF (par%master) WRITE(0,*) '   Reading PD observed climate data from file ', TRIM(PD_obs%netcdf%filename), '...'
    IF (par%master) CALL read_PD_obs_global_climate_file( PD_obs)
    CALL sync

    ! Determine process domains
    CALL partition_list( PD_obs%nlon, par%i, par%n, PD_obs%i1, PD_obs%i2)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_PD_obs_global
  SUBROUTINE initialise_climate_snapshot_global( snapshot, name, nc_filename, CO2, orbit_time)
    ! Allocate shared memory for the data fields of a GCM snapshot (stored in the climate matrix),
    ! read them from the specified NetCDF file (latter only done by master process).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_climate_snapshot_global), INTENT(INOUT) :: snapshot
    CHARACTER(LEN=*),                   INTENT(IN)    :: name
    CHARACTER(LEN=*),                   INTENT(IN)    :: nc_filename
    REAL(dp),                           INTENT(IN)    :: CO2
    REAL(dp),                           INTENT(IN)    :: orbit_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_snapshot_global'
    INTEGER                                           :: i,j,m
    REAL(dp), PARAMETER                               :: Precip_minval = 1E-5_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Metadata
    snapshot%name            = name
    snapshot%netcdf%filename = nc_filename

    ! General forcing info
    CALL allocate_shared_dp_0D( snapshot%CO2,        snapshot%wCO2       )
    CALL allocate_shared_dp_0D( snapshot%orbit_time, snapshot%worbit_time)
    CALL allocate_shared_dp_0D( snapshot%orbit_ecc,  snapshot%worbit_ecc )
    CALL allocate_shared_dp_0D( snapshot%orbit_obl,  snapshot%worbit_obl )
    CALL allocate_shared_dp_0D( snapshot%orbit_pre,  snapshot%worbit_pre )

    snapshot%CO2        = CO2
    snapshot%orbit_time = orbit_time

    ! Inquire if all required variables are present in the NetCDF file, and read the grid size.
    CALL allocate_shared_int_0D( snapshot%nlon, snapshot%wnlon)
    CALL allocate_shared_int_0D( snapshot%nlat, snapshot%wnlat)
    IF (par%master) CALL inquire_GCM_global_climate_file( snapshot)
    CALL sync

    ! Allocate memory
    CALL allocate_shared_dp_1D( snapshot%nlon,                    snapshot%lon,     snapshot%wlon    )
    CALL allocate_shared_dp_1D(                snapshot%nlat,     snapshot%lat,     snapshot%wlat    )
    CALL allocate_shared_dp_2D( snapshot%nlon, snapshot%nlat,     snapshot%Hs,      snapshot%wHs     )
    CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, 12, snapshot%T2m,     snapshot%wT2m    )
    CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, 12, snapshot%Precip,  snapshot%wPrecip )
    CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, 12, snapshot%Wind_WE, snapshot%wWind_WE)
    CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, 12, snapshot%Wind_SN, snapshot%wWind_SN)

    ! Read data from the NetCDF file
    IF (par%master) WRITE(0,*) '   Reading GCM snapshot ', TRIM(snapshot%name), ' from file ', TRIM(snapshot%netcdf%filename), '...'
    IF (par%master) CALL read_GCM_global_climate_file( snapshot)
    CALL sync

    ! Determine process domains
    CALL partition_list( snapshot%nlon, par%i, par%n, snapshot%i1, snapshot%i2)

    ! Very rarely zero precipitation can occur in GCM snapshots, which gives problems with the matrix interpolation. Fix this.
    DO i = snapshot%i1, snapshot%i2
    DO j = 1, snapshot%nlat
    DO m = 1, 12
      snapshot%Precip( i,j,m) = MAX( Precip_minval, snapshot%Precip( i,j,m))
    END DO
    END DO
    END DO
    CALL sync

    ! Safety checks
    DO i = snapshot%i1, snapshot%i2
    DO j = 1, snapshot%nlat
    DO m = 1, 12
      ! Temperature errors
      IF     (snapshot%T2m( i,j,m) < 0._dp) THEN
        CALL crash('negative temperature detected in snapshot ' // TRIM(snapshot%name) // ' at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      ELSEIF (snapshot%T2m( i,j,m) /= snapshot%T2m( i,j,m)) THEN
        CALL crash('NaN temperature detected in snapshot ' // TRIM(snapshot%name) // ' at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      END IF
      ! Precipitation errors
      IF     (snapshot%Precip( i,j,m) <= 0._dp) THEN
        CALL crash('zero/negative precipitation detected in snapshot ' // TRIM(snapshot%name) // ' at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      ELSEIF (snapshot%Precip( i,j,m) /= snapshot%Precip( i,j,m)) THEN
        CALL crash('NaN precipitation detected in snapshot ' // TRIM(snapshot%name) // ' at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      END IF
      ! Temperature warnings
      IF     (snapshot%T2m( i,j,m) < 150._dp) THEN
        CALL warning('excessively low temperature (< 150 K) detected in snapshot ' // TRIM(snapshot%name) // ' at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      ELSEIF (snapshot%T2m( i,j,m) > 350._dp) THEN
        CALL warning('excessively high temperature (< 350 K) detected in snapshot ' // TRIM(snapshot%name) // ' at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      END IF
      ! Precipitation warnings
      IF     (snapshot%Precip( i,j,m) > 10._dp) THEN
        CALL warning('excessively high precipitation (> 10 m/month) detected in snapshot ' // TRIM(snapshot%name) // ' at [{int_01},{int_02},{int_03}]', int_01 = m, int_02 = j, int_03 = i)
      END IF

    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_snapshot_global

  ! Initialising the region-specific climate model, containing all the subclimates
  ! (PD observations, GCM snapshots and the applied climate) on the model grid
  SUBROUTINE initialise_climate_matrix_regional( region, climate_matrix_global)
    ! Allocate shared memory for the regional climate models, containing the PD observed,
    ! GCM snapshots and applied climates as "subclimates"

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_climate_matrix_global),    INTENT(IN)    :: climate_matrix_global

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_matrix_regional'
    INTEGER                                            :: i,j,m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory for the regional ERA40 climate and the final applied climate
    CALL allocate_climate_snapshot_regional( region%grid, region%climate_matrix%PD_obs,   name = 'PD_obs'  )
    CALL allocate_climate_snapshot_regional( region%grid, region%climate_matrix%GCM_PI,   name = 'GCM_PI'  )
    CALL allocate_climate_snapshot_regional( region%grid, region%climate_matrix%GCM_warm, name = 'GCM_warm')
    CALL allocate_climate_snapshot_regional( region%grid, region%climate_matrix%GCM_cold, name = 'GCM_cold')
    CALL allocate_climate_snapshot_regional( region%grid, region%climate_matrix%applied,  name = 'applied' )

    ! Map the snapshots from global lat/lon-grid to model x/y-grid
    CALL map_snapshot_to_grid( region%grid, climate_matrix_global%PD_obs,   region%climate_matrix%PD_obs  )
    CALL map_snapshot_to_grid( region%grid, climate_matrix_global%GCM_PI,   region%climate_matrix%GCM_PI  )
    CALL map_snapshot_to_grid( region%grid, climate_matrix_global%GCM_warm, region%climate_matrix%GCM_warm)
    CALL map_snapshot_to_grid( region%grid, climate_matrix_global%GCM_cold, region%climate_matrix%GCM_cold)

    ! Right now, no wind is read from GCM output; just use PD observations everywhere
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
    DO m = 1, 12
      region%climate_matrix%GCM_warm%Wind_WE( m,j,i) = region%climate_matrix%PD_obs%Wind_WE( m,j,i)
      region%climate_matrix%GCM_warm%Wind_SN( m,j,i) = region%climate_matrix%PD_obs%Wind_SN( m,j,i)
      region%climate_matrix%GCM_cold%Wind_WE( m,j,i) = region%climate_matrix%PD_obs%Wind_WE( m,j,i)
      region%climate_matrix%GCM_cold%Wind_SN( m,j,i) = region%climate_matrix%PD_obs%Wind_SN( m,j,i)
    END DO
    END DO
    END DO
    CALL sync

    ! Calculate spatially variable lapse rate
    region%climate_matrix%GCM_warm%lambda( :,region%grid%i1:region%grid%i2) = 0.008_dp
    CALL sync
    IF     (region%name == 'NAM' .OR. region%name == 'EAS') THEN
      CALL initialise_snapshot_spatially_variable_lapserate( region%grid, region%climate_matrix%GCM_PI, region%climate_matrix%GCM_cold)
    ELSEIF (region%name == 'GLR' .OR. region%name == 'ANT') THEN
      region%climate_matrix%GCM_cold%lambda( :,region%grid%i1:region%grid%i2) = 0.008_dp
      CALL sync
    END IF

    ! Calculate and apply GCM bias correction
    CALL calculate_GCM_bias( region%grid, region%climate_matrix)
    CALL correct_GCM_bias(   region%grid, region%climate_matrix, region%climate_matrix%GCM_warm, do_correct_bias = C%climate_matrix_biascorrect_warm)
    CALL correct_GCM_bias(   region%grid, region%climate_matrix, region%climate_matrix%GCM_cold, do_correct_bias = C%climate_matrix_biascorrect_cold)

    ! Get reference absorbed insolation for the GCM snapshots
    CALL initialise_snapshot_absorbed_insolation( region%grid, region%climate_matrix%GCM_warm, region%name, region%mask_noice)
    CALL initialise_snapshot_absorbed_insolation( region%grid, region%climate_matrix%GCM_cold, region%name, region%mask_noice)

    ! Initialise applied climate with present-day observations
    DO i = region%grid%i2, region%grid%i2
    DO j = 1, region%grid%ny
    DO m = 1, 12
      region%climate_matrix%applied%T2m(     m,j,i) = region%climate_matrix%PD_obs%T2m(     m,j,i)
      region%climate_matrix%applied%Precip(  m,j,i) = region%climate_matrix%PD_obs%Precip(  m,j,i)
      region%climate_matrix%applied%Hs(        j,i) = region%climate_matrix%PD_obs%Hs(        j,i)
      region%climate_matrix%applied%Wind_LR( m,j,i) = region%climate_matrix%PD_obs%Wind_LR( m,j,i)
      region%climate_matrix%applied%Wind_DU( m,j,i) = region%climate_matrix%PD_obs%Wind_DU( m,j,i)
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_matrix_regional
  SUBROUTINE allocate_climate_snapshot_regional( grid, climate, name)
    ! Allocate shared memory for a "subclimate" (PD observed, GCM snapshot or applied climate) on the grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_climate_snapshot_regional), INTENT(INOUT) :: climate
    CHARACTER(LEN=*),                     INTENT(IN)    :: name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_climate_snapshot_regional'

    ! Add routine to path
    CALL init_routine( routine_name)

    climate%name = name

    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, climate%Hs,             climate%wHs            )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%T2m,            climate%wT2m           )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%Precip,         climate%wPrecip        )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%Wind_WE,        climate%wWind_WE       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%Wind_SN,        climate%wWind_SN       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%Wind_LR,        climate%wWind_LR       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%Wind_DU,        climate%wWind_DU       )

    CALL allocate_shared_dp_0D(                       climate%CO2,            climate%wCO2           )
    CALL allocate_shared_dp_0D(                       climate%orbit_time,     climate%worbit_time    )
    CALL allocate_shared_dp_0D(                       climate%orbit_ecc,      climate%worbit_ecc     )
    CALL allocate_shared_dp_0D(                       climate%orbit_obl,      climate%worbit_obl     )
    CALL allocate_shared_dp_0D(                       climate%orbit_pre,      climate%worbit_pre     )
    CALL allocate_shared_dp_0D(                       climate%sealevel,       climate%wsealevel      )

    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, climate%lambda,         climate%wlambda        )

    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%T2m_corr,       climate%wT2m_corr      )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%Precip_corr,    climate%wPrecip_corr   )

    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%Q_TOA,          climate%wQ_TOA         )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%Albedo,         climate%wAlbedo        )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, climate%I_abs,          climate%wI_abs         )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_climate_snapshot_regional
  SUBROUTINE calculate_GCM_bias( grid, climate_matrix)
    ! Calculate the GCM bias in temperature and precipitation
    !
    ! Account for the fact that the GCM PI snapshot has a lower resolution, and therefore
    ! a different surface elevation than the PD observed climatology!

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),  INTENT(INOUT) :: climate_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calculate_GCM_bias'
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: T2m_SL_GCM, T2m_SL_obs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%GCM_bias_T2m,    climate_matrix%wGCM_bias_T2m   )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%GCM_bias_Precip, climate_matrix%wGCM_bias_Precip)

    ! Calculate bias
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12

      ! Scale modelled and observed temperature to sea level using a constant lapse rate
      T2m_SL_obs = climate_matrix%PD_obs%T2m( m,j,i) + climate_matrix%PD_obs%Hs( j,i) * C%constant_lapserate
      T2m_SL_GCM = climate_matrix%GCM_PI%T2m( m,j,i) + climate_matrix%GCM_PI%Hs( j,i) * C%constant_lapserate

      ! Calculate temperature bias
      climate_matrix%GCM_bias_T2m(    m,j,i) = T2m_SL_GCM - T2m_SL_obs

      ! Calculate precipitation bias
      climate_matrix%GCM_bias_Precip( m,j,i) = climate_matrix%GCM_PI%Precip( m,j,i) / climate_matrix%PD_obs%Precip( m,j,i)

    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_GCM_bias
  SUBROUTINE correct_GCM_bias( grid, climate_matrix, climate, do_correct_bias)
    ! Calculate bias-corrected climate for this snapshot

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),   INTENT(IN)    :: climate_matrix
    TYPE(type_climate_snapshot_regional), INTENT(INOUT) :: climate
    LOGICAL,                              INTENT(IN)    :: do_correct_bias

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'correct_GCM_bias'
    INTEGER                                             :: i,j,m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no bias correction should be applied, set T2m_corr = T2m and Precip_corr = Precip
    IF (.NOT. do_correct_bias) THEN
      climate%T2m_corr(    :,:,grid%i1:grid%i2) = climate%T2m(    :,:,grid%i1:grid%i2)
      climate%Precip_corr( :,:,grid%i1:grid%i2) = climate%Precip( :,:,grid%i1:grid%i2)
      CALL sync
      RETURN
    END IF

    ! Apply bias correction
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12

      ! Temperature
      climate%T2m_corr(    m,j,i) = climate%T2m(    m,j,i) - climate_matrix%GCM_bias_T2m(    m,j,i)

      ! Precipitation
      climate%Precip_corr( m,j,i) = climate%Precip( m,j,i) / climate_matrix%GCM_bias_Precip( m,j,i)

    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE correct_GCM_bias
  SUBROUTINE initialise_snapshot_spatially_variable_lapserate( grid, climate_PI, climate)
    ! Calculate the spatially variable lapse-rate (for non-PI GCM climates; see Berends et al., 2018)
    ! Only meaningful for climates where there is ice (LGM, M2_Medium, M2_Large),
    ! and only intended for North America and Eurasia

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_climate_snapshot_regional), INTENT(IN)    :: climate_PI
    TYPE(type_climate_snapshot_regional), INTENT(INOUT) :: climate

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_snapshot_spatially_variable_lapserate'
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

      IF (climate%Hs( j,i) > climate_PI%Hs( j,i) + 100._dp) THEN
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
        dT_mean_nonice = dT_mean_nonice + climate%T2m( m,j,i) - climate_PI%T2m( m,j,i)
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
          climate%lambda( j,i) = climate%lambda( j,i) + 1/12._dp * MAX(lambda_min, MIN(lambda_max, &
            -(climate%T2m( m,j,i) - (climate_PI%T2m( m,j,i) + dT_mean_nonice)) / (climate%Hs( j,i) - climate_PI%Hs( j,i))))
        END DO

        lambda_mean_ice = lambda_mean_ice + climate%lambda( j,i)
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
      IF (mask_calc_lambda( j,i) == 0) climate%lambda( j,i) = lambda_mean_ice
    END DO
    END DO
    CALL sync

    ! Smooth the lapse rate field with a 160 km Gaussian filter
    CALL smooth_Gaussian_2D( grid, climate%lambda, 160000._dp)

    ! Normalise the entire region to a mean lapse rate of 8 K /km
    climate%lambda( :,grid%i1:grid%i2) = climate%lambda( :,grid%i1:grid%i2) * (0.008_dp / lambda_mean_ice)

    ! Clean up after yourself
    CALl deallocate_shared( wmask_calc_lambda)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_snapshot_spatially_variable_lapserate
  SUBROUTINE initialise_snapshot_absorbed_insolation( grid, climate, region_name, mask_noice)
    ! Calculate the yearly absorbed insolation for this (regional) GCM snapshot, to be used in the matrix interpolation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_climate_snapshot_regional), INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name
    INTEGER,  DIMENSION(:,:  ),           INTENT(IN)    :: mask_noice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_snapshot_absorbed_insolation'
    INTEGER                                             :: i,j,m
    TYPE(type_ice_model)                                :: ice_dummy
    TYPE(type_climate_matrix_regional)                  :: climate_dummy
    TYPE(type_SMB_model)                                :: SMB_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get insolation at the desired time from the insolation NetCDF file
    ! ==================================================================

    CALL get_insolation_at_time( grid, climate%orbit_time, climate%Q_TOA)

    ! Create temporary "dummy" climate, ice & SMB data structures,
    ! so we can run the SMB model and determine the reference albedo field
    ! ====================================================================

    ! Climate
    ! =======

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_dummy%applied%T2m,    climate_dummy%applied%wT2m)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_dummy%applied%Precip, climate_dummy%applied%wPrecip)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_dummy%applied%Q_TOA,  climate_dummy%applied%wQ_TOA)

    ! Copy climate fields
    climate_dummy%applied%T2m(    :,:,grid%i1:grid%i2) = climate%T2m_corr(    :,:,grid%i1:grid%i2)
    climate_dummy%applied%Precip( :,:,grid%i1:grid%i2) = climate%Precip_corr( :,:,grid%i1:grid%i2)
    climate_dummy%applied%Q_TOA(  :,:,grid%i1:grid%i2) = climate%Q_TOA(       :,:,grid%i1:grid%i2)

    ! Ice
    ! ===

    CALL allocate_shared_int_2D(    grid%ny, grid%nx, ice_dummy%mask_ocean_a   , ice_dummy%wmask_ocean_a   )
    CALL allocate_shared_int_2D(    grid%ny, grid%nx, ice_dummy%mask_ice_a     , ice_dummy%wmask_ice_a     )
    CALL allocate_shared_int_2D(    grid%ny, grid%nx, ice_dummy%mask_shelf_a   , ice_dummy%wmask_shelf_a   )

    ! Fill in masks for the SMB model
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (climate%Hs( j,i) == MINVAL(climate%Hs)) THEN
        ice_dummy%mask_ocean_a( j,i) = 1
      ELSE
        ice_dummy%mask_ocean_a( j,i) = 0
      END IF

      IF (climate%Hs( j,i) > 100._dp .AND. SUM(climate%T2m( :,j,i)) / 12._dp < 0._dp) THEN
        ice_dummy%mask_ice_a(   j,i) = 1
      ELSE
        ice_dummy%mask_ice_a(   j,i) = 0
      END IF

      ! mask_shelf is used in the SMB model only to find open ocean; since mask_ocean
      ! in this case already marks only open ocean, no need to look for shelves
      ice_dummy%mask_shelf_a( j,i) = 0

    END DO
    END DO
    CALL sync

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

    ! Run the SMB model for 10 years for this particular climate
    ! (experimentally determined to be long enough to converge)
    DO i = 1, 10
      CALL run_SMB_model( grid, ice_dummy, climate_dummy, 0._dp, SMB_dummy, mask_noice)
    END DO
    CALL sync

    ! Copy the resulting albedo to the climate
    climate%Albedo( :,:,grid%i1:grid%i2) = SMB_dummy%Albedo( :,:,grid%i1:grid%i2)
    CALL sync

    ! Calculate yearly total absorbed insolation
    climate%I_abs( :,grid%i1:grid%i2) = 0._dp
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      climate%I_abs( j,i) = climate%I_abs( j,i) + climate%Q_TOA( m,j,i) * (1._dp - climate%Albedo( m,j,i))
    END DO
    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( ice_dummy%wmask_ocean_a)
    CALL deallocate_shared( ice_dummy%wmask_ice_a)
    CALL deallocate_shared( ice_dummy%wmask_shelf_a)
    CALL deallocate_shared( climate_dummy%applied%wT2m)
    CALL deallocate_shared( climate_dummy%applied%wPrecip)
    CALL deallocate_shared( climate_dummy%applied%wQ_TOA)
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

  END SUBROUTINE initialise_snapshot_absorbed_insolation

! == Directly prescribed global climate
! =====================================

  SUBROUTINE run_climate_model_direct_climate_global( grid, clim_glob, climate_matrix, time)
    ! Run the regional climate model
    !
    ! Use a directly prescribed global climate

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_direct_climate_forcing_global), INTENT(INOUT) :: clim_glob
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix
    REAL(dp),                                 INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_direct_climate_global'
    REAL(dp)                                           :: wt0, wt1
    INTEGER                                            :: i,j,m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_global') THEN
      CALL crash('choice_climate_model should be "direct_global"!')
    END IF

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time < climate_matrix%direct%t0 .OR. time > climate_matrix%direct%t1) THEN

      ! Find and read the two global time frames
      CALL update_direct_global_climate_timeframes_from_file( clim_glob, time)

      ! Map the result to the regional grid
      CALL map_glob_to_grid_2D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%Hs0,      climate_matrix%direct%Hs0     )
      CALL map_glob_to_grid_3D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%T2m0,     climate_matrix%direct%T2m0    )
      CALL map_glob_to_grid_3D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%Precip0,  climate_matrix%direct%Precip0 )
      CALL map_glob_to_grid_3D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%Wind_WE0, climate_matrix%direct%Wind_WE0)
      CALL map_glob_to_grid_3D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%Wind_SN0, climate_matrix%direct%Wind_SN0)

      CALL map_glob_to_grid_2D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%Hs1,      climate_matrix%direct%Hs1     )
      CALL map_glob_to_grid_3D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%T2m1,     climate_matrix%direct%T2m1    )
      CALL map_glob_to_grid_3D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%Precip1,  climate_matrix%direct%Precip1 )
      CALL map_glob_to_grid_3D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%Wind_WE1, climate_matrix%direct%Wind_WE1)
      CALL map_glob_to_grid_3D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%Wind_SN1, climate_matrix%direct%Wind_SN1)

      ! Rotate the wind fields
      CALL rotate_wind_to_model_grid( grid, climate_matrix%direct%Wind_WE0, climate_matrix%direct%Wind_SN0, climate_matrix%direct%Wind_LR0, climate_matrix%direct%Wind_DU0)
      CALL rotate_wind_to_model_grid( grid, climate_matrix%direct%Wind_WE1, climate_matrix%direct%Wind_SN1, climate_matrix%direct%Wind_LR1, climate_matrix%direct%Wind_DU1)

      ! Update time stamps of regional climate forcing
      IF (par%master) THEN
        climate_matrix%direct%t0 = clim_glob%t0
        climate_matrix%direct%t1 = clim_glob%t1
      END IF
      CALL sync

    END IF ! IF (time >= climate_matrix%direct%t0 .AND. time <= climate_matrix%direct%t1) THEN

    ! Interpolate the two timeframes in time
    wt0 = (climate_matrix%direct%t1 - time) / (climate_matrix%direct%t1 - climate_matrix%direct%t0)
    wt1 = 1._dp - wt0

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      climate_matrix%applied%Hs(        j,i) = (wt0 * climate_matrix%direct%Hs0(        j,i)) + (wt1 * climate_matrix%direct%Hs1(        j,i))
      climate_matrix%applied%T2m(     m,j,i) = (wt0 * climate_matrix%direct%T2m0(     m,j,i)) + (wt1 * climate_matrix%direct%T2m1(     m,j,i))
      climate_matrix%applied%Precip(  m,j,i) = (wt0 * climate_matrix%direct%Precip0(  m,j,i)) + (wt1 * climate_matrix%direct%Precip1(  m,j,i))
      climate_matrix%applied%Wind_LR( m,j,i) = (wt0 * climate_matrix%direct%Wind_LR0( m,j,i)) + (wt1 * climate_matrix%direct%Wind_LR1( m,j,i))
      climate_matrix%applied%Wind_DU( m,j,i) = (wt0 * climate_matrix%direct%Wind_DU0( m,j,i)) + (wt1 * climate_matrix%direct%Wind_DU1( m,j,i))
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_direct_climate_global
  SUBROUTINE update_direct_global_climate_timeframes_from_file( clim_glob, time)
    ! Read the NetCDF file containing the global climate forcing data. Only read the time
    ! frames enveloping the current coupling timestep to save on memory usage.

    IMPLICIT NONE

    TYPE(type_direct_climate_forcing_global), INTENT(INOUT) :: clim_glob
    REAL(dp),                                 INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_direct_global_climate_timeframes_from_file'
    INTEGER                                            :: ti0, ti1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_global') THEN
      CALL crash('choice_climate_model should be "direct_global"!')
    END IF

    ! Check if data for model time is available
    IF (time < clim_glob%time(1)) THEN
      CALL crash('queried time out of range of direct transient climate forcing!')
    END IF

    ! Find time indices to be read
    IF (par%master) THEN

      IF     (time < clim_glob%time( 1)) THEN

        CALL warning('using constant start-of-record climate when extrapolating!')
        ti0 = 1
        ti1 = 1
        clim_glob%t0 = clim_glob%time( ti0) - 1._dp
        clim_glob%t1 = clim_glob%time( ti1)

      ELSEIF (time <= clim_glob%time( clim_glob%nyears)) THEN

        ti1 = 1
        DO WHILE (clim_glob%time( ti1) < time)
          ti1 = ti1 + 1
        END DO
        ti0 = ti1 - 1

        IF (ti0 == 0) THEN
          ti0 = 1
          ti1 = 2
        ELSEIF (ti1 == clim_glob%nyears) THEN
          ti0 = clim_glob%nyears - 1
          ti1 = ti0 + 1
        END IF

        clim_glob%t0 = clim_glob%time( ti0)
        clim_glob%t1 = clim_glob%time( ti1)

      ELSE ! IF     (time < clim_glob%time( 1)) THEN

        CALL warning('using constant end-of-record climate when extrapolating!')
        ti0 = clim_glob%nyears
        ti1 = clim_glob%nyears
        clim_glob%t0 = clim_glob%time( ti0) - 1._dp
        clim_glob%t1 = clim_glob%time( ti1)

      END IF ! IF     (time < clim_glob%time( 1)) THEN

    END IF ! IF (par%master) THEN

    ! Read new global climate fields from the NetCDF file
    IF (par%master) CALL read_direct_global_climate_file_timeframes( clim_glob, ti0, ti1)
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_direct_global_climate_timeframes_from_file
  SUBROUTINE initialise_climate_model_direct_climate_global( clim_glob)
    ! Initialise the global climate model
    !
    ! Use a directly prescribed global climate

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_direct_climate_forcing_global), INTENT(INOUT) :: clim_glob

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_model_direct_climate_global'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_global') THEN
      CALL crash('choice_climate_model should be "direct_global"!')
    END IF

    ! The times at which we have climate fields from input, between which we'll interpolate
    ! to find the climate at model time (t0 <= model_time <= t1)

    CALL allocate_shared_dp_0D( clim_glob%t0, clim_glob%wt0)
    CALL allocate_shared_dp_0D( clim_glob%t1, clim_glob%wt1)

    IF (par%master) THEN
      ! Give impossible values to timeframes, so that the first call to run_climate_model_direct_climate_global
      ! is guaranteed to first read two new timeframes from the NetCDF file
      clim_glob%t0 = C%start_time_of_run - 100._dp
      clim_glob%t1 = C%start_time_of_run - 90._dp
    END IF ! IF (par%master) THEN
    CALL sync

    IF (par%master) WRITE(0,*) '  Initialising direct global climate forcing from ', TRIM( C%filename_direct_global_climate), '...'

    ! Inquire into the direct global cliamte forcing netcdf file
    CALL allocate_shared_int_0D( clim_glob%nyears, clim_glob%wnyears)
    CALL allocate_shared_int_0D( clim_glob%nlon,   clim_glob%wnlon  )
    CALL allocate_shared_int_0D( clim_glob%nlat,   clim_glob%wnlat  )

    clim_glob%netcdf%filename = C%filename_direct_global_climate

    IF (par%master) CALL inquire_direct_global_climate_forcing_file( clim_glob)
    CALL sync

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( clim_glob%nyears, clim_glob%time, clim_glob%wtime)
    CALL allocate_shared_dp_1D( clim_glob%nlon,   clim_glob%lon,  clim_glob%wlon )
    CALL allocate_shared_dp_1D( clim_glob%nlat,   clim_glob%lat,  clim_glob%wlat )

    CALL allocate_shared_dp_2D( clim_glob%nlon, clim_glob%nlat,     clim_glob%Hs0,      clim_glob%wHs0     )
    CALL allocate_shared_dp_2D( clim_glob%nlon, clim_glob%nlat,     clim_glob%Hs1,      clim_glob%wHs1     )
    CALL allocate_shared_dp_3D( clim_glob%nlon, clim_glob%nlat, 12, clim_glob%T2m0,     clim_glob%wT2m0    )
    CALL allocate_shared_dp_3D( clim_glob%nlon, clim_glob%nlat, 12, clim_glob%T2m1,     clim_glob%wT2m1    )
    CALL allocate_shared_dp_3D( clim_glob%nlon, clim_glob%nlat, 12, clim_glob%Precip0,  clim_glob%wPrecip0 )
    CALL allocate_shared_dp_3D( clim_glob%nlon, clim_glob%nlat, 12, clim_glob%Precip1,  clim_glob%wPrecip1 )
    CALL allocate_shared_dp_3D( clim_glob%nlon, clim_glob%nlat, 12, clim_glob%Wind_WE0, clim_glob%wWind_WE0)
    CALL allocate_shared_dp_3D( clim_glob%nlon, clim_glob%nlat, 12, clim_glob%Wind_WE1, clim_glob%wWind_WE1)
    CALL allocate_shared_dp_3D( clim_glob%nlon, clim_glob%nlat, 12, clim_glob%Wind_SN0, clim_glob%wWind_SN0)
    CALL allocate_shared_dp_3D( clim_glob%nlon, clim_glob%nlat, 12, clim_glob%Wind_SN1, clim_glob%wWind_SN1)

    ! Read time and grid data
    IF (par%master) CALL read_direct_global_climate_file_time_latlon( clim_glob)
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_direct_climate_global
  SUBROUTINE initialise_climate_model_direct_climate_global_regional( grid, climate_matrix)
    ! Initialise the regional climate model
    !
    ! Use a directly prescribed global climate

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_model_direct_climate_global_regional'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_global') THEN
      CALL crash('choice_climate_model should be "direct_global"!')
    END IF

    ! The times at which we have climate fields from input, between which we'll interpolate
    ! to find the climate at model time (t0 <= model_time <= t1)

    CALL allocate_shared_dp_0D( climate_matrix%direct%t0, climate_matrix%direct%wt0)
    CALL allocate_shared_dp_0D( climate_matrix%direct%t1, climate_matrix%direct%wt1)

    IF (par%master) THEN
      ! Give impossible values to timeframes, so that the first call to run_climate_model_direct_climate_global
      ! is guaranteed to first read two new timeframes from the NetCDF file
      climate_matrix%direct%t0 = C%start_time_of_run - 100._dp
      climate_matrix%direct%t1 = C%start_time_of_run - 90._dp
    END IF ! IF (par%master) THEN
    CALL sync

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, climate_matrix%direct%Hs0,      climate_matrix%direct%wHs0     )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, climate_matrix%direct%Hs1,      climate_matrix%direct%wHs1     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%T2m0,     climate_matrix%direct%wT2m0    )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%T2m1,     climate_matrix%direct%wT2m1    )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%Precip0,  climate_matrix%direct%wPrecip0 )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%Precip1,  climate_matrix%direct%wPrecip1 )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%Wind_WE0, climate_matrix%direct%wWind_WE0)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%Wind_WE1, climate_matrix%direct%wWind_WE1)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%Wind_SN0, climate_matrix%direct%wWind_SN0)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%Wind_SN1, climate_matrix%direct%wWind_SN1)

    ! Lastly, allocate memory for the "applied" snapshot
    CALL allocate_climate_snapshot_regional( grid, climate_matrix%applied, name = 'applied')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_direct_climate_global_regional

! == Directly prescribed regional climate
! =======================================

  SUBROUTINE run_climate_model_direct_climate_regional( grid, climate_matrix, time)
    ! Run the regional climate model
    !
    ! Use a directly prescribed regional climate

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix
    REAL(dp),                                 INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_direct_climate_regional'
    REAL(dp)                                           :: wt0, wt1
    INTEGER                                            :: i,j,m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_regional') THEN
      CALL crash('choice_climate_model should be "direct_regional"!')
    END IF

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time < climate_matrix%direct%t0 .OR. time > climate_matrix%direct%t1) THEN

      ! Find and read the two global time frames
      CALL update_direct_regional_climate_timeframes_from_file( grid, climate_matrix%direct, time)

      ! Rotate the wind fields
      CALL rotate_wind_to_model_grid( grid, climate_matrix%direct%Wind_WE0, climate_matrix%direct%Wind_SN0, climate_matrix%direct%Wind_LR0, climate_matrix%direct%Wind_DU0)
      CALL rotate_wind_to_model_grid( grid, climate_matrix%direct%Wind_WE1, climate_matrix%direct%Wind_SN1, climate_matrix%direct%Wind_LR1, climate_matrix%direct%Wind_DU1)

    END IF ! IF (time >= climate_matrix%direct%t0 .AND. time <= climate_matrix%direct%t1) THEN

    ! Interpolate the two timeframes in time
    wt0 = (climate_matrix%direct%t1 - time) / (climate_matrix%direct%t1 - climate_matrix%direct%t0)
    wt1 = 1._dp - wt0

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      climate_matrix%applied%Hs(        j,i) = (wt0 * climate_matrix%direct%Hs0(        j,i)) + (wt1 * climate_matrix%direct%Hs1(        j,i))
      climate_matrix%applied%T2m(     m,j,i) = (wt0 * climate_matrix%direct%T2m0(     m,j,i)) + (wt1 * climate_matrix%direct%T2m1(     m,j,i))
      climate_matrix%applied%Precip(  m,j,i) = (wt0 * climate_matrix%direct%Precip0(  m,j,i)) + (wt1 * climate_matrix%direct%Precip1(  m,j,i))
      climate_matrix%applied%Wind_LR( m,j,i) = (wt0 * climate_matrix%direct%Wind_LR0( m,j,i)) + (wt1 * climate_matrix%direct%Wind_LR1( m,j,i))
      climate_matrix%applied%Wind_DU( m,j,i) = (wt0 * climate_matrix%direct%Wind_DU0( m,j,i)) + (wt1 * climate_matrix%direct%Wind_DU1( m,j,i))
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_direct_climate_regional
  SUBROUTINE update_direct_regional_climate_timeframes_from_file( grid, clim_reg, time)
    ! Read the NetCDF file containing the regional climate forcing data. Only read the time
    ! frames enveloping the current coupling timestep to save on memory usage.

    IMPLICIT NONE

    TYPE(type_grid),                            INTENT(IN)    :: grid
    TYPE(type_direct_climate_forcing_regional), INTENT(INOUT) :: clim_reg
    REAL(dp),                                   INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_direct_regional_climate_timeframes_from_file'
    INTEGER                                            :: ti0, ti1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_regional') THEN
      CALL crash('choice_climate_model should be "direct_regional"!')
    END IF

    ! Check if data for model time is available
    IF (time < clim_reg%time(1)) THEN
      CALL crash('queried time out of range of direct transient climate forcing!')
    END IF

    ! Find time indices to be read
    IF (par%master) THEN

      IF     (time < clim_reg%time( 1)) THEN

        CALL warning('using constant start-of-record climate when extrapolating!')
        ti0 = 1
        ti1 = 1
        clim_reg%t0 = clim_reg%time( ti0) - 1._dp
        clim_reg%t1 = clim_reg%time( ti1)

      ELSEIF (time <= clim_reg%time( clim_reg%nyears)) THEN

        ti1 = 1
        DO WHILE (clim_reg%time( ti1) < time)
          ti1 = ti1 + 1
        END DO
        ti0 = ti1 - 1

        IF (ti0 == 0) THEN
          ti0 = 1
          ti1 = 2
        ELSEIF (ti1 == clim_reg%nyears) THEN
          ti0 = clim_reg%nyears - 1
          ti1 = ti0 + 1
        END IF

        clim_reg%t0 = clim_reg%time( ti0)
        clim_reg%t1 = clim_reg%time( ti1)

      ELSE ! IF     (time < clim_reg%time( 1)) THEN

        CALL warning('using constant end-of-record climate when extrapolating!')
        ti0 = clim_reg%nyears
        ti1 = clim_reg%nyears
        clim_reg%t0 = clim_reg%time( ti0) - 1._dp
        clim_reg%t1 = clim_reg%time( ti1)

      END IF ! IF     (time < clim_reg%time( 1)) THEN

    END IF ! IF (par%master) THEN

    ! Read new regional climate fields from the NetCDF file
    IF (par%master) CALL read_direct_regional_climate_file_timeframes( clim_reg, ti0, ti1)
    CALL sync

    ! Map the newly read data to the model grid
    CALL map_square_to_square_cons_2nd_order_2D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%Hs0_raw,      clim_reg%Hs0     )
    CALL map_square_to_square_cons_2nd_order_2D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%Hs1_raw,      clim_reg%Hs1     )
    CALL map_square_to_square_cons_2nd_order_3D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%T2m0_raw,     clim_reg%T2m0    )
    CALL map_square_to_square_cons_2nd_order_3D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%T2m1_raw,     clim_reg%T2m1    )
    CALL map_square_to_square_cons_2nd_order_3D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%Precip0_raw,  clim_reg%Precip0 )
    CALL map_square_to_square_cons_2nd_order_3D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%Precip1_raw,  clim_reg%Precip1 )
    CALL map_square_to_square_cons_2nd_order_3D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%Wind_WE0_raw, clim_reg%Wind_WE0)
    CALL map_square_to_square_cons_2nd_order_3D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%Wind_WE1_raw, clim_reg%Wind_WE1)
    CALL map_square_to_square_cons_2nd_order_3D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%Wind_SN0_raw, clim_reg%Wind_SN0)
    CALL map_square_to_square_cons_2nd_order_3D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%Wind_SN1_raw, clim_reg%Wind_SN1)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_direct_regional_climate_timeframes_from_file
  SUBROUTINE initialise_climate_model_direct_climate_regional( grid, climate_matrix, region_name)
    ! Initialise the regional climate model
    !
    ! Use a directly prescribed regional climate

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix
    CHARACTER(LEN=3),                         INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_model_direct_climate_regional'
    INTEGER                                            :: nx, ny

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_regional') THEN
      CALL crash('choice_climate_model should be "direct_regional"!')
    END IF

    ! The times at which we have climate fields from input, between which we'll interpolate
    ! to find the climate at model time (t0 <= model_time <= t1)

    CALL allocate_shared_dp_0D( climate_matrix%direct%t0, climate_matrix%direct%wt0)
    CALL allocate_shared_dp_0D( climate_matrix%direct%t1, climate_matrix%direct%wt1)

    IF (par%master) THEN
      ! Give impossible values to timeframes, so that the first call to run_climate_model_direct_climate_global
      ! is guaranteed to first read two new timeframes from the NetCDF file
      climate_matrix%direct%t0 = C%start_time_of_run - 100._dp
      climate_matrix%direct%t1 = C%start_time_of_run - 90._dp
    END IF ! IF (par%master) THEN
    CALL sync

    ! Inquire into the direct global cliamte forcing netcdf file
    CALL allocate_shared_int_0D( climate_matrix%direct%nyears, climate_matrix%direct%wnyears)
    CALL allocate_shared_int_0D( climate_matrix%direct%nx_raw, climate_matrix%direct%wnx_raw)
    CALL allocate_shared_int_0D( climate_matrix%direct%ny_raw, climate_matrix%direct%wny_raw)

    ! Determine name of file to read data from
    IF     (region_name == 'NAM') THEN
      climate_matrix%direct%netcdf%filename = C%filename_direct_regional_climate_NAM
    ELSEIF (region_name == 'EAS') THEN
      climate_matrix%direct%netcdf%filename = C%filename_direct_regional_climate_EAS
    ELSEIF (region_name == 'GRL') THEN
      climate_matrix%direct%netcdf%filename = C%filename_direct_regional_climate_GRL
    ELSEIF (region_name == 'ANT') THEN
      climate_matrix%direct%netcdf%filename = C%filename_direct_regional_climate_ANT
    END IF

    IF (par%master) WRITE(0,*) '  Initialising direct regional climate forcing from ', TRIM( climate_matrix%direct%netcdf%filename), '...'

    IF (par%master) CALL inquire_direct_regional_climate_forcing_file( climate_matrix%direct)
    CALL sync

    ! Abbreviations for shorter code
    nx = climate_matrix%direct%nx_raw
    ny = climate_matrix%direct%ny_raw

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( climate_matrix%direct%nyears, climate_matrix%direct%time, climate_matrix%direct%wtime        )
    CALL allocate_shared_dp_1D(                   nx, climate_matrix%direct%x_raw,        climate_matrix%direct%wx_raw       )
    CALL allocate_shared_dp_1D(          ny,          climate_matrix%direct%y_raw,        climate_matrix%direct%wy_raw       )

    CALL allocate_shared_dp_2D(          ny,      nx, climate_matrix%direct%Hs0_raw,      climate_matrix%direct%wHs0_raw     )
    CALL allocate_shared_dp_2D(          ny,      nx, climate_matrix%direct%Hs1_raw,      climate_matrix%direct%wHs1_raw     )
    CALL allocate_shared_dp_3D( 12,      ny,      nx, climate_matrix%direct%T2m0_raw,     climate_matrix%direct%wT2m0_raw    )
    CALL allocate_shared_dp_3D( 12,      ny,      nx, climate_matrix%direct%T2m1_raw,     climate_matrix%direct%wT2m1_raw    )
    CALL allocate_shared_dp_3D( 12,      ny,      nx, climate_matrix%direct%Precip0_raw,  climate_matrix%direct%wPrecip0_raw )
    CALL allocate_shared_dp_3D( 12,      ny,      nx, climate_matrix%direct%Precip1_raw,  climate_matrix%direct%wPrecip1_raw )
    CALL allocate_shared_dp_3D( 12,      ny,      nx, climate_matrix%direct%Wind_WE0_raw, climate_matrix%direct%wWind_WE0_raw)
    CALL allocate_shared_dp_3D( 12,      ny,      nx, climate_matrix%direct%Wind_WE1_raw, climate_matrix%direct%wWind_WE1_raw)
    CALL allocate_shared_dp_3D( 12,      ny,      nx, climate_matrix%direct%Wind_SN0_raw, climate_matrix%direct%wWind_SN0_raw)
    CALL allocate_shared_dp_3D( 12,      ny,      nx, climate_matrix%direct%Wind_SN1_raw, climate_matrix%direct%wWind_SN1_raw)

    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, climate_matrix%direct%Hs0,          climate_matrix%direct%wHs0         )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, climate_matrix%direct%Hs1,          climate_matrix%direct%wHs1         )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%T2m0,         climate_matrix%direct%wT2m0        )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%T2m1,         climate_matrix%direct%wT2m1        )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%Precip0,      climate_matrix%direct%wPrecip0     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%Precip1,      climate_matrix%direct%wPrecip1     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%Wind_WE0,     climate_matrix%direct%wWind_WE0    )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%Wind_WE1,     climate_matrix%direct%wWind_WE1    )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%Wind_SN0,     climate_matrix%direct%wWind_SN0    )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate_matrix%direct%Wind_SN1,     climate_matrix%direct%wWind_SN1    )

    ! Read time and grid data
    IF (par%master) CALL read_direct_regional_climate_file_time_xy( climate_matrix%direct)
    CALL sync

    ! Lastly, allocate memory for the "applied" snapshot
    CALL allocate_climate_snapshot_regional( grid, climate_matrix%applied, name = 'applied')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_direct_climate_regional

! == Directly prescribed global SMB
! =================================

  SUBROUTINE run_climate_model_direct_SMB_global( grid, clim_glob, climate_matrix, time)
    ! Run the regional climate model
    !
    ! Use a directly prescribed global SMB

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_direct_SMB_forcing_global),     INTENT(INOUT) :: clim_glob
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix
    REAL(dp),                                 INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_climate_model_direct_SMB_global'
    REAL(dp)                                           :: wt0, wt1
    INTEGER                                            :: i,j,m

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_global') THEN
      CALL crash('choice_SMB_model should be "direct_global"!')
    END IF

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time < climate_matrix%SMB_direct%t0 .OR. time > climate_matrix%SMB_direct%t1) THEN

      ! Find and read the two global time frames
      CALL update_direct_global_SMB_timeframes_from_file( clim_glob, time)

      ! Map the result to the regional grid
      CALL map_glob_to_grid_2D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%T2m_year0, climate_matrix%SMB_direct%T2m_year0)
      CALL map_glob_to_grid_2D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%SMB_year0, climate_matrix%SMB_direct%SMB_year0)

      CALL map_glob_to_grid_2D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%T2m_year1, climate_matrix%SMB_direct%T2m_year1)
      CALL map_glob_to_grid_2D( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%SMB_year1, climate_matrix%SMB_direct%SMB_year1)

      ! Update time stamps of regional climate forcing
      IF (par%master) THEN
        climate_matrix%SMB_direct%t0 = clim_glob%t0
        climate_matrix%SMB_direct%t1 = clim_glob%t1
      END IF
      CALL sync

    END IF ! IF (time >= climate_matrix%SMB_direct%t0 .AND. time <= climate_matrix%SMB_direct%t1) THEN

    ! Interpolate the two timeframes in time
    wt0 = (climate_matrix%SMB_direct%t1 - time) / (climate_matrix%SMB_direct%t1 - climate_matrix%SMB_direct%t0)
    wt1 = 1._dp - wt0

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      climate_matrix%applied%T2m( m,j,i) = (wt0 * climate_matrix%SMB_direct%T2m_year0( j,i)) + (wt1 * climate_matrix%SMB_direct%T2m_year1( j,i))
    END DO
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_direct_SMB_global
  SUBROUTINE update_direct_global_SMB_timeframes_from_file( clim_glob, time)
    ! Read the NetCDF file containing the global SMB forcing data. Only read the time
    ! frames enveloping the current coupling timestep to save on memory usage.

    IMPLICIT NONE

    TYPE(type_direct_SMB_forcing_global),     INTENT(INOUT) :: clim_glob
    REAL(dp),                                 INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_direct_global_SMB_timeframes_from_file'
    INTEGER                                            :: ti0, ti1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_global') THEN
      CALL crash('choice_SMB_model should be "direct_global"!')
    END IF

    ! Check if data for model time is available
    IF (time < clim_glob%time(1)) THEN
      CALL crash('queried time out of range of direct transient SMB forcing!')
    END IF

    ! Find time indices to be read
    IF (par%master) THEN

      IF     (time < clim_glob%time( 1)) THEN

        CALL warning('using constant start-of-record SMB when extrapolating!')
        ti0 = 1
        ti1 = 1
        clim_glob%t0 = clim_glob%time( ti0) - 1._dp
        clim_glob%t1 = clim_glob%time( ti1)

      ELSEIF (time <= clim_glob%time( clim_glob%nyears)) THEN

        ti1 = 1
        DO WHILE (clim_glob%time( ti1) < time)
          ti1 = ti1 + 1
        END DO
        ti0 = ti1 - 1

        IF (ti0 == 0) THEN
          ti0 = 1
          ti1 = 2
        ELSEIF (ti1 == clim_glob%nyears) THEN
          ti0 = clim_glob%nyears - 1
          ti1 = ti0 + 1
        END IF

        clim_glob%t0 = clim_glob%time( ti0)
        clim_glob%t1 = clim_glob%time( ti1)

      ELSE ! IF     (time < clim_glob%time( 1)) THEN

        CALL warning('using constant end-of-record SMB when extrapolating!')
        ti0 = clim_glob%nyears
        ti1 = clim_glob%nyears
        clim_glob%t0 = clim_glob%time( ti0) - 1._dp
        clim_glob%t1 = clim_glob%time( ti1)

      END IF ! IF     (time < clim_glob%time( 1)) THEN

    END IF ! IF (par%master) THEN

    ! Read new global climate fields from the NetCDF file
    IF (par%master) CALL read_direct_global_SMB_file_timeframes( clim_glob, ti0, ti1)
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_direct_global_SMB_timeframes_from_file
  SUBROUTINE initialise_climate_model_direct_SMB_global( clim_glob)
    ! Initialise the global climate model
    !
    ! Use a directly prescribed global SMB

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_direct_SMB_forcing_global), INTENT(INOUT) :: clim_glob

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_model_direct_SMB_global'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_global') THEN
      CALL crash('choice_SMB_model should be "direct_global"!')
    END IF

    ! The times at which we have climate fields from input, between which we'll interpolate
    ! to find the climate at model time (t0 <= model_time <= t1)

    CALL allocate_shared_dp_0D( clim_glob%t0, clim_glob%wt0)
    CALL allocate_shared_dp_0D( clim_glob%t1, clim_glob%wt1)

    IF (par%master) THEN
      ! Give impossible values to timeframes, so that the first call to run_climate_model_direct_climate_global
      ! is guaranteed to first read two new timeframes from the NetCDF file
      clim_glob%t0 = C%start_time_of_run - 100._dp
      clim_glob%t1 = C%start_time_of_run - 90._dp
    END IF ! IF (par%master) THEN
    CALL sync

    IF (par%master) WRITE(0,*) '  Initialising direct global SMB forcing from ', TRIM( C%filename_direct_global_SMB), '...'

    ! Inquire into the direct global cliamte forcing netcdf file
    CALL allocate_shared_int_0D( clim_glob%nyears, clim_glob%wnyears)
    CALL allocate_shared_int_0D( clim_glob%nlon,   clim_glob%wnlon  )
    CALL allocate_shared_int_0D( clim_glob%nlat,   clim_glob%wnlat  )

    clim_glob%netcdf%filename = C%filename_direct_global_SMB

    IF (par%master) CALL inquire_direct_global_SMB_forcing_file( clim_glob)
    CALL sync

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( clim_glob%nyears, clim_glob%time, clim_glob%wtime)
    CALL allocate_shared_dp_1D( clim_glob%nlon,   clim_glob%lon,  clim_glob%wlon )
    CALL allocate_shared_dp_1D( clim_glob%nlat,   clim_glob%lat,  clim_glob%wlat )

    CALL allocate_shared_dp_2D( clim_glob%nlon, clim_glob%nlat, clim_glob%T2m_year0, clim_glob%wT2m_year0)
    CALL allocate_shared_dp_2D( clim_glob%nlon, clim_glob%nlat, clim_glob%T2m_year1, clim_glob%wT2m_year1)
    CALL allocate_shared_dp_2D( clim_glob%nlon, clim_glob%nlat, clim_glob%SMB_year0, clim_glob%wSMB_year0)
    CALL allocate_shared_dp_2D( clim_glob%nlon, clim_glob%nlat, clim_glob%SMB_year1, clim_glob%wSMB_year1)

    ! Read time and grid data
    IF (par%master) CALL read_direct_global_SMB_file_time_latlon( clim_glob)
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_direct_SMB_global
  SUBROUTINE initialise_climate_model_direct_SMB_global_regional( grid, climate_matrix)
    ! Initialise the regional climate model
    !
    ! Use a directly prescribed global SMB

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_climate_model_direct_SMB_global_regional'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_global') THEN
      CALL crash('choice_SMB_model should be "direct_global"!')
    END IF

    ! The times at which we have climate fields from input, between which we'll interpolate
    ! to find the climate at model time (t0 <= model_time <= t1)

    CALL allocate_shared_dp_0D( climate_matrix%SMB_direct%t0, climate_matrix%SMB_direct%wt0)
    CALL allocate_shared_dp_0D( climate_matrix%SMB_direct%t1, climate_matrix%SMB_direct%wt1)

    IF (par%master) THEN
      ! Give impossible values to timeframes, so that the first call to run_climate_model_direct_climate_global
      ! is guaranteed to first read two new timeframes from the NetCDF file
      climate_matrix%SMB_direct%t0 = C%start_time_of_run - 100._dp
      climate_matrix%SMB_direct%t1 = C%start_time_of_run - 90._dp
    END IF ! IF (par%master) THEN
    CALL sync

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%SMB_direct%T2m_year0, climate_matrix%SMB_direct%wT2m_year0)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%SMB_direct%T2m_year1, climate_matrix%SMB_direct%wT2m_year1)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%SMB_direct%SMB_year0, climate_matrix%SMB_direct%wSMB_year0)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%SMB_direct%SMB_year1, climate_matrix%SMB_direct%wSMB_year1)

    ! Lastly, allocate memory for the "applied" snapshot
    CALL allocate_climate_snapshot_regional( grid, climate_matrix%applied, name = 'applied')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_direct_SMB_global_regional

! == Directly prescribed regional SMB
! ===================================

  SUBROUTINE run_climate_model_direct_SMB_regional( grid, climate_matrix, time)
    ! Run the regional climate model
    !
    ! Use a directly prescribed regional SMB

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix
    REAL(dp),                                 INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'run_climate_model_direct_SMB_regional'
    REAL(dp)                                                :: wt0, wt1
    INTEGER                                                 :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
      CALL crash('choice_SMB_model should be "direct_regional"!')
    END IF

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time < climate_matrix%SMB_direct%t0 .OR. time > climate_matrix%SMB_direct%t1) THEN

      ! Find and read the two global time frames
      CALL sync
      CALL update_direct_regional_SMB_timeframes_from_file( grid, climate_matrix%SMB_direct, time)

    END IF ! IF (time >= climate_matrix%SMB_direct%t0 .AND. time <= climate_matrix%SMB_direct%t1) THEN

    ! Interpolate the two timeframes in time
    wt0 = (climate_matrix%SMB_direct%t1 - time) / (climate_matrix%SMB_direct%t1 - climate_matrix%SMB_direct%t0)
    wt1 = 1._dp - wt0

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      climate_matrix%applied%T2m( :,j,i) = (wt0 * climate_matrix%SMB_direct%T2m_year0( j,i)) + (wt1 * climate_matrix%SMB_direct%T2m_year1( j,i))
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_direct_SMB_regional
  SUBROUTINE update_direct_regional_SMB_timeframes_from_file( grid, clim_reg, time)
    ! Read the NetCDF file containing the regional climate forcing data. Only read the time
    ! frames enveloping the current coupling timestep to save on memory usage.

    IMPLICIT NONE

    TYPE(type_grid),                            INTENT(IN)    :: grid
    TYPE(type_direct_SMB_forcing_regional),     INTENT(INOUT) :: clim_reg
    REAL(dp),                                   INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_direct_regional_SMB_timeframes_from_file'
    INTEGER                                            :: ti0, ti1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
      CALL crash('choice_SMB_model should be "direct_regional"!')
    END IF

    ! Find time indices to be read
    IF (par%master) THEN

      IF     (time < clim_reg%time( 1)) THEN

        CALL warning('using constant start-of-record SMB when extrapolating!')
        ti0 = 1
        ti1 = 1
        clim_reg%t0 = clim_reg%time( ti0) - 1._dp
        clim_reg%t1 = clim_reg%time( ti1)

      ELSEIF (time <= clim_reg%time( clim_reg%nyears)) THEN

        ti1 = 1
        DO WHILE (clim_reg%time( ti1) < time)
          ti1 = ti1 + 1
        END DO
        ti0 = ti1 - 1

        IF (ti0 == 0) THEN
          ti0 = 1
          ti1 = 2
        ELSEIF (ti1 == clim_reg%nyears) THEN
          ti0 = clim_reg%nyears - 1
          ti1 = ti0 + 1
        END IF

        clim_reg%t0 = clim_reg%time( ti0)
        clim_reg%t1 = clim_reg%time( ti1)

      ELSE ! IF     (time < clim_reg%time( 1)) THEN

        CALL warning('using constant end-of-record SMB when extrapolating!')
        ti0 = clim_reg%nyears
        ti1 = clim_reg%nyears
        clim_reg%t0 = clim_reg%time( ti0) - 1._dp
        clim_reg%t1 = clim_reg%time( ti1)

      END IF ! IF     (time < clim_reg%time( 1)) THEN

    END IF ! IF (par%master) THEN

    ! Read new regional climate fields from the NetCDF file
    IF (par%master) CALL read_direct_regional_SMB_file_timeframes( clim_reg, ti0, ti1)
    CALL sync

    ! Map the newly read data to the model grid
    CALL map_square_to_square_cons_2nd_order_2D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%T2m_year0_raw, clim_reg%T2m_year0)
    CALL map_square_to_square_cons_2nd_order_2D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%T2m_year1_raw, clim_reg%T2m_year1)
    CALL map_square_to_square_cons_2nd_order_2D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%SMB_year0_raw, clim_reg%SMB_year0)
    CALL map_square_to_square_cons_2nd_order_2D( clim_reg%nx_raw, clim_reg%ny_raw, clim_reg%x_raw, clim_reg%y_raw, grid%nx, grid%ny, grid%x, grid%y, clim_reg%SMB_year1_raw, clim_reg%SMB_year1)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_direct_regional_SMB_timeframes_from_file
  SUBROUTINE initialise_climate_model_direct_SMB_regional( grid, climate_matrix, region_name)
    ! Initialise the regional climate model
    !
    ! Use a directly prescribed regional SMB

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix
    CHARACTER(LEN=3),                         INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'initialise_climate_model_direct_SMB_regional'
    INTEGER                                                 :: nx, ny

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
      CALL crash('choice_SMB_model should be "direct_regional"!')
    END IF

    ! The times at which we have climate fields from input, between which we'll interpolate
    ! to find the climate at model time (t0 <= model_time <= t1)

    CALL allocate_shared_dp_0D( climate_matrix%SMB_direct%t0, climate_matrix%SMB_direct%wt0)
    CALL allocate_shared_dp_0D( climate_matrix%SMB_direct%t1, climate_matrix%SMB_direct%wt1)

    IF (par%master) THEN
      ! Give impossible values to timeframes, so that the first call to run_climate_model_direct_climate_global
      ! is guaranteed to first read two new timeframes from the NetCDF file
      climate_matrix%SMB_direct%t0 = C%start_time_of_run - 100._dp
      climate_matrix%SMB_direct%t1 = C%start_time_of_run - 90._dp
    END IF ! IF (par%master) THEN
    CALL sync

    ! Inquire into the direct global cliamte forcing netcdf file
    CALL allocate_shared_int_0D( climate_matrix%SMB_direct%nyears, climate_matrix%SMB_direct%wnyears)
    CALL allocate_shared_int_0D( climate_matrix%SMB_direct%nx_raw, climate_matrix%SMB_direct%wnx_raw)
    CALL allocate_shared_int_0D( climate_matrix%SMB_direct%ny_raw, climate_matrix%SMB_direct%wny_raw)

    ! Determine name of file to read data from
    IF     (region_name == 'NAM') THEN
      climate_matrix%SMB_direct%netcdf%filename = C%filename_direct_regional_SMB_NAM
    ELSEIF (region_name == 'EAS') THEN
      climate_matrix%SMB_direct%netcdf%filename = C%filename_direct_regional_SMB_EAS
    ELSEIF (region_name == 'GRL') THEN
      climate_matrix%SMB_direct%netcdf%filename = C%filename_direct_regional_SMB_GRL
    ELSEIF (region_name == 'ANT') THEN
      climate_matrix%SMB_direct%netcdf%filename = C%filename_direct_regional_SMB_ANT
    END IF

    IF (par%master) WRITE(0,*) '  Initialising direct regional SMB forcing from ', TRIM( climate_matrix%SMB_direct%netcdf%filename), '...'

    IF (par%master) CALL inquire_direct_regional_SMB_forcing_file( climate_matrix%SMB_direct)
    CALL sync

    ! Abbreviations for shorter code
    nx = climate_matrix%SMB_direct%nx_raw
    ny = climate_matrix%SMB_direct%ny_raw

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( climate_matrix%SMB_direct%nyears, climate_matrix%SMB_direct%time, climate_matrix%SMB_direct%wtime  )
    CALL allocate_shared_dp_1D(               nx, climate_matrix%SMB_direct%x_raw,         climate_matrix%SMB_direct%wx_raw        )
    CALL allocate_shared_dp_1D(      ny,          climate_matrix%SMB_direct%y_raw,         climate_matrix%SMB_direct%wy_raw        )

    CALL allocate_shared_dp_2D(      ny,      nx, climate_matrix%SMB_direct%T2m_year0_raw, climate_matrix%SMB_direct%wT2m_year0_raw)
    CALL allocate_shared_dp_2D(      ny,      nx, climate_matrix%SMB_direct%T2m_year1_raw, climate_matrix%SMB_direct%wT2m_year1_raw)
    CALL allocate_shared_dp_2D(      ny,      nx, climate_matrix%SMB_direct%SMB_year0_raw, climate_matrix%SMB_direct%wSMB_year0_raw)
    CALL allocate_shared_dp_2D(      ny,      nx, climate_matrix%SMB_direct%SMB_year1_raw, climate_matrix%SMB_direct%wSMB_year1_raw)

    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%SMB_direct%T2m_year0,     climate_matrix%SMB_direct%wT2m_year0    )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%SMB_direct%T2m_year1,     climate_matrix%SMB_direct%wT2m_year1    )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%SMB_direct%SMB_year0,     climate_matrix%SMB_direct%wSMB_year0    )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%SMB_direct%SMB_year1,     climate_matrix%SMB_direct%wSMB_year1    )

    ! Read time and grid data
    IF (par%master) CALL read_direct_regional_SMB_file_time_xy( climate_matrix%SMB_direct)
    CALL sync

    ! Lastly, allocate memory for the "applied" snapshot
    CALL allocate_climate_snapshot_regional( grid, climate_matrix%applied, name = 'applied')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_direct_SMB_regional

! == ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing
! =================================================================

  SUBROUTINE run_climate_model_ISMIP_forcing( grid, climate_matrix, time, refgeo_PD, ice)
    ! Run the regional climate model
    !
    ! Use the ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix
    REAL(dp),                                 INTENT(IN)    :: time
    TYPE(type_reference_geometry),            INTENT(IN)    :: refgeo_PD
    TYPE(type_ice_model),                     INTENT(IN)    :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'run_climate_model_ISMIP_forcing'
    REAL(dp)                                                :: wt0, wt1
    INTEGER                                                 :: i,j
    REAL(dp)                                                :: dz

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time < climate_matrix%ISMIP_forcing%t0 .OR. time > climate_matrix%ISMIP_forcing%t1) THEN

      ! Find and read the two global time frames
      CALL sync
      CALL update_ISMIP_forcing_timeframes( grid, climate_matrix, time)

    END IF ! IF (time >= climate_matrix%SMB_direct%t0 .AND. time <= climate_matrix%SMB_direct%t1) THEN

    ! Interpolate the two timeframes in time
    wt0 = (climate_matrix%ISMIP_forcing%t1 - time) / (climate_matrix%ISMIP_forcing%t1 - climate_matrix%ISMIP_forcing%t0)
    wt1 = 1._dp - wt0

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      climate_matrix%ISMIP_forcing%aSMB(   j,i) = (wt0 * climate_matrix%ISMIP_forcing%aSMB0(   j,i)) + &
                                                  (wt1 * climate_matrix%ISMIP_forcing%aSMB1(   j,i))
      climate_matrix%ISMIP_forcing%dSMBdz( j,i) = (wt0 * climate_matrix%ISMIP_forcing%dSMBdz0( j,i)) + &
                                                  (wt1 * climate_matrix%ISMIP_forcing%dSMBdz1( j,i))
      climate_matrix%ISMIP_forcing%aST(    j,i) = (wt0 * climate_matrix%ISMIP_forcing%aST0(    j,i)) + &
                                                  (wt1 * climate_matrix%ISMIP_forcing%aST1(    j,i))
      climate_matrix%ISMIP_forcing%dSTdz(  j,i) = (wt0 * climate_matrix%ISMIP_forcing%dSTdz0(  j,i)) + &
                                                  (wt1 * climate_matrix%ISMIP_forcing%dSTdz1(  j,i))

    END DO
    END DO
    CALL sync

    ! Apply the anomaly and elevation correction to calculate the applied SMB and temperature
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      dz = ice%Hs_a( j,i) - refgeo_PD%Hs( j,i)

      climate_matrix%ISMIP_forcing%SMB( j,i) = climate_matrix%ISMIP_forcing%SMB_baseline( j,i) + &
                                               climate_matrix%ISMIP_forcing%aSMB(         j,i) + &
                                               climate_matrix%ISMIP_forcing%dSMBdz(       j,i) * dz

      climate_matrix%ISMIP_forcing%ST(  j,i) = climate_matrix%ISMIP_forcing%ST_baseline(  j,i) + &
                                               climate_matrix%ISMIP_forcing%aST(          j,i) + &
                                               climate_matrix%ISMIP_forcing%dSTdz(        j,i) * dz

    END DO
    END DO
    CALL sync

    ! Set the final values in the "applied" climate snapshot
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      climate_matrix%applied%T2m( :,j,i) = climate_matrix%ISMIP_forcing%ST( j,i)

    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_ISMIP_forcing
  SUBROUTINE update_ISMIP_forcing_timeframes( grid, climate_matrix, time)
    ! Update the two timeframes of the ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing data
    !
    ! This is where we deal with the fact that ISMIP supplies all these fields in separate NetCDF files
    ! for each year.

    USE parameters_module, ONLY: sec_per_year, ice_density

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix
    REAL(dp),                                 INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'update_ISMIP_forcing_timeframes'
    INTEGER                                                 :: i,j,i1,i2
    INTEGER                                                 :: year0, year1
    CHARACTER(LEN=4)                                        :: year0str, year1str
    TYPE(type_grid)                                         :: grid_raw
    REAL(dp), DIMENSION(:,:  ), POINTER                     ::  aSMB_raw,  dSMBdz_raw,  aST_raw,  dSTdz_raw
    INTEGER                                                 :: waSMB_raw, wdSMBdz_raw, waST_raw, wdSTdz_raw
    REAL(dp), DIMENSION(:,:,:), POINTER                     ::  data_combi_raw,  data_combi_mapped
    INTEGER                                                 :: wdata_combi_raw, wdata_combi_mapped

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Since ISMIP supplies NetCDF files for each individual year, determining
    ! which two files we need to read is simple enough.

    year0 = FLOOR(        time)
    year1 = MAX( CEILING( time), year0+1)

    WRITE( year0str,'(I4)') year0
    WRITE( year1str,'(I4)') year1

    ! Update timestamps
    climate_matrix%ISMIP_forcing%t0 = REAL( year0,dp)
    climate_matrix%ISMIP_forcing%t1 = REAL( year1,dp)

  ! ===== aSMB =====
  ! ================

    ! Timeframe 0
    ! ===========

    ! Determine the name of the file containing the timeframe
    climate_matrix%ISMIP_forcing%netcdf_aSMB0%filename = TRIM( C%ISMIP_forcing_foldername_aSMB) // '/' // TRIM( C%ISMIP_forcing_basefilename_aSMB) // year0str // '.nc'

    ! Inquire if everything we need is present in the file
    CALL allocate_shared_int_0D( grid_raw%nx, grid_raw%wnx)
    CALL allocate_shared_int_0D( grid_raw%ny, grid_raw%wny)
    CALL inquire_ISMIP_forcing_aSMB_file( climate_matrix%ISMIP_forcing%netcdf_aSMB0, grid_raw%nx, grid_raw%ny)
    CALL sync

    ! Allocate memory to store all raw data in a single array (to make remapping faster)
    CALL allocate_shared_dp_3D( grid_raw%nx, grid_raw%ny, 8, data_combi_raw, wdata_combi_raw)

    ! Allocate memory, read the data
    CALL allocate_shared_dp_0D(                           grid_raw%dx, grid_raw%wdx)
    CALL allocate_shared_dp_1D( grid_raw%nx,              grid_raw%x , grid_raw%wx )
    CALL allocate_shared_dp_1D(              grid_raw%ny, grid_raw%y , grid_raw%wy )
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, aSMB_raw   , waSMB_raw   )
    CALL read_ISMIP_forcing_aSMB_file( climate_matrix%ISMIP_forcing%netcdf_aSMB0, grid_raw%x, grid_raw%y, aSMB_raw)
    IF (par%master) grid_raw%dx = grid_raw%x( 2) - grid_raw%x( 1)
    CALL sync

    ! Clean up the raw data a bit
    CALL partition_list( grid_raw%nx, par%i, par%n, i1, i2)
    DO i = i1, i2
    DO j = 1, grid_raw%ny

      ! Convert the data from SI units (kg m-2 s-1) to ice model units (m.i.e./yr)
      aSMB_raw( i,j) = aSMB_raw( i,j) * sec_per_year / ice_density

      ! Replace NaN by zero
      IF (ISNAN( aSMB_raw( i,j))) THEN
        aSMB_raw( i,j) = 0._dp
      END IF

    END DO
    END DO
    CALL sync

    ! Store in the big 3-D array
    data_combi_raw( i1:i2,:,1) = aSMB_raw( i1:i2,:)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( grid_raw%wnx)
    CALL deallocate_shared( grid_raw%wny)
    CALL deallocate_shared( grid_raw%wdx)
    CALL deallocate_shared( grid_raw%wx )
    CALL deallocate_shared( grid_raw%wy )
    CALL deallocate_shared( waSMB_raw   )

    ! Timeframe 1
    ! ===========

    ! Determine the name of the file containing the timeframe
    climate_matrix%ISMIP_forcing%netcdf_aSMB1%filename = TRIM( C%ISMIP_forcing_foldername_aSMB) // '/' // TRIM( C%ISMIP_forcing_basefilename_aSMB) // year1str // '.nc'

    ! Inquire if everything we need is present in the file
    CALL allocate_shared_int_0D( grid_raw%nx, grid_raw%wnx)
    CALL allocate_shared_int_0D( grid_raw%ny, grid_raw%wny)
    CALL inquire_ISMIP_forcing_aSMB_file( climate_matrix%ISMIP_forcing%netcdf_aSMB1, grid_raw%nx, grid_raw%ny)
    CALL sync

    ! Allocate memory, read the data
    CALL allocate_shared_dp_0D(                           grid_raw%dx, grid_raw%wdx)
    CALL allocate_shared_dp_1D( grid_raw%nx,              grid_raw%x , grid_raw%wx )
    CALL allocate_shared_dp_1D(              grid_raw%ny, grid_raw%y , grid_raw%wy )
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, aSMB_raw   , waSMB_raw   )
    CALL read_ISMIP_forcing_aSMB_file( climate_matrix%ISMIP_forcing%netcdf_aSMB1, grid_raw%x, grid_raw%y, aSMB_raw)
    IF (par%master) grid_raw%dx = grid_raw%x( 2) - grid_raw%x( 1)
    CALL sync

    ! Clean up the raw data a bit
    CALL partition_list( grid_raw%nx, par%i, par%n, i1, i2)
    DO i = i1, i2
    DO j = 1, grid_raw%ny

      ! Convert the data from SI units (kg m-2 s-1) to ice model units (m.i.e./yr)
      aSMB_raw( i,j) = aSMB_raw( i,j) * sec_per_year / ice_density

      ! Replace NaN by zero
      IF (ISNAN( aSMB_raw( i,j))) THEN
        aSMB_raw( i,j) = 0._dp
      END IF

    END DO
    END DO
    CALL sync

    ! Store in the big 3-D array
    data_combi_raw( i1:i2,:,2) = aSMB_raw( i1:i2,:)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( grid_raw%wnx)
    CALL deallocate_shared( grid_raw%wny)
    CALL deallocate_shared( grid_raw%wdx)
    CALL deallocate_shared( grid_raw%wx )
    CALL deallocate_shared( grid_raw%wy )
    CALL deallocate_shared( waSMB_raw   )

  ! ===== dSMBdz =====
  ! ================

    ! Timeframe 0
    ! ===========

    ! Determine the name of the file containing the timeframe
    climate_matrix%ISMIP_forcing%netcdf_dSMBdz0%filename = TRIM( C%ISMIP_forcing_foldername_dSMBdz) // '/' // TRIM( C%ISMIP_forcing_basefilename_dSMBdz) // year0str // '.nc'

    ! Inquire if everything we need is present in the file
    CALL allocate_shared_int_0D( grid_raw%nx, grid_raw%wnx)
    CALL allocate_shared_int_0D( grid_raw%ny, grid_raw%wny)
    CALL inquire_ISMIP_forcing_dSMBdz_file( climate_matrix%ISMIP_forcing%netcdf_dSMBdz0, grid_raw%nx, grid_raw%ny)
    CALL sync

    ! Allocate memory, read the data
    CALL allocate_shared_dp_0D(                           grid_raw%dx, grid_raw%wdx)
    CALL allocate_shared_dp_1D( grid_raw%nx,              grid_raw%x , grid_raw%wx )
    CALL allocate_shared_dp_1D(              grid_raw%ny, grid_raw%y , grid_raw%wy )
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, dSMBdz_raw , wdSMBdz_raw )
    CALL read_ISMIP_forcing_dSMBdz_file( climate_matrix%ISMIP_forcing%netcdf_dSMBdz0, grid_raw%x, grid_raw%y, dSMBdz_raw)
    IF (par%master) grid_raw%dx = grid_raw%x( 2) - grid_raw%x( 1)
    CALL sync

    ! Clean up the raw data a bit
    CALL partition_list( grid_raw%nx, par%i, par%n, i1, i2)
    DO i = i1, i2
    DO j = 1, grid_raw%ny

      ! Convert the data from SI units (kg m-2 s-1) to ice model units (m.i.e./yr)
      dSMBdz_raw( i,j) = dSMBdz_raw( i,j) * sec_per_year / ice_density

      ! Replace NaN by zero
      IF (ISNAN( dSMBdz_raw( i,j))) THEN
        dSMBdz_raw( i,j) = 0._dp
      END IF

    END DO
    END DO
    CALL sync

    ! Store in the big 3-D array
    data_combi_raw( i1:i2,:,3) = dSMBdz_raw( i1:i2,:)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( grid_raw%wnx)
    CALL deallocate_shared( grid_raw%wny)
    CALL deallocate_shared( grid_raw%wdx)
    CALL deallocate_shared( grid_raw%wx )
    CALL deallocate_shared( grid_raw%wy )
    CALL deallocate_shared( wdSMBdz_raw )

    ! Timeframe 1
    ! ===========

    ! Determine the name of the file containing the timeframe
    climate_matrix%ISMIP_forcing%netcdf_dSMBdz1%filename = TRIM( C%ISMIP_forcing_foldername_dSMBdz) // '/' // TRIM( C%ISMIP_forcing_basefilename_dSMBdz) // year1str // '.nc'

    ! Inquire if everything we need is present in the file
    CALL allocate_shared_int_0D( grid_raw%nx, grid_raw%wnx)
    CALL allocate_shared_int_0D( grid_raw%ny, grid_raw%wny)
    CALL inquire_ISMIP_forcing_dSMBdz_file( climate_matrix%ISMIP_forcing%netcdf_dSMBdz1, grid_raw%nx, grid_raw%ny)
    CALL sync

    ! Allocate memory, read the data
    CALL allocate_shared_dp_0D(                           grid_raw%dx, grid_raw%wdx)
    CALL allocate_shared_dp_1D( grid_raw%nx,              grid_raw%x , grid_raw%wx )
    CALL allocate_shared_dp_1D(              grid_raw%ny, grid_raw%y , grid_raw%wy )
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, dSMBdz_raw , wdSMBdz_raw )
    CALL read_ISMIP_forcing_dSMBdz_file( climate_matrix%ISMIP_forcing%netcdf_dSMBdz1, grid_raw%x, grid_raw%y, dSMBdz_raw)
    IF (par%master) grid_raw%dx = grid_raw%x( 2) - grid_raw%x( 1)
    CALL sync

    ! Clean up the raw data a bit
    CALL partition_list( grid_raw%nx, par%i, par%n, i1, i2)
    DO i = i1, i2
    DO j = 1, grid_raw%ny

      ! Convert the data from SI units (kg m-2 s-1) to ice model units (m.i.e./yr)
      dSMBdz_raw( i,j) = dSMBdz_raw( i,j) * sec_per_year / ice_density

      ! Replace NaN by zero
      IF (ISNAN( dSMBdz_raw( i,j))) THEN
        dSMBdz_raw( i,j) = 0._dp
      END IF

    END DO
    END DO
    CALL sync

    ! Store in the big 3-D array
    data_combi_raw( i1:i2,:,4) = dSMBdz_raw( i1:i2,:)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( grid_raw%wnx)
    CALL deallocate_shared( grid_raw%wny)
    CALL deallocate_shared( grid_raw%wdx)
    CALL deallocate_shared( grid_raw%wx )
    CALL deallocate_shared( grid_raw%wy )
    CALL deallocate_shared( wdSMBdz_raw )

  ! ===== aST =====
  ! ================

    ! Timeframe 0
    ! ===========

    ! Determine the name of the file containing the timeframe
    climate_matrix%ISMIP_forcing%netcdf_aST0%filename = TRIM( C%ISMIP_forcing_foldername_aST) // '/' // TRIM( C%ISMIP_forcing_basefilename_aST) // year0str // '.nc'

    ! Inquire if everything we need is present in the file
    CALL allocate_shared_int_0D( grid_raw%nx, grid_raw%wnx)
    CALL allocate_shared_int_0D( grid_raw%ny, grid_raw%wny)
    CALL inquire_ISMIP_forcing_aST_file( climate_matrix%ISMIP_forcing%netcdf_aST0, grid_raw%nx, grid_raw%ny)
    CALL sync

    ! Allocate memory, read the data
    CALL allocate_shared_dp_0D(                           grid_raw%dx, grid_raw%wdx)
    CALL allocate_shared_dp_1D( grid_raw%nx,              grid_raw%x , grid_raw%wx )
    CALL allocate_shared_dp_1D(              grid_raw%ny, grid_raw%y , grid_raw%wy )
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, aST_raw    , waST_raw    )
    CALL read_ISMIP_forcing_aST_file( climate_matrix%ISMIP_forcing%netcdf_aST0, grid_raw%x, grid_raw%y, aST_raw)
    IF (par%master) grid_raw%dx = grid_raw%x( 2) - grid_raw%x( 1)
    CALL sync

    ! Clean up the raw data a bit
    CALL partition_list( grid_raw%nx, par%i, par%n, i1, i2)
    DO i = i1, i2
    DO j = 1, grid_raw%ny

      ! Replace NaN by zero
      IF (ISNAN( aST_raw( i,j))) THEN
        aST_raw( i,j) = 0._dp
      END IF

    END DO
    END DO
    CALL sync

    ! Store in the big 3-D array
    data_combi_raw( i1:i2,:,5) = aST_raw( i1:i2,:)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( grid_raw%wnx)
    CALL deallocate_shared( grid_raw%wny)
    CALL deallocate_shared( grid_raw%wdx)
    CALL deallocate_shared( grid_raw%wx )
    CALL deallocate_shared( grid_raw%wy )
    CALL deallocate_shared( waST_raw    )

    ! Timeframe 1
    ! ===========

    ! Determine the name of the file containing the timeframe
    climate_matrix%ISMIP_forcing%netcdf_aST1%filename = TRIM( C%ISMIP_forcing_foldername_aST) // '/' // TRIM( C%ISMIP_forcing_basefilename_aST) // year1str // '.nc'

    ! Inquire if everything we need is present in the file
    CALL allocate_shared_int_0D( grid_raw%nx, grid_raw%wnx)
    CALL allocate_shared_int_0D( grid_raw%ny, grid_raw%wny)
    CALL inquire_ISMIP_forcing_aST_file( climate_matrix%ISMIP_forcing%netcdf_aST1, grid_raw%nx, grid_raw%ny)
    CALL sync

    ! Allocate memory, read the data
    CALL allocate_shared_dp_0D(                           grid_raw%dx, grid_raw%wdx)
    CALL allocate_shared_dp_1D( grid_raw%nx,              grid_raw%x , grid_raw%wx )
    CALL allocate_shared_dp_1D(              grid_raw%ny, grid_raw%y , grid_raw%wy )
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, aST_raw    , waST_raw    )
    CALL read_ISMIP_forcing_aST_file( climate_matrix%ISMIP_forcing%netcdf_aST1, grid_raw%x, grid_raw%y, aST_raw)
    IF (par%master) grid_raw%dx = grid_raw%x( 2) - grid_raw%x( 1)
    CALL sync

    ! Clean up the raw data a bit
    CALL partition_list( grid_raw%nx, par%i, par%n, i1, i2)
    DO i = i1, i2
    DO j = 1, grid_raw%ny

      ! Replace NaN by zero
      IF (ISNAN( aST_raw( i,j))) THEN
        aST_raw( i,j) = 0._dp
      END IF

    END DO
    END DO
    CALL sync

    ! Store in the big 3-D array
    data_combi_raw( i1:i2,:,6) = aST_raw( i1:i2,:)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( grid_raw%wnx)
    CALL deallocate_shared( grid_raw%wny)
    CALL deallocate_shared( grid_raw%wdx)
    CALL deallocate_shared( grid_raw%wx )
    CALL deallocate_shared( grid_raw%wy )
    CALL deallocate_shared( waST_raw    )

  ! ===== dSTdz =====
  ! ================

    ! Timeframe 0
    ! ===========

    ! Determine the name of the file containing the timeframe
    climate_matrix%ISMIP_forcing%netcdf_dSTdz0%filename = TRIM( C%ISMIP_forcing_foldername_dSTdz) // '/' // TRIM( C%ISMIP_forcing_basefilename_dSTdz) // year0str // '.nc'

    ! Inquire if everything we need is present in the file
    CALL allocate_shared_int_0D( grid_raw%nx, grid_raw%wnx)
    CALL allocate_shared_int_0D( grid_raw%ny, grid_raw%wny)
    CALL inquire_ISMIP_forcing_dSTdz_file( climate_matrix%ISMIP_forcing%netcdf_dSTdz0, grid_raw%nx, grid_raw%ny)
    CALL sync

    ! Allocate memory, read the data
    CALL allocate_shared_dp_0D(                           grid_raw%dx, grid_raw%wdx)
    CALL allocate_shared_dp_1D( grid_raw%nx,              grid_raw%x , grid_raw%wx )
    CALL allocate_shared_dp_1D(              grid_raw%ny, grid_raw%y , grid_raw%wy )
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, dSTdz_raw  , wdSTdz_raw  )
    CALL read_ISMIP_forcing_dSTdz_file( climate_matrix%ISMIP_forcing%netcdf_dSTdz0, grid_raw%x, grid_raw%y, dSTdz_raw)
    IF (par%master) grid_raw%dx = grid_raw%x( 2) - grid_raw%x( 1)
    CALL sync

    ! Clean up the raw data a bit
    CALL partition_list( grid_raw%nx, par%i, par%n, i1, i2)
    DO i = i1, i2
    DO j = 1, grid_raw%ny

      ! Replace NaN by zero
      IF (ISNAN( dSTdz_raw( i,j))) THEN
        dSTdz_raw( i,j) = 0._dp
      END IF

    END DO
    END DO
    CALL sync

    ! Store in the big 3-D array
    data_combi_raw( i1:i2,:,7) = dSTdz_raw( i1:i2,:)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( grid_raw%wnx)
    CALL deallocate_shared( grid_raw%wny)
    CALL deallocate_shared( grid_raw%wdx)
    CALL deallocate_shared( grid_raw%wx )
    CALL deallocate_shared( grid_raw%wy )
    CALL deallocate_shared( wdSTdz_raw  )

    ! Timeframe 1
    ! ===========

    ! Determine the name of the file containing the timeframe
    climate_matrix%ISMIP_forcing%netcdf_dSTdz1%filename = TRIM( C%ISMIP_forcing_foldername_dSTdz) // '/' // TRIM( C%ISMIP_forcing_basefilename_dSTdz) // year1str // '.nc'

    ! Inquire if everything we need is present in the file
    CALL allocate_shared_int_0D( grid_raw%nx, grid_raw%wnx)
    CALL allocate_shared_int_0D( grid_raw%ny, grid_raw%wny)
    CALL inquire_ISMIP_forcing_dSTdz_file( climate_matrix%ISMIP_forcing%netcdf_dSTdz1, grid_raw%nx, grid_raw%ny)
    CALL sync

    ! Allocate memory, read the data
    CALL allocate_shared_dp_0D(                           grid_raw%dx, grid_raw%wdx)
    CALL allocate_shared_dp_1D( grid_raw%nx,              grid_raw%x , grid_raw%wx )
    CALL allocate_shared_dp_1D(              grid_raw%ny, grid_raw%y , grid_raw%wy )
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, dSTdz_raw  , wdSTdz_raw  )
    CALL read_ISMIP_forcing_dSTdz_file( climate_matrix%ISMIP_forcing%netcdf_dSTdz1, grid_raw%x, grid_raw%y, dSTdz_raw)
    IF (par%master) grid_raw%dx = grid_raw%x( 2) - grid_raw%x( 1)
    CALL sync

    ! Clean up the raw data a bit
    CALL partition_list( grid_raw%nx, par%i, par%n, i1, i2)
    DO i = i1, i2
    DO j = 1, grid_raw%ny

      ! Replace NaN by zero
      IF (ISNAN( dSTdz_raw( i,j))) THEN
        dSTdz_raw( i,j) = 0._dp
      END IF

    END DO
    END DO
    CALL sync

    ! Store in the big 3-D array
    data_combi_raw( i1:i2,:,8) = dSTdz_raw( i1:i2,:)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wdSTdz_raw  )

    ! Transpose the data from [i,j,k] to [k,j,i] indexing
    CALL transpose_dp_3D( data_combi_raw, wdata_combi_raw)

    ! Map the data to the model grid
    CALL allocate_shared_dp_3D( 8, grid%ny, grid%nx, data_combi_mapped, wdata_combi_mapped)
    CALL map_square_to_square_cons_2nd_order_3D( grid_raw%nx, grid_raw%ny, grid_raw%x, grid_raw%y, &
      grid%nx, grid%ny, grid%x, grid%y, data_combi_raw, data_combi_mapped)

    ! Get mapped data back from the 3-D combined array
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      climate_matrix%ISMIP_forcing%aSMB0(   j,i) = data_combi_mapped( 1,j,i)
      climate_matrix%ISMIP_forcing%aSMB1(   j,i) = data_combi_mapped( 2,j,i)
      climate_matrix%ISMIP_forcing%dSMBdz0( j,i) = data_combi_mapped( 3,j,i)
      climate_matrix%ISMIP_forcing%dSMBdz1( j,i) = data_combi_mapped( 4,j,i)
      climate_matrix%ISMIP_forcing%aST0(    j,i) = data_combi_mapped( 5,j,i)
      climate_matrix%ISMIP_forcing%aST1(    j,i) = data_combi_mapped( 6,j,i)
      climate_matrix%ISMIP_forcing%dSTdz0(  j,i) = data_combi_mapped( 7,j,i)
      climate_matrix%ISMIP_forcing%dSTdz1(  j,i) = data_combi_mapped( 8,j,i)
    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( grid_raw%wnx      )
    CALL deallocate_shared( grid_raw%wny      )
    CALL deallocate_shared( grid_raw%wdx      )
    CALL deallocate_shared( grid_raw%wx       )
    CALL deallocate_shared( grid_raw%wy       )
    CALL deallocate_shared( wdata_combi_raw   )
    CALL deallocate_shared( wdata_combi_mapped)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_ISMIP_forcing_timeframes
  SUBROUTINE initialise_climate_model_ISMIP_forcing( grid, climate_matrix)
    ! Initialise the regional climate model
    !
    ! Use the ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'initialise_climate_model_ISMIP_forcing'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise the baseline SMB and temperature
    CALL initialise_climate_model_ISMIP_forcing_SMB_baseline( grid, climate_matrix)
    CALL initialise_climate_model_ISMIP_forcing_ST_baseline(  grid, climate_matrix)

    ! Allocate memory for the timestamps of the two timeframes
    CALL allocate_shared_dp_0D( climate_matrix%ISMIP_forcing%t0, climate_matrix%ISMIP_forcing%wt0)
    CALL allocate_shared_dp_0D( climate_matrix%ISMIP_forcing%t1, climate_matrix%ISMIP_forcing%wt1)

    IF (par%master) THEN
      ! Give impossible values to timeframes, so that the first call to run_climate_model_ISMIP_forcing
      ! is guaranteed to first read two new timeframes from the NetCDF file
      climate_matrix%ISMIP_forcing%t0 = C%start_time_of_run - 100._dp
      climate_matrix%ISMIP_forcing%t1 = C%start_time_of_run - 90._dp
    END IF ! IF (par%master) THEN
    CALL sync

    ! Allocate memory for the two timeframes
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%aSMB0  , climate_matrix%ISMIP_forcing%waSMB0  )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%dSMBdz0, climate_matrix%ISMIP_forcing%wdSMBdz0)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%aST0   , climate_matrix%ISMIP_forcing%waST0   )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%dSTdz0 , climate_matrix%ISMIP_forcing%wdSTdz0 )

    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%aSMB1  , climate_matrix%ISMIP_forcing%waSMB1  )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%dSMBdz1, climate_matrix%ISMIP_forcing%wdSMBdz1)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%aST1   , climate_matrix%ISMIP_forcing%waST1   )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%dSTdz1 , climate_matrix%ISMIP_forcing%wdSTdz1 )

    ! Allocate memory for the time-interpolated values of aSMB, dSMBdz, ST, and dSTdz
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%aSMB   , climate_matrix%ISMIP_forcing%waSMB   )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%dSMBdz , climate_matrix%ISMIP_forcing%wdSMBdz )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%aST    , climate_matrix%ISMIP_forcing%waST    )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%dSTdz  , climate_matrix%ISMIP_forcing%wdSTdz  )

    ! Allocate memory for the applied values of SMB and ST (i.e. after applying the anomaly and elevation correction)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%SMB    , climate_matrix%ISMIP_forcing%wSMB    )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%ST     , climate_matrix%ISMIP_forcing%wST     )

    ! Lastly, allocate memory for the "applied" snapshot
    CALL allocate_climate_snapshot_regional( grid, climate_matrix%applied, name = 'applied')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_ISMIP_forcing
  SUBROUTINE initialise_climate_model_ISMIP_forcing_SMB_baseline( grid, climate_matrix)
    ! Use the ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing
    ! Initialise the baseline SMB

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'initialise_climate_model_ISMIP_forcing_SMB_baseline'
    TYPE(type_grid)                                         :: grid_raw
    REAL(dp), DIMENSION(:,:  ), POINTER                     ::  SMB_raw
    INTEGER                                                 :: wSMB_raw

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire if everything we need is present in the file
    climate_matrix%ISMIP_forcing%netcdf_SMB_baseline%filename = C%ISMIP_forcing_filename_SMB_baseline
    CALL allocate_shared_int_0D( grid_raw%nx, grid_raw%wnx)
    CALL allocate_shared_int_0D( grid_raw%ny, grid_raw%wny)
    CALL inquire_ISMIP_forcing_SMB_baseline_file( climate_matrix%ISMIP_forcing%netcdf_SMB_baseline, grid_raw%nx, grid_raw%ny)
    CALL sync

    ! Allocate memory, read the data
    CALL allocate_shared_dp_0D(                           grid_raw%dx, grid_raw%wdx)
    CALL allocate_shared_dp_1D( grid_raw%nx,              grid_raw%x , grid_raw%wx )
    CALL allocate_shared_dp_1D(              grid_raw%ny, grid_raw%y , grid_raw%wy )
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, SMB_raw    , wSMB_raw    )
    CALL read_ISMIP_forcing_SMB_baseline_file( climate_matrix%ISMIP_forcing%netcdf_SMB_baseline, grid_raw%x, grid_raw%y, SMB_raw)
    IF (par%master) grid_raw%dx = grid_raw%x( 2) - grid_raw%x( 1)
    CALL sync

    ! Transpose the data to [j,i] indexing
    CALL transpose_dp_2D( SMB_raw, wSMB_raw)

    ! Map the data to the model grid
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%SMB_baseline, climate_matrix%ISMIP_forcing%wSMB_baseline)
    CALL map_square_to_square_cons_2nd_order_2D( grid_raw%nx, grid_raw%ny, grid_raw%x, grid_raw%y, &
      grid%nx, grid%ny, grid%x, grid%y, SMB_raw, climate_matrix%ISMIP_forcing%SMB_baseline)

    ! Clean up after yourself
    CALL deallocate_shared( grid_raw%wnx)
    CALL deallocate_shared( grid_raw%wny)
    CALL deallocate_shared( grid_raw%wdx)
    CALL deallocate_shared( grid_raw%wx )
    CALL deallocate_shared( grid_raw%wy )
    CALL deallocate_shared( wSMB_raw    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_ISMIP_forcing_SMB_baseline
  SUBROUTINE initialise_climate_model_ISMIP_forcing_ST_baseline( grid, climate_matrix)
    ! Use the ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing
    ! Initialise the baseline temperature

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                          INTENT(IN)    :: grid
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'initialise_climate_model_ISMIP_forcing_ST_baseline'
    TYPE(type_grid)                                         :: grid_raw
    REAL(dp), DIMENSION(:,:  ), POINTER                     ::  ST_raw
    INTEGER                                                 :: wST_raw

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire if everything we need is present in the file
    climate_matrix%ISMIP_forcing%netcdf_ST_baseline%filename = C%ISMIP_forcing_filename_ST_baseline
    CALL allocate_shared_int_0D( grid_raw%nx, grid_raw%wnx)
    CALL allocate_shared_int_0D( grid_raw%ny, grid_raw%wny)
    CALL inquire_ISMIP_forcing_ST_baseline_file( climate_matrix%ISMIP_forcing%netcdf_ST_baseline, grid_raw%nx, grid_raw%ny)
    CALL sync

    ! Allocate memory, read the data
    CALL allocate_shared_dp_0D(                           grid_raw%dx, grid_raw%wdx)
    CALL allocate_shared_dp_1D( grid_raw%nx,              grid_raw%x , grid_raw%wx )
    CALL allocate_shared_dp_1D(              grid_raw%ny, grid_raw%y , grid_raw%wy )
    CALL allocate_shared_dp_2D( grid_raw%nx, grid_raw%ny, ST_raw     , wST_raw     )
    CALL read_ISMIP_forcing_ST_baseline_file( climate_matrix%ISMIP_forcing%netcdf_ST_baseline, grid_raw%x, grid_raw%y, ST_raw)
    IF (par%master) grid_raw%dx = grid_raw%x( 2) - grid_raw%x( 1)
    CALL sync

    ! Transpose the data to [j,i] indexing
    CALL transpose_dp_2D( ST_raw, wST_raw)

    ! Map the data to the model grid
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, climate_matrix%ISMIP_forcing%ST_baseline, climate_matrix%ISMIP_forcing%wST_baseline)
    CALL map_square_to_square_cons_2nd_order_2D( grid_raw%nx, grid_raw%ny, grid_raw%x, grid_raw%y, &
      grid%nx, grid%ny, grid%x, grid%y, ST_raw, climate_matrix%ISMIP_forcing%ST_baseline)

    ! Clean up after yourself
    CALL deallocate_shared( grid_raw%wnx)
    CALL deallocate_shared( grid_raw%wny)
    CALL deallocate_shared( grid_raw%wdx)
    CALL deallocate_shared( grid_raw%wx )
    CALL deallocate_shared( grid_raw%wy )
    CALL deallocate_shared( wST_raw    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_ISMIP_forcing_ST_baseline

! == Some generally useful tools
! ==============================

  ! Map a global climate snapshot to the regional model grid
  SUBROUTINE map_snapshot_to_grid( grid, cglob, creg)
    ! Map data from a global climate snapshot to the regional grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_climate_snapshot_global),   INTENT(IN)    :: cglob  ! Global climate
    TYPE(type_climate_snapshot_regional), INTENT(INOUT) :: creg   ! grid   climate

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'map_snapshot_to_grid'
    INTEGER                                             :: i,j,m
    REAL(dp), DIMENSION(:,:,:), POINTER                 ::  logPrecip
    INTEGER                                             :: wlogPrecip

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To make sure extrapolation never results in negative precipitation, take the logarithm and inter/extrapolate that
    CALL allocate_shared_dp_3D( cglob%nlon, cglob%nlat, 12, logPrecip, wlogPrecip)
    DO i = cglob%i1, cglob%i2
    DO j = 1, cglob%nlat
    DO m = 1, 12
      logPrecip( i,j,m) = log(cglob%Precip( i,j,m))
    END DO
    END DO
    END DO
    CALL sync

    ! Perform the mapping operations
    CALL map_glob_to_grid_2D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%Hs,      creg%Hs     )
    CALL map_glob_to_grid_3D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%T2m,     creg%T2m    )
    CALL map_glob_to_grid_3D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, logPrecip,     creg%Precip )
    CALL map_glob_to_grid_3D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%Wind_WE, creg%Wind_WE)
    CALL map_glob_to_grid_3D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%Wind_SN, creg%Wind_SN)

    ! Get precipitation back from logarithm
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      creg%Precip( m,j,i) = EXP( creg%Precip( m,j,i))
    END DO
    END DO
    END DO
    CALL sync

    ! Rotate zonal/meridional wind to x,y wind
    CALL rotate_wind_to_model_grid( grid, creg%wind_WE, creg%wind_SN, creg%wind_LR, creg%wind_DU)

    ! Copy snapshot metadata
    creg%CO2        = cglob%CO2
    creg%orbit_time = cglob%orbit_time
    creg%orbit_ecc  = cglob%orbit_ecc
    creg%orbit_obl  = cglob%orbit_obl
    creg%orbit_pre  = cglob%orbit_pre
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wlogPrecip)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_snapshot_to_grid

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


END MODULE climate_module
