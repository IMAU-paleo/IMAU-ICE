MODULE forcing_module

  ! Contains all the routines for reading and calculating the model forcing (CO2, d18O, insolation, sea level),
  ! as well as the "forcing" structure which stores all the results from these routines (so that they
  ! can be accessed from all four ice-sheet models and the coupling routine).

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE data_types_module,               ONLY: type_forcing_data, type_model_region, type_grid, type_grid_lonlat, type_ice_model
  USE netcdf_extra_module,             ONLY: inquire_insolation_data_file, read_insolation_data_file_time_lat, read_insolation_data_file_timeframes
  USE netcdf_input_module,             ONLY: read_field_from_lonlat_file_2D, read_field_from_file_2D, read_field_from_file_history_1D
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             map_glob_to_grid_2D
  USE netcdf_debug_module,             ONLY: save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, save_variable_as_netcdf_int_3D, &
                                             save_variable_as_netcdf_dp_1D,  save_variable_as_netcdf_dp_2D,  save_variable_as_netcdf_dp_3D

  IMPLICIT NONE

  ! The data structure containing model forcing data - CO2 record, d18O record, (global) insolation record
  ! Updated at every coupling time step.
  TYPE(type_forcing_data), SAVE :: forcing

CONTAINS

! == Main routines that are called from IMAU_ICE_program
  SUBROUTINE update_global_forcing( NAM, EAS, GRL, ANT, time)
    ! Update global forcing data (d18O, CO2, insolation, geothermal heat flux)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(IN)    :: NAM, EAS, GRL, ANT
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_global_forcing'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Climate forcing stuff: CO2, d18O, inverse routine data
    IF     (C%choice_forcing_method == 'none') THEN
      ! Nothing needed; climate is either parameterised, or prescribed directly

    ELSEIF (C%choice_forcing_method == 'CO2_direct') THEN
      ! The global climate is calculated based on a prescribed CO2 record (e.g. from ice cores),
      ! either using a glacial-index method or a climate-matrix method, following Berends et al. (2018)
      
      IF (C%do_calculate_benthic_d18O) THEN
        CALL calculate_modelled_d18O( NAM, EAS, GRL, ANT)
      END IF

    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      ! The global climate is calculated using the observed present-day climate plus a global
      ! temperature offset, which is calculated using the inverse routine, following de Boer et al. (2014)

      CALL update_d18O_at_model_time( time)
      CALL update_global_mean_temperature_change_history( NAM, EAS, GRL, ANT)
      CALL calculate_modelled_d18O( NAM, EAS, GRL, ANT)
      CALL inverse_routine_global_temperature_offset

    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! The global climate is calculated based on modelled CO2, which follows from the inverse routine
      ! (following Berends et al., 2019). The climate itself can then be calculated using either a
      ! glacial-index method or a climate-matrix method

      CALL update_d18O_at_model_time( time)
      CALL update_global_mean_temperature_change_history( NAM, EAS, GRL, ANT)
      CALL calculate_modelled_d18O( NAM, EAS, GRL, ANT)
      CALL inverse_routine_CO2
    ELSE
      CALL crash('unknown choice_forcing_method "' // TRIM(C%choice_forcing_method) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_global_forcing
  SUBROUTINE initialise_global_forcing
    ! Initialise global forcing data (d18O, CO2, insolation, geothermal heat flux)

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_global_forcing'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise benthic d18O
     IF (C%do_calculate_benthic_d18O) THEN
       CALL initialise_modelled_benthic_d18O_data   
     END IF

    ! Climate forcing stuff: CO2, d18O, inverse routine data
    IF     (C%choice_forcing_method == 'none') THEN
      ! Nothing needed; climate is either parameterised, or prescribed directly

    ELSEIF (C%choice_forcing_method == 'CO2_direct') THEN
      ! The global climate is calculated based on a prescribed CO2 record (e.g. from ice cores),
      ! either using a glacial-index method or a climate-matrix method, following Berends et al. (2018)

      CALL initialise_CO2_record

    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      ! The global climate is calculated using the observed present-day climate plus a global
      ! temperature offset, which is calculated using the inverse routine, following de Boer et al. (2014)
      CALL initialise_d18O_record
      CALL initialise_inverse_routine_data

    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! The global climate is calculated based on modelled CO2, which follows from the inverse routine
      ! (following Berends et al., 2019). The climate itself can then be calculated using either a
      ! glacial-index method or a climate-matrix method
      CALL initialise_d18O_record      
      CALL initialise_inverse_routine_data

    ELSE
      CALL crash('unknown choice_forcing_method "' // TRIM(C%choice_forcing_method) // '"!')
    END IF

    ! Insolation
    CALL initialise_insolation_data

    ! Sea level
    IF (C%choice_sealevel_model == 'prescribed') THEN
      ! The global sea level is calculated based on a prescribed sea-level record (e.g. from stacks),
      CALL initialise_sealevel_record
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_global_forcing

! == Modelled benthic d18O
  SUBROUTINE calculate_modelled_d18O( NAM, EAS, GRL, ANT)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_model_region),             INTENT(IN)    :: NAM, EAS, GRL, ANT

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calculate_modelled_d18O'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%do_calculate_benthic_d18O) THEN
       CALL crash('should only be called when do_calculate_benthic_d18O = .TRUE.!')
    END IF
    IF (C%choice_ice_isotopes_model == 'none') THEN
      CALL crash('choice_ice_isotopes_model = none; cannot calculate d18O contribution of ice sheets when no englacial isotope calculation is being done!')
    END IF

    IF (par%master) THEN

      ! Determine contributions to benthic d18O from the four ice sheets
      forcing%d18O_NAM = 0._dp
      forcing%d18O_EAS = 0._dp
      forcing%d18O_GRL = 0._dp
      forcing%d18O_ANT = 0._dp
      IF (C%do_NAM) forcing%d18O_NAM = NAM%d18O_contribution - NAM%d18O_contribution_PD
      IF (C%do_EAS) forcing%d18O_EAS = EAS%d18O_contribution - EAS%d18O_contribution_PD
      IF (C%do_GRL) forcing%d18O_GRL = GRL%d18O_contribution - GRL%d18O_contribution_PD
      IF (C%do_ANT) forcing%d18O_ANT = ANT%d18O_contribution - ANT%d18O_contribution_PD
      forcing%d18O_from_ice_volume_mod = forcing%d18O_NAM + forcing%d18O_EAS + forcing%d18O_GRL + forcing%d18O_ANT

      ! Determine contributions to benthic d18O from ocean temperature,
      ! which we assume scales with global mean surface temperature
      forcing%d18O_from_temperature_mod = forcing%dT_deepwater * C%d18O_dT_deepwater_ratio

      ! Determine total change in modelled benthic d18O
      forcing%d18O_mod = forcing%d18O_obs_PD + forcing%d18O_from_ice_volume_mod + forcing%d18O_from_temperature_mod

    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_modelled_d18O
  SUBROUTINE update_global_mean_temperature_change_history( NAM, EAS, GRL, ANT)
    ! Calculate the annual mean surface temperature change w.r.t PD for all
    ! model regions, and calculate a weighted average to get a "global" (high-latitude) anomaly.
    ! Add this value to the global temperature anomaly history, and average of the moving
    ! time window (and scale with a factor) to get the deep-water temperature anomaly
    ! (which is needed to calculate benthic d18O).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(IN)    :: NAM, EAS, GRL, ANT

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_global_mean_temperature_change_history'
    REAL(dp)                                           ::  dT_NAM,  dT_EAS,  dT_GRL,  dT_ANT
    REAL(dp)                                           :: dT_glob_average_over_window
    REAL(dp)                                           :: n_glob

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine annual mean surface temperature change for all model regions
    IF (C%do_NAM) CALL calculate_mean_temperature_change_region( NAM, dT_NAM)
    IF (C%do_EAS) CALL calculate_mean_temperature_change_region( EAS, dT_EAS)
    IF (C%do_GRL) CALL calculate_mean_temperature_change_region( GRL, dT_GRL)
    IF (C%do_ANT) CALL calculate_mean_temperature_change_region( ANT, dT_ANT)

    IF (par%master) THEN

      ! Weighted average of all model regions
      forcing%dT_glob = 0._dp
      n_glob          = 0

      ! Calculate the mean global temperature using the weighted average (across AREA) of the domains
      IF (C%do_NAM) THEN
        forcing%dT_glob = forcing%dT_glob + (dT_NAM * NAM%grid%nx * NAM%grid%ny * NAM%grid%dx * NAM%grid%dx)
        n_glob  = n_glob  + NAM%grid%nx * NAM%grid%ny * NAM%grid%dx * NAM%grid%dx
      END IF
      IF (C%do_EAS) THEN
        forcing%dT_glob = forcing%dT_glob + (dT_EAS * EAS%grid%nx * EAS%grid%ny * EAS%grid%dx * EAS%grid%dx)
        n_glob  = n_glob  + EAS%grid%nx * EAS%grid%ny * EAS%grid%dx * EAS%grid%dx
      END IF
      IF (C%do_GRL) THEN
        forcing%dT_glob = forcing%dT_glob + (dT_GRL * GRL%grid%nx * GRL%grid%ny * GRL%grid%dx * GRL%grid%dx)
        n_glob  = n_glob  + GRL%grid%nx * GRL%grid%ny * GRL%grid%dx * GRL%grid%dx
      END IF
      IF (C%do_ANT) THEN
        forcing%dT_glob = forcing%dT_glob + (dT_ANT * ANT%grid%nx * ANT%grid%ny * ANT%grid%dx * ANT%grid%dx)
        n_glob  = n_glob  + ANT%grid%nx * ANT%grid%ny * ANT%grid%dx * ANT%grid%dx
      END IF

      forcing%dT_glob = forcing%dT_glob / n_glob

      ! Update the global mean temperature change history
      ! 1st entry is the current value, 2nd is 1*dt_coupling ago, 3d is 2*dt_coupling ago, etc.
      forcing%dT_glob_history( 2:forcing%ndT_glob_history) = forcing%dT_glob_history( 1:forcing%ndT_glob_history-1)
      forcing%dT_glob_history( 1) = forcing%dT_glob

      ! Calculate deep-water temperature anomaly by averaging over the specified time window (default value 3000 years)
      ! and scaling with the specified factor (default value 0.25)
      dT_glob_average_over_window = SUM( forcing%dT_glob_history) / REAL( forcing%ndT_glob_history,dp)
      forcing%dT_deepwater = C%dT_deepwater_dT_surf_ratio * dT_glob_average_over_window

    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_global_mean_temperature_change_history
  SUBROUTINE calculate_mean_temperature_change_region( region, dT)
    ! Calculate the annual mean surface temperature change w.r.t PD (corrected for elevation changes) over a model region

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(IN)    :: region
    REAL(dp),                            INTENT(OUT)   :: dT

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calculate_mean_temperature_change_region'
    INTEGER                                            :: i,j,m
    REAL(dp), DIMENSION(:,:  ), POINTER                :: Hs_PD
    REAL(dp), DIMENSION(:,:,:), POINTER                :: T2m_PD
    INTEGER                                            :: wHs_PD, wT2m_PD
    REAL(dp)                                           :: dT_lapse_mod, dT_lapse_PD, T_pot_mod, T_pot_PD
    INTEGER                                            :: n_dT

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     region%grid%ny, region%grid%nx, Hs_PD , wHs_PD )
    CALL allocate_shared_dp_3D( 12, region%grid%ny, region%grid%nx, T2m_PD, wT2m_PD)

    ! Obtain present-day orography and temperature from the climate model
    IF     (C%choice_climate_model == 'PD_obs') THEN
      Hs_PD(    :,region%grid%i1:region%grid%i2) = region%climate%PD_obs%snapshot%Hs(    :,region%grid%i1:region%grid%i2)
      T2m_PD( :,:,region%grid%i1:region%grid%i2) = region%climate%PD_obs%snapshot%T2m( :,:,region%grid%i1:region%grid%i2)
    ELSEIF (C%choice_climate_model == 'matrix') THEN
      Hs_PD(    :,region%grid%i1:region%grid%i2) = region%climate%matrix%PD_obs%Hs(    :,region%grid%i1:region%grid%i2)
      T2m_PD( :,:,region%grid%i1:region%grid%i2) = region%climate%matrix%PD_obs%T2m( :,:,region%grid%i1:region%grid%i2)
    ELSE
      CALL crash('present-day temperature not known for choice_climate_model "'//TRIM( C%choice_climate_model)//'"!')
    END IF

    dT   = 0._dp
    n_dT = 0

    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny      
    ! Do not calculate dT_glob over no_ice regions (assume difference is 0)
      IF (region%mask_noice( j,i) == 0) THEN

        ! Correct to sea level to prevent temperature differences from topography
        dT_lapse_mod = region%ice%Hs_a(          j,i) * C%constant_lapserate
        dT_lapse_PD  = Hs_PD(                    j,i) * C%constant_lapserate
        DO m = 1, 12
          n_dT = n_dT + 1

          T_pot_mod = region%climate%T2m(         m,j,i) + dT_lapse_mod
          T_pot_PD  = T2m_PD(                     m,j,i) + dT_lapse_PD
          dT = dT + T_pot_mod - T_pot_PD
        END DO
      END IF
    END DO
    END DO
    CALL sync

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dT,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_dT, 1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate the mean temperature difference
    IF (par%master) THEN
      IF (n_dT == 0) THEN ! We don't want to divide by 0!
        dT = 0 
      ELSE
        dT = dT / REAL(n_dT,dp)
      END IF
    END IF
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wHs_PD )
    CALL deallocate_shared( wT2m_PD)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_mean_temperature_change_region
  SUBROUTINE initialise_modelled_benthic_d18O_data
    ! Allocate shared memory for the d18O and global temperature variables.

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_modelled_benthic_d18O_data'
    CHARACTER(LEN=256)                                 :: filename
    CHARACTER(LEN=256)                                 :: choice_d18O_inverse_init
    REAL(dp)                                           :: time_to_restart_from

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. C%do_calculate_benthic_d18O) THEN  
       CALL crash('should only be called when do_calculate_benthic_d18O = .TRUE.!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_0D( forcing%dT_glob,                   forcing%wdT_glob                  )
    CALL allocate_shared_dp_0D( forcing%dT_deepwater,              forcing%wdT_deepwater             )
    CALL allocate_shared_dp_0D( forcing%d18O_NAM,                  forcing%wd18O_NAM                 )
    CALL allocate_shared_dp_0D( forcing%d18O_EAS,                  forcing%wd18O_EAS                 )
    CALL allocate_shared_dp_0D( forcing%d18O_GRL,                  forcing%wd18O_GRL                 )
    CALL allocate_shared_dp_0D( forcing%d18O_ANT,                  forcing%wd18O_ANT                 )
    CALL allocate_shared_dp_0D( forcing%d18O_from_ice_volume_mod,  forcing%wd18O_from_ice_volume_mod )
    CALL allocate_shared_dp_0D( forcing%d18O_from_temperature_mod, forcing%wd18O_from_temperature_mod)
    CALL allocate_shared_dp_0D( forcing%d18O_obs,                  forcing%wd18O_obs                 )
    CALL allocate_shared_dp_0D( forcing%d18O_obs_PD,               forcing%wd18O_obs_PD              )
    CALL allocate_shared_dp_0D( forcing%d18O_mod,                  forcing%wd18O_mod                 )

    ! Determine which type of d18O initialization, and which filename
    IF  (C%do_NAM) THEN
      choice_d18O_inverse_init = C%choice_d18O_inverse_init_NAM
      filename                 = C%filename_d18O_inverse_init_NAM
      time_to_restart_from     = C%time_to_restart_from_NAM
    ELSEIF (C%do_EAS) THEN
      choice_d18O_inverse_init = C%choice_d18O_inverse_init_EAS
      filename                 = C%filename_d18O_inverse_init_EAS
      time_to_restart_from     = C%time_to_restart_from_EAS
    ELSEIF (C%do_GRL) THEN
      choice_d18O_inverse_init = C%choice_d18O_inverse_init_GRL
      filename                 = C%filename_d18O_inverse_init_GRL
      time_to_restart_from     = C%time_to_restart_from_GRL
    ELSEIF (C%do_ANT) THEN
      choice_d18O_inverse_init = C%choice_d18O_inverse_init_ANT
      filename                 = C%filename_d18O_inverse_init_ANT
      time_to_restart_from     = C%time_to_restart_from_ANT
    END IF

    ! Determine number of entries in the global mean temperature change history
    CALL allocate_shared_int_0D( forcing%ndT_glob_history, forcing%wndT_glob_history)
    IF (par%master) forcing%ndT_glob_history = CEILING( C%dT_deepwater_averaging_window / C%dt_coupling)
    CALL sync

    IF (par%master) forcing%d18O_obs_PD = 3.23_dp
    CALL sync

    IF (choice_d18O_inverse_init == 'init') THEN
       ! Starting from a clean run

       ! Allocate memory for the global mean temperature change history
       CALL allocate_shared_dp_1D( forcing%ndT_glob_history, forcing%dT_glob_history, forcing%wdT_glob_history)
        
       IF (par%master) forcing%dT_glob_history = 0._dp
       
       CALL sync

    ELSEIF (choice_d18O_inverse_init == 'restart') THEN
        ! If starting from a restart file
        
        ! Load from restart file
        CALL read_field_from_file_history_1D(         filename, 'dT_glob_history',    'time_dT_glob_history',     forcing%dT_glob_history     , forcing%wdT_glob_history,     time_to_restart_from )

        ! Crash if the length of the inverse history does not match between restart and model 
        IF (SIZE(forcing%dT_glob_history) /= forcing%ndT_glob_history) THEN
            CALL crash('Length of dT_glob_history is not the same between current model set-up and restart. Check dt_coupling and CO2_averaging_window')
        END IF
        
        CALL sync

    ELSE
      CALL crash('unknown choice_forcing_method "' // TRIM(choice_d18O_inverse_init) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_modelled_benthic_d18O_data

! == Inverse forward routine
  SUBROUTINE inverse_routine_global_temperature_offset
    ! Use the inverse routine to calculate a global temperature offset
    ! (i.e. the method used in de Boer et al., 2013)
    ! (de Boer, B., van de Wal, R., Lourens, L. J., Bintanja, R., and Reerink, T. J.:
    ! A continuous simulation of global ice volume over the past 1 million years with 3-D ice-sheet models, Climate Dynamics 41, 1365-1384, 2013)

    ! Local variables:
    REAL(dp)                                           :: dT_glob_inverse_average_over_window

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inverse_routine_global_temperature_offset'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! The inverse routine might not work properly when not all ice sheets are simulated
    IF (((.NOT. C%do_NAM) .OR. (.NOT. C%do_EAS) .OR. (.NOT. C%do_GRL) .OR. (.NOT. C%do_ANT)) .AND. par%master) THEN
      CALL warning('the inverse routine only works properly when all four ice sheets are simulated! ' // &
                  'Leaving one out means you will miss that contribution to the d18O, which the ' // &
                  'routine will try to compensate for by making the world colder. ' // &
                  'If you really want to simulate only one ice sheet, consider using direct CO2 forcing.')
    END IF

    IF (par%master) THEN

      ! Average dT_glob_inverse over the moving time window
      dT_glob_inverse_average_over_window = SUM( forcing%dT_glob_inverse_history) / REAL(forcing%ndT_glob_inverse_history,dp)

      ! Update dT_glob_inverse based on the difference between modelled and observed d18O
      forcing%dT_glob_inverse = dT_glob_inverse_average_over_window + (forcing%d18O_mod - forcing%d18O_obs) * C%inverse_d18O_to_dT_glob_scaling

      ! Update the moving time window
      forcing%dT_glob_inverse_history( 2:forcing%ndT_glob_inverse_history) = forcing%dT_glob_inverse_history( 1:forcing%ndT_glob_inverse_history-1)
      forcing%dT_glob_inverse_history( 1) = forcing%dT_glob_inverse

    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inverse_routine_global_temperature_offset
  SUBROUTINE inverse_routine_CO2
    ! Use the inverse routine to calculate modelled CO2
    ! (i.e. the method used in Berends et al., 2019)
    ! (Berends, C. J., de Boer, B., Dolan, A. M., Hill, D. J., and van de Wal, R. S. W.: Modelling ice sheet evolution and atmospheric CO2 during the Late Pliocene, Climate of the Past 15, 1603-1619, 2019)

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inverse_routine_CO2'
    REAL(dp)                                           :: CO2_inverse_average_over_window

    ! Add routine to path
    CALL init_routine( routine_name)

    ! The inverse routine might not work properly when not all ice sheets are simulated
    IF (((.NOT. C%do_NAM) .OR. (.NOT. C%do_EAS) .OR. (.NOT. C%do_GRL) .OR. (.NOT. C%do_ANT)) .AND. par%master) THEN
      CALL warning('the inverse routine only works properly when all four ice sheets are simulated! ' // &
                  'Leaving one out means you will miss that contribution to the d18O, which the ' // &
                  'routine will try to compensate for by making the world colder. ' // &
                  'If you really want to simulate only one ice sheet, consider using direct CO2 forcing.')
    END IF

    IF (par%master) THEN

      ! Average CO2_inverse over the moving time window
      CO2_inverse_average_over_window = SUM( forcing%CO2_inverse_history) / REAL(forcing%nCO2_inverse_history,dp)

      ! Update CO2_inverse based on the difference between modelled and observed d18O
      forcing%CO2_inverse = CO2_inverse_average_over_window + (forcing%d18O_mod - forcing%d18O_obs) * C%inverse_d18O_to_CO2_scaling

      ! Update the moving time window
      forcing%CO2_inverse_history( 2:forcing%nCO2_inverse_history) = forcing%CO2_inverse_history( 1:forcing%nCO2_inverse_history-1)
      forcing%CO2_inverse_history( 1) = forcing%CO2_inverse

      forcing%CO2_mod = forcing%CO2_inverse

      ! Keep track of the new CO2 concentration
      WRITE(0,'(A,F9.3,A)') '  - CO2 concentration has updated to ', forcing%CO2_mod, ' ppm'


    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inverse_routine_CO2
  SUBROUTINE initialise_inverse_routine_data
    ! Allocate shared memory for the moving time windows used in the inverse routine

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_inverse_routine_data'
    CHARACTER(LEN=256)                                 :: filename
    CHARACTER(LEN=256)                                 :: choice_d18O_inverse_init
    REAL(dp)                                           :: time_to_restart_from
    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN

      CALL crash('need to fix the inverse routine stuff to cope with restarting!')

!      CALL allocate_shared_dp_0D( forcing%dT_glob_inverse, forcing%wdT_glob_inverse)
!      ! Determine number of entries in the history
!      CALL allocate_shared_int_0D( forcing%ndT_glob_inverse_history, forcing%wndT_glob_inverse_history)
!      IF (par%master) forcing%ndT_glob_inverse_history = CEILING( C%dT_glob_inverse_averaging_window / C%dt_coupling)
!      CALL sync
!      ! Allocate memory for the global mean temperature change history
!      CALL allocate_shared_dp_1D( forcing%ndT_glob_inverse_history, forcing%dT_glob_inverse_history, forcing%wdT_glob_inverse_history)
!      IF (par%master) forcing%dT_glob_inverse_history = 0._dp
!      IF (par%master) forcing%dT_glob_inverse         = 0._dp
!      IF (par%master) forcing%CO2_mod                 = 0._dp
!
!      ! If we're restarting a previous run, read inverse routine history from one of the restart files
!      IF (C%is_restart) THEN
!        IF (C%do_NAM) THEN
!          filename = C%filename_init_NAM
    ! Allocate 
!        ELSEIF (C%do_EAS) THEN
!          filename = C%filename_init_EAS
!        ELSEIF (C%do_GRL) THEN
!          filename = C%filename_init_GRL
!        ELSEIF (C%do_ANT) THEN
!          filename = C%filename_init_ANT
!        END IF
!        IF (par%master) CALL read_inverse_routine_history_dT_glob(         forcing, C%filename_init_NAM)
!        IF (par%master) CALL read_inverse_routine_history_dT_glob_inverse( forcing, C%filename_init_NAM)
!        IF (par%master) forcing%dT_glob_inverse = forcing%dT_glob_inverse_history(1)
!        CALL sync
!      END IF

    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN

      ! Determine which type of d18O initialization, and which filename
      IF  (C%do_NAM) THEN
        choice_d18O_inverse_init = C%choice_d18O_inverse_init_NAM
        filename                 = C%filename_d18O_inverse_init_NAM
        time_to_restart_from     = C%time_to_restart_from_NAM
      ELSEIF (C%do_EAS) THEN
        choice_d18O_inverse_init = C%choice_d18O_inverse_init_EAS
        filename                 = C%filename_d18O_inverse_init_EAS
        time_to_restart_from     = C%time_to_restart_from_EAS
      ELSEIF (C%do_GRL) THEN
        choice_d18O_inverse_init = C%choice_d18O_inverse_init_GRL
        filename                 = C%filename_d18O_inverse_init_GRL
        time_to_restart_from     = C%time_to_restart_from_GRL
      ELSEIF (C%do_ANT) THEN
        choice_d18O_inverse_init = C%choice_d18O_inverse_init_ANT
        filename                 = C%filename_d18O_inverse_init_ANT
        time_to_restart_from     = C%time_to_restart_from_ANT
      END IF

      ! Allocate 
      CALL allocate_shared_dp_0D(  forcing%CO2_inverse,              forcing%wCO2_inverse)
      CALL allocate_shared_dp_0D(  forcing%CO2_mod,                  forcing%wCO2_mod)
      CALL allocate_shared_int_0D( forcing%nCO2_inverse_history,     forcing%wnCO2_inverse_history)
      CALL allocate_shared_dp_0D(  forcing%dT_glob_inverse,          forcing%wdT_glob_inverse)
      CALL allocate_shared_int_0D( forcing%ndT_glob_inverse_history, forcing%wndT_glob_inverse_history)

      ! Determine number of entries in the history
      IF (par%master) forcing%nCO2_inverse_history     = CEILING( C%CO2_inverse_averaging_window / C%dt_coupling)
      IF (par%master) forcing%ndT_glob_inverse_history = CEILING( C%dT_glob_inverse_averaging_window / C%dt_coupling)
      CALL sync

      IF (choice_d18O_inverse_init == 'init') THEN
        ! Starting from the beginning, so we use a constant value

        ! Allocate
        CALL allocate_shared_dp_1D(  forcing%nCO2_inverse_history, forcing%CO2_inverse_history, forcing%wCO2_inverse_history)
        CALL allocate_shared_dp_1D(  forcing%ndT_glob_inverse_history, forcing%dT_glob_inverse_history, forcing%wdT_glob_inverse_history)

        ! Add a uniform value for the inverse history
        IF (par%master) forcing%CO2_inverse_history = C%inverse_d18O_to_CO2_initial_CO2 
        IF (par%master) forcing%CO2_inverse         = C%inverse_d18O_to_CO2_initial_CO2
        IF (par%master) forcing%CO2_mod             = C%inverse_d18O_to_CO2_initial_CO2

        CALL sync

      ELSEIF (choice_d18O_inverse_init == 'restart') THEN
        ! If starting from a restart file
        ! Read dT_glob
        CALL read_field_from_file_history_1D(         filename, 'CO2_inverse_history',    'time_CO2_inverse_history',     forcing%CO2_inverse_history     , forcing%wCO2_inverse_history,     time_to_restart_from )
       
        ! Crash if the length of the inverse history does not match between restart and model 
        IF (SIZE(forcing%CO2_inverse_history) /= forcing%nCO2_inverse_history) THEN
            CALL crash('Length of CO2_inverse_history is not the same between current model set-up and restart. Check dt_coupling and CO2_averaging_window')
        END IF       
    
        ! Get the current CO2_inverse and CO2_mod from the restart file
        IF (par%master) forcing%CO2_inverse = forcing%CO2_inverse_history(1)
        IF (par%master) forcing%CO2_mod     = forcing%CO2_inverse 
        
        CALL sync
      
      ELSE
          CALL crash('unknown choice_inverse_clim_init "' // TRIM(choice_d18O_inverse_init) // '"!')
      END IF

    ELSE
      CALL crash('unknown choice_forcing_method "' // TRIM(C%choice_forcing_method) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_inverse_routine_data

! == Prescribed CO2 record
  SUBROUTINE update_CO2_at_model_time( time)
    ! Interpolate the data in forcing%CO2 to find the value at the queried time.
    ! If time lies outside the range of forcing%CO2_time, return the first/last value
    !
    ! NOTE: assumes time is listed in yr BP, so LGM would be -21000.0, and 0.0 corresponds to January 1st 1900.
    !
    ! NOTE: calculates average value over the preceding 30 years. For paleo this doesn't matter
    !       in the least, but for the historical period this makes everything more smooth.

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_CO2_at_model_time'
    INTEGER                                            :: ti1, ti2, til, tiu
    REAL(dp)                                           :: a, b, tl, tu, intCO2, dintCO2
    REAL(dp), PARAMETER                                :: dt_smooth = 60._dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      ! Observed CO2 is needed for these forcing methods.
    ELSE
      CALL crash('should only be called when choice_forcing_method = "CO2_direct"!')
    END IF

    IF (par%master) THEN

      IF     (time < MINVAL( forcing%CO2_time)) THEN ! times 1000 because forcing%CO2_time is in kyr
        ! Model time before start of CO2 record; using constant extrapolation
        forcing%CO2_obs = forcing%CO2_record( 1)
      ELSEIF (time > MAXVAL( forcing%CO2_time)) THEN
        ! Model time beyond end of CO2 record; using constant extrapolation
        forcing%CO2_obs = forcing%CO2_record( C%CO2_record_length)
      ELSE

        ! Find range of raw time frames enveloping model time
        ti1 = 1
        DO WHILE (forcing%CO2_time( ti1) < time - dt_smooth .AND. ti1 < C%CO2_record_length)
          ti1 = ti1 + 1
        END DO
        ti1 = MAX( 1, ti1 - 1)

        ti2 = 2
        DO WHILE (forcing%CO2_time( ti2) < time             .AND. ti2 < C%CO2_record_length)
          ti2 = ti2 + 1
        END DO

        ! Calculate conservatively-remapped time-averaged CO2
        intCO2 = 0._dp
        DO til = ti1, ti2 - 1
          tiu = til + 1

          ! Linear interpolation between til and tiu: CO2( t) = a + b*t
          b = (forcing%CO2_record( tiu) - forcing%CO2_record( til)) / (forcing%CO2_time( tiu) - forcing%CO2_time( til))
          a = forcing%CO2_record( til) - b*forcing%CO2_time( til)

          ! Window of overlap between [til,tiu] and [t - dt_smooth, t]
          tl = MAX( forcing%CO2_time( til), time - dt_smooth)
          tu = MIN( forcing%CO2_time( tiu), time            )
          dintCO2 = (tu - tl) * (a + b * (tl + tu) / 2._dp)
          intCO2 = intCO2 + dintCO2
        END DO
        forcing%CO2_obs = intCO2 / dt_smooth

      END IF

    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_CO2_at_model_time
  SUBROUTINE initialise_CO2_record
    ! Read the CO2 record specified in C%filename_CO2_record. Assumes this is an ASCII text file with at least two columns (time in kyr and CO2 in ppmv)
    ! and the number of rows being equal to C%CO2_record_length

    ! NOTE: assumes time is listed in kyr BP (so LGM would be -21.0); converts to yr after reading!

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_CO2_record'
    INTEGER                                            :: i,ios

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      ! Observed CO2 is needed for these forcing methods.
    ELSE
      CALL crash('should only be called when choice_forcing_method = "CO2_direct"!')
    END IF

    ! Allocate shared memory to take the data
    CALL allocate_shared_dp_1D( C%CO2_record_length, forcing%CO2_time,   forcing%wCO2_time  )
    CALL allocate_shared_dp_1D( C%CO2_record_length, forcing%CO2_record, forcing%wCO2_record)
    CALL allocate_shared_dp_0D(                      forcing%CO2_obs,    forcing%wCO2_obs   )

    ! Read CO2 record (time and values) from specified text file
    IF (par%master) THEN

      WRITE(0,*) ' Reading CO2 record from ', TRIM(C%filename_CO2_record), '...'

      OPEN(   UNIT = 1337, FILE=C%filename_CO2_record, ACTION='READ')

      DO i = 1, C%CO2_record_length
        READ( UNIT = 1337, FMT=*, IOSTAT=ios) forcing%CO2_time( i), forcing%CO2_record( i)
        IF (ios /= 0) THEN
          CALL crash('length of text file "' // TRIM(C%filename_CO2_record) // '" does not match C%CO2_record_length!')
        END IF
      END DO

      CLOSE( UNIT  = 1337)

      IF (C%start_time_of_run/1000._dp < forcing%CO2_time(1)) THEN
         CALL warning(' Model time starts before start of CO2 record; constant extrapolation will be used in that case!')
      END IF
      IF (C%end_time_of_run/1000._dp > forcing%CO2_time(C%CO2_record_length)) THEN
         CALL warning(' Model time will reach beyond end of CO2 record; constant extrapolation will be used in that case!')
      END IF

      ! Convert from kyr to yr
      forcing%CO2_time = forcing%CO2_time * 1000._dp

    END IF ! IF (par%master)
    CALL sync

    ! Set the value for the current (starting) model time
    CALL update_CO2_at_model_time( C%start_time_of_run)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_CO2_record

! == Prescribed d18O record
  SUBROUTINE update_d18O_at_model_time( time)
    ! Interpolate the data in forcing%d18O to find the value at the queried time.
    ! If time lies outside the range of forcing%d18O_time, return the first/last value
    ! NOTE: assumes time is listed in years! (so LGM would be -21000.0)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_d18O_at_model_time'
    INTEGER                                            :: il, iu
    REAL(dp)                                           :: wl, wu

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF     (C%choice_forcing_method == 'd18O_inverse_dT_glob' .OR. &
            C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! Observed d18O is needed for these forcing methods.
    ELSE
      CALL crash('should only be called when choice_forcing_method = "d18O_inverse_dT_glob" or "d18O_inverse_CO2"!')
    END IF

    IF (par%master) THEN

      IF     (time < MINVAL( forcing%d18O_time)) THEN ! times 1000 because forcing%d18O_time is in kyr
        CALL warning('model time before start of d18O record; using constant extrapolation!')
        forcing%d18O_obs = forcing%d18O_record( 1)
      ELSEIF (time > MAXVAL( forcing%d18O_time)) THEN
        CALL warning('model time beyond end of d18O record; using constant extrapolation!')
        forcing%d18O_obs = forcing%d18O_record( C%d18O_record_length)
      ELSE
        iu = 1
        DO WHILE (forcing%d18O_time(iu) < time)
          iu = iu+1
        END DO
        il = iu-1

        wl = (forcing%d18O_time(iu) - time) / (forcing%d18O_time(iu)-forcing%d18O_time(il))
        wu = (time - forcing%d18O_time(il)) / (forcing%d18O_time(iu)-forcing%d18O_time(il))

        forcing%d18O_obs = forcing%d18O_record(il) * wl + forcing%d18O_record(iu) * wu

      END IF

    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_d18O_at_model_time
  SUBROUTINE initialise_d18O_record
    ! Read the d18O record specified in C%filename_d18O_record. Assumes this is an ASCII text file with at least two columns (time in yr and d18O in per mil)
    ! and the number of rows being equal to C%d18O_record_length
    ! NOTE: assumes time is listed in years! (so LGM would be -21000.0)

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_d18O_record'
    INTEGER                                            :: i,ios

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF     (C%choice_forcing_method == 'd18O_inverse_dT_glob' .OR. &
             C%choice_forcing_method == 'd18O_inverse_CO2') THEN
       ! Observed d18O is needed for these forcing methods.
    ELSE
       CALL crash('should only be called when choice_forcing_method = "d18O_inverse_dT_glob" or "d18O_inverse_CO2"!')
    END IF

    IF (par%master) WRITE(0,*) ' Reading d18O record from ', TRIM(C%filename_d18O_record), '...'

    ! Allocate shared memory to take the data
    CALL allocate_shared_dp_1D( C%d18O_record_length, forcing%d18O_time,   forcing%wd18O_time  )
    CALL allocate_shared_dp_1D( C%d18O_record_length, forcing%d18O_record, forcing%wd18O_record)
    CALL allocate_shared_dp_0D(                       forcing%d18O_obs,    forcing%wd18O_obs   )

    ! Read d18O record (time and values) from specified text file
    IF (par%master) THEN

      OPEN(   UNIT = 1337, FILE=C%filename_d18O_record, ACTION='READ')
      DO i = 1, C%d18O_record_length
        READ( UNIT = 1337, FMT=*, IOSTAT=ios) forcing%d18O_time(i), forcing%d18O_record(i)
        IF (ios /= 0) THEN
          CALL crash('length of text file "' // TRIM(C%filename_d18O_record) // '" does not match C%d18O_record_length!')
        END IF
      END DO

      CLOSE( UNIT  = 1337)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Set the PD value
    CALL update_d18O_at_model_time( 0._dp)
    forcing%d18O_obs_PD = forcing%d18O_obs

    ! Set the value for the current (starting) model time
    CALL update_d18O_at_model_time( C%start_time_of_run)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_d18O_record

! == Insolation

SUBROUTINE get_insolation_at_time( grid, time, Q_TOA)
    ! Get monthly insolation at time t on the regional grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp),                            INTENT(IN)    :: time
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: Q_TOA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'get_insolation_at_time'
    REAL(dp)                                           :: time_applied
    INTEGER                                            :: i,j,m,ilat_l,ilat_u
    REAL(dp)                                           :: wt0, wt1, wlat_l, wlat_u
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  Q_TOA_int
    INTEGER                                            :: wQ_TOA_int

    ! Add routine to path
    CALL init_routine( routine_name)

    time_applied = 0._dp

    ! Safety
    IF     (C%choice_insolation_forcing == 'none') THEN
      CALL crash('insolation should not be used when choice_insolation_forcing = "none"!')
    ELSEIF (C%choice_insolation_forcing == 'static') THEN
      time_applied = C%static_insolation_time
    ELSEIF (C%choice_insolation_forcing == 'realistic') THEN
      time_applied = time
    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time_applied < forcing%ins_t0 .OR. time_applied > forcing%ins_t1) THEN
      CALL sync
      CALL update_insolation_timeframes_from_file( time_applied)
    END IF

    ! Allocate shared memory for timeframe-interpolated lat-month-only insolation
    CALL allocate_shared_dp_2D( forcing%ins_nlat, 12, Q_TOA_int, wQ_TOA_int)

    ! Calculate timeframe interpolation weights
    wt0 = (forcing%ins_t1 - time_applied) / (forcing%ins_t1 - forcing%ins_t0)
    wt1 = 1._dp - wt0

    ! Interpolate the two timeframes
    IF (par%master) THEN
      Q_TOA_int = wt0 * forcing%ins_Q_TOA0 + wt1 * forcing%ins_Q_TOA1
    END IF
    CALL sync

    ! Map the timeframe-interpolated lat-month-only insolation to the model grid
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ilat_l = FLOOR(grid%lat(j,i) + 91)
      ilat_u = ilat_l + 1

      wlat_l = forcing%ins_lat(ilat_u) - grid%lat(j,i)
      wlat_u = 1._dp - wlat_l

      DO m = 1, 12
        Q_TOA( m,j,i) = wlat_l * Q_TOA_int( ilat_l,m) + wlat_u * Q_TOA_int( ilat_u,m)
      END DO

    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wQ_TOA_int)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_insolation_at_time
  SUBROUTINE get_insolation_at_time_month_and_lat( time, month, lat, Q_TOA)
    ! Get monthly insolation at time t, month m and latitude l on the regional grid

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time
    INTEGER,                             INTENT(IN)    :: month
    REAL(dp),                            INTENT(IN)    :: lat
    REAL(dp),                            INTENT(OUT)   :: Q_TOA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'get_insolation_at_time_month_and_lat'
    REAL(dp)                                           :: time_applied
    INTEGER                                            :: ilat_l,ilat_u
    REAL(dp)                                           :: wt0, wt1, wlat_l, wlat_u

    ! Add routine to path
    CALL init_routine( routine_name)

    time_applied = 0._dp

    ! Safety
    IF     (C%choice_insolation_forcing == 'none') THEN
      CALL crash('insolation should not be used when choice_insolation_forcing = "none"!')
    ELSEIF (C%choice_insolation_forcing == 'static') THEN
      time_applied = C%static_insolation_time
    ELSEIF (C%choice_insolation_forcing == 'realistic') THEN
      time_applied = time
    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time_applied < forcing%ins_t0 .AND. time_applied > forcing%ins_t1) THEN
      CALL update_insolation_timeframes_from_file( time_applied)
    END IF

    ! Calculate timeframe interpolation weights
    wt0 = (forcing%ins_t1 - time) / (forcing%ins_t1 - forcing%ins_t0)
    wt1 = 1._dp - wt0

    ! Get value at month m and latitude l
    ilat_l = FLOOR(lat + 91)
    ilat_u = ilat_l + 1

    wlat_l = (forcing%ins_lat( ilat_u) - lat) / (forcing%ins_lat( ilat_u) - forcing%ins_lat( ilat_l))
    wlat_u = 1._dp - wlat_l

    IF (par%master) Q_TOA = wt0 * wlat_l * forcing%ins_Q_TOA0( ilat_l,month) + &
                            wt0 * wlat_u * forcing%ins_Q_TOA0( ilat_u,month) + &
                            wt1 * wlat_l * forcing%ins_Q_TOA1( ilat_l,month) + &
                            wt1 * wlat_u * forcing%ins_Q_TOA1( ilat_u,month)
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_insolation_at_time_month_and_lat
  
  SUBROUTINE get_summer_insolation_at_time( time)
    ! Obtain the summer insolation in the Northern (65N) and Southern (80S) Hemispheres

    IMPLICIT NONE

    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'get_summer_insolation_at_time'
    REAL(dp)                                           :: Q_TOA_JJA_65N_raw, Q_TOA_DJF_80S_raw
    INTEGER                                            :: m
    
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Northern Hemisphere (June, July, August)
    forcing%Q_TOA_JJA_65N = 0._dp
  
    DO m = 6,8
      CALL get_insolation_at_time_month_and_lat( time, m, 65._dp,  Q_TOA_JJA_65N_raw)
      IF (par%master) forcing%Q_TOA_JJA_65N = forcing%Q_TOA_JJA_65N + Q_TOA_JJA_65N_raw
    END DO
  
    ! Southern Hemisphere (December, January, February)
    forcing%Q_TOA_DJF_80S = 0._dp
   
    DO m = 1,12
      IF ((m == 1) .OR. (m == 2) .OR. (m == 12)) THEN
        CALL get_insolation_at_time_month_and_lat( time, m, -80._dp, Q_TOA_DJF_80S_raw)
        IF (par%master) forcing%Q_TOA_DJF_80S = forcing%Q_TOA_DJF_80S + Q_TOA_DJF_80S_raw
      END IF
    END DO 

    IF (par%master) THEN
      forcing%Q_TOA_DJF_80S = forcing%Q_TOA_DJF_80S/3._dp
      forcing%Q_TOA_JJA_65N = forcing%Q_TOA_JJA_65N/3._dp
    END IF
    CALL sync
   
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE
   
  SUBROUTINE update_insolation_timeframes_from_file( time)
    ! Read the NetCDF file containing the insolation forcing data. Only read the time frames enveloping the current
    ! coupling timestep to save on memory usage. Only done by master.

    ! NOTE: assumes time in forcing file is in kyr

    IMPLICIT NONE

    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_insolation_timeframes_from_file'
    INTEGER, POINTER                                   :: ti0, ti1
    INTEGER                                            :: wti0, wti1
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate ti0 and ti1
    CALL allocate_shared_int_0D( ti0, wti0)
    CALL allocate_shared_int_0D( ti1, wti1)
    
    IF     (C%choice_insolation_forcing == 'none') THEN
      CALL crash('insolation should not be used when choice_insolation_forcing = "none"!')
    ELSEIF (C%choice_insolation_forcing == 'static' .OR. &
            C%choice_insolation_forcing == 'realistic') THEN
      ! Update insolation

      ! Check if data for model time is available
      IF (time < forcing%ins_time(1)) THEN
        CALL crash('queried time before start of insolation record!')
      END IF

      ! Find time indices to be read
      IF (par%master) THEN
       IF (time <= forcing%ins_time( forcing%ins_nyears)) THEN
          ti1 = 1
          DO WHILE (forcing%ins_time(ti1) < time)
            ti1 = ti1 + 1
          END DO
          ti0 = ti1 - 1
 
          forcing%ins_t0 = forcing%ins_time(ti0)
          forcing%ins_t1 = forcing%ins_time(ti1)
        ELSE
          ! Constant PD insolation for future projections
          ti0 = forcing%ins_nyears
          ti1 = forcing%ins_nyears
 
          forcing%ins_t0 = forcing%ins_time(ti0) - 1._dp
          forcing%ins_t1 = forcing%ins_time(ti1)
       END IF
      END IF
      CALL sync

      ! Read new insolation fields from the NetCDF file
      CALL read_insolation_data_file_timeframes( forcing, ti0, ti1, forcing%ins_Q_TOA0, forcing%ins_Q_TOA1)
      CALL sync

    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Clean up after yourself 
    CALL deallocate_shared(wti0)
    CALL deallocate_shared(wti1)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_insolation_timeframes_from_file
  SUBROUTINE initialise_insolation_data
    ! Allocate shared memory for the forcing data fields

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_insolation_data'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_insolation_forcing == 'none') THEN
      ! No insolation included, likely because we're running an idealised-geometry experiment
    ELSEIF (C%choice_insolation_forcing == 'static' .OR. &
            C%choice_insolation_forcing == 'realistic') THEN
      ! Initialise insolation

      ! The times at which we have insolation fields from Laskar, between which we'll interpolate
      ! to find the insolation at model time (assuming that t0 <= model_time <= t1)

      CALL allocate_shared_dp_0D( forcing%ins_t0, forcing%wins_t0)
      CALL allocate_shared_dp_0D( forcing%ins_t1, forcing%wins_t1)

      IF (par%master) THEN
        ! Give impossible values to timeframes, so that the first call to get_insolation_at_time
        ! is guaranteed to first read two new timeframes from the NetCDF file
        forcing%ins_t0 = C%start_time_of_run - 100._dp
        forcing%ins_t1 = C%start_time_of_run - 90._dp
      END IF ! IF (par%master) THEN
      CALL sync

      IF (par%master) WRITE(0,*) '  Initialising insolation data from ', TRIM(C%filename_insolation), '...'

      ! Inquire into the insolation forcing netcdf file
      CALL allocate_shared_int_0D( forcing%ins_nyears, forcing%wins_nyears)
      CALL allocate_shared_int_0D( forcing%ins_nlat,   forcing%wins_nlat  )

      forcing%netcdf_ins%filename = C%filename_insolation

      CALL inquire_insolation_data_file( forcing)
      CALL sync

      ! Insolation
      CALL allocate_shared_dp_1D( forcing%ins_nyears,   forcing%ins_time,      forcing%wins_time     )
      CALL allocate_shared_dp_1D( forcing%ins_nlat,     forcing%ins_lat,       forcing%wins_lat      )
      CALL allocate_shared_dp_2D( forcing%ins_nlat, 12, forcing%ins_Q_TOA0,    forcing%wins_Q_TOA0   )
      CALL allocate_shared_dp_2D( forcing%ins_nlat, 12, forcing%ins_Q_TOA1,    forcing%wins_Q_TOA1   )
      CALL allocate_shared_dp_0D( forcing%Q_TOA_DJF_80S, forcing%wQ_TOA_DJF_80S)
      CALL allocate_shared_dp_0D( forcing%Q_TOA_JJA_65N, forcing%wQ_TOA_JJA_65N)

      ! Read time and latitude data
      CALL read_insolation_data_file_time_lat( forcing)

      IF (par%master) THEN
        IF (C%start_time_of_run < forcing%ins_time(1)) THEN
          CALL warning(' Model time starts before start of insolation record; the model will crash lol')
        END IF
        IF (C%end_time_of_run > forcing%ins_time(forcing%ins_nyears)) THEN
          CALL warning(' Model time will reach beyond end of insolation record; constant extrapolation will be used in that case!')
        END IF
      END IF
      CALL sync

    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_insolation_data

! == Geothermal heat flux

  SUBROUTINE initialise_geothermal_heat_flux_regional( grid, ice, region_name)
    ! Calculate the flow factor A in Glen's flow law

    USE parameters_module, ONLY: sec_per_year

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_geothermal_heat_flux_regional'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_geothermal_heat_flux == 'constant') THEN
      ice%GHF_a( :,grid%i1:grid%i2) = C%constant_geothermal_heat_flux
      CALL sync
    ELSEIF (C%choice_geothermal_heat_flux == 'spatial') THEN
      IF (par%master) WRITE(0,*) '  Initialising geothermal heat flux data from ', TRIM(C%filename_geothermal_heat_flux), '...'
      
      CALL read_field_from_file_2D(   C%filename_geothermal_heat_flux, 'hflux', grid, ice%GHF_a,  region_name)
      ice%GHF_a( :,grid%i1:grid%i2)  = ice%GHF_a( :,grid%i1:grid%i2)  * sec_per_year
    ELSE
      CALL crash('unknown choice_geothermal_heat_flux "' // TRIM( C%choice_geothermal_heat_flux) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_geothermal_heat_flux_regional

! ===== Sea level records =====
! =============================

  ! Read in a sea level record from a file
  SUBROUTINE initialise_sealevel_record
    ! Read the sea level record specified in C%filename_sealevel_record. Assumes
    ! this is an ASCII text file with at least two columns (time in kyr and sea level in m)
    ! and the number of rows being equal to C%sealevel_record_length
    ! NOTE: assumes time is listed in YEARS (so LGM would be -21000.0)

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_sealevel_record'
    INTEGER                                            :: i,ios

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF     (C%choice_sealevel_model == 'prescribed') THEN
      ! Observed sea level is needed for these methods.
    ELSE
      CALL crash('should only be called when choice_sealevel_method == "prescribed"!')
    END IF

    ! Allocate shared memory to take the data
    CALL allocate_shared_dp_1D( C%sealevel_record_length, forcing%sealevel_time,   forcing%wsealevel_time  )
    CALL allocate_shared_dp_1D( C%sealevel_record_length, forcing%sealevel_record, forcing%wsealevel_record)
    CALL allocate_shared_dp_0D(                           forcing%sealevel_obs,    forcing%wsealevel_obs   )

    ! Read sea level record (time and values) from specified text file
    IF (par%master) THEN

        WRITE(0,*) ' Reading sea level record from ', TRIM(C%filename_sealevel_record), '...'

        OPEN(   UNIT = 1337, FILE=C%filename_sealevel_record, ACTION='READ')
        DO i = 1, C%sealevel_record_length
          READ( UNIT = 1337, FMT=*, IOSTAT=ios) forcing%sealevel_time(i), forcing%sealevel_record(i)
          IF (ios /= 0) THEN
            CALL crash('length of text file "' // TRIM(C%filename_sealevel_record) // '" does not match C%sealevel_record_length!')
          END IF
        END DO
        CLOSE( UNIT  = 1337)

        IF (C%start_time_of_run < forcing%sealevel_time(1)) THEN
          CALL warning(' Model time starts before start of sea level record; constant extrapolation will be used in that case!')
        END IF
        IF (C%end_time_of_run > forcing%sealevel_time(C%sealevel_record_length)) THEN
          CALL warning(' Model time will reach beyond end of sea level record; constant extrapolation will be used in that case!')
        END IF

    END IF ! IF (par%master) THEN
    CALL sync

    ! Set the value for the current (starting) model time
    CALL update_sealevel_record_at_model_time( C%start_time_of_run)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_sealevel_record

  SUBROUTINE update_sealevel_record_at_model_time( time)
    ! Interpolate the data in forcing%sealevel to find the value at the queried time.
    ! If time lies outside the range of forcing%sealevel_time, return the first/last value
    !
    ! NOTE: assumes time is listed in YEARS (so LGM would be -21000.0)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_sealevel_record_at_model_time'
    INTEGER                                            :: il, iu
    REAL(dp)                                           :: wl, wu

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF     (C%choice_sealevel_model == 'prescribed') THEN
      ! Observed sea level is needed for these forcing methods.
    ELSE
      CALL crash('should only be called when choice_sealevel_model = "prescribed"!')
    END IF

    IF (par%master) THEN
      IF     (time < MINVAL( forcing%sealevel_time)) THEN ! Remember: forcing%sealevel_time is in YEARS, not kyr
        ! Model time before start of sea level record; using constant extrapolation
        forcing%sealevel_obs = forcing%sealevel_record( 1)
      ELSEIF (time > MAXVAL( forcing%sealevel_time)) THEN
        ! Model time beyond end of sea level record; using constant extrapolation
        forcing%sealevel_obs = forcing%sealevel_record( C%sealevel_record_length)
      ELSE
        iu = 1
        DO WHILE (forcing%sealevel_time(iu) < time)
          iu = iu+1
        END DO
        il = iu-1

        wl = (forcing%sealevel_time(iu) - time) / ((forcing%sealevel_time(iu)-forcing%sealevel_time(il)))
        wu = 1._dp - wl

        forcing%sealevel_obs = forcing%sealevel_record(il) * wl + forcing%sealevel_record(iu) * wu

      END IF
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_sealevel_record_at_model_time

END MODULE forcing_module
