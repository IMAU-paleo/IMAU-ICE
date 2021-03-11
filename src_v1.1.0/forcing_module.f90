MODULE forcing_module
 
  USE mpi
  USE configuration_module,        ONLY: dp, C
  USE parallel_module,             ONLY: par, sync, &
                                         allocate_shared_int_0D, allocate_shared_dp_0D, &
                                         allocate_shared_int_1D, allocate_shared_dp_1D, &
                                         allocate_shared_int_2D, allocate_shared_dp_2D, &
                                         allocate_shared_int_3D, allocate_shared_dp_3D, &
                                         deallocate_shared
  USE data_types_module,           ONLY: type_forcing_data, type_model_region, type_grid
  USE netcdf_module,               ONLY: inquire_insolation_data_file, read_insolation_data_file_time_lat, read_insolation_data_file, &
                                         read_inverse_routine_history_dT_glob, read_inverse_routine_history_dT_glob_inverse, read_inverse_routine_history_CO2_inverse, &
                                         inquire_geothermal_heat_flux_file, read_geothermal_heat_flux_file

  IMPLICIT NONE
  
  ! The data structure containing model forcing data - CO2 record, d18O record, (global) insolation record
  ! Updated at every coupling time step. 
  TYPE(type_forcing_data), SAVE :: forcing
    
CONTAINS
  
  ! =======================
  ! Inverse forward routine
  ! =======================
  
  SUBROUTINE inverse_routine_global_temperature_offset
    ! Use the inverse routine to calculate a global temperature offset
    ! (i.e. the method used in de Boer et al., 2013)
    ! (de Boer, B., van de Wal, R., Lourens, L. J., Bintanja, R., and Reerink, T. J.:
    ! A continuous simulation of global ice volume over the past 1 million years with 3-D ice-sheet models, Climate Dynamics 41, 1365-1384, 2013)
    
    ! Local variables:
    INTEGER                                            :: cerr,ierr
    REAL(dp)                                           :: dT_glob_inverse_average_over_window
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in inverse_routine_global_temperature_offset!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! The inverse routine might not work properly when not all ice sheets are simulated
    IF ((.NOT. C%do_NAM) .OR. (.NOT. C%do_EAS) .OR. (.NOT. C%do_GRL) .OR. (.NOT. C%do_ANT)) THEN
      IF (par%master) WRITE(0,*) '  WARNING: The inverse routine only works properly when all four ice sheets are simulated!'
      IF (par%master) WRITE(0,*) '           Leaving one out means you will miss that contribution to the d18O, which the'
      IF (par%master) WRITE(0,*) '           routine will try to compensate for by making the world colder.'
      IF (par%master) WRITE(0,*) '           If you really want to simulate only one ice sheet, consider using direct CO2 forcing.'
    END IF
    
    IF (par%master) THEN
    
      ! Average dT_glob_inverse over the moving time window
      dT_glob_inverse_average_over_window = SUM( forcing%dT_glob_inverse_history) / REAL(forcing%ndT_glob_inverse_history,dp)
      
      ! Update dT_glob_inverse based on the difference between modelled and observed d18O
      forcing%dT_glob_inverse = dT_glob_inverse_average_over_window + (forcing%d18O_mod - forcing%d18O_obs) * C%inverse_d18O_to_dT_glob_scaling
      
      ! Update the moving time window
      forcing%dT_glob_inverse_history( 2:forcing%ndT_glob_inverse_history) = forcing%dT_glob_inverse_history( 1:forcing%ndT_glob_inverse_history-1)
      forcing%dT_glob_inverse_history( 1) = forcing%dT_glob_inverse
      
      !WRITE(0,*) ' dT_glob_inverse_history = ', forcing%dT_glob_inverse_history
      
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE inverse_routine_global_temperature_offset
  SUBROUTINE inverse_routine_CO2
    ! Use the inverse routine to calculate modelled CO2
    ! (i.e. the method used in Berends et al., 2019)
    ! (Berends, C. J., de Boer, B., Dolan, A. M., Hill, D. J., and van de Wal, R. S. W.: Modelling ice sheet evolution and atmospheric CO2 during the Late Pliocene, Climate of the Past 15, 1603-1619, 2019)
    
    ! Local variables:
    INTEGER                                            :: cerr,ierr
    REAL(dp)                                           :: CO2_inverse_average_over_window
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in inverse_routine_global_temperature_offset!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    IF (par%master) THEN
    
      ! Average CO2_inverse over the moving time window
      CO2_inverse_average_over_window = SUM( forcing%CO2_inverse_history) / REAL(forcing%nCO2_inverse_history,dp)
      
      ! Update CO2_inverse based on the difference between modelled and observed d18O
      forcing%CO2_inverse = CO2_inverse_average_over_window + (forcing%d18O_mod - forcing%d18O_obs) * C%inverse_d18O_to_CO2_scaling
      
      ! Update the moving time window
      forcing%CO2_inverse_history( 2:forcing%nCO2_inverse_history) = forcing%CO2_inverse_history( 1:forcing%nCO2_inverse_history-1)
      forcing%CO2_inverse_history( 1) = forcing%CO2_inverse
      
      forcing%CO2_mod = forcing%CO2_inverse
      
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE inverse_routine_CO2
  SUBROUTINE calculate_modelled_d18O( NAM, EAS, GRL, ANT)
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_model_region),             INTENT(IN)    :: NAM, EAS, GRL, ANT
    
    ! Local variables
    INTEGER                                            :: cerr,ierr
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in calculate_modelled_d18O!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
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
    INTEGER                                            :: cerr,ierr
    REAL(dp)                                           :: dT_NAM, dT_EAS, dT_GRL, dT_ANT, dT_glob_average_over_window
    INTEGER                                            :: n_glob
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in update_global_mean_temperature_change_history!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Determine annual mean surface temperature change for all model regions
    dT_NAM = 0._dp
    dT_EAS = 0._dp
    dT_GRL = 0._dp
    dT_ANT = 0._dp
    IF (C%do_NAM) CALL calculate_mean_temperature_change_region( NAM, dT_NAM)
    IF (C%do_EAS) CALL calculate_mean_temperature_change_region( EAS, dT_EAS)
    IF (C%do_GRL) CALL calculate_mean_temperature_change_region( GRL, dT_GRL)
    IF (C%do_ANT) CALL calculate_mean_temperature_change_region( ANT, dT_ANT)
    
    IF (par%master) THEN
    
      ! Weighted average of all model regions
      forcing%dT_glob = 0._dp
      n_glob          = 0
      
      IF (C%do_NAM) THEN
        forcing%dT_glob = forcing%dT_glob + (dT_NAM * NAM%grid%nx * NAM%grid%ny)
        n_glob  = n_glob  + NAM%grid%nx * NAM%grid%ny
      END IF
      IF (C%do_EAS) THEN
        forcing%dT_glob = forcing%dT_glob + (dT_EAS * EAS%grid%nx * EAS%grid%ny)
        n_glob  = n_glob  + EAS%grid%nx * EAS%grid%ny
      END IF
      IF (C%do_GRL) THEN
        forcing%dT_glob = forcing%dT_glob + (dT_GRL * GRL%grid%nx * GRL%grid%ny)
        n_glob  = n_glob  + GRL%grid%nx * GRL%grid%ny
      END IF
      IF (C%do_ANT) THEN
        forcing%dT_glob = forcing%dT_glob + (dT_ANT * ANT%grid%nx * ANT%grid%ny)
        n_glob  = n_glob  + ANT%grid%nx * ANT%grid%ny
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
    
  END SUBROUTINE update_global_mean_temperature_change_history
  SUBROUTINE calculate_mean_temperature_change_region( region, dT)
    ! Calculate the annual mean surface temperature change w.r.t PD (corrected for elevation changes) over a model region
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_model_region),             INTENT(IN)    :: region
    REAL(dp),                            INTENT(OUT)   :: dT
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: dT_lapse_mod, dT_lapse_PD, T_pot_mod, T_pot_PD
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in calculate_mean_temperature_change_region!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    dT = 0._dp
    
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
      dT_lapse_mod = region%ice%Hs_Aa(         j,i) * C%constant_lapserate
      dT_lapse_PD  = region%climate%PD_obs%Hs( j,i) * C%constant_lapserate
      DO m = 1, 12
        T_pot_mod = region%climate%applied%T2m( m,j,i) - dT_lapse_mod
        T_pot_PD  = region%climate%PD_obs%T2m(  m,j,i) - dT_lapse_PD
        dT = dT + T_pot_mod - T_pot_PD
      END DO
    END DO
    END DO
    
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dT, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    dT = dT / (region%grid%nx * region%grid%ny * 12._dp)
    
  END SUBROUTINE calculate_mean_temperature_change_region
  SUBROUTINE initialise_d18O_data
    ! Allocate shared memory for the d18O and global temperature variables.
    
    IMPLICIT NONE
    
    ! Determine number of entries in the global mean temperature change history
    CALL allocate_shared_int_0D( forcing%ndT_glob_history, forcing%wndT_glob_history)
    IF (par%master) forcing%ndT_glob_history = CEILING( C%dT_deepwater_averaging_window / C%dt_coupling)
    CALL sync
    ! Allocate memory for the global mean temperature change history
    CALL allocate_shared_dp_1D( forcing%ndT_glob_history, forcing%dT_glob_history, forcing%wdT_glob_history)
    
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
    CALL allocate_shared_dp_0D( forcing%CO2_obs,                   forcing%wCO2_obs                  )
    CALL allocate_shared_dp_0D( forcing%CO2_mod,                   forcing%wCO2_mod                  )
    
    forcing%d18O_obs_PD = 3.23_dp
    
  END SUBROUTINE initialise_d18O_data
  SUBROUTINE initialise_inverse_routine_data
    ! Allocate shared memory for the moving time windows used in the inverse routine
    
    IMPLICIT NONE
    
    ! Local variables
    INTEGER                                            :: cerr,ierr
    CHARACTER(LEN=256)                                 :: filename
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_inverse_routine_data!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Not everything is needed for all forcing methods
    IF (C%choice_forcing_method == 'CO2_direct') THEN
    
      IF (par%master) forcing%CO2_mod = forcing%CO2_obs
      
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
    
      CALL allocate_shared_dp_0D( forcing%dT_glob_inverse, forcing%wdT_glob_inverse)
      ! Determine number of entries in the history
      CALL allocate_shared_int_0D( forcing%ndT_glob_inverse_history, forcing%wndT_glob_inverse_history)
      IF (par%master) forcing%ndT_glob_inverse_history = CEILING( C%dT_glob_inverse_averaging_window / C%dt_coupling)
      CALL sync
      ! Allocate memory for the global mean temperature change history
      CALL allocate_shared_dp_1D( forcing%ndT_glob_inverse_history, forcing%dT_glob_inverse_history, forcing%wdT_glob_inverse_history)
      IF (par%master) forcing%dT_glob_inverse_history = 0._dp
      IF (par%master) forcing%dT_glob_inverse         = 0._dp
      IF (par%master) forcing%CO2_mod                 = 0._dp
      
      ! If we're restarting a previous run, read inverse routine history from one of the restart files
      IF (C%is_restart) THEN
        IF (C%do_NAM) THEN
          filename = C%filename_init_NAM
        ELSEIF (C%do_EAS) THEN
          filename = C%filename_init_EAS
        ELSEIF (C%do_GRL) THEN
          filename = C%filename_init_GRL
        ELSEIF (C%do_ANT) THEN
          filename = C%filename_init_ANT
        END IF
        IF (par%master) CALL read_inverse_routine_history_dT_glob(         forcing, C%filename_init_NAM)
        IF (par%master) CALL read_inverse_routine_history_dT_glob_inverse( forcing, C%filename_init_NAM)
        IF (par%master) forcing%dT_glob_inverse = forcing%dT_glob_inverse_history(1)
        CALL sync
      END IF
      
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
    
      CALL allocate_shared_dp_0D( forcing%CO2_inverse, forcing%wCO2_inverse)
      ! Determine number of entries in the history
      CALL allocate_shared_int_0D( forcing%nCO2_inverse_history, forcing%wnCO2_inverse_history)
      IF (par%master) forcing%nCO2_inverse_history = CEILING( C%CO2_inverse_averaging_window / C%dt_coupling)
      CALL sync
      ! Allocate memory for the global mean temperature change history
      CALL allocate_shared_dp_1D( forcing%nCO2_inverse_history, forcing%CO2_inverse_history, forcing%wCO2_inverse_history)
      IF (par%master) forcing%CO2_inverse_history = C%inverse_d18O_to_CO2_initial_CO2
      IF (par%master) forcing%CO2_inverse         = C%inverse_d18O_to_CO2_initial_CO2
      IF (par%master) forcing%CO2_mod             = C%inverse_d18O_to_CO2_initial_CO2
      
      ! If we're restarting a previous run, read inverse routine history from one of the restart files
      IF (C%is_restart) THEN
        IF (C%do_NAM) THEN
          filename = C%filename_init_NAM
        ELSEIF (C%do_EAS) THEN
          filename = C%filename_init_EAS
        ELSEIF (C%do_GRL) THEN
          filename = C%filename_init_GRL
        ELSEIF (C%do_ANT) THEN
          filename = C%filename_init_ANT
        END IF
        IF (par%master) CALL read_inverse_routine_history_dT_glob(     forcing, C%filename_init_NAM)
        IF (par%master) CALL read_inverse_routine_history_CO2_inverse( forcing, C%filename_init_NAM)
        IF (par%master) forcing%CO2_inverse = forcing%CO2_inverse_history(1)
        IF (par%master) forcing%CO2_mod     = forcing%CO2_inverse
        CALL sync
      END IF
      
    ELSE
      WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in initialise_inverse_routine_data!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE initialise_inverse_routine_data

  ! =========================================
  ! Read and update forcing data from records
  ! =========================================

  ! CO2
  SUBROUTINE update_CO2_at_model_time( time)
    ! Interpolate the data in forcing%CO2 to find the value at the queried time.
    ! If time lies outside the range of forcing%CO2_time, return the first/last value
    ! NOTE: assumes time is listed in kyr (so LGM would be -21000.0)
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables
    INTEGER                                            :: cerr,ierr
    INTEGER                                            :: il, iu
    REAL(dp)                                           :: wl, wu
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in update_CO2_at_model_time!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Not needed for all forcing methods
    IF (C%choice_forcing_method == 'CO2_direct') THEN
      ! Observed CO2 is needed for these forcing methods.
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob' .OR. &
            C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! Observed CO2 is not needed for these forcing methods
      RETURN
    ELSE
      WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in update_CO2_at_model_time!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    IF (par%master) THEN
      IF     (time < MINVAL( forcing%CO2_time) * 1000._dp) THEN ! times 1000 because forcing%CO2_time is in kyr
        IF (par%master) WRITE(0,*) '  WARNING: model time before start of CO2 record, using constant extrapolation!'
        forcing%CO2_obs = forcing%CO2_record( 1)
      ELSEIF (time > MAXVAL( forcing%CO2_time) * 1000._dp) THEN
        IF (par%master) WRITE(0,*) '  WARNING: model time beyond end of CO2 record, using constant extrapolation!'
        forcing%CO2_obs = forcing%CO2_record( C%CO2_record_length)
      ELSE
        iu = 1
        DO WHILE (forcing%CO2_time(iu) * 1000._dp < time)
          iu = iu+1
        END DO
        il = iu-1
        
        wl = (forcing%CO2_time(iu)*1000._dp - time) / ((forcing%CO2_time(iu)-forcing%CO2_time(il))*1000._dp)
        wu = (time - forcing%CO2_time(il)*1000._dp) / ((forcing%CO2_time(iu)-forcing%CO2_time(il))*1000._dp)
        
        forcing%CO2_obs = forcing%CO2_record(il) * wl + forcing%CO2_record(iu) * wu
      END IF
    END IF
    CALL sync
    
  END SUBROUTINE update_CO2_at_model_time
  SUBROUTINE initialise_CO2_record
    ! Read the CO2 record specified in C%filename_CO2_record. Assumes this is an ASCII text file with at least two columns (time in kyr and CO2 in ppmv)
    ! and the number of rows being equal to C%CO2_record_length
    ! NOTE: assumes time is listed in kyr (so LGM would be -21.0)
    
    IMPLICIT NONE
    
    ! Local variables
    INTEGER                                            :: i,ios,cerr,ierr
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_CO2_record!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Not needed for all forcing methods
    IF (C%choice_forcing_method == 'CO2_direct') THEN
      ! Observed CO2 is needed for these forcing methods.
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob' .OR. &
            C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! Observed CO2 is not needed for these forcing methods
      RETURN
    ELSE
      WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in initialise_CO2_record!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Reading CO2 record from ', TRIM(C%filename_CO2_record), '...'
    
    ! Allocate shared memory to take the data
    CALL allocate_shared_dp_1D( C%CO2_record_length, forcing%CO2_time,   forcing%wCO2_time  )
    CALL allocate_shared_dp_1D( C%CO2_record_length, forcing%CO2_record, forcing%wCO2_record)
    
    ! Read CO2 record (time and values) from specified text file
    IF (par%master) THEN
      
      OPEN(   UNIT = 1337, FILE=C%filename_CO2_record, ACTION='READ')
      DO i = 1, C%CO2_record_length
        READ( UNIT = 1337, FMT=*, IOSTAT=ios) forcing%CO2_time(i), forcing%CO2_record(i) 
        IF (ios /= 0) THEN
          WRITE(0,*) ' read_CO2_record - ERROR: length of text file "', TRIM(C%filename_CO2_record), '" does not match C%CO2_record_length = ', C%CO2_record_length
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END DO
    
      CLOSE( UNIT  = 1337)
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Set the value for the current (starting) model time
    CALL update_CO2_at_model_time( C%start_time_of_run)
    
  END SUBROUTINE initialise_CO2_record
  
  ! d18O
  SUBROUTINE update_d18O_at_model_time( time)
    ! Interpolate the data in forcing%d18O to find the value at the queried time.
    ! If time lies outside the range of forcing%d18O_time, return the first/last value
    ! NOTE: assumes time is listed in years! (so LGM would be -21000.0)
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables
    INTEGER                                            :: cerr,ierr
    INTEGER                                            :: il, iu
    REAL(dp)                                           :: wl, wu
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in update_d18O_at_model_time!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Not needed for all forcing methods
    IF (C%choice_forcing_method == 'd18O_inverse_dT_glob' .OR. &
        C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! Observed d18O is needed for these forcing methods.
    ELSEIF (C%choice_forcing_method == 'CO2_direct') THEN
      ! Observed d18O is not needed for these forcing methods
      RETURN
    ELSE
      WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in update_d18O_at_model_time!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    IF (par%master) THEN
  
      IF     (time < MINVAL( forcing%d18O_time)) THEN ! times 1000 because forcing%d18O_time is in kyr
        IF (par%master) WRITE(0,*) '  WARNING: model time before start of d18O record, using constant extrapolation!'
        forcing%d18O_obs = forcing%d18O_record( 1)
      ELSEIF (time > MAXVAL( forcing%d18O_time)) THEN
        IF (par%master) WRITE(0,*) '  WARNING: model time beyond end of d18O record, using constant extrapolation!'
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
    
  END SUBROUTINE update_d18O_at_model_time
  SUBROUTINE initialise_d18O_record
    ! Read the d18O record specified in C%filename_d18O_record. Assumes this is an ASCII text file with at least two columns (time in yr and d18O in per mil)
    ! and the number of rows being equal to C%d18O_record_length
    ! NOTE: assumes time is listed in years! (so LGM would be -21000.0)
    
    IMPLICIT NONE
    
    ! Local variables
    INTEGER                                            :: i,ios,cerr,ierr
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_d18O_record!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Not needed for all forcing methods
    IF (C%choice_forcing_method == 'd18O_inverse_dT_glob' .OR. &
        C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! Observed d18O is needed for these forcing methods.
    ELSEIF (C%choice_forcing_method == 'CO2_direct') THEN
      ! Observed d18O is not needed for these forcing methods
      RETURN
    ELSE
      WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in initialise_d18O_record!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    IF (par%master) WRITE(0,*) ''
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
          WRITE(0,*) ' read_d18O_record - ERROR: length of text file "', TRIM(C%filename_d18O_record), '" does not match C%d18O_record_length = ', C%d18O_record_length
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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
    
  END SUBROUTINE initialise_d18O_record

  ! Insolation
  SUBROUTINE update_insolation_data( t_coupling)
    ! Read the NetCDF file containing the insolation forcing data. Only read the time frames enveloping the current
    ! coupling timestep to save on memory usage. Only done by master.
    
    ! NOTE: assumes time in forcing file is in kyr
    
    IMPLICIT NONE

    REAL(dp),                            INTENT(IN)    :: t_coupling
    
    ! Local variables
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: ti0, ti1
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in update_insolation_data!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Updating insolation data...'
    
    ! Initialise at zero
    IF (par%master) THEN
      forcing%ins_Q_TOA0 = 0._dp
      forcing%ins_Q_TOA1 = 0._dp
    END IF
    
    ! Check if data for model time is available
    IF (t_coupling < forcing%ins_time(1)) THEN
      WRITE(0,*) '  update_insolation_data - ERROR: insolation data only available between ', MINVAL(forcing%ins_time), ' y and ', MAXVAL(forcing%ins_time), ' y'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Find time indices to be read
    IF (par%master) THEN
      IF (t_coupling <= forcing%ins_time( forcing%ins_nyears)) THEN
        ti1 = 1
        DO WHILE (forcing%ins_time(ti1) < t_coupling)
          ti1 = ti1 + 1
        END DO
        ti0 = ti1 - 1
        
        forcing%ins_t0 = forcing%ins_time(ti0)
        forcing%ins_t1 = forcing%ins_time(ti1)
      ELSE
        IF (par%master) WRITE(0,*) '  WARNING: using constant PD insolation for future projections!'
        ti0 = forcing%ins_nyears
        ti1 = forcing%ins_nyears
        
        forcing%ins_t0 = forcing%ins_time(ti0) - 1._dp
        forcing%ins_t1 = forcing%ins_time(ti1)
      END IF
    END IF ! IF (par%master) THEN
        
    ! Read new insolation fields from the NetCDF file
    IF (par%master) CALL read_insolation_data_file( forcing, ti0, ti1, forcing%ins_Q_TOA0, forcing%ins_Q_TOA1)
    CALL sync
    
  END SUBROUTINE update_insolation_data
  SUBROUTINE map_insolation_to_grid( grid, ins_t0, ins_t1, Q_TOA0, Q_TOA1, time, Q_TOA)
    ! Interpolate two insolation timeframes to the desired time, and then map it to the model grid.
      
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp),                            INTENT(IN)    :: ins_t0, ins_t1
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Q_TOA0, Q_TOA1
    REAL(dp),                            INTENT(IN)    :: time
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: Q_TOA
    
    ! Local variables:
    INTEGER                                            :: i,j,m, ilat_l, ilat_u
    REAL(dp)                                           :: wt0, wt1, wlat_l, wlat_u
    
    ! Calculate time interpolation weights
    wt0 = (ins_t1 - time) / (ins_t1 - ins_t0)
    wt1 = 1._dp - wt0
        
    ! Interpolate on the grid
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
     
      ilat_l = FLOOR(grid%lat(j,i) + 91)
      ilat_u = ilat_l + 1
      
      wlat_l = forcing%ins_lat(ilat_u) - grid%lat(j,i)
      wlat_u = 1._dp - wlat_l
      
      DO m = 1, 12
        Q_TOA( m,j,i) = (wt0 * wlat_l * Q_TOA0( ilat_l,m)) + &
                        (wt0 * wlat_u * Q_TOA0( ilat_u,m)) + &
                        (wt1 * wlat_l * Q_TOA1( ilat_l,m)) + &
                        (wt1 * wlat_u * Q_TOA1( ilat_u,m))
      END DO    
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_insolation_to_grid
  SUBROUTINE initialise_insolation_data
    ! Allocate shared memory for the forcing data fields
    
    IMPLICIT NONE
    
    INTEGER :: cerr, ierr
    
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Initialising insolation data from ', TRIM(C%filename_insolation), '...'
        
    ! The times at which we have insolation fields from Laskar, between which we'll interpolate
    ! to find the insolation at model time (ins_t0 < model_time < ins_t1)
    
    CALL allocate_shared_dp_0D( forcing%ins_t0, forcing%wins_t0)
    CALL allocate_shared_dp_0D( forcing%ins_t1, forcing%wins_t1)
    
    IF (par%master) THEN
      forcing%ins_t0 = C%start_time_of_run
      forcing%ins_t1 = C%end_time_of_run
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_insolation_data!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Inquire into the insolation forcing netcdf file    
    CALL allocate_shared_int_0D( forcing%ins_nyears, forcing%wins_nyears)
    CALL allocate_shared_int_0D( forcing%ins_nlat,   forcing%wins_nlat  )
    
    forcing%netcdf%filename = C%filename_insolation
    
    IF (par%master) CALL inquire_insolation_data_file( forcing)
    CALL sync
    
    ! Insolation    
    CALL allocate_shared_dp_1D( forcing%ins_nyears,   forcing%ins_time,    forcing%wins_time   )
    CALL allocate_shared_dp_1D( forcing%ins_nlat,     forcing%ins_lat,     forcing%wins_lat    )
    CALL allocate_shared_dp_2D( forcing%ins_nlat, 12, forcing%ins_Q_TOA0,  forcing%wins_Q_TOA0 )
    CALL allocate_shared_dp_2D( forcing%ins_nlat, 12, forcing%ins_Q_TOA1,  forcing%wins_Q_TOA1 )
    
    ! Read time and latitude data
    IF (par%master) CALL read_insolation_data_file_time_lat( forcing)
    CALL sync
    
    ! Read insolation data
    CALL update_insolation_data( C%start_time_of_run)
    
  END SUBROUTINE initialise_insolation_data
  SUBROUTINE initialise_geothermal_heat_flux

    IMPLICIT NONE

    INTEGER :: cerr, ierr

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Initialising geothermal heat flux data from ', TRIM(C%filename_geothermal_heat_flux), '...'

    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_geothermal_heat_flux!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN

    ! Inquire into the insolation forcing netcdf file
    CALL allocate_shared_int_0D( forcing%ghf_nlat,   forcing%wghf_nlat  )
    CALL allocate_shared_int_0D( forcing%ghf_nlon,   forcing%wghf_nlon  )

    forcing%netcdf_ghf%filename = C%filename_geothermal_heat_flux

    IF (par%master) CALL inquire_geothermal_heat_flux_file( forcing)
    CALL sync
  
    ! Geothermal heat flux
    CALL allocate_shared_dp_1D( forcing%ghf_nlon,     forcing%ghf_lon,     forcing%wghf_lon    )
    CALL allocate_shared_dp_1D( forcing%ghf_nlat,     forcing%ghf_lat,     forcing%wghf_lat    )
    CALL allocate_shared_dp_2D( forcing%ghf_nlon,     forcing%ghf_nlat,    forcing%ghf_ghf,  forcing%wghf_ghf )

    IF (par%master) CALL read_geothermal_heat_flux_file( forcing, forcing%ghf_ghf)
    CALL sync

  END SUBROUTINE initialise_geothermal_heat_flux

END MODULE forcing_module
