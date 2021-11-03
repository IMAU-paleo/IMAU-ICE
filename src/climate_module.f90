MODULE climate_module

  ! Contains all the routines for calculating the climate forcing.

  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_climate_model, type_SMB_model, &
                                             type_climate_matrix, type_subclimate_global, type_subclimate_region
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, inquire_PD_obs_data_file, read_PD_obs_data_file, inquire_GCM_snapshot, read_GCM_snapshot, &
                                             read_insolation_data_file
  USE forcing_module,                  ONLY: forcing, map_insolation_to_grid
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             error_function, smooth_Gaussian_2D, smooth_Shepard_2D, &
                                             map_glob_to_grid_2D, map_glob_to_grid_3D, &
                                             map_square_to_square_cons_2nd_order_2D, map_square_to_square_cons_2nd_order_3D
  USE derivatives_and_grids_module,    ONLY: ddx_a_to_a_2D, ddy_a_to_a_2D
  USE SMB_module,                      ONLY: run_SMB_model
  USE ocean_module,                    ONLY: initialise_PD_obs_ocean_fields

  IMPLICIT NONE
    
CONTAINS

  ! Run the climate model on a region grid
  SUBROUTINE run_climate_model( grid, ice, SMB, climate, region_name, time)
    ! Run the regional climate model
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice 
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables
    REAL(dp)                                           :: wt0, wt1
    INTEGER                                            :: m
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  T2m_year
    INTEGER                                            :: wT2m_year
    
  ! ================================================
  ! ===== Exceptions for benchmark experiments =====
  ! ================================================
  
    IF (C%do_benchmark_experiment) THEN
      IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6') THEN
          
        CALL EISMINT_climate( grid, ice, climate, time)
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'Bueler' .OR. &
              C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
              C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_F' .OR. &
              C%choice_benchmark_experiment == 'MISMIPplus') THEN
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_climate_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  ! =======================================================
  ! ===== End of exceptions for benchmark experiments =====
  ! =======================================================
    
    ! Different kinds of climate forcing for realistic experiments
    IF (C%choice_forcing_method == 'SMB_direct') THEN
      ! Yearly SMB and surface temperature are prescribed as forcing

      ! Interpolate the two timeframes in time
      IF (par%master) THEN
        wt0 = (forcing%clim_t1 - time) / (forcing%clim_t1 - forcing%clim_t0)
        wt1 = 1._dp - wt0
        forcing%clim_SMB2      = (wt0 * forcing%clim_SMB0)      + (wt1 * forcing%clim_SMB1     )
        forcing%clim_T2m_year2 = (wt0 * forcing%clim_T2m_year0) + (wt1 * forcing%clim_T2m_year1)
      END IF
      CALL sync
      
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, T2m_year, wT2m_year)

      ! Map the interpolated timeframe to the model grid
      IF     (C%domain_climate_forcing == 'global') THEN     
        CALL map_glob_to_grid_2D(                    forcing%clim_nlat, forcing%clim_nlon, forcing%clim_lat, forcing%clim_lon, grid,                             forcing%clim_SMB2,      SMB%SMB_year)
        CALL map_glob_to_grid_2D(                    forcing%clim_nlat, forcing%clim_nlon, forcing%clim_lat, forcing%clim_lon, grid,                             forcing%clim_T2m_year2, T2m_year    )
      ELSEIF (C%domain_climate_forcing == 'regional') THEN
        CALL map_square_to_square_cons_2nd_order_2D( forcing%clim_nx,   forcing%clim_ny,   forcing%clim_x,   forcing%clim_y,   grid%nx, grid%ny, grid%x, grid%y, forcing%clim_SMB2,      SMB%SMB_year)
        CALL map_square_to_square_cons_2nd_order_2D( forcing%clim_nx,   forcing%clim_ny,   forcing%clim_x,   forcing%clim_y,   grid%nx, grid%ny, grid%x, grid%y, forcing%clim_T2m_year2, T2m_year    ) 
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: domain_climate_forcing "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_climate_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
      ! Set monthly temperatures equal to annual mean
      DO m = 1, 12
        climate%applied%T2m( m,:,grid%i1:grid%i2) = T2m_year( :,grid%i1:grid%i2)
      END DO
      
      ! Assume the prescribed climate is present-day (because we need a PD reference climate in the isotopes module)
      climate%PD_obs%T2m( :,:,grid%i1:grid%i2) = climate%applied%T2m( :,:,grid%i1:grid%i2)
      
      ! Clean up after yourself
      CALL deallocate_shared( wT2m_year)

    ELSEIF (C%choice_forcing_method == 'climate_direct') THEN
      ! Monthly climate fields are prescribed as forcing

      ! Interpolate the two timeframes in time
      IF (par%master) THEN
        wt0 = (forcing%clim_t1 - time) / (forcing%clim_t1 - forcing%clim_t0)
        wt1 = 1._dp - wt0
        forcing%clim_T2m2    = (wt0 * forcing%clim_T2m0)    + (wt1 * forcing%clim_T2m1   )   
        forcing%clim_Precip2 = (wt0 * forcing%clim_Precip0) + (wt1 * forcing%clim_Precip1)
      END IF
      CALL sync

      ! Map the interpolated timeframe to the model grid
      IF     (C%domain_climate_forcing == 'global') THEN
        CALL map_glob_to_grid_3D(                    forcing%clim_nlat, forcing%clim_nlon, forcing%clim_lat, forcing%clim_lon, grid,                             forcing%clim_T2m2,    climate%applied%T2m,    12)
        CALL map_glob_to_grid_3D(                    forcing%clim_nlat, forcing%clim_nlon, forcing%clim_lat, forcing%clim_lon, grid,                             forcing%clim_Precip2, climate%applied%Precip, 12) 
      ELSEIF (C%domain_climate_forcing == 'regional') THEN
        CALL map_square_to_square_cons_2nd_order_3D( forcing%clim_nx,   forcing%clim_ny,   forcing%clim_x,   forcing%clim_y,   grid%nx, grid%ny, grid%x, grid%y, forcing%clim_T2m2, climate%applied%T2m, 12)
        CALL map_square_to_square_cons_2nd_order_3D( forcing%clim_nx,   forcing%clim_ny,   forcing%clim_x,   forcing%clim_y,   grid%nx, grid%ny, grid%x, grid%y, forcing%clim_Precip2, climate%applied%Precip, 12)      
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: domain_climate_forcing "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_climate_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
      ! Assume the prescribed climate is present-day (because we need a PD reference climate in the isotopes module)
      climate%PD_obs%T2m(    :,:,grid%i1:grid%i2) = climate%applied%T2m(    :,:,grid%i1:grid%i2)
      climate%PD_obs%Precip( :,:,grid%i1:grid%i2) = climate%applied%Precip( :,:,grid%i1:grid%i2)
      CALL sync

    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      ! Use the global temperature offset as calculated by the inverse routine
      
      CALL run_climate_model_dT_glob( grid, ice, climate, region_name)
      
    ELSEIF (C%choice_forcing_method == 'CO2_direct' .OR. C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! Use CO2 (either taken directly from the specified record, or as calculated by the inverse routine) to force the climate matrix
      
      IF (C%choice_climate_matrix == 'warm_cold') THEN
        ! Use the two-snapshot climate matrix
        CALL run_climate_model_matrix_warm_cold( grid, ice, SMB, climate, region_name, time)
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_climate_matrix "', TRIM(C%choice_climate_matrix), '" not implemented in run_climate_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: forcing method "', TRIM(C%choice_forcing_method), '" not implemented in run_climate_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE run_climate_model
  
  ! Semi-parameterised climate (ERA40 + global temperature offset) from de Boer et al., 2013
  SUBROUTINE run_climate_model_dT_glob( grid, ice, climate, region_name)
    ! Use the climate parameterisation from de Boer et al., 2013 (global temperature offset calculated with the inverse routine,
    ! plus a precipitation correction based on temperature + orography changes (NAM & EAS; Roe & Lindzen model), or only temperature (GRL & ANT).
    ! (for more details, see de Boer, B., van de Wal, R., Lourens, L. J., Bintanja, R., and Reerink, T. J.:
    ! A continuous simulation of global ice volume over the past 1 million years with 3-D ice-sheet models, Climate Dynamics 41, 1365-1384, 2013)
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice 
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables:
    INTEGER                                            :: i,j,m
    
    REAL(dp), DIMENSION(:,:  ), POINTER                :: dHs_dx_ref, dHs_dy_ref
    REAL(dp)                                           :: dT_lapse
    REAL(dp), DIMENSION(:,:,:), POINTER                :: Precip_RL_ref, Precip_RL_mod, dPrecip_RL
    INTEGER                                            :: wdHs_dx_ref, wdHs_dy_ref, wPrecip_RL_ref, wPrecip_RL_mod, wdPrecip_RL
    
    REAL(dp), PARAMETER                                :: P_offset = 0.008_dp       ! Normalisation term in precipitation anomaly to avoid divide-by-nearly-zero
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, dHs_dx_ref   , wdHs_dx_ref   )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, dHs_dy_ref   , wdHs_dy_ref   )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, Precip_RL_ref, wPrecip_RL_ref)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, Precip_RL_mod, wPrecip_RL_mod)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, dPrecip_RL   , wdPrecip_RL   )
    
    ! Get surface slopes for the PD_obs reference orography
    CALL ddx_a_to_a_2D( grid, climate%PD_obs%Hs, dHs_dx_ref)
    CALL ddy_a_to_a_2D( grid, climate%PD_obs%Hs, dHs_dy_ref)

    ! Temperature: constant lapse rate plus global offset
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      dT_lapse = (ice%Hs_a( j,i) - climate%PD_obs%Hs( j,i)) * C%constant_lapserate
      DO m = 1, 12
        climate%applied%T2m( m,j,i) = climate%PD_obs%T2m( m,j,i) + dT_lapse + forcing%dT_glob_inverse
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
        
        CALL precipitation_model_Roe( climate%PD_obs%T2m(  m,j,i), dHs_dx_ref(   j,i), dHs_dy_ref(   j,i), climate%PD_obs%Wind_LR( m,j,i), climate%PD_obs%Wind_DU( m,j,i), Precip_RL_ref( m,j,i))
        CALL precipitation_model_Roe( climate%applied%T2m( m,j,i), ice%dHs_dx_a( j,i), ice%dHs_dy_a( j,i), climate%PD_obs%Wind_LR( m,j,i), climate%PD_obs%Wind_DU( m,j,i), Precip_RL_mod( m,j,i))
        dPrecip_RL( m,j,i) = MAX(0.01_dp, MIN( 2._dp, Precip_RL_mod( m,j,i) / Precip_RL_ref( m,j,i) ))
        
        climate%applied%Precip( m,j,i) = climate%PD_obs%Precip( m,j,i) * dPrecip_RL( m,j,i)
        
      END DO
      END DO
      END DO
      CALL sync
    
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
    
      CALL adapt_precip_CC(  grid, ice%Hs_a, climate%PD_obs%Hs, climate%PD_obs%T2m, climate%PD_obs%Precip, climate%applied%Precip, region_name)
    
    END IF
    
    ! Clean up after yourself
    CALL deallocate_shared( wdHs_dx_ref)
    CALL deallocate_shared( wdHs_dy_ref)
    CALL deallocate_shared( wPrecip_RL_ref)
    CALL deallocate_shared( wPrecip_RL_mod)
    CALL deallocate_shared( wdPrecip_RL)
    
    ! Safety
    CALL check_for_NaN_dp_3D( climate%applied%T2m   , 'climate%applied%T2m'   , 'run_climate_model_dT_glob')
    CALL check_for_NaN_dp_3D( climate%applied%Precip, 'climate%applied%Precip', 'run_climate_model_dT_glob')
    
  END SUBROUTINE run_climate_model_dT_glob
  
  ! Climate matrix with PI + LGM snapshots, forced with CO2 (from record or from inverse routine) from Berends et al., 2018
  ! Generalised for different timeframes, L.B. Stap (2021)
  SUBROUTINE run_climate_model_matrix_warm_cold( grid, ice, SMB, climate, region_name, time)
    ! Use CO2 (either prescribed or inversely modelled) to force the 2-snapshot (PI-LGM) climate matrix (Berends et al., 2018)
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables:
    INTEGER                                            :: i,j,m
    
    ! Update insolation forcing at model time
    CALL map_insolation_to_grid( grid, forcing%ins_t0, forcing%ins_t1, forcing%ins_Q_TOA0, forcing%ins_Q_TOA1, time, climate%applied%Q_TOA, climate%applied%Q_TOA_jun_65N, climate%applied%Q_TOA_jan_80S)
    
    ! Use the (CO2 + absorbed insolation)-based interpolation scheme for temperature
    CALL run_climate_model_matrix_warm_cold_temperature( grid, ice, SMB, climate, region_name)
    
    ! Use the (CO2 + ice-sheet geometry)-based interpolation scheme for precipitation
    CALL run_climate_model_matrix_warm_cold_precipitation( grid, ice, climate, region_name)
      
    ! Safety checks
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      IF (climate%applied%T2m( m,j,i) < 150._dp) THEN
        WRITE(0,*) ' WARNING - run_climate_model_matrix_warm_cold: excessively low temperatures (<150K) detected!'
      ELSEIF (climate%applied%T2m( m,j,i) < 0._dp) THEN
        WRITE(0,*) ' ERROR - run_climate_model_matrix_warm_cold: negative temperatures (<0K) detected!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (climate%applied%T2m( m,j,i) /= climate%applied%T2m( m,j,i)) THEN
        WRITE(0,*) ' ERROR - run_climate_model_matrix_warm_cold: NaN temperatures  detected!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (climate%applied%Precip( m,j,i) <= 0._dp) THEN
        WRITE(0,*) ' ERROR - run_climate_model_matrix_warm_cold: zero/negative precipitation detected!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (climate%applied%Precip( m,j,i) /= climate%applied%Precip( m,j,i)) THEN
        WRITE(0,*) ' ERROR - run_climate_model_matrix_warm_cold: NaN precipitation  detected!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE run_climate_model_matrix_warm_cold
  SUBROUTINE run_climate_model_matrix_warm_cold_temperature( grid, ice, SMB, climate, region_name)
    ! The (CO2 + absorbed insolation)-based matrix interpolation for temperature, from Berends et al. (2018)
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables:
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
    
    IF (C%choice_forcing_method == 'CO2_direct') THEN
      CO2 = forcing%CO2_obs
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CO2 = forcing%CO2_mod
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      CO2 = 0._dp
      WRITE(0,*) '  ERROR - run_climate_model_matrix_warm_cold must only be called with the correct forcing method, check your code!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSE
      CO2 = 0._dp
      WRITE(0,*) '  ERROR - choice_forcing_method "', C%choice_forcing_method, '" not implemented in run_climate_model_matrix_warm_cold!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    w_CO2 = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (CO2 - C%climate_matrix_low_CO2_level) / (C%climate_matrix_high_CO2_level - C%climate_matrix_low_CO2_level) ))   ! Berends et al., 2018 - Eq. 1
    
    ! Find the interpolation weights based on absorbed insolation
    ! ===========================================================
    
    ! Calculate modelled absorbed insolation
    climate%applied%I_abs( :,grid%i1:grid%i2) = 0._dp
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      climate%applied%I_abs( j,i) = climate%applied%I_abs( j,i) + climate%applied%Q_TOA( m,j,i) * (1._dp - SMB%Albedo( m,j,i))  ! Berends et al., 2018 - Eq. 2
    END DO
    END DO
    END DO
    CALL sync
    
    ! Calculate weighting field
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      w_ins( j,i) = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (    climate%applied%I_abs(  j,i) -     climate%GCM_cold%I_abs( j,i)) / &  ! Berends et al., 2018 - Eq. 3
                                                           (    climate%GCM_warm%I_abs( j,i) -     climate%GCM_cold%I_abs( j,i)) ))
    END DO
    END DO
    CALL sync
    w_ins_av      = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (SUM(climate%applied%I_abs )      - SUM(climate%GCM_cold%I_abs)     ) / &
                                                           (SUM(climate%GCM_warm%I_abs)      - SUM(climate%GCM_cold%I_abs)     ) ))
   
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
    IF     (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      w_tot( :,grid%i1:grid%i2) = (        w_CO2 + w_ice( :,grid%i1:grid%i2)) / 2._dp  ! Berends et al., 2018 - Eq. 5
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      w_tot( :,grid%i1:grid%i2) = (3._dp * w_CO2 + w_ice( :,grid%i1:grid%i2)) / 4._dp  ! Berends et al., 2018 - Eq. 9
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
    ! TODO Change bias correction, at least for warm_cold climate matrix forcing!
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! Find matrix-interpolated orography, lapse rate, and temperature
      Hs_GCM(       j,i) = (w_tot( j,i) * climate%GCM_warm%Hs(         j,i)) + ((1._dp - w_tot( j,i)) * climate%GCM_cold%Hs(         j,i))  ! Berends et al., 2018 - Eq. 8
      lambda_GCM(   j,i) = (w_tot( j,i) * climate%GCM_warm%lambda(     j,i)) + ((1._dp - w_tot( j,i)) * climate%GCM_cold%lambda(     j,i))  ! Not listed in the article, shame on me!
      T_ref_GCM(  :,j,i) = (w_tot( j,i) * climate%GCM_warm%T2m_corr( :,j,i)) + ((1._dp - w_tot( j,i)) * climate%GCM_cold%T2m_corr( :,j,i))  ! Berends et al., 2018 - Eq. 6
    
      ! Adapt temperature to model orography using matrix-derived lapse-rate
      DO m = 1, 12
        climate%applied%T2m( m,j,i) = T_ref_GCM( m,j,i) - lambda_GCM( j,i) * (ice%Hs_a( j,i) - Hs_GCM( j,i))  ! Berends et al., 2018 - Eq. 11
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
    CALL check_for_NaN_dp_3D( climate%applied%T2m   , 'climate%applied%T2m'   , 'run_climate_model_matrix_warm_cold_temperature')
    
  END SUBROUTINE run_climate_model_matrix_warm_cold_temperature
  SUBROUTINE run_climate_model_matrix_warm_cold_precipitation( grid, ice, climate, region_name)
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
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  w_warm,  w_cold
    INTEGER                                            :: ww_warm, ww_cold
    REAL(dp)                                           :: w_tot
    REAL(dp), DIMENSION(:,:,:), POINTER                :: T_ref_GCM, P_ref_GCM
    REAL(dp), DIMENSION(:,:  ), POINTER                :: Hs_GCM
    INTEGER                                            :: wT_ref_GCM, wP_ref_GCM, wHs_GCM
    
    REAL(dp), PARAMETER                                :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_warm,         ww_warm        )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_cold,         ww_cold        )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, T_ref_GCM,      wT_ref_GCM     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, P_ref_GCM,      wP_ref_GCM     )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, Hs_GCM,         wHs_GCM        )
    
    ! Calculate interpolation weights based on ice geometry
    ! =====================================================
    
    ! First calculate the total ice volume term (second term in the equation)
    w_tot = MAX(-w_cutoff, MIN(1._dp + w_cutoff, (SUM(ice%Hs_a) - SUM(climate%GCM_warm%Hs)) / (SUM(climate%GCM_cold%Hs) - SUM(climate%GCM_warm%Hs)) ))
    
    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      ! Combine total + local ice thicness; Berends et al., 2018, Eq. 12
    
      ! Then the local ice thickness term
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        
        IF (climate%GCM_warm%Hs( j,i) < climate%GCM_PI%Hs( j,i) + 50._dp) THEN
          IF (climate%GCM_cold%Hs( j,i) < climate%GCM_PI%Hs( j,i) + 50._dp) THEN
            ! No ice in any GCM state. Use only total ice volume.
            w_cold( j,i) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, w_tot ))
            w_warm( j,i) = 1._dp - w_cold( j,i)
          ELSE
            ! No ice in warm climate, ice in cold climate. Linear inter- / extrapolation.
            w_cold( j,i) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, ((ice%Hs_a( j,i) - climate%GCM_PI%Hs( j,i)) / (climate%GCM_cold%Hs( j,i) - climate%GCM_PI%Hs( j,i))) * w_tot ))
            w_warm( j,i)  = 1._dp - w_cold( j,i)  
          END IF
        ELSE
          ! Ice in both GCM states.  Linear inter- / extrapolation
          w_cold( j,i) = MAX(-w_cutoff, MIN(1._dp + w_cutoff, ((ice%Hs_a( j,i) - climate%GCM_PI%Hs( j,i)) / (climate%GCM_cold%Hs( j,i) - climate%GCM_PI%Hs( j,i))) * w_tot ))
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
      w_tot = 1._dp - (MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (forcing%CO2_obs - C%climate_matrix_low_CO2_level) / (C%climate_matrix_high_CO2_level - C%climate_matrix_low_CO2_level) )) )
      w_cold( :,grid%i1:grid%i2) = w_tot
      w_warm( :,grid%i1:grid%i2) = 1._dp - w_cold( :,grid%i1:grid%i2)
    END IF     
        
    ! Interpolate the GCM snapshots
    ! =============================
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      T_ref_GCM(  :,j,i) =      (w_warm( j,i) *      climate%GCM_warm%T2m_corr(    :,j,i))  + (w_cold( j,i) *     climate%GCM_cold%T2m_corr(    :,j,i))   ! Berends et al., 2018 - Eq. 6
      P_ref_GCM(  :,j,i) = EXP( (w_warm( j,i) *  LOG(climate%GCM_warm%Precip_corr( :,j,i))) + (w_cold( j,i) * LOG(climate%GCM_cold%Precip_corr( :,j,i)))) ! Berends et al., 2018 - Eq. 7
      Hs_GCM(       j,i) =      (w_warm( j,i) *      climate%GCM_warm%Hs(            j,i))  + (w_cold( j,i) *     climate%GCM_cold%Hs(            j,i))   ! Berends et al., 2018 - Eq. 8
      
    END DO
    END DO
    CALL sync
    
    ! Downscale precipitation from the coarse-resolution reference
    ! GCM orography to the fine-resolution ice-model orography
    ! ========================================================
    
    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      ! Use the Roe&Lindzen precipitation model to do this; Berends et al., 2018, Eqs. A3-A7
      CALL adapt_precip_Roe( grid, Hs_GCM,   T_ref_GCM,           climate%PD_obs%Wind_LR, climate%PD_obs%Wind_DU, P_ref_GCM, &
                                   ice%Hs_a, climate%applied%T2m, climate%PD_obs%Wind_LR, climate%PD_obs%Wind_DU, climate%applied%Precip)
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use a simpler temperature-based correction; Berends et al., 2018, Eq. 14
      CALL adapt_precip_CC( grid, ice%Hs_a, Hs_GCM, T_ref_GCM, P_ref_GCM, climate%applied%Precip, region_name)
    END IF
   
    ! Clean up after yourself
    CALL deallocate_shared( ww_warm)
    CALL deallocate_shared( ww_cold)
    CALL deallocate_shared( wT_ref_GCM)
    CALL deallocate_shared( wP_ref_GCM)
    CALL deallocate_shared( wHs_GCM)
    
    ! Safety
    CALL check_for_NaN_dp_3D( climate%applied%Precip, 'climate%applied%Precip', 'run_climate_model_matrix_warm_cold_precipitation')
    
  END SUBROUTINE run_climate_model_matrix_warm_cold_precipitation
  
  ! Two different parameterised precipitation models:
  ! - a simply Clausius-Clapeyron-based method, used for GRL and ANT
  ! - the Roe & Lindzen temperature/orography-based model, used for NAM and EAS
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
    
    ! Local variables
    INTEGER                                            :: i,j,m
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  T_inv,  T_inv_ref
    INTEGER                                            :: wT_inv, wT_inv_ref
    
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
      IF (par%master) WRITE(0,*) '  ERROR - adapt_precip_CC should only be used for Greenland and Antarctica!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Clean up after yourself
    CALL deallocate_shared( wT_inv)
    CALL deallocate_shared( wT_inv_ref)
    
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
    INTEGER                                            :: i,j,m
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dHs_dx1,  dHs_dx2
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dHs_dy1,  dHs_dy2
    INTEGER                                            :: wdHs_dx1, wdHs_dx2
    INTEGER                                            :: wdHs_dy1, wdHs_dy2
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  Precip_RL1,  Precip_RL2,  dPrecip_RL
    INTEGER                                            :: wPrecip_RL1, wPrecip_RL2, wdPrecip_RL
    
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
        
  END SUBROUTINE precipitation_model_Roe
  
  ! Climate parameterisation for the EISMINT experiments
  SUBROUTINE EISMINT_climate( grid, ice, climate, time)
    ! Simple lapse-rate temperature parameterisation
    
    USe parameters_module,           ONLY: T0, pi
    
    IMPLICIT NONE

    TYPE(type_grid),                     INTENT(IN)    :: grid   
    TYPE(type_ice_model),                INTENT(IN)    :: ice 
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    REAL(dp),                            INTENT(IN)    :: time
    
    REAL(dp), PARAMETER                                :: lambda = -0.010_dp
    
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: dT_lapse, d, dT
    
    ! Set precipitation to zero - SMB is parameterised anyway...
    climate%applied%Precip( :,:,grid%i1:grid%i2) = 0._dp
    
    ! Surface temperature for fixed or moving margin experiments
    IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_3') THEN
      ! Moving margin
          
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
      
        dT_lapse = ice%Hs_a(j,i) * lambda
          
        DO m = 1, 12
          climate%applied%T2m(m,j,i) = 270._dp + dT_lapse
        END DO
        
      END DO
      END DO
      
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_6') THEN
      ! Fixed margin
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
         
        d = MAX( ABS(grid%x(i)/1000._dp), ABS(grid%y(j)/1000._dp))
    
        DO m = 1, 12
          climate%applied%T2m(m,j,i) = 239._dp + (8.0E-08_dp * d**3)
        END DO
        
      END DO
      END DO
      
    END IF
    CALL sync
    
    ! Glacial cycles
    IF     (C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_5') THEN
      IF (time > 0._dp) THEN
        dT = 10._dp * SIN(2._dp * pi * time / 20000._dp)
        climate%applied%T2m( :,:,grid%i1:grid%i2) = climate%applied%T2m( :,:,grid%i1:grid%i2) + dT
      END IF
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_6') THEN
      IF (time > 0._dp) THEN
        dT = 10._dp * SIN(2._dp * pi * time / 40000._dp)
        climate%applied%T2m( :,:,grid%i1:grid%i2) = climate%applied%T2m( :,:,grid%i1:grid%i2) + dT
      END IF
    END IF
    CALL sync
    
  END SUBROUTINE EISMINT_climate

  ! Initialising the climate matrix, containing all the global subclimates
  ! (PD observations and GCM snapshots) on their own lat-lon grids
  SUBROUTINE initialise_climate_matrix( matrix)
    ! Allocate shared memory for the global climate matrix
  
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_climate_matrix),      INTENT(INOUT) :: matrix
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_F' .OR. &
          C%choice_benchmark_experiment == 'MISMIPplus') THEN
        RETURN
      ELSE 
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_climate_matrix!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Initialising the climate matrix...'
    
    ! The global ERA40 climate
    CALL initialise_PD_obs_data_fields( matrix%PD_obs, 'ERA40')
    
    ! The global WOA18 ocean, if needed
    CALL initialise_PD_obs_ocean_fields ( matrix%PD_obs_ocean, 'WOA18')
    
    ! The differenct GCM snapshots 
    IF ((C%choice_forcing_method == 'd18O_inverse_dT_glob') .OR. (C%choice_forcing_method == 'SMB_direct') .OR. (C%choice_forcing_method == 'climate_direct')) THEN
      ! These choices of forcing don't use any GCM (snapshot) data
      RETURN
    ELSEIF (C%choice_forcing_method == 'CO2_direct' .OR. C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! These two choices use the climate matrix
      
      IF (C%choice_climate_matrix == 'warm_cold') THEN
      
        ! Initialise the GCM snapshots
        CALL initialise_snapshot( matrix%GCM_PI,   name = 'ref_PI', nc_filename = C%filename_GCM_snapshot_PI,   CO2 = 280._dp,                         orbit_time =       0._dp)
        CALL initialise_snapshot( matrix%GCM_warm, name = 'warm',   nc_filename = C%filename_GCM_snapshot_warm, CO2 = C%climate_matrix_high_CO2_level, orbit_time = C%climate_matrix_warm_orbit_time)
        CALL initialise_snapshot( matrix%GCM_cold, name = 'cold',   nc_filename = C%filename_GCM_snapshot_cold, CO2 = C%climate_matrix_low_CO2_level,  orbit_time = C%climate_matrix_cold_orbit_time)
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_climate_matrix "', TRIM(C%choice_climate_matrix), '" not implemented in initialise_climate_matrix!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in initialise_climate_matrix!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE initialise_climate_matrix  
  SUBROUTINE initialise_PD_obs_data_fields( PD_obs, name)
    ! Allocate shared memory for the global PD observed climate data fields (stored in the climate matrix),
    ! read them from the specified NetCDF file (latter only done by master process).
     
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_subclimate_global),   INTENT(INOUT) :: PD_obs
    CHARACTER(LEN=*),               INTENT(IN)    :: name
    
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
    IF (par%master) CALL inquire_PD_obs_data_file(PD_obs)
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
    IF (par%master) CALL read_PD_obs_data_file(PD_obs)
    CALL sync
      
    ! Determine process domains
    CALL partition_list( PD_obs%nlon, par%i, par%n, PD_obs%i1, PD_obs%i2)
    
  END SUBROUTINE initialise_PD_obs_data_fields   
  SUBROUTINE initialise_snapshot( snapshot, name, nc_filename, CO2, orbit_time)
    ! Allocate shared memory for the data fields of a GCM snapshot (stored in the climate matrix),
    ! read them from the specified NetCDF file (latter only done by master process).
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_subclimate_global),   INTENT(INOUT) :: snapshot
    CHARACTER(LEN=*),               INTENT(IN)    :: name
    CHARACTER(LEN=*),               INTENT(IN)    :: nc_filename
    REAL(dp),                       INTENT(IN)    :: CO2
    REAL(dp),                       INTENT(IN)    :: orbit_time
    
    ! Local variables:
    INTEGER                                       :: i,j,m
    REAL(dp), PARAMETER                           :: Precip_minval = 1E-5_dp
    
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
    IF (par%master) CALL inquire_GCM_snapshot( snapshot)
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
    IF (par%master) CALL read_GCM_snapshot( snapshot)
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
      IF (snapshot%T2m( i,j,m) < 150._dp) THEN
        WRITE(0,*) ' WARNING - initialise_snapshot: excessively low temperatures (<150K) detected in snapshot ', snapshot%name, '!'
      ELSEIF (snapshot%T2m( i,j,m) < 0._dp) THEN
        WRITE(0,*) ' ERROR - initialise_snapshot: negative temperatures (<0K) detected in snapshot ', snapshot%name, '!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (snapshot%T2m( i,j,m) /= snapshot%T2m( i,j,m)) THEN
        WRITE(0,*) ' ERROR - initialise_snapshot: NaN temperatures  detected in snapshot ', snapshot%name, '!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (snapshot%Precip( i,j,m) <= 0._dp) THEN
        WRITE(0,*) ' ERROR - initialise_snapshot: zero/negative precipitation detected in snapshot ', snapshot%name, '!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (snapshot%Precip( i,j,m) /= snapshot%Precip( i,j,m)) THEN
        WRITE(0,*) ' ERROR - initialise_snapshot: NaN precipitation  detected in snapshot ', snapshot%name, '!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE initialise_snapshot
  
  ! Initialising the region-specific climate model, containing all the subclimates
  ! (PD observations, GCM snapshots and the applied climate) on the model grid
  SUBROUTINE initialise_climate_model( grid, climate, matrix, region_name, mask_noice)
    ! Allocate shared memory for the regional climate models, containing the PD observed,
    ! GCM snapshots and applied climates as "subclimates"
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    TYPE(type_climate_matrix),           INTENT(IN)    :: matrix
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_noice
    
    IF (par%master) WRITE (0,*) '  Initialising climate model...'
    
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
          C%choice_benchmark_experiment == 'ISMIP_HOM_F' .OR. &
          C%choice_benchmark_experiment == 'MISMIPplus') THEN
          
        ! Entirely parameterised climate, no ocean
        CALL initialise_subclimate( grid, climate%applied, 'applied')
        RETURN
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_climate_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
        
    ! Initialise data structures for the regional ERA40 climate and the final applied climate
    CALL initialise_subclimate( grid, climate%PD_obs,   'ERA40'  )
    CALL initialise_subclimate( grid, climate%applied,  'applied')
    
    ! Map these subclimates from global grid to model grid
    CALL map_subclimate_to_grid( grid,  matrix%PD_obs,  climate%PD_obs)
   
    ! The differenct GCM snapshots
    IF (C%choice_forcing_method == 'd18O_inverse_dT_glob' .OR. C%choice_forcing_method == 'SMB_direct' .OR. C%choice_forcing_method == 'climate_direct') THEN
      ! These choices of forcing don't use any GCM (snapshot) data
      RETURN
    ELSEIF (C%choice_forcing_method == 'CO2_direct' .OR. C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! These two choices use the climate matrix
      
      IF (C%choice_climate_matrix == 'warm_cold') THEN
      
        ! Initialise data structures for the GCM snapshots
        CALL initialise_subclimate( grid, climate%GCM_PI,   'ref_PI' )
        CALL initialise_subclimate( grid, climate%GCM_warm, 'Warm' )
        CALL initialise_subclimate( grid, climate%GCM_cold, 'Cold')
        
        ! Map these subclimates from global grid to model grid
        CALL map_subclimate_to_grid( grid,  matrix%GCM_PI,   climate%GCM_PI  )
        CALL map_subclimate_to_grid( grid,  matrix%GCM_warm, climate%GCM_warm)
        CALL map_subclimate_to_grid( grid,  matrix%GCM_cold, climate%GCM_cold)
        
        ! Right now, no wind is read from GCM output; just use PD observations everywhere
        climate%GCM_warm%Wind_WE( :,:,grid%i1:grid%i2) = climate%PD_obs%Wind_WE( :,:,grid%i1:grid%i2)
        climate%GCM_warm%Wind_SN( :,:,grid%i1:grid%i2) = climate%PD_obs%Wind_SN( :,:,grid%i1:grid%i2)
        climate%GCM_cold%Wind_WE( :,:,grid%i1:grid%i2) = climate%PD_obs%Wind_WE( :,:,grid%i1:grid%i2)
        climate%GCM_cold%Wind_SN( :,:,grid%i1:grid%i2) = climate%PD_obs%Wind_SN( :,:,grid%i1:grid%i2)
  
        ! Calculate spatially variable lapse rate
        climate%GCM_warm%lambda = 0.008_dp
        IF     (region_name == 'NAM' .OR. region_name == 'EAS') THEN
          CALL initialise_subclimate_spatially_variable_lapserate( grid, climate%GCM_PI, climate%GCM_cold)
        ELSEIF (region_name == 'GLR' .OR. region_name == 'ANT') THEN
          climate%GCM_cold%lambda = 0.008_dp
        END IF
        
        ! Calculate and apply GCM bias correction
        CALL calculate_GCM_bias( grid, climate)
        CALL correct_GCM_bias(   grid, climate, climate%GCM_warm, do_correct_bias = C%climate_matrix_biascorrect_warm)
        CALL correct_GCM_bias(   grid, climate, climate%GCM_cold, do_correct_bias = C%climate_matrix_biascorrect_cold)
        
        ! Get reference absorbed insolation for the GCM snapshots
        CALL initialise_subclimate_absorbed_insolation( grid, climate%GCM_warm, region_name, mask_noice)
        CALL initialise_subclimate_absorbed_insolation( grid, climate%GCM_cold, region_name, mask_noice)
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_climate_matrix "', TRIM(C%choice_climate_matrix), '" not implemented in initialise_climate_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in initialise_climate_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
  
    ! Initialise applied climate with present-day observations
    IF (par%master) THEN
      climate%applied%T2m     = climate%PD_obs%T2m
      climate%applied%Precip  = climate%PD_obs%Precip
      climate%applied%Hs      = climate%PD_obs%Hs
      climate%applied%Wind_LR = climate%PD_obs%Wind_LR
      climate%applied%Wind_DU = climate%PD_obs%Wind_DU
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE initialise_climate_model  
  SUBROUTINE initialise_subclimate( grid, subclimate, name)
    ! Allocate shared memory for a "subclimate" (PD observed, GCM snapshot or applied climate) on the grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_subclimate_region),        INTENT(INOUT) :: subclimate
    CHARACTER(LEN=*),                    INTENT(IN)    :: name
    
    subclimate%name = name
    
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, subclimate%Hs,             subclimate%wHs            )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%T2m,            subclimate%wT2m           )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Precip,         subclimate%wPrecip        )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Wind_WE,        subclimate%wWind_WE       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Wind_SN,        subclimate%wWind_SN       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Wind_LR,        subclimate%wWind_LR       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Wind_DU,        subclimate%wWind_DU       )
      
    CALL allocate_shared_dp_0D(                       subclimate%CO2,            subclimate%wCO2           )
    CALL allocate_shared_dp_0D(                       subclimate%orbit_time,     subclimate%worbit_time    )
    CALL allocate_shared_dp_0D(                       subclimate%orbit_ecc,      subclimate%worbit_ecc     )
    CALL allocate_shared_dp_0D(                       subclimate%orbit_obl,      subclimate%worbit_obl     )
    CALL allocate_shared_dp_0D(                       subclimate%orbit_pre,      subclimate%worbit_pre     )
    CALL allocate_shared_dp_0D(                       subclimate%sealevel,       subclimate%wsealevel      )
    
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, subclimate%lambda,         subclimate%wlambda        )
    
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%T2m_corr,       subclimate%wT2m_corr      )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Precip_corr,    subclimate%wPrecip_corr   )
    
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Q_TOA,          subclimate%wQ_TOA         )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Albedo,         subclimate%wAlbedo        )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, subclimate%I_abs,          subclimate%wI_abs         )
    CALL allocate_shared_dp_0D(                       subclimate%Q_TOA_jun_65N,  subclimate%wQ_TOA_jun_65N )
    CALL allocate_shared_dp_0D(                       subclimate%Q_TOA_jan_80S,  subclimate%wQ_TOA_jan_80S )
    
    CALL allocate_shared_dp_0D(                       subclimate%T_ocean_mean,   subclimate%wT_ocean_mean  )
    
  END SUBROUTINE initialise_subclimate
  SUBROUTINE calculate_GCM_bias( grid, climate)
    ! Calculate the GCM bias in temperature and precipitation
    !
    ! Account for the fact that the GCM PI snapshot has a lower resolution, and therefore
    ! a different surface elevation than the PD observed climatology!
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    
    ! Local variables:
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: T2m_SL_GCM, T2m_SL_obs
    
    ! Allocate shared memory
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%GCM_bias_T2m,    climate%wGCM_bias_T2m   )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, climate%GCM_bias_Precip, climate%wGCM_bias_Precip)
    
    ! Calculate bias
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
    
      ! Scale modelled and observed temperature to sea level using a constant lapse rate
      T2m_SL_obs = climate%PD_obs%T2m( m,j,i) + climate%PD_obs%Hs( j,i) * C%constant_lapserate
      T2m_SL_GCM = climate%GCM_PI%T2m( m,j,i) + climate%GCM_PI%Hs( j,i) * C%constant_lapserate
      
      ! Calculate temperature bias
      climate%GCM_bias_T2m(    m,j,i) = T2m_SL_GCM - T2m_SL_obs
      
      ! Calculate precipitation bias
      climate%GCM_bias_Precip( m,j,i) = climate%GCM_PI%Precip( m,j,i) / climate%PD_obs%Precip( m,j,i)
      
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calculate_GCM_bias
  SUBROUTINE correct_GCM_bias( grid, climate, subclimate, do_correct_bias)
    ! Calculate bias-corrected climate for this snapshot
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_climate_model),            INTENT(IN)    :: climate
    TYPE(type_subclimate_region),        INTENT(INOUT) :: subclimate
    LOGICAL,                             INTENT(IN)    :: do_correct_bias
    
    ! Local variables:
    INTEGER                                            :: i,j,m
    
    ! If no bias correction should be applied, set T2m_corr = T2m and Precip_corr = Precip
    IF (.NOT. do_correct_bias) THEN
      subclimate%T2m_corr(    :,:,grid%i1:grid%i2) = subclimate%T2m(    :,:,grid%i1:grid%i2)
      subclimate%Precip_corr( :,:,grid%i1:grid%i2) = subclimate%Precip( :,:,grid%i1:grid%i2)
      CALL sync
      RETURN
    END IF
    
    ! Apply bias correction
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      
      ! Temperature
      subclimate%T2m_corr(    m,j,i) = subclimate%T2m(    m,j,i) - climate%GCM_bias_T2m(    m,j,i)
      
      ! Precipitation
      subclimate%Precip_corr( m,j,i) = subclimate%Precip( m,j,i) / climate%GCM_bias_Precip( m,j,i)
      
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE correct_GCM_bias
  SUBROUTINE initialise_subclimate_spatially_variable_lapserate( grid, snapshot_PI, snapshot)
    ! Calculate the spatially variable lapse-rate (for non-PI GCM snapshots; see Berends et al., 2018)
    ! Only meaningful for snapshots where there is ice (LGM, M2_Medium, M2_Large),
    ! and only intended for North America and Eurasia
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_subclimate_region),        INTENT(IN)    :: snapshot_PI 
    TYPE(type_subclimate_region),        INTENT(INOUT) :: snapshot
    
    ! Local variables:
    INTEGER                                            :: i,j,m
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  mask_calc_lambda
    INTEGER                                            :: wmask_calc_lambda
    REAL(dp)                                           :: dT_mean_nonice
    INTEGER                                            :: n_nonice, n_ice
    REAL(dp)                                           :: lambda_mean_ice
    
    REAL(dp), PARAMETER                                :: lambda_min = 0.002_dp
    REAL(dp), PARAMETER                                :: lambda_max = 0.05_dp
    
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
    snapshot%lambda( :,grid%i1:grid%i2) = snapshot%lambda( :,grid%i1:grid%i2) * (0.008_dp / lambda_mean_ice)
    
    ! Clean up after yourself
    CALl deallocate_shared( wmask_calc_lambda)
    
  END SUBROUTINE initialise_subclimate_spatially_variable_lapserate
  SUBROUTINE initialise_subclimate_absorbed_insolation( grid, snapshot, region_name, mask_noice)
    ! Calculate the yearly absorbed insolation for this (regional) GCM snapshot, to be used in the matrix interpolation
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: snapshot
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_noice
    
    ! Local variables:
    INTEGER                                            :: i,j,m
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  Q_TOA0,  Q_TOA1
    INTEGER                                            :: wQ_TOA0, wQ_TOA1
    INTEGER                                            :: ti0, ti1
    REAL(dp)                                           :: ins_t0, ins_t1
    REAL(dp)                                           :: wt0, wt1
    INTEGER                                            :: ilat_l, ilat_u
    REAL(dp)                                           :: wlat_l, wlat_u
    
    TYPE(type_ice_model)                               :: ice_dummy
    TYPE(type_subclimate_region)                       :: climate_dummy
    TYPE(type_SMB_model)                               :: SMB_dummy
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( forcing%ins_nlat, 12, Q_TOA0, wQ_TOA0)
    CALL allocate_shared_dp_2D( forcing%ins_nlat, 12, Q_TOA1, wQ_TOA1)
    
    ! Get insolation at the desired time from the insolation NetCDF file
    ! ==================================================================
    
    ins_t0 = 0
    ins_t1 = 0
    
    ! Find time indices to be read
    IF (snapshot%orbit_time >= MINVAL(forcing%ins_time) .AND. snapshot%orbit_time <= MAXVAL(forcing%ins_time)) THEN
      ti1 = 1
      DO WHILE (forcing%ins_time( ti1) < snapshot%orbit_time)
        ti1 = ti1 + 1
      END DO
      ti0 = ti1 - 1
      
      ins_t0 = forcing%ins_time( ti0)
      ins_t1 = forcing%ins_time( ti1)
    ELSE
      WRITE(0,*) '  ERROR - orbit_time ', snapshot%orbit_time, ' for snapshot ', TRIM(snapshot%name), ' outside of range of insolation solution file "', TRIM(forcing%netcdf_ins%filename), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    CALL sync
    
    ! Read insolation time frames enveloping desired time from netcdf file
    IF (par%master) CALL read_insolation_data_file( forcing, ti0, ti1, Q_TOA0, Q_TOA1)
    CALL sync
    
    ! Map monthly insolation at the top of the atmosphere to the region grid
    ! Calculate time interpolation weights
    wt0 = (ins_t1 - snapshot%orbit_time) / (ins_t1 - ins_t0)
    wt1 = 1._dp - wt0
        
    ! Interpolate on the grid
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
     
      ilat_l = FLOOR(grid%lat(j,i) + 91)
      ilat_u = ilat_l + 1
      
      wlat_l = forcing%ins_lat(ilat_u) - grid%lat(j,i)
      wlat_u = 1._dp - wlat_l
      
      DO m = 1, 12
        snapshot%Q_TOA( m,j,i) = (wt0 * wlat_l * Q_TOA0( ilat_l,m)) + &
                                 (wt0 * wlat_u * Q_TOA0( ilat_u,m)) + &
                                 (wt1 * wlat_l * Q_TOA1( ilat_l,m)) + &
                                 (wt1 * wlat_u * Q_TOA1( ilat_u,m))
      END DO    
    END DO
    END DO
    CALL sync
    
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
    climate_dummy%T2m(    :,:,grid%i1:grid%i2) = snapshot%T2m_corr(    :,:,grid%i1:grid%i2)
    climate_dummy%Precip( :,:,grid%i1:grid%i2) = snapshot%Precip_corr( :,:,grid%i1:grid%i2)
    climate_dummy%Q_TOA(  :,:,grid%i1:grid%i2) = snapshot%Q_TOA(       :,:,grid%i1:grid%i2)
    
    ! Ice
    ! ===
    
    CALL allocate_shared_int_2D(    grid%ny, grid%nx, ice_dummy%mask_ocean_a   , ice_dummy%wmask_ocean_a   )
    CALL allocate_shared_int_2D(    grid%ny, grid%nx, ice_dummy%mask_ice_a     , ice_dummy%wmask_ice_a     )
    CALL allocate_shared_int_2D(    grid%ny, grid%nx, ice_dummy%mask_shelf_a   , ice_dummy%wmask_shelf_a   )
    
    ! Fill in masks for the SMB model
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      IF (snapshot%Hs( j,i) == MINVAL(snapshot%Hs)) THEN
        ice_dummy%mask_ocean_a( j,i) = 1
      ELSE
        ice_dummy%mask_ocean_a( j,i) = 0
      END IF
      
      IF (snapshot%Hs( j,i) > 100._dp .AND. SUM(snapshot%T2m( :,j,i)) / 12._dp < 0._dp) THEN
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
    
    ! Run the SMB model for 10 years for this particular snapshot
    ! (experimentally determined to be long enough to converge)
    DO i = 1, 10
      CALL run_SMB_model( grid, ice_dummy, climate_dummy, 0._dp, SMB_dummy, mask_noice)
    END DO
    CALL sync
    
    ! Copy the resulting albedo to the snapshot
    snapshot%Albedo( :,:,grid%i1:grid%i2) = SMB_dummy%Albedo( :,:,grid%i1:grid%i2)
    CALL sync
    
    ! Calculate yearly total absorbed insolation
    snapshot%I_abs( :,grid%i1:grid%i2) = 0._dp
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      snapshot%I_abs( j,i) = snapshot%I_abs( j,i) + snapshot%Q_TOA( m,j,i) * (1._dp - snapshot%Albedo( m,j,i))
    END DO
    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wQ_TOA0)
    CALL deallocate_shared( wQ_TOA1)
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
    
  END SUBROUTINE initialise_subclimate_absorbed_insolation
  
  ! Map a global subclimate from the matrix (PD observed or GCM snapshot) to a region grid
  SUBROUTINE map_subclimate_to_grid( grid,  cglob, creg)
    ! Map data from a global "subclimate" (PD observed or GCM snapshot) to the grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_subclimate_global),        INTENT(IN)    :: cglob  ! Global climate
    TYPE(type_subclimate_region),        INTENT(INOUT) :: creg   ! grid   climate
    
    ! Local variables:
    INTEGER                                            :: i,j,m
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  logPrecip
    INTEGER                                            :: wlogPrecip
    
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
    CALL map_glob_to_grid_2D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%Hs,      creg%Hs         )
    CALL map_glob_to_grid_3D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%T2m,     creg%T2m    , 12)
    CALL map_glob_to_grid_3D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, logPrecip,     creg%Precip , 12)
    CALL map_glob_to_grid_3D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%Wind_WE, creg%Wind_WE, 12)
    CALL map_glob_to_grid_3D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%Wind_SN, creg%Wind_SN, 12)
    
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
  
  END SUBROUTINE map_subclimate_to_grid
  
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
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: longitude_start, Uwind_x, Uwind_y, Vwind_x, Vwind_y

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
    
  END SUBROUTINE rotate_wind_to_model_grid


END MODULE climate_module
