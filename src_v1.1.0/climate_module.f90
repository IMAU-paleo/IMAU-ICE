MODULE climate_module

  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parallel_module,                 ONLY: par, sync, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_climate_model, type_SMB_model, &
                                             type_climate_matrix, type_subclimate_global, type_subclimate_region
  USE netcdf_module,                   ONLY: inquire_PD_obs_data_file, read_PD_obs_data_file, inquire_GCM_snapshot, read_GCM_snapshot,&
                                             read_insolation_data_file
  USE forcing_module,                  ONLY: forcing, map_insolation_to_grid
  USE utilities_module,                ONLY: error_function, smooth_Gaussian_2D
  USE derivatives_and_grids_module,    ONLY: ddx_Aa_to_Aa_2D, ddy_Aa_to_Aa_2D
  USE SMB_module,                      ONLY: run_SMB_model

  IMPLICIT NONE
    
CONTAINS

  ! Run the climate model on a region grid
  SUBROUTINE run_climate_model( grid, ice, SMB, climate, region_name, time)
    ! Run the regional climate model
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT)    :: ice 
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    
    ! Check if we need to apply any special benchmark experiment climate
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
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
              C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_climate_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Update insolation forcing at model time
    CALL map_insolation_to_grid( grid, forcing%ins_t0, forcing%ins_t1, forcing%ins_Q_TOA0, forcing%ins_Q_TOA1, time, climate%applied%Q_TOA)
        
    ! Different kinds of climate forcing for realistic experiments
    IF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      ! Use the global temperature offset as calculated by the inverse routine
      
      CALL run_climate_model_dT_glob( grid, ice, climate, region_name)
      
    ELSEIF (C%choice_forcing_method == 'CO2_direct' .OR. C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! Use CO2 (either taken directly from the specified record, or as calculated by the inverse routine) to force the climate matrix
      
      IF (C%choice_climate_matrix == 'PI_LGM') THEN
        ! Use the two-snapshot climate matrix
        CALL run_climate_model_matrix_PI_LGM( grid, ice, SMB, climate, region_name)
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_climate_matrix "', TRIM(C%choice_climate_matrix), '" not implemented in run_climate_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: forcing method "', TRIM(C%choice_forcing_method), '" not implemented in run_climate_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE run_climate_model
  
  ! Parameterised climate (ERA40 + global temperature offset) from de Boer et al., 2013
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
    CALL ddx_Aa_to_Aa_2D( grid, climate%PD_obs%Hs, dHs_dx_ref)
    CALL ddy_Aa_to_Aa_2D( grid, climate%PD_obs%Hs, dHs_dy_ref)

    ! Temperature: constant lapse rate plus global offset
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      dT_lapse = (ice%Hs_Aa( j,i) - climate%PD_obs%Hs( j,i)) * C%constant_lapserate
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
        
        CALL precipitation_model_Roe( climate%PD_obs%T2m(  m,j,i), dHs_dx_ref(    j,i), dHs_dy_ref(    j,i), climate%PD_obs%Wind_LR( m,j,i), climate%PD_obs%Wind_DU( m,j,i), Precip_RL_ref( m,j,i))
        CALL precipitation_model_Roe( climate%applied%T2m( m,j,i), ice%dHs_dx_Aa( j,i), ice%dHs_dy_Aa( j,i), climate%PD_obs%Wind_LR( m,j,i), climate%PD_obs%Wind_DU( m,j,i), Precip_RL_mod( m,j,i))
        dPrecip_RL( m,j,i) = MAX(0.01_dp, MIN( 2._dp, (Precip_RL_mod( m,j,i) + P_offset) / (Precip_RL_ref( m,j,i) + P_offset) ))
        
        climate%applied%Precip( m,j,i) = climate%PD_obs%Precip( m,j,i) * dPrecip_RL( m,j,i)
        
      END DO
      END DO
      END DO
      CALL sync
    
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
      DO m = 1, 12
        
        CALL adjust_precipitation_for_temperature_change( climate%PD_obs%T2m( m,j,i), climate%applied%T2m( m,j,i), climate%PD_obs%Precip( m,j,i), climate%applied%Precip( m,j,i))
        
      END DO
      END DO
      END DO
      CALL sync
    
    END IF
    
    ! Clean up after yourself
    CALL deallocate_shared( wdHs_dx_ref)
    CALL deallocate_shared( wdHs_dy_ref)
    CALL deallocate_shared( wPrecip_RL_ref)
    CALL deallocate_shared( wPrecip_RL_mod)
    CALL deallocate_shared( wdPrecip_RL)
    
  END SUBROUTINE run_climate_model_dT_glob
  
  ! Climate matrix with PI + LGM snapshots, forced with CO2 (from record or from inverse routine) from Berends et al., 2018
  SUBROUTINE run_climate_model_matrix_PI_LGM( grid, ice, SMB, climate, region_name)
    ! Use CO2 (either prescribed or inversely modelled) to force the 2-snapshot (PI-LGM) climate matrix (Berends et al., 2018)
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT)    :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Use the (CO2 + absorbed insolation)-based interpolation scheme for temperature
    CALL run_climate_model_matrix_PI_LGM_temperature( grid, ice, SMB, climate, region_name)
    
    ! Use the (CO2 + ice-sheet geometry)-based interpolation scheme for precipitation
    CALL run_climate_model_matrix_PI_LGM_precipitation( grid, ice, climate, region_name)
    
  END SUBROUTINE run_climate_model_matrix_PI_LGM
  SUBROUTINE run_climate_model_matrix_PI_LGM_temperature( grid, ice, SMB, climate, region_name)
    ! The (CO2 + absorbed insolation)-based matrix interpolation for temperature, from Berends et al. (2018)
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT)    :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: CO2, w_CO2
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  w_ins,  w_ins_smooth,  w_ice,  w_tot
    INTEGER                                            :: ww_ins, ww_ins_smooth, ww_ice, ww_tot
    REAL(dp)                                           :: w_ins_av
    REAL(dp), DIMENSION(:,: ,:), POINTER                :: T_pot_ref_GCM
    REAL(dp), DIMENSION(:,:  ), POINTER                :: Hs_ref_GCM, lambda_ref_GCM
    INTEGER                                            :: wT_pot_ref_GCM, wHs_ref_GCM, wlambda_ref_GCM
    
    REAL(dp), PARAMETER                                :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]
    REAL(dp), PARAMETER                                :: P_offset = 0.008_dp       ! Normalisation term in precipitation anomaly to avoid divide-by-nearly-zero
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_ins,          ww_ins         )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_ins_smooth,   ww_ins_smooth  )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_ice,          ww_ice         )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_tot,          ww_tot         )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, T_pot_ref_GCM,  wT_pot_ref_GCM )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, Hs_ref_GCM,     wHs_ref_GCM    )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, lambda_ref_GCM, wlambda_ref_GCM)
    
    ! Use either prescribed or modelled CO2
    CO2 = 0._dp
    IF (C%choice_forcing_method == 'CO2_direct') THEN
      CO2 = forcing%CO2_obs
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CO2 = forcing%CO2_mod
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      WRITE(0,*) '  ERROR - run_climate_model_matrix_PI_LGM must only be called with the correct forcing method, check your code!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSE
      WRITE(0,*) '  ERROR - choice_forcing_method "', C%choice_forcing_method, '" not implemented in run_climate_model_matrix_PI_LGM!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Find CO2 interpolation weights
    w_CO2 = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (CO2 - 190._dp) / (280._dp - 190._dp) ))   ! Berends et al., 2018 - Eq. 1
    
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
    
    ! Find the interpolation weights based on absorbed insolation (used to determine temperature)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      w_ins( j,i) = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (    climate%applied%I_abs( j,i) -     climate%GCM_LGM%I_abs( j,i)) / &  ! Berends et al., 2018 - Eq. 3
                                                           (    climate%GCM_PI%I_abs(  j,i) -     climate%GCM_LGM%I_abs( j,i)) ))
    END DO
    END DO
    CALL sync
    w_ins_av      = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (SUM(climate%applied%I_abs)      - SUM(climate%GCM_LGM%I_abs)     ) / &
                                                           (SUM(climate%GCM_PI%I_abs )      - SUM(climate%GCM_LGM%I_abs)     ) ))
   
    ! Smooth the weighting field (Berends et al., 2018, Eq. 4)
    w_ins_smooth( :,grid%i1:grid%i2) = w_ins( :,grid%i1:grid%i2)
    CALL smooth_Gaussian_2D( grid, w_ins_smooth, 5 * CEILING(40000._dp / grid%dx))
    w_ice( :,grid%i1:grid%i2) = (1._dp * w_ins(        :,grid%i1:grid%i2) + &
                                 3._dp * w_ins_smooth( :,grid%i1:grid%i2) + &
                                 3._dp * w_ins_av) / 7._dp
     
    ! Combine interpolation weights from absorbed insolation and CO2 into the final weights fields
    IF     (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      w_tot( :,grid%i1:grid%i2) = (        w_CO2 + w_ice( :,grid%i1:grid%i2)) / 2._dp  ! Berends et al., 2018 - Eq. 5
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      w_tot( :,grid%i1:grid%i2) = (3._dp * w_CO2 + w_ice( :,grid%i1:grid%i2)) / 4._dp  ! Berends et al., 2018 - Eq. 9
    END IF

    ! Find matrix-interpolated temperature
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! Find matrix-interpolated orography, lapse rate, and temperature
      Hs_ref_GCM(       j,i) = (w_tot( j,i) * climate%GCM_PI%Hs(             j,i)) + ((1._dp - w_tot( j,i)) * climate%GCM_LGM%Hs(             j,i))  ! Berends et al., 2018 - Eq. 8
      lambda_ref_GCM(   j,i) = (w_tot( j,i) * climate%GCM_PI%lambda(         j,i)) + ((1._dp - w_tot( j,i)) * climate%GCM_LGM%lambda(         j,i))  ! Not listed in the article, shame on me!
      T_pot_ref_GCM(  :,j,i) = (w_tot( j,i) * climate%GCM_PI%T2m_pot_corr( :,j,i)) + ((1._dp - w_tot( j,i)) * climate%GCM_LGM%T2m_pot_corr( :,j,i))  ! Berends et al., 2018 - Eq. 6
    
      ! Adapt temperature to model orography using matrix-derived lapse-rate
      DO m = 1, 12
        climate%applied%T2m( m,j,i) = T_pot_ref_GCM( m,j,i) + lambda_ref_GCM( j,i) * ice%Hs_Aa( j,i)  ! Berends et al., 2018 - Eq. 11
      END DO
    
    END DO
    END DO
    CALL sync
   
    ! Clean up after yourself
    CALL deallocate_shared( ww_ins)
    CALL deallocate_shared( ww_ins_smooth)
    CALL deallocate_shared( ww_ice)
    CALL deallocate_shared( ww_tot)
    CALL deallocate_shared( wT_pot_ref_GCM)
    CALL deallocate_shared( wHs_ref_GCM)
    CALL deallocate_shared( wlambda_ref_GCM)
    
  END SUBROUTINE run_climate_model_matrix_PI_LGM_temperature
  SUBROUTINE run_climate_model_matrix_PI_LGM_precipitation( grid, ice, climate, region_name)
    ! The (CO2 + ice geometry)-based matrix interpolation for precipitation, from Berends et al. (2018)
    ! This was developed specifically for North America and Eurasia, where glacial-interglacial
    ! changes in ice-sheet geometry are much more dramatic than in Greenland and Antarctica.
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: CO2, w_CO2
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  w_ice,  w_tot
    INTEGER                                            :: ww_ice, ww_tot
    REAL(dp)                                           :: w_ice_tot
    REAL(dp), DIMENSION(:,:,:), POINTER                :: T_pot_ref_GCM, T_ref_GCM, P_ref_GCM
    REAL(dp), DIMENSION(:,:  ), POINTER                :: lambda_GCM, Hs_ref_GCM
    INTEGER                                            :: wT_pot_ref_GCM, wlambda_GCM, wT_ref_GCM, wP_ref_GCM, wHs_ref_GCM
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dHs_dx_ref_GCM,  dHs_dy_ref_GCM
    INTEGER                                            :: wdHs_dx_ref_GCM, wdHs_dy_ref_GCM
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  P_RL_ref_GCM,  P_RL_mod,  dP_RL
    INTEGER                                            :: wP_RL_ref_GCM, wP_RL_mod, wdP_RL
    
    REAL(dp), PARAMETER                                :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]
    REAL(dp), PARAMETER                                :: P_offset = 0.008_dp       ! Normalisation term in precipitation anomaly to avoid divide-by-nearly-zero
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_ice,          ww_ice         )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, w_tot,          ww_tot         )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, T_pot_ref_GCM,  wT_pot_ref_GCM )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, lambda_GCM,     wlambda_GCM    )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, T_ref_GCM,      wT_ref_GCM     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, P_ref_GCM,      wP_ref_GCM     )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, Hs_ref_GCM,     wHs_ref_GCM    )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, dHs_dx_ref_GCM, wdHs_dx_ref_GCM)
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, dHs_dy_ref_GCM, wdHs_dy_ref_GCM)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, P_RL_ref_GCM,   wP_RL_ref_GCM  )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, P_RL_mod,       wP_RL_mod      )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, dP_RL,          wdP_RL         )
    
    ! Use either prescribed or modelled CO2
    CO2 = 0._dp
    IF (C%choice_forcing_method == 'CO2_direct') THEN
      CO2 = forcing%CO2_obs
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CO2 = forcing%CO2_mod
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      WRITE(0,*) '  ERROR - run_climate_model_matrix_PI_LGM must only be called with the correct forcing method, check your code!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSE
      WRITE(0,*) '  ERROR - choice_forcing_method "', C%choice_forcing_method, '" not implemented in run_climate_model_matrix_PI_LGM!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Find CO2 interpolation weights
    w_CO2 = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (CO2 - 190._dp) / (280._dp - 190._dp) ))   ! Berends et al., 2018 - Eq. 1
    
    ! Calculate interpolation weights based on ice geometry
    
    ! First calculate the total ice volume term (second term in the equation)
    w_ice_tot = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (SUM( climate%GCM_LGM%Hi ) - SUM( climate%GCM_PI%Hi - ice%Hi_Aa)) / (SUM( climate%GCM_LGM%Hi ) - SUM( climate%GCM_PI%Hi )) ))
    
    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      ! Combine total + local ice thicness; Berends et al., 2018, Eq. 12
    
      ! Then the local ice thickness term
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
       
        IF (climate%GCM_PI%Hi( j,i) == 0._dp) THEN
          IF (climate%GCM_LGM%Hi( j,i) == 0._dp) THEN
            ! This pixel has no ice in any GCM state. Use only total ice volume.
            w_ice( j,i) = 1._dp
          ELSE
            ! No ice at PD, ice at LGM. Linear inter- / extrapolation.
            w_ice( j,i) = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (climate%GCM_LGM%Hi( j,i) - ice%Hi_Aa( j,i)) / climate%GCM_LGM%Hi( j,i) ))
          END IF
        ELSE
          ! Ice in both GCM states.  Linear inter- / extrapolation (not possible in NAM or EAS?)
          w_ice(   j,i) = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (climate%GCM_LGM%Hi( j,i) - ice%Hi_Aa( j,i)) / (climate%GCM_LGM%Hi( j,i) - climate%GCM_PI%Hi( j,i)) ))
        END IF 
        
        ! Multiply with total ice volume term
        w_ice( j,i) = w_ice( j,i) * w_ice_tot
        
      END DO
      END DO
      CALL sync
      
      ! Smooth the weighting field
      CALL smooth_Gaussian_2D( grid, w_ice, 5 * CEILING(40000._dp / grid%dx))    
      
      ! Combine with the CO2 weight to get the final interpolation weights field
      w_tot( :,grid%i1:grid%i2) = (w_CO2 + w_ice( :,grid%i1:grid%i2)) / 2._dp
      
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use only local ice thickness; Berendset al., 2018, Eq. 13
      
      w_tot( :,grid%i1:grid%i2) = (w_CO2 + w_ice_tot) / 2._dp
      
    END IF
        
    ! Interpolate the GCM snapshots
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      T_pot_ref_GCM(  :,j,i) =      (w_tot( j,i) *     climate%GCM_PI%T2m_pot_corr( :,j,i)            ) + ((1._dp - w_tot( j,i)) *     climate%GCM_LGM%T2m_pot_corr( :,j,i)            )   ! Berends et al., 2018 - Eq. 6
      lambda_GCM(       j,i) =      (w_tot( j,i) *     climate%GCM_PI%lambda(         j,i)            ) + ((1._dp - w_tot( j,i)) *     climate%GCM_LGM%lambda(         j,i)            )
      Hs_ref_GCM(       j,i) =      (w_tot( j,i) *     climate%GCM_PI%Hs(             j,i)            ) + ((1._dp - w_tot( j,i)) *     climate%GCM_LGM%Hs(             j,i)            )   ! Berends et al., 2018 - Eq. 8
      T_ref_GCM(      :,j,i) = T_pot_ref_GCM( :,j,i) + lambda_GCM( j,i) * Hs_ref_GCM( j,i)
      P_ref_GCM(      :,j,i) = EXP( (w_tot( j,i) * LOG(climate%GCM_PI%Precip_corr(  :,j,i) + P_offset)) + ((1._dp - w_tot( j,i)) * LOG(climate%GCM_LGM%Precip_corr(  :,j,i) + P_offset)))  ! Berends et al., 2018 - Eq. 7
      
      
    END DO
    END DO
    CALL sync
    
    ! Downscale precipitation from the coarse-resolution reference GCM orography to the fine-resolution ice-model orography
    
    IF (region_name == 'NAM' .OR. region_name == 'EAS') THEN
      ! Use the Roe&Lindzen precipitation model to do this; Berends et al., 2018, Eqs. A3-A7
    
      ! Get the surface slopes of the coarse-resolution reference GCM orography
      CALL ddx_Aa_to_Aa_2D( grid, Hs_ref_GCM, dHs_dx_ref_GCM)
      CALL ddy_Aa_to_Aa_2D( grid, Hs_ref_GCM, dHs_dy_ref_GCM)
      
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
      DO m = 1, 12
      
        ! Calculate RL precipitation for the matrix-interpolated GCM reference state
        CALL precipitation_model_Roe( T_ref_GCM(           m,j,i), dHs_dx_ref_GCM( j,i), dHs_dy_ref_GCM( j,i), climate%PD_obs%Wind_LR( m,j,i), climate%PD_obs%Wind_DU( m,j,i), P_RL_ref_GCM( m,j,i))
        
        ! Calculate RL precipitation for the actual ice model state
        CALL precipitation_model_Roe( climate%applied%T2m( m,j,i), ice%dHs_dx_Aa(  j,i), ice%dHs_dy_Aa(  j,i), climate%PD_obs%Wind_LR( m,j,i), climate%PD_obs%Wind_DU( m,j,i), P_RL_mod(     m,j,i))
        
        ! Ratio between those two
        dP_RL( m,j,i) = MAX(0.01_dp, MIN( 2._dp, (P_RL_mod( m,j,i) + P_offset) / (P_RL_ref_GCM( m,j,i) + P_offset) ))
        
        ! Applied model precipitation = (matrix-interpolated GCM reference precipitation) * RL ratio
        climate%applied%Precip( m,j,i) = P_ref_GCM( m,j,i) * dP_RL( m,j,i)
        
      END DO
      END DO
      END DO
      CALL sync
    
    ELSEIF (region_name == 'GRL' .OR. region_name == 'ANT') THEN
      ! Use a simpler temperature-based correction; Berends et al., 2018, Eq. 14
      
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
      DO m = 1, 12
      
        CALL adjust_precipitation_for_temperature_change( T_ref_GCM( m,j,i), climate%applied%T2m( m,j,i), P_ref_GCM( m,j,i), climate%applied%Precip( m,j,i))
        
      END DO
      END DO
      END DO
      CALL sync
      
    END IF
   
    ! Clean up after yourself
    CALL deallocate_shared( ww_ice)
    CALL deallocate_shared( ww_tot)
    CALL deallocate_shared( wT_pot_ref_GCM)
    CALL deallocate_shared( wlambda_GCM)
    CALL deallocate_shared( wT_ref_GCM)
    CALL deallocate_shared( wP_ref_GCM)
    CALL deallocate_shared( wHs_ref_GCM)
    CALL deallocate_shared( wdHs_dx_ref_GCM)
    CALL deallocate_shared( wdHs_dy_ref_GCM)
    CALL deallocate_shared( wP_RL_ref_GCM)
    CALL deallocate_shared( wP_RL_mod)
    CALL deallocate_shared( wdP_RL)
    
  END SUBROUTINE run_climate_model_matrix_PI_LGM_precipitation
  
  ! The Roe & Lindzen precipitation model
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
  
  ! Adjust precipitation only based on temperature change (Jouzel & Merlivat, 1984, and also Huybrects, 1992), used in ANT and GRL
  SUBROUTINE adjust_precipitation_for_temperature_change( T_old, T_new, P_old, P_new)
    ! Precipitation model of Roe (J. Glac, 2002), integration from Roe and Lindzen (J. Clim. 2001)
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: T_old                ! old 2-m air temperature [K]
    REAL(dp),                            INTENT(IN)    :: T_new                ! new 2-m air temperature [K]
    REAL(dp),                            INTENT(IN)    :: P_old                ! old total monthly precipitation [m]
    REAL(dp),                            INTENT(OUT)   :: P_new                ! new total monthly precipitation [m]
    
    ! Local variables:
    REAL(dp)                                           :: T_inv_old, T_inv_new
    
    ! Adjust monthly temperature field and calculate monthly free atmospheric temperature
    ! Method of Jouzel and Merlivat (1984), see equation (4.82) in Huybrechts (1992): 
    T_inv_old = 88.9_dp + 0.67_dp * T_old
    T_inv_new = 88.9_dp + 0.67_dp * T_new

    ! Adjust precipation as function of Free atmospheric temperature (Tinv)
    P_new = P_old * 1.04_dp**(T_inv_new - T_inv_old)
    
  END SUBROUTINE adjust_precipitation_for_temperature_change
  
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
      
        dT_lapse = ice%Hs_Aa(j,i) * lambda
          
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
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        ! Entirely parameterised climate, no need to read anything here
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
    
    ! The differenct GCM snapshots 
    IF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      ! This choice of forcing doesn't use any GCM data
      RETURN
    ELSEIF (C%choice_forcing_method == 'CO2_direct' .OR. C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! These two choices use the climate matrix
      
      IF (C%choice_climate_matrix == 'PI_LGM') THEN
        CALL initialise_snapshot( matrix%GCM_PI,  name='HadCM3_PI',  nc_filename=C%filename_GCM_snapshot_PI,  CO2 = 280._dp, orbit_time =       0._dp)
        CALL initialise_snapshot( matrix%GCM_LGM, name='HadCM3_LGM', nc_filename=C%filename_GCM_snapshot_LGM, CO2 = 190._dp, orbit_time = -120000._dp)
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
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    
    PD_obs%name = name 
    
    ! Allocate memory for general forcing info (not relevant for PD but required for compatibility with GCM snapshots)
    CALL allocate_shared_dp_0D( PD_obs%CO2,        PD_obs%wCO2       )
    CALL allocate_shared_dp_0D( PD_obs%orbit_time, PD_obs%worbit_time)
    CALL allocate_shared_dp_0D( PD_obs%orbit_ecc,  PD_obs%worbit_ecc )
    CALL allocate_shared_dp_0D( PD_obs%orbit_obl,  PD_obs%worbit_obl )
    CALL allocate_shared_dp_0D( PD_obs%orbit_pre,  PD_obs%worbit_pre )
        
    PD_obs%netcdf%filename   = C%filename_PD_obs_climate 
    
    ! Determine size of Cartesian grid, so that memory can be allocated accordingly.
    CALL allocate_shared_int_0D( PD_obs%nlon, PD_obs%wnlon)
    CALL allocate_shared_int_0D( PD_obs%nlat, PD_obs%wnlat)
    
    ! For the benchmark experiments, use dummy input data.
    ! For realistic experiments, read the provided input file.
    IF (par%master) THEN
    
      IF (C%do_benchmark_experiment) THEN
        IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
            C%choice_benchmark_experiment == 'Halfar' .OR. &
            C%choice_benchmark_experiment == 'Bueler' .OR. &
            C%choice_benchmark_experiment == 'MISMIP_mod'.OR. &
            C%choice_benchmark_experiment == 'mesh_generation_test'.OR. &
            C%choice_benchmark_experiment == 'SSA_icestream') THEN
      
          PD_obs%nlon = 10
          PD_obs%nlat = 10
          
        ELSE
          IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_obs_data_fields!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
      ELSE ! IF (C%do_benchmark_experiment) THEN
      
        ! Read data from input file
        CALL inquire_PD_obs_data_file(PD_obs)
        
      END IF ! IF (C%do_benchmark_experiment) THEN
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Allocate memory  
    CALL allocate_shared_dp_1D( PD_obs%nlon,                  PD_obs%lon,         PD_obs%wlon        )
    CALL allocate_shared_dp_1D(              PD_obs%nlat,     PD_obs%lat,         PD_obs%wlat        )
    CALL allocate_shared_dp_2D( PD_obs%nlon, PD_obs%nlat,     PD_obs%Hi,          PD_obs%wHi         )
    CALL allocate_shared_dp_2D( PD_obs%nlon, PD_obs%nlat,     PD_obs%Hs,          PD_obs%wHs         )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%T2m,         PD_obs%wT2m        )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%Precip,      PD_obs%wPrecip     )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%Wind_WE,     PD_obs%wWind_WE    )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%Wind_SN,     PD_obs%wWind_SN    )
    
    ! For the benchmark experiments, use dummy input data.
    ! For realistic experiments, read the provided input file.
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod'.OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test'.OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        ! No need to read any climate data
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_obs_data_fields!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSE ! IF (C%do_benchmark_experiment) THEN
    
      ! Read data from input file
      IF (par%master) WRITE(0,*) '   Reading PD observed climate data from file ', TRIM(PD_obs%netcdf%filename), '...'
      IF (par%master) CALL read_PD_obs_data_file(PD_obs)
      CALL sync
      
    END IF ! IF (C%do_benchmark_experiment) THEN
      
    ! Determine process domains
    PD_obs%i1 = MAX(1,           FLOOR(REAL(PD_obs%nlon *  par%i      / par%n)) + 1)
    PD_obs%i2 = MIN(PD_obs%nlon, FLOOR(REAL(PD_obs%nlon * (par%i + 1) / par%n)))
    
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
    INTEGER                                       :: cerr, ierr
    
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
    ! At some point in the future, "orbit_time" will tell the model to read orbital parameters from
    ! an external file.
    
    ! Determine size of Cartesian grid, so that memory can be allocated accordingly.
    CALL allocate_shared_int_0D( snapshot%nlon, snapshot%wnlon)
    CALL allocate_shared_int_0D( snapshot%nlat, snapshot%wnlat)
    
    ! For the benchmark experiments, use dummy input data.
    ! For realistic experiments, read the provided input file.
    IF (par%master) THEN
    
      IF (C%do_benchmark_experiment) THEN
        IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
            C%choice_benchmark_experiment == 'Halfar' .OR. &
            C%choice_benchmark_experiment == 'Bueler' .OR. &
            C%choice_benchmark_experiment == 'MISMIP_mod'.OR. &
            C%choice_benchmark_experiment == 'mesh_generation_test'.OR. &
            C%choice_benchmark_experiment == 'SSA_icestream') THEN
      
          snapshot%nlon = 10
          snapshot%nlat = 10
          
        ELSE
          IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_obs_data_fields!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
      ELSE ! IF (C%do_benchmark_experiment) THEN
      
        ! Read data from input file
        CALL inquire_GCM_snapshot( snapshot)
        
      END IF ! IF (C%do_benchmark_experiment) THEN
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Allocate memory  
    CALL allocate_shared_dp_1D( snapshot%nlon,                    snapshot%lon,         snapshot%wlon        )
    CALL allocate_shared_dp_1D(                snapshot%nlat,     snapshot%lat,         snapshot%wlat        )
    CALL allocate_shared_dp_2D( snapshot%nlon, snapshot%nlat,     snapshot%Hi,          snapshot%wHi         )
    CALL allocate_shared_dp_2D( snapshot%nlon, snapshot%nlat,     snapshot%Hs,          snapshot%wHs         )
    CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, 12, snapshot%T2m,         snapshot%wT2m        )
    CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, 12, snapshot%Precip,      snapshot%wPrecip     )
    CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, 12, snapshot%Wind_WE,     snapshot%wWind_WE    )
    CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, 12, snapshot%Wind_SN,     snapshot%wWind_SN    )
    
    ! For the benchmark experiments, use dummy input data.
    ! For realistic experiments, read the provided input file.
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod'.OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test'.OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        ! No need to read any climate data
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_obs_data_fields!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSE ! IF (C%do_benchmark_experiment) THEN
    
      ! Read data from input file
      IF (par%master) WRITE(0,*) '   Reading GCM snapshot ', TRIM(snapshot%name), ' from file ', TRIM(snapshot%netcdf%filename), '...'
      IF (par%master) CALL read_GCM_snapshot( snapshot)
      CALL sync
      
    END IF ! IF (C%do_benchmark_experiment) THEN
      
    ! Determine process domains
    CALL partition_list( snapshot%nlon, par%i, par%n, snapshot%i1, snapshot%i2)
    
  END SUBROUTINE initialise_snapshot
  
  ! Initialising the region-specific climate model, containing all the subclimates
  ! (PD observations, GCM snapshots and the applied climate) on the model grid
  SUBROUTINE initialise_climate_model( grid, climate, matrix, region_name)
    ! Allocate shared memory for the regional climate models, containing the PD observed,
    ! GCM snapshots and applied climates as "subclimates"
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    TYPE(type_climate_matrix),           INTENT(IN)    :: matrix
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: i,j
    
    IF (par%master) WRITE (0,*) '  Initialising climate model...'
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod') THEN
          
        ! Entirely parameterised climate
        CALL initialise_subclimate( grid, climate%applied, 'applied')
        RETURN
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_climate_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Initialise data structures for the regional ERA40 climate and the final applied climate
    CALL initialise_subclimate( grid, climate%PD_obs,   'ERA40')
    CALL initialise_subclimate( grid, climate%applied,  'applied')
    
    ! Map these subclimates from global grid to model grid
    CALL map_subclimate_to_grid( grid,  matrix%PD_obs,  climate%PD_obs )
    
    ! Get PD_obs potential temperature
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      climate%PD_obs%T2m_pot( :,j,i) = climate%PD_obs%T2m( :,j,i) + climate%PD_obs%Hs( j,i) * 0.008_dp
    END DO
    END DO
    CALL sync
    
    ! The differenct GCM snapshots 
    IF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      ! This choice of forcing doesn't use any GCM data
      RETURN
    ELSEIF (C%choice_forcing_method == 'CO2_direct' .OR. C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! These two choices use the climate matrix
      
      IF (C%choice_climate_matrix == 'PI_LGM') THEN
      
        ! Initialise data structures for the GCM snapshots
        CALL initialise_subclimate( grid, climate%GCM_PI,   'HadCM3_PI')
        CALL initialise_subclimate( grid, climate%GCM_LGM,  'HadCM3_LGM')
        
        ! Map these subclimates from global grid to model grid
        CALL map_subclimate_to_grid( grid,  matrix%GCM_PI,  climate%GCM_PI )
        CALL map_subclimate_to_grid( grid,  matrix%GCM_LGM, climate%GCM_LGM)
        
        ! Calculate spatially variable lapse rate
        climate%GCM_PI%lambda = -0.008_dp
        IF     (region_name == 'NAM' .OR. region_name == 'EAS') THEN
          CALL initialise_subclimate_spatially_variable_lapserate( grid, climate%GCM_LGM, climate%GCM_PI)
        ELSEIF (region_name == 'GLR' .OR. region_name == 'ANT') THEN
          climate%GCM_LGM%lambda = -0.008_dp
        END IF
    
        ! Get bias-corrected climate for the GCM snapshots
        CALL correct_snapshot_PD_bias( grid, climate%GCM_PI,  climate%GCM_PI, climate%PD_obs)
        CALL correct_snapshot_PD_bias( grid, climate%GCM_LGM, climate%GCM_PI, climate%PD_obs)
        
        ! Get reference absorbed insolation for the GCM snapshots
        CALL initialise_subclimate_absorbed_insolation( grid, climate%GCM_PI , region_name)
        CALL initialise_subclimate_absorbed_insolation( grid, climate%GCM_LGM, region_name)
        
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
    
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: subclimate
    CHARACTER(LEN=*),                    INTENT(IN)    :: name
    
    subclimate%name = name
    
    CALL allocate_shared_dp_0D(                       subclimate%CO2,          subclimate%wCO2         )
    CALL allocate_shared_dp_0D(                       subclimate%orbit_time,   subclimate%worbit_time  )
    CALL allocate_shared_dp_0D(                       subclimate%orbit_ecc,    subclimate%worbit_ecc   )
    CALL allocate_shared_dp_0D(                       subclimate%orbit_obl,    subclimate%worbit_obl   )
    CALL allocate_shared_dp_0D(                       subclimate%orbit_pre,    subclimate%worbit_pre   )
    
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%T2m,          subclimate%wT2m         )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Precip,       subclimate%wPrecip      )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, subclimate%Hi,           subclimate%wHi          )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, subclimate%Hs,           subclimate%wHs          )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Wind_WE,      subclimate%wWind_WE     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Wind_SN,      subclimate%wWind_SN     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Wind_LR,      subclimate%wWind_LR     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Wind_DU,      subclimate%wWind_DU     )
    
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%T2m_pot,      subclimate%wT2m_pot     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%T2m_pot_corr, subclimate%wT2m_pot_corr)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Precip_corr,  subclimate%wPrecip_corr )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, subclimate%Hs_corr,      subclimate%wHs_corr     )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, subclimate%lambda,       subclimate%wlambda      )
    
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Q_TOA,        subclimate%wQ_TOA       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, subclimate%Albedo,       subclimate%wAlbedo      )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, subclimate%I_abs,        subclimate%wI_abs       )
    
    CALL allocate_shared_dp_0D(                       subclimate%T_ocean_mean, subclimate%wT_ocean_mean)
    
  END SUBROUTINE initialise_subclimate
  SUBROUTINE initialise_subclimate_absorbed_insolation( grid, snapshot, region_name)
    ! Calculate the yearly absorbed insolation for this (regional) GCM snapshot, to be used in the matrix interpolation
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: snapshot
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: i,j,m
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  Q_TOA0,  Q_TOA1
    INTEGER                                            :: wQ_TOA0, wQ_TOA1
    INTEGER                                            :: ti0, ti1
    REAL(dp)                                           :: ins_t0, ins_t1
    REAL(dp)                                           :: wt0, wt1
    INTEGER                                            :: ilat_l, ilat_u
    REAL(dp)                                           :: wlat_l, wlat_u
    
    TYPE(type_ice_model)                               :: ice_dummy
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
      DO WHILE (forcing%ins_time(ti1) < snapshot%orbit_time)
        ti1 = ti1 + 1
      END DO
      ti0 = ti1 - 1
      
      ins_t0 = forcing%ins_time(ti0)
      ins_t1 = forcing%ins_time(ti1)
    ELSE
      WRITE(0,*) '  ERROR - orbit_time ', snapshot%orbit_time, ' for snapshot ', TRIM(snapshot%name), ' outside of range of insolation solution file "', TRIM(forcing%netcdf%filename), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
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
    
    ! Create a temporary "dummy" ice & SMB data structure, so we can run the SMB model and determine the reference albedo field
    CALL allocate_shared_int_2D(    grid%ny, grid%nx, ice_dummy%mask_ocean_Aa   , ice_dummy%wmask_ocean_Aa   )
    CALL allocate_shared_int_2D(    grid%ny, grid%nx, ice_dummy%mask_ice_Aa     , ice_dummy%wmask_ice_Aa     )
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
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB_dummy%SMB             , SMB_dummy%wSMB             )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB_dummy%SMB_year        , SMB_dummy%wSMB_year        )
    
    ! Tuning parameters
    CALL allocate_shared_dp_0D( SMB_dummy%C_abl_constant, SMB_dummy%wC_abl_constant)
    CALL allocate_shared_dp_0D( SMB_dummy%C_abl_Ts,       SMB_dummy%wC_abl_Ts      )
    CALL allocate_shared_dp_0D( SMB_dummy%C_abl_Q,        SMB_dummy%wC_abl_Q       )
    CALL allocate_shared_dp_0D( SMB_dummy%C_refr,         SMB_dummy%wC_refr        )
    
    IF (par%master) THEN
      IF     (region_name == 'NAM') THEN
        SMB_dummy%C_abl_constant = C%C_abl_constant_NAM
        SMB_dummy%C_abl_Ts       = C%C_abl_Ts_NAM
        SMB_dummy%C_abl_Q        = C%C_abl_Q_NAM
        SMB_dummy%C_refr         = C%C_refr_NAM
      ELSEIF (region_name == 'EAS') THEN
        SMB_dummy%C_abl_constant = C%C_abl_constant_EAS
        SMB_dummy%C_abl_Ts       = C%C_abl_Ts_EAS
        SMB_dummy%C_abl_Q        = C%C_abl_Q_EAS
        SMB_dummy%C_refr         = C%C_refr_EAS
      ELSEIF (region_name == 'GRL') THEN
        SMB_dummy%C_abl_constant = C%C_abl_constant_GRL
        SMB_dummy%C_abl_Ts       = C%C_abl_Ts_GRL
        SMB_dummy%C_abl_Q        = C%C_abl_Q_GRL
        SMB_dummy%C_refr         = C%C_refr_GRL
      ELSEIF (region_name == 'ANT') THEN
        SMB_dummy%C_abl_constant = C%C_abl_constant_ANT
        SMB_dummy%C_abl_Ts       = C%C_abl_Ts_ANT
        SMB_dummy%C_abl_Q        = C%C_abl_Q_ANT
        SMB_dummy%C_refr         = C%C_refr_ANT
      END IF
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Fill in the dummy ice and ocean mask (the only two variables the SMB model needs) with the ice thickness and surface elevation
    ! data read from the GCM snapshot netcdf file
    ice_dummy%mask_ocean_Aa( :,grid%i1:grid%i2) = 0
    ice_dummy%mask_ice_Aa(   :,grid%i1:grid%i2) = 0
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (snapshot%Hs( j,i) == MINVAL(snapshot%Hs)) ice_dummy%mask_ocean_Aa( j,i) = 1
      IF (snapshot%Hi( j,i) > 0._dp) ice_dummy%mask_ice_Aa( j,i) = 1
    END DO
    END DO
    CALL sync
    
    ! Run the SMB model for 10 years for this particular snapshot
    ! (experimentally determined to be long enough to converge)
    DO i = 1, 10
      CALL run_SMB_model( grid, ice_dummy, snapshot, 0._dp, SMB_dummy)
    END DO
    
    ! Copy the resulting albedo to the snapshot
    snapshot%Albedo( :,:,grid%i1:grid%i2) = SMB_dummy%Albedo( :,:,grid%i1:grid%i2)
    
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
    CALL deallocate_shared( ice_dummy%wmask_ocean_Aa)
    CALL deallocate_shared( ice_dummy%wmask_ice_Aa)
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
    CALL deallocate_shared( SMB_dummy%wSMB)
    CALL deallocate_shared( SMB_dummy%wSMB_year)
    CALL deallocate_shared( SMB_dummy%wC_abl_constant)
    CALL deallocate_shared( SMB_dummy%wC_abl_Ts)
    CALL deallocate_shared( SMB_dummy%wC_abl_Q)
    CALL deallocate_shared( SMB_dummy%wC_refr)
    
  END SUBROUTINE initialise_subclimate_absorbed_insolation
  SUBROUTINE initialise_subclimate_spatially_variable_lapserate( grid, snapshot, snapshot_PI)
    ! Calculate the spatially variable lapse-rate (for non-PI GCM snapshots; see Berends et al., 2018)
    ! Only meaningful for snapshots where there is ice (LGM, M2_Medium, M2_Large),
    ! and only intended for North America and Eurasia
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: snapshot
    TYPE(type_subclimate_region),        INTENT(IN)    :: snapshot_PI
    
    ! Local variables:
    INTEGER                                            :: ierr
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: dT_mean_nonice
    INTEGER                                            :: n_nonice, n_ice
    REAL(dp)                                           :: lambda_mean_ice
    
    REAL(dp), PARAMETER                                :: lambda_min = 0.002_dp
    REAL(dp), PARAMETER                                :: lambda_max = 0.05_dp
      
    ! Calculate the regional average temperature change outside of the ice sheet.
    ! ===========================================================================
    
    dT_mean_nonice = 0._dp
    n_nonice              = 0
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO m = 1, 12
      IF (snapshot%Hi( j,i) == 0._dp) THEN
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
    
      snapshot%lambda( j,i) = 0._dp
      
      IF (snapshot%Hi( j,i) > 100._dp .AND. snapshot%Hs( j,i) > snapshot_PI%Hs( j,i)) THEN
      
        DO m = 1, 12
          snapshot%lambda( j,i) = snapshot%lambda( j,i) + 1/12._dp * MAX(lambda_min, MIN(lambda_max, &                        ! Berends et al., 2018 - Eq. 10
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
      IF (.NOT. (snapshot%Hi( j,i) > 100._dp .AND. snapshot%Hs( j,i) > snapshot_PI%Hs( j,i))) THEN
        snapshot%lambda( j,i) = lambda_mean_ice
      END IF
    END DO
    END DO
    CALL sync
    
    ! Normalise the entire region to a mean lapse rate of 8 K /km
    snapshot%lambda( :,grid%i1:grid%i2) = snapshot%lambda( :,grid%i1:grid%i2) * (-0.008_dp / lambda_mean_ice)
    
    ! Smooth the lapse rate field with a 160 km Gaussian filter
    CALL smooth_Gaussian_2D( grid, snapshot%lambda, 4 * CEILING(40000._dp / grid%dx))
    
  END SUBROUTINE initialise_subclimate_spatially_variable_lapserate
  SUBROUTINE correct_snapshot_PD_bias( grid, snapshot, snapshot_PI, PD_obs)
    ! Get bias-corrected climate for the GCM snapshots
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: snapshot
    TYPE(type_subclimate_region),        INTENT(IN)    :: snapshot_PI
    TYPE(type_subclimate_region),        INTENT(IN)    :: PD_obs
    
    ! Local variables:
    INTEGER                                       :: i,j
    REAL(dp), PARAMETER                           :: P_offset = 1E-4_dp       ! Normalisation term in precipitation anomaly to avoid divide-by-nearly-zero
    
    ! Calculate potential (= elevation-corrected) temperature
    ! =======================================================
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      snapshot%T2m_pot( :,j,i) = snapshot%T2m( :,j,i) - snapshot%Hs( j,i) * snapshot%lambda( j,i)
    END DO
    END DO
    CALL sync
    
    ! Calculate bias-corrected temperature
    ! ====================================
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      snapshot%T2m_pot_corr( :,j,i) = snapshot%T2m_pot( :,j,i) - ( snapshot_PI%T2m_pot( :,j,i)             -  PD_obs%T2m_pot( :,j,i)            )
      snapshot%Precip_corr(  :,j,i) = snapshot%Precip(  :,j,i) / ((snapshot_PI%Precip(  :,j,i) + P_offset) / (PD_obs%Precip(  :,j,i) + P_offset))
      snapshot%Hs_corr(        j,i) = snapshot%Hs(        j,i) - ( snapshot_PI%Hs(        j,i)             -  PD_obs%Hs(        j,i)            )
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE correct_snapshot_PD_bias
  
  ! Map a global subclimate from the matrix (PD observed or GCM snapshot) to a region grid
  SUBROUTINE map_subclimate_to_grid( grid,  cglob, creg)
    ! Map data from a global "subclimate" (PD observed or GCM snapshot) to the grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_subclimate_global),        INTENT(IN)    :: cglob  ! Global climate
    TYPE(type_subclimate_region),        INTENT(INOUT) :: creg   ! grid   climate
      
    CALL map_glob_to_grid_2D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%Hi,          creg%Hi             )
    CALL map_glob_to_grid_2D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%Hs,          creg%Hs             )
    CALL map_glob_to_grid_3D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%T2m,         creg%T2m    , 12    )
    CALL map_glob_to_grid_3D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%Precip,      creg%Precip , 12    )
    CALL map_glob_to_grid_3D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%Wind_WE,     creg%Wind_WE, 12    )
    CALL map_glob_to_grid_3D( cglob%nlat, cglob%nlon, cglob%lat, cglob%lon, grid, cglob%Wind_SN,     creg%Wind_SN, 12    )
    
    ! Copy snapshot metadata
    creg%CO2        = cglob%CO2
    creg%orbit_time = cglob%orbit_time
    creg%orbit_ecc  = cglob%orbit_ecc
    creg%orbit_obl  = cglob%orbit_obl
    creg%orbit_pre  = cglob%orbit_pre
    
    ! Rotate zonal/meridional wind to x,y wind
    CALL rotate_wind_to_model_grid( grid, creg%wind_WE, creg%wind_SN, creg%wind_LR, creg%wind_DU)
  
  END SUBROUTINE map_subclimate_to_grid
  SUBROUTINE map_glob_to_grid_2D( nlat, nlon, lat, lon, grid, d_glob, d_grid)
    ! Map a data field from a global lat-lon grid to the regional square grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    INTEGER,                         INTENT(IN)  :: nlat
    INTEGER,                         INTENT(IN)  :: nlon
    REAL(dp), DIMENSION(nlat),       INTENT(IN)  :: lat
    REAL(dp), DIMENSION(nlon),       INTENT(IN)  :: lon
    TYPE(type_grid),                 INTENT(IN)  :: grid
    REAL(dp), DIMENSION(:,:  ),      INTENT(IN)  :: d_glob
    REAL(dp), DIMENSION(:,:  ),      INTENT(OUT) :: d_grid
    
    ! Local variables:
    INTEGER                                                :: i, j, il, iu, jl, ju
    REAL(dp)                                               :: wil, wiu, wjl, wju
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Find enveloping lat-lon indices
      il  = MAX(1,MIN(nlon-1, 1 + FLOOR((grid%lon(j,i)-MINVAL(lon)) / (lon(2)-lon(1)))))
      iu  = il+1        
      wil = (lon(iu) - grid%lon(j,i))/(lon(2)-lon(1))
      wiu = 1-wil
      
      ! Exception for pixels near the zero meridian
      IF (grid%lon(j,i) < MINVAL(lon)) THEN
        il = nlon
        iu = 1      
        wil = (lon(iu) - grid%lon(j,i))/(lon(2)-lon(1))
        wiu = 1-wil
      ELSEIF (grid%lon(j,i) > MAXVAL(lon)) THEN
        il = nlon
        iu = 1
        wiu = (grid%lon(j,i) - lon(il))/(lon(2)-lon(1))
        wil = 1-wiu
      END IF
          
      jl  = MAX(1,MIN(nlat-1, 1 + FLOOR((grid%lat(j,i)-MINVAL(lat)) / (lat(2)-lat(1)))))
      ju  = jl+1        
      wjl = (lat(ju) - grid%lat(j,i))/(lat(2)-lat(1))
      wju = 1-wjl
      
      ! Interpolate data
      d_grid( j,i) = (d_glob( il,jl) * wil * wjl) + &
                     (d_glob( il,ju) * wil * wju) + &
                     (d_glob( iu,jl) * wiu * wjl) + &
                     (d_glob( iu,ju) * wiu * wju)
    
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_glob_to_grid_2D
  SUBROUTINE map_glob_to_grid_3D( nlat, nlon, lat, lon, grid, d_glob, d_grid, nz)
    ! Map a data field from a global lat-lon grid to the regional square grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    INTEGER,                         INTENT(IN)  :: nlat
    INTEGER,                         INTENT(IN)  :: nlon
    REAL(dp), DIMENSION(nlat),       INTENT(IN)  :: lat
    REAL(dp), DIMENSION(nlon),       INTENT(IN)  :: lon
    TYPE(type_grid),                 INTENT(IN)  :: grid
    REAL(dp), DIMENSION(:,:,:),      INTENT(IN)  :: d_glob
    REAL(dp), DIMENSION(:,:,:),      INTENT(OUT) :: d_grid
    INTEGER,                         INTENT(IN)  :: nz
    
    ! Local variables:
    INTEGER                                                :: i, j, il, iu, jl, ju, k
    REAL(dp)                                               :: wil, wiu, wjl, wju
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Find enveloping lat-lon indices
      il  = MAX(1,MIN(nlon-1, 1 + FLOOR((grid%lon(j,i)-MINVAL(lon)) / (lon(2)-lon(1)))))
      iu  = il+1        
      wil = (lon(iu) - grid%lon(j,i))/(lon(2)-lon(1))
      wiu = 1-wil
      
      ! Exception for pixels near the zero meridian
      IF (grid%lon(j,i) < MINVAL(lon)) THEN
        il = nlon
        iu = 1      
        wil = (lon(iu) - grid%lon(j,i))/(lon(2)-lon(1))
        wiu = 1-wil
      ELSEIF (grid%lon(j,i) > MAXVAL(lon)) THEN
        il = nlon
        iu = 1
        wiu = (grid%lon(j,i) - lon(il))/(lon(2)-lon(1))
        wil = 1-wiu
      END IF
          
      jl  = MAX(1,MIN(nlat-1, 1 + FLOOR((grid%lat(j,i)-MINVAL(lat)) / (lat(2)-lat(1)))))
      ju  = jl+1        
      wjl = (lat(ju) - grid%lat(j,i))/(lat(2)-lat(1))
      wju = 1-wjl
      
      ! Interpolate data
      DO k = 1, nz
        d_grid( k,j,i) = (d_glob( il,jl,k) * wil * wjl) + &
                         (d_glob( il,ju,k) * wil * wju) + &
                         (d_glob( iu,jl,k) * wiu * wjl) + &
                         (d_glob( iu,ju,k) * wiu * wju)
      END DO
    
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_glob_to_grid_3D
  
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
