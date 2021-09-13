MODULE SMB_module

  ! Contains all the routines for calculating the surface mass balance for the current climate.

  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_subclimate_region, type_init_data_fields, type_SMB_model
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  USE forcing_module,                  ONLY: forcing
  USE parameters_module,               ONLY: T0, L_fusion, sec_per_year, pi, ice_density

  IMPLICIT NONE
  
  REAL(dp), PARAMETER :: albedo_water        = 0.1_dp
  REAL(dp), PARAMETER :: albedo_soil         = 0.2_dp
  REAL(dp), PARAMETER :: albedo_ice          = 0.5_dp
  REAL(dp), PARAMETER :: albedo_snow         = 0.85_dp
  REAL(dp), PARAMETER :: initial_snow_depth  = 0.1_dp
    
CONTAINS

  ! Run the SMB model on the region grid
  SUBROUTINE run_SMB_model(            grid, ice, climate, time, SMB, mask_noice)
    ! Run the IMAU-ITM SMB model. Old version, exactly as it was in ANICE2.1 (so with the "wrong" refreezing)
    
    ! NOTE: all the SMB components and the total are in meters of water equivalent
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice 
    TYPE(type_subclimate_region),        INTENT(IN)    :: climate
    REAL(dp),                            INTENT(IN)    :: time
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_noice
    
    ! Local variables:
    INTEGER                                            :: i,j,m
    INTEGER                                            :: mprev
    REAL(dp)                                           :: snowfrac, liquid_water, sup_imp_wat
    
    ! Check if we need to apply any special benchmark experiment SMB
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6') THEN
          
        CALL EISMINT_SMB( grid, time, SMB)
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
        SMB%SMB_year(   :,grid%i1:grid%i2) = 0._dp
        SMB%SMB(      :,:,grid%i1:grid%i2) = 0._dp
        CALL sync
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'Bueler') THEN
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          SMB%SMB_year( j,i) = Bueler_solution_MB( grid%x(i), grid%y(j), time)
          SMB%SMB(    :,j,i) = SMB%SMB_year( j,i) / 12._dp
        END DO
        END DO
        CALL sync
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        SMB%SMB_year(   :,grid%i1:grid%i2) = 0.3_dp
        SMB%SMB(      :,:,grid%i1:grid%i2) = 0.3_dp / 12._dp
        CALL sync
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_SMB_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
        
    ! Calculate SMB components with IMAU_ITM
    ! ======================================
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! "Background" albedo (= surface without any firn, so either ice, soil, or water)
      SMB%AlbedoSurf( j,i) = albedo_soil
      IF ((ice%mask_ocean_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 0) .OR. mask_noice( j,i) == 1) SMB%AlbedoSurf( j,i) = albedo_water
      IF (ice%mask_ice_a(    j,i) == 1) SMB%AlbedoSurf( j,i) = albedo_ice
    
      DO m = 1, 12  ! Month loop
        
        mprev = m - 1
        IF (mprev == 0) mprev = 12
        
        SMB%Albedo( m,j,i) = MIN(albedo_snow, MAX( SMB%AlbedoSurf( j,i), albedo_snow - (albedo_snow - SMB%AlbedoSurf( j,i))  * &
                             EXP(-15._dp * SMB%FirnDepth( mprev,j,i)) - 0.015_dp * SMB%MeltPreviousYear( j,i)))
        IF ((ice%mask_ocean_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 0) .OR. mask_noice( j,i) == 1) SMB%Albedo( m,j,i) = albedo_water
               
        ! Determine ablation as function af surface temperature and albedo/insolation
        ! according to Bintanja et al. (2002) 
    
        SMB%Melt( m,j,i) = MAX(0._dp, ( SMB%C_abl_Ts         * (climate%T2m( m,j,i) - T0) + &
                                        SMB%C_abl_Q          * (1.0_dp - SMB%Albedo( m,j,i)) * climate%Q_TOA( m,j,i) - &
                                        SMB%C_abl_constant)  * sec_per_year / (L_fusion * 1000._dp * 12._dp))
                
        ! Determine accumulation with snow/rain fraction from Ohmura et al. (1999),
        ! liquid water content (rain and melt water) and snowdepth
        snowfrac = MAX(0._dp, MIN(1._dp, 0.5_dp   * (1 - ATAN((climate%T2m( m,j,i) - T0) / 3.5_dp)  / 1.25664_dp)))
    
        SMB%Snowfall( m,j,i) = climate%Precip( m,j,i) *          snowfrac
        SMB%Rainfall( m,j,i) = climate%Precip( m,j,i) * (1._dp - snowfrac)
    
        ! Refreezing, according to Janssens & Huybrechts, 2000)
        ! The refreezing (=effective retention) is the minimum value of the amount of super imposed 
        ! water and the available liquid water, with a maximum value of the total precipitation.
        ! (see also Huybrechts & de Wolde, 1999)
        sup_imp_wat  = 0.012_dp * MAX(0._dp, T0 - climate%T2m( m,j,i))
        liquid_water = SMB%Rainfall( m,j,i) + SMB%Melt( m,j,i)
        SMB%Refreezing( m,j,i) = MIN( MIN( sup_imp_wat, liquid_water), climate%Precip( m,j,i))
        IF (ice%mask_ice_a( j,i) == 0 .OR. mask_noice( j,i) == 1) SMB%Refreezing( m,j,i) = 0._dp
        
        ! Calculate runoff and total SMB
        SMB%Runoff( m,j,i) = SMB%Melt(     m,j,i) + SMB%Rainfall(   m,j,i) - SMB%Refreezing( m,j,i)
        SMB%SMB(    m,j,i) = SMB%Snowfall( m,j,i) + SMB%Refreezing( m,j,i) - SMB%Melt(       m,j,i)
    
        ! Add this month's snow accumulation to next month's initial snow depth.
        SMB%AddedFirn( m,j,i) = SMB%Snowfall( m,j,i) - SMB%Melt( m,j,i)
        SMB%FirnDepth( m,j,i) = MIN(10._dp, MAX(0._dp, SMB%FirnDepth( mprev,j,i) + SMB%AddedFirn( m,j,i) ))
    
      END DO ! DO m = 1, 12
      
      ! Calculate total SMB for the entire year
      SMB%SMB_year( j,i) = SUM(SMB%SMB( :,j,i))
      
      ! Calculate total melt over this year, to be used for determining next year's albedo
      SMB%MeltPreviousYear( j,i) = SUM(SMB%Melt( :,j,i))
      
      ! Calculate yearly mean albedo (diagnostic only)
      SMB%Albedo_year( j,i) = SUM(SMB%Albedo( :,j,i)) / 12._dp
      
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( SMB%AlbedoSurf      , 'SMB%AlbedoSurf'      , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%Albedo          , 'SMB%Albedo'          , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%Melt            , 'SMB%Melt'            , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%Snowfall        , 'SMB%Snowfall'        , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%Rainfall        , 'SMB%Rainfall'        , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%Refreezing      , 'SMB%Refreezing'      , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%Runoff          , 'SMB%Runoff'          , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%SMB             , 'SMB%SMB'             , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%AddedFirn       , 'SMB%AddedFirn'       , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%FirnDepth       , 'SMB%FirnDepth'       , 'run_SMB_model')
    CALL check_for_NaN_dp_2D( SMB%SMB_year        , 'SMB%SMB_year'        , 'run_SMB_model')
    CALL check_for_NaN_dp_2D( SMB%MeltPreviousYear, 'SMB%MeltPreviousYear', 'run_SMB_model')
    CALL check_for_NaN_dp_2D( SMB%Albedo_year     , 'SMB%Albedo_year'     , 'run_SMB_model')
          
  END SUBROUTINE run_SMB_model
  SUBROUTINE run_SMB_model_refr_fixed( grid, ice, climate, time, SMB, mask_noice)
    ! Run the IMAU-ITM SMB model. Based on the one from ANICE.
    
    ! NOTE: all the SMB components and the total are in meters of water equivalent
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice 
    TYPE(type_subclimate_region),        INTENT(IN)    :: climate
    REAL(dp),                            INTENT(IN)    :: time
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_noice
    
    ! Local variables:
    INTEGER                                            :: i,j,m
    INTEGER                                            :: mprev
    REAL(dp)                                           :: snowfrac, liquid_water, sup_imp_wat
    
    ! Check if we need to apply any special benchmark experiment SMB
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6') THEN
          
        CALL EISMINT_SMB( grid, time, SMB)
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'SSA_icestream') THEN
        SMB%SMB_year(   :,grid%i1:grid%i2) = 0._dp
        SMB%SMB(      :,:,grid%i1:grid%i2) = 0._dp
        CALL sync
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'Bueler') THEN
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          SMB%SMB_year( j,i) = Bueler_solution_MB( grid%x(i), grid%y(j), time)
          SMB%SMB(    :,j,i) = SMB%SMB_year( j,i) / 12._dp
        END DO
        END DO
        CALL sync
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        SMB%SMB_year(   :,grid%i1:grid%i2) = 0.3_dp
        SMB%SMB(      :,:,grid%i1:grid%i2) = 0.3_dp / 12._dp
        CALL sync
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_SMB_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Calculate SMB components with IMAU_ITM
    ! ======================================
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! Background albedo
      SMB%AlbedoSurf( j,i) = albedo_soil
      IF ((ice%mask_ocean_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 0) .OR. mask_noice( j,i) == 1) SMB%AlbedoSurf( j,i) = albedo_water
      IF (ice%mask_ice_a(    j,i) == 1) SMB%AlbedoSurf( j,i) = albedo_ice
    
      DO m = 1, 12  ! Month loop
        
        mprev = m - 1
        IF (mprev==0) mprev = 12
        
        SMB%Albedo( m,j,i) = MIN(albedo_snow, MAX( SMB%AlbedoSurf( j,i), albedo_snow - (albedo_snow - SMB%AlbedoSurf( j,i))  * &
                             EXP(-15._dp * SMB%FirnDepth( mprev,j,i)) - 0.015_dp * SMB%MeltPreviousYear( j,i)))
        IF ((ice%mask_ocean_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 0) .OR. mask_noice( j,i) == 1) SMB%Albedo( m,j,i) = albedo_water
               
        ! Determine albation as function af surface temperature and albedo/insolation
        ! according to Bintanja et al. (2002) 
    
        SMB%Melt( m,j,i) = MAX(0._dp, ( SMB%C_abl_Ts         * (climate%T2m( m,j,i) - T0) + &
                                        SMB%C_abl_Q          * (1.0_dp - SMB%Albedo( m,j,i)) * climate%Q_TOA( m,j,i) - &
                                        SMB%C_abl_constant)  * sec_per_year / (L_fusion * 1000._dp * 12._dp))
                
        ! Determine accumulation with snow/rain fraction from Ohmura et al. (1999),
        ! liquid water content (rain and melt water) and snowdepth
    
        ! NOTE: commented version is the old ANICE version, supposedly based on "physics" (which we cant check), but 
        !       the new version was tuned to RACMO output and produced significantly better snow fractions...
    
  !      snowfrac = MAX(0._dp, MIN(1._dp, 0.5_dp   * (1 - ATAN((climate%T2m(vi,m) - T0) / 3.5_dp)  / 1.25664_dp)))
        snowfrac = MAX(0._dp, MIN(1._dp, 0.725_dp * (1 - ATAN((climate%T2m( m,j,i) - T0) / 5.95_dp) / 1.8566_dp)))
    
        SMB%Snowfall( m,j,i) = climate%Precip( m,j,i) *          snowfrac
        SMB%Rainfall( m,j,i) = climate%Precip( m,j,i) * (1._dp - snowfrac)
    
        ! Refreezing, according to Janssens & Huybrechts, 2000)
        ! The refreezing (=effective retention) is the minimum value of the amount of super imposed 
        ! water and the available liquid water, with a maximum value of the total precipitation.
        ! (see also Huybrechts & de Wolde, 1999)
    
        ! Add this month's snow accumulation to next month's initial snow depth.
        SMB%AddedFirn( m,j,i) = SMB%Snowfall( m,j,i) - SMB%Melt( m,j,i)
        SMB%FirnDepth( m,j,i) = MIN(10._dp, MAX(0._dp, SMB%FirnDepth( mprev,j,i) + SMB%AddedFirn( m,j,i) ))
    
      END DO ! DO m = 1, 12
    
      ! Calculate refrezzing for the whole year, divide equally over the 12 months, then calculate resulting runoff and SMB.
      ! This resolves the problem with refreezing, where liquid water is mostly available in summer
      ! but "refreezing potential" mostly in winter, and there is no proper meltwater retention.
      
      sup_imp_wat  = SMB%C_refr * MAX(0._dp, T0 - SUM(climate%T2m( :,j,i))/12._dp)
      liquid_water = SUM(SMB%Rainfall( :,j,i)) + SUM(SMB%Melt( :,j,i))
      
      SMB%Refreezing_year( j,i) = MIN( MIN( sup_imp_wat, liquid_water), SUM(climate%Precip( :,j,i)))
      IF (ice%mask_ice_a( j,i)==0) SMB%Refreezing_year( j,i) = 0._dp
  
      DO m = 1, 12
        SMB%Refreezing( m,j,i) = SMB%Refreezing_year( j,i) / 12._dp
        SMB%Runoff(     m,j,i) = SMB%Melt( m,j,i) + SMB%Rainfall( m,j,i) - SMB%Refreezing( m,j,i)
        SMB%SMB(        m,j,i) = SMB%Snowfall( m,j,i) + SMB%Refreezing( m,j,i) - SMB%Melt( m,j,i)
      END DO
      
      SMB%SMB_year( j,i) = SUM(SMB%SMB( :,j,i))
      
      ! Calculate total melt over this year, to be used for determining next year's albedo
      SMB%MeltPreviousYear( j,i) = SUM(SMB%Melt( :,j,i))
      
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( SMB%AlbedoSurf      , 'SMB%AlbedoSurf'      , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%Albedo          , 'SMB%Albedo'          , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%Melt            , 'SMB%Melt'            , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%Snowfall        , 'SMB%Snowfall'        , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%Rainfall        , 'SMB%Rainfall'        , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%Refreezing      , 'SMB%Refreezing'      , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%Runoff          , 'SMB%Runoff'          , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%SMB             , 'SMB%SMB'             , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%AddedFirn       , 'SMB%AddedFirn'       , 'run_SMB_model')
    CALL check_for_NaN_dp_3D( SMB%FirnDepth       , 'SMB%FirnDepth'       , 'run_SMB_model')
    CALL check_for_NaN_dp_2D( SMB%SMB_year        , 'SMB%SMB_year'        , 'run_SMB_model')
    CALL check_for_NaN_dp_2D( SMB%MeltPreviousYear, 'SMB%MeltPreviousYear', 'run_SMB_model')
    CALL check_for_NaN_dp_2D( SMB%Albedo_year     , 'SMB%Albedo_year'     , 'run_SMB_model')
          
  END SUBROUTINE run_SMB_model_refr_fixed
  
  ! The EISMINT SMB parameterisations
  SUBROUTINE EISMINT_SMB( grid, time, SMB)
    ! Run the IMAU-ITM SMB model. Based on the one from ANICE.
    
    ! NOTE: all the SMB components and the total are in meters of water equivalent
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    REAL(dp),                            INTENT(IN)    :: time
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    REAL(dp)                                           :: E               ! Radius of circle where accumulation is M_max
    REAL(dp)                                           :: dist            ! distance to centre of circle
    REAL(dp)                                           :: S_b             ! Gradient of accumulation-rate change with horizontal distance
    REAL(dp)                                           :: M_max           ! Maximum accumulation rate 
    
    ! Default EISMINT configuration
    E         = 450000._dp
    S_b       = 0.01_dp / 1000._dp 
    M_max     = 0.5_dp
    
    IF     (C%choice_benchmark_experiment == 'EISMINT_1') THEN ! Moving margin, steady state
      ! No changes
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_2') THEN ! Moving margin, 20 kyr
      IF (time < 0._dp) THEN
        ! No changes; first 120 kyr are initialised with EISMINT_1
      ELSE
        E         = 450000._dp + 100000._dp * SIN( 2._dp * pi * time / 20000._dp)
      END IF
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_3') THEN ! Moving margin, 40 kyr
      IF (time < 0._dp) THEN
        ! No changes; first 120 kyr are initialised with EISMINT_1
      ELSE
        E         = 450000._dp + 100000._dp * SIN( 2._dp * pi * time / 40000._dp)
      END IF
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_4') THEN ! Fixed margin, steady state
      M_max       = 0.3_dp       
      E           = 999000._dp
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_5') THEN ! Fixed margin, 20 kyr
      IF (time < 0._dp) THEN
        M_max     = 0.3_dp
        E         = 999000._dp 
      ELSE
        M_max     = 0.3_dp + 0.2_dp * SIN( 2._dp * pi * time / 20000._dp)
        E         = 999000._dp 
      END IF
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_6') THEN ! Fixed margin, 40 kyr
      IF (time < 0._dp) THEN
        M_max     = 0.3_dp
        E         = 999000._dp 
      ELSE
        M_max     = 0.3_dp + 0.2_dp * SIN( 2._dp * pi * time / 40000._dp)
        E         = 999000._dp 
      END IF
    END IF

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      dist = SQRT(grid%x(i)**2+grid%y(j)**2)
      SMB%SMB_year( j,i) = MIN( M_max, S_b * (E - dist))
      SMB%SMB(    :,j,i) = SMB%SMB_year( j,i) / 12._dp
    END DO
    END DO
    CALL sync
          
  END SUBROUTINE EISMINT_SMB
  FUNCTION Bueler_solution_MB( x, y, t) RESULT(M)
    ! Describes an ice-sheet at time t (in years) conforming to the Bueler solution
    ! with dome thickness H0 and margin radius R0 at t0, with a surface mass balance
    ! determined by lambda.
    
    ! Input variables
    REAL(dp), INTENT(IN) :: x       ! x coordinate [m]
    REAL(dp), INTENT(IN) :: y       ! y coordinate [m]
    REAL(dp), INTENT(IN) :: t       ! Time from t0 [years]
    
    ! Result
    REAL(dp)             :: M
    
    ! Local variables
    REAL(dp) :: A_flow, rho, g, n, alpha, beta, Gamma, f1, f2, t0, tp, f3, f4, H
    
    REAL(dp), PARAMETER :: H0     = 3000._dp    ! Ice dome thickness at t=0 [m]
    REAL(dp), PARAMETER :: R0     = 500000._dp  ! Ice margin radius  at t=0 [m]
    REAL(dp), PARAMETER :: lambda = 5.0_dp      ! Mass balance parameter
  
    A_flow  = 1E-16_dp
    rho     = 910._dp
    g       = 9.81_dp
    n       = 3._dp
    
    alpha = (2._dp - (n+1._dp)*lambda) / ((5._dp*n)+3._dp)
    beta  = (1._dp + ((2._dp*n)+1._dp)*lambda) / ((5._dp*n)+3._dp)
    Gamma = 2._dp/5._dp * (A_flow/sec_per_year) * (rho * g)**n
    
    f1 = ((2._dp*n)+1)/(n+1._dp)
    f2 = (R0**(n+1._dp))/(H0**((2._dp*n)+1._dp))
    t0 = (beta / Gamma) * (f1**n) * f2 
    
    !tp = (t * sec_per_year) + t0; % Acutal equation needs t in seconds from zero , but we want to supply t in years from t0
    tp = t * sec_per_year
    
    f1 = (tp / t0)**(-alpha)
    f2 = (tp / t0)**(-beta)
    f3 = SQRT( (x**2._dp) + (y**2._dp) )/R0
    f4 = MAX(0._dp, 1._dp - (f2*f3)**((n+1._dp)/n))
    H = H0 * f1 * f4**(n/((2._dp*n)+1._dp))
    
    M = (lambda / tp) * H * sec_per_year
  
  END FUNCTION Bueler_solution_MB
  
  ! Initialise the SMB model (allocating shared memory)
  SUBROUTINE initialise_SMB_model( grid, init, SMB, region_name)
    ! Allocate memory for the data fields of the SMB model.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_init_data_fields),         INTENT(IN)    :: init
    TYPE(type_SMB_model),                INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables
    INTEGER                                            :: i,j
    
    IF (par%master) WRITE (0,*) '  Initialising SMB model...'
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB%AlbedoSurf      , SMB%wAlbedoSurf      )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB%MeltPreviousYear, SMB%wMeltPreviousYear)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB%FirnDepth       , SMB%wFirnDepth       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB%Rainfall        , SMB%wRainfall        )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB%Snowfall        , SMB%wSnowfall        )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB%AddedFirn       , SMB%wAddedFirn       )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB%Melt            , SMB%wMelt            )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB%Refreezing      , SMB%wRefreezing      )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB%Refreezing_year , SMB%wRefreezing_year )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB%Runoff          , SMB%wRunoff          )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB%Albedo          , SMB%wAlbedo          )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB%Albedo_year     , SMB%wAlbedo_year     )
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, SMB%SMB             , SMB%wSMB             )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB%SMB_year        , SMB%wSMB_year        )
    
    ! DENK DROM
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB%Snow_prev       , SMB%wSnow_prev       )
    CALL allocate_shared_dp_2D(     grid%ny, grid%nx, SMB%abla_prev       , SMB%wabla_prev       )
    
    ! Tuning parameters
    CALL allocate_shared_dp_0D( SMB%C_abl_constant, SMB%wC_abl_constant)
    CALL allocate_shared_dp_0D( SMB%C_abl_Ts,       SMB%wC_abl_Ts      )
    CALL allocate_shared_dp_0D( SMB%C_abl_Q,        SMB%wC_abl_Q       )
    CALL allocate_shared_dp_0D( SMB%C_refr,         SMB%wC_refr        )
    
    IF (par%master) THEN
      IF     (region_name == 'NAM') THEN
        SMB%C_abl_constant = C%C_abl_constant_NAM
        SMB%C_abl_Ts       = C%C_abl_Ts_NAM
        SMB%C_abl_Q        = C%C_abl_Q_NAM
        SMB%C_refr         = C%C_refr_NAM
      ELSEIF (region_name == 'EAS') THEN
        SMB%C_abl_constant = C%C_abl_constant_EAS
        SMB%C_abl_Ts       = C%C_abl_Ts_EAS
        SMB%C_abl_Q        = C%C_abl_Q_EAS
        SMB%C_refr         = C%C_refr_EAS
      ELSEIF (region_name == 'GRL') THEN
        SMB%C_abl_constant = C%C_abl_constant_GRL
        SMB%C_abl_Ts       = C%C_abl_Ts_GRL
        SMB%C_abl_Q        = C%C_abl_Q_GRL
        SMB%C_refr         = C%C_refr_GRL
      ELSEIF (region_name == 'ANT') THEN
        SMB%C_abl_constant = C%C_abl_constant_ANT
        SMB%C_abl_Ts       = C%C_abl_Ts_ANT
        SMB%C_abl_Q        = C%C_abl_Q_ANT
        SMB%C_refr         = C%C_refr_ANT
      END IF
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Initialise albedo to background albedo
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! Background albedo
      IF (init%Hb( j,i) < 0._dp) THEN
        SMB%AlbedoSurf( j,i) = albedo_water
      ELSE
        SMB%AlbedoSurf( j,i) = albedo_soil
      END IF
      
      IF (init%Hi( j,i) > 0._dp) THEN
        SMB%AlbedoSurf(  j,i) = albedo_snow
        SMB%FirnDepth( :,j,i) = initial_snow_depth   
      END IF
      
      SMB%Albedo( :,j,i) = SMB%AlbedoSurf( j,i)
      
    END DO
    END DO
    CALL sync
    
    IF (C%is_restart) THEN
      SMB%FirnDepth(        :,:,grid%i1:grid%i2) = init%FirnDepth(        :,:,grid%i1:grid%i2)
      SMB%MeltPreviousYear(   :,grid%i1:grid%i2) = init%MeltPreviousYear(   :,grid%i1:grid%i2)
      CALL sync
    END IF
  
  END SUBROUTINE initialise_SMB_model

END MODULE SMB_module
