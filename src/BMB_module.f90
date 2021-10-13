MODULE BMB_module

  ! Contains all the routines for calculating the basal mass balance.

  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_subclimate_region, type_BMB_model
  USE netcdf_module,                   ONLY: debug, write_to_debug_file 
  USE parameters_module,               ONLY: T0, L_fusion, seawater_density, ice_density, sec_per_year
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  USE forcing_module,                  ONLY: forcing

  IMPLICIT NONE
    
CONTAINS

  ! The main routine that is called from the IMAU_ICE_main_model
  SUBROUTINE run_BMB_model( grid, ice, climate, BMB, region_name)
    ! Run the selected BMB model
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables:
    INTEGER                                            :: i,j
    
  ! ================================================
  ! ===== Exceptions for benchmark experiments =====
  ! ================================================
  
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
        BMB%BMB( :,grid%i1:grid%i2) = 0._dp
        CALL sync
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'MISMIPplus') THEN
        ! No exception; use the actual basal melt model
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_BMB_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  ! =======================================================
  ! ===== End of exceptions for benchmark experiments =====
  ! =======================================================
    
    ! Run the selected shelf BMB model
    IF     (C%choice_BMB_shelf_model == 'uniform') THEN
      BMB%BMB_shelf( :,grid%i1:grid%i2) = C%BMB_shelf_uniform
      CALL sync
    ELSEIF (C%choice_BMB_shelf_model == 'ANICE_legacy') THEN
      CALL run_BMB_model_ANICE_legacy( grid, ice, climate, BMB, region_name)
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_shelf_model "', TRIM(C%choice_BMB_shelf_model), '" not implemented in run_BMB_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Run the selected sheet BMB model
    IF     (C%choice_BMB_sheet_model == 'uniform') THEN
      BMB%BMB_sheet( :,grid%i1:grid%i2) = C%BMB_sheet_uniform
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_sheet_model "', TRIM(C%choice_BMB_sheet_model), '" not implemented in run_BMB_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Add sheet and shelf melt rates together, applying the selected scheme for sub-grid shelf melt
    ! (see Leguy et al. 2021 for explanations of the three schemes)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! No sub-grid scaling for sub-sheet melt yet
      BMB%BMB( j,i) = 0._dp
      IF (ice%mask_sheet_a( j,i) == 1._dp) BMB%BMB( j,i) = BMB%BMB_sheet( j,i)
      
      ! Different sub-grid schemes for sub-shelf melt
      IF     (C%choice_BMB_subgrid == 'FCMP') THEN
        IF (ice%mask_shelf_a( j,i) == 1) BMB%BMB( j,i) = BMB%BMB( j,i) + BMB%BMB_shelf( j,i)
      ELSEIF (C%choice_BMB_subgrid == 'PMP') THEN
        BMB%BMB( j,i) = BMB%BMB( j,i) + (1._dp - ice%f_grnd_a( j,i)) * BMB%BMB_shelf( j,i)
      ELSEIF (C%choice_BMB_subgrid == 'NMP') THEN
        IF (ice%f_grnd_a( j,i) == 0._dp) BMB%BMB( j,i) = BMB%BMB( j,i) + BMB%BMB_shelf( j,i)
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_subgrid "', TRIM(C%choice_BMB_subgrid), '" not implemented in run_BMB_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE run_BMB_model

  ! The ANICE_legacy sub-shelf melt model
  SUBROUTINE run_BMB_model_ANICE_legacy( grid, ice, climate, BMB, region_name)
    ! Calculate sub-shelf melt with the ANICE_legacy model, which is based on the glacial-interglacial
    ! parameterisation by Pollard & DeConto (2012), the distance-to-open-ocean + subtended-angle parameterisation
    ! by Pollard & DeConto (2012), and the linear temperature-melt relation by Martin et al. (2011).
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables
    INTEGER                                            :: i,j
    REAL(dp)                                           :: BMB_shelf                             ! Sub-shelf melt rate for non-exposed shelf  [m/year]
    REAL(dp)                                           :: BMB_shelf_exposed                     ! Sub-shelf melt rate for exposed shelf      [m/year]
    REAL(dp)                                           :: BMB_deepocean                         ! Sub-shelf melt rate for deep-ocean areas   [m/year]
    REAL(dp)                                           :: w_ins, w_PD, w_warm, w_cold, w_deep, w_expo
    REAL(dp)                                           :: T_freeze                              ! Freezing temperature at the base of the shelf (Celcius)
    REAL(dp)                                           :: water_depth
    REAL(dp), PARAMETER                                :: cp0        = 3974._dp                 ! specific heat capacity of the ocean mixed layer (J kg-1 K-1) 
    REAL(dp), PARAMETER                                :: gamma_T    = 1.0E-04_dp               ! Thermal exchange velocity (m s-1)
      
    ! Initialise everything at zero
    BMB%BMB(       :,grid%i1:grid%i2) = 0._dp
    BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
    BMB%BMB_shelf( :,grid%i1:grid%i2) = 0._dp
    BMB%sub_angle( :,grid%i1:grid%i2) = 360._dp
    BMB%dist_open( :,grid%i1:grid%i2) = 0._dp
    w_ins                             = 0._dp
    w_PD                              = 0._dp
    w_warm                            = 0._dp
    w_cold                            = 0._dp
    w_deep                            = 0._dp
    w_expo                            = 0._dp
    BMB_shelf                         = 0._dp
    BMB_shelf_exposed                 = 0._dp
    BMB_deepocean                     = 0._dp
    CALL sync
    
    ! Find the "subtended angle" and distance-to-open-ocean of all shelf pixels
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%mask_shelf_a( j,i) == 1) THEN
        BMB%sub_angle( j,i) = subtended_angle(     grid, i, j, ice%mask_land_a, ice%mask_ocean_a, ice%mask_ice_a, ice%mask_sheet_a, ice%mask_shelf_a)
        BMB%dist_open( j,i) = distance_open_ocean( grid, i, j, ice%mask_land_a, ice%mask_ocean_a, ice%mask_ice_a,                   ice%mask_shelf_a)
      END IF
    END DO
    END DO
    CALL sync
    
    ! Find the weight from insolation
    IF (region_name == 'NAM' .OR. region_name == 'EAS' .OR. region_name == 'GRL') THEN
      w_ins = MAX(0._dp, (climate%Q_TOA_jun_65N - 462.29_dp) / 40._dp)
    ELSEIF (region_name == 'ANT') THEN
      w_ins = MAX(0._dp, (climate%Q_TOA_jan_80S - 532.19_dp) / 40._dp)
    END IF
    
    ! Determine mean ocean temperature and basal melt rates for deep ocean and exposed shelves
    IF (C%choice_ocean_temperature_model == 'fixed') THEN
      ! Use present-day values
      
      climate%T_ocean_mean = BMB%T_ocean_mean_PD
      BMB_deepocean        = BMB%BMB_deepocean_PD
      BMB_shelf_exposed    = BMB%BMB_shelf_exposed_PD
      
    ELSEIF (C%choice_ocean_temperature_model == 'scaled') THEN
      ! Scale between config values of mean ocean temperature and basal melt rates for PD, cold, and warm climates.
    
      ! Determine weight for scaling between different ocean temperatures
      IF (C%choice_forcing_method == 'CO2_direct') THEN
      
        ! Use the prescribed CO2 record as a glacial index
        IF (forcing%CO2_obs > 280._dp) THEN
          ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 400 ppmv
          w_PD   = MIN(1.25_dp, MAX(-0.25_dp, (400._dp - forcing%CO2_obs) / (400._dp - 280._dp) )) - w_ins
          w_warm = 1._dp - w_PD
          w_cold = 0._dp
        ELSE
          ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
          w_PD   = MIN(1.25_dp, MAX(-0.25_dp, (forcing%CO2_obs - 190._dp) / (280._dp - 190._dp) )) + w_ins
          w_cold = 1._dp - w_PD
          w_warm = 0._dp
        END IF
        
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      
        ! Use modelled CO2 as a glacial index
        IF (forcing%CO2_mod > 280._dp) THEN
          ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 400 ppmv
          w_PD   = MIN(1.25_dp, MAX(-0.25_dp, (400._dp - forcing%CO2_mod) / (400._dp - 280._dp) )) - w_ins
          w_warm = 1._dp - w_PD
          w_cold = 0._dp
        ELSE
          ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
          w_PD   = MIN(1.25_dp, MAX(-0.25_dp, (forcing%CO2_mod - 190._dp) / (280._dp - 190._dp) )) + w_ins
          w_cold = 1._dp - w_PD
          w_warm = 0._dp
        END IF
        
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      
        ! Use modelled global mean annual temperature change as a glacial index
        IF (forcing%dT_glob_inverse > 0._dp) THEN
          ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 405 ppmv
          w_warm = MIN(1.25_dp, MAX(-0.25_dp, forcing%dT_glob_inverse / 12._dp )) - w_ins
          w_PD   = 1._dp - w_warm
          w_cold = 0._dp
        ELSE
          ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
          w_cold = MIN(1.25_dp, MAX(-0.25_dp, -forcing%dT_glob_inverse / 12._dp )) + w_ins
          w_PD   = 1._dp - w_cold
          w_warm = 0._dp
        END IF
        
      ELSEIF (C%choice_forcing_method == 'climate_direct' .OR. C%choice_forcing_method == 'SMB_direct') THEN
        ! In this case, no CO2/d18O forcing is used; just assume PD weights
        
        w_warm = 0._dp
        w_cold = 0._dp
        w_PD   = 1._dp
        
      ELSE ! IF (C%choice_forcing_method == 'CO2_direct') THEN
        WRITE(0,*) '  ERROR: forcing method "', TRIM(C%choice_forcing_method), '" not implemented in run_BMB_model_ANICE_legacy!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF ! IF (C%choice_forcing_method == 'CO2_direct') THEN
      
      climate%T_ocean_mean = w_PD * BMB%T_ocean_mean_PD      + w_warm * BMB%T_ocean_mean_warm      + w_cold * BMB%T_ocean_mean_cold
      BMB_deepocean        = w_PD * BMB%BMB_deepocean_PD     + w_warm * BMB%BMB_deepocean_warm     + w_cold * BMB%BMB_deepocean_cold
      BMB_shelf_exposed    = w_PD * BMB%BMB_shelf_exposed_PD + w_warm * BMB%BMB_shelf_exposed_warm + w_cold * BMB%BMB_shelf_exposed_cold
      CALL sync
      
    ELSE ! IF (C%choice_ocean_temperature_model == 'fixed') THEN
      WRITE(0,*) '  ERROR: choice_ocean_temperature_model "', TRIM(C%choice_ocean_temperature_model), '" not implemented in run_BMB_model_ANICE_legacy!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF ! IF (C%choice_ocean_temperature_model == 'fixed') THEN
    
    ! Use the (interpolated, spatially uniform) ocean temperature and the subtended angle + distance-to-open-ocean
    ! to calculate sub-shelf melt rates using the parametrisation from Martin et al., 2011
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      IF (ice%mask_shelf_a( j,i) == 1) THEN
        ! Sub-shelf melt

        ! Freezing temperature at the bottom of the ice shelves, scaling with depth below water level
        T_freeze = 0.0939_dp - 0.057_dp * 35._dp - 7.64E-04_dp * ice%Hi_a( j,i) * ice_density / seawater_density

        ! Sub-shelf melt rate for non-exposed shelves (Martin, TC, 2011) - melt values, when T_ocean > T_freeze.
        BMB_shelf   = seawater_density * cp0 * sec_per_year * gamma_T * BMB%subshelf_melt_factor * &
                   (climate%T_ocean_mean - T_freeze) / (L_fusion * ice_density)

      ELSE
        BMB_shelf = 0._dp
      END IF

      IF (ice%mask_shelf_a( j,i) == 1 .OR. ice%mask_ocean_a( j,i) == 1) THEN
      
        water_depth = ice%SL_a( j,i) - ice%Hb_a( j,i)
        w_deep = MAX(0._dp, MIN(1._dp, (water_depth - BMB%deep_ocean_threshold_depth) / 200._dp))
        w_expo = MAX(0._dp, MIN(1._dp, (BMB%sub_angle( j,i) - 80._dp)/30._dp)) * EXP(-BMB%dist_open( j,i) / 100000._dp)
        
        BMB%BMB_shelf( j,i) = (w_deep * BMB_deepocean) + (1._dp - w_deep) * (w_expo * BMB_shelf_exposed + (1._dp - w_expo) * BMB_shelf)
        
      ELSE  
        BMB%BMB_shelf( j,i) = 0._dp
      END IF

    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( BMB%sub_angle, 'BMB%sub_angle', 'run_BMB_model_ANICE_legacy')
    CALL check_for_NaN_dp_2D( BMB%dist_open, 'BMB%dist_open', 'run_BMB_model_ANICE_legacy')
    CALL check_for_NaN_dp_2D( BMB%BMB_shelf, 'BMB%BMB_shelf', 'run_BMB_model_ANICE_legacy')
          
  END SUBROUTINE run_BMB_model_ANICE_legacy
  
  ! The distance-to-open-ocean and subtended-angle functions for the Pollard&DeConto (2012) sub-shelf melt parameterisation
  FUNCTION distance_open_ocean( grid, i_shelf, j_shelf, mask_land, mask_ocean, mask_ice, mask_shelf) RESULT( open_distance)
    ! Determine the distance to the open ocean (expressed in number of grid cells). A  solution is found for all 16 directions
    ! also when encountering land points, keep searching to find a grid point. 
    ! The minimum value is the solution

    ! The 16 directions
    !
    !              12 13 14
    !             11 ---- 15
    !            10 --  -- 16
    !           09 -- ij -- 01
    !            08 --  -- 02 
    !             07 ---- 03
    !              06 05 04 
    !

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    INTEGER,                             INTENT(IN)    :: j_shelf, i_shelf
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_land
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_ocean
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_ice
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_shelf
    REAL(dp)                                           :: open_distance
    
    ! Local variables:
    REAL(dp), DIMENSION(16)                            :: distance             ! list of all distances for 16 directions (m)
    INTEGER                                            :: iloop,imax           ! loop integer and maximum value (max of NX/NY)
    INTEGER                                            :: ii,jj
    
    imax          = MAX( grid%nx, grid%ny)
    distance      = 0
    
    IF (mask_shelf( j_shelf,i_shelf) == 0) THEN
      IF     (mask_land(  j_shelf,i_shelf) == 1) THEN
        open_distance = REAL(grid%nx+grid%ny,dp) * grid%dx
      ELSE
        open_distance = 0._dp
      END IF
      RETURN
    END IF
   
    ! direction 01: (should work fine with imax, no possibility that jj or ii will exceed nx or ny)
    ! but what about when jj or ii are lower than 1..
    DO iloop = 1,imax
      ii = i_shelf
      jj = j_shelf + iloop
       
        ! first check if the boundaries of the grid are reached
        IF (jj > grid%ny) THEN
          distance(1) = 100._dp
          EXIT
        END IF  
        
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(1) = real(jj) - real(j_shelf)
          EXIT
        END IF    
      END DO

      ! direction 02:
      DO iloop = 1,imax   
        ii = i_shelf - iloop
        jj = j_shelf + 2*iloop
        IF (jj > grid%ny .OR. ii < 1) THEN
          distance(2) = 100._dp
          EXIT
        END IF  
        
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(2) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 03:
      DO iloop = 1,imax   
        ii = i_shelf - iloop
        jj = j_shelf + iloop
        IF (jj > grid%ny .OR. ii < 1) THEN
          distance(3) = 100._dp
          EXIT
        END IF  

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(3) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 04:
      DO iloop = 1,imax   
        ii = i_shelf - 2*iloop
        jj = j_shelf + iloop
        IF (jj > grid%ny .OR. ii < 1) THEN
          distance(4) = 100._dp
          EXIT
        END IF  
        
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(4) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 05:
      DO iloop = 1,imax   
        ii = i_shelf - iloop
        jj = j_shelf
        IF ( ii < 1) THEN
          distance(5) = 100._dp
          EXIT
        END IF  

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(5) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 06:
      DO iloop = 1,imax   
        ii = i_shelf - 2*iloop
        jj = j_shelf - iloop
        IF (jj < 1 .OR. ii < 1) THEN
          distance(6) = 100._dp
          EXIT
        END IF
        
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(6) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 07:
      DO iloop = 1,imax   
        ii = i_shelf - iloop
        jj = j_shelf - iloop
        IF (jj < 1 .OR. ii < 1) THEN
          distance(7) = 100._dp
          EXIT
        END IF 
                 
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(7) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 08:
      DO iloop = 1,imax   
        ii = i_shelf - iloop
        jj = j_shelf - 2*iloop
        IF (jj < 1 .OR. ii < 1) THEN
          distance(8) = 100._dp
          EXIT
        END IF
                  
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(8) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 09:
      DO iloop = 1,imax   
        ii = i_shelf 
        jj = j_shelf - iloop
        IF (jj < 1) THEN
          distance(9) = 100._dp
          EXIT
        END IF
                  
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(9) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 10:
      DO iloop = 1,imax   
        ii = i_shelf + iloop 
        jj = j_shelf - 2*iloop
        IF (jj < 1 .OR. ii > grid%nx) THEN
          distance(10) = 100._dp
          EXIT
        END IF
                  
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(10) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 11:
      DO iloop = 1,imax   
        ii = i_shelf + iloop 
        jj = j_shelf - iloop
        IF (jj < 1 .OR. ii > grid%nx) THEN
          distance(11) = 100._dp
          EXIT
        END IF

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(11) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 12:
      DO iloop = 1,imax   
        ii = i_shelf + 2*iloop 
        jj = j_shelf - iloop
        IF (jj < 1 .OR. ii > grid%nx) THEN
          distance(12) = 100._dp
          EXIT
        END IF

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(12) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 13:
      DO iloop = 1,imax   
        ii = i_shelf + iloop 
        jj = j_shelf
        IF (ii > grid%nx) THEN
          distance(13) = 100._dp
          EXIT
        END IF

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(13) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 14:
      DO iloop = 1,imax   
        ii = i_shelf + 2*iloop 
        jj = j_shelf + iloop
        IF (jj > grid%ny .OR. ii > grid%nx) THEN
          distance(14) = 100._dp
          EXIT
        END IF

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(14) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 15:
      DO iloop = 1,imax   
        ii = i_shelf + iloop 
        jj = j_shelf + iloop
        IF (jj > grid%ny .OR. ii > grid%nx) THEN
          distance(15) = 100._dp
          EXIT
        END IF

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(15) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 16:
      DO iloop = 1,imax   
        ii = i_shelf + iloop 
        jj = j_shelf + 2*iloop
        IF (jj > grid%ny .OR. ii > grid%nx) THEN
          distance(16) = 100._dp
          EXIT
        END IF

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(16) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! Calculate minimum distance
      open_distance = MINVAL(distance) * grid%dx

    END FUNCTION distance_open_ocean
  FUNCTION subtended_angle( grid, i_shelf, j_shelf, mask_land, mask_ocean, mask_ice, mask_sheet, mask_shelf) RESULT(angle_sub)
    ! Determine the subtended angle to the open ocean. When encountering land points
    ! the angle is zet to zero (i.e. a loss of 1/16th of the fully possible 360). 
    ! The minimum value is the solution
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    INTEGER,                             INTENT(IN)    :: j_shelf, i_shelf
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_land
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_ocean
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_ice
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_sheet
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_shelf
    REAL(dp)                                           :: angle_sub
    
    ! Local variables:
    REAL(dp), DIMENSION(16)                            :: angle                !  list of open-ocean angles for 16 directions (degrees)
    INTEGER                                            :: iloop,imax           ! loop integer and maximum value (max of NX/NY)
    INTEGER                                            :: ii,jj

    imax      = MAX(grid%nx,grid%ny)
    angle     = 1._dp
    angle_sub = 360._dp
    
    IF (mask_shelf( j_shelf,i_shelf) == 0) THEN
      IF     (mask_land(  j_shelf,i_shelf) == 1) THEN
        angle_sub = 0._dp
      ELSE
        angle_sub = 360._dp
      END IF
      RETURN
    END IF
    
    ! direction 01:
    DO iloop = 1,imax   
      ii = i_shelf
      jj = j_shelf + iloop
      
      ! first check if the boundaries of the grid are reached
      IF (jj > grid%ny) THEN
        angle(1) = 0._dp
        EXIT
      END IF  
      
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(1) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(1) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 02:
    DO iloop = 1,imax   
      ii = i_shelf - iloop
      jj = j_shelf + 2*iloop
      IF (jj > grid%ny .OR. ii < 1) THEN
        angle(2) = 0._dp
        EXIT
      END IF  
      
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(2) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(2) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 03:
    DO iloop = 1,imax   
      ii = i_shelf - iloop
      jj = j_shelf + iloop
      IF (jj > grid%ny .OR. ii < 1) THEN
        angle(3) = 0._dp
        EXIT
      END IF  

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(3) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(3) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 04:
    DO iloop = 1,imax   
      ii = i_shelf - 2*iloop
      jj = j_shelf + iloop
      IF (jj > grid%ny .OR. ii < 1) THEN
        angle(4) = 0._dp
        EXIT
      END IF  
      
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(4) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(4) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 05:
    DO iloop = 1,imax   
      ii = i_shelf - iloop
      jj = j_shelf
      IF ( ii < 1) THEN
        angle(5) = 0._dp
        EXIT
      END IF  

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(5) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(5) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 06:
    DO iloop = 1,imax   
      ii = i_shelf - 2*iloop
      jj = j_shelf - iloop
      IF (jj < 1 .OR. ii < 1) THEN
        angle(6) = 0._dp
        EXIT
      END IF
      
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(6) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(6) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 07:
    DO iloop = 1,imax   
      ii = i_shelf - iloop
      jj = j_shelf - iloop
      IF (jj < 1 .OR. ii < 1) THEN
        angle(7) = 0._dp
        EXIT
      END IF 
               
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(7) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(7) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 08:
    DO iloop = 1,imax   
      ii = i_shelf - iloop
      jj = j_shelf - 2*iloop
      IF (jj < 1 .OR. ii < 1) THEN
        angle(8) = 0._dp
        EXIT
      END IF
                
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(8) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(8) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 09:
    DO iloop = 1,imax   
      ii = i_shelf 
      jj = j_shelf - iloop
      IF (jj < 1) THEN
        angle(9) = 0._dp
        EXIT
      END IF
                
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(9) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(9) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 10:
    DO iloop = 1,imax   
      ii = i_shelf + iloop 
      jj = j_shelf - 2*iloop
      IF (jj < 1 .OR. ii > grid%nx) THEN
        angle(10) = 0._dp
        EXIT
      END IF
                
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(10) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(10) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 11:
    DO iloop = 1,imax   
      ii = i_shelf + iloop 
      jj = j_shelf - iloop
      IF (jj < 1 .OR. ii > grid%nx) THEN
        angle(11) = 0._dp
        EXIT
      END IF

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(11) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(11) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 12:
    DO iloop = 1,imax   
      ii = i_shelf + 2*iloop 
      jj = j_shelf - iloop
      IF (jj < 1 .OR. ii > grid%nx) THEN
        angle(12) = 0._dp
        EXIT
      END IF

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(12) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(12) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 13:
    DO iloop = 1,imax   
      ii = i_shelf + iloop 
      jj = j_shelf
      IF (ii > grid%nx) THEN
        angle(13) = 0._dp
        EXIT
      END IF

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(13) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(13) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 14:
    DO iloop = 1,imax   
      ii = i_shelf + 2*iloop 
      jj = j_shelf + iloop
      IF (jj > grid%ny .OR. ii > grid%nx) THEN
        angle(14) = 0._dp
        EXIT
      END IF

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(14) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(14) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 15:
    DO iloop = 1,imax   
      ii = i_shelf + iloop 
      jj = j_shelf + iloop
      IF (jj > grid%ny .OR. ii > grid%nx) THEN
        angle(15) = 0._dp
        EXIT
      END IF

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(15) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(15) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 16:
    DO iloop = 1,imax   
      ii = i_shelf + iloop 
      jj = j_shelf + 2*iloop
      IF (jj > grid%ny .OR. ii > grid%nx) THEN
        angle(16) = 0._dp
        EXIT
      END IF

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(16) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(16) = 1._dp
        EXIT
      END IF    
    END DO

    ! Calculate total subtanded angle
    angle_sub = SUM(angle(1:16)) * 360._dp / 16._dp

  END FUNCTION subtended_angle
  
  ! Administration: allocation and initialisation
  SUBROUTINE initialise_BMB_model( grid, BMB, region_name)
    ! Allocate memory for the data fields of the SMB model.
    
    IMPLICIT NONE
    
    !In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    IF (par%master) WRITE (0,*) '  Initialising BMB model...'
    
    ! Total
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%BMB      , BMB%wBMB      )
    
    ! Shelf
    IF     (C%choice_BMB_shelf_model == 'uniform') THEN
    
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%BMB_shelf, BMB%wBMB_shelf)
      
    ELSEIF (C%choice_BMB_shelf_model == 'ANICE_legacy') THEN
    
      ! Variables
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%BMB_shelf, BMB%wBMB_shelf)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%sub_angle, BMB%wsub_angle)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%dist_open, BMB%wdist_open)
      
      ! Tuning parameters
      CALL allocate_shared_dp_0D( BMB%T_ocean_mean_PD,            BMB%wT_ocean_mean_PD           )
      CALL allocate_shared_dp_0D( BMB%T_ocean_mean_cold,          BMB%wT_ocean_mean_cold         )
      CALL allocate_shared_dp_0D( BMB%T_ocean_mean_warm,          BMB%wT_ocean_mean_warm         )
      CALL allocate_shared_dp_0D( BMB%BMB_deepocean_PD,           BMB%wBMB_deepocean_PD          )
      CALL allocate_shared_dp_0D( BMB%BMB_deepocean_cold,         BMB%wBMB_deepocean_cold        )
      CALL allocate_shared_dp_0D( BMB%BMB_deepocean_warm,         BMB%wBMB_deepocean_warm        )
      CALL allocate_shared_dp_0D( BMB%BMB_shelf_exposed_PD,       BMB%wBMB_shelf_exposed_PD      )
      CALL allocate_shared_dp_0D( BMB%BMB_shelf_exposed_cold,     BMB%wBMB_shelf_exposed_cold    )
      CALL allocate_shared_dp_0D( BMB%BMB_shelf_exposed_warm,     BMB%wBMB_shelf_exposed_warm    )
      CALL allocate_shared_dp_0D( BMB%subshelf_melt_factor,       BMB%wsubshelf_melt_factor      )
      CALL allocate_shared_dp_0D( BMB%deep_ocean_threshold_depth, BMB%wdeep_ocean_threshold_depth)
      
      IF (region_name == 'NAM') THEN
        BMB%T_ocean_mean_PD            = C%T_ocean_mean_PD_NAM
        BMB%T_ocean_mean_cold          = C%T_ocean_mean_cold_NAM
        BMB%T_ocean_mean_warm          = C%T_ocean_mean_warm_NAM
        BMB%BMB_deepocean_PD           = C%BMB_deepocean_PD_NAM
        BMB%BMB_deepocean_cold         = C%BMB_deepocean_cold_NAM
        BMB%BMB_deepocean_warm         = C%BMB_deepocean_warm_NAM
        BMB%BMB_shelf_exposed_PD       = C%BMB_shelf_exposed_PD_NAM
        BMB%BMB_shelf_exposed_cold     = C%BMB_shelf_exposed_cold_NAM
        BMB%BMB_shelf_exposed_warm     = C%BMB_shelf_exposed_warm_NAM
        BMB%subshelf_melt_factor       = C%subshelf_melt_factor_NAM
        BMB%deep_ocean_threshold_depth = C%deep_ocean_threshold_depth_NAM
      ELSEIF (region_name == 'EAS') THEN
        BMB%T_ocean_mean_PD            = C%T_ocean_mean_PD_EAS
        BMB%T_ocean_mean_cold          = C%T_ocean_mean_cold_EAS
        BMB%T_ocean_mean_warm          = C%T_ocean_mean_warm_EAS
        BMB%BMB_deepocean_PD           = C%BMB_deepocean_PD_EAS
        BMB%BMB_deepocean_cold         = C%BMB_deepocean_cold_EAS
        BMB%BMB_deepocean_warm         = C%BMB_deepocean_warm_EAS
        BMB%BMB_shelf_exposed_PD       = C%BMB_shelf_exposed_PD_EAS
        BMB%BMB_shelf_exposed_cold     = C%BMB_shelf_exposed_cold_EAS
        BMB%BMB_shelf_exposed_warm     = C%BMB_shelf_exposed_warm_EAS
        BMB%subshelf_melt_factor       = C%subshelf_melt_factor_EAS
        BMB%deep_ocean_threshold_depth = C%deep_ocean_threshold_depth_EAS
      ELSEIF (region_name == 'GRL') THEN
        BMB%T_ocean_mean_PD            = C%T_ocean_mean_PD_GRL
        BMB%T_ocean_mean_cold          = C%T_ocean_mean_cold_GRL
        BMB%T_ocean_mean_warm          = C%T_ocean_mean_warm_GRL
        BMB%BMB_deepocean_PD           = C%BMB_deepocean_PD_GRL
        BMB%BMB_deepocean_cold         = C%BMB_deepocean_cold_GRL
        BMB%BMB_deepocean_warm         = C%BMB_deepocean_warm_GRL
        BMB%BMB_shelf_exposed_PD       = C%BMB_shelf_exposed_PD_GRL
        BMB%BMB_shelf_exposed_cold     = C%BMB_shelf_exposed_cold_GRL
        BMB%BMB_shelf_exposed_warm     = C%BMB_shelf_exposed_warm_GRL
        BMB%subshelf_melt_factor       = C%subshelf_melt_factor_GRL
        BMB%deep_ocean_threshold_depth = C%deep_ocean_threshold_depth_GRL
      ELSEIF (region_name == 'ANT') THEN
        BMB%T_ocean_mean_PD            = C%T_ocean_mean_PD_ANT
        BMB%T_ocean_mean_cold          = C%T_ocean_mean_cold_ANT
        BMB%T_ocean_mean_warm          = C%T_ocean_mean_warm_ANT
        BMB%BMB_deepocean_PD           = C%BMB_deepocean_PD_ANT
        BMB%BMB_deepocean_cold         = C%BMB_deepocean_cold_ANT
        BMB%BMB_deepocean_warm         = C%BMB_deepocean_warm_ANT
        BMB%BMB_shelf_exposed_PD       = C%BMB_shelf_exposed_PD_ANT
        BMB%BMB_shelf_exposed_cold     = C%BMB_shelf_exposed_cold_ANT
        BMB%BMB_shelf_exposed_warm     = C%BMB_shelf_exposed_warm_ANT
        BMB%subshelf_melt_factor       = C%subshelf_melt_factor_ANT
        BMB%deep_ocean_threshold_depth = C%deep_ocean_threshold_depth_ANT
      END IF
    
    ELSE ! IF     (C%choice_BMB_shelf_model == 'uniform') THEN
      IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_shelf_model "', TRIM(C%choice_BMB_shelf_model), '" not implemented in initialise_BMB_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Sheet
    IF     (C%choice_BMB_sheet_model == 'uniform') THEN
    
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%BMB_sheet, BMB%wBMB_sheet)
    
    ELSE ! IF     (C%choice_BMB_sheet_model == 'uniform') THEN
      IF (par%master) WRITE(0,*) '  ERROR: choice_BMB_sheet_model "', TRIM(C%choice_BMB_sheet_model), '" not implemented in initialise_BMB_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
      
  END SUBROUTINE initialise_BMB_model
  
END MODULE BMB_module
