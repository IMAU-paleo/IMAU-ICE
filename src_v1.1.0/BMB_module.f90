MODULE BMB_module

  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parallel_module,                 ONLY: par, sync, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_subclimate_region, type_BMB_model
  USE parameters_module,               ONLY: T0, L_fusion, seawater_density, ice_density, sec_per_year
  USE forcing_module,                  ONLY: forcing

  IMPLICIT NONE
    
CONTAINS

  ! Run the SMB model on the region grid
  SUBROUTINE run_BMB_model( grid, ice, climate, BMB)
    ! Calculate mean ocean temperature (saved in "climate") and basal mass balance
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: i,j
    
    REAL(dp)                                           :: BMB_shelf                             ! Sub-shelf melt rate for non-exposed shelf  [m/year]
    REAL(dp)                                           :: BMB_shelf_exposed                     ! Sub-shelf melt rate for exposed shelf      [m/year]
    REAL(dp)                                           :: BMB_deepocean                         ! Sub-shelf melt rate for deep-ocean areas   [m/year]
    REAL(dp)                                           :: w_PD, w_warm, w_cold, w_deep, w_expo
    
    REAL(dp)                                           :: T_freeze                              ! Freezing temperature at the base of the shelf (Celcius)
    
    REAL(dp), PARAMETER                                :: cp0        = 3974._dp                 ! specific heat capacity of the ocean mixed layer (J kg-1 K-1) 
    REAL(dp), PARAMETER                                :: F_melt     = 5.0E-03_dp               ! Melt parameter (m s-1) (PISM-PIK original is 5E-03)
    REAL(dp), PARAMETER                                :: gamma_T    = 1.0E-04_dp               ! Thermal exchange velocity (m s-1)
    
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
          
        BMB%BMB( :,grid%i1:grid%i2) = 0._dp
        CALL sync
        RETURN
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_BMB_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
      
    ! Initialise everything at zero
    BMB%BMB(       :,grid%i1:grid%i2) = 0._dp
    BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
    BMB%BMB_shelf( :,grid%i1:grid%i2) = 0._dp
    BMB%sub_angle( :,grid%i1:grid%i2) = 360._dp
    BMB%dist_open( :,grid%i1:grid%i2) = 0._dp
    w_PD                              = 0._dp
    w_warm                            = 0._dp
    w_cold                            = 0._dp
    w_deep                            = 0._dp
    w_expo                            = 0._dp
    BMB_deepocean                     = 0._dp
    BMB_shelf_exposed                 = 0._dp
    
    ! Find the "subtended angle" and distance-to-open-ocean of all shelf pixels
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%mask_shelf_Aa( j,i) == 1) THEN
        BMB%sub_angle( j,i) = subtended_angle(     grid, i, j, ice%mask_land_Aa, ice%mask_sheet_Aa, ice%mask_shelf_Aa, ice%mask_ocean_Aa)
        BMB%dist_open( j,i) = distance_open_ocean( grid, i, j,                                      ice%mask_shelf_Aa, ice%mask_ocean_Aa)
      END IF
    END DO
    END DO
    
    ! Determine mean ocean temperature and basal melt rates for deep ocean and exposed shelves
    IF (C%choice_ocean_temperature_model == 'fixed') THEN
      ! Use present-day values
      
      climate%T_ocean_mean = C%ocean_temperature_PD
      BMB_deepocean        = C%BMB_deepocean_PD
      BMB_shelf_exposed    = C%BMB_shelf_exposed_PD
      
    ELSEIF (C%choice_ocean_temperature_model == 'scaled') THEN
      ! Scale between config values of mean ocean temperature and basal melt rates for PD, cold, and warm climates.
    
      ! Determine weight for scaling between different ocean temperatures
      IF (C%choice_forcing_method == 'CO2_direct') THEN
      
        ! Use the prescribed CO2 record as a glacial index
        IF (forcing%CO2_obs > 280._dp) THEN
          ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 405 ppmv
          w_PD   = MIN(1._dp, MAX(0._dp, (405._dp - forcing%CO2_obs) / (405._dp - 280._dp) ))
          w_warm = 1._dp - w_PD
          w_cold = 0._dp
        ELSE
          ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
          w_PD   = MIN(1._dp, MAX(0._dp, (forcing%CO2_obs - 190._dp) / (280._dp - 190._dp) ))
          w_cold = 1._dp - w_PD
          w_warm = 0._dp
        END IF
        
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      
        ! Use modelled CO2 as a glacial index
        IF (forcing%CO2_mod > 280._dp) THEN
          ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 405 ppmv
          w_PD   = MIN(1._dp, MAX(0._dp, (405._dp - forcing%CO2_mod) / (405._dp - 280._dp) ))
          w_warm = 1._dp - w_PD
          w_cold = 0._dp
        ELSE
          ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
          w_PD   = MIN(1._dp, MAX(0._dp, (forcing%CO2_mod - 190._dp) / (280._dp - 190._dp) ))
          w_cold = 1._dp - w_PD
          w_warm = 0._dp
        END IF
        
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      
        ! Use modelled global mean annual temperature change as a glacial index
        IF (forcing%dT_glob_inverse > 0._dp) THEN
          ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 405 ppmv
          w_warm = MIN(1._dp, MAX(0._dp, forcing%dT_glob_inverse / 12._dp ))
          w_PD   = 1._dp - w_warm
          w_cold = 0._dp
        ELSE
          ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
          w_cold = MIN(1._dp, MAX(0._dp, -forcing%dT_glob_inverse / 12._dp ))
          w_PD   = 1._dp - w_cold
          w_warm = 0._dp
        END IF
        
      ELSE ! IF (C%choice_forcing_method == 'CO2_direct') THEN
        WRITE(0,*) '  ERROR: forcing method "', TRIM(C%choice_forcing_method), '" not implemented in run_BMB_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF ! IF (C%choice_forcing_method == 'CO2_direct') THEN
      
      climate%T_ocean_mean = w_PD * C%ocean_temperature_PD + w_warm * C%ocean_temperature_warm + w_cold * C%ocean_temperature_cold
      BMB_deepocean        = w_PD * C%BMB_deepocean_PD     + w_warm * C%BMB_deepocean_warm     + w_cold * C%BMB_deepocean_cold
      BMB_shelf_exposed    = w_PD * C%BMB_shelf_exposed_PD + w_warm * C%BMB_shelf_exposed_warm + w_cold * C%BMB_shelf_exposed_cold
      
    ELSE ! IF (C%choice_ocean_temperature_model == 'fixed') THEN
      WRITE(0,*) '  ERROR: choice_ocean_temperature_model "', TRIM(C%choice_ocean_temperature_model), '" not implemented in run_BMB_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF ! IF (C%choice_ocean_temperature_model == 'fixed') THEN
    
    ! Use the (interpolated, spatially uniform) ocean temperature and the subtended angle + distance-to-open-ocean
    ! to calculate sub-shelf melt rates using the parametrisation from Martin et al., 2011
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      IF (ice%mask_shelf_Aa( j,i) == 1) THEN
        ! Sub-shelf melt

        ! Freezing temperature at the bottom of the ice shelves, scaling with depth below water level
        T_freeze = T0 + 0.0939_dp - 0.057_dp * 35._dp - 7.64E-04_dp * ice%Hi_Aa( j,i) * ice_density / seawater_density

        ! Sub-shelf melt rate for non-exposed shelves (Martin, TC, 2011) - melt values, when T_ocean > T_freeze.
        BMB_shelf   = seawater_density * cp0 * sec_per_year * gamma_T * C%subshelf_melt_factor * &
                   (climate%T_ocean_mean - T_freeze) / (L_fusion * ice_density)

      ELSE
        BMB_shelf = 0._dp
      END IF

      IF (ice%mask_shelf_Aa( j,i) == 1 .OR. ice%mask_ocean_Aa( j,i) == 1) THEN
      
        w_deep = MAX(0._dp, MIN(1._dp, (ice%SL_Aa( j,i) - ice%Hb_Aa( j,i) - 1800._dp)/200._dp))
        w_expo = MAX(0._dp, MIN(1._dp, (BMB%sub_angle( j,i) - 80._dp)/30._dp)) * EXP(-BMB%dist_open( j,i)* grid%dx/100000._dp)
        
        BMB%BMB_shelf( j,i) = (w_deep * BMB_deepocean) + (1._dp - w_deep) * (w_expo * BMB_shelf_exposed + (1._dp - w_expo) * BMB_shelf)
        
      ELSE  
        BMB%BMB_shelf( j,i) = 0._dp
      END IF

    END DO
    END DO
    CALL sync
    
    ! Add sheet and shelf melt together
    BMB%BMB( :,grid%i1:grid%i2) = BMB%BMB_sheet( :,grid%i1:grid%i2) + BMB%BMB_shelf( :,grid%i1:grid%i2)
    CALL sync
          
  END SUBROUTINE run_BMB_model
  
  FUNCTION distance_open_ocean( grid, i_shelf, j_shelf, mask_shelf, mask_ocean) RESULT( open_distance)
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
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_shelf
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_ocean
    REAL(dp)                                           :: open_distance
    
    ! Local variables:
    REAL(dp), DIMENSION(16)                            :: distance             ! list of all distances for 16 directions (m)
    INTEGER                                            :: iloop,imax           ! loop integer and maximum value (max of NX/NY)
    INTEGER                                            :: ii,jj
    
    imax          = MAX( grid%nx, grid%ny)
    distance      = 0
    open_distance = 0._dp
    
    IF (mask_shelf( j_shelf,i_shelf) == 0) RETURN
   
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
        
        IF (mask_ocean( jj,ii) == 1) THEN
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
        
        IF (mask_ocean( jj,ii) == 1) THEN
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

        IF (mask_ocean( jj,ii) == 1) THEN
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
        
        IF (mask_ocean( jj,ii) == 1) THEN
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

        IF (mask_ocean( jj,ii) == 1) THEN
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
        
        IF (mask_ocean( jj,ii) == 1) THEN
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
                 
        IF (mask_ocean( jj,ii) == 1) THEN
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
                  
        IF (mask_ocean( jj,ii) == 1) THEN
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
                  
        IF (mask_ocean( jj,ii) == 1) THEN
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
                  
        IF (mask_ocean( jj,ii) == 1) THEN
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

        IF (mask_ocean( jj,ii) == 1) THEN
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

        IF (mask_ocean( jj,ii) == 1) THEN
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

        IF (mask_ocean( jj,ii) == 1) THEN
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

        IF (mask_ocean( jj,ii) == 1) THEN
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

        IF (mask_ocean( jj,ii) == 1) THEN
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

        IF (mask_ocean( jj,ii) == 1) THEN
          distance(16) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! Calculate minimum distance
      open_distance = MINVAL(distance(1:16))

    END FUNCTION distance_open_ocean
  FUNCTION subtended_angle( grid, i_shelf, j_shelf, mask_land, mask_sheet, mask_shelf, mask_ocean) RESULT(angle_sub)
    ! Determine the subtended angle to the open ocean. When encountering land points
    ! the angle is zet to zero (i.e. a loss of 1/16th of the fully possible 360). 
    ! The minimum value is the solution
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    INTEGER,                             INTENT(IN)    :: j_shelf, i_shelf
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_land
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_sheet
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_shelf
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_ocean
    REAL(dp)                                           :: angle_sub
    
    ! Local variables:
    REAL(dp), DIMENSION(16)                            :: angle                !  list of open-ocean angles for 16 directions (degrees)
    INTEGER                                            :: iloop,imax           ! loop integer and maximum value (max of NX/NY)
    INTEGER                                            :: ii,jj

    imax      = MAX(grid%nx,grid%ny)
    angle     = 1._dp
    angle_sub = 360._dp
    
    IF (mask_shelf( j_shelf,i_shelf) == 0) RETURN
    
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
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
      ELSE IF (mask_ocean( jj,ii) == 1) THEN
        angle(16) = 1._dp
        EXIT
      END IF    
    END DO

    ! Calculate total subtanded angle
    angle_sub = SUM(angle(1:16)) * 360._dp / 16._dp

  END FUNCTION subtended_angle
  
  ! Administration: allocation, initialisation, and remapping
  SUBROUTINE initialise_BMB_model( grid, BMB)
    ! Allocate memory for the data fields of the SMB model.
    
    IMPLICIT NONE
    
    !In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    IF (par%master) WRITE (0,*) '  Initialising BMB model...'
    
    ! Allocate memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%BMB      , BMB%wBMB      )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%BMB_sheet, BMB%wBMB_sheet)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%BMB_shelf, BMB%wBMB_shelf)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%sub_angle, BMB%wsub_angle)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%dist_open, BMB%wdist_open)
      
  END SUBROUTINE initialise_BMB_model
  
END MODULE BMB_module
