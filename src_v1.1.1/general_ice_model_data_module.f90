MODULE general_ice_model_data_module

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parallel_module,                 ONLY: par, sync, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE data_types_module,               ONLY: type_grid, type_ice_model
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE parameters_module,               ONLY: seawater_density, ice_density
  USE utilities_module,                ONLY: vertical_average
  USE derivatives_and_grids_module,    ONLY: map_Aa_to_Acx_2D, map_Aa_to_Acy_2D, ddx_Aa_to_Acx_2D, ddy_Aa_to_Acy_2D, &
                                             ddy_Aa_to_Acx_2D, ddx_Aa_to_Acy_2D, map_Aa_to_Acx_3D, map_Aa_to_Acy_3D, &
                                             ddx_Aa_to_Aa_2D, ddy_Aa_to_Aa_2D, map_Aa_to_Ab_2D

  IMPLICIT NONE

CONTAINS
  
  ! Routines for calculating general ice model data - Hs, masks, ice physical properties
  SUBROUTINE update_general_ice_model_data( grid, ice, time)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables
    INTEGER                                            :: i,j
    
    ! Update Hs
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%Hs_Aa( j,i) = ice%Hi_Aa( j,i) + MAX( ice%SL_Aa( j,i) - ice_density / seawater_density * ice%Hi_Aa( j,i), ice%Hb_Aa( j,i))
    END DO
    END DO
    CALL sync
    
    ! Determine masks
    CALL determine_masks( grid, ice)
 
    ! Get ice thickness, flow factor, and surface slopes on the required grids
    CALL map_Aa_to_Acx_2D( grid, ice%Hi_Aa,          ice%Hi_Acx    )
    CALL map_Aa_to_Acy_2D( grid, ice%Hi_Aa,          ice%Hi_Acy    )
    CALL map_Aa_to_Acx_3D( grid, ice%A_flow_Aa,      ice%A_flow_Acx)
    CALL map_Aa_to_Acy_3D( grid, ice%A_flow_Aa,      ice%A_flow_Acy)
    CALL map_Aa_to_Acx_2D( grid, ice%A_flow_mean_Aa, ice%A_flow_mean_Acx)
    CALL map_Aa_to_Acy_2D( grid, ice%A_flow_mean_Aa, ice%A_flow_mean_Acy)
    
    CALL ddx_Aa_to_Aa_2D(  grid, ice%Hi_Aa,          ice%dHi_dx_Aa )
    CALL ddy_Aa_to_Aa_2D(  grid, ice%Hi_Aa,          ice%dHi_dy_Aa )
    
    CALL ddx_Aa_to_Aa_2D(  grid, ice%Hs_Aa,          ice%dHs_dx_Aa )
    CALL ddy_Aa_to_Aa_2D(  grid, ice%Hs_Aa,          ice%dHs_dy_Aa )
    CALL ddx_Aa_to_Acx_2D( grid, ice%Hs_Aa,          ice%dHs_dx_Acx)
    CALL ddy_Aa_to_Acx_2D( grid, ice%Hs_Aa,          ice%dHs_dy_Acx)
    CALL ddx_Aa_to_Acy_2D( grid, ice%Hs_Aa,          ice%dHs_dx_Acy)
    CALL ddy_Aa_to_Acy_2D( grid, ice%Hs_Aa,          ice%dHs_dy_Acy)
            
    ! Calculate physical properties (flow factor, heat capacity, thermal conductivity, pressure melting point)
    CALL ice_physical_properties( grid, ice, time)
    
  END SUBROUTINE update_general_ice_model_data
  SUBROUTINE determine_masks( grid, ice)
    ! Determine the different masks, on both the Aa and the Ac grid
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(INOUT) :: ice 
  
    INTEGER                                            :: i,j
    
    ! Start out with land everywhere, fill in the rest based on input.
    IF (par%master) THEN
      ice%mask_land_Aa    = 1
      ice%mask_ocean_Aa   = 0
      ice%mask_lake_Aa    = 0
      ice%mask_ice_Aa     = 0
      ice%mask_sheet_Aa   = 0
      ice%mask_sheet_Acx  = 0
      ice%mask_sheet_Acy  = 0
      ice%mask_shelf_Aa   = 0
      ice%mask_shelf_Acx  = 0
      ice%mask_shelf_Acy  = 0
      ice%mask_coast_Aa   = 0
      ice%mask_coast_Acx  = 0
      ice%mask_coast_Acy  = 0
      ice%mask_margin_Aa  = 0
      ice%mask_margin_Acx = 0
      ice%mask_margin_Acy = 0
      ice%mask_gl_Aa      = 0
      ice%mask_gl_Acx     = 0
      ice%mask_gl_Acy     = 0
      ice%mask_cf_Aa      = 0
      ice%mask_cf_Acx     = 0
      ice%mask_cf_Acy     = 0
      ice%mask_Aa         = C%type_land
      ice%mask_Acx        = C%type_land
      ice%mask_Acy        = C%type_land
    END IF
    CALL sync
  
    ! First on the Aa grid
    ! ====================
    
    IF (C%do_ocean_floodfill) THEN
      CALL ocean_floodfill( grid, ice)
    ELSE
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (is_floating( ice%Hi_Aa( j,i), ice%Hb_Aa( j,i), ice%SL_Aa( j,i))) THEN
          ice%mask_land_Aa(  j,i) = 0
          ice%mask_ocean_Aa( j,i) = 1
          ice%mask_Aa(       j,i) = C%type_ocean
        END IF
      END DO
      END DO
      CALL sync
    
    END IF ! IF (C%do_ocean_floodfill) THEN
  
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Ice
      IF (ice%Hi_Aa( j,i) > 0._dp) THEN
        ice%mask_ice_Aa(  j,i) = 1
      END IF
      
      ! Sheet
      IF (ice%mask_ice_Aa( j,i) == 1 .AND. ice%mask_land_Aa( j,i) == 1) THEN
        ice%mask_sheet_Aa( j,i) = 1
        ice%mask_Aa(       j,i) = C%type_sheet
      END IF
    
      ! Shelf
      IF (ice%mask_ice_Aa( j,i) == 1 .AND. ice%mask_ocean_Aa( j,i) == 1) THEN
        ice%mask_shelf_Aa( j,i) = 1
        ice%mask_Aa(       j,i) = C%type_shelf
      END IF
      
    END DO
    END DO
    CALL sync
  
    ! Determine coast, grounding line and calving front
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    
      ! Ice-free land bordering ocean equals coastline
      IF (ice%mask_land_Aa( j,i) == 1 .AND. ice%mask_ice_Aa( j,i) == 0) THEN  
        IF (ice%mask_ocean_Aa( j-1,i  ) == 1 .OR. &
            ice%mask_ocean_Aa( j+1,i  ) == 1 .OR. &
            ice%mask_ocean_Aa( j  ,i-1) == 1 .OR. &
            ice%mask_ocean_Aa( j  ,i+1) == 1) THEN
            ice%mask_Aa(       j,i) = C%type_coast
            ice%mask_coast_Aa( j,i) = 1
        END IF
      END IF
    
      ! Ice bordering non-ice equals margin
      IF (ice%mask_ice_Aa( j,i) == 1) THEN  
        IF (ice%mask_ice_Aa( j-1,i  ) == 0 .OR. &
            ice%mask_ice_Aa( j+1,i  ) == 0 .OR. &
            ice%mask_ice_Aa( j  ,i-1) == 0 .OR. &
            ice%mask_ice_Aa( j  ,i+1) == 0) THEN
            ice%mask_Aa(        j,i) = C%type_margin
            ice%mask_margin_Aa( j,i) = 1
        END IF
      END IF
    
      ! Sheet bordering shelf equals grounding line
      IF (ice%mask_sheet_Aa( j,i) == 1) THEN  
        IF (ice%mask_shelf_Aa( j-1,i  ) == 1 .OR. &
            ice%mask_shelf_Aa( j+1,i  ) == 1 .OR. &
            ice%mask_shelf_Aa( j  ,i-1) == 1 .OR. &
            ice%mask_shelf_Aa( j  ,i+1) == 1) THEN
            ice%mask_Aa(    j,i) = C%type_groundingline
            ice%mask_gl_Aa( j,i) = 1
        END IF
      END IF
  
      ! Ice (sheet or shelf) bordering open ocean equals calvingfront
      IF (ice%mask_ice_Aa( j,i) == 1) THEN  
        IF ((ice%mask_ocean_Aa( j-1,i  ) == 1 .AND. ice%mask_ice_Aa( j-1,i  ) == 0) .OR. &
            (ice%mask_ocean_Aa( j+1,i  ) == 1 .AND. ice%mask_ice_Aa( j+1,i  ) == 0) .OR. &
            (ice%mask_ocean_Aa( j  ,i-1) == 1 .AND. ice%mask_ice_Aa( j  ,i-1) == 0) .OR. &
            (ice%mask_ocean_Aa( j  ,i+1) == 1 .AND. ice%mask_ice_Aa( j  ,i+1) == 0)) THEN
            ice%mask_Aa(        j,i) = C%type_calvingfront
            ice%mask_cf_Aa( j,i) = 1
        END IF
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Then on the Ac grids
    ! ====================
  
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      
      IF (ice%mask_sheet_Aa( j,i) == 1 .OR.  ice%mask_sheet_Aa( j,i+1) == 1) ice%mask_sheet_Acx(  j,i) = 1
      IF (ice%mask_shelf_Aa( j,i) == 1 .OR.  ice%mask_shelf_Aa( j,i+1) == 1) ice%mask_shelf_Acx(  j,i) = 1
      IF (ice%mask_ocean_Aa( j,i) == 1 .OR.  ice%mask_ocean_Aa( j,i+1) == 1) ice%mask_ocean_Acx(  j,i) = 1
      IF (ice%mask_sheet_Aa( j,i) == 1 .AND. ice%mask_shelf_Aa( j,i+1) == 1) ice%mask_gl_Acx(     j,i) = 1
      IF (ice%mask_shelf_Aa( j,i) == 1 .AND. ice%mask_sheet_Aa( j,i+1) == 1) ice%mask_gl_Acx(     j,i) = 1
      IF (ice%mask_ice_Aa(   j,i) == 1 .AND. ice%mask_ice_Aa(   j,i+1) == 0) ice%mask_margin_Acx( j,i) = 1
      IF (ice%mask_ice_Aa(   j,i) == 0 .AND. ice%mask_ice_Aa(   j,i+1) == 1) ice%mask_margin_Acx( j,i) = 1
      IF (ice%mask_ice_Aa(   j,i) == 1 .AND. ice%mask_ocean_Aa( j,i+1) == 1) ice%mask_cf_Acx(     j,i) = 1
      IF (ice%mask_ocean_Aa( j,i) == 1 .AND. ice%mask_ice_Aa(   j,i+1) == 1) ice%mask_cf_Acx(     j,i) = 1
      
    END DO
    END DO
    CALL sync
  
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      
      IF (ice%mask_sheet_Aa( j,i) == 1 .OR.  ice%mask_sheet_Aa( j+1,i) == 1) ice%mask_sheet_Acy(  j,i) = 1
      IF (ice%mask_shelf_Aa( j,i) == 1 .OR.  ice%mask_shelf_Aa( j+1,i) == 1) ice%mask_shelf_Acy(  j,i) = 1
      IF (ice%mask_ocean_Aa( j,i) == 1 .OR.  ice%mask_ocean_Aa( j+1,i) == 1) ice%mask_ocean_Acy(  j,i) = 1
      IF (ice%mask_sheet_Aa( j,i) == 1 .AND. ice%mask_shelf_Aa( j+1,i) == 1) ice%mask_gl_Acy(     j,i) = 1
      IF (ice%mask_shelf_Aa( j,i) == 1 .AND. ice%mask_sheet_Aa( j+1,i) == 1) ice%mask_gl_Acy(     j,i) = 1
      IF (ice%mask_ice_Aa(   j,i) == 1 .AND. ice%mask_ice_Aa(   j+1,i) == 0) ice%mask_margin_Acy( j,i) = 1
      IF (ice%mask_ice_Aa(   j,i) == 0 .AND. ice%mask_ice_Aa(   j+1,i) == 1) ice%mask_margin_Acy( j,i) = 1
      IF (ice%mask_ice_Aa(   j,i) == 1 .AND. ice%mask_ocean_Aa( j+1,i) == 1) ice%mask_cf_Acy(     j,i) = 1
      IF (ice%mask_ocean_Aa( j,i) == 1 .AND. ice%mask_ice_Aa(   j+1,i) == 1) ice%mask_cf_Acy(     j,i) = 1
      
    END DO
    END DO
    CALL sync
  
  END SUBROUTINE determine_masks
  SUBROUTINE ice_physical_properties( grid, ice, time)
    ! Calculate the pressure melting point, flow parameter, specific heat and thermal conductivity of the ice.
      
    USE parameters_module, ONLY: T0, CC, n_flow, sec_per_year
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: time
  
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: Ti_mean     ! Mean ice temperature at the shelf [K]
    REAL(dp), DIMENSION(C%nZ)                          :: prof
    
    REAL(dp)                                           :: A_flow_MISMIP
    
    REAL(dp), PARAMETER                                :: A_low_temp  = 1.14E-05_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: A_high_temp = 5.47E+10_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_low_temp  = 6.0E+04_dp    ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_high_temp = 13.9E+04_dp   ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: R_gas       = 8.314_dp      ! Gas constant [J mol^-1 K^-1]
    
    ! If we're doing one of the EISMINT experiments, use fixed values
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler') THEN
        
        ice%A_flow_Aa(      :,:,grid%i1:grid%i2) = 1.0E-16_dp
        ice%A_flow_mean_Aa(   :,grid%i1:grid%i2) = 1.0E-16_dp
        ice%Ki_Aa(          :,:,grid%i1:grid%i2) = 2.1_dp * sec_per_year
        ice%Cpi_Aa(         :,:,grid%i1:grid%i2) = 2009._dp
        
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
        DO k = 1, C%nZ
          ice%Ti_pmp_Aa( k,j,i) = T0 - (C%zeta(k) * ice%Hi_Aa( j,i) * 8.7E-04_dp)
        END DO
        END DO
        END DO
        CALL sync
        
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
      
        A_flow_MISMIP = 1.0E-16_dp
        IF     (time < 25000._dp) THEN
          A_flow_MISMIP = 1.0E-16_dp
        ELSEIF (time < 50000._dp) THEN
          A_flow_MISMIP = 1.0E-17_dp
        ELSEIF (time < 75000._dp) THEN
          A_flow_MISMIP = 1.0E-16_dp
        END IF
        
        ice%A_flow_Aa(      :,:,grid%i1:grid%i2) = A_flow_MISMIP
        ice%A_flow_mean_Aa(   :,grid%i1:grid%i2) = A_flow_MISMIP
        CALL sync
        
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
      
        ice%A_flow_Aa(      :,:,grid%i1:grid%i2) = (3.7E8_dp ** (-n_flow)) * sec_per_year
        ice%A_flow_mean_Aa(   :,grid%i1:grid%i2) = (3.7E8_dp ** (-n_flow)) * sec_per_year
        CALL sync
        
        RETURN
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in ice_physical_properties!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
  
    ! First on the Aa grid
    ! ====================
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Calculate the pressure melting point temperature (= the maximum temperature) for each depth (see equation (11.2)):
      ice%Ti_pmp_Aa( :,j,i) = T0 - CC * ice%Hi_Aa( j,i) * C%zeta
      
      DO k = 1, C%nZ 

        ! Calculation of the flow parameter at the sheet and groundline as a function of the ice temperature 
        ! the Arrhenius relationship (see equation (11.10), Huybrechts (4.6)):
        IF (ice%Ti_Aa( k,j,i) < 263.15_dp) THEN
          ice%A_flow_Aa( k,j,i) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * ice%Ti_Aa( k,j,i)))  
        ELSE
          ice%A_flow_Aa( k,j,i) = A_high_temp * EXP(-Q_high_temp / (R_gas * ice%Ti_Aa( k,j,i)))  
        END IF
           
        ! Calculation of the parameterization of the specific heat capacity of ice, based on Pounder (1965):
        ice%Cpi_Aa( k,j,i) = 2115.3_dp + 7.79293_dp * (ice%Ti_Aa( k,j,i) - T0)  ! See equation (11.9)
           
        ! Calculation of the parameterization of the thermal conductivity of ice, based on Ritz (1987):
        ice%Ki_Aa( k,j,i)  = 3.101E+08_dp * EXP(-0.0057_dp * ice%Ti_Aa( k,j,i)) ! See equation (11.5), Huybrechts (4.40)
      END DO 

      IF (ice%mask_sheet_Aa( j,i) == 1) THEN
        ! Calculation of the vertical average flow parameter at the sheet and groundline
        prof = ice%A_flow_Aa( :,j,i)
        ice%A_flow_mean_Aa( j,i) = vertical_average( prof)
      ELSE
        ! Calculation of the flow parameter at the shelf as a function of the ice temperature 
        ! the Arrhenius relationship (see equation (11.10), Huybrechts (4.6)):
        prof = ice%Ti_Aa( :,j,i)
        Ti_mean = vertical_average( prof)
        IF (Ti_mean < 263.15_dp) THEN
          ice%A_flow_mean_Aa( j,i) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * Ti_mean))  
        ELSE
          ice%A_flow_mean_Aa( j,i) = A_high_temp * EXP(-Q_high_temp / (R_gas * Ti_mean))  
        END IF
      END IF  
          
      IF (ice%A_flow_mean_Aa( j,i) /= ice%A_flow_mean_Aa( j,i)) THEN
        WRITE(0,*) ' NaN values in A_flow_mean_Aa (ice_physical_properties)!'
        WRITE(0,*) '   A_flow_mean = ', ice%A_flow_mean_Aa( j,i), ', Ti = ', ice%Ti_Aa( :,j,i)
        CALl MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ice_physical_properties 
  
  FUNCTION is_floating( Hi, Hb, SL) RESULT( isso)
    ! The flotation criterion
      
    IMPLICIT NONE
    
    REAL(dp),                            INTENT(IN)    :: Hi, Hb, SL
    LOGICAL                                            :: isso
    
    isso = .FALSE.
    IF (Hi < (SL - Hb) * seawater_density/ice_density) isso = .TRUE.
    
  END FUNCTION is_floating
  SUBROUTINE determine_grounded_fraction( grid, ice)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER                :: T, T_B
    INTEGER                                            :: wT, wT_B
    REAL(dp)                                           :: T_nw, T_ne, T_sw, T_se
    REAL(dp)                                           :: xr_n, xr_s, yr_w, yr_e
    REAL(dp)                                           :: A_gr
    
    ! Calculate thickness above flotation
    CALL allocate_shared_dp_2D( grid%ny  , grid%nx  , T  , wT  )
    CALL allocate_shared_dp_2D( grid%ny-1, grid%nx-1, T_B, wT_B)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      T( j,i)  = ice%Hi_Aa( j,i) - ((ice%SL_Aa( j,i) - ice%Hb_Aa( j,i)) * (seawater_density / ice_density))
    END DO
    END DO
    CALL sync
    CALL map_Aa_to_Ab_2D( grid, T, T_B)

    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1

      ice%grounded_fraction( j,i) = 0._dp

      IF ((T_B( j-1,i-1) <  0._dp .OR. T_B( j-1,i  ) <  0._dp .OR. T_B( j  ,i-1) <  0._dp .OR. T_B( j  ,i  ) <  0._dp) .AND. &
          (T_B( j-1,i-1) >= 0._dp .OR. T_B( j-1,i  ) >= 0._dp .OR. T_B( j  ,i-1) >= 0._dp .OR. T_B( j  ,i  ) >= 0._dp)) THEN

        T_nw = T_B( j  ,i-1)
        T_ne = T_B( j  ,i  )
        T_sw = T_B( j-1,i-1)
        T_se = T_B( j-1,i  )

        xr_n = 0._dp
        xr_s = 0._dp
        yr_w = 0._dp
        yr_e = 0._dp
        
        IF (T_nw * T_ne < 0._dp) THEN
          xr_n = MAX(0._dp, MIN(1._dp, T_nw / (T_nw - T_ne) ))
        END IF
        IF (T_sw * T_se < 0._dp) THEN
          xr_s = MAX(0._dp, MIN(1._dp, T_sw / (T_sw - T_se) ))
        END IF
        IF (T_sw * T_nw < 0._dp) THEN
          yr_w = MAX(0._dp, MIN(1._dp, T_sw / (T_sw - T_nw) ))
        END IF
        IF (T_se * T_ne < 0._dp) THEN
          yr_e = MAX(0._dp, MIN(1._dp, T_se / (T_se - T_ne) ))
        END IF

        A_gr = 0._dp

        IF (T_nw >= 0._dp) THEN
          ! NW is grounded

          IF (T_ne < 0._dp .AND. T_sw < 0._dp .AND. T_se < 0._dp) THEN
            ! Only NW is grounded
            A_gr = 0.5_dp * (xr_n * (1._dp - yr_w))
          ELSEIF (T_ne < 0._dp .AND. T_sw >= 0._dp .AND. T_se < 0._dp) THEN
            ! Both NW and SW are grounded
            A_gr = 0.5_dp * (xr_n + xr_s)
          ELSEIF (T_ne >= 0._dp .AND. T_sw < 0._dp .AND. T_se < 0._dp) THEN
            ! Both NW and NE are grounded
            A_gr =  0.5_dp * ((1._dp - yr_w) + (1._dp - yr_e))
          ELSEIF (T_ne >= 0._dp .AND. T_sw >= 0._dp .AND. T_se < 0._dp) THEN
            ! Only SE is floating
            A_gr = 1._dp - 0.5_dp * ((1._dp - xr_s) * yr_e)
          ELSEIF (T_ne >= 0._dp .AND. T_sw < 0._dp .AND. T_se >= 0._dp) THEN
            ! Only SW is floating
            A_gr = 1._dp - 0.5_dp * (xr_s * yr_w)
          ELSEIF (T_ne < 0._dp .AND. T_sw >= 0._dp .AND. T_se >= 0._dp) THEN
            ! Only NE is floating
            A_gr = 1._dp - 0.5_dp * ((1._dp - xr_n) * (1._dp - yr_e))
          ELSEIF (T_ne < 0._dp .AND. T_sw < 0._dp .AND. T_se >= 0._dp) THEN
            ! Both NE and SW are floating
            A_gr = 1._dp - 0.5_dp * (xr_s * yr_w) - 0.5_dp * ((1._dp - xr_n) * (1._dp - yr_e))
          END IF

        ELSE ! IF (T_nw >= 0._dp)
          ! NW is floating

          IF (T_ne >= 0._dp .AND. T_sw >= 0._dp .AND. T_se >= 0._dp) THEN
            ! Only NW is floating
            A_gr = 1._dp - 0.5_dp * (xr_n * (1._dp - yr_w))
          ELSEIF (T_ne >= 0._dp .AND. T_sw < 0._dp .AND. T_se >= 0._dp) THEN
            ! Both NW and SW are floating
            A_gr = 0.5_dp * ((1._dp - xr_n) + (1._dp - xr_s))
          ELSEIF (T_ne < 0._dp .AND. T_sw >= 0._dp .AND. T_se >= 0._dp) THEN
            ! Both NW and NE are floating
            A_gr = 0.5_dp * (yr_w + yr_e)
          ELSEIF (T_ne < 0._dp .AND. T_sw < 0._dp .AND. T_se >= 0._dp) THEN
            ! Only SE is grounded
            A_gr = 0.5_dp * (1._dp - xr_s) * yr_e
          ELSEIF (T_ne < 0._dp .AND. T_sw >= 0._dp .AND. T_se < 0._dp) THEN
            ! Only SW is grounded
            A_gr = 0.5_dp * (xr_s * yr_w)
          ELSEIF (T_ne >= 0._dp .AND. T_se < 0._dp .AND. T_se < 0._dp) THEN
            ! Only NE is grounded
            A_gr = 0.5_dp * ((1._dp - xr_n) * (1._dp - yr_e))
          ELSEIF (T_ne >= 0._dp .AND. T_sw >= 0._dp .AND. T_se < 0._dp) THEN
            ! Both NE and SW are grounded
            A_gr = 0.5_dp * (xr_s * yr_w) + 0.5_dp * ((1._dp - xr_n) * (1._dp - yr_e))
          END IF

        END IF

      ELSEIF (T( j,i) >= 0._dp) THEN
        ! This entire grid cell is grounded
        A_gr = 1._dp

      ELSE
        ! This entire grid cell is floating
        A_gr = 0._dp

      END IF ! IF ((T_B( j-1,i-1) <  0._dp .OR. T_B( j-1,i  ) <  0._dp .OR. T_B( j  ,i-1) <  0._dp .OR. T_B( j  ,i  ) <  0._dp) .AND. &
             !     (T_B( j-1,i-1) >= 0._dp .OR. T_B( j-1,i  ) >= 0._dp .OR. T_B( j  ,i-1) >= 0._dp .OR. T_B( j  ,i  ) >= 0._dp))

      ice%grounded_fraction( j,i) = A_gr

    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wT)
    CALL deallocate_shared( wT_B)
    
  END SUBROUTINE determine_grounded_fraction
  SUBROUTINE ocean_floodfill( grid, ice)
    ! Use a simple floodfill algorithm to determine the ocean mask,
    ! to prevent the formation of (pro-/sub-glacial) lakes
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(INOUT) :: ice 
  
    INTEGER                                            :: i,j,ii,jj
    INTEGER, DIMENSION( :,:  ), POINTER                :: map
    INTEGER, DIMENSION( :,:  ), POINTER                :: stack
    INTEGER                                            :: wmap, wstack
    INTEGER                                            :: stackN
    
    ! Allocate shared memory for the map, because even though the flood-fill is done only
    ! by the Master, the other processes need access to the resulting filled map.
    CALL allocate_shared_int_2D( grid%ny, grid%nx,   map,   wmap  )
    CALL allocate_shared_int_2D( grid%ny*grid%nx, 2, stack, wstack)
    
    ! No easy way to parallelise flood-fill, just let the Master do it
    IF (par%master) THEN
    
      ice%mask_land_Aa  = 1
      ice%mask_ocean_Aa = 0
      ice%mask_Aa       = C%type_land
      
      map    = 0
      stack  = 0
      stackN = 0
      
      ! Let the ocean flow in from the domain edges
      DO i = 1, grid%nx
        j = 1
        stackN = stackN + 1
        stack( stackN,:) = [i,j]
        map( j,i) = 1
        
        j = grid%ny
        stackN = stackN + 1
        stack( stackN,:) = [i,j]
        map( j,i) = 1
      END DO
      DO j = 1, grid%ny
        i = 1
        stackN = stackN + 1
        stack( stackN,:) = [i,j]
        map( j,i) = 1
        
        i = grid%nx
        stackN = stackN + 1
        stack( stackN,:) = [i,j]
        map( j,i) = 1
      END DO
      
      ! Flood fill
      DO WHILE (stackN > 0)
        
        ! Inspect the last element of the stack
        i = stack( stackN,1)
        j = stack( stackN,2)
        
        ! Remove it from the stack
        stack( stackN,:) = [0,0]
        stackN = stackN-1
        
        IF (is_floating( ice%Hi_Aa( j,i), ice%Hb_Aa( j,i), ice%SL_Aa( j,i))) THEN
          ! This element is ocean
          
          ! Mark it as such on the map
          map( j,i) = 2
          ice%mask_ocean_Aa( j,i) = 1
          ice%mask_land_Aa(  j,i) = 0
          ice%mask_Aa(       j,i) = C%type_ocean
          
          ! Add its (valid) neighbours to the stack
          DO ii = i-1,i+1
          DO jj = j-1,j+1
            IF (ii>=1 .AND. ii<= grid%nx .AND. jj>=1 .AND. jj<=grid%ny .AND. (.NOT. (ii==i .AND. jj==j))) THEN
              ! This neighbour exists
              
              IF (map( jj,ii) == 0) THEN
                ! This neighbour isn't yet mapped or stacked
                
                map( jj,ii) = 1
                stackN = stackN + 1
                stack( stackN,:) = [ii,jj]
                
              END IF
              
            END IF
          END DO
          END DO
          
        ELSE ! IF (is_floating( ice%Hi_Aa( j,i), ice%Hb_Aa( j,i), ice%SL_Aa( j,i))) THEN
          ! This element is land
        END IF ! IF (is_floating( ice%Hi_Aa( j,i), ice%Hb_Aa( j,i), ice%SL_Aa( j,i))) THEN
        
      END DO ! DO WHILE (stackN > 0)
    
    END IF ! IF (par%master) THEN
    CALL sync
      
    ! Clean up after yourself
    CALL deallocate_shared( wmap)
    CALL deallocate_shared( wstack)
    
  END SUBROUTINE ocean_floodfill

END MODULE general_ice_model_data_module
