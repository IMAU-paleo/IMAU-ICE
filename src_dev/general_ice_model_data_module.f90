MODULE general_ice_model_data_module
  ! Contains routines for calculating "general" ice model data (surface slopes, staggered
  ! data fields, masks, physical properties, etc.)

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE data_types_module,               ONLY: type_grid, type_ice_model
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE parameters_module
  USE utilities_module,                ONLY: vertical_average
  USE derivatives_and_grids_module,    ONLY: map_a_to_cx_2D, map_a_to_cy_2D, ddx_a_to_cx_2D, ddy_a_to_cy_2D, &
                                             ddy_a_to_cx_2D, ddx_a_to_cy_2D, map_a_to_cx_3D, map_a_to_cy_3D, &
                                             ddx_a_to_a_2D, ddy_a_to_a_2D, map_a_to_b_2D
  USE basal_conditions_module,         ONLY: calc_basal_yield_stress

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
      ice%Hs_a( j,i) = ice%Hi_a( j,i) + MAX( ice%SL_a( j,i) - ice_density / seawater_density * ice%Hi_a( j,i), ice%Hb_a( j,i))
    END DO
    END DO
    CALL sync
    
    ! Calculate thickness above flotation
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%TAF_a( j,i) = ice%Hi_a( j,i) - MAX(0._dp, (ice%SL_a( j,i) - ice%Hb_a( j,i)) * (seawater_density / ice_density))
    END DO
    END DO
    CALL sync
    CALL map_a_to_b_2D( grid, ice%TAF_a, ice%TAF_b)
    
    ! Determine masks
    CALL determine_masks(                    grid, ice)
    CALL determine_grounded_fraction(        grid, ice)
    CALL determine_floating_margin_fraction( grid, ice)
 
    ! Get ice thickness, flow factor, and surface slopes on the required grids
    CALL map_a_to_cx_2D( grid, ice%Hi_a,          ice%Hi_cx    )
    CALL map_a_to_cy_2D( grid, ice%Hi_a,          ice%Hi_cy    )
    CALL map_a_to_b_2D(  grid, ice%Hi_a,          ice%Hi_b     )
    CALL map_a_to_cx_3D( grid, ice%A_flow_3D_a,   ice%A_flow_3D_cx)
    CALL map_a_to_cy_3D( grid, ice%A_flow_3D_a,   ice%A_flow_3D_cy)
    CALL map_a_to_cx_2D( grid, ice%A_flow_vav_a,  ice%A_flow_vav_cx)
    CALL map_a_to_cy_2D( grid, ice%A_flow_vav_a,  ice%A_flow_vav_cy)
    
    CALL ddx_a_to_a_2D(  grid, ice%Hi_a,          ice%dHi_dx_a )
    CALL ddy_a_to_a_2D(  grid, ice%Hi_a,          ice%dHi_dy_a )
    
    CALL ddx_a_to_a_2D(  grid, ice%Hs_a,          ice%dHs_dx_a )
    CALL ddy_a_to_a_2D(  grid, ice%Hs_a,          ice%dHs_dy_a )
    CALL ddx_a_to_cx_2D( grid, ice%Hs_a,          ice%dHs_dx_cx)
    CALL ddy_a_to_cx_2D( grid, ice%Hs_a,          ice%dHs_dy_cx)
    CALL ddx_a_to_cy_2D( grid, ice%Hs_a,          ice%dHs_dx_cy)
    CALL ddy_a_to_cy_2D( grid, ice%Hs_a,          ice%dHs_dy_cy)
            
    ! Calculate physical properties (flow factor, heat capacity, thermal conductivity, pressure melting point)
    CALL ice_physical_properties( grid, ice, time)
    
    ! Calculate the basal yield stress tau_c
    CALL calc_basal_yield_stress( grid, ice)
    
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
      ice%mask_land_a    = 1
      ice%mask_ocean_a   = 0
      ice%mask_lake_a    = 0
      ice%mask_ice_a     = 0
      ice%mask_sheet_a   = 0
      ice%mask_shelf_a   = 0
      ice%mask_coast_a   = 0
      ice%mask_coast_cx  = 0
      ice%mask_coast_cy  = 0
      ice%mask_margin_a  = 0
      ice%mask_margin_cx = 0
      ice%mask_margin_cy = 0
      ice%mask_gl_a      = 0
      ice%mask_gl_cx     = 0
      ice%mask_gl_cy     = 0
      ice%mask_cf_a      = 0
      ice%mask_cf_cx     = 0
      ice%mask_cf_cy     = 0
      ice%mask_a         = C%type_land
    END IF
    CALL sync
  
    ! First on the Aa grid
    ! ====================
    
    IF (C%do_ocean_floodfill) THEN
      CALL ocean_floodfill( grid, ice)
    ELSE
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (is_floating( ice%Hi_a( j,i), ice%Hb_a( j,i), ice%SL_a( j,i))) THEN
          ice%mask_land_a(  j,i) = 0
          ice%mask_ocean_a( j,i) = 1
          ice%mask_a(       j,i) = C%type_ocean
        END IF
      END DO
      END DO
      CALL sync
    
    END IF ! IF (C%do_ocean_floodfill) THEN
  
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Ice
      IF (ice%Hi_a( j,i) > 0._dp) THEN
        ice%mask_ice_a(  j,i) = 1
      END IF
      
      ! Sheet
      IF (ice%mask_ice_a( j,i) == 1 .AND. ice%mask_land_a( j,i) == 1) THEN
        ice%mask_sheet_a( j,i) = 1
        ice%mask_a(       j,i) = C%type_sheet
      END IF
    
      ! Shelf
      IF (ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ocean_a( j,i) == 1) THEN
        ice%mask_shelf_a( j,i) = 1
        ice%mask_a(       j,i) = C%type_shelf
      END IF
      
    END DO
    END DO
    CALL sync
  
    ! Determine coast, grounding line and calving front
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    
      ! Ice-free land bordering ocean equals coastline
      IF (ice%mask_land_a( j,i) == 1 .AND. ice%mask_ice_a( j,i) == 0) THEN  
        IF (ice%mask_ocean_a( j-1,i  ) == 1 .OR. &
            ice%mask_ocean_a( j+1,i  ) == 1 .OR. &
            ice%mask_ocean_a( j  ,i-1) == 1 .OR. &
            ice%mask_ocean_a( j  ,i+1) == 1) THEN
            ice%mask_a(       j,i) = C%type_coast
            ice%mask_coast_a( j,i) = 1
        END IF
      END IF
    
      ! Ice bordering non-ice equals margin
      IF (ice%mask_ice_a( j,i) == 1) THEN  
        IF (ice%mask_ice_a( j-1,i  ) == 0 .OR. &
            ice%mask_ice_a( j+1,i  ) == 0 .OR. &
            ice%mask_ice_a( j  ,i-1) == 0 .OR. &
            ice%mask_ice_a( j  ,i+1) == 0) THEN
            ice%mask_a(        j,i) = C%type_margin
            ice%mask_margin_a( j,i) = 1
        END IF
      END IF
    
      ! Sheet bordering shelf equals grounding line
      IF (ice%mask_sheet_a( j,i) == 1) THEN  
        IF (ice%mask_shelf_a( j-1,i  ) == 1 .OR. &
            ice%mask_shelf_a( j+1,i  ) == 1 .OR. &
            ice%mask_shelf_a( j  ,i-1) == 1 .OR. &
            ice%mask_shelf_a( j  ,i+1) == 1) THEN
            ice%mask_a(    j,i) = C%type_groundingline
            ice%mask_gl_a( j,i) = 1
        END IF
      END IF
  
      ! Ice (sheet or shelf) bordering open ocean equals calvingfront
      IF (ice%mask_ice_a( j,i) == 1) THEN  
        IF ((ice%mask_ocean_a( j-1,i  ) == 1 .AND. ice%mask_ice_a( j-1,i  ) == 0) .OR. &
            (ice%mask_ocean_a( j+1,i  ) == 1 .AND. ice%mask_ice_a( j+1,i  ) == 0) .OR. &
            (ice%mask_ocean_a( j  ,i-1) == 1 .AND. ice%mask_ice_a( j  ,i-1) == 0) .OR. &
            (ice%mask_ocean_a( j  ,i+1) == 1 .AND. ice%mask_ice_a( j  ,i+1) == 0)) THEN
            ice%mask_a(        j,i) = C%type_calvingfront
            ice%mask_cf_a( j,i) = 1
        END IF
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Then on the Ac grids
    ! ====================
  
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      
      IF (ice%mask_sheet_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i+1) == 1) ice%mask_gl_cx(     j,i) = 1
      IF (ice%mask_shelf_a( j,i) == 1 .AND. ice%mask_sheet_a( j,i+1) == 1) ice%mask_gl_cx(     j,i) = 1
      IF (ice%mask_ice_a(   j,i) == 1 .AND. ice%mask_ice_a(   j,i+1) == 0) ice%mask_margin_cx( j,i) = 1
      IF (ice%mask_ice_a(   j,i) == 0 .AND. ice%mask_ice_a(   j,i+1) == 1) ice%mask_margin_cx( j,i) = 1
      IF (ice%mask_ice_a(   j,i) == 1 .AND. ice%mask_ocean_a( j,i+1) == 1) ice%mask_cf_cx(     j,i) = 1
      IF (ice%mask_ocean_a( j,i) == 1 .AND. ice%mask_ice_a(   j,i+1) == 1) ice%mask_cf_cx(     j,i) = 1
      
    END DO
    END DO
    CALL sync
  
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      
      IF (ice%mask_sheet_a( j,i) == 1 .AND. ice%mask_shelf_a( j+1,i) == 1) ice%mask_gl_cy(     j,i) = 1
      IF (ice%mask_shelf_a( j,i) == 1 .AND. ice%mask_sheet_a( j+1,i) == 1) ice%mask_gl_cy(     j,i) = 1
      IF (ice%mask_ice_a(   j,i) == 1 .AND. ice%mask_ice_a(   j+1,i) == 0) ice%mask_margin_cy( j,i) = 1
      IF (ice%mask_ice_a(   j,i) == 0 .AND. ice%mask_ice_a(   j+1,i) == 1) ice%mask_margin_cy( j,i) = 1
      IF (ice%mask_ice_a(   j,i) == 1 .AND. ice%mask_ocean_a( j+1,i) == 1) ice%mask_cf_cy(     j,i) = 1
      IF (ice%mask_ocean_a( j,i) == 1 .AND. ice%mask_ice_a(   j+1,i) == 1) ice%mask_cf_cy(     j,i) = 1
      
    END DO
    END DO
    CALL sync
  
  END SUBROUTINE determine_masks
  SUBROUTINE ice_physical_properties( grid, ice, time)
    ! Calculate the flow factor, pressure melting point, specific heat, and thermal conductivity of the ice.
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: time
  
    ! Local variables:
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
      IF (C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_E') THEN
        
        ice%A_flow_3D_a(  :,:,grid%i1:grid%i2) = 1.0E-16_dp
        ice%A_flow_vav_a(   :,grid%i1:grid%i2) = 1.0E-16_dp
        ice%Ki_a(         :,:,grid%i1:grid%i2) = 2.1_dp * sec_per_year
        ice%Cpi_a(        :,:,grid%i1:grid%i2) = 2009._dp
        
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
        DO k = 1, C%nZ
          ice%Ti_pmp_a( k,j,i) = T0 - (C%zeta(k) * ice%Hi_a( j,i) * 8.7E-04_dp)
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
        
        ice%A_flow_3D_a(  :,:,grid%i1:grid%i2) = A_flow_MISMIP
        ice%A_flow_vav_a(   :,grid%i1:grid%i2) = A_flow_MISMIP
        CALL sync
        
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
      
        ice%A_flow_3D_a(  :,:,grid%i1:grid%i2) = (3.7E8_dp ** (-C%n_flow)) * sec_per_year
        ice%A_flow_vav_a(   :,grid%i1:grid%i2) = (3.7E8_dp ** (-C%n_flow)) * sec_per_year
        CALL sync
        
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
      
        ice%A_flow_3D_a(  :,:,grid%i1:grid%i2) = 2.140373E-7_dp
        ice%A_flow_vav_a(   :,grid%i1:grid%i2) = 2.140373E-7_dp
        CALL sync
        
        RETURN
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in ice_physical_properties!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Calculate the pressure melting point temperature (= the maximum temperature) for each depth (see equation (11.2)):
      ice%Ti_pmp_a( :,j,i) = T0 - CC * ice%Hi_a( j,i) * C%zeta
      
      DO k = 1, C%nZ 

        ! Calculation of the flow parameter at the sheet and groundline as a function of the ice temperature 
        ! the Arrhenius relationship (see equation (11.10), Huybrechts (4.6)):
        IF (ice%mask_ice_a( j,i) == 1) THEN
          IF (ice%Ti_a( k,j,i) < 263.15_dp) THEN
            ice%A_flow_3D_a( k,j,i) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * ice%Ti_a( k,j,i)))  
          ELSE
            ice%A_flow_3D_a( k,j,i) = A_high_temp * EXP(-Q_high_temp / (R_gas * ice%Ti_a( k,j,i)))  
          END IF
        ELSE
          IF (C%DIVA_choice_ice_margin == 'BC') THEN
            ice%A_flow_3D_a( k,j,i) = 0._dp
          ELSEIF (C%DIVA_choice_ice_margin == 'infinite_slab') THEN
            ! In the "infinite slab" case, calculate effective viscosity everywhere
            ! (even when there's technically no ice present)
            ice%A_flow_3D_a( k,j,i) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * 263.15_dp))  
          ELSE
            IF (par%master) WRITE(0,*) '  ERROR: DIVA_choice_ice_margin "', TRIM(C%DIVA_choice_ice_margin), '" not implemented in calc_effective_viscosity!'
            CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
          END IF
        END IF
           
        ! Calculation of the parameterization of the specific heat capacity of ice, based on Pounder (1965):
        ice%Cpi_a( k,j,i) = 2115.3_dp + 7.79293_dp * (ice%Ti_a( k,j,i) - T0)  ! See equation (11.9)
           
        ! Calculation of the parameterization of the thermal conductivity of ice, based on Ritz (1987):
        ice%Ki_a( k,j,i)  = 3.101E+08_dp * EXP(-0.0057_dp * ice%Ti_a( k,j,i)) ! See equation (11.5), Huybrechts (4.40)
      END DO 

      IF (ice%mask_sheet_a( j,i) == 1) THEN
        ! Calculation of the vertical average flow parameter at the sheet and groundline
        prof = ice%A_flow_3D_a( :,j,i)
        ice%A_flow_vav_a( j,i) = vertical_average( prof)
      ELSE
        ! Calculation of the flow parameter at the shelf as a function of the ice temperature 
        ! the Arrhenius relationship (see equation (11.10), Huybrechts (4.6)):
        prof = ice%Ti_a( :,j,i)
        Ti_mean = vertical_average( prof)
        IF (Ti_mean < 263.15_dp) THEN
          ice%A_flow_vav_a( j,i) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * Ti_mean))  
        ELSE
          ice%A_flow_vav_a( j,i) = A_high_temp * EXP(-Q_high_temp / (R_gas * Ti_mean))  
        END IF
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
    ! Determine the grounded fraction of next-to-grounding-line pixels
    ! (used for determining basal friction in the DIVA)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: T_nw, T_ne, T_sw, T_se
    REAL(dp)                                           :: xr_n, xr_s, yr_w, yr_e
    REAL(dp)                                           :: A_gr

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ice%f_grnd_a( j,i) = 0._dp
      
      IF (ice%Hi_a( j,i) == 0._dp) CYCLE

      ! Get TAF at surrounding b-grid cells
      IF (i == 1) THEN
        IF (j == 1) THEN
          T_nw = 0.5_dp * (ice%TAF_a( 1,1) + ice%TAF_a( 2,1))
          T_ne = ice%TAF_b( 1,1)
          T_sw = ice%TAF_a( 1,1)
          T_se = 0.5_dp * (ice%TAF_a( 1,1) + ice%TAF_a( 1,2))
        ELSEIF (j == grid%ny) THEN
          T_nw = ice%TAF_a( grid%ny,1)
          T_ne = 0.5_dp * (ice%TAF_a( grid%ny-1,1) + ice%TAF_a( grid%ny-1,2))
          T_sw = 0.5_dp * (ice%TAF_a( grid%ny-1,1) + ice%TAF_a( grid%ny-2,1))
          T_se = ice%TAF_b( grid%ny-1,1)
        ELSE
          T_nw = 0.5_dp * (ice%TAF_a( j,1) + ice%TAF_a( j+1,1))
          T_ne = ice%TAF_b( j,1)
          T_sw = 0.5_dp * (ice%TAF_a( j,1) + ice%TAF_a( j-1,1))
          T_se = ice%TAF_b( j-1,1)
        END IF
      ELSEIF (i == grid%nx) THEN
        IF (j == 1) THEN
          T_nw = ice%TAF_b( 1,grid%nx-1)
          T_ne = 0.5_dp * (ice%TAF_a( 1,grid%nx) + ice%TAF_a( 2,grid%nx))
          T_sw = 0.5_dp * (ice%TAF_a( 1,grid%nx) + ice%TAF_a( 1,grid%nx-1))
          T_se = ice%TAF_a( 1,grid%nx)
        ELSEIF (j == grid%ny) THEN
          T_nw = 0.5_dp * (ice%TAF_a( grid%ny,grid%nx) + ice%TAF_a( grid%ny,grid%nx-1))
          T_ne = ice%TAF_a( grid%ny,grid%nx)
          T_sw = ice%TAF_b( grid%ny-1,grid%nx-1)
          T_se = 0.5_dp * (ice%TAF_a( grid%ny,grid%nx) + ice%TAF_a( grid%ny-1,grid%nx))
        ELSE
          T_nw = 0.5_dp * (ice%TAF_a( j,grid%nx) + ice%TAF_a( j+1,grid%nx))
          T_ne = ice%TAF_b( j,grid%nx-1)
          T_sw = 0.5_dp * (ice%TAF_a( j,grid%nx) + ice%TAF_a( j-1,grid%nx))
          T_se = ice%TAF_b( j-1,grid%nx-1)
        END IF
      ELSE
        IF (j == 1) THEN
          T_nw = ice%TAF_b( 1,i-1)
          T_ne = ice%TAF_b( 1,i)
          T_sw = 0.5_dp * (ice%TAF_a( 1,i) + ice%TAF_a( 1,i-1))
          T_se = 0.5_dp * (ice%TAF_a( 1,i) + ice%TAF_a( 1,i+1))
        ELSEIF (j == grid%ny) THEN
          T_nw = ice%TAF_b( grid%ny-1,i-1)
          T_ne = ice%TAF_b( grid%ny-1,i)
          T_sw = 0.5_dp * (ice%TAF_a( grid%ny,i) + ice%TAF_a( grid%ny,i-1))
          T_se = 0.5_dp * (ice%TAF_a( grid%ny,i) + ice%TAF_a( grid%ny,i+1))
        ELSE
          T_nw = ice%TAF_b( j  ,i-1)
          T_ne = ice%TAF_b( j  ,i  )
          T_sw = ice%TAF_b( j-1,i-1)
          T_se = ice%TAF_b( j-1,i  )
        END IF
      END IF
 
      IF ((T_nw <  0._dp .OR. T_ne <  0._dp .OR. T_sw <  0._dp .OR. T_se <  0._dp) .AND. &
          (T_nw >= 0._dp .OR. T_ne >= 0._dp .OR. T_sw >= 0._dp .OR. T_se >= 0._dp)) THEN
        ! At least one of the four corners is floating and at least one is grounded;
        ! the grounded fraction is between 0 and 1 and should be calculated

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

      ELSEIF (ice%TAF_a( j,i) >= 0._dp) THEN
        ! This entire grid cell is grounded
        A_gr = 1._dp

      ELSE
        ! This entire grid cell is floating
        A_gr = 0._dp

      END IF ! IF ((T_sw <  0._dp .OR. T_se <  0._dp .OR. T_nw <  0._dp .OR. T_ne <  0._dp) .AND. &
             !     (T_sw >= 0._dp .OR. T_se >= 0._dp .OR. T_nw >= 0._dp .OR. T_ne >= 0._dp)) THEN

      ice%f_grnd_a( j,i) = A_gr

    END DO
    END DO
    CALL sync
    
    ! Get the grounded fraction on the staggered cx/cy grids
    CALL map_a_to_cx_2D( grid, ice%f_grnd_a, ice%f_grnd_cx)
    CALL map_a_to_cy_2D( grid, ice%f_grnd_a, ice%f_grnd_cy)
    
  END SUBROUTINE determine_grounded_fraction
  SUBROUTINE determine_floating_margin_fraction( grid, ice)
    ! Determine the ice-filled fraction of floating margin pixels
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j,ii,jj
    REAL(dp)                                           :: Hi_min

    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1

      ice%float_margin_frac_a( j,i) = 0._dp
      ice%Hi_actual_cf_a(      j,i) = 0._dp
      
      IF (ice%mask_cf_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 1) THEN
        ! Actual ice thickness is equal to that of the thinnest neighbouring non-margin grid cell
        
        Hi_min = 10000._dp
        
        DO ii = i-1,i+1
        DO jj = j-1,j+1
          IF (ice%mask_ice_a( jj,ii) == 1 .AND. ice%mask_cf_a( jj,ii) == 0) THEN
            Hi_min = MIN( Hi_min, ice%Hi_a( jj,ii))
          END IF
        END DO
        END DO
        
        ! Calculate ice-filled fraction
        IF (Hi_min < 10000._dp) THEN
          ice%float_margin_frac_a( j,i) = ice%Hi_a( j,i) / Hi_min
          ice%Hi_actual_cf_a(      j,i) = Hi_min
        ELSE
          ice%float_margin_frac_a( j,i) = 1._dp
          ice%Hi_actual_cf_a(      j,i) = ice%Hi_a( j,i)
        END IF
        
      END IF

    END DO
    END DO
    CALL sync
    
  END SUBROUTINE determine_floating_margin_fraction
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
    
      ice%mask_land_a  = 1
      ice%mask_ocean_a = 0
      ice%mask_a       = C%type_land
      
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
        
        IF (is_floating( ice%Hi_a( j,i), ice%Hb_a( j,i), ice%SL_a( j,i))) THEN
          ! This element is ocean
          
          ! Mark it as such on the map
          map( j,i) = 2
          ice%mask_ocean_a( j,i) = 1
          ice%mask_land_a(  j,i) = 0
          ice%mask_a(       j,i) = C%type_ocean
          
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
          
        ELSE ! IF (is_floating( ice%Hi_a( j,i), ice%Hb_a( j,i), ice%SL_a( j,i))) THEN
          ! This element is land
        END IF ! IF (is_floating( ice%Hi_a( j,i), ice%Hb_a( j,i), ice%SL_a( j,i))) THEN
        
      END DO ! DO WHILE (stackN > 0)
    
    END IF ! IF (par%master) THEN
    CALL sync
      
    ! Clean up after yourself
    CALL deallocate_shared( wmap)
    CALL deallocate_shared( wstack)
    
  END SUBROUTINE ocean_floodfill

END MODULE general_ice_model_data_module
