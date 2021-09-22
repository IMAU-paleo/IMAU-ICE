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
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             is_floating, surface_elevation, thickness_above_floatation, &
                                             vertical_average, interp_bilin_2D
  USE derivatives_and_grids_module,    ONLY: map_a_to_cx_2D, map_a_to_cy_2D, ddx_a_to_cx_2D, ddy_a_to_cy_2D, &
                                             ddy_a_to_cx_2D, ddx_a_to_cy_2D, map_a_to_cx_3D, map_a_to_cy_3D, &
                                             ddx_a_to_a_2D, ddy_a_to_a_2D, map_a_to_b_2D
  USE basal_conditions_module,         ONLY: calc_basal_conditions

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
    
    ! Calculate surface elevation and thickness above floatation
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%Hs_a( j,i) = surface_elevation( ice%Hi_a( j,i), ice%Hb_a( j,i), ice%SL_a( j,i))
      ice%TAF_a( j,i) = thickness_above_floatation( ice%Hi_a( j,i), ice%Hb_a( j,i), ice%SL_a( j,i))
    END DO
    END DO
    CALL sync
    
    ! Determine masks
    CALL determine_masks(                    grid, ice)
    CALL determine_grounded_fractions(       grid, ice)
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
    CALL calc_basal_conditions( grid, ice)
    
  END SUBROUTINE update_general_ice_model_data
  SUBROUTINE determine_masks( grid, ice)
    ! Determine the different masks, on both the Aa and the Ac grid
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(INOUT) :: ice 
  
    INTEGER                                            :: i,j
    
    ! Save the previous ice mask, for use in thermodynamics
    ice%mask_ice_a_prev( :,grid%i1:grid%i2) = ice%mask_ice_a( :,grid%i1:grid%i2)
    CALL sync
    
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
      
      IF (ice%mask_ice_a(   j,i  ) == 1 .AND. &
         (ice%mask_ocean_a( j,i+1) == 1 .AND. ice%mask_ice_a( j,i+1) == 0)) ice%mask_cf_cx(     j,i) = 1
      IF (ice%mask_ice_a(   j,i+1) == 1 .AND. &
         (ice%mask_ocean_a( j,i  ) == 1 .AND. ice%mask_ice_a( j,i  ) == 0)) ice%mask_cf_cx(     j,i) = 1
      
    END DO
    END DO
    CALL sync
  
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      
      IF (ice%mask_sheet_a( j,i) == 1 .AND. ice%mask_shelf_a( j+1,i) == 1) ice%mask_gl_cy(     j,i) = 1
      IF (ice%mask_shelf_a( j,i) == 1 .AND. ice%mask_sheet_a( j+1,i) == 1) ice%mask_gl_cy(     j,i) = 1
      IF (ice%mask_ice_a(   j,i) == 1 .AND. ice%mask_ice_a(   j+1,i) == 0) ice%mask_margin_cy( j,i) = 1
      IF (ice%mask_ice_a(   j,i) == 0 .AND. ice%mask_ice_a(   j+1,i) == 1) ice%mask_margin_cy( j,i) = 1
      
      IF (ice%mask_ice_a(   j  ,i) == 1 .AND. &
         (ice%mask_ocean_a( j+1,i) == 1 .AND. ice%mask_ice_a( j+1,i) == 0)) ice%mask_cf_cy(     j,i) = 1
      IF (ice%mask_ice_a(   j+1,i) == 1 .AND. &
         (ice%mask_ocean_a( j  ,i) == 1 .AND. ice%mask_ice_a( j  ,i) == 0)) ice%mask_cf_cy(     j,i) = 1
      
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
    REAL(dp), DIMENSION(C%nZ)                          :: prof
    REAL(dp)                                           :: A_flow_MISMIP
    
    REAL(dp), PARAMETER                                :: A_low_temp  = 1.14E-05_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: A_high_temp = 5.47E+10_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_low_temp  = 6.0E+04_dp    ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_high_temp = 13.9E+04_dp   ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    
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
    
  ! =======================================================
  ! ===== End of exceptions for benchmark experiments =====
  ! =======================================================
    
    ! Calculate the ice flow factor as a function of the ice temperature according to the Arrhenius relationship  (Huybrechts, 1992)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      DO k = 1, C%nz
        IF (ice%mask_ice_a( j,i) == 1) THEN
          IF (ice%Ti_a( k,j,i) < 263.15_dp) THEN
            ice%A_flow_3D_a( k,j,i) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * ice%Ti_a( k,j,i)))  
          ELSE
            ice%A_flow_3D_a( k,j,i) = A_high_temp * EXP(-Q_high_temp / (R_gas * ice%Ti_a( k,j,i)))  
          END IF
        ELSE
          IF (C%choice_ice_margin == 'BC') THEN
            ice%A_flow_3D_a( k,j,i) = 0._dp
          ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
            ! In the "infinite slab" case, calculate effective viscosity everywhere
            ! (even when there's technically no ice present)
            ice%A_flow_3D_a( k,j,i) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * 263.15_dp))  
          ELSE
            IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in calc_effective_viscosity!'
            CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
          END IF
        END IF
      END DO ! DO k = 1, C%nz
        
      ! Apply the flow enhancement factors
      IF (ice%mask_sheet_a( j,i) == 1) THEN
        ice%A_flow_3D_a( :,j,i) = ice%A_flow_3D_a( :,j,i) * C%m_enh_sheet
      ELSE
        ice%A_flow_3D_a( :,j,i) = ice%A_flow_3D_a( :,j,i) * C%m_enh_shelf
      END IF

      ! Vertical average
      prof = ice%A_flow_3D_a( :,j,i)
      ice%A_flow_vav_a( j,i) = vertical_average( prof)
       
    END DO
    END DO
    CALL sync
    
    ! Calculate the pressure melting point (Huybrechts, 1992), heat capacity (Pounder, 1965), and thermal conductivity (Ritz, 1987) of ice
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%Ti_pmp_a( :,j,i) = T0 - CC * ice%Hi_a( j,i) * C%zeta
      ice%Cpi_a(    :,j,i) = 2115.3_dp + 7.79293_dp * (ice%Ti_a( :,j,i) - T0)
      ice%Ki_a(     :,j,i) = 3.101E+08_dp * EXP(-0.0057_dp * ice%Ti_a( :,j,i))
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_3D( ice%A_flow_3D_a , 'ice%A_flow_3D_a' , 'ice_physical_properties')
    CALL check_for_NaN_dp_2D( ice%A_flow_vav_a, 'ice%A_flow_vav_a', 'ice_physical_properties')
    CALL check_for_NaN_dp_3D( ice%Ti_pmp_a    , 'ice%Ti_pmp_a'    , 'ice_physical_properties')
    CALL check_for_NaN_dp_3D( ice%Cpi_a       , 'ice%Cpi_a'       , 'ice_physical_properties')
    CALL check_for_NaN_dp_3D( ice%Ki_a        , 'ice%Ki_a'        , 'ice_physical_properties')
    
  END SUBROUTINE ice_physical_properties 
  
  SUBROUTINE determine_floating_margin_fraction( grid, ice)
    ! Determine the ice-filled fraction of floating margin pixels
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j,ii,jj
    LOGICAL                                            :: has_noncf_neighbours
    REAL(dp)                                           :: Hi_neighbour_max

    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1

      ice%float_margin_frac_a( j,i) = 1._dp
      ice%Hi_actual_cf_a(      j,i) = ice%Hi_a( j,i)
      
      IF (ice%mask_cf_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 1) THEN
        
        ! First check if any non-calving-front neighbours actually exist
        has_noncf_neighbours = .FALSE.
        DO ii = i-1,i+1
        DO jj = j-1,j+1
          IF (ice%mask_ice_a( jj,ii) == 1 .AND. ice%mask_cf_a( jj,ii) == 0) has_noncf_neighbours = .TRUE.
        END DO
        END DO
        
        ! If not, then the floating fraction is defined as 1
        IF (.NOT. has_noncf_neighbours) THEN
          ice%float_margin_frac_a( j,i) = 1._dp
          ice%Hi_actual_cf_a(      j,i) = ice%Hi_a( j,i)
          CYCLE
        END IF
        
        ! If so, find the ice thickness the thickest non-calving-front neighbour
        Hi_neighbour_max = 0._dp
        DO ii = i-1,i+1
        DO jj = j-1,j+1
          IF (ice%mask_ice_a( jj,ii) == 1 .AND. ice%mask_cf_a( jj,ii) == 0) THEN
            Hi_neighbour_max = MAX( Hi_neighbour_max, ice%Hi_a( jj,ii))
          END IF
        END DO
        END DO
        
        ! If the thickest non-calving-front neighbour has thinner ice, define the fraction as 1
        IF (Hi_neighbour_max < ice%Hi_a( j,i)) THEN
          ice%float_margin_frac_a( j,i) = 1._dp
          ice%Hi_actual_cf_a(      j,i) = ice%Hi_a( j,i)
          CYCLE
        END IF
        
        ! Calculate ice-filled fraction
        ice%float_margin_frac_a( j,i) = ice%Hi_a( j,i) / Hi_neighbour_max
        ice%Hi_actual_cf_a(      j,i) = Hi_neighbour_max
        
      END IF

    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%float_margin_frac_a, 'ice%float_margin_frac_a', 'determine_floating_margin_fraction')
    CALL check_for_NaN_dp_2D( ice%Hi_actual_cf_a     , 'ice%Hi_actual_cf_a'     , 'determine_floating_margin_fraction')
    
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
              
            END IF ! IF (neighbour exists) THEN
            
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
  
! == Routines for determining the grounded fraction on all four grids
  SUBROUTINE determine_grounded_fractions( grid, ice)
    ! Determine the grounded fraction of next-to-grounding-line pixels
    ! (used for determining basal friction in the DIVA)
    !
    ! Uses the bilinear interpolation scheme (with analytical solutions) from CISM (Leguy et al., 2021)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  f_grnd_a_NW,  f_grnd_a_NE,  f_grnd_a_SW,  f_grnd_a_SE
    INTEGER                                            :: wf_grnd_a_NW, wf_grnd_a_NE, wf_grnd_a_SW, wf_grnd_a_SE
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny  , grid%nx  , f_grnd_a_NW, wf_grnd_a_NW)
    CALL allocate_shared_dp_2D( grid%ny  , grid%nx  , f_grnd_a_NE, wf_grnd_a_NE)
    CALL allocate_shared_dp_2D( grid%ny  , grid%nx  , f_grnd_a_SW, wf_grnd_a_SW)
    CALL allocate_shared_dp_2D( grid%ny  , grid%nx  , f_grnd_a_SE, wf_grnd_a_SE)
    
    ! Calculate grounded fractions of all four quadrants of each a-grid cell
    CALL determine_grounded_fractions_CISM_quads( grid, ice, f_grnd_a_NW, f_grnd_a_NE, f_grnd_a_SW, f_grnd_a_SE)
    
    ! Get grounded fractions on all four grids by averaging over the quadrants
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
  
      ! a-grid
      ice%f_grnd_a( j,i) = 0.25_dp * (f_grnd_a_NW( j,i) + f_grnd_a_NE( j,i) + f_grnd_a_SW( j,i) + f_grnd_a_SE( j,i))
      
      ! cx-grid
      IF (i < grid%nx) THEN
        ice%f_grnd_cx( j,i) = 0.25_dp * (f_grnd_a_NE( j,i) + f_grnd_a_SE( j,i) + f_grnd_a_NW( j,i+1) + f_grnd_a_SW( j,i+1))
      END IF
      
      ! cy-grid
      IF (j < grid%ny) THEN
        ice%f_grnd_cy( j,i) = 0.25_dp * (f_grnd_a_NE( j,i) + f_grnd_a_NW( j,i) + f_grnd_a_SE( j+1,i) + f_grnd_a_SW( j+1,i))
      END IF
      
      ! b-grid
      IF (i < grid%nx .AND. j < grid%ny) THEN
        ice%f_grnd_b( j,i) = 0.25_dp * (f_grnd_a_NE( j,i) + f_grnd_a_NW( j,i+1) + f_grnd_a_SE( j+1,i) + f_grnd_a_SW( j+1,i+1))
      END IF
  
    END DO
    END DO
    
    ! Clean up after yourself
    CALL deallocate_shared( wf_grnd_a_NW)
    CALL deallocate_shared( wf_grnd_a_NE)
    CALL deallocate_shared( wf_grnd_a_SW)
    CALL deallocate_shared( wf_grnd_a_SE)
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%f_grnd_a , 'ice%f_grnd_a' , 'determine_grounded_fractions')
    CALL check_for_NaN_dp_2D( ice%f_grnd_cx, 'ice%f_grnd_cx', 'determine_grounded_fractions')
    CALL check_for_NaN_dp_2D( ice%f_grnd_cy, 'ice%f_grnd_cy', 'determine_grounded_fractions')
    CALL check_for_NaN_dp_2D( ice%f_grnd_b , 'ice%f_grnd_b' , 'determine_grounded_fractions')
    
  END SUBROUTINE determine_grounded_fractions
  SUBROUTINE determine_grounded_fractions_CISM_quads( grid, ice, f_grnd_a_NW, f_grnd_a_NE, f_grnd_a_SW, f_grnd_a_SE)
    ! Calculate grounded fractions of all four quadrants of each a-grid cell
    ! (using the approach from CISM, where grounded fractions are calculated
    !  based on analytical solutions to the bilinear interpolation)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: f_grnd_a_NW, f_grnd_a_NE, f_grnd_a_SW, f_grnd_a_SE
    
    ! Local variables:
    INTEGER                                            :: i,j,ii,jj
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  TAF_a_ext
    INTEGER                                            :: wTAF_a_ext
    REAL(dp)                                           :: f_NW, f_N, f_NE, f_W, f_m, f_E, f_SW, f_S, f_SE
    REAL(dp)                                           :: fq_NW, fq_NE, fq_SW, fq_SE
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny+2, grid%nx+2, TAF_a_ext, wTAF_a_ext)
    
    ! Fill in "extended" thickness above flotation (one extra row of pixels around the domain)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ii = i+1
      jj = j+1
      TAF_a_ext(jj,ii) = ice%TAF_a( j,i)
    END DO
    END DO
    IF (par%master) THEN
      TAF_a_ext( 2:grid%ny+1,1          ) = ice%TAF_a( :,1)
      TAF_a_ext( 2:grid%ny+1,  grid%nx+2) = ice%TAF_a( :,grid%nx)
      TAF_a_ext( 1          ,2:grid%nx+1) = ice%TAF_a( 1,:)
      TAF_a_ext(   grid%ny+2,2:grid%nx+1) = ice%TAF_a( grid%ny,:)
      TAF_a_ext( 1          ,1          ) = ice%TAF_a( 1,1)
      TAF_a_ext( 1          ,  grid%nx+2) = ice%TAF_a( 1,grid%nx)
      TAF_a_ext(   grid%ny+2,1          ) = ice%TAF_a( grid%ny,1)
      TAF_a_ext(   grid%ny+2,  grid%nx+2) = ice%TAF_a( grid%ny,grid%nx)
    END IF
    CALL sync
    
    ! Calculate grounded fractions of all four quadrants of each a-grid cell
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ii = i+1
      jj = j+1
      
      f_NW = 0.25_dp * (TAF_a_ext( jj+1,ii-1) + TAF_a_ext( jj+1,ii  ) + TAF_a_ext( jj  ,ii-1) + TAF_a_ext( jj  ,ii  ))
      f_N  = 0.5_dp  * (TAF_a_ext( jj+1,ii  ) + TAF_a_ext( jj  ,ii  ))
      f_NE = 0.25_dp * (TAF_a_ext( jj+1,ii  ) + TAF_a_ext( jj+1,ii+1) + TAF_a_ext( jj  ,ii  ) + TAF_a_ext( jj  ,ii+1))
      f_W  = 0.5_dp  * (TAF_a_ext( jj  ,ii-1) + TAF_a_ext( jj  ,ii  ))
      f_m  = TAF_a_ext( jj,ii)
      f_E  = 0.5_dp  * (TAF_a_ext( jj  ,ii  ) + TAF_a_ext( jj  ,ii+1))
      f_SW = 0.25_dp * (TAF_a_ext( jj  ,ii-1) + TAF_a_ext( jj  ,ii  ) + TAF_a_ext( jj-1,ii-1) + TAF_a_ext( jj-1,ii  ))
      f_S  = 0.5_dp  * (TAF_a_ext( jj  ,ii  ) + TAF_a_ext( jj-1,ii  ))
      f_SE = 0.25_dp * (TAF_a_ext( jj  ,ii  ) + TAF_a_ext( jj  ,ii+1) + TAF_a_ext( jj-1,ii  ) + TAF_a_ext( jj-1,ii+1))
      
      ! NW
      fq_NW = f_NW
      fq_NE = f_N
      fq_SW = f_W
      fq_SE = f_m
      CALL calc_fraction_above_zero( fq_NW, fq_NE,  fq_SW,  fq_SE,  f_grnd_a_NW( j,i))
      
      ! NE
      fq_NW = f_N
      fq_NE = f_NE
      fq_SW = f_m
      fq_SE = f_E
      CALL calc_fraction_above_zero( fq_NW, fq_NE,  fq_SW,  fq_SE,  f_grnd_a_NE( j,i))
      
      ! SW
      fq_NW = f_W
      fq_NE = f_m
      fq_SW = f_SW
      fq_SE = f_S
      CALL calc_fraction_above_zero( fq_NW, fq_NE,  fq_SW,  fq_SE,  f_grnd_a_SW( j,i))
      
      ! SE
      fq_NW = f_m
      fq_NE = f_E
      fq_SW = f_S
      fq_SE = f_SE
      CALL calc_fraction_above_zero( fq_NW, fq_NE,  fq_SW,  fq_SE,  f_grnd_a_SE( j,i))
      
    END DO
    END DO
    
    ! Clean up after yourself
    CALL deallocate_shared( wTAF_a_ext)
    
  END SUBROUTINE determine_grounded_fractions_CISM_quads
  SUBROUTINE calc_fraction_above_zero( f_NW, f_NE, f_SW, f_SE, phi)
    ! Given a square with function values at the four corners,
    ! calculate the fraction phi of the square where the function is larger than zero.
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: f_NW, f_NE, f_SW, f_SE
    REAL(dp),                            INTENT(OUT)   :: phi
    
    ! Local variables:
    REAL(dp)                                           :: f_NWp, f_NEp, f_SWp, f_SEp
    REAL(dp)                                           :: aa,bb,cc,dd,x,f1,f2
    INTEGER                                            :: scen
    
    ! The analytical solutions sometime give problems when one or more of the corner
    ! values is VERY close to zero; avoid this.
    IF (f_NW == 0._dp) THEN
      f_NWp = 1E-8_dp
    ELSEIF (f_NW > 0._dp) THEN
      f_NWp = MAX(  1E-8_dp, f_NW)
    ELSEIF (f_NW < 0._dp) THEN
      f_NWp = MIN( -1E-8_dp, f_NW)
    ELSE
      f_NWp = f_NW
    END IF
    IF (f_NE == 0._dp) THEN
      f_NEp = 1E-8_dp
    ELSEIF (f_NE > 0._dp) THEN
      f_NEp = MAX(  1E-8_dp, f_NE)
    ELSEIF (f_NE < 0._dp) THEN
      f_NEp = MIN( -1E-8_dp, f_NE)
    ELSE
      f_NEp = f_NE
    END IF
    IF (f_SW == 0._dp) THEN
      f_SWp = 1E-8_dp
    ELSEIF (f_SW > 0._dp) THEN
      f_SWp = MAX(  1E-8_dp, f_SW)
    ELSEIF (f_SW < 0._dp) THEN
      f_SWp = MIN( -1E-8_dp, f_SW)
    ELSE
      f_SWp = f_SW
    END IF
    IF (f_SE == 0._dp) THEN
      f_SEp = 1E-8_dp
    ELSEIF (f_SE > 0._dp) THEN
      f_SEp = MAX(  1E-8_dp, f_SE)
    ELSEIF (f_SE < 0._dp) THEN
      f_SEp = MIN( -1E-8_dp, f_SE)
    ELSE
      f_SEp = f_SE
    END IF
    
    IF     (f_NWp > 0._dp .AND. f_NEp > 0._dp .AND. f_SWp > 0._dp .AND. f_SEp > 0._dp) THEN
      ! All four corners are grounded.
      
      phi = 1._dp
      
    ELSEIF (f_NWp < 0._dp .AND. f_NEp < 0._dp .AND. f_SWp < 0._dp .AND. f_SEp < 0._dp) THEN
      ! All four corners are floating
      
      phi = 0._dp
      
    ELSE
      ! At least one corner is grounded and at least one is floating;
      ! the grounding line must pass through this square!
      
      ! Only four "scenarios" exist (with rotational symmetries):
      ! 1: SW grounded, rest floating
      ! 2: SW floating, rest grounded
      ! 3: south grounded, north floating
      ! 4: SW & NE grounded, SE & NW floating
      ! Rotate the four-corner world until it matches one of these scenarios.
      CALL rotate_quad_until_match( f_NWp, f_NEp, f_SWp, f_SEp, scen)
    
      IF (scen == 1) THEN
        ! 1: SW grounded, rest floating
        aa  = f_SWp
        bb  = f_SEp - f_SWp
        cc  = f_NWp - f_SWp
        dd  = f_NEp + f_SWp - f_NWp - f_SEp
        phi =         ((bb*cc - aa*dd) * LOG(ABS(1._dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
      ELSEIF (scen == 2) THEN
        ! 2: SW floating, rest grounded
        aa  = -(f_SWp)
        bb  = -(f_SEp - f_SWp)
        cc  = -(f_NWp - f_SWp)
        dd  = -(f_NEp + f_SWp - f_NWp - f_SEp)
        phi = 1._dp - ((bb*cc - aa*dd) * LOG(ABS(1._dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
      ELSEIF (scen == 3) THEN
        ! 3: south grounded, north floating
        aa  = f_SWp
        bb  = f_SEp - f_SWp
        cc  = f_NWp - f_SWp
        dd  = f_NEp + f_SWp - f_NWp - f_SEp
        x   = 0._dp
        f1  = ((bb*cc - aa*dd) * LOG(ABS(cc+dd*x)) - bb*dd*x) / (dd**2)
        x   = 1._dp
        f2  = ((bb*cc - aa*dd) * LOG(ABS(cc+dd*x)) - bb*dd*x) / (dd**2)
        phi = f2-f1
      ELSEIF (scen == 4) THEN
        ! 4: SW & NE grounded, SE & NW floating
        ! SW corner
        aa  = f_SWp
        bb  = f_SEp - f_SWp
        cc  = f_NWp - f_SWp
        dd  = f_NEp + f_SWp - f_NWp - f_SEp
        phi = ((bb*cc - aa*dd) * LOG(ABS(1._dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
        ! NE corner
        CALL rotate_quad( f_NWp, f_NEp, f_SWp, f_SEp)
        CALL rotate_quad( f_NWp, f_NEp, f_SWp, f_SEp)
        aa  = f_SWp
        bb  = f_SEp - f_SWp
        cc  = f_NWp - f_SWp
        dd  = f_NEp + f_SWp - f_NWp - f_SEp
        phi = phi + ((bb*cc - aa*dd) * LOG(ABS(1._dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
      ELSE
        WRITE(0,*) 'determine_grounded_fractions_CISM_quads - calc_fraction_above_zero - ERROR: unknown scenario [', scen, ']!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END IF
    
  END SUBROUTINE calc_fraction_above_zero
  SUBROUTINE rotate_quad_until_match( f_NW, f_NE, f_SW, f_SE, scen)
    ! Rotate the four corners until one of the four possible scenarios is found.
    ! 1: SW grounded, rest floating
    ! 2: SW floating, rest grounded
    ! 3: south grounded, north floating
    ! 4: SW & NE grounded, SE & NW floating
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(INOUT) :: f_NW, f_NE, f_SW, f_SE
    INTEGER,                             INTENT(OUT)   :: scen
    
    ! Local variables:
    LOGICAL                                            :: found_match
    INTEGER                                            :: nit
    
    found_match = .FALSE.
    scen        = 0
    nit         = 0
    
    DO WHILE (.NOT. found_match)
      
      nit = nit+1
      
      CALL rotate_quad( f_NW, f_NE, f_SW, f_SE)
      
      IF     (f_SW > 0._dp .AND. f_SE < 0._dp .AND. f_NE < 0._dp .AND. f_NW < 0._dp) THEN
        ! 1: SW grounded, rest floating
        scen = 1
        found_match = .TRUE.
      ELSEIF (f_SW < 0._dp .AND. f_SE > 0._dp .AND. f_NE > 0._dp .AND. f_NW > 0._dp) THEN
        ! 2: SW floating, rest grounded
        scen = 2
        found_match = .TRUE.
      ELSEIF (f_SW > 0._dp .AND. f_SE > 0._dp .AND. f_NE < 0._dp .AND. f_NW < 0._dp) THEN
        ! 3: south grounded, north floating
        scen = 3
        found_match = .TRUE.
      ELSEIF (f_SW > 0._dp .AND. f_SE < 0._dp .AND. f_NE > 0._dp .AND. f_NW < 0._dp) THEN
        ! 4: SW & NE grounded, SE & NW floating
        scen = 4
        found_match = .TRUE.
      END IF
      
      IF (nit > 4) THEN
        IF (par%master) WRITE(0,*) 'determine_grounded_fractions_CISM_quads - rotate_quad_until_match - ERROR: couldnt find matching scenario!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END DO
    
  END SUBROUTINE rotate_quad_until_match
  SUBROUTINE rotate_quad( f_NW, f_NE, f_SW, f_SE)
    ! Rotate the four corners anticlockwise by 90 degrees
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(INOUT) :: f_NW, f_NE, f_SW, f_SE
    
    ! Local variables:
    REAL(dp), DIMENSION(4)                             :: fvals
    
    fvals = [f_NW,f_NE,f_SE,f_SW]
    f_NW = fvals( 2)
    f_NE = fvals( 3)
    f_SE = fvals( 4)
    f_SW = fvals( 1)
    
  END SUBROUTINE rotate_quad

END MODULE general_ice_model_data_module
