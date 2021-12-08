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
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_model_region
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE parameters_module
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             is_floating, surface_elevation, thickness_above_floatation, &
                                             is_in_polygon, oblique_sg_projection
  USE derivatives_and_grids_module,    ONLY: map_a_to_cx_2D, map_a_to_cy_2D, ddx_a_to_cx_2D, ddy_a_to_cy_2D, &
                                             ddy_a_to_cx_2D, ddx_a_to_cy_2D, map_a_to_cx_3D, map_a_to_cy_3D, &
                                             ddx_a_to_a_2D, ddy_a_to_a_2D, map_a_to_b_2D

  IMPLICIT NONE

CONTAINS
  
! == Update general ice model data - Hs, masks, ice physical properties (the main routine that's called from the main ice model)
  SUBROUTINE update_general_ice_model_data( grid, ice)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables
    INTEGER                                            :: i,j
    
    ! Calculate surface elevation and thickness above floatation
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%Hs_a(  j,i) = surface_elevation( ice%Hi_a( j,i), ice%Hb_a( j,i), ice%SL_a( j,i))
      ice%TAF_a( j,i) = thickness_above_floatation( ice%Hi_a( j,i), ice%Hb_a( j,i), ice%SL_a( j,i))
    END DO
    END DO
    CALL sync
    
    ! Determine masks
    CALL determine_masks( grid, ice)
 
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
    
  END SUBROUTINE update_general_ice_model_data
  
! == Routines for determining the model masks
  SUBROUTINE determine_masks( grid, ice)
    ! Determine the different masks
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(INOUT) :: ice 
    
    CALL determine_masks_landocean(   grid, ice)
    CALL determine_masks_ice(         grid, ice)
    CALL determine_masks_transitions( grid, ice)
  
  END SUBROUTINE determine_masks
  SUBROUTINE determine_masks_landocean( grid, ice)
    ! Determine the different masks
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(INOUT) :: ice 
  
    INTEGER                                            :: i,j
    
    ! Initialise: start with land everywhere.
    ice%mask_land_a(   :,grid%i1:grid%i2) = 1
    ice%mask_ocean_a(  :,grid%i1:grid%i2) = 0
    ice%mask_a(        :,grid%i1:grid%i2) = C%type_land
    CALL sync
    
    ! Determine land/ocean masks
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
  
  END SUBROUTINE determine_masks_landocean
  SUBROUTINE determine_masks_ice( grid, ice)
    ! Determine the different masks
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(INOUT) :: ice 
  
    INTEGER                                            :: i,j
    
    CALL determine_masks_landocean( grid, ice)
    
    ! Initialise: start with land everywhere.
    ice%mask_ice_a(    :,grid%i1:grid%i2) = 0
    ice%mask_sheet_a(  :,grid%i1:grid%i2) = 0
    ice%mask_shelf_a(  :,grid%i1:grid%i2) = 0
    CALL sync
    
    ! Determine ice/sheet/shelf masks
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Ice
      IF (ice%Hi_a( j,i) > TINY(ice%Hi_a( 1,1))) THEN
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
  
  END SUBROUTINE determine_masks_ice
  SUBROUTINE determine_masks_transitions( grid, ice)
    ! Determine the different masks
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(INOUT) :: ice 
  
    INTEGER                                            :: i,j
    
    CALL determine_masks_landocean( grid, ice)
    
    ! Initialise
    ice%mask_coast_a(  :,grid%i1:grid%i2) = 0
    ice%mask_margin_a( :,grid%i1:grid%i2) = 0
    ice%mask_gl_a(     :,grid%i1:grid%i2) = 0
    ice%mask_cf_a(     :,grid%i1:grid%i2) = 0
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
  
  END SUBROUTINE determine_masks_transitions
  
! == Routines for calculating sub-grid grounded fractions
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
    REAL(dp), PARAMETER                                :: ftol = 1E-4_dp
    
    ! The analytical solutions sometime give problems when one or more of the corner
    ! values is VERY close to zero; avoid this.
    IF (f_NW == 0._dp) THEN
      f_NWp = ftol * 1.1_dp
    ELSEIF (f_NW > 0._dp) THEN
      f_NWp = MAX(  ftol * 1.1_dp, f_NW)
    ELSEIF (f_NW < 0._dp) THEN
      f_NWp = MIN( -ftol * 1.1_dp, f_NW)
    ELSE
      f_NWp = f_NW
    END IF
    IF (f_NE == 0._dp) THEN
      f_NEp = ftol * 1.21_dp
    ELSEIF (f_NE > 0._dp) THEN
      f_NEp = MAX(  ftol * 1.21_dp, f_NE)
    ELSEIF (f_NE < 0._dp) THEN
      f_NEp = MIN( -ftol * 1.21_dp, f_NE)
    ELSE
      f_NEp = f_NE
    END IF
    IF (f_SW == 0._dp) THEN
      f_SWp = ftol * 1.13_dp
    ELSEIF (f_SW > 0._dp) THEN
      f_SWp = MAX(  ftol * 1.13_dp, f_SW)
    ELSEIF (f_SW < 0._dp) THEN
      f_SWp = MIN( -ftol * 1.13_dp, f_SW)
    ELSE
      f_SWp = f_SW
    END IF
    IF (f_SE == 0._dp) THEN
      f_SEp = ftol * 0.97_dp
    ELSEIF (f_SE > 0._dp) THEN
      f_SEp = MAX(  ftol * 0.97_dp, f_SE)
    ELSEIF (f_SE < 0._dp) THEN
      f_SEp = MIN( -ftol * 0.97_dp, f_SE)
    ELSE
      f_SEp = f_SE
    END IF
    
    IF     (f_NWp >= 0._dp .AND. f_NEp >= 0._dp .AND. f_SWp >= 0._dp .AND. f_SEp >= 0._dp) THEN
      ! All four corners are grounded.
      
      phi = 1._dp
      
    ELSEIF (f_NWp <= 0._dp .AND. f_NEp <= 0._dp .AND. f_SWp <= 0._dp .AND. f_SEp <= 0._dp) THEN
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
          
        ! Exception for when d=0
        IF (ABS(dd) < 1E-4_dp) THEN
          IF (f_SWp > 0._dp) THEN
            f_SWp = f_SWp + 0.1_dp
          ELSE
            f_SWp = f_SWp - 0.1_dp
          END IF
          aa  = f_SWp
          bb  = f_SEp - f_SWp
          cc  = f_NWp - f_SWp
          dd  = f_NEp + f_SWp - f_NWp - f_SEp
          phi =         ((bb*cc - aa*dd) * LOG(ABS(1._dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
        END IF
        
      ELSEIF (scen == 2) THEN
        ! 2: SW floating, rest grounded
        
        aa  = -(f_SWp)
        bb  = -(f_SEp - f_SWp)
        cc  = -(f_NWp - f_SWp)
        dd  = -(f_NEp + f_SWp - f_NWp - f_SEp)
        phi = 1._dp - ((bb*cc - aa*dd) * LOG(ABS(1._dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
          
        ! Exception for when d=0
        IF (ABS(dd) < 1E-4_dp) THEN
          IF (f_SWp > 0._dp) THEN
            f_SWp = f_SWp + 0.1_dp
          ELSE
            f_SWp = f_SWp - 0.1_dp
          END IF
          aa  = -(f_SWp)
          bb  = -(f_SEp - f_SWp)
          cc  = -(f_NWp - f_SWp)
          dd  = -(f_NEp + f_SWp - f_NWp - f_SEp)
          phi = 1._dp - ((bb*cc - aa*dd) * LOG(ABS(1._dp - (aa*dd)/(bb*cc))) + aa*dd) / (dd**2)
        END IF
        
      ELSEIF (scen == 3) THEN
        ! 3: south grounded, north floating
        
        ! Exception for when the GL runs parallel to the x-axis
        IF (ABS( 1._dp - f_NWp/f_NEp) < 1E-6_dp .AND. ABS( 1._dp - f_SWp/f_SEp) < 1E-6_dp) THEN
          
          phi = f_SWp / (f_SWp - f_NWp)
          
        ELSE
        
          aa  = f_SWp
          bb  = f_SEp - f_SWp
          cc  = f_NWp - f_SWp
          dd  = f_NEp + f_SWp - f_NWp - f_SEp
          x   = 0._dp
          f1  = ((bb*cc - aa*dd) * LOG(ABS(cc+dd*x)) - bb*dd*x) / (dd**2)
          x   = 1._dp
          f2  = ((bb*cc - aa*dd) * LOG(ABS(cc+dd*x)) - bb*dd*x) / (dd**2)
          phi = f2-f1
          
          ! Exception for when d=0
          IF (ABS(dd) < 1E-4_dp) THEN
            IF (f_SWp > 0._dp) THEN
              f_SWp = f_SWp + 0.1_dp
            ELSE
              f_SWp = f_SWp - 0.1_dp
            END IF
            aa  = f_SWp
            bb  = f_SEp - f_SWp
            cc  = f_NWp - f_SWp
            dd  = f_NEp + f_SWp - f_NWp - f_SEp
            x   = 0._dp
            f1  = ((bb*cc - aa*dd) * LOG(ABS(cc+dd*x)) - bb*dd*x) / (dd**2)
            x   = 1._dp
            f2  = ((bb*cc - aa*dd) * LOG(ABS(cc+dd*x)) - bb*dd*x) / (dd**2)
            phi = f2-f1
          END IF
        
        END IF
        
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
    
    IF (pi < -0.01_dp .OR. phi > 1.01_dp .OR. phi /= phi) THEN
      WRITE(0,*) 'calc_fraction_above_zero - ERROR: phi = ', phi
      WRITE(0,*) 'f = [', f_NWp, ',', f_NEp, ',', f_SWp, ',', f_SEp, ']'
      WRITE(0,*) 'scen = ', scen
      WRITE(0,*) 'aa = ', aa, ', bb = ', bb, ', cc = ', cc, ', dd = ', dd, ', f1 = ', f1, ',f2 = ', f2
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    phi = MAX( 0._dp, MIN( 1._dp, phi))
    
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
  
! == Routines for calculating the sub-grid ice-filled fraction and effective ice thickness at the calving front
  SUBROUTINE determine_floating_margin_fraction( grid, ice)
    ! Determine the ice-filled fraction and effective ice thickness of floating margin pixels
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j,ii,jj
    LOGICAL                                            :: has_noncf_neighbours
    REAL(dp)                                           :: Hi_neighbour_max

    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1

      ! Initialise
      IF (ice%mask_ice_a( j,i) == 1) THEN
        ice%float_margin_frac_a( j,i) = 1._dp
        ice%Hi_eff_cf_a(         j,i) = ice%Hi_a( j,i)
      ELSE
        ice%float_margin_frac_a( j,i) = 0._dp
        ice%Hi_eff_cf_a(         j,i) = 0._dp
      END IF
      
      IF (ice%mask_cf_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 1) THEN
        
        ! First check if any non-calving-front neighbours actually exist
        has_noncf_neighbours = .FALSE.
        DO ii = MAX( 1, i-1), MIN( grid%nx, i+1)
        DO jj = MAX( 1, j-1), MIN( grid%ny, j+1)
          IF (ice%mask_ice_a( jj,ii) == 1 .AND. ice%mask_cf_a( jj,ii) == 0) has_noncf_neighbours = .TRUE.
        END DO
        END DO
        
        ! If not, then the floating fraction is defined as 1
        IF (.NOT. has_noncf_neighbours) THEN
          ice%float_margin_frac_a( j,i) = 1._dp
          ice%Hi_eff_cf_a(         j,i) = ice%Hi_a( j,i)
          CYCLE
        END IF
        
        ! If so, find the ice thickness the thickest non-calving-front neighbour
        Hi_neighbour_max = 0._dp
        DO ii = MAX( 1, i-1), MIN( grid%nx, i+1)
        DO jj = MAX( 1, j-1), MIN( grid%ny, j+1)
          IF (ice%mask_ice_a( jj,ii) == 1 .AND. ice%mask_cf_a( jj,ii) == 0) THEN
            Hi_neighbour_max = MAX( Hi_neighbour_max, ice%Hi_a( jj,ii))
          END IF
        END DO
        END DO
        
        ! If the thickest non-calving-front neighbour has thinner ice, define the fraction as 1
        IF (Hi_neighbour_max < ice%Hi_a( j,i)) THEN
          ice%float_margin_frac_a( j,i) = 1._dp
          ice%Hi_eff_cf_a(         j,i) = ice%Hi_a( j,i)
          CYCLE
        END IF
        
        ! Calculate ice-filled fraction
        ice%float_margin_frac_a( j,i) = ice%Hi_a( j,i) / Hi_neighbour_max
        ice%Hi_eff_cf_a(         j,i) = Hi_neighbour_max
        
      END IF

    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%float_margin_frac_a, 'ice%float_margin_frac_a', 'determine_floating_margin_fraction')
    CALL check_for_NaN_dp_2D( ice%Hi_eff_cf_a        , 'ice%Hi_eff_cf_a'        , 'determine_floating_margin_fraction')
    
  END SUBROUTINE determine_floating_margin_fraction
  
! == A simple flood-fill algorithm for determining the ocean mask without creating inland seas / lakes
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
  
! == Routines for defining ice drainage basins from an external polygon file
  SUBROUTINE initialise_basins( grid, basin_ID, nbasins, region_name)
    ! Define the ice basins mask from an external text file
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    INTEGER,  DIMENSION(:,:  ),          INTENT(INOUT) :: basin_ID
    INTEGER,                             INTENT(INOUT) :: nbasins
    CHARACTER(LEN=3)                                   :: region_name
    
    ! Local variables:
    INTEGER                                            :: i,j,vi,vj,bi
    CHARACTER(LEN=256)                                 :: choice_basin_scheme
    CHARACTER(LEN=256)                                 :: filename_basins
    LOGICAL                                            :: recognised_file
    INTEGER                                            :: n_skip
    INTEGER                                            :: n_header_lines, n_vertices
    REAL(dp), DIMENSION(:    ), POINTER                :: Vlat, Vlon, Vx, Vy
    INTEGER,  DIMENSION(:    ), POINTER                :: Vid
    INTEGER                                            :: wVlat, wVlon, wVx, wVy, wVid
    REAL(dp)                                           :: VID_dp
    INTEGER                                            :: ios, dummy
    INTEGER                                            :: vi1, vi2
    REAL(dp), DIMENSION(:,:  ), POINTER                :: poly_bi
    INTEGER                                            :: wpoly_bi
    REAL(dp)                                           :: xmin,xmax,ymin,ymax
    REAL(dp), DIMENSION(2)                             :: p
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  basin_ID_loc,  basin_ID_ext_loc
    INTEGER                                            :: wbasin_ID_loc, wbasin_ID_ext_loc
    INTEGER                                            :: n_front
    INTEGER,  DIMENSION(:,:  ), POINTER                :: ij_front
    INTEGER                                            :: wij_front
    INTEGER                                            :: ii,jj,k,i_nearest,j_nearest
    REAL(dp)                                           :: dist,dist_min
    
    ! Determine what to do for this region
    choice_basin_scheme = 'none'
    filename_basins     = ''
    IF (region_name == 'NAM') THEN
      choice_basin_scheme = C%choice_basin_scheme_NAM
      filename_basins     = C%filename_basins_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_basin_scheme = C%choice_basin_scheme_EAS
      filename_basins     = C%filename_basins_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_basin_scheme = C%choice_basin_scheme_GRL
      filename_basins     = C%filename_basins_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_basin_scheme = C%choice_basin_scheme_ANT
      filename_basins     = C%filename_basins_ANT
    END IF
    
    IF     (choice_basin_scheme == 'none') THEN
      ! No basins are defined (i.e. the whole region is one big, big basin)
      
      basin_ID( :,grid%i1:grid%i2) = 1
      IF (par%master) nbasins = 1
      CALL sync
      
    ELSEIF (choice_basin_scheme == 'file') THEN
      ! Define basins from an external text file describing the polygons
      
      IF (par%master) WRITE(0,*) '  Reading basins for ', TRIM(region_name), ' from file "', TRIM(filename_basins), '"'
      
      IF (region_name == 'ANT') THEN
        ! Antarctica: ant_full_drainagesystem_polygons.txt
        ! Can be downloaded from: https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems
        ! A text file with 7 header lines, followed by three columns of data:
        !    Lat, Lon, basin ID
        
      ! ===== Check if this is really the file we're reading =====
      ! ==========================================================
        
        recognised_file = .FALSE.
        n_header_lines  = 0
        n_vertices      = 0
        n_skip          = 0
        
        DO i = 1, 256-36
          IF (filename_basins(i:i+35) == 'ant_full_drainagesystem_polygons.txt') THEN
            recognised_file = .TRUE.
            n_header_lines  = 7
            n_skip          = 5 ! Since the polygon file is at a ridiculously high resolution, downscale it a bit for efficiency
            n_vertices      = CEILING( REAL(901322,dp) / REAL(n_skip,dp))
          END IF
        END DO
        
        IF ((.NOT. recognised_file) .AND. par%master) THEN
          WRITE(0,*) ''
          WRITE(0,*) ' ===== '
          WRITE(0,*) 'initialise_basins - WARNING: for Antarctica, we expect the file "ant_full_drainagesystem_polygons.txt"'
          WRITE(0,*) '                             This can be downloaded from https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems'
          WRITE(0,*) '                             If another file is used, make sure that it has 7 header lines, or alternatively just change the code in initialise_basins'
          WRITE(0,*) ' ===== '
          WRITE(0,*) ''
        END IF
        
        ! Allocate shared memory
        CALL allocate_shared_dp_1D(  n_vertices, Vlat, wVlat)
        CALL allocate_shared_dp_1D(  n_vertices, Vlon, wVlon)
        CALL allocate_shared_dp_1D(  n_vertices, Vx,   wVx  )
        CALL allocate_shared_dp_1D(  n_vertices, Vy,   wVy  )
        CALL allocate_shared_int_1D( n_vertices, VID,  wVID )
        
      ! ===== Read the file =====
      ! =========================
        
        IF (par%master) THEN
          
          ! Open the file
          OPEN( UNIT = 1337, FILE=filename_basins, ACTION='READ')
          
          ! Skip the header lines
          DO i = 1, n_header_lines
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
          END DO
          
          ! Read the actual data
          DO vi = 1, n_vertices
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) Vlat( vi), Vlon( vi), VID( vi)
            IF (ios /= 0) THEN
              WRITE(0,*) ' initialise_basins - ERROR: length of text file "', TRIM( filename_basins), '" does not match n_vertices = ', n_vertices
              CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
            END IF
            
            DO vj = 1, n_skip-1
              READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
            END DO
          END DO
          
          ! Close the file
          CLOSE( UNIT = 1337)
          
        END IF ! IF (par%master) THEN
        CALL sync
        
        ! Project [lat,lon] to [x,y]
        CALL partition_list( n_vertices, par%i, par%n, vi1, vi2)
        DO vi = vi1, vi2
          CALL oblique_sg_projection( Vlon( vi), Vlat( vi), grid%lambda_M, grid%phi_M, grid%alpha_stereo, Vx( vi), Vy( vi))
        END DO
        
      ELSEIF (region_name == 'GRL') THEN
        ! Greenland: grndrainagesystems_ekholm.txt
        ! Can be downloaded from: https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems
        ! A text file with 7 header lines, followed by three columns of data:
        !    basin ID, Lat, Lon   (NOTE: here basin ID is a floating-point rather than an integer!)
        
      ! ===== Check if this is really the file we're reading =====
      ! ==========================================================
        
        recognised_file = .FALSE.
        n_header_lines  = 0
        n_vertices      = 0
        n_skip          = 0
        
        DO i = 1, 256-29
          IF (filename_basins(i:i+28) == 'grndrainagesystems_ekholm.txt') THEN
            recognised_file = .TRUE.
            n_header_lines  = 7
            n_skip          = 5 ! Since the polygon file is at a ridiculously high resolution, downscale it a bit for efficiency
            n_vertices      = CEILING( REAL(272965,dp) / REAL(n_skip,dp))
          END IF
        END DO
        
        IF ((.NOT. recognised_file) .AND. par%master) THEN
          WRITE(0,*) ''
          WRITE(0,*) ' ===== '
          WRITE(0,*) 'initialise_basins - WARNING: for Greenland, we expect the file "grndrainagesystems_ekholm.txt"'
          WRITE(0,*) '                             This can be downloaded from https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems'
          WRITE(0,*) '                             If another file is used, make sure that it has 7 header lines, or alternatively just change the code in initialise_basins'
          WRITE(0,*) ' ===== '
          WRITE(0,*) ''
        END IF
        
        ! Allocate shared memory
        CALL allocate_shared_dp_1D(  n_vertices, Vlat,   wVlat  )
        CALL allocate_shared_dp_1D(  n_vertices, Vlon,   wVlon  )
        CALL allocate_shared_dp_1D(  n_vertices, Vx,     wVx    )
        CALL allocate_shared_dp_1D(  n_vertices, Vy,     wVy    )
        CALL allocate_shared_int_1D( n_vertices, VID,    wVID   )
        
      ! ===== Read the file =====
      ! =========================
        
        IF (par%master) THEN
          
          ! Open the file
          OPEN( UNIT = 1337, FILE=filename_basins, ACTION='READ')
          
          ! Skip the header lines
          DO i = 1, n_header_lines
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
          END DO
          
          ! Read the actual data
          DO vi = 1, n_vertices
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) VID_dp, Vlat( vi), Vlon( vi)
            IF (ios /= 0) THEN
              WRITE(0,*) ' initialise_basins - ERROR: length of text file "', TRIM( filename_basins), '" does not match n_vertices = ', n_vertices
              CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
            END IF
            
            DO vj = 1, n_skip-1
              READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy
            END DO
            
            ! Convert basin ID from floating point to integer
            IF     ( viD_dp == 1.1_dp) THEN
              VID( vi) = 1
            ELSEIF ( viD_dp == 1.2_dp) THEN
              VID( vi) = 2
            ELSEIF ( viD_dp == 1.3_dp) THEN
              VID( vi) = 3
            ELSEIF ( viD_dp == 1.4_dp) THEN
              VID( vi) = 4
            ELSEIF ( viD_dp == 2.1_dp) THEN
              VID( vi) = 5
            ELSEIF ( viD_dp == 2.2_dp) THEN
              VID( vi) = 6
            ELSEIF ( viD_dp == 3.1_dp) THEN
              VID( vi) = 7
            ELSEIF ( viD_dp == 3.2_dp) THEN
              VID( vi) = 8
            ELSEIF ( viD_dp == 3.3_dp) THEN
              VID( vi) = 9
            ELSEIF ( viD_dp == 4.1_dp) THEN
              VID( vi) = 10
            ELSEIF ( viD_dp == 4.2_dp) THEN
              VID( vi) = 11
            ELSEIF ( viD_dp == 4.3_dp) THEN
              VID( vi) = 12
            ELSEIF ( viD_dp == 5.0_dp) THEN
              VID( vi) = 13
            ELSEIF ( viD_dp == 6.1_dp) THEN
              VID( vi) = 14
            ELSEIF ( viD_dp == 6.2_dp) THEN
              VID( vi) = 15
            ELSEIF ( viD_dp == 7.1_dp) THEN
              VID( vi) = 16
            ELSEIF ( viD_dp == 7.2_dp) THEN
              VID( vi) = 17
            ELSEIF ( viD_dp == 8.1_dp) THEN
              VID( vi) = 18
            ELSEIF ( viD_dp == 8.2_dp) THEN
              VID( vi) = 19
            ELSE
              WRITE(0,*) 'initialise_basins - ERROR: unrecognised floating-point basin ID in file "', TRIM( filename_basins), '"!'
              CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
            END IF
            
          END DO
          
          ! Close the file
          CLOSE( UNIT = 1337)
          
        END IF ! IF (par%master) THEN
        CALL sync
        
        ! Project [lat,lon] to [x,y]
        CALL partition_list( n_vertices, par%i, par%n, vi1, vi2)
        DO vi = vi1, vi2
          CALL oblique_sg_projection( Vlon( vi), Vlat( vi), grid%lambda_M, grid%phi_M, grid%alpha_stereo, Vx( vi), Vy( vi))
        END DO
        
      END IF ! IF (region_name == 'ANT') THEN
      
    ! ===== Fill in the basins =====
    ! ==============================
    
      ! Allocate shared memory
      CALL allocate_shared_int_2D( grid%ny, grid%nx, basin_ID_loc, wbasin_ID_loc)
      
      ! Determine number of basins
      IF (par%master) nbasins = MAXVAL( VID)
      CALL sync
      
      DO bi = 1, nbasins
        
        ! Find range of vertices [vi1,vi2] for this basin
        vi1 = 1
        DO WHILE (VID( vi1) /= bi)
          vi1 = vi1 + 1
        END DO
        vi2 = vi1
        DO WHILE (VID( vi2) == bi .AND. vi2 < n_vertices)
          vi2 = vi2 + 1
        END DO
        vi2 = vi2 - 1
        
        ! Copy these vertices to a single array
        CALL allocate_shared_dp_2D( vi2+1-vi1, 2, poly_bi, wpoly_bi)
        IF (par%master) THEN
          poly_bi( :,1) = Vx( vi1:vi2)
          poly_bi( :,2) = Vy( vi1:vi2)
        END IF
        CALL sync
  
        ! Determine maximum polygon extent, for quick checks
        xmin = MINVAL( poly_bi(:,1))
        xmax = MAXVAL( poly_bi(:,1))
        ymin = MINVAL( poly_bi(:,2))
        ymax = MAXVAL( poly_bi(:,2))
        
        ! Check which grid cells lie inside the polygon spanned by these vertices
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          p = [grid%x( i), grid%y( j)]
  
          ! Quick test
          IF (p(1) < xmin .OR. p(1) > xmax .OR. &
              p(2) < ymin .OR. p(2) > ymax) THEN
            ! p cannot lie in the polygon, don't bother checking
          ElSE
            IF (is_in_polygon( poly_bi, p)) basin_ID_loc( j,i) = bi
          END IF
          
        END DO
        END DO
        CALL sync
        
        ! Clean up this basin's polygon
        CALL deallocate_shared( wpoly_bi)
        
      END DO ! DO bi = 1, nbasins
      
    ! ===== Extend basins into the ocean =====
    ! ========================================
    
      ! Allocate shared memory
      CALL allocate_shared_int_2D( grid%ny, grid%nx, basin_ID_ext_loc, wbasin_ID_ext_loc)
      
      ! Copy data
      basin_ID_ext_loc( :,grid%i1:grid%i2) = basin_ID_loc( :,grid%i1:grid%i2)
      CALL sync
      
      ! Compile list of ice-front pixels and their basin ID's
      IF (par%master) THEN
        n_front = 0
        DO i = 2, grid%nx-1
        DO j = 2, grid%ny-1
          IF (basin_ID_loc( j,i) > 0) THEN
            IF (MINVAL( basin_ID_loc( j-1:j+1,i-1:i+1)) == 0) THEN
              n_front = n_front + 1
            END IF
          END IF
        END DO
        END DO
      END IF
      CALL MPI_BCAST( n_front, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      
      CALL allocate_shared_int_2D( n_front, 2, ij_front, wij_front)
      
      IF (par%master) THEN
        k = 0
        DO i = 2, grid%nx-1
        DO j = 2, grid%ny-1
          IF (basin_ID_loc( j,i) > 0) THEN
            IF (MINVAL( basin_ID_loc( j-1:j+1,i-1:i+1)) == 0) THEN
              k = k + 1
              ij_front( k,:) = [i,j]
            END IF
          END IF
        END DO
        END DO
      END IF
      CALL sync
      
      ! For all non-assigned grid cells, find the nearest front cell and copy that basin ID
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        
        IF (basin_ID_ext_loc( j,i) == 0) THEN
          
          ! Go over all front cells, find the nearest
          dist_min  = REAL(MAX(grid%nx,grid%ny),dp) * grid%dx
          i_nearest = 0
          j_nearest = 0
          DO k = 1, n_front
            ii = ij_front( k,1)
            jj = ij_front( k,2)
            dist = SQRT( REAL(ii-i,dp)**2 + REAL(jj-j,dp)**2) * grid%dx
            IF (dist < dist_min) THEN
              dist_min = dist
              i_nearest = ii
              j_nearest = jj
            END IF
          END DO ! DO k = 1, n_front
          
          ! Assign basin ID of nearest front cell
          basin_ID_ext_loc( j,i) = basin_ID_loc( j_nearest, i_nearest)
          
        END IF
        
      END DO
      END DO
      CALL sync
      
    ! ===== Copy final result to the ice structure =====
    ! ==================================================
      
      basin_ID( :,grid%i1:grid%i2) = basin_ID_ext_loc( :,grid%i1:grid%i2)
      CALL sync
      
      ! Clean up after yourself
      CALL deallocate_shared( wVlat            )
      CALL deallocate_shared( wVlon            )
      CALL deallocate_shared( wVx              )
      CALL deallocate_shared( wVy              )
      CALL deallocate_shared( wVID             )
      CALL deallocate_shared( wbasin_ID_loc    )
      CALL deallocate_shared( wbasin_ID_ext_loc)
      CALL deallocate_shared( wij_front        )
      
    ! ==== If so specified, merge certain ice basins =====
    ! ====================================================
      
      IF     (region_name == 'ANT' .AND. C%do_merge_basins_ANT) THEN
        CALL merge_basins_ANT( grid, basin_ID, nbasins)
      ELSEIF (region_name == 'GRL' .AND. C%do_merge_basins_GRL) THEN
        CALL merge_basins_GRL( grid, basin_ID, nbasins)
      END IF
      
    ELSE
      IF (par%master) WRITE(0,*) 'initialise_basins - ERROR: unknown choice_basin_scheme "', TRIM(choice_basin_scheme), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE initialise_basins
  SUBROUTINE merge_basins_ANT(  grid, basin_ID, nbasins)
    ! Merge some of the Antarctic basins to match the more common definitions
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    INTEGER,  DIMENSION(:,:  ),          INTENT(INOUT) :: basin_ID
    INTEGER,                             INTENT(INOUT) :: nbasins
    
    ! Local variables
    INTEGER                                            :: i,j,k
    INTEGER,  DIMENSION( 27)                           :: T
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  basin_ID_old
    INTEGER                                            :: wbasin_ID_old
    
    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx, basin_ID_old, wbasin_ID_old)
    
    ! Copy data
    basin_ID_old( :,grid%i1:grid%i2) = basin_ID( :,grid%i1:grid%i2)
    CALL sync
    
    ! Define the merging table
    T(  1) = 1
    T(  2) = 1
    T(  3) = 1
    T(  4) = 2
    T(  5) = 3
    T(  6) = 4
    T(  7) = 5
    T(  8) = 6
    T(  9) = 6
    T( 10) = 6
    T( 11) = 6
    T( 12) = 7
    T( 13) = 8
    T( 14) = 9
    T( 15) = 10
    T( 16) = 11
    T( 17) = 11
    T( 18) = 11
    T( 19) = 11
    T( 20) = 12
    T( 21) = 13
    T( 22) = 13
    T( 23) = 13
    T( 24) = 14
    T( 25) = 15
    T( 26) = 16
    T( 27) = 17
    
    ! Perform the merge
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      DO k = 1, 27
        IF (basin_ID_old( j,i) == k) basin_ID( j,i) = T( k)
      END DO
    END DO
    END DO
    IF (par%master) nbasins = 17
    CALL sync
    
    
    ! Clean up after yourself
    CALL deallocate_shared( wbasin_ID_old)
    
  END SUBROUTINE merge_basins_ANT
  SUBROUTINE merge_basins_GRL(  grid, basin_ID, nbasins)
    ! Merge some of the Greenland basins to match the more common definitions
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    INTEGER,  DIMENSION(:,:  ),          INTENT(INOUT) :: basin_ID
    INTEGER,                             INTENT(INOUT) :: nbasins
    
    ! Local variables
    INTEGER                                            :: i,j,k
    INTEGER,  DIMENSION( 19)                           :: T
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  basin_ID_old
    INTEGER                                            :: wbasin_ID_old
    
    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx, basin_ID_old, wbasin_ID_old)
    
    ! Copy data
    basin_ID_old( :,grid%i1:grid%i2) = basin_ID( :,grid%i1:grid%i2)
    CALL sync
    
    ! Define the merging table
    T(  1) = 1
    T(  2) = 1
    T(  3) = 1
    T(  4) = 1
    T(  5) = 2
    T(  6) = 2
    T(  7) = 3
    T(  8) = 3
    T(  9) = 3
    T( 10) = 4
    T( 11) = 4
    T( 12) = 4
    T( 13) = 5
    T( 14) = 6
    T( 15) = 6
    T( 16) = 7
    T( 17) = 7
    T( 18) = 8
    T( 19) = 8
    
    ! Perform the merge
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      DO k = 1, 19
        IF (basin_ID_old( j,i) == k) basin_ID( j,i) = T( k)
      END DO
    END DO
    END DO
    IF (par%master) nbasins = 8
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wbasin_ID_old)
    
  END SUBROUTINE merge_basins_GRL
  
! == The no-ice mask, to prevent ice growth in certain areas
  SUBROUTINE initialise_mask_noice( region)
    ! Mask a certain area where no ice is allowed to grow. This is used to "remove"
    ! Greenland from NAM and EAS, and Ellesmere Island from GRL.
    ! 
    ! Also used to define calving fronts in certain idealised-geometry experiments
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_model_region),             INTENT(INOUT) :: region
    
    ! Allocate shared memory
    CALL allocate_shared_int_2D( region%grid%ny, region%grid%nx, region%mask_noice, region%wmask_noice)
    
    ! Initialise
    region%mask_noice( :,region%grid%i1:region%grid%i2) = 0
    CALL sync
    
    IF     (region%name == 'NAM') THEN
      ! Define a no-ice mask for North America
    
      IF     (C%choice_mask_noice_NAM == 'none') THEN
        ! No no-ice mask is defined for North America
      ELSEIF (C%choice_mask_noice_NAM == 'NAM_remove_GRL') THEN
        ! Prevent ice growth in the Greenlandic part of the North America domain
        CALL initialise_mask_noice_NAM_remove_GRL( region%grid, region%mask_noice)
      ELSE
        IF (par%master) WRITE(0,*) 'initialise_mask_noice - ERROR: unknown choice_mask_noice_NAM "', TRIM(C%choice_mask_noice_NAM), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSEIF (region%name == 'EAS') THEN
      ! Define a no-ice mask for Eurasia
    
      IF     (C%choice_mask_noice_EAS == 'none') THEN
        ! No no-ice mask is defined for Eurasia
      ELSEIF (C%choice_mask_noice_EAS == 'EAS_remove_GRL') THEN
        ! Prevent ice growth in the Greenlandic part of the Eurasia domain
        CALL initialise_mask_noice_EAS_remove_GRL( region%grid, region%mask_noice)
      ELSE
        IF (par%master) WRITE(0,*) 'initialise_mask_noice - ERROR: unknown choice_mask_noice_EAS "', TRIM(C%choice_mask_noice_EAS), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSEIF (region%name == 'GRL') THEN
      ! Define a no-ice mask for Greenland
    
      IF     (C%choice_mask_noice_GRL == 'none') THEN
        ! No no-ice mask is defined for Greenland
      ELSEIF (C%choice_mask_noice_GRL == 'GRL_remove_Ellesmere') THEN
        ! Prevent ice growth in the Ellesmere Island part of the Greenland domain
        CALL initialise_mask_noice_GRL_remove_Ellesmere( region%grid, region%mask_noice)
      ELSE
        IF (par%master) WRITE(0,*) 'initialise_mask_noice - ERROR: unknown choice_mask_noice_GRL "', TRIM(C%choice_mask_noice_GRL), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSEIF (region%name == 'ANT') THEN
      ! Define a no-ice mask for Antarctica, or for an idealised-geometry experiment
    
      IF     (C%choice_mask_noice_ANT == 'none') THEN
        ! No no-ice mask is defined for Antarctica
      ELSEIF (C%choice_mask_noice_ANT == 'MISMIP_mod') THEN
        ! Confine ice to the circular shelf around the cone-shaped island of the MISMIP_mod idealised geometry
        CALL initialise_mask_noice_MISMIP_mod( region%grid, region%mask_noice)
      ELSEIF (C%choice_mask_noice_ANT == 'MISMIP+') THEN
        ! Enforce the static calving front at x = 640 km in the MISMIP+ idealised geometry
        CALL initialise_mask_noice_MISMIPplus( region%grid, region%mask_noice)
      ELSE
        IF (par%master) WRITE(0,*) 'initialise_mask_noice - ERROR: unknown choice_mask_noice_ANT "', TRIM(C%choice_mask_noice_ANT), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END IF
  
  END SUBROUTINE initialise_mask_noice
  SUBROUTINE initialise_mask_noice_NAM_remove_GRL( grid, mask_noice)
    ! Prevent ice growth in the Greenlandic part of the North America domain
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    INTEGER,  DIMENSION(:,:  ),          INTENT(OUT)   :: mask_noice
  
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(2)                             :: pa, pb
    REAL(dp)                                           :: yl_ab
      
    pa = [ 490000._dp, 1530000._dp]
    pb = [2030000._dp,  570000._dp]
    
    DO i = grid%i1, grid%i2
      yl_ab = pa(2) + (grid%x(i) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
      DO j = 1, grid%ny
        IF (grid%y(j) > yl_ab .AND. grid%x(i) > pa(1) .AND. grid%y(j) > pb(2)) THEN
          mask_noice( j,i) = 1
        END IF
      END DO
    END DO
    CALL sync
  
  END SUBROUTINE initialise_mask_noice_NAM_remove_GRL
  SUBROUTINE initialise_mask_noice_EAS_remove_GRL( grid, mask_noice)
    ! Prevent ice growth in the Greenlandic part of the Eurasia domain
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    INTEGER,  DIMENSION(:,:  ),          INTENT(OUT)   :: mask_noice
  
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc, pd
    REAL(dp)                                           :: yl_ab, yl_bc, yl_cd
    
    pa = [-2900000._dp, 1300000._dp]
    pb = [-1895000._dp,  900000._dp]
    pc = [ -835000._dp, 1135000._dp]
    pd = [ -400000._dp, 1855000._dp]
    
    DO i = grid%i1, grid%i2
      yl_ab = pa(2) + (grid%x(i) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
      yl_bc = pb(2) + (grid%x(i) - pb(1))*(pc(2)-pb(2))/(pc(1)-pb(1))
      yl_cd = pc(2) + (grid%x(i) - pc(1))*(pd(2)-pc(2))/(pd(1)-pc(1))
      DO j = 1, grid%ny
        IF ((grid%x(i) <  pa(1) .AND. grid%y(j) > pa(2)) .OR. &
            (grid%x(i) >= pa(1) .AND. grid%x(i) < pb(1) .AND. grid%y(j) > yl_ab) .OR. &
            (grid%x(i) >= pb(1) .AND. grid%x(i) < pc(1) .AND. grid%y(j) > yl_bc) .OR. &
            (grid%x(i) >= pc(1) .AND. grid%x(i) < pd(1) .AND. grid%y(j) > yl_cd)) THEN
          mask_noice( j,i) = 1
        END IF
      END DO
    END DO
    CALL sync
  
  END SUBROUTINE initialise_mask_noice_EAS_remove_GRL
  SUBROUTINE initialise_mask_noice_GRL_remove_Ellesmere( grid, mask_noice)
    ! Prevent ice growth in the Ellesmere Island part of the Greenland domain
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    INTEGER,  DIMENSION(:,:  ),          INTENT(OUT)   :: mask_noice
  
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(2)                             :: pa, pb
    REAL(dp)                                           :: yl_ab
      
    pa = [-750000._dp,  900000._dp]
    pb = [-250000._dp, 1250000._dp]
    
    DO i = grid%i1, grid%i2
      yl_ab = pa(2) + (grid%x(i) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
      DO j = 1, grid%ny
        IF (grid%y(j) > pa(2) .AND. grid%y(j) > yl_ab .AND. grid%x(i) < pb(1)) THEN
          mask_noice( j,i) = 1
        END IF
      END DO
    END DO
    CALL sync
  
  END SUBROUTINE initialise_mask_noice_GRL_remove_Ellesmere
  SUBROUTINE initialise_mask_noice_MISMIP_mod( grid, mask_noice)
    ! Confine ice to the circular shelf around the cone-shaped island of the MISMIP_mod idealised-geometry experiment
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    INTEGER,  DIMENSION(:,:  ),          INTENT(OUT)   :: mask_noice
  
    ! Local variables:
    INTEGER                                            :: i,j
        
    ! Create a nice circular ice shelf
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (SQRT(grid%x(i)**2+grid%y(j)**2) > grid%xmax * 0.95_dp) THEN
        mask_noice( j,i) = 1
      END IF
    END DO
    END DO
    CALL sync
  
  END SUBROUTINE initialise_mask_noice_MISMIP_mod
  SUBROUTINE initialise_mask_noice_MISMIPplus( grid, mask_noice)
    ! Enforce the static calving front at x = 640 km in the MISMIP+ idealised geometry
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    INTEGER,  DIMENSION(:,:  ),          INTENT(OUT)   :: mask_noice
  
    ! Local variables:
    INTEGER                                            :: i,j
        
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ! NOTE: because IMAU-ICE wants to centre the domain at x=0, the front now lies at x = 240 km
      IF (grid%x( i) > 240000._dp) THEN
        mask_noice( j,i) = 1
      END IF
    END DO
    END DO
    CALL sync
  
  END SUBROUTINE initialise_mask_noice_MISMIPplus
  
! == Automatically tuning the ice flow factor A for the grounding-line position in the MISMIPplus experiment
  SUBROUTINE MISMIPplus_adapt_flow_factor( grid, ice)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,jmid
    REAL(dp)                                           :: TAF1,TAF2,lambda_GL, x_GL
    REAL(dp)                                           :: A_flow_old, f, A_flow_new
    
    IF (par%master) THEN
      
      ! Determine mid-channel grounding-line position
      jmid = CEILING( REAL(grid%ny,dp) / 2._dp)
      i = 1
      DO WHILE (ice%mask_sheet_a( jmid,i) == 1 .AND. i < grid%nx)
        i = i+1
      END DO
      
      TAF1 = ice%TAF_a( jmid,i-1)
      TAF2 = ice%TAF_a( jmid,i  )
      lambda_GL = TAF1 / (TAF1 - TAF2)
      x_GL = lambda_GL * grid%x( i) + (1._dp - lambda_GL) * grid%x( i-1)
      
      ! Adjust for the fact that the IMAU-ICE coordinate system is different than the one used in MISMIPplus
      x_GL = x_GL + 400000._dp
      
      ! Adjust the flow factor
      A_flow_old = ice%A_flow_vav_a( 1,1)
      f = 2._dp ** ((x_GL - C%MISMIPplus_xGL_target) / 80000._dp)
      A_flow_new = A_flow_old * f
      
      C%uniform_flow_factor = A_flow_new
      
      IF (par%master) WRITE(0,*) '    MISMIPplus_adapt_flow_factor: x_GL = ', x_GL/1E3, ' km; changed flow factor from ', A_flow_old, ' to ', A_flow_new
      
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE MISMIPplus_adapt_flow_factor

END MODULE general_ice_model_data_module
