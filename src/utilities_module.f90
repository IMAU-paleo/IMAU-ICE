MODULE utilities_module

  ! Contains some general "utilities" that are used at multiple points in the model.

  USE mpi
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, partition_list, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE data_types_module,               ONLY: type_grid, type_sparse_matrix_CSR
  USE petsc_module,                    ONLY: solve_matrix_equation_CSR_PETSc

CONTAINS

! == Some operations on the scaled vertical coordinate
  FUNCTION vertical_integration_from_bottom_to_zeta( f) RESULT( integral_f)
    ! This subroutine integrates f from the bottom level at C%zeta(k=C%NZ) = 1 up to the level C%zeta(k):
    !  See Eq. (12.1)
    ! If the integrand f is positive (our case) the integral is negative because the integration is in
    ! the opposite zeta direction. A 1D array which contains for each k-layer the integrated value from
    ! the bottom up to that k-layer is returned. The value of the integrand f at some integration step k
    ! is the average of f(k+1) and f(k):
    !  integral_f(k) = integral_f(k+1) + 0.5 * (f(k+1) + f(k)) * (-dzeta)
    ! with dzeta = C%zeta(k+1) - C%zeta(k). So for f > 0  integral_f < 0.
    
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%NZ), INTENT(IN)  :: f

    ! Output variables:
    REAL(dp), DIMENSION(C%NZ)              :: integral_f
    
    ! Local variables:
    INTEGER                                :: k

    integral_f(C%NZ) = 0._dp
    DO k = C%NZ-1, 1, -1
      integral_f(k) = integral_f(k+1) - 0.5_dp * (f(k+1) + f(k)) * (C%zeta(k+1) - C%zeta(k))
    END DO
  END FUNCTION vertical_integration_from_bottom_to_zeta
  FUNCTION vertical_integration_from_top_to_zeta(    f) RESULT( integral_f)
    ! This subroutine integrates f from the top level at C%zeta(k=1) = 0 down to the level C%zeta(k): Eq. (12.2)
    ! Similar to Eq. (12.1) but in the other direction.
    ! If the integrand f is positive (our case) the integral is positive because the integration is in
    ! the zeta direction. A 1D array which contains for each k-layer the integrated value from
    ! the top down to that k-layer is returned. The value of the integrand f at some integration step k
    ! is the average of f(k) and f(k-1):
    ! integral_f(k) = integral_f(k-1) + 0.5 * (f(k) + f(k-1)) * (dzeta); with dzeta = C%zeta(k+1) - C%zeta(k). 
    ! Heiko Goelzer (h.goelzer@uu.nl) Jan 2016

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%NZ), INTENT(IN)  :: f

    ! Output variables:
    REAL(dp), DIMENSION(C%NZ)              :: integral_f
    
    ! Local variables:
    INTEGER                                :: k

    integral_f(1) = 0._dp
    DO k = 2, C%NZ, 1
      integral_f(k) = integral_f(k-1) + 0.5_dp * (f(k) + f(k-1)) * (C%zeta(k) - C%zeta(k-1))
    END DO
  END FUNCTION vertical_integration_from_top_to_zeta
  FUNCTION vertical_integrate(                       f) RESULT( integral_f)
    ! Integrate f over the ice column (from the base to the surface)
    
    IMPLICIT NONE

    ! Input variable:
    REAL(dp), DIMENSION(C%nz), INTENT(IN) :: f

    ! Output variable:
    REAL(dp)                              :: integral_f

    ! Local variable:
    INTEGER                               :: k

    ! Initial value is zero
    integral_f = 0.0_dp 

    ! Intermediate values include sum of all previous values 
    ! Take current value as average between points
    DO k = 2, C%nz
       integral_f = integral_f + 0.5_dp*(f(k)+f(k-1))*(C%zeta(k) - C%zeta(k-1))
    END DO

  END FUNCTION vertical_integrate
  FUNCTION vertical_average(                         f) RESULT( average_f)
    ! Calculate the vertical average of any given function f defined at the vertical zeta grid.
    !  See Eq. (11.3) in DOCUMENTATION/icedyn-documentation/icedyn-documentation.pdf.
    ! The integration is in the direction of the positive zeta-axis from C%zeta(k=1) = 0 up to C%zeta(k=C%NZ) = 1.
    ! Numerically: de average between layer k and k+1 is calculated and multiplied by the distance between those 
    ! layers k and k+1, which is imediately the weight factor for this contribution because de total layer distance 
    ! is scaled to 1. The sum of all the weighted contribution gives average_f the vertical average of f.
    
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%NZ), INTENT(IN) :: f

    ! Result variables:
    REAL(dp)                              :: average_f
    
    ! Local variables:
    INTEGER                               :: k

    !  See Eq. (11.4) in DOCUMENTATION/icedyn-documentation/icedyn-documentation.pdf
    average_f = 0._dp
    DO k = 1, C%NZ-1
      average_f = average_f + 0.5_dp * (f(k+1) + f(k)) * (C%zeta(k+1) - C%zeta(k))
    END DO

    RETURN
  END FUNCTION vertical_average
  
! == Floatation criterion, surface elevation, and thickness above floatation
  FUNCTION is_floating( Hi, Hb, SL) RESULT( isso)
    ! The flotation criterion
      
    IMPLICIT NONE
    
    REAL(dp),                            INTENT(IN)    :: Hi, Hb, SL
    LOGICAL                                            :: isso
    
    isso = .FALSE.
    IF (Hi < (SL - Hb) * seawater_density/ice_density) isso = .TRUE.
    
  END FUNCTION is_floating
  FUNCTION surface_elevation( Hi, Hb, SL) RESULT( Hs)
    ! The surface elevation equation
      
    IMPLICIT NONE
    
    REAL(dp),                            INTENT(IN)    :: Hi, Hb, SL
    REAL(dp)                                           :: Hs
    
    Hs = Hi + MAX( SL - ice_density / seawater_density * Hi, Hb)
    
  END FUNCTION surface_elevation
  FUNCTION thickness_above_floatation( Hi, Hb, SL) RESULT( TAF)
    ! The thickness-above-floatation equation
      
    IMPLICIT NONE
    
    REAL(dp),                            INTENT(IN)    :: Hi, Hb, SL
    REAL(dp)                                           :: TAF
    
    TAF = Hi - MAX(0._dp, (SL - Hb) * (seawater_density / ice_density))
  
  END FUNCTION thickness_above_floatation
  FUNCTION Hi_from_Hb_and_Hs( Hb, Hs, SL) RESULT( Hi)
    ! Given Hb, Hs, and SL, calculate Hi
      
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: Hb, Hs, SL
    REAL(dp)                                           :: Hi
    
    ! Local variables:
    REAL(dp)                                           :: Hs_float
    
    ! Calculate surface elevation if the ice were exactly at floating thickness
    Hs_float = Hb + (SL - Hb) * seawater_density / ice_density
    
    IF (Hs > Hs_float) THEN
      ! Grounded ice
      Hi = Hs - Hb
    ELSE
      ! Floating ice
      Hi = (Hs - SL) / (1._dp - ice_density / seawater_density)
    END IF
    
  END FUNCTION Hi_from_Hb_and_Hs
  
! == The error function (used in the Roe&Lindzen precipitation model)
  SUBROUTINE error_function(X, ERR)
    ! Purpose: Compute error function erf(x)
    ! Input:   x   --- Argument of erf(x)
    ! Output:  ERR --- erf(x)
    
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: X
    
    ! Output variables:
    REAL(dp), INTENT(OUT) :: ERR
    
    ! Local variables:
    REAL(dp)              :: EPS
    REAL(dp)              :: X2
    REAL(dp)              :: ER
    REAL(dp)              :: R
    REAL(dp)              :: C0
    INTEGER               :: k
    
    EPS = 1.0E-15_dp
    X2  = X * X
    IF(ABS(X) < 3.5_dp) THEN
     ER = 1.0_dp
     R  = 1.0_dp
     DO k = 1, 50
       R  = R * X2 / (REAL(k, dp) + 0.5_dp)
       ER = ER+R
       IF(ABS(R) < ABS(ER) * EPS) THEN
        C0  = 2.0_dp / SQRT(pi) * X * EXP(-X2)
        ERR = C0 * ER
        EXIT
       END IF
     END DO
    ELSE
     ER = 1.0_dp
     R  = 1.0_dp
     DO k = 1, 12
       R  = -R * (REAL(k, dp) - 0.5_dp) / X2
       ER = ER + R
       C0  = EXP(-X2) / (ABS(X) * SQRT(pi))
       ERR = 1.0_dp - C0 * ER
       IF(X < 0.0_dp) ERR = -ERR
     END DO
    ENDIF

    RETURN
  END SUBROUTINE error_function
  
! == The oblique stereographic projection
  SUBROUTINE oblique_sg_projection( lambda, phi, lambda_M_deg, phi_M_deg, beta_deg, x_IM_P_prime, y_IM_P_prime, k_P)
    ! This subroutine projects with an oblique stereographic projection the longitude-latitude
    ! coordinates to a rectangular coordinate system, with coordinates (x,y).
    !
    ! For more information about M, beta_deg, the center of projection and the used
    ! projection method see: Reerink et al. (2010), Mapping technique of climate fields
    ! between GCM's and ice models, GMD

    ! For North and South Pole: lambda_M_deg = 0._dp, to generate the correct coordinate
    ! system, see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)            :: lambda        ! in degrees
    REAL(dp), INTENT(IN)            :: phi           ! in degrees

    ! Polar stereographic projection parameters
    REAL(dp), INTENT(IN)            :: lambda_M_deg  ! in degrees
    REAL(dp), INTENT(IN)            :: phi_M_deg     ! in degrees
    REAL(dp), INTENT(IN)            :: beta_deg      ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT)           :: x_IM_P_prime  ! in metres
    REAL(dp), INTENT(OUT)           :: y_IM_P_prime  ! in metres
    REAL(dp), INTENT(OUT), OPTIONAL :: k_P           ! Length scale factor [-],  k in Snyder (1987)

    ! Local variables:
    REAL(dp)                        :: alpha_deg     ! in degrees
    REAL(dp)                        :: phi_P         ! in radians
    REAL(dp)                        :: lambda_P      ! in radians
    REAL(dp)                        :: t_P_prime
    REAL(dp)                        :: lambda_M, phi_M, alpha
    
    ! Convert beta to alpha
    alpha_deg = 90._dp - beta_deg

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = (pi / 180._dp) * phi
    lambda_P = (pi / 180._dp) * lambda

    ! Convert projection parameters to radians:
    lambda_M = (pi / 180._dp) * lambda_M_deg
    phi_M    = (pi / 180._dp) * phi_M_deg
    alpha    = (pi / 180._dp) * alpha_deg

    ! See equation (2.6) or equation (A.56) in Reerink et al. (2010):
    t_P_prime = (1._dp + COS(alpha)) / (1._dp + COS(phi_P) * COS(phi_M) * COS(lambda_P - lambda_M) + SIN(phi_P) * SIN(phi_M))

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010):
    x_IM_P_prime =  earth_radius * (COS(phi_P) * SIN(lambda_P - lambda_M)) * t_P_prime
    y_IM_P_prime =  earth_radius * (SIN(phi_P) * COS(phi_M) - (COS(phi_P) * SIN(phi_M)) * COS(lambda_P - lambda_M)) * t_P_prime

    ! See equation (21-4) on page 157 in Snyder (1987):
    IF(PRESENT(k_P)) k_P = (1._dp + COS(alpha)) / (1._dp + SIN(phi_M) * SIN(phi_P) + COS(phi_M) * COS(phi_P) * COS(lambda_P - lambda_M))

  END SUBROUTINE oblique_sg_projection
  SUBROUTINE inverse_oblique_sg_projection( x_IM_P_prime, y_IM_P_prime, lambda_M_deg, phi_M_deg, beta_deg, lambda_P, phi_P)
    ! This subroutine projects with an inverse oblique stereographic projection the
    ! (x,y) coordinates to a longitude-latitude coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, alpha_deg, the center of projection and the used
    ! projection method see: Reerink et al. (2010), Mapping technique of climate fields
    ! between GCM's and ice models, GMD

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in metres
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in metres

    ! Polar stereographic projection parameters
    REAL(dp), INTENT(IN)  :: lambda_M_deg  ! in degrees
    REAL(dp), INTENT(IN)  :: phi_M_deg     ! in degrees
    REAL(dp), INTENT(IN)  :: beta_deg      ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees

    ! Local variables:
    REAL(dp)              :: alpha_deg     ! in degrees
    REAL(dp)              :: x_3D_P_prime  ! in metres
    REAL(dp)              :: y_3D_P_prime  ! in metres
    REAL(dp)              :: z_3D_P_prime  ! in metres
    REAL(dp)              :: a
    REAL(dp)              :: t_P
    REAL(dp)              :: x_3D_P        ! in metres
    REAL(dp)              :: y_3D_P        ! in metres
    REAL(dp)              :: z_3D_P        ! in metres
    REAL(dp)              :: lambda_M, phi_M, alpha
    
    ! Convert beta to alpha
    alpha_deg = 90._dp - beta_deg

    ! Convert projection parameters to radians:
    lambda_M = (pi / 180._dp) * lambda_M_deg
    phi_M    = (pi / 180._dp) * phi_M_deg
    alpha    = (pi / 180._dp) * alpha_deg

    ! See equations (2.14-2.16) or equations (B.21-B.23) in Reerink et al. (2010):
    x_3D_P_prime = earth_radius * COS(alpha) * COS(lambda_M) * COS(phi_M) - SIN(lambda_M) * x_IM_P_prime - COS(lambda_M) * SIN(phi_M) * y_IM_P_prime
    y_3D_P_prime = earth_radius * COS(alpha) * SIN(lambda_M) * COS(phi_M) + COS(lambda_M) * x_IM_P_prime - SIN(lambda_M) * SIN(phi_M) * y_IM_P_prime
    z_3D_P_prime = earth_radius * COS(alpha) *                 SIN(phi_M)                                +                 COS(phi_M) * y_IM_P_prime

    ! See equation (2.13) or equation (B.20) in Reerink et al. (2010):
    a = COS(lambda_M) * COS(phi_M) * x_3D_P_prime  +  SIN(lambda_M) * COS(phi_M) * y_3D_P_prime  +  SIN(phi_M) * z_3D_P_prime

    ! See equation (2.12) or equation (B.19) in Reerink et al. (2010):
    t_P = (2._dp * earth_radius**2 + 2._dp * earth_radius * a) / (earth_radius**2 + 2._dp * earth_radius * a + x_3D_P_prime**2 + y_3D_P_prime**2 + z_3D_P_prime**2)

    ! See equations (2.9-2.11) or equations (B.16-B.18) in Reerink et al. (2010):
    x_3D_P =  earth_radius * COS(lambda_M) * COS(phi_M) * (t_P - 1._dp) + x_3D_P_prime * t_P
    y_3D_P =  earth_radius * SIN(lambda_M) * COS(phi_M) * (t_P - 1._dp) + y_3D_P_prime * t_P
    z_3D_P =  earth_radius *                 SIN(phi_M) * (t_P - 1._dp) + z_3D_P_prime * t_P

    ! See equation (2.7) or equation (B.24) in Reerink et al. (2010):
    IF(x_3D_P <  0._dp                      ) THEN
      lambda_P = 180._dp + (180._dp / pi) * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P >= 0._dp) THEN
      lambda_P =           (180._dp / pi) * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P <  0._dp) THEN
      lambda_P = 360._dp + (180._dp / pi) * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P >  0._dp) THEN
      lambda_P =  90._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P <  0._dp) THEN
      lambda_P = 270._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P == 0._dp) THEN
      lambda_P =   0._dp
    END IF

    ! See equation (2.8) or equation (B.25) in Reerink et al. (2010):
    IF(x_3D_P /= 0._dp .OR. y_3D_P /= 0._dp) THEN
      phi_P = (180._dp / pi) * ATAN(z_3D_P / sqrt(x_3D_P**2 + y_3D_P**2))
    ELSE IF(z_3D_P >  0._dp) THEN
      phi_P =   90._dp
    ELSE IF(z_3D_P <  0._dp) THEN
      phi_P =  -90._dp
    END IF

  END SUBROUTINE inverse_oblique_sg_projection

! == Smoothing operations
  SUBROUTINE smooth_Gaussian_2D( grid, d, r)
    ! Apply a Gaussian smoothing filter of with sigma = n*dx to the 2D data field d
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d
    REAL(dp),                            INTENT(IN)    :: r      ! Smoothing radius in m
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Gaussian_2D'
    INTEGER                                            :: i,j,n,k,ii,jj
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_smooth
    INTEGER                                            :: wd_smooth
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: f
    
    n = CEILING( r / grid%dx) * 3  ! Number of cells to extend the data by (3 standard deviations is enough to capture the tails of the normal distribution)
    
    ! Fill in the smoothing filter
    ALLOCATE( f( -n:n))
    f = 0._dp
    DO i = -n, n
      f(i) = EXP( -0.5_dp * (REAL(i,dp) * grid%dx/r)**2)
    END DO
    f = f / SUM(f)
    
    ! Allocate temporary shared memory for the smoothed data field
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, d_smooth, wd_smooth)
    d_smooth( :,grid%i1:grid%i2) = 0._dp
    CALL sync
    
    ! Smooth in x-direction
    DO i = grid%i1, grid%i2
    DO j = 1,       grid%ny
      DO k = -n, n
        ii = MAX( 1, MIN( grid%nx, i + k))
        d_smooth( j,i) = d_smooth( j,i) + d( j,ii) * f(k)
      END DO
    END DO
    END DO
    CALL sync
    
    ! Copy x-smoothed data
    d( :,grid%i1:grid%i2) = d_smooth( :,grid%i1:grid%i2)
    CALL sync
    
    ! Reset smoothed data field
    d_smooth( :,grid%i1:grid%i2) = 0._dp
    CALL sync
    
    ! Smooth in y-direction
    DO i = grid%i1, grid%i2
    DO j = 1,       grid%ny
      DO k = -n, n
        jj = MAX( 1, MIN( grid%ny, j + k))
        d_smooth( j,i) = d_smooth( j,i) + d( jj,i) * f(k)
      END DO
    END DO
    END DO
    CALL sync
    
    ! Copy smoothed data
    d( :,grid%i1:grid%i2) = d_smooth( :,grid%i1:grid%i2)
    CALL sync
    
    ! Clean up after yourself
    DEALLOCATE( f)
    CALL deallocate_shared( wd_smooth)
    
  END SUBROUTINE smooth_Gaussian_2D
  SUBROUTINE smooth_Gaussian_3D( grid, d, r, nz)
    ! Apply a Gaussian smoothing filter of with sigma = n*dx to the 3D data field d
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d
    REAL(dp),                            INTENT(IN)    :: r      ! Smoothing radius in km
    INTEGER,                             INTENT(IN)    :: nz
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Gaussian_3D'
    INTEGER                                            :: k
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_2D
    INTEGER                                            :: wd_2D
    
    ! Allocate temporary shared memory for the extended and smoothed data fields
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, d_2D, wd_2D)
    
    DO k = 1, nz
      d_2D( :,grid%i1:grid%i2) = d( k,:,grid%i1:grid%i2)
      CALL smooth_Gaussian_2D( grid, d_2D, r)
      d( k,:,grid%i1:grid%i2) = d_2D( :,grid%i1:grid%i2)
    END DO
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_2D)
    
  END SUBROUTINE smooth_Gaussian_3D
  SUBROUTINE smooth_Shepard_2D( grid, d, r)
    ! Apply a Shepard smoothing filter of with sigma = n*dx to the 2D data field d
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d
    REAL(dp),                            INTENT(IN)    :: r      ! Smoothing radius in m
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Shepard_2D'
    INTEGER                                            :: i,j,k,l,n
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_ext,  d_ext_smooth
    INTEGER                                            :: wd_ext, wd_ext_smooth
    REAL(dp)                                           :: ShepNumSum     ! The sum of the numerators in the Shepard weighting
    REAL(dp)                                           :: ShepDenSum     ! The sum of the denumerators in the Shepard weighting
    REAL(dp)                                           :: distance       ! in gridsize units
    REAL(dp)                                           :: smooth_radius  ! in gridsize units
    REAL(dp)                                           :: exponent       ! The distance weighting exponent in the Shepard weighting
    
    n = CEILING( r / grid%dx) ! Number of cells to extend the data by (3 standard deviations is enough to capture the tails of the normal distribution)
    
    smooth_radius = REAL(n,dp)
    exponent      = 2._dp
    
    ! Allocate temporary shared memory for the extended and smoothed data fields
    CALL allocate_shared_dp_2D( grid%ny + 2*n, grid%nx + 2*n, d_ext,        wd_ext       )
    CALL allocate_shared_dp_2D( grid%ny + 2*n, grid%nx + 2*n, d_ext_smooth, wd_ext_smooth)
    
    ! Copy data to the extended array and fill in the margins
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      d_ext( j+n,i+n) = d( j,i)
    END DO
    END DO
    IF (par%master) THEN
      ! West
      d_ext( n+1:n+grid%ny, 1            ) = d( :      ,1      )
      ! East
      d_ext( n+1:n+grid%ny, grid%nx+2*n  ) = d( :      ,grid%nx)
      ! South
      d_ext( 1            , n+1:n+grid%nx) = d( 1      ,:      )
      ! North
      d_ext( grid%ny+2*n  , n+1:n+grid%nx) = d( grid%ny,:      )
      ! Corners
      d_ext( 1:n,                     1:n                    ) = d( 1      ,1      )
      d_ext( 1:n,                     grid%nx+n+1:grid%nx+2*n) = d( 1      ,grid%nx)
      d_ext( grid%ny+n+1:grid%ny+2*n, 1:n                    ) = d( grid%ny,1      )
      d_ext( grid%ny+n+1:grid%ny+2*n, grid%nx+n+1:grid%nx+2*n) = d( grid%ny,grid%nx)
    END IF
    CALL sync
    
    ! Convolute extended data with the smoothing filter
    d_ext_smooth( :,grid%i1+n:grid%i2+n) = 0._dp
    CALL sync
    
    DO i = grid%i1, grid%i2
    DO j = 1,       grid%ny
      ShepNumSum = 0._dp
      ShepDenSum = 0._dp
      DO k = -n, n
      DO l = -n, n
        distance = SQRT(REAL(k,dp)**2 + REAL(l,dp)**2)
        IF (distance <= smooth_radius .AND. distance > 0._dp) THEN
          ShepNumSum = ShepNumSum + (d_ext( j+n+l,i+n+k)/(distance**exponent))
          ShepDenSum = ShepDenSum + (1._dp/(distance**exponent))
        END IF
      END DO
      END DO
      d_ext_smooth( j+n,i+n) = ShepNumSum/ShepDenSum
    END DO
    END DO
    CALL sync
    
    ! Copy data back
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      d( j,i) = d_ext_smooth( j+n, i+n)
    END DO
    END DO
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_ext)
    CALL deallocate_shared( wd_ext_smooth)
    
  END SUBROUTINE smooth_Shepard_2D
  SUBROUTINE smooth_Shepard_3D( grid, d, r, nz)
    ! Apply a Shepard smoothing filter of with sigma = n*dx to the 3D data field d
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d
    REAL(dp),                            INTENT(IN)    :: r      ! Smoothing radius in km
    INTEGER,                             INTENT(IN)    :: nz
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Shepard_3D'
    INTEGER                                            :: k
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_2D
    INTEGER                                            :: wd_2D
    
    ! Allocate temporary shared memory for the extended and smoothed data fields
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, d_2D, wd_2D)
    
    DO k = 1, nz
      d_2D( :,grid%i1:grid%i2) = d( k,:,grid%i1:grid%i2)
      CALL smooth_Shepard_2D( grid, d_2D, r)
      d( k,:,grid%i1:grid%i2) = d_2D( :,grid%i1:grid%i2)
    END DO
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_2D)
    
  END SUBROUTINE smooth_Shepard_3D
  
! == Solve matrix equations in CSR format
  SUBROUTINE solve_matrix_equation_CSR( CSR, choice_matrix_solver, SOR_nit, SOR_tol, SOR_omega, PETSc_rtol, PETSc_abstol)
    ! Solve the matrix equation Ax = b
    ! The matrix A is provided in Compressed Sparse Row (CSR) format
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: CSR
    CHARACTER(LEN=256),                  INTENT(IN)    :: choice_matrix_solver
    INTEGER,                             INTENT(IN)    :: SOR_nit
    REAL(dp),                            INTENT(IN)    :: SOR_tol
    REAL(dp),                            INTENT(IN)    :: SOR_omega
    REAL(dp),                            INTENT(IN)    :: PETSc_rtol
    REAL(dp),                            INTENT(IN)    :: PETSc_abstol
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_matrix_equation_CSR'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Safety
    CALL check_for_NaN_dp_1D(  CSR%A_val, 'CSR%A_val')
    CALL check_for_NaN_dp_1D(  CSR%b,     'CSR%b'    )
    CALL check_for_NaN_dp_1D(  CSR%x,     'CSR%x'    )
    
    IF (choice_matrix_solver == 'SOR') THEN
      ! Use the old simple SOR solver
      
      CALL solve_matrix_equation_CSR_SOR( CSR, SOR_nit, SOR_tol, SOR_omega)
      
    ELSEIF (choice_matrix_solver == 'PETSc') THEN
      ! Use the PETSc solver (much preferred, this is way faster and more stable!)
    
      CALL solve_matrix_equation_CSR_PETSc( CSR, PETSc_rtol, PETSc_abstol)
    
    ELSE
      CALL crash('unknown choice_matrix_solver "' // TRIM( choice_matrix_solver) // '"!')
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE solve_matrix_equation_CSR
  SUBROUTINE solve_matrix_equation_CSR_SOR( CSR, nit, tol, omega)
    ! Solve the matrix equation Ax = b using successive over-relaxation (SOR)
    ! The matrix A is provided in Compressed Sparse Row (CSR) format
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: CSR
    INTEGER,                             INTENT(IN)    :: nit
    REAL(dp),                            INTENT(IN)    :: tol
    REAL(dp),                            INTENT(IN)    :: omega
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_matrix_equation_CSR_SOR'
    INTEGER                                            :: i,j,k,it,i1,i2
    REAL(dp)                                           :: lhs, res, cij, res_max, omega_dyn
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Partition equations over the processors
    CALL partition_list( CSR%m, par%i, par%n, i1, i2)
    
    omega_dyn = omega
    
    res_max = tol * 2._dp
    it = 0
    SOR_iterate: DO WHILE (res_max > tol .AND. it < nit)
      it = it+1
      
      res_max = 0._dp

      DO i = i1, i2
      
        lhs = 0._dp
        cij = 0._dp
        DO k = CSR%A_ptr( i), CSR%A_ptr( i+1)-1
          j = CSR%A_index( k)
          lhs = lhs + CSR%A_val( k) * CSR%x( j)
          IF (j == i) cij = CSR%A_val( k)
        END DO
        
        res = (lhs - CSR%b( i)) / cij
        res_max = MAX( res_max, ABS(res))
        
        CSR%x( i) = CSR%x( i) - omega_dyn * res
        
      END DO ! DO i = i1, i2
      CALL sync
      
      ! Check if we've reached a stable solution
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, res_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      
      !IF (par%master) WRITE(0,*) '      SOR iteration ', it, ': res_max = ', res_max
      
      IF (it > 100 .AND. res_max > 1E3_dp ) THEN
        
        ! Divergence detected - decrease omega, reset solution to zero, restart SOR.
        IF (par%master) WRITE(0,*) '  solve_matrix_equation_CSR_SOR - divergence detected; decrease omega, reset solution to zero, restart SOR'
        omega_dyn = omega_dyn - 0.1_dp
        it = 0
        CSR%x( i1:i2) = 0._dp
        CALL sync
        
        IF (omega_dyn <= 0.1_dp) THEN
          CALL crash('divergence detected even with extremely low relaxation parameter!')
        END IF
      END IF
      
    END DO SOR_iterate
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE solve_matrix_equation_CSR_SOR
  SUBROUTINE initialise_matrix_equation_CSR( CSR, m, n, nnz_per_row_max)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: CSR
    INTEGER,                             INTENT(IN)    :: m,n,nnz_per_row_max
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_matrix_equation_CSR'
    
    CALL allocate_shared_int_0D( CSR%m,               CSR%wm              )
    CALL allocate_shared_int_0D( CSR%n,               CSR%wn              )
    CALL allocate_shared_int_0D( CSR%nnz_per_row_max, CSR%wnnz_per_row_max)
    CALL allocate_shared_int_0D( CSR%nnz_max,         CSR%wnnz_max        )
    
    IF (par%master) THEN
      CSR%m               = m
      CSR%n               = n
      CSR%nnz_per_row_max = nnz_per_row_max
      CSR%nnz_max         = nnz_per_row_max * m
    END IF
    CALL sync
    
    CALL allocate_shared_int_1D( CSR%m+1,     CSR%A_ptr,   CSR%wA_ptr  )
    CALL allocate_shared_int_1D( CSR%nnz_max, CSR%A_index, CSR%wA_index)
    CALL allocate_shared_dp_1D(  CSR%nnz_max, CSR%A_val,   CSR%wA_val  )
    CALL allocate_shared_dp_1D(  CSR%m,       CSR%b,       CSR%wb      )
    CALL allocate_shared_dp_1D(  CSR%n,       CSR%x,       CSR%wx      )
    
  END SUBROUTINE initialise_matrix_equation_CSR
  SUBROUTINE check_CSR_for_double_entries( CSR)
    ! Check a CSR matrix representation for double entries
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: CSR
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_CSR_for_double_entries'
    INTEGER                                            :: i,j,k,i1,i2,k2,j2
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Partition equations over the processors
    CALL partition_list( CSR%m, par%i, par%n, i1, i2)

    DO i = i1, i2  
    
      DO k = CSR%A_ptr( i), CSR%A_ptr( i+1)-1
        j = CSR%A_index( k)
        
        DO k2 = CSR%A_ptr( i), CSR%A_ptr( i+1)-1
          IF (k2 == k) CYCLE
          j2 = CSR%A_index( k2)
          IF (j2 == j) THEN
            CALL crash('double entry detected!')
          END IF
        END DO
        
      END DO
      
    END DO
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE check_CSR_for_double_entries
  
! == Solve a tridiagonal matrix equation using LAPACK
  FUNCTION tridiagonal_solve( ldiag, diag, udiag, rhs) RESULT(x)
    ! Lapack tridiagnal solver (in double precision):
    ! Matrix system solver for tridiagonal matrices. 
    ! Used e.g. in solving the ADI scheme. 
    ! ldiag = lower diagonal elements (j,j-1) of the matrix
    ! diag  = diagonal elements (j,j) of the matrix
    ! udiag = upper diagonal elements (j,j+1) of the matrix
    ! rhs   = right hand side of the matrix equation in the ADI scheme
    USE configuration_module, ONLY: dp
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(:),            INTENT(IN) :: diag
    REAL(dp), DIMENSION(SIZE(diag)-1), INTENT(IN) :: udiag, ldiag
    REAL(dp), DIMENSION(SIZE(diag)),   INTENT(IN) :: rhs

    ! Result variable:
    REAL(dp), DIMENSION(SIZE(diag))               :: x
    
    ! Local variables:     
    INTEGER                                       :: info
    REAL(dp), DIMENSION(SIZE(diag))               :: diag_copy
    REAL(dp), DIMENSION(SIZE(udiag))              :: udiag_copy, ldiag_copy

    ! External subroutines:      
    EXTERNAL DGTSV ! Lapack routine that solves tridiagonal systems (in double precision).

    ! The LAPACK solver will overwrite the rhs with the solution x. Therefore we 
    ! first copy the rhs in the solution vector x:
    x = rhs

    ! The LAPACK solver will change the elements in the matrix, therefore we copy them:
    diag_copy  =  diag
    udiag_copy = udiag
    ldiag_copy = ldiag

    CALL DGTSV(SIZE(diag), 1, ldiag_copy, diag_copy, udiag_copy, x, SIZE(diag), info)
    
    IF (info /= 0) THEN
      CALL crash('tridiagonal_solve - ERROR: LAPACK solver DGTSV returned error message info = {int_01}', int_01 = info)
    END IF
    
  END FUNCTION tridiagonal_solve 
  
! == Analytical solution by Schoof 2006 for the "SSA_icestream" benchmark experiment
  SUBROUTINE SSA_Schoof2006_analytical_solution( tantheta, h0, A_flow, y, U, tauc)
      
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: tantheta   ! Surface slope in the x-direction
    REAL(dp),                            INTENT(IN)    :: h0         ! Ice thickness
    REAL(dp),                            INTENT(IN)    :: A_flow     ! Ice flow factor
    REAL(dp),                            INTENT(IN)    :: y          ! y-coordinate
    REAL(dp),                            INTENT(OUT)   :: U          ! Ice velocity in the x-direction
    REAL(dp),                            INTENT(OUT)   :: tauc       ! Till yield stress
    
    ! Local variables:
    REAL(dp)                                           :: B, m, L, f, W, ua, ub, uc, ud, ue
    
    ! Constants
    B = A_flow**(-1._dp/C%n_flow)
    m = C%SSA_icestream_m
    L = C%SSA_icestream_L
    
    ! Calculate the gravitational driving stress f
    f = ice_density * grav * h0 * tantheta    
    
    ! Calculate the "ice stream half-width" W
    W = L * (m+1._dp)**(1._dp/m)
    
    ! Calculate the till yield stress across the stream
    tauc = f * ABS(y/L)**m
    
    ! Calculate the analytical solution for u
    ua = -2._dp * f**3 * L**4 / (B**3 * h0**3)
    ub = ( 1._dp / 4._dp                           ) * (   (y/L)**     4._dp  - (m+1._dp)**(       4._dp/m) )
    uc = (-3._dp / ((m+1._dp)    * (      m+4._dp))) * (ABS(y/L)**(  m+4._dp) - (m+1._dp)**(1._dp+(4._dp/m)))
    ud = ( 3._dp / ((m+1._dp)**2 * (2._dp*m+4._dp))) * (ABS(y/L)**(2*m+4._dp) - (m+1._dp)**(2._dp+(4._dp/m)))
    ue = (-1._dp / ((m+1._dp)**3 * (3._dp*m+4._dp))) * (ABS(y/L)**(3*m+4._dp) - (m+1._dp)**(3._dp+(4._dp/m)))
    u = ua * (ub + uc + ud + ue)
    
    ! Outside the ice-stream, velocity is zero
    IF (ABS(y) > w) u = 0._dp
    
  END SUBROUTINE SSA_Schoof2006_analytical_solution
  
! == Map data between two square grids using 2nd-order conservative remapping
  SUBROUTINE map_square_to_square_cons_2nd_order_2D( nx_src, ny_src, x_src, y_src, nx_dst, ny_dst, x_dst, y_dst, d_src, d_dst)
    ! Map data from one square grid to another (e.g. PD ice thickness from the square grid in the input file to the model square grid)
    
    IMPLICIT NONE
  
    ! Input and output variables
    INTEGER,                            INTENT(IN)    :: nx_src
    INTEGER,                            INTENT(IN)    :: ny_src
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: x_src
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: y_src
    INTEGER,                            INTENT(IN)    :: nx_dst
    INTEGER,                            INTENT(IN)    :: ny_dst
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: x_dst
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: y_dst
    REAL(dp), DIMENSION(:,:  ),         INTENT(IN)    :: d_src
    REAL(dp), DIMENSION(:,:  ),         INTENT(OUT)   :: d_dst
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'map_square_to_square_cons_2nd_order_2D'
    INTEGER                                           :: i,j,i_src,j_src,i_dst,j_dst,i1,i2,igmin,igmax,jgmin,jgmax,j1,j2
    REAL(dp)                                          :: dx_src, dy_src, dx_dst, dy_dst
    INTEGER,  DIMENSION(nx_dst)                       :: il_src, iu_src
    INTEGER,  DIMENSION(ny_dst)                       :: jl_src, ju_src
    REAL(dp)                                          :: xomin, xomax, yomin, yomax, w0, w1x, w1y
    REAL(dp)                                          :: Ad, Asd, Asum
    REAL(dp), DIMENSION(:,:  ), POINTER               ::  ddx_src,  ddy_src
    INTEGER                                           :: wddx_src, wddy_src
    INTEGER,  DIMENSION(:,:  ), POINTER               ::  mask_dst_outside_src
    INTEGER                                           :: wmask_dst_outside_src
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D(  ny_src, nx_src, ddx_src,              wddx_src             )
    CALL allocate_shared_dp_2D(  ny_src, nx_src, ddy_src,              wddy_src             )
    CALL allocate_shared_int_2D( ny_dst, nx_dst, mask_dst_outside_src, wmask_dst_outside_src)
    
    ! Find grid spacings
    dx_src = x_src(2) - x_src(1)
    dy_src = y_src(2) - y_src(1)
    dx_dst = x_dst(2) - x_dst(1)
    dy_dst = y_dst(2) - y_dst(1)
    Ad = dx_dst * dy_dst
    
    ! If the grids are equal, the solution is trivial; just copy the data
    IF (dx_src == dx_dst .AND. dy_src == dy_dst .AND. nx_src == nx_dst .AND. ny_src == ny_dst) THEN
      CALL partition_list( nx_dst, par%i, par%n, i1, i2)
      d_dst( :,i1:i2) = d_src( :,i1:i2)
      CALL sync
      CALL deallocate_shared(  wddx_src             )
      CALL deallocate_shared(  wddy_src             )
      CALL deallocate_shared(  wmask_dst_outside_src)
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Find overlaps between grids
    ! dst cell [i,j] overlaps with src cells [il_src( i):iu_src( i), jl_src( j):ju_src( j)]
    
    ! il
    DO i_dst = 1, nx_dst
  
      IF (i_dst == 1) THEN
        il_src( i_dst) = 1
      ELSE
        il_src( i_dst) = il_src( i_dst-1)
      END IF
  
      DO WHILE (x_src( il_src( i_dst)) + dx_src / 2._dp < x_dst( i_dst) - dx_dst / 2._dp)
        il_src( i_dst) = il_src( i_dst) + 1
        IF (il_src( i_dst) >= nx_src) THEN
          EXIT
        END IF
      END DO
      il_src( i_dst) = MAX( 1, MIN( nx_src, il_src( i_dst) ))
  
    END DO ! DO i_dst = 1, nx_dst
    
    il_src = il_src - 1
    il_src = MAX( 1, MIN( nx_src, il_src))
    
    ! iu
    DO i_dst = nx_dst, 1, -1
  
      IF (i_dst == nx_dst) THEN
        iu_src( i_dst) = nx_src
      ELSE
        iu_src( i_dst) = il_src( i_dst+1)
      END IF
  
      DO WHILE (x_src( iu_src( i_dst)) - dx_src / 2._dp > x_dst( i_dst) + dx_dst / 2._dp)
        iu_src( i_dst) = iu_src( i_dst) - 1
        IF (iu_src( i_dst) <= 1) THEN
          EXIT
        END IF
      END DO
      iu_src( i_dst) = MAX( 1, MIN( nx_src, iu_src( i_dst) ))
  
    END DO ! DO i_dst = nx_dst, 1, -1
    
    iu_src = iu_src + 1
    iu_src = MAX( 1, MIN( nx_src, iu_src))
    
    ! jl
    DO j_dst = 1, ny_dst
  
      IF (j_dst == 1) THEN
        jl_src( j_dst) = 1
      ELSE
        jl_src( j_dst) = jl_src( j_dst-1)
      END IF
  
      DO WHILE (y_src( jl_src( j_dst)) + dx_src / 2._dp < y_dst( j_dst) - dx_dst / 2._dp)
        jl_src( j_dst) = jl_src( j_dst) + 1
        IF (jl_src( j_dst) >= ny_src) THEN
          EXIT
        END IF
      END DO
      jl_src( j_dst) = MAX( 1, MIN( ny_src, jl_src( j_dst) ))
  
    END DO ! DO j_dst = 1, ny_dst
    
    jl_src = jl_src - 1
    jl_src = MAX( 1, MIN( ny_src, jl_src))
    
    ! ju
    DO j_dst = ny_dst, 1, -1
  
      IF (j_dst == ny_dst) THEN
        ju_src( j_dst) = ny_src
      ELSE
        ju_src( j_dst) = jl_src( j_dst+1)
      END IF
  
      DO WHILE (y_src( ju_src( j_dst)) - dx_src / 2._dp > y_dst( j_dst) + dx_dst / 2._dp)
        ju_src( j_dst) = ju_src( j_dst) - 1
        IF (ju_src( j_dst) <= 1) THEN
          EXIT
        END IF
      END DO
      ju_src( j_dst) = MAX( 1, MIN( ny_src, ju_src( j_dst) ))
  
    END DO ! DO j_dst = ny_dst, 1, -1
    
    ju_src = ju_src + 1
    ju_src = MAX( 1, MIN( ny_src, ju_src))
    
    ! Get derivatives of d_src
    CALL partition_list( nx_src, par%i, par%n, i1, i2)
    DO i = MAX(2,i1), MIN(nx_src-1,i2)
    DO j = 2, ny_src-1
      ddx_src( j,i) = (d_src( j,i+1) - d_src( j,i-1)) / (2._dp * dx_src)
      ddy_src( j,i) = (d_src( j+1,i) - d_src( j-1,i)) / (2._dp * dy_src)
    END DO
    END DO
    CALL sync
    
    ! Find parallelisation domains
    CALL partition_list( nx_dst, par%i, par%n, i1, i2)
    CALL partition_list( ny_dst, par%i, par%n, j1, j2)
    
    DO i = i1, i2
    DO j = 1, ny_dst
      
      d_dst( j,i) = 0._dp
      Asum        = 0._dp
      
      ! If this dst cell lies (partly) outside of the src grid, mark it as such;
      ! in that case, use nearest-neighbour extrapolation instead of conservative remapping
      IF (x_dst( i) - dx_dst/2._dp < MINVAL( x_src) - dx_src/2._dp .OR. &
          x_dst( i) + dx_dst/2._dp > MAXVAL( x_src) + dx_src/2._dp .OR. &
          y_dst( j) - dx_dst/2._dp < MINVAL( y_src) - dx_src/2._dp .OR. &
          y_dst( j) + dx_dst/2._dp > MAXVAL( y_src) + dx_src/2._dp) THEN
        mask_dst_outside_src( j,i) = 1
        CYCLE
      ELSE
        mask_dst_outside_src( j,i) = 0
      END IF
      
      DO i_src = il_src( i), iu_src( i)
      DO j_src = jl_src( j), ju_src( j)
        
        xomin = MAX( x_dst( i) - dx_dst/2._dp, x_src( i_src) - dx_src/2._dp)
        xomax = MIN( x_dst( i) + dx_dst/2._dp, x_src( i_src) + dx_src/2._dp)
        yomin = MAX( y_dst( j) - dy_dst/2._dp, y_src( j_src) - dy_src/2._dp)
        yomax = MIN( y_dst( j) + dy_dst/2._dp, y_src( j_src) + dy_src/2._dp)
        
        IF (xomax <= xomin .OR. yomax <= yomin) CYCLE
        
        Asd  = (xomax - xomin) * (yomax - yomin)
        Asum = Asum + Asd
        
        w0  = Asd / Ad
        w1x = 1._dp / Ad * (line_integral_mxydx( [xomin,yomin], [xomax,yomin], 1E-9_dp) + &
                            line_integral_mxydx( [xomax,yomin], [xomax,yomax], 1E-9_dp) + &
                            line_integral_mxydx( [xomax,yomax], [xomin,yomax], 1E-9_dp) + &
                            line_integral_mxydx( [xomin,yomax], [xomin,yomin], 1E-9_dp)) - w0 * x_src( i_src)
        w1y = 1._dp / Ad * (line_integral_xydy(  [xomin,yomin], [xomax,yomin], 1E-9_dp) + &
                            line_integral_xydy(  [xomax,yomin], [xomax,yomax], 1E-9_dp) + &
                            line_integral_xydy(  [xomax,yomax], [xomin,yomax], 1E-9_dp) + &
                            line_integral_xydy(  [xomin,yomax], [xomin,yomin], 1E-9_dp)) - w0 * y_src( j_src)
        
        d_dst( j,i) = d_dst( j,i) + w0  * d_src(   j_src,i_src) + &
                                    w1x * ddx_src( j_src,i_src) + &
                                    w1y * ddy_src( j_src,i_src)
        
      END DO ! DO j_src = jr_src( j,1), jr_src( j,2)
      END DO ! DO i_src = ir_src( i,1), ir_src( i,2)
      
      ! Safety
      IF (ABS( 1._dp - Asum / Ad) > 1E-4_dp) THEN
        CALL crash('dst grid cell [{int_01},{int_02}] couldnt be completely filled! Asum = {dp_01}, Ad = {dp_02}', int_01 = j, int_02 = i, dp_01 = Asum, dp_02 = Ad)
      END IF
      
    END DO ! DO j = 1, ny_dst
    END DO ! DO i = i1, i2
    CALL sync
    
    ! Use nearest-neighbour extrapolation for dst cells outside of the src grid
    ! =========================================================================
    
    ! Find the range of grid cells that were mapped correctly
    igmin = 0
    igmax = 0
    jgmin = 0
    jgmax = 0
    
    j = INT( REAL(ny_dst,dp)/2._dp)
    DO i = 1, nx_dst
      IF (mask_dst_outside_src( j,i) == 0) THEN
        igmin = i
        EXIT
      END IF
    END DO
    DO i = nx_dst, 1, -1
      IF (mask_dst_outside_src( j,i) == 0) THEN
        igmax = i
        EXIT
      END IF
    END DO
    
    i = INT( REAL(nx_dst,dp)/2._dp)
    DO j = 1, ny_dst
      IF (mask_dst_outside_src( j,i) == 0) THEN
        jgmin = j
        EXIT
      END IF
    END DO
    DO j = ny_dst, 1, -1
      IF (mask_dst_outside_src( j,i) == 0) THEN
        jgmax = j
        EXIT
      END IF
    END DO
    
    ! Corners
    IF (par%master) THEN
      ! Southwest
      d_dst( 1      :jgmin-1 ,1      :igmin-1) = d_dst( jgmin,igmin)
      ! Southeast
      d_dst( 1      :jgmin-1 ,igmax+1:nx_dst ) = d_dst( jgmin,igmax)
      ! Northwest
      d_dst( jgmax+1:ny_dst  ,1      :igmin-1) = d_dst( jgmax,igmin)
      ! Northeast
      d_dst( jgmax+1:ny_dst  ,igmax+1:nx_dst ) = d_dst( jgmax,igmax)
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Borders
    DO i = MAX(i1,igmin), MIN(i2,igmax)
      ! South
      d_dst( 1      :jgmin-1,i) = d_dst( jgmin,i)
      ! North
      d_dst( jgmax+1:ny_dst ,i) = d_dst( jgmax,i)
    END DO
    DO j = MAX(j1,jgmin), MIN(j2,jgmax)
      ! West
      d_dst( j,1      :igmin-1) = d_dst( j,igmin)
      ! East
      d_dst( j,igmax+1:nx_dst ) = d_dst( j,igmax)
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wddx_src             )
    CALL deallocate_shared( wddy_src             )
    CALL deallocate_shared( wmask_dst_outside_src)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE map_square_to_square_cons_2nd_order_2D
  SUBROUTINE map_square_to_square_cons_2nd_order_3D( nx_src, ny_src, x_src, y_src, nx_dst, ny_dst, x_dst, y_dst, d_src, d_dst)
    ! Map data from one square grid to another (e.g. PD ice thickness from the square grid in the input file to the model square grid)
    
    IMPLICIT NONE
  
    ! Input and output variables
    INTEGER,                            INTENT(IN)    :: nx_src
    INTEGER,                            INTENT(IN)    :: ny_src
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: x_src
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: y_src
    INTEGER,                            INTENT(IN)    :: nx_dst
    INTEGER,                            INTENT(IN)    :: ny_dst
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: x_dst
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: y_dst
    REAL(dp), DIMENSION(:,:,:),         INTENT(IN)    :: d_src
    REAL(dp), DIMENSION(:,:,:),         INTENT(OUT)   :: d_dst
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'map_square_to_square_cons_2nd_order_3D'
    INTEGER                                           :: i,j,i_src,j_src,i_dst,j_dst,i1,i2,igmin,igmax,jgmin,jgmax,j1,j2,k
    REAL(dp)                                          :: dx_src, dy_src, dx_dst, dy_dst
    INTEGER,  DIMENSION(nx_dst)                       :: il_src, iu_src
    INTEGER,  DIMENSION(ny_dst)                       :: jl_src, ju_src
    REAL(dp)                                          :: xomin, xomax, yomin, yomax, w0, w1x, w1y
    REAL(dp)                                          :: Ad, Asd, Asum
    REAL(dp), DIMENSION(:,:,:), POINTER               ::  ddx_src,  ddy_src
    INTEGER                                           :: wddx_src, wddy_src
    INTEGER,  DIMENSION(:,:  ), POINTER               ::  mask_dst_outside_src
    INTEGER                                           :: wmask_dst_outside_src
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Safety
    IF (SIZE( d_src,1) /= SIZE( d_dst,1)) THEN
      CALL crash('data field sizes do not match!')
    END IF
    
    ! Allocate shared memory
    CALL allocate_shared_dp_3D(  SIZE( d_src,1), ny_src, nx_src, ddx_src,              wddx_src             )
    CALL allocate_shared_dp_3D(  SIZE( d_src,1), ny_src, nx_src, ddy_src,              wddy_src             )
    CALL allocate_shared_int_2D( ny_dst, nx_dst, mask_dst_outside_src, wmask_dst_outside_src)
    
    ! Find grid spacings
    dx_src = x_src(2) - x_src(1)
    dy_src = y_src(2) - y_src(1)
    dx_dst = x_dst(2) - x_dst(1)
    dy_dst = y_dst(2) - y_dst(1)
    Ad = dx_dst * dy_dst
    
    ! If the grids are equal, the solution is trivial; just copy the data
    IF (dx_src == dx_dst .AND. dy_src == dy_dst .AND. nx_src == nx_dst .AND. ny_src == ny_dst) THEN
      CALL partition_list( nx_dst, par%i, par%n, i1, i2)
      d_dst( :,:,i1:i2) = d_src( :,:,i1:i2)
      CALL sync
      CALL deallocate_shared(  wddx_src             )
      CALL deallocate_shared(  wddy_src             )
      CALL deallocate_shared(  wmask_dst_outside_src)
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Find overlaps between grids
    ! dst cell [i,j] overlaps with src cells [il_src( i):iu_src( i), jl_src( j):ju_src( j)]
    
    ! il
    DO i_dst = 1, nx_dst
  
      IF (i_dst == 1) THEN
        il_src( i_dst) = 1
      ELSE
        il_src( i_dst) = il_src( i_dst-1)
      END IF
  
      DO WHILE (x_src( il_src( i_dst)) + dx_src / 2._dp < x_dst( i_dst) - dx_dst / 2._dp)
        il_src( i_dst) = il_src( i_dst) + 1
        IF (il_src( i_dst) >= nx_src) THEN
          EXIT
        END IF
      END DO
      il_src( i_dst) = MAX( 1, MIN( nx_src, il_src( i_dst) ))
  
    END DO ! DO i_dst = 1, nx_dst
    
    il_src = il_src - 1
    il_src = MAX( 1, MIN( nx_src, il_src))
    
    ! iu
    DO i_dst = nx_dst, 1, -1
  
      IF (i_dst == nx_dst) THEN
        iu_src( i_dst) = nx_src
      ELSE
        iu_src( i_dst) = il_src( i_dst+1)
      END IF
  
      DO WHILE (x_src( iu_src( i_dst)) - dx_src / 2._dp > x_dst( i_dst) + dx_dst / 2._dp)
        iu_src( i_dst) = iu_src( i_dst) - 1
        IF (iu_src( i_dst) <= 1) THEN
          EXIT
        END IF
      END DO
      iu_src( i_dst) = MAX( 1, MIN( nx_src, iu_src( i_dst) ))
  
    END DO ! DO i_dst = nx_dst, 1, -1
    
    iu_src = iu_src + 1
    iu_src = MAX( 1, MIN( nx_src, iu_src))
    
    ! jl
    DO j_dst = 1, ny_dst
  
      IF (j_dst == 1) THEN
        jl_src( j_dst) = 1
      ELSE
        jl_src( j_dst) = jl_src( j_dst-1)
      END IF
  
      DO WHILE (y_src( jl_src( j_dst)) + dx_src / 2._dp < y_dst( j_dst) - dx_dst / 2._dp)
        jl_src( j_dst) = jl_src( j_dst) + 1
        IF (jl_src( j_dst) >= ny_src) THEN
          EXIT
        END IF
      END DO
      jl_src( j_dst) = MAX( 1, MIN( ny_src, jl_src( j_dst) ))
  
    END DO ! DO j_dst = 1, ny_dst
    
    jl_src = jl_src - 1
    jl_src = MAX( 1, MIN( ny_src, jl_src))
    
    ! ju
    DO j_dst = ny_dst, 1, -1
  
      IF (j_dst == ny_dst) THEN
        ju_src( j_dst) = ny_src
      ELSE
        ju_src( j_dst) = jl_src( j_dst+1)
      END IF
  
      DO WHILE (y_src( ju_src( j_dst)) - dx_src / 2._dp > y_dst( j_dst) + dx_dst / 2._dp)
        ju_src( j_dst) = ju_src( j_dst) - 1
        IF (ju_src( j_dst) <= 1) THEN
          EXIT
        END IF
      END DO
      ju_src( j_dst) = MAX( 1, MIN( ny_src, ju_src( j_dst) ))
  
    END DO ! DO j_dst = ny_dst, 1, -1
    
    ju_src = ju_src + 1
    ju_src = MAX( 1, MIN( ny_src, ju_src))
    
    ! Get derivatives of d_src
    CALL partition_list( nx_src, par%i, par%n, i1, i2)
    DO i = MAX(2,i1), MIN(nx_src-1,i2)
    DO j = 2, ny_src-1
      ddx_src( :,j,i) = (d_src( :,j,i+1) - d_src( :,j,i-1)) / (2._dp * dx_src)
      ddy_src( :,j,i) = (d_src( :,j+1,i) - d_src( :,j-1,i)) / (2._dp * dy_src)
    END DO
    END DO
    CALL sync
    
    ! Find parallelisation domains
    CALL partition_list( nx_dst, par%i, par%n, i1, i2)
    CALL partition_list( ny_dst, par%i, par%n, j1, j2)
    
    DO i = i1, i2
    DO j = 1, ny_dst
      
      d_dst( :,j,i) = 0._dp
      Asum          = 0._dp
      
      ! If this dst cell lies (partly) outside of the src grid, mark it as such;
      ! in that case, use nearest-neighbour extrapolation instead of conservative remapping
      IF (x_dst( i) - dx_dst/2._dp < MINVAL( x_src) - dx_src/2._dp .OR. &
          x_dst( i) + dx_dst/2._dp > MAXVAL( x_src) + dx_src/2._dp .OR. &
          y_dst( j) - dx_dst/2._dp < MINVAL( y_src) - dx_src/2._dp .OR. &
          y_dst( j) + dx_dst/2._dp > MAXVAL( y_src) + dx_src/2._dp) THEN
        mask_dst_outside_src( j,i) = 1
        CYCLE
      ELSE
        mask_dst_outside_src( j,i) = 0
      END IF
      
      DO i_src = il_src( i), iu_src( i)
      DO j_src = jl_src( j), ju_src( j)
        
        xomin = MAX( x_dst( i) - dx_dst/2._dp, x_src( i_src) - dx_src/2._dp)
        xomax = MIN( x_dst( i) + dx_dst/2._dp, x_src( i_src) + dx_src/2._dp)
        yomin = MAX( y_dst( j) - dy_dst/2._dp, y_src( j_src) - dy_src/2._dp)
        yomax = MIN( y_dst( j) + dy_dst/2._dp, y_src( j_src) + dy_src/2._dp)
        
        IF (xomax <= xomin .OR. yomax <= yomin) CYCLE
        
        Asd  = (xomax - xomin) * (yomax - yomin)
        Asum = Asum + Asd
        
        w0  = Asd / Ad
        w1x = 1._dp / Ad * (line_integral_mxydx( [xomin,yomin], [xomax,yomin], 1E-9_dp) + &
                            line_integral_mxydx( [xomax,yomin], [xomax,yomax], 1E-9_dp) + &
                            line_integral_mxydx( [xomax,yomax], [xomin,yomax], 1E-9_dp) + &
                            line_integral_mxydx( [xomin,yomax], [xomin,yomin], 1E-9_dp)) - w0 * x_src( i_src)
        w1y = 1._dp / Ad * (line_integral_xydy(  [xomin,yomin], [xomax,yomin], 1E-9_dp) + &
                            line_integral_xydy(  [xomax,yomin], [xomax,yomax], 1E-9_dp) + &
                            line_integral_xydy(  [xomax,yomax], [xomin,yomax], 1E-9_dp) + &
                            line_integral_xydy(  [xomin,yomax], [xomin,yomin], 1E-9_dp)) - w0 * y_src( j_src)
        
        d_dst( :,j,i) = d_dst( :,j,i) + w0  * d_src(   :,j_src,i_src) + &
                                        w1x * ddx_src( :,j_src,i_src) + &
                                        w1y * ddy_src( :,j_src,i_src)
        
      END DO ! DO j_src = jr_src( j,1), jr_src( j,2)
      END DO ! DO i_src = ir_src( i,1), ir_src( i,2)
      
      ! Safety
      IF (ABS( 1._dp - Asum / Ad) > 1E-4_dp) THEN
        CALL crash('dst grid cell [{int_01},{int_02}] couldnt be completely filled! Asum = {dp_01}, Ad = {dp_02}', int_01 = j, int_02 = i, dp_01 = Asum, dp_02 = Ad)
      END IF
      
    END DO ! DO j = 1, ny_dst
    END DO ! DO i = i1, i2
    CALL sync
    
    ! Use nearest-neighbour extrapolation for dst cells outside of the src grid
    ! =========================================================================
    
    ! Find the range of grid cells that were mapped correctly
    igmin = 0
    igmax = 0
    jgmin = 0
    jgmax = 0
    
    j = INT( REAL(ny_dst,dp)/2._dp)
    DO i = 1, nx_dst
      IF (mask_dst_outside_src( j,i) == 0) THEN
        igmin = i
        EXIT
      END IF
    END DO
    DO i = nx_dst, 1, -1
      IF (mask_dst_outside_src( j,i) == 0) THEN
        igmax = i
        EXIT
      END IF
    END DO
    
    i = INT( REAL(nx_dst,dp)/2._dp)
    DO j = 1, ny_dst
      IF (mask_dst_outside_src( j,i) == 0) THEN
        jgmin = j
        EXIT
      END IF
    END DO
    DO j = ny_dst, 1, -1
      IF (mask_dst_outside_src( j,i) == 0) THEN
        jgmax = j
        EXIT
      END IF
    END DO
    
    ! Corners
    IF (par%master) THEN
      DO k = 1, SIZE( d_src,1)
        ! Southwest
        d_dst( k,1      :jgmin-1 ,1      :igmin-1) = d_dst( k,jgmin,igmin)
        ! Southeast
        d_dst( k,1      :jgmin-1 ,igmax+1:nx_dst ) = d_dst( k,jgmin,igmax)
        ! Northwest
        d_dst( k,jgmax+1:ny_dst  ,1      :igmin-1) = d_dst( k,jgmax,igmin)
        ! Northeast
        d_dst( k,jgmax+1:ny_dst  ,igmax+1:nx_dst ) = d_dst( k,jgmax,igmax)
      END DO
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Borders
    DO k = 1, SIZE( d_src,1)
      DO i = MAX(i1,igmin), MIN(i2,igmax)
        ! South
        d_dst( k,1      :jgmin-1,i) = d_dst( k,jgmin,i)
        ! North
        d_dst( k,jgmax+1:ny_dst ,i) = d_dst( k,jgmax,i)
      END DO
      DO j = MAX(j1,jgmin), MIN(j2,jgmax)
        ! West
        d_dst( k,j,1      :igmin-1) = d_dst( k,j,igmin)
        ! East
        d_dst( k,j,igmax+1:nx_dst ) = d_dst( k,j,igmax)
      END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wddx_src             )
    CALL deallocate_shared( wddy_src             )
    CALL deallocate_shared( wmask_dst_outside_src)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE map_square_to_square_cons_2nd_order_3D
  FUNCTION line_integral_xdy(   p, q, tol_dist) RESULT( I_pq)
    ! Calculate the line integral x dy from p to q
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p, q
    REAL(dp),                                INTENT(IN)    :: tol_dist
    REAL(dp)                                               :: I_pq
    
    ! Local variables:
    REAL(dp)                                               :: xp, yp, xq, yq, dx, dy
        
    xp = p(1)
    yp = p(2)
    xq = q(1)
    yq = q(2)
    
    IF (ABS(yp-yq) < tol_dist) THEN
      I_pq = 0._dp
      RETURN
    END IF
    
    dx = q(1)-p(1)
    dy = q(2)-p(2)
    
    I_pq = xp*dy - yp*dx + (dx / (2._dp*dy)) * (yq**2 - yp**2)
    
  END FUNCTION line_integral_xdy
  FUNCTION line_integral_mxydx( p, q, tol_dist) RESULT( I_pq)
    ! Calculate the line integral -xy dx from p to q    
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p, q
    REAL(dp),                                INTENT(IN)    :: tol_dist
    REAL(dp)                                               :: I_pq
    
    ! Local variables:
    REAL(dp)                                               :: xp, yp, xq, yq, dx, dy
        
    xp = p(1)
    yp = p(2)
    xq = q(1)
    yq = q(2)
    
    IF (ABS(xp-xq) < tol_dist) THEN
      I_pq = 0._dp
      RETURN
    END IF
    
    dx = q(1)-p(1)
    dy = q(2)-p(2)
    
    I_pq = (1._dp/2._dp * (xp*dy/dx - yp) * (xq**2-xp**2)) - (1._dp/3._dp * dy/dx * (xq**3-xp**3))
    
  END FUNCTION line_integral_mxydx
  FUNCTION line_integral_xydy(  p, q, tol_dist) RESULT( I_pq)
    ! Calculate the line integral xy dy from p to q    
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p, q
    REAL(dp),                                INTENT(IN)    :: tol_dist
    REAL(dp)                                               :: I_pq
    
    ! Local variables:
    REAL(dp)                                               :: xp, yp, xq, yq, dx, dy
        
    xp = p(1)
    yp = p(2)
    xq = q(1)
    yq = q(2)
    
    IF (ABS(yp-yq) < tol_dist) THEN
      I_pq = 0._dp
      RETURN
    END IF
    
    dx = q(1)-p(1)
    dy = q(2)-p(2)
    
    I_pq = (1._dp/2._dp * (xp - yp*dx/dy) * (yq**2-yp**2)) + (1._dp/3._dp * dx/dy * (yq**3-yp**3))
    
  END FUNCTION line_integral_xydy

! == Map data from a global lon/lat-grid to the model y/x-grid
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
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'map_glob_to_grid_2D'
    INTEGER                                                :: i, j, il, iu, jl, ju
    REAL(dp)                                               :: wil, wiu, wjl, wju
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Find enveloping lat-lon indices
      il  = MAX(1,MIN(nlon-1, 1 + FLOOR((grid%lon(j,i)-MINVAL(lon)) / (lon(2)-lon(1)))))
      iu  = il+1        
      wil = MAX( 0._dp, MIN( 1._dp, (lon(iu) - grid%lon(j,i))/(lon(2)-lon(1)) ))
      wiu = 1._dp - wil
      
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
      wjl = MAX( 0._dp, MIN( 1._dp, (lat(ju) - grid%lat(j,i))/(lat(2)-lat(1)) ))
      wju = 1._dp - wjl
      
      ! Interpolate data
      d_grid( j,i) = (d_glob( il,jl) * wil * wjl) + &
                     (d_glob( il,ju) * wil * wju) + &
                     (d_glob( iu,jl) * wiu * wjl) + &
                     (d_glob( iu,ju) * wiu * wju)
    
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE map_glob_to_grid_2D
  SUBROUTINE map_glob_to_grid_3D( nlat, nlon, lat, lon, grid, d_glob, d_grid)
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
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'map_glob_to_grid_3D'
    INTEGER                                      :: i, j, il, iu, jl, ju, k, nz
    REAL(dp)                                     :: wil, wiu, wjl, wju
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    nz = SIZE( d_glob,3)
    
    ! Safety
    IF (SIZE(d_grid,1) /= nz) THEN
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Find enveloping lat-lon indices
      il  = MAX(1,MIN(nlon-1, 1 + FLOOR((grid%lon(j,i)-MINVAL(lon)) / (lon(2)-lon(1)))))
      iu  = il+1        
      wil = MAX( 0._dp, MIN( 1._dp, (lon(iu) - grid%lon(j,i))/(lon(2)-lon(1)) ))
      wiu = 1._dp - wil
      
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
      wjl = MAX( 0._dp, MIN( 1._dp, (lat(ju) - grid%lat(j,i))/(lat(2)-lat(1)) ))
      wju = 1._dp - wjl
      
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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE map_glob_to_grid_3D
  
! == Bilinear and bicubic interpolation
  FUNCTION interp_bilin_2D( d, x, y, xq, yq) RESULT( dq)
    ! Simple bilinear interpolation
    
    IMPLICIT NONE
    
    ! In- and output variables
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: x, y
    REAL(dp),                            INTENT(IN)    :: xq, yq
    REAL(dp)                                           :: dq
    
    ! Local variables:
    INTEGER                                            :: nx,ny
    REAL(dp)                                           :: dx
    INTEGER                                            :: il,iu,jl,ju
    REAL(dp)                                           :: wil,wiu,wjl,wju
    
    nx = SIZE(x)
    ny = SIZE(y)
    dx = x(2)-x(1)
    
    il  = MAX( 1, MIN( nx-1, 1 + FLOOR((xq-MINVAL(x)) / dx)))
    iu  = il + 1
    wil = (x(iu) - xq) / dx
    wiu = 1._dp - wil
    
    jl  = MAX( 1, MIN( ny-1, 1 + FLOOR((yq-MINVAL(y)) / dx)))
    ju  = jl + 1
    wjl = (y(ju) - yq) / dx
    wju = 1._dp - wjl
    
    ! interpolate
    dq = (wil * wjl * d( jl,il)) + &
         (wil * wju * d( ju,il)) + &
         (wiu * wjl * d( jl,iu)) + &
         (wiu * wju * d( ju,iu))
        
  END FUNCTION interp_bilin_2D
  
! == Check if a point lies inside a polygon (used in basin definition)
  FUNCTION is_in_polygon( Vpoly, p) RESULT( isso)
    ! Use the ray-casting algorithm to check if the point p = [px,py]
    ! lies inside the polygon spanned by poly = [x1,y1; x2,y2; ...; xn,yn]
    
    IMPLICIT NONE
    
    ! In- and output variables
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Vpoly
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    LOGICAL                                            :: isso
    
    ! Local variables:
    REAL(dp), DIMENSION(2)                             :: q,r,s,llis
    REAL(dp)                                           :: xmin,xmax,ymin,ymax
    INTEGER                                            :: n_intersects
    INTEGER                                            :: vi,vj,n_vertices
    LOGICAL                                            :: do_cross
    REAL(dp), PARAMETER                                :: tol_dist = 1E-5_dp
    
    isso = .FALSE.
    
    xmin = MINVAL( Vpoly(:,1))
    xmax = MAXVAL( Vpoly(:,1))
    ymin = MINVAL( Vpoly(:,2))
    ymax = MAXVAL( Vpoly(:,2))
    
    ! Quick test
    IF (p(1) < xmin .OR. p(1) > xmax .OR. &
        p(2) < ymin .OR. p(2) > ymax) THEN
      isso = .FALSE.
      RETURN
    END IF
    
    ! Define the endpoint of the east-pointing ray
    q = [xmax + (xmax - xmin) / 10._dp, p(2)]
    
    ! Determine how often the ray intersects the polygon
    
    n_vertices   = SIZE( Vpoly,1)
    n_intersects = 0
    
    DO vi = 1, n_vertices
    
      ! Find vertices spanning a polygon line section
      IF (vi < n_vertices) THEN
        vj = vi + 1
      ELSE
        vj = 1
      END IF
      
      ! Define the line section
      r = Vpoly( vi,:)
      s = Vpoly( vj,:)
      
      ! Determine if the ray intersects the line section
      IF ((r(2) < p(2) .AND. s(2) < p(2)) .OR. (r(2) > p(2) .AND. s(2) > p(2))) THEN
        do_cross = .FALSE.
      ELSE
        CALL segment_intersection( p, q, r, s, llis, do_cross, tol_dist)
      END IF
      
      IF (do_cross) n_intersects = n_intersects + 1
      
    END DO ! DO vi = 1, n_vertices
    
    ! If the number of intersections is odd, p lies inside the polygon
    IF (MOD( n_intersects,2) == 1) THEN
      isso = .TRUE.
    ELSE
      isso = .FALSE.
    END IF
  
  END FUNCTION is_in_polygon
  SUBROUTINE segment_intersection( p, q, r, s, llis, do_cross, tol_dist)
    ! Find out if the line segments [pq] and [rs] intersect. If so, return
    ! the coordinates of the point of intersection    
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(2),   INTENT(IN)          :: p, q, r, s
    REAL(dp), DIMENSION(2),   INTENT(OUT)         :: llis
    LOGICAL,                  INTENT(OUT)         :: do_cross
    REAL(dp),                 INTENT(IN)          :: tol_dist
    
    ! Local variables:
    REAL(dp), DIMENSION(2,2)                      :: A
    REAL(dp), DIMENSION(2)                        :: x, b
    INTEGER,  DIMENSION(2)                        :: IPIV
    INTEGER                                       :: info

    ! External subroutines:      
    EXTERNAL DGESV ! Lapack routine that solves matrix equation Ax=b for x (in double precision)
    
    ! If pq and rs are colinear, define them as not intersecting
    IF ( ABS( cross2( [q(1)-p(1), q(2)-p(2)], [s(1)-r(1), s(2)-r(2)] )) < tol_dist) THEN
      llis = [0._dp, 0._dp]
      do_cross = .FALSE.
      RETURN
    END IF
    
    A(1,:) = [(p(1)-q(1)), (r(1)-s(1))]
    A(2,:) = [(p(2)-q(2)), (r(2)-s(2))]
    b = [(r(1)-q(1)), (r(2)-q(2))]
    
    ! The LAPACK solver will overwrite the right-hand side b with the solution x. Therefore we 
    ! first copy the rhs in the solution vector x:
    x = b
    
    ! Solve Ax = b using LAPACK
    CALL DGESV( 2, 1, A, 2, IPIV, x, 2, info)
    
    llis = [q(1) + x(1) * (p(1)-q(1)), q(2) + x(1) * (p(2)-q(2))]
    
    IF (x(1)>0._dp .AND. x(1)<1._dp .AND. x(2)>0._dp .AND. x(2)<1._dp) THEN
      do_cross = .TRUE.
    ELSE
      do_cross = .FALSE.
    END IF
    
  END SUBROUTINE segment_intersection
  FUNCTION cross2( a,b) RESULT(z)
    ! Vector product z between 2-dimensional vectors a and b    
    
    IMPLICIT NONE
    
    REAL(dp), DIMENSION(2),     INTENT(IN)        :: a, b
    REAL(dp)                                      :: z

    z = (a(1)*b(2)) - (a(2)*b(1))

  END FUNCTION cross2
  
! == Transpose a data field (i.e. go from [i,j] to [j,i] indexing or the other way round)
  SUBROUTINE transpose_dp_2D( d, wd)
    ! Transpose a data field (i.e. go from [i,j] to [j,i] indexing or the other way round)
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: d
    INTEGER,                             INTENT(INOUT) :: wd
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'transpose_dp_2D'
    INTEGER                                      :: i,j,nx,ny,i1,i2
    REAL(dp), DIMENSION(:,:  ), POINTER          ::  d_temp
    INTEGER                                      :: wd_temp
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    nx = SIZE( d,1)
    ny = SIZE( d,2)
    CALL partition_list( nx, par%i, par%n, i1, i2)
    
    ! Allocate temporary memory
    CALL allocate_shared_dp_2D( nx, ny, d_temp, wd_temp)
    
    ! Copy data to temporary memory
    DO i = i1,i2
    DO j = 1, ny
      d_temp( i,j) = d( i,j)
    END DO
    END DO
    CALL sync
    
    ! Deallocate memory
    CALL deallocate_shared( wd)
    
    ! Reallocate transposed memory
    CALL allocate_shared_dp_2D( ny, nx, d, wd)
    
    ! Copy and transpose data from temporary memory
    DO i = i1, i2
    DO j = 1, ny
      d( j,i) = d_temp( i,j)
    END DO
    END DO
    CALL sync
    
    ! Deallocate temporary memory
    CALL deallocate_shared( wd_temp)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE transpose_dp_2D
  SUBROUTINE transpose_dp_3D( d, wd)
    ! Transpose a data field (i.e. go from [i,j] to [j,i] indexing or the other way round)
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: d
    INTEGER,                             INTENT(INOUT) :: wd
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'transpose_dp_3D'
    INTEGER                                      :: i,j,k,nx,ny,nz,i1,i2
    REAL(dp), DIMENSION(:,:,:), POINTER          ::  d_temp
    INTEGER                                      :: wd_temp
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    nx = SIZE( d,1)
    ny = SIZE( d,2)
    nz = SIZE( d,3)
    CALL partition_list( nx, par%i, par%n, i1, i2)
    
    ! Allocate temporary memory
    CALL allocate_shared_dp_3D( nx, ny, nz, d_temp, wd_temp)
    
    ! Copy data to temporary memory
    DO i = i1,i2
    DO j = 1, ny
    DO k = 1, nz
      d_temp( i,j,k) = d( i,j,k)
    END DO
    END DO
    END DO
    CALL sync
    
    ! Deallocate memory
    CALL deallocate_shared( wd)
    
    ! Reallocate transposed memory
    CALL allocate_shared_dp_3D( nz, ny, nx, d, wd)
    
    ! Copy and transpose data from temporary memory
    DO i = i1, i2
    DO j = 1, ny
    DO k = 1, nz
      d( k,j,i) = d_temp( i,j,k)
    END DO
    END DO
    END DO
    CALL sync
    
    ! Deallocate temporary memory
    CALL deallocate_shared( wd_temp)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE transpose_dp_3D
  SUBROUTINE transpose_int_2D( d, wd)
    ! Transpose a data field (i.e. go from [i,j] to [j,i] indexing or the other way round)
    
    IMPLICIT NONE
    
    ! In/output variables:
    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: d
    INTEGER,                             INTENT(INOUT) :: wd
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'transpose_int_2D'
    INTEGER                                      :: i,j,nx,ny,i1,i2
    INTEGER,  DIMENSION(:,:  ), POINTER          ::  d_temp
    INTEGER                                      :: wd_temp
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    nx = SIZE( d,1)
    ny = SIZE( d,2)
    CALL partition_list( nx, par%i, par%n, i1, i2)
    
    ! Allocate temporary memory
    CALL allocate_shared_int_2D( nx, ny, d_temp, wd_temp)
    
    ! Copy data to temporary memory
    DO i = i1,i2
    DO j = 1, ny
      d_temp( i,j) = d( i,j)
    END DO
    END DO
    CALL sync
    
    ! Deallocate memory
    CALL deallocate_shared( wd)
    
    ! Reallocate transposed memory
    CALL allocate_shared_int_2D( ny, nx, d, wd)
    
    ! Copy and transpose data from temporary memory
    DO i = i1, i2
    DO j = 1, ny
      d( j,i) = d_temp( i,j)
    END DO
    END DO
    CALL sync
    
    ! Deallocate temporary memory
    CALL deallocate_shared( wd_temp)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE transpose_int_2D
  SUBROUTINE transpose_int_3D( d, wd)
    ! Transpose a data field (i.e. go from [i,j] to [j,i] indexing or the other way round)
    
    IMPLICIT NONE
    
    ! In/output variables:
    INTEGER,  DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: d
    INTEGER,                             INTENT(INOUT) :: wd
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'transpose_int_3D'
    INTEGER                                      :: i,j,k,nx,ny,nz,i1,i2
    INTEGER,  DIMENSION(:,:,:), POINTER          ::  d_temp
    INTEGER                                      :: wd_temp
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    nx = SIZE( d,1)
    ny = SIZE( d,2)
    nz = SIZE( d,3)
    CALL partition_list( nx, par%i, par%n, i1, i2)
    
    ! Allocate temporary memory
    CALL allocate_shared_int_3D( nx, ny, nz, d_temp, wd_temp)
    
    ! Copy data to temporary memory
    DO i = i1,i2
    DO j = 1, ny
    DO k = 1, nz
      d_temp( i,j,k) = d( i,j,k)
    END DO
    END DO
    END DO
    CALL sync
    
    ! Deallocate memory
    CALL deallocate_shared( wd)
    
    ! Reallocate transposed memory
    CALL allocate_shared_int_3D( nz, ny, nx, d, wd)
    
    ! Copy and transpose data from temporary memory
    DO i = i1, i2
    DO j = 1, ny
    DO k = 1, nz
      d( k,j,i) = d_temp( i,j,k)
    END DO
    END DO
    END DO
    CALL sync
    
    ! Deallocate temporary memory
    CALL deallocate_shared( wd_temp)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE transpose_int_3D
  
! == Interpolate ocean column data to a queried depth
  SUBROUTINE interpolate_ocean_depth( nz_ocean, z_ocean, f_ocean, z_query, f_query)
    ! Interpolate ocean column data to a queried depth using a simple bisection method.
    
    IMPLICIT NONE
    
    ! In/output variables:
    INTEGER,                             INTENT(IN)    :: nz_ocean    ! Number of vertical layers
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: z_ocean     ! Depth of layers (assumed to be monotonically increasing, does not need to be regularly spaced)
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: f_ocean     ! Value of whatever function we want to interpolate
    REAL(dp),                            INTENT(IN)    :: z_query     ! Depth at which we want to know the function
    REAL(dp),                            INTENT(OUT)   :: f_query     ! Interpolated function value
    
    ! Local variables:
    INTEGER                                            :: k_lo,k_hi,k_mid
    LOGICAL                                            :: foundit
    REAL(dp)                                           :: w
    
    ! Safety
    IF (z_query < 0._dp) THEN
      CALL crash('interpolate_ocean_depth - ERROR: z_query = {dp_01} < 0; cannot extrapolate above the sea surface, obviously!', dp_01 = z_query)
    ELSEIF (z_query > 12000._dp) THEN
      CALL crash('interpolate_ocean_depth - ERROR: z_query = {dp_01} > 12 km;  the ocean is not that deep!', dp_01 = z_query)
    ELSEIF (SIZE(z_ocean,1) /= nz_ocean) THEN
      CALL crash('interpolate_ocean_depth - ERROR: SIZE( z_ocean,1) = {int_01} /= nz_ocean = {int_02}!', int_01 = SIZE( z_ocean,1), int_02 = nz_ocean)
    ELSEIF (SIZE(f_ocean,1) /= nz_ocean) THEN
      CALL crash('interpolate_ocean_depth - ERROR: SIZE( f_ocean,1) = {int_01} /= nz_ocean = {int_02}!', int_01 = SIZE( f_ocean,1), int_02 = nz_ocean)
    ELSEIF (z_query > MAXVAL( z_ocean)) THEN
      ! Nearest-neighbour extrapolation when querying data beneath the end of the ocean data column
      f_query = f_ocean( nz_ocean)
      RETURN
    END IF
    
    ! Exception for when z_query = 0 (the World Ocean Atlas depth starts at 1.25...)
    IF (z_query < MINVAL( z_ocean)) THEN
      f_query = f_ocean(1)
      RETURN
    END IF
    
    ! Bisection method
    k_lo  = 1
    k_hi  = nz_ocean
    k_mid = INT( REAL(k_lo + k_hi,dp) / 2._dp)
    
    ! Exceptions
    IF     (ABS( z_query - z_ocean( k_lo )) < 1E-4_dp) THEN
      f_query = f_ocean( k_lo)
      RETURN
    ELSEIF (ABS( z_query - z_ocean( k_hi )) < 1E-4_dp) THEN
      f_query = f_ocean( k_hi)
      RETURN
    ELSEIF (ABS( z_query - z_ocean( k_mid)) < 1E-4_dp) THEN
      f_query = f_ocean( k_mid)
      RETURN
    END IF
    
    ! Bisection method
    foundit = .FALSE.
    DO WHILE (.NOT. foundit)
    
      IF (ABS( z_query - z_ocean( k_mid)) < 1E-4_dp) THEN
        ! Exception for when the queried depth is exactly at the midpoint index depth
        f_query = f_ocean( k_mid)
        RETURN
      ELSEIF ( z_query > z_ocean( k_mid)) THEN
        ! Queried depth lies to the right of the midpoint
        k_lo = k_mid
        k_mid = INT( REAL(k_lo + k_hi,dp) / 2._dp)
      ELSE
        ! Queried depth lies to the left of the midpoint
        k_hi = k_mid
        k_mid = INT( REAL(k_lo + k_hi,dp) / 2._dp)
      END IF
      
      ! Stop iterating when endpoints lie next to each other; then just do linear interpolation between those two.
      IF (k_hi == k_lo+1) foundit = .TRUE.
      
    END DO ! DO WHILE (.NOT. foundit)
    
    ! Linear interpolation between nearest layers
    w = (z_query - z_ocean( k_lo)) / (z_ocean( k_hi) - z_ocean( k_lo))
    f_query = w * f_ocean( k_hi) + (1._dp - w) * f_ocean( k_lo)
    
  END SUBROUTINE interpolate_ocean_depth

! == 2nd-order conservative remapping of a 1-D variable
  SUBROUTINE remap_cons_2nd_order_1D( z_src, mask_src, d_src, z_dst, mask_dst, d_dst)
    ! 2nd-order conservative remapping of a 1-D variable
    ! 
    ! Used to remap ocean data from the provided vertical grid to the IMAU-ICE ocean vertical grid
    !
    ! Both z_src and z_dst can be irregular.
    ! 
    ! Both the src and dst data have a mask, with 0 indicating grid points where no data is defined.
    !
    ! This subroutine is serial, as it will be applied to single grid cells when remapping 3-D data fields,
    !   with the parallelisation being done by distributing the 2-D grid cells over the processes.
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: z_src
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_src
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_src
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: z_dst
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_dst
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_dst
    
    ! Local variables:
    LOGICAL                                            :: all_are_masked
    INTEGER                                            :: nz_src, nz_dst
    INTEGER                                            :: k
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: ddz_src
    INTEGER                                            :: k_src, k_dst
    REAL(dp)                                           :: zl_src, zu_src, zl_dst, zu_dst, z_lo, z_hi, z, d
    REAL(dp)                                           :: dz_overlap, dz_overlap_tot, d_int, d_int_tot
    REAL(dp)                                           :: dist_to_dst, dist_to_dst_min, max_dist
    INTEGER                                            :: k_src_nearest_to_dst
    REAL(dp)                                           :: z_src_lowest, z_src_highest, NaN
    
    NaN = -1._dp
    NaN = SQRT( NaN)
    
    ! Initialise
    d_dst = 0._dp
    
    ! Sizes
    nz_src = SIZE( z_src,1)
    nz_dst = SIZE( z_dst,1)
    
    ! Maximum distance on combined grids
    max_dist = MAXVAL([ ABS( z_src( nz_src) - z_src( 1)), &
                        ABS( z_dst( nz_dst) - z_dst( 1)), &
                        ABS( z_src( nz_src) - z_dst( 1)), &
                        ABS( z_dst( nz_dst) - z_src( 1))])
    
    ! Exception for when the entire src field is masked
    all_are_masked = .TRUE.
    DO k = 1, nz_src
      IF (mask_src( k) == 1) all_are_masked = .FALSE.
    END DO
    IF (all_are_masked) RETURN
    
    ! Exception for when the entire dst field is masked
    all_are_masked = .TRUE.
    DO k = 1, nz_dst
      IF (mask_dst( k) == 1) all_are_masked = .FALSE.
    END DO
    IF (all_are_masked) RETURN
    
    ! Calculate derivative d_src/dz (one-sided differencing at the boundary or in cells with masked neighbours, central differencing everywhere else)
    ALLOCATE( ddz_src( nz_src))
    DO k = 2, nz_src-1
      IF (mask_src( k-1) == 1 .AND. mask_src( k) == 1 .AND. mask_src( k+1) == 1) THEN
        ddz_src( k    ) = (d_src( k+1   ) - d_src( k-1     )) / (z_src( k+1   ) - z_src( k-1     ))
      ELSEIF (mask_src( k) == 0) THEN
        ddz_src( k    ) = 0._dp
      ELSEIF (mask_src( k-1) == 0) THEN
        ddz_src( k    ) = (d_src( k+1   ) - d_src( k       )) / (z_src( k+1   ) - z_src( k       ))
      ELSEIF (mask_src( k+1) == 0) THEN
        ddz_src( k    ) = (d_src( k     ) - d_src( k-1     )) / (z_src( k     ) - z_src( k-1     ))
      ELSE
        CALL crash('whaa!')
      END IF
    END DO
    ddz_src(    1     ) = (d_src( 2     ) - d_src( 1       )) / (z_src( 2     ) - z_src( 1       ))
    ddz_src(    nz_src) = (d_src( nz_src) - d_src( nz_src-1)) / (z_src( nz_src) - z_src( nz_src-1))
    
    ! Find lowest and highest z value with source data
    k_src = 1
    DO WHILE (mask_src( k_src) == 0)
      k_src = k_src + 1
    END DO
    z_src_lowest = z_src( k_src)
    
    k_src = nz_src
    DO WHILE (mask_src( k_src) == 0)
      k_src = k_src - 1
    END DO
    z_src_highest = z_src( k_src)
    
    ! Perform conservative remapping by finding regions of overlap
    ! between source and destination grid cells
    
    DO k_dst = 1, nz_dst
      
      ! Skip masked grid cells
      IF (mask_dst( k_dst) == 0) THEN
        d_dst( k_dst) = 0._dp
        CYCLE
      END IF
      
      ! Find z range covered by this dst grid cell
      IF (k_dst > 1) THEN
        zl_dst = 0.5_dp * (z_dst( k_dst - 1) + z_dst( k_dst))
      ELSE
        zl_dst = z_dst( 1) - 0.5_dp * (z_dst( 2) - z_dst( 1))
      END IF
      IF (k_dst < nz_dst) THEN
        zu_dst = 0.5_dp * (z_dst( k_dst + 1) + z_dst( k_dst))
      ELSE
        zu_dst = z_dst( nz_dst) + 0.5_dp * (z_dst( nz_dst) - z_dst( nz_dst-1))
      END IF
      
      ! Find all overlapping src grid cells
      d_int_tot      = 0._dp
      dz_overlap_tot = 0._dp
      DO k_src = 1, nz_src
      
        ! Skip masked grid cells
        IF (mask_src( k_src) == 0) CYCLE
      
        ! Find z range covered by this src grid cell
        IF (k_src > 1) THEN
          zl_src = 0.5_dp * (z_src( k_src - 1) + z_src( k_src))
        ELSE
          zl_src = z_src( 1) - 0.5_dp * (z_src( 2) - z_src( 1))
        END IF
        IF (k_src < nz_src) THEN
          zu_src = 0.5_dp * (z_src( k_src + 1) + z_src( k_src))
        ELSE
          zu_src = z_src( nz_src) + 0.5_dp * (z_src( nz_src) - z_src( nz_src-1))
        END IF
        
        ! Find region of overlap
        z_lo = MAX( zl_src, zl_dst)
        z_hi = MIN( zu_src, zu_dst)
        dz_overlap = MAX( 0._dp, z_hi - z_lo)
        
        ! Calculate integral over region of overlap and add to sum
        IF (dz_overlap > 0._dp) THEN
          z = 0.5_dp * (z_lo + z_hi)
          d = d_src( k_src) + ddz_src( k_src) * (z - z_src( k_src))
          d_int = d * dz_overlap
          
          d_int_tot      = d_int_tot      + d_int
          dz_overlap_tot = dz_overlap_tot + dz_overlap
        END IF
        
      END DO ! DO k_src = 1, nz_src
      
      IF (dz_overlap_tot > 0._dp) THEN
        ! Calculate dst value
        d_dst( k_dst) = d_int_tot / dz_overlap_tot
      ELSE
        ! Exception for when no overlapping src grid cells were found; use nearest-neighbour extrapolation
        
        k_src_nearest_to_dst = 0._dp
        dist_to_dst_min      = max_dist
        DO k_src = 1, nz_src
          IF (mask_src( k_src) == 1) THEN
            dist_to_dst = ABS( z_src( k_src) - z_dst( k_dst))
            IF (dist_to_dst < dist_to_dst_min) THEN
              dist_to_dst_min      = dist_to_dst
              k_src_nearest_to_dst = k_src
            END IF
          END IF
        END DO
        
        ! Safety
        IF (k_src_nearest_to_dst == 0) THEN
          CALL crash('remap_cons_2nd_order_1D: couldnt find nearest neighbour on source grid!')
        END IF
        
        d_dst( k_dst) = d_src( k_src_nearest_to_dst)
        
      END IF ! IF (dz_overlap_tot > 0._dp) THEN
      
      ! Forbid extrapolation
      IF     (z_dst( k_dst) < z_src_lowest) THEN
        d_dst( k_dst) = NaN
      ELSEIF (z_dst( k_dst) > z_src_highest) THEN
        d_dst( k_dst) = NaN
      END IF
        
    END DO ! DO k_dst = 1, nz_dst
    
    ! Clean up after yourself
    DEALLOCATE( ddz_src)
    
  END SUBROUTINE remap_cons_2nd_order_1D
  
! == Remove Lake Vostok from Antarctic input geometry data
  SUBROUTINE remove_Lake_Vostok( x, y, Hi, Hb, Hs)
    ! Remove Lake Vostok from Antarctic input geometry data
    ! by manually increasing ice thickness so that Hi = Hs - Hb
    !
    ! NOTE: since IMAU-ICE doesn't consider subglacial lakes, Vostok simply shows
    !       up as a "dip" in the initial geometry. The model will run fine, the dip
    !       fills up in a few centuries, but it slows down the model for a while and
    !       it looks ugly, so we just remove it right away.
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: x,y
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: Hi
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Hb
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Hs
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'remove_Lake_Vostok'
    INTEGER                                       :: i,j,nx,ny
    REAL(dp), PARAMETER                           :: lake_Vostok_xmin = 1164250.0
    REAL(dp), PARAMETER                           :: lake_Vostok_xmax = 1514250.0
    REAL(dp), PARAMETER                           :: lake_Vostok_ymin = -470750.0
    REAL(dp), PARAMETER                           :: lake_Vostok_ymax = -220750.0
    INTEGER                                       :: il,iu,jl,ju
    
    IF (par%master) THEN
      
      nx = SIZE( Hi,2)
      ny = SIZE( Hi,1)
      
      il = 1
      DO WHILE (x( il) < lake_Vostok_xmin)
        il = il+1
      END DO
      iu = nx
      DO WHILE (x( iu) > lake_Vostok_xmax)
        iu = iu-1
      END DO
      jl = 1
      DO WHILE (y( jl) < lake_Vostok_ymin)
        jl = jl+1
      END DO
      ju = ny
      DO WHILE (y( ju) > lake_Vostok_ymax)
        ju = ju-1
      END DO
        
      DO i = il, iu
      DO j = jl, ju
        Hi( j,i) = Hs( j,i) - Hb( j,i)
      END DO
      END DO
      
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE remove_Lake_Vostok
  
! == Gaussian extrapolation (used for e.g. ocean data and basal roughness inversion)
  SUBROUTINE extrapolate_Gaussian_floodfill( grid, mask, d, sigma, mask_filled)
    ! Extrapolate the data field d into the area designated by the mask,
    ! using Gaussian extrapolation of sigma
    ! 
    ! Note about the mask:
    !    2 = data provided
    !    1 = no data provided, fill allowed
    !    0 = no fill allowed
    ! (so basically this routine extrapolates data from the area
    !  where mask == 2 into the area where mask == 1)
    ! 
    ! NOTE: Not parallelised! When using this routine to extrapolate ocean data,
    !       parallelisation can be achieved by dividing the vertical layers over
    !       the processes. Otherwise, be careful to only let the Master process
    !       call this routine!
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d
    REAL(dp),                            INTENT(IN)    :: sigma
    INTEGER,  DIMENSION(:,:  ),          INTENT(OUT)   :: mask_filled   ! 1 = successfully filled, 2 = failed to fill (region of to-be-filled pixels not connected to source data)
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extrapolate_Gaussian_floodfill'
    INTEGER                                            :: i,j,k,ii,jj,it
    INTEGER                                            :: stackN1, stackN2
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE            :: stack1, stack2
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE            :: map
    INTEGER                                            :: n_search
    LOGICAL                                            :: has_filled_neighbours
    INTEGER                                            :: n
    REAL(dp)                                           :: sum_d, w, sum_w
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    n_search = 1 + CEILING( 2._dp * sigma / grid%dx)
    
    ! Allocate map and stacks
    ALLOCATE( map(    grid%ny,  grid%nx   ))
    ALLOCATE( stack1( grid%ny * grid%nx, 2))
    ALLOCATE( stack2( grid%ny * grid%nx, 2))
      
    map         = 0
    stack1      = 0
    stack2      = 0
    stackN1     = 0
    stackN2     = 0
    mask_filled = 0
  
    ! Initialise the map from the mask
    DO i = 1, grid%nx
    DO j = 1, grid%ny
      IF (mask( j,i) == 2) THEN
        map( j,i) = 2
      END IF
    END DO
    END DO
    
    ! Initialise the stack with all empty-next-to-filled grid cells
    DO i = 1, grid%nx
    DO j = 1, grid%ny
        
      IF (mask( j,i) == 1) THEN
        ! This grid cell is empty and should be filled
        
        has_filled_neighbours = .FALSE.
        DO ii = MAX(1 ,i-1), MIN(grid%nx,i+1)
        DO jj = MAX(1 ,j-1), MIN(grid%ny,j+1)
          IF (mask( jj,ii) == 2) THEN
            has_filled_neighbours = .TRUE.
            EXIT
          END IF
        END DO
        IF (has_filled_neighbours) EXIT
        END DO
        
        IF (has_filled_neighbours) THEN
          ! Add this empty-with-filled-neighbours grid cell to the stack,
          ! and mark it as stacked on the map
          map( j,i) = 1
          stackN2 = stackN2 + 1
          stack2( stackN2,:) = [i,j]
        END IF
      
      END IF ! IF (map( i,j) == 0) THEN
        
    END DO
    END DO
    
    ! Perform the flood-fill
    it = 0
    DO WHILE (stackN2 > 0)
      
      it = it + 1
      
      ! Go over all the stacked empty-next-to-filled grid cells, perform the
      ! Gaussian-kernel extrapolation to fill, and mark them as as such on the map
      DO k = 1, stackN2
        
        ! Get grid cell indices
        i = stack2( k,1)
        j = stack2( k,2)
        
        ! Find Gaussian-weighted average value over nearby filled pixels within the basin
        n     = 0
        sum_d = 0._dp
        sum_w = 0._dp
        
        DO ii = MAX( 1 ,i - n_search), MIN( grid%nx,i + n_search)
        DO jj = MAX( 1 ,j - n_search), MIN( grid%ny,j + n_search)
        
          IF (map( jj,ii) == 2) THEN
            n     = n + 1
            w     = EXP( -0.5_dp * (SQRT(REAL(ii-i,dp)**2 + REAL(jj-j,dp)**2) / sigma)**2)
            sum_w = sum_w + w
            sum_d = sum_d + w * d( jj,ii)
          END IF
          
        END DO
        END DO
        
        ! Fill in averaged value
        d( j,i) = sum_d / sum_w
        
      END DO ! DO k = 1, stackN2
      
      ! Mark all newly-filled grid cells as such
      DO k = 1, stackN2
        
        ! Get grid cell indices
        i = stack2( k,1)
        j = stack2( k,2)
        
        ! Mark grid cell as filled
        map( j,i) = 2
        mask_filled( j,i) = 1
        
      END DO ! DO k = 1, stackN2
      
      ! Cycle the stacks
      stack1  = stack2
      stackN1 = stackN2
      stack2  = 0
      stackN2 = 0
      
      ! List new empty-next-to-filled grid cells
      DO k = 1, stackN1
        
        ! Get grid cell indices
        i = stack1( k,1)
        j = stack1( k,2)
        
        ! Find empty neighbours; if unlisted, list them and mark them on the map
        DO ii = MAX( 1, i-1), MIN( grid%nx, i+1)
        DO jj = MAX( 1, j-1), MIN( grid%ny, j+1)
        
          IF (map( jj,ii) == 0 .AND. mask( jj,ii) == 1) THEN
            map( jj,ii) = 1
            stackN2 = stackN2 + 1
            stack2( stackN2,:) = [ii,jj]
          END IF
          
        END DO
        END DO
        
      END DO ! DO k = 1, stackN1
      
      ! Safety
      IF (it > 2 * MAX( grid%ny, grid%nx)) THEN
        CALL crash('flood-fill got stuck!')
      END IF
              
    END DO ! DO WHILE (stackN2 > 0)
    
    ! Mark grid cells that could not be filled
    DO i = 1, grid%nx
    DO j = 1, grid%ny
      IF (mask_filled( j,i) == 0 .AND. mask( j,i) == 1) THEN
        mask_filled( j,i) = 2
      END IF
    END DO
    END DO
    
    ! Clean up after yourself
    DEALLOCATE( map   )
    DEALLOCATE( stack1)
    DEALLOCATE( stack2)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE extrapolate_Gaussian_floodfill
  
! == Debugging
  SUBROUTINE check_for_NaN_dp_1D(  d, d_name)
    ! Check if NaN values occur in the 1-D dp data field d
    ! NOTE: parallelised!
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:    ),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name
    
    ! Local variables:
    INTEGER                                                :: nx,i,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc
    
    ! Only do this when so specified in the config
    IF (.NOT. C%do_check_for_NaN) RETURN
    
    ! Field size
    nx = SIZE(d,1)
    
    ! Parallelisation range
    CALL partition_list( nx, par%i, par%n, i1, i2)
    
    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF
    
    ! Inspect data field
    DO i = i1, i2
    
      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...
      
      IF     (d( i) /= d( i)) THEN
        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01}]' ), int_01 = i)
      ELSEIF (d( i) > HUGE( d( i))) THEN
        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]' ), int_01 = i)
      ELSEIF (d( i)  < -HUGE( d( i))) THEN
        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]'), int_01 = i)
      END IF
      
    END DO
    CALL sync
    
  END SUBROUTINE check_for_NaN_dp_1D
  SUBROUTINE check_for_NaN_dp_2D(  d, d_name)
    ! Check if NaN values occur in the 2-D dp data field d
    ! NOTE: parallelised!
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name
    
    ! Local variables:
    INTEGER                                                :: nx,ny,i,j,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc
    
    ! Only do this when so specified in the config
    IF (.NOT. C%do_check_for_NaN) RETURN
    
    ! Field size
    nx = SIZE(d,2)
    ny = SIZE(d,1)
    
    ! Parallelisation range
    CALL partition_list( nx, par%i, par%n, i1, i2)
    
    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF
    
    ! Inspect data field
    DO i = i1, i2
    DO j = 1, ny
    
      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...
      
      IF     (d( j,i) /= d( j,i)) THEN
        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
      ELSEIF (d( j,i) > HUGE( d( j,i))) THEN
        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
      ELSEIF (d( j,i) < -HUGE( d( j,i))) THEN
        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
      END IF
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE check_for_NaN_dp_2D
  SUBROUTINE check_for_NaN_dp_3D(  d, d_name)
    ! Check if NaN values occur in the 3-D dp data field d
    ! NOTE: parallelised!
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name
    
    ! Local variables:
    INTEGER                                                :: nx,ny,nz,i,j,k,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc
    
    ! Only do this when so specified in the config
    IF (.NOT. C%do_check_for_NaN) RETURN
    
    ! Field size
    nx = SIZE(d,3)
    ny = SIZE(d,2)
    nz = SIZE(d,1)
    
    ! Parallelisation range
    CALL partition_list( nx, par%i, par%n, i1, i2)
    
    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF
    
    ! Inspect data field
    DO i = i1, i2
    DO j = 1, ny
    DO k = 1, nz
    
      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...
      
      IF     (d( k,j,i) /= d( k,j,i)) THEN
        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (d( k,j,i) > HUGE( d( k,j,i))) THEN
        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (d( k,j,i) < -HUGE( d( k,j,i))) THEN
        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
      END IF
      
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE check_for_NaN_dp_3D
  SUBROUTINE check_for_NaN_int_1D( d, d_name)
    ! Check if NaN values occur in the 1-D int data field d
    ! NOTE: parallelised!
    
    IMPLICIT NONE
    
    ! In/output variables:
    INTEGER,  DIMENSION(:    ),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name
    
    ! Local variables:
    INTEGER                                                :: nx,i,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc
    
    ! Only do this when so specified in the config
    IF (.NOT. C%do_check_for_NaN) RETURN
    
    ! Field size
    nx = SIZE(d,1)
    
    ! Parallelisation range
    CALL partition_list( nx, par%i, par%n, i1, i2)
    
    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF
    
    ! Inspect data field
    DO i = i1, i2
    
      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...
      
      IF     (d( i) /= d( i)) THEN
        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01}]' ), int_01 = i)
      ELSEIF (d( i) > HUGE( d( i))) THEN
        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]' ), int_01 = i)
      ELSEIF (d( i) < -HUGE( d( i))) THEN
        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]'), int_01 = i)
      END IF
      
    END DO
    CALL sync
    
  END SUBROUTINE check_for_NaN_int_1D
  SUBROUTINE check_for_NaN_int_2D( d, d_name)
    ! Check if NaN values occur in the 2-D int data field d
    ! NOTE: parallelised!
    
    IMPLICIT NONE
    
    ! In/output variables:
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name
    
    ! Local variables:
    INTEGER                                                :: nx,ny,i,j,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc
    
    ! Only do this when so specified in the config
    IF (.NOT. C%do_check_for_NaN) RETURN
    
    ! Field size
    nx = SIZE(d,2)
    ny = SIZE(d,1)
    
    ! Parallelisation range
    CALL partition_list( nx, par%i, par%n, i1, i2)
    
    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF
    
    ! Inspect data field
    DO i = i1, i2
    DO j = 1, ny
    
      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...
      
      IF     (d( j,i) /= d( j,i)) THEN
        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
      ELSEIF (d( j,i) > HUGE( d( j,i))) THEN
        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
      ELSEIF (d( j,i) < -HUGE( d( j,i))) THEN
        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
      END IF
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE check_for_NaN_int_2D
  SUBROUTINE check_for_NaN_int_3D( d, d_name)
    ! Check if NaN values occur in the 3-D int data field d
    ! NOTE: parallelised!
    
    IMPLICIT NONE
    
    ! In/output variables:
    INTEGER,  DIMENSION(:,:,:),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name
    
    ! Local variables:
    INTEGER                                                :: nx,ny,nz,i,j,k,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc
    
    ! Only do this when so specified in the config
    IF (.NOT. C%do_check_for_NaN) RETURN
    
    ! Field size
    nx = SIZE(d,3)
    ny = SIZE(d,2)
    nz = SIZE(d,1)
    
    ! Parallelisation range
    CALL partition_list( nx, par%i, par%n, i1, i2)
    
    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF
    
    ! Inspect data field
    DO i = i1, i2
    DO j = 1, ny
    DO k = 1, nz
    
      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...
      
      IF     (d( k,j,i) /= d( k,j,i)) THEN
        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (d( k,j,i) > HUGE( d( k,j,i))) THEN
        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (d( k,j,i) < -HUGE( d( k,j,i))) THEN
        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
      END IF
      
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE check_for_NaN_int_3D

END MODULE utilities_module
