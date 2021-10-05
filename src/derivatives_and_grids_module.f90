MODULE derivatives_and_grids_module

  ! Contains all the routines involved in the Arakawa grids, for mapping and
  ! calculating derivatives between all the different grids.

  USE mpi
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, partition_list
  USE configuration_module,            ONLY: dp, C
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_zeta_coefficients
  USE netcdf_module,                   ONLY: debug
  
  ! The vertical scaled coordinate zeta transformation coefficients
  TYPE(type_zeta_coefficients) :: zeta

CONTAINS

! ======================
! ==== Derivatives =====
! ======================
  
  ! Aa to Aa
  
  ! 2D
  SUBROUTINE ddx_a_to_a_2D(  grid, d_a, dx_a)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dx_a
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      dx_a( j,i) = (d_a( j,i+1) - d_a( j,i-1)) / (2 * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    dx_a( grid%j1:grid%j2,1      ) = (d_a( grid%j1:grid%j2,2      ) - d_a( grid%j1:grid%j2,1        )) / grid%dx
    dx_a( grid%j1:grid%j2,grid%nx) = (d_a( grid%j1:grid%j2,grid%nx) - d_a( grid%j1:grid%j2,grid%nx-1)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddx_a_to_a_2D
  SUBROUTINE ddy_a_to_a_2D(  grid, d_a, dy_a)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dy_a
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
      dy_a( j,i) = (d_a( j+1,i) - d_a( j-1,i)) / (2 * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    dy_a( 1      ,grid%i1:grid%i2) = (d_a( 2      ,grid%i1:grid%i2) - d_a( 1        ,grid%i1:grid%i2)) / grid%dx
    dy_a( grid%ny,grid%i1:grid%i2) = (d_a( grid%ny,grid%i1:grid%i2) - d_a( grid%ny-1,grid%i1:grid%i2)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddy_a_to_a_2D
  SUBROUTINE ddxx_a_to_a_2D( grid, d_a, dxx_a)
    ! Input:  scalar on the Aa grid
    ! Output: its xx-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dxx_a
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      dxx_a( j,i) = (d_a( j,i+1) + d_a( j,i-1) - 2._dp * d_a( j,i)) / grid%dx**2
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    dxx_a( grid%j1:grid%j2,1      ) = (d_a( grid%j1:grid%j2,3      ) + d_a( grid%j1:grid%j2,1        ) - 2._dp * d_a( grid%j1:grid%j2,2        )) / grid%dx**2
    dxx_a( grid%j1:grid%j2,grid%nx) = (d_a( grid%j1:grid%j2,grid%nx) + d_a( grid%j1:grid%j2,grid%nx-2) - 2._dp * d_a( grid%j1:grid%j2,grid%nx-1)) / grid%dx**2
    CALL sync
    
  END SUBROUTINE ddxx_a_to_a_2D
  SUBROUTINE ddyy_a_to_a_2D( grid, d_a, dyy_a)
    ! Input:  scalar on the Aa grid
    ! Output: its yy-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dyy_a
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
      dyy_a( j,i) = (d_a( j+1,i) + d_a( j-1,i) - 2._dp * d_a( j,i)) / grid%dx**2
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    dyy_a( 1      ,grid%i1:grid%i2) = (d_a( 3      ,grid%i1:grid%i2) + d_a( 1        ,grid%i1:grid%i2) - 2._dp * d_a( 2        ,grid%i1:grid%i2)) / grid%dx**2
    dyy_a( grid%ny,grid%i1:grid%i2) = (d_a( grid%ny,grid%i1:grid%i2) + d_a( grid%ny-2,grid%i1:grid%i2) - 2._dp * d_a( grid%ny-1,grid%i1:grid%i2)) / grid%dx**2
    CALL sync
    
  END SUBROUTINE ddyy_a_to_a_2D
  SUBROUTINE ddxy_a_to_a_2D( grid, d_a, dxy_a)
    ! Input:  scalar on the Aa grid
    ! Output: its xy-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dxy_a
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dxy_a( j,i) = (d_a( j+1,i+1) + d_a( j-1,i-1) - d_a( j+1,i-1) - d_a( j-1,i+1)) / (4._dp * grid%dx * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    ! NO IDEA HOW TO DO THIS...
    dxy_a( 1              ,grid%i1:grid%i2) = 0._dp
    dxy_a( grid%ny        ,grid%i1:grid%i2) = 0._dp
    dxy_a( grid%j1:grid%j2,1              ) = 0._dp
    dxy_a( grid%j1:grid%j2,grid%nx        ) = 0._dp
    CALL sync
    
  END SUBROUTINE ddxy_a_to_a_2D
  ! 3D
  SUBROUTINE ddx_a_to_a_3D(  grid, d_a, dx_a)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dx_a
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      dx_a( k,j,i) = (d_a( k,j,i+1) - d_a( k,j,i-1)) / (2 * grid%dx)
    END DO
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    dx_a( :,grid%j1:grid%j2,1      ) = (d_a( :,grid%j1:grid%j2,2      ) - d_a( :,grid%j1:grid%j2,1        )) / grid%dx
    dx_a( :,grid%j1:grid%j2,grid%nx) = (d_a( :,grid%j1:grid%j2,grid%nx) - d_a( :,grid%j1:grid%j2,grid%nx-1)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddx_a_to_a_3D
  SUBROUTINE ddy_a_to_a_3D(  grid, d_a, dy_a)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dy_a
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      dy_a( k,j,i) = (d_a( k,j+1,i) - d_a( k,j-1,i)) / (2 * grid%dx)
    END DO
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    dy_a( :,1      ,grid%i1:grid%i2) = (d_a( :,2      ,grid%i1:grid%i2) - d_a( :,1        ,grid%i1:grid%i2)) / grid%dx
    dy_a( :,grid%ny,grid%i1:grid%i2) = (d_a( :,grid%ny,grid%i1:grid%i2) - d_a( :,grid%ny-1,grid%i1:grid%i2)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddy_a_to_a_3D
  SUBROUTINE ddxx_a_to_a_3D( grid, d_a, dxx_a)
    ! Input:  scalar on the Aa grid
    ! Output: its xx-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dxx_a
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      dxx_a( k,j,i) = (d_a( k,j,i+1) + d_a( k,j,i-1) - 2._dp * d_a( k,j,i)) / grid%dx**2
    END DO
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    dxx_a( :,grid%j1:grid%j2,1      ) = (d_a( :,grid%j1:grid%j2,3      ) + d_a( :,grid%j1:grid%j2,1        ) - 2._dp * d_a( :,grid%j1:grid%j2,2        )) / grid%dx**2
    dxx_a( :,grid%j1:grid%j2,grid%nx) = (d_a( :,grid%j1:grid%j2,grid%nx) + d_a( :,grid%j1:grid%j2,grid%nx-2) - 2._dp * d_a( :,grid%j1:grid%j2,grid%nx-1)) / grid%dx**2
    CALL sync
    
  END SUBROUTINE ddxx_a_to_a_3D
  SUBROUTINE ddyy_a_to_a_3D( grid, d_a, dyy_a)
    ! Input:  scalar on the Aa grid
    ! Output: its yy-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dyy_a
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      dyy_a( k,j,i) = (d_a( k,j+1,i) + d_a( k,j-1,i) - 2._dp * d_a( k,j,i)) / grid%dx**2
    END DO
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    dyy_a( :,1      ,grid%i1:grid%i2) = (d_a( :,3      ,grid%i1:grid%i2) + d_a( :,1        ,grid%i1:grid%i2) - 2._dp * d_a( :,2        ,grid%i1:grid%i2)) / grid%dx**2
    dyy_a( :,grid%ny,grid%i1:grid%i2) = (d_a( :,grid%ny,grid%i1:grid%i2) + d_a( :,grid%ny-2,grid%i1:grid%i2) - 2._dp * d_a( :,grid%ny-1,grid%i1:grid%i2)) / grid%dx**2
    CALL sync
    
  END SUBROUTINE ddyy_a_to_a_3D
  SUBROUTINE ddxy_a_to_a_3D( grid, d_a, dxy_a)
    ! Input:  scalar on the Aa grid
    ! Output: its xy-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dxy_a
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      dxy_a( k,j,i) = (d_a( k,j+1,i+1) + d_a( k,j-1,i-1) - d_a( k,j+1,i-1) - d_a( k,j-1,i+1)) / (4._dp * grid%dx * grid%dx)
    END DO
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    ! NO IDEA HOW TO DO THIS...
    dxy_a( :,1              ,grid%i1:grid%i2) = 0._dp
    dxy_a( :,grid%ny        ,grid%i1:grid%i2) = 0._dp
    dxy_a( :,grid%j1:grid%j2,1              ) = 0._dp
    dxy_a( :,grid%j1:grid%j2,grid%nx        ) = 0._dp
    CALL sync
    
  END SUBROUTINE ddxy_a_to_a_3D
  ! 3D upwind, for thermodynamics
  SUBROUTINE ddx_a_to_a_3D_upwind( grid, d_a, dx_a, U_3D_a)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Aa grid, using upwind one-sided differencing
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dx_a
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: U_3D_a
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Upwind one-sided differencing
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      IF (U_3D_a( k,j,i) > 0._dp) THEN
        dx_a( k,j,i) = (d_a( k,j,i  ) - d_a( k,j,i-1)) / grid%dx
      ELSE
        dx_a( k,j,i) = (d_a( k,j,i+1) - d_a( k,j,i  )) / grid%dx
      END IF
    END DO
    END DO
    END DO
    CALL sync
    
    dx_a( :,grid%j1:grid%j2,1      ) = 0._dp
    dx_a( :,grid%j1:grid%j2,grid%nx) = 0._dp
    CALL sync
    
  END SUBROUTINE ddx_a_to_a_3D_upwind
  SUBROUTINE ddy_a_to_a_3D_upwind( grid, d_a, dy_a, V_3D_a)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Aa grid, using upwind one-sided differencing
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dy_a
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: V_3D_a
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Upwind one-sided differencing
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      IF (V_3D_a( k,j,i) > 0._dp) THEN
        dy_a( k,j,i) = (d_a( k,j  ,i) - d_a( k,j-1,i)) / grid%dx
      ELSE
        dy_a( k,j,i) = (d_a( k,j+1,i) - d_a( k,j  ,i)) / grid%dx
      END IF
    END DO
    END DO
    END DO
    CALL sync
    
    dy_a( :,1      ,grid%i1:grid%i2) = 0._dp
    dy_a( :,grid%ny,grid%i1:grid%i2) = 0._dp
    CALL sync
    
  END SUBROUTINE ddy_a_to_a_3D_upwind
  
  ! Aa to Acx/Acy
  
  ! 2D
  SUBROUTINE ddx_a_to_cx_2D( grid, d_a, dx_cx)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: dx_cx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      dx_cx( j,i) = (d_a( j,i+1) - d_a( j,i)) / grid%dx
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddx_a_to_cx_2D
  SUBROUTINE ddy_a_to_cy_2D( grid, d_a, dy_cy)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: dy_cy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      dy_cy( j,i) = (d_a( j+1,i) - d_a( j,i)) / grid%dx
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddy_a_to_cy_2D
  SUBROUTINE ddx_a_to_cy_2D( grid, d_a, dx_cy)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: dx_cy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      dx_cy( j,i) = (d_a( j,i+1) + d_a( j+1,i+1) - d_a( j,i-1) - d_a( j+1,i-1)) / (4._dp * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundary
    DO j = grid%j1, MIN(grid%ny-1,grid%j2)
      dx_cy( j,1      ) = (d_a( j,2      ) + d_a( j+1,2      ) - d_a( j,1        ) - d_a( j+1,1        )) / (2._dp * grid%dx)
      dx_cy( j,grid%nx) = (d_a( j,grid%nx) + d_a( j+1,grid%nx) - d_a( j,grid%nx-1) - d_a( j+1,grid%nx-1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddx_a_to_cy_2D
  SUBROUTINE ddy_a_to_cx_2D( grid, d_a, dy_cx)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: dy_cx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dy_cx( j,i) = (d_a( j+1,i) + d_a( j+1,i+1) - d_a( j-1,i) - d_a( j-1,i+1)) / (4._dp * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundary
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      dy_cx( 1      ,i) = (d_a( 2,      i) + d_a( 2,      i+1) - d_a( 1,        i) - d_a( 1,        i+1)) / (2._dp * grid%dx)
      dy_cx( grid%ny,i) = (d_a( grid%ny,i) + d_a( grid%ny,i+1) - d_a( grid%ny-1,i) - d_a( grid%ny-1,i+1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddy_a_to_cx_2D
  ! 3D
  SUBROUTINE ddx_a_to_cx_3D( grid, d_a, dx_cx)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), INTENT(OUT)   :: dx_cx
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      dx_cx( k,j,i) = (d_a( k,j,i+1) - d_a( k,j,i)) / grid%dx
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddx_a_to_cx_3D
  SUBROUTINE ddy_a_to_cy_3D( grid, d_a, dy_cy)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), INTENT(OUT)   :: dy_cy
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
    DO k = 1, C%nZ
      dy_cy( k,j,i) = (d_a( k,j+1,i) - d_a( k,j,i)) / grid%dx
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddy_a_to_cy_3D
  SUBROUTINE ddx_a_to_cy_3D( grid, d_a, dx_cy)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), INTENT(OUT)   :: dx_cy
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
    DO k = 1, C%nZ
      dx_cy( k,j,i) = (d_a( k,j,i+1) + d_a( k,j+1,i+1) - d_a( k,j,i-1) - d_a( k,j+1,i-1)) / (4._dp * grid%dx)
    END DO
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundary
    DO j = grid%j1, MIN(grid%ny-1,grid%j2)
      dx_cy( :,j,1      ) = (d_a( :,j,2      ) + d_a( :,j+1,2      ) - d_a( :,j,1        ) - d_a( :,j+1,1        )) / (2._dp * grid%dx)
      dx_cy( :,j,grid%nx) = (d_a( :,j,grid%nx) + d_a( :,j+1,grid%nx) - d_a( :,j,grid%nx-1) - d_a( :,j+1,grid%nx-1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddx_a_to_cy_3D
  SUBROUTINE ddy_a_to_cx_3D( grid, d_a, dy_cx)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), INTENT(OUT)   :: dy_cx
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      dy_cx( k,j,i) = (d_a( k,j+1,i) + d_a( k,j+1,i+1) - d_a( k,j-1,i) - d_a( k,j-1,i+1)) / (4._dp * grid%dx)
    END DO
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundary
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      dy_cx( :,1      ,i) = (d_a( :,2,      i) + d_a( :,2,      i+1) - d_a( :,1,        i) - d_a( :,1,        i+1)) / (2._dp * grid%dx)
      dy_cx( :,grid%nx,i) = (d_a( :,grid%nx,i) + d_a( :,grid%nx,i+1) - d_a( :,grid%nx-1,i) - d_a( :,grid%nx-1,i+1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddy_a_to_cx_3D
  
  ! Acx/Acy to Aa
  
  ! 2D
  SUBROUTINE ddx_cx_to_a_2D( grid, d_cx, dx_a)
    ! Input:  scalar on the Acx grid
    ! Output: its x-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dx_a
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      dx_a( j,i) = (d_cx( j,i) - d_cx( j,i-1)) / grid%dx
    END DO
    END DO
    CALL sync
    
    dx_a( grid%j1:grid%j2,1      ) = dx_a( grid%j1:grid%j2,2        )
    dx_a( grid%j1:grid%j2,grid%nx) = dx_a( grid%j1:grid%j2,grid%nx-1)
    CALL sync
    
  END SUBROUTINE ddx_cx_to_a_2D
  SUBROUTINE ddy_cy_to_a_2D( grid, d_cy, dy_a)
    ! Input:  scalar on the Acy grid
    ! Output: its y-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dy_a
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
      dy_a( j,i) = (d_cy( j,i) - d_cy( j-1,i)) / grid%dx
    END DO
    END DO
    CALL sync
    
    dy_a( 1      ,grid%i1:grid%i2) = dy_a( 2        ,grid%i1:grid%i2)
    dy_a( grid%ny,grid%i1:grid%i2) = dy_a( grid%ny-1,grid%i1:grid%i2)
    CALL sync
    
  END SUBROUTINE ddy_cy_to_a_2D
  SUBROUTINE ddy_cx_to_a_2D( grid, d_cx, dy_a)
    ! Input:  scalar on the Acx grid
    ! Output: its y-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dy_a
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dy_a( j,i) = (d_cx( j+1,i-1) + d_cx( j+1,i) - d_cx( j-1,i-1) - d_cx( j-1,i)) / (4._dp * grid%dx)
    END DO
    END DO
    CALL sync
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
      ! South ex. corners
      j = 1
      dy_a( j,i) = (d_cx( j+1,i-1) + d_cx( j+1,i) - d_cx( j  ,i-1) - d_cx( j  ,i)) / (4._dp * grid%dx)
      ! North ex. corners
      j = grid%ny
      dy_a( j,i) = (d_cx( j  ,i-1) + d_cx( j  ,i) - d_cx( j-1,i-1) - d_cx( j-1,i)) / (4._dp * grid%dx)
    END DO
    CALL sync
    
    DO j = MAX(2,grid%j1), MIN(grid%ny-1,grid%j2)
      ! West ex. corners
      i = 1
      dy_a( j,i) = (d_cx( j+1,i  ) - d_cx( j-1,i  )) / (2._dp * grid%dx)
      ! East ex. corners
      i = grid%nx
      dy_a( j,i) = (d_cx( j+1,i-1) - d_cx( j-1,i-1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
    ! Corners
    IF (par%master) THEN
    dy_a( 1      ,1      ) = (d_cx( 2      ,1        ) - d_cx( 1        ,1        )) / grid%dx
    dy_a( 1      ,grid%nx) = (d_cx( 2      ,grid%nx-1) - d_cx( 1        ,grid%nx-1)) / grid%dx
    dy_a( grid%ny,1      ) = (d_cx( grid%ny,1        ) - d_cx( grid%ny-1,1        )) / grid%dx
    dy_a( grid%ny,grid%nx) = (d_cx( grid%ny,grid%nx-1) - d_cx( grid%ny-1,grid%nx-1)) / grid%dx
    END IF
    CALL sync
    
  END SUBROUTINE ddy_cx_to_a_2D
  SUBROUTINE ddx_cy_to_a_2D( grid, d_cy, dx_a)
    ! Input:  scalar on the Acy grid
    ! Output: its x-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dx_a
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dx_a( j,i) = (d_cy( j-1,i+1) + d_cy( j,i+1) - d_cy( j-1,i-1) - d_cy( j,i-1)) / (4._dp * grid%dx)
    END DO
    END DO
    CALL sync
    
    DO j = MAX(2,grid%j1), MIN(grid%ny-1,grid%j2)
      ! West ex. corners
      i = 1
      dx_a( j,i) = (d_cy( j-1,i+1) + d_cy( j  ,i+1) - d_cy( j-1,i  ) - d_cy( j  ,i  )) / (4._dp * grid%dx)
      ! East ex. corners
      i = grid%nx
      dx_a( j,i) = (d_cy( j-1,i  ) + d_cy( j  ,i  ) - d_cy( j-1,i-1) - d_cy( j  ,i-1)) / (4._dp * grid%dx)
    END DO
    CALL sync
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
      ! South ex. corners
      j = 1
      dx_a( j,i) = (d_cy( j  ,i+1) - d_cy( j  ,i-1)) / (2._dp * grid%dx)
      ! North ex. corners
      j = grid%ny
      dx_a( j,i) = (d_cy( j-1,i+1) - d_cy( j-1,i-1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
    ! Corners
    IF (par%master) THEN
    dx_a( 1      ,      1) = (d_cy( 1        ,2      ) - d_cy( 1        ,1        )) / grid%dx
    dx_a( 1      ,grid%nx) = (d_cy( 1        ,grid%nx) - d_cy( 1        ,grid%nx-1)) / grid%dx
    dx_a( grid%ny,1      ) = (d_cy( grid%ny-1,2      ) - d_cy( grid%ny-1,1        )) / grid%dx
    dx_a( grid%ny,grid%nx) = (d_cy( grid%ny-1,grid%nx) - d_cy( grid%ny-1,grid%nx-1)) / grid%dx
    END IF
    CALL sync
    
  END SUBROUTINE ddx_cy_to_a_2D
  ! 3D
  SUBROUTINE ddx_cx_to_a_3D( grid, d_cx, dx_a)
    ! Input:  scalar on the Acx grid
    ! Output: its x-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dx_a
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      dx_a( k,j,i) = (d_cx( k,j,i) - d_cx( k,j,i-1)) / grid%dx
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddx_cx_to_a_3D
  SUBROUTINE ddy_cy_to_a_3D( grid, d_cy, dy_a)
    ! Input:  scalar on the Acy grid
    ! Output: its y-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dy_a
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      dy_a( k,j,i) = (d_cy( k,j,i) - d_cy( k,j-1,i)) / grid%dx
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddy_cy_to_a_3D
  
  ! Acx/Acy to Acx/Acy
  SUBROUTINE ddx_cx_to_cx_2D( grid, d_cx, dx_cx)
    ! Input:  scalar on the Acx grid
    ! Output: its x-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: dx_cx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = 2, grid%nx-2
    DO j = 1, grid%ny
      dx_cx( j,i) = (d_cx( j,i+1) - d_cx( j,i-1)) / (2 * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    dx_cx( grid%j1:grid%j2,1        ) = (d_cx( grid%j1:grid%j2,2        ) - d_cx( grid%j1:grid%j2,1        )) / grid%dx
    dx_cx( grid%j1:grid%j2,grid%nx-1) = (d_cx( grid%j1:grid%j2,grid%nx-1) - d_cx( grid%j1:grid%j2,grid%nx-2)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddx_cx_to_cx_2D
  SUBROUTINE ddy_cx_to_cx_2D( grid, d_cx, dy_cx)
    ! Input:  scalar on the Acx grid
    ! Output: its y-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: dy_cx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dy_cx( j,i) = (d_cx( j+1,i) - d_cx( j-1,i)) / (2 * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    dy_cx( 1      ,grid%i1:grid%i2) = (d_cx( 2      ,grid%i1:grid%i2) - d_cx( 1        ,grid%i1:grid%i2)) / grid%dx
    dy_cx( grid%ny,grid%i1:grid%i2) = (d_cx( grid%ny,grid%i1:grid%i2) - d_cx( grid%ny-1,grid%i1:grid%i2)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddy_cx_to_cx_2D
  SUBROUTINE ddx_cy_to_cy_2D( grid, d_cy, dx_cy)
    ! Input:  scalar on the Acy grid
    ! Output: its x-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: dx_cy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      dx_cy( j,i) = (d_cy( j,i+1) - d_cy( j,i-1)) / (2 * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    dx_cy( grid%j1:MIN(grid%j2,grid%ny-1),1      ) = (d_cy( grid%j1:MIN(grid%j2,grid%ny-1),2      ) - d_cy( grid%j1:MIN(grid%j2,grid%ny-1),1        )) / grid%dx
    dx_cy( grid%j1:MIN(grid%j2,grid%ny-1),grid%nx) = (d_cy( grid%j1:MIN(grid%j2,grid%ny-1),grid%nx) - d_cy( grid%j1:MIN(grid%j2,grid%ny-1),grid%nx-1)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddx_cy_to_cy_2D
  SUBROUTINE ddy_cy_to_cy_2D( grid, d_cy, dy_cy)
    ! Input:  scalar on the Acy grid
    ! Output: its y-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: dy_cy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-2
      dy_cy( j,i) = (d_cy( j+1,i) - d_cy( j-1,i)) / (2 * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! One-sided differencing on the boundaries
    dy_cy( 1        ,grid%i1:grid%i2) = (d_cy( 2        ,grid%i1:grid%i2) - d_cy( 1        ,grid%i1:grid%i2)) / grid%dx
    dy_cy( grid%ny-1,grid%i1:grid%i2) = (d_cy( grid%ny-1,grid%i1:grid%i2) - d_cy( grid%ny-2,grid%i1:grid%i2)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddy_cy_to_cy_2D
  SUBROUTINE ddx_cx_to_cy_2D( grid, d_cx, dx_cy)
    ! Input:  scalar on the Acx grid
    ! Output: its x-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: dx_cy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      dx_cy( j,i) = (d_cx( j,i  ) + d_cx( j+1,i  ) - d_cx( j,i-1) - d_cx( j+1,i-1)) / (2._dp * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! Boundaries
    DO j = grid%j1, MIN(grid%ny-1,grid%j2)
      ! West
      dx_cy( j,1      ) = dx_cy( j,2        )
      ! East
      dx_cy( j,grid%nx) = dx_cy( j,grid%nx-1)
    END DO
    CALL sync
    
  END SUBROUTINE ddx_cx_to_cy_2D
  SUBROUTINE ddy_cx_to_cy_2D( grid, d_cx, dy_cy)
    ! Input:  scalar on the Acx grid
    ! Output: its y-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: dy_cy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      dy_cy( j,i) = (d_cx( j+1,i-1) + d_cx( j+1,i  ) - d_cx( j  ,i-1) - d_cx( j  ,i  )) / (2._dp * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! Boundaries
    DO j = grid%j1, MIN(grid%ny-1,grid%j2)
      ! West
      i = 1
      dy_cy( j,i) = (d_cx( j,i  ) - d_cx( j,i  )) / grid%dx
      ! East
      i = grid%nx
      dy_cy( j,i) = (d_cx( j,i-1) - d_cx( j,i-1)) / grid%dx
    END DO
    CALL sync
    
  END SUBROUTINE ddy_cx_to_cy_2D
  SUBROUTINE ddx_cy_to_cx_2D( grid, d_cy, dx_cx)
    ! Input:  scalar on the Acy grid
    ! Output: its x-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: dx_cx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dx_cx( j,i) = (d_cy( j-1,i+1) + d_cy( j  ,i+1) - d_cy( j-1,i  ) - d_cy( j  ,i  )) / (2._dp * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! Boundaries
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      ! South
      j = 1
      dx_cx( j,i) = (d_cy( j  ,i+1) - d_cy( j  ,i  )) / grid%dx
      ! North
      j = grid%ny
      dx_cx( j,i) = (d_cy( j-1,i+1) - d_cy( j-1,i  )) / grid%dx
    END DO
    CALL sync
    
  END SUBROUTINE ddx_cy_to_cx_2D
  SUBROUTINE ddy_cy_to_cx_2D( grid, d_cy, dy_cx)
    ! Input:  scalar on the Acy grid
    ! Output: its y-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: dy_cx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dy_cx( j,i) = (d_cy( j  ,i  ) + d_cy( j  ,i+1) - d_cy( j-1,i  ) - d_cy( j-1,i+1)) / (2._dp * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! Boundaries
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      ! South
      j = 1
      dy_cx( j,i) = (d_cy( j+1,i  ) + d_cy(j+1,i+1) - d_cy( j  ,i  ) - d_cy( j  ,i+1)) / (2._dp * grid%dx)
      ! North
      j = grid%ny
      dy_cx( j,i) = (d_cy( j-1,i  ) + d_cy(j-1,i+1) - d_cy( j-2,i  ) - d_cy( j-2,i+1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddy_cy_to_cx_2D
  
  ! Acx to Ab
  SUBROUTINE ddx_cx_to_b_2D( grid, d_cx, dx_b)
    ! Input:  scalar on the Acx grid
    ! Output: its x-derivative on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: dx_b
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-2,grid%i2)
    DO j = 1, grid%ny-1
      dx_b( j,i) = (d_cx( j+1,i+1) + d_cx( j  ,i+1) - d_cx( j+1,i-1) - d_cx( j  ,i-1)) / (4._dp * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! Boundaries
    DO j = grid%j1, MIN(grid%ny-1,grid%j2)
      i = 1
      dx_b( j,i) = (d_cx( j+1,i+1) + d_cx( j  ,i+1) - d_cx( j+1,i  ) - d_cx( j  ,i  )) / (2._dp * grid%dx)
      i = grid%nx-1
      dx_b( j,i) = (d_cx( j+1,i  ) + d_cx( j  ,i  ) - d_cx( j+1,i-1) - d_cx( j  ,i-1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddx_cx_to_b_2D
  SUBROUTINE ddy_cy_to_b_2D( grid, d_cy, dy_b)
    ! Input:  scalar on the Acy grid
    ! Output: its y-derivative on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: dy_b
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-2,grid%i2)
    DO j = 2, grid%ny-2
      dy_b( j,i) = (d_cy( j+1,i+1) + d_cy( j+1,i  ) - d_cy( j-1,i+1) - d_cy( j-1,i  )) / (4._dp * grid%dx)
    END DO
    END DO
    CALL sync
    
    ! Boundaries
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      j = 1
      dy_b( j,i) = (d_cy( j+1,i+1) + d_cy( j+1,i  ) - d_cy( j  ,i+1) - d_cy( j  ,i  )) / (2._dp * grid%dx)
      j = grid%ny-1
      dy_b( j,i) = (d_cy( j  ,i+1) + d_cy( j  ,i  ) - d_cy( j-1,i+1) - d_cy( j-1,i  )) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddy_cy_to_b_2D
  SUBROUTINE ddx_cy_to_b_2D( grid, d_cy, dx_b)
    ! Input:  scalar on the Acy grid
    ! Output: its x-derivative on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: dx_b
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      dx_b( j,i) = (d_cy( j,i+1) - d_cy( j,i)) / grid%dx
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddx_cy_to_b_2D
  SUBROUTINE ddy_cx_to_b_2D( grid, d_cx, dy_b)
    ! Input:  scalar on the Acx grid
    ! Output: its y-derivative on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: dy_b
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      dy_b( j,i) = (d_cx( j+1,i) - d_cx( j,i)) / grid%dx
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddy_cx_to_b_2D
  
! =============================================
! ===== Mapping between (staggered) grids =====
! =============================================

  ! Aa to Acx/Acy
  
  ! 2D
  SUBROUTINE map_a_to_cx_2D( grid, d_a, d_cx)
    ! Input:  scalar on the Aa grid
    ! Output: the same on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: d_cx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      d_cx( j,i) = (d_a( j,i) + d_a( j,i+1)) / 2._dp
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_a_to_cx_2D
  SUBROUTINE map_a_to_cy_2D( grid, d_a, d_cy)
    ! Input:  scalar on the Aa grid
    ! Output: the same on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: d_cy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      d_cy( j,i) = (d_a( j,i) + d_a( j+1,i)) / 2._dp
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_a_to_cy_2D
  ! 3D
  SUBROUTINE map_a_to_cx_3D( grid, d_a, d_cx)
    ! Input:  scalar on the Aa grid
    ! Output: the same on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), INTENT(OUT)   :: d_cx
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      d_cx( k,j,i) = (d_a( k,j,i) + d_a( k,j,i+1)) / 2._dp
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_a_to_cx_3D
  SUBROUTINE map_a_to_cy_3D( grid, d_a, d_cy)
    ! Input:  scalar on the Aa grid
    ! Output: the same on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), INTENT(OUT)   :: d_cy
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
    DO k = 1, C%nZ
      d_cy( k,j,i) = (d_a( k,j,i) + d_a( k,j+1,i)) / 2._dp
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_a_to_cy_3D
  
  ! Acx/Acy to Aa
  
  ! 2D
  SUBROUTINE map_cx_to_a_2D( grid, d_cx, d_a)
    ! Input:  scalar on the Acx grid
    ! Output: the same on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: d_a
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      d_a( j,i) = (d_cx( j,i-1) + d_cx( j,i)) / 2._dp
    END DO
    END DO
    CALL sync
    
    d_a( grid%j1:grid%j2,1      ) = d_cx( grid%j1:grid%j2,1        )
    d_a( grid%j1:grid%j2,grid%nx) = d_cx( grid%j1:grid%j2,grid%nx-1)
    CALL sync
    
  END SUBROUTINE map_cx_to_a_2D
  SUBROUTINE map_cy_to_a_2D( grid, d_cy, d_a)
    ! Input:  scalar on the Acy grid
    ! Output: the same on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: d_a
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
      d_a( j,i) = (d_cy( j-1,i) + d_cy( j,i)) / 2._dp
    END DO
    END DO
    CALL sync
    
    d_a( 1      ,grid%i1:grid%i2) = d_cy( 1        ,grid%i1:grid%i2)
    d_a( grid%ny,grid%i1:grid%i2) = d_cy( grid%ny-1,grid%i1:grid%i2)
    CALL sync
    
  END SUBROUTINE map_cy_to_a_2D
  ! 3D
  SUBROUTINE map_cx_to_a_3D( grid, d_cx, d_a)
    ! Input:  scalar on the Acx grid
    ! Output: the same on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: d_a
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      d_a( k,j,i) = (d_cx( k,j,i-1) + d_cx( k,j,i)) / 2._dp
    END DO
    END DO
    END DO
    CALL sync
    
    d_a( :,grid%j1:grid%j2,1      ) = d_cx( :,grid%j1:grid%j2,1        )
    d_a( :,grid%j1:grid%j2,grid%nx) = d_cx( :,grid%j1:grid%j2,grid%nx-1)
    CALL sync
    
  END SUBROUTINE map_cx_to_a_3D
  SUBROUTINE map_cy_to_a_3D( grid, d_cy, d_a)
    ! Input:  scalar on the Acy grid
    ! Output: the same on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: d_a
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      d_a( k,j,i) = (d_cy( k,j-1,i) + d_cy( k,j,i)) / 2._dp
    END DO
    END DO
    END DO
    CALL sync
    
    d_a( :,1      ,grid%i1:grid%i2) = d_cy( :,1        ,grid%i1:grid%i2)
    d_a( :,grid%ny,grid%i1:grid%i2) = d_cy( :,grid%ny-1,grid%i1:grid%i2)
    CALL sync
    
  END SUBROUTINE map_cy_to_a_3D
  
  ! Acx/Acy to Acy/Acx
  SUBROUTINE map_cx_to_cy_2D( grid, d_cx, d_cy)
    ! Input:  scalar on the Acx grid
    ! Output: the same on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: d_cy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      d_cy( j,i) = (d_cx( j  ,i-1) + d_cx( j  ,i  ) + d_cx( j+1,i-1) + d_cx( j+1,i  )) / 4._dp
    END DO
    END DO
    CALL sync
    
    ! Boundaries
    DO j = grid%j1, MIN(grid%ny-1,grid%j2)
      i = 1
      d_cy( j,i) = (d_cx( j  ,i  ) + d_cx( j+1,i  )) / 2._dp
      i = grid%nx
      d_cy( j,i) = (d_cx( j  ,i-1) + d_cx( j+1,i-1)) / 2._dp
    END DO
    CALL sync
    
  END SUBROUTINE map_cx_to_cy_2D
  SUBROUTINE map_cy_to_cx_2D( grid, d_cy, d_cx)
    ! Input:  scalar on the Acy grid
    ! Output: the same on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: d_cx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      d_cx( j,i) = (d_cy( j-1,i  ) + d_cy( j-1,i+1) + d_cy( j  ,i  ) + d_cy( j  ,i+1)) / 4._dp
    END DO
    END DO
    CALL sync
    
    ! Boundaries
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      j = 1
      d_cx( j,i) = (d_cy( j  ,i  ) + d_cy( j  ,i+1)) / 2._dp
      j = grid%ny
      d_cx( j,i) = (d_cy( j-1,i  ) + d_cy( j-1,i+1)) / 2._dp
    END DO
    CALL sync
    
  END SUBROUTINE map_cy_to_cx_2D
  
  ! Aa to Ab
  SUBROUTINE map_a_to_b_2D( grid, d_a, d_b)
    ! Input:  scalar on the Aa grid
    ! Output: the same on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_a
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: d_b
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      d_b( j,i) = (d_a( j,i) + d_a( j,i+1) + d_a( j+1,i) + d_a( j+1,i+1)) / 4._dp
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_a_to_b_2D
  
  ! Acx/Acy to Ab
  SUBROUTINE map_cx_to_b_2D( grid, d_cx, d_b)
    ! Input:  scalar on the Acx grid
    ! Output: the same on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_cx
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: d_b
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      d_b( j,i) = (d_cx( j,i) + d_cx( j+1,i)) / 2._dp
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_cx_to_b_2D
  SUBROUTINE map_cy_to_b_2D( grid, d_cy, d_b)
    ! Input:  scalar on the Acy grid
    ! Output: the same on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_cy
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: d_b
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      d_b( j,i) = (d_cy( j,i) + d_cy( j,i+1)) / 2._dp
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_cy_to_b_2D
  
! ============================
! ===== Zeta derivatives =====
! ============================

  SUBROUTINE calculate_zeta_derivatives( grid, ice)
    ! Calculate the spatial and temporal derivatives of the scaled vertical coordinate zeta,
    ! required for the calculation of 3D spatial derivatives in the thermodynamics routine
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    TYPE(type_ice_model),                             INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    REAL(dp)                                                        :: inverse_Hi

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      inverse_Hi  = 1._dp / ice%Hi_a( j,i)

      ice%dzeta_dz_a( j,i)  = -inverse_Hi
      DO k = 1, C%nZ
        ice%dzeta_dt_a(  k,j,i) = inverse_Hi * (ice%dHs_dt_a(  j,i) - C%zeta(k) * ice%dHi_dt_a(  j,i))
        ice%dzeta_dx_a(  k,j,i) = inverse_Hi * (ice%dHs_dx_a(  j,i) - C%zeta(k) * ice%dHi_dx_a(  j,i))
        ice%dzeta_dy_a(  k,j,i) = inverse_Hi * (ice%dHs_dy_a(  j,i) - C%zeta(k) * ice%dHi_dy_a(  j,i))
      END DO
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calculate_zeta_derivatives
  SUBROUTINE initialise_zeta_discretisation
    ! Initialise the coefficients for calculating df/dzeta
    ! (used in thermodynamics)
  
    IMPLICIT NONE

    ! Local variables:
    INTEGER :: k

    ALLOCATE( zeta%a_k(          2:C%nz  ))
    ALLOCATE( zeta%b_k(          1:C%nz-1))
    ALLOCATE( zeta%c_k(          3:C%nz  ))
    ALLOCATE( zeta%d_k(          1:C%nz-2))
    ALLOCATE( zeta%a_zeta(       2:C%nz-1))
    ALLOCATE( zeta%b_zeta(       2:C%nz-1))
    ALLOCATE( zeta%c_zeta(       2:C%nz-1))
    ALLOCATE( zeta%a_zetazeta(   2:C%nz-1))
    ALLOCATE( zeta%b_zetazeta(   2:C%nz-1))
    ALLOCATE( zeta%c_zetazeta(   2:C%nz-1))

    ALLOCATE( zeta%z_zeta_minus( 3:C%nz  ))
    ALLOCATE( zeta%a_zeta_minus( 3:C%nz  ))
    ALLOCATE( zeta%b_zeta_minus( 3:C%nz  ))
    ALLOCATE( zeta%b_zeta_plus(  1:C%nz-2))
    ALLOCATE( zeta%c_zeta_plus(  1:C%nz-2))
    ALLOCATE( zeta%d_zeta_plus(  1:C%nz-2))

    DO k = 2, C%nz
      zeta%a_k( k) = C%zeta( k  ) - C%zeta( k-1)
    END DO
    DO k = 1, C%nz-1
      zeta%b_k( k) = C%zeta( k+1) - C%zeta( k  )
    END DO
    DO k = 3, C%nz
      zeta%c_k( k) = C%zeta( k  ) - C%zeta( k-2)
    END DO
    DO k = 1, C%nz-2
      zeta%d_k( k) = C%zeta( k+2) - C%zeta( k  )
    END DO

    DO k = 2, C%nz-1
      zeta%a_zeta( k)     =                - zeta%b_k( k)  / ( zeta%a_k( k) * ( zeta%a_k( k) + zeta%b_k( k)))
      zeta%b_zeta( k)     = ( zeta%b_k( k) - zeta%a_k( k)) / ( zeta%a_k( k) *                  zeta%b_k( k) )
      zeta%c_zeta( k)     =   zeta%a_k( k)                 / ( zeta%b_k( k) * ( zeta%a_k( k) + zeta%b_k( k)))
      zeta%a_zetazeta( k) =                      2.0_dp    / ( zeta%a_k( k) * ( zeta%a_k( k) + zeta%b_k( k)))
      zeta%b_zetazeta( k) =                     -2.0_dp    / ( zeta%a_k( k) *                  zeta%b_k( k) )
      zeta%c_zetazeta( k) =                      2.0_dp    / ( zeta%b_k( k) * ( zeta%a_k( k) + zeta%b_k( k)))
    END DO

    ! Not all of these are in use:
    DO k = 1, C%nz-2
      zeta%b_zeta_plus( k) = -( zeta%b_k( k) + zeta%d_k( k)) / ( zeta%b_k( k) *   zeta%d_k( k)                )
      zeta%c_zeta_plus( k) =                   zeta%d_k( k)  / ( zeta%b_k( k) * ( zeta%d_k( k) - zeta%b_k( k)))
      zeta%d_zeta_plus( k) =                   zeta%b_k( k)  / ( zeta%d_k( k) * ( zeta%b_k( k) - zeta%d_k( k)))
    END DO

    ! Not all of these are in use:
    DO k = 3, C%nz
      zeta%z_zeta_minus( k) =                  zeta%a_k( k)  / ( zeta%c_k( k) * ( zeta%c_k( k) - zeta%a_k( k)))
      zeta%a_zeta_minus( k) =                  zeta%c_k( k)  / ( zeta%a_k( k) * ( zeta%a_k( k) - zeta%c_k( k)))
      zeta%b_zeta_minus( k) = ( zeta%a_k( k) + zeta%c_k( k)) / ( zeta%a_k( k) *   zeta%c_k( k)                )
    END DO
    
  END SUBROUTINE initialise_zeta_discretisation
  
! ======================================
! ===== Neumann boundary condition =====
! ======================================

  SUBROUTINE Neumann_BC_a_2D( grid, d_a)
    ! Input:  scalar on the Aa grid
    ! Output: the same, with boundary values set to next-to-boundary values,
    !         so that d/dx=0 on the west&east boundaries, and d/dy-0 on north&south
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(INOUT) :: d_a
    
    d_a( 1              ,grid%i1:grid%i2) = d_a( 2              ,grid%i1:grid%i2)
    d_a( grid%ny        ,grid%i1:grid%i2) = d_a( grid%ny-1      ,grid%i1:grid%i2)
    CALL sync
    d_a( grid%j1:grid%j2,1              ) = d_a( grid%j1:grid%j2,2              )
    d_a( grid%j1:grid%j2,grid%nx        ) = d_a( grid%j1:grid%j2,grid%nx-1      )
    CALL sync
    
  END SUBROUTINE Neumann_BC_a_2D
  SUBROUTINE Neumann_BC_a_3D( grid, d_a)
    ! Input:  scalar on the Aa grid
    ! Output: the same, with boundary values set to next-to-boundary values,
    !         so that d/dx=0 on the west&east boundaries, and d/dy-0 on north&south
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(INOUT) :: d_a
    
    d_a( :,1              ,grid%i1:grid%i2) = d_a( :,2              ,grid%i1:grid%i2)
    d_a( :,grid%ny        ,grid%i1:grid%i2) = d_a( :,grid%ny-1      ,grid%i1:grid%i2)
    CALL sync
    d_a( :,grid%j1:grid%j2,1              ) = d_a( :,grid%j1:grid%j2,2              )
    d_a( :,grid%j1:grid%j2,grid%nx        ) = d_a( :,grid%j1:grid%j2,grid%nx-1      )
    CALL sync
    
  END SUBROUTINE Neumann_BC_a_3D
  SUBROUTINE Neumann_BC_cx_2D( grid, d_cx)
    ! Input:  scalar on the Aa grid
    ! Output: the same, with boundary values set to next-to-boundary values,
    !         so that d/dx=0 on the west&east boundaries, and d/dy-0 on north&south
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(INOUT) :: d_cx
    
    d_cx( 1              ,grid%i1:MIN(grid%nx-1,grid%i2)) = d_cx( 2              ,grid%i1:MIN(grid%nx-1,grid%i2))
    d_cx( grid%ny        ,grid%i1:MIN(grid%nx-1,grid%i2)) = d_cx( grid%ny-1      ,grid%i1:MIN(grid%nx-1,grid%i2))
    CALL sync
    d_cx( grid%j1:grid%j2,1              ) = d_cx( grid%j1:grid%j2,2              )
    d_cx( grid%j1:grid%j2,grid%nx        ) = d_cx( grid%j1:grid%j2,grid%nx-1      )
    CALL sync
    
  END SUBROUTINE Neumann_BC_cx_2D
  SUBROUTINE Neumann_BC_cy_2D( grid, d_cy)
    ! Input:  scalar on the Aa grid
    ! Output: the same, with boundary values set to next-to-boundary values,
    !         so that d/dx=0 on the west&east boundaries, and d/dy-0 on north&south
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(INOUT) :: d_cy
    
    d_cy( 1              ,grid%i1:grid%i2) = d_cy( 2              ,grid%i1:grid%i2)
    d_cy( grid%ny        ,grid%i1:grid%i2) = d_cy( grid%ny-1      ,grid%i1:grid%i2)
    CALL sync
    d_cy( grid%j1:MIN(grid%ny-1,grid%j2),1              ) = d_cy( grid%j1:MIN(grid%ny-1,grid%j2),2              )
    d_cy( grid%j1:MIN(grid%ny-1,grid%j2),grid%nx        ) = d_cy( grid%j1:MIN(grid%ny-1,grid%j2),grid%nx-1      )
    CALL sync
    
  END SUBROUTINE Neumann_BC_cy_2D

END MODULE derivatives_and_grids_module



