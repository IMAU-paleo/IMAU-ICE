! File name: derivatives_and_grids_module.f90
!
! This file is part of the IMAUICE model
!
! IMAU, Utrecht University, The Netherlands

MODULE derivatives_and_grids_module

  USE mpi
  USE parallel_module,                 ONLY: par, sync
  USE configuration_module,            ONLY: dp, C
  USE data_types_module,               ONLY: type_grid, type_ice_model

CONTAINS

! ======================
! ==== Derivatives =====
! ======================
  
  ! Aa to Aa
  
  ! 2D
  SUBROUTINE ddx_Aa_to_Aa_2D(  grid, d_Aa, dx_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dx_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      dx_Aa( j,i) = (d_Aa( j,i+1) - d_Aa( j,i-1)) / (2 * grid%dx)
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    dx_Aa( grid%j1:grid%j2,1      ) = (d_Aa( grid%j1:grid%j2,2      ) - d_Aa( grid%j1:grid%j2,1        )) / grid%dx
    dx_Aa( grid%j1:grid%j2,grid%nx) = (d_Aa( grid%j1:grid%j2,grid%nx) - d_Aa( grid%j1:grid%j2,grid%nx-1)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddx_Aa_to_Aa_2D
  SUBROUTINE ddy_Aa_to_Aa_2D(  grid, d_Aa, dy_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dy_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
      dy_Aa( j,i) = (d_Aa( j+1,i) - d_Aa( j-1,i)) / (2 * grid%dx)
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    dy_Aa( 1      ,grid%i1:grid%i2) = (d_Aa( 2      ,grid%i1:grid%i2) - d_Aa( 1        ,grid%i1:grid%i2)) / grid%dx
    dy_Aa( grid%ny,grid%i1:grid%i2) = (d_Aa( grid%ny,grid%i1:grid%i2) - d_Aa( grid%ny-1,grid%i1:grid%i2)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddy_Aa_to_Aa_2D
  SUBROUTINE ddxx_Aa_to_Aa_2D( grid, d_Aa, dxx_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: its xx-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dxx_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      dxx_Aa( j,i) = (d_Aa( j,i+1) + d_Aa( j,i-1) - 2._dp * d_Aa( j,i)) / grid%dx**2
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    dxx_Aa( grid%j1:grid%j2,1      ) = (d_Aa( grid%j1:grid%j2,3      ) + d_Aa( grid%j1:grid%j2,1        ) - 2._dp * d_Aa( grid%j1:grid%j2,2        )) / grid%dx**2
    dxx_Aa( grid%j1:grid%j2,grid%nx) = (d_Aa( grid%j1:grid%j2,grid%nx) + d_Aa( grid%j1:grid%j2,grid%nx-2) - 2._dp * d_Aa( grid%j1:grid%j2,grid%nx-1)) / grid%dx**2
    CALL sync
    
  END SUBROUTINE ddxx_Aa_to_Aa_2D
  SUBROUTINE ddyy_Aa_to_Aa_2D( grid, d_Aa, dyy_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: its yy-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dyy_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
      dyy_Aa( j,i) = (d_Aa( j+1,i) + d_Aa( j-1,i) - 2._dp * d_Aa( j,i)) / grid%dx**2
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    dyy_Aa( 1      ,grid%i1:grid%i2) = (d_Aa( 3      ,grid%i1:grid%i2) + d_Aa( 1        ,grid%i1:grid%i2) - 2._dp * d_Aa( 2        ,grid%i1:grid%i2)) / grid%dx**2
    dyy_Aa( grid%ny,grid%i1:grid%i2) = (d_Aa( grid%ny,grid%i1:grid%i2) + d_Aa( grid%ny-2,grid%i1:grid%i2) - 2._dp * d_Aa( grid%ny-1,grid%i1:grid%i2)) / grid%dx**2
    CALL sync
    
  END SUBROUTINE ddyy_Aa_to_Aa_2D
  SUBROUTINE ddxy_Aa_to_Aa_2D( grid, d_Aa, dxy_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: its xy-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dxy_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dxy_Aa( j,i) = (d_Aa( j+1,i+1) + d_Aa( j-1,i-1) - d_Aa( j+1,i-1) - d_Aa( j-1,i+1)) / (4._dp * grid%dx * grid%dx)
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    ! NO IDEA HOW TO DO THIS...
    dxy_Aa( 1              ,grid%i1:grid%i2) = 0._dp
    dxy_Aa( grid%ny        ,grid%i1:grid%i2) = 0._dp
    dxy_Aa( grid%j1:grid%j2,1              ) = 0._dp
    dxy_Aa( grid%j1:grid%j2,grid%nx        ) = 0._dp
    CALL sync
    
  END SUBROUTINE ddxy_Aa_to_Aa_2D
  ! 3D
  SUBROUTINE ddx_Aa_to_Aa_3D(  grid, d_Aa, dx_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dx_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      dx_Aa( k,j,i) = (d_Aa( k,j,i+1) - d_Aa( k,j,i-1)) / (2 * grid%dx)
    END DO
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    dx_Aa( :,grid%j1:grid%j2,1      ) = (d_Aa( :,grid%j1:grid%j2,2      ) - d_Aa( :,grid%j1:grid%j2,1        )) / grid%dx
    dx_Aa( :,grid%j1:grid%j2,grid%nx) = (d_Aa( :,grid%j1:grid%j2,grid%nx) - d_Aa( :,grid%j1:grid%j2,grid%nx-1)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddx_Aa_to_Aa_3D
  SUBROUTINE ddy_Aa_to_Aa_3D(  grid, d_Aa, dy_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dy_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      dy_Aa( k,j,i) = (d_Aa( k,j+1,i) - d_Aa( k,j-1,i)) / (2 * grid%dx)
    END DO
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    dy_Aa( :,1      ,grid%i1:grid%i2) = (d_Aa( :,2      ,grid%i1:grid%i2) - d_Aa( :,1        ,grid%i1:grid%i2)) / grid%dx
    dy_Aa( :,grid%ny,grid%i1:grid%i2) = (d_Aa( :,grid%ny,grid%i1:grid%i2) - d_Aa( :,grid%ny-1,grid%i1:grid%i2)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddy_Aa_to_Aa_3D
  SUBROUTINE ddxx_Aa_to_Aa_3D( grid, d_Aa, dxx_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: its xx-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dxx_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      dxx_Aa( k,j,i) = (d_Aa( k,j,i+1) + d_Aa( k,j,i-1) - 2._dp * d_Aa( k,j,i)) / grid%dx**2
    END DO
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    dxx_Aa( :,grid%j1:grid%j2,1      ) = (d_Aa( :,grid%j1:grid%j2,3      ) + d_Aa( :,grid%j1:grid%j2,1        ) - 2._dp * d_Aa( :,grid%j1:grid%j2,2        )) / grid%dx**2
    dxx_Aa( :,grid%j1:grid%j2,grid%nx) = (d_Aa( :,grid%j1:grid%j2,grid%nx) + d_Aa( :,grid%j1:grid%j2,grid%nx-2) - 2._dp * d_Aa( :,grid%j1:grid%j2,grid%nx-1)) / grid%dx**2
    CALL sync
    
  END SUBROUTINE ddxx_Aa_to_Aa_3D
  SUBROUTINE ddyy_Aa_to_Aa_3D( grid, d_Aa, dyy_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: its yy-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dyy_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      dyy_Aa( k,j,i) = (d_Aa( k,j+1,i) + d_Aa( k,j-1,i) - 2._dp * d_Aa( k,j,i)) / grid%dx**2
    END DO
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    dyy_Aa( :,1      ,grid%i1:grid%i2) = (d_Aa( :,3      ,grid%i1:grid%i2) + d_Aa( :,1        ,grid%i1:grid%i2) - 2._dp * d_Aa( :,2        ,grid%i1:grid%i2)) / grid%dx**2
    dyy_Aa( :,grid%ny,grid%i1:grid%i2) = (d_Aa( :,grid%ny,grid%i1:grid%i2) + d_Aa( :,grid%ny-2,grid%i1:grid%i2) - 2._dp * d_Aa( :,grid%ny-1,grid%i1:grid%i2)) / grid%dx**2
    CALL sync
    
  END SUBROUTINE ddyy_Aa_to_Aa_3D
  SUBROUTINE ddxy_Aa_to_Aa_3D( grid, d_Aa, dxy_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: its xy-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dxy_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      dxy_Aa( k,j,i) = (d_Aa( k,j+1,i+1) + d_Aa( k,j-1,i-1) - d_Aa( k,j+1,i-1) - d_Aa( k,j-1,i+1)) / (4._dp * grid%dx * grid%dx)
    END DO
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    ! NO IDEA HOW TO DO THIS...
    dxy_Aa( :,1              ,grid%i1:grid%i2) = 0._dp
    dxy_Aa( :,grid%ny        ,grid%i1:grid%i2) = 0._dp
    dxy_Aa( :,grid%j1:grid%j2,1              ) = 0._dp
    dxy_Aa( :,grid%j1:grid%j2,grid%nx        ) = 0._dp
    CALL sync
    
  END SUBROUTINE ddxy_Aa_to_Aa_3D
  ! 3D upwind, for thermodynamics
  SUBROUTINE ddx_Aa_to_Aa_3D_upwind( grid, d_Aa, dx_Aa, U_3D_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Aa grid, using upwind one-sided differencing
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dx_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: U_3D_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Upwind one-sided differencing
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      IF (U_3D_Aa( k,j,i) > 0._dp) THEN
        dx_Aa( k,j,i) = (d_Aa( k,j,i  ) - d_Aa( k,j,i-1)) / grid%dx
      ELSE
        dx_Aa( k,j,i) = (d_Aa( k,j,i+1) - d_Aa( k,j,i  )) / grid%dx
      END IF
    END DO
    END DO
    END DO
    
    dx_Aa( :,grid%j1:grid%j2,1      ) = 0._dp
    dx_Aa( :,grid%j1:grid%j2,grid%nx) = 0._dp
    CALL sync
    
  END SUBROUTINE ddx_Aa_to_Aa_3D_upwind
  SUBROUTINE ddy_Aa_to_Aa_3D_upwind( grid, d_Aa, dy_Aa, V_3D_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Aa grid, using upwind one-sided differencing
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dy_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: V_3D_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Upwind one-sided differencing
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      IF (V_3D_Aa( k,j,i) > 0._dp) THEN
        dy_Aa( k,j,i) = (d_Aa( k,j  ,i) - d_Aa( k,j-1,i)) / grid%dx
      ELSE
        dy_Aa( k,j,i) = (d_Aa( k,j+1,i) - d_Aa( k,j  ,i)) / grid%dx
      END IF
    END DO
    END DO
    END DO
    
    dy_Aa( :,1      ,grid%i1:grid%i2) = 0._dp
    dy_Aa( :,grid%ny,grid%i1:grid%i2) = 0._dp
    CALL sync
    
  END SUBROUTINE ddy_Aa_to_Aa_3D_upwind
  
  ! Aa to Acx/Acy
  
  ! 2D
  SUBROUTINE ddx_Aa_to_Acx_2D( grid, d_Aa, dx_Acx)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: dx_Acx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      dx_Acx( j,i) = (d_Aa( j,i+1) - d_Aa( j,i)) / grid%dx
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddx_Aa_to_Acx_2D
  SUBROUTINE ddy_Aa_to_Acy_2D( grid, d_Aa, dy_Acy)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: dy_Acy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      dy_Acy( j,i) = (d_Aa( j+1,i) - d_Aa( j,i)) / grid%dx
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddy_Aa_to_Acy_2D
  SUBROUTINE ddx_Aa_to_Acy_2D( grid, d_Aa, dx_Acy)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: dx_Acy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      dx_Acy( j,i) = (d_Aa( j,i+1) + d_Aa( j+1,i+1) - d_Aa( j,i-1) - d_Aa( j+1,i-1)) / (4._dp * grid%dx)
    END DO
    END DO
    
    ! One-sided differencing on the boundary
    DO j = grid%j1, MIN(grid%ny-1,grid%j2)
      dx_Acy( j,1      ) = (d_Aa( j,2      ) + d_Aa( j+1,2      ) - d_Aa( j,1        ) - d_Aa( j+1,1        )) / (2._dp * grid%dx)
      dx_Acy( j,grid%nx) = (d_Aa( j,grid%nx) + d_Aa( j+1,grid%nx) - d_Aa( j,grid%nx-1) - d_Aa( j+1,grid%nx-1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddx_Aa_to_Acy_2D
  SUBROUTINE ddy_Aa_to_Acx_2D( grid, d_Aa, dy_Acx)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: dy_Acx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dy_Acx( j,i) = (d_Aa( j+1,i) + d_Aa( j+1,i+1) - d_Aa( j-1,i) - d_Aa( j-1,i+1)) / (4._dp * grid%dx)
    END DO
    END DO
    
    ! One-sided differencing on the boundary
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      dy_Acx( 1      ,i) = (d_Aa( 2,      i) + d_Aa( 2,      i+1) - d_Aa( 1,        i) - d_Aa( 1,        i+1)) / (2._dp * grid%dx)
      dy_Acx( grid%ny,i) = (d_Aa( grid%ny,i) + d_Aa( grid%ny,i+1) - d_Aa( grid%ny-1,i) - d_Aa( grid%ny-1,i+1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddy_Aa_to_Acx_2D
  ! 3D
  SUBROUTINE ddx_Aa_to_Acx_3D( grid, d_Aa, dx_Acx)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), INTENT(OUT)   :: dx_Acx
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      dx_Acx( k,j,i) = (d_Aa( k,j,i+1) - d_Aa( k,j,i)) / grid%dx
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddx_Aa_to_Acx_3D
  SUBROUTINE ddy_Aa_to_Acy_3D( grid, d_Aa, dy_Acy)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), INTENT(OUT)   :: dy_Acy
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
    DO k = 1, C%nZ
      dy_Acy( k,j,i) = (d_Aa( k,j+1,i) - d_Aa( k,j,i)) / grid%dx
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddy_Aa_to_Acy_3D
  SUBROUTINE ddx_Aa_to_Acy_3D( grid, d_Aa, dx_Acy)
    ! Input:  scalar on the Aa grid
    ! Output: its x-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), INTENT(OUT)   :: dx_Acy
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
    DO k = 1, C%nZ
      dx_Acy( k,j,i) = (d_Aa( k,j,i+1) + d_Aa( k,j+1,i+1) - d_Aa( k,j,i-1) - d_Aa( k,j+1,i-1)) / (4._dp * grid%dx)
    END DO
    END DO
    END DO
    
    ! One-sided differencing on the boundary
    DO j = grid%j1, MIN(grid%ny-1,grid%j2)
      dx_Acy( :,j,1      ) = (d_Aa( :,j,2      ) + d_Aa( :,j+1,2      ) - d_Aa( :,j,1        ) - d_Aa( :,j+1,1        )) / (2._dp * grid%dx)
      dx_Acy( :,j,grid%nx) = (d_Aa( :,j,grid%nx) + d_Aa( :,j+1,grid%nx) - d_Aa( :,j,grid%nx-1) - d_Aa( :,j+1,grid%nx-1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddx_Aa_to_Acy_3D
  SUBROUTINE ddy_Aa_to_Acx_3D( grid, d_Aa, dy_Acx)
    ! Input:  scalar on the Aa grid
    ! Output: its y-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), INTENT(OUT)   :: dy_Acx
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    ! Central differencing in the interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      dy_Acx( k,j,i) = (d_Aa( k,j+1,i) + d_Aa( k,j+1,i+1) - d_Aa( k,j-1,i) - d_Aa( k,j-1,i+1)) / (4._dp * grid%dx)
    END DO
    END DO
    END DO
    
    ! One-sided differencing on the boundary
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      dy_Acx( :,1      ,i) = (d_Aa( :,2,      i) + d_Aa( :,2,      i+1) - d_Aa( :,1,        i) - d_Aa( :,1,        i+1)) / (2._dp * grid%dx)
      dy_Acx( :,grid%nx,i) = (d_Aa( :,grid%nx,i) + d_Aa( :,grid%nx,i+1) - d_Aa( :,grid%nx-1,i) - d_Aa( :,grid%nx-1,i+1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddy_Aa_to_Acx_3D
  
  ! Acx/Acy to Aa
  
  ! 2D
  SUBROUTINE ddx_Acx_to_Aa_2D( grid, d_Acx, dx_Aa)
    ! Input:  scalar on the Acx grid
    ! Output: its x-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dx_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      dx_Aa( j,i) = (d_Acx( j,i) - d_Acx( j,i-1)) / grid%dx
    END DO
    END DO
    
    dx_Aa( grid%j1:grid%j2,1      ) = dx_Aa( grid%j1:grid%j2,2        )
    dx_Aa( grid%j1:grid%j2,grid%nx) = dx_Aa( grid%j1:grid%j2,grid%nx-1)
    CALL sync
    
  END SUBROUTINE ddx_Acx_to_Aa_2D
  SUBROUTINE ddy_Acy_to_Aa_2D( grid, d_Acy, dy_Aa)
    ! Input:  scalar on the Acy grid
    ! Output: its y-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dy_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
      dy_Aa( j,i) = (d_Acy( j,i) - d_Acy( j-1,i)) / grid%dx
    END DO
    END DO
    
    dy_Aa( 1      ,grid%i1:grid%i2) = dy_Aa( 2        ,grid%i1:grid%i2)
    dy_Aa( grid%ny,grid%i1:grid%i2) = dy_Aa( grid%ny-1,grid%i1:grid%i2)
    CALL sync
    
  END SUBROUTINE ddy_Acy_to_Aa_2D
  SUBROUTINE ddy_Acx_to_Aa_2D( grid, d_Acx, dy_Aa)
    ! Input:  scalar on the Acx grid
    ! Output: its y-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dy_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dy_Aa( j,i) = (d_Acx( j+1,i-1) + d_Acx( j+1,i) - d_Acx( j-1,i-1) - d_Acx( j-1,i)) / (4._dp * grid%dx)
    END DO
    END DO
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
      ! South ex. corners
      j = 1
      dy_Aa( j,i) = (d_Acx( j+1,i-1) + d_Acx( j+1,i) - d_Acx( j  ,i-1) - d_Acx( j  ,i)) / (4._dp * grid%dx)
      ! North ex. corners
      j = grid%ny
      dy_Aa( j,i) = (d_Acx( j  ,i-1) + d_Acx( j  ,i) - d_Acx( j-1,i-1) - d_Acx( j-1,i)) / (4._dp * grid%dx)
    END DO
    
    DO j = MAX(2,grid%j1), MIN(grid%ny-1,grid%j2)
      ! West ex. corners
      i = 1
      dy_Aa( j,i) = (d_Acx( j+1,i  ) - d_Acx( j-1,i  )) / (2._dp * grid%dx)
      ! East ex. corners
      i = grid%nx
      dy_Aa( j,i) = (d_Acx( j+1,i-1) - d_Acx( j-1,i-1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
    ! Corners
    IF (par%master) THEN
    dy_Aa( 1      ,1      ) = (d_Acx( 2      ,1        ) - d_Acx( 1        ,1        )) / grid%dx
    dy_Aa( 1      ,grid%nx) = (d_Acx( 2      ,grid%nx-1) - d_Acx( 1        ,grid%nx-1)) / grid%dx
    dy_Aa( grid%ny,1      ) = (d_Acx( grid%ny,1        ) - d_Acx( grid%ny-1,1        )) / grid%dx
    dy_Aa( grid%ny,grid%nx) = (d_Acx( grid%ny,grid%nx-1) - d_Acx( grid%ny-1,grid%nx-1)) / grid%dx
    END IF
    CALL sync
    
  END SUBROUTINE ddy_Acx_to_Aa_2D
  SUBROUTINE ddx_Acy_to_Aa_2D( grid, d_Acy, dx_Aa)
    ! Input:  scalar on the Acy grid
    ! Output: its x-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: dx_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dx_Aa( j,i) = (d_Acy( j-1,i+1) + d_Acy( j,i+1) - d_Acy( j-1,i-1) - d_Acy( j,i-1)) / (4._dp * grid%dx)
    END DO
    END DO
    
    DO j = MAX(2,grid%j1), MIN(grid%ny-1,grid%j2)
      ! West ex. corners
      i = 1
      dx_Aa( j,i) = (d_Acy( j-1,i+1) + d_Acy( j  ,i+1) - d_Acy( j-1,i  ) - d_Acy( j  ,i  )) / (4._dp * grid%dx)
      ! East ex. corners
      i = grid%nx
      dx_Aa( j,i) = (d_Acy( j-1,i  ) + d_Acy( j  ,i  ) - d_Acy( j-1,i-1) - d_Acy( j  ,i-1)) / (4._dp * grid%dx)
    END DO
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
      ! South ex. corners
      j = 1
      dx_Aa( j,i) = (d_Acy( j  ,i+1) - d_Acy( j  ,i-1)) / (2._dp * grid%dx)
      ! North ex. corners
      j = grid%ny
      dx_Aa( j,i) = (d_Acy( j-1,i+1) - d_Acy( j-1,i-1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
    ! Corners
    IF (par%master) THEN
    dx_Aa( 1      ,      1) = (d_Acy( 1        ,2      ) - d_Acy( 1        ,1        )) / grid%dx
    dx_Aa( 1      ,grid%nx) = (d_Acy( 1        ,grid%nx) - d_Acy( 1        ,grid%nx-1)) / grid%dx
    dx_Aa( grid%ny,1      ) = (d_Acy( grid%ny-1,2      ) - d_Acy( grid%ny-1,1        )) / grid%dx
    dx_Aa( grid%ny,grid%nx) = (d_Acy( grid%ny-1,grid%nx) - d_Acy( grid%ny-1,grid%nx-1)) / grid%dx
    END IF
    CALL sync
    
  END SUBROUTINE ddx_Acy_to_Aa_2D
  ! 3D
  SUBROUTINE ddx_Acx_to_Aa_3D( grid, d_Acx, dx_Aa)
    ! Input:  scalar on the Acx grid
    ! Output: its x-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dx_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      dx_Aa( k,j,i) = (d_Acx( k,j,i) - d_Acx( k,j,i-1)) / grid%dx
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddx_Acx_to_Aa_3D
  SUBROUTINE ddy_Acy_to_Aa_3D( grid, d_Acy, dy_Aa)
    ! Input:  scalar on the Acy grid
    ! Output: its y-derivative on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: dy_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      dy_Aa( k,j,i) = (d_Acy( k,j,i) - d_Acy( k,j-1,i)) / grid%dx
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddy_Acy_to_Aa_3D
  
  ! Acx/Acy to Acx/Acy
  SUBROUTINE ddx_Acx_to_Acx_2D( grid, d_Acx, dx_Acx)
    ! Input:  scalar on the Acx grid
    ! Output: its x-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: dx_Acx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = 2, grid%nx-2
    DO j = 1, grid%ny
      dx_Acx( j,i) = (d_Acx( j,i+1) - d_Acx( j,i-1)) / (2 * grid%dx)
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    dx_Acx( grid%j1:grid%j2,1        ) = (d_Acx( grid%j1:grid%j2,2        ) - d_Acx( grid%j1:grid%j2,1        )) / grid%dx
    dx_Acx( grid%j1:grid%j2,grid%nx-1) = (d_Acx( grid%j1:grid%j2,grid%nx-1) - d_Acx( grid%j1:grid%j2,grid%nx-2)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddx_Acx_to_Acx_2D
  SUBROUTINE ddy_Acx_to_Acx_2D( grid, d_Acx, dy_Acx)
    ! Input:  scalar on the Acx grid
    ! Output: its y-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: dy_Acx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dy_Acx( j,i) = (d_Acx( j+1,i) - d_Acx( j-1,i)) / (2 * grid%dx)
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    dy_Acx( 1      ,grid%i1:grid%i2) = (d_Acx( 2      ,grid%i1:grid%i2) - d_Acx( 1        ,grid%i1:grid%i2)) / grid%dx
    dy_Acx( grid%ny,grid%i1:grid%i2) = (d_Acx( grid%ny,grid%i1:grid%i2) - d_Acx( grid%ny-1,grid%i1:grid%i2)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddy_Acx_to_Acx_2D
  SUBROUTINE ddx_Acy_to_Acy_2D( grid, d_Acy, dx_Acy)
    ! Input:  scalar on the Acy grid
    ! Output: its x-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: dx_Acy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      dx_Acy( j,i) = (d_Acy( j,i+1) - d_Acy( j,i-1)) / (2 * grid%dx)
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    dx_Acy( grid%j1:grid%j2,1      ) = (d_Acy( grid%j1:grid%j2,2      ) - d_Acy( grid%j1:grid%j2,1        )) / grid%dx
    dx_Acy( grid%j1:grid%j2,grid%nx) = (d_Acy( grid%j1:grid%j2,grid%nx) - d_Acy( grid%j1:grid%j2,grid%nx-1)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddx_Acy_to_Acy_2D
  SUBROUTINE ddy_Acy_to_Acy_2D( grid, d_Acy, dy_Acy)
    ! Input:  scalar on the Acy grid
    ! Output: its y-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: dy_Acy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Central differencing in the interior
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-2
      dy_Acy( j,i) = (d_Acy( j+1,i) - d_Acy( j-1,i)) / (2 * grid%dx)
    END DO
    END DO
    
    ! One-sided differencing on the boundaries
    dy_Acy( 1        ,grid%i1:grid%i2) = (d_Acy( 2        ,grid%i1:grid%i2) - d_Acy( 1        ,grid%i1:grid%i2)) / grid%dx
    dy_Acy( grid%ny-1,grid%i1:grid%i2) = (d_Acy( grid%ny-1,grid%i1:grid%i2) - d_Acy( grid%ny-2,grid%i1:grid%i2)) / grid%dx
    CALL sync
    
  END SUBROUTINE ddy_Acy_to_Acy_2D
  SUBROUTINE ddx_Acx_to_Acy_2D( grid, d_Acx, dx_Acy)
    ! Input:  scalar on the Acx grid
    ! Output: its x-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: dx_Acy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      dx_Acy( j,i) = (d_Acx( j,i  ) + d_Acx( j+1,i  ) - d_Acx( j,i-1) - d_Acx( j+1,i-1)) / (2._dp * grid%dx)
    END DO
    END DO
    
    ! Boundaries
    DO j = grid%j1, MIN(grid%ny-1,grid%j2)
      ! West
      i = 1
      dx_Acy( j,i) = (d_Acx( j,i+1) + d_Acx( j+1,i+1) - d_Acx( j,i  ) - d_Acx( j+1,i  )) / (2._dp * grid%dx)
      ! East
      i = grid%nx
      dx_Acy( j,i) = (d_Acx( j,i-1) + d_Acx( j+1,i-1) - d_Acx( j,i-2) - d_Acx( j+1,i-2)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddx_Acx_to_Acy_2D
  SUBROUTINE ddy_Acx_to_Acy_2D( grid, d_Acx, dy_Acy)
    ! Input:  scalar on the Acx grid
    ! Output: its y-derivative on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: dy_Acy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      dy_Acy( j,i) = (d_Acx( j+1,i-1) + d_Acx( j+1,i  ) - d_Acx( j  ,i-1) - d_Acx( j  ,i  )) / (2._dp * grid%dx)
    END DO
    END DO
    
    ! Boundaries
    DO j = grid%j1, MIN(grid%ny-1,grid%j2)
      ! West
      i = 1
      dy_Acy( j,i) = (d_Acx( j,i  ) - d_Acx( j,i  )) / grid%dx
      ! East
      i = grid%nx
      dy_Acy( j,i) = (d_Acx( j,i-1) - d_Acx( j,i-1)) / grid%dx
    END DO
    CALL sync
    
  END SUBROUTINE ddy_Acx_to_Acy_2D
  SUBROUTINE ddx_Acy_to_Acx_2D( grid, d_Acy, dx_Acx)
    ! Input:  scalar on the Acy grid
    ! Output: its x-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: dx_Acx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dx_Acx( j,i) = (d_Acy( j-1,i+1) + d_Acy( j  ,i+1) - d_Acy( j-1,i  ) - d_Acy( j  ,i  )) / (2._dp * grid%dx)
    END DO
    END DO
    
    ! Boundaries
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      ! South
      j = 1
      dx_Acx( j,i) = (d_Acy( j  ,i+1) - d_Acy( j  ,i  )) / grid%dx
      ! North
      j = grid%ny
      dx_Acx( j,i) = (d_Acy( j-1,i+1) - d_Acy( j-1,i  )) / grid%dx
    END DO
    CALL sync
    
  END SUBROUTINE ddx_Acy_to_Acx_2D
  SUBROUTINE ddy_Acy_to_Acx_2D( grid, d_Acy, dy_Acx)
    ! Input:  scalar on the Acy grid
    ! Output: its y-derivative on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: dy_Acx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      dy_Acx( j,i) = (d_Acy( j  ,i  ) + d_Acy( j  ,i+1) - d_Acy( j-1,i  ) - d_Acy( j-1,i+1)) / (2._dp * grid%dx)
    END DO
    END DO
    
    ! Boundaries
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      ! South
      j = 1
      dy_Acx( j,i) = (d_Acy( j+1,i  ) + d_Acy(j+1,i+1) - d_Acy( j  ,i  ) - d_Acy( j  ,i+1)) / (2._dp * grid%dx)
      ! North
      j = grid%ny
      dy_Acx( j,i) = (d_Acy( j-1,i  ) + d_Acy(j-1,i+1) - d_Acy( j-2,i  ) - d_Acy( j-2,i+1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddy_Acy_to_Acx_2D
  
  ! Acx to Ab
  SUBROUTINE ddx_Acx_to_Ab_2D( grid, d_Acx, dx_Ab)
    ! Input:  scalar on the Acx grid
    ! Output: its x-derivative on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: dx_Ab
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = 2, grid%nx-2
    DO j = 1, grid%ny-1
      dx_Ab( j,i) = (d_Acx( j+1,i+1) + d_Acx( j  ,i+1) - d_Acx( j+1,i-1) - d_Acx( j  ,i-1)) / (4._dp * grid%dx)
    END DO
    END DO
    
    ! Boundaries
    DO j = grid%j1, MIN(grid%ny-1,grid%j2)
      i = 1
      dx_Ab( j,i) = (d_Acx( j+1,i+1) + d_Acx( j  ,i+1) - d_Acx( j+1,i  ) - d_Acx( j  ,i  )) / (2._dp * grid%dx)
      i = grid%nx-1
      dx_Ab( j,i) = (d_Acx( j+1,i  ) + d_Acx( j  ,i  ) - d_Acx( j+1,i-1) - d_Acx( j  ,i-1)) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddx_Acx_to_Ab_2D
  SUBROUTINE ddy_Acy_to_Ab_2D( grid, d_Acy, dy_Ab)
    ! Input:  scalar on the Acy grid
    ! Output: its y-derivative on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: dy_Ab
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = grid%i1, grid%i2-2
    DO j = 2, grid%ny-2
      dy_Ab( j,i) = (d_Acy( j+1,i+1) + d_Acy( j+1,i  ) - d_Acy( j-1,i+1) - d_Acy( j-1,i  )) / (4._dp * grid%dx)
    END DO
    END DO
    
    ! Boundaries
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      j = 1
      dy_Ab( j,i) = (d_Acy( j+1,i+1) + d_Acy( j+1,i  ) - d_Acy( j  ,i+1) - d_Acy( j  ,i  )) / (2._dp * grid%dx)
      j = grid%ny-1
      dy_Ab( j,i) = (d_Acy( j  ,i+1) + d_Acy( j  ,i  ) - d_Acy( j-1,i+1) - d_Acy( j-1,i  )) / (2._dp * grid%dx)
    END DO
    CALL sync
    
  END SUBROUTINE ddy_Acy_to_Ab_2D
  SUBROUTINE ddx_Acy_to_Ab_2D( grid, d_Acy, dx_Ab)
    ! Input:  scalar on the Acy grid
    ! Output: its x-derivative on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: dx_Ab
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      dx_Ab( j,i) = (d_Acy( j,i+1) - d_Acy( j,i)) / grid%dx
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddx_Acy_to_Ab_2D
  SUBROUTINE ddy_Acx_to_Ab_2D( grid, d_Acx, dy_Ab)
    ! Input:  scalar on the Acx grid
    ! Output: its y-derivative on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: dy_Ab
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      dy_Ab( j,i) = (d_Acx( j+1,i) - d_Acx( j,i)) / grid%dx
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE ddy_Acx_to_Ab_2D
  
! =============================================
! ===== Mapping between (staggered) grids =====
! =============================================

  ! Aa to Acx/Acy
  
  ! 2D
  SUBROUTINE map_Aa_to_Acx_2D( grid, d_Aa, d_Acx)
    ! Input:  scalar on the Aa grid
    ! Output: the same on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: d_Acx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      d_Acx( j,i) = (d_Aa( j,i) + d_Aa( j,i+1)) / 2._dp
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_Aa_to_Acx_2D
  SUBROUTINE map_Aa_to_Acy_2D( grid, d_Aa, d_Acy)
    ! Input:  scalar on the Aa grid
    ! Output: the same on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: d_Acy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      d_Acy( j,i) = (d_Aa( j,i) + d_Aa( j+1,i)) / 2._dp
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_Aa_to_Acy_2D
  ! 3D
  SUBROUTINE map_Aa_to_Acx_3D( grid, d_Aa, d_Acx)
    ! Input:  scalar on the Aa grid
    ! Output: the same on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), INTENT(OUT)   :: d_Acx
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      d_Acx( k,j,i) = (d_Aa( k,j,i) + d_Aa( k,j,i+1)) / 2._dp
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_Aa_to_Acx_3D
  SUBROUTINE map_Aa_to_Acy_3D( grid, d_Aa, d_Acy)
    ! Input:  scalar on the Aa grid
    ! Output: the same on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), INTENT(OUT)   :: d_Acy
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
    DO k = 1, C%nZ
      d_Acy( k,j,i) = (d_Aa( k,j,i) + d_Aa( k,j+1,i)) / 2._dp
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_Aa_to_Acy_3D
  
  ! Acx/Acy to Aa
  
  ! 2D
  SUBROUTINE map_Acx_to_Aa_2D( grid, d_Acx, d_Aa)
    ! Input:  scalar on the Acx grid
    ! Output: the same on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: d_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      d_Aa( j,i) = (d_Acx( j,i-1) + d_Acx( j,i)) / 2._dp
    END DO
    END DO
    
    d_Aa( grid%j1:grid%j2,1      ) = d_Acx( grid%j1:grid%j2,1        )
    d_Aa( grid%j1:grid%j2,grid%nx) = d_Acx( grid%j1:grid%j2,grid%nx-1)
    CALL sync
    
  END SUBROUTINE map_Acx_to_Aa_2D
  SUBROUTINE map_Acy_to_Aa_2D( grid, d_Acy, d_Aa)
    ! Input:  scalar on the Acy grid
    ! Output: the same on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(OUT)   :: d_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
      d_Aa( j,i) = (d_Acy( j-1,i) + d_Acy( j,i)) / 2._dp
    END DO
    END DO
    
    d_Aa( 1      ,grid%i1:grid%i2) = d_Acy( 1        ,grid%i1:grid%i2)
    d_Aa( grid%ny,grid%i1:grid%i2) = d_Acy( grid%ny-1,grid%i1:grid%i2)
    CALL sync
    
  END SUBROUTINE map_Acy_to_Aa_2D
  ! 3D
  SUBROUTINE map_Acx_to_Aa_3D( grid, d_Acx, d_Aa)
    ! Input:  scalar on the Acx grid
    ! Output: the same on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: d_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    DO k = 1, C%nZ
      d_Aa( k,j,i) = (d_Acx( k,j,i-1) + d_Acx( k,j,i)) / 2._dp
    END DO
    END DO
    END DO
    
    d_Aa( :,grid%j1:grid%j2,1      ) = d_Acx( :,grid%j1:grid%j2,1        )
    d_Aa( :,grid%j1:grid%j2,grid%nx) = d_Acx( :,grid%j1:grid%j2,grid%nx-1)
    CALL sync
    
  END SUBROUTINE map_Acx_to_Aa_3D
  SUBROUTINE map_Acy_to_Aa_3D( grid, d_Acy, d_Aa)
    ! Input:  scalar on the Acy grid
    ! Output: the same on the Aa grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(OUT)   :: d_Aa
    
    ! Local variables:
    INTEGER                                                         :: i,j,k
    
    DO i = grid%i1, grid%i2
    DO j = 2, grid%ny-1
    DO k = 1, C%nZ
      d_Aa( k,j,i) = (d_Acy( k,j-1,i) + d_Acy( k,j,i)) / 2._dp
    END DO
    END DO
    END DO
    
    d_Aa( :,1      ,grid%i1:grid%i2) = d_Acy( :,1        ,grid%i1:grid%i2)
    d_Aa( :,grid%ny,grid%i1:grid%i2) = d_Acy( :,grid%ny-1,grid%i1:grid%i2)
    CALL sync
    
  END SUBROUTINE map_Acy_to_Aa_3D
  
  ! Acx/Acy to Acy/Acx
  SUBROUTINE map_Acx_to_Acy_2D( grid, d_Acx, d_Acy)
    ! Input:  scalar on the Acx grid
    ! Output: the same on the Acy grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(OUT)   :: d_Acy
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      d_Acy( j,i) = (d_Acx( j  ,i-1) + d_Acx( j  ,i  ) + d_Acx( j+1,i-1) + d_Acx( j+1,i  )) / 4._dp
    END DO
    END DO
    
    ! Boundaries
    DO j = grid%j1, MIN(grid%ny-1,grid%j2)
      i = 1
      d_Acy( j,i) = (d_Acx( j  ,i  ) + d_Acx( j+1,i  )) / 2._dp
      i = grid%nx
      d_Acy( j,i) = (d_Acx( j  ,i-1) + d_Acx( j+1,i-1)) / 2._dp
    END DO
    CALL sync
    
  END SUBROUTINE map_Acx_to_Acy_2D
  SUBROUTINE map_Acy_to_Acx_2D( grid, d_Acy, d_Acx)
    ! Input:  scalar on the Acy grid
    ! Output: the same on the Acx grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(OUT)   :: d_Acx
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    ! Interior
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      d_Acx( j,i) = (d_Acy( j-1,i  ) + d_Acy( j-1,i+1) + d_Acy( j  ,i  ) + d_Acy( j  ,i+1)) / 4._dp
    END DO
    END DO
    
    ! Boundaries
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      j = 1
      d_Acx( j,i) = (d_Acy( j  ,i  ) + d_Acy( j  ,i+1)) / 2._dp
      j = grid%ny
      d_Acx( j,i) = (d_Acy( j-1,i  ) + d_Acy( j-1,i+1)) / 2._dp
    END DO
    CALL sync
    
  END SUBROUTINE map_Acy_to_Acx_2D
  
  ! Aa to Ab
  SUBROUTINE map_Aa_to_Ab_2D( grid, d_Aa, d_Ab)
    ! Input:  scalar on the Aa grid
    ! Output: the same on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(IN)    :: d_Aa
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: d_Ab
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      d_Ab( j,i) = (d_Aa( j,i) + d_Aa( j,i+1) + d_Aa( j+1,i) + d_Aa( j+1,i+1)) / 4._dp
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_Aa_to_Ab_2D
  
  ! Acx/Acy to Ab
  SUBROUTINE map_Acx_to_Ab_2D( grid, d_Acx, d_Ab)
    ! Input:  scalar on the Acx grid
    ! Output: the same on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(IN)    :: d_Acx
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: d_Ab
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      d_Ab( j,i) = (d_Acx( j,i) + d_Acx( j+1,i)) / 2._dp
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_Acx_to_Ab_2D
  SUBROUTINE map_Acy_to_Ab_2D( grid, d_Acy, d_Ab)
    ! Input:  scalar on the Acy grid
    ! Output: the same on the Ab grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(IN)    :: d_Acy
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx-1), INTENT(OUT)   :: d_Ab
    
    ! Local variables:
    INTEGER                                                         :: i,j
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
      d_Ab( j,i) = (d_Acy( j,i) + d_Acy( j,i+1)) / 2._dp
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_Acy_to_Ab_2D
  
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
!    REAL(dp)                                                        :: inverse_Hi2

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      inverse_Hi  = 1._dp / ice%Hi_Aa( j,i)
!      inverse_Hi2 = inverse_Hi * inverse_Hi

      ice%dzeta_dz_Aa( j,i)  = -inverse_Hi
!      ice%dzeta_dxz(j,i) =  inverse_Hi2 * dHi_dx_Ab(j,i)
!      ice%dzeta_dyz(j,i) =  inverse_Hi2 * dHi_dy_Ab(j,i)
      DO k = 1, C%nZ
        ice%dzeta_dt_Aa(  k,j,i) = inverse_Hi * (ice%dHs_dt_Aa(  j,i) - C%zeta(k) * ice%dHi_dt_Aa(  j,i))
        ice%dzeta_dx_Aa(  k,j,i) = inverse_Hi * (ice%dHs_dx_Aa(  j,i) - C%zeta(k) * ice%dHi_dx_Aa(  j,i))
        ice%dzeta_dy_Aa(  k,j,i) = inverse_Hi * (ice%dHs_dy_Aa(  j,i) - C%zeta(k) * ice%dHi_dy_Aa(  j,i))
!        ice%dzeta_dxx_Aa( k,j,i) = inverse_Hi * (ice%dHs_dxx_Ab( j,i) - C%zeta(k) * ice%dHi_dxx_Ab( j,i) - 2._dp * ice%dzeta_dx( k,j,i) * dHi_dx_Ab( j,i))
!        ice%dzeta_dxy_Aa( k,j,i) = inverse_Hi * (ice%dHs_dxy_Ab( j,i) - C%zeta(k) * ice%dHi_dxy_Ab( j,i) -         ice%dzeta_dy( k,j,i) * dHi_dx_Ab( j,i) - ice%dzeta_dx( k,j,i) * dHi_dy_Ab( j,i))
!        ice%dzeta_dyy_Aa( k,j,i) = inverse_Hi * (ice%dHs_dyy_Ab( j,i) - C%zeta(k) * ice%dHi_dyy_Ab( j,i) - 2._dp * ice%dzeta_dy( k,j,i) * dHi_dy_Ab( j,i))
      END DO
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calculate_zeta_derivatives
  
! ======================================
! ===== Neumann boundary condition =====
! ======================================

  SUBROUTINE Neumann_BC_Aa_2D( grid, d_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: the same, with boundary values set to next-to-boundary values,
    !         so that d/dx=0 on the west&east boundaries, and d/dy-0 on north&south
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx  ), INTENT(INOUT) :: d_Aa
    
    d_Aa( 1              ,grid%i1:grid%i2) = d_Aa( 2              ,grid%i1:grid%i2)
    d_Aa( grid%ny        ,grid%i1:grid%i2) = d_Aa( grid%ny-1      ,grid%i1:grid%i2)
    CALL sync
    d_Aa( grid%j1:grid%j2,1              ) = d_Aa( grid%j1:grid%j2,2              )
    d_Aa( grid%j1:grid%j2,grid%nx        ) = d_Aa( grid%j1:grid%j2,grid%nx-1      )
    CALL sync
    
  END SUBROUTINE Neumann_BC_Aa_2D
  SUBROUTINE Neumann_BC_Aa_3D( grid, d_Aa)
    ! Input:  scalar on the Aa grid
    ! Output: the same, with boundary values set to next-to-boundary values,
    !         so that d/dx=0 on the west&east boundaries, and d/dy-0 on north&south
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION( C%nZ, grid%ny  , grid%nx  ), INTENT(INOUT) :: d_Aa
    
    d_Aa( :,1              ,grid%i1:grid%i2) = d_Aa( :,2              ,grid%i1:grid%i2)
    d_Aa( :,grid%ny        ,grid%i1:grid%i2) = d_Aa( :,grid%ny-1      ,grid%i1:grid%i2)
    CALL sync
    d_Aa( :,grid%j1:grid%j2,1              ) = d_Aa( :,grid%j1:grid%j2,2              )
    d_Aa( :,grid%j1:grid%j2,grid%nx        ) = d_Aa( :,grid%j1:grid%j2,grid%nx-1      )
    CALL sync
    
  END SUBROUTINE Neumann_BC_Aa_3D
  SUBROUTINE Neumann_BC_Acx_2D( grid, d_Acx)
    ! Input:  scalar on the Aa grid
    ! Output: the same, with boundary values set to next-to-boundary values,
    !         so that d/dx=0 on the west&east boundaries, and d/dy-0 on north&south
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny  , grid%nx-1), INTENT(INOUT) :: d_Acx
    
    d_Acx( 1              ,grid%i1:MIN(grid%nx-1,grid%i2)) = d_Acx( 2              ,grid%i1:MIN(grid%nx-1,grid%i2))
    d_Acx( grid%ny        ,grid%i1:MIN(grid%nx-1,grid%i2)) = d_Acx( grid%ny-1      ,grid%i1:MIN(grid%nx-1,grid%i2))
    CALL sync
    d_Acx( grid%j1:grid%j2,1              ) = d_Acx( grid%j1:grid%j2,2              )
    d_Acx( grid%j1:grid%j2,grid%nx        ) = d_Acx( grid%j1:grid%j2,grid%nx-1      )
    CALL sync
    
  END SUBROUTINE Neumann_BC_Acx_2D
  SUBROUTINE Neumann_BC_Acy_2D( grid, d_Acy)
    ! Input:  scalar on the Aa grid
    ! Output: the same, with boundary values set to next-to-boundary values,
    !         so that d/dx=0 on the west&east boundaries, and d/dy-0 on north&south
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                                  INTENT(IN)    :: grid
    REAL(dp), DIMENSION(       grid%ny-1, grid%nx  ), INTENT(INOUT) :: d_Acy
    
    d_Acy( 1              ,grid%i1:grid%i2) = d_Acy( 2              ,grid%i1:grid%i2)
    d_Acy( grid%ny        ,grid%i1:grid%i2) = d_Acy( grid%ny-1      ,grid%i1:grid%i2)
    CALL sync
    d_Acy( grid%j1:MIN(grid%ny-1,grid%j2),1              ) = d_Acy( grid%j1:MIN(grid%ny-1,grid%j2),2              )
    d_Acy( grid%j1:MIN(grid%ny-1,grid%j2),grid%nx        ) = d_Acy( grid%j1:MIN(grid%ny-1,grid%j2),grid%nx-1      )
    CALL sync
    
  END SUBROUTINE Neumann_BC_Acy_2D

END MODULE derivatives_and_grids_module



