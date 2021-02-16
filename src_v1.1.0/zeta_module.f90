MODULE zeta_module

  USE mpi
  USE parallel_module,                 ONLY: par, sync
  USE configuration_module,            ONLY: dp, C
  USE data_types_module,               ONLY: type_grid, type_ice_model
  
  IMPLICIT NONE
  
  TYPE type_zeta_parameters
    ! Some useful constants for the scaled vertical coordinate transformation
    
    REAL(dp), DIMENSION(:), ALLOCATABLE :: a_zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: b_zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: c_zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: a_zetazeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: b_zetazeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: c_zetazeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: a_zeta_minus
    REAL(dp), DIMENSION(:), ALLOCATABLE :: b_zeta_minus
    REAL(dp), DIMENSION(:), ALLOCATABLE :: z_zeta_minus
    
    REAL(dp) :: a_N
    REAL(dp) :: b_1
    REAL(dp) :: b_t
    REAL(dp) :: s_t
  
  END TYPE type_zeta_parameters
  
  TYPE( type_zeta_parameters), SAVE :: p_zeta
  
CONTAINS

  FUNCTION vertical_average(f) RESULT(average_f)
    ! Calculate the vertical average of any given function f defined at the vertical zeta grid.
    !  average_f = INTEGRAL_{zeta=0}^{zeta=1} (f(zeta))dzeta
    ! So the the integration is in the direction of the positive zeta-axis from C%zeta(k=1) = 0 up to C%zeta(k=C%nZ) = 1.
    ! Numerically: de average between layer k and k+1 is calculated and multiplied by the distance between those 
    ! layers k and k+1, which is immediately the weight factor for this contribution because de total layer distance 
    ! is scaled to 1. The sum of all the weighted contribution gives the vertical average average_f of f.    
    
    IMPLICIT NONE

    ! Input variable:
    REAL(dp), DIMENSION(C%nZ), INTENT(IN) :: f

    ! Result variable:
    REAL(dp)                              :: average_f
    
    ! Local variable:
    INTEGER                               :: k

    ! See equation (11.14):
    average_f = 0._dp
    DO k = 1, C%nZ-1
      average_f = average_f + 0.5_dp * (f(k+1) + f(k)) * (C%zeta(k+1) - C%zeta(k))
    END DO
    
  END FUNCTION vertical_average
  FUNCTION vertical_integrate(f) RESULT(int_f)
    ! This subroutine calcualtes the integral int_f:
    !  int_f(k) = INTEGRAL_bottom[C%zeta(k=C%nZ)=1]^zeta[C%zeta(k)] f(zeta) dzeta
    ! In case the integrant f is positive (our cases) our result must be negative because we
    ! integrate from C%zeta(k=C%nZ) = 1 opposite to the zeta-axis up to C%zeta(k). (Our dzeta's are
    ! negative).
    ! This subroutine returns the integral for each layer k from bottom to this layer, so
    ! inf_f is an array with length C%nZ: int_f(k=1:C%nZ)
    ! The value of the integrant f at some integration step k is the average of f(k+1) and f(k).
    !  int_f(k) = int_f(k+1) + 0.5*(f(k+1) + f(k))*(-dzeta)
    ! So for f > 0  int_f < 0.    
    
    IMPLICIT NONE

    ! Input variable:
    REAL(dp), DIMENSION(C%nZ), INTENT(IN) :: f

    ! Output variable:
    REAL(dp), DIMENSION(C%nZ)             :: int_f

    ! Local variable:
    INTEGER                               :: k

    int_f(C%nZ) = 0._dp
    DO k = C%nZ-1, 1, -1
      int_f(k) = int_f(k+1) - 0.5_dp * (f(k+1) + f(k)) * (C%zeta(k+1) - C%zeta(k))
    END DO
    
  END FUNCTION vertical_integrate  
  SUBROUTINE calculate_zeta_derivatives( grid, ice)
    ! This subroutine calculates the global struct which contains the derivatives of 
    ! zeta, which are used in the transformation to the t,x,y,zeta coordinates. 
    ! This zeta derivatives are the Jacobians.    
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: inverse_Hi                 ! Contains the inverse of Hi
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      inverse_Hi  = 1._dp / MAX(0.1_dp, ice%Hi_Aa( j,i))
      ice%dzeta_dz_Aa( j,i)  = -inverse_Hi                                                      
      DO k = 1, C%nZ
        ice%dzeta_dt_Aa( k,j,i)  =  inverse_Hi * (ice%dHs_dt_Aa( j,i)  - C%zeta(k) * ice%dHi_dt_Aa( j,i))                                   
        ice%dzeta_dx_Aa( k,j,i)  =  inverse_Hi * (ice%dHs_dx_Aa( j,i)  - C%zeta(k) * ice%dHi_dx_Aa( j,i))                                   
        ice%dzeta_dy_Aa( k,j,i)  =  inverse_Hi * (ice%dHs_dy_Aa( j,i)  - C%zeta(k) * ice%dHi_dy_Aa( j,i))                                   
      END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calculate_zeta_derivatives
  SUBROUTINE initialize_zeta_discretization   
    
    IMPLICIT NONE
    
    ! Local variables:
    INTEGER                       :: k
    REAL(dp), DIMENSION(2:C%nZ)   :: a_k
    REAL(dp), DIMENSION(1:C%nZ-1) :: b_k
    REAL(dp), DIMENSION(3:C%nZ)   :: c_k
    REAL(dp), DIMENSION(1:C%nZ-2) :: d_k

    ALLOCATE(p_zeta%a_zeta(2:C%nZ-1))
    ALLOCATE(p_zeta%b_zeta(2:C%nZ-1))
    ALLOCATE(p_zeta%c_zeta(2:C%nZ-1))
    ALLOCATE(p_zeta%a_zetazeta(2:C%nZ-1))
    ALLOCATE(p_zeta%b_zetazeta(2:C%nZ-1))
    ALLOCATE(p_zeta%c_zetazeta(2:C%nZ-1))
    
    ALLOCATE(p_zeta%z_zeta_minus(3:C%nZ))
    ALLOCATE(p_zeta%a_zeta_minus(3:C%nZ))
    ALLOCATE(p_zeta%b_zeta_minus(3:C%nZ))
    
    DO k = 2, C%nZ
     a_k(k) = C%zeta(k)   - C%zeta(k-1)
    END DO
    DO k = 1, C%nZ-1
     b_k(k) = C%zeta(k+1) - C%zeta(k)
    END DO
    DO k = 3, C%nZ
     c_k(k) = C%zeta(k)   - C%zeta(k-2)
    END DO
    DO k = 1, C%nZ-2
     d_k(k) = C%zeta(k+2) - C%zeta(k)
    END DO
    p_zeta%a_N = a_k(C%nZ)
    p_zeta%b_1 = b_k(1)
    
    DO k = 2, C%nZ-1
     p_zeta%a_zeta(k)     =         - b_k(k)  / (a_k(k) * (a_k(k) + b_k(k))) 
     p_zeta%b_zeta(k)     = (b_k(k) - a_k(k)) / (a_k(k) *  b_k(k))  
     p_zeta%c_zeta(k)     =           a_k(k)  / (b_k(k) * (a_k(k) + b_k(k)))
     p_zeta%a_zetazeta(k) =           2._dp   / (a_k(k) * (a_k(k) + b_k(k)))
     p_zeta%b_zetazeta(k) =         - 2._dp   / (a_k(k) *  b_k(k))       
     p_zeta%c_zetazeta(k) =           2._dp   / (b_k(k) * (a_k(k) + b_k(k)))
    END DO

    ! Not all of these are in use:
    DO k = 3, C%nZ
      p_zeta%z_zeta_minus(k) =           a_k(k)  / (c_k(k) * (c_k(k) - a_k(k)))   
      p_zeta%a_zeta_minus(k) =           c_k(k)  / (a_k(k) * (a_k(k) - c_k(k)))   
      p_zeta%b_zeta_minus(k) = (a_k(k) + c_k(k)) / (a_k(k) * c_k(k))          
    END DO
    
    ! Discretization coefficients in time
    p_zeta%b_t = -1._dp / C%dt_thermo                      ! See equations 14.38
    p_zeta%s_t =  1._dp / C%dt_thermo                      ! See equations 14.39
    
  END SUBROUTINE initialize_zeta_discretization

END MODULE zeta_module
