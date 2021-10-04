  MODULE common_taboo

  USE configuration_module,  ONLY: dp
  USE FMZM
  IMPLICIT NONE
!
! This module declares variables common to all of the three TABOO Tasks 
! Revised
!
 INTEGER, PARAMETER  :: i4b = SELECTED_INT_KIND(9)
 INTEGER, PARAMETER  :: sp  = KIND(1.0)
 INTEGER, PARAMETER  :: qp  = KIND(1.0Q0)
! SAVE   ! ???
!
 Real(dp), parameter :: raggio   = 6.371D6     ! Earth radius (m)
 Real(dp), parameter :: emass    = 5.97D24     ! Earth mass (kg)
 Real(dp), parameter :: rhoea    = 5.51157D3   ! Average density (kg/m3)
 Real(dp), parameter :: rhoice   = 1.D3        ! Ice density (kg/m3)
!
 Integer(i4b), parameter :: nv_max     =  24   ! maximum allowed number of v.e. layers 
 Integer(i4b), parameter :: nroots_max =  96   ! maximum allowed number of modes  
 Integer(i4b), parameter :: llmin      =   0   ! minimum allowed degree 
 Integer(i4b), parameter :: llmax      = 256   ! maximum allowed degree 
 Integer     :: vec (llmin:llmax,nroots_max)   ! Marker for fake modes 
! 
 Integer(i4b) :: IV           ! Verbose (1) or Silent (0) mode
 Integer(i4b) :: lmin, lmax2     ! Current values of lmin e lmax, taken from TASK_*.DAT    (LMAX is stored in C_S
 Integer(i4b) :: NROOTS         ! NROOTS is the number of viscoelastic modes 
 Integer(i4b) :: Only_Elastic   ! Switch for elastic fields 
!
 TYPE(FM)  :: r  (0:nv_max+1)	 ! Radii 
 TYPE(FM)  :: rho(0:nv_max+1)	 ! Density 
 TYPE(FM)  :: rmu(0:nv_max+1)	 ! Shear Moduli 
 TYPE(FM)  :: vis(0:nv_max+1)	 ! Viscosity
 TYPE(FM)  :: pon(0:nv_max+1)	 ! Gravity  
!
 REAL(sp)  :: s_vis(0:nv_max+1)	 ! Single precision Viscosity 
!
 TYPE(FM)  :: h_e(llmin:llmax), l_e(llmin:llmax), k_e(llmin:llmax)  ! Elastic modes
 TYPE(FM)  :: h_f(llmin:llmax), l_f(llmin:llmax), k_f(llmin:llmax)  ! Fluid modes 
!
 type(fm) :: h_v (llmin:llmax,1:nroots_max) ! Viscoelastic residues h
 type(fm) :: l_v (llmin:llmax,1:nroots_max) !	    "	       "    l 
 type(fm) :: k_v (llmin:llmax,1:nroots_max) !	    "	       "    k 
 type(fm) :: s (llmin:llmax,0:nroots_max)   ! Roots of the secular equation, in kyr**(-1)  
!
 Real(dp) :: r_h (llmin:llmax,1:nroots_max)   ! Normalized residue for h
 Real(dp) :: r_l (llmin:llmax,1:nroots_max)   ! Normalized residue for l
 Real(dp) :: r_k (llmin:llmax,1:nroots_max)   ! Normalized residue for k
 Real(dp) :: tek (llmin:llmax,1:nroots_max)   ! Relaxation time in k-yrs
!
  END MODULE common_taboo
!
!
!
!
 MODULE COMMON_FOR_SPECTRUM
  USE configuration_module, ONLY: dp
!
! Variables in use by Sbr. 'spectrum' and in the routines depending on it  
! Revised GS May 2010 for g95 implementation 
!
 USE common_taboo
 IMPLICIT NONE
!
!
 type(fm) :: pivo	 ! Largest coeff. of the secular polynomial 
 type(fm) :: rubens	 ! For the elastic part of the solution 
 type(fm) :: ggg	 ! Newton constant
 type(fm) :: pi_c    ! pi for in taboo
 type(fm) :: xmass	 ! Earth mass 
!
 type(fm) :: aco (0:nv_max+1) ! <<A>> coefficients (bottom to top)  
 type(fm) :: g (0:nv_max+1)   ! gravity at the interfaces (bottom to top) 
!
 type(fm) :: a (6, 6), b (6, 6), c (6, 6), d (6, 6)  ! 6*6 propagators 
 type(fm) :: ac(6, 6), ad(6, 6), bc(6, 6), bd(6, 6)  ! Outputs of MATPROD 
 type(fm) :: k0 (6, 6), k1 (6, 6), k2 (6, 6)	      ! 6*6 propagators 
 type(fm) :: matela (6, 6)			      ! Elastic product 
 type(fm) :: sinist_2 (3, 6), sinist_1 (3, 6)        ! Left products 
 type(fm) :: rr (0:2 * nv_max, 3, 3)		      ! Matrix to be inverted 
 type(fm) :: co (0:3 * 2 * nv_max), aa(nroots_max + 1), op(nroots_max+1) ! Polynomial coefficints 
!
 type(fm) :: rad (llmin:llmax,nroots_max)  ! Roots in years 
 type(fm) :: cc (0:2, 1:nv_max, 6, 6)      ! Propagator 
 type(fm) :: coefmat (0:2 * nv_max, 6, 6)  ! Propagator 
 type(fm) :: ctmp (0:3 * 2 * nv_max - 1)   ! Derivative of the secular poly. 
 type(fm) :: rrrr (0:2 * nv_max, 3, 6), qqqq (0:2 * nv_max, 3, 6)  ! Left products 
 type(fm) :: cmb(6,3), bcs (3) 	    ! CMB and surface Boundary conditions 
 type(fm) :: qq (0:2 * nv_max, 3, 3)	    ! Q matrix 
 type(fm) :: r_r (3, 3, 0:nroots_max), q_q (3, 3, 0:nroots_max)  ! Ausilium di R e Q  
 type(fm) :: qr (3, 3, 0:nroots_max)	    ! Product Q*R 
 type(fm) :: raggiu (3, 3, 0:nroots_max)   ! Adjoint 
 type(fm) :: derpo (1:nroots_max)	    ! Derivative of the secular poly. 
!
 type(fm) :: rt1 (1:nroots_max), rt2 (1:nroots_max)  ! Real and Imag. parts of the roots, in kyears  
 type(fm) ::  x_el (3), xx (3), xr (3, nroots_max)   ! Solution in vector form 
 type(fm) :: ded			      ! Useful to compute the determinant
 EXTERNAL DED
!
 END MODULE COMMON_FOR_SPECTRUM


