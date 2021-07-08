! Modified by Bas de Boer 21 Januari 2014: combined some loops for (hopefully) speed up of computations
!
MODULE SELEN_spline_interpolation_module

  USE configuration_module, ONLY: dp, C
  IMPLICIT NONE

CONTAINS

!! SPLINE fit
SUBROUTINE spline(x,y,yp1,ypn,y2)
! Given arrays x(1:C%SELEN_irreg_time_n) and y(1:C%SELEN_irreg_time_n) containing a tabulated function, i.e., yi = f(xi), with
! x1 < x2 < ::: < xN, and given values yp1 and ypn for the rst derivative of the interpolating
! function at points 1 and C%SELEN_irreg_time_n, respectively, this routine returns an array y2(1:C%SELEN_irreg_time_n) of
! length C%SELEN_irreg_time_n which contains the second derivatives of the interpolating function at the tabulated
! points xi. If yp1 and/or ypn are equal to 10E30 or larger, the routine is signaled to set
! the corresponding boundary condition for a natural spline, with zero second derivative on
! that boundary.

  INTEGER                                                        :: i,k
  REAL(dp),   DIMENSION(             0:C%SELEN_irreg_time_n), INTENT(IN)  :: x
  COMPLEX*16, DIMENSION(  C%SELEN_jmax,0:C%SELEN_irreg_time_n), INTENT(IN)  :: y ! = S
  COMPLEX*16, DIMENSION(  C%SELEN_jmax              ), INTENT(IN)  :: yp1
  COMPLEX*16, DIMENSION(  C%SELEN_jmax              ), INTENT(IN)  :: ypn

  COMPLEX*16, DIMENSION(  C%SELEN_jmax,0:C%SELEN_irreg_time_n), INTENT(OUT) :: y2

  COMPLEX*16, DIMENSION(  C%SELEN_jmax,0:C%SELEN_irreg_time_n)              :: u
  COMPLEX*16, DIMENSION(  C%SELEN_jmax              )              :: p,un
  REAL(dp)                                                       :: sig,qn


  write (0,*) 'spline interpolation: n = ',C%SELEN_irreg_time_n

!if (yp1 == 0._dp) then                                               ! The lower boundary condition is set either to be
!  y2(1) = 0._dp                                                      ! natural
!  u(1)  = 0._dp
!else                                                                 ! or else to have a specied rst derivative.
!  y2(1) = -0.5_dp
!  u(1)  = (3._dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
!endif

   y2(:,0) = -0.5
   u(:,0)  = (3./(x(1)-x(0)))*((y(:,1)-y(:,0))/(x(1)-x(0))-yp1(:))

do i = 1,C%SELEN_irreg_time_n-1                                               ! This is the decomposition loop of the tridiagonal
  sig     = (x(i)-x(i-1))/(x(i+1)-x(i-1))                            ! algorithm. y2 and u are used for temporary
  p(:)    = sig * y2(:,i-1) + 2._dp                                  ! storage of the decomposed factors.
  y2(:,i) = (sig - 1.) / p(:)
  u(:,i)  = (6._dp*((y(:,i+1)-y(:,i))/(x(i+1)-x(i))-(y(:,i)-y(:,i-1)) / & 
            (x(i)-x(i-1)))/(x(i+1)-x(i-1)) - sig * u(:,i-1)) / p(:)
end do

!if (ypn == 0._dp) then                                              ! The upper boundary condition is set either to be
!  qn = 0._dp                                                        ! natural
!  un = 0._dp
!else                                                                ! or else to have a specied rst derivative.
!  qn = 0.5_dp
!  un = (3._dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
!end if

  qn    = 0.5
  un(:) = (3./(x(C%SELEN_irreg_time_n)-x(C%SELEN_irreg_time_n-1)))*(ypn(:) - &
          (y(:,C%SELEN_irreg_time_n)-y(:,C%SELEN_irreg_time_n-1))/(x(C%SELEN_irreg_time_n)-x(C%SELEN_irreg_time_n-1)))

  y2(:,C%SELEN_irreg_time_n) = (un(:)-qn*u(:,C%SELEN_irreg_time_n-1))/(qn*y2(:,C%SELEN_irreg_time_n-1)+1._dp)

do k = C%SELEN_irreg_time_n-1,1,-1                                           ! This is the backsubstitution loop of the tridiagonal algorithm.
  y2(:,k) = y2(:,k)*y2(:,k+1)+u(:,k)
end do

END SUBROUTINE



!! SPLINE Interpolation
SUBROUTINE splint(xa,ya,y2a,x,y)
! Given the arrays xa(1:C%SELEN_irreg_time_n) and ya(1:C%SELEN_irreg_time_n) of length C%SELEN_irreg_time_n, which tabulate a function 
! (with the xai 's in order), and given the array y2a(1:(with the), which is the output from spline above,
! and given a value of x, this routine returns a cubic-spline interpolated value y.

  REAL(dp),   DIMENSION(               C%SELEN_irreg_time_n), INTENT(IN)  :: xa   ! time points
  COMPLEX*16, DIMENSION(  C%SELEN_jmax,0:C%SELEN_irreg_time_n), INTENT(IN)  :: ya   ! S
  COMPLEX*16, DIMENSION(  C%SELEN_jmax,0:C%SELEN_irreg_time_n), INTENT(IN)  :: y2a  ! spline coordinates

  REAL(dp),                                          INTENT(IN)  :: x    ! time points for interpolation
  COMPLEX*16, DIMENSION(  C%SELEN_jmax              ), INTENT(OUT) :: y    ! interpolated function at time point x

  INTEGER, SAVE                       :: k1
  INTEGER                             :: k,khi,klo,i
  REAL(dp)                            :: a,b,h

! We will find the right place in the table by means of bisection.
! This is optimal if sequential calls to this routine are at random
! values of x. If sequential calls are in order, and closely
! spaced, one would do better to store previous values of
! klo and khi and test if they remain appropriate on the next call.

klo = 0
khi = C%SELEN_irreg_time_n

  do while(khi-klo > 1)
    k = (khi+klo) / 2
    if (xa(k) > x) then
      khi = k
    else
      klo = k
    end if
  end do

! save
!k1 = klo

!write (0,'(3f12.3,2i5)') x,xa(klo),xa(khi), klo,khi

h = xa(khi) - xa(klo)

!if (h.eq.0.)  pause 'bad xa input in splint'           ! The xa's must be distinct.

a = (xa(khi) - x) / h                                   ! Cubic spline polynomial is now evaluated.
b = (x - xa(klo)) / h

y(:) = a*ya(:,klo) + b*ya(:,khi) + ((a**3-a)*y2a(:,klo) + (b**3-b)*y2a(:,khi))*(h**2) / 6._dp

END SUBROUTINE

END MODULE SELEN_spline_interpolation_module
