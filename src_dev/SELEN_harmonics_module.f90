MODULE SELEN_harmonics_module

  USE configuration_module,              ONLY: dp
  USE data_types_module,                 ONLY: type_SELEN_global
  
  IMPLICIT NONE
!
! This is the include file "HARMONICS.F90"
!
! HARMONICS.F90 is a set of Fortran utilities for Spherical Harmonics (SH) 
! analysis, and some other useful stuff, listed below. 
!
! *** Modified by GS & FC 06-13-2008  = INTEL PORT = V 2.6 
! *** Romberg integration program added July 09, 2008. 
! *** On July 19, 2008 I have added to_real10 --- GS 
! *** Updated with new routines on August 2008, for v. 2.7
! *** Some typos corrected on September 29, 2008. 
! *** On october 2, I have moved here the pixelization routines 
! *** Reviewed GS & FC July 2009 -  "Varying coastlines" & ALMA coupling
! *** Reviewed GS & FC November 2009 - Porting under gfortran 
! *** Revised GS May 2010 - g95 - Ice breaker routine included 
!
!
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
! Copyright (C) 2008 Giorgio Spada, Florence Colleoni, and Paolo Stocchi 
!
! This file is part of SELEN. 
!  
! SELEN is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or at your option) any later 
! version. 
!
! SELEN is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details. 
! 
! You should have received a copy of the GNU General Public License along 
! with SELEN.  If not, see <http://www.gnu.org/licenses/>.
! +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! This file contains several units *mainly* related to spherical harmonics:  
! 
! - Subroutine PLegendre (*)          Legendre polynomials 
! - Subroutine PlmBar (*)             Normalized Legendre functions
! - Subroutine PLMBAR_MOD             (Modified) normalized Legendre functions
! - Subroutine Pleg                   Legendre polynomials              
! - Function J_INDEX                  J=l*(l+1)/2+m+1 
! - Function DOM                      Delta(0,m) 
! - Function MJ                       Retrives m (order) from J
! - Function LJ                       Retrives l (degree) from J
! - Subroutine HARMO                  Complex spherical harmonics (CSH)
! - Subroutine Scan_string            Scans for valid strings
! - Subroutine STOP_CONFIG            Stop configuration if errors are detected
! - Subroutine CHAR100_2_REAL         Converts from CHARACTER*100 to FLOATING POINT 
! - Function SINDD                    Sine with argument in degrees 
! - Function COSDD                    Cosine with argument in degrees 
! - SUBROUTINE READ_JUNK              Read "I" lines from Unit "J"
!
!   (*)  adapted from SHTOOLS, Copyright (c) 2005, Mark A. Wieczorek
!   (**) implemented by the ACM TOMS algorithm #TOMS351. The "early" algorithm 
!        used here is by John Burkardt (see http://people.scs.fsu.edu/~burkardt). 
!        No modifications have been made on the code available from http://
!        people.scs.fsu.edu/~burkardt/f77_src/toms351/toms351.html.  
! ----------------------------------------------------------------------------
! 
!
!
!
!
CONTAINS
    FUNCTION SINDD(ALFA)
    implicit NONE
!
! --- Assuming that ALFA is given in DEGREES, it returns the SINE of ALFA
!     *********************** GS and FC 11-02-2009 **********************
!
    REAL(dp),  PARAMETER :: PI=3.14159265358979323840 
    REAL(dp) SINDD, ALFA, ALFARAD              
    ALFARAD=ALFA*PI/180.
    SINDD=SIN(ALFARAD)

    END FUNCTION SINDD 
!
!
!
!
!
!
    FUNCTION COSDD(ALFA)
    implicit NONE
!
! --- Assuming that ALFA is given in DEGREES, it returns the COSINE of ALFA
!     ************************ GS and FC 11-02-2009 ***********************
!    
    REAL(dp),  PARAMETER :: PI=3.14159265358979323840 
    REAL(dp) COSDD, ALFA, ALFARAD               
    ALFARAD=ALFA*PI/180.    
    COSDD=COS(ALFARAD)    
!
    END FUNCTION COSDD
!
!
!
 SUBROUTINE SCAN_STRING (S, N, SS, M)
 IMPLICIT NONE 
 INTEGER I, J, N, M, APO(2*N) 
 CHARACTER*200 S ; CHARACTER*100 SS(N)
 CHARACTER*1, PARAMETER :: PRIME="'" 
!
! Given the input string s, this routine determines the sub-strings ss of s 
! which are delimited by consecutive pairs of primes ('). n (input) is the 
! expected number of substrings, and m (output) is the effective number 
! found. Apo is a pointer to the primes. The execution is stopped whenever 
! i) n=0, ii) m=0, iii) m/=n, or iv) if m is odd. The maximum lenght of the 
! substrings is 20, that of the input string is 200.  == GS Nov 17 2007 == 
!
! write(*,'(a100)') s
!
! ----- Exits for n=0
!
!
 if(n==0) then 
      Write(88,*) "SCAN_STRING has nothing to do" 
      Write(88,*) "The program will STOP ----------------"
      Write(*, *) "SCAN_STRING has nothing to do" 
      Write(*, *) "The program will STOP ----------------"     
          call stop_config 
      Stop
 endif
!
! ----- Looking for primes (', not ") within the input string  
 i=0
 do j=1, 200
    if(s(j:j)==prime) then 
        i=i+1 ; apo(i)=j
    endif
 enddo 
 m=i/2 
!
! ----- Exits if no primes are found 
 if(i==0) then 
      Write(*,*) "SCAN_STRING has found NO primes" 
      Write(*,*) "The program will STOP ----------------"
      Write(*, *) "SCAN_STRING has found NO primes" 
      Write(*,'(a100)') s 
      Write(*, *) "The program will STOP ----------------"     
          call stop_config 
      Stop
 endif 
!
! ----- Exits if an odd number of primes is found
 if(mod(i,2)/=0) then 
      Write(*,*) "SCAN_STRING has found an even number of primes:", i
      Write(*,*) "The program will STOP ----------------"
      Write(*, *) "SCAN_STRING has found an even number of primes:", i
      Write(*, *) "The program will STOP ----------------"     
          call stop_config 
      Stop
 endif
!
! ----- Exits if m/=n, otherwise determines the substrings
 if(m/=n) then 
      Write(*,*) "SCAN_STRING has found ", m, " substrings"
          Write(*,*) "SCAN_STRING  expected ", n, " substrings"
      Write(*,*) "The program will STOP ----------------"   
      Write(*, *) "SCAN_STRING has found ", m, " substrings"
          Write(*, *) "SCAN_STRING  expected ", n, " substrings"
      Write(*, *) "The program will STOP ----------------"
          call stop_config 
      Stop
 else
!
!        Write(*,*) "SCAN_STRING has found ", m, " substrings"
!        Write(*,*) "SCAN_STRING  expected ", n, " substrings"   
!
         do i=1, n 
            ss(i)=s(apo(2*i-1)+1:apo(2*i)-1)
         enddo
 endif  
! 
 end subroutine scan_string 
!
!
!
  Subroutine STOP_CONFIG 
!
! --- Creates a STOP file that stops the execution of makeseles.sh 
!     ******************* GS and FC 13-12-2007 *******************
!
    open (11,file='stop.dat',status='new')
    Write(11,*) "Errors have been detected - SELEN will stop"
    close(11)
!
    write(0 ,*) 'STOP_CONFIG: SELEN will STOP!' 
    write(88,*) 'STOP_CONFIG: SELEN will STOP!' 
    STOP 
!
    End subroutine STOP_CONFIG     
!
!
!
 SUBROUTINE CHAR100_2_REAL(STRING, REAL_VAL)
 IMPLICIT NONE
 INTEGER, PARAMETER :: IU=30
 CHARACTER*100 STRING 
 REAL REAL_VAL  
!
! Given the input CHARACTER*100 string STRING, this routine converts it in Floating 
! point value by copying it on a file and reading it back. *** GS Oct 19 2007 *** 
! Revised GS July 27, 2009. 
!
 open(iu,file='junk.dat',status='unknown') ; write(iu,*)  string ; close(iu) 
 open(iu,file='junk.dat',status='unknown') ; read(iu,*) real_val;  close(iu)  
!
 end SUBROUTINE CHAR100_2_REAL
!
!
!
subroutine PLegendre(p, lmax, z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    This subroutine evalutates all of the unnormalized Legendre polynomials 
!    up to degree lmax. 
!
!    Calling Parameters:
!        Out
!            p:    A vector of all unnormalized Legendgre polynomials evaluated at 
!                z up to lmax. The lenght must by greater or equal to (lmax+1).
!        IN
!            lmax:    Maximum degree to compute.
!            z:    Value within [-1, 1], cos(colatitude) or sin(latitude).
!
!    Notes:
!    
!    1.    The integral of Pl**2 over (-1,1) is 2/(2l+1).
!    2.    Values are calculated accoring to the following recursion scheme:
!            P_0(z) = 1.0, P_1(z) = z, and 
!            P_l(z) = (2l-1) * z * P_{l-1}(z) / l - (l-1) * P_{l-2}(z) / l
!
!    Dependencies:    None
!
!    Written by Mark Wieczorek June 2004
!
! ----> Modified to SINGLE PRECISION by Giorgio Spada 2007 
! ----> Also modified for the management of ERROR conditions 
!
!    Copyright (c) 2005, Mark A. Wieczorek
!    All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    integer, intent(in) ::    lmax
    REAL(dp), intent(out) ::    p(:)
           REAL(dp), intent(in) ::    z
           REAL(dp) ::    pm2, pm1, pl
          integer ::    l

    if(size(p) < lmax+1) then
             Write(88,*) "Error --- PlegendreL" 
        Write(88,*) "P must be dimensioned as (LMAX+1) where LMAX is ", lmax
        Write(88,*) "Input array is dimensioned ", size(p)
         Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error --- PlegendreL" 
        write(0 ,*) "P must be dimensioned as (LMAX+1) where LMAX is ", lmax
        write(0 ,*) "Input array is dimensioned ", size(p)
        write(0 ,*) "The program will STOP ----------------"       
                call stop_config 
            Stop
         elseif (lmax < 0) then 
            Write(88,*) "Error --- PlegendreL" 
        Write(88,*) "LMAX must be greater than or equal to 0."
        Write(88,*) "Input value is ", lmax
         Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error --- PlegendreL" 
        write(0 ,*) "LMAX must be greater than or equal to 0."
        write(0 ,*) "Input value is ", lmax
        write(0 ,*) "The program will STOP ----------------"       
                call stop_config 
            Stop
         elseif(abs(z) > 1.) then
            Write(88,*) "Error --- PlegendreL" 
        Write(88,*) "ABS(Z) must be less than or equal to 1."
        Write(88,*) "Input value is ", z
         Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error --- PlegendreL" 
        write(0 ,*) "ABS(Z) must be less than or equal to 1."
        write(0 ,*) "Input value is ", z
        write(0 ,*) "The program will STOP ----------------"       
                call stop_config 
            Stop
         endif
          
       pm2  = 1.
          p(1) = 1.
          
          pm1  = z
          p(2) = pm1
          
          do l = 2, lmax
             pl = (dble(2*l-1) * z * pm1 - dble(l-1) * pm2)/dble(l)
             p(l+1) = pl
             pm2  = pm1
             pm1  = pl
          enddo
!
end subroutine PLegendre
!
!
!
!
!
!
  SUBROUTINE PLEG (LMAX, LAT, LEGP)

  INTEGER, PARAMETER :: LLMAX=256
  INTEGER L, LMAX 
  REAL(dp) LAT, LEGP(0:LLMAX), P(1:LLMAX+1)
!
! Given latitude LAT (in degrees), computes *all* the Legendre polynomials P_{l}(x), 
! with x=cos(colatitude), to degree LMAX. Uses the SHTOOLS subroutine PLegendre in
! REAL(dp) precision - Note that here the degree 0 is at place '0' in the array PLEG. 
! Last modified by GS 10-9-2007 -  
!
  If(abs(lat).gt.90._dp) then 
             Write(88,*) "Error in Sbr. PLEG: Latitude is out of range" 
         Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error in Sbr. PLEG: Latitude is out of range" 
         write(0 ,*) "The program will STOP ----------------"
                call stop_config 
            Stop 1
  elseif(lmax>llmax) then 
        Write(88,*) "Error in Sbr. PLEG: The degree exceeds 256", lmax 
         Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error in Sbr. PLEG: The degree exceeds 256", lmax 
         write(0 ,*) "The program will STOP ----------------"
                call stop_config 
            Stop 2        
  Endif  
!
!
! --- Computes all the the P_l's to degree LMAX using Mark Wieczorek's code 
!     Note that p(1)=Leg(0), p(2)=Leg(1), ..., p(LMAX+1)=Leg(LMAX) 
!
      call PLegendre(p, lmax, cosdd(90.-lat))
!  
!
! --- Shifts according to our conventions (degree 0 is in place 0 in the array)
!
    do l=0, lmax
     legp(l)=p(l+1)
    end do 
! 
  END SUBROUTINE PLEG
!
!
!
FUNCTION J_INDEX(L,M) 

INTEGER J_INDEX, L, M 
!
! given l (degree) and m(order), this function computes 
! the j-index "J", with J=l(l+1)/2+m+1 - GS 14-09-2007
!
 If(l<0.or.m<0.or.m>l) then 
      Write(88,*) "J_INDEX: The degree and order are out of bounds" 
       Write(88,*) "The program will STOP -------------------------"
       write(0, *) "J_INDEX: The degree and order are out of bounds"  
       write(0, *) "The program will STOP -------------------------"       
          call stop_config 
      Stop
 endif
!
  J_INDEX = l*(l+1)/2+m+1
!
END FUNCTION J_INDEX 
!
!
 INTEGER FUNCTION MJ (J)

   INTEGER J, L, K
   REAL(dp) Z
!
! Given the harmonic index j=l(l+1)/2+m+1, returns m (order) 
! -GS 10.09.07 
!
 If(j<0) then 
      Write(88,*) "MJ: The J-index is out of bounds" 
       Write(88,*) "The program will STOP ----------"
       write(0, *) "MJ: The J-index is out of bounds"  
       write(0, *) "The program will STOP ----------"       
          call stop_config 
      Stop
 endif
!
   mj=0
   k=j 
10 z=(-1.0+sqrt(8.*dble(k)-7.))/2.
   l = int(z) 
   if(z-int(z)==0.) return 
   k=k-1
   mj=mj+1
   goto 10
   return 
   end FUNCTION MJ 

!
!
!
  SUBROUTINE PLMBAR_MOD ( SELEN, LMAX, LAT, PLM)

  TYPE(type_SELEN_global),             INTENT(IN)    :: SELEN
  INTEGER, PARAMETER :: LLMAX=256, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
  INTEGER J, LMAX
  REAL(dp) Z, LAT, PLM(JJMAX)
!
! Given latitude LAT (in degrees), computes *all* the fully normalized 
! associated Legendre functions \bar P_{lm}(x), with x=cos(colatitude), 
! to degree LMAX. Uses the SHTOOLS Legendre functions by PlmBar, and 
! rescales by sqrt(2-delta(0,m)).   Last modified by GS 10-9-2007 - 
!
!
! ---- Tests the latitude bounds 
!
  If(abs(lat).gt.90.) then 
             Write(88,*) "Error in Sbr. PLMBAR_MOD: Latitude is out of range" 
         Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error in Sbr. PLMBAR_MOD: Latitude is out of range" 
         write(0 ,*) "The program will STOP ----------------"
                call stop_config 
            Stop 1
  elseif(lmax>llmax) then 
        Write(88,*) "Error in Sbr. PLMBAR_MOD: The degree exceeds 256", lmax 
         Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error in Sbr. PLMBAR_MOD: The degree exceeds 256", lmax 
         write(0 ,*) "The program will STOP ----------------"
                call stop_config 
            Stop 2        
  Endif  
!
!
!
! ---- "z" is cos(theta), theta is co-latitude 
!
    z=cosdd(90.-lat) 
!
! ---- Builds the SHTOOLS Legendre functions including the Condon-Shortley phase 
!
      call PlmBar(plm, lmax, z, -1)  
!
!
! ---- Scales to obtain our Legendre functions, dividing by sqrt(2-delta(0,m))
!  
      do j=1, j_index(lmax,lmax)
             if(SELEN%MJ_VAL(j)/=0) plm(j)=plm(j)/sqrt(2.)
!             if(mj(j)/=0) plm(j)=plm(j)/sqrt(2.)
      end do
!
  END SUBROUTINE PLMBAR_MOD 
!
!
!
   INTEGER FUNCTION DOM ( SELEN, J)

   TYPE(type_SELEN_global),             INTENT(IN)    :: SELEN
   INTEGER J
!
! Given the index j=l(l+1)/2+m+1, returns 
! delta(0,m), with m=order - GS 14.09.07 
!
 If(j<0) then 
      Write(88,*) "DOM: The J-index is out of bounds" 
       Write(88,*) "The program will STOP -----------"
       write(0, *) "DOM: The J-index is out of bounds"  
       write(0, *) "The program will STOP -----------"       
          call stop_config 
      Stop
 endif
!
   DOM=0 
   If(SELEN%MJ_VAL(j)==0) DOM=1
!   If(mj(j)==0) DOM=1
!
   END FUNCTION DOM
!
!
!
!
!
!
   INTEGER FUNCTION LJ (J)

   INTEGER J, M, K  
   REAL(dp) Z
!
! Given the harmonic index j=l(l+1)/2+m+1, returns l (degree) 
! -GS 10.09.07 
!
 If(j<0) then 
      Write(88,*) "LJ: The J-index is out of bounds" 
       Write(88,*) "The program will STOP ----------"
       write(0, *) "LJ: The J-index is out of bounds"  
       write(0, *) "The program will STOP ----------"       
          call stop_config 
      Stop
 endif
!
   m=0
   k=j 
10 z=(-1.0+sqrt(8.*dble(k)-7.))/2.
   lj=int(z) 
   if(z-int(z)==0.) return 
   k=k-1
   m=m+1
   goto 10
   return 
   end function LJ 
!   
!
!   
!
!
!   
 SUBROUTINE HARMO( SELEN, LMAX, LON, LAT, ARMOY)

  TYPE(type_SELEN_global),             INTENT(IN)    :: SELEN
 INTEGER, PARAMETER :: LLMAX=256, JJMAX=(LLMAX+1)*(LLMAX+2)/2  
 INTEGER J, LMAX  
 COMPLEX*16 ARMOY(JJMAX) 
 REAL(dp) LON, LAT, PLM(JJMAX) 
!
!
! Given longitude LON and latitude LAT  - in degrees - this routine computes 
! *all* the 4-pi normalized complex harmonics $\cal Y_{lm}(\theta,\lambda)$, 
! with \theta=colatitude and \lambda=longitude, to degree LMAX given. Uses  
! the (modified) SHTOOLS Legendre functions by PlmBar (see Sbr PLMBAR_MOD). 
! - Last modified by GS 12-14-2007 - 
!
!
! ---- Tests the longitude and latitude bounds 
!
  If(abs(lat)>90.) then 
             Write(88,*) "Error in Sbr. HARMO: Latitude is out of range" 
         Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error in Sbr. HARMO: Latitude is out of range" 
         write(0 ,*) "The program will STOP ----------------"
                call stop_config 
            Stop 1        
  elseif (lon<0.or.lon>360.) then
             Write(88,*) "Error in Sbr. HARMO: Longitude is out of range" 
         Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error in Sbr. HARMO: Longitude is out of range" 
         write(0 ,*) "The program will STOP ----------------"
                call stop_config 
            Stop 2    
  elseif(lmax>256) then 
        Write(88,*) "Error in Sbr. HARMO: The degree exceeds 256", lmax 
         Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error in Sbr. HARMO: The degree exceeds 256", lmax 
         write(0 ,*) "The program will STOP ----------------"
                call stop_config 
  Endif   
!  
!
! ---- Builds the SHTOOLS Legendre functions, with Condon-Shortley phase 
!
    CALL PLMBAR_MOD( SELEN, LMAX, LAT, PLM)  
!
! 
! ---- Computes the 4-pi normalized Spherical Harmonics up to degree LMAX
!
    do j=1, j_index(lmax,lmax)
           armoy(j) = plm(j)*cmplx(cosdd(dble(SELEN%MJ_VAL(j))*lon),sindd(dble(SELEN%MJ_VAL(j))*lon)) 
!           armoy(j) = plm(j)*cmplx(cosdd(dble(mj(j))*lon),sindd(dble(mj(j))*lon)) 
    end do
!
END SUBROUTINE HARMO
!
!
!
!
!
!
SUBROUTINE PLMBAR(P, LMAX, Z, CSPHASE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    This function evalutates all of the normalized associated Legendre
!    functions up to degree lmax. The functions are initially scaled by 
!    10^280 sin^m in order to minimize the effects of underflow at large m 
!    near the poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299). 
!    On a Mac OSX system with a maximum allowable double precision value of 
!    2.225073858507203E-308 the scaled portion of the algorithm will not overflow 
!    for degrees less than or equal to 2800.
!
!    For each value of m, the rescaling factor is computed as rescalem=rescalem*sin(theta), 
!    with the intial value of rescalem being equal to 1/scalef (which is here equal 
!    to 10^280). This will gradually reduce this huge number to a tiny number, and will 
!    ultimately underflow. In order to prevent this underflow, when rescalem becomes less than
!    10^(-280), the subsequent rescaling factors of sin(theta) will be directly applied to Plm, and then this
!    number will be multipled by the old value of rescalem.
!
!    Temporary variables in saved in an allocated array. In order to explicitly deallocate this
!    memory, call this routine with a spherical harmonic degree of -1.
!
!    Calling Parameters:
!        OUT
!            p:        A vector of all associated Legendgre polynomials evaluated at 
!                    z up to lmax. The length must by greater or equal to (lmax+1)*(lmax+2)/2.
!        OPTIONAL (IN)
!            csphase:    1: Do not include the phase factor of (-1)^m
!                    -1: Apply the phase factor of (-1)^m.
!        IN
!            lmax:        Maximum spherical harmonic degree to compute.
!            z:        cos(colatitude) or sin(latitude).
!
!    Notes:
!    
!    1.    The employed normalization is the "geophysical convention." The integral of
!        (plm*cos(m theta))**2 or (plm*sin (m theta))**2 over all space is 4 pi.
!    2.    The integral of plm**2 over (-1,1) is 2 * (2 - delta(0,m))
!    3.    The index of the array p corresponds to l*(l+1)/2 + m + 1. As such
!        the array p should be dimensioned as (lmax+1)*(lmax+2)/2 in the 
!        calling routine.
!    4.     The default is to exclude the Condon-Shortley phase of (-1)^m.
!
!
!    Dependencies:    CSPHASE_DEFAULT
!
!    Written by Mark Wieczorek September 25, 2005.
!
!    Copyright (c) 2005, Mark A. Wieczorek
!    All rights reserved.
!
!
! ----> Modified by GS December 14, 2007 - error messages 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!
    integer, intent(in) ::    lmax
    REAL(dp), intent(out) ::    p(:)
           REAL(dp), intent(in) ::    z
           integer, intent(in), optional :: csphase
           REAL(dp) ::    pmm, rescalem, phase, u, scalef
          REAL(dp), save, allocatable ::    f1(:), f2(:), sqr(:)
          integer ::    k, kstart, m, l, astat(3)
          integer, save ::    lmax_old  = 0
!    
        integer CSPHASE_DEFAULT
!
    ! The default is to EXCLUDE the CONDON-SHORTLEY phase of (-1)^m
    ! in front of the Legendre functions.
    ! To use this phase function, set CSPHASE_DEFAULT = -1
    CSPHASE_DEFAULT = 1

    if (lmax == -1) then
        if (allocated(sqr)) deallocate(sqr)
        if (allocated(f1)) deallocate(f1)
        if (allocated(f2)) deallocate(f2)
        lmax_old = 0
        return
    endif

    if (size(p) < (lmax+1)*(lmax+2)/2) then 
             Write(88,*) "Error --- PlmBar"
             Write(88,*) "P must be dimensioned as (LMAX+1)*(LMAX+2)/2 where LMAX is ", lmax
             Write(88,*) "Input array is dimensioned ", size(p)
         Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error --- PlmBar"
             write(0 ,*) "P must be dimensioned as (LMAX+1)*(LMAX+2)/2 where LMAX is ", lmax
             write(0 ,*) "Input array is dimensioned ", size(p)
         write(0 ,*) "The program will STOP ----------------"
                call stop_config 
            Stop         
    elseif (lmax < 0) then 
             Write(88,*) "Error --- PlmBar"
             Write(88,*) "LMAX must be greater than or equal to 0."
             Write(88,*) "Input value is ", lmax
        Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error --- PlmBar"
             write(0 ,*) "LMAX must be greater than or equal to 0."
             write(0 ,*) "Input value is ", lmax
        write(0 ,*) "The program will STOP ----------------"
                call stop_config 
            Stop     
    elseif(abs(z) > 1.00) then
             Write(88,*) "Error --- PlmBar"
             Write(88,*) "ABS(Z) must be less than or equal to 1."
             Write(88,*) "Input value is ", z
        Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "Error --- PlmBar"
             write(0 ,*) "ABS(Z) must be less than or equal to 1."
             write(0 ,*) "Input value is ", z
        write(0, *) "The program will STOP ----------------"
                call stop_config 
            Stop     
    endif         
!         
         if (present(csphase)) then
             if (csphase == -1) then
                 phase = -1.00
             elseif (csphase == 1) then
                 phase = 1.00
             else
             Write(88,*) "PlmBar --- Error"
             Write(88,*) "CSPHASE must be 1 (exclude) or -1 (include)."
             Write(88,*) "Input value is ", csphase
        Write(88,*) "The program will STOP ----------------"
             write(0 ,*) "PlmBar --- Error"
             write(0 ,*) "CSPHASE must be 1 (exclude) or -1 (include)."
             write(0 ,*) "Input value is ", csphase
        write(0 ,*) "The program will STOP ----------------"
                call stop_config 
            Stop 
             endif
         else
             phase = dble(CSPHASE_DEFAULT)
         endif
             
    scalef = 1.0e-20
    
    
    if (lmax /= lmax_old) then
        
        if (allocated(sqr)) deallocate(sqr)
        if (allocated(f1)) deallocate(f1)
        if (allocated(f2)) deallocate(f2)
        
        allocate(sqr(2*lmax+1), stat=astat(1))
        allocate(f1((lmax+1)*(lmax+2)/2), stat=astat(2))
        allocate(f2((lmax+1)*(lmax+2)/2), stat=astat(3))
        
        if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0) then
             write(0 ,*) "PlmBar --- Error"
             write(0 ,*) "Problem allocating arrays SQR, F1 and F2", astat(1), astat(2), astat(3)
                call stop_config 
            Stop 
        endif
            
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        !    Precompute square roots of integers that are used several times.
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        do l=1, 2*lmax+1
            sqr(l) = sqrt(dble(l))   !sqrt(dble(l))
        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        !     Precompute multiplicative factors used in recursion relationships
        !         Plmbar(l,m) = x*f1(l,m)*Plmbar(l-1,m) - Plmbar(l-2,m)*f2(l,m)
        !        k = l*(l+1)/2 + m + 1
        !    Note that prefactors are not used for the case when m=l and m=l-1,
        !    as a different recursion is used for these two values.
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        k = 3
    
        do l=2, lmax, 1
            k = k + 1
            f1(k) = sqr(2*l-1) * sqr(2*l+1) / dble(l)   !dble(l)
            f2(k) = dble(l-1) * sqr(2*l+1) / sqr(2*l-3) / dble(l)   !dble(l)
            do m=1, l-2
                k = k+1
                f1(k) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                        f2(k) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
                               / sqr(2*l-3) / sqr(l+m) / sqr(l-m) 
            enddo
            k = k + 2
        enddo
    
        lmax_old = lmax
    
    endif
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    
    !    Calculate P(l,0). These are not scaled.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    u = sqrt((1.00-z)*(1.00+z)) ! sin(theta)
    
          p(1) = 1.00
          
          if (lmax == 0) return
          
          p(2)  = sqr(3)*z
          
          k = 2

          do l = 2, lmax, 1
             k = k+l
             p(k) = f1(k)*z*p(k-l)-f2(k)*p(k-2*l+1)
          enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !    Calculate P(m,m), P(m+1,m), and P(l,m)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
          pmm = sqr(2)*scalef
          rescalem = 1./scalef
          kstart = 1

          do m = 1, lmax - 1, 1
          
              rescalem = rescalem * u
          
        ! Calculate P(m,m)
            kstart = kstart+m+1
         
             pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
            p(kstart) = pmm

        ! Calculate P(m+1,m)
        k = kstart+m+1
           p(k) = z * sqr(2*m+3) * pmm

        ! Calculate P(l,m)
                   do l = m+2, lmax, 1
                       k = k+l
                      p(k) = z*f1(k)*p(k-l)-f2(k)*p(k-2*l+1)
                      p(k-2*l+1) = p(k-2*l+1) * rescalem
                   enddo
                   
                   p(k) = p(k) * rescalem
                   p(k-lmax) = p(k-lmax) * rescalem
                                 
          enddo
          
          ! Calculate P(lmax,lmax)
          
          rescalem = rescalem * u
                
        kstart = kstart+m+1
        p(kstart) = phase * pmm * sqr(2*lmax+1) / sqr(2*lmax) * rescalem
              
end subroutine PlmBar
!
!
!
    SUBROUTINE READ_JUNK(I,J)
    IMPLICIT NONE 
    INTEGER I, J, L  
    CHARACTER*50 JUNK
!   
! ///////////////////////////////////////////////
!
! Reads "i" junk 50-characters lines from unit "j"
! Added by GS on August 2, 2010 
!
! ///////////////////////////////////////////////
!   
    DO 1 L=1, I 
             READ(J,'(A50)')JUNK 
1       CONTINUE    
!
    END Subroutine read_junk
!

!
!
!
!
END MODULE SELEN_harmonics_module
