MODULE SELEN_sealevel_equation_module

  ! Contains the main routine that solves the sea-level equation (with some small help routines).
  ! This is called from "run_SELEN" in the SELEN_main_module.

  USE mpi
  USE configuration_module,              ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE parallel_module,                   ONLY: par, sync, cerr, ierr, &
                                               allocate_shared_int_0D, allocate_shared_dp_0D, &
                                               allocate_shared_int_1D, allocate_shared_dp_1D, &
                                               allocate_shared_int_2D, allocate_shared_dp_2D, &
                                               allocate_shared_int_3D, allocate_shared_dp_3D, &
                                               allocate_shared_complex_1D, allocate_shared_complex_2D, allocate_shared_complex_3D, &
                                               deallocate_shared, partition_list
  USE data_types_module,                 ONLY: type_model_region, type_SELEN_global
  USE utilities_module,                  ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                               check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  USE SELEN_spline_interpolation_module, ONLY: spline, splint
  USE SELEN_taboo_hp_module,             ONLY: BETAS, BETAU, ES, EU
  
  IMPLICIT NONE
 
CONTAINS

  SUBROUTINE solve_SLE( SELEN, NAM, EAS, GRl, ANT, time, ocean_area, ocean_depth)
    ! Solve the sea-level equation for the current ice loading history 
    
    IMPLICIT NONE
 
    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(INOUT) :: SELEN
    TYPE(type_model_region),             INTENT(INOUT) :: NAM, EAS, GRL, ANT
    REAL(dp),                            INTENT(IN)    :: time
    REAL(dp),                            INTENT(OUT)   :: ocean_area   ! Area of the Worlds oceans (m^2) to be used in Anice
    REAL(dp),                            INTENT(OUT)   :: ocean_depth  ! Averaged depth of the Worlds oceans (m) to be used in Anice
 
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_SLE'
    INTEGER                               :: i,j,k,p,iters,tdof_iteration
    REAL(dp)                              :: resh,imsh
    REAL(dp)                              :: RHOI_O_RHOE_X3, RHOW_O_RHOE_X3, RHOI_O_RHOW, NP_INV
    
    INTEGER,    DIMENSION(:    ), POINTER :: DM                  ! INTEGER,    DIMENSION( C%SELEN_jmax                ) :: DM
    REAL(dp),   DIMENSION(:,:  ), POINTER :: X                   ! REAL(dp),   DIMENSION( SELEN%mesh%nV,  0:C%SELEN_irreg_time_n  ) :: X
    REAL(dp),   DIMENSION(            0:C%SELEN_irreg_time_n  ) :: AAVV
    REAL(dp),   DIMENSION(            0:C%SELEN_irreg_time_n  ) :: BBVV
   
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: ZE                  ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: ZE                  ! Eustatic Z array for water loading
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: SE                  ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: SE                  ! Eustatic S array for ice loading
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: AAAA                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: AAAA
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: AAAA_MOD            ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: AAAA_MOD
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: BBBB                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: BBBB
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: BBBB_MOD            ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: BBBB_MOD
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: HHHH                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: HHHH
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: KKKK                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: KKKK
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: TTTT                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: TTTT                ! Ice loading in SH
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: IIII                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: IIII                ! Relative ice loading in SH
    
    ! New arrays for the rotational feedback
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: load_ice_and_water  ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: load_ice_and_water  ! Ice and water loading in SH for rotational feedback
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: ZROT                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: ZROT
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: ZROT_MOD            ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: ZROT_MOD
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: SROT                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: SROT
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: SROT_MOD            ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: SROT_MOD
    COMPLEX*16, DIMENSION(:    ), POINTER :: ZROTVV              ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: ZROTVV
    COMPLEX*16, DIMENSION(:    ), POINTER :: SROTVV              ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: SROTVV
 
    ! Results
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: S                   ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: S                   ! relative sea level
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: U                   ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: U                   ! Solid Earth deformation (N = S+U) 
    COMPLEX*16, DIMENSION(:,:,:), POINTER :: Z                   ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n,0:C%SELEN_n_recursion_iterations) :: Z
 
    ! New arrays for tdof
    REAL(dp),   DIMENSION(:,:  ), POINTER :: newalf              ! REAL(dp),   DIMENSION( C%SELEN_jmax,  C_SLE%NANCH  ) :: newalf              ! New values of alf
    REAL(dp),   DIMENSION(:,:  ), POINTER :: slc                 ! REAL(dp),   DIMENSION( SELEN%mesh%nV,  0:C%SELEN_irreg_time_n  ) :: slc                 ! sea level forward in time for all elements
    INTEGER,    DIMENSION(:,:  ), POINTER :: WET                 ! INTEGER,    DIMENSION( SELEN%mesh%nV,  0:C%SELEN_irreg_time_n  ) :: WET                 ! new ocean function (0 or 1) forward in time
    REAL(dp),   DIMENSION(:,:  ), POINTER :: newtopo             ! REAL(dp),   DIMENSION( SELEN%mesh%nV,  0:C%SELEN_irreg_time_n  ) :: newtopo             ! new topograhpy forward in time
 
    INTEGER,    DIMENSION(:,:  ), POINTER :: fj                  ! INTEGER,    DIMENSION( SELEN%mesh%nV,  0:C%SELEN_irreg_time_n  ) :: fj                  ! floating function for ice thickness (float or not)
    INTEGER,    DIMENSION(:,:  ), POINTER :: xcf                 ! INTEGER,    DIMENSION( SELEN%mesh%nV,  0:C%SELEN_irreg_time_n+1) :: xcf
 
    INTEGER,    DIMENSION(:    ), POINTER :: newet2              ! INTEGER,    DIMENSION( SELEN%mesh%nV                  ) :: newet2              ! new ocean function after iteration
    REAL(dp),   DIMENSION(:    ), POINTER :: newtopo2            ! REAL(dp),   DIMENSION( SELEN%mesh%nV                  ) :: newtopo2            ! new topograhpy after iteration
 
    REAL(dp),   DIMENSION(:    ), POINTER :: S_global, U_global  ! REAL(dp),   DIMENSION( SELEN%mesh%nV                  ) :: S_global, U_global  ! Global fields of S and U for output
 
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: MSint               ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: MSint               ! Interpolated memory of relative sea level
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: MUint               ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: MUint               ! Interpolationd memory of solid earth deformation (N is geoid)
 
    COMPLEX*16, DIMENSION(:    ), POINTER :: MSread              ! COMPLEX*16, DIMENSION( C%SELEN_jmax                ) :: MSread              ! Initial memory of relative sea level
 
    ! New variables for spline interpolation
    REAL(dp),   DIMENSION(:    ), POINTER :: int_time            ! REAL(dp),   DIMENSION(            0:C%SELEN_irreg_time_n  ) :: int_time            ! time points (depending on DELTA) for interpolation
    COMPLEX*16, DIMENSION(:    ), POINTER :: sp1, spn, preSint   ! COMPLEX*16, DIMENSION( C%SELEN_jmax)                 :: sp1, spn, preSint   ! spline in- and output for S
    COMPLEX*16, DIMENSION(:    ), POINTER :: up1, upn, preUint   ! COMPLEX*16, DIMENSION( C%SELEN_jmax)                 :: up1, upn, preUint   ! spline in- and output for U
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: s2,u2               ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: s2,u2               ! spline output
    REAL(dp)                              :: xi                  ! REAL(dp)                                           :: xi                  ! time point of interpolation
 
    ! for the interpolation
    REAL(dp)                              :: Tbeg, Tend, Tint    ! interpolation time points
 
    REAL(dp),   DIMENSION(0:C%SELEN_irreg_time_n)  :: Seust               ! eustatic rsl change for each time step
    REAL(dp)                              :: Seust2              ! total eustatic rsl change
    REAL(dp)                              :: Tw1,Tw2             ! changes in ice loading in meters eustatic sea level, area of the ocean
 
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: varpreoc            ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: varpreoc            ! pre determined SH coefficients ocean function
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: varoc               ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: varoc               ! variable SH coefficient of ocean function
    COMPLEX*16, DIMENSION(:    ), POINTER :: varoc_inv           ! COMPLEX*16, DIMENSION(            0:C%SELEN_irreg_time_n  ) :: varoc_inv           ! 1 / varoc(1,:) - the first degree

    ! Saved memory of S and U to add all viscous responses
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: MEM_S               ! ALLOCATE( MEM_S( C%SELEN_jmax,0:C%SELEN_irreg_time_n))
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: MEM_U               ! ALLOCATE( MEM_U( C%SELEN_jmax,0:C%SELEN_irreg_time_n))
    REAL(dp),   DIMENSION(:,:  ), POINTER :: icethick            ! ALLOCATE( icethick( SELEN%mesh%nV, 0:C%SELEN_irreg_time_n ))
    
    ! Data that used to be read from other modules (now stored in shared memory in SELEN structure):
    INTEGER,    DIMENSION(:    ), POINTER :: MJ_VAL              ! ALLOCATE(MJ_VAL(C%SELEN_jmax))
    INTEGER,    DIMENSION(:    ), POINTER :: LJ_VAL              ! ALLOCATE(LJ_VAL(C%SELEN_jmax))
    INTEGER,    DIMENSION(:    ), POINTER :: SELEN_anc           ! ALLOCATE(SELEN_anc (SELEN%mesh%nV))
    REAL(dp),   DIMENSION(:    ), POINTER :: INIT_TOPO           ! ALLOCATE(INIT_TOPO (SELEN%mesh%nV))
    REAL(dp),   DIMENSION(:,:  ), POINTER :: ALF                 ! ALLOCATE(ALF       (C%SELEN_jmax,C_SLE%NANCH))
    COMPLEX*16, DIMENSION(:,:  ), POINTER :: LONG_TABLE          ! ALLOCATE(LONG_TABLE(0:C_SLE%LMAX,SELEN%mesh%nV))
    REAL(dp),   DIMENSION(:    ), POINTER :: INIT_ICE            ! ALLOCATE(INIT_ICE  (SELEN%mesh%nV))
    
    ! Track progress of large loops for printing to terminal
    INTEGER                               :: progprev, prog
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Allocate memory and bind pointers
    CALL allocate_memory_and_bind_pointers( SELEN, DM, X, ZE, SE, AAAA, AAAA_MOD, BBBB, BBBB_MOD, HHHH, KKKK, TTTT, IIII, &
      load_ice_and_water, ZROT, ZROT_MOD, SROT, SROT_MOD, ZROTVV, SROTVV, S, U, Z, newalf, slc, WET, newtopo, fj, xcf, &
      newet2, newtopo2, S_global, U_global, MSint, MUint, MSread, int_time, sp1, spn, preSint, up1, upn, preUint, s2,u2, &
      varpreoc, varoc, varoc_inv, MEM_S, MEM_U, icethick, MJ_VAL, LJ_VAL, SELEN_anc, INIT_TOPO, ALF, LONG_TABLE, INIT_ICE)
    
    ! Some often-used derived constants
    RHOI_O_RHOE_X3 = 3._dp * ice_density / earth_density
    RHOW_O_RHOE_X3 = 3._dp * seawater_density / earth_density
    RHOI_O_RHOW    =         ice_density / seawater_density
    NP_INV         = 1._dp / DBLE(SELEN%mesh%nV)
    
    ! Share data from TABOO among the processes (since TABOO isnt parallelised and doesnt use shared memory)
    IF (.NOT. par%master) THEN
      ALLOCATE( BETAS( 0:C%SELEN_n_harmonics, 0:C%SELEN_reg_time_n))
      ALLOCATE( BETAU( 0:C%SELEN_n_harmonics, 0:C%SELEN_reg_time_n))
      ALLOCATE( ES   ( 0:C%SELEN_n_harmonics                      ))
      ALLOCATE( EU   ( 0:C%SELEN_n_harmonics                      ))
    END IF
    CALL MPI_BCAST( BETAS, (C%SELEN_n_harmonics+1)*(C%SELEN_reg_time_n+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( BETAU, (C%SELEN_n_harmonics+1)*(C%SELEN_reg_time_n+1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( ES,    (C%SELEN_n_harmonics+1)                       , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( EU,    (C%SELEN_n_harmonics+1)                       , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
 
    ! Pre-computing 'll' and 'mm' corresponding to degree 'J'
    newalf( C%SELEN_j1:C%SELEN_j2,:) = 0._dp
    DO j = C%SELEN_j1, C%SELEN_j2
      dm(j) = 2
      IF (MJ_VAL(j) == 0) dm(j) = 1
      newalf(j  ,:) = ALF(j  ,:) * dm(j  )
    END DO
    CALL sync
    
    ! Initialize Sea level geoid
    S( C%SELEN_j1:C%SELEN_j2,:) = (0.,0.)
    U( C%SELEN_j1:C%SELEN_j2,:) = (0.,0.)
    CALL sync

    ! ===================================
    ! == Start of TDOF-NFWD iterations ==
    ! ===================================

    DO tdof_iteration = 1, C%SELEN_n_TDOF_iterations
    
      IF (par%master) WRITE (0,'(A,I3,A,I3)') '   TDOF iteration ',tdof_iteration, ' of ', C%SELEN_n_TDOF_iterations
      
      ! For each new iteration calculate slc change from MEM_S at each pixel (global mesh) and in time
      IF (par%master) WRITE (0,'(A)') '    Calculating sea-level change from MEM_S...'
      progprev = -1
      
      DO i = C%SELEN_i1, C%SELEN_i2
            
        prog = FLOOR(100*REAL(i)/REAL(SELEN%mesh%nV))
        IF (par%master .AND. prog>progprev .AND. C%SELEN_display_progress) THEN
          progprev = prog
          WRITE(0,'(A,I3,A)',ADVANCE='NO') CHAR(13)//'      Progress: ', prog, ' %'
        END IF
      
        slc(i,:) = 0._dp
        DO j = 1, C%SELEN_jmax
          ! here S is summed to the global sea-level change
          slc(i,:) = slc(i,:) + newalf(j,SELEN_anc(i))*dble(MEM_S(j,:) * LONG_TABLE(MJ_VAL(j),i))
        END DO
      END DO ! DO i = 1, SELEN%mesh%nV 
      IF (par%master .AND. C%SELEN_display_progress) WRITE(0,'(A)',ADVANCE='NO') CHAR(13)
    
      ! Generate (0:C%SELEN_irreg_time_n) topographies from slc that is calculated at the previous CALL
      CALL futopos( SELEN, icethick, INIT_TOPO, slc, WET, fj, newtopo)
      
      ! Update Ocean functions coefficients (C%SELEN_irreg_time_n ocean function sh tables are computed)
      ! First calculate the "continental function" for wet=0 and compute the OF taking into account that OF + CF = 1
      IF (par%master) WRITE (0,'(A)') '    Updating ocean functions coefficients...'
      progprev = -1
      
      DO k = 0, C%SELEN_irreg_time_n
            
        prog = FLOOR(100*REAL(k)/REAL(C%SELEN_irreg_time_n))
        IF (par%master .AND. prog>progprev .AND. C%SELEN_display_progress) THEN
          progprev = prog
          WRITE(0,'(A,I3,A)',ADVANCE='NO') CHAR(13)//'      Progress: ', prog, ' %'
        END IF
        
        DO j = C%SELEN_j1, C%SELEN_j2
          varpreoc(j,k) = 0.
          DO i = 1, SELEN%mesh%nV 
            IF( WET( i,k)==0) varpreoc(j,k) = varpreoc(j,k) + NP_INV * ALF(j,SELEN_anc(i)) * conjg(LONG_TABLE(MJ_VAL(j),i))
          END DO ! DO i = 1, SELEN%mesh%nV 
  
          IF(j==1) THEN
            resh       = 1._dp - dble(varpreoc(j,k))
          ELSE
            resh       = -dble(varpreoc(j,k))
          END IF
          imsh       = -dimag(varpreoc(j,k))   ! or aimag..
          varoc(j,k) = cmplx(resh, imsh)
        END DO ! DO j = 1, C%SELEN_jmax
        CALL sync
        varoc_inv(k) = 1._dp / varoc(1,k)
      END DO ! DO k = 0, C%SELEN_irreg_time_n
      IF (par%master .AND. C%SELEN_display_progress) WRITE(0,'(A)',ADVANCE='NO') CHAR(13)
      
      ! Compute the spherical harmonic transform IIII of the ice loading history
      IF (par%master) WRITE (0,'(A)') '    Computing spherical harmonics transform of ice loading history...'
      CALL update_ice_sh_and_load( SELEN, NAM, EAS, GRL, ANT, icethick, TTTT, IIII)
    
      ! Computing the eustatic Z and S arrays (always from 0 to C%SELEN_irreg_time_n, read IIII for the respective time steps)
      IF (par%master) WRITE (0,'(A)') '    Computing the eustatic Z array...'
      ZE( C%SELEN_j1:C%SELEN_j2,:) = (0.,0.)
      SE( C%SELEN_j1:C%SELEN_j2,:) = (0.,0.)
      CALL sync
      DO k = 0,C%SELEN_irreg_time_n
        ZE(  C%SELEN_j1:C%SELEN_j2,k) = - RHOI_O_RHOW * (IIII(1,k) * varoc_inv(k)) * varoc(C%SELEN_j1:C%SELEN_j2,k)
        SE(  1                    ,k) = - RHOI_O_RHOW * (IIII(1,k) * varoc_inv(k)) 
        xcf( C%SELEN_j1:C%SELEN_j2,k) =   WET( C%SELEN_j1:C%SELEN_j2,k)
      END DO ! DO k = 0,C%SELEN_irreg_time_n
      xcf(:,C%SELEN_irreg_time_n+1) = xcf(:,C%SELEN_irreg_time_n)
      CALL sync
    
      IF (par%master) WRITE (0,'(A)') '    Computing the A array...'
      progprev = -1
      
      ! I(p) - I(p-1), should we take into account the time step length here?? Compute A relative to the previous time step, but..
      ! time step lenght is that an issue here... Ok, is taken into account with BetaS..
      AAAA( C%SELEN_j1:C%SELEN_j2,:) = (0.,0.)
      DO k = 0, C%SELEN_irreg_time_n   ! only use TCONV, the last time step including all history for AAAA
          
        prog = FLOOR(100*REAL(k)/REAL(C%SELEN_irreg_time_n))
        IF (par%master .AND. prog>progprev .AND. C%SELEN_display_progress) THEN
          progprev = prog
          WRITE(0,'(A,I3,A)',ADVANCE='NO') CHAR(13)//'      Progress: ', prog, ' %'
        END IF
      
        ! Computing the ocean average of A 
        AAVV(k) = 0.      
        DO j = C%SELEN_j1, C%SELEN_j2
          AAAA(j,k) = ES(LJ_VAL(j)) * IIII(j,k)
  
          DO p = 0, k
            IF     (p == 0 .AND. k == 0) THEN
              AAAA( j,k) = AAAA( j,k) - (IIII( j,p)               ) * BETAS( LJ_VAL( j),0)
            ELSEIF (p == 0 .AND. k /= 0) THEN
              AAAA( j,k) = AAAA( j,k) - (IIII( j,p)               ) * BETAS( LJ_VAL( j),INT(SUM(C%SELEN_irreg_time_window(p+1:k))))
            ELSEIF (p >  0 .AND. p <  k) THEN
              AAAA( j,k) = AAAA( j,k) - (IIII( j,p) - IIII( j,p-1)) * BETAS( LJ_VAL( j),INT(SUM(C%SELEN_irreg_time_window(p+1:k))))
            ELSE ! p==k
              AAAA( j,k) = AAAA( j,k) - (IIII( j,p) - IIII( j,p-1)) * BETAS( LJ_VAL( j),0)
            END IF
            ! For values of Beta see table in taboo module, where betas/u are computed..
          END DO ! DO p = 0, k
          AAAA( j,k) = RHOI_o_RHOE_x3 * AAAA( j,k)
          AAVV(   k) = AAVV( k) + dm( j) * DBLE( varoc( j,k) * CONJG( AAAA( j,k))) * varoc_inv( k)
        END DO ! DO j = 1, C%SELEN_jmax
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, AAVV( k), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      END DO ! DO k = 0, C%SELEN_irreg_time_n
      IF (par%master .AND. C%SELEN_display_progress) WRITE(0,'(A)',ADVANCE='NO') CHAR(13)
      
      ! Computing the modified ocean average of A 
      AAAA_MOD( C%SELEN_j1:C%SELEN_j2,:) = AAAA( C%SELEN_j1:C%SELEN_j2,:)
      IF (par%master) AAAA_MOD(1,:) = AAAA(1,:) - AAVV
      CALL sync
      
      IF (par%master) WRITE (0,'(A)') '    Computing the R array...'
      progprev = -1
      
      DO k = 0, C%SELEN_irreg_time_n
          
        prog = FLOOR(50*REAL(k)/REAL(C%SELEN_irreg_time_n))
        IF (par%master .AND. prog>progprev .AND. C%SELEN_display_progress) THEN
          progprev = prog
          WRITE(0,'(A,I3,A)',ADVANCE='NO') CHAR(13)//'      Progress: ', prog, ' %'
        END IF
        
        DO i = C%SELEN_i1, C%SELEN_i2
          X(i,k) = 0.
          DO j = 1, C%SELEN_jmax
            X( i,k) = X( i,k) + ALF( j,SELEN_anc( i)) * dm( j) * DBLE( AAAA_MOD( j,k) * LONG_TABLE( MJ_VAL( j),i))
          END DO ! DO j = 1, C%SELEN_jmax
        END DO ! DO i = 1, SELEN%mesh%nV
        CALL sync
      END DO ! DO k = 0, C%SELEN_irreg_time_n
      
      DO k = 0, C%SELEN_irreg_time_n
          
        prog = 50 + FLOOR(50*REAL(k)/REAL(C%SELEN_irreg_time_n))
        IF (par%master .AND. prog>progprev .AND. C%SELEN_display_progress) THEN
          progprev = prog
          WRITE(0,'(A,I3,A)',ADVANCE='NO') CHAR(13)//'      Progress: ', prog, ' %'
        END IF
        
        DO j = C%SELEN_j1, C%SELEN_j2
          HHHH( j,k) = (0.,0.)
          DO i = 1, SELEN%mesh%nV ! fusion
            HHHH( j,k) = HHHH( j,k) + (X( i,k) * xcf( i,k+1) + newtopo( i,k) * (xcf( i,k) - xcf( i,k+1))) * &
                         ALF( j,SELEN_anc( i)) * CONJG( LONG_TABLE( MJ_VAL( j),i)) * NP_INV
          END DO ! DO i = 1, SELEN%mesh%nV
        END DO ! DO j = 1, C%SELEN_jmax
        CALL sync 
      END DO ! DO k = 0, C%SELEN_irreg_time_n
      IF (par%master .AND. C%SELEN_display_progress) WRITE(0,'(A)',ADVANCE='NO') CHAR(13)

    ! =========================
    ! == Rotational feedback ==
    ! =========================
    
!      IF (C%SELEN_use_rotational_feedback) THEN
!      
!        IF (par%master) WRITE(0,*) 'ERROR - rotational feedback isnt parallelised and also it doesnt work!'
!        STOP
!      
!        ! Call the polar motion transfer function; generates a file which is used by the s_rotaz
!        CALL PMTF()
!
!        ! set total load of ice and ocean
!        load_ice_and_water = ice_density * IIII + seawater_density * ZE
!
!        ! calculate rotational feedback Srot
!        CALL s_rotaz(load_ice_and_water,srot)
!        
!        ! Computing S^{rot} at all pixels [Note: Using variable "x" in order to save memory]
!        DO k = 0, C%SELEN_irreg_time_n
!          DO i = 1, SELEN%mesh%nV
!            x(i,k)=0d0
!            DO j = 1, C%SELEN_jmax
!              x(i,k) = x(i,k) + ALF(j,SELEN_anc(i))*dm(j)*dble(SROT(j,k)*long_table(MJ_VAL(j),i))
!            END DO ! DO j = 1, C%SELEN_jmax
!          END DO ! DO i = 1, SELEN%mesh%nV
!        END DO ! DO k = 0, C%SELEN_irreg_time_n
!        
!        ! Computing the harmonics of Z^{rot} by integration over wet pixels 
!        ! Computing the ocean average of S^{rot}  
!        ! Computing the ocean average of Z^{rot}  
!
!        zrot = 0d0
!        DO k = 0, C%SELEN_irreg_time_n
!          srotvv(k) = 0.
!          zrotvv(k) = 0.
!          DO j = 1, C%SELEN_jmax 
!            zrot(j,k) = 0d0
!            DO i = 1, SELEN%mesh%nV 
!              IF (WET(i,k)==1) zrot(j,k) = zrot(j,k) + x(i,k)*alf(j,SELEN_anc(i))*conjg(long_table(MJ_VAL(j),i))  
!            END DO ! DO i = 1, SELEN%mesh%nV 
!            srotvv(k) = srotvv(k) + dm(j)*dble(varoc(j,k)*conjg(srot(j,k))) * varoc_inv(k)
!            zrotvv(k) = zrotvv(k) + dm(j)*dble(varoc(j,k)*conjg(zrot(j,k))) * varoc_inv(k)
!          END DO ! DO j = 1, C%SELEN_jmax 
!        END DO ! DO k = 0, C%SELEN_irreg_time_n
!        zrot = zrot * NP_INV
!        
!        ! Computing the modified array S^{rot} 
!        srot_mod = srot
!        srot_mod(1,:) = srot(1,:) - srotvv(:) 
!        
!        ! Computing the modified array Z^{rot} 
!        zrot_mod      = zrot
!        zrot_mod(1,:) = zrot(1,:) - zrotvv(:) 
!        
!      END IF ! IF (C%SELEN_use_rotational_feedback) THEN
    
    ! ===============================
    ! == Rotational feedback - end ==
    ! ===============================
    
      IF (par%master) WRITE (0,'(A)') '    Initializing the Z and S arrays...'
      ! Initializing the Z and S arrays 
      Z( C%SELEN_j1:C%SELEN_j2,:,0) = ZE( C%SELEN_j1:C%SELEN_j2,:)
      S( C%SELEN_j1:C%SELEN_j2,:  ) = SE( C%SELEN_j1:C%SELEN_j2,:)
      CALL sync
    
    ! ===============
    ! == Recursion ==
    ! ===============
      
      IF ( C%SELEN_n_recursion_iterations /= 0) THEN

        IF (par%master) WRITE (0,'(A)') '    Starting the recursion...'
        
        DO iters = 1, C%SELEN_n_recursion_iterations
        
          IF (par%master) WRITE (0,'(A,I3,A,I3)') '     Iteration ', iters, '/', C%SELEN_n_recursion_iterations
          
          ! Computing the 'B' array
          IF (par%master) WRITE (0,'(A)') '      Computing the B array...'
          progprev = -1
          
          BBBB( C%SELEN_j1:C%SELEN_j2,:) = (0.,0.)
          CALL sync
          DO k = 0,C%SELEN_irreg_time_n
            
            prog = FLOOR(100*REAL(k)/REAL(C%SELEN_irreg_time_n))
            IF (par%master .AND. prog>progprev .AND. C%SELEN_display_progress) THEN
              progprev = prog
              WRITE(0,'(A,I3,A)',ADVANCE='NO') CHAR(13)//'      Progress: ', prog, ' %'
            END IF
        
            ! Computing the ocean-average of array B array
            BBVV( k) = 0.
            DO j = C%SELEN_j1, C%SELEN_j2
              BBBB( j,k) = ES( LJ_VAL( j)) * Z( j,k,iters-1)
              DO p = 0,k
                IF     (p == 0 .AND. k == 0) THEN
                  BBBB( j,k) = BBBB( j,k) -(Z( j,p,iters-1) - 0.              ) * BETAS( LJ_VAL( j),0)
                ELSEIF (p == 0 .AND. k /= 0) THEN
                  BBBB( j,k) = BBBB( j,k) -(Z( j,p,iters-1) - 0.              ) * BETAS( LJ_VAL( j),INT(SUM(C%SELEN_irreg_time_window(p+1:k))))
                ELSEIF (p >  0 .AND. p <  k) THEN
                  BBBB( j,k) = BBBB( j,k) -(Z( j,p,iters-1) - Z(j,p-1,iters-1)) * BETAS( LJ_VAL( j),INT(SUM(C%SELEN_irreg_time_window(p+1:k))))
                ELSE ! p==k
                  BBBB( j,k) = BBBB( j,k) -(Z( j,p,iters-1) - Z(j,p-1,iters-1)) * BETAS( LJ_VAL( j),0) 
                END IF
              END DO ! DO p = 0,k
              BBBB( j,k) = RHOW_o_RHOE_X3 * BBBB( j,k)
              BBVV(   k) = BBVV(k) + dm(j) * DBLE( varoc( j,k) * CONJG( BBBB( j,k))) * varoc_inv( k)
            END DO ! DO j = 1, C%SELEN_jmax 
            CALL MPI_ALLREDUCE( MPI_IN_PLACE, BBVV( k), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
          END DO ! DO k = 0,C%SELEN_irreg_time_n
          IF (par%master .AND. C%SELEN_display_progress) WRITE(0,'(A)',ADVANCE='NO') CHAR(13)
          
          ! Computing modified 'B' array
          BBBB_MOD( C%SELEN_j1:C%SELEN_j2,:) = BBBB( C%SELEN_j1:C%SELEN_j2,:)
          IF (par%master) BBBB_MOD(1,:) = BBBB(1,:) - BBVV
          CALL sync
          
          ! Computing array K...
          IF (par%master) WRITE (0,'(A)') '      Computing the K array...'
          progprev = -1
          
          DO k = 0, C%SELEN_irreg_time_n
            
            prog = FLOOR(50*REAL(k)/REAL(C%SELEN_irreg_time_n))
            IF (par%master .AND. prog>progprev .AND. C%SELEN_display_progress) THEN
              progprev = prog
              WRITE(0,'(A,I3,A)',ADVANCE='NO') CHAR(13)//'      Progress: ', prog, ' %'
            END IF
            
            DO i = C%SELEN_i1, C%SELEN_i2
              X( i,k) = 0._dp
              DO j = 1, C%SELEN_jmax
                X( i,k) = X( i,k) + ALF( j,SELEN_anc( i)) * dm( j) * DBLE( BBBB_MOD( j,k) * LONG_TABLE( MJ_VAL( j),i))
              END DO ! DO j = 1, C%SELEN_jmax
            END DO ! DO i = 1, SELEN%mesh%nV 
          END DO ! DO k = 0, C%SELEN_irreg_time_n
          
          DO k = 0, C%SELEN_irreg_time_n
            
            prog = 50+FLOOR(50*REAL(k)/REAL(C%SELEN_irreg_time_n))
            IF (par%master .AND. prog>progprev .AND. C%SELEN_display_progress) THEN
              progprev = prog
              WRITE(0,'(A,I3,A)',ADVANCE='NO') CHAR(13)//'      Progress: ', prog, ' %'
            END IF
            
            DO j = C%SELEN_j1, C%SELEN_j2
              KKKK( j,k) = (0.,0.)
              DO i = 1, SELEN%mesh%nV ! fusion
                KKKK( j,k) = KKKK( j,k) + ( X( i,k) * xcf( i,k+1) + newtopo( i,k) * (xcf( i,k) - xcf(i,k+1))) * &
                             ALF( j,SELEN_anc( i)) * CONJG( LONG_TABLE( MJ_VAL( j),i)) ! * NP_INV 
              END DO ! DO i = 1, SELEN%mesh%nV
            END DO ! DO j = 1, C%SELEN_jmax
            CALL sync
          END DO ! DO k = 0, C%SELEN_irreg_time_n
          IF (par%master .AND. C%SELEN_display_progress) WRITE(0,'(A)',ADVANCE='NO') CHAR(13)
          
          KKKK( C%SELEN_j1:C%SELEN_j2,:) = KKKK( C%SELEN_j1:C%SELEN_j2,:) * NP_INV
          CALL sync
          
          ! Solving for arrays 'Z' and 'S' 
          ! Only comput S(k=TCONV), solve S(k=TCONV) similat to the real solution, interpolation is done
          IF (C%SELEN_use_rotational_feedback) THEN
            DO j = C%SELEN_j1, C%SELEN_j2
              Z( j,:,iters) = ZE( j,:) + HHHH(     j,:) + KKKK(     j,:) + ZROT_MOD( j,:) ! added rotational feedback terms
              S( j,:      ) = SE( j,:) + AAAA_MOD( j,:) + BBBB_MOD( j,:) + SROT_MOD( j,:) ! added rotational feedback terms
            END DO
            CALL sync
          ELSE
            DO j = C%SELEN_j1, C%SELEN_j2
              Z( j,:,iters) = ZE( j,:) + HHHH(     j,:) + KKKK(     j,:)
              S( j,:      ) = SE( j,:) + AAAA_MOD( j,:) + BBBB_MOD( j,:)
            END DO
            CALL sync
          END IF ! IF (C%SELEN_use_rotational_feedback) THEN
          
!          ! Additional Srot and Zrot calculation needed when the number of iterations is larger then 1
!          IF (iters > 1 .AND. C%SELEN_use_rotational_feedback) THEN
!      
!            IF (par%master) WRITE(0,*) 'ERROR - rotational feedback isnt parallelised and also it doesnt work!'
!            STOP
!            
!            ! Update the total loading of ice and water
!            load_ice_and_water = ice_density * IIII + seawater_density * Z(:,:,iters)
!            CALL s_rotaz(load_ice_and_water,srot)
!            
!            ! Computing S^{rot} at all pixels [Note: Using variable "x" in order to save memory] 
!            DO i = 1, SELEN%mesh%nV
!              DO k = 0, C%SELEN_irreg_time_n
!                x(i,k) = 0d0
!                DO j = 1, C%SELEN_jmax
!                  x(i,k) = x(i,k) + ALF(j,SELEN_anc(i))*dm(j)*dble(SROT(j,k)*long_table(MJ_VAL(j),i))
!                END DO ! DO j = 1, C%SELEN_jmax
!              END DO ! DO k = 0, C%SELEN_irreg_time_n
!            END DO ! DO i = 1, SELEN%mesh%nV
!            
!            ! Computing the harmonics of Z^{rot} by integration over wet pixels 
!            ! Computing the ocean average of S^{rot}
!            ! Computing the ocean average of Z^{rot}
!            zrot(:,:)=0d0
!            DO k = 0, C%SELEN_irreg_time_n
!              srotvv(k) = 0
!              zrotvv(k) = 0.
!              DO j = 1, C%SELEN_jmax
!                zrot(j,k)=0d0
!                DO i = 1, SELEN%mesh%nV
!                  IF (WET(i,k)==1) zrot(j,k) = zrot(j,k) + x(i,k)*alf(j,SELEN_anc(i))*conjg(long_table(MJ_VAL(j),i))
!                END DO ! DO i = 1, SELEN%mesh%nV
!                srotvv(k) = srotvv(k) + dm(j)*dble(varoc(j,k)*conjg(srot(j,k))) * varoc_inv(k)
!                zrotvv(k) = zrotvv(k) + dm(j)*dble(varoc(j,k)*conjg(zrot(j,k))) * varoc_inv(k)
!              END DO ! DO j = 1, C%SELEN_jmax
!            END DO ! DO k = 0, C%SELEN_irreg_time_n
!            zrot(:,:)=zrot(:,:) * NP_INV
!            
!            ! Computing the modified array S^{rot}
!            srot_mod      = srot
!            srot_mod(1,:) = srot(1,:) - srotvv(:)
!            
!            ! Computing the modified array Z^{rot}
!            zrot_mod      = zrot
!            zrot_mod(1,:) = zrot(1,:) - zrotvv(:)
!            
!          END IF ! IF (iters > 1 .AND. C%SELEN_use_rotational_feedback) THEN
    
        END DO ! DO iters = 1, C%SELEN_n_recursion_iterations

        IF (par%master) WRITE (0,'(A)') '    Finished the recursion'
    
    ! =====================
    ! == Recursion - end ==
    ! =====================

      END IF ! IF ( (.NOT. C%use_selen_with_eustatic) .OR. C%SELEN_n_recursion_iterations /= 0) THEN

      ! Array "B" for vertical displacement from water loading
      ! Array "A" for vertical displacement from ice loading
      IF (par%master) WRITE (0,'(A)') '    Computing vertical displacement...'
      progprev = -1
        
      BBBB( C%SELEN_j1:C%SELEN_j2,:) = (0.,0.)
      AAAA( C%SELEN_j1:C%SELEN_j2,:) = (0.,0.)
      DO k = 0,C%SELEN_irreg_time_n
            
        prog = FLOOR(100*REAL(k)/REAL(C%SELEN_irreg_time_n))
        IF (par%master .AND. prog>progprev .AND. C%SELEN_display_progress) THEN
          progprev = prog
          WRITE(0,'(A,I3,A)',ADVANCE='NO') CHAR(13)//'      Progress: ', prog, ' %'
        END IF
            
        DO j = C%SELEN_j1, C%SELEN_j2
          BBBB(j,k) = EU(LJ_VAL(j)) * Z(j,k,C%SELEN_n_recursion_iterations)
          AAAA(j,k) = EU(LJ_VAL(j)) * IIII(j,k)
          DO p = 0,k
            IF     (p == 0 .AND. k == 0) THEN
              BBBB( j,k) = BBBB( j,k) - (Z  (  j,p,C%SELEN_n_recursion_iterations)                                              ) * BETAU( LJ_VAL( j),0)
              AAAA( j,k) = AAAA( j,k) - (IIII( j,p)                                                                             ) * BETAU( LJ_VAL( j),0)
            ELSEIF (p == 0 .AND. k /= 0) THEN
              BBBB( j,k) = BBBB( j,k) - (Z(    j,p,C%SELEN_n_recursion_iterations)                                              ) * BETAU( LJ_VAL( j),INT(SUM(C%SELEN_irreg_time_window(p+1:k))))
              AAAA( j,k) = AAAA( j,k) - (IIII( j,p)                                                                             ) * BETAU( LJ_VAL( j),INT(SUM(C%SELEN_irreg_time_window(p+1:k))))
            ELSEIF (p >  0 .AND. p <  k) THEN
              BBBB( j,k) = BBBB( j,k) - (Z(    j,p,C%SELEN_n_recursion_iterations) - Z(    j,p-1,C%SELEN_n_recursion_iterations)) * BETAU( LJ_VAL( j),INT(SUM(C%SELEN_irreg_time_window(p+1:k))))
              AAAA( j,k) = AAAA( j,k) - (IIII( j,p)                                - IIII( j,p-1)                               ) * BETAU( LJ_VAL( j),INT(SUM(C%SELEN_irreg_time_window(p+1:k))))
            ELSE ! p==k
              BBBB( j,k) = BBBB( j,k) - (Z(    j,p,C%SELEN_n_recursion_iterations) - Z(    j,p-1,C%SELEN_n_recursion_iterations)) * BETAU( LJ_VAL( j),0) 
              AAAA( j,k) = AAAA( j,k) - (IIII( j,p)                                - IIII( j,p-1)                               ) * BETAU( LJ_VAL( j),0)
            END IF
          END DO ! DO p = 0,k
          BBBB( j,k) = RHOW_o_RHOE_X3 * BBBB( j,k)
          AAAA( j,k) = RHOI_o_RHOE_X3 * AAAA( j,k)
        END DO ! DO j = 1, C%SELEN_jmax
        CALL sync
      END DO ! DO k = 0,C%SELEN_irreg_time_n
      IF (par%master .AND. C%SELEN_display_progress) WRITE(0,'(A)',ADVANCE='NO') CHAR(13)
      
      ! Vertical displacement
      U( C%SELEN_j1:C%SELEN_j2,:) = AAAA( C%SELEN_j1:C%SELEN_j2,:) + BBBB( C%SELEN_j1:C%SELEN_j2,:)
      
      ! For each new iteration calcualte set S in MEM_S
      IF (C%SELEN_n_TDOF_iterations > 1 .AND. tdof_iteration /= C%SELEN_n_TDOF_iterations) THEN
        MEM_S( C%SELEN_j1:C%SELEN_j2,:) = S( C%SELEN_j1:C%SELEN_j2,:)
      END IF
      
    END DO ! DO tdof_iteration = 1, C%SELEN_n_TDOF_iterations

  ! =================================
  ! == End of TDOF-NFWD iterations ==
  ! =================================
  
    IF (par%master) WRITE(0,'(A,I3,A)') '    SLE solved in ', C%SELEN_n_TDOF_iterations, ' iterations'
        
    ! Interpolate MEM_S and MEM_U to be used in the next time step, although the last stays the current one.. the rest shifts 1 interval
    
    ! Interpolation as a function of the CALL time step (C%dt_bedrock in years)
    MSint( C%SELEN_j1:C%SELEN_j2,0) = S( C%SELEN_j1:C%SELEN_j2,0)
    MUint( C%SELEN_j1:C%SELEN_j2,0) = U( C%SELEN_j1:C%SELEN_j2,0)
    
    DO k = 1, C%SELEN_irreg_time_n-1
    
      ! find the corresponding range of the original data for the next time step k+1
      Tbeg = -SUM(C%SELEN_irreg_time_window(k  :C%SELEN_irreg_time_n-1))    ! time of 1st time point  
      Tend = -SUM(C%SELEN_irreg_time_window(k+1:C%SELEN_irreg_time_n-1))    ! time of 2nd time point
      Tint = Tbeg + (C%dt_SELEN / 1000._dp)        ! time of interpolated time point, as function of the CALL interval
      
      ! The previous time point to interpolate from (if DELTA = dt_bedrock, set the 1st and 2nd time step the same)
      MSint( C%SELEN_j1:C%SELEN_j2,k) = S( C%SELEN_j1:C%SELEN_j2,k) + (Tint - Tbeg) * (S( C%SELEN_j1:C%SELEN_j2,k+1) - S( C%SELEN_j1:C%SELEN_j2,k)) / C%SELEN_irreg_time_window(k)
      MUint( C%SELEN_j1:C%SELEN_j2,k) = U( C%SELEN_j1:C%SELEN_j2,k) + (Tint - Tbeg) * (U( C%SELEN_j1:C%SELEN_j2,k+1) - U( C%SELEN_j1:C%SELEN_j2,k)) / C%SELEN_irreg_time_window(k)
      
      ! Forward extrapolation when past time is not the same..
      MSint( C%SELEN_j1:C%SELEN_j2,k) = S( C%SELEN_j1:C%SELEN_j2,k) + (Tint - Tbeg) * (S( C%SELEN_j1:C%SELEN_j2,k+1) - S( C%SELEN_j1:C%SELEN_j2,k)) / C%SELEN_irreg_time_window(k)
      MUint( C%SELEN_j1:C%SELEN_j2,k) = U( C%SELEN_j1:C%SELEN_j2,k) + (Tint - Tbeg) * (U( C%SELEN_j1:C%SELEN_j2,k+1) - U( C%SELEN_j1:C%SELEN_j2,k)) / C%SELEN_irreg_time_window(k)
      
    END DO ! DO k = 1, C%SELEN_irreg_time_n-1
  
    MSint( C%SELEN_j1:C%SELEN_j2,C%SELEN_irreg_time_n) = S( C%SELEN_j1:C%SELEN_j2,C%SELEN_irreg_time_n)
    MUint( C%SELEN_j1:C%SELEN_j2,C%SELEN_irreg_time_n) = U( C%SELEN_j1:C%SELEN_j2,C%SELEN_irreg_time_n)
    
    !  Store interpolated S and U in memory
    MEM_S( C%SELEN_j1:C%SELEN_j2,:) = MSint( C%SELEN_j1:C%SELEN_j2,:)
    MEM_U( C%SELEN_j1:C%SELEN_j2,:) = MUint( C%SELEN_j1:C%SELEN_j2,:)
    CALL sync
    
    ! Output for global relative sea-level S, new bedrock topography and ocean function
    IF (par%master) WRITE (0,'(A)') '    Calculating global fields...'
    S_global(          C%SELEN_i1:C%SELEN_i2)  = 0._dp
    U_global(          C%SELEN_i1:C%SELEN_i2)  = 0._dp
    SELEN%Hi_glob(     C%SELEN_i1:C%SELEN_i2) = 0._dp
    SELEN%Hi_rel_glob( C%SELEN_i1:C%SELEN_i2) = 0._dp
    SELEN%U_glob(      C%SELEN_i1:C%SELEN_i2) = 0._dp
    SELEN%N_glob(      C%SELEN_i1:C%SELEN_i2) = 0._dp
    SELEN%of_glob(     C%SELEN_i1:C%SELEN_i2) = 0._dp
    
    DO i = C%SELEN_i1, C%SELEN_i2
      DO j = 1, C%SELEN_jmax
        S_global(          i) = S_global(          i)  + newalf( j,SELEN_anc( i)) * dble( MEM_S( j,C%SELEN_irreg_time_n) * LONG_TABLE( MJ_VAL( j),i))
        U_global(          i) = U_global(          i)  + newalf( j,SELEN_anc( i)) * dble( MEM_U( j,C%SELEN_irreg_time_n) * LONG_TABLE( MJ_VAL( j),i))
        SELEN%Hi_glob(     i) = SELEN%Hi_glob(     i)  + newalf( j,SELEN_anc( i)) * dble( TTTT(  j,C%SELEN_irreg_time_n) * LONG_TABLE( MJ_VAL( j),i))
        SELEN%Hi_rel_glob( i) = SELEN%Hi_rel_glob( i)  + newalf( j,SELEN_anc( i)) * dble( IIII(  j,C%SELEN_irreg_time_n) * LONG_TABLE( MJ_VAL( j),i))
      END DO ! DO j=1, C%SELEN_jmax
      
      SELEN%N_glob(  i) = S_global( i) + U_global( i)
      SELEN%U_glob(  i) = U_global( i)
      SELEN%of_glob( i) = WET( i, C%SELEN_irreg_time_n)
      
    END DO
    CALL sync
    
    DO k = 0, C%SELEN_irreg_time_n
      Seust( k) = 0._dp
      DO i = C%SELEN_i1, C%SELEN_i2
        IF (WET(i,k) == 1) THEN
          DO j = 1, C%SELEN_jmax
            Seust( k) = Seust( k) + newalf(j,SELEN_anc(i)) * dble( MEM_S(j,k) * LONG_TABLE(MJ_VAL(j),i))
          END DO ! DO j = 1, C%SELEN_jmax
        END IF ! IF (WET(i,k) == 1) THEN
      END DO ! DO i = 1, SELEN%mesh%nV
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, Seust( k), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    END DO ! DO k = 0, C%SELEN_irreg_time_n
    
    IF (par%master) THEN
    DO k = 0, C%SELEN_irreg_time_n
      Seust(k) = Seust(k) / DBLE(COUNT(WET(:,k)==1))
    END DO ! DO k = 0, C%SELEN_irreg_time_n
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Update px-table and save new topography and new land-ocean mask (newet2).
    CALL nextopo( SELEN, icethick, INIT_TOPO, S_global, newtopo2, newet2)
    
    ! Current change in ice loading, eustatic change
    Tw1    = -RHOI_O_RHOW * SUM(icethick(1:SELEN%mesh%nV,C%SELEN_irreg_time_n) - icethick(1:SELEN%mesh%nV,C%SELEN_irreg_time_n-1)) / DBLE(COUNT(newet2 == 1))    
    ! Total Eustatic sea level change
    Seust2 = SUM(S_global, MASK = (newet2 == 1)) / DBLE(COUNT(newet2 == 1))
    ! Total change in ice loading, eustatic change
    Tw2    = -RHOI_O_RHOW * SUM(icethick(1:SELEN%mesh%nV,C%SELEN_irreg_time_n) - INIT_ICE(1:SELEN%mesh%nV)  ) / DBLE(COUNT(newet2 == 1))

    ocean_area  = COUNT(newet2 == 1) * 2._dp * pi * earth_radius**2 * (1._dp - COS(C%SELEN_alfa))
    ocean_depth = (SUM(S_global - INIT_TOPO, MASK = (newet2 == 1)) * 2._dp * pi * earth_radius**2 * (1._dp - COS(C%SELEN_alfa))) / ocean_area
  
    CALL MPI_BCAST( ocean_area,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( ocean_depth, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Clean up after yourself
    CALL deallocate_memory_SELEN( SELEN)
    IF (.NOT. par%master) THEN
      DEALLOCATE( BETAS)
      DEALLOCATE( BETAU)
      DEALLOCATE( ES   )
      DEALLOCATE( EU   )
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_SLE
  
  SUBROUTINE futopos( SELEN, icethick, topo, slc, newet, fj, newtopo)
    ! Generate (0:C%SELEN_irreg_time_n) topographies from slc that is calculated at the previous CALL
    
    IMPLICIT NONE
 
    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(INOUT) :: SELEN
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(IN)    :: icethick
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(IN)    :: topo
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(IN)    :: slc
    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: newet
    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: fj
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: newtopo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'futopos'
    INTEGER                                            :: i,k,oj
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    DO i = C%SELEN_i1, C%SELEN_i2
    DO k = 0, C%SELEN_irreg_time_n
    
      newtopo( i,k) = topo( i) - slc( i,k)
      
      IF (newtopo( i,k) <= 0) THEN
        oj = 1
      ELSE
        oj = 0
      END IF
      
      IF (-newtopo( i,k) * seawater_density >  icethick( i,k) * ice_density) THEN
        fj( i,k) = 1
      ELSE
        fj( i,k) = 0
      END IF
      
      newet( i,k) = oj * fj( i,k)
      
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE futopos
  SUBROUTINE nextopo( SELEN, icethick, topo, slc, newtopo2, newet2)
    !  Same as futopos buth only for the current time step
    
    IMPLICIT NONE
 
    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(INOUT) :: SELEN
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(IN)    :: icethick
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(IN)    :: topo
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(IN)    :: slc
    INTEGER,  DIMENSION(:    ), POINTER, INTENT(INOUT) :: newet2
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(INOUT) :: newtopo2
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'nextopo'
    INTEGER                                            :: i,oj,fj
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    DO i = C%SELEN_i1, C%SELEN_i2
      
      newtopo2( i) = topo( i) - slc( i)
      
      IF (newtopo2( i) <= 0._dp) THEN
        oj = 1
      ELSE
        oj = 0
      END IF
      
      IF (-newtopo2( i) * seawater_density > icethick( i,C%SELEN_irreg_time_n) * ice_density) THEN
        fj = 1
      ELSE
        fj = 0
      END IF
      
      newet2( i) = oj * fj
      
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE nextopo
  
  SUBROUTINE update_ice_sh_and_load( SELEN, NAM, EAS, GRL, ANT, icethick, TTTT, IIII)
    ! Computes the spherical harmonics transform IIII of the ice loading history.
    
    IMPLICIT NONE
 
    ! In/output variables:
    TYPE(type_SELEN_global),               INTENT(INOUT) :: SELEN
    TYPE(type_model_region),               INTENT(INOUT) :: NAM, EAS, GRL, ANT
    REAL(dp),   DIMENSION(:,:  ), POINTER, INTENT(IN)    :: icethick
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: TTTT
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: IIII
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_ice_sh_and_load'
    INTEGER                                              :: ki 
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Compute spherical harmonics transformation of the ice loading history
    CALL sh_transform_all_regions_3D( NAM, EAS, GRL, ANT, icethick, TTTT)
    CALL sync
    
    ! Calculate "relative ice loading history"
    IIII( C%SELEN_j1:C%SELEN_j2,:) = (0.,0.)
    DO ki = 0, C%SELEN_irreg_time_n
      IF (ki < C%SELEN_irreg_time_n) THEN
        IIII( C%SELEN_j1:C%SELEN_j2,ki) = TTTT( C%SELEN_j1:C%SELEN_j2,ki+1) - TTTT( C%SELEN_j1:C%SELEN_j2,0) ! Is new time points minus reference (increasing k with older points
      ELSE
        IIII( C%SELEN_j1:C%SELEN_j2,ki) = TTTT( C%SELEN_j1:C%SELEN_j2,ki  ) - TTTT( C%SELEN_j1:C%SELEN_j2,0)
      END IF
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE update_ice_sh_and_load
  
  SUBROUTINE sh_transform_all_regions_2D( NAM, EAS, GRL, ANT, d, sh)
    ! Computes the spherical harmonics transform sh of a data field d (defined on the
    ! SELEN global grid, but assumed to have non-zero values only on pixels lying
    ! inside UFEMISM regions).
    ! Can be parallelised, is not right now.
    
    IMPLICIT NONE
 
    ! In/output variables:
    TYPE(type_model_region),               INTENT(INOUT) :: NAM, EAS, GRL, ANT
    REAL(dp),   DIMENSION(:    ), POINTER, INTENT(IN)    :: d
    COMPLEX*16, DIMENSION(:    ), POINTER, INTENT(INOUT) :: sh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'sh_transform_all_regions_2D'
    COMPLEX*16, DIMENSION( C%SELEN_jmax)                   :: sh_core
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    sh_core = (0.,0.)
    
    IF (C%do_NAM) CALL sh_transform_single_region_2D( NAM, d, sh_core)
    IF (C%do_EAS) CALL sh_transform_single_region_2D( EAS, d, sh_core)
    IF (C%do_GRL) CALL sh_transform_single_region_2D( GRL, d, sh_core)
    IF (C%do_ANT) CALL sh_transform_single_region_2D( ANT, d, sh_core)
    
    CALL MPI_REDUCE( sh_core, sh, C%SELEN_jmax, MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE sh_transform_all_regions_2D
  SUBROUTINE sh_transform_all_regions_3D( NAM, EAS, GRL, ANT, d, sh)
    ! Computes the spherical harmonics transform sh of a data field d (defined on the
    ! SELEN global grid, but assumed to have non-zero values only on pixels lying
    ! inside UFEMISM regions).
    ! Can be parallelised, is not right now.
    
    IMPLICIT NONE
 
    ! In/output variables:
    TYPE(type_model_region),               INTENT(INOUT) :: NAM, EAS, GRL, ANT
    REAL(dp),   DIMENSION(:,:  ), POINTER, INTENT(IN)    :: d
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: sh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'sh_transform_all_regions_3D'
    COMPLEX*16, DIMENSION( C%SELEN_jmax, 0:C%SELEN_irreg_time_n)    :: sh_core
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    sh_core = (0.,0.)
    
    IF (C%do_NAM) CALL sh_transform_single_region_3D( NAM, d, sh_core)
    IF (C%do_EAS) CALL sh_transform_single_region_3D( EAS, d, sh_core)
    IF (C%do_GRL) CALL sh_transform_single_region_3D( GRL, d, sh_core)
    IF (C%do_ANT) CALL sh_transform_single_region_3D( ANT, d, sh_core)
    
    CALL MPI_REDUCE( sh_core, sh, C%SELEN_jmax * (C%SELEN_irreg_time_n+1), MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE sh_transform_all_regions_3D
  SUBROUTINE sh_transform_single_region_2D( region, d, sh_core)
    ! Here, the spherical harmonics are read from the external files and the actual transform is done.
    
    IMPLICIT NONE
 
    ! In/output variables:
    TYPE(type_model_region),               INTENT(INOUT) :: region
    REAL(dp),   DIMENSION(:    ), POINTER, INTENT(IN)    :: d
    COMPLEX*16, DIMENSION( C%SELEN_jmax),  INTENT(INOUT) :: sh_core
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'sh_transform_single_region_2D'
    CHARACTER(LEN=256)                                   :: sh_filename
    COMPLEX*16, DIMENSION( C%SELEN_jmax)                 :: Ypx
    INTEGER                                              :: nblocks, bi, b1, b2, ir, ir1, ir2, i
    
    ! Add routine to path
    CALL init_routine( routine_name)
  
    ! Determine how many blocks of 100 pixels there are, divide those over the processes
    nblocks = CEILING( REAL(region%SELEN%nr) / 100._dp)
    CALL partition_list( nblocks, par%i, par%n, b1, b2)
    
    ! Go over all the blocks assigned to this process
    DO bi = b1, b2
    
      ! Determine file name
      sh_filename = ''
      IF     (bi < 10)    THEN
        WRITE( sh_filename,'(A,A,I1,A)') TRIM(region%SELEN%sh_foldername), '/sh_glob_000', bi, '.bin'
      ELSEIF (bi < 100)   THEN
        WRITE( sh_filename,'(A,A,I2,A)') TRIM(region%SELEN%sh_foldername), '/sh_glob_00',  bi, '.bin'
      ELSEIF (bi < 1000)  THEN
        WRITE( sh_filename,'(A,A,I3,A)') TRIM(region%SELEN%sh_foldername), '/sh_glob_0',   bi, '.bin'
      ELSEIF (bi < 10000) THEN
        WRITE( sh_filename,'(A,A,I4,A)') TRIM(region%SELEN%sh_foldername), '/sh_glob_',    bi, '.bin'
      END IF
      
      ! Open the binary file
      OPEN( UNIT = 1337+par%i, FILE = sh_filename, STATUS = 'OLD', FORM = 'UNFORMATTED')
    
      ! Calculate and write spherical harmonics for all pixels in this block
      ir1 = (bi-1)*100+1
      ir2 = MIN( bi*100, region%SELEN%nr)
      DO ir = ir1, ir2
      
        ! Look up this pixel's global index
        i = region%SELEN%map_isl_region2glob( ir)
      
        ! Read
        READ( UNIT = 1337+par%i) Ypx
        
        ! Calculate
        sh_core = sh_core + d( i) * Ypx
        
      END DO ! DO ir = ir1, ir2
      
      ! Close the binary file
      CLOSE( UNIT = 1337+par%i)
    
    END DO ! DO bi = b1, b2
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE sh_transform_single_region_2D
  SUBROUTINE sh_transform_single_region_3D( region, d, sh_core)
    ! Here, the spherical harmonics are read from the external files and the actual transform is done.
    
    IMPLICIT NONE
 
    ! In/output variables:
    TYPE(type_model_region),               INTENT(INOUT) :: region
    REAL(dp),   DIMENSION(:,:  ), POINTER, INTENT(IN)    :: d
    COMPLEX*16, DIMENSION( C%SELEN_jmax, 0:C%SELEN_irreg_time_n), INTENT(INOUT)   :: sh_core
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'sh_transform_single_region_3D'
    INTEGER                                              :: ki
    CHARACTER(LEN=256)                                   :: sh_filename
    COMPLEX*16, DIMENSION( C%SELEN_jmax)                 :: Ypx
    INTEGER                                              :: nblocks, bi, b1, b2, ir, ir1, ir2, i
    
    ! Add routine to path
    CALL init_routine( routine_name)
  
    ! Determine how many blocks of 100 pixels there are, divide those over the processes
    nblocks = CEILING( REAL(region%SELEN%nr) / 100._dp)
    CALL partition_list( nblocks, par%i, par%n, b1, b2)
    
    ! Go over all the blocks assigned to this process
    DO bi = b1, b2
    
      ! Determine file name
      sh_filename = ''
      IF     (bi < 10)    THEN
        WRITE( sh_filename,'(A,A,I1,A)') TRIM(region%SELEN%sh_foldername), '/sh_glob_000', bi, '.bin'
      ELSEIF (bi < 100)   THEN
        WRITE( sh_filename,'(A,A,I2,A)') TRIM(region%SELEN%sh_foldername), '/sh_glob_00',  bi, '.bin'
      ELSEIF (bi < 1000)  THEN
        WRITE( sh_filename,'(A,A,I3,A)') TRIM(region%SELEN%sh_foldername), '/sh_glob_0',   bi, '.bin'
      ELSEIF (bi < 10000) THEN
        WRITE( sh_filename,'(A,A,I4,A)') TRIM(region%SELEN%sh_foldername), '/sh_glob_',    bi, '.bin'
      END IF
      
      ! Open the binary file
      OPEN( UNIT = 1337+par%i, FILE = sh_filename, STATUS = 'OLD', FORM = 'UNFORMATTED')
    
      ! Calculate and write spherical harmonics for all pixels in this block
      ir1 = (bi-1)*100+1
      ir2 = MIN( bi*100, region%SELEN%nr)
      DO ir = ir1, ir2
      
        ! Look up this pixel's global index
        i = region%SELEN%map_isl_region2glob( ir)
      
        ! Read
        READ( UNIT = 1337+par%i) Ypx
        
        ! Calculate
        DO ki = 0, C%SELEN_irreg_time_n
          sh_core( :,ki) = sh_core( :,ki) + d( i,ki) * Ypx
        END DO
        
      END DO ! DO ir = ir1, ir2
      
      ! Close the binary file
      CLOSE( UNIT = 1337+par%i)
    
    END DO ! DO bi = b1, b2
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE sh_transform_single_region_3D
  
  SUBROUTINE inverse_sh_transform_single_region_2D( region, sh, d)
    ! Performs the inverse spherical harmonical transformation, converting
    ! the spherical harmonics data fields sh back to the gridded data field d
    ! (on a regional square model grid)
    
    IMPLICIT NONE
 
    ! In/output variables:
    TYPE(type_model_region),               INTENT(INOUT) :: region
    COMPLEX*16, DIMENSION(:    ), POINTER, INTENT(IN)    :: sh
    REAL(dp),   DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: d
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inverse_sh_transform_single_region_2D'
    INTEGER                                              :: i,j
    CHARACTER(LEN=256)                                   :: sh_filename
    COMPLEX*16, DIMENSION( C%SELEN_jmax)                 :: Ypx
    INTEGER                                              :: li
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Go over all the columns assigned to this process
    DO i = region%grid_GIA%i1, region%grid_GIA%i2
    
      ! Determine file name
      sh_filename = ''
      IF     (i < 10)    THEN
        WRITE( sh_filename,'(A,A,I1,A)') TRIM(region%SELEN%sh_foldername), '/sh_square_000', i, '.bin'
      ELSEIF (i < 100)   THEN
        WRITE( sh_filename,'(A,A,I2,A)') TRIM(region%SELEN%sh_foldername), '/sh_square_00',  i, '.bin'
      ELSEIF (i < 1000)  THEN
        WRITE( sh_filename,'(A,A,I3,A)') TRIM(region%SELEN%sh_foldername), '/sh_square_0',   i, '.bin'
      ELSEIF (i < 10000) THEN
        WRITE( sh_filename,'(A,A,I4,A)') TRIM(region%SELEN%sh_foldername), '/sh_square_',    i, '.bin'
      END IF
      
      ! Open the binary file
      OPEN( UNIT = 1337+par%i, FILE = sh_filename, STATUS = 'OLD', FORM = 'UNFORMATTED')
      
      ! Calculate and write spherical harmonics for all rows in this column
      DO j = 1, region%grid_GIA%ny
           
        ! Read spherical harmonics for this GIA grid pixel
        READ( UNIT = 1337+par%i) Ypx
      
        ! Compute inverse transform
        d( j,i) = 0._dp
        DO li = 1, C%SELEN_jmax
          d( j,i) = d( j,i) + dble( sh( li) * Ypx( li)) 
        END DO ! DO li = 1, C%SELEN_jmax
        
      END DO ! DO j = 1, region%grid_GIA%ny
      
      ! Close the binary file
      CLOSE( UNIT = 1337+par%i)
    
    END DO ! DO i = region%grid_GIA%i1, region%grid_GIA%i2
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inverse_sh_transform_single_region_2D
  
  SUBROUTINE allocate_memory_and_bind_pointers( SELEN, DM, X, ZE, SE, AAAA, AAAA_MOD, BBBB, BBBB_MOD, HHHH, KKKK, TTTT, IIII, &
    load_ice_and_water, ZROT, ZROT_MOD, SROT, SROT_MOD, ZROTVV, SROTVV, S, U, Z, newalf, slc, WET, newtopo, fj, xcf, &
    newet2, newtopo2, S_global, U_global, MSint, MUint, MSread, int_time, sp1, spn, preSint, up1, upn, preUint, s2,u2, &
    varpreoc, varoc, varoc_inv, MEM_S, MEM_U, icethick, MJ_VAL, LJ_VAL, SELEN_anc, INIT_TOPO, ALF, LONG_TABLE, INIT_ICE)
    ! Allocate shared memory for SELEN internal data and bind local pointers to it.
  
    IMPLICIT NONE  
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_memory_and_bind_pointers'
    TYPE(type_SELEN_global),             INTENT(INOUT) :: SELEN
    INTEGER,    DIMENSION(:    ), POINTER, INTENT(OUT) :: DM                  ! INTEGER,    DIMENSION( C%SELEN_jmax                ) :: DM
    REAL(dp),   DIMENSION(:,:  ), POINTER, INTENT(OUT) :: X                   ! REAL(dp),   DIMENSION( SELEN%mesh%nV,  0:C%SELEN_irreg_time_n  ) :: X
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: ZE                  ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: ZE                  ! Eustatic Z array for water loading
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: SE                  ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: SE                  ! Eustatic S array for ice loading
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: AAAA                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: AAAA
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: AAAA_MOD            ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: AAAA_MOD
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: BBBB                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: BBBB
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: BBBB_MOD            ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: BBBB_MOD
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: HHHH                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: HHHH
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: KKKK                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: KKKK
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: TTTT                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: TTTT                ! Relative ice loading in SH
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: IIII                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: IIII                ! Relative ice loading in SH
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: load_ice_and_water  ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: load_ice_and_water  ! Ice and water loading in SH for rotational feedback
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: ZROT                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: ZROT
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: ZROT_MOD            ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: ZROT_MOD
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: SROT                ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: SROT
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: SROT_MOD            ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: SROT_MOD
    COMPLEX*16, DIMENSION(:    ), POINTER, INTENT(OUT) :: ZROTVV              ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: ZROTVV
    COMPLEX*16, DIMENSION(:    ), POINTER, INTENT(OUT) :: SROTVV              ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: SROTVV
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: S                   ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: S                   ! relative sea level
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: U                   ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: U                   ! Solid Earth deformation (N = S+U) 
    COMPLEX*16, DIMENSION(:,:,:), POINTER, INTENT(OUT) :: Z                   ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n,0:C%SELEN_n_recursion_iterations) :: Z
    REAL(dp),   DIMENSION(:,:  ), POINTER, INTENT(OUT) :: newalf              ! REAL(dp),   DIMENSION( C%SELEN_jmax,  SELEN%mesh%nanc  ) :: newalf              ! New values of alf
    REAL(dp),   DIMENSION(:,:  ), POINTER, INTENT(OUT) :: slc                 ! REAL(dp),   DIMENSION( SELEN%mesh%nV,  0:C%SELEN_irreg_time_n  ) :: slc                 ! sea level forward in time for all elements
    INTEGER,    DIMENSION(:,:  ), POINTER, INTENT(OUT) :: WET                 ! INTEGER,    DIMENSION( SELEN%mesh%nV,  0:C%SELEN_irreg_time_n  ) :: WET                 ! new ocean function (0 or 1) forward in time
    REAL(dp),   DIMENSION(:,:  ), POINTER, INTENT(OUT) :: newtopo             ! REAL(dp),   DIMENSION( SELEN%mesh%nV,  0:C%SELEN_irreg_time_n  ) :: newtopo             ! new topograhpy forward in time
    INTEGER,    DIMENSION(:,:  ), POINTER, INTENT(OUT) :: fj                  ! INTEGER,    DIMENSION( SELEN%mesh%nV,  0:C%SELEN_irreg_time_n  ) :: fj                  ! floating function for ice thickness (float or not)
    INTEGER,    DIMENSION(:,:  ), POINTER, INTENT(OUT) :: xcf                 ! INTEGER,    DIMENSION( SELEN%mesh%nV,  0:C%SELEN_irreg_time_n+1) :: xcf
    INTEGER,    DIMENSION(:    ), POINTER, INTENT(OUT) :: newet2              ! INTEGER,    DIMENSION( SELEN%mesh%nV                  ) :: newet2              ! new ocean function after iteration
    REAL(dp),   DIMENSION(:    ), POINTER, INTENT(OUT) :: newtopo2            ! REAL(dp),   DIMENSION( SELEN%mesh%nV                  ) :: newtopo2            ! new topograhpy after iteration
    REAL(dp),   DIMENSION(:    ), POINTER, INTENT(OUT) :: S_global, U_global  ! REAL(dp),   DIMENSION( SELEN%mesh%nV                  ) :: S_global, U_global  ! Global fields of S and U for output
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: MSint               ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: MSint               ! Interpolated memory of relative sea level
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: MUint               ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: MUint               ! Interpolationd memory of solid earth deformation (N is geoid)
    COMPLEX*16, DIMENSION(:    ), POINTER, INTENT(OUT) :: MSread              ! COMPLEX*16, DIMENSION( C%SELEN_jmax                ) :: MSread              ! Initial memory of relative sea level
    REAL(dp),   DIMENSION(:    ), POINTER, INTENT(OUT) :: int_time            ! REAL(dp),   DIMENSION(            0:C%SELEN_irreg_time_n  ) :: int_time            ! time points (depending on DELTA) for interpolation
    COMPLEX*16, DIMENSION(:    ), POINTER, INTENT(OUT) :: sp1, spn, preSint   ! COMPLEX*16, DIMENSION( C%SELEN_jmax)                 :: sp1, spn, preSint   ! spline in- and output for S
    COMPLEX*16, DIMENSION(:    ), POINTER, INTENT(OUT) :: up1, upn, preUint   ! COMPLEX*16, DIMENSION( C%SELEN_jmax)                 :: up1, upn, preUint   ! spline in- and output for U
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: s2,u2               ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: s2,u2               ! spline output
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: varpreoc            ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: varpreoc            ! pre determined SH coefficients ocean function
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: varoc               ! COMPLEX*16, DIMENSION( C%SELEN_jmax,0:C%SELEN_irreg_time_n  ) :: varoc               ! variable SH coefficient of ocean function
    COMPLEX*16, DIMENSION(:    ), POINTER, INTENT(OUT) :: varoc_inv           ! COMPLEX*16, DIMENSION(            0:C%SELEN_irreg_time_n  ) :: varoc_inv           ! 1 / varoc(1,:) - the first degree
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: MEM_S               ! COMPLEX*16, DIMENSION(:,:), ALLOCATABLE, SAVE :: MEM_S
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: MEM_U               ! COMPLEX*16, DIMENSION(:,:), ALLOCATABLE       :: MEM_U
    REAL(dp),   DIMENSION(:,:  ), POINTER, INTENT(OUT) :: icethick            ! REAL(dp),   DIMENSION(:,:), ALLOCATABLE       :: icethick
    INTEGER,    DIMENSION(:    ), POINTER, INTENT(OUT) :: MJ_VAL              ! ALLOCATE(MJ_VAL(C%SELEN_jmax))
    INTEGER,    DIMENSION(:    ), POINTER, INTENT(OUT) :: LJ_VAL              ! ALLOCATE(LJ_VAL(C%SELEN_jmax))
    INTEGER,    DIMENSION(:    ), POINTER, INTENT(OUT) :: SELEN_anc           ! INTEGER,    DIMENSION(:  ), ALLOCATABLE, SAVE :: SELEN_anc
    REAL(dp),   DIMENSION(:    ), POINTER, INTENT(OUT) :: INIT_TOPO           ! REAL(dp),   DIMENSION(:  ), ALLOCATABLE, SAVE :: INIT_TOPO
    REAL(dp),   DIMENSION(:,:  ), POINTER, INTENT(OUT) :: ALF                 ! REAL(dp),   DIMENSION(:,:), ALLOCATABLE, SAVE :: ALF
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT) :: LONG_TABLE          ! COMPLEX*16, DIMENSION(:,:), ALLOCATABLE, SAVE :: LONG_TABLE
    REAL(dp),   DIMENSION(:    ), POINTER, INTENT(OUT) :: INIT_ICE            ! REAL(dp),   DIMENSION(:  ), ALLOCATABLE, SAVE :: INIT_ICE
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! First, shared memory is allocated in the SELEN structure.
    ! Then, the pointers (which are variables local to the sle_model) are associated
    ! with this memory. For arrays that contain time memory, which are indexed from 0 in SELEN,
    ! indexing is altered here to make sure the code is unaffected.
    ! For reasons unknown, this only works when the shared memory with which the pointer is
    ! associated is 1-D.
    
    CALL allocate_shared_int_1D(     C%SELEN_jmax                  , SELEN%DM,                 SELEN%wDM                )
    CALL allocate_shared_dp_1D(      SELEN%mesh%nV   * (C%SELEN_irreg_time_n+1), SELEN%X,                  SELEN%wX                 )
    
    DM(                 1:C%SELEN_jmax                 ) => SELEN%DM
    X(                  1:SELEN%mesh%nV,   0:C%SELEN_irreg_time_n  ) => SELEN%X
    
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%ZE,                 SELEN%wZE                )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%SE,                 SELEN%wSE                )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%AAAA,               SELEN%wAAAA              )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%AAAA_MOD,           SELEN%wAAAA_MOD          )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%BBBB,               SELEN%wBBBB              )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%BBBB_MOD,           SELEN%wBBBB_MOD          )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%HHHH,               SELEN%wHHHH              )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%KKKK,               SELEN%wKKKK              )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%TTTT,               SELEN%wTTTT              )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%IIII,               SELEN%wIIII              )
    
    ZE(                 1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%ZE
    SE(                 1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%SE
    AAAA(               1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%AAAA
    AAAA_MOD(           1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%AAAA_MOD
    BBBB(               1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%BBBB
    BBBB_MOD(           1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%BBBB_MOD
    HHHH(               1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%HHHH
    KKKK(               1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%KKKK
    TTTT(               1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%TTTT 
    IIII(               1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%IIII 
    
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%load_ice_and_water, SELEN%wload_ice_and_water)
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%ZROT,               SELEN%wZROT              )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%ZROT_MOD,           SELEN%wZROT_MOD          )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%SROT,               SELEN%wSROT              )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%SROT_MOD,           SELEN%wSROT_MOD          )
    CALL allocate_shared_complex_1D(              (C%SELEN_irreg_time_n+1), SELEN%ZROTVV,             SELEN%wZROTVV            )
    CALL allocate_shared_complex_1D(              (C%SELEN_irreg_time_n+1), SELEN%SROTVV,             SELEN%wSROTVV            )
    
    load_ice_and_water( 1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%load_ice_and_water
    ZROT(               1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%ZROT
    ZROT_MOD(           1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%ZROT_MOD
    SROT(               1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%SROT
    SROT_MOD(           1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%SROT_MOD
    ZROTVV(                           0:C%SELEN_irreg_time_n  ) => SELEN%ZROTVV
    SROTVV(                           0:C%SELEN_irreg_time_n  ) => SELEN%SROTVV
    
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%S,                  SELEN%wS                 )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%U,                  SELEN%wU                 )
    
    S(                  1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%S
    U(                  1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%U
    
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1) * (C%SELEN_n_recursion_iterations+1), SELEN%Z, SELEN%wZ      )
    Z( 1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n, 0:C%SELEN_n_recursion_iterations) => SELEN%Z
    
    CALL allocate_shared_dp_1D(      C%SELEN_jmax *  SELEN%mesh%nanc   , SELEN%newalf,             SELEN%wnewalf            )
    CALL allocate_shared_dp_1D(      SELEN%mesh%nV   * (C%SELEN_irreg_time_n+1), SELEN%slc,                SELEN%wslc               )
    CALL allocate_shared_int_1D(     SELEN%mesh%nV   * (C%SELEN_irreg_time_n+1), SELEN%WET,                SELEN%wWET               )
    CALL allocate_shared_dp_1D(      SELEN%mesh%nV   * (C%SELEN_irreg_time_n+1), SELEN%newtopo,            SELEN%wnewtopo           )
    CALL allocate_shared_int_1D(     SELEN%mesh%nV   * (C%SELEN_irreg_time_n+1), SELEN%fj,                 SELEN%wfj                )
    CALL allocate_shared_int_1D(     SELEN%mesh%nV   * (C%SELEN_irreg_time_n+2), SELEN%xcf,                SELEN%wxcf               )

    newalf(             1:C%SELEN_jmax, 1:SELEN%mesh%nanc  ) => SELEN%newalf
    slc(                1:SELEN%mesh%nV,   0:C%SELEN_irreg_time_n  ) => SELEN%slc
    WET(                1:SELEN%mesh%nV,   0:C%SELEN_irreg_time_n  ) => SELEN%WET
    newtopo(            1:SELEN%mesh%nV,   0:C%SELEN_irreg_time_n  ) => SELEN%newtopo
    fj(                 1:SELEN%mesh%nV,   0:C%SELEN_irreg_time_n  ) => SELEN%fj
    xcf(                1:SELEN%mesh%nV,   0:C%SELEN_irreg_time_n+1) => SELEN%xcf
    
    CALL allocate_shared_int_1D(     SELEN%mesh%nV                    , SELEN%newet2,             SELEN%wnewet2            )
    CALL allocate_shared_dp_1D(      SELEN%mesh%nV                    , SELEN%newtopo2,           SELEN%wnewtopo2          )
    CALL allocate_shared_dp_1D(      SELEN%mesh%nV                    , SELEN%S_global,           SELEN%wS_global          )
    CALL allocate_shared_dp_1D(      SELEN%mesh%nV                    , SELEN%U_global,           SELEN%wU_global          )
    
    newet2(             1:SELEN%mesh%nV                   ) => SELEN%newet2
    newtopo2(           1:SELEN%mesh%nV                   ) => SELEN%newtopo2
    S_global(           1:SELEN%mesh%nV                   ) => SELEN%S_global
    U_global(           1:SELEN%mesh%nV                   ) => SELEN%U_global
    
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%MSint,              SELEN%wMSint             )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%MUint,              SELEN%wMUint             )
    CALL allocate_shared_complex_1D( C%SELEN_jmax                  , SELEN%MSread,             SELEN%wMSread            )
    
    MSint(              1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%MSint
    MUint(              1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%MUint
    MSread(             1:C%SELEN_jmax                 ) => SELEN%MSread
    
    CALL allocate_shared_dp_1D(                   (C%SELEN_irreg_time_n+1), SELEN%int_time,           SELEN%wint_time          )
    CALL allocate_shared_complex_1D( C%SELEN_jmax                  , SELEN%sp1,                SELEN%wsp1               )
    CALL allocate_shared_complex_1D( C%SELEN_jmax                  , SELEN%spn,                SELEN%wspn               )
    CALL allocate_shared_complex_1D( C%SELEN_jmax                  , SELEN%preSint,            SELEN%wpreSint           )
    CALL allocate_shared_complex_1D( C%SELEN_jmax                  , SELEN%up1,                SELEN%wup1               )
    CALL allocate_shared_complex_1D( C%SELEN_jmax                  , SELEN%upn,                SELEN%wupn               )
    CALL allocate_shared_complex_1D( C%SELEN_jmax                  , SELEN%preUint,            SELEN%wpreUint           )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%s2,                 SELEN%ws2                )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%u2,                 SELEN%wu2                )
    
    int_time(                         0:C%SELEN_irreg_time_n  ) => SELEN%int_time
    sp1(                1:C%SELEN_jmax                 ) => SELEN%sp1
    spn(                1:C%SELEN_jmax                 ) => SELEN%spn
    preSint(            1:C%SELEN_jmax                 ) => SELEN%preSint
    up1(                1:C%SELEN_jmax                 ) => SELEN%up1
    upn(                1:C%SELEN_jmax                 ) => SELEN%upn
    preUint(            1:C%SELEN_jmax                 ) => SELEN%preUint
    s2(                 1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%s2
    u2(                 1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%u2
    
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%varpreoc,           SELEN%wvarpreoc          )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%varoc,              SELEN%wvaroc             )
    CALL allocate_shared_complex_1D(                (C%SELEN_irreg_time_n+1), SELEN%varoc_inv,          SELEN%wvaroc_inv         )
    
    varpreoc(           1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%varpreoc
    varoc(              1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%varoc
    varoc_inv(                          0:C%SELEN_irreg_time_n  ) => SELEN%varoc_inv
    
    ! Some external variables
    MEM_S(              1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%MEM_S
    MEM_U(              1:C%SELEN_jmax, 0:C%SELEN_irreg_time_n  ) => SELEN%MEM_U
    icethick(           1:SELEN%mesh%nV,0:C%SELEN_irreg_time_n  ) => SELEN%dice_loading_history_irreg_glob
    MJ_VAL(             1:C%SELEN_jmax                          ) => SELEN%MJ_VAL
    LJ_VAL(             1:C%SELEN_jmax                          ) => SELEN%LJ_VAL
    SELEN_anc(          1:SELEN%mesh%nV                         ) => SELEN%mesh%ianc
    INIT_TOPO(          1:SELEN%mesh%nV                         ) => SELEN%topo_ref
    ALF                                                           => SELEN%ALF
    LONG_TABLE(         0:C%SELEN_n_harmonics, 1:SELEN%mesh%nV  ) => SELEN%dLONG_TABLE
    INIT_ICE(           1:SELEN%mesh%nV                         ) => SELEN%load_ref
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE allocate_memory_and_bind_pointers
  SUBROUTINE deallocate_memory_SELEN( SELEN)
    ! Deallocate shared memory for SELEN internal data
  
    IMPLICIT NONE  
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_memory_SELEN'
    TYPE(type_SELEN_global),            INTENT(INOUT) :: SELEN
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    CALL deallocate_shared( SELEN%wDM)
    CALL deallocate_shared( SELEN%wX)
    CALL deallocate_shared( SELEN%wZE)
    CALL deallocate_shared( SELEN%wSE)
    CALL deallocate_shared( SELEN%wAAAA)
    CALL deallocate_shared( SELEN%wAAAA_MOD)
    CALL deallocate_shared( SELEN%wBBBB)
    CALL deallocate_shared( SELEN%wBBBB_MOD)
    CALL deallocate_shared( SELEN%wHHHH)
    CALL deallocate_shared( SELEN%wKKKK)
    CALL deallocate_shared( SELEN%wIIII)
    CALL deallocate_shared( SELEN%wload_ice_and_water)
    CALL deallocate_shared( SELEN%wZROT)
    CALL deallocate_shared( SELEN%wZROT_MOD)
    CALL deallocate_shared( SELEN%wSROT)
    CALL deallocate_shared( SELEN%wSROT_MOD)
    CALL deallocate_shared( SELEN%wZROTVV)
    CALL deallocate_shared( SELEN%wSROTVV)
    CALL deallocate_shared( SELEN%wS)
    CALL deallocate_shared( SELEN%wU)
    CALL deallocate_shared( SELEN%wZ)
    CALL deallocate_shared( SELEN%wnewalf)
    CALL deallocate_shared( SELEN%wslc)
    CALL deallocate_shared( SELEN%wWET)
    CALL deallocate_shared( SELEN%wnewtopo)
    CALL deallocate_shared( SELEN%wfj)
    CALL deallocate_shared( SELEN%wxcf)
    CALL deallocate_shared( SELEN%wnewet2)
    CALL deallocate_shared( SELEN%wnewtopo2)
    CALL deallocate_shared( SELEN%wS_global)
    CALL deallocate_shared( SELEN%wU_global)
    CALL deallocate_shared( SELEN%wMSint)
    CALL deallocate_shared( SELEN%wMUint)
    CALL deallocate_shared( SELEN%wMSread)
    CALL deallocate_shared( SELEN%wint_time)
    CALL deallocate_shared( SELEN%wsp1)
    CALL deallocate_shared( SELEN%wspn)
    CALL deallocate_shared( SELEN%wpreSint)
    CALL deallocate_shared( SELEN%wup1)
    CALL deallocate_shared( SELEN%wupn)
    CALL deallocate_shared( SELEN%wpreUint)
    CALL deallocate_shared( SELEN%ws2)
    CALL deallocate_shared( SELEN%wu2)
    CALL deallocate_shared( SELEN%wvarpreoc)
    CALL deallocate_shared( SELEN%wvaroc)
    CALL deallocate_shared( SELEN%wvaroc_inv)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE deallocate_memory_SELEN

END MODULE SELEN_sealevel_equation_module