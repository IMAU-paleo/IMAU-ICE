MODULE basal_conditions_and_sliding_module

  ! Contains all the routines for calculating the basal conditions underneath the ice.

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_ice_model
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             SSA_Schoof2006_analytical_solution

  IMPLICIT NONE
  
CONTAINS

! == Basal conditions
! ===================

  SUBROUTINE calc_basal_conditions( grid, ice)
    ! Determine the basal conditions underneath the ice

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! If no sliding is allowed, no point in determining basal conditions
    IF (C%choice_sliding_law == 'no_sliding') RETURN
    
    IF     (C%choice_basal_conditions == 'idealised') THEN
      ! Idealised cases
      
      CALL calc_basal_conditions_idealised( grid, ice)
      
    ELSEIF (C%choice_basal_conditions == 'Martin2011') THEN
      ! The Martin et al. (2011) till model
      
      CALL calc_basal_conditions_Martin2011( grid, ice)
      
    ELSE
      IF (par%master) WRITE(0,*) 'calc_basal_conditions - ERROR: unknown choice_basal_conditions "', TRIM(C%choice_basal_conditions), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_basal_conditions
  
  ! The Martin et al. (2011) till model
  SUBROUTINE calc_basal_conditions_Martin2011( grid, ice)
    ! Calculate the till friction angle phi_fric and basal yield stress tauc,
    ! using the till model by Martin et al. (2011).
    ! 
    ! Only applicable when choice_sliding_law = "Coulomb" or "Coulomb_regularised"

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: lambda_p
    REAL(dp)                                           :: pore_water_pressure
    REAL(dp)                                           :: w_Hb
    
    ! Safety
    IF (.NOT. (C%choice_sliding_law == 'Coulomb' .OR. C%choice_sliding_law == 'Coulomb_regularised')) THEN
      IF (par%master) WRITE(0,*) 'calc_basal_conditions_Martin2011 - ERROR: only applicable when choice_sliding_law = "Coulomb" or "Coulomb_regularised"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    IF (C%choice_basal_conditions /= 'Martin2011') THEN
      IF (par%master) WRITE(0,*) 'calc_basal_conditions_Martin2011 - ERROR: choice_basal_conditions should be "Martin2011"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
  
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Pore water pressure scaling factor (Martin et al., 2011, Eq. 12)
      lambda_p = MIN( 1._dp, MAX( 0._dp, 1._dp - (ice%Hb_a( j,i) - ice%SL_a( j,i) - C%Martin2011till_pwp_Hb_min) / (C%Martin2011till_pwp_Hb_max - C%Martin2011till_pwp_Hb_min) ))
  
      ! Pore water pressure (Martin et al., 2011, Eq. 11)
      pore_water_pressure = 0.96_dp * ice_density * grav * MAX(0.1_dp, ice%Hi_a( j,i)) * lambda_p 
  
      ! Till friction angle (Martin et al., 2011, Eq. 10)
      w_Hb = MIN( 1._dp, MAX( 0._dp, (ice%Hb_a( j,i) - C%Martin2011till_phi_Hb_min) / (C%Martin2011till_phi_Hb_max - C%Martin2011till_phi_Hb_min) ))
      ice%phi_fric_a( j,i) = (1._dp - w_Hb) * C%Martin2011till_phi_min + w_Hb * C%Martin2011till_phi_max
  
      ! Till yield stress (Martin et al., 2011, Eq. 9)
      ice%tauc_a( j,i) = TAN((pi / 180._dp) * ice%phi_fric_a( j,i)) * (ice_density * grav * MAX(0.1_dp, ice%Hi_a( j,i)) - pore_water_pressure)
    
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%phi_fric_a, 'ice%phi_fric_a', 'Martin_2011_till_model')
    CALL check_for_NaN_dp_2D( ice%tauc_a    , 'ice%tauc_a'    , 'Martin_2011_till_model')
    
  END SUBROUTINE calc_basal_conditions_Martin2011
  
  ! Idealised cases
  SUBROUTINE calc_basal_conditions_idealised( grid, ice)
    ! Determine the basal conditions underneath the ice
    ! 
    ! Idealised cases

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    IF     (C%choice_idealised_basal_conditions == 'SSA_icestream') THEN
      ! Idealised case: SSA_icestream (i.e. the Schoof 2006 analytical solution)
      
      CALL calc_basal_conditions_idealised_SSA_icestream( grid, ice)
      
    ELSEIF (C%choice_idealised_basal_conditions == 'MISMIP+') THEN
      ! Idealised case: MISMIP+ (see Asay-Davis et al., 2016)
      
      CALL calc_basal_conditions_idealised_MISMIPplus( grid, ice)
      
    ELSE
      IF (par%master) WRITE(0,*) 'calc_basal_conditions_idealised - ERROR: unknown choice_idealised_basal_conditions "', TRIM(C%choice_idealised_basal_conditions), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_basal_conditions_idealised
  SUBROUTINE calc_basal_conditions_idealised_SSA_icestream( grid, ice)
    ! Determine the basal conditions underneath the ice
    ! 
    ! Idealised case: SSA_icestream (i.e. the Schoof 2006 analytical solution)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: dummy1
  
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      CALL SSA_Schoof2006_analytical_solution( ABS(ice%dHs_dx_a( j,i)), ice%Hi_a( j,i), ice%A_flow_vav_a( j,i), grid%y(j), dummy1, ice%tauc_a( j,i))
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calc_basal_conditions_idealised_SSA_icestream
  SUBROUTINE calc_basal_conditions_idealised_MISMIPplus( grid, ice)
    ! Determine the basal conditions underneath the ice
    ! 
    ! Idealised case: MISMIP+ (see Asay-Davis et al., 2016)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), PARAMETER                                :: MISMIPplus_alpha_sq = 0.5_dp   ! Coulomb-law friction coefficient [unitless];         see Asay-Davis et al., 2016
    REAL(dp), PARAMETER                                :: MISMIPplus_beta_sq  = 1.0E4_dp ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3]; idem dito
  
    IF     (C%choice_sliding_law == 'Weertman') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the first (Weertman) sliding law option
       
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%beta_sq_a( j,i) = MISMIPplus_beta_sq
      END DO
      END DO
      CALL sync
          
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the second (Tsai et al., 2015) sliding law option
       
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%alpha_sq_a( j,i) = MISMIPplus_alpha_sq
        ice%beta_sq_a(  j,i) = MISMIPplus_beta_sq
      END DO
      END DO
      CALL sync
      
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the third (Schoof, 2005) sliding law option
       
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%alpha_sq_a( j,i) = MISMIPplus_alpha_sq
        ice%beta_sq_a(  j,i) = MISMIPplus_beta_sq
      END DO
      END DO
      CALL sync
      
    ELSE
      IF (par%master) WRITE(0,*) 'calc_basal_conditions_idealised_MISMIPplus - ERROR: only defined when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_basal_conditions_idealised_MISMIPplus
  
! == Sliding laws
! ===============

  SUBROUTINE calc_sliding_law( grid, ice, u_a, v_a, beta_a)
    ! Calculate the sliding term beta in the SSA/DIVA using the specified sliding law
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
      
    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No sliding allowed
      
      ice%beta_a( :,grid%i1:grid%i2) = 0._dp
      CALL sync
      
    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! Sliding laws for some idealised experiments
      
      CALL calc_sliding_law_idealised( grid, ice, u_a, v_a, beta_a)
      
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      
      CALL calc_sliding_law_Weertman( grid, ice, u_a, v_a, beta_a)
      
    ELSEIF (C%choice_sliding_law == 'Coulomb') THEN
      ! Coulomb-type sliding law
      
      CALL calc_sliding_law_Coulomb( grid, ice, u_a, v_a, beta_a)
      
    ELSEIF (C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Regularised Coulomb-type sliding law
      
      CALL calc_sliding_law_Coulomb_regularised( grid, ice, u_a, v_a, beta_a)
      
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      
      CALL calc_sliding_law_Tsai2015( grid, ice, u_a, v_a, beta_a)
      
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      
      CALL calc_sliding_law_Schoof2005( grid, ice, u_a, v_a, beta_a)
      
    ELSE
      IF (par%master) WRITE(0,*) 'calc_sliding_law - ERROR: unknown choice_sliding_law "', TRIM(C%choice_sliding_law), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_sliding_law

  SUBROUTINE calc_sliding_law_Weertman( grid, ice, u_a, v_a, beta_a)
    ! Weertman-type ("power law") sliding law
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)
      beta_a( j,i) = ice%beta_sq_a( j,i) * uabs_a ** (1._dp / C%slid_Weertman_m - 1._dp) ! Asay-Davis et al. (2016), Eq. 6
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_Weertman')
    
  END SUBROUTINE calc_sliding_law_Weertman
  SUBROUTINE calc_sliding_law_Coulomb( grid, ice, u_a, v_a, beta_a)
    ! Coulomb-type sliding law
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)
      beta_a( j,i) = ice%tauc_a( j,i) / uabs_a
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_Coulomb')
    
  END SUBROUTINE calc_sliding_law_Coulomb
  SUBROUTINE calc_sliding_law_Coulomb_regularised( grid, ice, u_a, v_a, beta_a)
    ! Regularised Coulomb-type sliding law
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)
      beta_a( j,i) = ice%tauc_a( j,i) * uabs_a ** (C%slid_Coulomb_reg_q_plastic - 1._dp) / (C%slid_Coulomb_reg_u_threshold ** C%slid_Coulomb_reg_q_plastic)
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_Coulomb_regularised')
    
  END SUBROUTINE calc_sliding_law_Coulomb_regularised
  SUBROUTINE calc_sliding_law_Tsai2015( grid, ice, u_a, v_a, beta_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Davis et al., 2016)
    ! 
    ! Asay-Davis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    ! 
    ! Tsai et al.: Marine ice-sheet profiles and stability under Coulomb basal conditions,
    ! Journal of Glaciology 61, 205–215, doi:10.3189/2015JoG14J221, 2015.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a, hf, N
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)
      hf = MAX( 0._dp, -seawater_density * ice%Hb_a( j,i) / ice_density)                                                              ! Asay-Davis et al. (2016), Eq. 10
      N  = MAX( 0._dp, ice_density * grav * (ice%Hi_a( j,i) - hf) )                                                                   ! Asay-Davis et al. (2016), Eq. 9
      beta_a( j,i) = MIN( ice%alpha_sq_a( j,i) * N, ice%beta_sq_a( j,i) * uabs_a ** (1._dp / C%slid_Weertman_m)) * uabs_a**(-1._dp)   ! Asay-Davis et al. (2016), Eq. 7
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_Tsai2015')
    
  END SUBROUTINE calc_sliding_law_Tsai2015
  SUBROUTINE calc_sliding_law_Schoof2005( grid, ice, u_a, v_a, beta_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Davis et al., 2016)
    ! 
    ! Asay-Davis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    ! 
    ! Schoof: The effect of cavitation on glacier sliding, P. Roy. Soc. A-Math. Phy., 461, 609–627, doi:10.1098/rspa.2004.1350, 2005
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a, hf, N, alpha_sq, beta_sq
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)
      hf = MAX( 0._dp, -seawater_density * ice%Hb_a( j,i) / ice_density)                                                              ! Asay-Davis et al. (2016), Eq. 10
      N  = MAX( 0._dp, ice_density * grav * (ice%Hi_a( j,i) - hf) )                                                                   ! Asay-Davis et al. (2016), Eq. 9
      alpha_sq = ice%alpha_sq_a( j,i)
      beta_sq  = ice%beta_sq_a(  j,i)
      beta_a( j,i) = ((beta_sq * uabs_a**(1._dp / C%slid_Weertman_m) * alpha_sq * N) / &                                              ! Asay-Davis et al. (2016), Eq. 11
        ((beta_sq**C%slid_Weertman_m * uabs_a + (alpha_sq * N)**C%slid_Weertman_m)**(1._dp / C%slid_Weertman_m))) * uabs_a**(-1._dp)
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_Schoof2005')
    
  END SUBROUTINE calc_sliding_law_Schoof2005
  
  SUBROUTINE calc_sliding_law_idealised( grid, ice, u_a, v_a, beta_a)
    ! Sliding laws for some idealised experiments
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    REAL(dp) :: dummy_dp
    
    ! To prevent compiler warnings...
    dummy_dp = u_a( 1,1)
    dummy_dp = v_a( 1,1)
      
    IF     (C%choice_idealised_sliding_law == 'ISMIP_HOM_C') THEN
      ! ISMIP-HOM experiment C
      
      CALL calc_sliding_law_idealised_ISMIP_HOM_C( grid, beta_a)
      
    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_D') THEN
      ! ISMIP-HOM experiment D
      
      CALL calc_sliding_law_idealised_ISMIP_HOM_D( grid, beta_a)
      
    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_E') THEN
      ! ISMIP-HOM experiment E
      
      CALL calc_sliding_law_idealised_ISMIP_HOM_E( grid, beta_a)
      
    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_F') THEN
      ! ISMIP-HOM experiment F
      
      CALL calc_sliding_law_idealised_ISMIP_HOM_F( grid, ice, beta_a)
      
    ELSE
      IF (par%master) WRITE(0,*) 'calc_sliding_law_idealised - ERROR: unknown choice_idealised_sliding_law "', TRIM(C%choice_idealised_sliding_law), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_sliding_law_idealised
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C( grid, beta_a)
    ! Sliding laws for some idealised experiments
    ! 
    ! ISMIP-HOM experiment C
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      beta_a( j,i) = 1000._dp + 1000._dp * SIN( 2._dp * pi * grid%x( i) / C%ISMIP_HOM_L) * SIN( 2._dp * pi * grid%y( j) / C%ISMIP_HOM_L)
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_idealised_ISMIP_HOM_C')
    
  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D( grid, beta_a)
    ! Sliding laws for some idealised experiments
    ! 
    ! ISMIP-HOM experiment D
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      beta_a( j,i) = 1000._dp + 1000._dp * SIN( 2._dp * pi * grid%x( i) / C%ISMIP_HOM_L)
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_idealised_ISMIP_HOM_D')
    
  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_E( grid, beta_a)
    ! Sliding laws for some idealised experiments
    ! 
    ! ISMIP-HOM experiment E: use the externally prescribed slip zone in Haut Glacier d'Arolla
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: dummy1, dummy2, dummy3
    INTEGER                                            :: ios,slides
    
    IF (par%master) THEN
      OPEN(   UNIT = 1337, FILE=C%ISMIP_HOM_E_Arolla_filename, ACTION='READ')
      DO i = 1, 51
        READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy1, dummy2, dummy3, slides
        IF (slides == 1) THEN
          DO j = 1, grid%ny
            beta_a( j,i) = 0._dp
          END DO
        ELSE
          DO j = 1, grid%ny
            beta_a( j,i) = 1.0E30_dp
          END DO
        END IF
        IF (ios /= 0) THEN
          IF (par%master) WRITE(0,*) ' calc_sliding_law_idealised_ISMIP_HOM_E - ERROR: length of text file "', TRIM(C%ISMIP_HOM_E_Arolla_filename), '" should be 51 lines!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END DO
      CLOSE( UNIT  = 1337)
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_idealised_ISMIP_HOM_E')
    
  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_E
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F( grid, ice, beta_a)
    ! Sliding laws for some idealised experiments
    ! 
    ! ISMIP-HOM experiment F
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      beta_a( j,i) = (ice%A_flow_vav_a( j,i) * 1000._dp)**(-1._dp)
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_idealised_ISMIP_HOM_F')
    
  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F
  
END MODULE basal_conditions_and_sliding_module
