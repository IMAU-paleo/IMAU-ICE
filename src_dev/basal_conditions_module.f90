MODULE basal_conditions_module

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

  ! The main routine that is called from update_general_ice_model_data
  SUBROUTINE calc_basal_conditions( grid, ice)
    ! Determine the basal conditions underneath the ice
    ! When choice_sliding_law = 'Weertman', calculate the sliding factor A_slid
    ! When choice_sliding_law = 'Coulomb' or 'Coulomb_regularised', calculate the till friction angle phi_fric and the till yield stress tauc.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
  ! ================================================
  ! ===== Exceptions for benchmark experiments =====
  ! ================================================
  
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6') THEN
        ! No sliding included here anyway
        IF (.NOT. C%no_sliding) THEN
          IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment ', C%choice_benchmark_experiment, ' should have C%no_sliding = .TRUE.!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSE
          RETURN
        END IF
      ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
        IF (.NOT. C%choice_sliding_law == 'Coulomb') THEN
          IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "SSA_icestream" should have C%choice_sliding_law = "Coulomb"!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSE
          ! Don't do anything; the sliding term beta is calculated separately here (using the analytical value for tauc)
          RETURN
        END IF
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        ! No exception; use the same parameterisation as in our realistic experiments
      ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_B') THEN
        ! No sliding included here
        IF (.NOT. C%no_sliding) THEN
          IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment ', C%choice_benchmark_experiment, ' should have C%no_sliding = .TRUE.!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSE
          RETURN
        END IF
      ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
        ! Don't do anything; the sliding term beta is calculated separately in these experiments
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_E') THEN
      
        IF (.NOT. C%choice_sliding_law == 'Coulomb') THEN
          IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "ISMIP_HOM_E" should have C%choice_sliding_law = "Coulomb"!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSE
          ! Don't do anything; the sliding term beta is calculated separately here
          RETURN
        END IF
      
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in calc_basal_conditions!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  ! =======================================================
  ! ===== End of exceptions for benchmark experiments =====
  ! =======================================================
    
    IF     (C%choice_sliding_law == 'Weertman') THEN
    
      CALL calc_Weertman_sliding_factor( grid, ice)
      
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
            
      CALL Martin_2011_till_model( grid, ice)
      
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_sliding_law "', TRIM(C%choice_sliding_law), '" not implemented in calc_basal_conditions!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_basal_conditions
  
  SUBROUTINE calc_Weertman_sliding_factor( grid, ice)
    ! Calculate the sliding factor A_slid (used when choice_sliding_law = 'Weertman')
    
    ! DENK DROM - no spatial variation included yet!

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    INTEGER                                            :: i,j

    ! Use a uniform sliding parameter
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%A_slid_a( j,i) = 3.0E-14_dp
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%A_slid_a, 'ice%A_slid_a', 'calc_Weertman_sliding_factor')
    
  END SUBROUTINE calc_Weertman_sliding_factor
  SUBROUTINE Martin_2011_till_model( grid, ice)
    ! Calculate the till friction angle phi_fric and basal yield stress tauc,
    ! using the till model by Martin et al. (2011).
    !
    ! NOTE: applied on the b-grid, following the approach from f.ETISh (Lars Zipf, personal communication, August 2021)
  
    USE derivatives_and_grids_module,    ONLY: map_a_to_b_2D

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: lambda_p
    REAL(dp)                                           :: pore_water_pressure
    REAL(dp)                                           :: w_Hb    
  
    ! Apply Matin et al. (2011) till model
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
    
  END SUBROUTINE Martin_2011_till_model
  
END MODULE basal_conditions_module
