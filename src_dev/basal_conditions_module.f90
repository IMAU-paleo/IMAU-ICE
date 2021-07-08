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
  USE utilities_module,                ONLY: SSA_Schoof2006_analytical_solution

  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE calc_basal_yield_stress( grid, ice)
    ! Calculate the basal yield stress tau_c

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: lambda_p               ! scaling of the pore water pressure
    REAL(dp)                                           :: pore_water_pressure    ! pore water pressure [Pa]
    
    REAL(dp), PARAMETER                                :: scaling_Hb_minimal_friction_angle = -1000._dp
    REAL(dp), PARAMETER                                :: scaling_Hb_maximal_friction_angle = 0._dp
    REAL(dp), PARAMETER                                :: minimal_friction_angle            = 5._dp
    REAL(dp), PARAMETER                                :: maximal_friction_angle            = 20._dp 
    REAL(dp)                                           :: x,Hs,Hb
    INTEGER                                            :: ios,slides
    
    ! Exceptions for benchmark experiments
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
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
        ! Use the Schoof solution
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          CALL SSA_Schoof2006_analytical_solution( ABS(ice%dHs_dx_a( j,i)), ice%Hi_a( j,i), ice%A_flow_vav_a( j,i), grid%y(j), lambda_p, ice%tauc_a( j,i))
        END DO
        END DO
        CALL sync
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        ! No exception; use the same parameterisation as in our realistic experiments
      ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_B') THEN
        ! No sliding included here
        IF (.NOT. C%no_sliding) THEN
          IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment ISMIP_HOM_A should have C%no_sliding = .TRUE.!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
        ! Beta is calculated separately in these experiments
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_E') THEN
        ! Use tau_c as a guide for the slip zone in Haut Glacier d'Arolla
      
        IF (par%master) THEN
          OPEN(   UNIT = 1337, FILE=C%ISMIP_HOM_E_Arolla_filename, ACTION='READ')
          DO i = 1, 51
            READ( UNIT = 1337, FMT=*, IOSTAT=ios) x, Hs, Hb, slides
            DO j = 1, grid%ny
              ice%tauc_a( j,i) = REAL(slides,dp)
            END DO
            IF (ios /= 0) THEN
              WRITE(0,*) ' initialise_initial_model_state_schematic_benchmarks - ERROR: length of text file "', TRIM(C%ISMIP_HOM_E_Arolla_filename), '" should be 51 lines!'
              CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
            END IF
          END DO
          CLOSE( UNIT  = 1337)
        END IF
        CALL sync
        RETURN
      
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in calc_basal_yield_stress!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    
  
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! The pore water pressure is scaled with a bedrock height dependend parameterisation
      ! Equation (12) in Martin et al. (2011)
      IF     ((ice%Hb_a( j,i) - ice%SL_a( j,i)) <= 0._dp) THEN
       lambda_p = 1._dp
      ELSEIF ((ice%Hb_a( j,i) - ice%SL_a( j,i)) >= 1000._dp) THEN
       lambda_p = 0._dp
      ELSE ! between 0 and 1000
       lambda_p = 1._dp - (ice%Hb_a( j,i) - ice%SL_a( j,i)) / 1000._dp
      END IF
  
      ! The pore water pressure, equation (11) in Martin et al. (2011)
      pore_water_pressure = 0.96_dp * ice_density * grav * MAX(0.1_dp, ice%Hi_a( j,i)) * lambda_p 
  
      ! The friction angle, used for the yield stress, equation (10) in Martin et al. (2011)
      IF     (ice%Hb_a( j,i) <= scaling_Hb_minimal_friction_angle) THEN
        ice%phi_fric_a( j,i) = minimal_friction_angle
      ELSEIF (ice%Hb_a( j,i) >= scaling_Hb_maximal_friction_angle) THEN
        ice%phi_fric_a( j,i) = maximal_friction_angle
      ELSE ! between C%scaling_Hb_maximal_friction_angle and C%scaling_Hb_minimal_friction_angle
        ice%phi_fric_a( j,i) = minimal_friction_angle + (maximal_friction_angle - minimal_friction_angle) * (1._dp &
                    + (ice%Hb_a( j,i) - scaling_Hb_maximal_friction_angle) / (scaling_Hb_maximal_friction_angle - scaling_Hb_minimal_friction_angle))
      END IF
  
      ! calculate yield stress everywhere, restrictions to sheet applied elsewhere, equation (9) in Martin et al. (2011)
      ice%tauc_a( j,i) = TAN((pi / 180._dp) * ice%phi_fric_a( j,i)) * (ice_density * grav * MAX(0.1_dp, ice%Hi_a( j,i)) - pore_water_pressure)
    
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calc_basal_yield_stress
  
END MODULE basal_conditions_module
