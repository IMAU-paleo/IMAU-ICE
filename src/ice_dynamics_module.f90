MODULE ice_dynamics_module

  ! Contains all the routines needed to calculate ice-sheet geometry at the next time
  ! step, including routines to determine said time step.
  ! NOTE: routines for calculating ice velocities             have been moved to the ice_velocity_module;
  !       routines for integrating the ice thickness equation have been moved to the ice_thickness_module.

  USE mpi
  USE configuration_module,                ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE parallel_module,                     ONLY: par, sync, cerr, ierr, &
                                                 allocate_shared_int_0D, allocate_shared_dp_0D, &
                                                 allocate_shared_int_1D, allocate_shared_dp_1D, &
                                                 allocate_shared_int_2D, allocate_shared_dp_2D, &
                                                 allocate_shared_int_3D, allocate_shared_dp_3D, &
                                                 deallocate_shared, partition_list
  USE data_types_module,                   ONLY: type_model_region, type_grid, type_ice_model, type_reference_geometry, type_restart_data
  USE utilities_module,                    ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                                 check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                                 SSA_Schoof2006_analytical_solution, vertical_average, surface_elevation, is_floating
  USE ice_velocity_module,                 ONLY: initialise_SSADIVA_solution_matrix, solve_SIA, solve_SSA, solve_DIVA, &
                                                 initialise_ice_velocity_ISMIP_HOM, initialise_velocities_from_restart_file
  USE ice_thickness_module,                ONLY: calc_dHi_dt, initialise_implicit_ice_thickness_matrix_tables, apply_ice_thickness_BC, &
                                                 remove_unconnected_shelves, initialise_target_dHi_dt
  USE general_ice_model_data_module,       ONLY: update_general_ice_model_data, determine_floating_margin_fraction, determine_masks_ice, &
                                                 determine_masks_transitions
  USE basal_conditions_and_sliding_module, ONLY: initialise_basal_conditions
  USE calving_module,                      ONLY: apply_calving_law

  USE netcdf_debug_module,                 ONLY: save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, save_variable_as_netcdf_int_3D, &
                                                 save_variable_as_netcdf_dp_1D,  save_variable_as_netcdf_dp_2D,  save_variable_as_netcdf_dp_3D
  USE netcdf_input_module,                 ONLY: read_field_from_file_2D

  IMPLICIT NONE

CONTAINS

! ===== Run ice dynamics =====
! ============================

  SUBROUTINE run_ice_model( region, t_end)
    ! Calculate ice velocities and the resulting change in ice geometry

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ice_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_timestepping == 'direct') THEN
      CALL run_ice_dynamics_direct( region, t_end)
    ELSEIF (C%choice_timestepping == 'pc') THEN
        CALL run_ice_dynamics_pc( region, t_end)
    ELSE
      CALL crash('unknown choice_timestepping "' // TRIM(C%choice_timestepping) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_model

  SUBROUTINE run_ice_dynamics_direct( region, t_end)
    ! Ice dynamics and time-stepping with the "direct" method (as was the default up to IMAU-ICE v1.1.1)
    ! NOTE: this does not work for DIVA ice dynamics!

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ice_dynamics_direct'
    INTEGER                                            :: i1,i2
    REAL(dp)                                           :: dt_crit_SIA, dt_crit_SSA, r_solver_acc, dt_max

    ! Add routine to path
    CALL init_routine( routine_name)

    i1 = region%grid%i1
    i2 = region%grid%i2

    ! Start-up phase
    ! ==============

    ! Get a more accurate velocity solution during the start-up phase to prevent initialisation "bumps"
    IF (C%dt_startup_phase > 0) THEN
      IF (region%time <= C%start_time_of_run + C%dt_startup_phase) THEN
        r_solver_acc = 0.01_dp + 0.99_dp * (region%time - C%start_time_of_run) / C%dt_startup_phase
      ELSEIF (region%time >= C%end_time_of_run   - C%dt_startup_phase) THEN
        r_solver_acc = 0.01_dp + 0.99_dp * (C%end_time_of_run - region%time) / C%dt_startup_phase
      ELSE
        r_solver_acc = 1._dp
      END IF
    ELSE
      r_solver_acc = 1._dp
    END IF

    region%ice%DIVA_SOR_nit      = C%DIVA_SOR_nit      * CEILING( 1._dp / r_solver_acc)
    region%ice%DIVA_SOR_tol      = C%DIVA_SOR_tol      * r_solver_acc
    region%ice%DIVA_SOR_omega    = C%DIVA_SOR_omega
    region%ice%DIVA_PETSc_rtol   = C%DIVA_PETSc_rtol   * r_solver_acc
    region%ice%DIVA_PETSc_abstol = C%DIVA_PETSc_abstol * r_solver_acc

    ! Reduce the time-step during the start-up phase
    IF (C%dt_startup_phase > 0) THEN
      IF     (region%time <= C%start_time_of_run + C%dt_startup_phase) THEN
        dt_max = C%dt_min + (C%dt_max - C%dt_min) * ((region%time - C%start_time_of_run) / C%dt_startup_phase)**2
      ELSEIF (region%time >= C%end_time_of_run   - C%dt_startup_phase) THEN
        dt_max = C%dt_min + (C%dt_max - C%dt_min) * ((C%end_time_of_run - region%time  ) / C%dt_startup_phase)**2
      ELSE
        dt_max = C%dt_max
      END IF
    ELSE
      dt_max = C%dt_max
    END IF

    ! Calculate ice velocities with the selected ice-dynamical approximation
    ! ======================================================================
    IF (C%do_read_velocities_from_restart) THEN
    ! Do nothing, velocities are read from the restart file
    ELSE
      IF     (C%choice_ice_dynamics == 'none') THEN
        ! Fixed ice geometry

      ELSEIF (C%choice_ice_dynamics == 'SIA') THEN
        ! Shallow ice approximation

        IF (region%time == region%t_next_SIA) THEN

          ! Calculate new ice velocities
          CALL solve_SIA( region%grid, region%ice)

          ! Calculate critical time step
          CALL calc_critical_timestep_SIA( region%grid, region%ice, dt_crit_SIA)

          IF (par%master) THEN

            ! Apply conditions to the time step
            dt_crit_SIA = MAX( C%dt_min, MIN( dt_max, dt_crit_SIA))

            ! Update timer
            region%dt_crit_SIA = dt_crit_SIA
            region%t_last_SIA  = region%time
            region%t_next_SIA  = region%time + region%dt_crit_SIA

          END IF
          CALL sync

        END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN

      ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
        ! Shallow shelf approximation

        IF (region%time == region%t_next_SSA) THEN

          ! Calculate new ice velocities
          CALL solve_SSA( region%grid, region%ice)

          ! Calculate critical time step
          CALL calc_critical_timestep_adv( region%grid, region%ice, dt_crit_SSA)

          IF (par%master) THEN

            ! Apply conditions to the time step
            dt_crit_SSA = MAX( C%dt_min, MIN( dt_max, dt_crit_SSA))

            ! Update timer
            region%dt_crit_SSA = dt_crit_SSA
            region%t_last_SSA  = region%time
            region%t_next_SSA  = region%time + region%dt_crit_SSA

          END IF
          CALL sync

        END IF ! IF (ABS(region%time - region%t_next_SSA) < dt_tol) THEN

      ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
        ! Hybrid SIA/SSA (Bueler and Brown, 2009)

        IF (region%time == region%t_next_SIA) THEN

          ! Calculate new ice velocities
          CALL solve_SIA( region%grid, region%ice)


          ! Calculate critical time step
          CALL calc_critical_timestep_SIA( region%grid, region%ice, dt_crit_SIA)

          IF (par%master) THEN

            ! Apply conditions to the time step
            dt_crit_SIA = MAX( C%dt_min, MIN( dt_max, dt_crit_SIA))

            ! Update timer
            region%dt_crit_SIA = dt_crit_SIA
            region%t_last_SIA  = region%time
            region%t_next_SIA  = region%time + region%dt_crit_SIA

          END IF
          CALL sync

        END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN

        IF (region%time == region%t_next_SSA) THEN
          ! Calculate new ice velocities

          CALL solve_SSA( region%grid, region%ice)

          ! Calculate critical time step
          CALL calc_critical_timestep_adv( region%grid, region%ice, dt_crit_SSA)

          IF (par%master) THEN

            ! Apply conditions to the time step
            dt_crit_SSA = MAX( C%dt_min, MIN( dt_max, dt_crit_SSA))

            ! Update timer
            region%dt_crit_SSA = dt_crit_SSA
            region%t_last_SSA  = region%time
            region%t_next_SSA  = region%time + region%dt_crit_SSA

          END IF
          CALL sync

        END IF ! IF (ABS(region%time - region%t_next_SSA) < dt_tol) THEN

      ELSE ! IF     (C%choice_ice_dynamics == 'SIA') THEN
        CALL crash('"direct" time stepping works only with SIA, SSA, or SIA/SSA ice dynamics, not with DIVA!')
      END IF ! IF     (C%choice_ice_dynamics == 'SIA') THEN
    END IF !(do_update_ice_velocity) THEN

    ! Adjust the time step to prevent overshooting other model components (thermodynamics, SMB, output, etc.)
    CALL determine_timesteps( region, t_end)

    !IF (par%master) WRITE(0,'(A,F7.4,A,F7.4,A,F7.4)') 'dt_crit_SIA = ', dt_crit_SIA, ', dt_crit_SSA = ', dt_crit_SSA, ', dt = ', region%dt

    ! Calculate new ice geometry
    ! ==========================

    IF (C%choice_ice_dynamics == 'none') THEN
      ! Fixed ice geometry
      region%ice%dHi_dt_a(     :,i1:i2) = 0._dp
      region%ice%Hi_tplusdt_a( :,i1:i2) = region%ice%Hi_a( :,i1:i2)
      CALL sync
    ELSE
      CALL calc_dHi_dt( region%grid, region%ice, region%SMB, region%BMB, region%dt, region%time)
      region%ice%Hi_tplusdt_a( :,i1:i2) = region%ice%Hi_a( :,i1:i2) + region%dt * region%ice%dHi_dt_a( :,i1:i2)
      CALL sync
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_dynamics_direct

  SUBROUTINE run_ice_dynamics_pc( region, t_end)
    ! Ice dynamics and time-stepping with the predictor/correct method
    ! (adopted from Yelmo, originally based on Cheng et al., 2017)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ice_dynamics_pc'
    INTEGER                                            :: i1,i2
    LOGICAL                                            :: do_update_ice_velocity
    REAL(dp)                                           :: dt_from_pc, dt_crit_adv, r_solver_acc, dt_max

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Abbreviations for cleaner code
    i1 = region%grid%i1
    i2 = region%grid%i2

    ! Determine whether or not we need to update ice velocities
    do_update_ice_velocity = .FALSE.

    IF (C%do_read_velocities_from_restart) THEN
      ! If restart, do not update ice_velocity but only the timers. Velocities are read from the restart file.
      IF     (C%choice_ice_dynamics == 'SIA') THEN
        IF (par%master) region%t_last_SIA = region%time
        ! IF (par%master) region%t_next_SIA = region%time + region%dt_crit_ice
        ! CALL sync
      ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
        IF (par%master) region%t_last_SSA = region%time
        ! IF (par%master) region%t_next_SSA = region%time + region%dt_crit_ice
        ! CALL sync
      ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
        IF (par%master) region%t_last_SIA = region%time
        IF (par%master) region%t_last_SSA = region%time
        ! IF (par%master) region%t_next_SIA = region%time + region%dt_crit_ice
        ! IF (par%master) region%t_next_SSA = region%time + region%dt_crit_ice
        ! CALL sync
      ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
        IF (par%master) region%t_last_DIVA = region%time
        ! IF (par%master) region%t_next_DIVA = region%time + region%dt_crit_ice
        ! CALL sync
      ELSE
        CALL crash('unknown choice_ice_dynamics "' // TRIM(C%choice_ice_dynamics) // '"!')
      END IF
    ELSE
      IF     (C%choice_ice_dynamics == 'none') THEN
        region%ice%dHi_dt_a(     :,i1:i2) = 0._dp
        region%ice%Hi_corr(      :,i1:i2) = region%ice%Hi_a( :,i1:i2)
        region%ice%Hi_tplusdt_a( :,i1:i2) = region%ice%Hi_a( :,i1:i2)
        CALL sync
      ELSEIF (C%choice_ice_dynamics == 'SIA') THEN
        IF (region%time == region%t_next_SIA ) do_update_ice_velocity = .TRUE.
      ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
        IF (region%time == region%t_next_SSA ) do_update_ice_velocity = .TRUE.
      ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
        IF (region%time == region%t_next_SIA ) do_update_ice_velocity = .TRUE.
      ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
        IF (region%time == region%t_next_DIVA) do_update_ice_velocity = .TRUE.
      ELSE
        CALL crash('unknown choice_ice_dynamics "' // TRIM(C%choice_ice_dynamics) // '"!')
      END IF
    END IF

    ! Start-up phase
    ! ==============

    ! Get a more accurate velocity solution during the start-up phase to prevent initialisation "bumps"
    IF (C%dt_startup_phase > 0) THEN
      IF (region%time <= C%start_time_of_run + C%dt_startup_phase) THEN
        r_solver_acc = 0.01_dp + 0.99_dp * (region%time - C%start_time_of_run) / C%dt_startup_phase
      ELSEIF (region%time >= C%end_time_of_run   - C%dt_startup_phase) THEN
        r_solver_acc = 0.01_dp + 0.99_dp * (C%end_time_of_run - region%time) / C%dt_startup_phase
      ELSE
        r_solver_acc = 1._dp
      END IF
    ELSE
      r_solver_acc = 1._dp
    END IF

    region%ice%DIVA_SOR_nit      = C%DIVA_SOR_nit      * CEILING( 1._dp / r_solver_acc)
    region%ice%DIVA_SOR_tol      = C%DIVA_SOR_tol      * r_solver_acc
    region%ice%DIVA_SOR_omega    = C%DIVA_SOR_omega
    region%ice%DIVA_PETSc_rtol   = C%DIVA_PETSc_rtol   * r_solver_acc
    region%ice%DIVA_PETSc_abstol = C%DIVA_PETSc_abstol * r_solver_acc

    ! Reduce the time-step during the start-up phase
    IF (C%dt_startup_phase > 0) THEN
      IF     (region%time <= C%start_time_of_run + C%dt_startup_phase) THEN
        dt_max = C%dt_min + (C%dt_max - C%dt_min) * ((region%time - C%start_time_of_run) / C%dt_startup_phase)**2
      ELSEIF (region%time >= C%end_time_of_run   - C%dt_startup_phase) THEN
        dt_max = C%dt_min + (C%dt_max - C%dt_min) * ((C%end_time_of_run - region%time  ) / C%dt_startup_phase)**2
      ELSE
        dt_max = C%dt_max
      END IF
    ELSE
      dt_max = C%dt_max
    END IF

    IF (do_update_ice_velocity) THEN

      ! Calculate time step based on the truncation error in ice thickness (Robinson et al., 2020, Eq. 33)
      CALL calc_critical_timestep_adv( region%grid, region%ice, dt_crit_adv)
      IF (par%master) THEN

        ! Calculate critical time step
        region%dt_crit_ice_prev = region%dt_crit_ice
        dt_from_pc              = (C%pc_epsilon / region%ice%pc_eta)**(C%pc_k_I + C%pc_k_p) * (C%pc_epsilon / region%ice%pc_eta_prev)**(-C%pc_k_p) * region%dt
        region%dt_crit_ice      = MAX( C%dt_min, MAX( 0.5_dp * region%dt_crit_ice_prev, MINVAL([ C%dt_max, 2._dp * region%dt_crit_ice_prev, dt_crit_adv, dt_from_pc])))

        ! Apply conditions to the time step
        region%dt_crit_ice = MAX( C%dt_min, MIN( dt_max, region%dt_crit_ice))

        ! Calculate zeta
        region%ice%pc_zeta      = region%dt_crit_ice / region%dt_crit_ice_prev

      END IF
      CALL sync

      ! Predictor step
      ! ==============

      ! Calculate new ice geometry
      region%ice%dHidt_Hnm1_unm1( :,i1:i2) = region%ice%dHidt_Hn_un( :,i1:i2)

      CALL calc_dHi_dt( region%grid, region%ice, region%SMB, region%BMB, region%dt_crit_ice, region%time)
      region%ice%dHidt_Hn_un( :,i1:i2) = region%ice%dHi_dt_a( :,i1:i2)

      ! Robinson et al. (2020), Eq. 30)
      region%ice%Hi_pred( :,i1:i2) = MAX(0._dp, region%ice%Hi_a( :,i1:i2) + region%dt_crit_ice * &
        ((1._dp + region%ice%pc_zeta / 2._dp) * region%ice%dHidt_Hn_un( :,i1:i2) - (region%ice%pc_zeta / 2._dp) * region%ice%dHidt_Hnm1_unm1( :,i1:i2)))
      CALL sync

      ! Update step
      ! ===========

      ! Calculate velocities for predicted geometry
      region%ice%Hi_old( :,i1:i2) = region%ice%Hi_a(    :,i1:i2)
      region%ice%Hi_a(   :,i1:i2) = region%ice%Hi_pred( :,i1:i2)
      CALL update_general_ice_model_data( region%grid, region%ice)

      IF     (C%choice_ice_dynamics == 'SIA') THEN

        ! Calculate velocities
        CALL solve_SIA(  region%grid, region%ice)

        ! Update timer
        IF (par%master) region%t_last_SIA = region%time
        IF (par%master) region%t_next_SIA = region%time + region%dt_crit_ice
        CALL sync

      ELSEIF (C%choice_ice_dynamics == 'SSA') THEN

        ! Calculate velocities
        CALL solve_SSA(  region%grid, region%ice)

        ! Update timer
        IF (par%master) region%t_last_SSA = region%time
        IF (par%master) region%t_next_SSA = region%time + region%dt_crit_ice
        CALL sync

      ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN

        ! Calculate velocities
        CALL solve_SIA(  region%grid, region%ice)
        CALL solve_SSA(  region%grid, region%ice)

        ! Update timer
        IF (par%master) region%t_last_SIA = region%time
        IF (par%master) region%t_last_SSA = region%time
        IF (par%master) region%t_next_SIA = region%time + region%dt_crit_ice
        IF (par%master) region%t_next_SSA = region%time + region%dt_crit_ice
        CALL sync

      ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN

        ! Calculate velocities
        CALL solve_DIVA( region%grid, region%ice)

        ! Update timer
        IF (par%master) region%t_last_DIVA = region%time
        IF (par%master) region%t_next_DIVA = region%time + region%dt_crit_ice
        CALL sync

      ELSE
        CALL crash('unknown choice_ice_dynamics "' // TRIM(C%choice_ice_dynamics) // '"!')
      END IF

      ! Corrector step
      ! ==============

      ! Calculate dHi_dt for the predicted ice thickness and updated velocity
      CALL calc_dHi_dt( region%grid, region%ice, region%SMB, region%BMB, region%dt_crit_ice, region%time)
      region%ice%dHidt_Hstarnp1_unp1( :,i1:i2) = region%ice%dHi_dt_a( :,i1:i2)

      ! Go back to old ice thickness. Run all the other modules (climate, SMB, BMB, thermodynamics, etc.)
      ! and only go to new (corrected) ice thickness at the end of this time loop.
      region%ice%Hi_a(     :,i1:i2) = region%ice%Hi_old(   :,i1:i2)
      CALL update_general_ice_model_data( region%grid, region%ice)

      ! Calculate "corrected" ice thickness (Robinson et al. (2020), Eq. 31)
      region%ice%Hi_corr( :,i1:i2) = MAX(0._dp, region%ice%Hi_old( :,i1:i2) + 0.5_dp * region%dt_crit_ice * &
        (region%ice%dHidt_Hn_un( :,i1:i2) + region%ice%dHidt_Hstarnp1_unp1( :,i1:i2)))

      ! Calculate applied ice thickness rate of change
      region%ice%dHi_dt_a( :,i1:i2) = (region%ice%Hi_corr( :,i1:i2) - region%ice%Hi_a( :,i1:i2)) / region%dt_crit_ice

      ! Determine truncation error
      CALL calc_pc_truncation_error( region%grid, region%ice, region%dt_crit_ice)

    END IF ! IF (do_update_ice_velocity) THEN

    ! Update region%dt
    CALL determine_timesteps( region, t_end)

    ! Calculate ice thickness at the end of the model time loop
    region%ice%Hi_tplusdt_a( :,i1:i2) = MAX( 0._dp, region%ice%Hi_a( :,i1:i2) + region%dt * region%ice%dHi_dt_a( :,i1:i2))

    !IF (par%master) WRITE(0,'(A,F7.4,A,F7.4,A,F7.4)') 'dt_crit_adv = ', dt_crit_adv, ', dt_from_pc = ', dt_from_pc, ', dt = ', region%dt

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_dynamics_pc

! ===== Ice thickness update =====
! ================================

  SUBROUTINE update_ice_thickness( grid, ice, mask_noice, refgeo_PD, refgeo_GIAeq, time)
    ! Update the ice thickness at the end of a model time loop

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_noice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_GIAeq
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_ice_thickness'
    INTEGER                                            :: i,j,ii,jj
    LOGICAL                                            :: is_shelf_or_GL, is_sheet_or_GL, is_GL

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Save the previous ice mask, for use in thermodynamics
    ice%mask_ice_a_prev( :,grid%i1:grid%i2) = ice%mask_ice_a( :,grid%i1:grid%i2)
    CALL sync

    CALL alter_ice_thickness( grid, ice, refgeo_PD, time)

    ! Set ice thickness to new value
    ice%Hi_a( :,grid%i1:grid%i2) = MAX( 0._dp, ice%Hi_tplusdt_a( :,grid%i1:grid%i2))
    CALL sync

    ! Apply calving law
    !
    ! NOTE: done twice, so that the calving front is also allowed to retreat
    IF (C%choice_calving_law == 'none') THEN
      CALL determine_masks_ice(                grid, ice)
      CALL determine_masks_transitions(        grid, ice)
      CALL determine_floating_margin_fraction( grid, ice)

      ! Remove unconnected shelves
      CALL determine_masks_ice(                grid, ice)
      CALL determine_masks_transitions(        grid, ice)
      CALL remove_unconnected_shelves(         grid, ice)
    ELSE
      CALL determine_masks_ice(                grid, ice)
      CALL determine_masks_transitions(        grid, ice)
      CALL determine_floating_margin_fraction( grid, ice)
      CALL apply_calving_law(                  grid, ice)

      CALL determine_masks_ice(                grid, ice)
      CALL determine_masks_transitions(        grid, ice)
      CALL determine_floating_margin_fraction( grid, ice)
      CALL apply_calving_law(                  grid, ice)

      ! Remove unconnected shelves
      CALL determine_masks_ice(                grid, ice)
      CALL determine_masks_transitions(        grid, ice)
      CALL remove_unconnected_shelves(         grid, ice)
    END IF

    ! Remove ice in areas where no ice is allowed (e.g. Greenland in NAM and EAS, and Ellesmere Island in GRL)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (mask_noice( j,i) == 1) THEN
        ice%Hi_a( j,i) = 0._dp
      END IF
    END DO
    END DO
    CALL sync

    ! If so specified, remove all floating ice
    IF (C%do_remove_shelves) THEN
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (is_floating( ice%Hi_a( j,i), ice%Hb_a( j,i), ice%SL_a( j,i))) THEN
        ice%Hi_a( j,i) = 0._dp
        END IF
      END DO
      END DO
      CALL sync
    END IF ! IF (C%do_remove_shelves) THEN

    ! If so specified, remove all floating ice beyond the present-day calving front
    IF (C%remove_shelves_larger_than_PD) THEN
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (refgeo_PD%Hi( j,i) == 0._dp .AND. refgeo_PD%Hb( j,i) < 0._dp) THEN
          ice%Hi_a( j,i) = 0._dp
        END IF
      END DO
      END DO
      CALL sync
    END IF ! IF (C%remove_shelves_larger_than_PD) THEN

    ! If so specified, remove all floating ice crossing the continental shelf edge
    IF (C%continental_shelf_calving) THEN
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (refgeo_GIAeq%Hi( j,i) == 0._dp .AND. refgeo_GIAeq%Hb( j,i) < C%continental_shelf_min_height) THEN
          ice%Hi_a( j,i) = 0._dp
        END IF
      END DO
      END DO
      CALL sync
    END IF ! IF (C%continental_shelf_calving) THEN


    ! Update the masks, slopes, etc.
    CALL update_general_ice_model_data( grid, ice)

    ! If wanted, relax a bit the target geometry for a while by
    ! making it equal to the current model geometry. Idea is that
    ! during this time, no SMB or BMB is applied, so the ice sheet
    ! can just flow and reshape itself. Ideal during initial years
    ! of spinup.
    IF (time >= C%relax_thick_t_start .AND. time <= C%relax_thick_t_end) THEN
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ! Make sure we are not spreading beyond observed margins
        IF (refgeo_PD%Hi( j,i) > 0._dp) THEN
          refgeo_PD%Hi( j,i) = ice%Hi_a( j,i)
          refgeo_PD%Hs( j,i) = ice%Hs_a( j,i)
          refgeo_PD%Hb( j,i) = ice%Hb_a( j,i)
        END IF
      END DO
      END DO
    END IF
    CALL sync

    ! Update deviations w.r.t. PD
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%dHi_a( j,i) = ice%Hi_a( j,i) - refgeo_PD%Hi( j,i)
      ice%dHs_a( j,i) = ice%Hs_a( j,i) - surface_elevation( refgeo_PD%Hi( j,i), refgeo_PD%Hb( j,i), 0._dp)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_ice_thickness

  SUBROUTINE alter_ice_thickness( grid, ice, refgeo_PD, time)
    ! Check if we want to keep the ice thickness fixed in time for specific areas,
    ! or delay its evolution by a given percentage each time step, both options
    ! including the possibility to reduce this constraint over time.

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),               INTENT(IN)    :: grid
    TYPE(type_ice_model),          INTENT(INOUT) :: ice
    TYPE(type_reference_geometry), INTENT(IN)    :: refgeo_PD
    REAL(dp),                      INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'alter_ice_thickness'
    INTEGER                                      :: i,j
    REAL(dp)                                     :: decay_start, decay_end, fixiness

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    CALL init_routine( routine_name)

    ! === Fixiness ===
    ! ================

    ! Intial value
    fixiness = 1._dp

    ! Make sure that the start and end times make sense
    decay_start = C%fixed_decay_t_start
    decay_end   = C%fixed_decay_t_end

    ! Compute decaying fixiness
    IF (decay_start >= decay_end) THEN
      ! This makes no sense
      fixiness = 0._dp
    ELSEIF (time <= decay_start) THEN
      ! Apply full fix/delay
      fixiness = 1._dp
    ELSEIF (time >= decay_end) THEN
      ! Remove any fix/delay
      fixiness = 0._dp
    ELSE
      ! Fixiness decreases with time
      fixiness = 1._dp - (time - decay_start) / (decay_end - decay_start)
    END IF

    ! Just in case
    fixiness = MIN( 1._dp, MAX( 0._dp, fixiness))

    IF (time >= C%relax_thick_t_start .AND. time <= C%relax_thick_t_end) THEN
      fixiness = 0._dp
    END IF

    ! === Fix, delay, limit ====
    ! ==========================

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Grounding line (grounded side)
      ! ==============================

      IF (ice%mask_gl_a( j,i) == 1) THEN

        ! Compute new ice thickness
        ice%Hi_tplusdt_a( j,i) = ice%Hi_a( j,i)         *          C%fixed_grounding_line_g * fixiness + &
                                 ice%Hi_tplusdt_a( j,i) * (1._dp - C%fixed_grounding_line_g * fixiness)

      ! Grounding line (floating side)
      ! ==============================

      ELSEIF (ice%mask_glf_a( j,i) == 1) THEN

        ! Compute new ice thickness
        ice%Hi_tplusdt_a( j,i) = ice%Hi_a( j,i)         *          C%fixed_grounding_line_f * fixiness + &
                                 ice%Hi_tplusdt_a( j,i) * (1._dp - C%fixed_grounding_line_f * fixiness)

      ! Grounded ice
      ! ============

      ELSEIF (ice%mask_sheet_a( j,i) == 1) THEN

        ! Compute new ice thickness
        ice%Hi_tplusdt_a( j,i) = ice%Hi_a( j,i)         *          C%fixed_sheet_geometry * fixiness + &
                                 ice%Hi_tplusdt_a( j,i) * (1._dp - C%fixed_sheet_geometry * fixiness)

      ! Floating ice
      ! ============

      ELSEIF (ice%mask_shelf_a( j,i) == 1) THEN

        ! Compute new ice thickness
        ice%Hi_tplusdt_a( j,i) = ice%Hi_a( j,i)         *          C%fixed_shelf_geometry * fixiness + &
                                 ice%Hi_tplusdt_a( j,i) * (1._dp - C%fixed_shelf_geometry * fixiness)

      END IF

    END DO
    END DO
    CALL sync

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE alter_ice_thickness

! ===== Time stepping =====
! =========================

  SUBROUTINE calc_critical_timestep_SIA( grid, ice, dt_crit_SIA)
    ! Calculate the critical time step for advective ice flow (CFL criterion)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp),                            INTENT(OUT)   :: dt_crit_SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_critical_timestep_SIA'
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: u_SIA_vav, v_SIA_vav, D_SIA, dt, u_a, v_a
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

    ! Add routine to path
    CALL init_routine( routine_name)

    dt_crit_SIA = C%dt_max

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Find the SIA ice diffusivity
      D_SIA = 1E-9_dp
      IF (i < grid%nx) THEN
        IF (ice%mask_sheet_a( j,i) == 1 .AND. ice%mask_sheet_a( j,i+1) == 1) THEN
          u_SIA_vav = vertical_average( ice%u_3D_SIA_cx( :,j,i))
          D_SIA = MAX( D_SIA, ABS( u_SIA_vav * ice%Hi_cx( j,i) / ice%dHs_dx_cx( j,i) ))
        END IF
      END IF
      IF (j < grid%ny) THEN
        IF (ice%mask_sheet_a( j,i) == 1 .AND. ice%mask_sheet_a( j+1,i) == 1) THEN
          v_SIA_vav = vertical_average( ice%v_3D_SIA_cy( :,j,i))
          D_SIA = MAX( D_SIA, ABS( v_SIA_vav * ice%Hi_cy( j,i) / ice%dHs_dy_cy( j,i) ))
        END IF
      END IF

      dt = grid%dx**2 / (6._dp * D_SIA)
      dt_crit_SIA = MIN(dt_crit_SIA, dt)

      ! Also check the 3D advective time step
      DO k = 1, C%nZ

        IF     (i == 1) THEN
          u_a = ABS(ice%u_3D_cx( k,j,1))
        ELSEIF (i == grid%nx) THEN
          u_a = ABS(ice%u_3D_cx( k,j,grid%nx-1))
        ELSE
          u_a = MAX( ABS(ice%u_3D_cx( k,j,i-1)), ABS(ice%u_3D_cx( k,j,i)) )
        END IF
        IF     (j == 1) THEN
          v_a = ABS(ice%v_3D_cy( k,1,i))
        ELSEIF (j == grid%ny) THEN
          v_a = ABS(ice%v_3D_cy( k,grid%ny-1,i))
        ELSE
          v_a = MAX( ABS(ice%v_3D_cy( k,j-1,i)), ABS(ice%v_3D_cy( k,j,i)) )
        END IF

        IF (u_a + v_a > 0._dp) THEN
          dt = grid%dx / (u_a + v_a)
          dt_crit_SIA = MIN(dt_crit_SIA, dt)
        END IF

      END DO

    END DO
    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_SIA, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_SIA = dt_crit_SIA * dt_correction_factor

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_critical_timestep_SIA

  SUBROUTINE calc_critical_timestep_adv( grid, ice, dt_crit_adv)
    ! Calculate the critical time step for advective ice flow (CFL criterion)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp),                            INTENT(OUT)   :: dt_crit_adv

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_critical_timestep_adv'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: dt
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

    ! Add routine to path
    CALL init_routine( routine_name)

    dt_crit_adv = 2._dp * C%dt_max

    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1

      IF (ABS(ice%u_vav_cx( j,i)) + ABS(ice%v_vav_cy( j,i)) > 0._dp ) THEN
        dt = grid%dx / (ABS(ice%u_vav_cx( j,i)) + ABS(ice%v_vav_cy( j,i)))
        dt_crit_adv = MIN( dt_crit_adv, dt)
      END IF

    END DO
    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_adv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_adv = MIN(C%dt_max, dt_crit_adv * dt_correction_factor)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_critical_timestep_adv

  SUBROUTINE calc_pc_truncation_error( grid, ice, dt)
    ! Calculate the truncation error in the ice thickness rate of change (Robinson et al., 2020, Eq. 32)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pc_truncation_error'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: eta_proc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Save previous value of eta
    ice%pc_eta_prev  = ice%pc_eta

    ! Find maximum truncation error
    eta_proc = C%pc_eta_min

    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1

      ! Calculate truncation error (Robinson et al., 2020, Eq. 32)
      ice%pc_tau( j,i) = ABS( ice%pc_zeta * (ice%Hi_corr( j,i) - ice%Hi_pred( j,i)) / ((3._dp * ice%pc_zeta + 3._dp) * dt))

      ! IF (ice%mask_sheet_a( j,i) == 1 .AND. ice%mask_gl_a( j,i) == 0) THEN
      IF (ice%mask_sheet_a( j,i) == 1 .AND. &
          ice%mask_gl_a( j+1,i-1) == 0 .AND. &
          ice%mask_gl_a( j+1,i  ) == 0 .AND. &
          ice%mask_gl_a( j+1,i+1) == 0 .AND. &
          ice%mask_gl_a( j  ,i-1) == 0 .AND. &
          ice%mask_gl_a( j  ,i  ) == 0 .AND. &
          ice%mask_gl_a( j  ,i+1) == 0 .AND. &
          ice%mask_gl_a( j-1,i-1) == 0 .AND. &
          ice%mask_gl_a( j-1,i  ) == 0 .AND. &
          ice%mask_gl_a( j-1,i+1) == 0) THEN

        eta_proc = MAX( eta_proc, ice%pc_tau( j,i))

      END IF

    END DO
    END DO
    CALL MPI_REDUCE( eta_proc, ice%pc_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pc_truncation_error

  SUBROUTINE determine_timesteps( region, t_end)
    ! Determine how long we can run just ice dynamics before another "action" (thermodynamics,
    ! GIA, output writing, inverse routine, etc.) has to be performed, and adjust the time step accordingly.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_timesteps'
    REAL(dp)                                           :: t_next

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN

      ! Determine when each model components should be updated

      t_next = region%time + C%dt_max

      ! First the ice dynamics
      ! ======================

      IF     (C%choice_ice_dynamics == 'none') THEN
        ! Just stick to the maximum time step
      ELSEIF (C%choice_ice_dynamics == 'SIA') THEN
        t_next = MIN( t_next, region%t_next_SIA)
      ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
        t_next = MIN( t_next, region%t_next_SSA)
      ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
        t_next = MIN( t_next, region%t_next_SIA)
        t_next = MIN( t_next, region%t_next_SSA)
      ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
        t_next = MIN( t_next, region%t_next_DIVA)
      ELSE
        CALL crash('unknown choice_ice_dynamics "' // TRIM(C%choice_ice_dynamics) // '"!')
      END IF ! IF (C%choice_ice_dynamics == 'SIA') THEN

      ! Then the other model components
      ! ===============================

      IF (region%time == region%t_next_thermo) THEN
        region%t_last_thermo  = region%time
        region%t_next_thermo  = region%t_last_thermo + C%dt_thermo
      END IF
      t_next = MIN( t_next, region%t_next_thermo)

      IF (region%time == region%t_next_climate) THEN
        region%t_last_climate = region%time
        region%t_next_climate = region%t_last_climate + C%dt_climate
      END IF
      t_next = MIN( t_next, region%t_next_climate)

      IF (region%time == region%t_next_ocean) THEN
        region%t_last_ocean   = region%time
        region%t_next_ocean   = region%t_last_ocean + C%dt_ocean
      END IF
      t_next = MIN( t_next, region%t_next_ocean)

      IF (region%time == region%t_next_SMB) THEN
        region%t_last_SMB     = region%time
        region%t_next_SMB     = region%t_last_SMB + C%dt_SMB
      END IF
      t_next = MIN( t_next, region%t_next_SMB)

      IF (region%time == region%t_next_BMB) THEN
        region%t_last_BMB     = region%time
        region%t_next_BMB     = region%t_last_BMB + C%dt_BMB
      END IF
      t_next = MIN( t_next, region%t_next_BMB)

      IF (region%time == region%t_next_ELRA) THEN
        region%t_last_ELRA    = region%time
        region%t_next_ELRA    = region%t_last_ELRA + C%dt_bedrock_ELRA
      END IF
      t_next = MIN( t_next, region%t_next_ELRA)

      IF (C%do_BIVgeo) THEN
        IF (region%time == region%t_next_BIV) THEN
          region%t_last_BIV     = region%time
          region%t_next_BIV     = region%t_last_BIV + C%BIVgeo_dt
        END IF
        t_next = MIN( t_next, region%t_next_BIV)
      END IF ! IF (C%do_BIVgeo) THEN

      IF (region%time == region%t_next_output) THEN
        region%t_last_output  = region%time
        region%t_next_output  = region%t_last_output + C%dt_output
      END IF
      t_next = MIN( t_next, region%t_next_output)

      IF (region%time == region%t_next_output_restart) THEN
        region%t_last_output_restart  = region%time
        region%t_next_output_restart  = region%t_last_output_restart + C%dt_output_restart
      END IF
      t_next = MIN( t_next, region%t_next_output_restart)

      IF (region%time == region%t_next_output_regional_scalar) THEN
        region%t_last_output_regional_scalar  = region%time
        region%t_next_output_regional_scalar  = region%t_last_output_regional_scalar + C%dt_output_regional_scalar
      END IF
      t_next = MIN( t_next, region%t_next_output_regional_scalar)

      ! Set time step so that we move forward to the next action
      IF (C%do_read_velocities_from_restart) THEN
        ! do nothing, dt is read from restart file
      ELSE
        region%dt = t_next - region%time
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_timesteps

  SUBROUTINE determine_actions( region)
    ! Determine how long we can run just ice dynamics before another "action" (thermodynamics,
    ! GIA, output writing, inverse routine, etc.) has to be performed, and adjust the time step accordingly.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_actions'
    REAL(dp)                                           :: t_next

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN

      ! Determine when each model components should be updated

      region%do_thermo  = .FALSE.
      IF (region%time == region%t_next_thermo) THEN
        region%do_thermo      = .TRUE.
      END IF

      region%do_climate = .FALSE.
      IF (region%time == region%t_next_climate) THEN
        region%do_climate     = .TRUE.
      END IF

      region%do_ocean   = .FALSE.
      IF (region%time == region%t_next_ocean) THEN
        region%do_ocean       = .TRUE.
      END IF

      region%do_SMB     = .FALSE.
      IF (region%time == region%t_next_SMB) THEN
        region%do_SMB         = .TRUE.
      END IF

      region%do_BMB     = .FALSE.
      IF (region%time == region%t_next_BMB) THEN
        region%do_BMB         = .TRUE.
      END IF


      region%do_ELRA    = .FALSE.
      IF (region%time == region%t_next_ELRA) THEN
        region%do_ELRA        = .TRUE.
      END IF

      region%do_BIV     = .FALSE.
      IF (C%do_BIVgeo) THEN
        IF (region%time == region%t_next_BIV) THEN
          region%do_BIV         = .TRUE.
        END IF
      END IF

      region%do_output  = .FALSE.
      IF (region%time == region%t_next_output) THEN
        region%do_output      = .TRUE.
      END IF

      region%do_output_restart  = .FALSE.
      IF (region%time == region%t_next_output_restart) THEN
        region%do_output_restart      = .TRUE.
      END IF

      region%do_output_regional_scalar  = .FALSE.
      IF (region%time == region%t_next_output_regional_scalar) THEN
        region%do_output_regional_scalar      = .TRUE.
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_actions

! ===== Administration: allocation and initialisation =====
! =========================================================

  SUBROUTINE initialise_ice_model( region, grid, ice, refgeo_init, refgeo_PD, region_name)
    ! Allocate shared memory for all the data fields of the ice dynamical module, and
    ! initialise some of them

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    !TYPE(type_restart_data),             INTENT(IN)    :: restart

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_model'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: tauc_analytical
    LOGICAL                                            :: is_ISMIP_HOM

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) '  Initialising ice dynamics model...'

    ! Allocate shared memory
    CALL allocate_ice_model( grid, ice)

    ! Initialise the ice-sheet geometry
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%Hi_a( j,i) = refgeo_init%Hi( j,i)
      ice%Hb_a( j,i) = refgeo_init%Hb( j,i)
      ice%Hs_a( j,i) = surface_elevation( ice%Hi_a( j,i), ice%Hb_a( j,i), 0._dp)
    END DO
    END DO
    CALL sync

    ! Differences w.r.t. present day
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%dHi_a( j,i) = ice%Hi_a( j,i) - refgeo_PD%Hi( j,i)
      ice%dHs_a( j,i) = ice%Hs_a( j,i) - surface_elevation( refgeo_PD%Hi( j,i), refgeo_PD%Hb( j,i), 0._dp)
    END DO
    END DO
    CALL sync

    ! if restart, intialise dHi_dt and dHb_dt from restart file.
    IF (C%do_read_velocities_from_restart) THEN
      CALL initialise_changing_geometry_from_restart_file( grid, ice, region_name)
    END IF

    ! Make sure we already start with correct boundary conditions
    CALL apply_ice_thickness_BC(        grid, ice, C%dt_min)
    CALL update_general_ice_model_data( grid, ice)
    CALL remove_unconnected_shelves(    grid, ice)
    CALL update_general_ice_model_data( grid, ice)

    ! Allocate and initialise basal conditions
    CALL initialise_basal_conditions( grid, ice, region_name)

    ! Initialise the "previous ice mask", so that the first call to thermodynamics works correctly
    ice%mask_ice_a_prev( :,grid%i1:grid%i2) = ice%mask_ice_a( :,grid%i1:grid%i2)
    CALL sync

    ! Initialise some numbers for the predictor/corrector ice thickness update method
    IF (par%master) THEN
      ice%pc_zeta        = 1._dp
      ice%pc_eta         = C%pc_epsilon
      ! ice%pc_eta_prev    = C%pc_epsilon
    END IF
    CALL sync

    IF (C%choice_refgeo_init_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'SSA_icestream') THEN
      ! Initialise ice velocity with the exact solution for the SSA icestream benchmark
      ! (since the solver converges very slowly for this one)
      DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      DO j = 1, grid%ny
        CALL SSA_Schoof2006_analytical_solution( C%SSA_icestream_tantheta, C%SSA_icestream_H, C%SSA_icestream_A, grid%y(j), ice%u_vav_cx( j,i), tauc_analytical)
        ice%u_vav_cx(  j,i) = ice%u_vav_cx( j,i) * 0.5_dp
        ice%u_base_cx( j,i) = ice%u_vav_cx( j,i)
        ice%u_SSA_cx(  j,i) = ice%u_vav_cx( j,i)
      END DO
      END DO
      CALL sync
    END IF

    ! Initialise the ISMIP-HOM experiments for faster convergence
    is_ISMIP_HOM = .FALSE.
    IF (C%choice_refgeo_init_ANT == 'idealised' .AND. &
       (C%choice_refgeo_init_idealised == 'ISMIP_HOM_A' .OR. &
        C%choice_refgeo_init_idealised == 'ISMIP_HOM_B' .OR. &
        C%choice_refgeo_init_idealised == 'ISMIP_HOM_C' .OR. &
        C%choice_refgeo_init_idealised == 'ISMIP_HOM_D')) THEN
      is_ISMIP_HOM = .TRUE.
    END IF
    IF (is_ISMIP_HOM) CALL initialise_ice_velocity_ISMIP_HOM( grid, ice)

    ! Read velocity fields from restart data
    IF (C%do_read_velocities_from_restart) THEN
      CALL initialise_velocities_from_restart_file( region, grid, ice, region_name)
    END IF

    ! Read target dHi_dt from external file
    IF (C%do_target_dhdt) THEN
      CALL initialise_target_dHi_dt( grid, ice, region_name)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_model
  SUBROUTINE initialise_changing_geometry_from_restart_file( grid, ice, region_name)
    ! Initialise dHi_dt and dHb_dt with data from a previous simulation's restart file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_ice_model),           INTENT(INOUT) :: ice
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_changing_geometry_from_restart_file'
    CHARACTER(LEN=256)                            :: filename_restart
    REAL(dp)                                      :: time_to_restart_from

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Select filename and time to restart from
    IF     (region_name == 'NAM') THEN
      filename_restart     = C%filename_refgeo_init_NAM
      time_to_restart_from = C%time_to_restart_from_NAM
    ELSEIF (region_name == 'EAS') THEN
      filename_restart     = C%filename_refgeo_init_EAS
      time_to_restart_from = C%time_to_restart_from_EAS
    ELSEIF (region_name == 'GRL') THEN
      filename_restart     = C%filename_refgeo_init_GRL
      time_to_restart_from = C%time_to_restart_from_GRL
    ELSEIF (region_name == 'ANT') THEN
      filename_restart     = C%filename_refgeo_init_ANT
      time_to_restart_from = C%time_to_restart_from_ANT
    END IF

    CALL read_field_from_file_2D(   filename_restart, 'dHi_dt', grid,  ice%dHi_dt_a,  region_name, time_to_restart_from)
    CALL read_field_from_file_2D(   filename_restart, 'dHb_dt', grid,  ice%dHb_dt_a,  region_name, time_to_restart_from)
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( ice%dHi_dt_a, 'ice%dHi_dt')
    CALL check_for_NaN_dp_2D( ice%dHb_dt_a, 'ice%dHb_dt')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_changing_geometry_from_restart_file
  SUBROUTINE allocate_ice_model( grid, ice)

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_ice_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ! ===============

    ! Basic data - ice thickness, bedrock & surface elevation, sea level (geoid elevation), englacial temperature
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_a                 , ice%wHi_a                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%Hi_cx                , ice%wHi_cx                )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%Hi_cy                , ice%wHi_cy                )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%Hi_b                 , ice%wHi_b                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_tplusdt_a         , ice%wHi_tplusdt_a         )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hb_a                 , ice%wHb_a                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%Hb_cx                , ice%wHb_cx                )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%Hb_cy                , ice%wHb_cy                )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%Hb_b                 , ice%wHb_b                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hs_a                 , ice%wHs_a                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%Hs_cx                , ice%wHs_cx                )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%Hs_cy                , ice%wHs_cy                )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%Hs_b                 , ice%wHs_b                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%SL_a                 , ice%wSL_a                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%SL_cx                , ice%wSL_cx                )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%SL_cy                , ice%wSL_cy                )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%SL_b                 , ice%wSL_b                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%TAF_a                , ice%wTAF_a                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%TAF_cx               , ice%wTAF_cx               )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%TAF_cy               , ice%wTAF_cy               )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%TAF_b                , ice%wTAF_b                )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%Ti_a                 , ice%wTi_a                 )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx-1, ice%Ti_cx                , ice%wTi_cx                )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny-1, grid%nx  , ice%Ti_cy                , ice%wTi_cy                )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny-1, grid%nx-1, ice%Ti_b                 , ice%wTi_b                 )

    ! Ice velocities
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%u_3D_a               , ice%wu_3D_a               )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%v_3D_a               , ice%wv_3D_a               )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx-1, ice%u_3D_cx              , ice%wu_3D_cx              )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny-1, grid%nx  , ice%v_3D_cy              , ice%wv_3D_cy              )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%w_3D_a               , ice%ww_3D_a               )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%u_vav_a              , ice%wu_vav_a              )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%v_vav_a              , ice%wv_vav_a              )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%u_vav_cx             , ice%wu_vav_cx             )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%v_vav_cy             , ice%wv_vav_cy             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%u_vav_cx_forrestart  , ice%wu_vav_cx_forrestart  )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%v_vav_cy_forrestart  , ice%wv_vav_cy_forrestart  )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%uabs_vav_a           , ice%wuabs_vav_a           )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%u_surf_a             , ice%wu_surf_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%v_surf_a             , ice%wv_surf_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%w_surf_a             , ice%ww_surf_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%u_surf_cx            , ice%wu_surf_cx            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%v_surf_cy            , ice%wv_surf_cy            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%uabs_surf_a          , ice%wuabs_surf_a          )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%u_base_a             , ice%wu_base_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%v_base_a             , ice%wv_base_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%w_base_a             , ice%ww_base_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%u_base_cx            , ice%wu_base_cx            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%v_base_cy            , ice%wv_base_cy            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%uabs_base_a          , ice%wuabs_base_a          )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx-1, ice%u_3D_SIA_cx          , ice%wu_3D_SIA_cx          )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny-1, grid%nx  , ice%v_3D_SIA_cy          , ice%wv_3D_SIA_cy          )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%u_SSA_cx             , ice%wu_SSA_cx             )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%v_SSA_cy             , ice%wv_SSA_cy             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%R_shear              , ice%wR_shear              )

    ! Different masks
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_land_a          , ice%wmask_land_a          )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_ocean_a         , ice%wmask_ocean_a         )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_lake_a          , ice%wmask_lake_a          )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_ice_a           , ice%wmask_ice_a           )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_sheet_a         , ice%wmask_sheet_a         )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_shelf_a         , ice%wmask_shelf_a         )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_coast_a         , ice%wmask_coast_a         )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_margin_a        , ice%wmask_margin_a        )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_gl_a            , ice%wmask_gl_a            )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_glf_a           , ice%wmask_glf_a           )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_cf_a            , ice%wmask_cf_a            )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_a               , ice%wmask_a               )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%f_grnd_a             , ice%wf_grnd_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%f_grnd_cx            , ice%wf_grnd_cx            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%f_grnd_cy            , ice%wf_grnd_cy            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%f_grnd_b             , ice%wf_grnd_b             )

    ! Ice physical properties
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%A_flow_3D_a          , ice%wA_flow_3D_a          )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx-1, ice%A_flow_3D_cx         , ice%wA_flow_3D_cx         )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny-1, grid%nx  , ice%A_flow_3D_cy         , ice%wA_flow_3D_cy         )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%A_flow_vav_a         , ice%wA_flow_vav_a         )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%A_flow_vav_cx        , ice%wA_flow_vav_cx        )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%A_flow_vav_cy        , ice%wA_flow_vav_cy        )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%A_flow_vav_b         , ice%wA_flow_vav_b         )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%Ti_pmp_a             , ice%wTi_pmp_a             )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%Cpi_a                , ice%wCpi_a                )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%Ki_a                 , ice%wKi_a                 )

    ! Spatial and temporal derivatives
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHi_dt_a             , ice%wdHi_dt_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHi_dx_a             , ice%wdHi_dx_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%dHi_dx_cx            , ice%wdHi_dx_cx            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHi_dy_a             , ice%wdHi_dy_a             )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%dHi_dy_cy            , ice%wdHi_dy_cy            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHb_a                , ice%wdHb_a                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHb_dt_a             , ice%wdHb_dt_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHs_dt_a             , ice%wdHs_dt_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHs_dx_a             , ice%wdHs_dx_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHs_dy_a             , ice%wdHs_dy_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%dHs_dx_cx            , ice%wdHs_dx_cx            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%dHs_dy_cx            , ice%wdHs_dy_cx            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%dHs_dx_cy            , ice%wdHs_dx_cy            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%dHs_dy_cy            , ice%wdHs_dy_cy            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dSL_dt_a             , ice%wdSL_dt_a             )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%dTi_dx_a             , ice%wdTi_dx_a             )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%dTi_dy_a             , ice%wdTi_dy_a             )

    ! Zeta derivatives
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%dzeta_dt_a           , ice%wdzeta_dt_a           )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%dzeta_dx_a           , ice%wdzeta_dx_a           )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%dzeta_dy_a           , ice%wdzeta_dy_a           )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dzeta_dz_a           , ice%wdzeta_dz_a           )

    ! Ice dynamics - driving stress
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%taudx_cx             , ice%wtaudx_cx             )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%taudy_cy             , ice%wtaudy_cy             )

    ! Ice dynamics - physical terms in the SSA/DIVA
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%du_dx_b              , ice%wdu_dx_b              )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%du_dy_b              , ice%wdu_dy_b              )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%dv_dx_b              , ice%wdv_dx_b              )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%dv_dy_b              , ice%wdv_dy_b              )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx-1, ice%du_dz_3D_cx          , ice%wdu_dz_3D_cx          )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny-1, grid%nx  , ice%dv_dz_3D_cy          , ice%wdv_dz_3D_cy          )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx-1, ice%du_dz_3D_cx_forrestart, ice%wdu_dz_3D_cx_forrestart)
    CALL allocate_shared_dp_3D(  C%nz, grid%ny-1, grid%nx  , ice%dv_dz_3D_cy_forrestart, ice%wdv_dz_3D_cy_forrestart)
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%visc_eff_3D_a        , ice%wvisc_eff_3D_a        )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny-1, grid%nx-1, ice%visc_eff_3D_b        , ice%wvisc_eff_3D_b        )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%visc_eff_int_a       , ice%wvisc_eff_int_a       )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%visc_eff_int_b       , ice%wvisc_eff_int_b       )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%N_a                  , ice%wN_a                  )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%N_cx                 , ice%wN_cx                 )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%N_cy                 , ice%wN_cy                 )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%N_b                  , ice%wN_b                  )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%beta_a               , ice%wbeta_a               )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%F2_a                 , ice%wF2_a                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%beta_eff_a           , ice%wbeta_eff_a           )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%beta_eff_cx          , ice%wbeta_eff_cx          )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%beta_eff_cy          , ice%wbeta_eff_cy          )
    CALL allocate_shared_dp_3D( C%nz,  grid%ny  , grid%nx  , ice%F1_3D_a              , ice%wF1_3D_a              )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%taub_cx              , ice%wtaub_cx              )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%taub_cy              , ice%wtaub_cy              )

    ! Ice dynamics - additional solver fields for the SSA/DIVA
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%u_cx_prev            , ice%wu_cx_prev            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%v_cy_prev            , ice%wv_cy_prev            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%DIVA_err_cx          , ice%wDIVA_err_cx          )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%DIVA_err_cy          , ice%wDIVA_err_cy          )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%DIVA_mask_combi      , ice%wDIVA_mask_combi      )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%DIVA_isfront_inner   , ice%wDIVA_isfront_inner   )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%DIVA_isfront_outer   , ice%wDIVA_isfront_outer   )
    CALL initialise_SSADIVA_solution_matrix( grid, ice)

    ! Ice dynamics - ice thickness calculation
    IF     (C%choice_ice_integration_method == 'none') THEN
    ELSEIF (C%choice_ice_integration_method == 'explicit') THEN
      CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%Qx_cx                , ice%wQx_cx                )
      CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%Qy_cy                , ice%wQy_cy                )
    ELSEIF (C%choice_ice_integration_method == 'semi-implicit') THEN
      CALL initialise_implicit_ice_thickness_matrix_tables( grid, ice)
    ELSE
      CALL crash('unknown choice_ice_integration_method "' // TRIM(C%choice_ice_integration_method) // '"!')
    END IF

    ! Ice dynamics - calving
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%float_margin_frac_a  , ice%wfloat_margin_frac_a  )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_eff_cf_a          , ice%wHi_eff_cf_a          )

    ! Ice dynamics - predictor/corrector ice thickness update
    CALL allocate_shared_dp_0D(                              ice%pc_zeta              , ice%wpc_zeta              )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%pc_tau               , ice%wpc_tau               )
    CALL allocate_shared_dp_0D(                              ice%pc_eta               , ice%wpc_eta               )
    CALL allocate_shared_dp_0D(                              ice%pc_eta_prev          , ice%wpc_eta_prev          )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHidt_Hnm1_unm1      , ice%wdHidt_Hnm1_unm1      )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHidt_Hn_un          , ice%wdHidt_Hn_un          )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHidt_Hn_un_forrestart, ice%wdHidt_Hn_un_forrestart)
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHidt_Hstarnp1_unp1  , ice%wdHidt_Hstarnp1_unp1  )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_old               , ice%wHi_old               )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_pred              , ice%wHi_pred              )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_corr              , ice%wHi_corr              )

    ! Thermodynamics
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_ice_a_prev      , ice%wmask_ice_a_prev      )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%frictional_heating_a , ice%wfrictional_heating_a )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%GHF_a                , ice%wGHF_a                )

    ! Isotopes
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_a_prev            , ice%wHi_a_prev            )

    ! Useful stuff
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHi_a                , ice%wdHi_a                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHs_a                , ice%wdHs_a                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHi_dt_target        , ice%wdHi_dt_target        )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_ice_model

END MODULE ice_dynamics_module
