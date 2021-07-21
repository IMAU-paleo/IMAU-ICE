MODULE ice_dynamics_module
  ! Contains all the routines needed to calculate ice-sheet geometry at the next time
  ! step, including routines to determine said time step.
  ! NOTE: routines for calculating ice velocities have been moved to the ice_velocity_module.

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_model_region, type_grid, type_ice_model, type_init_data_fields, type_SMB_model, type_BMB_model
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE utilities_module,                ONLY: SSA_Schoof2006_analytical_solution, vertical_average,  initialise_matrix_equation_CSR, &
                                             solve_matrix_equation_CSR, check_CSR_for_double_entries
  USE ice_velocity_module,             ONLY: initialise_DIVA_matrix_tables, solve_SIA, solve_SSA, solve_DIVA
  USE general_ice_model_data_module,   ONLY: update_general_ice_model_data
  USE calving_module,                  ONLY: calculate_calving_flux

  IMPLICIT NONE
  
CONTAINS

  ! The main ice dynamics routine
  SUBROUTINE run_ice_model( region, t_end)
    ! Calculate ice velocities and the resulting change in ice geometry
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end
    
    IF (C%choice_timestepping == 'direct') THEN
      CALL run_ice_dynamics_direct( region, t_end)
    ELSEIF (C%choice_timestepping == 'pc') THEN
      CALL run_ice_dynamics_pc( region, t_end)
    ELSE
      IF (par%master) WRITE(0,*) 'run_ice_model - ERROR: unknown choice_timestepping "', C%choice_timestepping, '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE run_ice_model
  SUBROUTINE run_ice_dynamics_direct( region, t_end)
    ! Ice dynamics and time-stepping with the "direct" method (as was the default up to IMAU-ICE v1.1.1)
    ! NOTE: this does not work for DIVA ice dynamics!
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end
    
    ! Local variables:
    INTEGER                                            :: i1,i2
    REAL(dp)                                           :: dt_crit_SIA, dt_crit_SSA
    
    i1 = region%grid%i1
    i2 = region%grid%i2
    
    IF     (C%choice_ice_dynamics == 'SIA') THEN
    
      IF (region%time == region%t_next_SIA) THEN
    
        ! Calculate new ice velocities
        CALL solve_SIA( region%grid, region%ice)
    
        ! Calculate critical time step
        CALL calculate_critical_timestep_SIA( region%grid, region%ice, dt_crit_SIA)
        IF (par%master) region%dt_crit_SIA = dt_crit_SIA
        
        ! Update timer
        IF (par%master) region%t_last_SIA = region%time
        IF (par%master) region%t_next_SIA = region%time + region%dt_crit_SIA
        CALL sync
        
      END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN
      
    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
    
      IF (region%time == region%t_next_SSA) THEN
    
        ! Calculate new ice velocities
        CALL solve_SSA( region%grid, region%ice)
    
        ! Calculate critical time step
        CALL calculate_critical_timestep_adv( region%grid, region%ice, dt_crit_SSA)
        IF (par%master) region%dt_crit_SSA = dt_crit_SSA
        
        ! Update timer
        IF (par%master) region%t_last_SSA = region%time
        IF (par%master) region%t_next_SSA = region%time + region%dt_crit_SSA
        CALL sync
        
      END IF ! IF (ABS(region%time - region%t_next_SSA) < dt_tol) THEN
      
    ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
    
      IF (region%time == region%t_next_SIA) THEN
    
        ! Calculate new ice velocities
        CALL solve_SIA( region%grid, region%ice)
    
        ! Calculate critical time step
        CALL calculate_critical_timestep_SIA( region%grid, region%ice, dt_crit_SIA)
        IF (par%master) region%dt_crit_SIA = dt_crit_SIA
        
        ! Update timer
        IF (par%master) region%t_last_SIA = region%time
        IF (par%master) region%t_next_SIA = region%time + region%dt_crit_SIA
        CALL sync
        
      END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN
      
      IF (region%time == region%t_next_SSA) THEN
    
        ! Calculate new ice velocities
        CALL solve_SSA( region%grid, region%ice)
    
        ! Calculate critical time step
        CALL calculate_critical_timestep_adv( region%grid, region%ice, dt_crit_SSA)
        IF (par%master) region%dt_crit_SSA = dt_crit_SSA
        
        ! Update timer
        IF (par%master) region%t_last_SSA = region%time
        IF (par%master) region%t_next_SSA = region%time + region%dt_crit_SSA
        CALL sync
        
      END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN
      
    ELSE ! IF     (C%choice_ice_dynamics == 'SIA') THEN
    
      IF (par%master) WRITE(0,*) 'run_ice_dynamics_direct - ERROR: "direct" time stepping works only with SIA, SSA, or SIA/SSA ice dynamics, not with DIVA!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      
    END IF ! IF     (C%choice_ice_dynamics == 'SIA') THEN
    
    ! Adjust the time step to prevent overshooting other model components (thermodynamics, SMB, output, etc.)
    CALL determine_timesteps_and_actions( region, t_end)
      
    IF (par%master) WRITE(0,'(A,F7.4,A,F7.4,A,F7.4)') 'dt_crit_SIA = ', dt_crit_SIA, ', dt_crit_SSA = ', dt_crit_SSA, ', dt = ', region%dt
    
    ! Calculate new ice geometry
    CALL calculate_dHi_dt( region%grid, region%ice, region%SMB, region%BMB, region%dt, region%mask_noice)
    region%ice%Hi_a_new( :,i1:i2) = region%ice%Hi_a( :,i1:i2) + region%dt * region%ice%dHi_dt_a( :,i1:i2)
    CALL sync
    
  END SUBROUTINE run_ice_dynamics_direct
  SUBROUTINE run_ice_dynamics_pc( region, t_end)
    ! Ice dynamics and time-stepping with the predictor/correct method (adopted from Yelmo,
    ! originally based on Cheng et al., 2017)
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end
    
    ! Local variables:
    INTEGER                                            :: i1,i2
    LOGICAL                                            :: do_update_ice_velocity
    REAL(dp)                                           :: dt_from_pc, dt_crit_adv
    
    ! Abbreviations for cleaner code
    i1 = region%grid%i1
    i2 = region%grid%i2
    
    ! Determine whether or not we need to update ice velocities
    do_update_ice_velocity = .FALSE.
    IF     (C%choice_ice_dynamics == 'SIA') THEN
      IF (region%time == region%t_next_SIA ) do_update_ice_velocity = .TRUE.
    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
      IF (region%time == region%t_next_SSA ) do_update_ice_velocity = .TRUE.
    ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
      IF (region%time == region%t_next_SIA ) do_update_ice_velocity = .TRUE.
    ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
      IF (region%time == region%t_next_DIVA) do_update_ice_velocity = .TRUE.
    ELSE
      IF (par%master) WRITE(0,*) 'run_ice_dynamics_pc - ERROR: unknown choice_ice_dynamics "', C%choice_ice_dynamics, '"'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    IF (do_update_ice_velocity) THEN
    
      ! Calculate time step based on the truncation error in ice thickness (Robinson et al., 2020, Eq. 33)
      CALL calculate_critical_timestep_adv( region%grid, region%ice, dt_crit_adv)
      IF (par%master) THEN
        region%dt_crit_ice_prev = region%dt_crit_ice
        region%ice%pc_eta_prev  = region%ice%pc_eta
        dt_from_pc              = (C%pc_epsilon / region%ice%pc_eta)**(C%pc_k_I + C%pc_k_p) * (C%pc_epsilon / region%ice%pc_eta_prev)**(-C%pc_k_p) * region%dt
        region%dt_crit_ice      = MAX(C%dt_min, MINVAL([ C%dt_max, 2._dp * region%dt_crit_ice_prev, dt_crit_adv, dt_from_pc]))
        region%ice%pc_zeta      = region%dt_crit_ice / region%dt_crit_ice_prev
        region%ice%pc_beta1     = 1._dp + region%ice%pc_zeta / 2._dp
        region%ice%pc_beta2     =       - region%ice%pc_zeta / 2._dp
      END IF
      CALL sync
      
      ! Predictor step
      ! ==============
      
      ! Calculate new ice geometry
      region%ice%pc_f2(   :,i1:i2) = region%ice%dHi_dt_a( :,i1:i2)
      CALL calculate_dHi_dt( region%grid, region%ice, region%SMB, region%BMB, region%dt_crit_ice, region%mask_noice)
      region%ice%pc_f1(   :,i1:i2) = region%ice%dHi_dt_a( :,i1:i2)
      region%ice%Hi_pred( :,i1:i2) = region%ice%Hi_a(     :,i1:i2) + region%dt_crit_ice * region%ice%dHi_dt_a( :,i1:i2)
      CALL sync
      
      ! Update step
      ! ===========
  
      ! Calculate velocities for predicted geometry
      region%ice%Hi_old( :,i1:i2) = region%ice%Hi_a(    :,i1:i2)
      region%ice%Hi_a(   :,i1:i2) = region%ice%Hi_pred( :,i1:i2)
      CALL update_general_ice_model_data( region%grid, region%ice, region%time)
      
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
        IF (par%master) WRITE(0,*) 'run_ice_dynamics_pc - ERROR: unknown choice_ice_dynamics "', C%choice_ice_dynamics, '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    
      ! Corrector step
      ! ==============
      
      ! Calculate "corrected" ice thickness based on new velocities
      CALL calculate_dHi_dt( region%grid, region%ice, region%SMB, region%BMB, region%dt_crit_ice, region%mask_noice)
      region%ice%pc_f3(   :,i1:i2) = region%ice%dHi_dt_a( :,i1:i2)
      region%ice%pc_f4(   :,i1:i2) = region%ice%pc_f1(    :,i1:i2)
      region%ice%Hi_corr( :,i1:i2) = region%ice%Hi_old(   :,i1:i2) + 0.5_dp * region%dt_crit_ice * (region%ice%pc_f3( :,i1:i2) + region%ice%pc_f4( :,i1:i2))
  
      ! Determine truncation error
      CALL calculate_pc_truncation_error( region%grid, region%ice, region%dt_crit_ice, region%dt_prev)
    
      ! Go back to old ice thickness. Run all the other modules (climate, SMB, BMB, thermodynamics, etc.)
      ! and only go to new (corrected) ice thickness at the end of this time loop.
      region%ice%Hi_a(     :,i1:i2) = region%ice%Hi_old(  :,i1:i2)
      region%ice%dHi_dt_a( :,i1:i2) = (region%ice%Hi_pred( :,i1:i2) - region%ice%Hi_a( :,i1:i2)) / region%dt_crit_ice
    
    END IF ! IF (do_update_ice_velocity) THEN
      
    ! Adjust the time step to prevent overshooting other model components (thermodynamics, SMB, output, etc.)
    CALL determine_timesteps_and_actions( region, t_end)
    
    ! Calculate ice thickness at the end of this model loop
    region%ice%Hi_a_new( :,i1:i2) = region%ice%Hi_a( :,i1:i2) + region%dt * region%ice%dHi_dt_a( :,i1:i2)
    
    !IF (par%master) WRITE(0,'(A,F7.4,A,F7.4,A,F7.4)') 'dt_crit_adv = ', dt_crit_adv, ', dt_from_pc = ', dt_from_pc, ', dt = ', region%dt
    
  END SUBROUTINE run_ice_dynamics_pc
  
  ! Time stepping
  SUBROUTINE calculate_critical_timestep_SIA( grid, ice, dt_crit_SIA)
    ! Calculate the critical time step for advective ice flow (CFL criterion)
    
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp),                            INTENT(OUT)   :: dt_crit_SIA
    
    ! Local variables:
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: u_SIA_vav, v_SIA_vav, D_SIA, dt, u_a, v_a
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

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
        
        dt = grid%dx / (u_a + v_a)
        dt_crit_SIA = MIN(dt_crit_SIA, dt)
      
      END DO
       
    END DO
    END DO
    
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_SIA, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_SIA = dt_crit_SIA * dt_correction_factor
        
  END SUBROUTINE calculate_critical_timestep_SIA
  SUBROUTINE calculate_critical_timestep_adv( grid, ice, dt_crit_adv)
    ! Calculate the critical time step for advective ice flow (CFL criterion)
    
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp),                            INTENT(OUT)   :: dt_crit_adv
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: dt 
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

    dt_crit_adv = 2._dp * C%dt_max
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      
      dt = grid%dx / (ABS(ice%u_vav_cx( j,i)) + ABS(ice%v_vav_cy( j,i)))
      dt_crit_adv = MIN( dt_crit_adv, dt)
       
    END DO
    END DO
    
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_adv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_adv = MIN(C%dt_max, dt_crit_adv * dt_correction_factor)
        
  END SUBROUTINE calculate_critical_timestep_adv
  SUBROUTINE calculate_pc_truncation_error( grid, ice, dt, dt_prev)
    ! Calculate the truncation error in the ice thickness rate of change (Robinson et al., 2020, Eq. 32)
    
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: dt, dt_prev
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: zeta, eta_proc
    
    ! Ratio of time steps
    zeta = dt / dt_prev
        
    ! Calculate truncation error (Robinson et al., 2020, Eq. 32)
    ice%pc_tau(:,grid%i1:grid%i2) = zeta * (ice%Hi_corr(:,grid%i1:grid%i2) - ice%Hi_pred(:,grid%i1:grid%i2)) / ((3._dp * zeta + 3._dp) * dt)
    
    ! Find maximum truncation error
    
    eta_proc = C%pc_eta_min
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%mask_sheet_a( j,i) == 1 .AND. ice%mask_gl_a( j,i) == 0) THEN
        eta_proc = MAX( eta_proc, ABS(ice%pc_tau( j,i)))
      END IF
    END DO
    END DO
    
    CALL MPI_REDUCE( eta_proc, ice%pc_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
        
  END SUBROUTINE calculate_pc_truncation_error
  SUBROUTINE determine_timesteps_and_actions( region, t_end)
    ! Determine how long we can run just ice dynamics before another "action" (thermodynamics,
    ! GIA, output writing, inverse routine, etc.) has to be performed, and adjust the time step accordingly.
    
    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end
    
    ! Local variables:
    REAL(dp)                                           :: t_next
    
    IF (par%master) THEN
      
      ! Determine when each model components should be updated
      
      t_next = MIN(t_end, region%time + C%dt_max)
      
      ! First the ice dynamics
      ! ======================
      
      IF     (C%choice_ice_dynamics == 'SIA') THEN
        t_next = MIN( t_next, region%t_next_SIA)
      ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
        t_next = MIN( t_next, region%t_next_SSA)
      ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
        t_next = MIN( t_next, region%t_next_SIA)
        t_next = MIN( t_next, region%t_next_SSA)
      ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
        t_next = MIN( t_next, region%t_next_DIVA)
      ELSE
        WRITE(0,*) 'determine_timesteps_and_actions_direct - ERROR: unknown choice_ice_dynamics "', C%choice_ice_dynamics, '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF ! IF (C%choice_ice_dynamics == 'SIA') THEN
      
      ! Then the other model components
      ! ===============================
      
      region%do_thermo  = .FALSE.
      IF (region%time == region%t_next_thermo) THEN
        region%do_thermo      = .TRUE.
        region%t_last_thermo  = region%time
        region%t_next_thermo  = region%t_last_thermo + C%dt_thermo
      END IF
      t_next = MIN( t_next, region%t_next_thermo)
      
      region%do_climate = .FALSE.
      IF (region%time == region%t_next_climate) THEN
        region%do_climate     = .TRUE.
        region%t_last_climate = region%time
        region%t_next_climate = region%t_last_climate + C%dt_climate
      END IF
      t_next = MIN( t_next, region%t_next_climate)
      
      region%do_SMB     = .FALSE.
      IF (region%time == region%t_next_SMB) THEN
        region%do_SMB         = .TRUE.
        region%t_last_SMB     = region%time
        region%t_next_SMB     = region%t_last_SMB + C%dt_SMB
      END IF
      t_next = MIN( t_next, region%t_next_SMB)
      
      region%do_BMB     = .FALSE.
      IF (region%time == region%t_next_BMB) THEN
        region%do_BMB         = .TRUE.
        region%t_last_BMB     = region%time
        region%t_next_BMB     = region%t_last_BMB + C%dt_BMB
      END IF
      t_next = MIN( t_next, region%t_next_BMB)
      
      region%do_ELRA    = .FALSE.
      IF (region%time == region%t_next_ELRA) THEN
        region%do_ELRA        = .TRUE.
        region%t_last_ELRA    = region%time
        region%t_next_ELRA    = region%t_last_ELRA + C%dt_bedrock_ELRA
      END IF
      t_next = MIN( t_next, region%t_next_ELRA)
      
      region%do_output  = .FALSE.
      IF (region%time == region%t_next_output) THEN
        region%do_output      = .TRUE.
        region%t_last_output  = region%time
        region%t_next_output  = region%t_last_output + C%dt_output
      END IF
      t_next = MIN( t_next, region%t_next_output)
      
      ! Set time step so that we move forward to the next action
      region%dt = t_next - region%time
    
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE determine_timesteps_and_actions
  SUBROUTINE calc_checkerboard_factor( grid, ice)
    ! Calculate the local "checkerboard factor" of a data field
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j,ii,jj
    LOGICAL                                            :: skip_this_one
    REAL(dp)                                           :: d1, d2
    
    ice%pc_fcb( :,grid%i1:grid%i2) = 0._dp
    
    DO i = MAX(3,grid%i1), MIN(grid%nx-2,grid%i2)
    DO j = 3, grid%ny-2
    
      skip_this_one = .FALSE.
      DO ii = i-2,i+2
      DO jj = j-2,j+2
        IF (ice%mask_ice_a( jj,ii) == 0) skip_this_one = .TRUE.
        IF (ice%mask_shelf_a( jj,ii) == 0) skip_this_one = .TRUE.
        IF (ice%mask_sheet_a( j,i) == 1 .AND. ice%mask_sheet_a( jj,ii) == 0) skip_this_one = .TRUE.
        IF (ice%mask_shelf_a( j,i) == 1 .AND. ice%mask_shelf_a( jj,ii) == 0) skip_this_one = .TRUE.
      END DO
      END DO
      IF (skip_this_one) CYCLE
      
      d1 = (ice%Hi_a( j  ,i-1) + ice%Hi_a( j  ,i+1) + ice%Hi_a( j-1,i  ) + ice%Hi_a( j+1,i  )                ) / 4._dp
      d2 = (ice%Hi_a( j  ,i-2) + ice%Hi_a( j  ,i+2) + ice%Hi_a( j-2,i  ) + ice%Hi_a( j+2,i  ) + ice%Hi_a(j,i)) / 5._dp
      
      ice%pc_fcb( j,i) = ABS(d1 - d2) / (d1 + d2)
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calc_checkerboard_factor

  ! Calculate dHi_dt
  SUBROUTINE calculate_dHi_dt( grid, ice, SMB, BMB, dt, mask_noice)
    ! Use the total ice velocities to update the ice thickness
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    REAL(dp),                            INTENT(IN)    :: dt
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_noice
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    ! Exceptions for benchmark experiments with no time evolution
    IF (C%do_benchmark_experiment) THEN
      IF     (C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'Bueler' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
              C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
        ! No exceptions here; these experiments have evolving ice geometry
      ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_E') THEN
        ! These experiments have no time evolution; don't change ice thickness
        ice%dHi_dt_a( :,grid%i1:grid%i2) = 0._dp
        CALL sync
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in calculate_dHi_dt!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Use the specified time integration method to calculate the ice thickness at t+dt
    IF     (C%choice_ice_integration_method == 'explicit') THEN
      CALL calculate_dHi_dt_explicit( grid, ice, SMB, BMB, dt)
    ELSEIF (C%choice_ice_integration_method == 'implicit') THEN
      CALL calculate_dHi_dt_implicit( grid, ice, SMB, BMB, dt)
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_ice_integration_method "', TRIM(C%choice_ice_integration_method), '" not implemented in calculate_dHi_dt!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Remove ice in areas where no ice is allowed (i.e. Greenland in NAM and EAS, and Ellesmere Island in GRL)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (mask_noice( j,i) == 1) THEN
        ice%Hi_a(     j,i) = 0._dp
        ice%dHi_dt_a( j,i) = 0._dp
      END IF
    END DO
    END DO
    CALL sync
    
    ! Apply boundary conditions
    CALL apply_ice_thickness_BC( grid, ice)
    
    ! Remove free-floating shelves not connected to any grounded ice
    CALL remove_unconnected_shelves( grid, ice)
    
    ! Calculate the calving flux
    CALL calculate_calving_flux( grid, ice, dt)
    ice%dHi_dt_a( :,grid%i1:grid%i2) = ice%dHi_dt_a( :,grid%i1:grid%i2) + ice%dHi_dt_calving_a( :,grid%i1:grid%i2)
    CALL sync
    
  END SUBROUTINE calculate_dHi_dt
  SUBROUTINE calculate_dHi_dt_explicit( grid, ice, SMB, BMB, dt)
    ! Use the total ice velocities to update the ice thickness
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    REAL(dp),                            INTENT(IN)    :: dt
    
    ! Local variables
    INTEGER                                            :: i,j
    REAL(dp)                                           :: dVi_in, dVi_out, Vi_available, rescale_factor
    REAL(dp), DIMENSION(:,:), POINTER                  :: dVi_MB
    INTEGER                                            :: wdVi_MB
        
    ! Calculate ice fluxes on the Ac grids, with ice thickness
    ! defined on the Aa grid in the upwind direction 
    ! ========================================================
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      
      ! Advective flow, so use upwind ice thickness
      IF (ice%u_vav_cx( j,i) > 0._dp) THEN
        ! Ice moves from left to right, use left-hand ice thickness
        ice%Qx_cx( j,i) = ice%u_vav_cx( j,i) * ice%Hi_a( j,i  ) * grid%dx * dt
      ELSE
        ! Ice moves from right to left, use right-hand ice thickness
        ice%Qx_cx( j,i) = ice%u_vav_cx( j,i) * ice%Hi_a( j,i+1) * grid%dx * dt
      END IF
      
      ! Flow from floating ice to open ocean is only allowed once the floating pixel is completely filled
      IF     (ice%mask_shelf_a( j  ,i  ) == 1 .AND. ice%mask_ocean_a( j  ,i+1) == 1 .AND. ice%mask_ice_a( j  ,i+1) == 0) THEN
        IF (ice%float_margin_frac_a( j  ,i  ) < 1._dp) ice%Qx_cx( j  ,i  ) = 0._dp
      ELSEIF (ice%mask_shelf_a( j  ,i+1) == 1 .AND. ice%mask_ocean_a( j  ,i  ) == 1 .AND. ice%mask_ice_a( j  ,i  ) == 0) THEN
        IF (ice%float_margin_frac_a( j  ,i+1) < 1._dp) ice%Qx_cx( j  ,i  ) = 0._dp
      END IF
      
    END DO
    END DO
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      
      ! Advective flow, so use upwind ice thickness
      IF (ice%v_vav_cy( j,i) > 0._dp) THEN
        ! Ice moves from left to right, use left-hand ice thickness
        ice%Qy_cy( j,i) = ice%v_vav_cy( j,i) * ice%Hi_a( j,i  ) * grid%dx * dt
      ELSE
        ! Ice moves from right to left, use right-hand ice thickness
        ice%Qy_cy( j,i) = ice%v_vav_cy( j,i) * ice%Hi_a( j+1,i) * grid%dx * dt
      END IF
      
      ! Flow from floating ice to open ocean is only allowed once the floating pixel is completely filled
      IF     (ice%mask_shelf_a( j  ,i  ) == 1 .AND. ice%mask_ocean_a( j+1,i  ) == 1 .AND. ice%mask_ice_a( j+1,i  ) == 0) THEN
        IF (ice%float_margin_frac_a( j  ,i  ) < 1._dp) ice%Qy_cy( j  ,i  ) = 0._dp
      ELSEIF (ice%mask_shelf_a( j+1,i  ) == 1 .AND. ice%mask_ocean_a( j  ,i  ) == 1 .AND. ice%mask_ice_a( j  ,i  ) == 0) THEN
        IF (ice%float_margin_frac_a( j+1,i  ) < 1._dp) ice%Qy_cy( j  ,i  ) = 0._dp
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Correct outfluxes for possible resulting negative ice thicknesses
    ! =================================================================
    
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dVi_MB, wdVi_MB)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Ice volume added to each grid cell through the (surface + basal) mass balance
      ! => With an exception for the calving front, where we only apply
      !    the mass balance to the floating fraction
      IF (ice%mask_cf_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 1) THEN
        dVi_MB( j,i) = (SMB%SMB_year( j,i) + BMB%BMB( j,i)) * grid%dx * grid%dx * dt * ice%float_margin_frac_a( j,i)
      ELSE
        dVi_MB( j,i) = (SMB%SMB_year( j,i) + BMB%BMB( j,i)) * grid%dx * grid%dx * dt
      END IF
    
      ! Check how much ice is available for melting or removing (in m^3)
      Vi_available = ice%Hi_a( j,i) * grid%dx * grid%dx
      
      ! Check how much ice is moving into this grid cell from all four (existing) neighbours
      dVi_in   = 0._dp
      dVi_out  = 0._dp
      IF (i > 1) THEN
        IF (ice%Qx_cx( j  ,i-1) > 0._dp) THEN
          dVi_in  = dVi_in  + ice%Qx_cx( j  ,i-1)
        ELSE
          dVi_out = dVi_out - ice%Qx_cx( j  ,i-1)
        END IF
      END IF
      IF (i < grid%nx) THEN
        IF (ice%Qx_cx( j  ,i  ) < 0._dp) THEN
          dVi_in  = dVi_in  - ice%Qx_cx( j  ,i  )
        ELSE
          dVi_out = dVi_out + ice%Qx_cx( j  ,i  )
        END IF
      END IF
      IF (j > 1) THEN
        IF (ice%Qy_cy( j-1,i  ) > 0._dp) THEN
          dVi_in  = dVi_in  + ice%Qy_cy( j-1,i  )
        ELSE
          dVi_out = dVi_out - ice%Qy_cy( j-1,i  )
        END IF
      END IF
      IF (j < grid%ny) THEN
        IF (ice%Qy_cy( j  ,i  ) < 0._dp) THEN
          dVi_in  = dVi_in  - ice%Qy_cy( j  ,i  )
        ELSE
          dVi_out = dVi_out + ice%Qy_cy( j  ,i  )
        END IF
      END IF
      
      rescale_factor = 1._dp      
      
      ! If all the ice already present melts away, there can be no outflux.
      IF (-dVi_MB( j,i) >= Vi_available) THEN      
        ! All ice in this vertex melts, nothing remains to move around. Rescale outfluxes to zero.
        dVi_MB( j,i) = -Vi_available
        rescale_factor = 0._dp
      END IF
            
      ! If the total outflux exceeds the available ice plus SMB plus total influx, rescale outfluxes
      IF (dVi_out > Vi_available + dVi_MB( j,i)) THEN      
        ! Total outflux plus melt exceeds available ice volume. Rescale outfluxes to correct for this.
        rescale_factor = (Vi_available + dVi_MB( j,i)) / dVi_out        
      END IF
      
      ! Rescale ice outfluxes out of this grid cell
      IF (rescale_factor < 1._dp) THEN
        IF (i > 1) THEN
          IF (ice%Qx_cx( j  ,i-1) < 0._dp) ice%Qx_cx( j  ,i-1) = ice%Qx_cx( j  ,i-1) * rescale_factor
        END IF
        IF (i < grid%nx) THEN
          IF (ice%Qx_cx( j  ,i  ) > 0._dp) ice%Qx_cx( j  ,i  ) = ice%Qx_cx( j  ,i  ) * rescale_factor
        END IF
        IF (j > 1) THEN
          IF (ice%Qy_cy( j-1,i  ) < 0._dp) ice%Qy_cy( j-1,i  ) = ice%Qy_cy( j-1,i  ) * rescale_factor
        END IF
        IF (j < grid%ny) THEN
          IF (ice%Qy_cy( j  ,i  ) > 0._dp) ice%Qy_cy( j  ,i  ) = ice%Qy_cy( j  ,i  ) * rescale_factor
        END IF
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Calculate change in ice thickness over time at every vertex
    ! ===========================================================
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%dHi_dt_a( j,i) = dVi_MB( j,i) / (grid%dx * grid%dx * dt) ! m/y
      IF (i > 1      ) ice%dHi_dt_a( j,i) = ice%dHi_dt_a( j,i) + ice%Qx_cx( j  ,i-1) / (grid%dx * grid%dx * dt)
      IF (i < grid%nx) ice%dHi_dt_a( j,i) = ice%dHi_dt_a( j,i) - ice%Qx_cx( j  ,i  ) / (grid%dx * grid%dx * dt)
      IF (j > 1      ) ice%dHi_dt_a( j,i) = ice%dHi_dt_a( j,i) + ice%Qy_cy( j-1,i  ) / (grid%dx * grid%dx * dt)
      IF (j < grid%ny) ice%dHi_dt_a( j,i) = ice%dHi_dt_a( j,i) - ice%Qy_cy( j  ,i  ) / (grid%dx * grid%dx * dt)
    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wdVi_MB)
    
  END SUBROUTINE calculate_dHi_dt_explicit
  SUBROUTINE calculate_dHi_dt_implicit( grid, ice, SMB, BMB, dt)
    ! Calculate the ice thickness at time t+dt using an implicit matrix-based method
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    REAL(dp),                            INTENT(IN)    :: dt
    
    ! Local variables
    INTEGER                                            :: i,j,n,k
    REAL(dp)                                           :: Hi_new
    
    ! Fill the sparse matrix
    ! ======================
    
    IF (par%master) THEN
      
      n = 0
      k = 0
      ice%dHi_m%A_ptr(1) = 1
      
      DO i = 1, grid%nx
      DO j = 1, grid%ny
      
        n = n+1
        
        IF (i == 1 .AND. j == 1) THEN
          ! Exception for the southwestern corner
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j,i)
          ice%dHi_m%A_val(   k) = 1._dp / dt + (MAX( ice%u_vav_cx( j,i),0._dp)                                   ) / grid%dx + &
                                               (MAX( ice%v_vav_cy( j,i),0._dp)                                   ) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i+1)
          ice%dHi_m%A_val(   k) =  MIN( ice%u_vav_cx( j  ,i  ),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j+1,i  )
          ice%dHi_m%A_val(   k) =  MIN( ice%v_vav_cy( j  ,i  ),0._dp) / grid%dx
          
        ELSEIF (i == grid%nx .AND. j == 1) THEN
          ! Exception for the southeastern corner
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j,i)
          ice%dHi_m%A_val(   k) = 1._dp / dt + (                                -MIN( ice%u_vav_cx( j,i-1),0._dp)) / grid%dx + &
                                               (MAX( ice%v_vav_cy( j,i),0._dp)                                   ) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i-1)
          ice%dHi_m%A_val(   k) = -MAX( ice%u_vav_cx( j  ,i-1),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j+1,i  )
          ice%dHi_m%A_val(   k) =  MIN( ice%v_vav_cy( j  ,i  ),0._dp) / grid%dx
          
        ELSEIF (i == 1 .AND. j == grid%ny) THEN
          ! Exception for the northwestern corner
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j,i)
          ice%dHi_m%A_val(   k) = 1._dp / dt + (MAX( ice%u_vav_cx( j,i),0._dp)                                   ) / grid%dx + &
                                               (                                -MIN( ice%v_vav_cy( j-1,i),0._dp)) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i+1)
          ice%dHi_m%A_val(   k) =  MIN( ice%u_vav_cx( j  ,i  ),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j-1,i  )
          ice%dHi_m%A_val(   k) = -MAX( ice%v_vav_cy( j-1,i  ),0._dp) / grid%dx
          
        ELSEIF (i == grid%nx .AND. j == grid%ny) THEN
          ! Exception for the northeastern corner
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j,i)
          ice%dHi_m%A_val(   k) = 1._dp / dt + (                                -MIN( ice%u_vav_cx( j,i-1),0._dp)) / grid%dx + &
                                               (                                -MIN( ice%v_vav_cy( j-1,i),0._dp)) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i-1)
          ice%dHi_m%A_val(   k) = -MAX( ice%u_vav_cx( j  ,i-1),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j-1,i  )
          ice%dHi_m%A_val(   k) = -MAX( ice%v_vav_cy( j-1,i  ),0._dp) / grid%dx
          
        ELSEIF (i == 1) THEN
          ! Exception for the western boundary
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j,i)
          ice%dHi_m%A_val(   k) = 1._dp / dt + (MAX( ice%u_vav_cx( j,i),0._dp)                                   ) / grid%dx + &
                                               (MAX( ice%v_vav_cy( j,i),0._dp) - MIN( ice%v_vav_cy( j-1,i),0._dp)) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i+1)
          ice%dHi_m%A_val(   k) =  MIN( ice%u_vav_cx( j  ,i  ),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j+1,i  )
          ice%dHi_m%A_val(   k) =  MIN( ice%v_vav_cy( j  ,i  ),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j-1,i  )
          ice%dHi_m%A_val(   k) = -MAX( ice%v_vav_cy( j-1,i  ),0._dp) / grid%dx
          
        ELSEIF (i == grid%nx) THEN
          ! Exception for the eastern boundary
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j,i)
          ice%dHi_m%A_val(   k) = 1._dp / dt + (                                -MIN( ice%u_vav_cx( j,i-1),0._dp)) / grid%dx + &
                                               (MAX( ice%v_vav_cy( j,i),0._dp) - MIN( ice%v_vav_cy( j-1,i),0._dp)) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i-1)
          ice%dHi_m%A_val(   k) = -MAX( ice%u_vav_cx( j  ,i-1),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j+1,i  )
          ice%dHi_m%A_val(   k) =  MIN( ice%v_vav_cy( j  ,i  ),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j-1,i  )
          ice%dHi_m%A_val(   k) = -MAX( ice%v_vav_cy( j-1,i  ),0._dp) / grid%dx
          
        ELSEIF (j == 1) THEN
          ! Exception for the southern boundary
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j,i)
          ice%dHi_m%A_val(   k) = 1._dp / dt + (MAX( ice%u_vav_cx( j,i),0._dp) - MIN( ice%u_vav_cx( j,i-1),0._dp)) / grid%dx + &
                                               (MAX( ice%v_vav_cy( j,i),0._dp)                                   ) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i+1)
          ice%dHi_m%A_val(   k) =  MIN( ice%u_vav_cx( j  ,i  ),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i-1)
          ice%dHi_m%A_val(   k) = -MAX( ice%u_vav_cx( j  ,i-1),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j+1,i  )
          ice%dHi_m%A_val(   k) =  MIN( ice%v_vav_cy( j  ,i  ),0._dp) / grid%dx
          
        ELSEIF (j == grid%ny) THEN
          ! Exception for the northern boundary
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j,i)
          ice%dHi_m%A_val(   k) = 1._dp / dt + (MAX( ice%u_vav_cx( j,i),0._dp) - MIN( ice%u_vav_cx( j,i-1),0._dp)) / grid%dx + &
                                               (                                -MIN( ice%v_vav_cy( j-1,i),0._dp)) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i+1)
          ice%dHi_m%A_val(   k) =  MIN( ice%u_vav_cx( j  ,i  ),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i-1)
          ice%dHi_m%A_val(   k) = -MAX( ice%u_vav_cx( j  ,i-1),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j-1,i  )
          ice%dHi_m%A_val(   k) = -MAX( ice%v_vav_cy( j-1,i  ),0._dp) / grid%dx
          
        ELSE ! IF (i == 1 .AND. j == 1) THEN
          ! General equation for mass continuity
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j,i)
          ice%dHi_m%A_val(   k) = 1._dp / dt + (MAX( ice%u_vav_cx( j,i),0._dp) - MIN( ice%u_vav_cx( j,i-1),0._dp)) / grid%dx + &
                                               (MAX( ice%v_vav_cy( j,i),0._dp) - MIN( ice%v_vav_cy( j-1,i),0._dp)) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i+1)
          ice%dHi_m%A_val(   k) =  MIN( ice%u_vav_cx( j  ,i  ),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i-1)
          ice%dHi_m%A_val(   k) = -MAX( ice%u_vav_cx( j  ,i-1),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j+1,i  )
          ice%dHi_m%A_val(   k) =  MIN( ice%v_vav_cy( j  ,i  ),0._dp) / grid%dx
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j-1,i  )
          ice%dHi_m%A_val(   k) = -MAX( ice%v_vav_cy( j-1,i  ),0._dp) / grid%dx
          
        END IF ! IF (i == 1 .AND. j == 1) THEN
          
        ! Right-hand side and initial guess
        ice%dHi_m%b( n) = MAX(0._dp, SMB%SMB_year( j,i) + BMB%BMB( j,i) + ice%Hi_a( j,i) / dt)
        ice%dHi_m%x( n) = ice%Hi_a( j,i)
        
        ! Next equation
        ice%dHi_m%A_ptr( n+1) = k+1
        
      END DO
      END DO
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Solve the matrix equation
    ! =========================
    
    CALL solve_matrix_equation_CSR( ice%dHi_m, C%dHi_choice_matrix_solver, &
      SOR_nit = C%dHi_SOR_nit, SOR_tol = C%dHi_SOR_tol, SOR_omega = C%dHi_SOR_omega, &
      PETSc_rtol = C%dHi_PETSc_rtol, PETSc_abstol = C%dHi_PETSc_abstol)
    
    ! Map the results back to the model grid
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      n = ice%dHi_ij2n( j,i)
      Hi_new = ice%dHi_m%x( n)
      ice%dHi_dt_a( j,i) = (Hi_new - ice%Hi_a( j,i)) / dt
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calculate_dHi_dt_implicit
  SUBROUTINE initialise_implicit_ice_thickness_matrix_tables( grid, ice)
    ! Allocate and initialise the translation tables and sparse matrix lists used 
    ! in the implicit ice thickness update method.
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j,n,neq,nnz_per_row_max
    
    ! Initialise the sparse matrix
    neq             = grid%nx * grid%ny
    nnz_per_row_max = 5
    
    CALL initialise_matrix_equation_CSR( ice%dHi_m, neq, neq, nnz_per_row_max)
    
    ! Initialise and fill the matrix-vector table
    CALL allocate_shared_int_2D( grid%ny, grid%nx, ice%dHi_ij2n, ice%wdHi_ij2n)
    IF (par%master) THEN
      n = 0
      DO i = 1, grid%nx
      DO j = 1, grid%ny
        n = n+1
        ice%dHi_ij2n( j,i) = n
      END DO
      END DO
    END IF
    CALL sync
    
  END SUBROUTINE initialise_implicit_ice_thickness_matrix_tables
  SUBROUTINE apply_ice_thickness_BC( grid, ice)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j

    IF (.NOT. C%do_benchmark_experiment) THEN
          
      ! For realistic experiments set ice thickness to zero at the domain boundary
      ice%Hi_a(     1              ,grid%i1:grid%i2) = 0._dp
      ice%Hi_a(     grid%ny        ,grid%i1:grid%i2) = 0._dp
      ice%Hi_a(     grid%j1:grid%j2,1              ) = 0._dp
      ice%Hi_a(     grid%j1:grid%j2,grid%nx        ) = 0._dp
      ice%dHi_dt_a( 1              ,grid%i1:grid%i2) = 0._dp
      ice%dHi_dt_a( grid%ny        ,grid%i1:grid%i2) = 0._dp
      ice%dHi_dt_a( grid%j1:grid%j2,1              ) = 0._dp
      ice%dHi_dt_a( grid%j1:grid%j2,grid%nx        ) = 0._dp
      CALL sync
      
    ELSE
        
      IF     (C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_E') THEN
        ! No exceptions here; in these cases, ice thickness is not updated anyway
        ! (so we should never reach this point!)
        
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" should not change ice thickness!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        
      ELSEIF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
              C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'Bueler') THEN
          
        ! Apply boundary conditions: set ice thickness to zero at the domain boundary
        ice%Hi_a(     1              ,grid%i1:grid%i2) = 0._dp
        ice%Hi_a(     grid%ny        ,grid%i1:grid%i2) = 0._dp
        ice%Hi_a(     grid%j1:grid%j2,1              ) = 0._dp
        ice%Hi_a(     grid%j1:grid%j2,grid%nx        ) = 0._dp
        ice%dHi_dt_a( 1              ,grid%i1:grid%i2) = 0._dp
        ice%dHi_dt_a( grid%ny        ,grid%i1:grid%i2) = 0._dp
        ice%dHi_dt_a( grid%j1:grid%j2,1              ) = 0._dp
        ice%dHi_dt_a( grid%j1:grid%j2,grid%nx        ) = 0._dp
        CALL sync
        
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        
        ! Create a nice circular ice shelf
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          IF (SQRT(grid%x(i)**2+grid%y(j)**2) > grid%xmax * 0.95_dp) THEN
            ice%Hi_a(     j,i) = 0._dp
            ice%dHi_dt_a( j,i) = 0._dp
          END IF
        END DO
        END DO
        CALL sync
        
      ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
      
        ! Fixed ice thickness at the boundary
        ice%Hi_a(     1              ,grid%i1:grid%i2) = 1000._dp
        ice%Hi_a(     grid%ny        ,grid%i1:grid%i2) = 1000._dp
        ice%Hi_a(     grid%j1:grid%j2,1              ) = 1000._dp
        ice%Hi_a(     grid%j1:grid%j2,grid%nx        ) = 1000._dp
        ice%dHi_dt_a( 1              ,grid%i1:grid%i2) = 0._dp
        ice%dHi_dt_a( grid%ny        ,grid%i1:grid%i2) = 0._dp
        ice%dHi_dt_a( grid%j1:grid%j2,1              ) = 0._dp
        ice%dHi_dt_a( grid%j1:grid%j2,grid%nx        ) = 0._dp
        CALL sync
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in apply_ice_thickness_BC!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END IF ! IF (.NOT. C%do_benchmark_experiment) THEN
    
  END SUBROUTINE apply_ice_thickness_BC
  SUBROUTINE remove_unconnected_shelves( grid, ice)
    ! Use a flood-fill algorithm to find all shelves connected to sheets.
    ! Remove all other shelves.
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables
    INTEGER                                            :: i,j
    INTEGER, DIMENSION(:,:  ), ALLOCATABLE             :: map
    INTEGER, DIMENSION(:,:  ), ALLOCATABLE             :: stack
    INTEGER                                            :: stackN
    
    IF (par%master) THEN
      
      ! Allocate memory for map and stack
      ALLOCATE( map(   grid%ny,  grid%nx   ))
      ALLOCATE( stack( grid%ny * grid%nx, 2))
      
      map    = 0
      stack  = 0
      stackN = 0
      
      ! Fill the stack with all shelf-next-to-sheet points
      DO i = 2, grid%nx-1
      DO j = 2, grid%ny-1
        IF (ice%mask_shelf_a( j,i) == 1) THEN
          IF (ice%mask_sheet_a( j  ,i+1) == 1 .OR. &
              ice%mask_sheet_a( j  ,i-1) == 1 .OR. &
              ice%mask_sheet_a( j+1,i  ) == 1 .OR. &
              ice%mask_sheet_a( j-1,i  ) == 1) THEN
            stackN = stackN + 1
            stack( stackN,:) = [i,j]
            map( j,i) = 1
          END IF
        END IF
      END DO
      END DO
      
      ! Mark all connected shelf points on the map using a flood fill
      DO WHILE (stackN > 0)
      
        ! Remove the last element from the stack
        i = stack( stackN,1)
        j = stack( stackN,2)
        stackN = stackN - 1
        
        ! Check if this element is shelf. If so, mark it on the map,
        ! and add its neighbours to the stack (if they're not in it yet).
        ! If not, do nothing
        IF (ice%mask_shelf_a( j,i) == 1) THEN
          ! This element is shelf. Mark it on the map
          map( j,i) = 2
          ! Add its four non-stacked neighbours to the stack
          IF (i>1) THEN
            IF (map(j,i-1)==0) THEN
              map( j,i-1) = 1
              stackN = stackN + 1
              stack( stackN,:) = [i-1,j]
            END IF
          END IF
          IF (i<grid%nx) THEN
            IF (map(j,i+1)==0) THEN
              map( j,i+1) = 1
              stackN = stackN + 1
              stack( stackN,:) = [i+1,j]
            END IF
          END IF
          IF (j>1) THEN
            IF (map(j-1,i)==0) THEN
              map( j-1,i) = 1
              stackN = stackN + 1
              stack( stackN,:) = [i,j-1]
            END IF
          END IF
          IF (j<grid%ny) THEN
            IF (map(j+1,i)==0) THEN
              map( j+1,i) = 1
              stackN = stackN + 1
              stack( stackN,:) = [i,j+1]
            END IF
          END IF
        END IF
      
      END DO ! DO WHILE (stackN > 0)
      
      ! Remove ice for all unconnected shelves
      DO i = 1, grid%nx
      DO j = 1, grid%ny
        IF (ice%mask_shelf_a( j,i) == 1 .AND. map( j,i) == 0) THEN
          ice%Hi_a( j,i) = 0._dp
        END IF
      END DO
      END DO
      
      ! Clean up after yourself
      DEALLOCATE( map)
      DEALLOCATE( stack)
      
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE remove_unconnected_shelves
  
  ! Administration: allocation and initialisation
  SUBROUTINE initialise_ice_model( grid, ice, init)
    ! Allocate shared memory for all the data fields of the ice dynamical module, and
    ! initialise some of them
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_init_data_fields),         INTENT(IN)    :: init
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: tauc_analytical
    
    IF (par%master) WRITE(0,*) '  Initialising ice dynamics model...'
    
    ! Allocate shared memory
    CALL allocate_ice_model( grid, ice)
        
    ! Initialise with data from initial file or restart file
    IF (.NOT. C%is_restart) THEN
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%Hi_a( j,i) = init%Hi( j,i)
        ice%Hb_a( j,i) = init%Hb( j,i)
        ice%Hs_a( j,i) = MAX( ice%SL_a( j,i), ice%Hb_a( j,i) + ice%Hi_a( j,i))
      END DO
      END DO
      CALL sync
    
    ELSE
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%Hi_a( j,i) = init%Hi( j,i)
        ice%Hb_a( j,i) = init%Hb( j,i)
        ice%Hs_a( j,i) = MAX( ice%SL_a( j,i), ice%Hb_a( j,i) + ice%Hi_a( j,i))
      END DO
      END DO
      CALL sync
      
    END IF ! IF (.NOT. C%is_restart) THEN
    
    ! Make sure we already start with correct boundary conditions
    CALL apply_ice_thickness_BC( grid, ice)
    
    ! Allocate and initialise grid-cell-to-matrix-row translation tables used in the DIVA solver
    CALL initialise_DIVA_matrix_tables( grid, ice)
    
    ! Initialise some numbers for the predictor/corrector ice thickness update method
    IF (par%master) THEN
      ice%pc_zeta        = 1._dp
      ice%pc_eta         = C%pc_epsilon
      ice%pc_eta_prev    = C%pc_epsilon
      ice%pc_beta1       =  1.5_dp
      ice%pc_beta2       = -0.5_dp
      ice%pc_beta3       =  0.5_dp
      ice%pc_beta4       =  0.5_dp
    END IF
    CALL sync
    
    IF (C%do_benchmark_experiment .AND. C%choice_benchmark_experiment == 'SSA_icestream') THEN
      ! Initialise ice velocity with the exact solution for the SSA icestream benchmark
      ! (since the solver converges very slowly for this one)
      DO i = grid%i1, MIN(grid%nx-1,grid%i2)
      DO j = 1, grid%ny
        CALL SSA_Schoof2006_analytical_solution( 0.001_dp, 2000._dp, (3.7E8_dp ** (-C%n_flow)) * sec_per_year, grid%y(j), ice%u_vav_cx( j,i), tauc_analytical)
        ice%u_vav_cx(  j,i) = ice%u_vav_cx( j,i) * (1._dp - 1E-7_dp)
        ice%u_base_cx( j,i) = ice%u_vav_cx( j,i)
      END DO
      END DO
      CALL sync
    END IF
    
  END SUBROUTINE initialise_ice_model
  SUBROUTINE allocate_ice_model( grid, ice)  
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
      
    ! Allocate memory
    ! ===============
    
    ! Basic data - ice thickness, bedrock & surface elevation, sea level (geoid elevation), englacial temperature
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_a                 , ice%wHi_a                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%Hi_cx                , ice%wHi_cx                )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%Hi_cy                , ice%wHi_cy                )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%Hi_b                 , ice%wHi_b                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_a_new             , ice%wHi_a_new             )
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
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%uabs_vav_a           , ice%wuabs_vav_a           )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%u_surf_a             , ice%wu_surf_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%v_surf_a             , ice%wv_surf_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%u_surf_cx            , ice%wu_surf_cx            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%v_surf_cy            , ice%wv_surf_cy            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%uabs_surf_a          , ice%wuabs_surf_a          )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%u_base_a             , ice%wu_base_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%v_base_a             , ice%wv_base_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%u_base_cx            , ice%wu_base_cx            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%v_base_cy            , ice%wv_base_cy            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%uabs_base_a          , ice%wuabs_base_a          )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx-1, ice%u_3D_SIA_cx          , ice%wu_3D_SIA_cx          )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny-1, grid%nx  , ice%v_3D_SIA_cy          , ice%wv_3D_SIA_cy          )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%u_SSA_cx             , ice%wu_SSA_cx             )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%v_SSA_cy             , ice%wv_SSA_cy             )
    
    ! Different masks
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_land_a          , ice%wmask_land_a          )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_ocean_a         , ice%wmask_ocean_a         )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_lake_a          , ice%wmask_lake_a          )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_ice_a           , ice%wmask_ice_a           )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_sheet_a         , ice%wmask_sheet_a         )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_shelf_a         , ice%wmask_shelf_a         )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_coast_a         , ice%wmask_coast_a         )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%mask_coast_cx        , ice%wmask_coast_cx        )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%mask_coast_cy        , ice%wmask_coast_cy        )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_margin_a        , ice%wmask_margin_a        )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%mask_margin_cx       , ice%wmask_margin_cx       )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%mask_margin_cy       , ice%wmask_margin_cy       )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_gl_a            , ice%wmask_gl_a            )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%mask_gl_cx           , ice%wmask_gl_cx           )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%mask_gl_cy           , ice%wmask_gl_cy           )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_cf_a            , ice%wmask_cf_a            )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%mask_cf_cx           , ice%wmask_cf_cx           )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%mask_cf_cy           , ice%wmask_cf_cy           )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_a               , ice%wmask_a               )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%f_grnd_a             , ice%wf_grnd_a             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%f_grnd_cx            , ice%wf_grnd_cx            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%f_grnd_cy            , ice%wf_grnd_cy            )
    
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
    
    ! Ice dynamics - DIVA - physical terms
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%phi_fric_a           , ice%wphi_fric_a           )    
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%tauc_a               , ice%wtauc_a               )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%taudx_cx             , ice%wtaudx_cx             )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%taudy_cy             , ice%wtaudy_cy             )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%du_dx_b              , ice%wdu_dx_b              )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%du_dy_b              , ice%wdu_dy_b              )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%dv_dx_b              , ice%wdv_dx_b              )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%dv_dy_b              , ice%wdv_dy_b              )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx-1, ice%du_dz_3D_cx          , ice%wdu_dz_3D_cx          )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny-1, grid%nx  , ice%dv_dz_3D_cy          , ice%wdv_dz_3D_cy          )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny  , grid%nx  , ice%visc_eff_3D_a        , ice%wvisc_eff_3D_a        )
    CALL allocate_shared_dp_3D(  C%nz, grid%ny-1, grid%nx-1, ice%visc_eff_3D_b        , ice%wvisc_eff_3D_b        )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%visc_eff_int_a       , ice%wvisc_eff_int_a       )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%visc_eff_int_b       , ice%wvisc_eff_int_b       )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%N_a                  , ice%wN_a                  )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%N_b                  , ice%wN_b                  )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%beta_a               , ice%wbeta_a               )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%beta_cx              , ice%wbeta_cx              )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%beta_cy              , ice%wbeta_cy              )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%F2_a                 , ice%wF2_a                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%beta_eff_a           , ice%wbeta_eff_a           )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%beta_eff_cx          , ice%wbeta_eff_cx          )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%beta_eff_cy          , ice%wbeta_eff_cy          )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%taub_cx              , ice%wtaub_cx              )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%taub_cy              , ice%wtaub_cy              )
    CALL allocate_shared_dp_3D( C%nz,  grid%ny  , grid%nx  , ice%F1_3D_a              , ice%wF1_3D_a              )
        
    ! Ice dynamics - DIVA - solver data
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%u_cx_prev            , ice%wu_cx_prev            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%v_cy_prev            , ice%wv_cy_prev            )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%DIVA_solve_mask_cx   , ice%wDIVA_solve_mask_cx   )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%DIVA_solve_mask_cy   , ice%wDIVA_solve_mask_cy   )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%DIVA_err_cx          , ice%wDIVA_err_cx          )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%DIVA_err_cy          , ice%wDIVA_err_cy          )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%DIVA_mask_combi      , ice%wDIVA_mask_combi      )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%DIVA_isfront_inner   , ice%wDIVA_isfront_inner   )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%DIVA_isfront_outer   , ice%wDIVA_isfront_outer   )
    ! NOTE - the matrix solver arrays are allocated in initialise_DIVA_matrix_tables
    
    ! Ice dynamics - GL flux
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%Qx_GL_cx             , ice%wQx_GL_cx             )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%Qy_GL_cy             , ice%wQy_GL_cy             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%u_GL_cx              , ice%wu_GL_cx              )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%v_GL_cy              , ice%wv_GL_cy              )
    
    ! Ice dynamics - ice thickness calculation
    IF (C%choice_ice_integration_method == 'explicit') THEN
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%Qx_cx                , ice%wQx_cx                )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%Qy_cy                , ice%wQy_cy                )
    ELSEIF (C%choice_ice_integration_method == 'implicit') THEN
      CALL initialise_implicit_ice_thickness_matrix_tables( grid, ice)
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_ice_integration_methodt "', TRIM(C%choice_ice_integration_method), '" not implemented in allocate_ice_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Ice dynamics - calving
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%float_margin_frac_a  , ice%wfloat_margin_frac_a  )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_actual_cf_a       , ice%wHi_actual_cf_a       )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHi_dt_calving_a     , ice%wdHi_dt_calving_a     )
    
    ! Ice dynamics - predictor/corrector ice thickness update
    CALL allocate_shared_dp_0D(                              ice%pc_zeta              , ice%wpc_zeta              )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%pc_tau               , ice%wpc_tau               )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%pc_fcb               , ice%wpc_fcb               )
    CALL allocate_shared_dp_0D(                              ice%pc_eta               , ice%wpc_eta               )
    CALL allocate_shared_dp_0D(                              ice%pc_eta_prev          , ice%wpc_eta_prev          )
    CALL allocate_shared_dp_0D(                              ice%pc_beta1             , ice%wpc_beta1             )
    CALL allocate_shared_dp_0D(                              ice%pc_beta2             , ice%wpc_beta2             )
    CALL allocate_shared_dp_0D(                              ice%pc_beta3             , ice%wpc_beta3             )
    CALL allocate_shared_dp_0D(                              ice%pc_beta4             , ice%wpc_beta4             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%pc_f1                , ice%wpc_f1                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%pc_f2                , ice%wpc_f2                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%pc_f3                , ice%wpc_f3                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%pc_f4                , ice%wpc_f4                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_old               , ice%wHi_old               )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_pred              , ice%wHi_pred              )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_corr              , ice%wHi_corr              )
    
    ! Thermodynamics
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%frictional_heating_a , ice%wfrictional_heating_a )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%GHF_a                , ice%wGHF_a                )
    
    ! Isotopes
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_a_prev            , ice%wHi_a_prev            )
    
  END SUBROUTINE allocate_ice_model
  
END MODULE ice_dynamics_module
