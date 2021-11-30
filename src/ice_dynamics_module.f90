MODULE ice_dynamics_module

  ! Contains all the routines needed to calculate ice-sheet geometry at the next time
  ! step, including routines to determine said time step.
  ! NOTE: routines for calculating ice velocities             have been moved to the ice_velocity_module;
  !       routines for integrating the ice thickness equation have been moved to the ice_thickness_module.

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_model_region, type_grid, type_ice_model, type_reference_geometry
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             SSA_Schoof2006_analytical_solution, vertical_average, surface_elevation
  USE ice_velocity_module,             ONLY: initialise_SSADIVA_solution_matrix, solve_SIA, solve_SSA, solve_DIVA
  USE ice_thickness_module,            ONLY: calc_dHi_dt, initialise_implicit_ice_thickness_matrix_tables, apply_ice_thickness_BC, &
                                             remove_unconnected_shelves
  USE general_ice_model_data_module,   ONLY: update_general_ice_model_data

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
    
    ! Calculate ice velocities with the selected ice-dynamical approximation
    ! ======================================================================
    
    IF     (C%choice_ice_dynamics == 'none') THEN
      ! Fixed ice geometry
    
    ELSEIF (C%choice_ice_dynamics == 'SIA') THEN
      ! Shallow ice approximation
    
      IF (region%time == region%t_next_SIA) THEN
    
        ! Calculate new ice velocities
        CALL solve_SIA( region%grid, region%ice)
    
        ! Calculate critical time step
        CALL calc_critical_timestep_SIA( region%grid, region%ice, dt_crit_SIA)
        IF (par%master) region%dt_crit_SIA = dt_crit_SIA
        
        ! Update timer
        IF (par%master) region%t_last_SIA = region%time
        IF (par%master) region%t_next_SIA = region%time + region%dt_crit_SIA
        CALL sync
        
      END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN
      
    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
      ! Shallow shelf approximation
    
      IF (region%time == region%t_next_SSA) THEN
    
        ! Calculate new ice velocities
        CALL solve_SSA( region%grid, region%ice)
    
        ! Calculate critical time step
        CALL calc_critical_timestep_adv( region%grid, region%ice, dt_crit_SSA)
        IF (par%master) region%dt_crit_SSA = dt_crit_SSA
        
        ! Update timer
        IF (par%master) region%t_last_SSA = region%time
        IF (par%master) region%t_next_SSA = region%time + region%dt_crit_SSA
        CALL sync
        
      END IF ! IF (ABS(region%time - region%t_next_SSA) < dt_tol) THEN
      
    ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
      ! Hybrid SIA/SSA (Bueler and Brown, 2009)
    
      IF (region%time == region%t_next_SIA) THEN
    
        ! Calculate new ice velocities
        CALL solve_SIA( region%grid, region%ice)
    
        ! Calculate critical time step
        CALL calc_critical_timestep_SIA( region%grid, region%ice, dt_crit_SIA)
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
        CALL calc_critical_timestep_adv( region%grid, region%ice, dt_crit_SSA)
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
      
    !IF (par%master) WRITE(0,'(A,F7.4,A,F7.4,A,F7.4)') 'dt_crit_SIA = ', dt_crit_SIA, ', dt_crit_SSA = ', dt_crit_SSA, ', dt = ', region%dt
    
    ! Calculate new ice geometry
    ! ==========================
    
    IF (C%choice_ice_dynamics == 'none') THEN
      ! Fixed ice geometry
      region%ice%dHi_dt_a( :,i1:i2) = 0._dp
      CALL sync
    ELSE
      CALL calc_dHi_dt( region%grid, region%ice, region%SMB, region%BMB, region%dt, region%mask_noice)
    END IF
    
  END SUBROUTINE run_ice_dynamics_direct
  SUBROUTINE run_ice_dynamics_pc( region, t_end)
    ! Ice dynamics and time-stepping with the predictor/correct method
    ! (adopted from Yelmo, originally based on Cheng et al., 2017)
    
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
    IF     (C%choice_ice_dynamics == 'none') THEN
      region%ice%dHi_dt_a( :,i1:i2) = 0._dp
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
      IF (par%master) WRITE(0,*) 'run_ice_dynamics_pc - ERROR: unknown choice_ice_dynamics "', C%choice_ice_dynamics, '"'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    IF (do_update_ice_velocity) THEN
    
      ! Calculate time step based on the truncation error in ice thickness (Robinson et al., 2020, Eq. 33)
      CALL calc_critical_timestep_adv( region%grid, region%ice, dt_crit_adv)
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
      CALL calc_dHi_dt( region%grid, region%ice, region%SMB, region%BMB, region%dt_crit_ice, region%mask_noice)
      region%ice%pc_f1(   :,i1:i2) = region%ice%dHi_dt_a( :,i1:i2)
      region%ice%Hi_pred( :,i1:i2) = MAX(0._dp, region%ice%Hi_a(     :,i1:i2) + region%dt_crit_ice * region%ice%dHi_dt_a( :,i1:i2))
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
        IF (par%master) WRITE(0,*) 'run_ice_dynamics_pc - ERROR: unknown choice_ice_dynamics "', C%choice_ice_dynamics, '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    
      ! Corrector step
      ! ==============
      
      ! Calculate "corrected" ice thickness based on new velocities
      CALL calc_dHi_dt( region%grid, region%ice, region%SMB, region%BMB, region%dt_crit_ice, region%mask_noice)
      region%ice%pc_f3(   :,i1:i2) = region%ice%dHi_dt_a( :,i1:i2)
      region%ice%pc_f4(   :,i1:i2) = region%ice%pc_f1(    :,i1:i2)
      region%ice%Hi_corr( :,i1:i2) = MAX(0._dp, region%ice%Hi_old(   :,i1:i2) + 0.5_dp * region%dt_crit_ice * (region%ice%pc_f3( :,i1:i2) + region%ice%pc_f4( :,i1:i2))      )
      CALL sync
  
      ! Determine truncation error
      CALL calc_pc_truncation_error( region%grid, region%ice, region%dt_crit_ice, region%dt_prev)
    
      ! Go back to old ice thickness. Run all the other modules (climate, SMB, BMB, thermodynamics, etc.)
      ! and only go to new (corrected) ice thickness at the end of this time loop.
      region%ice%Hi_a(     :,i1:i2) = region%ice%Hi_old(   :,i1:i2)
      region%ice%dHi_dt_a( :,i1:i2) = (region%ice%Hi_corr( :,i1:i2) - region%ice%Hi_a( :,i1:i2)) / region%dt_crit_ice
      CALL sync
    
    END IF ! IF (do_update_ice_velocity) THEN
      
    ! Adjust the time step to prevent overshooting other model components (thermodynamics, SMB, output, etc.)
    CALL determine_timesteps_and_actions( region, t_end)
    
    ! Calculate ice thickness at the end of this model loop
    region%ice%Hi_tplusdt_a( :,i1:i2) = region%ice%Hi_a( :,i1:i2) + region%dt * region%ice%dHi_dt_a( :,i1:i2)
    
    !IF (par%master) WRITE(0,'(A,F7.4,A,F7.4,A,F7.4)') 'dt_crit_adv = ', dt_crit_adv, ', dt_from_pc = ', dt_from_pc, ', dt = ', region%dt
    
  END SUBROUTINE run_ice_dynamics_pc
  
  ! Time stepping
  SUBROUTINE calc_critical_timestep_SIA( grid, ice, dt_crit_SIA)
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
        
  END SUBROUTINE calc_critical_timestep_SIA
  SUBROUTINE calc_critical_timestep_adv( grid, ice, dt_crit_adv)
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
        
  END SUBROUTINE calc_critical_timestep_adv
  SUBROUTINE calc_pc_truncation_error( grid, ice, dt, dt_prev)
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
        
    ! Find maximum truncation error
    eta_proc = C%pc_eta_min
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    
      ! Calculate truncation error (Robinson et al., 2020, Eq. 32)
      ice%pc_tau( j,i) = ABS( zeta * (ice%Hi_corr( j,i) - ice%Hi_pred( j,i)) / ((3._dp * zeta + 3._dp) * dt))
    
!      IF (ice%mask_sheet_a( j,i) == 1 .AND. ice%mask_gl_a( j,i) == 0) THEN
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
        
  END SUBROUTINE calc_pc_truncation_error
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
      
      region%do_ocean   = .FALSE.
      IF (region%time == region%t_next_ocean) THEN
        region%do_ocean       = .TRUE.
        region%t_last_ocean   = region%time
        region%t_next_ocean   = region%t_last_ocean + C%dt_ocean
      END IF
      t_next = MIN( t_next, region%t_next_ocean)
      
      region%do_SMB     = .FALSE.
      IF (region%time == region%t_next_SMB) THEN
        region%do_SMB         = .TRUE.
        region%t_last_SMB     = region%time
        region%t_next_SMB     = region%t_last_SMB + C%dt_SMB
      END IF
      t_next = MIN( t_next, region%t_next_SMB)
      
      region%do_BMB     = .FALSE.
      IF (C%do_asynchronous_BMB) THEN
        IF (region%time == region%t_next_BMB) THEN
          region%do_BMB         = .TRUE.
          region%t_last_BMB     = region%time
          region%t_next_BMB     = region%t_last_BMB + C%dt_BMB
        END IF
        t_next = MIN( t_next, region%t_next_BMB)
      ELSE
        ! Don't use separate timestepping for the BMB; just run it in every ice dynamics time step
        region%do_BMB = .TRUE.
      END IF
      
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
  
  ! Administration: allocation and initialisation
  SUBROUTINE initialise_ice_model( grid, ice, refgeo_init)
    ! Allocate shared memory for all the data fields of the ice dynamical module, and
    ! initialise some of them
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: tauc_analytical
    
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
    
    ! Make sure we already start with correct boundary conditions
    CALL apply_ice_thickness_BC(        grid, ice, C%dt_min)
    CALL update_general_ice_model_data( grid, ice)
    CALL remove_unconnected_shelves(    grid, ice, C%dt_min)
    CALL update_general_ice_model_data( grid, ice)
    
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
    
    IF (C%do_benchmark_experiment .AND. ( &
          C%choice_benchmark_experiment == 'MISMIPplus' .OR. &
          C%choice_benchmark_experiment == 'MISOMIP1')) THEN
      ! Set the ice flow factor only during initialisation; don't allow the "ice_physical_properties" routine
      ! to change it, but instead let the tune-for-GL-position routine do that
      ice%A_flow_3D_a(  :,:,grid%i1:grid%i2) = C%MISMIPplus_A_flow_initial
      ice%A_flow_vav_a(   :,grid%i1:grid%i2) = C%MISMIPplus_A_flow_initial
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
    
    ! Ice dynamics - basal conditions and driving stress
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%taudx_cx             , ice%wtaudx_cx             )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%taudy_cy             , ice%wtaudy_cy             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%phi_fric_a           , ice%wphi_fric_a           )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%tauc_a               , ice%wtauc_a               )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%A_slid_a             , ice%wA_slid_a             )
        
    ! Ice dynamics - physical terms in the SSA/DIVA
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
    IF     (C%choice_ice_integration_method == 'explicit') THEN
      CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%Qx_cx                , ice%wQx_cx                )
      CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%Qy_cy                , ice%wQy_cy                )
    ELSEIF (C%choice_ice_integration_method == 'semi-implicit') THEN
      CALL initialise_implicit_ice_thickness_matrix_tables( grid, ice)
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_ice_integration_methodt "', TRIM(C%choice_ice_integration_method), '" not implemented in allocate_ice_model!'
     CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Ice dynamics - calving
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%float_margin_frac_a  , ice%wfloat_margin_frac_a  )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_actual_cf_a       , ice%wHi_actual_cf_a       )
    
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
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_ice_a_prev      , ice%wmask_ice_a_prev      )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%frictional_heating_a , ice%wfrictional_heating_a )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%GHF_a                , ice%wGHF_a                )
    
    ! Isotopes
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_a_prev            , ice%wHi_a_prev            )
    
  END SUBROUTINE allocate_ice_model
  
END MODULE ice_dynamics_module
