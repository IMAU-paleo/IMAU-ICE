MODULE ice_velocity_module
  ! Contains all the routines needed to calculate instantaneous ice velocities for the
  ! current modelled ice-sheet geometry.

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_model_region, type_grid, type_ice_model
  USE netcdf_module,                   ONLY: debug, write_to_debug_file 
  USE utilities_module,                ONLY: vertical_integration_from_bottom_to_zeta, vertical_average, vertical_integrate, &
                                             SSA_Schoof2006_analytical_solution, initialise_matrix_equation_CSR, &
                                             solve_matrix_equation_CSR, check_CSR_for_double_entries
  USE derivatives_and_grids_module,    ONLY: ddx_cx_to_b_2D, ddy_cx_to_b_2D, ddy_cy_to_b_2D, ddx_cy_to_b_2D, &
                                             map_cx_to_a_2D, map_cy_to_a_2D, map_cx_to_a_3D, map_cy_to_a_3D
  USE general_ice_model_data_module,   ONLY: is_floating

  IMPLICIT NONE
  
CONTAINS
  
  ! The different ice-dynamical solvers that can be called from run_ice_model
  SUBROUTINE solve_SIA( grid, ice)
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
      
    ! Local variables:    
    INTEGER                                            :: i,j
    REAL(dp)                                           :: D_0
    REAL(dp), DIMENSION(C%nZ)                          :: D_deformation
    REAL(dp), DIMENSION(C%nZ)                          :: D_SIA_3D
    REAL(dp), PARAMETER                                :: D_uv_3D_cutoff = -1E5_dp
    
    ! Check that this routine is called correctly
    IF (.NOT. (C%choice_ice_dynamics == 'SIA' .OR. C%choice_ice_dynamics == 'SIA/SSA')) THEN
      IF (par%master) WRITE(0,*) 'ERROR - solve_SIA should only be called when choice_ice_dynamics is set to SIA or SIA/SSA!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Calculate 3D horizontal velocities
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! u
      IF (i < grid%nx) THEN
      
        IF (ice%mask_sheet_a( j,i) == 0 .AND. ice%mask_sheet_a( j,i+1) == 0) THEN
    
          ice%u_3D_SIA_cx( :,j,i) = 0._dp
        
        ELSE ! IF (ice%mask_sheet_a( j,i) == 0 .AND. ice%mask_sheet_a( j,i+1) == 0) THEN
          
          D_0           = (ice_density * grav * ice%Hi_cx( j,i))**C%n_flow * ((ice%dHs_dx_cx( j,i)**2 + ice%dHs_dy_cx( j,i)**2))**((C%n_flow - 1._dp) / 2._dp)
          D_deformation = vertical_integration_from_bottom_to_zeta( ice%A_flow_3D_cx( :,j,i) * C%zeta**C%n_flow)
          D_deformation = 2._dp * ice%Hi_cx( j,i) * D_deformation
          D_SIA_3D      = MAX(D_0 * D_deformation, D_uv_3D_cutoff)
         
          ice%u_3D_SIA_cx( :,j,i) = D_SIA_3D * ice%dHs_dx_cx( j,i)
        
        END IF ! IF (ice%mask_sheet_a( j,i) == 0 .AND. ice%mask_sheet_a( j,i+1) == 0) THEN
      
      END IF ! IF (i < grid%nx) THEN
      
      ! v
      IF (j < grid%ny) THEN
      
        IF (ice%mask_sheet_a( j,i) == 0 .AND. ice%mask_sheet_a( j+1,i) == 0) THEN
    
          ice%v_3D_SIA_cy( :,j,i) = 0._dp
        
        ELSE
          
          D_0           = (ice_density * grav * ice%Hi_cy( j,i))**C%n_flow * ((ice%dHs_dx_cy( j,i)**2 + ice%dHs_dy_cy( j,i)**2))**((C%n_flow - 1._dp) / 2._dp)
          D_deformation = vertical_integration_from_bottom_to_zeta( ice%A_flow_3D_cy( :,j,i) * C%zeta**C%n_flow)
          D_deformation = 2._dp * ice%Hi_cy( j,i) * D_deformation
          D_SIA_3D      = MAX(D_0 * D_deformation, D_uv_3D_cutoff)
         
          ice%v_3D_SIA_cy( :,j,i) = D_SIA_3D * ice%dHs_dy_cy( j,i)
        
        END IF ! IF (ice%mask_sheet_a( j,i) == 0 .AND. ice%mask_sheet_a( j+1,i) == 0) THEN
        
      END IF ! IF (j < grid%ny) THEN
      
    END DO
    END DO
    CALL sync
    
    ! Calculate secondary velocities (surface, base, etc.)
    CALL calculate_secondary_velocities( grid, ice)

  END SUBROUTINE solve_SIA
  SUBROUTINE solve_SSA( grid, ice)
    ! Calculate ice velocities using the SSA. Velocities are calculated on the staggered
    ! cx/cy-grids, using the discretisation scheme adopted from Yelmo/SICOPOLIS.
    !
    ! At the ice margin (both grounded and floating), stress boundary conditions are prescribed,
    ! instead of the old PISM-based approach of "very thin ice everywhere".
    !
    ! Also, instead of the old "grid" SOR solver, we now write the equation explicitly as
    ! matrix equation (in compressed sparse row [CSR] format) and solve that. Due to the
    ! different ordering of the equations, this solver is more stable (again, adopted from
    ! Yelmo/SICOPOLIS). This also makes it possible to use an external solver in the future.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    LOGICAL                                            :: set_velocities_to_zero
    LOGICAL                                            :: has_converged
    INTEGER                                            :: viscosity_iteration_i
    REAL(dp)                                           :: resid_UV
    REAL(dp)                                           :: umax_analytical, tauc_analytical
    REAL(dp)                                           :: t0
    
    t0 = MPI_WTIME()
    
    ! Check that this routine is called correctly
    IF (.NOT. (C%choice_ice_dynamics == 'SSA' .OR. C%choice_ice_dynamics == 'SIA/SSA')) THEN
      IF (par%master) WRITE(0,*) 'ERROR - solve_SSA should only be called when choice_ice_dynamics is set to SSA or SIA/SSA!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! If there's no grounded ice anywhere, don't bother
    set_velocities_to_zero = .FALSE.
    IF (SUM( ice%mask_sheet_a) == 0) set_velocities_to_zero = .TRUE.
    
    ! If we're prescribing no sliding, set velocities to zero
    IF (C%no_sliding) set_velocities_to_zero = .TRUE.
    
    IF (set_velocities_to_zero) THEN
      ice%u_SSA_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 0._dp
      ice%v_SSA_cy( :,grid%i1:              grid%i2 ) = 0._dp
      CALL sync
      RETURN
    END IF
    
    ! Find analytical solution for the SSA icestream experiment (used only to print numerical error to screen)
    CALL SSA_Schoof2006_analytical_solution( 0.001_dp, 2000._dp, ice%A_flow_vav_a( 1,1), 0._dp, umax_analytical, tauc_analytical)
    
    ! Calculate the driving stresses taudx, taudy
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (i < grid%nx) ice%taudx_cx( j,i) = ice_density * grav * ice%Hi_cx( j,i) * ice%dHs_dx_cx( j,i)
      IF (j < grid%ny) ice%taudy_cy( j,i) = ice_density * grav * ice%Hi_cy( j,i) * ice%dHs_dy_cy( j,i)
    END DO
    END DO
    
    ! Initialise DIVA solving masks
    CALL initialise_DIVA_solve_masks( grid, ice)
            
    ! Initially set error very high 
    ice%DIVA_err_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 1E5_dp
    ice%DIVA_err_cy( :,grid%i1:              grid%i2 ) = 1E5_dp
    CALL sync
    
    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged = .FALSE.
    viscosity_iteration: DO WHILE (.NOT. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1
      
      ! Calculate the effective viscosity and the product term N = eta * H
      CALL calc_effective_viscosity( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy)
            
      ! Calculate the sliding term beta (on both the A and Cx/Cy grids)
      CALL calc_sliding_term_beta( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy)
      
      ! Set beta_eff equal to beta, this turns the DIVA into the SSA
      ice%beta_eff_a( :,grid%i1:grid%i2) = ice%beta_a( :,grid%i1:grid%i2)
      CALL sync
      
      ! Get beta_eff on the staggered grids
      CALL stagger_beta_eff( grid, ice)
      
      ! Update the SICOPOLIS DIVA solving masks
      IF (viscosity_iteration_i > 1) THEN
        CALL update_DIVA_solve_masks( grid, ice)
      END IF
      
      ! Store the previous solution so we can check for convergence later
      ice%u_cx_prev( :,grid%i1:MIN(grid%nx-1,grid%i2)) = ice%u_SSA_cx( :,grid%i1:MIN(grid%nx-1,grid%i2))
      ice%v_cy_prev( :,grid%i1:              grid%i2 ) = ice%v_SSA_cy( :,grid%i1:              grid%i2 )
      CALL sync
      
      ! Solve the linearised DIVA with the SICOPOLIS solver
      CALL solve_DIVA_stag_linearised( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy)

      ! Apply velocity limits (both overflow and underflow) for improved stability
      CALL apply_velocity_limits( grid, ice%u_SSA_cx)
      CALL apply_velocity_limits( grid, ice%v_SSA_cy)
      
      ! "relax" subsequent viscosity iterations for improved stability
      CALL relax_DIVA_visc_iterations( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy, C%DIVA_visc_it_relax)
      
      ! Check if the viscosity iteration has converged
      CALL calculate_visc_iter_UV_resid( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy, resid_UV)
      !IF (par%master) WRITE(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ': resid_UV = ', resid_UV, ', u = [', MINVAL(ice%u_SSA_cx), ' - ', MAXVAL(ice%u_SSA_cx), ']'
      
      IF (par%master .AND. C%do_benchmark_experiment .AND. C%choice_benchmark_experiment == 'SSA_icestream') &
        WRITE(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ': err = ', ABS(1._dp - MAXVAL(ice%u_vav_cx) / umax_analytical), ': resid_UV = ', resid_UV

      has_converged = .FALSE.
      IF     (resid_UV < C%DIVA_visc_it_norm_dUV_tol) THEN
        has_converged = .TRUE.
      ELSEIF (viscosity_iteration_i >= C%DIVA_visc_it_nit) THEN
        has_converged = .TRUE.
      END IF
      
      ! If needed, estimate the error in the velocity fields so we can
      ! update the SICOPOLIS-style DIVA solving masks (for better efficiency)
      IF (.NOT. has_converged) CALL estimate_visc_iter_UV_errors( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy)
      
    END DO viscosity_iteration
    
    ! Calculate secondary velocities (surface, base, etc.)
    CALL calculate_secondary_velocities( grid, ice)
    
    !IF (par%master) WRITE(0,*) '   ', TRIM(C%DIVA_choice_matrix_solver), ' solved ', TRIM(C%choice_ice_dynamics), ' in ', MPI_WTIME()-t0, ' seconds.'
    
  END SUBROUTINE solve_SSA
  SUBROUTINE solve_DIVA( grid, ice)
    ! Calculate ice velocities using the DIVA. Velocities are calculated on the staggered
    ! cx/cy-grids, using the discretisation scheme adopted from Yelmo/SICOPOLIS.
    !
    ! At the ice margin (both grounded and floating), stress boundary conditions are prescribed,
    ! instead of the old PISM-based approach of "very thin ice everywhere".
    !
    ! Also, instead of the old "grid" SOR solver, we now write the equation explicitly as
    ! matrix equation (in compressed sparse row [CSR] format) and solve that. Due to the
    ! different ordering of the equations, this solver is more stable (again, adopted from
    ! Yelmo/SICOPOLIS). This also makes it possible to use an external solver in the future.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    LOGICAL                                            :: set_velocities_to_zero
    LOGICAL                                            :: has_converged
    INTEGER                                            :: viscosity_iteration_i
    REAL(dp)                                           :: resid_UV
    REAL(dp)                                           :: t0
    
    t0 = MPI_WTIME()
    
    ! Check that this routine is called correctly
    IF (.NOT. C%choice_ice_dynamics == 'DIVA') THEN
      IF (par%master) WRITE(0,*) 'ERROR - solve_DIVA should only be called when choice_ice_dynamics is set to DIVA!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! If there's no grounded ice anywhere, don't bother
    set_velocities_to_zero = .FALSE.
    IF (SUM( ice%mask_sheet_a) == 0) set_velocities_to_zero = .TRUE.
    
    IF (set_velocities_to_zero) THEN
      ice%u_vav_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 0._dp
      ice%v_vav_cy( :,grid%i1:              grid%i2 ) = 0._dp
      CALL sync
      RETURN
    END IF
    
    ! Calculate the driving stresses taudx, taudy
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (i < grid%nx) ice%taudx_cx( j,i) = ice_density * grav * ice%Hi_cx( j,i) * ice%dHs_dx_cx( j,i)
      IF (j < grid%ny) ice%taudy_cy( j,i) = ice_density * grav * ice%Hi_cy( j,i) * ice%dHs_dy_cy( j,i)
    END DO
    END DO
    
    ! If needed, calculate analytical GL flux
    IF (C%use_analytical_GL_flux) CALL calculate_GL_flux( grid, ice)
    
    ! Initialise DIVA solving masks
    CALL initialise_DIVA_solve_masks( grid, ice)
            
    ! Initially set error very high 
    ice%DIVA_err_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 1E5_dp
    ice%DIVA_err_cy( :,grid%i1:              grid%i2 ) = 1E5_dp
    CALL sync
    
    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged = .FALSE.
    viscosity_iteration: DO WHILE (.NOT. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1
      
      ! Calculate the vertical shear strain rates
      CALL calc_vertical_shear_strain_rates( grid, ice)
      
      ! Calculate the effective viscosity and the product term N = eta * H
      CALL calc_effective_viscosity( grid, ice, ice%u_vav_cx, ice%v_vav_cy)
            
      ! Calculate the sliding term beta (on both the A and Cx/Cy grids)
      CALL calc_sliding_term_beta( grid, ice, ice%u_vav_cx, ice%v_vav_cy)
      
      ! Calculate the F-integral F2
      CALL calc_F_integral( grid, ice, n = 2._dp)
      
      ! Calculate beta_eff
      CALL calc_beta_eff( grid, ice)
      
      ! Get beta_eff on the staggered grids
      CALL stagger_beta_eff( grid, ice)
      
      ! Update the SICOPOLIS DIVA solving masks
      IF (viscosity_iteration_i > 1) THEN
        CALL update_DIVA_solve_masks( grid, ice)
      END IF
      
      ! Store the previous solution so we can check for convergence later
      ice%u_cx_prev( :,grid%i1:MIN(grid%nx-1,grid%i2)) = ice%u_vav_cx( :,grid%i1:MIN(grid%nx-1,grid%i2))
      ice%v_cy_prev( :,grid%i1:              grid%i2 ) = ice%v_vav_cy( :,grid%i1:              grid%i2 )
      CALL sync
      
      ! Solve the linearised DIVA with the SICOPOLIS solver
      CALL solve_DIVA_stag_linearised( grid, ice, ice%u_vav_cx, ice%v_vav_cy)

      ! Apply velocity limits (both overflow and underflow) for improved stability
      CALL apply_velocity_limits( grid, ice%u_vav_cx)
      CALL apply_velocity_limits( grid, ice%v_vav_cy)
      
      ! "relax" subsequent viscosity iterations for improved stability
      CALL relax_DIVA_visc_iterations( grid, ice, ice%u_vav_cx, ice%v_vav_cy, C%DIVA_visc_it_relax)
      
      ! Check if the viscosity iteration has converged
      CALL calculate_visc_iter_UV_resid( grid, ice, ice%u_vav_cx, ice%v_vav_cy, resid_UV)
     ! IF (par%master) WRITE(0,*) '    DIVA - viscosity iteration ', viscosity_iteration_i, ': resid_UV = ', resid_UV

      has_converged = .FALSE.
      IF     (resid_UV < C%DIVA_visc_it_norm_dUV_tol) THEN
        has_converged = .TRUE.
      ELSEIF (viscosity_iteration_i >= C%DIVA_visc_it_nit) THEN
        has_converged = .TRUE.
      END IF
      
      ! If needed, estimate the error in the velocity fields so we can
      ! update the SICOPOLIS-style DIVA solving masks (for better efficiency)
      IF (.NOT. has_converged) CALL estimate_visc_iter_UV_errors( grid, ice, ice%u_vav_cx, ice%v_vav_cy)

      ! Calculate basal stress 
      CALL calc_basal_stress( grid, ice)
      
      ! Calculate basal velocity from depth-averaged solution and basal stress
      CALL calc_basal_velocities( grid, ice)
      
    END DO viscosity_iteration
    
    ! Calculate full 3D velocities
    CALL calc_3D_horizontal_velocities_DIVA( grid, ice)
    
    ! Calculate secondary velocities (surface, base, etc.)
    CALL calculate_secondary_velocities( grid, ice)
    
    !IF (par%master) WRITE(0,*) '   ', TRIM(C%DIVA_choice_matrix_solver), ' solved ', TRIM(C%choice_ice_dynamics), ' in ', MPI_WTIME()-t0, ' seconds.'
    
  END SUBROUTINE solve_DIVA
  SUBROUTINE calculate_secondary_velocities( grid, ice)
    ! Calculate secondary velocities
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j,k
    
    IF (C%choice_ice_dynamics == 'SIA') THEN
      ! No SSA or sliding, just the SIA
      
      ! Set 3D velocity field equal to SIA answer
      ice%u_3D_cx( :,:,grid%i1:MIN(grid%nx-1,grid%i2)) = ice%u_3D_SIA_cx( :,:,grid%i1:MIN(grid%nx-1,grid%i2))
      ice%v_3D_cy( :,:,grid%i1:              grid%i2 ) = ice%v_3D_SIA_cy( :,:,grid%i1:              grid%i2 )
      CALL sync
      
      ! Basal velocity is zero
      ice%u_base_cx(   :,grid%i1:MIN(grid%nx-1,grid%i2)) = 0._dp
      ice%v_base_cy(   :,grid%i1:              grid%i2 ) = 0._dp
      CALL sync
    
      ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
      CALL calc_3D_vertical_velocities( grid, ice)
      
      ! Get surface velocity from the 3D fields
      ice%u_surf_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = ice%u_3D_cx( 1,:,grid%i1:MIN(grid%nx-1,grid%i2))
      ice%v_surf_cy( :,grid%i1:              grid%i2 ) = ice%v_3D_cy( 1,:,grid%i1:              grid%i2 )
      CALL sync
      
      ! Get vertically averaged velocities
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (i < grid%nx) ice%u_vav_cx( j,i) = vertical_average( ice%u_3D_cx( :,j,i))
        IF (j < grid%ny) ice%v_vav_cy( j,i) = vertical_average( ice%v_3D_cy( :,j,i))
      END DO
      END DO
      CALL sync
      
    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
      ! No SIA, just the SSA
      
      ! Set basal velocity equal to SSA answer
      ice%u_base_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = ice%u_SSA_cx( :,grid%i1:MIN(grid%nx-1,grid%i2))
      ice%v_base_cy( :,grid%i1:              grid%i2 ) = ice%v_SSA_cy( :,grid%i1:              grid%i2 )
      CALL sync
      
      ! No vertical variations in velocity
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (i < grid%nx) ice%u_vav_cx(    j,i) = ice%u_SSA_cx( j,i)
        IF (j < grid%ny) ice%v_vav_cy(    j,i) = ice%v_SSA_cy( j,i)
        IF (i < grid%nx) ice%u_surf_cx(   j,i) = ice%u_SSA_cx( j,i)
        IF (j < grid%ny) ice%v_surf_cy(   j,i) = ice%v_SSA_cy( j,i)
        IF (i < grid%nx) ice%u_3D_cx(   :,j,i) = ice%u_SSA_cx( j,i)
        IF (j < grid%ny) ice%v_3D_cy(   :,j,i) = ice%v_SSA_cy( j,i)
      END DO
      END DO
      CALL sync
      
    ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
      ! Hybrid SIA/SSA
      
      ! Add the two fields together
      DO k = 1, C%nz
        ice%u_3D_cx( k,:,grid%i1:MIN(grid%nx-1,grid%i2)) = ice%u_3D_SIA_cx( k,:,grid%i1:MIN(grid%nx-1,grid%i2)) + ice%u_SSA_cx( :,grid%i1:MIN(grid%nx-1,grid%i2))
        ice%v_3D_cy( k,:,grid%i1:              grid%i2 ) = ice%v_3D_SIA_cy( k,:,grid%i1:              grid%i2 ) + ice%v_SSA_cy( :,grid%i1:              grid%i2 )
      END DO
      CALL sync
      
      ! Set basal velocity equal to SSA answer
      ice%u_base_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = ice%u_SSA_cx( :,grid%i1:MIN(grid%nx-1,grid%i2))
      ice%v_base_cy( :,grid%i1:              grid%i2 ) = ice%v_SSA_cy( :,grid%i1:              grid%i2 )
      CALL sync
      
      ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
      CALL calc_3D_vertical_velocities( grid, ice)
      
      ! Get surface velocity from the 3D fields
      ice%u_surf_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = ice%u_3D_cx( 1,:,grid%i1:MIN(grid%nx-1,grid%i2))
      ice%v_surf_cy( :,grid%i1:              grid%i2 ) = ice%v_3D_cy( 1,:,grid%i1:              grid%i2 )
      CALL sync
      
      ! Get vertically averaged velocities
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (i < grid%nx) ice%u_vav_cx( j,i) = vertical_average( ice%u_3D_cx( :,j,i))
        IF (j < grid%ny) ice%v_vav_cy( j,i) = vertical_average( ice%v_3D_cy( :,j,i))
      END DO
      END DO
      CALL sync
    
    ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
      ! DIVA
      
      CALL calc_3D_vertical_velocities( grid, ice)
      
      ! Get surface velocity from the 3D fields
      ice%u_surf_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = ice%u_3D_cx( 1,:,grid%i1:MIN(grid%nx-1,grid%i2))
      ice%v_surf_cy( :,grid%i1:              grid%i2 ) = ice%v_3D_cy( 1,:,grid%i1:              grid%i2 )
      CALL sync
      
    ELSE ! IF (C%choice_ice_dynamics == 'SIA')
    
      IF (par%master) WRITE(0,*) 'calculate_secondary_velocities - ERROR: unknown choice_ice_dynamics "', C%choice_ice_dynamics, '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      
    END IF ! IF (C%choice_ice_dynamics == 'SIA')
      
    ! Get velocity components on the A-grid (only used for writing to output)
    CALL map_velocity_from_cx_to_a_3D( grid, ice, ice%u_3D_cx,   ice%u_3D_a  )
    CALL map_velocity_from_cy_to_a_3D( grid, ice, ice%v_3D_cy,   ice%v_3D_a  )
    CALL map_velocity_from_cx_to_a_2D( grid, ice, ice%u_vav_cx,  ice%u_vav_a )
    CALL map_velocity_from_cy_to_a_2D( grid, ice, ice%v_vav_cy,  ice%v_vav_a )
    CALL map_velocity_from_cx_to_a_2D( grid, ice, ice%u_surf_cx, ice%u_surf_a)
    CALL map_velocity_from_cy_to_a_2D( grid, ice, ice%v_surf_cy, ice%v_surf_a)
    CALL map_velocity_from_cx_to_a_2D( grid, ice, ice%u_base_cx, ice%u_base_a)
    CALL map_velocity_from_cy_to_a_2D( grid, ice, ice%v_base_cy, ice%v_base_a)
    
    ! Get absolute velocities on the A-grid (only used for writing to output)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%uabs_vav_a(  j,i) = SQRT( ice%u_vav_a(  j,i)**2 + ice%v_vav_a(  j,i)**2)
      ice%uabs_surf_a( j,i) = SQRT( ice%u_surf_a( j,i)**2 + ice%v_surf_a( j,i)**2)
      ice%uabs_base_a( j,i) = SQRT( ice%u_base_a( j,i)**2 + ice%v_base_a( j,i)**2)
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calculate_secondary_velocities

  ! Calculating some physical terms (basal yield stress, effective viscosity, etc.)
  SUBROUTINE calc_vertical_shear_strain_rates( grid, ice)
    ! Calculate vertical shear rates (Lipscomb et al. 2019, Eq. 36)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: visc_eff_cx, visc_eff_cy
    REAL(dp), PARAMETER                                :: visc_min = 1E3_dp
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! u
      IF (i < grid%nx) THEN
        IF (ice%mask_ice_a( j,i) == 0 .AND. ice%mask_ice_a( j,i+1) == 0) THEN
          ! No ice
          ice%du_dz_3D_cx( :,j,i) = 0._dp
        ELSE
          DO k = 1, C%nZ 
            visc_eff_cx  = calc_staggered_margin( ice%visc_eff_3D_a( k,j,i), ice%visc_eff_3D_a( k,j,i+1), ice%Hi_a( j,i), ice%Hi_a( j,i+1))
            visc_eff_cx  = MAX( visc_eff_cx, visc_min)    ! For safety 
            ice%du_dz_3D_cx( k,j,i) = (ice%taub_cx( j,i) / visc_eff_cx) * C%zeta( k)
          END DO
        END IF
      END IF
      
      ! v
      IF (j < grid%ny) THEN
        IF (ice%mask_ice_a( j,i) == 0 .AND. ice%mask_ice_a( j+1,i) == 0) THEN
          ! No ice
          ice%dv_dz_3D_cy( :,j,i) = 0._dp
        ELSE
          DO k = 1, C%nZ 
            visc_eff_cy  = calc_staggered_margin( ice%visc_eff_3D_a( k,j,i), ice%visc_eff_3D_a( k,j+1,i), ice%Hi_a( j,i), ice%Hi_a( j+1,i))
            visc_eff_cy  = MAX( visc_eff_cy, visc_min)    ! For safety 
            ice%dv_dz_3D_cy( k,j,i) = (ice%taub_cy( j,i) / visc_eff_cy) * C%zeta( k)
          END DO
        END IF
      END IF
    
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calc_vertical_shear_strain_rates
  SUBROUTINE calc_effective_viscosity( grid, ice, u_cx, v_cy)
    ! Calculate 3D effective viscosity following Lipscomb et al. (2019), Eq. 2

    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_cx
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_cy
    
    ! Local variables
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: dudz_b, dvdz_b
    REAL(dp), PARAMETER                                :: epsilon_sq_0 = 1E-15_dp  
    REAL(dp)                                           :: eps_sq, A_flow_b, w
    REAL(dp), DIMENSION(C%nz)                          :: visc_eff_3D_a
    
    ! Calculate effective strain components from horizontal stretching on ab-nodes
    CALL ddx_cx_to_b_2D( grid, u_cx, ice%du_dx_b)
    CALL ddy_cx_to_b_2D( grid, u_cx, ice%du_dy_b)
    CALL ddx_cy_to_b_2D( grid, v_cy, ice%dv_dx_b)
    CALL ddy_cy_to_b_2D( grid, v_cy, ice%dv_dy_b)

    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
    
      IF (C%choice_ice_margin == 'BC') THEN
        IF (ice%mask_ice_a( j  ,i  ) == 0 .AND. &
            ice%mask_ice_a( j  ,i+1) == 0 .AND. &
            ice%mask_ice_a( j+1,i  ) == 0 .AND. &
            ice%mask_ice_a( j+1,i+1) == 0) THEN
          ! No ice
          ice%visc_eff_3D_b( :,j,i) = 0._dp
          CYCLE
        END IF
      ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
        ! In the "infinite slab" case, calculate effective viscosity everywhere
        ! (even when there's technically no ice present)
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in calc_effective_viscosity!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

      DO k = 1, C%nz

        ! Un-stagger shear terms to central aa-nodes in horizontal
        dudz_b = calc_staggered_margin( ice%du_dz_3D_cx( k,j,i), ice%du_dz_3D_cx( k,j+1,i), ice%Hi_a( j,i), ice%Hi_a( j+1,i))
        dvdz_b = calc_staggered_margin( ice%dv_dz_3D_cy( k,j,i), ice%dv_dz_3D_cy( k,j,i+1), ice%Hi_a( j,i), ice%Hi_a( j,i+1))
    
        ! Calculate the total effective strain rate from L19, Eq. 21 
        eps_sq = ice%du_dx_b( j,i)**2 + &
                 ice%dv_dy_b( j,i)**2 + &
                 ice%du_dx_b( j,i) * ice%dv_dy_b( j,i) + &
                 0.25_dp * (ice%du_dy_b( j,i) + ice%dv_dx_b( j,i))**2 + &
                 0.25_dp * (dudz_b**2 + dvdz_b**2) + &
                 epsilon_sq_0
        
        ! Get flow factor on the B-grid
        A_flow_b = 0._dp
        IF (C%choice_ice_margin == 'BC') THEN
        
          A_flow_b = 0._dp
          w        = 0._dp
          
          IF (ice%mask_ice_a( j,i) == 1) THEN
            A_flow_b = A_flow_b + ice%A_flow_3D_a( k,j,i)
            w = w + 1._dp
          END IF
          IF (ice%mask_ice_a( j,i+1) == 1) THEN
            A_flow_b = A_flow_b + ice%A_flow_3D_a( k,j,i+1)
            w = w + 1._dp
          END IF
          IF (ice%mask_ice_a( j+1,i) == 1) THEN
            A_flow_b = A_flow_b + ice%A_flow_3D_a( k,j+1,i)
            w = w + 1._dp
          END IF
          IF (ice%mask_ice_a( j+1,i+1) == 1) THEN
            A_flow_b = A_flow_b + ice%A_flow_3D_a( k,j+1,i+1)
            w = w + 1._dp
          END IF
          
          A_flow_b = A_flow_b / w
          
        ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
          
          A_flow_b = (ice%A_flow_3D_a( k,j,i) + ice%A_flow_3D_a( k,j,i+1) + ice%A_flow_3D_a( k,j+1,i) + ice%A_flow_3D_a( k,j+1,i+1)) / 4._dp
          
        ELSE
          IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in calc_effective_viscosity!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
        ! Calculate effective viscosity on ab-nodes
        ice%visc_eff_3D_b( k,j,i) = 0.5_dp * A_flow_b**(-1._dp/C%n_flow) * (eps_sq)**((1._dp - C%n_flow)/(2._dp*C%n_flow))

      END DO ! DO k = 1, C%nz
      
      ! Vertical integral
      ice%visc_eff_int_b( j,i) = vertical_integrate( ice%visc_eff_3D_b( :,j,i))
      
      ! Product term N = eta * H
      IF (C%choice_ice_margin == 'BC') THEN
        ice%N_b( j,i)  = ice%visc_eff_int_b( j,i) * ice%Hi_b( j,i)
      ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
        ice%N_b( j,i)  = ice%visc_eff_int_b( j,i) * MAX(0.1_dp,ice%Hi_b( j,i))
      END IF

    END DO  
    END DO
    CALL sync
    
    ! Unstagger back from B-grid to A-grid
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    
      IF (C%choice_ice_margin == 'BC') THEN
      
        IF (ice%mask_ice_a( j,i) == 0) THEN
          ! No ice
          ice%visc_eff_3D_a( :,j,i) = 0._dp
          CYCLE
        END IF
        
        visc_eff_3D_a = 0._dp
        w             = 0._dp
        
        IF (ice%visc_eff_3D_b( 1,j-1,i-1) > 0._dp) THEN
          visc_eff_3D_a = visc_eff_3D_a + ice%visc_eff_3D_b( :,j-1,i-1)
          w = w + 1._dp
        END IF
        IF (ice%visc_eff_3D_b( 1,j-1,i  ) > 0._dp) THEN
          visc_eff_3D_a = visc_eff_3D_a + ice%visc_eff_3D_b( :,j-1,i  )
          w = w + 1._dp
        END IF
        IF (ice%visc_eff_3D_b( 1,j  ,i-1) > 0._dp) THEN
          visc_eff_3D_a = visc_eff_3D_a + ice%visc_eff_3D_b( :,j  ,i-1)
          w = w + 1._dp
        END IF
        IF (ice%visc_eff_3D_b( 1,j  ,i  ) > 0._dp) THEN
          visc_eff_3D_a = visc_eff_3D_a + ice%visc_eff_3D_b( :,j  ,i  )
          w = w + 1._dp
        END IF
        
        ice%visc_eff_3D_a( :,j,i) = visc_eff_3D_a / w
        
      ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
        ! In the "infinite slab" case, calculate effective viscosity everywhere
        ! (even when there's technically no ice present)
        
        ice%visc_eff_3D_a( :,j,i) = (ice%visc_eff_3D_b( :,j-1,i-1) + ice%visc_eff_3D_b( :,j-1,i  ) + ice%visc_eff_3D_b( :,j  ,i-1) + ice%visc_eff_3D_b( :,j  ,i  )) / 4._dp
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in calc_effective_viscosity!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
      ! Vertical integral
      ice%visc_eff_int_a( j,i) = vertical_integrate( ice%visc_eff_3D_a( :,j,i))
      
      ! Product term N = eta * H
      IF (C%choice_ice_margin == 'BC') THEN
        ice%N_a( j,i) = ice%visc_eff_int_a( j,i) * ice%Hi_a( j,i)
      ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
        ice%N_a( j,i) = ice%visc_eff_int_a( j,i) * MAX(0.1_dp,ice%Hi_a( j,i))
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Domain boundaries
    IF (par%master) THEN
      ice%N_a(            :      ,1      ) = ice%N_a(            :        ,2        )
      ice%N_a(            :      ,grid%nx) = ice%N_a(            :        ,grid%nx-1)
      ice%N_a(            1      ,:      ) = ice%N_a(            2        ,:        )
      ice%N_a(            grid%ny,:      ) = ice%N_a(            grid%ny-1,:        )
    END IF
    CALL sync

  END SUBROUTINE calc_effective_viscosity
  SUBROUTINE stagger_effective_viscosity_a_b( grid, ice)
    ! Get N_b from N_a

    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables
    INTEGER                                            :: i,j
    REAL(dp)                                           :: w

    ! Stagger viscosity only using contributions from neighbors that have ice  
    DO i = 1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1

      w = 0._dp
      ice%N_b( j,i) = 0._dp
  
      IF (ice%mask_ice_a( j,i) == 1) THEN
        w = w + 1._dp
        ice%N_b( j,i) = ice%N_b( j,i) + ice%N_a( j,i)
      END IF
  
      IF (ice%mask_ice_a( j,i+1) == 1) THEN
        w = w + 1._dp
        ice%N_b( j,i) = ice%N_b( j,i) + ice%N_a( j,i+1)
      END IF
  
      IF (ice%mask_ice_a( j+1,i) == 1) THEN
        w = w + 1._dp
        ice%N_b( j,i) = ice%N_b( j,i) + ice%N_a( j+1,i)
      END IF
  
      IF (ice%mask_ice_a( j+1,i+1) == 1) THEN
        w = w + 1._dp
        ice%N_b( j,i) = ice%N_b( j,i) + ice%N_a( j+1,i+1)
      END IF
  
      IF (w > 0._dp) ice%N_b( j,i) = ice%N_b( j,i) / w

    END DO
    END DO
    CALL sync

  END SUBROUTINE stagger_effective_viscosity_a_b
  SUBROUTINE calc_sliding_term_beta( grid, ice, u_cx, v_cy)
    ! Calculate the sliding term beta = tau_c / norm([U,V]) in the DIVA
    ! First calculate it on the regular A-grid, then stagger it to the cx/cy-grids
    ! using a "sub-grid" interpolation scheme at the grounding line (adopted from Yelmo)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_cx
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_cy
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: u_a, v_a
    REAL(dp), PARAMETER                                :: delta_v = 1E-3_dp                    ! Normalisation parameter to prevent errors when velocity is zero
    REAL(dp), PARAMETER                                :: q_plastic            = 0.30_dp       ! Parameter used for basal stress (inverse of m_flow)
    REAL(dp), PARAMETER                                :: u_threshold          = 100._dp       ! scaling of tau_yield to get the correct unit (function of q_plastic)
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        ! No exceptions here; either no sliding is included, or the default model sliding law is used.
        
      ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_C') THEN
        
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          ice%beta_a( j,i) = 1000._dp + 1000._dp * SIN( 2._dp * pi * grid%x( i) / C%ISMIP_HOM_L) * SIN( 2._dp * pi * grid%y( j) / C%ISMIP_HOM_L)
        END DO
        END DO
        CALL sync
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_D') THEN
        
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          ice%beta_a( j,i) = 1000._dp + 1000._dp * SIN( 2._dp * pi * grid%x( i) / C%ISMIP_HOM_L)
        END DO
        END DO
        CALL sync
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_E') THEN
        ! Use tau_c as a guide for the slip zone in Haut Glacier d'Arolla
        
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          IF (ice%tauc_a( j,i) == 0._dp) THEN
            ice%beta_a( j,i) = 1E30_dp
          ELSE
            ice%beta_a( j,i) = 0._dp
          END IF
        END DO
        END DO
        CALL sync
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
        ! Use tau_c as a guide for the slip zone in Haut Glacier d'Arolla
        
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          ice%beta_a( j,i) = (ice%A_flow_vav_a( j,i) * 1000._dp)**(-1._dp)
        END DO
        END DO
        CALL sync
        RETURN
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in calc_sliding_term_beta!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF
      
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      IF (C%choice_ice_margin == 'BC') THEN
        IF (ice%mask_ice_a( j,i) == 0) THEN
          ! No ice
          ice%beta_a( j,i) = 0._dp
          CYCLE
        END IF
      ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
        ! In the "infinite slab" case, calculate effective viscosity everywhere
        ! (even when there's technically no ice present)
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in calc_sliding_term_beta!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    
      ! Get u,v on the A grid
      IF (i == 1) THEN
        u_a = u_cx( j,1)
      ELSEIF (i == grid%nx) THEN
        u_a = u_cx( j,grid%nx-1)
      ELSE
        u_a = 0.5_dp * (u_cx( j,i-1) + u_cx( j,i))
      END IF
      IF (j == 1) THEN
        v_a = v_cy( 1,i)
      ELSEIF (j == grid%ny) THEN
        v_a = v_cy( grid%ny-1,i)
      ELSE
        v_a = 0.5_dp * (v_cy( j-1,i) + v_cy( j,i))
      END IF
    
      IF (C%choice_sliding_law == 'Coulomb') THEN
        ! Coulomb-type sliding law
      
        IF (.NOT. is_floating( ice%Hi_a( j,i), ice%Hb_a( j,i), ice%SL_a( j,i))) THEN
          ! The sliding term beta, with a normalisation term (following Bueler & Brown, 2009) to prevent divide-by-zero errors.
          ice%beta_a( j,i) = ice%tauc_a( j,i) / SQRT(delta_v**2 + u_a**2 + v_a**2)
        ELSE
          ice%beta_a( j,i) = 0._dp
        END IF
      
      ELSEIF (C%choice_sliding_law == 'Coulomb_regularised') THEN
        ! Regularised Coulomb-type sliding law
      
        IF (.NOT. is_floating( ice%Hi_a( j,i), ice%Hb_a( j,i), ice%SL_a( j,i))) THEN
          ! The sliding term beta, with a normalisation term (following Bueler & Brown, 2009) to prevent divide-by-zero errors.
          ice%beta_a( j,i) = ice%tauc_a( j,i) * ( (delta_v**2 + u_a**2 + v_a**2)**(0.5_dp * (q_plastic-1._dp)) ) / (u_threshold**q_plastic)
        ELSE
          ice%beta_a( j,i) = 0._dp
        END IF
      
      ELSE ! IF (C%choice_sliding_law == 'Coulomb') THEN
        WRITE(0,*) ' ERROR: choice_sliding_law "', C%choice_sliding_law,'" not implemented in DIVA_sliding_term!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
        
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calc_sliding_term_beta
  SUBROUTINE calc_F_integral( grid, ice, n)
    ! Calculate the integral F2 in Lipscomb et al. (2019), Eq. 30
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: n
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: F_int_min
    REAL(dp), PARAMETER                                :: visc_min = 1E5_dp

    ! Determine the minimum value of F_int, to assign when H_ice == 0,
    ! since F_int should be nonzero everywhere for numerics
    F_int_min = vertical_integrate( (1._dp / visc_min) * C%zeta**n)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      IF (ice%mask_ice_a( j,i) == 1) THEN
      
        ice%F2_a( j,i) = MAX( F_int_min, vertical_integrate( (ice%Hi_a( j,i) / ice%visc_eff_3D_a( :,j,i)) * C%zeta**n))
        
      ELSE
      
        IF (C%choice_ice_margin == 'BC') THEN
          ice%F2_a( j,i) = 0._dp
        ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
          ! In the "infinite slab" case, calculate effective viscosity everywhere
          ! (even when there's technically no ice present)
          ice%F2_a( j,i) = F_int_min
        ELSE
          IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in calc_F_integral!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
      END IF

    END DO
    END DO
    CALL sync

  END SUBROUTINE calc_F_integral
  SUBROUTINE calc_beta_eff( grid, ice)
    ! Calculate the "effective basal friction" beta_eff, used in the DIVA
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j

    IF (C%no_sliding) THEN
      ! No basal sliding allowed, impose beta_eff derived from viscosity 

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ! Lipscomb et al., 2019, Eq. 35
        ice%beta_eff_a( j,i) = 1._dp / ice%F2_a( j,i)
      END DO
      END DO
      CALL sync

    ELSE
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ! Lipscomb et al., 2019, Eq. 33
        ice%beta_eff_a( j,i) = ice%beta_a( j,i) / (1._dp + ice%beta_a( j,i) * ice%F2_a( j,i))
      END DO
      END DO
      CALL sync

    END IF

  END SUBROUTINE calc_beta_eff
  SUBROUTINE stagger_beta_eff( grid, ice)
    ! "sub-grid" staggering of beta from the A to the Cx/Cy grids (adapted from Yelmo)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j

    ! acx-nodes
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny

      IF (C%choice_ice_margin == 'BC') THEN
      
        IF     (ice%Hi_a( j,i) >  0._dp .AND. ice%Hi_a( j,i+1) == 0._dp) THEN
          ! Ice to the left
          ice%beta_eff_cx( j,i) = ice%beta_eff_a( j,i)
        ELSEIF (ice%Hi_a( j,i) == 0._dp .AND. ice%Hi_a( j,i+1) >  0._dp) THEN
          ! Ice to the right
          ice%beta_eff_cx( j,i) = ice%beta_eff_a( j,i+1)
        ELSE
          ! Ice on both sides (or on neither)
          IF     (ice%f_grnd_a( j,i) >  0._dp .AND. ice%f_grnd_a( j,i+1) == 0._dp) THEN
            ! Floating to the right 
            ice%beta_eff_cx( j,i) = ice%f_grnd_cx( j,i) * ice%beta_eff_a( j,i) + (1._dp - ice%f_grnd_cx( j,i)) * ice%beta_eff_a( j,i+1)
          ELSEIF (ice%f_grnd_a( j,i) == 0._dp .AND. ice%f_grnd_a( j,i+1) >  0._dp) THEN 
            ! Floating to the left 
            ice%beta_eff_cx( j,i) = (1._dp - ice%f_grnd_cx( j,i)) * ice%beta_eff_a( j,i) + ice%f_grnd_cx( j,i) * ice%beta_eff_a( j,i+1)
          ELSEIF (ice%f_grnd_a( j,i) >  0._dp .AND. ice%f_grnd_a( j,i+1) >  0._dp) THEN 
            ! Fully grounded, simple staggering
            ice%beta_eff_cx( j,i) = 0.5_dp * (ice%beta_eff_a( j,i) + ice%beta_eff_a( j,i+1))
          ELSE 
            ! Fully floating
            ice%beta_eff_cx( j,i) = 0._dp 
          END IF
        END IF
        
      ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
        ! In the "infinite slab" case, calculate stuff everywhere
        ! (even when there's technically no ice present)
        
        IF     (ice%f_grnd_a( j,i) >  0._dp .AND. ice%f_grnd_a( j,i+1) == 0._dp) THEN
          ! Floating to the right 
          ice%beta_eff_cx( j,i) = ice%f_grnd_cx( j,i) * ice%beta_eff_a( j,i) + (1._dp - ice%f_grnd_cx( j,i)) * ice%beta_eff_a( j,i+1)
        ELSEIF (ice%f_grnd_a( j,i) == 0._dp .AND. ice%f_grnd_a( j,i+1) >  0._dp) THEN 
          ! Floating to the left 
          ice%beta_eff_cx( j,i) = (1._dp - ice%f_grnd_cx( j,i)) * ice%beta_eff_a( j,i) + ice%f_grnd_cx( j,i) * ice%beta_eff_a( j,i+1)
        ELSEIF (ice%f_grnd_a( j,i) >  0._dp .AND. ice%f_grnd_a( j,i+1) >  0._dp) THEN 
          ! Fully grounded, simple staggering
          ice%beta_eff_cx( j,i) = 0.5_dp * (ice%beta_eff_a( j,i) + ice%beta_eff_a( j,i+1))
        ELSE 
          ! Fully floating
          ice%beta_eff_cx( j,i) = 0._dp 
        END IF
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in stagger_beta_eff!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

    END DO 
    END DO
    CALL sync
    
    ! acy-nodes 
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1

      IF (C%choice_ice_margin == 'BC') THEN

        IF     (ice%Hi_a( j,i) >  0._dp .AND. ice%Hi_a( j+1,i) == 0._dp) THEN
          ! Ice to the bottom
          ice%beta_eff_cy( j,i) = ice%beta_eff_a( j,i)
        ELSEIF (ice%Hi_a( j,i) == 0._dp .AND. ice%Hi_a( j+1,i) >  0._dp) THEN
          ! Ice to the top
          ice%beta_eff_cy( j,i) = ice%beta_eff_a( j+1,i)
        ELSE
          ! Ice on both sides (or on neither)
          IF     (ice%f_grnd_a( j,i) >  0._dp .AND. ice%f_grnd_a( j+1,i) == 0._dp) THEN 
            ! Floating to the top 
            ice%beta_eff_cy( j,i) = ice%f_grnd_cy( j,i) * ice%beta_eff_a( j,i) + (1._dp - ice%f_grnd_cy( j,i)) * ice%beta_eff_a( j+1,i)
          ELSEIF (ice%f_grnd_a( j,i) == 0._dp .AND. ice%f_grnd_a( j+1,i) >  0._dp) THEN 
            ! Floating to the bottom 
            ice%beta_eff_cy( j,i) = (1._dp - ice%f_grnd_cy( j,i)) * ice%beta_eff_a( j,i) + ice%f_grnd_cy( j,i) * ice%beta_eff_a( j+1,i)
          ELSEIF (ice%f_grnd_a( j,i) >  0._dp .AND. ice%f_grnd_a( j+1,i) >  0._dp) THEN 
            ! Fully grounded, simple staggering 
            ice%beta_eff_cy( j,i) = 0.5_dp *(ice%beta_eff_a( j,i) + ice%beta_eff_a( j+1,i))
          ELSE 
            ! Fully floating 
            ice%beta_eff_cy( j,i) = 0._dp 
          END IF
        END IF
        
      ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
        ! In the "infinite slab" case, calculate stuff everywhere
        ! (even when there's technically no ice present)
        
        IF     (ice%f_grnd_a( j,i) >  0._dp .AND. ice%f_grnd_a( j+1,i) == 0._dp) THEN 
          ! Floating to the top 
          ice%beta_eff_cy( j,i) = ice%f_grnd_cy( j,i) * ice%beta_eff_a( j,i) + (1._dp - ice%f_grnd_cy( j,i)) * ice%beta_eff_a( j+1,i)
        ELSEIF (ice%f_grnd_a( j,i) == 0._dp .AND. ice%f_grnd_a( j+1,i) >  0._dp) THEN 
          ! Floating to the bottom 
          ice%beta_eff_cy( j,i) = (1._dp - ice%f_grnd_cy( j,i)) * ice%beta_eff_a( j,i) + ice%f_grnd_cy( j,i) * ice%beta_eff_a( j+1,i)
        ELSEIF (ice%f_grnd_a( j,i) >  0._dp .AND. ice%f_grnd_a( j+1,i) >  0._dp) THEN 
          ! Fully grounded, simple staggering 
          ice%beta_eff_cy( j,i) = 0.5_dp *(ice%beta_eff_a( j,i) + ice%beta_eff_a( j+1,i))
        ELSE 
          ! Fully floating 
          ice%beta_eff_cy( j,i) = 0._dp 
        END IF
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in stagger_beta_eff!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

    END DO 
    END DO 
    CALL sync
    
  END SUBROUTINE stagger_beta_eff
  SUBROUTINE calc_basal_stress( grid, ice)
    ! Calculate the basal stress resulting from sliding (friction times velocity)
    ! Note: calculated on ac-nodes.
    ! taub = beta*u (here defined with taub in the same direction as u)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (i < grid%nx) ice%taub_cx( j,i) = ice%beta_eff_cx( j,i) * ice%u_vav_cx( j,i)
      IF (j < grid%ny) ice%taub_cy( j,i) = ice%beta_eff_cy( j,i) * ice%v_vav_cy( j,i)

    END DO
    END DO
    CALL sync

  END SUBROUTINE calc_basal_stress
  SUBROUTINE calc_basal_velocities( grid, ice)
    ! Calculate basal sliding following Goldberg (2011), Eq. 34
    ! (or it can also be obtained from L19, Eq. 32 given ub*beta=taub)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: F2_stag

    IF (C%no_sliding) THEN
      ! Set basal velocities to zero 
      ! (this comes out naturally more or less with beta_eff set as above, 
      !  but ensuring basal velocity is zero adds stability)
      ice%u_base_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 0._dp
      ice%v_base_cy( :,grid%i1:              grid%i2 ) = 0._dp
      CALL sync
      RETURN
    END IF
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! x-direction
      IF (i < grid%nx) THEN

        ! Stagger the F2 integral to the cx-grid
        F2_stag = calc_staggered_margin( ice%F2_a( j,i), ice%F2_a( j,i+1), ice%Hi_a( j,i), ice%Hi_a( j,i+1))
  
        ! Calculate basal velocity component 
        ice%u_base_cx( j,i) = ice%u_vav_cx( j,i) - ice%taub_cx( j,i) * F2_stag
      
      END IF

      ! y-direction
      IF (j < grid%ny) THEN
      
        ! Stagger the F2 integral to the cy-grid
        F2_stag = calc_staggered_margin( ice%F2_a( j,i), ice%F2_a( j+1,i), ice%Hi_a( j,i), ice%Hi_a( j+1,i))
            
        ! Calculate basal velocity component 
        ice%v_base_cy( j,i) = ice%v_vav_cy( j,i) - ice%taub_cy( j,i) * F2_stag
      
      END IF

    END DO
    END DO
    CALL sync
        
  END SUBROUTINE calc_basal_velocities
  SUBROUTINE calc_3D_horizontal_velocities_DIVA( grid, ice)
    ! Caluculate the 3D horizontal velocity field (following Lipscomb et al., 2019, Eq. 29)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION( C%nz)                         :: F1_stag

    ! First calculate F1 array on aa-nodes 
    ! (performing integral before staggering seems to improve result slightly)
    ! Note: L19 define the F1 integral as purely going from the base to the surface,
    ! whereas here F1 is calculated from the base to each point in the vertical. So, 
    ! it is not technically "F1" as defined by L19, Eq. 30, except at the surface.
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
        ice%F1_3D_a( :,j,i) = vertical_integration_from_bottom_to_zeta( (-ice%Hi_a( j,i) / ice%visc_eff_3D_a( :,j,i)) * C%zeta)
    END DO
    END DO
    CALL sync

    ! Next calculate 3D horizontal velocity components 
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! x direction
      IF (i < grid%nx) THEN
      
        ! Stagger F1 column to ac-nodes
        IF (C%choice_ice_margin == 'BC') THEN
          IF     (ice%Hi_a( j,i) >  0._dp .AND. ice%Hi_a( j,i+1) == 0._dp) THEN 
            F1_stag = ice%F1_3D_a( :,j,i) 
          ELSEIF (ice%Hi_a( j,i) == 0._dp .AND. ice%Hi_a( j,i+1) >  0._dp) THEN
            F1_stag = ice%F1_3D_a( :,j,i+1)
          ELSEIF (ice%Hi_a( j,i) >  0._dp .AND. ice%Hi_a( j,i+1) >  0._dp) THEN
            F1_stag = 0.5_dp * (ice%F1_3D_a( :,j,i) + ice%F1_3D_a( :,j,i+1))
          ELSE
            F1_stag = 0._dp
          END IF 
        ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
          ! In the "infinite slab" case, calculate effective viscosity everywhere
          ! (even when there's technically no ice present)
          F1_stag = 0.5_dp * (ice%F1_3D_a( :,j,i) + ice%F1_3D_a( :,j,i+1))
        ELSE
          IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in calc_effective_viscosity!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF

        ! Calculate velocity column 
        ice%u_3D_cx( :,j,i) = ice%u_base_cx( j,i) + ice%taub_cx( j,i) * F1_stag
        
      END IF

      ! y direction
      IF (j < grid%ny) THEN
      
        ! Stagger F1 column to ac-nodes 
        IF (C%choice_ice_margin == 'BC') THEN
          IF     (ice%Hi_a( j,i) >  0._dp .AND. ice%Hi_a( j+1,i) == 0._dp) THEN 
            F1_stag = ice%F1_3D_a( :,j,i) 
          ELSEIF (ice%Hi_a( j,i) == 0._dp .AND. ice%Hi_a( j+1,i) >  0._dp) THEN
            F1_stag = ice%F1_3D_a( :,j+1,i)
          ELSEIF (ice%Hi_a( j,i) >  0._dp .AND. ice%Hi_a( j+1,i) >  0._dp) THEN
            F1_stag = 0.5_dp * (ice%F1_3D_a(:,j,i) + ice%F1_3D_a(:,j+1,i))
          ELSE
            F1_stag = 0._dp
          END IF 
        ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
          ! In the "infinite slab" case, calculate effective viscosity everywhere
          ! (even when there's technically no ice present)
          F1_stag = 0.5_dp * (ice%F1_3D_a(:,j,i) + ice%F1_3D_a(:,j+1,i))
        ELSE
          IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in calc_effective_viscosity!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF

        ! Calculate velocity column
        ice%v_3D_cy( :,j,i) = ice%v_base_cy( j,i) + ice%taub_cy( j,i) * F1_stag
        
      END IF

    END DO
    END DO
    CALL sync

  END SUBROUTINE calc_3D_horizontal_velocities_DIVA
  SUBROUTINE calc_3D_vertical_velocities( grid, ice)
    ! Use simple conservation of mass to calculate the vertical velocity w_3D
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: dHb_dx, dHb_dy, u_base_a, v_base_a
    REAL(dp)                                           :: du_dx_k,   dv_dy_k
    REAL(dp)                                           :: du_dx_kp1, dv_dy_kp1
    REAL(dp)                                           :: w1, w2, w3, w4

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! Base-of-ice slopes (not the same as bedrock slope when ice is floating!)
      dHb_dx = ice%dHs_dx_a( j,i) - ice%dHi_dx_a( j,i)
      dHb_dy = ice%dHs_dy_a( j,i) - ice%dHi_dy_a( j,i)
      
      ! Horizontal basal velocities
      IF     (i == 1) THEN
        u_base_a = ice%u_3D_cx( C%nz,j,1)
      ELSEIF (i == grid%nx) THEN
        u_base_a = ice%u_3D_cx( C%nz,j,grid%nx-1)
      ELSE
        u_base_a = 0.5_dp * (ice%u_3D_cx( C%nz,j,i-1) + ice%u_3D_cx( C%nz,j,i))
      END IF
      IF     (j == 1) THEN
        v_base_a = ice%v_3D_cy( C%nz,1,i)
      ELSEIF (j == grid%ny) THEN
        v_base_a = ice%v_3D_cy( C%nz,grid%ny-1,i)
      ELSE
        v_base_a = 0.5_dp * (ice%v_3D_cy( C%nz,j-1,i) + ice%v_3D_cy( C%nz,j,i))
      END IF
      
      ! Vertical velocity at the base
      ice%w_3D_a( C%nz,j,i) = ice%dHb_dt_a( j,i) + u_base_a * dHb_dx + v_base_a * dHb_dy 
                           
      ! The integrand is calculated half way the layer of integration at k+1/2. This integrant is multiplied with the layer thickness and added to the integral
      ! of all layers below, giving the integral up to and including this layer:
      DO k = C%nz-1, 1, -1
        
        ! Flow divergence
        IF     (i == 1) THEN
          du_dx_k   = (ice%u_3D_cx( k  ,j,2) - ice%u_3D_cx( k  ,j,1)) / grid%dx
          du_dx_kp1 = (ice%u_3D_cx( k+1,j,2) - ice%u_3D_cx( k+1,j,1)) / grid%dx
        ELSEIF (i == grid%nx) THEN
          du_dx_k   = (ice%u_3D_cx( k  ,j,grid%nx-1) - ice%u_3D_cx(k  ,j,grid%nx-2)) / grid%dx
          du_dx_kp1 = (ice%u_3D_cx( k+1,j,grid%nx-1) - ice%u_3D_cx(k+1,j,grid%nx-2)) / grid%dx
        ELSE
          du_dx_k   = (ice%u_3D_cx( k  ,j,i) - ice%u_3D_cx( k  ,j,i-1)) / grid%dx
          du_dx_kp1 = (ice%u_3D_cx( k+1,j,i) - ice%u_3D_cx( k+1,j,i-1)) / grid%dx
        END IF
        IF     (j == 1) THEN
          dv_dy_k   = (ice%v_3D_cy( k  ,2,i) - ice%v_3D_cy( k  ,1,i)) / grid%dx
          dv_dy_kp1 = (ice%v_3D_cy( k+1,2,i) - ice%v_3D_cy( k+1,1,i)) / grid%dx
        ELSEIF (j == grid%ny) THEN
          dv_dy_k   = (ice%v_3D_cy( k  ,grid%ny-1,i) - ice%v_3D_cy( k  ,grid%ny-2,i)) / grid%dx
          dv_dy_kp1 = (ice%v_3D_cy( k+1,grid%ny-1,i) - ice%v_3D_cy( k+1,grid%ny-2,i)) / grid%dx
        ELSE
          dv_dy_k   = (ice%v_3D_cy( k  ,j,i) - ice%v_3D_cy( k  ,j-1,i)) / grid%dx
          dv_dy_kp1 = (ice%v_3D_cy( k+1,j,i) - ice%v_3D_cy( k+1,j-1,i)) / grid%dx
        END IF
        
        w1 = (du_dx_k + du_dx_kp1) / 2._dp
        w2 = (dv_dy_k + dv_dy_kp1) / 2._dp   
        
        w3 = ((ice%dHs_dx_a( j,i) - 0.5_dp * (C%zeta(k+1) + C%zeta(k)) * ice%dHi_dx_a( j,i)) / MAX(0.1_dp, ice%Hi_a( j,i))) * &
             ((ice%u_3D_a( k+1,j,i) - ice%u_3D_a( k,j,i)) / (C%zeta(k+1) - C%zeta(k)))
        w4 = ((ice%dHs_dy_a( j,i) - 0.5_dp * (C%zeta(k+1) + C%zeta(k)) * ice%dHi_dy_a( j,i)) / MAX(0.1_dp, ice%Hi_a( j,i))) * &
             ((ice%v_3D_a( k+1,j,i) - ice%v_3D_a( k,j,i)) / (C%zeta(k+1) - C%zeta(k)))

        ice%w_3D_a( k,j,i) = ice%w_3D_a( k+1,j,i) - ice%Hi_a( j,i) * (w1 + w2 + w3 + w4) * (C%zeta(k+1) - C%zeta(k))
        
      END DO ! DO k = C%nZ - 1, 1, -1

    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calc_3D_vertical_velocities
  
  FUNCTION calc_staggered_margin( var0, var1, Hi0, Hi1) result(var_mid)
    ! Calculate a staggered point but taking upstream point at the margin 

    IMPLICIT NONE
    
    REAL(dp),                            INTENT(IN)    :: var0
    REAL(dp),                            INTENT(IN)    :: var1
    REAL(dp),                            INTENT(IN)    :: Hi0
    REAL(dp),                            INTENT(IN)    :: Hi1
    REAL(dp)                                           :: var_mid

    IF     (Hi0 >  0._dp .AND. Hi1 == 0._dp) THEN
      var_mid = var0 
    ELSEIF (Hi0 == 0._dp .AND. Hi1 >  0._dp) THEN
      var_mid = var1
    ELSE 
      var_mid = 0.5_dp * (var0 + var1) 
    END IF

  END FUNCTION calc_staggered_margin
    
  ! Some "administration" routines that help speed up and stabilise the solver
  SUBROUTINE initialise_DIVA_solve_masks( grid, ice)
    ! Approach copied from Yelmo/SICOPOLIS: for grid cells where no ice is present,
    ! or where the viscosity iteration has already converged, don't solve the DIVA.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables
    INTEGER                                            :: i,j
    
    IF (C%do_benchmark_experiment .AND. C%choice_benchmark_experiment == 'SSA_icestream') THEN
      ice%DIVA_solve_mask_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 1
      ice%DIVA_solve_mask_cy( :,grid%i1:              grid%i2 ) = 1
      CALL sync
      RETURN
    END IF
    
    IF (C%choice_ice_margin == 'BC') THEN
      ! Update the solve masks
    ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
      ! Just solve the equations everywhere
      ice%DIVA_solve_mask_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 1
      ice%DIVA_solve_mask_cy( :,grid%i1:              grid%i2 ) = 1
      RETURN
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in initialise_DIVA_solve_masks!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Initially no active ssa points
    ice%DIVA_solve_mask_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 0
    ice%DIVA_solve_mask_cy( :,grid%i1:              grid%i2 ) = 0

    ! x-direction
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny

      IF (ice%Hi_a( j,i) > 0._dp .or. ice%Hi_a( j,i+1) > 0._dp) THEN 
        ! Ice is present on ac-node
            
        IF (ice%f_grnd_cx( j,i) > 0._dp) THEN 
          ! Grounded ice or grounding line (ie, shelfy-stream)
          ice%DIVA_solve_mask_cx( j,i) = 1
        ELSE 
          ! Shelf ice 
          ice%DIVA_solve_mask_cx( j,i) = 2
        END IF 

        ! Deactivate if dragging is too high and away from grounding line
        IF ( ice%beta_cx( j,i) >= C%DIVA_beta_max .AND. ice%f_grnd_cx( j,i) == 1.0 ) ice%DIVA_solve_mask_cx( j,i) = 0 
            
      END IF

    END DO 
    END DO

    ! y-direction
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1

      IF (ice%Hi_a( j,i) > 0._dp .or. ice%Hi_a( j+1,i) > 0._dp) THEN 
        ! Ice is present on ac-node
            
        IF (ice%f_grnd_cy( j,i) > 0._dp) THEN 
          ! Grounded ice or grounding line (ie, shelfy-stream)
          ice%DIVA_solve_mask_cy( j,i) = 1
        ELSE 
          ! Shelf ice 
          ice%DIVA_solve_mask_cy( j,i) = 2
        END IF 

        ! Deactivate if dragging is too high and away from grounding line
        IF ( ice%beta_cy( j,i) >= C%DIVA_beta_max .AND. ice%f_grnd_cy( j,i) == 1.0 ) ice%DIVA_solve_mask_cy( j,i) = 0 
            
      END IF
         
    END DO 
    END DO
    CALL sync

    ! Final check on both masks to avoid isolated non-ssa points
    IF (par%master) THEN
      DO i = 2, grid%nx-1
      DO j = 2, grid%ny-1
  
        ! acx-nodes 
        IF (i < grid%nx-1) THEN
          IF (  ice%DIVA_solve_mask_cx( j,i) == 0 .AND. &
            ice%DIVA_solve_mask_cx( j,i+1) > 0 .AND. ice%DIVA_solve_mask_cx( j,i-1) > 0 .AND.  &
            ice%DIVA_solve_mask_cx( j+1,i) > 0 .AND. ice%DIVA_solve_mask_cx( j-1,i) > 0 ) THEN 
    
            IF (ice%f_grnd_cx( j,i) > 0._dp) THEN 
              ice%DIVA_solve_mask_cx( j,i) = 1
            ELSE 
              ice%DIVA_solve_mask_cx( j,i) = 2
            END IF 
    
          END IF 
        END IF ! IF (i < grid%nx-1) THEN
  
        ! acy-nodes 
        IF (j < grid%ny-1) THEN
          IF (  ice%DIVA_solve_mask_cy( j,i) == 0 .AND. &
            ice%DIVA_solve_mask_cy( j,i+1) > 0 .AND. ice%DIVA_solve_mask_cy( j,i-1) > 0 .AND.  &
            ice%DIVA_solve_mask_cy( j+1,i) > 0 .AND. ice%DIVA_solve_mask_cy( j-1,i) > 0 ) THEN 
    
            IF (ice%f_grnd_cy( j,i) > 0._dp) THEN
              ice%DIVA_solve_mask_cy( j,i) = 1
            ELSE
              ice%DIVA_solve_mask_cy( j,i) = 2
            END IF
    
          END IF
        END IF ! IF (j < grid%ny-1) THEN
  
      END DO 
      END DO
    END IF ! IF (par%master) THEN
    CALL sync
        
  END SUBROUTINE initialise_DIVA_solve_masks
  SUBROUTINE update_DIVA_solve_masks( grid, ice)
    ! Approach copied from Yelmo/SICOPOLIS: for grid cells where no ice is present,
    ! or where the viscosity iteration has already converged, don't solve the DIVA.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables
    INTEGER                                            :: i,j,ii,jj
    
    IF (C%do_benchmark_experiment .AND. C%choice_benchmark_experiment == 'SSA_icestream') THEN
      ice%DIVA_solve_mask_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 1
      ice%DIVA_solve_mask_cy( :,grid%i1:              grid%i2 ) = 1
      CALL sync
      RETURN
    END IF

    ! Initially set candidate 'converged' points to -2
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      IF (ice%DIVA_solve_mask_cx( j,i) == 1 .AND. ice%DIVA_err_cx( j,i) < C%DIVA_err_lim) ice%DIVA_solve_mask_cx( j,i) = -2
    END DO
    END DO
    CALL sync

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      IF (ice%DIVA_solve_mask_cy( j,i) == 1 .AND. ice%DIVA_err_cy( j,i) < C%DIVA_err_lim) ice%DIVA_solve_mask_cy( j,i) = -2
    END DO
    END DO
    CALL sync
    
    ! Fill in neighbours of points that are still ssa (mask=1) to keep things clean 
    IF (par%master) THEN
      ! cx
      DO i = 2, grid%nx - 2
      DO j = 2, grid%ny - 1
          IF (ice%DIVA_solve_mask_cx(j,i) == 1) THEN
            DO ii = i-1, i+1
            DO jj = j-1, j+1
              IF (ice%DIVA_solve_mask_cx( jj,ii) == -2) ice%DIVA_solve_mask_cx( jj,ii) = 1
            END DO
            END DO
          END IF 
      END DO
      END DO
      ! cy
      DO i = 2, grid%nx - 1
      DO j = 2, grid%ny - 2
          IF (ice%DIVA_solve_mask_cy(j,i) == 1) THEN
            DO ii = i-1, i+1
            DO jj = j-1, j+1
              IF (ice%DIVA_solve_mask_cy( jj,ii) == -2) ice%DIVA_solve_mask_cy( jj,ii) = 1
            END DO
            END DO
          END IF 
      END DO
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Finally, replace temporary -2 values with -1 to prescribe ssa vel here 
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      IF (ice%DIVA_solve_mask_cx( j,i) == -2) ice%DIVA_solve_mask_cx( j,i) = -1
    END DO
    END DO
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      IF (ice%DIVA_solve_mask_cy( j,i) == -2) ice%DIVA_solve_mask_cy( j,i) = -1
    END DO
    END DO
    CALL sync

  END SUBROUTINE update_DIVA_solve_masks
  SUBROUTINE set_sico_masks( grid, ice)
    ! mask = 0: Grounded ice 
    ! mask = 1: Ice-free land 
    ! mask = 2: Open ocean  
    ! mask = 3: Ice shelf 

    ! Note: this mask is defined on central Aa nodes 
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
        
    ! Local variables:
    INTEGER                                            :: i, j, i1, i2, j1, j2 
    LOGICAL                                            :: is_float 
    LOGICAL, PARAMETER                                 :: disable_grounded_fronts = .TRUE. 

    ! First determine general ice coverage mask 
    ! (land==1,ocean==2,floating_ice==3,grounded_ice==0)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
        
        ! Check IF this point would be floating
        is_float = ice%TAF_a( j,i) < 0._dp 

        IF (ice%Hi_a( j,i) == 0._dp) THEN 
            ! Ice-free point 

            IF (is_float) THEN 
                ! Open ocean 
                ice%DIVA_mask_combi( j,i) = 2
            ELSE 
                ! Ice-free land 
                ice%DIVA_mask_combi( j,i) = 1 
            END IF 

        ELSE 
            ! Ice-covered point

            IF (is_float) THEN 
                ! Ice shelf 
                ice%DIVA_mask_combi( j,i) = 3
            ELSE
                ! Grounded ice 
                ice%DIVA_mask_combi( j,i) = 0 
            END IF 
            
        END IF 

    END DO 
    END DO
    CALL sync
    
    !-------- Detection of the grounding line and the calving front --------

    ice%DIVA_isfront_inner( :,grid%i1:MIN(grid%ny-1,grid%i2)) = 0
    ice%DIVA_isfront_outer( :,grid%i1:              grid%i2 ) = 0
    CALL sync

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
     
      i1 = max(i-1,1)
      i2 = min(i+1,grid%nx)
      j1 = max( j-1,1)
      j2 = min( j+1,grid%ny)
      
      IF (ice%Hi_a( j,i) > 0._dp .AND. &
         (ice%Hi_a( j,i1) == 0._dp .OR. &
          ice%Hi_a( j,i2) == 0._dp .OR. &
          ice%Hi_a( j1,i) == 0._dp .OR. &
          ice%Hi_a( j2,i) == 0._dp)) THEN 

        ice%DIVA_isfront_inner( j,i) = 1

      END IF 

      IF (ice%Hi_a( j,i) == 0._dp .AND. &
         (ice%Hi_a( j,i1) > 0._dp .OR. &
          ice%Hi_a( j,i2) > 0._dp .OR. &
          ice%Hi_a( j1,i) > 0._dp .OR. &
          ice%Hi_a( j2,i) > 0._dp)) THEN 

        ice%DIVA_isfront_outer( j,i) = 1

      END IF 

      IF (disable_grounded_fronts) THEN 
        ! Disable detection of grounded fronts for now,
        ! because it is more stable this way...

        IF ( ice%DIVA_isfront_inner( j,i) == 1 .AND. ice%DIVA_mask_combi( j,i) == 0 ) ice%DIVA_isfront_inner( j,i) = 0

      END IF 

    END DO
    END DO
    CALL sync
        
  END SUBROUTINE set_sico_masks
  SUBROUTINE apply_velocity_limits( grid, u)
    ! Apply a velocity limit (for stability)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: u

    ! Local variables:
    INTEGER                                            :: i,j
    INTEGER                                            :: nx, ny
    
    nx = SIZE( u,2)
    ny = SIZE( u,1)
    
    DO i = grid%i1, MIN(nx-1,grid%i2)
    DO j = 1, ny
      IF     (u( j,i) >  C%DIVA_vel_max) THEN
        u( j,i) =  C%DIVA_vel_max
      ELSEIF (u( j,i) < -C%DIVA_vel_max) THEN
        u( j,i) = -C%DIVA_vel_max
      ELSEIF (ABS(u( j,i)) < C%DIVA_vel_min) THEN
        u( j,i) = 0._dp
      END IF
    END DO
    END DO
    CALL sync

  END SUBROUTINE apply_velocity_limits
  SUBROUTINE relax_DIVA_visc_iterations( grid, ice, u_cx, v_cy, rel)
    ! Relax velocity solution with results of the previous viscosity iteration 
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: u_cx
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: v_cy
    REAL(dp),                            INTENT(IN)    :: rel
        
    ! Local variables:
    INTEGER                                            :: i,j

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (i < grid%nx) u_cx( j,i) = rel * u_cx( j,i) + (1._dp - rel) * ice%u_cx_prev( j,i)
      IF (j < grid%ny) v_cy( j,i) = rel * v_cy( j,i) + (1._dp - rel) * ice%v_cy_prev( j,i)
    END DO
    END DO
    CALL sync

  END SUBROUTINE relax_DIVA_visc_iterations
  SUBROUTINE calculate_visc_iter_UV_resid( grid, ice, u_cx, v_cy, resid_UV)
    ! Check if the viscosity iteration has converged to a stable solution
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_cx
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_cy
    REAL(dp),                            INTENT(OUT)   :: resid_UV
    
    ! Local variables:
    INTEGER                                            :: ierr
    INTEGER                                            :: i,j,ncx,ncy
    REAL(dp)                                           :: res1, res2
    REAL(dp), PARAMETER                                :: DIVA_vel_tolerance = 1e-6   ! [m/a] only consider points with velocity above this tolerance limit

    ! Calculate the L2 norm based on velocity solution between previous
    ! and current viscosity iteration (as in Yelmo/SICOPOLIS)
    
    ncx  = 0
    ncy  = 0
    res1 = 0._dp
    res2 = 0._dp
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      IF (i < grid%nx) THEN
        IF (ABS(u_cx( j,i)) > DIVA_vel_tolerance .AND. ice%DIVA_solve_mask_cx( j,i) == 1) THEN
          ncx = ncx + 1
          res1 = res1 + (u_cx( j,i) - ice%u_cx_prev( j,i))**2._dp
          res2 = res2 + (u_cx( j,i) + ice%u_cx_prev( j,i))**2._dp
        END IF
      END IF
      
      IF (j < grid%ny) THEN
        IF (ABS(v_cy( j,i)) > DIVA_vel_tolerance .AND. ice%DIVA_solve_mask_cy( j,i) == 1) THEN
          ncy = ncy + 1
          res1 = res1 + (v_cy( j,i) - ice%v_cy_prev( j,i))**2._dp
          res2 = res2 + (v_cy( j,i) + ice%v_cy_prev( j,i))**2._dp
        END IF
      END IF
      
    END DO
    END DO
    
    ! Combine results from all processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, ncx,  1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, ncy,  1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    IF (ncx + ncy > 0) THEN
      res1 = SQRT( res1)
      res2 = SQRT( res2 )
      res2 = MAX(  res2,1E-8_dp)
      resid_UV = 2._dp * res1 / res2 
    ELSE 
      ! No points available for comparison, set residual equal to zero 
      resid_UV = 0._dp
    END IF

  END SUBROUTINE calculate_visc_iter_UV_resid
  SUBROUTINE estimate_visc_iter_UV_errors( grid, ice, u_cx, v_cy)
    ! Estimate the error in the velocity fields so we can
    ! update the SICOPOLIS-style DIVA solving masks (for better efficiency)
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_cx
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_cy

    ! Local variables
    INTEGER                                            :: i,j
    REAL(dp), PARAMETER                                :: DIVA_vel_tolerance = 1e-6   ! [m/a] only consider points with velocity above this tolerance limit
    REAL(dp), PARAMETER                                :: tol = 1E-5_dp

    ! Error in x-direction
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      IF (i < grid%nx) THEN
        IF (ABS(ice%u_vav_cx( j,i)) > DIVA_vel_tolerance) THEN
          ice%DIVA_err_cx( j,i) = 2._dp * ABS(u_cx( j,i) - ice%u_cx_prev( j,i)) / &
                                          ABS(u_cx( j,i) + ice%u_cx_prev( j,i) + tol)
        ELSE
          ice%DIVA_err_cx( j,i) = 0._dp
        END IF
      END IF
      
      IF (j < grid%ny) THEN
        IF (ABS(ice%v_vav_cy( j,i)) > DIVA_vel_tolerance) THEN
          ice%DIVA_err_cy( j,i) = 2._dp * ABS(v_cy( j,i) - ice%v_cy_prev( j,i)) / &
                                          ABS(v_cy( j,i) + ice%v_cy_prev( j,i) + tol)
        ELSE
          ice%DIVA_err_cy( j,i) = 0._dp
        END IF
      END IF
      
    END DO
    END DO
    CALL sync

  END SUBROUTINE estimate_visc_iter_UV_errors
  
  ! Routines for solving the "linearised" DIVA by writing and solving it as a matrix equation
  ! (this is where the discretisation happens)
  SUBROUTINE solve_DIVA_stag_linearised( grid, ice, u, v)
    ! Solve the "linearised" version of the DIVA (i.e. assuming viscosity and basal stress are
    ! constant rather than functions of velocity).
    ! Adapted from Yelmo, which was in turn copied from SICOPOLIS.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: u
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: v

    ! Local variables
    INTEGER                                            :: i, j, k, n
    REAL(dp)                                           :: t0

    ! Fill in the old SICOPOLIS masks (used to determine
    ! where to apply ice margin boundary conditions)
    CALL set_sico_masks( grid, ice)

    ! Calculate and list (in compressed sparse row [CSR] format) all the matrix coefficients
    ! ======================================================================================
    
    IF (par%master) THEN
  
      ice%DIVA_m%A_ptr   = 0
      ice%DIVA_m%A_index = 0
      ice%DIVA_m%A_val   = 0._dp
      ice%DIVA_m%b       = 0._dp
      ice%DIVA_m%x       = 0._dp
  
      k = 0
      ice%DIVA_m%A_ptr(1) = 1
      
      IF (C%choice_ice_margin == 'BC') THEN
        ! Apply velocity boundary conditions at the ice margin, don't solve
        ! the equations for non-ice pixels
  
        DO n = 1, ice%DIVA_m%m
          IF (ice%DIVA_m_n2ij_uv( n,1) > 0) THEN
            ! Equation 1
            CALL list_DIVA_matrix_coefficients_eq_1_with_BC( grid, ice, u, k, n)
          ELSE
            ! Equation 2
            CALL list_DIVA_matrix_coefficients_eq_2_with_BC( grid, ice, v, k, n)
          END IF
          ice%DIVA_m%A_ptr( n+1) = k+1
        END DO ! DO n = 1, ice%DIVA_m_neq
        
      ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
        ! Use the "infinite slab" approach from PISM and ANICE
  
        DO n = 1, ice%DIVA_m%m
          IF (ice%DIVA_m_n2ij_uv( n,1) > 0) THEN
            ! Equation 1
            CALL list_DIVA_matrix_coefficients_eq_1_infinite_slab( grid, ice, u, k, n)
          ELSE
            ! Equation 2
            CALL list_DIVA_matrix_coefficients_eq_2_infinite_slab( grid, ice, v, k, n)
          END IF
          ice%DIVA_m%A_ptr( n+1) = k+1
        END DO ! DO n = 1, ice%DIVA_m_neq
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_ice_margin "', TRIM(C%choice_ice_margin), '" not implemented in solve_DIVA_stag_linearised!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Solve the matrix equation using SOR
    ! ===================================
    
    t0 = MPI_WTIME()
    CALL solve_matrix_equation_CSR( ice%DIVA_m, C%DIVA_choice_matrix_solver, &
      SOR_nit = C%DIVA_SOR_nit, SOR_tol = C%DIVA_SOR_tol, SOR_omega = C%DIVA_SOR_omega, &
      PETSc_rtol = C%DIVA_PETSc_rtol, PETSc_abstol = C%DIVA_PETSc_abstol)
    !IF (par%master) WRITE(0,*) '   ', TRIM(C%DIVA_choice_matrix_solver), ' solved linearised ', TRIM(C%choice_ice_dynamics), ' in ', MPI_WTIME()-t0, ' seconds.'

    ! Map the solution back from vector format to the model grids
    ! ===========================================================
    
    ! u
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      n = ice%DIVA_m_ij2n_u( j,i)
      u( j,i) = ice%DIVA_m%x( n)
    END DO
    END DO
    
    ! v
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      n = ice%DIVA_m_ij2n_v( j,i)
      v( j,i) = ice%DIVA_m%x( n)
    END DO
    END DO
    CALL sync

  END SUBROUTINE solve_DIVA_stag_linearised
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_with_BC(          grid, ice, u, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: u
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n
    
    ! Local variables:
    INTEGER                                            :: i,j

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,1)
    j = ice%DIVA_m_n2ij_uv( n,2)

    ! Exception for boundary conditions for the domain boundary
    IF (i == 1 .OR. i == grid%nx-1 .OR. j == 1 .OR. j == grid%ny) THEN
      CALL list_DIVA_matrix_coefficients_eq_1_boundary( grid, ice, u, k, n)
      RETURN
    END IF
    
    ! Exception for the grounding line
    IF (C%use_analytical_GL_flux .AND. ice%u_GL_cx( j,i) /= 0._dp) THEN
      k = k+1
      ice%DIVA_m%A_index( k) = n
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = u( j,i) + 0.3_dp * (ice%u_GL_cx( j,i) - u( j,i)) ! Include a little "relaxation" for better numerical stability
      ice%DIVA_m%x( n) = u( j,i) + 0.3_dp * (ice%u_GL_cx( j,i) - u( j,i))
      RETURN
    END IF

    ! Exception for locations where the viscosity iteration has already converged
    IF (ice%DIVA_solve_mask_cx( j,i) == -1) THEN
      k = k+1
      ice%DIVA_m%A_index( k) = n
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = u( j,i)
      ice%DIVA_m%x( n) = u( j,i)
      RETURN
    END IF

    ! Exception for grid cells where neither neighbour has ice (set velocity to zero)
    IF (ice%DIVA_solve_mask_cx( j,i) <= 0) THEN
      k = k+1
      ice%DIVA_m%A_index( k) = n
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = 0._dp
      ice%DIVA_m%x( n) = 0._dp
      RETURN
    END IF
        
    ! Exception for floating ice margin bordering ice-free land; set velocity to zero
    IF (((ice%DIVA_mask_combi( j,i) == 3) .AND. (ice%DIVA_mask_combi( j,i+1) == 1)) .OR. &
        ((ice%DIVA_mask_combi( j,i) == 1) .AND. (ice%DIVA_mask_combi( j,i+1) == 3))) THEN
      k = k+1
      ice%DIVA_m%A_index( k) = n
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = 0._dp
      ice%DIVA_m%x( n) = 0._dp
      RETURN
    END IF
     
    ! Exception for the ice margin (stress boundary condition)
    IF ((ice%DIVA_isfront_inner( j,i) == 1 .AND. ice%DIVA_isfront_outer( j,i+1) == 1) .OR. &
        (ice%DIVA_isfront_outer( j,i) == 1 .AND. ice%DIVA_isfront_inner( j,i+1) == 1)) THEN
      ! One neighbour is ice-covered and the other is ice-free (calving front, grounded ice front)
      CALL list_DIVA_matrix_coefficients_eq_1_margin( grid, ice, u, k, n)
      RETURN
    END IF

    ! No more exceptions; solve the complete, free DIVA
    CALL list_DIVA_matrix_coefficients_eq_1_free( grid, ice, u, k, n)
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_with_BC
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_with_BC(          grid, ice, v, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: v
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n
    
    ! Local variables:
    INTEGER                                            :: i,j

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,3)
    j = ice%DIVA_m_n2ij_uv( n,4)
    
    ! Exception for boundary conditions for the domain boundary
    IF (i == 1 .OR. i == grid%nx .OR. j == 1 .OR. j == grid%ny-1) THEN
      CALL list_DIVA_matrix_coefficients_eq_2_boundary( grid, ice, v, k, n)
      RETURN
    END IF
    
    ! Exception for the grounding line
    IF (C%use_analytical_GL_flux .AND. ice%v_GL_cy( j,i) /= 0._dp) THEN
      k = k+1
      ice%DIVA_m%A_index( k) = n
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = v( j,i) + 0.3_dp * (ice%v_GL_cy( j,i) - v( j,i)) ! Include a little "relaxation" for better numerical stability
      ice%DIVA_m%x( n) = v( j,i) + 0.3_dp * (ice%v_GL_cy( j,i) - v( j,i))
      RETURN
    END IF
    
    ! Exception for locations where the viscosity iteration has already converged
    IF (ice%DIVA_solve_mask_cy( j,i) == -1) THEN
      k = k+1
      ice%DIVA_m%A_index( k)  = n
      ice%DIVA_m%A_val(   k)  = 1._dp
      ice%DIVA_m%b( n) = v( j,i)
      ice%DIVA_m%x( n) = v( j,i)
      RETURN
    END IF

    ! Exception for grid cells where neither neighbour has ice (set velocity to zero)
    IF (ice%DIVA_solve_mask_cy(j,i) <= 0) THEN
      k = k+1
      ice%DIVA_m%A_index( k) = n
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = 0._dp
      ice%DIVA_m%x( n) = 0._dp
      RETURN
    END IF

    ! Exception for floating ice margin bordering ice-free land; set velocity to zero
    IF (((ice%DIVA_mask_combi( j,i) == 3) .AND. (ice%DIVA_mask_combi( j+1,i) == 1)) .OR. &
        ((ice%DIVA_mask_combi( j,i) == 1) .AND. (ice%DIVA_mask_combi( j+1,i) == 3))) THEN
      k = k+1
      ice%DIVA_m%A_index( k) = n
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = 0._dp
      ice%DIVA_m%x( n) = 0._dp
      RETURN
    END IF
    
    ! Exception for the ice margin (stress boundary condition)
    IF ((ice%DIVA_isfront_inner( j,i) == 1 .AND. ice%DIVA_isfront_outer( j+1,i) == 1) .OR. &
        (ice%DIVA_isfront_outer( j,i) == 1 .AND. ice%DIVA_isfront_inner( j+1,i) == 1)) THEN
      ! One neighbour is ice-covered and the other is ice-free (calving front, grounded ice front)
      CALL list_DIVA_matrix_coefficients_eq_2_margin( grid, ice, v, k, n)
      RETURN
    END IF

    ! No more exceptions; solve the complete, free DIVA
    CALL list_DIVA_matrix_coefficients_eq_2_free( grid, ice, v, k, n)
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_with_BC
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_infinite_slab(    grid, ice, u, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: u
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n
    
    ! Local variables:
    INTEGER                                            :: i,j

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,1)
    j = ice%DIVA_m_n2ij_uv( n,2)

    ! Exception for boundary conditions for the domain boundary
    IF (i == 1 .OR. i == grid%nx-1 .OR. j == 1 .OR. j == grid%ny) THEN
      CALL list_DIVA_matrix_coefficients_eq_1_boundary( grid, ice, u, k, n)
      RETURN
    END IF
    
    ! Exception for the grounding line
    IF (C%use_analytical_GL_flux .AND. ice%u_GL_cx( j,i) /= 0._dp) THEN
      k = k+1
      ice%DIVA_m%A_index( k) = n
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = u( j,i) + 0.3_dp * (ice%u_GL_cx( j,i) - u( j,i)) ! Include a little "relaxation" for better numerical stability
      ice%DIVA_m%x( n) = u( j,i) + 0.3_dp * (ice%u_GL_cx( j,i) - u( j,i))
      RETURN
    END IF

    ! Exception for locations where the viscosity iteration has already converged
    IF (ice%DIVA_solve_mask_cx( j,i) == -1) THEN
      k = k+1
      ice%DIVA_m%A_index( k) = n
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = u( j,i)
      ice%DIVA_m%x( n) = u( j,i)
      RETURN
    END IF

    ! No more exceptions; solve the complete, free DIVA
    CALL list_DIVA_matrix_coefficients_eq_1_free( grid, ice, u, k, n)
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_infinite_slab
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_infinite_slab(    grid, ice, v, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: v
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n
    
    ! Local variables:
    INTEGER                                            :: i,j

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,3)
    j = ice%DIVA_m_n2ij_uv( n,4)
    
    ! Exception for boundary conditions for the domain boundary
    IF (i == 1 .OR. i == grid%nx .OR. j == 1 .OR. j == grid%ny-1) THEN
      CALL list_DIVA_matrix_coefficients_eq_2_boundary( grid, ice, v, k, n)
      RETURN
    END IF
    
    ! Exception for the grounding line
    IF (C%use_analytical_GL_flux .AND. ice%v_GL_cy( j,i) /= 0._dp) THEN
      k = k+1
      ice%DIVA_m%A_index( k) = n
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = v( j,i) + 0.3_dp * (ice%v_GL_cy( j,i) - v( j,i)) ! Include a little "relaxation" for better numerical stability
      ice%DIVA_m%x( n) = v( j,i) + 0.3_dp * (ice%v_GL_cy( j,i) - v( j,i))
      RETURN
    END IF
    
    ! Exception for locations where the viscosity iteration has already converged
    IF (ice%DIVA_solve_mask_cy( j,i) == -1) THEN
      k = k+1
      ice%DIVA_m%A_index( k)  = n
      ice%DIVA_m%A_val(   k)  = 1._dp
      ice%DIVA_m%b( n) = v( j,i)
      ice%DIVA_m%x( n) = v( j,i)
      RETURN
    END IF

    ! No more exceptions; solve the complete, free DIVA
    CALL list_DIVA_matrix_coefficients_eq_2_free( grid, ice, v, k, n)
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_infinite_slab
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_boundary( grid, ice, u, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: u
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n
    
    ! Local variables:
    INTEGER                                            :: i,j

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,1)
    j = ice%DIVA_m_n2ij_uv( n,2)

    ! Exception for boundary conditions for the domain boundary
    IF (i == 1) THEN ! Western domain boundary
      IF     (C%DIVA_boundary_BC_u_west == 'infinite') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = u( j,2)
        ice%DIVA_m%x( n) = u( j,2)
      ELSEIF (C%DIVA_boundary_BC_u_west == 'zero') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = 0._dp
        ice%DIVA_m%x( n) = 0._dp
      ELSEIF (C%DIVA_boundary_BC_u_west == 'periodic') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = u( j,grid%nx-2)
        ice%DIVA_m%x( n) = u( j,grid%nx-2)
      ELSE
        WRITE(0,*) ' solve_DIVA_stag_linearised - ERROR: unknown DIVA_boundary_BC_u_west "', C%DIVA_boundary_BC_u_west, '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF 
    ELSEIF (i == grid%nx-1) THEN ! Eastern domain boundary
      IF     (C%DIVA_boundary_BC_u_east == 'infinite') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = u( j,grid%nx-2)
        ice%DIVA_m%x( n) = u( j,grid%nx-2)
      ELSEIF (C%DIVA_boundary_BC_u_east == 'zero') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = 0._dp
        ice%DIVA_m%x( n) = 0._dp
      ELSEIF (C%DIVA_boundary_BC_u_east == 'periodic') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = u( j,2)
        ice%DIVA_m%x( n) = u( j,2)
      ELSE
        WRITE(0,*) ' solve_DIVA_stag_linearised - ERROR: unknown DIVA_boundary_BC_u_east "', C%DIVA_boundary_BC_u_east, '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSEIF (j == 1) THEN ! Southern domain boundary
      IF     (C%DIVA_boundary_BC_u_south == 'infinite') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = u( 2,i)
        ice%DIVA_m%x( n) = u( 2,i)
      ELSEIF (C%DIVA_boundary_BC_u_south == 'zero') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = 0._dp
        ice%DIVA_m%x( n) = 0._dp
      ELSEIF (C%DIVA_boundary_BC_u_south == 'periodic') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = u( grid%ny-1,i)
        ice%DIVA_m%x( n) = u( grid%ny-1,i)
      ELSE
        WRITE(0,*) ' solve_DIVA_stag_linearised - ERROR: unknown DIVA_boundary_BC_u_south "', C%DIVA_boundary_BC_u_south, '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF 
    ELSEIF (j == grid%ny) THEN ! Northern domain boundary
      IF     (C%DIVA_boundary_BC_u_north == 'infinite') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = u( grid%ny-1,i)
        ice%DIVA_m%x( n) = u( grid%ny-1,i)
      ELSEIF (C%DIVA_boundary_BC_u_north == 'zero') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = 0._dp
        ice%DIVA_m%x( n) = 0._dp
      ELSEIF (C%DIVA_boundary_BC_u_north == 'periodic') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = u( 2,i)
        ice%DIVA_m%x( n) = u( 2,i)
      ELSE
        WRITE(0,*) ' solve_DIVA_stag_linearised - ERROR: unknown DIVA_boundary_BC_u_north "', C%DIVA_boundary_BC_u_north, '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSE
      WRITE(0,*) ' solve_DIVA_stag_linearised - ERROR: grid cell [',i,',',j,'] doesnt lie on the domain boundary!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_boundary
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_boundary( grid, ice, v, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: v
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n
    
    ! Local variables:
    INTEGER                                            :: i,j

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,3)
    j = ice%DIVA_m_n2ij_uv( n,4)
    
    ! Exception for boundary conditions for the domain boundary
    IF (i == 1) THEN ! Western domain boundary
      IF     (C%DIVA_boundary_BC_v_west == 'infinite') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = v( j,2)
        ice%DIVA_m%x( n) = v( j,2)
      ELSEIF (C%DIVA_boundary_BC_v_west == 'zero') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = 0._dp
        ice%DIVA_m%x( n) = 0._dp
      ELSEIF (C%DIVA_boundary_BC_v_west == 'periodic') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = v( j,grid%nx-1)
        ice%DIVA_m%x( n) = v( j,grid%nx-1)
      ELSE
        WRITE(0,*) ' solve_DIVA_stag_linearised - ERROR: unknown DIVA_boundary_BC_v_west "', C%DIVA_boundary_BC_v_west, '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF 
    ELSEIF (i == grid%nx) THEN ! Eastern domain boundary
      IF     (C%DIVA_boundary_BC_v_east == 'infinite') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = v( j,grid%nx-1)
        ice%DIVA_m%x( n) = v( j,grid%nx-1)
      ELSEIF (C%DIVA_boundary_BC_v_east == 'zero') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = 0._dp
        ice%DIVA_m%x( n) = 0._dp
      ELSEIF (C%DIVA_boundary_BC_v_east == 'periodic') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = v( j,2)
        ice%DIVA_m%x( n) = v( j,2)
      ELSE
        WRITE(0,*) ' solve_DIVA_stag_linearised - ERROR: unknown DIVA_boundary_BC_v_east "', C%DIVA_boundary_BC_v_east, '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSEIF (j == 1) THEN ! Southern domain boundary
      IF     (C%DIVA_boundary_BC_v_south == 'infinite') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = v( 2,i)
        ice%DIVA_m%x( n) = v( 2,i)
      ELSEIF (C%DIVA_boundary_BC_v_south == 'zero') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = 0._dp
        ice%DIVA_m%x( n) = 0._dp
      ELSEIF (C%DIVA_boundary_BC_v_south == 'periodic') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = v( grid%ny-2,i)
        ice%DIVA_m%x( n) = v( grid%ny-2,i)
      ELSE
        WRITE(0,*) ' solve_DIVA_stag_linearised - ERROR: unknown DIVA_boundary_BC_v_south "', C%DIVA_boundary_BC_v_south, '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF 
    ELSEIF (j == grid%ny-1) THEN ! Northern domain boundary
      IF     (C%DIVA_boundary_BC_v_north == 'infinite') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = v( grid%ny-2,i)
        ice%DIVA_m%x( n) = v( grid%ny-2,i)
      ELSEIF (C%DIVA_boundary_BC_v_north == 'zero') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = 0._dp
        ice%DIVA_m%x( n) = 0._dp
      ELSEIF (C%DIVA_boundary_BC_v_north == 'periodic') THEN
        k = k+1
        ice%DIVA_m%A_index( k) = n
        ice%DIVA_m%A_val(   k) = 1._dp
        ice%DIVA_m%b( n) = v( 2,i)
        ice%DIVA_m%x( n) = v( 2,i)
      ELSE
        WRITE(0,*) ' solve_DIVA_stag_linearised - ERROR: unknown DIVA_boundary_BC_v_north "', C%DIVA_boundary_BC_v_north, '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSE
      WRITE(0,*) ' solve_DIVA_stag_linearised - ERROR: grid cell [',i,',',j,'] doesnt lie on the domain boundary!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_boundary
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_margin(   grid, ice, u, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: u
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n
    
    ! Local variables:
    INTEGER                                            :: i,j,i1
    REAL(dp)                                           :: inv_dx
    REAL(dp)                                           :: H_ocn, H_ice

    ! Abbreviations of common factors in the equations
    inv_dx  = 1._dp /  grid%dx

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,1)
    j = ice%DIVA_m_n2ij_uv( n,2)

    ! Get the index of the neighbouring A-grid cell that has ice
    IF (ice%DIVA_isfront_inner( j,i) == 1) THEN
        i1 = i
    ELSE
        i1 = i+1
    END IF

    ! Exception for a single-pixel "ice ridge" (in the x-direction), set velocity to zero
    IF (ice%DIVA_isfront_outer( j,i1-1) == 1 .AND. ice%DIVA_isfront_outer( j,i1+1) == 1) THEN
    
      k = k+1
      ice%DIVA_m%A_index( k) = n
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = 0._dp
      ice%DIVA_m%x( n) = 0._dp
  
    ELSE

      k = k+1
      ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j,i1-1)
      ice%DIVA_m%A_val(   k) = -4._dp * inv_dx * ice%N_a( j,i1)
      
      k = k+1
      ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j,i1)
      ice%DIVA_m%A_val(   k) =  4._dp * inv_dx * ice%N_a( j,i1)

      k = k+1
      ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j-1,i1)
      ice%DIVA_m%A_val(   k) = -2._dp * inv_dx * ice%N_a( j,i1)

      k = k+1
      ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j,i1)
      ice%DIVA_m%A_val(   k) =  2._dp * inv_dx * ice%N_a( j,i1)

      ! Generalized solution for all ice fronts (floating and grounded)
      H_ice = ice%Hi_a( j,i1)
      IF (ice%SL_a( j,i1) > ice%Hb_a( j,i1)) THEN ! Bed below sea level 
        H_ocn = MIN(ice_density / seawater_density * ice%Hi_a( j,i1), & ! Flotation depth 
                    ice%SL_a( j,i1) - ice%Hb_a( j,i1))                  ! Grounded depth 
      ELSE
        ! Bed above sea level 
        H_ocn = 0._dp
      END IF

      ice%DIVA_m%b( n) = 0.5_dp * ice_density      * grav * H_ice * H_ice &
                       - 0.5_dp * seawater_density * grav * H_ocn * H_ocn  
      ice%DIVA_m%x( n) = u( j,i)

    END IF
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_margin
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_margin(   grid, ice, v, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: v
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n
    
    ! Local variables:
    INTEGER                                            :: i,j,j1
    REAL(dp)                                           :: inv_dx
    REAL(dp)                                           :: H_ocn, H_ice

    ! Abbreviations of common factors in the equations
    inv_dx  = 1._dp /  grid%dx

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,3)
    j = ice%DIVA_m_n2ij_uv( n,4)

    ! Get the index of the neighbouring A-grid cell that has ice
    IF (ice%DIVA_isfront_inner( j,i) == 1) THEN
        j1 = j
    ELSE
        j1 = j+1
    END IF

    ! Exception for a single-pixel "ice ridge" (in the x-direction), set velocity to zero
    IF (ice%DIVA_isfront_outer( j1-1,i) == 1 .AND. ice%DIVA_isfront_outer( j1+1,i) == 1) THEN
    
      k = k+1
      ice%DIVA_m%A_index( k) = n
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = 0._dp
      ice%DIVA_m%x( n) = 0._dp
  
    ELSE

      k = k+1
      ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j1-1,i)
      ice%DIVA_m%A_val(   k) = -4._dp * inv_dx * ice%N_a( j1,i)
      
      k = k+1
      ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j1,i)
      ice%DIVA_m%A_val(   k) =  4._dp * inv_dx * ice%N_a( j1,i)

      k = k+1
      ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j1,i-1)
      ice%DIVA_m%A_val(   k) = -2._dp * inv_dx * ice%N_a( j1,i)

      k = k+1
      ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j1,i)
      ice%DIVA_m%A_val(   k) =  2._dp * inv_dx * ice%N_a( j1,i)

      ! Generalized solution for all ice fronts (floating and grounded)
      H_ice = ice%Hi_a( j1,i)
      IF (ice%SL_a( j1,i) > ice%Hb_a( j1,i)) THEN ! Bed below sea level 
        H_ocn = MIN(ice_density / seawater_density * ice%Hi_a( j1,i), & ! Flotation depth 
                    ice%SL_a( j1,i) - ice%Hb_a( j1,i))                  ! Grounded depth 
      ELSE
        ! Bed above sea level 
        H_ocn = 0._dp
      END IF

      ice%DIVA_m%b( n) = 0.5_dp * ice_density      * grav * H_ice * H_ice &
                       - 0.5_dp * seawater_density * grav * H_ocn * H_ocn  
      ice%DIVA_m%x( n) = v( j,i)

    END IF
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_margin
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_free(     grid, ice, u, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: u
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n
    
    ! Local variables:
    INTEGER                                            :: i,j

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,1)
    j = ice%DIVA_m_n2ij_uv( n,2)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j,i+1)
    ice%DIVA_m%A_val(   k) = 4._dp * ice%N_a( j,i+1)
    
    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j,i-1)
    ice%DIVA_m%A_val(   k) = 4._dp * ice%N_a( j,i)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j+1,i)
    ice%DIVA_m%A_val(   k) = ice%N_b( j,i)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j-1,i)
    ice%DIVA_m%A_val(   k) = ice%N_b( j-1,i)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j,i)
    ice%DIVA_m%A_val(   k) = -4._dp * (ice%N_a( j,i+1) + ice%N_a( j  ,i)) + &
                             -1._dp * (ice%N_b( j,i  ) + ice%N_b( j-1,i)) - &
                             ice%beta_eff_cx( j,i) * grid%dx**2

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j,i)
    ice%DIVA_m%A_val(   k) = -2._dp * ice%N_a( j,i  ) - ice%N_b( j  ,i)

    k  = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j,i+1)
    ice%DIVA_m%A_val(   k) =  2._dp * ice%N_a( j,i+1) + ice%N_b( j  ,i)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j-1,i)
    ice%DIVA_m%A_val(   k) =  2._dp * ice%N_a( j,i  ) + ice%N_b( j-1,i)

    k  = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j-1,i+1)
    ice%DIVA_m%A_val(   k) = -2._dp * ice%N_a( j,i+1) - ice%N_b( j-1,i)

    ice%DIVA_m%b( n) = ice%taudx_cx( j,i) * grid%dx**2
    ice%DIVA_m%x( n) = u( j,i)
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_free
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_free(     grid, ice, v, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: v
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: inv_dx2

    ! Abbreviations of common factors in the equations
    inv_dx2 = 1._dp / (grid%dx * grid%dx)

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,3)
    j = ice%DIVA_m_n2ij_uv( n,4)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j+1,i)
    ice%DIVA_m%A_val(   k) = 4._dp * ice%N_a( j+1,i)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j-1,i)
    ice%DIVA_m%A_val(   k) = 4._dp * ice%N_a( j,i)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j,i+1)
    ice%DIVA_m%A_val(   k) = ice%N_b( j,i)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j,i-1)
    ice%DIVA_m%A_val(   k) = ice%N_b( j,i-1)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j,i)
    ice%DIVA_m%A_val(   k) = -4._dp * (ice%N_a( j+1,i) + ice%N_a( j,i  )) + &
                             -1._dp * (ice%N_b( j  ,i) + ice%N_b( j,i-1)) - &
                             ice%beta_eff_cy( j,i) * grid%dx**2

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j+1,i-1)
    ice%DIVA_m%A_val(   k) = -2._dp * ice%N_a( j+1,i) - ice%N_b( j,i-1)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j+1,i)
    ice%DIVA_m%A_val(   k) =  2._dp * ice%N_a( j+1,i) + ice%N_b( j,i)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j,i-1)
    ice%DIVA_m%A_val(   k) =  2._dp * ice%N_a( j,i) + ice%N_b( j,i-1)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j,i)
    ice%DIVA_m%A_val(   k) = -2._dp * ice%N_a( j,i) - ice%N_b( j,i)

    ice%DIVA_m%b( n) = ice%taudy_cy( j,i) * grid%dx**2
    ice%DIVA_m%x( n) = v( j,i)
    
  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_free
  SUBROUTINE initialise_DIVA_matrix_tables( grid, ice)
    ! Allocate and initialise the translation tables and sparse matrix lists used
    ! in the matrix-based SSA/DIVA solver.
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j,neq,nnz_per_row_max,n
    
    neq             = (grid%nx - 1) * grid%ny + grid%nx * (grid%ny - 1)
    nnz_per_row_max = 9
    
    CALL initialise_matrix_equation_CSR( ice%DIVA_m, neq, neq, nnz_per_row_max)
    
    CALL allocate_shared_int_2D( grid%ny  , grid%nx-1, ice%DIVA_m_ij2n_u  , ice%wDIVA_m_ij2n_u )
    CALL allocate_shared_int_2D( grid%ny-1, grid%nx  , ice%DIVA_m_ij2n_v  , ice%wDIVA_m_ij2n_v )
    CALL allocate_shared_int_2D( neq,       4        , ice%DIVA_m_n2ij_uv , ice%wDIVA_m_n2ij_uv)

    ! Alternate equations 1 and 2 in the matrix rows for better stability
    n = 0
    DO i = 1, grid%nx
    DO j = 1, grid%ny
      IF (i < grid%nx) THEN
        n = n + 1
        ice%DIVA_m_ij2n_u( j,i) = n
        ice%DIVA_m_n2ij_uv( n,1:2) = [i,j]
      END IF
      IF (j < grid%ny) THEN
        n = n+1
        ice%DIVA_m_ij2n_v( j,i) = n
        ice%DIVA_m_n2ij_uv( n,3:4) = [i,j]
      END IF
    END DO
    END DO
    
  END SUBROUTINE initialise_DIVA_matrix_tables
  
  ! Map velocity to the A-grid for writing to output (diagnostic only)
  SUBROUTINE map_velocity_from_cx_to_a_2D( grid, ice, u_cx, u_a)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_cx
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: u_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (i == 1) THEN
        IF (ice%mask_ice_a( j,i) == 1) THEN
          u_a( j,i) = u_cx( j,i)
        ELSE
          u_a( j,i) = 0._dp
        END IF
      ELSEIF (i == grid%nx) THEN
        IF (ice%mask_ice_a( j,i) == 1) THEN
          u_a( j,i) = u_cx( j,i-1)
        ELSE
          u_a( j,i) = 0._dp
        END IF
      ELSE
        IF     (ice%mask_ice_a( j,i-1) == 1 .AND. ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a( j,i+1) == 1) THEN
          u_a( j,i) = 0.5_dp * (u_cx( j,i-1) + u_cx( j,i))
        ELSEIF (ice%mask_ice_a( j,i-1) == 0 .AND. ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a( j,i+1) == 1) THEN
          u_a( j,i) = u_cx( j,i)
        ELSEIF (ice%mask_ice_a( j,i-1) == 1 .AND. ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a( j,i+1) == 0) THEN
          u_a( j,i) = u_cx( j,i-1)
        ELSE
          u_a( j,i) = 0._dp
        END IF
      END IF
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_velocity_from_cx_to_a_2D
  SUBROUTINE map_velocity_from_cx_to_a_3D( grid, ice, u_cx, u_a)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: u_cx
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: u_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (i == 1) THEN
        IF (ice%mask_ice_a( j,i) == 1) THEN
          u_a( :,j,i) = u_cx( :,j,i)
        ELSE
          u_a( :,j,i) = 0._dp
        END IF
      ELSEIF (i == grid%nx) THEN
        IF (ice%mask_ice_a( j,i) == 1) THEN
          u_a( :,j,i) = u_cx( :,j,i-1)
        ELSE
          u_a( :,j,i) = 0._dp
        END IF
      ELSE
        IF     (ice%mask_ice_a( j,i-1) == 1 .AND. ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a( j,i+1) == 1) THEN
          u_a( :,j,i) = 0.5_dp * (u_cx( :,j,i-1) + u_cx( :,j,i))
        ELSEIF (ice%mask_ice_a( j,i-1) == 0 .AND. ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a( j,i+1) == 1) THEN
          u_a( :,j,i) = u_cx( :,j,i)
        ELSEIF (ice%mask_ice_a( j,i-1) == 1 .AND. ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a( j,i+1) == 0) THEN
          u_a( :,j,i) = u_cx( :,j,i-1)
        ELSE
          u_a( :,j,i) = 0._dp
        END IF
      END IF
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_velocity_from_cx_to_a_3D
  SUBROUTINE map_velocity_from_cy_to_a_2D( grid, ice, v_cy, v_a)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_cy
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: v_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (j == 1) THEN
        IF (ice%mask_ice_a( j,i) == 1) THEN
          v_a( j,i) = v_cy( j,i)
        ELSE
          v_a( j,i) = 0._dp
        END IF
      ELSEIF (j == grid%ny) THEN
        IF (ice%mask_ice_a( j,i) == 1) THEN
          v_a( j,i) = v_cy( j-1,i)
        ELSE
          v_a( j,i) = 0._dp
        END IF
      ELSE
        IF     (ice%mask_ice_a( j-1,i) == 1 .AND. ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a( j+1,i) == 1) THEN
          v_a( j,i) = 0.5_dp * (v_cy( j-1,i) + v_cy( j,i))
        ELSEIF (ice%mask_ice_a( j-1,i) == 0 .AND. ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a( j+1,i) == 1) THEN
          v_a( j,i) = v_cy( j,i)
        ELSEIF (ice%mask_ice_a( j-1,i) == 1 .AND. ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a( j+1,i) == 0) THEN
          v_a( j,i) = v_cy( j-1,i)
        ELSE
          v_a( j,i) = 0._dp
        END IF
      END IF
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_velocity_from_cy_to_a_2D
  SUBROUTINE map_velocity_from_cy_to_a_3D( grid, ice, v_cy, v_a)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: v_cy
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: v_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (j == 1) THEN
        IF (ice%mask_ice_a( j,i) == 1) THEN
          v_a( :,j,i) = v_cy( :,j,i)
        ELSE
          v_a( :,j,i) = 0._dp
        END IF
      ELSEIF (j == grid%ny) THEN
        IF (ice%mask_ice_a( j,i) == 1) THEN
          v_a( :,j,i) = v_cy( :,j-1,i)
        ELSE
          v_a( :,j,i) = 0._dp
        END IF
      ELSE
        IF     (ice%mask_ice_a( j-1,i) == 1 .AND. ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a( j+1,i) == 1) THEN
          v_a( :,j,i) = 0.5_dp * (v_cy( :,j-1,i) + v_cy( :,j,i))
        ELSEIF (ice%mask_ice_a( j-1,i) == 0 .AND. ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a( j+1,i) == 1) THEN
          v_a( :,j,i) = v_cy( :,j,i)
        ELSEIF (ice%mask_ice_a( j-1,i) == 1 .AND. ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a( j+1,i) == 0) THEN
          v_a( :,j,i) = v_cy( :,j-1,i)
        ELSE
          v_a( :,j,i) = 0._dp
        END IF
      END IF
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE map_velocity_from_cy_to_a_3D
  
  ! Some routines for applying the semi-analytical GL flux solution
  ! as a boundary condition to the SSA/DIVA solver
  SUBROUTINE calculate_GL_flux( grid, ice)
    ! Use the semi-analytical solutions by Schoof (2007, in case of a Weertman-type sliding law),
    ! or by Tsai et al. (2012, in case of a Coulomb-type sliding law) to determine
    ! the grounding-line flux.
      
    IMPLICIT NONE
    
    ! In/output variables: 
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    LOGICAL                                            :: is_GL
    REAL(dp)                                           :: TAF1, TAF2, lambda_GL, Hi_GL, phi_fric_GL, A_flow_GL
    REAL(dp)                                           :: exponent_Schoof, C_sliding, factor_Schoof, factor_Tsai
    REAL(dp), PARAMETER                                :: Q0 = 0.61_dp
    REAL(dp)                                           :: Q_GL, u_GL
    REAL(dp)                                           :: dTAF_dx, dTAF_dy, abs_grad_TAF
    REAL(dp)                                           :: TAF3, TAF4, lambda_GL_w, lambda_GL_e, lambda_GL_s, lambda_GL_n
      
    ! Some constant terms for the Schoof grounding line flux
    exponent_Schoof  = (C%m_sliding + C%n_flow + 3._dp) / (C%m_sliding + 1._dp)
    C_sliding        = C%C_sliding / (sec_per_year**C%m_sliding)
    
    ! cx-grid
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    
      ! Initialise at zero
      ice%u_GL_cx( j,i) = 0._dp
      is_GL = .FALSE.
    
      IF ((ice%mask_sheet_a( j,i  ) == 1 .AND. ice%mask_shelf_a( j,i+1) == 1) .OR. &
          (ice%mask_sheet_a( j,i+1) == 1 .AND. ice%mask_shelf_a( j,i  ) == 1)) THEN
        ! "Proper" staggered GL
        
        is_GL = .TRUE.
        
        ! Interpolate ice thickness to the grounding-line position
        TAF1 = ice%TAF_a( j,i)
        TAF2 = ice%TAF_a( j,i+1)
        lambda_GL = TAF1 / (TAF1 - TAF2)
        Hi_GL = (ice%Hi_a( j,i) * (1._dp - lambda_GL)) + (ice%Hi_a( j,i+1) * lambda_GL)
        
        ! Take the friction angle and flow factor only from the grounded side
        IF (ice%mask_sheet_a( j,i) == 1) THEN
          phi_fric_GL = ice%phi_fric_a(   j,i)
          A_flow_GL   = ice%A_flow_vav_a( j,i)
        ELSE
          phi_fric_GL = ice%phi_fric_a(   j,i+1)
          A_flow_GL   = ice%A_flow_vav_a( j,i+1)
        END IF
        
        ! Calculate thickness-above-flotation gradient (for flow direction)
        dTAF_dx = (ice%TAF_a( j,i+1) - ice%TAF_a( j,i)) / grid%dx
        dTAF_dy = (ice%TAF_a( j+1,i+1) + ice%TAF_a( j+1,i) - ice%TAF_a( j-1,i+1) - ice%TAF_a( j-1,i)) / (4._dp * grid%dx)
        abs_grad_TAF = SQRT( dTAF_dx**2 + dTAF_dy**2)
        
      ELSEIF (ice%mask_sheet_a( j  ,i  ) == 1 .AND. ice%mask_sheet_a( j  ,i+1) == 1 .AND. &
              ice%mask_sheet_a( j-1,i  ) == 1 .AND. ice%mask_sheet_a( j-1,i+1) == 1 .AND. &
              ice%mask_shelf_a( j+1,i  ) == 1 .AND. ice%mask_shelf_a( j+1,i+1) == 1) THEN
        ! Parallel GL to the north
        
        is_GL = .TRUE.
        
        ! Interpolate ice thickness to the grounding-line position
        TAF1 = ice%TAF_a( j  ,i  )
        TAF2 = ice%TAF_a( j+1,i  )
        TAF3 = ice%TAF_a( j  ,i+1)
        TAF4 = ice%TAF_a( j+1,i+1)
        lambda_GL_w = TAF1 / (TAF1 - TAF2)
        lambda_GL_e = TAF3 / (TAF3 - TAF4)
        Hi_GL = 0.5_dp * ((ice%Hi_a( j  ,i  ) * (1._dp - lambda_GL_w)) + (ice%Hi_a( j+1,i  ) * lambda_GL_w)) + &
                0.5_dp * ((ice%Hi_a( j  ,i+1) * (1._dp - lambda_GL_e)) + (ice%Hi_a( j+1,i+1) * lambda_GL_e))
        
        ! Take the friction angle and flow factor only from the grounded side
        phi_fric_GL = 0.5_dp * (ice%phi_fric_a(   j,i) + ice%phi_fric_a(   j,i+1))
        A_flow_GL   = 0.5_dp * (ice%A_flow_vav_a( j,i) + ice%A_flow_vav_a( j,i+1))
        
        ! Calculate thickness-above-flotation gradient (for flow direction)
        dTAF_dx = (ice%TAF_a( j,i+1) + ice%TAF_a( j+1,i+1) - ice%TAF_a( j,i) - ice%TAF_a( j+1,i)) / (2._dp * grid%dx)
        dTAF_dy = (ice%TAF_a( j+1,i) + ice%TAF_a( j+1,i+1) - ice%TAF_a( j,i) - ice%TAF_a( j,i+1)) / (2._dp * grid%dx)
        abs_grad_TAF = SQRT( dTAF_dx**2 + dTAF_dy**2)
        
      ELSEIF (ice%mask_sheet_a( j  ,i  ) == 1 .AND. ice%mask_sheet_a( j  ,i+1) == 1 .AND. &
              ice%mask_sheet_a( j+1,i  ) == 1 .AND. ice%mask_sheet_a( j+1,i+1) == 1 .AND. &
              ice%mask_shelf_a( j-1,i  ) == 1 .AND. ice%mask_shelf_a( j-1,i+1) == 1) THEN
        ! Parallel GL to the south
        
        is_GL = .TRUE.
        
        ! Interpolate ice thickness to the grounding-line position
        TAF1 = ice%TAF_a( j  ,i  )
        TAF2 = ice%TAF_a( j-1,i  )
        TAF3 = ice%TAF_a( j  ,i+1)
        TAF4 = ice%TAF_a( j-1,i+1)
        lambda_GL_w = TAF1 / (TAF1 - TAF2)
        lambda_GL_e = TAF3 / (TAF3 - TAF4)
        Hi_GL = 0.5_dp * ((ice%Hi_a( j  ,i  ) * (1._dp - lambda_GL_w)) + (ice%Hi_a( j-1,i  ) * lambda_GL_w)) + &
                0.5_dp * ((ice%Hi_a( j  ,i+1) * (1._dp - lambda_GL_e)) + (ice%Hi_a( j-1,i+1) * lambda_GL_e))
        
        ! Take the friction angle and flow factor only from the grounded side
        phi_fric_GL = 0.5_dp * (ice%phi_fric_a(   j,i) + ice%phi_fric_a(   j,i+1))
        A_flow_GL   = 0.5_dp * (ice%A_flow_vav_a( j,i) + ice%A_flow_vav_a( j,i+1))
        
        ! Calculate thickness-above-flotation gradient (for flow direction)
        dTAF_dx = (ice%TAF_a( j-1,i+1) + ice%TAF_a( j  ,i+1) - ice%TAF_a( j-1,i  ) - ice%TAF_a( j  ,i  )) / (2._dp * grid%dx)
        dTAF_dy = (ice%TAF_a( j  ,i  ) + ice%TAF_a( j  ,i+1) - ice%TAF_a( j-1,i  ) - ice%TAF_a( j-1,i+1)) / (2._dp * grid%dx)
        abs_grad_TAF = SQRT( dTAF_dx**2 + dTAF_dy**2)
        
      END IF ! IF (ice%mask_GL_cx( j,i) == 1) THEN
      
      IF (is_GL) THEN

        ! The constant factor in the Schoof (2007) grounding line flux solution:
        factor_Schoof = ( A_flow_GL &
                         * ((ice_density * grav) ** (C%n_flow + 1._dp)) &
                         * ((1._dp - (ice_density / seawater_density)) ** C%n_flow) &
                         / ((4._dp**C%n_flow) * C_sliding) ) &
                         ** (1._dp / (REAL(C%m_sliding,dp) + 1._dp))
  
        ! The constant factor in the Tsai et al. (2015) grounding line flux solution:
        factor_Tsai   = ( 8._dp * Q0 * A_flow_GL &
                         * ((ice_density * grav) ** C%n_flow) &
                         * ((1._dp - (ice_density / seawater_density)) ** (C%n_flow - 1._dp)) &
                         / (4._dp**C%n_flow) )

        ! Calculate grounding line flux
        Q_GL = 0._dp
        IF     (C%choice_sliding_law == 'Coulomb_regularised') THEN
          Q_GL = factor_Tsai   * (Hi_GL**(C%n_flow + 2._dp)) / TAN(phi_fric_GL * (pi / 180._dp))
        ELSEIF (C%choice_sliding_law == 'Weertman') THEN
          Q_GL = factor_Schoof * Hi_GL**exponent_Schoof
        ELSE
          IF (par%master) WRITE(0,*) 'ERROR: choice_sliding_law "'//TRIM(C%choice_sliding_law)//'" not implemented in calculate_GL_flux!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
        ! Calculate grounding line speed
        u_GL = Q_GL / Hi_GL
        
        ! Calculate GL velocity components (assuming the flow direction is perpendicular to the grounding line,
        ! so parallel to the gradient of the thickness above flotation
        ice%Qx_GL_cx( j,i) = -Q_GL * dTAF_dx / abs_grad_TAF
        ice%u_GL_cx(  j,i) = -u_GL * dTAF_dx / abs_grad_TAF
        
      END IF ! IF (is_GL) THEN
      
    END DO
    END DO
    CALL sync
    
    ! cy-grid
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1
    
      ! Initialise at zero
      ice%v_GL_cy( j,i) = 0._dp
      is_GL = .FALSE.
    
      IF ((ice%mask_sheet_a( j  ,i) == 1 .AND. ice%mask_shelf_a( j+1,i) == 1) .OR. &
          (ice%mask_sheet_a( j+1,i) == 1 .AND. ice%mask_shelf_a( j  ,i) == 1)) THEN
        ! "Proper" staggered GL
        
        is_GL = .TRUE.
        
        ! Interpolate ice thickness to the grounding-line position
        TAF1 = ice%TAF_a( j,i)
        TAF2 = ice%TAF_a( j+1,i)
        lambda_GL = TAF1 / (TAF1 - TAF2)
        Hi_GL = (ice%Hi_a( j,i) * (1._dp - lambda_GL)) + (ice%Hi_a( j+1,i) * lambda_GL)
        
        ! Take the friction angle and flow factor only from the grounded side
        IF (ice%mask_sheet_a( j,i) == 1) THEN
          phi_fric_GL = ice%phi_fric_a(   j,i)
          A_flow_GL   = ice%A_flow_vav_a( j,i)
        ELSE
          phi_fric_GL = ice%phi_fric_a(   j+1,i)
          A_flow_GL   = ice%A_flow_vav_a( j+1,i)
        END IF
        
        ! Calculate thickness-above-flotation gradient (for flow direction)
        dTAF_dx = (ice%TAF_a( j+1,i+1) + ice%TAF_a( j,i+1) - ice%TAF_a( j+1,i-1) - ice%TAF_a( j,i-1)) / (4._dp * grid%dx)
        dTAF_dy = (ice%TAF_a( j+1,i) - ice%TAF_a( j,i)) / grid%dx
        abs_grad_TAF = SQRT( dTAF_dx**2 + dTAF_dy**2)
        
      ELSEIF (ice%mask_sheet_a( j  ,i  ) == 1 .AND. ice%mask_sheet_a( j+1,i  ) == 1 .AND. &
              ice%mask_sheet_a( j  ,i-1) == 1 .AND. ice%mask_sheet_a( j+1,i-1) == 1 .AND. &
              ice%mask_shelf_a( j  ,i+1) == 1 .AND. ice%mask_shelf_a( j+1,i+1) == 1) THEN
        ! Parallel GL to the east
        
        is_GL = .TRUE.
        
        ! Interpolate ice thickness to the grounding-line position
        TAF1 = ice%TAF_a( j  ,i  )
        TAF2 = ice%TAF_a( j  ,i+1)
        TAF3 = ice%TAF_a( j+1,i  )
        TAF4 = ice%TAF_a( j+1,i+1)
        lambda_GL_s = TAF1 / (TAF1 - TAF2)
        lambda_GL_n = TAF3 / (TAF3 - TAF4)
        Hi_GL = 0.5_dp * ((ice%Hi_a( j  ,i  ) * (1._dp - lambda_GL_s)) + (ice%Hi_a( j  ,i+1) * lambda_GL_s)) + &
                0.5_dp * ((ice%Hi_a( j+1,i  ) * (1._dp - lambda_GL_n)) + (ice%Hi_a( j+1,i+1) * lambda_GL_n))
        
        ! Take the friction angle and flow factor only from the grounded side
        phi_fric_GL = 0.5_dp * (ice%phi_fric_a(   j,i) + ice%phi_fric_a(   j+1,i))
        A_flow_GL   = 0.5_dp * (ice%A_flow_vav_a( j,i) + ice%A_flow_vav_a( j+1,i))
        
        ! Calculate thickness-above-flotation gradient (for flow direction)
        dTAF_dx = (ice%TAF_a( j,i+1) + ice%TAF_a( j+1,i+1) - ice%TAF_a( j,i) - ice%TAF_a( j+1,i)) / (2._dp * grid%dx)
        dTAF_dy = (ice%TAF_a( j+1,i) + ice%TAF_a( j+1,i+1) - ice%TAF_a( j,i) - ice%TAF_a( j,i+1)) / (2._dp * grid%dx)
        abs_grad_TAF = SQRT( dTAF_dx**2 + dTAF_dy**2)
        
      ELSEIF (ice%mask_sheet_a( j  ,i  ) == 1 .AND. ice%mask_sheet_a( j+1,i  ) == 1 .AND. &
              ice%mask_sheet_a( j  ,i+1) == 1 .AND. ice%mask_sheet_a( j+1,i+1) == 1 .AND. &
              ice%mask_shelf_a( j  ,i-1) == 1 .AND. ice%mask_shelf_a( j+1,i-1) == 1) THEN
        ! Parallel GL to the west
        
        is_GL = .TRUE.
        
        ! Interpolate ice thickness to the grounding-line position
        TAF1 = ice%TAF_a( j  ,i  )
        TAF2 = ice%TAF_a( j  ,i-1)
        TAF3 = ice%TAF_a( j+1,i  )
        TAF4 = ice%TAF_a( j+1,i-1)
        lambda_GL_s = TAF1 / (TAF1 - TAF2)
        lambda_GL_n = TAF3 / (TAF3 - TAF4)
        Hi_GL = 0.5_dp * ((ice%Hi_a( j  ,i  ) * (1._dp - lambda_GL_s)) + (ice%Hi_a( j  ,i-1) * lambda_GL_s)) + &
                0.5_dp * ((ice%Hi_a( j+1,i  ) * (1._dp - lambda_GL_n)) + (ice%Hi_a( j+1,i-1) * lambda_GL_n))
        
        ! Take the friction angle and flow factor only from the grounded side
        phi_fric_GL = 0.5_dp * (ice%phi_fric_a(   j,i) + ice%phi_fric_a(   j+1,i))
        A_flow_GL   = 0.5_dp * (ice%A_flow_vav_a( j,i) + ice%A_flow_vav_a( j+1,i))
        
        ! Calculate thickness-above-flotation gradient (for flow direction)
        dTAF_dx = (ice%TAF_a( j+1,i  ) + ice%TAF_a( j  ,i  ) - ice%TAF_a( j+1,i-1) - ice%TAF_a( j  ,i-1)) / (2._dp * grid%dx)
        dTAF_dy = (ice%TAF_a( j+1,i-1) + ice%TAF_a( j+1,i  ) - ice%TAF_a( j  ,i-1) - ice%TAF_a( j  ,i  )) / (2._dp * grid%dx)
        abs_grad_TAF = SQRT( dTAF_dx**2 + dTAF_dy**2)
        
      END IF ! IF (ice%mask_GL_cx( j,i) == 1) THEN
      
      IF (is_GL) THEN

        ! The constant factor in the Schoof (2007) grounding line flux solution:
        factor_Schoof = ( A_flow_GL &
                         * ((ice_density * grav) ** (C%n_flow + 1._dp)) &
                         * ((1._dp - (ice_density / seawater_density)) ** C%n_flow) &
                         / ((4._dp**C%n_flow) * C_sliding) ) &
                         ** (1._dp / (REAL(C%m_sliding,dp) + 1._dp))
  
        ! The constant factor in the Tsai et al. (2015) grounding line flux solution:
        factor_Tsai   = ( 8._dp * Q0 * A_flow_GL &
                         * ((ice_density * grav) ** C%n_flow) &
                         * ((1._dp - (ice_density / seawater_density)) ** (C%n_flow - 1._dp)) &
                         / (4._dp**C%n_flow) )

        ! Calculate grounding line flux
        Q_GL = 0._dp
        IF     (C%choice_sliding_law == 'Coulomb_regularised') THEN
          Q_GL = factor_Tsai   * (Hi_GL**(C%n_flow + 2._dp)) / TAN(phi_fric_GL * (pi / 180._dp))
        ELSEIF (C%choice_sliding_law == 'Weertman') THEN
          Q_GL = factor_Schoof * Hi_GL**exponent_Schoof
        ELSE
          IF (par%master) WRITE(0,*) 'ERROR: choice_sliding_law "'//TRIM(C%choice_sliding_law)//'" not implemented in calculate_GL_flux!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
        ! Calculate grounding line speed
        u_GL = Q_GL / Hi_GL
        
        ! Calculate GL velocity components (assuming the flow direction is perpendicular to the grounding line,
        ! so parallel to the gradient of the thickness above flotation
        ice%Qy_GL_cy( j,i) = -Q_GL * dTAF_dy / abs_grad_TAF
        ice%v_GL_cy(  j,i) = -u_GL * dTAF_dy / abs_grad_TAF
        
      END IF ! IF (is_GL) THEN
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calculate_GL_flux
  
END MODULE ice_velocity_module
