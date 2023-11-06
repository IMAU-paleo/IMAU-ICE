MODULE ice_velocity_module

  ! Contains all the routines needed to calculate instantaneous ice velocities for the
  ! current modelled ice-sheet geometry.

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_model_region, type_grid, type_ice_model, type_restart_data
  USE netcdf_module,                   ONLY: write_CSR_matrix_to_NetCDF
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             vertical_integration_from_bottom_to_zeta, vertical_average, &
                                             vertical_integrate, SSA_Schoof2006_analytical_solution, initialise_matrix_equation_CSR, &
                                             solve_matrix_equation_CSR, check_CSR_for_double_entries, is_floating, transpose_dp_2D
  USE derivatives_and_grids_module,    ONLY: ddx_cx_to_b_2D, ddy_cx_to_b_2D, ddy_cy_to_b_2D, ddx_cy_to_b_2D, &
                                             map_cx_to_a_2D, map_cy_to_a_2D, map_cx_to_a_3D, map_cy_to_a_3D, &
                                             map_cx_to_cy_2D, map_cy_to_cx_2D, map_a_to_cx_2D, map_a_to_cy_2D, &
                                             ddx_cx_to_cx_2D, ddy_cy_to_cx_2D, ddx_cx_to_cy_2D, ddy_cy_to_cy_2D, &
                                             ddx_cx_to_a_2D, ddy_cx_to_a_2D, ddx_cy_to_a_2D, ddy_cy_to_a_2D, &
                                             ddx_a_to_a_2D, ddy_a_to_a_2D
  USE basal_conditions_and_sliding_module, ONLY: calc_basal_conditions, calc_sliding_law
  USE general_ice_model_data_module,   ONLY: determine_grounded_fractions

  USE netcdf_debug_module,             ONLY: save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, save_variable_as_netcdf_int_3D, &
                                             save_variable_as_netcdf_dp_1D,  save_variable_as_netcdf_dp_2D,  save_variable_as_netcdf_dp_3D
  USE netcdf_input_module,             ONLY: read_field_from_xy_file_2D, read_field_from_file_2D

  IMPLICIT NONE

CONTAINS

  ! The different velocity equation solvers that can be called from run_ice_model
  SUBROUTINE solve_SIA(  grid, ice)

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_SIA'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: D_0
    REAL(dp), DIMENSION(C%nZ)                          :: D_deformation
    REAL(dp), DIMENSION(C%nZ)                          :: D_SIA_3D
    REAL(dp), PARAMETER                                :: D_uv_3D_cutoff = -1E5_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check that this routine is called correctly
    IF (.NOT. (C%choice_ice_dynamics == 'SIA' .OR. C%choice_ice_dynamics == 'SIA/SSA')) THEN
      CALL crash('should only be called when choice_ice_dynamics is set to "SIA" or "SIA/SSA"!')
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
    CALL calc_secondary_velocities( grid, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_SIA
  SUBROUTINE solve_SSA(  grid, ice)
    ! Calculate ice velocities using the SSA. Velocities are calculated on the staggered
    ! cx/cy-grids, using the discretisation scheme adopted from Yelmo/SICOPOLIS.

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_SSA'
    INTEGER                                            :: i,j
    LOGICAL                                            :: set_velocities_to_zero
    LOGICAL                                            :: has_converged
    INTEGER                                            :: viscosity_iteration_i
    REAL(dp)                                           :: resid_UV
    REAL(dp)                                           :: umax_analytical, tauc_analytical

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check that this routine is called correctly
    IF (.NOT. (C%choice_ice_dynamics == 'SSA' .OR. C%choice_ice_dynamics == 'SIA/SSA')) THEN
      CALL crash('should only be called when choice_ice_dynamics is set to "SSA" or "SIA/SSA"!')
    END IF

    ! If there's no grounded ice anywhere, don't bother
    set_velocities_to_zero = .FALSE.
    IF (SUM( ice%mask_sheet_a) == 0) set_velocities_to_zero = .TRUE.

    ! If we're prescribing no sliding, set velocities to zero
    IF (C%choice_sliding_law == 'no_sliding') set_velocities_to_zero = .TRUE.

    IF (set_velocities_to_zero) THEN
      ice%u_SSA_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 0._dp
      ice%v_SSA_cy( :,grid%i1:              grid%i2 ) = 0._dp
      CALL sync
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Calculate the driving stresses taudx, taudy
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (i < grid%nx) ice%taudx_cx( j,i) = -ice_density * grav * ice%Hi_cx( j,i) * ice%dHs_dx_cx( j,i)
      IF (j < grid%ny) ice%taudy_cy( j,i) = -ice_density * grav * ice%Hi_cy( j,i) * ice%dHs_dy_cy( j,i)
    END DO
    END DO
    CALL sync

    IF (par%master) print*, 'taudx_cx = ',SUM(ice%taudx_cx) !CvC

    ! Calculate the basal yield stress tau_c
    CALL calc_basal_conditions( grid, ice)

    ! Determine sub-grid grounded fractions for scaling the basal friction
    CALL determine_grounded_fractions( grid, ice)

    ! Find analytical solution for the SSA icestream experiment (used only to print numerical error to screen)
    CALL SSA_Schoof2006_analytical_solution( C%SSA_icestream_tantheta, C%SSA_icestream_H, ice%A_flow_vav_a( 1,1), 0._dp, umax_analytical, tauc_analytical)

    ! Initially set error very high
    ice%DIVA_err_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 1E5_dp
    ice%DIVA_err_cy( :,grid%i1:              grid%i2 ) = 1E5_dp
    CALL sync

    IF (C%do_read_velocities_from_restart) THEN
      ! do nothing
    ELSE
      ! The viscosity iteration
      viscosity_iteration_i = 0
      has_converged         = .FALSE.
      viscosity_iteration: DO WHILE (.NOT. has_converged)
        viscosity_iteration_i = viscosity_iteration_i + 1

        ! Calculate the effective viscosity and the product term N = eta * H
        CALL calc_effective_viscosity( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy)

        ! Calculate the sliding term beta (on both the A and Cx/Cy grids)
        CALL calc_sliding_term_beta( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy)

        ! Set beta_eff equal to beta; this turns the DIVA into the SSA
        ice%beta_eff_a( :,grid%i1:grid%i2) = ice%beta_a( :,grid%i1:grid%i2)

        ! Map beta_eff from the a-grid to the cx/cy-grids
        CALL map_a_to_cx_2D( grid, ice%beta_eff_a, ice%beta_eff_cx)
        CALL map_a_to_cy_2D( grid, ice%beta_eff_a, ice%beta_eff_cy)

        ! Apply the sub-grid grounded fraction
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          IF (i < grid%nx) ice%beta_eff_cx( j,i) = ice%beta_eff_cx( j,i) * ice%f_grnd_cx( j,i)**C%subgrid_friction_exponent
          IF (j < grid%ny) ice%beta_eff_cy( j,i) = ice%beta_eff_cy( j,i) * ice%f_grnd_cy( j,i)**C%subgrid_friction_exponent
          ice%beta_eff_a( j,i) = ice%beta_eff_a( j,i) * ice%f_grnd_a( j,i)**C%subgrid_friction_exponent ! Just for output
        END DO
        END DO
        CALL sync
        IF (par%master) print*, 'beta_eff_a = ',SUM(ice%beta_eff_a) !CvC

        ! Store the previous solution so we can check for convergence later
        ice%u_cx_prev( :,grid%i1:MIN(grid%nx-1,grid%i2)) = ice%u_SSA_cx( :,grid%i1:MIN(grid%nx-1,grid%i2))
        ice%v_cy_prev( :,grid%i1:              grid%i2 ) = ice%v_SSA_cy( :,grid%i1:              grid%i2 )
        CALL sync

        ! Solve the linearised DIVA with the SICOPOLIS solver
        CALL solve_DIVA_stag_linearised( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy)

        ! Apply velocity limits (both overflow and underflow) for improved stability
        CALL apply_velocity_limits( grid, ice%u_SSA_cx, ice%v_SSA_cy)

        ! "relax" subsequent viscosity iterations for improved stability
        CALL relax_DIVA_visc_iterations( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy, C%DIVA_visc_it_relax)

        ! Check if the viscosity iteration has converged
        CALL calc_visc_iter_UV_resid( grid, ice, ice%u_SSA_cx, ice%v_SSA_cy, resid_UV)
        !IF (par%master) WRITE(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ': resid_UV = ', resid_UV, ', u = [', MINVAL(ice%u_SSA_cx), ' - ', MAXVAL(ice%u_SSA_cx), ']'

        IF (par%master .AND. C%choice_refgeo_init_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'SSA_icestream') &
          WRITE(0,*) '    SSA - viscosity iteration ', viscosity_iteration_i, ': err = ', ABS(1._dp - MAXVAL(ice%u_SSA_cx) / umax_analytical), ': resid_UV = ', resid_UV

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
    END IF! (C%do_read_velocities_from_restart) THEN

    ! Calculate secondary velocities (surface, base, etc.)
    CALL calc_secondary_velocities( grid, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_SSA
  SUBROUTINE solve_DIVA( grid, ice)
    ! Calculate ice velocities using the DIVA. Velocities are calculated on the staggered
    ! cx/cy-grids, using the discretisation scheme adopted from Yelmo/SICOPOLIS.

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_DIVA'
    INTEGER                                            :: i,j
    LOGICAL                                            :: set_velocities_to_zero
    LOGICAL                                            :: has_converged
    INTEGER                                            :: viscosity_iteration_i
    REAL(dp)                                           :: resid_UV

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check that this routine is called correctly
    IF (.NOT. C%choice_ice_dynamics == 'DIVA') THEN
      CALL crash('should only be called when choice_ice_dynamics is set to "DIVA"!')
    END IF

    ! If there's no grounded ice anywhere, don't bother
    set_velocities_to_zero = .FALSE.
    IF (SUM( ice%mask_sheet_a) == 0) set_velocities_to_zero = .TRUE.

    IF (set_velocities_to_zero) THEN
      ice%u_vav_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 0._dp
      ice%v_vav_cy( :,grid%i1:              grid%i2 ) = 0._dp
      CALL sync
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Calculate the driving stresses taudx, taudy
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (i < grid%nx) ice%taudx_cx( j,i) = -ice_density * grav * ice%Hi_cx( j,i) * ice%dHs_dx_cx( j,i)
      IF (j < grid%ny) ice%taudy_cy( j,i) = -ice_density * grav * ice%Hi_cy( j,i) * ice%dHs_dy_cy( j,i)
    END DO
    END DO
    CALL sync

    ! Calculate the basal yield stress tau_c
    CALL calc_basal_conditions( grid, ice)

    ! Determine sub-grid grounded fractions for scaling the basal friction
    CALL determine_grounded_fractions( grid, ice)

    ! Initially set error very high
    ice%DIVA_err_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 1E5_dp
    ice%DIVA_err_cy( :,grid%i1:              grid%i2 ) = 1E5_dp
    CALL sync

    IF (C%do_read_velocities_from_restart) THEN
      ! do nothing
    ELSE
      ! The viscosity iteration
      viscosity_iteration_i = 0
      has_converged         = .FALSE.
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

        ! Store the previous solution so we can check for convergence later
        ice%u_cx_prev( :,grid%i1:MIN(grid%nx-1,grid%i2)) = ice%u_vav_cx( :,grid%i1:MIN(grid%nx-1,grid%i2))
        ice%v_cy_prev( :,grid%i1:              grid%i2 ) = ice%v_vav_cy( :,grid%i1:              grid%i2 )
        CALL sync

        ! Solve the linearised DIVA with the SICOPOLIS solver
        CALL solve_DIVA_stag_linearised( grid, ice, ice%u_vav_cx, ice%v_vav_cy)

        ! Apply velocity limits (both overflow and underflow) for improved stability
        CALL apply_velocity_limits( grid, ice%u_vav_cx, ice%v_vav_cy)

        ! "relax" subsequent viscosity iterations for improved stability
        CALL relax_DIVA_visc_iterations( grid, ice, ice%u_vav_cx, ice%v_vav_cy, C%DIVA_visc_it_relax)

        ! Check if the viscosity iteration has converged
        CALL calc_visc_iter_UV_resid( grid, ice, ice%u_vav_cx, ice%v_vav_cy, resid_UV)
        !IF (par%master) WRITE(0,*) '   DIVA - viscosity iteration ', viscosity_iteration_i, ': resid_UV = ', resid_UV, ', u = [', MINVAL(ice%u_vav_cx), ' - ', MAXVAL(ice%u_vav_cx), ']'

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
    END IF! (C%do_read_velocities_from_restart) THEN

    ! Calculate full 3D velocities
    CALL calc_3D_horizontal_velocities_DIVA( grid, ice)

    ! Calculate secondary velocities (surface, base, etc.)
    CALL calc_secondary_velocities( grid, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_DIVA

  ! Calculate "secondary" velocities (surface, basal, vertically averaged, on the A-grid, etc.)
  SUBROUTINE calc_secondary_velocities( grid, ice)
    ! Calculate "secondary" velocities (surface, basal, vertically averaged, on the A-grid, etc.)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_secondary_velocities'
    INTEGER                                            :: i,j,k

    ! Add routine to path
    CALL init_routine( routine_name)

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
      print*, "u_3D_cx = ", SUM(ice%u_3D_cx) !CvC

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
      CALL crash('unknown choice_ice_dynamics "' // TRIM(C%choice_ice_dynamics) // '"!')
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

    ! Calculate shearing ratio
    IF     (C%choice_ice_dynamics == 'SIA') THEN
      ! SIA = pure vertical shear
      ice%R_shear( :,grid%i1:grid%i2) = 1._dp
    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
      ! SSA = no vertical shear
      ice%R_shear( :,grid%i1:grid%i2) = 0._dp
    ELSE
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (ice%mask_ice_a( j,i) == 1 .AND. ice%uabs_base_a( j,i) /= ice%uabs_surf_a( j,i) ) THEN
          ice%R_shear( j,i) = MAX( 0._dp, MIN( 1._dp, (1._dp - ice%uabs_base_a( j,i) / ice%uabs_surf_a( j,i)) ))
        ELSE
          ice%R_shear( j,i) = 0._dp
        END IF
      END DO
      END DO
      CALL sync
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_secondary_velocities
  SUBROUTINE calc_3D_vertical_velocities( grid, ice)
    ! Use simple conservation of mass to calculate the vertical velocity w_3D

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_3D_vertical_velocities'
    INTEGER                                            :: i,j,k
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dHs_dx_a,  dHs_dy_a,  dHi_dx_a,  dHi_dy_a
    INTEGER                                            :: wdHs_dx_a, wdHs_dy_a, wdHi_dx_a, wdHi_dy_a
    REAL(dp)                                           :: dHbase_dx, dHbase_dy, u_base_a, v_base_a
    REAL(dp)                                           :: du_dx_k,   dv_dy_k
    REAL(dp)                                           :: du_dx_kp1, dv_dy_kp1
    REAL(dp)                                           :: w1, w2, w3, w4

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dHs_dx_a, wdHs_dx_a)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dHs_dy_a, wdHs_dy_a)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dHi_dx_a, wdHi_dx_a)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dHi_dy_a, wdHi_dy_a)

    ! Calculate surface & ice thickness gradients
    CALL ddx_a_to_a_2D( grid, ice%Hs_a, dHs_dx_a)
    CALL ddy_a_to_a_2D( grid, ice%Hs_a, dHs_dy_a)
    CALL ddx_a_to_a_2D( grid, ice%Hi_a, dHi_dx_a)
    CALL ddy_a_to_a_2D( grid, ice%Hi_a, dHi_dy_a)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Calculate the ice basal surface slope (not the same as bedrock slope when ice is floating!)
      dHbase_dx = dHs_dx_a( j,i) - dHi_dx_a( j,i)
      dHbase_dy = dHs_dy_a( j,i) - dHi_dy_a( j,i)

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
      IF (ice%mask_sheet_a( j,i) == 1) THEN
        ice%w_3D_a( C%nz,j,i) = (u_base_a * dHbase_dx) + (v_base_a * dHbase_dy) + ice%dHb_dt_a( j,i)
      ELSE
        ice%w_3D_a( C%nz,j,i) = (u_base_a * dHbase_dx) + (v_base_a * dHbase_dy)  ! Should this include the ice thinning rate?
      END IF

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

    ! Surface and basal vertical velocities

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ice%w_surf_a( j,i) = ice%w_3D_a( 1   ,j,i)
      ice%w_base_a( j,i) = ice%w_3D_a( C%nz,j,i)

    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wdHs_dx_a)
    CALL deallocate_shared( wdHs_dy_a)
    CALL deallocate_shared( wdHi_dx_a)
    CALL deallocate_shared( wdHi_dy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_3D_vertical_velocities

  ! Calculating some physical terms (basal yield stress, effective viscosity, etc.)
  SUBROUTINE calc_vertical_shear_strain_rates( grid, ice)
    ! Calculate vertical shear rates (Lipscomb et al. 2019, Eq. 36)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_vertical_shear_strain_rates'
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: visc_eff_cx, visc_eff_cy

    ! Add routine to path
    CALL init_routine( routine_name)

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
            ice%dv_dz_3D_cy( k,j,i) = (ice%taub_cy( j,i) / visc_eff_cy) * C%zeta( k)
          END DO
        END IF
      END IF

    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_3D( ice%du_dz_3D_cx, 'ice%du_dz_3D_cx')
    CALL check_for_NaN_dp_3D( ice%dv_dz_3D_cy, 'ice%dv_dz_3D_cy')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_shear_strain_rates
  SUBROUTINE calc_effective_viscosity( grid, ice, u_cx, v_cy)
    ! Calculate 3D effective viscosity following Lipscomb et al. (2019), Eq. 2

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_cx
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_cy

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_effective_viscosity'
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: du_dz_a, dv_dz_a
    REAL(dp)                                           :: eps_sq, w
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  du_dx_a,  du_dy_a,  dv_dx_a,  dv_dy_a
    INTEGER                                            :: wdu_dx_a, wdu_dy_a, wdv_dx_a, wdv_dy_a

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, du_dx_a, wdu_dx_a)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, du_dy_a, wdu_dy_a)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dv_dx_a, wdv_dx_a)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dv_dy_a, wdv_dy_a)

    ! Calculate effective strain components from horizontal stretching on the regular grid
    CALL ddx_cx_to_a_2D( grid, u_cx, du_dx_a)
    CALL ddy_cx_to_a_2D( grid, u_cx, du_dy_a)
    CALL ddx_cy_to_a_2D( grid, v_cy, dv_dx_a)
    CALL ddy_cy_to_a_2D( grid, v_cy, dv_dy_a)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (C%choice_ice_margin == 'BC') THEN
        IF (ice%mask_ice_a( j,i) == 0) THEN
          ! No ice
          ice%visc_eff_3D_a( :,j,i) = 0._dp
          CYCLE
        END IF
      ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
        ! In the "infinite slab" case, calculate effective viscosity everywhere
        ! (even when there's technically no ice present)
      ELSE
        CALL crash('unknown choice_ice_margin "' // TRIM(C%choice_ice_margin) // '"!')
      END IF

      DO k = 1, C%nz

        ! Calculate vertical shear terms on the a-grid
        IF     (i == 1) THEN
          du_dz_a = ice%du_dz_3D_cx( k,j,1)
        ELSEIF (i == grid%nx) THEN
          du_dz_a = ice%du_dz_3D_cx( k,j,grid%nx-1)
        ELSE
          du_dz_a = calc_staggered_margin( ice%du_dz_3D_cx( k,j,i-1), ice%du_dz_3D_cx( k,j,i), ice%Hi_a( j,i-1), ice%Hi_a( j,i))
        END IF
        IF     (j == 1) THEN
          dv_dz_a = ice%dv_dz_3D_cy( k,1,i)
        ELSEIF (j == grid%ny) THEN
          dv_dz_a = ice%dv_dz_3D_cy( k,grid%ny-1,i)
        ELSE
          dv_dz_a = calc_staggered_margin( ice%dv_dz_3D_cy( k,j-1,i), ice%dv_dz_3D_cy( k,j,i), ice%Hi_a( j-1,i), ice%Hi_a( j,i))
        END IF

        ! Calculate the total effective strain rate from L19, Eq. 21
        eps_sq = du_dx_a( j,i)**2 + &
                 dv_dy_a( j,i)**2 + &
                 du_dx_a( j,i) * dv_dy_a( j,i) + &
                 0.25_dp * (du_dy_a( j,i) + dv_dx_a( j,i))**2 + &
                 0.25_dp * (du_dz_a**2 + dv_dz_a**2) + &
                 C%DIVA_epsilon_sq_0

        ! Calculate effective viscosity
        ice%visc_eff_3D_a( k,j,i) = 0.5_dp * ice%A_flow_3D_a( k,j,i)**(-1._dp/C%n_flow) * (eps_sq)**((1._dp - C%n_flow)/(2._dp*C%n_flow))

        ! Safety
        ice%visc_eff_3D_a( k,j,i) = MAX( ice%visc_eff_3D_a( k,j,i), C%DIVA_visc_eff_min)

      END DO ! DO k = 1, C%nz

      ! Vertical integral
      ice%visc_eff_int_a( j,i) = vertical_integrate( ice%visc_eff_3D_a( :,j,i))

      ! Calculate the product term N = eta * H
      IF (C%choice_ice_margin == 'BC') THEN
        ice%N_a( j,i)  = ice%visc_eff_int_a( j,i) * ice%Hi_a( j,i)
      ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
        ice%N_a( j,i)  = ice%visc_eff_int_a( j,i) * MAX(0.1_dp,ice%Hi_a( j,i))
      END IF

    END DO
    END DO
    CALL sync

    ! Stagger from A-grid to B-grid
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny-1

      IF (C%choice_ice_margin == 'BC') THEN

        ice%visc_eff_3D_b( :,j,i) = 0._dp
        w                         = 0._dp

        IF (ice%mask_ice_a( j  ,i  ) == 1) THEN
          ice%visc_eff_3D_b( :,j,i) = ice%visc_eff_3D_b( :,j,i) + ice%visc_eff_3D_a( :,j  ,i  )
          w = w + 1._dp
        END IF
        IF (ice%mask_ice_a( j  ,i+1) == 1) THEN
          ice%visc_eff_3D_b( :,j,i) = ice%visc_eff_3D_b( :,j,i) + ice%visc_eff_3D_a( :,j  ,i+1)
          w = w + 1._dp
        END IF
        IF (ice%mask_ice_a( j+1,i  ) == 1) THEN
          ice%visc_eff_3D_b( :,j,i) = ice%visc_eff_3D_b( :,j,i) + ice%visc_eff_3D_a( :,j+1,i  )
          w = w + 1._dp
        END IF
        IF (ice%mask_ice_a( j+1,i+1) == 1) THEN
          ice%visc_eff_3D_b( :,j,i) = ice%visc_eff_3D_b( :,j,i) + ice%visc_eff_3D_a( :,j+1,i+1)
          w = w + 1._dp
        END IF

        IF (w > 0._dp) ice%visc_eff_3D_b( :,j,i) = ice%visc_eff_3D_b( :,j,i) / w

      ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
        ! In the "infinite slab" case, calculate effective viscosity everywhere
        ! (even when there's technically no ice present)

        ice%visc_eff_3D_b( :,j,i) = &
          (ice%visc_eff_3D_a( :,j  ,i  ) + &
           ice%visc_eff_3D_a( :,j  ,i+1) + &
           ice%visc_eff_3D_a( :,j+1,i  ) + &
           ice%visc_eff_3D_a( :,j+1,i+1)) / 4._dp

      ELSE
        CALL crash('unknown choice_ice_margin "' // TRIM(C%choice_ice_margin) // '"!')
      END IF

      ! Vertical integral
      ice%visc_eff_int_b( j,i) = vertical_integrate( ice%visc_eff_3D_b( :,j,i))

      ! Product term N = eta * H
      IF (C%choice_ice_margin == 'BC') THEN
        ice%N_b( j,i) = ice%visc_eff_int_b( j,i) * ice%Hi_b( j,i)
      ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
        ice%N_b( j,i) = ice%visc_eff_int_b( j,i) * MAX(0.1_dp,ice%Hi_b( j,i))
      END IF

    END DO
    END DO
    CALL sync

    ! If we're using the sans-cross-terms SSA/DIVA solver, get N on the cx/cy-grids
    IF (.NOT. C%include_SSADIVA_crossterms) THEN
      CALL map_a_to_cx_2D( grid, ice%N_a, ice%N_cx)
      CALL map_a_to_cy_2D( grid, ice%N_a, ice%N_cy)
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wdu_dx_a)
    CALL deallocate_shared( wdu_dy_a)
    CALL deallocate_shared( wdv_dx_a)
    CALL deallocate_shared( wdv_dy_a)

    ! Safety
    CALL check_for_NaN_dp_3D( ice%visc_eff_3D_a,  'ice%visc_eff_3D_a' )
    CALL check_for_NaN_dp_3D( ice%visc_eff_3D_b,  'ice%visc_eff_3D_b' )
    CALL check_for_NaN_dp_2D( ice%visc_eff_int_a, 'ice%visc_eff_int_a')
    CALL check_for_NaN_dp_2D( ice%visc_eff_int_b, 'ice%visc_eff_int_b')
    CALL check_for_NaN_dp_2D( ice%N_a,            'ice%N_a'           )
    CALL check_for_NaN_dp_2D( ice%N_b,            'ice%N_b'           )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_effective_viscosity
  SUBROUTINE calc_sliding_term_beta( grid, ice, u_cx, v_cy)
    ! Calculate the sliding term beta in the SSA/DIVA using the specified sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_cx
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_cy

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_term_beta'
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  u_a,  v_a
    INTEGER                                            :: wu_a, wv_a
    REAL(dp)                                           :: beta_min

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, u_a, wu_a)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, v_a, wv_a)

    ! Get velocities on the a-grid
    CALL map_cx_to_a_2D( grid, u_cx, u_a)
    CALL map_cy_to_a_2D( grid, v_cy, v_a)

    ! Calculate the basal friction coefficients beta on the a-grid
    CALL calc_sliding_law( grid, ice, u_a, v_a, ice%beta_a)

    ! Clean up after yourself
    CALl deallocate_shared( wu_a)
    CALl deallocate_shared( wv_a)

    ! Limit minimum beta over thin-ice areas to prevent super-fast ice on steep slopes
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ! Multiply config-reference minimum beta by a number between 0 and 1,
      ! depending on how much thinner than a config-reference threshold
      ! the current model thickness is. If thicker, no limit is applied.
      beta_min = C%DIVA_beta_min_thin_ice * MIN( 1._dp, MAX( 0._dp, 1._dp - ice%Hi_a( j,i) / C%DIVA_beta_min_thin_ice_H0))
      ! Apply dynamic minimum limit to beta
      ice%beta_a( j,i) = MAX( ice%beta_a( j,i), beta_min )
    END DO
    END DO
    CALL sync

    ! Limit maximum beta to improve stability
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%beta_a( j,i) = MIN( ice%beta_a( j,i), C%DIVA_beta_max)
    END DO
    END DO
    CALL sync

    ! LEGACY - Apply the flotation mask (this is how we did it before we introduced the PISM/CISM-style sub-grid grounded fraction)
    IF (.NOT. C%do_GL_subgrid_friction) THEN
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (ice%mask_ocean_a( j,i) == 1) ice%beta_a( j,i) = 0._dp
      END DO
      END DO
      CALL sync
    END IF

    ! Safety
    CALL check_for_NaN_dp_2D( ice%beta_a, 'ice%beta_a')

    ! Finalise routine path
      CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_term_beta
  SUBROUTINE calc_F_integral( grid, ice, n)
    ! Calculate the integral F2 in Lipscomb et al. (2019), Eq. 30

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: n

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_F_integral'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: F_int_min

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine the minimum value of F_int, to assign when H_ice == 0,
    ! since F_int should be nonzero everywhere for numerics
    F_int_min = vertical_integrate( (1._dp / C%DIVA_visc_eff_min) * C%zeta**n)

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
          CALL crash('unknown choice_ice_margin "' // TRIM(C%choice_ice_margin) // '"!')
        END IF

      END IF

    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( ice%F2_a, 'ice%F2_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_F_integral
  SUBROUTINE calc_beta_eff( grid, ice)
    ! Calculate the "effective basal friction" beta_eff, used in the DIVA

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_beta_eff'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate beta_eff on the a-grid
    IF (C%choice_sliding_law == 'no_sliding') THEN
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

        ! Exception for basal freezing
        IF (C%include_basal_freezing) THEN
          IF (ice%Ti_pmp_a( C%nz,j,i) - ice%Ti_a( C%nz,j,i) > C%deltaT_basal_freezing) THEN
            ice%beta_eff_a( j,i) = 1._dp / ice%F2_a( j,i)
          END IF
        END IF

      END DO
      END DO
      CALL sync

    END IF

    ! Map beta_eff from the a-grid to the cx/cy-grids
    CALL map_a_to_cx_2D( grid, ice%beta_eff_a, ice%beta_eff_cx)
    CALL map_a_to_cy_2D( grid, ice%beta_eff_a, ice%beta_eff_cy)

    ! Apply the sub-grid grounded fraction
    IF (C%do_GL_subgrid_friction) THEN
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (i < grid%nx) ice%beta_eff_cx( j,i) = ice%beta_eff_cx( j,i) * ice%f_grnd_cx( j,i)**C%subgrid_friction_exponent
        IF (j < grid%ny) ice%beta_eff_cy( j,i) = ice%beta_eff_cy( j,i) * ice%f_grnd_cy( j,i)**C%subgrid_friction_exponent
      END DO
      END DO
      CALL sync
    END IF

    ! Safety
    CALL check_for_NaN_dp_2D( ice%beta_eff_a , 'ice%beta_eff_a' )
    CALL check_for_NaN_dp_2D( ice%beta_eff_cx, 'ice%beta_eff_cx')
    CALL check_for_NaN_dp_2D( ice%beta_eff_cy, 'ice%beta_eff_cy')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_beta_eff
  SUBROUTINE calc_basal_stress( grid, ice)
    ! Calculate the basal stress resulting from sliding (friction times velocity)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_stress'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ! Lipscomp et al., 2019, just above Eq. 33
      IF (i < grid%nx) ice%taub_cx( j,i) = ice%beta_eff_cx( j,i) * ice%u_vav_cx( j,i)
      IF (j < grid%ny) ice%taub_cy( j,i) = ice%beta_eff_cy( j,i) * ice%v_vav_cy( j,i)

    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( ice%taub_cx, 'ice%taub_cx')
    CALL check_for_NaN_dp_2D( ice%taub_cy, 'ice%taub_cy')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_stress
  SUBROUTINE calc_basal_velocities( grid, ice)
    ! Calculate basal sliding following Lipscomp et al., 2019, Eq. 32 
    ! (or it can also be obtained from Goldberg (2011), Eq. 34)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_velocities'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: F2_stag

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_sliding_law == 'no_sliding') THEN
      ! Exception for the case of no sliding: set basal velocities to zero
      ! (this comes out naturally more or less with beta_eff set as above,
      !  but ensuring basal velocity is zero adds stability)
      ice%u_base_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = 0._dp
      ice%v_base_cy( :,grid%i1:              grid%i2 ) = 0._dp
      CALL sync

      ! Finalise routine path and RETURN
      CALL finalise_routine( routine_name)
      RETURN

    END IF ! IF (C%choice_sliding_law == 'no_sliding') THEN

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! x-direction
      IF (i < grid%nx) THEN

        ! Stagger the F2 integral to the cx-grid
        F2_stag = calc_staggered_margin( ice%F2_a( j,i), ice%F2_a( j,i+1), ice%Hi_a( j,i), ice%Hi_a( j,i+1))

        ! Calculate basal velocity component (Lipscomb et al., 2019, Eq. 32)
        ice%u_base_cx( j,i) = ice%u_vav_cx( j,i) / (1._dp + ice%beta_a( j,i) * F2_stag)

        ! Exception for basal freezing
        IF (C%include_basal_freezing) THEN
          IF (ice%Ti_pmp_a( C%nz,j,i) - ice%Ti_a( C%nz,j,i) > C%deltaT_basal_freezing) THEN
            ice%u_base_cx( j,i) = 0._dp
          END IF
        END IF

      END IF ! IF (i < grid%nx) THEN

      ! y-direction
      IF (j < grid%ny) THEN

        ! Stagger the F2 integral to the cy-grid
        F2_stag = calc_staggered_margin( ice%F2_a( j,i), ice%F2_a( j+1,i), ice%Hi_a( j,i), ice%Hi_a( j+1,i))

        ! Calculate basal velocity component (Lipscomb et al., 2019, Eq. 32)
        ice%v_base_cy( j,i) = ice%v_vav_cy( j,i) / (1._dp + ice%beta_a( j,i) * F2_stag)

        ! Exception for basal freezing
        IF (C%include_basal_freezing) THEN
          IF (ice%Ti_pmp_a( C%nz,j,i) - ice%Ti_a( C%nz,j,i) > C%deltaT_basal_freezing) THEN
            ice%v_base_cy( j,i) = 0._dp
          END IF
        END IF

      END IF ! IF (j < grid%ny) THEN

    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( ice%u_base_cx, 'ice%u_base_cx')
    CALL check_for_NaN_dp_2D( ice%v_base_cy, 'ice%v_base_cy')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_velocities
  SUBROUTINE calc_3D_horizontal_velocities_DIVA( grid, ice)
    ! Caluculate the 3D horizontal velocity field (following Lipscomb et al., 2019, Eq. 29)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_3D_horizontal_velocities_DIVA'
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION( C%nz)                         :: F1_stag

    ! Add routine to path
    CALL init_routine( routine_name)

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
          CALL crash('unknown choice_ice_margin "' // TRIM(C%choice_ice_margin) // '"!')
        END IF

        IF (C%choice_sliding_law == 'no_sliding') THEN
        ! Exception for the case of no sliding
        
          ! Lipscomb et al., 2019, Eq. 29, and text between Eqs. 33 and 34
          ice%u_3D_cx( :,j,i) = ice%taub_cx( j,i) * F1_stag

        ELSE
        ! Calculate velocity column
        
          ! Lipscomb et al., 2019, Eq. 29
          ice%u_3D_cx( :,j,i) = ice%u_base_cx( j,i) * (1._dp + ice%beta_a( j,i) * F1_stag)
        
        END IF ! IF (C%choice_sliding_law == 'no_sliding') THEN

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
          CALL crash('unknown choice_ice_margin "' // TRIM(C%choice_ice_margin) // '"!')
        END IF

        IF (C%choice_sliding_law == 'no_sliding') THEN
        ! Exception for the case of no sliding

          ! Lipscomb et al., 2019, Eq. 29, and text between Eqs. 33 and 34 
          ice%v_3D_cy( :,j,i) = ice%taub_cy( j,i) * F1_stag

        ELSE
        ! Calculate velocity column

          ! Lipscomb et al., 2019, Eq. 29
          ice%v_3D_cy( :,j,i) = ice%v_base_cy( j,i) * (1._dp + ice%beta_a( j,i) * F1_stag)
        
        END IF ! IF (C%choice_sliding_law == 'no_sliding') THEN

      END IF

    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_3D( ice%u_3D_cx, 'ice%u_3D_cx')
    CALL check_for_NaN_dp_3D( ice%v_3D_cy, 'ice%v_3D_cy')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_3D_horizontal_velocities_DIVA

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
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'set_sico_masks'
    INTEGER                                            :: i, j, i1, i2, j1, j2
    LOGICAL                                            :: is_float
    LOGICAL, PARAMETER                                 :: disable_grounded_fronts = .TRUE.

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_sico_masks
  SUBROUTINE apply_velocity_limits( grid, u_cx, v_cy)
    ! Apply a velocity limit (for stability)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: u_cx, v_cy

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_velocity_limits'
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  u_cy,  v_cx,  uabs_cx,  uabs_cy
    INTEGER                                            :: wu_cy, wv_cx, wuabs_cx, wuabs_cy

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL allocate_shared_dp_2D( grid%ny  , grid%nx-1, v_cx,    wv_cx   )
    CALL allocate_shared_dp_2D( grid%ny-1, grid%nx  , u_cy,    wu_cy   )
    CALL allocate_shared_dp_2D( grid%ny  , grid%nx-1, uabs_cx, wuabs_cx)
    CALL allocate_shared_dp_2D( grid%ny-1, grid%nx  , uabs_cy, wuabs_cy)

    CALL map_cx_to_cy_2D( grid, u_cx, u_cy)
    CALL map_cy_to_cx_2D( grid, v_cy, v_cx)

    DO i = grid%i1, MIN(grid%nx-1, grid%i2)
    DO j = 1, grid%ny
      uabs_cx( j,i) = SQRT( u_cx( j,i)**2 + v_cx( j,i)**2)
    END DO
    END DO
    CALL sync

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      uabs_cy( j,i) = SQRT( u_cy( j,i)**2 + v_cy( j,i)**2)
    END DO
    END DO
    CALL sync

    DO i = grid%i1, MIN(grid%nx-1, grid%i2)
    DO j = 1, grid%ny
      IF (uabs_cx( j,i) > C%DIVA_vel_max) THEN
        u_cx( j,i) = u_cx( j,i) * C%DIVA_vel_max / uabs_cx( j,i)
      END IF
    END DO
    END DO
    CALL sync

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      IF (uabs_cy( j,i) > C%DIVA_vel_max) THEN
        v_cy( j,i) = v_cy( j,i) * C%DIVA_vel_max / uabs_cy( j,i)
      END IF
    END DO
    END DO
    CALL sync

    CALL deallocate_shared( wv_cx   )
    CALL deallocate_shared( wu_cy   )
    CALL deallocate_shared( wuabs_cx)
    CALL deallocate_shared( wuabs_cy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_velocity_limits
  SUBROUTINE relax_DIVA_visc_iterations( grid, ice, u_cx, v_cy, rel)
    ! Relax velocity solution with results of the previous viscosity iteration

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: u_cx
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: v_cy
    REAL(dp),                            INTENT(IN)    :: rel

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'relax_DIVA_visc_iterations'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (i < grid%nx) u_cx( j,i) = rel * u_cx( j,i) + (1._dp - rel) * ice%u_cx_prev( j,i)
      IF (j < grid%ny) v_cy( j,i) = rel * v_cy( j,i) + (1._dp - rel) * ice%v_cy_prev( j,i)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE relax_DIVA_visc_iterations
  SUBROUTINE calc_visc_iter_UV_resid( grid, ice, u_cx, v_cy, resid_UV)
    ! Check if the viscosity iteration has converged to a stable solution

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_cx
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_cy
    REAL(dp),                            INTENT(OUT)   :: resid_UV

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_visc_iter_UV_resid'
    INTEGER                                            :: ierr
    INTEGER                                            :: i,j,ncx,ncy
    REAL(dp)                                           :: res1, res2
    REAL(dp), PARAMETER                                :: DIVA_vel_tolerance = 1e-6   ! [m/a] only consider points with velocity above this tolerance limit

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the L2 norm based on velocity solution between previous
    ! and current viscosity iteration (as in Yelmo/SICOPOLIS)

    ncx  = 0
    ncy  = 0
    res1 = 0._dp
    res2 = 0._dp

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (i < grid%nx) THEN
        IF (ABS(u_cx( j,i)) > DIVA_vel_tolerance) THEN
          ncx = ncx + 1
          res1 = res1 + (u_cx( j,i) - ice%u_cx_prev( j,i))**2._dp
          res2 = res2 + (u_cx( j,i) + ice%u_cx_prev( j,i))**2._dp
        END IF
      END IF

      IF (j < grid%ny) THEN
        IF (ABS(v_cy( j,i)) > DIVA_vel_tolerance) THEN
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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_visc_iter_UV_resid
  SUBROUTINE estimate_visc_iter_UV_errors( grid, ice, u_cx, v_cy)
    ! Estimate the error in the velocity fields so we can
    ! update the SICOPOLIS-style DIVA solving masks (for better efficiency)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_cx
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_cy

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'estimate_visc_iter_UV_errors'
    INTEGER                                            :: i,j
    REAL(dp), PARAMETER                                :: DIVA_vel_tolerance = 1e-6   ! [m/a] only consider points with velocity above this tolerance limit
    REAL(dp), PARAMETER                                :: tol = 1E-5_dp

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

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

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_DIVA_stag_linearised'
    INTEGER                                            :: i, j, k, n

    ! Add routine to path
    CALL init_routine( routine_name)

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
        CALL crash('unknown choice_ice_margin "' // TRIM(C%choice_ice_margin) // '"!')
      END IF

    END IF ! IF (par%master) THEN
    CALL sync

    ! Solve the matrix equation using SOR
    ! ===================================

    CALL solve_matrix_equation_CSR( ice%DIVA_m, C%DIVA_choice_matrix_solver, &
      SOR_nit = ice%DIVA_SOR_nit, SOR_tol = ice%DIVA_SOR_tol, SOR_omega = ice%DIVA_SOR_omega, &
      PETSc_rtol = ice%DIVA_PETSc_rtol, PETSc_abstol = ice%DIVA_PETSc_abstol)

    ! Map the solution back from vector format to the model grids
    ! ===========================================================

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (i < grid%nx) THEN
        n = ice%DIVA_m_ij2n_u( j,i)
        u( j,i) = ice%DIVA_m%x( n)
      END IF
      IF (j < grid%ny) THEN
        n = ice%DIVA_m_ij2n_v( j,i)
        v( j,i) = ice%DIVA_m%x( n)
      END IF
    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( u, 'u')
    CALL check_for_NaN_dp_2D( v, 'v')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_DIVA_stag_linearised
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_with_BC(       grid, ice, u, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u
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
    IF (C%include_SSADIVA_crossterms) THEN
      CALL list_DIVA_matrix_coefficients_eq_1_free(      grid, ice, u, k, n)
    ELSE
      CALL list_DIVA_matrix_coefficients_eq_1_free_sans( grid, ice, u, k, n)
    END IF

  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_with_BC
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_with_BC(       grid, ice, v, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v
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
    IF (C%include_SSADIVA_crossterms) THEN
      CALL list_DIVA_matrix_coefficients_eq_2_free(      grid, ice, v, k, n)
    ELSE
      CALL list_DIVA_matrix_coefficients_eq_2_free_sans( grid, ice, v, k, n)
    END IF

  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_with_BC
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_infinite_slab( grid, ice, u, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u
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

    ! No more exceptions; solve the complete, free DIVA
    IF (C%include_SSADIVA_crossterms) THEN
      CALL list_DIVA_matrix_coefficients_eq_1_free(      grid, ice, u, k, n)
    ELSE
      CALL list_DIVA_matrix_coefficients_eq_1_free_sans( grid, ice, u, k, n)
    END IF

  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_infinite_slab
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_infinite_slab( grid, ice, v, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v
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

    ! No more exceptions; solve the complete, free DIVA
    IF (C%include_SSADIVA_crossterms) THEN
      CALL list_DIVA_matrix_coefficients_eq_2_free(      grid, ice, v, k, n)
    ELSE
      CALL list_DIVA_matrix_coefficients_eq_2_free_sans( grid, ice, v, k, n)
    END IF

  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_infinite_slab
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_boundary(      grid, ice, u, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists
    ! Exception for boundary conditions for the domain boundary

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n

    ! Local variables:
    INTEGER                                            :: i,j,i_in,i_opp,j_in,j_opp,nu,nu_in,nu_opp
    CHARACTER(LEN=256)                                 :: BC

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,1)
    j = ice%DIVA_m_n2ij_uv( n,2)

    ! Determine what kind of boundary conditions to apply
    IF (i == 1 .AND. j == 1) THEN
      ! Southwest corner
      IF     (C%DIVA_boundary_BC_u_south == 'zero' .OR. C%DIVA_boundary_BC_u_west == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_u_south == 'infinite' .OR. C%DIVA_boundary_BC_u_west == 'infinite') THEN
        BC = 'infinite'
      ELSE
        BC = 'periodic'
      END IF
    ELSEIF (i == grid%nx-1 .AND. j == 1) THEN
      ! Southeast corner
      IF     (C%DIVA_boundary_BC_u_south == 'zero' .OR. C%DIVA_boundary_BC_u_east == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_u_south == 'infinite' .OR. C%DIVA_boundary_BC_u_east == 'infinite') THEN
        BC = 'infinite'
      ELSE
        BC = 'periodic'
      END IF
    ELSEIF (i == 1 .AND. j == grid%ny) THEN
      ! Northwest corner
      IF     (C%DIVA_boundary_BC_u_north == 'zero' .OR. C%DIVA_boundary_BC_u_west == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_u_north == 'infinite' .OR. C%DIVA_boundary_BC_u_west == 'infinite') THEN
        BC = 'infinite'
      ELSE
        BC = 'periodic'
      END IF
    ELSEIF (i == grid%nx-1 .AND. j == grid%ny) THEN
      ! Northeast corner
      IF     (C%DIVA_boundary_BC_u_north == 'zero' .OR. C%DIVA_boundary_BC_u_east == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_u_north == 'infinite' .OR. C%DIVA_boundary_BC_u_east == 'infinite') THEN
        BC = 'infinite'
      ELSE
        BC = 'periodic'
      END IF
    ELSEIF (i == 1) THEN
      ! Western boundary
      BC = C%DIVA_boundary_BC_u_west
    ELSEIF (i == grid%nx-1) THEN
      ! Eastern boundary
      BC = C%DIVA_boundary_BC_u_east
    ELSEIF (j == 1) THEN
      ! Southern boundary
      BC = C%DIVA_boundary_BC_u_south
    ELSEIF (j == grid%ny) THEN
      ! Northern boundary
      BC = C%DIVA_boundary_BC_u_north
    ELSE
      CALL crash('list_DIVA_matrix_coefficients_eq_1_boundary: grid cell [{int_01},{int_02}] doesnt lie on the domain boundary!', int_01 = j, int_02 = i)
    END IF

    ! Find indices of relevant grid cells
    IF     (i == 1) THEN
      i_in  = 2
      i_opp = grid%nx-2
    ELSEIF (i == grid%nx-1) THEN
      i_in  = grid%nx-2
      i_opp = 2
    ELSE
      i_in  = i
      i_opp = i
    END IF
    IF     (j == 1) THEN
      j_in  = 2
      j_opp = grid%ny-1
    ELSEIF (j == grid%ny) THEN
      j_in  = grid%ny-1
      j_opp = 2
    ELSE
      j_in  = j
      j_opp = j
    END IF

    nu     = ice%DIVA_m_ij2n_u( j    ,i    )
    nu_in  = ice%DIVA_m_ij2n_u( j_in ,i_in )
    nu_opp = ice%DIVA_m_ij2n_u( j_opp,i_opp)

    ! Add entries to sparse matrix lists
    IF     (BC == 'infinite') THEN
      k = k+1
      ice%DIVA_m%A_index( k) = nu
      ice%DIVA_m%A_val(   k) = 1._dp
      k = k+1
      ice%DIVA_m%A_index( k) = nu_in
      ice%DIVA_m%A_val(   k) = -1._dp
      ice%DIVA_m%b( n) = 0._dp
      ice%DIVA_m%x( n) = u( j,i)
    ELSEIF (BC == 'zero') THEN
      k = k+1
      ice%DIVA_m%A_index( k) = nu
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = 0._dp
      ice%DIVA_m%x( n) = 0._dp
    ELSEIF (BC == 'periodic') THEN
      k = k+1
      ice%DIVA_m%A_index( k) = nu
      ice%DIVA_m%A_val(   k) = 1._dp
      k = k+1
      ice%DIVA_m%A_index( k) = nu_opp
      ice%DIVA_m%A_val(   k) = -1._dp
      ice%DIVA_m%b( n) = 0._dp
      ice%DIVA_m%x( n) = u( j,i)
    ELSE
      CALL crash('list_DIVA_matrix_coefficients_eq_1_boundary: unknown BC "' // TRIM(BC) // '"')
    END IF

  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_boundary
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_boundary(      grid, ice, v, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists
    ! Exception for boundary conditions for the domain boundary

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n

    ! Local variables:
    INTEGER                                            :: i,j,i_in,i_opp,j_in,j_opp,nu,nu_in,nu_opp
    CHARACTER(LEN=256)                                 :: BC

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,3)
    j = ice%DIVA_m_n2ij_uv( n,4)

    ! Determine what kind of boundary conditions to apply
    IF (i == 1 .AND. j == 1) THEN
      ! Southwest corner
      IF     (C%DIVA_boundary_BC_v_south == 'zero' .OR. C%DIVA_boundary_BC_v_west == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_v_south == 'infinite' .OR. C%DIVA_boundary_BC_v_west == 'infinite') THEN
        BC = 'infinite'
      ELSE
        BC = 'periodic'
      END IF
    ELSEIF (i == grid%nx .AND. j == 1) THEN
      ! Southeast corner
      IF     (C%DIVA_boundary_BC_v_south == 'zero' .OR. C%DIVA_boundary_BC_v_east == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_v_south == 'infinite' .OR. C%DIVA_boundary_BC_v_east == 'infinite') THEN
        BC = 'infinite'
      ELSE
        BC = 'periodic'
      END IF
    ELSEIF (i == 1 .AND. j == grid%ny-1) THEN
      ! Northwest corner
      IF     (C%DIVA_boundary_BC_v_north == 'zero' .OR. C%DIVA_boundary_BC_v_west == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_v_north == 'infinite' .OR. C%DIVA_boundary_BC_v_west == 'infinite') THEN
        BC = 'infinite'
      ELSE
        BC = 'periodic'
      END IF
    ELSEIF (i == grid%nx .AND. j == grid%ny-1) THEN
      ! Northeast corner
      IF     (C%DIVA_boundary_BC_v_north == 'zero' .OR. C%DIVA_boundary_BC_v_east == 'zero') THEN
        BC = 'zero'
      ELSEIF (C%DIVA_boundary_BC_v_north == 'infinite' .OR. C%DIVA_boundary_BC_v_east == 'infinite') THEN
        BC = 'infinite'
      ELSE
        BC = 'periodic'
      END IF
    ELSEIF (i == 1) THEN
      ! Western boundary
      BC = C%DIVA_boundary_BC_v_west
    ELSEIF (i == grid%nx) THEN
      ! Eastern boundary
      BC = C%DIVA_boundary_BC_v_east
    ELSEIF (j == 1) THEN
      ! Southern boundary
      BC = C%DIVA_boundary_BC_v_south
    ELSEIF (j == grid%ny-1) THEN
      ! Northern boundary
      BC = C%DIVA_boundary_BC_v_north
    ELSE
      CALL crash('list_DIVA_matrix_coefficients_eq_2_boundary: grid cell [{int_01},{int_02}] doesnt lie on the domain boundary!', int_01 = j, int_02 = i)
    END IF

    ! Find indices of relevant grid cells
    IF     (i == 1) THEN
      i_in  = 2
      i_opp = grid%nx-1
    ELSEIF (i == grid%nx) THEN
      i_in  = grid%nx-1
      i_opp = 2
    ELSE
      i_in  = i
      i_opp = i
    END IF
    IF     (j == 1) THEN
      j_in  = 2
      j_opp = grid%ny-2
    ELSEIF (j == grid%ny-1) THEN
      j_in  = grid%ny-2
      j_opp = 2
    ELSE
      j_in  = j
      j_opp = j
    END IF

    nu     = ice%DIVA_m_ij2n_v( j    ,i    )
    nu_in  = ice%DIVA_m_ij2n_v( j_in ,i_in )
    nu_opp = ice%DIVA_m_ij2n_v( j_opp,i_opp)

    ! Add entries to sparse matrix lists
    IF     (BC == 'infinite') THEN
      k = k+1
      ice%DIVA_m%A_index( k) = nu
      ice%DIVA_m%A_val(   k) = 1._dp
      k = k+1
      ice%DIVA_m%A_index( k) = nu_in
      ice%DIVA_m%A_val(   k) = -1._dp
      ice%DIVA_m%b( n) = 0._dp
      ice%DIVA_m%x( n) = v( j,i)
    ELSEIF (BC == 'zero') THEN
      k = k+1
      ice%DIVA_m%A_index( k) = nu
      ice%DIVA_m%A_val(   k) = 1._dp
      ice%DIVA_m%b( n) = 0._dp
      ice%DIVA_m%x( n) = 0._dp
    ELSEIF (BC == 'periodic') THEN
      k = k+1
      ice%DIVA_m%A_index( k) = nu
      ice%DIVA_m%A_val(   k) = 1._dp
      k = k+1
      ice%DIVA_m%A_index( k) = nu_opp
      ice%DIVA_m%A_val(   k) = -1._dp
      ice%DIVA_m%b( n) = 0._dp
      ice%DIVA_m%x( n) = v( j,i)
    ELSE
      CALL crash('list_DIVA_matrix_coefficients_eq_2_boundary: unknown BC "' // TRIM(BC) // '"')
    END IF

  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_boundary
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_margin(        grid, ice, u, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u
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
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_margin(        grid, ice, v, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v
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
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_free(          grid, ice, u, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n

    ! Local variables:
    INTEGER                                            :: i,j

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,1)
    j = ice%DIVA_m_n2ij_uv( n,2)

    ! Matrix coefficients
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

    ! Right-hand side and initial guess
    ice%DIVA_m%b( n) = -ice%taudx_cx( j,i) * grid%dx**2
    ice%DIVA_m%x( n) = u( j,i)

  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_free
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_free(          grid, ice, v, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n

    ! Local variables:
    INTEGER                                            :: i,j

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,3)
    j = ice%DIVA_m_n2ij_uv( n,4)

    ! Matrix coefficients
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

    ! Right-hand side and initial guess
    ice%DIVA_m%b( n) = -ice%taudy_cy( j,i) * grid%dx**2
    ice%DIVA_m%x( n) = v( j,i)

  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_free
  SUBROUTINE list_DIVA_matrix_coefficients_eq_1_free_sans(     grid, ice, u, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n

    ! Local variables:
    INTEGER                                            :: i,j

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,1)
    j = ice%DIVA_m_n2ij_uv( n,2)

    ! Matrix coefficients
    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j,i+1)
    ice%DIVA_m%A_val(   k) = 4._dp

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j,i-1)
    ice%DIVA_m%A_val(   k) = 4._dp

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j+1,i)
    ice%DIVA_m%A_val(   k) = 1._dp

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j-1,i)
    ice%DIVA_m%A_val(   k) = 1._dp

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j,i)
    ice%DIVA_m%A_val(   k) = -10._dp - grid%dx**2 * ice%beta_eff_cx( j,i) / ice%N_cx( j,i)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j,i)
    ice%DIVA_m%A_val(   k) = -3._dp

    k  = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j,i+1)
    ice%DIVA_m%A_val(   k) =  3._dp

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j-1,i)
    ice%DIVA_m%A_val(   k) =  3._dp

    k  = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j-1,i+1)
    ice%DIVA_m%A_val(   k) = -3._dp

    ! Right-hand side and initial guess
    ice%DIVA_m%b( n) = -ice%taudx_cx( j,i) * grid%dx**2 / ice%N_cx( j,i)
    ice%DIVA_m%x( n) = u( j,i)

  END SUBROUTINE list_DIVA_matrix_coefficients_eq_1_free_sans
  SUBROUTINE list_DIVA_matrix_coefficients_eq_2_free_sans(     grid, ice, v, k, n)
    ! Find the matrix coefficients of equation n and add them to the lists

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v
    INTEGER,                             INTENT(INOUT) :: k
    INTEGER,                             INTENT(IN)    :: n

    ! Local variables:
    INTEGER                                            :: i,j

    ! Grid indices of the grid cell represented by equation n
    i = ice%DIVA_m_n2ij_uv( n,3)
    j = ice%DIVA_m_n2ij_uv( n,4)

    ! Matrix coefficients
    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j+1,i)
    ice%DIVA_m%A_val(   k) = 4._dp

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j-1,i)
    ice%DIVA_m%A_val(   k) = 4._dp

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j,i+1)
    ice%DIVA_m%A_val(   k) = 1._dp

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j,i-1)
    ice%DIVA_m%A_val(   k) = 1._dp

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_v( j,i)
    ice%DIVA_m%A_val(   k) = -10._dp - grid%dx**2 * ice%beta_eff_cy( j,i) / ice%N_cy( j,i)

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j+1,i-1)
    ice%DIVA_m%A_val(   k) = -3._dp

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j+1,i)
    ice%DIVA_m%A_val(   k) =  3._dp

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j,i-1)
    ice%DIVA_m%A_val(   k) =  3._dp

    k = k+1
    ice%DIVA_m%A_index( k) = ice%DIVA_m_ij2n_u( j,i)
    ice%DIVA_m%A_val(   k) = -3._dp

    ! Right-hand side and initial guess
    ice%DIVA_m%b( n) = -ice%taudy_cy( j,i) * grid%dx**2 / ice%N_cy( j,i)
    ice%DIVA_m%x( n) = v( j,i)

  END SUBROUTINE list_DIVA_matrix_coefficients_eq_2_free_sans

  ! Map velocity components to the a-grid for writing to output (diagnostic only)
  SUBROUTINE map_velocity_from_cx_to_a_2D( grid, ice, u_cx, u_a)
    ! Map velocity components from the staggered cx/cy-grids to the regular a-grid.
    ! When both neighbours of the a-grid cell are ice-covered, take the average
    ! of the two staggered velocities; when a neighbour is ice-free, take only the
    ! velocity on the ice-covered neighbour side (to prevent spurious velocity
    ! drops at the ice margin).

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_cx
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: u_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocity_from_cx_to_a_2D'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocity_from_cx_to_a_2D
  SUBROUTINE map_velocity_from_cx_to_a_3D( grid, ice, u_cx, u_a)
    ! Map velocity components from the staggered cx/cy-grids to the regular a-grid.
    ! When both neighbours of the a-grid cell are ice-covered, take the average
    ! of the two staggered velocities; when a neighbour is ice-free, take only the
    ! velocity on the ice-covered neighbour side (to prevent spurious velocity
    ! drops at the ice margin).

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: u_cx
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: u_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocity_from_cx_to_a_3D'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocity_from_cx_to_a_3D
  SUBROUTINE map_velocity_from_cy_to_a_2D( grid, ice, v_cy, v_a)
    ! Map velocity components from the staggered cx/cy-grids to the regular a-grid.
    ! When both neighbours of the a-grid cell are ice-covered, take the average
    ! of the two staggered velocities; when a neighbour is ice-free, take only the
    ! velocity on the ice-covered neighbour side (to prevent spurious velocity
    ! drops at the ice margin).

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_cy
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocity_from_cy_to_a_2D'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocity_from_cy_to_a_2D
  SUBROUTINE map_velocity_from_cy_to_a_3D( grid, ice, v_cy, v_a)
    ! Map velocity components from the staggered cx/cy-grids to the regular a-grid.
    ! When both neighbours of the a-grid cell are ice-covered, take the average
    ! of the two staggered velocities; when a neighbour is ice-free, take only the
    ! velocity on the ice-covered neighbour side (to prevent spurious velocity
    ! drops at the ice margin).

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: v_cy
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: v_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocity_from_cy_to_a_3D'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocity_from_cy_to_a_3D

  ! Allocate and initialise the matrix-based SSA/DIVA solver
  SUBROUTINE initialise_SSADIVA_solution_matrix( grid, ice)
    ! Allocate and initialise the translation tables and sparse matrix lists used
    ! in the matrix-based SSA/DIVA solver.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_SSADIVA_solution_matrix'
    INTEGER                                            :: i,j,neq,nnz_per_row_max,n

    ! Add routine to path
    CALL init_routine( routine_name)

    neq             = (grid%nx - 1) * grid%ny + grid%nx * (grid%ny - 1)
    nnz_per_row_max = 9

    CALL initialise_matrix_equation_CSR( ice%DIVA_m, neq, neq, nnz_per_row_max)

    CALL allocate_shared_int_2D( grid%ny  , grid%nx-1, ice%DIVA_m_ij2n_u , ice%wDIVA_m_ij2n_u )
    CALL allocate_shared_int_2D( grid%ny-1, grid%nx  , ice%DIVA_m_ij2n_v , ice%wDIVA_m_ij2n_v )
    CALL allocate_shared_int_2D( neq,       4        , ice%DIVA_m_n2ij_uv, ice%wDIVA_m_n2ij_uv)

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SSADIVA_solution_matrix
  SUBROUTINE initialise_ice_velocity_ISMIP_HOM( grid, ice)
    ! Initialise ice velocity fields so that the ISMIP-HOM experiments converge faster

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_velocity_ISMIP_HOM'
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  u_ISMIP_HOM
    INTEGER                                            :: wu_ISMIP_HOM
    INTEGER                                            :: i,j
    REAL(dp)                                           :: x, y, umin, umax

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, u_ISMIP_HOM, wu_ISMIP_HOM)

    umin = 0._dp
    umax = 0._dp

    ! Calculate an approximation of the solution
    IF (C%choice_refgeo_init_idealised == 'ISMIP_HOM_A') THEN

      IF     (C%ISMIP_HOM_L == 160000._dp) THEN
        umin = 1.6_dp
        umax = 108.84_dp
      ELSEIF (C%ISMIP_HOM_L ==  80000._dp) THEN
        umin = 1.75_dp
        umax = 95.73_dp
      ELSEIF (C%ISMIP_HOM_L ==  40000._dp) THEN
        umin = 2.27_dp
        umax = 74.45_dp
      ELSEIF (C%ISMIP_HOM_L ==  20000._dp) THEN
        umin = 4.49_dp
        umax = 49.99_dp
      ELSEIF (C%ISMIP_HOM_L ==  10000._dp) THEN
        umin = 11.09_dp
        umax = 32.74_dp
      ELSEIF (C%ISMIP_HOM_L ==   5000._dp) THEN
        umin = 18.38_dp
        umax = 24.79_dp
      END IF

      DO i = grid%i1, MIN( grid%nx-1, grid%i2)
      DO j = 1, grid%ny
        x = (grid%x( i) + grid%x( i+1)) / 2._dp
        y = grid%y( j)
        u_ISMIP_HOM( j,i) = umin + (umax - umin) * ( (1._dp - (SIN( x) * SIN( y))) / 2._dp)**2
      END DO
      END DO
      CALL sync

    ELSEIF (C%choice_refgeo_init_idealised == 'ISMIP_HOM_B') THEN

      IF     (C%ISMIP_HOM_L == 160000._dp) THEN
        umin = 1.57_dp
        umax = 111.41_dp
      ELSEIF (C%ISMIP_HOM_L ==  80000._dp) THEN
        umin = 1.69_dp
        umax = 100.73_dp
      ELSEIF (C%ISMIP_HOM_L ==  40000._dp) THEN
        umin = 2.09_dp
        umax = 82.3_dp
      ELSEIF (C%ISMIP_HOM_L ==  20000._dp) THEN
        umin = 3.92_dp
        umax = 57.84_dp
      ELSEIF (C%ISMIP_HOM_L ==  10000._dp) THEN
        umin = 10.23_dp
        umax = 35.2_dp
      ELSEIF (C%ISMIP_HOM_L ==   5000._dp) THEN
        umin = 17.22_dp
        umax = 23.53_dp
      END IF

      DO i = grid%i1, MIN( grid%nx-1, grid%i2)
      DO j = 1, grid%ny
        x = (grid%x( i) + grid%x( i+1)) / 2._dp
        y = grid%y( j)
        u_ISMIP_HOM( j,i) = umin + (umax - umin) * ( (1._dp - (SIN( x) * SIN( y))) / 2._dp)**2
      END DO
      END DO
      CALL sync

    ELSEIF (C%choice_refgeo_init_idealised == 'ISMIP_HOM_C') THEN

      IF     (C%ISMIP_HOM_L == 160000._dp) THEN
        umin = 8.77_dp
        umax = 143.45_dp
      ELSEIF (C%ISMIP_HOM_L ==  80000._dp) THEN
        umin = 9.8_dp
        umax = 60.28_dp
      ELSEIF (C%ISMIP_HOM_L ==  40000._dp) THEN
        umin = 11.84_dp
        umax = 28.57_dp
      ELSEIF (C%ISMIP_HOM_L ==  20000._dp) THEN
        umin = 14.55_dp
        umax = 18.48_dp
      ELSEIF (C%ISMIP_HOM_L ==  10000._dp) THEN
        umin = 15.7_dp
        umax = 16.06_dp
      ELSEIF (C%ISMIP_HOM_L ==   5000._dp) THEN
        umin = 13.38_dp
        umax = 13.51_dp
      END IF

      DO i = grid%i1, MIN( grid%nx-1, grid%i2)
      DO j = 1, grid%ny
        x = (grid%x( i) + grid%x( i+1)) / 2._dp
        y = grid%y( j)
        u_ISMIP_HOM( j,i) = umin + (umax - umin) * ( (1._dp - (SIN( x) * SIN( y))) / 2._dp)**2
      END DO
      END DO
      CALL sync

    ELSEIF (C%choice_refgeo_init_idealised == 'ISMIP_HOM_D') THEN

      IF     (C%ISMIP_HOM_L == 160000._dp) THEN
        umin = 8.62_dp
        umax = 227.23_dp
      ELSEIF (C%ISMIP_HOM_L ==  80000._dp) THEN
        umin = 9.65_dp
        umax = 94.79_dp
      ELSEIF (C%ISMIP_HOM_L ==  40000._dp) THEN
        umin = 12.18_dp
        umax = 40.06_dp
      ELSEIF (C%ISMIP_HOM_L ==  20000._dp) THEN
        umin = 15.28_dp
        umax = 20.29_dp
      ELSEIF (C%ISMIP_HOM_L ==  10000._dp) THEN
        umin = 15.93_dp
        umax = 16.25_dp
      ELSEIF (C%ISMIP_HOM_L ==   5000._dp) THEN
        umin = 14.43_dp
        umax = 14.59_dp
      END IF

      DO i = grid%i1, MIN( grid%nx-1, grid%i2)
      DO j = 1, grid%ny
        x = (grid%x( i) + grid%x( i+1)) / 2._dp
        y = grid%y( j)
        u_ISMIP_HOM( j,i) = umin + (umax - umin) * ( (1._dp - (SIN( x) * SIN( y))) / 2._dp)**2
      END DO
      END DO
      CALL sync

    END IF

    ! Initialise velocity fields with the approximation
    IF     (C%choice_ice_dynamics == 'SIA/SSA') THEN
      ice%u_SSA_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = u_ISMIP_HOM( :,grid%i1:MIN(grid%nx-1,grid%i2))
      CALL sync
    ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
      ice%u_vav_cx( :,grid%i1:MIN(grid%nx-1,grid%i2)) = u_ISMIP_HOM( :,grid%i1:MIN(grid%nx-1,grid%i2))
      CALL sync
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wu_ISMIP_HOM)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_velocity_ISMIP_HOM

  SUBROUTINE initialise_velocities_from_restart_file( grid, ice, region_name)
    ! Initialise velocities with data from a previous simulation's restart file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_ice_model),           INTENT(INOUT) :: ice
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_velocities_from_restart_file'
    CHARACTER(LEN=256)                            :: filename_restart
    REAL(dp)                                      :: time_to_restart_from
    REAL(dp), DIMENSION( :,:), POINTER            :: u_SSA_cx_a, v_SSA_cy_a, u_vav_cx_a, v_vav_cy_a
    INTEGER                                       :: wu_SSA_cx_a, wv_SSA_cy_a, wu_vav_cx_a, wv_vav_cy_a

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate room for velocity fields
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, u_SSA_cx_a, wu_SSA_cx_a)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, v_SSA_cy_a, wv_SSA_cy_a)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, u_vav_cx_a, wu_vav_cx_a)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, v_vav_cy_a, wv_vav_cy_a)

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

    ! Check whether SIA/SSA or DIVA is used
    IF (C%choice_ice_dynamics == 'SIA/SSA' .OR. C%choice_ice_dynamics == 'SSA') THEN
      u_SSA_cx_a = 0._dp
      v_SSA_cy_a = 0._dp
      CALL read_field_from_file_2D(   filename_restart, 'u_SSA_cx_a', grid,  u_SSA_cx_a,  region_name, time_to_restart_from)
      CALL read_field_from_file_2D(   filename_restart, 'v_SSA_cy_a', grid,  v_SSA_cy_a,  region_name, time_to_restart_from)
      IF (par%master) THEN
        ice%u_SSA_cx = u_SSA_cx_a(:, 1:grid%nx-1)
        ice%v_SSA_cy = v_SSA_cy_a(1:grid%ny-1, :)
        WRITE(0,*) '  Initialising velocities from restart file...'
      END IF
      CALL sync
    END IF

    IF (C%choice_timestepping == 'pc' .OR. C%choice_ice_dynamics == 'DIVA') THEN

    ! ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
      u_vav_cx_a = 0._dp
      v_vav_cy_a = 0._dp
      CALL read_field_from_file_2D(   filename_restart, 'u_vav_cx_a', grid,  u_vav_cx_a,  region_name, time_to_restart_from)
      CALL read_field_from_file_2D(   filename_restart, 'v_vav_cy_a', grid,  v_vav_cy_a,  region_name, time_to_restart_from)
      IF (par%master) THEN
        ice%u_vav_cx = u_vav_cx_a(:, 1:grid%nx-1)
        ice%v_vav_cy = v_vav_cy_a(1:grid%ny-1, :)
        WRITE(0,*) '  Initialising velocities from restart file...'
      END IF
      CALL sync
    END IF

    ! Safety
    CALL check_for_NaN_dp_2D( ice%u_SSA_cx, 'ice%wu_SSA_cx')
    CALL check_for_NaN_dp_2D( ice%v_SSA_cy, 'ice%wv_SSA_cy')
    CALL check_for_NaN_dp_2D( ice%u_vav_cx, 'ice%wu_vav_cx')
    CALL check_for_NaN_dp_2D( ice%v_vav_cy, 'ice%wv_vav_cy')

    CALL deallocate_shared( wu_SSA_cx_a)
    CALL deallocate_shared( wv_SSA_cy_a)
    CALL deallocate_shared( wu_vav_cx_a)
    CALL deallocate_shared( wv_vav_cy_a)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_velocities_from_restart_file

END MODULE ice_velocity_module
