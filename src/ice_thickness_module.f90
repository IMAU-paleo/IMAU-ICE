MODULE ice_thickness_module

  ! Contains the routines for solving the ice thickness equation

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_SMB_model, type_BMB_model
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             initialise_matrix_equation_CSR, solve_matrix_equation_CSR, check_CSR_for_double_entries

  IMPLICIT NONE
  
CONTAINS

  ! The main routine that is called from "run_ice_model" in the ice_dynamics_module
  SUBROUTINE calc_dHi_dt( grid, ice, SMB, BMB, dt, mask_noice)
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
              C%choice_benchmark_experiment == 'ISMIP_HOM_F' .OR. &
              C%choice_benchmark_experiment == 'MISMIPplus') THEN
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
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in calc_dHi_dt!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Use the specified time integration method to calculate the ice thickness at t+dt
    IF     (C%choice_ice_integration_method == 'explicit') THEN
      CALL calc_dHi_dt_explicit(     grid, ice, SMB, BMB, dt)
    ELSEIF (C%choice_ice_integration_method == 'semi-implicit') THEN
      CALL calc_dHi_dt_semiimplicit( grid, ice, SMB, BMB, dt)
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_ice_integration_method "', TRIM(C%choice_ice_integration_method), '" not implemented in calc_dHi_dt!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Remove ice in areas where no ice is allowed (i.e. Greenland in NAM and EAS, and Ellesmere Island in GRL)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (mask_noice(     j,i) == 1) THEN
        ice%dHi_dt_a(     j,i) = -ice%Hi_a( j,i) / dt
        ice%Hi_tplusdt_a( j,i) = 0._dp
      END IF
    END DO
    END DO
    CALL sync
    
    ! Apply boundary conditions
    CALL apply_ice_thickness_BC( grid, ice, dt)
    
    ! Remove free-floating shelves not connected to any grounded ice
    CALL remove_unconnected_shelves( grid, ice, dt)
    
  END SUBROUTINE calc_dHi_dt
  
  ! Different solvers for the ice thickness equation (explicit & semi-implicit)
  SUBROUTINE calc_dHi_dt_explicit(     grid, ice, SMB, BMB, dt)
    ! The explicit solver for the ice thickness equation
    
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
      ! => And for ice-free ocean, where no accumulation is allowed
      IF (ice%mask_cf_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 1) THEN
        dVi_MB( j,i) = (SMB%SMB_year( j,i) + BMB%BMB( j,i)) * grid%dx * grid%dx * dt * ice%float_margin_frac_a( j,i)
      ELSEIF (ice%mask_ocean_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 0) THEN
        dVi_MB( j,i) = 0._dp
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
      
      ! Rate of change
      ice%dHi_dt_a( j,i) = dVi_MB( j,i) / (grid%dx * grid%dx * dt) ! m/y
      IF (i > 1      ) ice%dHi_dt_a( j,i) = ice%dHi_dt_a( j,i) + ice%Qx_cx( j  ,i-1) / (grid%dx * grid%dx * dt)
      IF (i < grid%nx) ice%dHi_dt_a( j,i) = ice%dHi_dt_a( j,i) - ice%Qx_cx( j  ,i  ) / (grid%dx * grid%dx * dt)
      IF (j > 1      ) ice%dHi_dt_a( j,i) = ice%dHi_dt_a( j,i) + ice%Qy_cy( j-1,i  ) / (grid%dx * grid%dx * dt)
      IF (j < grid%ny) ice%dHi_dt_a( j,i) = ice%dHi_dt_a( j,i) - ice%Qy_cy( j  ,i  ) / (grid%dx * grid%dx * dt)
      
      ! Don't allow negative ice thickness
      ice%dHi_dt_a( j,i) = MAX( ice%dHi_dt_a( j,i), -ice%Hi_a( j,i) / dt)
      
      ! Thickness at t+dt
      ice%Hi_tplusdt_a( j,i) = ice%Hi_a( j,i) + ice%dHi_dt_a( j,i) * dt
      
    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wdVi_MB)
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%Qx_cx   , 'ice%Qx_cx'   , 'calc_dHi_dt_explicit')
    CALL check_for_NaN_dp_2D( ice%Qy_cy   , 'ice%Qx_cx'   , 'calc_dHi_dt_explicit')
    CALL check_for_NaN_dp_2D( ice%dHi_dt_a, 'ice%dHi_dt_a', 'calc_dHi_dt_explicit')
    
  END SUBROUTINE calc_dHi_dt_explicit
  SUBROUTINE calc_dHi_dt_semiimplicit( grid, ice, SMB, BMB, dt)
    ! Calculate the ice thickness at time t+dt using an implicit matrix-based method
    ! (uses an explicit scheme around the calving front to improve stability)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    REAL(dp),                            INTENT(IN)    :: dt
    
    ! Local variables
    INTEGER                                            :: i,j,n,k,im1,ip1,jm1,jp1
    REAL(dp)                                           :: MB_net, Hi_tplusdt
    CHARACTER(LEN=256)                                 :: local_scheme_choice
    REAL(dp), PARAMETER                                :: fs = 2.5_dp           ! Semi-implicit scale factor (0 = explicit, 1 = implicit, 1/2 = Crank-Nicolson, >1 = over-implicit)
    
    ! Fill the sparse matrix
    ! ======================
    
    IF (par%master) THEN
      
      ice%dHi_m%A_ptr   = 0
      ice%dHi_m%A_index = 0
      ice%dHi_m%A_val   = 0._dp
      ice%dHi_m%b       = 0._dp
      ice%dHi_m%x       = 0._dp
      
      k = 0
      ice%dHi_m%A_ptr(1) = 1
      
      DO i = 1, grid%nx
      DO j = 1, grid%ny
        
        ! Equation number
        n = ice%dHi_ij2n( j,i)
      
      ! Determine which scheme to use locally
      ! (explicit around the calving front, over-implicit everywhere else)
      ! ==================================================================
        
        im1 = MAX( 1,      i-1)
        ip1 = MIN( grid%nx,i+1)
        jm1 = MAX( 1,      j-1)
        jp1 = MIN( grid%ny,j+1)
        
        IF ( ice%mask_ice_a( j  ,i  ) == 0 .AND. &
             ice%mask_ice_a( j  ,im1) == 0 .AND. &
             ice%mask_ice_a( j  ,ip1) == 0 .AND. &
             ice%mask_ice_a( jm1,i  ) == 0 .AND. &
             ice%mask_ice_a( jp1,i  ) == 0) THEN
          ! No ice in the local neighbourhood, only use net mass balance
          
          local_scheme_choice = 'MB-only'
          
        ELSEIF (ice%mask_shelf_a( j,i) == 1 .AND. &
               ((ice%mask_ocean_a( j  ,im1) == 1 .AND. ice%mask_ice_a( j  ,im1) == 0) .OR. &
                (ice%mask_ocean_a( j  ,ip1) == 1 .AND. ice%mask_ice_a( j  ,ip1) == 0) .OR. &
                (ice%mask_ocean_a( jm1,i  ) == 1 .AND. ice%mask_ice_a( jm1,i  ) == 0) .OR. &
                (ice%mask_ocean_a( jp1,i  ) == 1 .AND. ice%mask_ice_a( jp1,i  ) == 0))) THEN
          ! Shelf next to open ocean; use explicit scheme
          
          local_scheme_choice = 'explicit'
          
        ELSEIF (ice%mask_ocean_a( j,i) == 1 .AND. ice%mask_ice_a( j,i) == 0 .AND. &
               (ice%mask_shelf_a( j  ,im1) == 1 .OR. &
                ice%mask_shelf_a( j  ,ip1) == 1 .OR. &
                ice%mask_shelf_a( jm1,i  ) == 1 .OR. &
                ice%mask_shelf_a( jp1,i  ) == 1)) THEN
          ! Open ocean next to shelf; use explicit scheme
          
          local_scheme_choice = 'explicit'
        
        ELSE
          
          local_scheme_choice = 'semi-implicit'
          
        END IF
        
      ! Calculate the net local mass balance
      ! => With an exception for the calving front, where we only apply
      !    the mass balance to the floating fraction
      ! => And for ice-free ocean, where no accumulation is allowed
      ! ==================================================================
      
        IF (ice%mask_cf_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 1) THEN
          MB_net = (SMB%SMB_year( j,i) + BMB%BMB( j,i)) * ice%float_margin_frac_a( j,i)
        ELSEIF (ice%mask_ocean_a( j,i) == 1 .AND. ice%mask_shelf_a( j,i) == 0) THEN
          MB_net = 0._dp
        ELSE
          MB_net = (SMB%SMB_year( j,i) + BMB%BMB( j,i))
        END IF
        
      ! Fill in matrix coefficients for the chosen local integration scheme
      ! ===================================================================
        
        IF     (local_scheme_choice == 'MB-only') THEN
          ! Neglect advective flow; apply only the local mass balance (useful when no ice exists in the local neighbourhood)
          
          Hi_tplusdt = MAX( 0._dp, ice%Hi_a( j,i) + (MB_net * dt))
          
          k = k+1
          ice%dHi_m%A_index( k) = ice%dHi_ij2n( j,i)
          ice%dHi_m%A_val(   k) = 1._dp
          ice%dHi_m%b(       n) = Hi_tplusdt
          ice%dHi_m%x(       n) = Hi_tplusdt
          
        ELSEIF (local_scheme_choice == 'explicit') THEN
          ! Explicit scheme
          
          CALL calc_dHi_dt_semiimplicit_add_matrix_coefficients_explicit(     grid, ice, MB_net, dt, i,j, k)
          
        ELSEIF (local_scheme_choice == 'semi-implicit') THEN
          ! Semi-implicit scheme
          
          CALL calc_dHi_dt_semiimplicit_add_matrix_coefficients_semiimplicit( grid, ice, MB_net, dt, i,j, k, fs)
        
        END IF ! IF (local_scheme_choice == 'MB-only') THEN
        
        ! Cycle pointer
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
    
    ! Map the results back from vector form to the model grid
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      n = ice%dHi_ij2n( j,i)
      ice%Hi_tplusdt_a( j,i) = MAX( 0._dp, ice%dHi_m%x( n))
      ice%dHi_dt_a( j,i) = (ice%Hi_tplusdt_a( j,i) - ice%Hi_a( j,i)) / dt
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%dHi_dt_a, 'ice%dHi_dt_a', 'calc_dHi_dt_semiimplicit')
    
  END SUBROUTINE calc_dHi_dt_semiimplicit
  SUBROUTINE calc_dHi_dt_semiimplicit_add_matrix_coefficients_explicit(     grid, ice, MB_net, dt, i,j, k)
    ! Add matrix coefficients representing the explicit scheme
    ! (used around the calving front in the matrix-based semi-implicit scheme)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: MB_net
    REAL(dp),                            INTENT(IN)    :: dt
    INTEGER,                             INTENT(IN)    :: i,j
    INTEGER,                             INTENT(INOUT) :: k
    
    ! Local variables
    INTEGER                                            :: n
    REAL(dp)                                           :: Qw, Qe, Qs, Qn, Hi_tplusdt
        
    ! Equation number
    n = ice%dHi_ij2n( j,i)
      
    ! Calculate ice flux through all four cell faces
    ! ==============================================
    
    IF (i > 1) THEN
      IF (ice%u_vav_cx( j  ,i-1) > 0._dp) THEN
        Qw =  ice%u_vav_cx( j  ,i-1) * ice%Hi_a( j  ,i-1) / grid%dx
      ELSE
        Qw =  ice%u_vav_cx( j  ,i-1) * ice%Hi_a( j  ,i  ) / grid%dx
      END IF
    ELSE
      Qw = 0._dp
    END IF
    
    IF (i < grid%nx) THEN
      IF (ice%u_vav_cx( j  ,i  ) > 0._dp) THEN
        Qe = -ice%u_vav_cx( j  ,i  ) * ice%Hi_a( j  ,i  ) / grid%dx
      ELSE
        Qe = -ice%u_vav_cx( j  ,i  ) * ice%Hi_a( j  ,i+1) / grid%dx
      END IF
    ELSE
      Qe = 0._dp
    END IF
    
    IF (j > 1) THEN
      IF (ice%v_vav_cy( j-1,i  ) > 0._dp) THEN
        Qs =  ice%v_vav_cy( j-1,i  ) * ice%Hi_a( j-1,i  ) / grid%dx
      ELSE
        Qs =  ice%v_vav_cy( j-1,i  ) * ice%Hi_a( j  ,i  ) / grid%dx
      END IF
    ELSE
      Qs = 0._dp
    END IF
    
    IF (j < grid%ny) THEN
      IF (ice%v_vav_cy( j  ,i  ) > 0._dp) THEN
        Qn = -ice%v_vav_cy( j  ,i  ) * ice%Hi_a( j  ,i  ) / grid%dx
      ELSE
        Qn = -ice%v_vav_cy( j  ,i  ) * ice%Hi_a( j+1,i  ) / grid%dx
      END IF
    ELSE
      Qn = 0._dp
    END IF

    ! Flow from floating ice to open ocean is only
    ! allowed when the floating pixel is completely filled
    ! ====================================================
      
    IF (i > 1) THEN
      IF     (ice%mask_shelf_a( j  ,i  ) == 1 .AND. ice%mask_ocean_a( j  ,i-1) == 1 .AND. ice%mask_ice_a( j  ,i-1) == 0) THEN
        ! Shelf here, open ocean to the west
        IF (ice%float_margin_frac_a( j  ,i  ) < 1._dp) Qw = 0._dp
      ELSEIF (ice%mask_shelf_a( j  ,i-1) == 1 .AND. ice%mask_ocean_a( j  ,i  ) == 1 .AND. ice%mask_ice_a( j  ,i  ) == 0) THEN
        ! Open ocean here, shelf to the west
        IF (ice%float_margin_frac_a( j  ,i-1) < 1._dp) Qw = 0._dp
      END IF
    END IF
    IF (i < grid%nx) THEN
      IF     (ice%mask_shelf_a( j  ,i  ) == 1 .AND. ice%mask_ocean_a( j  ,i+1) == 1 .AND. ice%mask_ice_a( j  ,i+1) == 0) THEN
        ! Shelf here, open ocean to the east
        IF (ice%float_margin_frac_a( j  ,i  ) < 1._dp) Qe = 0._dp
      ELSEIF (ice%mask_shelf_a( j  ,i+1) == 1 .AND. ice%mask_ocean_a( j  ,i  ) == 1 .AND. ice%mask_ice_a( j  ,i  ) == 0) THEN
        ! Open ocean here, shelf to the east
        IF (ice%float_margin_frac_a( j  ,i+1) < 1._dp) Qe = 0._dp
      END IF
    END IF
    IF (j > 1) THEN
      IF     (ice%mask_shelf_a( j  ,i  ) == 1 .AND. ice%mask_ocean_a( j-1,i  ) == 1 .AND. ice%mask_ice_a( j-1,i  ) == 0) THEN
        ! Shelf here, open ocean to the south
        IF (ice%float_margin_frac_a( j  ,i  ) < 1._dp) Qs = 0._dp
      ELSEIF (ice%mask_shelf_a( j-1,i  ) == 1 .AND. ice%mask_ocean_a( j  ,i  ) == 1 .AND. ice%mask_ice_a( j  ,i  ) == 0) THEN
        ! Open ocean here, shelf to the south
        IF (ice%float_margin_frac_a( j-1,i  ) < 1._dp) Qs = 0._dp
      END IF
    END IF
    IF (j < grid%ny) THEN
      IF     (ice%mask_shelf_a( j  ,i  ) == 1 .AND. ice%mask_ocean_a( j+1,i  ) == 1 .AND. ice%mask_ice_a( j+1,i  ) == 0) THEN
        ! Shelf here, open ocean to the north
        IF (ice%float_margin_frac_a( j  ,i  ) < 1._dp) Qn = 0._dp
      ELSEIF (ice%mask_shelf_a( j+1,i  ) == 1 .AND. ice%mask_ocean_a( j  ,i  ) == 1 .AND. ice%mask_ice_a( j  ,i  ) == 0) THEN
        ! Open ocean here, shelf to the north
        IF (ice%float_margin_frac_a( j+1,i  ) < 1._dp) Qn = 0._dp
      END IF
    END IF
      
    ! Calculate Hi at next time step
    ! ==============================
    
    Hi_tplusdt = MAX(0._dp, ice%Hi_a( j,i) + dt * (MB_net + Qw + Qe + Qs + Qn))
      
    ! Fill in sparse matrix coefficients
    ! ==================================
    
    k = k+1
    ice%dHi_m%A_index( k) = ice%dHi_ij2n( j,i)
    ice%dHi_m%A_val(   k) = 1._dp
    ice%dHi_m%b( n) = Hi_tplusdt
    ice%dHi_m%x( n) = Hi_tplusdt
    
  END SUBROUTINE calc_dHi_dt_semiimplicit_add_matrix_coefficients_explicit
  SUBROUTINE calc_dHi_dt_semiimplicit_add_matrix_coefficients_semiimplicit( grid, ice, MB_net, dt, i,j, k, fs)
    ! Add matrix coefficients for the semi-implicit scheme
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: MB_net
    REAL(dp),                            INTENT(IN)    :: dt
    INTEGER,                             INTENT(IN)    :: i,j
    INTEGER,                             INTENT(INOUT) :: k
    REAL(dp),                            INTENT(IN)    :: fs
    
    ! Local variables
    INTEGER                                            :: n
    REAL(dp)                                           :: u_west, u_east, v_south, v_north
        
    ! Equation number
    n = ice%dHi_ij2n( j,i)
    
    ! Flow from floating ice to open ocean is only
    ! allowed when the floating pixel is completely filled
    ! ====================================================
      
    IF (i > 1) THEN
      u_west = ice%u_vav_cx( j  ,i-1)
      IF     (ice%mask_shelf_a( j  ,i  ) == 1 .AND. ice%mask_ocean_a( j  ,i-1) == 1 .AND. ice%mask_ice_a( j  ,i-1) == 0) THEN
        ! Shelf here, open ocean to the west
        IF (ice%float_margin_frac_a( j  ,i  ) < 1._dp) u_west = 0._dp
      ELSEIF (ice%mask_shelf_a( j  ,i-1) == 1 .AND. ice%mask_ocean_a( j  ,i  ) == 1 .AND. ice%mask_ice_a( j  ,i  ) == 0) THEN
        ! Open ocean here, shelf to the west
        IF (ice%float_margin_frac_a( j  ,i-1) < 1._dp) u_west = 0._dp
      END IF
    ELSE
      u_west = 0._dp
    END IF
    
    IF (i < grid%nx) THEN
      u_east = ice%u_vav_cx( j  ,i  )
      IF     (ice%mask_shelf_a( j  ,i  ) == 1 .AND. ice%mask_ocean_a( j  ,i+1) == 1 .AND. ice%mask_ice_a( j  ,i+1) == 0) THEN
        ! Shelf here, open ocean to the east
        IF (ice%float_margin_frac_a( j  ,i  ) < 1._dp) u_east = 0._dp
      ELSEIF (ice%mask_shelf_a( j  ,i+1) == 1 .AND. ice%mask_ocean_a( j  ,i  ) == 1 .AND. ice%mask_ice_a( j  ,i  ) == 0) THEN
        ! Open ocean here, shelf to the east
        IF (ice%float_margin_frac_a( j  ,i+1) < 1._dp) u_east = 0._dp
      END IF
    ELSE
      u_east = 0._dp
    END IF
    
    IF (j > 1) THEN
      v_south = ice%v_vav_cy( j-1,i  )
      IF     (ice%mask_shelf_a( j  ,i  ) == 1 .AND. ice%mask_ocean_a( j-1,i  ) == 1 .AND. ice%mask_ice_a( j-1,i  ) == 0) THEN
        ! Shelf here, open ocean to the south
        IF (ice%float_margin_frac_a( j  ,i  ) < 1._dp) v_south = 0._dp
      ELSEIF (ice%mask_shelf_a( j-1,i  ) == 1 .AND. ice%mask_ocean_a( j  ,i  ) == 1 .AND. ice%mask_ice_a( j  ,i  ) == 0) THEN
        ! Open ocean here, shelf to the south
        IF (ice%float_margin_frac_a( j-1,i  ) < 1._dp) v_south = 0._dp
      END IF
    ELSE
      v_south = 0._dp
    END IF
    
    IF (j < grid%ny) THEN
      v_north = ice%v_vav_cy( j  ,i  )
      IF     (ice%mask_shelf_a( j  ,i  ) == 1 .AND. ice%mask_ocean_a( j+1,i  ) == 1 .AND. ice%mask_ice_a( j+1,i  ) == 0) THEN
        ! Shelf here, open ocean to the north
        IF (ice%float_margin_frac_a( j  ,i  ) < 1._dp) v_north = 0._dp
      ELSEIF (ice%mask_shelf_a( j+1,i  ) == 1 .AND. ice%mask_ocean_a( j  ,i  ) == 1 .AND. ice%mask_ice_a( j  ,i  ) == 0) THEN
        ! Open ocean here, shelf to the north
        IF (ice%float_margin_frac_a( j+1,i  ) < 1._dp) v_north = 0._dp
      END IF
    ELSE
      v_north = 0._dp
    END IF
    
    ! Fill in sparse matrix coefficients
    ! ==================================
    
    ! Matrix coefficient for Hi(i,j,t+dt)
    k = k+1
    ice%dHi_m%A_index( k) = ice%dHi_ij2n( j,i)
    ice%dHi_m%A_val(   k) = 1._dp / dt
    IF (i > 1      ) ice%dHi_m%A_val( k) = ice%dHi_m%A_val( k) - (fs / grid%dx) * MIN( u_west, 0._dp)
    IF (i < grid%nx) ice%dHi_m%A_val( k) = ice%dHi_m%A_val( k) + (fs / grid%dx) * MAX( u_east, 0._dp)
    IF (j > 1      ) ice%dHi_m%A_val( k) = ice%dHi_m%A_val( k) - (fs / grid%dx) * MIN( v_south,0._dp)
    IF (j < grid%ny) ice%dHi_m%A_val( k) = ice%dHi_m%A_val( k) + (fs / grid%dx) * MAX( v_north,0._dp)
    
    ! Matrix coefficient for Hi( i-1,j  ,t+dt)
    IF (i > 1) THEN
      k = k+1
      ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i-1)
      ice%dHi_m%A_val(   k) = (-fs / grid%dx) * MAX( u_west,0._dp)
    END IF
    
    ! Matrix coefficient for Hi( i+1,j  ,t+dt)
    IF (i < grid%nx) THEN
      k = k+1
      ice%dHi_m%A_index( k) = ice%dHi_ij2n( j  ,i+1)
      ice%dHi_m%A_val(   k) = ( fs / grid%dx) * MIN( u_east,0._dp)
    END IF
      
    ! Matrix coefficient for Hi( i  ,j-1,t+dt)
    IF (j > 1) THEN
      k = k+1
      ice%dHi_m%A_index( k) = ice%dHi_ij2n( j-1,i  )
      ice%dHi_m%A_val(   k) = (-fs / grid%dx) * MAX( v_south,0._dp)
    END IF
      
    ! Matrix coefficient for Hi( i  ,j+1,t+dt)
    IF (j < grid%ny) THEN
      k = k+1
      ice%dHi_m%A_index( k) = ice%dHi_ij2n( j+1,i  )
      ice%dHi_m%A_val(   k) = ( fs / grid%dx) * MIN( v_north,0._dp)
    END IF
    
    ! Right-hand side
    ice%dHi_m%b( n) = MB_net + ice%Hi_a( j,i) / dt
    IF (i > 1      ) ice%dHi_m%b( n) = ice%dHi_m%b( n) + ice%Hi_a( j  ,i  ) * (1._dp - fs) / grid%dx * MIN( u_west, 0._dp) &
                                                       + ice%Hi_a( j  ,i-1) * (1._dp - fs) / grid%dx * MAX( u_west, 0._dp)
    IF (i < grid%nx) ice%dHi_m%b( n) = ice%dHi_m%b( n) - ice%Hi_a( j  ,i  ) * (1._dp - fs) / grid%dx * MAX( u_east, 0._dp) &
                                                       - ice%Hi_a( j  ,i+1) * (1._dp - fs) / grid%dx * MIN( u_east, 0._dp)
    IF (j > 1      ) ice%dHi_m%b( n) = ice%dHi_m%b( n) + ice%Hi_a( j  ,i  ) * (1._dp - fs) / grid%dx * MIN( v_south,0._dp) &
                                                       + ice%Hi_a( j-1,i  ) * (1._dp - fs) / grid%dx * MAX( v_south,0._dp)
    IF (j < grid%ny) ice%dHi_m%b( n) = ice%dHi_m%b( n) - ice%Hi_a( j  ,i  ) * (1._dp - fs) / grid%dx * MAX( v_north,0._dp) &
                                                       - ice%Hi_a( j+1,i  ) * (1._dp - fs) / grid%dx * MIN( v_north,0._dp)
      
    ! Initial guess (= previous solution)
    ice%dHi_m%x( n) = ice%Hi_a( j,i)
    
  END SUBROUTINE calc_dHi_dt_semiimplicit_add_matrix_coefficients_semiimplicit
    
  ! Some useful tools
  SUBROUTINE apply_ice_thickness_BC( grid, ice, dt)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: dt
    
    ! Local variables:
    INTEGER                                            :: i,j

    IF (.NOT. C%do_benchmark_experiment) THEN
          
      ! For realistic experiments set ice thickness to zero at the domain boundary
      ice%dHi_dt_a(     1              ,grid%i1:grid%i2) = -ice%Hi_a( 1              ,grid%i1:grid%i2) / dt
      ice%dHi_dt_a(     grid%ny        ,grid%i1:grid%i2) = -ice%Hi_a( grid%ny        ,grid%i1:grid%i2) / dt
      ice%dHi_dt_a(     grid%j1:grid%j2,1              ) = -ice%Hi_a( grid%j1:grid%j2,1              ) / dt
      ice%dHi_dt_a(     grid%j1:grid%j2,grid%nx        ) = -ice%Hi_a( grid%j1:grid%j2,grid%nx        ) / dt
      ice%Hi_tplusdt_a( 1              ,grid%i1:grid%i2) = 0._dp
      ice%Hi_tplusdt_a( grid%ny        ,grid%i1:grid%i2) = 0._dp
      ice%Hi_tplusdt_a( grid%j1:grid%j2,1              ) = 0._dp
      ice%Hi_tplusdt_a( grid%j1:grid%j2,grid%nx        ) = 0._dp
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
        ice%dHi_dt_a(     1              ,grid%i1:grid%i2) = -ice%Hi_a( 1              ,grid%i1:grid%i2) / dt
        ice%dHi_dt_a(     grid%ny        ,grid%i1:grid%i2) = -ice%Hi_a( grid%ny        ,grid%i1:grid%i2) / dt
        ice%dHi_dt_a(     grid%j1:grid%j2,1              ) = -ice%Hi_a( grid%j1:grid%j2,1              ) / dt
        ice%dHi_dt_a(     grid%j1:grid%j2,grid%nx        ) = -ice%Hi_a( grid%j1:grid%j2,grid%nx        ) / dt
        ice%Hi_tplusdt_a( 1              ,grid%i1:grid%i2) = 0._dp
        ice%Hi_tplusdt_a( grid%ny        ,grid%i1:grid%i2) = 0._dp
        ice%Hi_tplusdt_a( grid%j1:grid%j2,1              ) = 0._dp
        ice%Hi_tplusdt_a( grid%j1:grid%j2,grid%nx        ) = 0._dp
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
        ice%dHi_dt_a(     1              ,grid%i1:grid%i2) = (1000._dp - ice%Hi_a( 1              ,grid%i1:grid%i2)) / dt
        ice%dHi_dt_a(     grid%ny        ,grid%i1:grid%i2) = (1000._dp - ice%Hi_a( grid%ny        ,grid%i1:grid%i2)) / dt
        ice%dHi_dt_a(     grid%j1:grid%j2,1              ) = (1000._dp - ice%Hi_a( grid%j1:grid%j2,1              )) / dt
        ice%dHi_dt_a(     grid%j1:grid%j2,grid%nx        ) = (1000._dp - ice%Hi_a( grid%j1:grid%j2,grid%nx        )) / dt
        ice%Hi_tplusdt_a( 1              ,grid%i1:grid%i2) = 1000._dp
        ice%Hi_tplusdt_a( grid%ny        ,grid%i1:grid%i2) = 1000._dp
        ice%Hi_tplusdt_a( grid%j1:grid%j2,1              ) = 1000._dp
        ice%Hi_tplusdt_a( grid%j1:grid%j2,grid%nx        ) = 1000._dp
        CALL sync
        
      ELSEIF (C%choice_benchmark_experiment == 'MISMIPplus') THEN
        ! Calving front at x = 640 km
        
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          IF (grid%x( i) > 240000._dp) THEN ! x = 640 km in the original set-up, coordinates are centered around zero in IMAU-ICE
            ice%Hi_a(     j,i) = 0._dp
            ice%dHi_dt_a( j,i) = 0._dp
          END IF
        END DO
        END DO
        CALL sync
        
        ! Ice divides at west, south, and north boundaries
        ice%Hi_a(     grid%j1:grid%j2,1              ) = ice%Hi_a( grid%j1:grid%j2,2              )
        ice%Hi_a(     1,              grid%i1:grid%i2) = ice%Hi_a( 2,              grid%i1:grid%i2)
        ice%Hi_a(     grid%ny,        grid%i1:grid%i2) = ice%Hi_a( grid%ny-1,      grid%i1:grid%i2)
        ice%dHi_dt_a( grid%j1:grid%j2,1              ) = 0._dp
        ice%dHi_dt_a( 1,              grid%i1:grid%i2) = 0._dp
        ice%dHi_dt_a( grid%ny,        grid%i1:grid%i2) = 0._dp
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in apply_ice_thickness_BC!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END IF ! IF (.NOT. C%do_benchmark_experiment) THEN
    
  END SUBROUTINE apply_ice_thickness_BC
  SUBROUTINE remove_unconnected_shelves( grid, ice, dt)
    ! Use a flood-fill algorithm to find all shelves connected to sheets.
    ! Remove all other shelves.
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: dt
    
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
          ice%Hi_a(         j,i) = 0._dp
          ice%dHi_dt_a(     j,i) = -ice%Hi_a( j,i) / dt
          ice%Hi_tplusdt_a( j,i) = 0._dp
        END IF
      END DO
      END DO
      
      ! Clean up after yourself
      DEALLOCATE( map)
      DEALLOCATE( stack)
      
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE remove_unconnected_shelves
  
  ! Initialise the sparse matrix lists and grid2vector translation tables needed for the matrix-based implicit solving scheme
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
  
END MODULE ice_thickness_module
