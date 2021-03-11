MODULE ice_dynamics_module

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parallel_module,                 ONLY: par, sync, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_init_data_fields, type_SMB_model, type_BMB_model
  USE netcdf_module,                   ONLY: debug, write_to_debug_file 
  USE utilities_module,                ONLY: vertical_integration_from_bottom_to_zeta, vertical_average
  USE derivatives_and_grids_module,    ONLY: map_Acx_to_Aa_2D, map_Acy_to_Aa_2D, Neumann_BC_Aa_3D, ddx_Aa_to_Aa_3D, ddy_Aa_to_Aa_3D, &
                                             map_Aa_to_Acx_2D, map_Aa_to_Acy_2D, ddx_Aa_to_Aa_2D, ddy_Aa_to_Aa_2D, Neumann_BC_Aa_2D, &
                                             map_Acx_to_Aa_3D, map_Acy_to_Aa_3D

  IMPLICIT NONE
  
CONTAINS   

  ! Update ice thickness from ice dynamic changes and SMB
  SUBROUTINE calculate_ice_thickness_change( grid, ice, SMB, BMB, dt, mask_noice)
    ! Use the total ice velocities to update the ice thickness
    
    USE parameters_module, ONLY: ice_density
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    REAL(dp),                            INTENT(IN)    :: dt
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_noice
    
    ! Local variables
    INTEGER                                            :: i,j,cerr,ierr
    REAL(dp)                                           :: Q_SIA, Q_SSA, dVi_in, dVi_out, Vi_available, rescale_factor
    REAL(dp), DIMENSION(:,:), POINTER                  :: dVi_MB
    INTEGER                                            :: wdVi_MB
    
    IF (C%do_benchmark_experiment .AND. C%choice_benchmark_experiment == 'SSA_icestream') THEN
      ! Don't change ice thickness in this experiment, we're only interested in ice velocity
      RETURN
    END IF
        
    ! Calculate ice fluxes on the Ac grids, with ice thickness
    ! defined on the Aa grid in the upwind direction 
    ! ========================================================
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    
      ! Q_SIA: diffusive, so with staggered ice thickness
      Q_SIA = ice%U_vav_SIA_Acx( j,i) * ice%Hi_Acx( j,i) * grid%dx * dt
      
      ! Q_SSA: advective, so with upwind ice thickness
      IF (ice%U_vav_SSA_Acx( j,i) > 0._dp) THEN
        ! Ice moves from left to right, use left-hand ice thickness
        Q_SSA = ice%U_vav_SSA_Acx( j,i) * ice%Hi_Aa( j,i  ) * grid%dx * dt
      ELSE
        ! Ice moves from right to left, use right-hand ice thickness
        Q_SSA = ice%U_vav_SSA_Acx( j,i) * ice%Hi_Aa( j,i+1) * grid%dx * dt
      END IF
      
      ice%Qx_Acx( j,i) = Q_SIA + Q_SSA
      
    END DO
    END DO
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
    
      ! Q_SIA: diffusive, so with staggered ice thickness
      Q_SIA = ice%V_vav_SIA_Acy( j,i) * ice%Hi_Acy( j,i) * grid%dx * dt
      
      ! Q_SSA: advective, so with upwind ice thickness
      IF (ice%V_vav_SSA_Acy( j,i) > 0._dp) THEN
        ! Ice moves from left to right, use left-hand ice thickness
        Q_SSA = ice%V_vav_SSA_Acy( j,i) * ice%Hi_Aa( j,i  ) * grid%dx * dt
      ELSE
        ! Ice moves from right to left, use right-hand ice thickness
        Q_SSA = ice%V_vav_SSA_Acy( j,i) * ice%Hi_Aa( j+1,i) * grid%dx * dt
      END IF
      
      ice%Qy_Acy( j,i) = Q_SIA + Q_SSA
      
    END DO
    END DO
    CALL sync
    
    ! Correct outfluxes for possible resulting negative ice thicknesses
    ! =================================================================
    
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dVi_MB, wdVi_MB)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Ice volume added to each grid cell through the (surface + basal) mass balance
      dVi_MB( j,i) = (SMB%SMB_year( j,i) + BMB%BMB( j,i)) * grid%dx * grid%dx * dt
    
      ! Check how much ice is available for melting or removing (in m^3)
      Vi_available = ice%Hi_Aa( j,i) * grid%dx * grid%dx
      
      ! Check how much ice is moving into this grid cell from all four (existing) neighbours
      dVi_in   = 0._dp
      dVi_out  = 0._dp
      IF (i > 1) THEN
        IF (ice%Qx_Acx( j  ,i-1) > 0._dp) THEN
          dVi_in  = dVi_in  + ice%Qx_Acx( j  ,i-1)
        ELSE
          dVi_out = dVi_out - ice%Qx_Acx( j  ,i-1)
        END IF
      END IF
      IF (i < grid%nx) THEN
        IF (ice%Qx_Acx( j  ,i  ) < 0._dp) THEN
          dVi_in  = dVi_in  - ice%Qx_Acx( j  ,i  )
        ELSE
          dVi_out = dVi_out + ice%Qx_Acx( j  ,i  )
        END IF
      END IF
      IF (j > 1) THEN
        IF (ice%Qy_Acy( j-1,i  ) > 0._dp) THEN
          dVi_in  = dVi_in  + ice%Qy_Acy( j-1,i  )
        ELSE
          dVi_out = dVi_out - ice%Qy_Acy( j-1,i  )
        END IF
      END IF
      IF (j < grid%ny) THEN
        IF (ice%Qy_Acy( j  ,i  ) < 0._dp) THEN
          dVi_in  = dVi_in  - ice%Qy_Acy( j  ,i  )
        ELSE
          dVi_out = dVi_out + ice%Qy_Acy( j  ,i  )
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
          IF (ice%Qx_Acx( j  ,i-1) < 0._dp) ice%Qx_Acx( j  ,i-1) = ice%Qx_Acx( j  ,i-1) * rescale_factor
        END IF
        IF (i < grid%nx) THEN
          IF (ice%Qx_Acx( j  ,i  ) > 0._dp) ice%Qx_Acx( j  ,i  ) = ice%Qx_Acx( j  ,i  ) * rescale_factor
        END IF
        IF (j > 1) THEN
          IF (ice%Qy_Acy( j-1,i  ) < 0._dp) ice%Qy_Acy( j-1,i  ) = ice%Qy_Acy( j-1,i  ) * rescale_factor
        END IF
        IF (j < grid%ny) THEN
          IF (ice%Qy_Acy( j  ,i  ) > 0._dp) ice%Qy_Acy( j  ,i  ) = ice%Qy_Acy( j  ,i  ) * rescale_factor
        END IF
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Calculate change in ice thickness over time at every vertex
    ! ===========================================================
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%dHi_dt_Aa( j,i) = dVi_MB( j,i) / (grid%dx * grid%dx * dt) ! m/y
    END DO
    END DO
    CALL sync
    
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      ice%dHi_dt_Aa( j,i+1) = ice%dHi_dt_Aa( j,i+1) + ice%Qx_Acx( j,i) / (grid%dx * grid%dx * dt)
      ice%dHi_dt_Aa( j,i  ) = ice%dHi_dt_Aa( j,i  ) - ice%Qx_Acx( j,i) / (grid%dx * grid%dx * dt)
    END DO
    END DO
    CALL sync
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
      ice%dHi_dt_Aa( j+1,i) = ice%dHi_dt_Aa( j+1,i) + ice%Qy_Acy( j,i) / (grid%dx * grid%dx * dt)
      ice%dHi_dt_Aa( j  ,i) = ice%dHi_dt_Aa( j  ,i) - ice%Qy_Acy( j,i) / (grid%dx * grid%dx * dt)
    END DO
    END DO
    CALL sync
    
    ! Manually set to zero in first timestep (shouldnt be needed as velocities should be zero, but you never know...)
    IF (dt == 0._dp) ice%dHi_dt_Aa( :,grid%i1:grid%i2) = 0._dp
    
    ! Add dynamic ice thickness change and SMB to update the ice thickness
    ! ====================================================================
    
    ice%Hi_Aa_prev( :,grid%i1:grid%i2) = ice%Hi_Aa( :,grid%i1:grid%i2)
    ice%Hi_Aa(      :,grid%i1:grid%i2) = ice%Hi_Aa( :,grid%i1:grid%i2) + ice%dHi_dt_Aa( :,grid%i1:grid%i2) * dt

    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler') THEN
          
        ! Apply boundary conditions: set ice thickness to zero at the domain boundary
        ice%Hi_Aa( 1              ,grid%i1:grid%i2) = 0._dp
        ice%Hi_Aa( grid%ny        ,grid%i1:grid%i2) = 0._dp
        ice%Hi_Aa( grid%j1:grid%j2,1              ) = 0._dp
        ice%Hi_Aa( grid%j1:grid%j2,grid%nx        ) = 0._dp
        CALL sync
        
      ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
        ! No exception here, as we already exited the routine at the top.
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        
        ! Create a nice circular ice shelf
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          IF (SQRT(grid%x(i)**2+grid%y(j)**2) > grid%xmax * 0.95_dp) ice%Hi_Aa( j,i) = 0._dp
        END DO
        END DO
        CALL sync
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in calculate_ice_thickness_change!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSE
          
      ! Apply boundary conditions: set ice thickness to zero at the domain boundary
      ice%Hi_Aa( 1              ,grid%i1:grid%i2) = 0._dp
      ice%Hi_Aa( grid%ny        ,grid%i1:grid%i2) = 0._dp
      ice%Hi_Aa( grid%j1:grid%j2,1              ) = 0._dp
      ice%Hi_Aa( grid%j1:grid%j2,grid%nx        ) = 0._dp
      CALL sync
      
    END IF ! IF (C%do_benchmark_experiment) THEN

    ! Remove shelves unconnected to sheets
    CALL remove_unconnected_shelves( grid, ice)
    
    ! Remove ice in areas where no ice is allowed (i.e. Greenland in NAM and EAS, and Ellesmere Island in GRL)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (mask_noice( j,i) == 1) ice%Hi_Aa( j,i) = 0._dp
    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wdVi_MB)
    
  END SUBROUTINE calculate_ice_thickness_change
  
  ! Remove shelves unconnected to sheets
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
        IF (ice%mask_shelf_Aa( j,i) == 1) THEN
          IF (ice%mask_sheet_Aa( j  ,i+1) == 1 .OR. &
              ice%mask_sheet_Aa( j  ,i-1) == 1 .OR. &
              ice%mask_sheet_Aa( j+1,i  ) == 1 .OR. &
              ice%mask_sheet_Aa( j-1,i  ) == 1) THEN
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
        IF (ice%mask_shelf_Aa( j,i) == 1) THEN
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
        IF (ice%mask_shelf_Aa( j,i) == 1 .AND. map( j,i) == 0) THEN
          ice%Hi_Aa( j,i) = 0._dp
        END IF
      END DO
      END DO
      
      ! Clean up after yourself
      DEALLOCATE( map)
      DEALLOCATE( stack)
      
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE remove_unconnected_shelves
  
  ! SIA solver
  SUBROUTINE solve_SIA( grid, ice)
    ! Calculate ice velocities using the SIA
    
    USE parameters_module, ONLY: n_flow, ice_density, grav
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: D_0
    REAL(dp), DIMENSION(C%nZ)                          :: D_deformation
    REAL(dp)                                           :: D_uv_2D 
    REAL(dp), DIMENSION(C%nZ)                          :: D_SIA_3D_prof
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  D_SIA_Aa_from_Acx,  D_SIA_Aa_from_Acy
    INTEGER                                            :: wD_SIA_Aa_from_Acx, wD_SIA_Aa_from_Acy
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  D_SIA_3D_Aa_from_Acx,  D_SIA_3D_Aa_from_Acy
    INTEGER                                            :: wD_SIA_3D_Aa_from_Acx, wD_SIA_3D_Aa_from_Acy
    
    REAL(dp), PARAMETER                                :: D_uv_3D_cutoff = -1E5_dp
    
    ! Calculate 3D ice diffusivity
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    
      IF (ice%mask_sheet_Acx( j,i) == 1 .OR. ice%mask_gl_Acx( j,i) == 1) THEN
        
        D_0 = (ice_density * grav * ice%Hi_Acx( j,i))**n_flow * ((ice%dHs_dx_Acx( j,i)**2 + ice%dHs_dy_Acx( j,i)**2))**((n_flow - 1._dp) / 2._dp)       
        D_deformation = vertical_integration_from_bottom_to_zeta( C%m_enh_sia * ice%A_flow_Acx( :,j,i) * C%zeta**n_flow)
        D_deformation = 2._dp * ice%Hi_Acx( j,i) * D_deformation
       
        ice%D_SIA_3D_Acx( :,j,i) = D_0 * D_deformation
      
        ! Check for very large D_uv_3D's, causing very large velocities
        ice%D_SIA_3D_Acx( :,j,i) = MAX( ice%D_SIA_3D_Acx( :,j,i), D_uv_3D_cutoff)

      END IF

    END DO
    END DO
    CALL sync
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
    
      IF (ice%mask_sheet_Acy( j,i) == 1 .OR. ice%mask_gl_Acy( j,i) == 1) THEN
        
        D_0 = (ice_density * grav * ice%Hi_Acy( j,i))**n_flow * ((ice%dHs_dx_Acy( j,i)**2 + ice%dHs_dy_Acy( j,i)**2))**((n_flow - 1._dp) / 2._dp)       
        D_deformation = vertical_integration_from_bottom_to_zeta( C%m_enh_sia * ice%A_flow_Acy( :,j,i) * C%zeta**n_flow)
        D_deformation = 2._dp * ice%Hi_Acy( j,i) * D_deformation
       
        ice%D_SIA_3D_Acy( :,j,i) = D_0 * D_deformation
      
        ! Check for very large D_uv_3D's, causing very large velocities
        ice%D_SIA_3D_Acy( :,j,i) = MAX( ice%D_SIA_3D_Acy( :,j,i), D_uv_3D_cutoff)

      END IF

    END DO
    END DO
    CALL sync
    
    ! Calculate vertically averaged ice velocities on the Ac grids
    DO i = grid%i1, MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
    
      IF (ice%mask_sheet_Acx( j,i) == 1 .OR. ice%mask_gl_Acx( j,i) == 1) THEN
        D_SIA_3D_prof = ice%D_SIA_3D_Acx( :,j,i)
        D_uv_2D = vertical_average( D_SIA_3D_prof)
       
        ice%D_SIA_Acx(     j,i) = ice%Hi_Acx( j,i) * D_uv_2D
        ice%U_vav_SIA_Acx( j,i) = D_uv_2D * ice%dHs_dx_Acx( j,i)
      END IF
          
    END DO
    END DO
    CALL sync
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny-1
    
      IF (ice%mask_sheet_Acy( j,i) == 1 .OR. ice%mask_gl_Acy( j,i) == 1) THEN
        D_SIA_3D_prof = ice%D_SIA_3D_Acy( :,j,i)
        D_uv_2D = vertical_average( D_SIA_3D_prof)
       
        ice%D_SIA_Acy(     j,i) = ice%Hi_Acy( j,i) * D_uv_2D
        ice%V_vav_SIA_Acy( j,i) = D_uv_2D * ice%dHs_dy_Acy( j,i)
      END IF
          
    END DO
    END DO
    CALL sync
    
    ! Map data to Aa grid for writing to output (not actually used anywhere...)
    CALL allocate_shared_dp_2D(       grid%ny, grid%nx, D_SIA_Aa_from_Acx,    wD_SIA_Aa_from_Acx   )
    CALL allocate_shared_dp_2D(       grid%ny, grid%nx, D_SIA_Aa_from_Acy,    wD_SIA_Aa_from_Acy   )
    CALL allocate_shared_dp_3D( C%nZ, grid%ny, grid%nx, D_SIA_3D_Aa_from_Acx, wD_SIA_3D_Aa_from_Acx)
    CALL allocate_shared_dp_3D( C%nZ, grid%ny, grid%nx, D_SIA_3D_Aa_from_Acy, wD_SIA_3D_Aa_from_Acy)
    
    CALL map_Acx_to_Aa_2D( grid, ice%D_SIA_Acx,    D_SIA_Aa_from_Acx)
    CALL map_Acy_to_Aa_2D( grid, ice%D_SIA_Acy,    D_SIA_Aa_from_Acy)
    CALL map_Acx_to_Aa_3D( grid, ice%D_SIA_3D_Acx, D_SIA_3D_Aa_from_Acx)
    CALL map_Acy_to_Aa_3D( grid, ice%D_SIA_3D_Acy, D_SIA_3D_Aa_from_Acy)
    
    ice%D_SIA_Aa(      :,grid%i1:grid%i2) = (D_SIA_Aa_from_Acx(      :,grid%i1:grid%i2) + D_SIA_Aa_from_Acy(      :,grid%i1:grid%i2)) / 2._dp
    ice%D_SIA_3D_Aa( :,:,grid%i1:grid%i2) = (D_SIA_3D_Aa_from_Acx( :,:,grid%i1:grid%i2) + D_SIA_3D_Aa_from_Acy( :,:,grid%i1:grid%i2)) / 2._dp
    
    CALL deallocate_shared( wD_SIA_Aa_from_Acx)
    CALL deallocate_shared( wD_SIA_Aa_from_Acy)
    CALL deallocate_shared( wD_SIA_3D_Aa_from_Acx)
    CALL deallocate_shared( wD_SIA_3D_Aa_from_Acy)
    
    CALL map_Acx_to_Aa_2D( grid, ice%U_vav_SIA_Acx, ice%U_vav_SIA_Aa)
    CALL map_Acy_to_Aa_2D( grid, ice%V_vav_SIA_Acy, ice%V_vav_SIA_Aa)
      
  END SUBROUTINE solve_SIA

  ! 3D SIA solver (used only for thermodynamics)
  SUBROUTINE solve_SIA_3D( grid, ice)
    
    USE parameters_module, ONLY: n_flow, ice_density, grav
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
      
    ! Local variables:    
    REAL(dp)                                           :: D_0
    REAL(dp), DIMENSION(C%nZ)                          :: D_deformation
    REAL(dp), DIMENSION(C%nZ)                          :: D_SIA_3D
    REAL(dp), PARAMETER                                :: D_uv_3D_cutoff = -1E5_dp 
    REAL(dp)                                           :: dHb_dx
    REAL(dp)                                           :: dHb_dy
    INTEGER                                            :: i,j,k
    REAL(dp), DIMENSION(C%nZ,grid%ny,grid%nx)          :: dU_dx, dU_dy, dV_dx, dV_dy
    REAL(dp)                                           :: w1, w2, w3, w4
                                        
    ice%U_3D_Aa(:,:,grid%i1:grid%i2) = 0._dp
    ice%V_3D_Aa(:,:,grid%i1:grid%i2) = 0._dp
    ice%W_3D_Aa(:,:,grid%i1:grid%i2) = 0._dp
    CALL sync

    ! Direct calculation of U and V for the ice sheet:
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ice%U_3D_Aa( :,j,i) = ice%U_base_Aa( j,i)
      ice%V_3D_Aa( :,j,i) = ice%V_base_Aa( j,i)
    
      IF (ice%mask_shelf_Aa( j,i) == 1    ) CYCLE
      IF (ice%Hi_Aa(         j,i) == 0._dp) CYCLE
        
      D_0           = (ice_density * grav * ice%Hi_Aa( j,i))**n_flow * ((ice%dHs_dx_Aa( j,i)**2 + ice%dHs_dy_Aa( j,i)**2))**((n_flow - 1._dp) / 2._dp)
      D_deformation = vertical_integration_from_bottom_to_zeta( C%m_enh_sia * ice%A_flow_Aa( :,j,i) * C%zeta**n_flow)
      D_deformation = 2._dp * ice%Hi_Aa( j,i) * D_deformation
      D_SIA_3D      = MAX(D_0 * D_deformation, D_uv_3D_cutoff)
       
      ice%U_3D_Aa( :,j,i) = D_SIA_3D * ice%dHs_dx_Aa( j,i) + ice%U_base_Aa( j,i)
      ice%V_3D_Aa( :,j,i) = D_SIA_3D * ice%dHs_dy_Aa( j,i) + ice%V_base_Aa( j,i)
      
    END DO
    END DO
    CALL sync
        
    CALL Neumann_BC_Aa_3D( grid, ice%U_3D_Aa)
    CALL Neumann_BC_Aa_3D( grid, ice%V_3D_Aa)
    
    CALL ddx_Aa_to_Aa_3D( grid, ice%U_3D_Aa, dU_dx)
    CALL ddy_Aa_to_Aa_3D( grid, ice%U_3D_Aa, dU_dy)
    CALL ddx_Aa_to_Aa_3D( grid, ice%V_3D_Aa, dV_dx)
    CALL ddy_Aa_to_Aa_3D( grid, ice%V_3D_Aa, dV_dy)

    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    
      IF (ice%mask_sheet_Aa( j,i) == 0) CYCLE
      
      dHb_dx = ice%dHs_dx_Aa( j,i) - ice%dHi_dx_Aa( j,i)
      dHb_dy = ice%dHs_dy_Aa( j,i) - ice%dHi_dy_Aa( j,i)   
      
      ice%W_3D_Aa( C%nZ,j,i) = ice%dHb_dt_Aa( j,i) + ice%U_3D_Aa( C%nZ,j,i) * dHb_dx + ice%V_3D_Aa( C%nZ,j,i) * dHb_dy 
                           
      ! The integrant is calculated half way the layer of integration at k+1/2. This integrant is multiplied with the layer thickness and added to the integral
      ! of all layers below, giving the integral up to and including this layer:
      DO k = C%nZ - 1, 1, -1
        
        w1 = (dU_dx( k,j,i) + dU_dx( k+1,j,i)) / 2._dp
        w2 = (dV_dy( k,j,i) + dV_dy( k+1,j,i)) / 2._dp   
        
        w3              = ((ice%dHs_dx_Aa( j,i) - 0.5_dp * (C%zeta(k+1) + C%zeta(k)) * ice%dHi_dx_Aa( j,i)) / MAX(0.1_dp, ice%Hi_Aa( j,i))) *   &
                          ((ice%U_3D_Aa( k+1,j,i) - ice%U_3D_Aa( k,j,i)) / (C%zeta(k+1) - C%zeta(k)))
        w4              = ((ice%dHs_dy_Aa( j,i) - 0.5_dp * (C%zeta(k+1) + C%zeta(k)) * ice%dHi_dy_Aa( j,i)) / MAX(0.1_dp, ice%Hi_Aa( j,i))) *   &
                          ((ice%V_3D_Aa( k+1,j,i) - ice%V_3D_Aa( k,j,i)) / (C%zeta(k+1) - C%zeta(k)))

        ice%W_3D_Aa( k,j,i) = ice%W_3D_Aa( k+1,j,i) - ice%Hi_Aa( j,i) * (w1 + w2 + w3 + w4) * (C%zeta(k+1) - C%zeta(k))
        
      END DO ! DO k = C%nZ - 1, 1, -1

    END DO
    END DO
    CALL sync
    
    CALL Neumann_BC_Aa_3D( grid, ice%W_3D_Aa)

  END SUBROUTINE solve_SIA_3D
  
  ! SSA solver
  SUBROUTINE solve_SSA( grid, ice)
    ! Calculate ice velocities using the SSA
    
    USE parameters_module,          ONLY: ice_density, seawater_density
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j, cerr, ierr
    LOGICAL                                            :: set_SSA_velocities_to_zero
    LOGICAL                                            :: has_converged
    INTEGER                                            :: viscosity_iteration_i
    REAL(dp)                                           :: sum_DN_sq, sum_N_sq, RN
    
    ! Exceptions for benchmark experiments
    set_SSA_velocities_to_zero = .FALSE.
    IF (C%do_benchmark_experiment) THEN
      IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
              C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'Bueler') THEN
        ! The SSA is not solved in these experiments
        set_SSA_velocities_to_zero = .TRUE.
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
              C%choice_benchmark_experiment == 'SSA_icestream') THEN
        ! The SSA is solved in these experiments
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in solve_SSA!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! If there's no grounded ice anywhere, don't bother
    IF (SUM( ice%mask_sheet_Aa) == 0) set_SSA_velocities_to_zero = .TRUE.
    
    IF (set_SSA_velocities_to_zero) THEN
      ice%U_vav_SSA_Aa( :,grid%i1:grid%i2) = 0._dp
      ice%V_vav_SSA_Aa( :,grid%i1:grid%i2) = 0._dp
      ice%U_vav_SSA_Acx(:,grid%i1:MIN(grid%nx-1,grid%i2)) = 0._dp
      ice%V_vav_SSA_Acy(:,grid%i1:grid%i2) = 0._dp
      CALL sync
      RETURN
    END IF
    
    ! Calculate the basal yield stress tau_c
    CALL basal_yield_stress( grid, ice)
    
    ! One-sided differenced surface slope to get accurate grounding-line velocities
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%mask_sheet_Aa( j,i)==1) THEN
        ice%dHs_dx_shelf_Aa( j,i) = ice%dHs_dx_Aa( j,i)
        ice%dHs_dy_shelf_Aa( j,i) = ice%dHs_dy_Aa( j,i)
      ELSE
        ice%dHs_dx_shelf_Aa( j,i) = (1._dp - ice_density / seawater_density) * ice%dHi_dx_Aa( j,i)
        ice%dHs_dy_shelf_Aa( j,i) = (1._dp - ice_density / seawater_density) * ice%dHi_dy_Aa( j,i)
      END IF
    END DO
    END DO
    CALL sync
    
    ! The viscosity iteration
    has_converged = .FALSE.
    viscosity_iteration_i  = 0
    viscosity_iteration: DO WHILE ((.NOT. has_converged) .AND. (viscosity_iteration_i < C%SSA_max_outer_loops))
      viscosity_iteration_i = viscosity_iteration_i + 1
      
      ! Store the previous solution, so we can check for converge later
      ice%N_Aa_prev( 2:grid%ny-1,MAX(2,grid%i1):MIN(grid%nx-1,grid%i2))  = ice%N_Aa( 2:grid%ny-1,MAX(2,grid%i1):MIN(grid%nx-1,grid%i2))
      CALL sync
    
      ! Update the effective viscosity eta and the product term N = eta * H
      CALL SSA_effective_viscosity( grid, ice)
      
      ! Check if a stable solution has been reached
      sum_DN_sq = SUM( (ice%N_Aa(      2:grid%ny-1, MAX(2,grid%i1):MIN(grid%nx-1,grid%i2)) - &
                        ice%N_Aa_prev( 2:grid%ny-1, MAX(2,grid%i1):MIN(grid%nx-1,grid%i2)))**2)
      sum_N_sq  = SUM(  ice%N_Aa(      2:grid%ny-1, MAX(2,grid%i1):MIN(grid%nx-1,grid%i2))**2)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, sum_DN_sq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, sum_N_sq,  1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      RN = SQRT( sum_DN_sq / sum_N_sq)
      
      !IF (par%master) WRITE(0,'(A,I5,A,E10.4)') '   SSA - viscosity iteration ', viscosity_iteration_i, ', RN = ', RN
      
      IF (RN < C%SSA_RN_tol) THEN
        has_converged = .TRUE.
        EXIT viscosity_iteration
      END IF
    
      ! Calculate the sliding term S = tau_c / norm([U,V])
      CALL SSA_sliding_term( grid, ice)
      
      ! Solve the linearised SSA using SOR
      CALL solve_SSA_linearised( grid, ice)
            
    END DO viscosity_iteration ! DO WHILE ((.NOT. has_converged) .AND. (viscosity_iteration_i < C%SSA_max_outer_loops))
    
    ! Map velocities to the Acx grid for use in ice thickness updating
    CALL map_Aa_to_Acx_2D( grid, ice%U_vav_SSA_Aa, ice%U_vav_SSA_Acx)
    CALL map_Aa_to_Acy_2D( grid, ice%V_vav_SSA_Aa, ice%V_vav_SSA_Acy)
    
  END SUBROUTINE solve_SSA
  SUBROUTINE solve_SSA_linearised( grid, ice)
    ! Calculate ice velocities using the SSA
    
    USE parameters_module,          ONLY: ice_density, grav
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: cerr,ierr
    INTEGER                                            :: i,j,m,j0
    LOGICAL                                            :: has_converged
    INTEGER                                            :: inner_loop_i
    REAL(dp)                                           :: max_residual_UV
    
    ! Pointer to bind to ice%U_vav_SSA_Aa for easier notation
    REAL(dp), DIMENSION(:,:), POINTER                  :: U,V
    U => ice%U_vav_SSA_Aa
    V => ice%V_vav_SSA_Aa
       
    ! Calculate the right-hand sides of the equations (i.e. the gravitational driving stress)
    DO i = grid%i1,grid%i2
    DO j = 2, grid%ny-1
      ice%RHSx( j,i) = ice_density * grav * ice%dHs_dx_shelf_Aa( j,i) / ice%eta_Aa( j,i)
      ice%RHSy( j,i) = ice_density * grav * ice%dHs_dy_shelf_Aa( j,i) / ice%eta_Aa( j,i)
    END DO
    END DO
    CALL sync
    
    ! Calculate the centre coefficients
    DO i = grid%i1,grid%i2
    DO j = 2, grid%ny-1
      IF (ice%mask_ocean_Aa( j,i) == 0) THEN
        ice%eu_ij( j,i) = (-8._dp / grid%dx**2) + (-2._dp / grid%dx**2) - ice%S_Aa( j,i) / (MAX(0.1_dp,ice%Hi_Aa( j,i)) * ice%eta_Aa( j,i))
        ice%ev_ij( j,i) = (-2._dp / grid%dx**2) + (-8._dp / grid%dx**2) - ice%S_Aa( j,i) / (MAX(0.1_dp,ice%Hi_Aa( j,i)) * ice%eta_Aa( j,i))
      ELSE
        ice%eu_ij( j,i) = (-8._dp / grid%dx**2) + (-2._dp / grid%dx**2)
        ice%ev_ij( j,i) = (-2._dp / grid%dx**2) + (-8._dp / grid%dx**2)
      END IF
    END DO
    END DO
    CALL sync
    
    ! SOR iteration
    has_converged = .FALSE.
    inner_loop_i  = 0
    DO WHILE ((.NOT. has_converged) .AND. (inner_loop_i < C%SSA_max_inner_loops))
      inner_loop_i = inner_loop_i + 1
      
      max_residual_UV = 0._dp
      
      ! Red-black scheme
      DO m = 0,1
      
        DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
        
        j0 = MOD(i+m,2)
        
        DO j = j0+2, grid%ny-1, 2
          
          ! Calculate the left-hand side of the equation, using partially updated values
          ice%LHSx( j,i) = (4._dp * (U( j  ,i+1) + U( j  ,i-1)) / grid%dx**2) + &
                           (        (U( j+1,i  ) + U( j-1,i  )) / grid%dx**2) + &
                           (3._dp * (V( j+1,i+1) + V( j-1,i-1) - V( j+1,i-1) - V( j-1,i+1)) / (4._dp*grid%dx*grid%dx)) + (ice%eu_ij( j,i) * U( j,i))
          ice%LHSy( j,i) = (4._dp * (V( j+1,i  ) + V( j-1,i  )) / grid%dx**2) + &
                           (        (V( j  ,i+1) + V( j  ,i-1)) / grid%dx**2) + &
                           (3._dp * (U( j+1,i+1) + U( j-1,i-1) - U( j+1,i-1) - U( j-1,i+1)) / (4._dp*grid%dx*grid%dx)) + (ice%ev_ij( j,i) * V( j,i))
          
          ! Calculate the residuals
          ice%resU( j,i) = (ice%LHSx( j,i) - ice%RHSx( j,i)) / ice%eu_ij( j,i)
          ice%resV( j,i) = (ice%LHSy( j,i) - ice%RHSy( j,i)) / ice%ev_ij( j,i)
          
          max_residual_UV = MAX( max_residual_UV, ABS(ice%resU( j,i)))
          max_residual_UV = MAX( max_residual_UV, ABS(ice%resV( j,i)))
          
          ! Update velocities
          U( j,i) = U( j,i) - C%SSA_SOR_omega * ice%resU( j,i)
          V( j,i) = V( j,i) - C%SSA_SOR_omega * ice%resV( j,i)
          
          IF (U( j,i) /= U( j,i) .OR. V( j,i) /= V( j,i)) THEN
            WRITE(0,*) ' SSA - NaN ice velocities in SSA solver!'
            WRITE(0,*) '   LHSx = ', ice%LHSx( j,i), ', RHSx = ', ice%RHSx( j,i), ', eu_ij = ', ice%eu_ij( j,i), ', resU = ', ice%resU( j,i), ', dHs_dx_shelf = ', ice%dHs_dx_shelf_Aa( j,i)
            WRITE(0,*) '   LHSy = ', ice%LHSy( j,i), ', RHSy = ', ice%RHSy( j,i), ', ev_ij = ', ice%ev_ij( j,i), ', resV = ', ice%resV( j,i), ', dHs_dy_shelf = ', ice%dHs_dy_shelf_Aa( j,i)
            WRITE(0,*) '   dHs_dx = ', ice%dHs_dx_Aa( j,i), ', dHi_dx = ', ice%dHi_dx_Aa( j,i)
            WRITE(0,*) '   dHs_dy = ', ice%dHs_dy_Aa( j,i), ', dHi_dy = ', ice%dHi_dy_Aa( j,i)
            WRITE(0,*) '   eta = ', ice%eta_Aa( j,i), ', A    = ', ice%A_flow_mean_Aa( j,i), ', Ti = ', ice%Ti_Aa( :,j,i)
            CALl MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
          END IF
        
        END DO ! DO j = j0+2, grid%ny-1, 2
        END DO ! DO i = 2, grid%nx-1
        CALL sync
      
      END DO ! DO m = 0, 1
      
      ! Apply Neumann boundary conditions
      CALL Neumann_BC_Aa_2D( grid, U)
      CALL Neumann_BC_Aa_2D( grid, V)
      
      ! Check if we've reached a stable solution
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, max_residual_UV, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      
      !IF (par%master) WRITE(0,*) ' SSA -  inner loop ', inner_loop_i, ': largest residual = ', max_residual_UV
      IF (max_residual_UV < C%SSA_max_residual_UV) THEN
        has_converged = .TRUE.
      ELSEIF (max_residual_UV > 1E6_dp) THEN
        WRITE(0,*) ' ERROR - instability in SSA SOR solver!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (inner_loop_i == C%SSA_max_inner_loops) THEN
        WRITE(0,*) ' WARNING - SSA SOR solver doesnt converge!'
        !IF (par%master) WRITE(0,*) ' ERROR - SSA SOR solver doesnt converge!'
        !CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END DO ! DO WHILE ((.NOT. has_converged) .AND. (inner_loop_i < max_inner_loops))
    CALL sync
    
  END SUBROUTINE solve_SSA_linearised
  SUBROUTINE SSA_effective_viscosity( grid, ice)
    ! Calculate the effective viscosity eta, the product term N = eta*H, and the gradient of N in the SSA
    
    USE parameters_module, ONLY: n_flow
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    REAL(dp), PARAMETER                                :: epsilon_sq_0 = 1E-12_dp
    
    ! Calculate derivatives of U and V
    CALL ddx_Aa_to_Aa_2D( grid, ice%U_vav_SSA_Aa, ice%Ux_Aa)
    CALL ddy_Aa_to_Aa_2D( grid, ice%U_vav_SSA_Aa, ice%Uy_Aa)
    CALL ddx_Aa_to_Aa_2D( grid, ice%V_vav_SSA_Aa, ice%Vx_Aa)
    CALL ddy_Aa_to_Aa_2D( grid, ice%V_vav_SSA_Aa, ice%Vy_Aa)
    
    ! The effective viscosity eta on both grids, with a normalisation term (following Bueler & Brown, 2009) to prevent divide-by-zero errors,
    ! and the product term N = eta * H
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%eta_Aa( j,i) = 0.5_dp * (C%m_enh_SSA * ice%A_flow_mean_Aa( j,i))**(-1._dp/n_flow) * (&
        ice%Ux_Aa( j,i)**2 + ice%Vy_Aa( j,i)**2 + ice%Ux_Aa( j,i)*ice%Vy_Aa( j,i) + 0.25_dp*(ice%Uy_Aa( j,i)+ice%Vx_Aa( j,i))**2 + epsilon_sq_0) ** ((1._dp - n_flow) / (2._dp * n_flow))
      ice%N_Aa( j,i) = ice%eta_Aa( j,i) * MAX(0.1_dp, ice%Hi_Aa( j,i))
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE SSA_effective_viscosity
  SUBROUTINE SSA_sliding_term( grid, ice)
    ! Calculate the sliding term S = tau_c / norm([U,V]) in the SSA
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: i,j
    REAL(dp), PARAMETER                                :: delta_v = 1E-3_dp                    ! Normalisation parameter to prevent errors when velocity is zero
    REAL(dp), PARAMETER                                :: q_plastic            = 0.30_dp       ! Parameter used for basal stress (inverse of m_flow)
    REAL(dp), PARAMETER                                :: u_threshold          = 100._dp       ! scaling of tau_yield to get the correct unit (function of q_plastic)
    
    IF (C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Regularised Coulomb-type sliding law
      
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ! The sliding term S, with a normalisation term (following Bueler & Brown, 2009) to prevent divide-by-zero errors.
        ice%S_Aa( j,i) = ice%tau_c_Aa( j,i) * ( (delta_v**2 + ice%U_vav_SSA_Aa( j,i)**2 + ice%V_vav_SSA_Aa( j,i)**2)**(0.5_dp * (q_plastic-1._dp)) ) / (u_threshold**q_plastic)
      END DO
      END DO
      CALL sync
      
    ELSE
      WRITE(0,*) ' ERROR: choice_sliding_law "', C%choice_sliding_law,'" not implemented in SSA_sliding_term!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE SSA_sliding_term
  SUBROUTINE basal_yield_stress( grid, ice)
    ! Calculate the basal yield stress for a regularised Coulomb sliding law
    
    USE parameters_module,          ONLY: pi, n_flow, ice_density, seawater_density, grav

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
  
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! The pore water pressure is scaled with a bedrock height dependend parameterisation
      ! Equation (12) in Martin et al. (2011)
      IF     ((ice%Hb_Aa( j,i) - ice%SL_Aa( j,i)) <= 0._dp) THEN
       lambda_p = 1._dp
      ELSEIF ((ice%Hb_Aa( j,i) - ice%SL_Aa( j,i)) >= 1000._dp) THEN
       lambda_p = 0._dp
      ELSE ! between 0 and 1000
       lambda_p = 1._dp - (ice%Hb_Aa( j,i) - ice%SL_Aa( j,i)) / 1000._dp
      END IF
  
      ! The pore water pressure, equation (11) in Martin et al. (2011)
      pore_water_pressure = 0.96_dp * ice_density * grav * MAX(0.1_dp, ice%Hi_Aa( j,i)) * lambda_p 
  
      ! The friction angle, used for the yield stress, equation (10) in Martin et al. (2011)
      IF     (ice%Hb_Aa( j,i) <= scaling_Hb_minimal_friction_angle) THEN
        ice%phi_fric_Aa( j,i) = minimal_friction_angle
      ELSEIF (ice%Hb_Aa( j,i) >= scaling_Hb_maximal_friction_angle) THEN
        ice%phi_fric_Aa( j,i) = maximal_friction_angle
      ELSE ! between C%scaling_Hb_maximal_friction_angle and C%scaling_Hb_minimal_friction_angle
        ice%phi_fric_Aa( j,i) = minimal_friction_angle + (maximal_friction_angle - minimal_friction_angle) * (1._dp &
                    + (ice%Hb_Aa( j,i) - scaling_Hb_maximal_friction_angle) / (scaling_Hb_maximal_friction_angle - scaling_Hb_minimal_friction_angle))
      END IF
  
      ! calculate yield stress everywhere, restrictions to sheet applied elsewhere, equation (9) in Martin et al. (2011)
      ice%tau_c_Aa( j,i) = TAN((pi / 180._dp) * ice%phi_fric_Aa( j,i)) * (ice_density * grav * MAX(0.1_dp, ice%Hi_Aa( j,i)) - pore_water_pressure)
    
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE basal_yield_stress
  
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
    
    IF (par%master) WRITE(0,*) '  Initialising ice dynamics model...'
    
    ! Allocate shared memory
    CALL allocate_ice_model( grid, ice)
        
    ! Initialise with data from initial file
    IF (.NOT. C%is_restart) THEN
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%Hi_Aa( j,i) = init%Hi( j,i)
        ice%Hb_Aa( j,i) = init%Hb( j,i)
        ice%Hs_Aa( j,i) = MAX( ice%SL_Aa( j,i), ice%Hb_Aa( j,i) + ice%Hi_Aa( j,i))
      END DO
      END DO
      CALL sync
    
    ELSE
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%Hi_Aa(        j,i) = init%Hi( j,i)
        ice%Hb_Aa(        j,i) = init%Hb( j,i)
        ice%Hs_Aa(        j,i) = MAX( ice%SL_Aa( j,i), ice%Hb_Aa( j,i) + ice%Hi_Aa( j,i))
        ice%U_vav_SSA_Aa( j,i) = init%U_SSA( j,i)
        ice%V_vav_SSA_Aa( j,i) = init%V_SSA( j,i)
      END DO
      END DO
      CALL sync
      
    END IF ! IF (.NOT. C%is_restart) THEN
    
  END SUBROUTINE initialise_ice_model
  SUBROUTINE allocate_ice_model( grid, ice)  
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
      
    ! Allocate memory
    ! ===============
    
    ! Basic data - ice thickness, bedrock & surface elevation, sea level (geoid elevation), englacial temperature
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_Aa                , ice%wHi_Aa                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%Hi_Acx               , ice%wHi_Acx               )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%Hi_Acy               , ice%wHi_Acy               )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%Hi_Ab                , ice%wHi_Ab                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hb_Aa                , ice%wHb_Aa                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%Hb_Acx               , ice%wHb_Acx               )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%Hb_Acy               , ice%wHb_Acy               )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%Hb_Ab                , ice%wHb_Ab                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hs_Aa                , ice%wHs_Aa                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%Hs_Acx               , ice%wHs_Acx               )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%Hs_Acy               , ice%wHs_Acy               )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%Hs_Ab                , ice%wHs_Ab                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%SL_Aa                , ice%wSL_Aa                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%SL_Acx               , ice%wSL_Acx               )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%SL_Acy               , ice%wSL_Acy               )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%SL_Ab                , ice%wSL_Ab                )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%Ti_Aa                , ice%wTi_Aa                )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx-1, ice%Ti_Acx               , ice%wTi_Acx               )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny-1, grid%nx  , ice%Ti_Acy               , ice%wTi_Acy               )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny-1, grid%nx-1, ice%Ti_Ab                , ice%wTi_Ab                )
    
    ! Ice velocities
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%U_vav_SIA_Aa         , ice%wU_vav_SIA_Aa         )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%V_vav_SIA_Aa         , ice%wV_vav_SIA_Aa         )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%U_vav_SIA_Acx        , ice%wU_vav_SIA_Acx        )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%V_vav_SIA_Acy        , ice%wV_vav_SIA_Acy        )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%U_vav_SSA_Aa         , ice%wU_vav_SSA_Aa         )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%V_vav_SSA_Aa         , ice%wV_vav_SSA_Aa         )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%U_vav_SSA_Acx        , ice%wU_vav_SSA_Acx        )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%V_vav_SSA_Acy        , ice%wV_vav_SSA_Acy        )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%U_vav_Aa             , ice%wU_vav_Aa             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%V_vav_Aa             , ice%wV_vav_Aa             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%U_vav_Acx            , ice%wU_vav_Acx            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%V_vav_Acy            , ice%wV_vav_Acy            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%U_surf_Aa            , ice%wU_surf_Aa            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%V_surf_Aa            , ice%wV_surf_Aa            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%U_base_Aa            , ice%wU_base_Aa            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%V_base_Aa            , ice%wV_base_Aa            )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%U_3D_Aa              , ice%wU_3D_Aa              )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%V_3D_Aa              , ice%wV_3D_Aa              )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%W_3D_Aa              , ice%wW_3D_Aa              )
    
    ! Different masks
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_land_Aa         , ice%wmask_land_Aa         )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_ocean_Aa        , ice%wmask_ocean_Aa        )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%mask_ocean_Acx       , ice%wmask_ocean_Acx       )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%mask_ocean_Acy       , ice%wmask_ocean_Acy       )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_lake_Aa         , ice%wmask_lake_Aa         )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_ice_Aa          , ice%wmask_ice_Aa          )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_sheet_Aa        , ice%wmask_sheet_Aa        )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%mask_sheet_Acx       , ice%wmask_sheet_Acx       )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%mask_sheet_Acy       , ice%wmask_sheet_Acy       )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_shelf_Aa        , ice%wmask_shelf_Aa        )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%mask_shelf_Acx       , ice%wmask_shelf_Acx       )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%mask_shelf_Acy       , ice%wmask_shelf_Acy       )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_coast_Aa        , ice%wmask_coast_Aa        )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%mask_coast_Acx       , ice%wmask_coast_Acx       )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%mask_coast_Acy       , ice%wmask_coast_Acy       )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_margin_Aa       , ice%wmask_margin_Aa       )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%mask_margin_Acx      , ice%wmask_margin_Acx      )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%mask_margin_Acy      , ice%wmask_margin_Acy      )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_gl_Aa           , ice%wmask_gl_Aa           )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%mask_gl_Acx          , ice%wmask_gl_Acx          )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%mask_gl_Acy          , ice%wmask_gl_Acy          )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_cf_Aa           , ice%wmask_cf_Aa           )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%mask_cf_Acx          , ice%wmask_cf_Acx          )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%mask_cf_Acy          , ice%wmask_cf_Acy          )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx  , ice%mask_Aa              , ice%wmask_Aa              )
    CALL allocate_shared_int_2D(       grid%ny  , grid%nx-1, ice%mask_Acx             , ice%wmask_Acx             )
    CALL allocate_shared_int_2D(       grid%ny-1, grid%nx  , ice%mask_Acy             , ice%wmask_Acy             )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%grounded_fraction    , ice%wgrounded_fraction    )
    
    ! Ice physical properties
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%A_flow_Aa            , ice%wA_flow_Aa            )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx-1, ice%A_flow_Acx           , ice%wA_flow_Acx           )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny-1, grid%nx  , ice%A_flow_Acy           , ice%wA_flow_Acy           )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%A_flow_mean_Aa       , ice%wA_flow_mean_Aa       )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%A_flow_mean_Acx      , ice%wA_flow_mean_Acx      )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%A_flow_mean_Acy      , ice%wA_flow_mean_Acy      )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx-1, ice%A_flow_mean_Ab       , ice%wA_flow_mean_Ab       )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%Ti_pmp_Aa            , ice%wTi_pmp_Aa            )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%Cpi_Aa               , ice%wCpi_Aa               )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%Ki_Aa                , ice%wKi_Aa                )
    
    ! Spatial and temporal derivatives
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHi_dt_Aa            , ice%wdHi_dt_Aa            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHi_dx_Aa            , ice%wdHi_dx_Aa            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHi_dy_Aa            , ice%wdHi_dy_Aa            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHb_dt_Aa            , ice%wdHb_dt_Aa            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHs_dt_Aa            , ice%wdHs_dt_Aa            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHs_dx_Aa            , ice%wdHs_dx_Aa            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHs_dy_Aa            , ice%wdHs_dy_Aa            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%dHs_dx_Acx           , ice%wdHs_dx_Acx           )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%dHs_dy_Acx           , ice%wdHs_dy_Acx           )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%dHs_dx_Acy           , ice%wdHs_dx_Acy           )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%dHs_dy_Acy           , ice%wdHs_dy_Acy           )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%dTi_dx_Aa            , ice%wdTi_dx_Aa            )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%dTi_dy_Aa            , ice%wdTi_dy_Aa            )
    
    ! Zeta derivatives
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%dzeta_dt_Aa          , ice%wdzeta_dt_Aa          )    
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%dzeta_dx_Aa          , ice%wdzeta_dx_Aa          )    
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%dzeta_dy_Aa          , ice%wdzeta_dy_Aa          )    
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dzeta_dz_Aa          , ice%wdzeta_dz_Aa          )
        
    ! Ice dynamics - SIA
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx-1, ice%D_SIA_3D_Acx         , ice%wD_SIA_3D_Acx         )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny-1, grid%nx  , ice%D_SIA_3D_Acy         , ice%wD_SIA_3D_Acy         )
    CALL allocate_shared_dp_3D(  C%nZ, grid%ny  , grid%nx  , ice%D_SIA_3D_Aa          , ice%wD_SIA_3D_Aa          )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%D_SIA_Acx            , ice%wD_SIA_Acx            )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%D_SIA_Acy            , ice%wD_SIA_Acy            )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%D_SIA_Aa             , ice%wD_SIA_Aa             )
    
    ! Ice dynamics - SSA
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%phi_fric_Aa          , ice%wphi_fric_Aa          )    
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%tau_c_Aa             , ice%wtau_c_Aa             )    
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Ux_Aa                , ice%wUx_Aa                )    
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Uy_Aa                , ice%wUy_Aa                )   
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Vx_Aa                , ice%wVx_Aa                )    
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Vy_Aa                , ice%wVy_Aa                )    
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%eta_Aa               , ice%weta_Aa               )    
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%N_Aa                 , ice%wN_Aa                 )   
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%N_Aa_prev            , ice%wN_Aa_prev            )     
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%S_Aa                 , ice%wS_Aa                 )   
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHs_dx_shelf_Aa      , ice%wdHs_dx_shelf_Aa      ) 
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%dHs_dy_shelf_Aa      , ice%wdHs_dy_shelf_Aa      )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%LHSx                 , ice%wLHSx                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%LHSy                 , ice%wLHSy                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%RHSx                 , ice%wRHSx                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%RHSy                 , ice%wRHSy                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%resU                 , ice%wresU                 )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%resV                 , ice%wresV                 ) 
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%eu_ij                , ice%weu_ij                )
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%ev_ij                , ice%wev_ij                ) 
    
    ! Ice dynamics - ice thickness calculation
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx-1, ice%Qx_Acx               , ice%wQx_Acx               )
    CALL allocate_shared_dp_2D(        grid%ny-1, grid%nx  , ice%Qy_Acy               , ice%wQy_Acy               )
    
    ! Thermodynamics
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%frictional_heating_Aa, ice%wfrictional_heating_Aa)
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%GHF_Aa               , ice%wGHF_Aa               )
    
    ! Isotopes
    CALL allocate_shared_dp_2D(        grid%ny  , grid%nx  , ice%Hi_Aa_prev           , ice%wHi_Aa_prev           )
    
  END SUBROUTINE allocate_ice_model
  
  ! Analytical solution by Schoof 2006 for the "SSA_icestream" benchmark experiment
  SUBROUTINE SSA_Schoof2006_analytical_solution( tantheta, h0, A_flow, L, m, y, U)
    
    USE parameters_module, ONLY: n_flow, ice_density, grav
      
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: tantheta   ! Surface slope in the x-direction
    REAL(dp),                            INTENT(IN)    :: h0         ! Ice thickness
    REAL(dp),                            INTENT(IN)    :: A_flow     ! Ice flow factor
    REAL(dp),                            INTENT(IN)    :: L          ! Ice stream width
    REAL(dp),                            INTENT(IN)    :: m          ! Ice stream transition exponent (sort of like the abruptness of the transition from high to low basal yield stress)
    REAL(dp),                            INTENT(IN)    :: y          ! y-coordinate
    REAL(dp),                            INTENT(OUT)   :: U          ! Ice velocity in the x-direction
    
    ! Local variables:
    REAL(dp)                                           :: B, f, W, ua, ub, uc, ud, ue
    
    ! Calculate the gravitational driving stress f
    f = ice_density * grav * h0 * tantheta
    
    ! Calculate the ice hardness factor B
    B = A_flow**(-1._dp/n_flow)
    
    ! Calculate the "ice stream half-width" W
    W = L * (m+1._dp)**(1._dp/m)
    
    ! Calculate the analytical solution for u
    ua = -2._dp * f**3 * L**4 / (B**3 * h0**3)
    ub = ( 1._dp / 4._dp                           ) * (   (y/L)**     4._dp  - (m+1._dp)**(       4._dp/m) )
    uc = (-3._dp / ((m+1._dp)    * (      m+4._dp))) * (ABS(y/L)**(  m+4._dp) - (m+1._dp)**(1._dp+(4._dp/m)))
    ud = ( 3._dp / ((m+1._dp)**2 * (2._dp*m+4._dp))) * (ABS(y/L)**(2*m+4._dp) - (m+1._dp)**(2._dp+(4._dp/m)))
    ue = (-1._dp / ((m+1._dp)**3 * (3._dp*m+4._dp))) * (ABS(y/L)**(3*m+4._dp) - (m+1._dp)**(3._dp+(4._dp/m)))
    u = ua * (ub + uc + ud + ue)
    
    ! Outside the ice-stream, velocity is zero
    IF (ABS(y) > w) U = 0._dp
    
  END SUBROUTINE SSA_Schoof2006_analytical_solution
  
END MODULE ice_dynamics_module
