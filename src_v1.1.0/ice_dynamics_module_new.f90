MODULE ice_dynamics_module

  USE configuration_module,          ONLY: dp, C                                         
  USE data_types_module,             ONLY: type_grid, type_ice_model, type_init_data_fields, type_SMB_model, type_BMB_model
  USE utilities_module,              ONLY: vertical_integration_from_bottom_to_zeta, vertical_average
  USE derivatives_and_grids_module,  ONLY: map_Acx_to_Aa_2D, map_Acy_to_Aa_2D, Neumann_BC_Aa_3D, ddx_Aa_to_Aa_3D, ddy_Aa_to_Aa_3D, &
                                           ddx_Aa_to_Acx_2D, ddy_Aa_to_Acx_2D, ddx_Aa_to_Acy_2D, ddy_Aa_to_Acy_2D, &
                                           ddx_Acx_to_Aa_2D, ddy_Acy_to_Aa_2D, map_Aa_to_Acx_2D, map_Aa_to_Acy_2D, &
                                           ddx_Acx_to_Ab_2D, ddy_Acx_to_Ab_2D, ddx_Acy_to_Ab_2D, ddy_Acy_to_Ab_2D, &
                                           map_Acx_to_Ab_2D, map_Acy_to_Ab_2D, ddy_Acx_to_Aa_2D, ddx_Acy_to_Aa_2D
  USE general_ice_model_data_module, ONLY: basal_yield_stress

  IMPLICIT NONE
  
CONTAINS   

  ! Update ice thickness from ice dynamic changes and SMB
  SUBROUTINE calculate_ice_thickness_change( grid, ice, SMB, BMB, dt)
    ! Use the total ice velocities to update the ice thickness
    
    USE parameters_module, ONLY: ice_density
    
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
    REAL(dp), DIMENSION( grid%ny, grid%nx)             :: dVi_MB
    
    IF (C%do_benchmark_experiment .AND. C%choice_benchmark_experiment == 'SSA_icestream') THEN
      ! Don't change ice thickness in this experiment, we're only interested in ice velocity
      RETURN
    END IF
            
    ! Initialise at zero
    ice%Qx_Acx = 0._dp
    ice%Qy_Acy = 0._dp
        
    ! Calculate ice fluxes on the Ac grids, with ice thickness
    ! defined on the Aa grid in the upwind direction 
    ! ========================================================
    
    DO i = 1, grid%nx-1
    DO j = 1, grid%ny
      
      IF (ice%U_vav_Acx( j,i) > 0._dp) THEN
        ! Ice moves from left to right, use left-hand ice thickness
        ice%Qx_Acx( j,i) = ice%U_vav_Acx( j,i) * ice%Hi_Aa( j,i  ) * grid%dy * dt
      ELSE
        ! Ice moves from right to left, use right-hand ice thickness
        ice%Qx_Acx( j,i) = ice%U_vav_Acx( j,i) * ice%Hi_Aa( j,i+1) * grid%dy * dt
      END IF
      
    END DO
    END DO
    
    DO i = 1, grid%nx
    DO j = 1, grid%ny-1
      
      IF (ice%V_vav_Acy( j,i) > 0._dp) THEN
        ! Ice moves from bottom to top, use bottom ice thickness
        ice%Qy_Acy( j,i) = ice%V_vav_Acy( j,i) * ice%Hi_Aa( j,i  ) * grid%dx * dt
      ELSE
        ! Ice moves from top to bottom, use top ice thickness
        ice%Qy_Acy( j,i) = ice%V_vav_Acy( j,i) * ice%Hi_Aa( j+1,i) * grid%dx * dt
      END IF
      
    END DO
    END DO
    
    ! Correct outfluxes for possible resulting negative ice thicknesses
    ! =================================================================
    
    ! Ice volume added to each grid cell through the (surface + basal) mass balance
    dVi_MB = (SMB%SMB_year + BMB%BMB) * grid%dx * grid%dy * dt
    
    DO i = 1, grid%nx
    DO j = 1, grid%ny
    
      ! Check how much ice is available for melting or removing (in m^3)
      Vi_available = ice%Hi_Aa( j,i) * grid%dx * grid%dy
      
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
    
    ! Calculate change in ice thickness over time at every vertex
    ! ===========================================================
    
    ice%dHi_dt_Aa = dVi_MB / (grid%dx * grid%dy * dt) ! m/y
    
    DO i = 1, grid%nx-1
    DO j = 1, grid%ny
      ice%dHi_dt_Aa( j,i+1) = ice%dHi_dt_Aa( j,i+1) + ice%Qx_Acx( j,i) / (grid%dx * grid%dy * dt)
      ice%dHi_dt_Aa( j,i  ) = ice%dHi_dt_Aa( j,i  ) - ice%Qx_Acx( j,i) / (grid%dx * grid%dy * dt)
    END DO
    END DO
    DO i = 1, grid%nx
    DO j = 1, grid%ny-1
      ice%dHi_dt_Aa( j+1,i) = ice%dHi_dt_Aa( j+1,i) + ice%Qy_Acy( j,i) / (grid%dx * grid%dy * dt)
      ice%dHi_dt_Aa( j  ,i) = ice%dHi_dt_Aa( j  ,i) - ice%Qy_Acy( j,i) / (grid%dx * grid%dy * dt)
    END DO
    END DO
    
    ! Manually set to zero in first timestep (shouldnt be needed as velocities should be zero, but you never know...)
    IF (dt == 0._dp) ice%dHi_dt_Aa = 0._dp
    
    ! Add dynamic ice thickness change and SMB to update the ice thickness
    ! ====================================================================
    
    ice%Hi_Aa = ice%Hi_Aa + ice%dHi_dt_Aa * dt

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
        ice%Hi_Aa( 1      ,:      ) = 0._dp
        ice%Hi_Aa( grid%ny,:      ) = 0._dp
        ice%Hi_Aa( :      ,1      ) = 0._dp
        ice%Hi_Aa( :      ,grid%nx) = 0._dp
        
      ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
        ! No exception here, as we already exited the routine at the top.
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        
        ! Create a nice circular ice shelf
        DO i = 1, grid%nx
        DO j = 1, grid%ny
          IF (SQRT(grid%x(i)**2+grid%y(j)**2) > grid%xmax * 0.95_dp) ice%Hi_Aa( j,i) = 0._dp
        END DO
        END DO
        
      ELSE
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in calculate_ice_thickness_change!'
        STOP
      END IF
      
    ELSE
          
      ! Apply boundary conditions: set ice thickness to zero at the domain boundary
      ice%Hi_Aa( 1      ,:      ) = 0._dp
      ice%Hi_Aa( grid%ny,:      ) = 0._dp
      ice%Hi_Aa( :      ,1      ) = 0._dp
      ice%Hi_Aa( :      ,grid%nx) = 0._dp
      
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  END SUBROUTINE calculate_ice_thickness_change
  
  
  SUBROUTINE ice_thickness( grid, dt, SL_Aa, Hi_Aa, Hb_Aa, MB, Us_Acx, Vs_Acy, A_flow_mean_Aa, mask_GL_Aa, mask_shelf_Aa, &
      Hi_new, MB_applied, domain_leaving_flux, mass_conservation, buttressing_term_x, buttressing_term_y, dHi_D)
    ! Update the ice thickness for all points
    ! Heiko Goelzer (h.goelzer@uu.nl) 2016

   ! USE configuration_module, ONLY: dp, C, G
   ! USE utilities_module, ONLY: stag_2D_ab
   ! USE schoof_grounding_line_flux_module, ONLY: schoof_grounding_line_flux, shelf_buttressing
   ! USE schoof_grounding_line_flux_module, ONLY: schoof_grounding_line_flux_down_stream_buttressing_stress
    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp),                             INTENT(IN)   :: dt
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(IN)  :: SL_Aa                   ! In meter relative to present day eustatic sea level
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(IN)  :: Hi_Aa                      ! In meter
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(IN)  :: Hb_Aa                      ! In meter
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(IN)  :: MB                      ! In meter ice equivalent per year
    REAL(dp), DIMENSION(grid%ny,grid%nx-1), INTENT(IN)  :: Us_Acx                      ! Velocity in x-direction [m yr^-1] on staggered Ac-x grid
    REAL(dp), DIMENSION(grid%ny-1,grid%nx), INTENT(IN)  :: Vs_Acy                      ! Velocity in y-direction [m yr^-1] on staggered Ac-y grid
    ! needed for schoof
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(IN)  :: A_flow_mean_Aa             ! A_flow on grid Aa
    INTEGER, DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: mask_GL_Aa
    INTEGER, DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: mask_shelf_Aa

    ! Output variables:
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(OUT) :: Hi_new                  ! The new calculated ice thickness [m]
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(OUT) :: MB_applied              ! Applied mass balance [meter ice equivalent per year]
    REAL(dp),                       INTENT(OUT) :: domain_leaving_flux     ! residual claving flux
    REAL(dp),                       INTENT(OUT) :: mass_conservation       ! In case all ice volume is conserved this variable is zero
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(OUT) :: buttressing_term_x      ! Buttressing for schoof (Aa)
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(OUT) :: buttressing_term_y      ! Buttressing for schoof (Aa)
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(OUT) :: dHi_D                   ! Dynamic thickness change (Aa)

    ! Local variables:
    INTEGER                                     :: i, j
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: dHi                     ! Ice thickness change due to fluxes divergence and MB, in meter ice equivalent for this G%dt 
    REAL(dp)                                    :: dHi_available           ! In meter ice equivalent for this G%dt 
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: MB_remaining            ! Remaining mass balance at first stage [meter ice equivalent for this G%dt]

    ! MacAyeal
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: flux_xr                 ! The flux in the x-direction to the right
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: flux_xl                 ! The flux in the x-direction to the left
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: flux_yu                 ! The flux in the y-direction upwards
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: flux_yd                 ! The flux in the y-direction downwards
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: Hi_Ab                   ! Ice thickness on staggered grid Ab
    
    LOGICAL,  PARAMETER :: apply_remaining_mass_balance = .FALSE.
    REAL(dp), PARAMETER :: maximum_dHi_change_per_dt = 100._dp

    buttressing_term_x = 0._dp
    buttressing_term_y = 0._dp

    dHi = 0._dp
    dHi_D = 0._dp

    CALL stag_2D_ab( grid, Hi_Aa, Hi_Ab)

    ! Solve diffusion equation

    ! Formulate as function of velocity, flux_xr and flux_xl on Ac_x, flux_yu and flux_yd on Ac_y
!$OMP Parallel Do private(i,j)
    DO i = 2,grid%nx-1
    DO j = 2,grid%ny-1
      !flux_xr(j,i) = Us_Acx(j  ,i  ) * 0.5_dp * (Hi_Ab(j  ,i  ) + Hi_Ab(j-1,i  ))
      !flux_xl(j,i) = Us_Acx(j  ,i-1) * 0.5_dp * (Hi_Ab(j  ,i-1) + Hi_Ab(j-1,i-1))
      !flux_yu(j,i) = Vs_Acy(j  ,i  ) * 0.5_dp * (Hi_Ab(j  ,i  ) + Hi_Ab(j  ,i-1))
      !flux_yd(j,i) = Vs_Acy(j-1,i  ) * 0.5_dp * (Hi_Ab(j-1,i  ) + Hi_Ab(j-1,i-1))
      IF (Us_Acx( j  ,i  ) > 0._dp) THEN
        flux_xr( j,i) = Us_Acx( j  ,i  ) * Hi_Aa( j  ,i  )
      ELSE
        flux_xr( j,i) = Us_Acx( j  ,i  ) * Hi_Aa( j  ,i+1)
      END IF
      IF (Us_Acx( j  ,i-1) > 0._dp) THEN
        flux_xl( j,i) = Us_Acx( j  ,i-1) * Hi_Aa( j  ,i-1)
      ELSE
        flux_xl( j,i) = Us_Acx( j  ,i-1) * Hi_Aa( j  ,i  )
      END IF
      IF (Vs_Acy( j  ,i  ) > 0._dp) THEN
        flux_yu( j,i) = Vs_Acy( j  ,i  ) * Hi_Aa( j  ,i  )
      ELSE
        flux_yu( j,i) = Vs_Acy( j  ,i  ) * Hi_Aa( j+1,i  )
      END IF
      IF (Vs_Acy( j-1,i  ) > 0._dp) THEN
        flux_yd( j,i) = Vs_Acy( j-1,i  ) * Hi_Aa( j-1,i  )
      ELSE
        flux_yd( j,i) = Vs_Acy( j-1,i  ) * Hi_Aa( j  ,i  )
      END IF
    END DO
    END DO
!$OMP END PARALLEL DO

!    IF(C%apply_schoof_grounding_line_flux) THEN
!     ! Apply Schoof flux boundary condition on grounding line
!     
!     IF(C%use_down_stream_buttressing_stress) THEN
!      ! Calculating the velocity gradients which are used to compute the buttressing stress on the shelf near the grounding line points
!      ! in a down-stream manner (the staggering might be not conform the default staggering approach, so a check on that is needed):
!      ! In/Output: flux_xr, flux_xl, flux_yu, flux_yd; Output: buttressing_term_x, buttressing_term_y
!      CALL schoof_grounding_line_flux_down_stream_buttressing_stress(sealevel, Hi, Hb, A_flow_mean, mask, Us_Ac, Vs_Ac, flux_xr, flux_xl, flux_yu, flux_yd, buttressing_term_x, buttressing_term_y)
!      CALL schoof_grounding_line_flux_down_stream_buttressing_stress( grid, SL_Aa, Hi_Aa, Hb_Aa, A_flow_mean_Aa, mask_GL_Aa, mask_shelf_Aa, Us_Acx, Vs_Acy, flux_xr, flux_xl, flux_yu, flux_yd, buttressing_term_x, buttressing_term_y)
!     ELSE
!      ! Calculate buttressing term on Aa
!      CALL shelf_buttressing(Hi, A_flow_mean, mask, Us_Ac, Vs_Ac, buttressing_term_x, buttressing_term_y)

      !CALL shelf_buttressing( grid, Hi_Aa, A_flow_mean_Aa, mask_shelf_Aa, Us_Acx, Vs_Acy, buttressing_term_x, buttressing_term_y)
      
!
!      ! In/Output: flux_xr, flux_xl, flux_yu, flux_yd
!      CALL schoof_grounding_line_flux(sealevel, Hi, Hb, A_flow_mean, mask, buttressing_term_x, buttressing_term_y, flux_xr, flux_xl, flux_yu, flux_yd)
!      CALL schoof_grounding_line_flux( grid, SL_Aa, Hi_Aa, Hb_Aa, A_flow_mean_Aa, mask_GL_Aa, mask_shelf_Aa, buttressing_term_x, buttressing_term_y, flux_xr, flux_xl, flux_yu, flux_yd)
!     END IF
!    END IF

    ! Calculate flux divergence on grid Aa
!$OMP Parallel Do private(i,j)
    DO i = 2,grid%nx-1
    DO j = 2,grid%ny-1
      dHi_D(j,i) = (dt / grid%dx) * (flux_xl(j,i) - flux_xr(j,i)) + (dt / grid%dy) * (flux_yd(j,i) - flux_yu(j,i))
    END DO
    END DO
!$OMP END PARALLEL DO

!$OMP Parallel Do private(i,j,dHi_available)
    DO i = 1, grid%nx
    DO j = 1, grid%ny

      ! The maximum ice thickness (for this j,i point) available for melting
      dHi_available      = Hi_Aa(j,i)

      ! The mass balance
      MB_applied(j,i)   = MB(j,i) * dt

      ! The remaining mass balance, the part which is not directly applied (possibly in a later stage it is):
      MB_remaining(j,i) = MIN(dHi_available + MB_applied(j,i), 0._dp)

      ! Set the applied mass balance equal to the maximum available amount of ice (or allowed ice change):
      IF(MB_remaining(j,i) < 0._dp) MB_applied(j,i) = - dHi_available 

      ! continuity equation
      dHi(j,i) = MB_applied(j,i) + dHi_D(j,i)

    END DO
    END DO
!$OMP END PARALLEL DO
    
    IF(apply_remaining_mass_balance) THEN
     ! In some points for which the MB_remaining /= 0 there might be again some ice layer
     ! after the flux contributions are distributed from neighbour grid points. For those cases the option
     ! is available to apply the MB_remaining to this renewed available ice volume:
!$OMP Parallel Do private(i,j)
     DO i = 1, grid%nx
     DO j = 1, grid%ny
       IF(MB_remaining(j,i) < 0._dp) THEN
        ! The part of MB_remaining which can be applied is determined:
        MB_remaining(j,i) = MAX(MB_remaining(j,i), - MAX(dHi(j,i), 0._dp))
        dHi(j,i)          = dHi(j,i)        + MB_remaining(j,i)
        MB_applied(j,i)   = MB_applied(j,i) + MB_remaining(j,i)
       END IF
     END DO
     END DO
!$OMP END PARALLEL DO
    END IF

    ! Expressing the mass balance again in meter ice equivalent per year for diagnostics:
    MB_applied = MB_applied / dt
    
    ! Manually set to zero in first timestep (shouldnt be needed as velocities should be zero, but you never know...)
    IF (dt == 0._dp) dHi = 0._dp

    ! Update Hi_new:
    Hi_new = Hi_Aa + dHi

    ! Boundary conditions
    !    IF(C%eismint_case .AND. (C%choice_eismint_experiment == 'P' .OR. C%choice_eismint_experiment == 'Q' .OR. C%choice_eismint_experiment == 'R')) THEN
     ! set eismint 1 fixed margin exp boundary condition

!    IF(.NOT. C%include_shelf_dynamics) THEN
!     Hi_new(1   ,:) = C%Hi_min
!     Hi_new(grid%ny,:) = C%Hi_min
!     Hi_new(:,1   ) = C%Hi_min
!     Hi_new(:,grid%nx) = C%Hi_min
!    ELSE
     ! left
     Hi_new(2:grid%ny-1,1   ) = 2._dp * Hi_new(2:grid%ny-1,2)-Hi_new(2:grid%ny-1,3)
     ! right
     Hi_new(2:grid%ny-1,grid%nx) = 2._dp * Hi_new(2:grid%ny-1,grid%nx-1)-Hi_new(2:grid%ny-1,grid%nx-2)
     ! bottom
     Hi_new(1   ,2:grid%nx-1) = 2._dp * Hi_new(2,2:grid%nx-1)-Hi_new(3,2:grid%nx-1)
     ! top
     Hi_new(grid%ny,2:grid%nx-1) = 2._dp * Hi_new(grid%ny-1,2:grid%nx-1)-Hi_new(grid%ny-2,2:grid%nx-1)

     ! Additional check on gradient. No increase in northward direction
     Hi_new(2:grid%ny-1,1) = MIN(Hi_new(2:grid%ny-1,1),Hi_new(2:grid%ny-1,2))
     Hi_new(2:grid%ny-1,grid%nx) = MIN(Hi_new(2:grid%ny-1,grid%nx),Hi_new(2:grid%ny-1,grid%nx-1))
     Hi_new(1,2:grid%nx-1) = MIN(Hi_new(1,2:grid%nx-1),Hi_new(2,2:grid%nx-1))
     Hi_new(grid%ny,2:grid%nx-1) = MIN(Hi_new(grid%ny,2:grid%nx-1),Hi_new(grid%ny-1,2:grid%nx-1))
!    END IF

    !    END IF

    ! Check Hi won't be below C%Hi_min:
!    IF(C%choice_write_message) THEN
!     IF(ANY(Hi_new < C%Hi_min - 1.0E-15_dp)) WRITE(UNIT=C%stdlog, FMT='(2(A, F20.16))') 'ALARM: Some Hi < C%Hi_min, should not be the case!, Minimum Hi = ', MINVAL(Hi_new), ' at time ', G%dt
!    END IF
!    Hi_new = MAX(Hi_new, C%Hi_min)

    ! Corners
    ! BL
    Hi_new(1,1) = Hi_new(1,2)+Hi_new(2,1)-Hi_new(2,2)
    ! BR
    Hi_new(1,grid%nx) = Hi_new(1,grid%nx-1)+Hi_new(2,grid%nx)-Hi_new(2,grid%nx-1)
    ! TL
    Hi_new(grid%ny,1) = Hi_new(grid%ny,2)+Hi_new(grid%ny-1,1)-Hi_new(grid%ny-1,2)
    ! BL
    Hi_new(grid%ny,grid%nx) = Hi_new(grid%ny,grid%nx-1)+Hi_new(grid%ny-1,grid%nx)-Hi_new(grid%ny-1,grid%nx-1)

!    ! Check Hi won't be below C%Hi_min:
!    IF(C%choice_write_message) THEN
!     IF(ANY(Hi_new < C%Hi_min - 1.0E-15_dp)) WRITE(UNIT=C%stdlog, FMT='(2(A, F20.16))') 'ALARM: Some Hi < C%Hi_min, should not be the case!, Minimum Hi = ', MINVAL(Hi_new), ' at time ', G%dt
!    END IF
!    Hi_new = MAX(Hi_new, C%Hi_min)

    ! The total ice thickness change and applied mass balance are compared, to give the residual calving flux 
    domain_leaving_flux = SUM(Hi_new - Hi_Aa - MB_applied * dt) 
    mass_conservation = 0._dp

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
        Hi_new( 1      ,:      ) = 0._dp
        Hi_new( grid%ny,:      ) = 0._dp
        Hi_new( :      ,1      ) = 0._dp
        Hi_new( :      ,grid%nx) = 0._dp
        
      ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
        ! No exception here, as we already exited the routine at the top.
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        
        ! Create a nice circular ice shelf
        DO i = 1, grid%nx
        DO j = 1, grid%ny
          IF (SQRT(grid%x(i)**2+grid%y(j)**2) > grid%xmax * 0.95_dp) Hi_new( j,i) = 0._dp
        END DO
        END DO
        
      ELSE
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in calculate_ice_thickness_change!'
        STOP
      END IF
      
    ELSE
          
      ! Apply boundary conditions: set ice thickness to zero at the domain boundary
      Hi_new( 1      ,:      ) = 0._dp
      Hi_new( grid%ny,:      ) = 0._dp
      Hi_new( :      ,1      ) = 0._dp
      Hi_new( :      ,grid%nx) = 0._dp
      
    END IF ! IF (C%do_benchmark_experiment) THEN

  END SUBROUTINE ice_thickness
  SUBROUTINE stag_2D_ab( grid, var_2D_Aa, var_2D_Ab)
    ! This subroutine staggers a variable from grid Aa to Ab
    ! Heiko Goelzer (h.goelzer@uu.nl) Jan 2016

    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(IN)       :: var_2D_Aa
    ! Output variables:  
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(OUT)      :: var_2D_Ab
    ! Local variables:
    INTEGER                                          :: i, j
    
    var_2D_Ab = 0._dp

!$OMP Parallel Do private(i,j)
    DO i = 1, grid%NX-1
    DO j = 1, grid%NY-1
      ! averaging from Aa to grid Ab
      var_2D_Ab(j,i) = 0.25_dp * ( var_2D_Aa(j,i) + var_2D_Aa(j,i+1) + var_2D_Aa(j+1,i+1) + var_2D_Aa(j+1,i) )
    END DO
    END DO
!$OMP END PARALLEL DO
    
  END SUBROUTINE stag_2D_ab
  SUBROUTINE ice_stream_velocity_by_ssa( grid, Hi_Aa, Hi_Ab, Hb_Ab, SL_Ab, mask_sheet_Aa, mask_shelf_Aa, mask_land_Aa, mask_ocean_Aa, dHs_dx_Acx, dHs_dy_Acy, A_flow_mean_Aa, basal_yield_stress_Aa, U_ssa_Ac, V_ssa_Ac)
    ! Calculate staggered (Ac) U_ssa and V_ssa for the entire sheet-shelf system with a hybrid model approach (Bueler&Brown, 2009).
    ! The iterative solver uses the successive over-relaxation (SOR) scheme (See Numerical recepies eq (19.5.29)):
    !
    !   2 U_ssa_xx + 1/2 U_ssa_yy + 3/2 V_ssa_xy = rhs_x = C%ice_density g dHs_dx / C_uv
    !   2 V_ssa_yy + 1/2 V_ssa_xx + 3/2 U_ssa_xy = rhs_y = C%ice_density g dHs_dy / C_uv
    ! with:
    !                2         2                                      2               (1-n)/2n
    ! C_uv = ( dU_ssa_dx  + dV_ssa_dy  + dU_ssa_dx dV_ssa_dy + 1/4 (dU_ssa_dy+dV_ssa_dx) + epsilon_sq_0)        (see Huybrechts 4.28)
    !
    ! The Red-Black Gauss-Seidel + over relaxation iterative method is used. At entrance the velocity
    ! field of the previous time step serves as an initial guess for the velocity fields. At exit it is
    ! the calculated solution. At grounded points the values of U_ssa and V_ssa are not changed.
    ! Heiko Goelzer (h.goelzer@uu.nl) 2017
    
    USE parameters_module,          ONLY: n_flow, ice_density, seawater_density, grav
    
    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(IN)    :: Hi_Aa
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1), INTENT(IN)    :: Hi_Ab
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1), INTENT(IN)    :: Hb_Ab
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1), INTENT(IN)    :: SL_Ab
    INTEGER,  DIMENSION(grid%ny,grid%nx), INTENT(IN)    :: mask_sheet_Aa, mask_shelf_Aa, mask_land_Aa, mask_ocean_Aa
    REAL(dp), DIMENSION(grid%ny,grid%nx-1), INTENT(IN)    :: dHs_dx_Acx
    REAL(dp), DIMENSION(grid%ny-1,grid%nx), INTENT(IN)    :: dHs_dy_Acy
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(IN)    :: A_flow_mean_Aa                  ! Vertical mean of the flow parameter A
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(IN)    :: basal_yield_stress_Aa           ! Basal yield stress (Mohr-Coulomb model) (eqn 13, de Boer, Clim. Dyn., 2013)

    ! Output variables:
    REAL(dp), DIMENSION(grid%ny,grid%nx-1), INTENT(INOUT) :: U_ssa_Ac                        ! SSA velocity field on grid Ac_x
    REAL(dp), DIMENSION(grid%ny-1,grid%nx), INTENT(INOUT) :: V_ssa_Ac                        ! SSA velocity field on grid Ac_y
    
    REAL(dp), DIMENSION(grid%ny,grid%nx-1) :: mask_shelf_Acx
    REAL(dp), DIMENSION(grid%ny-1,grid%nx) :: mask_shelf_Acy
    
    REAL(dp), DIMENSION(grid%ny,grid%nx-1) :: residue_Ac_x                    ! The SOR residue of the shelf velocity equation in the x-direction: 4 U_ssa_xx + U_ssa_yy + 3 V_ssa_xy = 2 rhs_x  (See Numerical recepies eq (19.5.28))
    REAL(dp), DIMENSION(grid%ny-1,grid%nx) :: residue_Ac_y                    ! The SOR residue of the shelf velocity equation in the y-direction: 4 V_ssa_yy + V_ssa_xx + 3 U_ssa_xy = 2 rhs_y  (See Numerical recepies eq (19.5.28))
    INTEGER    :: velocity_iterations          ! Only for printing to recording file

    ! Local variables:
    INTEGER                                       :: i, j
    INTEGER                                       :: explicit_update_iteration    ! The counter for the outermost loop in which the fields are updated which are used explicitly
    INTEGER                                       :: iter_gauss_seidel            ! Counter over the Gauss-Seidel iterations
    INTEGER                                       :: sweep                        ! To interchange between the 'red' and the 'black' grid points in the Red-Black Gauss-Seidel sheme
    INTEGER                                       :: s0, s1, sd
    INTEGER                                       :: start_odd_or_even
    REAL(dp)                                      :: largest_rest
    REAL(dp), DIMENSION(grid%ny,grid%nx)          :: eu_ij                        ! The coefficient of U_ij
    REAL(dp), DIMENSION(grid%ny,grid%nx)          :: ev_ij                        ! The coefficient of V_ij
    REAL(dp)                                      :: U_ssa_xy
    REAL(dp)                                      :: V_ssa_xy
    REAL(dp), DIMENSION(grid%ny,grid%nx)          :: rhs_x                        ! right hand side of the velocity equation in the x-direction
    REAL(dp), DIMENSION(grid%ny,grid%nx)          :: rhs_y                        ! right hand side of the velocity equation in the y-direction
    REAL(dp)                                      :: dU_ssa_dx                    ! Derivative in the x-direction of the U_ssa field
    REAL(dp)                                      :: dU_ssa_dy                    ! Derivative in the y-direction of the U_ssa field
    REAL(dp)                                      :: dV_ssa_dx                    ! Derivative in the x-direction of the V_ssa field
    REAL(dp)                                      :: dV_ssa_dy                    ! Derivative in the y-direction of the V_ssa field
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1)                :: C_uv
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1)                :: beta_base                    ! basal shear stress coefficient (plastic till model) (eqn (27) B&B,09)
    
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1)                :: A_flow_mean_Ab                  ! Vertical mean of the flow parameter A
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1)                :: phi_fric_Ab
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1)                :: basal_yield_stress_Ab           ! Basal yield stress (Mohr-Coulomb model) (eqn 13, de Boer, Clim. Dyn., 2013)
    
    INTEGER,  PARAMETER :: number_of_explicit_shelf_iterations    = 100 
    INTEGER,  PARAMETER :: maximum_number_gauss_seidel_iterations = 5000
    REAL(dp), PARAMETER :: SOR_omega                              = 1.2_dp
    REAL(dp), PARAMETER :: max_residual_U                         = 2.5_dp
    REAL(dp), PARAMETER :: epsilon_sq_0         = 1E-12_dp
    REAL(dp), PARAMETER :: dummy_delta_v        = 1E-3_dp
    REAL(dp), PARAMETER :: q_plastic            = 0.30_dp       ! Parameter used for basal stress (inverse of m_flow)
    REAL(dp), PARAMETER :: u_threshold          = 100._dp       ! scaling of tau_yield to get the correct unit (function of q_plastic)
    
    ! Map A_flow to Ab grid
    CALL stag_2D_ab( grid, A_flow_mean_Aa, A_flow_mean_Ab)
    
    ! Calculate basal yield stress on Ab grid
    CALL calculate_basal_yield_stress_Ab( grid, Hi_Ab, Hb_Ab, SL_Ab, phi_fric_Ab, basal_yield_stress_Ab)
    !CALL stag_2D_ab( grid, basal_yield_stress_Aa, basal_yield_stress_Ab)
    
    ! Get "shelf mask" (according to Heiko's definition) on staggered grid
    mask_shelf_Acx = 0
    DO i = 1, grid%nx-1
    DO j = 1, grid%ny
      IF ((mask_shelf_Aa( j,i) == 1 .OR. mask_ocean_Aa( j,i) == 1) .AND. (mask_shelf_Aa( j,i+1) == 1 .OR. mask_ocean_Aa( j,i+1) == 1)) mask_shelf_Acx( j,i) = 1
    END DO
    END DO
    DO i = 1, grid%nx
    DO j = 1, grid%ny-1
      IF ((mask_shelf_Aa( j,i) == 1 .OR. mask_ocean_Aa( j,i) == 1) .AND. (mask_shelf_Aa( j+1,i) == 1 .OR. mask_ocean_Aa( j+1,i) == 1)) mask_shelf_Acy( j,i) = 1
    END DO
    END DO

    C_uv = 0.0_dp
    
    velocity_iterations = 0

    ! The fields U_ssa and V_ssa are used in an explicit way to calcute the expression for C_uv. In this
    ! loop the fields U_ssa and V_ssa are updated, and thereafter the kernel iteration is repeated:
    DO explicit_update_iteration = 1, number_of_explicit_shelf_iterations

      ! Calculate the non-linear term C_uv on grid Ab 'fixed' for the current solution {U_ssa, V_ssa}
      ! Could replace DO loops by WHERE statements with appropriate masks  
      DO i = 2, grid%nx - 2
        DO j = 2, grid%ny - 2
          ! Gradients 2 Ac to Ac, average to Ab
          dU_ssa_dx = (U_ssa_Ac(j  ,i+1) - U_ssa_Ac(j  ,i-1) + U_ssa_Ac(j+1,i+1) - U_ssa_Ac(j+1,i-1)) / (4._dp * grid%dx) 
          dV_ssa_dy = (V_ssa_Ac(j+1,i  ) - V_ssa_Ac(j-1,i  ) + V_ssa_Ac(j+1,i+1) - V_ssa_Ac(j-1,i+1)) / (4._dp * grid%dy) 
          ! Gradients from Ac to Ab
          dU_ssa_dy = (U_ssa_Ac(j+1,i  ) - U_ssa_Ac(j,i)) / (grid%dy)
          dV_ssa_dx = (V_ssa_Ac(j  ,i+1) - V_ssa_Ac(j,i)) / (grid%dx)
          ! Common term for both directions
          C_uv(j,i) = A_flow_mean_Ab(j,i)**(-1._dp / n_flow) * (dU_ssa_dx**2 + dV_ssa_dy**2 + dU_ssa_dx * dV_ssa_dy &
               + 0.25_dp * (dU_ssa_dy + dV_ssa_dx)**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow))
          ! Determine beta_base on grid Ab, for basal stress (beta in eqn (27) Bueler&Brown, 2009), here equivalent to the sliding parameter (A**-1/m)
          beta_base(j,i) = basal_yield_stress_Ab(j,i) * ( (dummy_delta_v**2 + (0.5_dp*(U_ssa_Ac(j,i)+U_ssa_Ac(j+1,i)))**2 + (0.5_dp*(V_ssa_Ac(j,i)+V_ssa_Ac(j,i+1)))**2)**(0.5_dp * (q_plastic-1._dp)) ) / (u_threshold**q_plastic)

       END DO
      END DO

      ! Calculate right hand sides 
      DO i = 2, grid%nx - 2
      DO j = 2, grid%ny - 1
        ! Calculate the right hand sides: rhs_x on Ac_x
        rhs_x(j,i) =  ice_density * grav * dHs_dx_Acx(j,i) / (0.5_dp * (C_uv(j,i)+C_uv(j-1,i)))
        ! update coefficients
        IF (mask_shelf_Acx( j,i) == 1) THEN
         ! beta_base is zero on shelf
         eu_ij(j,i) = -(8._dp / (grid%dx**2) + 2._dp / (grid%dy**2))
        ELSE
         ! On sheet, coefficients on Ac_x now dependent on beta_base/Hi, because the whole equation is divided by Hi.
         ! Factors 0.5 in numerator and denominator in beta_base/Hi from averaging beta_base (Ab) and Hi (Aa) to Ac falls out.  
         eu_ij(j,i) = -(8._dp / (grid%dx**2) + 2._dp / (grid%dy**2)) - 2._dp * (beta_base(j,i)+beta_base(j-1,i))/(Hi_Aa(j,i)+Hi_Aa(j,i+1)) / (0.5_dp * (C_uv(j,i)+C_uv(j-1,i)))
        END IF
      END DO
      END DO
      DO i = 2, grid%nx - 1
      DO j = 2, grid%ny - 2
        ! Calculate the right hand sides: rhs_y on Ac_y
        rhs_y(j,i) =  ice_density * grav * dHs_dy_Acy(j,i) / (0.5_dp * (C_uv(j,i)+C_uv(j,i-1))) 
        ! update coefficients
        IF (mask_shelf_Acy( j,i) == 1) THEN
         ! beta_base is zero on shelf
         ev_ij(j,i) = -(2._dp / (grid%dx**2) + 8._dp / (grid%dy**2))
         ELSE
         ! On sheet, coefficients on Ac_x now dependent on beta_base/Hi, because the whole equation is divided by Hi.
         ! Factors 0.5 in numerator and denominator in beta_base/Hi from averaging beta_base (Ab) and Hi (Aa) to Ac falls out.  
         ev_ij(j,i) = -(2._dp / (grid%dx**2) + 8._dp / (grid%dy**2)) - 2._dp * (beta_base(j,i)+beta_base(j,i-1))/(Hi_Aa(j,i)+Hi_Aa(j+1,i)) / (0.5_dp * (C_uv(j,i)+C_uv(j,i-1))) 
        END IF
      END DO
      END DO

      iterate: DO iter_gauss_seidel = 1, maximum_number_gauss_seidel_iterations
        ! Calculate the new U_ssa, and V_ssa at shelf points. First sweep: j,i starts at (2,2) and is even; second sweep: j,i starts at (3,3) and is odd:

        largest_rest = 0._dp

        residue_Ac_x = 0._dp
        residue_Ac_y = 0._dp

        s0 = 0
        s1 = 1
        sd = s1-s0

        ! Red black scheme:
        DO sweep = s0, s1, sd
          DO j = 2, grid%ny - 1
            start_odd_or_even = 2 + MOD(j + sweep, 2)
            DO i = start_odd_or_even, grid%nx - 1, 2

              ! On grid Ac staggered in x for residue_x 
              V_ssa_xy = ((V_ssa_Ac(j,i+1)- V_ssa_Ac(j,i)) - (V_ssa_Ac(j-1,i+1)-V_ssa_Ac(j-1,i))) / (4._dp * grid%dx * grid%dy) 
              ! On grid Ac staggered in y for residue_y 
              U_ssa_xy = ((U_ssa_Ac(j+1,i)- U_ssa_Ac(j,i)) - (U_ssa_Ac(j+1,i-1)-U_ssa_Ac(j,i-1))) / (4._dp * grid%dx * grid%dy) 
              
              ! The first four terms of residue_x and residue_y use 'old' values at the first sweep and 'new' values at the second sweep (the red-black effect)
              IF(i<=grid%nx-2) THEN ! Ac_x limited in x but not in y
               residue_Ac_x(j,i) = (4._dp * (U_ssa_Ac(j,i+1) + U_ssa_Ac(j,i-1)) / (grid%dx**2) + (U_ssa_Ac(j+1,i) + U_ssa_Ac(j-1,i)) / (grid%dy**2) + eu_ij(j,i) * U_ssa_Ac(j,i) + 3._dp * V_ssa_xy - 2._dp * rhs_x(j,i)) / eu_ij(j,i)  
               U_ssa_Ac(j,i) = U_ssa_Ac(j,i) - SOR_omega * residue_Ac_x(j,i)
              END IF
              IF(j<=grid%ny-2) THEN ! Ac_y limited in y but not in x
               residue_Ac_y(j,i) = (4._dp * (V_ssa_Ac(j+1,i) + V_ssa_Ac(j-1,i)) / (grid%dy**2) + (V_ssa_Ac(j,i+1) + V_ssa_Ac(j,i-1)) / (grid%dx**2) + ev_ij(j,i) * V_ssa_Ac(j,i) + 3._dp * U_ssa_xy - 2._dp * rhs_y(j,i)) / ev_ij(j,i)  
               V_ssa_Ac(j,i) = V_ssa_Ac(j,i) - SOR_omega * residue_Ac_y(j,i)
              END IF
              
              ! Keeping the largest rest for the iteration-stop-check:
              largest_rest = MAX(largest_rest, ABS(residue_Ac_x(j,i)), ABS(residue_Ac_y(j,i)))

            END DO
          END DO
        END DO ! end of the sweep loop

        ! Neumann BC conditions at the domain borders:
        ! On Ac_x grid j=1:NY, i=1:(NX-1); On Ac_y grid j=1:(NY-1), i=1:NX; 
        U_ssa_Ac(   1, 2:grid%nx - 2) = U_ssa_Ac(       2, 2:grid%nx - 2)
        U_ssa_Ac(grid%ny, 2:grid%nx - 2) = U_ssa_Ac(grid%ny - 1, 2:grid%nx - 2)
        V_ssa_Ac(   1, 2:grid%nx - 1) = V_ssa_Ac(       2, 2:grid%nx - 1)
        V_ssa_Ac(grid%ny - 1, 2:grid%nx - 1) = V_ssa_Ac(grid%ny - 2, 2:grid%nx - 1)

        U_ssa_Ac(:,    1) = U_ssa_Ac(:,        2)
        U_ssa_Ac(:, grid%nx - 1) = U_ssa_Ac(:, grid%nx - 2)
        V_ssa_Ac(1:grid%ny - 1,    1) = V_ssa_Ac(1:grid%ny - 1,        2)
        V_ssa_Ac(1:grid%ny - 1, grid%nx) = V_ssa_Ac(1:grid%ny - 1, grid%nx - 1)

!        ! Huybrechts boundary condition (cf. Huybrechts p99/100)
!        U_ssa_Ac(2:grid%ny-1,1     ) = U_ssa_Ac(2:grid%ny-1,2     ) - 3._dp * grid%dx * A_flow_mean_Ab(2:grid%ny-1,1     ) * (ice_density * &
!             grav * (1._dp - ice_density / seawater_density) * (Hi_Aa(2:grid%ny-1,1     ) + Hi_Aa(2:grid%ny-1,2     )) / 12._dp)**3._dp
!        U_ssa_Ac(2:grid%ny-1,grid%nx-1) = U_ssa_Ac(2:grid%ny-1,grid%nx-2) + 3._dp * grid%dx * A_flow_mean_Ab(2:grid%ny-1,grid%nx-1) * (ice_density * &
!             grav * (1._dp - ice_density / seawater_density) * (Hi_Aa(2:grid%ny-1,grid%nx-1) + Hi_Aa(2:grid%ny-1,grid%nx-2)) / 12._dp)**3._dp
!        V_ssa_Ac(     1,2:grid%nx-1) = V_ssa_Ac(     2,2:grid%nx-1) - 3._dp * grid%dy * A_flow_mean_Ab(     1,2:grid%nx-1) * (ice_density * &
!             grav * (1._dp - ice_density / seawater_density) * (Hi_Aa(     1,2:grid%nx-1) + Hi_Aa(     2,2:grid%nx-1)) / 12._dp)**3._dp
!        V_ssa_Ac(grid%ny-1,2:grid%nx-1) = V_ssa_Ac(grid%ny-2,2:grid%nx-1) + 3._dp * grid%dy * A_flow_mean_Ab(grid%ny-1,2:grid%nx-1) * (ice_density * &
!             grav * (1._dp - ice_density / seawater_density) * (Hi_Aa(grid%ny-1,2:grid%nx-1) + Hi_Aa(grid%ny-2,2:grid%nx-1)) / 12._dp)**3._dp
!
!        U_ssa_Ac(   1,2:grid%nx-2) = U_ssa_Ac(     2,2:grid%nx-2) + 0.25_dp * grid%dx / grid%dy * (V_ssa_Ac(   1,3:grid%nx-1) - V_ssa_Ac(   1,1:grid%nx-3) + V_ssa_Ac(     2,3:grid%nx-1) - V_ssa_Ac(     2,1:grid%nx-3))
!        U_ssa_Ac(grid%ny,2:grid%nx-2) = U_ssa_Ac(grid%ny-1,2:grid%nx-2) - 0.25_dp * grid%dx / grid%dy * (V_ssa_Ac(grid%ny,3:grid%nx-1) - V_ssa_Ac(grid%ny,1:grid%nx-3) + V_ssa_Ac(grid%ny-1,3:grid%nx-1) - V_ssa_Ac(grid%ny-1,1:grid%nx-3))
!        V_ssa_Ac(2:grid%ny-2,1   ) = V_ssa_Ac(2:grid%ny-2,2     ) + 0.25_dp * grid%dx / grid%dy * (U_ssa_Ac(3:grid%ny-1,1   ) - U_ssa_Ac(1:grid%ny-3,1   ) + U_ssa_Ac(3:grid%ny-1,2     ) - U_ssa_Ac(1:grid%ny-3,2     ))
!        V_ssa_Ac(2:grid%ny-2,grid%nx) = V_ssa_Ac(2:grid%ny-2,grid%nx-1) - 0.25_dp * grid%dx / grid%dy * (U_ssa_Ac(3:grid%ny-1,grid%nx) - U_ssa_Ac(1:grid%ny-3,grid%nx) + U_ssa_Ac(3:grid%ny-1,grid%nx-1) - U_ssa_Ac(1:grid%ny-3,grid%nx-1))

        ! Corner points from average
        U_ssa_Ac(1   ,1     )= U_ssa_Ac(     2,1     ) + U_ssa_Ac(   1,2     ) - U_ssa_Ac(     2,2     )
        U_ssa_Ac(1   ,grid%nx-1)= U_ssa_Ac(     2,grid%nx-1) + U_ssa_Ac(   1,grid%nx-2) - U_ssa_Ac(     2,grid%nx-2)
        U_ssa_Ac(grid%ny,1     )= U_ssa_Ac(grid%ny-1,1     ) + U_ssa_Ac(grid%ny,2     ) - U_ssa_Ac(grid%ny-1,2     )
        U_ssa_Ac(grid%ny,grid%nx-1)= U_ssa_Ac(grid%ny-1,grid%nx-1) + U_ssa_Ac(grid%ny,grid%nx-2) - U_ssa_Ac(grid%ny-1,grid%nx-2)

        V_ssa_Ac(1   ,1     )= V_ssa_Ac(     2,1   ) + V_ssa_Ac(     1,2     ) - V_ssa_Ac(     2,2     )
        V_ssa_Ac(1   ,grid%nx  )= V_ssa_Ac(     2,grid%nx) + V_ssa_Ac(     1,grid%nx-1) - V_ssa_Ac(     2,grid%nx-1)
        V_ssa_Ac(grid%ny-1,1   )= V_ssa_Ac(grid%ny-2,1   ) + V_ssa_Ac(grid%ny-1,2     ) - V_ssa_Ac(grid%ny-2,2     )
        V_ssa_Ac(grid%ny-1,grid%nx)= V_ssa_Ac(grid%ny-2,grid%nx) + V_ssa_Ac(grid%ny-1,grid%nx-1) - V_ssa_Ac(grid%ny-2,grid%nx-1)

        ! Check convergence:
        IF(largest_rest < max_residual_U) THEN
         velocity_iterations = velocity_iterations + 1
         EXIT iterate
        ELSE IF(iter_gauss_seidel == maximum_number_gauss_seidel_iterations) THEN
         !IF(C%choice_write_message) WRITE(C%stderr,FMT='(A, I4, A, F13.1)') 'From ice_stream_velocity_by_ssa_stag(): No convergence after iteration: ', maximum_number_gauss_seidel_iterations, ', at time: ', G%time
         WRITE(0,*) ' ERROR - SOR iteration wont converge!'
         STOP
        END IF
        velocity_iterations = velocity_iterations + 1

      END DO iterate

    END DO ! End: explicit_update_iteration

  END SUBROUTINE ice_stream_velocity_by_ssa
  SUBROUTINE calculate_basal_yield_stress_Ab( grid, Hi_Ab, Hb_Ab, sealevel_Ab, phi_fric_Ab, basal_yield_stress_Ab)
    ! Calculate the basal yield stress, which is used to determine the basal stress, used
    ! for sliding. Using the parameterisations given by Martin et al. (TC: 2011, Equ 9)
    ! And as described in de Boer et al. (Climate Dynamics, 2013).
    ! Heiko Goelzer (h.goelzer@uu.nl) 2017
    
    USE parameters_module,          ONLY: n_flow, ice_density, seawater_density, grav

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1), INTENT(IN)  :: Hi_Ab
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1), INTENT(IN)  :: Hb_Ab                  ! the bedrock height [m]
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1), INTENT(IN)  :: sealevel_Ab            ! sea level relative to PD [m]

    ! Output variables:
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1), INTENT(OUT) :: phi_fric_Ab            ! the basal yield stress [Pa]
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1), INTENT(OUT) :: basal_yield_stress_Ab  ! the basal yield stress [Pa]

    ! Local variables:
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1)              :: lambda_p               ! scaling of the pore water pressure
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1)              :: pore_water_pressure    ! pore water pressure [Pa]
    
    REAL(dp), PARAMETER                                 :: scaling_Hb_minimal_friction_angle   = -1000._dp
    REAL(dp), PARAMETER                                 :: scaling_Hb_maximal_friction_angle   = 0._dp
    REAL(dp), PARAMETER                                 :: minimal_friction_angle = 5._dp
    REAL(dp), PARAMETER                                 :: maximal_friction_angle = 20._dp 

    ! The pore water pressure is scaled with a bedrock height dependend parameterisation
    ! Equation (12) in Martin et al. (2011)
    WHERE ((Hb_Ab - sealevel_Ab) <= 0._dp)
     lambda_p = 1._dp
    ELSE WHERE ((Hb_Ab - sealevel_Ab) >= 1000._dp)
     lambda_p = 0._dp
    ELSE WHERE ! between 0 and 1000
     lambda_p = 1._dp - (Hb_Ab - sealevel_Ab) / 1000._dp
    END WHERE

    ! The pore water pressure, equation (11) in Martin et al. (2011)
    pore_water_pressure = 0.96_dp * ice_density * grav * Hi_Ab * lambda_p 

    ! The friction angle, used for the yield stress, equation (10) in Martin et al. (2011)
    WHERE (Hb_Ab <= scaling_Hb_minimal_friction_angle)
      phi_fric_Ab = minimal_friction_angle
    ELSE WHERE (Hb_Ab >= scaling_Hb_maximal_friction_angle)
      phi_fric_Ab = maximal_friction_angle
    ELSE WHERE ! between C%scaling_Hb_maximal_friction_angle and C%scaling_Hb_minimal_friction_angle
      phi_fric_Ab = minimal_friction_angle + (maximal_friction_angle - minimal_friction_angle) * (1._dp &
                  + (Hb_Ab - scaling_Hb_maximal_friction_angle) / (scaling_Hb_maximal_friction_angle - scaling_Hb_minimal_friction_angle))
    END WHERE

    ! calculate yield stress everywhere, restrictions to sheet applied elsewhere, equation (9) in Martin et al. (2011)
    basal_yield_stress_Ab = TAN(C%deg2rad * phi_fric_Ab) * (ice_density * grav * Hi_Ab - pore_water_pressure)
    
  END SUBROUTINE calculate_basal_yield_stress_Ab
  SUBROUTINE smooth_ice_thickness_shelf( grid, mask_shelf_Aa, Hi_Aa, rx, ry)
    ! Smooth the ice shelf thickness for numerical stability with the SOR solver (cf. Huybrechts P.126)
    ! Heiko Goelzer (h.goelzer@uu.nl) 2016
    
    IMPLICIT NONE
 
    ! In/Output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    INTEGER,  DIMENSION(grid%ny,grid%nx), INTENT(IN)     :: mask_shelf_Aa
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(INOUT)  :: Hi_Aa
    INTEGER                       , INTENT(IN)     :: rx, ry                    ! The smooth radius

    ! Local variables:
    INTEGER                                        :: i, j, ii, jj
    REAL(dp)                                       :: hnc, hnd, dist

    ! smoothing the shelf (cf. Huybrechts P.126)
    DO J = 2,grid%ny-1
      DO I = 2,grid%nx-1
        IF(mask_shelf_Aa(j,i) == 1) THEN
         hnc=0
         hnd=0
         DO jj = j-ry,j+ry
           DO ii = i-rx,i+rx
             IF(mask_shelf_Aa(jj,ii) == 1) THEN
              IF(ii.EQ.i.AND.jj.EQ.j) THEN
               hnc=hnc+Hi_Aa(jj,ii)
               hnd=hnd+1.
              ELSE
               dist=(i-ii)**2+(j-jj)**2
               hnc=hnc+0.05*Hi_Aa(jj,ii)/dist
               hnd=hnd+0.05/dist
              ENDIF
             ENDIF
           END DO
         END DO
         Hi_Aa(j,i)=hnc/hnd
        ENDIF
      END DO
    END DO
    
  END SUBROUTINE smooth_ice_thickness_shelf
  SUBROUTINE schoof_grounding_line_flux( grid, SL_Aa, Hi_Aa, Hb_Aa, A_flow_mean_Aa, mask_GL_Aa, mask_shelf_Aa, buttressing_term_x, buttressing_term_y, flux_xr, flux_xl, flux_yu, flux_yd)
    ! Adjusts the grounding line flux according to the semi-analytical solution of Schoof (2007, Eq. (16))
    ! This includes a theta term for the buttressing effect following Pollard & DeConto (2012, Eq. (8)) and Pattyn (2017)
    
    USE parameters_module, ONLY: n_flow, grav, ice_density, seawater_density
    
    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: SL_Aa           ! In meter relative to present day eustatic sea level
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: Hi_Aa
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: Hb_Aa
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: A_flow_mean_Aa
    INTEGER,  DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: mask_GL_Aa
    INTEGER,  DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: mask_shelf_Aa
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: buttressing_term_x ! Buttressing term in Schoof equation
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: buttressing_term_y ! Buttressing term in Schoof equation

    ! In/Output variables:
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(INOUT) :: flux_xr            ! The flux in the x-direction to the right coinciding with the Cx staggered grid 
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(INOUT) :: flux_xl            ! The flux in the x-direction to the left  coinciding with the Cx staggered grid 
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(INOUT) :: flux_yu            ! The flux in the y-direction upwards      coinciding with the Cy staggered grid 
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(INOUT) :: flux_yd            ! The flux in the y-direction downwards    coinciding with the Cy staggered grid 

    ! Local variables:
    INTEGER                                            :: i, j
    REAL(dp)                                           :: exponent           ! The exponent of the Hi_gl_subgrid term in the Schoof grounding line parameterisation
    REAL(dp), DIMENSION(     grid%ny,grid%nx)          :: q_gl               ! The Schoof grounding line flux at the subgrid grounding line position

    REAL(dp)                                           :: Hi_gl_subgrid
    REAL(dp)                                           :: factor             ! The factor in front of the Hi_gl_subgrid term in the Schoof grounding line parameterisation
    REAL(dp)                                           :: C_sliding
    
    REAL(dp), PARAMETER :: C_m_sliding = 1._dp / 3._dp
    REAL(dp), PARAMETER :: C_C_sliding = 7.624E+06_dp
    INTEGER,  PARAMETER :: choice_heuristic_flux_rule_at_grounding_line = 1

    exponent  = (C_m_sliding + n_flow + 3._dp) / (C_m_sliding + 1._dp)
    C_sliding = C_C_sliding

    DO i = 2, grid%nx-1
    DO j = 2, grid%ny-1

      ! At grounding line points for which the bottom of the sheet ends on a point above sea level (Hb > sea level), the schoof flux is not applied.
      IF (mask_GL_Aa( j,i) == 1 .AND. (Hb_Aa( j,i) < SL_Aa( j,i))) THEN  

       ! The factor in front of the Hi_gl_subgrid term in the Schoof grounding line parameterisation:
       factor = (( A_flow_mean_Aa(j,i) * (ice_density * grav)**(n_flow + 1._dp) &   ! It is possible here to interpolate A_flow_mean as well, but maybe not so important
                    * (1._dp - ice_density / seawater_density)**n_flow &
                 ) / (4._dp**n_flow * C_sliding) &
                )**(1._dp / (C_m_sliding + 1._dp))

       !IF(.NOT. C%one_dimensional_experiment == 'in-y-direction') THEN
        ! Flux in x-direction right
        IF (mask_shelf_Aa( j,i+1) == 1) THEN
         ! Output: Hi_gl_subgrid
         CALL ice_thickness_at_subgrid_grounding_line( SL_Aa( j,i), Hb_Aa( j,i), Hb_Aa( j,i+1), Hi_Aa( j,i), Hi_Aa( j,i+1), Hi_gl_subgrid)
         ! We combine here factor on the last grounded point with buttressing_term_x defined on the first shelf point, which seems OK given that Hi_gl_subgrid is somehwere in between
         q_gl( j,i) = factor * buttressing_term_x(j,i+1) * Hi_gl_subgrid**exponent
         ! Apply the heuristic rule for the flux boundary condition according to Pollard & DeConto (2012, Eq. (9)):
         IF(q_gl(j,i) > flux_xr(j,i)) THEN
          flux_xr(j,i) = q_gl(j,i)
         ELSE
          IF(choice_heuristic_flux_rule_at_grounding_line == 1) flux_xr(j,i+1) = q_gl(j,i)            ! Pollard & DeCconto (2012, GMD) heuristic rule
         END IF ! IF(q_gl(j,i) > flux_xr(j,i)) THEN
        END IF ! IF(mask(j,i+1) == C%type_shelf) THEN

        ! Flux in x-direction left
        IF (mask_shelf_Aa( j,i-1) == 1) THEN
         ! Output: Hi_gl_subgrid
         CALL ice_thickness_at_subgrid_grounding_line( SL_Aa(j,i), Hb_Aa(j,i), Hb_Aa(j,i-1), Hi_Aa(j,i), Hi_Aa(j,i-1), Hi_gl_subgrid)
         q_gl(j,i) = factor * buttressing_term_x(j,i-1) * Hi_gl_subgrid**exponent
         ! Apply the heuristic rule for the flux boundary condition according to Pollard & DeConto (2012, Eq. (9)):
         IF(-q_gl(j,i) < flux_xl(j,i)) THEN
          flux_xl(j,i) = -q_gl(j,i)
         ELSE
          IF(choice_heuristic_flux_rule_at_grounding_line == 1) flux_xl(j,i-1) = -q_gl(j,i)           ! Pollard & DeCconto (2012, GMD) heuristic rule
         END IF ! IF(-q_gl(j,i) < flux_xl(j,i)) THEN
        END IF ! IF(mask(j,i-1) == C%type_shelf) THEN
       !END IF

       !IF(.NOT. C%one_dimensional_experiment == 'in-x-direction') THEN
        ! Flux in y-direction up
        IF (mask_shelf_Aa( j+1,i) == 1) THEN
         ! Output: Hi_gl_subgrid
         CALL ice_thickness_at_subgrid_grounding_line( SL_Aa(j,i), Hb_Aa(j,i), Hb_Aa(j+1,i), Hi_Aa(j,i), Hi_Aa(j+1,i), Hi_gl_subgrid)
         q_gl(j,i) = factor * buttressing_term_y(j+1,i) * Hi_gl_subgrid**exponent
         ! Apply the heuristic rule for the flux boundary condition according to Pollard & DeConto (2012, Eq. (9)):
         IF(q_gl(j,i) > flux_yu(j,i)) THEN
          flux_yu(j,i) = q_gl(j,i)
         ELSE
          IF(choice_heuristic_flux_rule_at_grounding_line == 1) flux_yu(j+1,i) = q_gl(j,i)            ! Pollard & DeCconto (2012, GMD) heuristic rule
         END IF ! IF(q_gl(j,i) > flux_yu(j,i)) THEN
        END IF ! IF(mask(j+1,i) == C%type_shelf) THEN

        ! Flux in y-direction down
        IF (mask_shelf_Aa( j-1,i) == 1) THEN
         ! Output: Hi_gl_subgrid
         CALL ice_thickness_at_subgrid_grounding_line( SL_Aa(j,i), Hb_Aa(j,i), Hb_Aa(j-1,i), Hi_Aa(j,i), Hi_Aa(j-1,i), Hi_gl_subgrid)
         q_gl(j,i) = factor * buttressing_term_y(j-1,i) * Hi_gl_subgrid**exponent
         ! Apply the heuristic rule for the flux boundary condition according to Pollard & DeConto (2012, Eq. (9)):
         IF(-q_gl(j,i) < flux_yd(j,i)) THEN
          flux_yd(j,i) = -q_gl(j,i)
         ELSE
          IF(choice_heuristic_flux_rule_at_grounding_line == 1) flux_yd(j-1,i) = -q_gl(j,i)           ! Pollard & DeCconto (2012, GMD) heuristic rule
         END IF ! IF(-q_gl(j,i) < flux_yd(j,i)) THEN
        END IF ! IF(mask(j-1,i) == C%type_shelf) THEN
       !END IF ! IF(.NOT. C%one_dimensional_experiment == 'in-x-direction') THEN
     END IF ! IF(mask(j,i) == C%type_groundline .AND. (Hb(j,i) < sealevel(j,i))) THEN  

    END DO
    END DO
    
  END SUBROUTINE schoof_grounding_line_flux
  SUBROUTINE ice_thickness_at_subgrid_grounding_line( SL, Hb_gl, Hb_shelf, Hi_gl, Hi_shelf, Hi_gl_subgrid)
    ! This routine calculates the ice thickness at the subgrid grounding line. This is done in a similar
    ! way as in Gladstone et al. (2010, eqs. (1) -- (9)) and Pollard & DeConto (2012). First the subgrid grounding
    ! line position has to be found be interpolating the flotation condition.
    ! In Gladstone et al. (2010) the b definition is inconsequent, but we follow the same idea here, using:
    !  The floatation condition: rho_i * Hi = rho_w * (sealevel - Hb)
    !  Using a dimensionless scaled coordinate in the interpolation direction, e.g. in the x-direction this is: lambda = (x - x_gl) / dx
    !  and the linear interpolation function for Hi: Hi(lambda) = Hi_gl(1-lambda) + Hi_shelf * lambda
    !  and the linear interpolation function for Hb: Hb(lambda) = Hb_gl(1-lambda) + Hb_shelf * lambda
    ! The lambda_gl follows from: 
    !  rho_i * Hi(lambda_gl_subgrid) = rho_w * (sealevel - Hb(lambda_gl_subgrid))
    ! with
    !  Hb(lambda_gl_subgrid) = Hb_gl(1-lambda_gl_subgrid) + Hb_shelf * lambda_gl_subgrid
    !  Hi(lambda_gl_subgrid) = Hi_gl(1-lambda_gl_subgrid) + Hi_shelf * lambda_gl_subgrid
    ! Finally Hi(lambda_gl_subgrid) can be calculated by evaluating the this expression by substituting the determined lambda_gl_subgrid.
    
    USE parameters_module, ONLY: ice_density, seawater_density
    
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: SL                ! The sealevel [m] which is assumed to be locally of uniform level
    REAL(dp), INTENT(IN)  :: Hb_gl             ! The Hb [m] on the model grid for the grounding line  point
    REAL(dp), INTENT(IN)  :: Hb_shelf          ! The Hb [m] on the model grid for the neighbour shelf point
    REAL(dp), INTENT(IN)  :: Hi_gl             ! The Hi [m] on the model grid for the grounding line  point
    REAL(dp), INTENT(IN)  :: Hi_shelf          ! The Hi [m] on the model grid for the neighbour shelf point

    ! Output variables:
    REAL(dp), INTENT(OUT) :: Hi_gl_subgrid     ! The Hi on the subgrid for the grounding line point at the point were the interpolated Hi - Hb geometry equals the floatation condition

    ! Local variables:
    REAL(dp)              :: lambda_gl_subgrid ! The dimensionless scaled coordinate in the interpolation direction at the point of the subgrid grounding line

    lambda_gl_subgrid = (seawater_density * (Hb_gl - SL      ) + ice_density * Hi_gl) &
                      / (seawater_density * (Hb_gl - Hb_shelf) + ice_density * (Hi_gl - Hi_shelf))

    Hi_gl_subgrid = Hi_gl * (1._dp - lambda_gl_subgrid) + Hi_shelf * lambda_gl_subgrid
    
  END SUBROUTINE ice_thickness_at_subgrid_grounding_line
  SUBROUTINE shelf_buttressing( grid, Hi_Aa, A_flow_mean_Aa, mask_shelf_Aa, Us_Acx, Vs_Acy, buttressing_term_x, buttressing_term_y)
    ! Calculate the buttressing term which is used in the Schoof parameterisation.
    
    USE parameters_module, ONLY: n_flow, grav, ice_density, seawater_density
    
    IMPLICIT NONE

    ! Input variables: 
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(IN)  :: Hi_Aa                 ! Aa
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(IN)  :: A_flow_mean_aa        ! Aa
    INTEGER,  DIMENSION(grid%ny,grid%nx), INTENT(IN)  :: mask_shelf_Aa               ! Aa
    REAL(dp), DIMENSION(grid%ny,grid%nx-1), INTENT(IN)  :: Us_Acx
    REAL(dp), DIMENSION(grid%ny-1,grid%nx), INTENT(IN)  :: Vs_Acy

    ! Output variables:
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(OUT) :: buttressing_term_x ! Aa; shelf points only
    REAL(dp), DIMENSION(grid%ny,grid%nx), INTENT(OUT) :: buttressing_term_y ! Aa; shelf points only

    ! Local variables:
    INTEGER                                     :: i, j
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: dUs_dx             ! Aa
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: dVs_dy             ! Aa
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: dUs_dy             ! Aa
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: dVs_dx             ! Aa
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: tau_stress_free    ! Free stress without buttressing Aa
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: tau_xx             ! Stress with buttressing Aa
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: tau_yy             ! Stress with buttressing Aa
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: theta_x            ! Aa; shelf points only
    REAL(dp), DIMENSION(grid%ny,grid%nx)              :: theta_y            ! Aa; shelf points only
    
    REAL(dp), PARAMETER :: epsilon_sq_0 = 1E-12_dp
    REAL(dp), PARAMETER :: C_m_sliding = 1._dp / 3._dp
    INTEGER,  PARAMETER :: choice_regulation_of_buttressing_range = 1
    INTEGER,  PARAMETER :: choice_buttressing_method = 1
    REAL(dp), PARAMETER :: buttressing_enhancement_factor = 1._dp

!    DO i = 2, grid%nx-1
!    DO j = 2, grid%ny-1
!      ! Gradients from Ac to Aa
!      ! Just wondering if the velocities wHi_Aach are used to determine the gradients must be shelf velocities and not GL velocities because it concerns the shelf stress.
!      dUs_dx(j,i) = (Us_Ac(j,i) - Us_Ac(j  ,i-1)) / C%dx 
!      dVs_dy(j,i) = (Vs_Ac(j,i) - Vs_Ac(j-1,i  )) / C%dy 
!      ! ??: Gradients from Ac to Ab averaged to Aa    or    gradients 2Ac to Ac averaged to Aa
!     !dUs_dy = 0.25_dp / (C%dy) * ( Us_Ac(j+1,i) - Us_Ac(j,i) + Us_Ac(j,i) - Us_Ac(j-1,i) + Us_Ac(j+1,i-1) - Us_Ac(j,i-1) + Us_Ac(j,i-1) - Us_Ac(j-1,i-1) )
!      dUs_dy(j,i) = 0.25_dp / C%dy * ( Us_Ac(j+1,i) - Us_Ac(j-1,i) + Us_Ac(j+1,i-1) - Us_Ac(j-1,i-1) )
!      dVs_dx(j,i) = 0.25_dp / C%dx * ( Vs_Ac(j,i+1) - Vs_Ac(j,i-1) + Vs_Ac(j-1,i+1) - Vs_Ac(j-1,i-1) ) 
!    END DO
!    END DO
    
    CALL ddx_Acx_to_Aa_2D( grid, Us_Acx, dUs_dx)
    CALL ddy_Acx_to_Aa_2D( grid, Us_Acx, dUs_dy)
    CALL ddx_Acy_to_Aa_2D( grid, Vs_Acy, dVs_dx)
    CALL ddy_Acy_to_Aa_2D( grid, Vs_Acy, dVs_dy)

    DO i = 2, grid%nx-1
    DO j = 2, grid%ny-1
     IF (mask_shelf_Aa( j,i) == 1) THEN
      ! The free stress in absence of any buttressing (Eq 8, Pollard & DeConto, 2012): tau_f = 0.5 rho_i * g * h_g * (1- rho_i / rho_w). 
      SELECT CASE(choice_regulation_of_buttressing_range)
      CASE(1)
       tau_stress_free(j,i) = 0.5_dp * ice_density * grav *         Hi_Aa(j,i)                    * (1._dp - ice_density / seawater_density)
      CASE(2)
       tau_stress_free(j,i) = 0.5_dp * ice_density * grav * max(min(Hi_Aa(j,i),1000._dp),100._dp) * (1._dp - ice_density / seawater_density)   ! Regularisation on Hi_Aa (Pattyn, 2017, Eq. (9))
      END SELECT

      SELECT CASE(choice_buttressing_method)
      CASE(1,3)
       ! The longitudinal deviatoric stress (See Pollard & DeConto (2012, GMD) Eq 8):   ! Note: Pollard & DeConto, 2012 take (0.01_dp / C%dx)^2 instead of C%epsilon_sq_0
       tau_xx(j,i) = abs(A_flow_mean_Aa(j,i)**(-1._dp / n_flow) * ( dUs_dx(j,i)**2 + dVs_dy(j,i)**2 + dUs_dx(j,i) * dVs_dy(j,i) + 0.25_dp * (dUs_dy(j,i) + dVs_dx(j,i))**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow)) * dUs_dx(j,i))
       tau_yy(j,i) = abs(A_flow_mean_Aa(j,i)**(-1._dp / n_flow) * ( dUs_dx(j,i)**2 + dVs_dy(j,i)**2 + dUs_dx(j,i) * dVs_dy(j,i) + 0.25_dp * (dUs_dy(j,i) + dVs_dx(j,i))**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow)) * dVs_dy(j,i))
      CASE(2,4)
       ! The longitudinal deviatoric stress (See Pollard & DeConto (2012, GMD) Eq 8):
       ! A fixed A_flow_mean_Aa value is prescibed, since artefacts across the grounding line due to artefacts in Ti wHi_Aach in turn are due to artefacts in W (namely in W(zeta=1) at the bottom)
       tau_xx(j,i) = abs(       3.e-16_dp**(-1._dp / n_flow) * ( dUs_dx(j,i)**2 + dVs_dy(j,i)**2 + dUs_dx(j,i) * dVs_dy(j,i) + 0.25_dp * (dUs_dy(j,i) + dVs_dx(j,i))**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow)) * dUs_dx(j,i))
       tau_yy(j,i) = abs(       3.e-16_dp**(-1._dp / n_flow) * ( dUs_dx(j,i)**2 + dVs_dy(j,i)**2 + dUs_dx(j,i) * dVs_dy(j,i) + 0.25_dp * (dUs_dy(j,i) + dVs_dx(j,i))**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow)) * dVs_dy(j,i))
      CASE(5)
       ! The longitudinal stress for the given shelf geometry; (See Pattyn, 2017;  Cuffey and Paterson, Eq 8.129)
       ! A fixed A_flow_mean_Aa value is prescibed, since artefacts across the grounding line due to artefacts in Ti wHi_Aach in turn are due to artefacts in W (namely in W(zeta=1) at the bottom)
       tau_xx(j,i) = abs(A_flow_mean_Aa(j,i)**(-1._dp / n_flow) * ( dUs_dx(j,i)**2 + dVs_dy(j,i)**2 + dUs_dx(j,i) * dVs_dy(j,i) + 0.25_dp * (dUs_dy(j,i) + dVs_dx(j,i))**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow)) * (2._dp * dUs_dx(j,i) +         dVs_dy(j,i)))
       tau_yy(j,i) = abs(A_flow_mean_Aa(j,i)**(-1._dp / n_flow) * ( dUs_dx(j,i)**2 + dVs_dy(j,i)**2 + dUs_dx(j,i) * dVs_dy(j,i) + 0.25_dp * (dUs_dy(j,i) + dVs_dx(j,i))**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow)) * (        dUs_dx(j,i) + 2._dp * dVs_dy(j,i)))
      CASE(6)
       ! The longitudinal stress for the given shelf geometry; (See Pattyn, 2017;  Cuffey and Paterson, Eq 8.129)
       tau_xx(j,i) = abs(       3.e-16_dp**(-1._dp / n_flow) * ( dUs_dx(j,i)**2 + dVs_dy(j,i)**2 + dUs_dx(j,i) * dVs_dy(j,i) + 0.25_dp * (dUs_dy(j,i) + dVs_dx(j,i))**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow)) * (2._dp * dUs_dx(j,i) +         dVs_dy(j,i)))
       tau_yy(j,i) = abs(       3.e-16_dp**(-1._dp / n_flow) * ( dUs_dx(j,i)**2 + dVs_dy(j,i)**2 + dUs_dx(j,i) * dVs_dy(j,i) + 0.25_dp * (dUs_dy(j,i) + dVs_dx(j,i))**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow)) * (        dUs_dx(j,i) + 2._dp * dVs_dy(j,i)))
      END SELECT

      SELECT CASE(choice_buttressing_method)
      CASE(3,4)
       ! The longitudinal full stress, with help of the longitudinal deviatoric stress:
       tau_xx(j,i) = tau_xx(j,i) - 0.5_dp * ice_density * grav * Hi_Aa(j,i)  ! From deviatoric stress to full stress
       tau_yy(j,i) = tau_yy(j,i) - 0.5_dp * ice_density * grav * Hi_Aa(j,i)  ! From deviatoric stress to full stress
      END SELECT

      ! The ratio gives the buttressing term theta for both directions:
      theta_x(j,i) = buttressing_enhancement_factor * tau_xx(j,i) / tau_stress_free(j,i)
      theta_y(j,i) = buttressing_enhancement_factor * tau_yy(j,i) / tau_stress_free(j,i)
     ELSE
      theta_x(j,i) = 0.0_dp
      theta_y(j,i) = 0.0_dp
     END IF
    END DO
    END DO

    ! Set the grid margins to zero:
    theta_x(:,   1) = 0.0_dp
    theta_x(:,grid%nx) = 0.0_dp
    theta_x(1,   :) = 0.0_dp
    theta_x(grid%ny,:) = 0.0_dp

    theta_y(:,   1) = 0.0_dp
    theta_y(:,grid%nx) = 0.0_dp
    theta_y(1,   :) = 0.0_dp
    theta_y(grid%ny,:) = 0.0_dp

    ! Limit the theta values witHi_Aan a range from 0 to 1:
    theta_x = max(min(theta_x, 1.0_dp), 0.0_dp)
    theta_y = max(min(theta_y, 1.0_dp), 0.0_dp)

    buttressing_term_x = theta_x**(n_flow / (C_m_sliding + 1._dp))
    buttressing_term_y = theta_y**(n_flow / (C_m_sliding + 1._dp))
    
  END SUBROUTINE shelf_buttressing
  SUBROUTINE schoof_grounding_line_flux_down_stream_buttressing_stress( grid, SL_Aa, Hi_Aa, Hb_Aa, A_flow_mean_Aa, mask_GL_Aa, mask_shelf_Aa, Us, Vs, flux_xr, flux_xl, flux_yu, flux_yd, buttressing_term_x, buttressing_term_y)
    ! Adjusts the grounding line flux according to the semi-analytical solution of Schoof (2007, Eq. (16))
    ! This includes a theta term for the buttressing effect following Pollard & DeConto (2012, Eq. (8)) and Pattyn (2017)
    ! Calculating the velocity gradients which are used to compute the buttressing stress on the shelf near the grounding line points
    ! in a down-stream manner (the staggering might be not conform the default staggering approach, so a check on that is needed):
    
    USE parameters_module, ONLY: n_flow, grav, ice_density, seawater_density
    
    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: SL_Aa           ! In meter relative to present day eustatic sea level
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: Hi_Aa
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: Hb_Aa
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: A_flow_mean_Aa
    INTEGER,  DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: mask_GL_Aa
    INTEGER,  DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: mask_shelf_Aa
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: Us                 ! Us_Ac?
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(IN)    :: Vs                 ! Vs_Ac?

    ! In/Output variables:
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(INOUT) :: flux_xr            ! The flux in the x-direction to the right coinciding with the Cx staggered grid 
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(INOUT) :: flux_xl            ! The flux in the x-direction to the left  coinciding with the Cx staggered grid 
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(INOUT) :: flux_yu            ! The flux in the y-direction upwards      coinciding with the Cy staggered grid 
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(INOUT) :: flux_yd            ! The flux in the y-direction downwards    coinciding with the Cy staggered grid 
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(OUT)   :: buttressing_term_x ! Buttressing term in Schoof equation, output for monitoring the field
    REAL(dp), DIMENSION(     grid%ny,grid%nx), INTENT(OUT)   :: buttressing_term_y ! Buttressing term in Schoof equation, output for monitoring the field

    ! Local variables:
    INTEGER                                            :: i, j
    REAL(dp)                                           :: exponent           ! The exponent of the Hi_gl_subgrid term in the Schoof grounding line parameterisation
    REAL(dp), DIMENSION(     grid%ny,grid%nx)                :: q_gl               ! The Schoof grounding line flux at the subgrid grounding line position

    REAL(dp)                                           :: Hi_gl_subgrid
    REAL(dp)                                           :: factor             ! The factor in front of the Hi_gl_subgrid term in the Schoof grounding line parameterisation
    REAL(dp)                                           :: C_sliding
    REAL(dp)                                           :: buttress           ! Buttressing term in Schoof equation
    
    REAL(dp), PARAMETER :: C_m_sliding = 1._dp / 3._dp
    REAL(dp), PARAMETER :: C_C_sliding = 7.624E+06_dp
    INTEGER,  PARAMETER :: choice_heuristic_flux_rule_at_grounding_line = 1

    exponent  = (C_m_sliding + n_flow + 3._dp) / (C_m_sliding + 1._dp)
    C_sliding = C_C_sliding

    buttressing_term_x = 0._dp
    buttressing_term_y = 0._dp

    DO i = 2, grid%nx-1
    DO j = 2, grid%ny-1

      ! At grounding line points for which the bottom of the sheet ends on a point above sea level (Hb_Aa > sea level), the schoof flux is not applied.
      IF (mask_GL_Aa( j,i) == 1 .AND. (Hb_Aa(j,i) < SL_Aa(j,i))) THEN 

       ! The factor in front of the Hi_gl_subgrid term in the Schoof grounding line parameterisation:
       factor = (( A_flow_mean_Aa(j,i) * (ice_density * grav)**(n_flow + 1._dp) &   ! It is possible here to interpolate A_flow_mean_Aa as well, but maybe not so important
                    * (1._dp - ice_density / seawater_density)**n_flow &
                 ) / (4._dp**n_flow * C_sliding) &
                )**(1._dp / (C_m_sliding + 1._dp))

        ! Flux in x-direction right
        IF (mask_shelf_Aa( j,i+1) == 1) THEN
         ! Output: Hi_gl_subgrid
         CALL ice_thickness_at_subgrid_grounding_line(SL_Aa(j,i), Hb_Aa(j,i), Hb_Aa(j,i+1), Hi_Aa(j,i), Hi_Aa(j,i+1), Hi_gl_subgrid)
         ! Output: buttress
         CALL buttressing(Hi_Aa(j,i+1), A_flow_mean_Aa(j,i+1), dUs_dx = (Us(j  ,i+2) - Us(j  ,i+1)) / grid%dx, &
                                                         dUs_dy = (Us(j  ,i+1) - Us(j-1,i+1)) / grid%dy, &      !(0.25_dp / C%dy) * ( )
                                                         dVs_dx = (Vs(j  ,i+2) - Vs(j  ,i+1)) / grid%dx, &      !(0.25_dp / C%dx) * ( )
                                                         dVs_dy = (Vs(j  ,i+1) - Vs(j-1,i+1)) / grid%dy, &
                                                         buttressing_term_x = buttress)
         buttressing_term_x(j,i+1) = buttress
         ! We combine here factor on the last grounded point with buttressing_term_x defined on the first shelf point, which seems OK given that Hi_gl_subgrid is somehwere in between
         q_gl(j,i) = factor * buttress * Hi_gl_subgrid**exponent
         ! Apply the heuristic rule for the flux boundary condition according to Pollard & DeConto (2012, Eq. (9)):
         IF(q_gl(j,i) > flux_xr(j,i)) THEN
          flux_xr(j,i) = q_gl(j,i)
         ELSE
          IF(choice_heuristic_flux_rule_at_grounding_line == 1) flux_xr(j,i+1) = q_gl(j,i)            ! Pollard & DeCconto (2012, GMD) heuristic rule
          !IF(choice_heuristic_flux_rule_at_grounding_line == 1) flux_xl(j,i+1) = q_gl(j,i)            ! Pollard & DeCconto (2012, GMD) heuristic rule
         END IF
        END IF

        ! Flux in x-direction left
        IF (mask_shelf_Aa( j,i-1) == 1) THEN
         ! Output: Hi_gl_subgrid
         CALL ice_thickness_at_subgrid_grounding_line(SL_Aa(j,i), Hb_Aa(j,i), Hb_Aa(j,i-1), Hi_Aa(j,i), Hi_Aa(j,i-1), Hi_gl_subgrid)
         ! Output: buttress
         CALL buttressing(Hi_Aa(j,i-1), A_flow_mean_Aa(j,i-1), dUs_dx = (Us(j  ,i-1) - Us(j  ,i-2)) / grid%dx, &
                                                         dUs_dy = (Us(j  ,i-1) - Us(j-1,i-1)) / grid%dy, &      !(0.25_dp / C%dy) * ( )
                                                         dVs_dx = (Vs(j  ,i-1) - Vs(j  ,i-2)) / grid%dx, &      !(0.25_dp / C%dx) * ( )
                                                         dVs_dy = (Vs(j  ,i-1) - Vs(j-1,i-1)) / grid%dy, &
                                                         buttressing_term_x = buttress)
         buttressing_term_x(j,i-1) = buttress
         q_gl(j,i) = factor * buttress * Hi_gl_subgrid**exponent
         ! Apply the heuristic rule for the flux boundary condition according to Pollard & DeConto (2012, Eq. (9)):
         IF(-q_gl(j,i) < flux_xl(j,i)) THEN
          flux_xl(j,i) = -q_gl(j,i)
         ELSE
          IF(choice_heuristic_flux_rule_at_grounding_line == 1) flux_xl(j,i-1) = -q_gl(j,i)           ! Pollard & DeCconto (2012, GMD) heuristic rule
          !IF(choice_heuristic_flux_rule_at_grounding_line == 1) flux_xr(j,i-1) = -q_gl(j,i)           ! Pollard & DeCconto (2012, GMD) heuristic rule
         END IF
        END IF

        ! Flux in y-direction up
        IF (mask_shelf_Aa( j+1,i) == 1) THEN
         ! Output: Hi_gl_subgrid
         CALL ice_thickness_at_subgrid_grounding_line(SL_Aa(j,i), Hb_Aa(j,i), Hb_Aa(j+1,i), Hi_Aa(j,i), Hi_Aa(j+1,i), Hi_gl_subgrid)
         ! Output: buttress
         CALL buttressing(Hi_Aa(j+1,i), A_flow_mean_Aa(j+1,i), dUs_dx = (Us(j+1,i  ) - Us(j+1,i-1)) / grid%dx, &
                                                         dUs_dy = (Us(j+2,i  ) - Us(j+1,i  )) / grid%dy, &      !(0.25_dp / C%dy) * ( )
                                                         dVs_dx = (Vs(j+1,i  ) - Vs(j+1,i-1)) / grid%dx, &      !(0.25_dp / C%dx) * ( )
                                                         dVs_dy = (Vs(j+2,i  ) - Vs(j+1,i  )) / grid%dy, &
                                                         buttressing_term_y = buttress)
         buttressing_term_y(j+1,i) = buttress
         q_gl(j,i) = factor * buttress * Hi_gl_subgrid**exponent
         ! Apply the heuristic rule for the flux boundary condition according to Pollard & DeConto (2012, Eq. (9)):
         IF(q_gl(j,i) > flux_yu(j,i)) THEN
          flux_yu(j,i) = q_gl(j,i)
         ELSE
          IF(choice_heuristic_flux_rule_at_grounding_line == 1) flux_yu(j+1,i) = q_gl(j,i)            ! Pollard & DeCconto (2012, GMD) heuristic rule
          !IF(choice_heuristic_flux_rule_at_grounding_line == 1) flux_yd(j+1,i) = q_gl(j,i)            ! Pollard & DeCconto (2012, GMD) heuristic rule
         END IF
        END IF

        ! Flux in y-direction down
        IF (mask_shelf_Aa( j-1,i) == 1) THEN
         ! Output: Hi_gl_subgrid
         CALL ice_thickness_at_subgrid_grounding_line(SL_Aa(j,i), Hb_Aa(j,i), Hb_Aa(j-1,i), Hi_Aa(j,i), Hi_Aa(j-1,i), Hi_gl_subgrid)
         ! Output: buttress
         CALL buttressing(Hi_Aa(j-1,i), A_flow_mean_Aa(j-1,i), dUs_dx = (Us(j-1,i  ) - Us(j-1,i-1)) / grid%dx, &
                                                         dUs_dy = (Us(j-1,i  ) - Us(j-2,i  )) / grid%dy, &      !(0.25_dp / C%dy) * ( )
                                                         dVs_dx = (Vs(j-1,i  ) - Vs(j-1,i-1)) / grid%dx, &      !(0.25_dp / C%dx) * ( )
                                                         dVs_dy = (Vs(j-1,i  ) - Vs(j-2,i  )) / grid%dy, &
                                                         buttressing_term_y = buttress)
         buttressing_term_y(j-1,i) = buttress
         q_gl(j,i) = factor * buttress * Hi_gl_subgrid**exponent
         ! Apply the heuristic rule for the flux boundary condition according to Pollard & DeConto (2012, Eq. (9)):
         IF(-q_gl(j,i) < flux_yd(j,i)) THEN
          flux_yd(j,i) = -q_gl(j,i)
         ELSE
          IF(choice_heuristic_flux_rule_at_grounding_line == 1) flux_yd(j-1,i) = -q_gl(j,i)           ! Pollard & DeCconto (2012, GMD) heuristic rule
          !IF(choice_heuristic_flux_rule_at_grounding_line == 1) flux_yu(j-1,i) = -q_gl(j,i)           ! Pollard & DeCconto (2012, GMD) heuristic rule
         END IF
        END IF
        
     END IF

    END DO
    END DO
  END SUBROUTINE schoof_grounding_line_flux_down_stream_buttressing_stress
  SUBROUTINE buttressing(Hi, A_flow_mean, dUs_dx, dUs_dy, dVs_dx, dVs_dy, buttressing_term_x, buttressing_term_y)
    ! Calculate the buttressing term which is used in the Schoof parameterisation.
    
    USE parameters_module, ONLY: n_flow, grav, ice_density, seawater_density
    
    IMPLICIT NONE

    ! Input variables: 
    REAL(dp), INTENT(IN)            :: Hi
    REAL(dp), INTENT(IN)            :: A_flow_mean
    REAL(dp), INTENT(IN)            :: dUs_dx
    REAL(dp), INTENT(IN)            :: dUs_dy
    REAL(dp), INTENT(IN)            :: dVs_dx
    REAL(dp), INTENT(IN)            :: dVs_dy

    ! Output variables:
    REAL(dp), INTENT(OUT), OPTIONAL :: buttressing_term_x
    REAL(dp), INTENT(OUT), OPTIONAL :: buttressing_term_y

    ! Local variables:
    REAL(dp)                        :: tau_stress_free     ! Free stress without buttressing
    REAL(dp)                        :: common_tau_part     ! The common part in the stress in the x- and y-direction
    
    REAL(dp), PARAMETER :: epsilon_sq_0 = 1E-12_dp
    REAL(dp), PARAMETER :: C_m_sliding = 1._dp / 3._dp
    INTEGER,  PARAMETER :: choice_regulation_of_buttressing_range = 1
    INTEGER,  PARAMETER :: choice_buttressing_method = 1
    REAL(dp), PARAMETER :: buttressing_enhancement_factor = 1._dp

    ! The free stress in absence of any buttressing (Eq 8, Pollard & DeConto, 2012): tau_f = 0.5 rho_i * g * h_g * (1- rho_i / rho_w).
    SELECT CASE(choice_regulation_of_buttressing_range)
    CASE(1)
     tau_stress_free = 0.5_dp * ice_density * grav *         Hi                    * (1._dp - ice_density / seawater_density)
    CASE(2)
     tau_stress_free = 0.5_dp * ice_density * grav * max(min(Hi,1000._dp),100._dp) * (1._dp - ice_density / seawater_density)   ! Regularisation on Hi (Pattyn, 2017, Eq. (9))
    CASE DEFAULT
     STOP
    END SELECT

    SELECT CASE(choice_buttressing_method)
    CASE(1,3,5)
     ! The longitudinal deviatoric stress (See Pollard & DeConto (2012, GMD) Eq 8):   ! Note: Pollard & DeConto, 2012 take (0.01_dp / C%dx)^2 instead of C%epsilon_sq_0
     common_tau_part = abs(A_flow_mean**(-1._dp / n_flow) * ( dUs_dx**2 + dVs_dy**2 + dUs_dx * dVs_dy + 0.25_dp * (dUs_dy + dVs_dx)**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow)) )
    CASE(2,4,6)
     ! The longitudinal deviatoric stress (See Pollard & DeConto (2012, GMD) Eq 8):
     ! A fixed A_flow_mean value is prescibed, since artefacts across the grounding line due to artefacts in Ti which in turn are due to artefacts in W (namely in W(zeta=1) at the bottom)
     common_tau_part = abs(  3.e-16_dp**(-1._dp / n_flow) * ( dUs_dx**2 + dVs_dy**2 + dUs_dx * dVs_dy + 0.25_dp * (dUs_dy + dVs_dx)**2 + epsilon_sq_0)**((1._dp - n_flow) / (2._dp * n_flow)) )
    CASE DEFAULT
     STOP
    END SELECT

    SELECT CASE(choice_buttressing_method)
    CASE(1,2)
     ! The ratio gives the buttressing term theta for both directions:
     IF(PRESENT(buttressing_term_x)) buttressing_term_x = ((buttressing_enhancement_factor * common_tau_part * dUs_dx) / tau_stress_free)**(n_flow / (C_m_sliding + 1._dp))
     IF(PRESENT(buttressing_term_y)) buttressing_term_y = ((buttressing_enhancement_factor * common_tau_part * dVs_dy) / tau_stress_free)**(n_flow / (C_m_sliding + 1._dp))
    CASE(3,4)
     ! The ratio gives the buttressing term theta for both directions:
     IF(PRESENT(buttressing_term_x)) buttressing_term_x = ((buttressing_enhancement_factor * common_tau_part * dUs_dx - 0.5_dp * ice_density * grav * Hi) / tau_stress_free)**(n_flow / (C_m_sliding + 1._dp))  ! From deviatoric stress to full stress
     IF(PRESENT(buttressing_term_y)) buttressing_term_y = ((buttressing_enhancement_factor * common_tau_part * dVs_dy - 0.5_dp * ice_density * grav * Hi) / tau_stress_free)**(n_flow / (C_m_sliding + 1._dp))  ! From deviatoric stress to full stress
    CASE(5,6)
     IF(PRESENT(buttressing_term_x)) buttressing_term_x = ((buttressing_enhancement_factor * common_tau_part * (2._dp * dUs_dx +         dVs_dy)) / tau_stress_free)**(n_flow / (C_m_sliding + 1._dp))
     IF(PRESENT(buttressing_term_y)) buttressing_term_y = ((buttressing_enhancement_factor * common_tau_part * (        dUs_dx + 2._dp * dVs_dy)) / tau_stress_free)**(n_flow / (C_m_sliding + 1._dp))
    END SELECT

    ! Limit the buttressing_term values within a range from 0 to 1:
    IF(PRESENT(buttressing_term_x)) buttressing_term_x = max(min(buttressing_term_x, 1.0_dp), 0.0_dp)
    IF(PRESENT(buttressing_term_y)) buttressing_term_y = max(min(buttressing_term_y, 1.0_dp), 0.0_dp)
  END SUBROUTINE buttressing
  
  SUBROUTINE solve_SSA_ANICE_stag( grid, ice)
    ! Calculate ice velocities using the SSA
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    LOGICAL                                            :: set_SSA_velocities_to_zero
    LOGICAL                                            :: has_converged
    INTEGER                                            :: outer_loop_i
    REAL(dp)                                           :: norm_dN, norm_N, RN
    
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
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in solve_SSA!'
        STOP
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! If there's no grounded ice anywhere, don't bother
    IF (SUM( ice%mask_sheet_Aa) == 0) set_SSA_velocities_to_zero = .TRUE.
    
    IF (set_SSA_velocities_to_zero) THEN
      ice%U_vav_SSA_Acx  = 0._dp
      ice%V_vav_SSA_Acy  = 0._dp
      ice%U_vav_SSA_Aa   = 0._dp
      ice%V_vav_SSA_Aa   = 0._dp
      RETURN
    END IF
    
    ! Calculate the basal yield stress tau_c
    CALL calculate_basal_yield_stress_Ab( grid, ice%Hi_Ab, ice%Hb_Ab, ice%SL_Ab, ice%phi_fric_Ab, ice%tau_c_Ab)
    
    ! Analytical solution at the grounding line
    !CALL SSA_analytical_GL_velocity_stag( grid, ice)
    
    ! DENK DROM
    ice%dummy2D_02(:,1:grid%nx-1) = ice%mask_GL_Acx
    ice%dummy2D_03(1:grid%ny-1,:) = ice%mask_GL_Acy
    ice%dummy2D_04(:,1:grid%nx-1) = ice%U_vav_SSA_Acx
    ice%dummy2D_05(1:grid%ny-1,:) = ice%V_vav_SSA_Acy
    !RETURN
    
    ! The viscosity iteration
    has_converged = .FALSE.
    outer_loop_i  = 0
    DO WHILE ((.NOT. has_converged) .AND. (outer_loop_i < C%SSA_max_outer_loops))
      outer_loop_i = outer_loop_i + 1
      
      ! Store the previous solution, so we can check for converge later
      ice%N_Ab_prev  = ice%N_Ab
    
      ! Update the effective viscosity eta and the product term N = eta * H
      CALL SSA_effective_viscosity_stag( grid, ice)
    
      ! Calculate the sliding term S = tau_c / norm([U,V])
      CALL SSA_sliding_term_stag( grid, ice)
      
      ! Solve the linearised SSA using SOR
      CALL solve_SSA_linearised_ANICE_stag( grid, ice)
      
      ! Check if a stable solution has been reached
      norm_N  = SQRT( SUM( ice%N_Ab**2))
      norm_dN = SQRT( SUM( (ice%N_Ab - ice%N_Ab_prev)**2 ))
      RN = norm_dN / norm_N
      
      !WRITE(0,'(A,I5,A,E10.4)') '   SSA - viscosity iteration ', outer_loop_i, ', RN = ', RN
      
      IF (RN < C%SSA_RN_tol) THEN
        has_converged = .TRUE.
      END IF
            
    END DO ! DO WHILE ((.NOT. has_converged) .AND. (outer_loop_i < C%SSA_max_outer_loops))
    
    ! Finally, map velocities to the staggered grid for use in the mass continuity integration
    CALL map_Acx_to_Aa_2D( grid, ice%U_vav_SSA_Acx, ice%U_vav_SSA_Aa)
    CALL map_Acy_to_Aa_2D( grid, ice%V_vav_SSA_Acy, ice%V_vav_SSA_Aa)
    
    ! And save the absolute value of velocity in one of the dummy output fields for easy inspection
    ice%dummy2D_01 = SQRT( ice%U_vav_SSA_Aa**2 + ice%V_vav_SSA_Aa**2)
    
  END SUBROUTINE solve_SSA_ANICE_stag
  SUBROUTINE solve_SSA_linearised_ANICE_stag( grid, ice)
    ! Calculate ice velocities using the SSA
    
    USE parameters_module,          ONLY: ice_density, grav
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j,m,j0
    REAL(dp), DIMENSION( grid%ny  , grid%nx-1)         :: LHSx
    REAL(dp), DIMENSION( grid%ny-1, grid%nx  )         :: LHSy
    REAL(dp), DIMENSION( grid%ny  , grid%nx-1)         :: RHSx
    REAL(dp), DIMENSION( grid%ny-1, grid%nx  )         :: RHSy
    REAL(dp), DIMENSION( grid%ny  , grid%nx-1)         :: resU
    REAL(dp), DIMENSION( grid%ny-1, grid%nx  )         :: resV
    REAL(dp), DIMENSION( grid%ny  , grid%nx-1)         :: eu_ij
    REAL(dp), DIMENSION( grid%ny-1, grid%nx  )         :: ev_ij
    INTEGER,  DIMENSION( grid%ny  , grid%nx-1)         :: mask_sliding_Acx
    INTEGER,  DIMENSION( grid%ny-1, grid%nx  )         :: mask_sliding_Acy
    REAL(dp), DIMENSION( grid%ny  , grid%nx-1)         :: U
    REAL(dp), DIMENSION( grid%ny-1, grid%nx  )         :: V
    LOGICAL                                            :: has_converged
    INTEGER                                            :: inner_loop_i
    REAL(dp)                                           :: max_residual_UV
    
    C%SSA_SOR_omega = 1.4_dp
    
    ! Local copies of U and V, to make the equations more readable
    U = ice%U_vav_SSA_Acx
    V = ice%V_vav_SSA_Acy
       
    ! Calculate the right-hand sides of the equations (i.e. the gravitational driving stress) on Acx/Acy
    DO i = 2, grid%nx-2
    DO j = 2, grid%ny-1
      RHSx( j,i) = ice_density * grav * ice%dHs_dx_Acx( j,i) / (0.5_dp * (ice%eta_Ab( j-1,i) + ice%eta_Ab( j,i)))
    END DO
    END DO
    DO i = 2, grid%nx-1
    DO j = 2, grid%ny-2
      RHSy( j,i) = ice_density * grav * ice%dHs_dy_Acy( j,i) / (0.5_dp * (ice%eta_Ab( j,i-1) + ice%eta_Ab( j,i)))
    END DO
    END DO
    
    ! Get "basal sliding mask" (i.e. all grounded grid points where basal traction > 0) on Acx/Acy
    mask_sliding_Acx = 1
    DO i = 1, grid%nx-1
    DO j = 1, grid%ny
      IF ((ice%mask_shelf_Aa( j,i) == 1 .OR. ice%mask_ocean_Aa( j,i) == 1) .AND. (ice%mask_shelf_Aa( j,i+1) == 1 .OR. ice%mask_ocean_Aa( j,i+1) == 1)) THEN
        mask_sliding_Acx( j,i) = 0
      END IF
    END DO
    END DO
    mask_sliding_Acy = 1
    DO i = 1, grid%nx
    DO j = 1, grid%ny-1
      IF ((ice%mask_shelf_Aa( j,i) == 1 .OR. ice%mask_ocean_Aa( j,i) == 1) .AND. (ice%mask_shelf_Aa( j+1,i) == 1 .OR. ice%mask_ocean_Aa( j+1,i) == 1)) THEN
        mask_sliding_Acy( j,i) = 0
      END IF
    END DO
    END DO
    
    ! Calculate the centre coefficients on Acx/Acy
    DO i = 2, grid%nx-2
    DO j = 2, grid%ny-1
      IF (mask_sliding_Acx( j,i) == 1) THEN
        eu_ij( j,i) = (-8._dp / grid%dx**2) + (-2._dp / grid%dy**2) - (ice%S_Ab( j-1,i)+ice%S_Ab( j,i)) / (MAX(0.1_dp,ice%Hi_Acx( j,i)) * (ice%eta_Ab( j-1,i)+ice%eta_Ab( j,i)))
      ELSE
        eu_ij( j,i) = (-8._dp / grid%dx**2) + (-2._dp / grid%dy**2)
      END IF
    END DO
    END DO
    DO i = 2, grid%nx-1
    DO j = 2, grid%ny-2
      IF (mask_sliding_Acy( j,i) == 1) THEN
        ev_ij( j,i) = (-2._dp / grid%dx**2) + (-8._dp / grid%dy**2) - (ice%S_Ab( j,i-1)+ice%S_Ab( j,i)) / (MAX(0.1_dp,ice%Hi_Acy( j,i)) * (ice%eta_Ab( j,i-1)+ice%eta_Ab( j,i)))
      ELSE
        ev_ij( j,i) = (-2._dp / grid%dx**2) + (-8._dp / grid%dy**2)
      END IF
    END DO
    END DO
    
    ! SOR iteration
    has_converged = .FALSE.
    inner_loop_i  = 0
    DO WHILE ((.NOT. has_converged) .AND. (inner_loop_i < C%SSA_max_inner_loops))
      inner_loop_i = inner_loop_i + 1
      
      LHSx = 0._dp
      LHSy = 0._dp
      resU = 0._dp
      resV = 0._dp
      
      ! Red-black scheme
      DO m = 0, 1
      
        j0 = m
        DO i = 2, grid%nx-1
        j0 = 1 - j0
        
        DO j = j0+2, grid%ny-1, 2
        
          ! Equation 1
          IF ( i <= grid%nx-2) THEN
            
            !IF (ice%mask_GL_Acx( j,i) == 1) CYCLE
            
            ! Calculate the left-hand side of the equation, using partially updated values
            LHSx( j,i) = (4._dp * (U( j  ,i+1) + U( j  ,i-1)) / grid%dx**2) + &
                         (        (U( j+1,i  ) + U( j-1,i  )) / grid%dy**2) + &
                         (3._dp * (V( j  ,i+1) + V( j-1,i  ) - V( j  ,i  ) - V( j-1,i+1)) / (4._dp*grid%dx*grid%dy)) + (eu_ij( j,i) * U( j,i))
                         
            
            ! Calculate the residuals
            resU( j,i) = (LHSx( j,i) - RHSx( j,i)) / eu_ij( j,i)
            
            ! Update velocities
            U( j,i) = U( j,i) - C%SSA_SOR_omega * resU( j,i)
          END IF
        
          ! Equation 2
          IF ( j <= grid%ny-2) THEN
            
            !IF (ice%mask_GL_Acy( j,i) == 1) CYCLE
            
            ! Calculate the left-hand side of the equation, using partially updated values
            LHSy( j,i) = (4._dp * (V( j+1,i  ) + V( j-1,i  )) / grid%dy**2) + &
                         (        (V( j  ,i+1) + V( j  ,i-1)) / grid%dx**2) + &
                         (3._dp * (U( j+1,i  ) + U( j  ,i-1) - U( j+1,i-1) - U( j  ,i  )) / (4._dp*grid%dx*grid%dy)) + (ev_ij( j,i) * V( j,i))
            
            ! Calculate the residuals
            resV( j,i) = (LHSy( j,i) - RHSy( j,i)) / ev_ij( j,i)
            
            ! Update velocities
            V( j,i) = V( j,i) - C%SSA_SOR_omega * resV( j,i)
          END IF
        
        END DO ! DO j = j0+2, grid%ny-1, 2
        END DO ! DO i = 2, grid%nx-1
      
      END DO ! DO m = 0, 1
      
      ! Apply Neumann boundary conditions
      U( :        ,1        ) = U( :        ,2        )
      U( :        ,grid%nx-1) = U( :        ,grid%nx-2)
      U( 1        ,:        ) = U( 2        ,:        )
      U( grid%ny  ,:        ) = U( grid%ny-1,:        )
      V( :        ,1        ) = V( :        ,2        )
      V( :        ,grid%nx  ) = V( :        ,grid%nx-1)
      V( 1        ,:        ) = V( 2        ,:        )
      V( grid%ny-1,:        ) = V( grid%ny-2,:        )
      
      ! Check if we've reached a stable solution
      max_residual_UV = MAX( MAXVAL( ABS(resU(2:grid%ny-1,2:grid%nx-2))), MAXVAL( ABS(resV(2:grid%ny-2,2:grid%nx-1))))
      !WRITE(0,*) ' SSA -  inner loop ', inner_loop_i, ': largest residual = ', max_residual_UV
      IF (max_residual_UV < C%SSA_max_residual_UV) THEN
        has_converged = .TRUE.
      ELSEIF (max_residual_UV > 1E8_dp) THEN
        WRITE(0,*) ' ERROR - instability in SSA SOR solver!'
        STOP
      ELSEIF (inner_loop_i == C%SSA_max_inner_loops) THEN
        WRITE(0,*) ' ERROR - SSA SOR solver doesnt converge!'
        STOP
      END IF
      
!      ! DENK DROM
!      ice%dummy2D_01 = 0._dp
!      ice%dummy2D_02 = 0._dp
!      ice%dummy2D_03 = 0._dp
!      ice%dummy2D_04 = 0._dp
!      ice%dummy2D_05 = 0._dp
!      ice%dummy2D_06 = 0._dp
!      ice%dummy2D_07 = 0._dp
!      ice%dummy2D_08 = 0._dp
!      ice%dummy2D_09 = 0._dp
!      ice%dummy2D_10 = 0._dp
!      ice%dummy3D_01 = 0._dp
!      ice%dummy2D_01(   1:grid%ny-1,1:grid%nx-1) = ice%A_flow_mean_Ab
!      ice%dummy2D_02(   1:grid%ny-1,1:grid%nx-1) = ice%eta_Ab
!      ice%dummy2D_03(   1:grid%ny-1,1:grid%nx-1) = ice%tau_c_Ab
!      ice%dummy2D_04(   1:grid%ny-1,1:grid%nx-1) = ice%S_Ab
!      ice%dummy2D_05(   :          ,1:grid%nx-1) = RHSx
!      ice%dummy2D_06(   1:grid%ny-1,:          ) = RHSy
!      ice%dummy2D_07(   :          ,1:grid%nx-1) = eu_ij
!      ice%dummy2D_08(   1:grid%ny-1,:          ) = ev_ij
!      ice%dummy2D_09(   :          ,1:grid%nx-1) = LHSx
!      ice%dummy2D_10(   1:grid%ny-1,:          ) = LHSy
!      ice%dummy3D_01( 1,:          ,1:grid%nx-1) = resU
!      ice%dummy3D_01( 2,1:grid%ny-1,:          ) = resV
!      IF (inner_loop_i == 2) RETURN
      
    END DO ! DO WHILE ((.NOT. has_converged) .AND. (inner_loop_i < max_inner_loops))
    
    ! Copy data back from local copies to ice structure:
    ice%U_vav_SSA_Acx = U
    ice%V_vav_SSA_Acy = V
    
  END SUBROUTINE solve_SSA_linearised_ANICE_stag
  SUBROUTINE SSA_effective_viscosity_stag( grid, ice)
    ! Calculate the effective viscosity eta, the product term N = eta*H, and the gradient of N in the SSA
    
    USE parameters_module, ONLY: n_flow
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION( grid%ny-1, grid%nx-1)         :: Ux
    REAL(dp), DIMENSION( grid%ny-1, grid%nx-1)         :: Uy
    REAL(dp), DIMENSION( grid%ny-1, grid%nx-1)         :: Vx
    REAL(dp), DIMENSION( grid%ny-1, grid%nx-1)         :: Vy
    
    REAL(dp), PARAMETER                                :: epsilon_sq_0 = 1E-12_dp
    
    ! Calculate derivatives of U and V
    CALL ddx_Acx_to_Ab_2D( grid, ice%U_vav_SSA_Acx, Ux)
    CALL ddy_Acx_to_Ab_2D( grid, ice%U_vav_SSA_Acx, Uy)
    CALL ddx_Acy_to_Ab_2D( grid, ice%V_vav_SSA_Acy, Vx)
    CALL ddy_Acy_to_Ab_2D( grid, ice%V_vav_SSA_Acy, Vy)
    
    ! The effective viscosity eta on both grids, with a normalisation term (following Bueler & Brown, 2009) to prevent divide-by-zero errors,
    ! and the product term N = eta * H
    DO i = 1, grid%nx-1
    DO j = 1, grid%ny-1
      ice%eta_Ab( j,i) = 0.5_dp * ice%A_flow_mean_Ab( j,i)**(-1._dp/n_flow) * (&
        Ux( j,i)**2 + Vy( j,i)**2 + Ux( j,i)*Vy( j,i) + 0.25_dp*(Uy( j,i)+Vx( j,i))**2 + epsilon_sq_0) ** ((1._dp - n_flow) / (2._dp * n_flow))
      ice%N_Ab( j,i) = ice%eta_Ab( j,i) * MAX(0.1_dp, ice%Hi_Ab( j,i))
    END DO
    END DO
    
  END SUBROUTINE SSA_effective_viscosity_stag
  SUBROUTINE SSA_sliding_term_stag( grid, ice)
    ! Calculate the sliding term S = tau_c / norm([U,V]) in the SSA
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    REAL(dp), PARAMETER                                :: delta_v = 1E-3_dp
    REAL(dp), PARAMETER                                :: q_plastic            = 0.30_dp       ! Parameter used for basal stress (inverse of m_flow)
    REAL(dp), PARAMETER                                :: u_threshold          = 100._dp       ! scaling of tau_yield to get the correct unit (function of q_plastic)
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(grid%ny-1,grid%nx-1)           :: U_Ab, V_Ab
    
    ! Map ice velocities to the B grid
    CALL map_Acx_to_Ab_2D( grid, ice%U_vav_SSA_Acx, U_Ab)
    CALL map_Acy_to_Ab_2D( grid, ice%V_vav_SSA_Acy, V_Ab)
    
    ! Get beta_base on the B grid
    ice%S_Ab = 0._dp
    DO i = 1, grid%nx-1
    DO j = 1, grid%ny-1
      ! The sliding term S, with a normalisation term (following Bueler & Brown, 2009) to prevent divide-by-zero errors.
      ice%S_Ab( j,i) = ice%tau_c_Ab( j,i) * ( (delta_v**2 + U_Ab( j,i)**2 + V_Ab( j,i)**2)**(0.5_dp * (q_plastic-1._dp)) ) / (u_threshold**q_plastic)
    END DO
    END DO
    
  END SUBROUTINE SSA_sliding_term_stag
  SUBROUTINE SSA_analytical_GL_velocity_stag( grid, ice)
    ! Calculate the effective viscosity eta, the product term N = eta*H, and the gradient of N in the SSA
    
    USE parameters_module, ONLY: seawater_density, ice_density, grav, n_flow 
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                        INTENT(IN)    :: grid
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: factor_Tsai
    REAL(dp)                                           :: TAF_left
    REAL(dp)                                           :: TAF_right
    REAL(dp)                                           :: xp_gl
    REAL(dp)                                           :: Hi
    REAL(dp)                                           :: phi_fric
    REAL(dp)                                           :: A_flow
    REAL(dp)                                           :: Q_GL
    REAL(dp)                                           :: hx, hy
    REAL(dp)                                           :: U_GL, V_GL
    REAL(dp)                                           :: Fx, Fy, F

    ! Some constant terms for the Schoof grounding line flux
    REAL(dp), PARAMETER                                :: Q0 = 0.61_dp
    
    ! First in the x direction
    DO i = 1, grid%nx-1
    DO j = 1, grid%ny
    
      IF (ice%mask_GL_Acx( j,i) == 1) THEN
        
        ! Interpolate ice thickness, till friction angle and ice flow factor
        TAF_left  = ice%Hi_Aa( j,i  ) - ((ice%SL_Aa( j,i  ) - ice%Hb_Aa( j,i  )) * (seawater_density / ice_density))
        TAF_right = ice%Hi_Aa( j,i+1) - ((ice%SL_Aa( j,i+1) - ice%Hb_Aa( j,i+1)) * (seawater_density / ice_density))
        xp_gl = TAF_left / (TAF_left - TAF_right)
        
        Hi       = (ice%Hi_Aa(          j,i) * (1._dp - xp_gl)) + (ice%Hi_Aa(          j,i+1) * xp_gl)
        A_flow   = (ice%A_flow_mean_Aa( j,i) * (1._dp - xp_gl)) + (ice%A_flow_mean_Aa( j,i+1) * xp_gl)
        phi_fric = 0.5_dp * (ice%phi_fric_Ab( j-1,i) + ice%phi_fric_Ab( j,i))

        ! The constant factor in the Tsai et al. (2015) grounding line flux solution:
        factor_Tsai   = ( 8._dp * Q0 * A_flow &
                         * ((ice_density * grav) ** n_flow) &
                         * ((1._dp - (ice_density / seawater_density)) ** (n_flow-1._dp)) &
                         / (4._dp ** n_flow) );

        ! Calculate grounding line flux
        Q_GL = factor_Tsai * (Hi**(n_flow + 2._dp)) / TAN(phi_fric * C%deg2rad)

        ! The flux is outward-perpendicular to the grounding line, and therefore antiparallel
        ! to the gradient of the thickness-above-flotation
        Fx = -(ice%dHi_dx_Acx( j,i) - ((ice%dSL_dx_Acx( j,i) - ice%dHb_dx_Acx( j,i)) * (seawater_density / ice_density)))
        Fy = -(ice%dHi_dy_Acx( j,i) - ((ice%dSL_dy_Acx( j,i) - ice%dHb_dy_Acx( j,i)) * (seawater_density / ice_density)))
        F  = NORM2([Fx,Fy])
        Fx = Fx / F
        Fy = Fy / F
  
        ! Analytical ice velocities (antiparallel to TAF gradient)
        ice%U_vav_SSA_Acx( j,i) = Q_GL * Fx / Hi
        
      END IF ! IF (is grounding line in x direction)
    
    END DO
    END DO
    
    ! Then in the y direction
    DO i = 1, grid%nx
    DO j = 1, grid%ny-1
    
      IF (ice%mask_GL_Acy( j,i) == 1) THEN
        
        ! Interpolate ice thickness, till friction angle and ice flow factor
        TAF_left  = ice%Hi_Aa( j  ,i) - ((ice%SL_Aa( j  ,i) - ice%Hb_Aa( j  ,i)) * (seawater_density / ice_density))
        TAF_right = ice%Hi_Aa( j+1,i) - ((ice%SL_Aa( j+1,i) - ice%Hb_Aa( j+1,i)) * (seawater_density / ice_density))
        xp_gl = TAF_left / (TAF_left - TAF_right)
        
        Hi       = (ice%Hi_Aa(          j,i) * (1._dp - xp_gl)) + (ice%Hi_Aa(          j+1,i) * xp_gl)
        A_flow   = (ice%A_flow_mean_Aa( j,i) * (1._dp - xp_gl)) + (ice%A_flow_mean_Aa( j+1,i) * xp_gl)
        phi_fric = 0.5_dp * (ice%phi_fric_Ab( j,i-1) + ice%phi_fric_Ab( j,i))

        ! The constant factor in the Tsai et al. (2015) grounding line flux solution:
        factor_Tsai   = ( 8._dp * Q0 * A_flow &
                         * ((ice_density * grav) ** n_flow) &
                         * ((1._dp - (ice_density / seawater_density)) ** (n_flow-1._dp)) &
                         / (4._dp ** n_flow) );

        ! Calculate grounding line flux
        Q_GL = factor_Tsai * (Hi**(n_flow + 2._dp)) / TAN(phi_fric * C%deg2rad)

        ! The flux is outward-perpendicular to the grounding line, and therefore antiparallel
        ! to the gradient of the thickness-above-flotation
        Fx = -(ice%dHi_dx_Acy( j,i) - ((ice%dSL_dx_Acy( j,i) - ice%dHb_dx_Acy( j,i)) * (seawater_density / ice_density)))
        Fy = -(ice%dHi_dy_Acy( j,i) - ((ice%dSL_dy_Acy( j,i) - ice%dHb_dy_Acy( j,i)) * (seawater_density / ice_density)))
        F  = NORM2([Fx,Fy])
        Fx = Fx / F
        Fy = Fy / F
  
        ! Analytical ice velocities (antiparallel to TAF gradient)
        ice%V_vav_SSA_Acy( j,i) = Q_GL * Fy / Hi
        
      END IF ! IF (is grounding line in x direction)
    
    END DO
    END DO
    
    
  END SUBROUTINE SSA_analytical_GL_velocity_stag
  
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
    REAL(dp), DIMENSION(grid%ny,grid%nx)               :: D_SIA_Aa_from_Acx, D_SIA_Aa_from_Acy
    
    REAL(dp), PARAMETER                                :: D_uv_3D_cutoff = -1E5_dp  
        
    ice%D_SIA_3D_Acx = 0._dp
    ice%D_SIA_3D_Acy = 0._dp
    ice%D_SIA_Acx    = 0._dp
    ice%D_SIA_Acy    = 0._dp
    ice%D_SIA_Aa     = 0._dp
    
    ! Calculate 3D ice diffusivity
    DO i = 1, grid%nx-1
    DO j = 1, grid%ny-1
    
      IF (ice%mask_sheet_Acx( j,i) == 1 .OR. ice%mask_gl_Acx( j,i) == 1) THEN
        
        D_0 = (ice_density * grav * ice%Hi_Acx( j,i))**n_flow * ((ice%dHs_dx_Acx( j,i)**2 + ice%dHs_dy_Acx( j,i)**2))**((n_flow - 1._dp) / 2._dp)       
        D_deformation = vertical_integration_from_bottom_to_zeta( C%m_enh_sia * ice%A_flow_Acx( :,j,i) * C%zeta**n_flow)
        D_deformation = 2._dp * ice%Hi_Acx( j,i) * D_deformation
       
        ice%D_SIA_3D_Acx( :,j,i) = D_0 * D_deformation

      END IF
    
      IF (ice%mask_sheet_Acy( j,i) == 1 .OR. ice%mask_gl_Acy( j,i) == 1) THEN
        
        D_0 = (ice_density * grav * ice%Hi_Acy( j,i))**n_flow * ((ice%dHs_dx_Acy( j,i)**2 + ice%dHs_dy_Acy( j,i)**2))**((n_flow - 1._dp) / 2._dp)       
        D_deformation = vertical_integration_from_bottom_to_zeta( C%m_enh_sia * ice%A_flow_Acy( :,j,i) * C%zeta**n_flow)
        D_deformation = 2._dp * ice%Hi_Acy( j,i) * D_deformation
       
        ice%D_SIA_3D_Acy( :,j,i) = D_0 * D_deformation

      END IF

    END DO
    END DO
      
    ! Check for very large D_uv_3D's, causing very large velocities RESET --DIRTY
    ice%D_SIA_3D_Acx = MAX( ice%D_SIA_3D_Acx, D_uv_3D_cutoff)
    ice%D_SIA_3D_Acy = MAX( ice%D_SIA_3D_Acy, D_uv_3D_cutoff)
    
    ! Calculate vertically averaged ice velocities on the Ac grids
    DO i = 1, grid%nx-1
    DO j = 1, grid%ny-1
    
      IF (ice%mask_sheet_Acx( j,i) == 1 .OR. ice%mask_gl_Acx( j,i) == 1) THEN
        D_SIA_3D_prof = ice%D_SIA_3D_Acx( :,j,i)
        D_uv_2D = vertical_average( D_SIA_3D_prof)
       
        ice%D_SIA_Acx(     j,i) = ice%Hi_Acx( j,i) * D_uv_2D
        ice%U_vav_SIA_Acx( j,i) = D_uv_2D * ice%dHs_dx_Acx( j,i)
      END IF
    
      IF (ice%mask_sheet_Acy( j,i) == 1 .OR. ice%mask_gl_Acy( j,i) == 1) THEN
        D_SIA_3D_prof = ice%D_SIA_3D_Acy( :,j,i)
        D_uv_2D = vertical_average( D_SIA_3D_prof)
       
        ice%D_SIA_Acy(     j,i) = ice%Hi_Acy( j,i) * D_uv_2D
        ice%V_vav_SIA_Acy( j,i) = D_uv_2D * ice%dHs_dy_Acy( j,i)
      END IF
          
    END DO
    END DO
    
    ! Map data to Aa grid for writing to output (not actually used anywhere...)
    CALL map_Acx_to_Aa_2D( grid, ice%D_SIA_Acx, D_SIA_Aa_from_Acx)
    CALL map_Acy_to_Aa_2D( grid, ice%D_SIA_Acy, D_SIA_Aa_from_Acy)
    ice%D_SIA_Aa = (D_SIA_Aa_from_Acx + D_SIA_Aa_from_Acy) / 2._dp
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
                                        
    ice%U_3D_Aa = 0._dp
    ice%V_3D_Aa = 0._dp
    ice%W_3D_Aa = 0._dp

    ! Direct calculation of U and V for the ice sheet:
    DO i = 1, grid%nx
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
        
    CALL Neumann_BC_Aa_3D( grid, ice%U_3D_Aa)
    CALL Neumann_BC_Aa_3D( grid, ice%V_3D_Aa)
    
    CALL ddx_Aa_to_Aa_3D( grid, ice%U_3D_Aa, dU_dx)
    CALL ddy_Aa_to_Aa_3D( grid, ice%U_3D_Aa, dU_dy)
    CALL ddx_Aa_to_Aa_3D( grid, ice%V_3D_Aa, dV_dx)
    CALL ddy_Aa_to_Aa_3D( grid, ice%V_3D_Aa, dV_dy)

    DO i = 2, grid%nx-1
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
    
    CALL Neumann_BC_Aa_3D( grid, ice%W_3D_Aa)

  END SUBROUTINE solve_SIA_3D
  
  ! SSA solver
  SUBROUTINE solve_SSA( grid, ice)
    ! Calculate ice velocities using the SSA
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    LOGICAL                                            :: set_SSA_velocities_to_zero
    
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
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in solve_SSA!'
        STOP
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! If there's no grounded ice anywhere, don't bother
    IF (SUM( ice%mask_sheet_Aa) == 0) set_SSA_velocities_to_zero = .TRUE.
    
    IF (set_SSA_velocities_to_zero) THEN
      ice%U_vav_SSA_Acx  = 0._dp
      ice%V_vav_SSA_Acy  = 0._dp
      ice%U_vav_SSA_Aa   = 0._dp
      ice%V_vav_SSA_Aa   = 0._dp
      RETURN
    END IF
    
    ! Calculate the basal yield stress tau_c
    CALL basal_yield_stress( grid, ice)
    
    ! First, solve the non-linear, "unconstrained" SSA
    ice%mask_U_analytical = 0
    ice%mask_V_analytical = 0
    CALL solve_SSA_nonlinear( grid, ice, ice%mask_U_analytical, ice%mask_V_analytical, ice%U_SSA_free, ice%V_SSA_free)
    
   ! IF (.NOT. C%use_analytical_GL_flux) THEN
      ice%U_SSA_cons = ice%U_SSA_free
      ice%V_SSA_cons = ice%V_SSA_free
   ! ELSE
   !   ! First, calculate the semi-analytical grounding-line flux, using the P&dC heuristic for where to apply it
   !   CALL analytical_GL_flux( grid, ice, ice%U_SSA_free, ice%V_SSA_free, ice%mask_U_analytical, ice%mask_V_analytical, ice%U_SSA_cons, ice%V_SSA_cons)
   !   ! Solve the non-linear SSA again, with the analytical solution as a constraint
   !   CALL solve_SSA_nonlinear( grid, ice, ice%mask_U_analytical, ice%mask_V_analytical, ice%U_SSA_cons, ice%V_SSA_cons)
   ! END IF
    
!    ! DENK DROM
!    ice%dummy2D_02 = ice%U_SSA_free
!    ice%dummy2D_03 = ice%V_SSA_free
!    ice%dummy2D_04 = SQRT(ice%U_SSA_free**2+ice%V_SSA_free**2)
!    ice%dummy2D_05 = REAL(ice%mask_U_analytical,dp)
!    ice%dummy2D_06 = REAL(ice%mask_V_analytical,dp)
!    ice%dummy2D_07 = ice%U_SSA_cons
!    ice%dummy2D_08 = ice%V_SSA_cons
!    ice%dummy2D_09 = SQRT(ice%U_SSA_cons**2+ice%V_SSA_cons**2)
    
!    ice%U_vav_SSA_Aa = ice%U_SSA_free
!    ice%V_vav_SSA_Aa = ice%V_SSA_free
    ice%U_vav_SSA_Aa = ice%U_SSA_cons
    ice%V_vav_SSA_Aa = ice%V_SSA_cons
    
    ! Finally, map velocities to the staggered grid for use in the mass continuity integration
    CALL map_Aa_to_Acx_2D( grid, ice%U_vav_SSA_Aa, ice%U_vav_SSA_Acx)
    CALL map_Aa_to_Acy_2D( grid, ice%V_vav_SSA_Aa, ice%V_vav_SSA_Acy)
    
    ! And save the absolute value of velocity in one of the dummy output fields for easy inspection
    ice%dummy2D_01 = SQRT( ice%U_vav_SSA_Aa**2 + ice%V_vav_SSA_Aa**2)
    
  END SUBROUTINE solve_SSA
  
  SUBROUTINE solve_SSA_nonlinear( grid, ice, mask_U_analytical, mask_V_analytical, U, V)
    ! Calculate ice velocities using the SSA
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                        INTENT(IN)    :: grid
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    INTEGER,  DIMENSION( grid%ny, grid%nx), INTENT(IN   ) :: mask_U_analytical
    INTEGER,  DIMENSION( grid%ny, grid%nx), INTENT(IN   ) :: mask_V_analytical
    REAL(dp), DIMENSION( grid%ny, grid%nx), INTENT(INOUT) :: U, V
    
    ! Local variables:
    LOGICAL                                            :: has_converged
    INTEGER                                            :: outer_loop_i
    REAL(dp)                                           :: U_analytical
    REAL(dp)                                           :: norm_dN, norm_N, RN
    
    ! The viscosity iteration
    has_converged = .FALSE.
    outer_loop_i  = 0
    DO WHILE ((.NOT. has_converged) .AND. (outer_loop_i < C%SSA_max_outer_loops))
      outer_loop_i = outer_loop_i + 1
      
      ! Store the previous solution, so we can check for converge later
      ice%N_Aa_prev  = ice%N_Aa
    
      ! Update the effective viscosity eta and the product term N = eta * H
      CALL SSA_effective_viscosity( grid, ice, U, V)
    
      ! Calculate the sliding term S = tau_c / norm([U,V])
      CALL SSA_sliding_term( grid, ice, U, V)
      
      ! Solve the linearised SSA using SOR
      CALL solve_SSA_linearised_BB_SOR( grid, ice, mask_U_analytical, mask_V_analytical, U, V)
!      CALL solve_SSA_linearised_FP_SOR( grid, ice, mask_U_analytical, mask_V_analytical, U, V)
      
      ! Check if a stable solution has been reached
      norm_N  = SQRT( SUM( ice%N_Aa**2))
      norm_dN = SQRT( SUM( (ice%N_Aa - ice%N_Aa_prev)**2 ))
      RN = norm_dN / norm_N
      
      IF (C%do_benchmark_experiment .AND. C%choice_benchmark_experiment == 'SSA_icestream') THEN
        CALL SSA_Schoof2006_analytical_solution( 0.001_dp, 2000._dp, ice%A_flow_mean_Aa( 1,1), 40000._dp, C%m_SSA_icestream, 0._dp, U_analytical)
        WRITE(0,'(A,I5,A,F7.3,A,F7.3,A,E10.4,A,E10.4)') '   Benchmark experiment SSA_icestream - viscosity iteration ', outer_loop_i, &
          ': U_analytical(y=0) = ', U_analytical, ', U_modelled(y=0) = ', MAXVAL(U), ', err_rel = ', ABS(1._dp - MAXVAL(U) / U_analytical), ', RN = ', RN
      END IF
      
      !WRITE(0,'(A,I5,A,E10.4)') '   SSA - viscosity iteration ', outer_loop_i, ', RN = ', RN
      
      IF (RN < C%SSA_RN_tol) THEN
        has_converged = .TRUE.
      END IF
            
    END DO ! DO WHILE ((.NOT. has_converged) .AND. (outer_loop_i < C%SSA_max_outer_loops))
    
  END SUBROUTINE solve_SSA_nonlinear
  
  SUBROUTINE solve_SSA_linearised_BB_SOR( grid, ice, mask_U_analytical, mask_V_analytical, U, V)
    ! Use successive over-relaxation (SOR) to solve the system of linear equations that
    ! is the SSA when the non-linear terms N and S are assumed constant.
    ! Using the discretisation described in Bueler and Brown, 2009.
    
    USE parameters_module, ONLY: ice_density, grav
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                        INTENT(IN)    :: grid
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    INTEGER,  DIMENSION( grid%ny, grid%nx), INTENT(IN   ) :: mask_U_analytical
    INTEGER,  DIMENSION( grid%ny, grid%nx), INTENT(IN   ) :: mask_V_analytical
    REAL(dp), DIMENSION( grid%ny, grid%nx), INTENT(INOUT) :: U, V
    
    ! Local variables:
    INTEGER                                            :: i,j,m,j0
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: N_ip12
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: N_im12
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: N_jp12
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: N_jm12
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: LHSx
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: LHSy
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: RHSx
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: RHSy
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: resU
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: resV
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: eu_ij
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: ev_ij
    LOGICAL                                            :: has_converged
    INTEGER                                            :: inner_loop_i
    REAL(dp)                                           :: max_residual_UV
    
    !  N on the staggered grids
    N_ip12( :,1:grid%nx-1) = ice%N_Acx
    N_im12( :,2:grid%nx  ) = ice%N_Acx
    N_jp12( 1:grid%ny-1,:) = ice%N_Acy
    N_jm12( 2:grid%ny  ,:) = ice%N_Acy
       
    ! Calculate the right-hand sides of the equations (i.e. the gravitational driving stress)
    DO i = 2, grid%nx-1
    DO j = 2, grid%ny-1
      RHSx( j,i) = ice_density * grav * ice%Hi_Aa( j,i) * ice%dHs_dx_shelf_Aa( j,i)
      RHSy( j,i) = ice_density * grav * ice%Hi_Aa( j,i) * ice%dHs_dy_shelf_Aa( j,i)
      !RHSx( j,i) = ice_density * grav * ice%Hi_Aa( j,i) * ice%dHs_dx_Aa( j,i)
      !RHSy( j,i) = ice_density * grav * ice%Hi_Aa( j,i) * ice%dHs_dy_Aa( j,i)
    END DO
    END DO
    
    ! Calculate the centre coefficients
    DO i = 2, grid%nx-1
    DO j = 2, grid%ny-1
      eu_ij( j,i) = (-4._dp * (N_ip12( j,i) + N_im12( j,i)) / grid%dx**2) + (-1._dp * (N_jp12( j,i) + N_jm12( j,i)) / grid%dy**2) - ice%S_Aa( j,i)
      ev_ij( j,i) = (-1._dp * (N_ip12( j,i) + N_im12( j,i)) / grid%dx**2) + (-4._dp * (N_jp12( j,i) + N_jm12( j,i)) / grid%dy**2) - ice%S_Aa( j,i)
    END DO
    END DO
    
    ! SOR iteration
    has_converged = .FALSE.
    inner_loop_i  = 0
    DO WHILE ((.NOT. has_converged) .AND. (inner_loop_i < C%SSA_max_inner_loops))
      inner_loop_i = inner_loop_i + 1
      
      ! Red-black scheme
      DO m = 0, 1
      
        j0 = m
        DO i = 2, grid%nx-1
        j0 = 1 - j0
        
        DO j = j0+2, grid%ny-1, 2
        
          ! Keep the values from the analytical flux condition as a constraint
         ! IF (mask_Q_analytical( j,i) == 1) CYCLE
        
          ! Calculate the left-hand sides of the equations
          LHSx( j,i) = (2._dp * N_ip12( j,i) / grid%dx) * ((2._dp * (U( j  ,i+1) - U( j  ,i  )) / grid%dx) + ((V( j+1,i+1) - V( j-1,i+1) + V( j+1,i  ) - V( j-1,i  )) / (4._dp * grid%dy))) &
                     - (2._dp * N_im12( j,i) / grid%dx) * ((2._dp * (U( j  ,i  ) - U( j  ,i-1)) / grid%dx) + ((V( j+1,i  ) - V( j-1,i  ) + V( j+1,i-1) - V( j-1,i-1)) / (4._dp * grid%dy))) &
                     + (        N_jp12( j,i) / grid%dy) * ((        (U( j+1,i  ) - U( j  ,i  )) / grid%dy) + ((V( j+1,i+1) - V( j+1,i-1) + V( j  ,i+1) - V( j  ,i-1)) / (4._dp * grid%dx))) &
                     - (        N_jm12( j,i) / grid%dy) * ((        (U( j  ,i  ) - U( j-1,i  )) / grid%dy) + ((V( j  ,i+1) - V( j  ,i-1) + V( j-1,i+1) - V( j-1,i-1)) / (4._dp * grid%dx))) &
                     - (ice%S_Aa( j,i) * U( j,i))

          LHSy( j,i) = (2._dp * N_jp12( j,i) / grid%dy) * ((2._dp * (V( j+1,i  ) - V( j  ,i  )) / grid%dy) + ((U( j+1,i+1) - U( j+1,i-1) + U( j  ,i+1) - U( j  ,i-1)) / (4._dp * grid%dx))) &
                     - (2._dp * N_jm12( j,i) / grid%dy) * ((2._dp * (V( j  ,i  ) - V( j-1,i  )) / grid%dy) + ((U( j  ,i+1) - U( j  ,i-1) + U( j-1,i+1) - U( j-1,i-1)) / (4._dp * grid%dx))) &
                     + (        N_ip12( j,i) / grid%dx) * ((        (V( j  ,i+1) - V( j  ,i  )) / grid%dx) + ((U( j+1,i+1) - U( j-1,i+1) + U( j+1,i  ) - U( j-1,i  )) / (4._dp * grid%dy))) &
                     - (        N_im12( j,i) / grid%dx) * ((        (V( j  ,i  ) - V( j  ,i-1)) / grid%dx) + ((U( j+1,i  ) - U( j-1,i  ) + U( j+1,i-1) - U( j-1,i-1)) / (4._dp * grid%dy))) &
                     - (ice%S_Aa( j,i) * V( j,i))
          
          ! Calculate the residuals
          resU( j,i) = (LHSx( j,i) - RHSx( j,i)) / eu_ij( j,i)
          resV( j,i) = (LHSy( j,i) - RHSy( j,i)) / ev_ij( j,i)
          
          ! Update velocities
          U( j,i) = U( j,i) - C%SSA_SOR_omega * resU( j,i)
          V( j,i) = V( j,i) - C%SSA_SOR_omega * resV( j,i)
        
        END DO ! DO j = j0+2, grid%ny-1, 2
        END DO ! DO i = 2, grid%nx-1
      
      END DO ! DO m = 0, 1
      
      ! Apply Neumann boundary conditions
      U( :      ,1      ) = U( :        ,2        )
      U( :      ,grid%nx) = U( :        ,grid%nx-1)
      U( 1      ,:      ) = U( 2        ,:        )
      U( grid%ny,:      ) = U( grid%ny-1,:        )
      V( :      ,1      ) = V( :        ,2        )
      V( :      ,grid%nx) = V( :        ,grid%nx-1)
      V( 1      ,:      ) = V( 2        ,:        )
      V( grid%ny,:      ) = V( grid%ny-1,:        )
      
      ! Check if we've reached a stable solution
      max_residual_UV = MAX( MAXVAL( ABS(resU(2:grid%ny-1,2:grid%nx-1))), MAXVAL( ABS(resV(2:grid%ny-1,2:grid%nx-1))))
      !WRITE(0,*) ' SSA -  inner loop ', inner_loop_i, ': largest residual = ', max_residual_UV
      IF (max_residual_UV < C%SSA_max_residual_UV) THEN
        has_converged = .TRUE.
      ELSEIF (max_residual_UV > 1E8_dp) THEN
        WRITE(0,*) ' ERROR - instability in SSA SOR solver!'
        STOP
      END IF
      
    END DO ! DO WHILE ((.NOT. has_converged) .AND. (inner_loop_i < max_inner_loops))
    
  END SUBROUTINE solve_SSA_linearised_BB_SOR
  SUBROUTINE solve_SSA_linearised_FP_SOR( grid, ice, mask_U_analytical, mask_V_analytical, U, V)
    ! Use successive over-relaxation (SOR) to solve the system of linear equations that
    ! is the SSA when the non-linear terms N and S are assumed constant.
    ! Using the discretisation from Frank Pattyn's f.ETISh model.
    
    USE parameters_module, ONLY: ice_density, grav
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                        INTENT(IN)    :: grid
    TYPE(type_ice_model),                   INTENT(INOUT) :: ice
    INTEGER,  DIMENSION( grid%ny, grid%nx), INTENT(IN   ) :: mask_U_analytical
    INTEGER,  DIMENSION( grid%ny, grid%nx), INTENT(IN   ) :: mask_V_analytical
    REAL(dp), DIMENSION( grid%ny, grid%nx), INTENT(INOUT) :: U, V
    
    ! Local variables:
    INTEGER                                            :: i,j,m,j0
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: LHSx
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: LHSy
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: RHSx
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: RHSy
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: resU
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: resV
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: eu_ij
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: ev_ij
    LOGICAL                                            :: has_converged
    INTEGER                                            :: inner_loop_i
    REAL(dp)                                           :: max_residual_UV
    
    ! Calculate the right-hand sides of the equations (i.e. the gravitational driving stress)
    DO i = 2, grid%nx-1
    DO j = 2, grid%ny-1
      RHSx( j,i) = ice_density * grav * ice%Hi_Aa( j,i) * ice%dHs_dx_shelf_Aa( j,i)
      RHSy( j,i) = ice_density * grav * ice%Hi_Aa( j,i) * ice%dHs_dy_shelf_Aa( j,i)
    END DO
    END DO
    
    ! Calculate the centre coefficients
    DO i = 2, grid%nx-1
    DO j = 2, grid%ny-1
      eu_ij( j,i) = (-8._dp * ice%N_Aa( j,i) / grid%dx**2) + (-2._dp * ice%N_Aa( j,i) / grid%dy**2) - ice%S_Aa( j,i)
      ev_ij( j,i) = (-2._dp * ice%N_Aa( j,i) / grid%dx**2) + (-8._dp * ice%N_Aa( j,i) / grid%dy**2) - ice%S_Aa( j,i)
    END DO
    END DO
    
    ! SOR iteration
    has_converged = .FALSE.
    inner_loop_i  = 0
    DO WHILE ((.NOT. has_converged) .AND. (inner_loop_i < C%SSA_max_inner_loops))
      inner_loop_i = inner_loop_i + 1
      
      ! Red-black scheme
      DO m = 0, 1
      
        j0 = m
        DO i = 2, grid%nx-1
        j0 = 1 - j0
        
        DO j = j0+2, grid%ny-1, 2
        
          ! Keep the values from the analytical flux condition as a constraint
          IF (mask_U_analytical( j,i) == 0) THEN
          
            ! Calculate the left-hand sides of the equations
            LHSx( j,i) = (4._dp * ice%N_Aa(  j,i) * (U( j  ,i+1) + U( j  ,i-1) - 2._dp*U( j  ,i  )) /          grid%dx**2) &
                       + (4._dp * ice%Nx_Aa( j,i) * (U( j  ,i+1) - U( j  ,i-1)                    ) / (2._dp * grid%dx)  ) & 
                       + (        ice%N_Aa(  j,i) * (U( j+1,i  ) + U( j-1,i  ) - 2._dp*U( j  ,i  )) /          grid%dy**2) &
                       + (        ice%Ny_Aa( j,i) * (U( j+1,i  ) - U( j-1,i  )                    ) / (2._dp * grid%dy)  ) &
                       + (3._dp * ice%N_Aa(  j,i) * (V( j+1,i+1) + V( j-1,i-1) - V( j+1,i-1) - V( j-1,i+1)) / (4._dp * grid%dx * grid%dy)) &
                       + (2._dp * ice%Nx_Aa( j,i) * (V( j+1,i  ) - V( j-1,i  )) / (2._dp * grid%dy)) &
                       + (        ice%Ny_Aa( j,i) * (V( j  ,i+1) - V( j  ,i-1)) / (2._dp * grid%dx)) & 
                       - ice%S_Aa( j,i) * U( j,i)
            
            ! Calculate the residuals
            resU( j,i) = (LHSx( j,i) - RHSx( j,i)) / eu_ij( j,i)
            
            ! Update velocities
            U( j,i) = U( j,i) - C%SSA_SOR_omega * resU( j,i)
                       
          END IF
        
          ! Keep the values from the analytical flux condition as a constraint
          IF (mask_V_analytical( j,i) == 0) THEN
            
            LHSy( j,i) = (4._dp * ice%N_Aa(  j,i) * (V( j+1,i  ) + V( j-1,i  ) - 2._dp*V( j  ,i  )) /          grid%dy**2) &
                       + (4._dp * ice%Ny_Aa( j,i) * (V( j+1,i  ) - V( j-1,i  )                    ) / (2._dp * grid%dy)  ) & 
                       + (        ice%N_Aa(  j,i) * (V( j  ,i+1) + V( j  ,i-1) - 2._dp*V( j  ,i  )) /          grid%dx**2) &
                       + (        ice%Nx_Aa( j,i) * (V( j  ,i+1) - V( j  ,i-1)                    ) / (2._dp * grid%dx)  ) &
                       + (3._dp * ice%N_Aa(  j,i) * (U( j+1,i+1) + U( j-1,i-1) - U( j+1,i-1) - U( j-1,i+1)) / (4._dp * grid%dx * grid%dy)) &
                       + (2._dp * ice%Ny_Aa( j,i) * (U( j  ,i+1) - U( j  ,i-1)) / (2._dp * grid%dx)) &
                       + (        ice%Nx_Aa( j,i) * (U( j+1,i  ) - U( j-1,i  )) / (2._dp * grid%dy)) & 
                       - ice%S_Aa( j,i) * V( j,i)
                       
            ! Calculate the residuals
            resV( j,i) = (LHSy( j,i) - RHSy( j,i)) / ev_ij( j,i)
            
            ! Update velocities
            V( j,i) = V( j,i) - C%SSA_SOR_omega * resV( j,i)
                         
          END IF
        
        END DO ! DO j = j0+2, grid%ny-1, 2
        END DO ! DO i = 2, grid%nx-1
      
      END DO ! DO m = 0, 1
      
      ! Apply Neumann boundary conditions
      U( :      ,1      ) = U( :        ,2        )
      U( :      ,grid%nx) = U( :        ,grid%nx-1)
      U( 1      ,:      ) = U( 2        ,:        )
      U( grid%ny,:      ) = U( grid%ny-1,:        )
      V( :      ,1      ) = V( :        ,2        )
      V( :      ,grid%nx) = V( :        ,grid%nx-1)
      V( 1      ,:      ) = V( 2        ,:        )
      V( grid%ny,:      ) = V( grid%ny-1,:        )
      
      ! Check if we've reached a stable solution
      max_residual_UV = MAX( MAXVAL( ABS(resU(2:grid%ny-1,2:grid%nx-1))), MAXVAL( ABS(resV(2:grid%ny-1,2:grid%nx-1))))
      !WRITE(0,*) ' SSA -  inner loop ', inner_loop_i, ': largest residual = ', max_residual_UV
      IF (max_residual_UV < C%SSA_max_residual_UV .AND. (.NOT. (C%do_benchmark_experiment .AND. C%choice_benchmark_experiment == 'SSA_icestream'))) THEN
        has_converged = .TRUE.
      ELSEIF (max_residual_UV > 1E8_dp) THEN
        WRITE(0,*) ' ERROR - instability in SSA SOR solver!'
        STOP
      END IF
      
    END DO ! DO WHILE ((.NOT. has_converged) .AND. (inner_loop_i < max_inner_loops))
    
  END SUBROUTINE solve_SSA_linearised_FP_SOR
  
  SUBROUTINE SSA_effective_viscosity( grid, ice, U, V)
    ! Calculate the effective viscosity eta, the product term N = eta*H, and the gradient of N in the SSA
    
    USE parameters_module, ONLY: n_flow
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION( grid%ny  , grid%nx  ), INTENT(IN) :: U
    REAL(dp), DIMENSION( grid%ny  , grid%nx  ), INTENT(IN) :: V
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION( grid%ny  , grid%nx-1)         :: Ux_Cx
    REAL(dp), DIMENSION( grid%ny  , grid%nx-1)         :: Uy_Cx
    REAL(dp), DIMENSION( grid%ny  , grid%nx-1)         :: Vx_Cx
    REAL(dp), DIMENSION( grid%ny  , grid%nx-1)         :: Vy_Cx
    REAL(dp), DIMENSION( grid%ny-1, grid%nx  )         :: Ux_Cy
    REAL(dp), DIMENSION( grid%ny-1, grid%nx  )         :: Uy_Cy
    REAL(dp), DIMENSION( grid%ny-1, grid%nx  )         :: Vx_Cy
    REAL(dp), DIMENSION( grid%ny-1, grid%nx  )         :: Vy_Cy
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: N_Aa_from_Cx
    REAL(dp), DIMENSION( grid%ny  , grid%nx  )         :: N_Aa_from_Cy
    
    REAL(dp), PARAMETER                                :: epsilon_sq_0 = 1E-12_dp
    
    ! Calculate derivatives of U and V
    CALL ddx_Aa_to_Acx_2D( grid, U, Ux_Cx)
    CALL ddy_Aa_to_Acx_2D( grid, U, Uy_Cx)
    CALL ddx_Aa_to_Acx_2D( grid, V, Vx_Cx)
    CALL ddy_Aa_to_Acx_2D( grid, V, Vy_Cx)
    CALL ddx_Aa_to_Acy_2D( grid, U, Ux_Cy)
    CALL ddy_Aa_to_Acy_2D( grid, U, Uy_Cy)
    CALL ddx_Aa_to_Acy_2D( grid, V, Vx_Cy)
    CALL ddy_Aa_to_Acy_2D( grid, V, Vy_Cy)
    
    ! The effective viscosity eta on both grids, with a normalisation term (following Bueler & Brown, 2009) to prevent divide-by-zero errors,
    ! and the product term N = eta * H
    DO i = 1, grid%nx-1
    DO j = 1, grid%ny
      ice%eta_Acx( j,i) = 0.5_dp * ice%A_flow_mean_Acx( j,i)**(-1._dp/n_flow) * (&
        Ux_Cx( j,i)**2 + Vy_Cx( j,i)**2 + Ux_Cx( j,i)*Vy_Cx( j,i) + 0.25_dp*(Uy_Cx( j,i)+Vx_Cx( j,i))**2 + epsilon_sq_0) ** ((1._dp - n_flow) / (2._dp * n_flow))
      ice%N_Acx( j,i) = ice%eta_Acx( j,i) * MAX( 0.1_dp, ice%Hi_Acx( j,i))
    END DO
    END DO
    DO i = 1, grid%nx
    DO j = 1, grid%ny-1
      ice%eta_Acy( j,i) = 0.5_dp * ice%A_flow_mean_Acy( j,i)**(-1._dp/n_flow) * (&
        Ux_Cy( j,i)**2 + Vy_Cy( j,i)**2 + Ux_Cy( j,i)*Vy_Cy( j,i) + 0.25_dp*(Uy_Cy( j,i)+Vx_Cy( j,i))**2 + epsilon_sq_0) ** ((1._dp - n_flow) / (2._dp * n_flow))
      ice%N_Acy( j,i) = ice%eta_Acy( j,i) * MAX( 0.1_dp, ice%Hi_Acy( j,i))
    END DO
    END DO
    
    ! Map N to the A grid and get its derivatives there
    CALL map_Acx_to_Aa_2D( grid, ice%N_Acx, N_Aa_from_Cx)
    CALL map_Acy_to_Aa_2D( grid, ice%N_Acy, N_Aa_from_Cy)
    ice%N_Aa = 0.5_dp * (N_Aa_from_Cx + N_Aa_from_Cy)
    
    CALL ddx_Acx_to_Aa_2D( grid, ice%N_Acx, ice%Nx_Aa)
    CALL ddy_Acy_to_Aa_2D( grid, ice%N_Acy, ice%Ny_Aa)
    
    ! Crop the derivatives of N (tip from Frank Pattyn, personal communication, 2021)
    ice%Nx_Aa = MAX( -C%SSA_max_grad_N, MIN( C%SSA_max_grad_N, ice%Nx_Aa))
    ice%Ny_Aa = MAX( -C%SSA_max_grad_N, MIN( C%SSA_max_grad_N, ice%Ny_Aa))
    
  END SUBROUTINE SSA_effective_viscosity
  SUBROUTINE SSA_sliding_term( grid, ice, U, V)
    ! Calculate the sliding term S = tau_c / norm([U,V]) in the SSA
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION( grid%ny  , grid%nx  ), INTENT(IN) :: U
    REAL(dp), DIMENSION( grid%ny  , grid%nx  ), INTENT(IN) :: V
    
    REAL(dp), PARAMETER                                :: delta_v = 1E-3_dp
    REAL(dp), PARAMETER                                :: q_plastic            = 0.30_dp       ! Parameter used for basal stress (inverse of m_flow)
    REAL(dp), PARAMETER                                :: u_threshold          = 100._dp       ! scaling of tau_yield to get the correct unit (function of q_plastic)
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    ice%S_Aa = 0._dp
    
    DO i = 1, grid%nx
    DO j = 1, grid%ny
      IF (ice%mask_sheet_Aa( j,i) == 1 .OR. ice%mask_land_Aa( j,i) == 1) THEN    
        ! The sliding term S, with a normalisation term (following Bueler & Brown, 2009) to prevent divide-by-zero errors.
        !ice%S_Aa( j,i) = ice%tau_c_Aa( j,i) / SQRT( delta_v + U( j,i)**2 + V( j,i)**2)
        ice%S_Aa( j,i) = ice%tau_c_Aa( j,i) * ( (delta_v**2 + U( j,i)**2 + V( j,i)**2)**(0.5_dp * (q_plastic-1._dp)) ) / (u_threshold**q_plastic)
        
        !ice%S_Aa( j,i) = ice%S_Aa( j,i) * ice%grounded_fraction( j,i)
      END IF
    END DO
    END DO
    
  END SUBROUTINE SSA_sliding_term
  
  ! Administration: allocation and initialisation
  SUBROUTINE initialise_ice_model( grid, ice, init)
    ! Allocate shared memory for all the data fields of the ice dynamical module, and
    ! initialise some of them  
    
    USE parameters_module, ONLY: n_flow
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_init_data_fields),         INTENT(IN)    :: init
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    ! Allocate shared memory
    CALL allocate_ice_model( grid, ice)
        
    ! Initialise with data from initial file    
    ice%Hi_Aa = init%Hi
    ice%Hb_Aa = init%Hb
    ice%Hs_Aa = MAX( ice%SL_Aa, ice%Hb_Aa + ice%Hi_Aa)
    
    ice%U_vav_SSA_Aa  = 0._dp
    ice%V_vav_SSA_Aa  = 0._dp
    ice%U_vav_SSA_Acx = 0._dp
    ice%V_vav_SSA_Acy = 0._dp
    
  END SUBROUTINE initialise_ice_model
  SUBROUTINE allocate_ice_model( grid, ice)  
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
      
    ! Allocate memory
    ! ===============
    
    ! Basic data - ice thickness, bedrock & surface elevation, sea level (geoid elevation), englacial temperature, and ice velocities
    ALLOCATE( ice%Hi_Aa(                 grid%ny  , grid%nx   ))
    ALLOCATE( ice%Hi_Acx(                grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%Hi_Acy(                grid%ny-1, grid%nx   ))
    ALLOCATE( ice%Hi_Ab(                 grid%ny-1, grid%nx-1 ))
    ALLOCATE( ice%Hb_Aa(                 grid%ny  , grid%nx   ))
    ALLOCATE( ice%Hb_Acx(                grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%Hb_Acy(                grid%ny-1, grid%nx   ))
    ALLOCATE( ice%Hb_Ab(                 grid%ny-1, grid%nx-1 ))
    ALLOCATE( ice%Hs_Aa(                 grid%ny  , grid%nx   ))
    ALLOCATE( ice%Hs_Acx(                grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%Hs_Acy(                grid%ny-1, grid%nx   ))
    ALLOCATE( ice%Hs_Ab(                 grid%ny-1, grid%nx-1 ))
    ALLOCATE( ice%SL_Aa(                 grid%ny  , grid%nx   ))
    ALLOCATE( ice%SL_Acx(                grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%SL_Acy(                grid%ny-1, grid%nx   ))
    ALLOCATE( ice%SL_Ab(                 grid%ny-1, grid%nx-1 ))
    ALLOCATE( ice%Ti_Aa(           C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%Ti_Acx(          C%nZ, grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%Ti_Acy(          C%nZ, grid%ny-1, grid%nx   ))
    
    ALLOCATE( ice%U_vav_SIA_Aa(          grid%ny  , grid%nx   ))
    ALLOCATE( ice%V_vav_SIA_Aa(          grid%ny  , grid%nx   ))
    ALLOCATE( ice%U_vav_SIA_Acx(         grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%V_vav_SIA_Acy(         grid%ny-1, grid%nx   ))
    ALLOCATE( ice%U_vav_SSA_Aa(          grid%ny  , grid%nx   ))
    ALLOCATE( ice%V_vav_SSA_Aa(          grid%ny  , grid%nx   ))
    ALLOCATE( ice%U_vav_SSA_Acx(         grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%V_vav_SSA_Acy(         grid%ny-1, grid%nx   ))
    ALLOCATE( ice%U_vav_Aa(              grid%ny  , grid%nx   ))
    ALLOCATE( ice%V_vav_Aa(              grid%ny  , grid%nx   ))
    ALLOCATE( ice%U_vav_Acx(             grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%V_vav_Acy(             grid%ny-1, grid%nx   ))
    ALLOCATE( ice%U_surf_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%V_surf_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%U_base_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%V_base_Aa(             grid%ny  , grid%nx   ))
    
    ! Different masks
    ALLOCATE( ice%mask_land_Aa(          grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_ocean_Aa(         grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_ocean_Acx(        grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%mask_ocean_Acy(        grid%ny-1, grid%nx   ))
    ALLOCATE( ice%mask_lake_Aa(          grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_ice_Aa(           grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_sheet_Aa(         grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_sheet_Acx(        grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%mask_sheet_Acy(        grid%ny-1, grid%nx   ))
    ALLOCATE( ice%mask_shelf_Aa(         grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_shelf_Acx(        grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%mask_shelf_Acy(        grid%ny-1, grid%nx   ))
    ALLOCATE( ice%mask_coast_Aa(         grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_coast_Acx(        grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%mask_coast_Acy(        grid%ny-1, grid%nx   ))
    ALLOCATE( ice%mask_margin_Aa(        grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_margin_Acx(       grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%mask_margin_Acy(       grid%ny-1, grid%nx   ))
    ALLOCATE( ice%mask_gl_Aa(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_gl_Acx(           grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%mask_gl_Acy(           grid%ny-1, grid%nx   ))
    ALLOCATE( ice%mask_cf_Aa(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_cf_Acx(           grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%mask_cf_Acy(           grid%ny-1, grid%nx   ))
    ALLOCATE( ice%mask_Aa(               grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_Acx(              grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%mask_Acy(              grid%ny-1, grid%nx   ))
    ALLOCATE( ice%grounded_fraction(     grid%ny  , grid%nx   ))
    
    ! Ice physical properties
    ALLOCATE( ice%A_flow_Aa(       C%nZ, grid%ny  , grid%nx   ))     ! Flow parameter [Pa^-3 y^-1]
    ALLOCATE( ice%A_flow_Acx(      C%nZ, grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%A_flow_Acy(      C%nZ, grid%ny-1, grid%nx   ))
    ALLOCATE( ice%A_flow_mean_Aa(        grid%ny  , grid%nx   ))     ! Vertically averaged flow parameter [Pa^-3 y^-1]
    ALLOCATE( ice%A_flow_mean_Acx(       grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%A_flow_mean_Acy(       grid%ny-1, grid%nx   ))
    ALLOCATE( ice%A_flow_mean_Ab(        grid%ny-1, grid%nx-1 ))
    ALLOCATE( ice%Ti_pmp_Aa(       C%nZ, grid%ny  , grid%nx   ))     ! The pressure melting point temperature [K]
    ALLOCATE( ice%Cpi_Aa(          C%nZ, grid%ny  , grid%nx   ))     ! Specific heat capacity of ice [J kg^-1 K^-1]
    ALLOCATE( ice%Ki_Aa(           C%nZ, grid%ny  , grid%nx   ))     ! Conductivity of ice [J m^-1 K^-1 yr^-1]
    
    ! Spatial and temporal derivatives
    ALLOCATE( ice%dHi_dt_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%dHi_dx_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%dHi_dy_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%dHi_dx_Acx(            grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%dHi_dy_Acx(            grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%dHi_dx_Acy(            grid%ny-1, grid%nx   ))
    ALLOCATE( ice%dHi_dy_Acy(            grid%ny-1, grid%nx   ))
    ALLOCATE( ice%dHb_dt_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%dHb_dx_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%dHb_dy_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%dHb_dx_Acx(            grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%dHb_dy_Acx(            grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%dHb_dx_Acy(            grid%ny-1, grid%nx   ))
    ALLOCATE( ice%dHb_dy_Acy(            grid%ny-1, grid%nx   ))
    ALLOCATE( ice%dHs_dt_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%dHs_dx_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%dHs_dy_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%dHs_dx_Acx(            grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%dHs_dy_Acx(            grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%dHs_dx_Acy(            grid%ny-1, grid%nx   ))
    ALLOCATE( ice%dHs_dy_Acy(            grid%ny-1, grid%nx   ))
    ALLOCATE( ice%dSL_dx_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%dSL_dy_Aa(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%dSL_dx_Acx(            grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%dSL_dy_Acx(            grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%dSL_dx_Acy(            grid%ny-1, grid%nx   ))
    ALLOCATE( ice%dSL_dy_Acy(            grid%ny-1, grid%nx   ))
    ALLOCATE( ice%dTi_dx_Aa(       C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dTi_dy_Aa(       C%nZ, grid%ny  , grid%nx   ))
    
    ! Zeta derivatives
    ALLOCATE( ice%dzeta_dt_Aa(     C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dzeta_dx_Aa(     C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dzeta_dy_Aa(     C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dzeta_dz_Aa(           grid%ny  , grid%nx   ))
    ALLOCATE( ice%dzeta_dxx_Aa (   C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dzeta_dxy_Aa (   C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dzeta_dyy_Aa (   C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dzeta_dxz_Aa (         grid%ny  , grid%nx   ))
    ALLOCATE( ice%dzeta_dyz_Aa (         grid%ny  , grid%nx   ))
        
    ! Ice dynamics - SIA
    ALLOCATE( ice%D_SIA_3D_Acx(    C%nZ, grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%D_SIA_3D_Acy(    C%nZ, grid%ny-1, grid%nx   ))
    ALLOCATE( ice%D_SIA_Acx(             grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%D_SIA_Acy(             grid%ny-1, grid%nx   ))
    ALLOCATE( ice%D_SIA_Aa(              grid%ny  , grid%nx   ))
    
    ! Ice dynamics - SSA
    ALLOCATE( ice%phi_fric_Aa(           grid%ny  , grid%nx   ))
    ALLOCATE( ice%phi_fric_Ab(           grid%ny-1, grid%nx-1 ))
    ALLOCATE( ice%tau_c_Aa(              grid%ny  , grid%nx   ))
    ALLOCATE( ice%tau_c_Ab(              grid%ny-1, grid%nx-1 ))
    ALLOCATE( ice%eta_Acx(               grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%eta_Acy(               grid%ny-1, grid%nx   ))
    ALLOCATE( ice%eta_Ab(                grid%ny-1, grid%nx-1 ))
    ALLOCATE( ice%N_Acx(                 grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%N_Acy(                 grid%ny-1, grid%nx   ))
    ALLOCATE( ice%N_Aa(                  grid%ny  , grid%nx   ))
    ALLOCATE( ice%N_Aa_prev(             grid%ny  , grid%nx   ))
    ALLOCATE( ice%N_Ab(                  grid%ny-1, grid%nx-1 ))
    ALLOCATE( ice%N_Ab_prev(             grid%ny-1, grid%nx-1 ))
    ALLOCATE( ice%Nx_Aa(                 grid%ny  , grid%nx   ))
    ALLOCATE( ice%Ny_Aa(                 grid%ny  , grid%nx   ))
    ALLOCATE( ice%S_Aa(                  grid%ny  , grid%nx   ))
    ALLOCATE( ice%S_Acx(                 grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%S_Acy(                 grid%ny-1, grid%nx   ))
    ALLOCATE( ice%S_Ab(                  grid%ny-1, grid%nx-1 ))
    ALLOCATE( ice%dHs_dx_shelf_Aa(       grid%ny  , grid%nx   ))
    ALLOCATE( ice%dHs_dy_shelf_Aa(       grid%ny  , grid%nx   ))
    ALLOCATE( ice%U_SSA_free(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%V_SSA_free(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%U_SSA_cons(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%V_SSA_cons(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_U_analytical(     grid%ny  , grid%nx   ))
    ALLOCATE( ice%mask_V_analytical(     grid%ny  , grid%nx   ))
    
    ! Ice dynamics- SSA - semi-analytical GL flux
    ALLOCATE( ice%Qx_GL_Acx(             grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%Qy_GL_Acy(             grid%ny-1, grid%nx   ))
    
    ! Ice dynamics - ice thickness calculation
    ALLOCATE( ice%Qx_Acx(                grid%ny  , grid%nx-1 ))
    ALLOCATE( ice%Qy_Acy(                grid%ny-1, grid%nx   ))
    
    ! Thermodynamics
    ALLOCATE( ice%frictional_heating_Aa( grid%ny  , grid%nx   ))
    ALLOCATE( ice%Fr_Aa(                 grid%ny  , grid%nx   ))
    ALLOCATE( ice%U_3D_Aa(         C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%V_3D_Aa(         C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%W_3D_Aa(         C%nZ, grid%ny  , grid%nx   ))
        
    ! Dummy variables for debugging
    ALLOCATE( ice%dummy2D_01(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy2D_02(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy2D_03(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy2D_04(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy2D_05(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy2D_06(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy2D_07(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy2D_08(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy2D_09(            grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy2D_10(            grid%ny  , grid%nx   ))
    
    ALLOCATE( ice%dummy3D_01(      C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy3D_02(      C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy3D_03(      C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy3D_04(      C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy3D_05(      C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy3D_06(      C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy3D_07(      C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy3D_08(      C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy3D_09(      C%nZ, grid%ny  , grid%nx   ))
    ALLOCATE( ice%dummy3D_10(      C%nZ, grid%ny  , grid%nx   ))
    
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
