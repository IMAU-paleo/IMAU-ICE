MODULE thermodynamics_module

  ! Contains all the routines for calculating the englacial temperature field.

  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE netcdf_module,                   ONLY: debug, write_to_debug_file 
  USE parameters_module
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_subclimate_region, type_SMB_model, type_init_data_fields
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             tridiagonal_solve
  USE derivatives_and_grids_module,    ONLY: zeta, calculate_zeta_derivatives, ddx_a_to_a_2D, ddy_a_to_a_2D, Neumann_BC_a_3D
  USE general_ice_model_data_module,   ONLY: ice_physical_properties
  
  IMPLICIT NONE
  
CONTAINS
   
  SUBROUTINE run_thermo_model( grid, ice, climate, SMB, do_solve_heat_equation)
    ! Run the thermodynamics model. If so specified, solve the heat equation;
    ! if not, only prescribe a simple temperature profile to newly ice-covered grid cells.
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_subclimate_region),        INTENT(IN)    :: climate
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    LOGICAL,                             INTENT(IN)    :: do_solve_heat_equation

    ! Local variables:
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: T_surf_annual
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6') THEN
        ! Thermodynamics are included in these experiments
      ELSEIF (C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'Bueler' .OR. &
              C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
              C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_F' .OR. &
              C%choice_benchmark_experiment == 'MISMIPplus' .OR. &
              C%choice_benchmark_experiment == 'MISOMIP1') THEN
        ! Thermodynamics are not included in these experiments
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_thermo_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Prescribe a simple temperature profile to newly ice-covered grid cells.
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 1, grid%ny
      
      IF (ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a_prev( j,i) == 0) THEN
        ! This grid cell is newly ice-covered
        ! If one of its neighbours was already ice-covered, assume the temperature
        ! profile here is equal to the profile from the upstream neighbour (due to advection).
        ! If no neighbours were ice-covered, the new ice must come from accumulation;
        ! just set a simple linear profile instead.
        
        IF     (ice%u_vav_cx( j  ,i-1) > 0._dp .AND. ice%mask_ice_a_prev( j  ,i-1) == 1) THEN
          ! Ice probably came from the west
          ice%Ti_a( :,j,i) = ice%Ti_a( :,j  ,i-1)
        ELSEIF (ice%u_vav_cx( j  ,i  ) < 0._dp .AND. ice%mask_ice_a_prev( j  ,i+1) == 1) THEN
          ! Ice probably came from the east
          ice%Ti_a( :,j,i) = ice%Ti_a( :,j  ,i+1)
        ELSEIF (ice%v_vav_cy( j-1,i  ) > 0._dp .AND. ice%mask_ice_a_prev( j-1,i  ) == 1) THEN
          ! Ice probably came from the south
          ice%Ti_a( :,j,i) = ice%Ti_a( :,j-1,i  )
        ELSEIF (ice%v_vav_cy( j  ,i  ) < 0._dp .AND. ice%mask_ice_a_prev( j+1,i  ) == 1) THEN
          ! Ice probably came from the north
          ice%Ti_a( :,j,i) = ice%Ti_a( :,j+1,i  )
        ELSE
          ! Ice probably came from surface accumulation; initialise with a linear temperature profile.
          
          T_surf_annual = MIN( SUM( climate%T2m( :,j,i)) / 12._dp, T0)
          DO k = 1, C%nz
            ice%Ti_a( k,j,i) = T_surf_annual + C%zeta(k) * (ice%Ti_pmp_a( C%nz,j,i) - T_surf_annual)
          END DO
          
        END IF
        
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! If so specified, solve the heat equation
    IF (do_solve_heat_equation) CALL solve_heat_equation( grid, ice, climate, SMB)
    
    ! Safety
    CALL check_for_NaN_dp_3D( ice%Ti_a, 'ice%Ti_a', 'run_thermo_model')

  END SUBROUTINE run_thermo_model
  SUBROUTINE solve_heat_equation( grid, ice, climate, SMB)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_subclimate_region),        INTENT(IN)    :: climate
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB

    ! Local variables:
    INTEGER                                            :: i,j,k
    REAL(dp), DIMENSION(:,:,:), POINTER                :: Ti_new
    INTEGER                                            :: wTi_new
    REAL(dp)                                           :: internal_heating, u_times_dT_dx_upwind, v_times_dT_dy_upwind, f1, f2, f3
    REAL(dp), DIMENSION(2:C%nz)                        :: alpha
    REAL(dp), DIMENSION(C%nz)                          :: beta
    REAL(dp), DIMENSION(C%nz-1)                        :: gamma
    REAL(dp), DIMENSION(C%nz)                          :: delta
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  is_unstable
    INTEGER                                            :: wis_unstable
    LOGICAL                                            :: hasnan
    INTEGER                                            :: n_unstable
    
    ! Allocate shared memory
    CALL allocate_shared_int_2D(       grid%ny, grid%nx, is_unstable, wis_unstable)
    CALL allocate_shared_dp_3D(  C%nz, grid%ny, grid%nx, Ti_new,      wTi_new     )

    ! Calculate the required derivatives of Hi and Hs   
    CALL ddx_a_to_a_2D( grid, ice%Hi_a, ice%dHi_dx_a)
    CALL ddy_a_to_a_2D( grid, ice%Hi_a, ice%dHi_dy_a)
    CALL ddx_a_to_a_2D( grid, ice%Hs_a, ice%dHs_dx_a)
    CALL ddy_a_to_a_2D( grid, ice%Hs_a, ice%dHs_dy_a)
    
    ! (dHi_dt and dHb_dt are set in other routines)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%mask_shelf_a( j,i) == 0) THEN
        ice%dHs_dt_a( j,i) = ice%dHi_dt_a( j,i) + ice%dHb_dt_a( j,i)
      ELSE
        ice%dHs_dt_a( j,i) = ice%dHi_dt_a( j,i) * (1._dp - ice_density / seawater_density)
      END IF
    END DO
    END DO
    CALL sync
    
    ! Calculate zeta derivatives required for solving the heat equation
    CALL calculate_zeta_derivatives( grid, ice)
    
    ! Set ice surface temperature equal to annual mean 2m air temperature
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%Ti_a( 1,j,i) = MIN( T0, SUM( climate%T2m( :,j,i)) / 12._dp)
    END DO
    END DO
    CALL sync
    
    ! Calculate bottom frictional heating
    CALL bottom_frictional_heating(  grid, ice)
    
    ! Solve the heat equation for all interior grid cells
    is_unstable( :,grid%i1:grid%i2) = 0
    n_unstable                      = 0
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    
      ! Exceptions
      ! ==========
      
      ! Skip ice-less elements
      IF (ice%mask_ice_a( j,i) == 0) THEN
         Ti_new( :,j,i) = ice%Ti_a( 1,j,i)
         CYCLE
      END IF
      
      ! For very thin ice, just let the profile equal the surface temperature
      IF (ice%Hi_a( j,i) < 1._dp) THEN
        Ti_new( :,j,i) = MIN( T0, SUM( climate%T2m( :,j,i)) / 12._dp)
        CYCLE
      END IF
      
      ! Fill in the tridiagonal matrix equation
      ! ========================================
  
      DO k = 2, C%nz-1
        
        IF (ice%mask_sheet_a( j,i) == 1) THEN
          internal_heating = ((- grav * C%zeta(k)) / ice%Cpi_a( k,j,i)) * ( &
               (zeta%a_zeta(k) * ice%U_3D_a( k-1,j,i) + zeta%b_zeta(k) * ice%U_3D_a( k,j,i) + zeta%c_zeta(k) * ice%U_3D_a( k+1,j,i)) * ice%dHs_dx_a( j,i) + &
               (zeta%a_zeta(k) * ice%V_3D_a( k-1,j,i) + zeta%b_zeta(k) * ice%V_3D_a( k,j,i) + zeta%c_zeta(k) * ice%V_3D_a( k+1,j,i)) * ice%dHs_dy_a( j,i) )
        ELSE
          internal_heating = 0._dp
        END IF
        
        IF (ice%u_3D_a( k,j,i) > 0._dp) THEN
          u_times_dT_dx_upwind = ice%u_3D_cx( k,j,i-1) * (ice%Ti_a( k,j,i  ) - ice%Ti_a( k,j,i-1)) / grid%dx
        ELSE
          u_times_dT_dx_upwind = ice%u_3D_cx( k,j,i  ) * (ice%Ti_a( k,j,i+1) - ice%Ti_a( k,j,i  )) / grid%dx
        END IF
        IF (ice%v_3D_a( k,j,i) > 0._dp) THEN
          v_times_dT_dy_upwind = ice%v_3D_cy( k,j-1,i) * (ice%Ti_a( k,j  ,i) - ice%Ti_a( k,j-1,i)) / grid%dx
        ELSE
          v_times_dT_dy_upwind = ice%v_3D_cy( k,j,i  ) * (ice%Ti_a( k,j+1,i) - ice%Ti_a( k,j  ,i)) / grid%dx
        END IF

        f1 = (ice%Ki_a( k,j,i) * ice%dzeta_dz_a( j,i)**2) / (ice_density * ice%Cpi_a( k,j,i))

        f2 = ice%dzeta_dt_a( k,j,i) + ice%dzeta_dx_a( k,j,i) * ice%U_3D_a( k,j,i) + ice%dzeta_dy_a( k,j,i) * ice%V_3D_a( k,j,i) + ice%dzeta_dz_a( j,i) * ice%W_3D_a( k,j,i)

        f3 = internal_heating + (u_times_dT_dx_upwind + v_times_dT_dy_upwind) - ice%Ti_a( k,j,i) / C%dt_thermo

        alpha(k) = f1 * zeta%a_zetazeta(k) - f2 * zeta%a_zeta(k)
        beta (k) = f1 * zeta%b_zetazeta(k) - f2 * zeta%b_zeta(k) - 1._dp / C%dt_thermo
        gamma(k) = f1 * zeta%c_zetazeta(k) - f2 * zeta%c_zeta(k)
        delta(k) = f3
        
!        ! DENK DROM - only vertical diffusion
!        alpha(k) = (zeta%a_zeta(k) * ice%dzeta_dt(vi,k)) - ((zeta%a_zetazeta(k) * ice%Ki(vi,k)) / (ice_density * ice%Cpi(vi,k) * ice%Hi(vi)**2))
!        beta( k) = (zeta%b_zeta(k) * ice%dzeta_dt(vi,k)) - ((zeta%b_zetazeta(k) * ice%Ki(vi,k)) / (ice_density * ice%Cpi(vi,k) * ice%Hi(vi)**2)) + (1._dp / C%dt_thermo)
!        gamma(k) = (zeta%c_zeta(k) * ice%dzeta_dt(vi,k)) - ((zeta%c_zetazeta(k) * ice%Ki(vi,k)) / (ice_density * ice%Cpi(vi,k) * ice%Hi(vi)**2))
!        delta(k) = ice%Ti(vi,k) / C%dt_thermo
        
!        ! DENK DROM - only vertical diffusion + vertical advection
!        alpha(k) = (zeta%a_zeta(k) * (ice%dzeta_dt(vi,k) - ice%W_3D(vi,k) / ice%Hi(vi))) - ((zeta%a_zetazeta(k) * ice%Ki(vi,k)) / (ice_density * ice%Cpi(vi,k) * ice%Hi(vi)**2))
!        beta( k) = (zeta%b_zeta(k) * (ice%dzeta_dt(vi,k) - ice%W_3D(vi,k) / ice%Hi(vi))) - ((zeta%b_zetazeta(k) * ice%Ki(vi,k)) / (ice_density * ice%Cpi(vi,k) * ice%Hi(vi)**2)) + (1._dp / C%dt_thermo)
!        gamma(k) = (zeta%c_zeta(k) * (ice%dzeta_dt(vi,k) - ice%W_3D(vi,k) / ice%Hi(vi))) - ((zeta%c_zetazeta(k) * ice%Ki(vi,k)) / (ice_density * ice%Cpi(vi,k) * ice%Hi(vi)**2))
!        delta(k) = ice%Ti(vi,k) / C%dt_thermo

      END DO ! DO k = 2, C%nz-1
      
      ! Boundary conditions at the surface: set ice temperature equal to annual mean surface temperature
      beta(  1) = 1._dp
      gamma( 1) = 0._dp
      delta( 1) =  MIN( T0, SUM( climate%T2m( :,j,i)) / 12._dp)
      
      ! Boundary conditions at the base
      IF (ice%mask_shelf_a( j,i) == 1) THEN
        ! Set ice bottom temperature equal to seawater temperature (limited to the PMP)
        alpha( C%nz) = 0._dp
        beta ( C%nz) = 1._dp
        delta( C%nz) = MIN( T0, MIN( ice%Ti_pmp_a( C%nz,j,i), T0 + climate%T_ocean_mean ))
      ELSE
        IF (ice%Ti_a( C%nz,j,i) >= ice%Ti_pmp_a( C%nz,j,i)) THEN
          ! Ice is already at/above pressure melting point; set temperature equal to PMP
          alpha( C%nz) = 0._dp
          beta ( C%nz) = 1._dp
          delta( C%nz) = ice%Ti_pmp_a( C%nz,j,i)
        ELSE
          ! Set a Neumann BC so the temperature gradient at the base is equal to basal heating rate (= geothermal + friction)
          alpha( C%nz) = 1._dp
          beta ( C%nz) = -1._dp
          delta( C%nz) = (C%zeta(C%nz) - C%zeta(C%nz-1)) * (ice%GHF_a( j,i) + ice%frictional_heating_a( j,i)) / (ice%dzeta_dz_a( j,i) * ice%Ki_a( C%nz,j,i)) 
        END IF
      END IF ! IF (ice%mask_shelf(vi) == 1 .OR. ice%mask_groundingline(vi) == 1) THEN

      ! Solve the tridiagonal matrix equation representing the heat equation for this grid cell
      Ti_new( :,j,i) = tridiagonal_solve( alpha, beta, gamma, delta)
      
      ! Make sure ice temperature doesn't exceed pressure melting point
      DO k = 1, C%nz-1
        Ti_new( k,j,i) = MIN( Ti_new( k,j,i), ice%Ti_pmp_a( k,j,i))
      END DO
      
      IF (Ti_new( C%nz,j,i) >= ice%Ti_pmp_a( C%nz,j,i)) THEN
        Ti_new( C%nz,j,i) = MIN( ice%Ti_pmp_a( C%nz,j,i), ice%Ti_a( C%nz-1,j,i) - (C%zeta(C%nz) - C%zeta(C%nz-1)) * &
          (ice%GHF_a( j,i) + ice%frictional_heating_a( j,i)) / (ice%dzeta_dz_a( j,i) * ice%Ki_a( C%nz,j,i)))
      END IF
    
      ! Mark temperatures below 150 K or NaN as unstable, to be replaced with the Robin solution.
      hasnan = .FALSE.
      DO k = 1, C%nz
        IF (Ti_new( k,j,i) /= Ti_new( k,j,i)) THEN
          hasnan = .TRUE.
        END IF
      END DO
      IF (MINVAL(Ti_new( :,j,i)) < 150._dp .OR. hasnan) THEN
        is_unstable( j,i) = 1
        n_unstable  = n_unstable + 1
        !WRITE(0,*) 'instability detected; Hi = ', ice%Hi_a( j,i), ', dHi_dt = ', ice%dHi_dt_a( j,i)
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Apply Neumann boundary conditions to the temperature field
    CALL Neumann_BC_a_3D( grid, Ti_new)
    
    ! Cope with instability
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_unstable, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)    
    IF (n_unstable < CEILING(REAL(grid%nx*grid%ny) / 100._dp)) THEN
      ! Instability is limited to an acceptably small number (< 1%) of grid cells;
      ! replace the temperature profile in those cells with the Robin solution
      
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (is_unstable( j,i) == 1) CALL replace_Ti_with_robin_solution( ice, climate, SMB, Ti_new, i,j)
      END DO
      END DO
      CALL sync
      
    ELSE
      ! An unacceptably large number of grid cells was unstable; throw an error.
      
      IF (par%master) WRITE(0,*) '   solve_heat_equation - ERROR:  heat equation solver unstable for more than 1% of grid cells!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      
    END IF
    
    ! Move the new temperature field to the ice data structure
    ice%Ti_a( :,:,grid%i1:grid%i2) = Ti_new( :,:,grid%i1:grid%i2)
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wTi_new)
    CALL deallocate_shared( wis_unstable)
    
    ! Safety
    CALL check_for_NaN_dp_3D( ice%Ti_a, 'ice%Ti_a', 'solve_heat_equation')

  END SUBROUTINE solve_heat_equation
  
  SUBROUTINE replace_Ti_with_robin_solution( ice, climate, SMB, Ti, i,j)
    ! This function calculates for one horizontal grid point the temperature profiles
    ! using the surface temperature and the geothermal heat flux as boundary conditions.
    ! See Robin solution in: Cuffey & Paterson 2010, 4th ed, chapter 9, eq. (9.13) - (9.22).
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_subclimate_region),        INTENT(IN)    :: climate
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: Ti
    INTEGER,                             INTENT(IN)    :: i,j

    ! Local variables:
    INTEGER                                            :: k
    REAL(dp)                                           :: Ts
    REAL(dp)                                           :: thermal_length_scale
    REAL(dp)                                           :: distance_above_bed
    REAL(dp)                                           :: erf1
    REAL(dp)                                           :: erf2
    
    REAL(dp)                                           :: thermal_conductivity_robin
    REAL(dp)                                           :: thermal_diffusivity_robin
    REAL(dp)                                           :: bottom_temperature_gradient_robin
    
    REAL(dp), PARAMETER                                :: kappa_0_ice_conductivity     = 9.828_dp                      ! The linear constant in the thermal conductivity of ice [J m^-1 K^-1 s^-1], see equation (12.6), Ritz (1987), Cuffey & Paterson (2010, p. 400), Zwinger (2007)
    REAL(dp), PARAMETER                                :: kappa_e_ice_conductivity     = 0.0057_dp                     ! The exponent constant in the thermal conductivity of ice [K^-1], see equation (12.6), Ritz (1987), Cuffey & Paterson (2010, p. 400), Zwinger (2007)
    REAL(dp), PARAMETER                                :: c_0_specific_heat            = 2127.5_dp                     ! The constant in the specific heat capacity of ice [J kg^-1 K^-1], see equation (12.5), Zwinger (2007), Cuffey & Paterson (2010, p. 400)
    
    thermal_conductivity_robin        = kappa_0_ice_conductivity * sec_per_year * EXP(-kappa_e_ice_conductivity * T0)  ! Thermal conductivity            [J m^-1 K^-1 y^-1]
    thermal_diffusivity_robin         = thermal_conductivity_robin / (ice_density * c_0_specific_heat)                 ! Thermal diffusivity             [m^2 y^-1]
    bottom_temperature_gradient_robin = - ice%GHF_a( j,i) / thermal_conductivity_robin                                ! Temperature gradient at bedrock
    
    Ts = MIN( T0, SUM(climate%T2m( :,j,i)) / 12._dp)
    
    IF (ice%mask_sheet_a( j,i) == 1 ) THEN
    
      IF (SMB%SMB_year( j,i) > 0._dp) THEN    
        ! The Robin solution can be used to estimate the subsurface temperature profile in an accumulation area
        
        thermal_length_scale = SQRT(2._dp * thermal_diffusivity_robin * ice%Hi_a( j,i) / SMB%SMB_year( j,i))
        DO k = 1, C%nz
          distance_above_bed = (1._dp - C%zeta(k)) * ice%Hi_a( j,i)
          erf1 = erf( distance_above_bed / thermal_length_scale)
          erf2 = erf( ice%Hi_a( j,i) / thermal_length_scale)
          Ti( k,j,i) = Ts + SQRT(pi) / 2._dp * thermal_length_scale * bottom_temperature_gradient_robin * (erf1 - erf2)
        END DO
      
      ELSE
    
        ! Ablation area: use linear temperature profile from Ts to (offset below) T_pmp
        Ti( :,j,i) = Ts + ((T0 - CC * ice%Hi_a( j,i)) - Ts) * C%zeta(:)
      
      END IF
      
    ELSEIF( ice%mask_shelf_a( j,i) == 1) THEN
    
      ! Use a linear profile between Ts and seawater temperature:
      Ti( :,j,i) = Ts + C%zeta(:) * (T0 + climate%T_ocean_mean - Ts)
      
    ELSE
    
      ! No ice present: use Ts everywhere
      Ti( :,j,i) = Ts
      
    END IF

    ! Correct all temperatures above T_pmp:
    DO k = 1, C%nz
      Ti( k,j,i) = MIN( Ti( k,j,i), T0 - CC * ice%Hi_a( j,i) * C%zeta(k))
    END DO

  END SUBROUTINE replace_Ti_with_robin_solution
  SUBROUTINE bottom_frictional_heating( grid, ice)
    ! Calculation of the frictional heating at the bottom due to sliding at the sheet/Gl - bedrock interface.
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_grid),                     INTENT(IN)    :: grid

    ! Local variables
    INTEGER                                            :: i,j
    REAL(dp)                                           :: beta
    
    REAL(dp), PARAMETER                                :: delta_v = 1E-4_dp
    REAL(dp), PARAMETER                                :: q_plastic            = 0.30_dp       ! Parameter used for basal stress (inverse of m_flow)
    REAL(dp), PARAMETER                                :: u_threshold          = 100._dp       ! scaling of tau_yield to get the correct unit (function of q_plastic)
    
    ! Exception for when no sliding can pccur
    IF (C%choice_ice_dynamics == 'SIA' .OR. C%no_sliding) THEN
      ice%frictional_heating_a( :,grid%i1:grid%i2) = 0._dp
      RETURN
    END IF
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    
      ice%frictional_heating_a( j,i) = 0._dp
    
      IF (ice%mask_sheet_a( j,i) == 1) THEN
      
        beta = ice%tauc_a( j,i) * ( (delta_v**2 + ice%U_base_a( j,i)**2 + ice%V_base_a( j,i)**2)**(0.5_dp * (q_plastic-1._dp)) ) / (u_threshold**q_plastic)
        ice%frictional_heating_a( j,i) = beta * (ice%U_base_a( j,i)**2 + ice%V_base_a( j,i)**2)
        
      END IF 
           
    END DO
    END DO
    CALL sync

  END SUBROUTINE bottom_frictional_heating
  
  SUBROUTINE initialise_ice_temperature( grid, ice, climate, init, SMB)
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_subclimate_region),        INTENT(IN)    :: climate
    TYPE(type_init_data_fields),         INTENT(IN)    :: init
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    
    ! Local variables
    INTEGER                                            :: i,j
    
    ! Special cases for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3') THEN
              
        ice%Ti_a( :,:,grid%i1:grid%i2) = 270._dp
        CALL sync
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6') THEN
              
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          ice%Ti_a( :,j,i) = climate%T2m( 1,j,i)
        END DO
        END DO
        CALL sync
        RETURN
              
      ELSEIF (C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'Bueler' .OR. &
              C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
              C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_F' .OR. &
              C%choice_benchmark_experiment == 'MISMIPplus' .OR. &
              C%choice_benchmark_experiment == 'MISOMIP1') THEN
              
        ice%Ti_a( :,:,grid%i1:grid%i2) = 270._dp
        CALL sync
        RETURN
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_ice_temperature!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF
    
    ! Realistic experiments
    ! =====================
      
    IF (C%is_restart) THEN
      ! Initialise with temperatures read from restart file
      ice%Ti_a( :,:,grid%i1:grid%i2) = init%Ti( :,:,grid%i1:grid%i2)
      CALL sync
      RETURN
    END IF ! IF (C%is_restart) THEN
    
    ! Not a restart; initialise from scratch with a simple linear profile
    
    ! First set all temperatures to -10C so thermal properties can be determined
    ice%Ti_a( :,:,grid%i1:grid%i2) = 260._dp
    CALL sync
   
    ! Calculate Ti_pmp
    CALL ice_physical_properties( grid, ice, C%start_time_of_run)
      
    ! Initialise with the Robin solution
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      CALL replace_Ti_with_robin_solution( ice, climate, SMB, ice%Ti_a, i,j)    
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE initialise_ice_temperature 
  
END MODULE thermodynamics_module
