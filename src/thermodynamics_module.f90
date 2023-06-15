MODULE thermodynamics_module

  ! Contains all the routines for calculating the englacial temperature field.

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, inquire_restart_file_temperature, read_restart_file_temperature
  USE parameters_module
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_climate_model, type_SMB_model, &
                                             type_restart_data, type_ocean_snapshot_regional
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             tridiagonal_solve, map_square_to_square_cons_2nd_order_3D, transpose_dp_3D, &
                                             interpolate_ocean_depth, vertical_average
  USE derivatives_and_grids_module,    ONLY: zeta, calculate_zeta_derivatives, ddx_a_to_a_2D, ddy_a_to_a_2D, Neumann_BC_a_3D, &
                                             map_a_to_cx_2D, map_a_to_cy_2D, map_a_to_cx_3D, map_a_to_cy_3D

  IMPLICIT NONE

CONTAINS

! == Run the chosen thermodynamics model
  SUBROUTINE run_thermo_model( grid, ice, climate, ocean, SMB, time, do_solve_heat_equation)
    ! Run the thermodynamics model. If so specified, solve the heat equation;
    ! if not, only prescribe a vertically uniform temperature to newly ice-covered grid cells.

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_model)   ,          INTENT(IN)    :: climate
    TYPE(type_ocean_snapshot_regional),   INTENT(IN)    :: ocean
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB
    REAL(dp),                             INTENT(IN)    :: time
    LOGICAL,                              INTENT(IN)    :: do_solve_heat_equation

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'run_thermo_model'
    INTEGER                                             :: i,j
    REAL(dp)                                            :: T_surf_annual

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_thermo_model == 'none') THEN
      ! No need to do anything
      ! NOTE: choice_ice_rheology_model should be set to "uniform"!
    ELSEIF (C%choice_thermo_model == '3D_heat_equation') THEN
      ! Solve the 3-D heat equation

      ! NOTE: solved asynchronously from the ice dynamical equations.
      !       Since newly ice-covered pixels won't have a temperature assigned
      !       until the heat equation is solved again, treat these separately every time step.

      ! Prescribe a simple temperature profile to newly ice-covered grid cells.
      DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
      DO j = 2, grid%ny-1

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
            ! Ice probably came from surface accumulation; initialise with a vertically uniform temperature.

            T_surf_annual = MIN( SUM( climate%T2m( :,j,i)) / 12._dp, T0)
            ice%Ti_a( :,j,i) = T_surf_annual

          END IF

        END IF ! IF (ice%mask_ice_a( j,i) == 1 .AND. ice%mask_ice_a_prev( j,i) == 0) THEN

      END DO
      END DO
      CALL sync

      ! Calculate various physical terms
      CALL calc_bottom_frictional_heating( grid, ice)
      CALL calc_heat_capacity(             grid, ice)
      CALL calc_thermal_conductivity(      grid, ice)
      CALL calc_pressure_melting_point(    grid, ice)

      ! If so specified, solve the heat equation
      IF (do_solve_heat_equation) CALL solve_3D_heat_equation( grid, ice, climate, ocean, SMB)

      ! Safety
      CALL check_for_NaN_dp_3D( ice%Ti_a, 'ice%Ti_a')

    ELSE
      CALL crash('unknown choice_thermo_model "' // TRIM( C%choice_thermo_model) // '"!')
    END IF

    ! Calculate the ice flow factor for the new temperature solution
    CALL calc_ice_rheology( grid, ice, time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_thermo_model

! == Solve the 3-D heat equation
  SUBROUTINE solve_3D_heat_equation( grid, ice, climate, ocean, SMB)

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_model)   ,          INTENT(IN)    :: climate
    TYPE(type_ocean_snapshot_regional),   INTENT(IN)    :: ocean
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'solve_3D_heat_equation'
    INTEGER                                             :: i,j,k
    REAL(dp), DIMENSION(:,:,:), POINTER                 :: Ti_new
    INTEGER                                             :: wTi_new
    REAL(dp)                                            :: internal_heating, u_times_dT_dx_upwind, v_times_dT_dy_upwind, f1, f2, f3
    REAL(dp), DIMENSION(2:C%nz)                         :: alpha
    REAL(dp), DIMENSION(C%nz)                           :: beta
    REAL(dp), DIMENSION(C%nz-1)                         :: gamma
    REAL(dp), DIMENSION(C%nz)                           :: delta
    INTEGER,  DIMENSION(:,:  ), POINTER                 ::  is_unstable
    INTEGER                                             :: wis_unstable
    LOGICAL                                             :: hasnan
    INTEGER                                             :: n_unstable
    REAL(dp)                                            :: depth
    REAL(dp), DIMENSION(:,:  ), POINTER                 ::  T_ocean_at_shelf_base
    INTEGER                                             :: wT_ocean_at_shelf_base

    ! Add routine to path
    CALL init_routine( routine_name)

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

    ! Find ocean temperature at the shelf base
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, T_ocean_at_shelf_base, wT_ocean_at_shelf_base)

    IF (C%choice_ocean_model == 'none') THEN
      ! No ocean data available; use constant PD value

      T_ocean_at_shelf_base( :,grid%i1:grid%i2) = T0 + C%T_ocean_mean_PD_ANT
      CALL sync

    ELSE ! IF (C%choice_ocean_model == 'none') THEN
      ! Calculate shelf base temperature from ocean data

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny

        IF (ice%mask_shelf_a( j,i) == 1) THEN
          depth = MAX( 0.1_dp, ice%Hi_a( j,i) - ice%Hs_a( j,i))   ! Depth is positive when below the sea surface!
          CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( :,j,i), depth, T_ocean_at_shelf_base( j,i))
        ELSE
          T_ocean_at_shelf_base( j,i) = 0._dp
        END IF

        ! NOTE: ocean data gives temperature in Celsius, thermodynamics wants Kelvin!
        T_ocean_at_shelf_base( j,i) = T_ocean_at_shelf_base( j,i) + T0

      END DO
      END DO
      CALL sync

    END IF ! IF (C%choice_ocean_model == 'none') THEN

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

        f2 = ice%dzeta_dt_a( k,j,i) + ice%dzeta_dx_a( k,j,i) * ice%u_3D_a( k,j,i) + ice%dzeta_dy_a( k,j,i) * ice%v_3D_a( k,j,i) + ice%dzeta_dz_a( j,i) * ice%w_3D_a( k,j,i)

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
        delta( C%nz) = MIN( T0, MIN( ice%Ti_pmp_a( C%nz,j,i), T_ocean_at_shelf_base( j,i) ))
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
      IF (MINVAL( Ti_new( :,j,i)) < 150._dp .OR. hasnan) THEN
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
    IF (n_unstable < CEILING( REAL( grid%nx*grid%ny) / 100._dp)) THEN
      ! Instability is limited to an acceptably small number (< 1%) of grid cells;
      ! replace the temperature profile in those cells with the Robin solution

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (is_unstable( j,i) == 1) CALL replace_Ti_with_robin_solution( ice, climate, ocean, SMB, Ti_new, i,j)
      END DO
      END DO
      CALL sync

    ELSE
      ! An unacceptably large number of grid cells was unstable; throw an error.
      CALL crash('heat equation solver unstable for more than 1% of grid cells!')
    END IF

    ! Move the new temperature field to the ice data structure
    ice%Ti_a( :,:,grid%i1:grid%i2) = Ti_new( :,:,grid%i1:grid%i2)
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wTi_new)
    CALL deallocate_shared( wis_unstable)
    CALL deallocate_shared( wT_ocean_at_shelf_base)

    ! Safety
    CALL check_for_NaN_dp_3D( ice%Ti_a, 'ice%Ti_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_3D_heat_equation

! == The Robin temperature solution
  SUBROUTINE replace_Ti_with_robin_solution( ice, climate, ocean, SMB, Ti, i,j)
    ! This function calculates for one horizontal grid point the temperature profiles
    ! using the surface temperature and the geothermal heat flux as boundary conditions.
    ! See Robin solution in: Cuffey & Paterson 2010, 4th ed, chapter 9, eq. (9.13) - (9.22).

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_climate_model)   ,          INTENT(IN)    :: climate
    TYPE(type_ocean_snapshot_regional),   INTENT(IN)    :: ocean
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB
    REAL(dp), DIMENSION(:,:,:),           INTENT(INOUT) :: Ti
    INTEGER,                              INTENT(IN)    :: i,j

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'replace_Ti_with_robin_solution'
    INTEGER                                             :: k
    REAL(dp)                                            :: Ts
    REAL(dp)                                            :: thermal_length_scale
    REAL(dp)                                            :: distance_above_bed
    REAL(dp)                                            :: erf1
    REAL(dp)                                            :: erf2

    REAL(dp)                                            :: thermal_conductivity_robin
    REAL(dp)                                            :: thermal_diffusivity_robin
    REAL(dp)                                            :: bottom_temperature_gradient_robin

    REAL(dp), PARAMETER                                 :: kappa_0_ice_conductivity     = 9.828_dp  ! The linear constant in the thermal conductivity of ice [J m^-1 K^-1 s^-1], see equation (12.6), Ritz (1987), Cuffey & Paterson (2010, p. 400), Zwinger (2007)
    REAL(dp), PARAMETER                                 :: kappa_e_ice_conductivity     = 0.0057_dp ! The exponent constant in the thermal conductivity of ice [K^-1], see equation (12.6), Ritz (1987), Cuffey & Paterson (2010, p. 400), Zwinger (2007)
    REAL(dp), PARAMETER                                 :: c_0_specific_heat            = 2127.5_dp ! The constant in the specific heat capacity of ice [J kg^-1 K^-1], see equation (12.5), Zwinger (2007), Cuffey & Paterson (2010, p. 400)

    REAL(dp)                                            :: depth
    REAL(dp)                                            :: T_ocean_at_shelf_base

    ! Add routine to path
    CALL init_routine( routine_name)

    thermal_conductivity_robin        = kappa_0_ice_conductivity * sec_per_year * EXP(-kappa_e_ice_conductivity * T0)  ! Thermal conductivity            [J m^-1 K^-1 y^-1]
    thermal_diffusivity_robin         = thermal_conductivity_robin / (ice_density * c_0_specific_heat)                 ! Thermal diffusivity             [m^2 y^-1]
    bottom_temperature_gradient_robin = -ice%GHF_a( j,i) / thermal_conductivity_robin                                ! Temperature gradient at bedrock

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

      IF (C%choice_ocean_model == 'none') THEN
        ! No ocean data available; use constant PD value

        T_ocean_at_shelf_base = T0 + C%T_ocean_mean_PD_ANT

      ELSE ! IF (C%choice_ocean_model == 'none') THEN
        ! Calculate shelf base temperature from ocean data

        depth = MAX( 0.1_dp, ice%Hi_a( j,i) - ice%Hs_a( j,i))   ! Depth is positive when below the sea surface!
        CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( :,j,i), depth, T_ocean_at_shelf_base)

      END IF

      Ti( :,j,i) = Ts + C%zeta(:) * (T0 + T_ocean_at_shelf_base - Ts)

    ELSE

      ! No ice present: use Ts everywhere
      Ti( :,j,i) = Ts

    END IF

    ! Correct all temperatures above T_pmp:
    DO k = 1, C%nz
      Ti( k,j,i) = MIN( Ti( k,j,i), ice%Ti_pmp_a( k,j,i))
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE replace_Ti_with_robin_solution

! == Calculate various physical terms
  SUBROUTINE calc_bottom_frictional_heating( grid, ice)
    ! Calculation of the frictional heating at the bottom due to sliding at the sheet/Gl - bedrock interface.

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_grid),                     INTENT(IN)    :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bottom_frictional_heating'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Exception for when no sliding can occur
    IF (C%choice_ice_dynamics == 'SIA' .OR. C%choice_sliding_law == 'no_sliding') THEN
      ice%frictional_heating_a( :,grid%i1:grid%i2) = 0._dp
      CALL sync
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      IF (ice%mask_sheet_a( j,i) == 1) THEN
        ice%frictional_heating_a( j,i) = ice%beta_a( j,i) * (ice%U_base_a( j,i)**2 + ice%V_base_a( j,i)**2)
      ELSE
        ice%frictional_heating_a( j,i) = 0._dp
      END IF
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bottom_frictional_heating
  SUBROUTINE calc_heat_capacity( grid, ice)
    ! Calculate the heat capacity of the ice

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_heat_capacity'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_heat_capacity == 'uniform') THEN
      ! Apply a uniform value for the heat capacity

      ice%Cpi_a( :,:,grid%i1:grid%i2) = C%uniform_ice_heat_capacity
      CALL sync

    ELSEIF (C%choice_ice_heat_capacity == 'Pounder1965') THEN
      ! Calculate the heat capacity of ice according to Pounder: The Physics of Ice (1965)

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%Cpi_a( :,j,i) = 2115.3_dp + 7.79293_dp * (ice%Ti_a( :,j,i) - T0)
      END DO
      END DO
      CALL sync

    ELSE
      CALL crash('unknown choice_ice_heat_capacity "' // TRIM( C%choice_ice_heat_capacity) // '"!')
    END IF

    ! Safety
    CALL check_for_NaN_dp_3D( ice%Cpi_a, 'ice%Cpi_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_heat_capacity
  SUBROUTINE calc_thermal_conductivity( grid, ice)
    ! Calculate the thermal conductivity of the ice

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_thermal_conductivity'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_thermal_conductivity == 'uniform') THEN
      ! Apply a uniform value for the thermal conductivity

      ice%Ki_a( :,:,grid%i1:grid%i2) = C%uniform_ice_thermal_conductivity
      CALL sync

    ELSEIF (C%choice_ice_thermal_conductivity == 'Ritz1987') THEN
      ! Calculate the thermal conductivity of ice according to Ritz (1987)

      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%Ki_a( :,j,i) = 3.101E+08_dp * EXP(-0.0057_dp * ice%Ti_a( :,j,i))
      END DO
      END DO
      CALL sync

    ELSE
      CALL crash('unknown choice_ice_thermal_conductivity "' // TRIM( C%choice_ice_thermal_conductivity) // '"!')
    END IF

    ! Safety
    CALL check_for_NaN_dp_3D( ice%Ki_a, 'ice%Ki_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_thermal_conductivity
  SUBROUTINE calc_pressure_melting_point( grid, ice)
    ! Calculate the pressure melting point of the ice according to Huybrechts (1992)

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pressure_melting_point'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%Ti_pmp_a( :,j,i) = T0 - CC * ice%Hi_a( j,i) * C%zeta
    END DO
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_3D( ice%Ti_pmp_a, 'ice%Ti_pmp_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pressure_melting_point

! == Calculate the  flow factor A in Glen's flow law
  SUBROUTINE calc_ice_rheology( grid, ice, time)
    ! Calculate the flow factor A in Glen's flow law

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_ice_rheology'
    INTEGER                                            :: i,j,k
    REAL(dp), DIMENSION(C%nZ)                          :: prof
    REAL(dp), PARAMETER                                :: A_low_temp  = 1.14E-05_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: A_high_temp = 5.47E+10_dp   ! [Pa^-3 yr^-1] The constant a in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_low_temp  = 6.0E+04_dp    ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp), PARAMETER                                :: Q_high_temp = 13.9E+04_dp   ! [J mol^-1] Activation energy for creep in the Arrhenius relationship
    REAL(dp)                                           :: A_flow_MISMIP

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_rheology == 'uniform') THEN
      ! Apply a uniform value for the ice flow factor

      ice%A_flow_3D_a( :,:,grid%i1:grid%i2) = C%uniform_flow_factor
      CALL sync

    ELSEIF (C%choice_ice_rheology == 'Huybrechts1992') THEN

      ! Calculate the ice flow factor as a function of the ice temperature according to the Arrhenius relationship (Huybrechts, 1992)
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny

        DO k = 1, C%nz
          IF (ice%mask_ice_a( j,i) == 1) THEN
            IF (ice%Ti_a( k,j,i) < 263.15_dp) THEN
              ice%A_flow_3D_a( k,j,i) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * ice%Ti_a( k,j,i)))
            ELSE
              ice%A_flow_3D_a( k,j,i) = A_high_temp * EXP(-Q_high_temp / (R_gas * ice%Ti_a( k,j,i)))
            END IF
          ELSE
            IF (C%choice_ice_margin == 'BC') THEN
              ice%A_flow_3D_a( k,j,i) = 0._dp
            ELSEIF (C%choice_ice_margin == 'infinite_slab') THEN
              ! In the "infinite slab" case, calculate effective viscosity everywhere
              ! (even when there's technically no ice present)
              ice%A_flow_3D_a( k,j,i) = A_low_temp  * EXP(-Q_low_temp  / (R_gas * 263.15_dp))
            ELSE
              CALL crash('unknown choice_ice_margin "' // TRIM( C%choice_ice_margin) // '"!')
            END IF
          END IF
        END DO ! DO k = 1, C%nz

      END DO
      END DO
      CALL sync

    ELSEIF (C%choice_ice_rheology == 'MISMIP_mod') THEN
      ! The time-dependent, step-wise changing uniform flow factor in the MISMIP_mod experiment

      A_flow_MISMIP = 1.0E-16_dp
      IF     (time < 15000._dp) THEN
        A_flow_MISMIP = 1.0E-16_dp
      ELSEIF (time < 30000._dp) THEN
        A_flow_MISMIP = 1.0E-17_dp
      ELSEIF (time < 45000._dp) THEN
        A_flow_MISMIP = 1.0E-16_dp
      END IF

      ice%A_flow_3D_a( :,:,grid%i1:grid%i2) = A_flow_MISMIP
      CALL sync

    ELSE
      CALL crash('unknown choice_ice_rheology "' // TRIM( C%choice_ice_rheology) // '"!')
    END IF

    ! Apply the flow enhancement factors
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF     (ice%mask_sheet_a( j,i) == 1) THEN
        ice%A_flow_3D_a( :,j,i) = ice%A_flow_3D_a( :,j,i) * C%m_enh_sheet
      ELSEIF (ice%mask_shelf_a( j,i) == 1) THEN
        ice%A_flow_3D_a( :,j,i) = ice%A_flow_3D_a( :,j,i) * C%m_enh_shelf
      END IF
    END DO
    END DO
    CALL sync

    ! Calculate vertical average
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      prof = ice%A_flow_3D_a( :,j,i)
      ice%A_flow_vav_a( j,i) = vertical_average( prof)
    END DO
    END DO
    CALL sync

    ! Map flow factor to staggered grids
    CALL map_a_to_cx_3D( grid, ice%A_flow_3D_a,   ice%A_flow_3D_cx)
    CALL map_a_to_cy_3D( grid, ice%A_flow_3D_a,   ice%A_flow_3D_cy)
    CALL map_a_to_cx_2D( grid, ice%A_flow_vav_a,  ice%A_flow_vav_cx)
    CALL map_a_to_cy_2D( grid, ice%A_flow_vav_a,  ice%A_flow_vav_cy)

    ! Safety
    CALL check_for_NaN_dp_3D( ice%A_flow_3D_a , 'ice%A_flow_3D_a' )
    CALL check_for_NaN_dp_2D( ice%A_flow_vav_a, 'ice%A_flow_vav_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_ice_rheology

! == Initialise the englacial ice temperature at the start of a simulation
  SUBROUTINE initialise_ice_temperature( grid, ice, climate, ocean, SMB, region_name)
    ! Initialise the englacial ice temperature at the start of a simulation

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_model)   ,          INTENT(IN)    :: climate
    TYPE(type_ocean_snapshot_regional),   INTENT(IN)    :: ocean
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_temperature'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) '  Initialising ice temperature profile "', TRIM(C%choice_initial_ice_temperature), '"...'

    IF     (C%choice_initial_ice_temperature == 'uniform') THEN
      ! Simple uniform temperature
      CALL initialise_ice_temperature_uniform( grid, ice)
    ELSEIF (C%choice_initial_ice_temperature == 'linear') THEN
      ! Simple linear temperature profile
      CALL initialise_ice_temperature_linear( grid, ice, climate)
    ELSEIF (C%choice_initial_ice_temperature == 'Robin') THEN
      ! Initialise with the Robin solution
      CALL initialise_ice_temperature_Robin( grid, ice, climate, ocean, SMB)
    ELSEIF (C%choice_initial_ice_temperature == 'restart') THEN
      ! Initialise with the temperature field from the provided restart file
      CALL initialise_ice_temperature_restart( grid, ice, region_name)
    ELSE
      CALL crash('unknown choice_initial_ice_temperature "' // TRIM( C%choice_initial_ice_temperature) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature
  SUBROUTINE initialise_ice_temperature_uniform( grid, ice)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Simple uniform temperature

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_temperature_uniform'
    INTEGER                                            :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (ice%Hi_a( j,i) > 0._dp) THEN
        ice%Ti_a( :,j,i) = C%uniform_ice_temperature
      ELSE
        ice%Ti_a( :,j,i) = 0._dp
      END IF

    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature_uniform
  SUBROUTINE initialise_ice_temperature_linear( grid, ice, climate)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Simple linear temperature profile

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_model)   ,          INTENT(IN)    :: climate

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_temperature_linear'
    INTEGER                                             :: i,j
    REAL(dp)                                            :: T_surf_annual, T_PMP_base

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (ice%Hi_a( j,i) > 0._dp) THEN
        T_surf_annual = MIN( SUM( climate%T2m( :,j,i)) / 12._dp, T0)
        T_PMP_base    = T0 - (ice%Hi_a( j,i) * 8.7E-04_dp)
        ice%Ti_a( :,j,i) = T_surf_annual - C%zeta * (T_surf_annual - T_PMP_base)
      ELSE
        ice%Ti_a( :,j,i) = 0._dp
      END IF

    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature_linear
  SUBROUTINE initialise_ice_temperature_Robin( grid, ice, climate, ocean, SMB)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Initialise with the Robin solution

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                      INTENT(IN)    :: grid
    TYPE(type_ice_model),                 INTENT(INOUT) :: ice
    TYPE(type_climate_model)   ,          INTENT(IN)    :: climate
    TYPE(type_ocean_snapshot_regional),   INTENT(IN)    :: ocean
    TYPE(type_SMB_model),                 INTENT(IN)    :: SMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_temperature_Robin'
    INTEGER                                             :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate Ti_pmp
    CALL calc_pressure_melting_point( grid, ice)

    ! Initialise with the Robin solution
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      CALL replace_Ti_with_robin_solution( ice, climate, ocean, SMB, ice%Ti_a, i,j)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature_Robin
  SUBROUTINE initialise_ice_temperature_restart( grid, ice, region_name)
    ! Initialise the englacial ice temperature at the start of a simulation
    !
    ! Initialise with the temperature field from the provided restart file

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_temperature_restart'
    CHARACTER(LEN=256)                                 :: filename_restart
    REAL(dp)                                           :: time_to_restart_from
    TYPE(type_restart_data)                            :: restart

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Assume that temperature and geometry are read from the same restart file
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

    ! Inquire if all the required fields are present in the specified NetCDF file,
    ! and determine the dimensions of the memory to be allocated.
    CALL allocate_shared_int_0D( restart%nx, restart%wnx)
    CALL allocate_shared_int_0D( restart%ny, restart%wny)
    CALL allocate_shared_int_0D( restart%nz, restart%wnz)
    CALL allocate_shared_int_0D( restart%nt, restart%wnt)
    IF (par%master) THEN
      restart%netcdf%filename = filename_restart
      CALL inquire_restart_file_temperature( restart)
    END IF
    CALL sync

    ! Allocate memory for raw data
    CALL allocate_shared_dp_1D( restart%nx, restart%x,    restart%wx   )
    CALL allocate_shared_dp_1D( restart%ny, restart%y,    restart%wy   )
    CALL allocate_shared_dp_1D( restart%nz, restart%zeta, restart%wzeta)
    CALL allocate_shared_dp_1D( restart%nt, restart%time, restart%wtime)

    CALL allocate_shared_dp_3D( restart%nx, restart%ny, restart%nz, restart%Ti,               restart%wTi              )

    ! Read data from input file
    IF (par%master) CALL read_restart_file_temperature( restart, time_to_restart_from)
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_3D( restart%Ti, 'restart%Ti')

    ! Since we want data represented as [j,i] internally, transpose the data we just read.
    CALL transpose_dp_3D( restart%Ti, restart%wTi)

    ! Map (transposed) raw data to the model grid
    CALL map_square_to_square_cons_2nd_order_3D( restart%nx, restart%ny, restart%x, restart%y, grid%nx, grid%ny, grid%x, grid%y, restart%Ti, ice%Ti_a)

    ! Deallocate raw data
    CALL deallocate_shared( restart%wnx              )
    CALL deallocate_shared( restart%wny              )
    CALL deallocate_shared( restart%wnz              )
    CALL deallocate_shared( restart%wnt              )
    CALL deallocate_shared( restart%wx               )
    CALL deallocate_shared( restart%wy               )
    CALL deallocate_shared( restart%wzeta            )
    CALL deallocate_shared( restart%wtime            )
    CALL deallocate_shared( restart%wTi              )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ice_temperature_restart

END MODULE thermodynamics_module
