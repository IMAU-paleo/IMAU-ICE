MODULE IMAU_ICE_main_model

  ! Contains all the routines for initialising and running the IMAU-ICE regional ice-sheet model.

  USE, INTRINSIC :: ISO_C_BINDING,         ONLY: c_backspace
  USE mpi
  USE configuration_module,                ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,                     ONLY: par, sync, cerr, ierr, &
                                                 allocate_shared_int_0D, allocate_shared_dp_0D, &
                                                 allocate_shared_int_1D, allocate_shared_dp_1D, &
                                                 allocate_shared_int_2D, allocate_shared_dp_2D, &
                                                     allocate_shared_int_3D, allocate_shared_dp_3D, &
                                                 allocate_shared_bool_0D, &
                                                 deallocate_shared, partition_list
  USE data_types_module,                   ONLY: type_model_region, type_ice_model, type_reference_geometry, &
                                                 type_SMB_model, type_BMB_model, type_forcing_data, type_grid, &
                                                 type_ocean_matrix_global
  USE utilities_module,                    ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                                 check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                                 inverse_oblique_sg_projection, surface_elevation, is_floating
  USE parameters_module
  USE reference_fields_module,             ONLY: initialise_reference_geometries
  USE netcdf_output_module,                ONLY: write_to_restart_file_grid, write_to_help_fields_file_grid, create_restart_file_grid, &
                                                 create_help_fields_file_grid, create_regional_scalar_file
  USE forcing_module,                      ONLY: forcing, initialise_geothermal_heat_flux_regional, update_sealevel_record_at_model_time
  USE general_ice_model_data_module,       ONLY: initialise_basins, initialise_mask_noice
  USE ice_velocity_module,                 ONLY: solve_DIVA
  USE ice_dynamics_module,                 ONLY: initialise_ice_model,              run_ice_model, update_ice_thickness
  USE thermodynamics_module,               ONLY: initialise_ice_temperature,        run_thermo_model, calc_ice_rheology
  USE ocean_module,                        ONLY: initialise_ocean_model_regional,   run_ocean_model, ocean_temperature_inversion, write_inverted_ocean_temperature_to_file
  USE climate_module,                      ONLY: initialise_climate_model,          run_climate_model
  USE SMB_module,                          ONLY: initialise_SMB_model,              run_SMB_model
  USE BMB_module,                          ONLY: initialise_BMB_model,              run_BMB_model
  USE isotopes_module,                     ONLY: initialise_isotopes_model,         run_isotopes_model
  USE bedrock_ELRA_module,                 ONLY: initialise_ELRA_model,             run_ELRA_model
# if (defined(DO_SELEN))
  USE SELEN_main_module,                   ONLY: apply_SELEN_bed_geoid_deformation_rates
# endif
  USE scalar_data_output_module,           ONLY: write_regional_scalar_data
  USE basal_conditions_and_sliding_module, ONLY: basal_inversion_geo, write_inverted_bed_roughness_to_file

  USE netcdf_debug_module,                 ONLY: save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, save_variable_as_netcdf_int_3D, &
                                                 save_variable_as_netcdf_dp_1D,  save_variable_as_netcdf_dp_2D,  save_variable_as_netcdf_dp_3D

  IMPLICIT NONE

CONTAINS

  SUBROUTINE run_model( region, t_end)
    ! Run the model until t_end (usually a 100 years further)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end

    ! Local variables:
    !CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_model'
    CHARACTER(LEN=256)                                 :: routine_name
    REAL(dp)                                           :: tstart, tstop, t1, t2, dt_ave
    INTEGER                                            :: it
    CHARACTER(LEN=9)                                   :: r_time, r_step, r_adv


    ! Add routine to path
    routine_name = 'run_model('  //  region%name  //  ')'
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,'(A,A3,A2,A,A,F9.3,A,F9.3,A)') '  Running model region ', region%name, ' (', TRIM(region%long_name), &
                                                           ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' kyr'

    ! Computation time tracking
    region%tcomp_total          = 0._dp
    region%tcomp_ice            = 0._dp
    region%tcomp_thermo         = 0._dp
    region%tcomp_climate        = 0._dp
    region%tcomp_GIA            = 0._dp

    tstart = MPI_WTIME()

  ! ====================================
  ! ===== The main model time loop =====
  ! ====================================

    it = 0
    dt_ave = 0._dp
    DO WHILE (region%time < t_end)
      it = it + 1

    ! GIA
    ! ===

      t1 = MPI_WTIME()
      IF     (C%choice_GIA_model == 'none') THEN
        ! Nothing to be done
      ELSEIF (C%choice_GIA_model == 'ELRA') THEN
        CALL run_ELRA_model( region)
# if (defined(DO_SELEN))
      ELSEIF (C%choice_GIA_model == 'SELEN') THEN
        CALL apply_SELEN_bed_geoid_deformation_rates( region)
# endif
      ELSEIF (C%choice_GIA_model == '3DGIA') THEN       ! CvC
        CALL update_Hb_with_3D_GIA_model_output(region%grid, region%ice, region%time, region%refgeo_init)
      ELSE
        CALL crash('unknown choice_GIA_model "' // TRIM(C%choice_GIA_model) // '"!')
      END IF
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_GIA = region%tcomp_GIA + t2 - t1

    ! Ice dynamics
    ! ============

      ! Calculate ice velocities and the resulting change in ice geometry
      ! NOTE: geometry is not updated yet; this happens at the end of the time loop
      t1 = MPI_WTIME()
      CALL run_ice_model( region, t_end)
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_ice = region%tcomp_ice + t2 - t1

    ! == Time display
    ! ===============

      if (par%master .and. C%do_time_display) then
        if (region%time + region%dt < t_end) then
          r_adv = "no"
          write(r_time,"(F9.3)") min(region%time,t_end) / 1000._dp
          write(r_step,"(F6.3)") max(region%dt,0.001_dp)
          write(*,"(A)",advance=trim(r_adv)) repeat(c_backspace,999) // &
                  "   t = " // trim(r_time) // " kyr - dt = " // trim(r_step) // " yr"
        else
          r_adv = "yes"
          write(r_time,"(F9.3)") min(region%time,t_end) / 1000._dp
          write(r_step, "(F6.3)") dt_ave / real(it,dp)
          write(*,"(A)",advance=trim(r_adv)) repeat(c_backspace,999) // &
                "   t = " // trim(r_time) // " kyr - dt_ave = " // trim(r_step) // " yr"
        end if
        if (region%do_output .OR. region%do_output_restart) then
          r_adv = "no"
          write(*,"(A)",advance=trim(r_adv)) repeat(c_backspace,999)
        end if
      end if

    ! Climate , SMB and BMB
    ! =====================

      t1 = MPI_WTIME()

      ! Run the climate model
      IF (region%do_climate) THEN
        CALL run_climate_model( region, region%time)
      END IF

      ! Run the ocean model
      IF (region%do_ocean) THEN
        CALL run_ocean_model( region%grid, region%ice, region%ocean_matrix, region%climate, region%name, region%time, region%refgeo_PD)
      END IF

      ! Run the SMB model
      IF (region%do_SMB) THEN
        CALL run_SMB_model( region%grid, region%ice, region%climate, region%time, region%SMB, region%mask_noice)
      END IF

      ! Run the BMB model
      IF (region%do_BMB) THEN
        CALL run_BMB_model( region%grid, region%ice, region%ocean_matrix%applied, region%BMB, region%name, region%time, region%refgeo_PD)
      END IF

      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_climate = region%tcomp_climate + t2 - t1

    ! Thermodynamics
    ! ==============

      t1 = MPI_WTIME()
      CALL run_thermo_model( region%grid, region%ice, region%climate, region%ocean_matrix%applied, region%SMB, region%time, do_solve_heat_equation = region%do_thermo)
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_thermo = region%tcomp_thermo + t2 - t1

    ! Isotopes
    ! ========

      CALL run_isotopes_model( region)

    ! Geometry-based basal inversion
    ! ==============================

      IF (region%do_BIV) THEN
        IF (region%time > C%BIVgeo_t_start .AND. region%time < C%BIVgeo_t_end) THEN
          CALL basal_inversion_geo( region%grid, region%ice, region%refgeo_PD, C%BIVgeo_dt)
        END IF
      END IF

    ! Ocean temperature inversion
    ! ===========================

      IF (C%do_ocean_temperature_inversion .AND. region%do_BMB) THEN
        IF (region%time > C%ocean_temperature_inv_t_start .AND. region%time < C%ocean_temperature_inv_t_end) THEN
          ! Adjust ocean temperatures
          IF (C%do_asynchronous_BMB) THEN
            ! Use custom BMB time step
            CALL ocean_temperature_inversion( region%grid, region%ice, region%ocean_matrix%applied, region%refgeo_PD, C%dt_BMB)
          ELSE
            ! Use main model time step
            CALL ocean_temperature_inversion( region%grid, region%ice, region%ocean_matrix%applied, region%refgeo_PD, region%dt)
          END IF
        END IF
      END IF

    ! Time step and output
    ! ====================

      ! Write main output
      IF (region%do_output) THEN
        CALL write_to_help_fields_file_grid( region%help_fields_filename, region)
      END IF

      ! Write to restart file
      IF (region%do_output_restart) THEN
        CALL write_to_restart_file_grid( region%restart_filename, region)
      END IF

      ! Write scalar output
      CALL calculate_icesheet_volume_and_area(region)
      CALL write_regional_scalar_data( region, region%time)

      ! Update ice geometry and advance region time
      CALL update_ice_thickness( region%grid, region%ice, region%mask_noice, region%refgeo_PD, region%refgeo_GIAeq, region%time)
      IF (par%master) region%time = region%time + region%dt
      IF (par%master) dt_ave = dt_ave + region%dt
      CALL sync

      ! DENK DROM
      ! region%time = t_end

    END DO ! DO WHILE (region%time < t_end)

  ! ===========================================
  ! ===== End of the main model time loop =====
  ! ===========================================

    ! Write to NetCDF output one last time at the end of the simulation
    IF (region%time == C%end_time_of_run) THEN

      IF (C%choice_GIA_model == '3DGIA') THEN
        CALL calculate_output_for_3D_GIA_model(region%grid, region%ice, region%refgeo_GIAeq)
      END IF

      CALL write_to_restart_file_grid( region%restart_filename, region)
      CALL write_to_help_fields_file_grid( region%help_fields_filename, region)

      ! Write inverted bed roughness field to file
      IF (C%do_BIVgeo) THEN
        CALL write_inverted_bed_roughness_to_file( region%grid, region%ice)
      END IF
      ! Write inverted ocean temperature field to file
      IF (C%do_ocean_temperature_inversion) THEN
        CALL write_inverted_ocean_temperature_to_file( region%grid, region%ocean_matrix%applied)
      END IF

    END IF

    ! Determine total ice sheet area, volume, volume-above-flotation and GMSL contribution,
    ! used for writing to text output and in the inverse routine
    CALL calculate_icesheet_volume_and_area(region)

    ! Write to text output
    CALL write_regional_scalar_data( region, region%time)

    tstop = MPI_WTIME()
    region%tcomp_total = tstop - tstart

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_model

  ! 3D GIA subroutines
  SUBROUTINE read_dHb_3D_file( grid, ice)
  ! Caroline van Calcar 08/2022

    IMPLICIT NONE

    ! in/output variables:
    TYPE(type_grid),                  INTENT(IN)        :: grid
    TYPE(type_ice_model),             INTENT(INOUT)     :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'read_dHb_3D_file'
    INTEGER                                             :: i,j

    CALL init_routine( routine_name)

    IF (par%master) THEN
      ! open (91, file = TRIM( C%deformation_foldername) // '/dHb_3D.dat', status = 'old')
            ! open (91, file = '/home/caroline/HetGroteKoppelScript_IMAUICE/IMAU-ICE/Datasets/3DGIA/reference_surface_load.dat', status = 'old')

        open (91, file = '/home/caroline/HetGroteKoppelScript_IMAUICE/IMAU-ICE/Datasets/GIA_V2/dHb_3D.dat', status = 'old')
        do i = 1,grid%NX
          read(91,*) (ice%dHb_3D(j,i),j=1,grid%NY)
        enddo
      close(91)
    END IF
    CALL SYNC

    CALL finalise_routine( routine_name)

  END SUBROUTINE read_dHb_3D_file
  SUBROUTINE update_Hb_with_3D_GIA_model_output( grid, ice, time, refgeo_init)
    ! Caroline van Calcar 08/2022

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                  INTENT(IN)        :: grid
    TYPE(type_ice_model),             INTENT(INOUT)     :: ice
    TYPE(type_reference_geometry),    INTENT(IN)        :: refgeo_init
    REAL(dp),                         INTENT(IN)        :: time


    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'Update_Hb_with_3D_GIA_model_output'
    INTEGER                                             :: i,j

    CALL init_routine( routine_name)


    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%dHb_dt_a(j,i) = ice%dHb_3D(j,i) / (C%end_time_of_run - C%start_time_of_run) * C%dt_coupling
      ice%Hb_a(j,i) = refgeo_init%Hb(j,i) + ice%dHb_3D(j,i) * (time - C%start_time_of_run)/ (C%end_time_of_run - C%start_time_of_run)
      ice%dHb_a(j,i) = ice%Hb_a(j,i) - refgeo_init%Hb(j,i)
    END DO
    END DO

    CALL SYNC

  ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_Hb_with_3D_GIA_model_output
  SUBROUTINE calculate_output_for_3D_GIA_model(grid, ice, refgeo_GIAeq)
    ! Caroline van Calcar 08/2022

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_ice_model),             INTENT(INOUT)     :: ice
    TYPE(type_grid),                  INTENT(IN)        :: grid
    TYPE(type_reference_geometry),    INTENT(IN)    :: refgeo_GIAeq


    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'calculate_output_for_3D_GIA_model'
    INTEGER                                             :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (is_floating( ice%Hi_a( j,i), ice%Hb_a( j,i), 0._dp)) THEN
            IF (is_floating( refgeo_GIAeq%Hi( j,i), refgeo_GIAeq%Hb( j,i), 0._dp)) THEN
              ice%surface_load_rel( j,i) = 0; ! % dHi for the GIA model is zero (from floating to floating)
            ELSE
              ice%surface_load_rel( j,i) = -(refgeo_GIAeq%Hi( j,i)+refgeo_GIAeq%Hb( j,i)*seawater_density/ice_density); !%dHi for GIA model is minus the floating part of Hi,1
            END IF

      ELSE !% Hi,2 is grounded
            IF (is_floating( refgeo_GIAeq%Hi( j,i), refgeo_GIAeq%Hb( j,i), 0._dp)) THEN
                ice%surface_load_rel( j,i) = ice%Hi_a( j,i) +ice%Hb_a( j,i)*(seawater_density/ice_density);
            ELSE
                ice%surface_load_rel( j,i) = ice%Hi_a( j,i) - refgeo_GIAeq%Hi( j,i) + (ice%Hb_a( j,i) - refgeo_GIAeq%Hb( j,i))*(seawater_density/ice_density); ! % dHi for GIA model is the difference in ice thickness minus the difference in bedrock elevation*density compensation
            END IF
      END IF
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_output_for_3D_GIA_model

  ! Initialise the entire model region - read initial and PD data, initialise the ice dynamics, climate and SMB sub models
  SUBROUTINE initialise_model( region, name, ocean_matrix_global)

    USE climate_module, ONLY: map_glob_to_grid_2D

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),          INTENT(INOUT)     :: region
    CHARACTER(LEN=3),                 INTENT(IN)        :: name
    TYPE(type_ocean_matrix_global),   INTENT(INOUT)     :: ocean_matrix_global

    ! Local variables:
    !CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_model'
    CHARACTER(LEN=256)                                 :: routine_name
    CHARACTER(LEN=20)                                   :: short_filename
    INTEGER                                             :: n

    ! Add routine to path
    routine_name = 'initialise_model('  //  name  //  ')'
    CALL init_routine( routine_name)

    ! Region name
    region%name      = name
    IF     (region%name == 'NAM') THEN
      region%long_name = 'North America'
    ELSEIF (region%name == 'EAS') THEN
      region%long_name = 'Eurasia'
    ELSEIF (region%name == 'GRL') THEN
      region%long_name = 'Greenland'
    ELSEIF (region%name == 'ANT') THEN
      region%long_name = 'Antarctica'
    END IF

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Initialising model region ', region%name, ' (', TRIM(region%long_name), ')...'

    ! ===== Allocate memory for timers and scalars =====
    ! ==================================================

    CALL allocate_region_timers_and_scalars( region)

    ! ===== Initialise this region's grid =====
    ! =========================================

    CALL initialise_model_grid( region)
    IF (par%master) WRITE (0,'(A,F5.2,A,I4,A,I4,A)') '   Initialised model grid at ', region%grid%dx/1000._dp, ' km resolution: [', region%grid%nx, ' x ', region%grid%ny, '] pixels'

    ! ===== Set up output files ======
    ! =================================================

    ! Restart file
    short_filename = 'restart_NAM.nc'
    short_filename(9:11) = region%name
    DO n = 1, 256
      region%restart_filename(n:n) = ' '
    END DO
    region%restart_filename = TRIM(C%output_dir) // TRIM(short_filename)

    ! Help fields file
    short_filename = 'help_fields_NAM.nc'
    short_filename(13:15) = region%name
    DO n = 1, 256
      region%help_fields_filename(n:n) = ' '
    END DO
    region%help_fields_filename = TRIM(C%output_dir) // TRIM(short_filename)

    ! Create the (empty) NetCDF files
    CALL create_restart_file_grid(     region%restart_filename,     region%grid)
    CALL create_help_fields_file_grid( region%help_fields_filename, region%grid)

    ! ===== Initialise initial, present-day, and GIA equilibrium reference geometries =====
    ! =====================================================================================

    CALL initialise_reference_geometries( region%grid, region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq, region%name)

    ! ===== Define mask where no ice is allowed to form (i.e. Greenland in NAM and EAS, Ellesmere Island in GRL)
    ! ==========================================================================================================

    CALL initialise_mask_noice( region)

    ! ===== The ice dynamics model
    ! ============================

    CALL initialise_ice_model( region%grid, region%ice, region%refgeo_init, region%refgeo_PD, region%name)

    ! ===== Set sea level if prescribed externally =====
    ! ==================================================

    IF     (C%choice_sealevel_model == 'fixed') THEN
      region%ice%SL_a( :,region%grid%i1:region%grid%i2) = C%fixed_sealevel
    ELSEIF (C%choice_sealevel_model == 'eustatic' .OR. C%choice_sealevel_model == 'SELEN') THEN
      ! FIXME
    ELSEIF     (C%choice_sealevel_model == 'prescribed') THEN
      CALL update_sealevel_record_at_model_time( C%start_time_of_run)
      region%ice%SL_a( :,region%grid%i1:region%grid%i2) = forcing%sealevel_obs
    ELSE
      CALL crash('unknown choice_sealevel_model "' // TRIM(C%choice_sealevel_model) // '"!')
    END IF
    CALL sync

    ! ===== Define ice basins =====
    ! =============================

    ! Allocate shared memory
    CALL allocate_shared_int_2D( region%grid%ny, region%grid%nx, region%ice%basin_ID, region%ice%wbasin_ID)
    CALL allocate_shared_int_0D(                                 region%ice%nbasins,  region%ice%wnbasins )

    ! Define basins
    CALL initialise_basins( region%grid, region%ice%basin_ID, region%ice%nbasins, region%name)

    ! ===== The climate model =====
    ! =============================

    CALL initialise_climate_model( region)

    ! ===== The ocean model =====
    ! ===========================

    CALL initialise_ocean_model_regional( region, ocean_matrix_global)

    ! ===== The SMB model =====
    ! =========================

    CALL initialise_SMB_model( region%grid, region%ice, region%SMB, region%name)

    ! ===== The BMB model =====
    ! =========================

    CALL initialise_BMB_model( region%grid, region%ice, region%BMB, region%name)

    ! ===== The GIA model =====
    ! =========================

    IF     (C%choice_GIA_model == 'none') THEN
      ! Nothing to be done
    ELSEIF (C%choice_GIA_model == 'ELRA') THEN
      CALL initialise_GIA_model_grid( region)
      CALL initialise_ELRA_model( region%grid, region%grid_GIA, region%ice, region%refgeo_GIAeq)
# if (defined(DO_SELEN))
    ELSEIF (C%choice_GIA_model == 'SELEN') THEN
      CALL initialise_GIA_model_grid( region)
# endif
    ELSEIF (C%choice_GIA_model == '3DGIA') THEN
      CALL allocate_shared_dp_2D(        region%grid%ny  , region%grid%nx  , region%ice%dHb_3D                , region%ice%wdHb_3D)
      CALL allocate_shared_dp_2D(        region%grid%ny  , region%grid%nx  , region%ice%surface_load_rel , region%ice%wsurface_load_rel)
      CALL read_dHb_3D_file(region%grid, region%ice)
    ELSE
      CALL crash('unknown choice_GIA_model "' // TRIM(C%choice_GIA_model) // '"!')
    END IF

    ! ===== The isotopes model =====
    ! ==============================

    CALL initialise_isotopes_model( region)

    ! ===== Geothermal heat flux =====
    ! ================================

    CALL initialise_geothermal_heat_flux_regional( region%grid, region%ice, region%name)

    ! ===== Initialise the ice temperature profile =====
    ! ==================================================

    ! Run the climate and SMB models once, to get the correct surface temperature+SMB fields for the ice temperature initialisation
    CALL run_climate_model( region, C%start_time_of_run)
    CALL run_SMB_model( region%grid, region%ice, region%climate, C%start_time_of_run, region%SMB, region%mask_noice)

    ! Initialise the temperature field
    CALL initialise_ice_temperature( region%grid, region%ice, region%climate, region%ocean_matrix%applied, region%SMB, region%name)

    ! Initialise the rheology
    CALL calc_ice_rheology( region%grid, region%ice, C%start_time_of_run)

    ! ===== Scalar output (regionally integrated ice volume, SMB components, etc.)
    ! ============================================================================

    ! Create the file
    CALL create_regional_scalar_file( region%name, region%scalar_filename)

    ! Calculate and write the first entry (ice volume and area, GMSL contribution, isotope stuff)
    CALL calculate_PD_sealevel_contribution( region)
    CALL calculate_icesheet_volume_and_area(region)
    CALL write_regional_scalar_data( region, C%start_time_of_run)

    ! ===
    ! === If we're running with choice_ice_dynamics == "none", calculate a velocity field
    ! === once during initialisation (so that the thermodynamics are solved correctly)
    ! ===

    IF (C%choice_ice_dynamics == 'none') THEN
      C%choice_ice_dynamics = 'DIVA'
      CALL solve_DIVA( region%grid, region%ice)
      C%choice_ice_dynamics = 'none'
    END IF

    IF (par%master) WRITE (0,*) ' Finished initialising model region ', region%name, '.'

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_model
  SUBROUTINE allocate_region_timers_and_scalars( region)
    ! Allocate shared memory for this region's timers (used for the asynchronous coupling between the
    ! ice dynamics and the secondary model components), and for the scalars (integrated ice volume and
    ! area, SMB components, computation times, etc.)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_region_timers_and_scalars'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Timers and time steps
    ! =====================

    CALL allocate_shared_dp_0D(   region%time,             region%wtime            )
    CALL allocate_shared_dp_0D(   region%dt,               region%wdt              )
    CALL allocate_shared_dp_0D(   region%dt_prev,          region%wdt_prev         )
    CALL allocate_shared_dp_0D(   region%dt_crit_SIA,      region%wdt_crit_SIA     )
    CALL allocate_shared_dp_0D(   region%dt_crit_SSA,      region%wdt_crit_SSA     )
    CALL allocate_shared_dp_0D(   region%dt_crit_ice,      region%wdt_crit_ice     )
    CALL allocate_shared_dp_0D(   region%dt_crit_ice_prev, region%wdt_crit_ice_prev)

    region%dt               = C%dt_min
    region%dt_prev          = C%dt_min
    region%dt_crit_ice      = C%dt_min
    region%dt_crit_ice_prev = C%dt_min

    CALL allocate_shared_dp_0D(   region%t_last_SIA,            region%wt_last_SIA           )
    CALL allocate_shared_dp_0D(   region%t_next_SIA,            region%wt_next_SIA           )
    CALL allocate_shared_bool_0D( region%do_SIA,                region%wdo_SIA               )

    CALL allocate_shared_dp_0D(   region%t_last_SSA,            region%wt_last_SSA           )
    CALL allocate_shared_dp_0D(   region%t_next_SSA,            region%wt_next_SSA           )
    CALL allocate_shared_bool_0D( region%do_SSA,                region%wdo_SSA               )

    CALL allocate_shared_dp_0D(   region%t_last_DIVA,           region%wt_last_DIVA          )
    CALL allocate_shared_dp_0D(   region%t_next_DIVA,           region%wt_next_DIVA          )
    CALL allocate_shared_bool_0D( region%do_DIVA,               region%wdo_DIVA              )

    CALL allocate_shared_dp_0D(   region%t_last_thermo,         region%wt_last_thermo        )
    CALL allocate_shared_dp_0D(   region%t_next_thermo,         region%wt_next_thermo        )
    CALL allocate_shared_bool_0D( region%do_thermo,             region%wdo_thermo            )

    CALL allocate_shared_dp_0D(   region%t_last_output,         region%wt_last_output        )
    CALL allocate_shared_dp_0D(   region%t_next_output,         region%wt_next_output        )
    CALL allocate_shared_bool_0D( region%do_output,             region%wdo_output            )

    CALL allocate_shared_dp_0D(   region%t_last_output_restart, region%wt_last_output_restart)
    CALL allocate_shared_dp_0D(   region%t_next_output_restart, region%wt_next_output_restart)
    CALL allocate_shared_bool_0D( region%do_output_restart,     region%wdo_output_restart    )

    CALL allocate_shared_dp_0D(   region%t_last_climate,        region%wt_last_climate       )
    CALL allocate_shared_dp_0D(   region%t_next_climate,        region%wt_next_climate       )
    CALL allocate_shared_bool_0D( region%do_climate,            region%wdo_climate           )

    CALL allocate_shared_dp_0D(   region%t_last_ocean,          region%wt_last_ocean         )
    CALL allocate_shared_dp_0D(   region%t_next_ocean,          region%wt_next_ocean         )
    CALL allocate_shared_bool_0D( region%do_ocean,              region%wdo_ocean             )

    CALL allocate_shared_dp_0D(   region%t_last_SMB,            region%wt_last_SMB           )
    CALL allocate_shared_dp_0D(   region%t_next_SMB,            region%wt_next_SMB           )
    CALL allocate_shared_bool_0D( region%do_SMB,                region%wdo_SMB               )

    CALL allocate_shared_dp_0D(   region%t_last_BMB,            region%wt_last_BMB           )
    CALL allocate_shared_dp_0D(   region%t_next_BMB,            region%wt_next_BMB           )
    CALL allocate_shared_bool_0D( region%do_BMB,                region%wdo_BMB               )

    CALL allocate_shared_dp_0D(   region%t_last_ELRA,           region%wt_last_ELRA          )
    CALL allocate_shared_dp_0D(   region%t_next_ELRA,           region%wt_next_ELRA          )
    CALL allocate_shared_bool_0D( region%do_ELRA,               region%wdo_ELRA              )

    CALL allocate_shared_dp_0D(   region%t_last_BIV,            region%wt_last_BIV           )
    CALL allocate_shared_dp_0D(   region%t_next_BIV,            region%wt_next_BIV           )
    CALL allocate_shared_bool_0D( region%do_BIV,                region%wdo_BIV               )

    IF (par%master) THEN
      region%time                  = C%start_time_of_run
      region%dt                    = C%dt_min
      region%dt_prev               = C%dt_min

      region%t_last_SIA            = C%start_time_of_run
      region%t_next_SIA            = C%start_time_of_run
      region%do_SIA                = .TRUE.

      region%t_last_SSA            = C%start_time_of_run
      region%t_next_SSA            = C%start_time_of_run
      region%do_SSA                = .TRUE.

      region%t_last_DIVA           = C%start_time_of_run
      region%t_next_DIVA           = C%start_time_of_run
      region%do_DIVA               = .TRUE.

      region%t_last_thermo         = C%start_time_of_run
      region%t_next_thermo         = C%start_time_of_run + C%dt_thermo
      region%do_thermo             = .FALSE.

      region%t_last_climate        = C%start_time_of_run
      region%t_next_climate        = C%start_time_of_run
      region%do_climate            = .TRUE.

      region%t_last_ocean          = C%start_time_of_run
      region%t_next_ocean          = C%start_time_of_run
      region%do_ocean              = .TRUE.

      region%t_last_SMB            = C%start_time_of_run
      region%t_next_SMB            = C%start_time_of_run
      region%do_SMB                = .TRUE.

      region%t_last_BMB            = C%start_time_of_run
      region%t_next_BMB            = C%start_time_of_run
      region%do_BMB                = .TRUE.

      region%t_last_ELRA           = C%start_time_of_run
      region%t_next_ELRA           = C%start_time_of_run
      region%do_ELRA               = .TRUE.

      region%t_last_BIV            = C%start_time_of_run
      region%t_next_BIV            = C%start_time_of_run + C%BIVgeo_dt
      region%do_BIV                = .FALSE.

      region%t_last_output         = C%start_time_of_run
      region%t_next_output         = C%start_time_of_run
      region%do_output             = .TRUE.

      region%t_last_output_restart = C%start_time_of_run
      region%t_next_output_restart = C%start_time_of_run
      region%do_output_restart     = .TRUE.
    END IF

    ! ===== Scalars =====
    ! ===================

    ! Ice-sheet volume and area
    CALL allocate_shared_dp_0D( region%ice_area                     , region%wice_area                     )
    CALL allocate_shared_dp_0D( region%ice_volume                   , region%wice_volume                   )
    CALL allocate_shared_dp_0D( region%ice_volume_PD                , region%wice_volume_PD                )
    CALL allocate_shared_dp_0D( region%ice_volume_above_flotation   , region%wice_volume_above_flotation   )
    CALL allocate_shared_dp_0D( region%ice_volume_above_flotation_PD, region%wice_volume_above_flotation_PD)

    ! Regionally integrated SMB components
    CALL allocate_shared_dp_0D( region%int_T2m                      , region%wint_T2m                      )
    CALL allocate_shared_dp_0D( region%int_snowfall                 , region%wint_snowfall                 )
    CALL allocate_shared_dp_0D( region%int_rainfall                 , region%wint_rainfall                 )
    CALL allocate_shared_dp_0D( region%int_melt                     , region%wint_melt                     )
    CALL allocate_shared_dp_0D( region%int_refreezing               , region%wint_refreezing               )
    CALL allocate_shared_dp_0D( region%int_runoff                   , region%wint_runoff                   )
    CALL allocate_shared_dp_0D( region%int_SMB                      , region%wint_SMB                      )
    CALL allocate_shared_dp_0D( region%int_BMB                      , region%wint_BMB                      )
    CALL allocate_shared_dp_0D( region%int_MB                       , region%wint_MB                       )

    ! Englacial isotope content
    CALL allocate_shared_dp_0D( region%GMSL_contribution            , region%wGMSL_contribution            )
    CALL allocate_shared_dp_0D( region%mean_isotope_content         , region%wmean_isotope_content         )
    CALL allocate_shared_dp_0D( region%mean_isotope_content_PD      , region%wmean_isotope_content_PD      )
    CALL allocate_shared_dp_0D( region%d18O_contribution            , region%wd18O_contribution            )
    CALL allocate_shared_dp_0D( region%d18O_contribution_PD         , region%wd18O_contribution_PD         )

    ! Computation times
    CALL allocate_shared_dp_0D( region%tcomp_total   , region%wtcomp_total   )
    CALL allocate_shared_dp_0D( region%tcomp_ice     , region%wtcomp_ice     )
    CALL allocate_shared_dp_0D( region%tcomp_thermo  , region%wtcomp_thermo  )
    CALL allocate_shared_dp_0D( region%tcomp_climate , region%wtcomp_climate )
    CALL allocate_shared_dp_0D( region%tcomp_GIA     , region%wtcomp_GIA     )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_region_timers_and_scalars
  SUBROUTINE initialise_model_grid( region)
    ! Initialise the square grid for this model region

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_model_grid'
    INTEGER                                            :: nsx, nsy, i, j
    REAL(dp)                                           :: xmin, xmax, ymin, ymax, xmid, ymid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Assign dummy values to suppress compiler warnings
    xmin = 0._dp
    xmax = 0._dp
    ymin = 0._dp
    ymax = 0._dp
    xmid = 0._dp
    ymid = 0._dp
    nsx  = 0
    nsy  = 0

    ! Allocate shared memory
    CALL allocate_shared_int_0D( region%grid%nx,           region%grid%wnx          )
    CALL allocate_shared_int_0D( region%grid%ny,           region%grid%wny          )
    CALL allocate_shared_dp_0D(  region%grid%dx,           region%grid%wdx          )
    CALL allocate_shared_dp_0D(  region%grid%xmin,         region%grid%wxmin        )
    CALL allocate_shared_dp_0D(  region%grid%xmax,         region%grid%wxmax        )
    CALL allocate_shared_dp_0D(  region%grid%ymin,         region%grid%wymin        )
    CALL allocate_shared_dp_0D(  region%grid%ymax,         region%grid%wymax        )
    CALL allocate_shared_dp_0D(  region%grid%lambda_M,     region%grid%wlambda_M    )
    CALL allocate_shared_dp_0D(  region%grid%phi_M,        region%grid%wphi_M       )
    CALL allocate_shared_dp_0D(  region%grid%beta_stereo,  region%grid%wbeta_stereo )

    ! Let the Master do the work
    IF (par%master) THEN

      ! Resolution, domain size, and projection parameters for this region are determined from the config
      IF     (region%name == 'NAM') THEN
        xmin                     = C%xmin_NAM
        xmax                     = C%xmax_NAM
        ymin                     = C%ymin_NAM
        ymax                     = C%ymax_NAM
        region%grid%dx           = C%dx_NAM
        region%grid%lambda_M     = C%lambda_M_NAM
        region%grid%phi_M        = C%phi_M_NAM
        region%grid%beta_stereo  = C%beta_stereo_NAM
      ELSEIF (region%name == 'EAS') THEN
        xmin                     = C%xmin_EAS
        xmax                     = C%xmax_EAS
        ymin                     = C%ymin_EAS
        ymax                     = C%ymax_EAS
        region%grid%dx           = C%dx_EAS
        region%grid%lambda_M     = C%lambda_M_EAS
        region%grid%phi_M        = C%phi_M_EAS
        region%grid%beta_stereo  = C%beta_stereo_EAS
      ELSEIF (region%name == 'GRL') THEN
        xmin                     = C%xmin_GRL
        xmax                     = C%xmax_GRL
        ymin                     = C%ymin_GRL
        ymax                     = C%ymax_GRL
        region%grid%dx           = C%dx_GRL
        region%grid%lambda_M     = C%lambda_M_GRL
        region%grid%phi_M        = C%phi_M_GRL
        region%grid%beta_stereo  = C%beta_stereo_GRL
      ELSEIF (region%name == 'ANT') THEN
        xmin                     = C%xmin_ANT
        xmax                     = C%xmax_ANT
        ymin                     = C%ymin_ANT
        ymax                     = C%ymax_ANT
        region%grid%dx           = C%dx_ANT
        region%grid%lambda_M     = C%lambda_M_ANT
        region%grid%phi_M        = C%phi_M_ANT
        region%grid%beta_stereo  = C%beta_stereo_ANT
      END IF

      ! Determine the number of grid cells we can fit in this domain
      xmid = (xmax + xmin) / 2._dp
      ymid = (ymax + ymin) / 2._dp
      nsx = FLOOR( (xmax - xmid) / region%grid%dx)
      nsy = FLOOR( (ymax - ymid) / region%grid%dx)

      ! Small exceptions for very weird benchmark experiments
      IF (C%choice_refgeo_init_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'SSA_icestream') nsx = 3
      IF (C%choice_refgeo_init_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'ISMIP_HOM_E')   nsx = 25

      region%grid%nx = 1 + 2*nsx
      region%grid%ny = 1 + 2*nsy

    END IF ! IF (par%master) THEN
    CALL sync

    ! Assign range to each processor
    CALL partition_list( region%grid%nx, par%i, par%n, region%grid%i1, region%grid%i2)
    CALL partition_list( region%grid%ny, par%i, par%n, region%grid%j1, region%grid%j2)

    ! Allocate shared memory for x and y
    CALL allocate_shared_dp_1D( region%grid%nx, region%grid%x, region%grid%wx)
    CALL allocate_shared_dp_1D( region%grid%ny, region%grid%y, region%grid%wy)

    ! Fill in x and y
    IF (par%master) THEN

      ! x
      region%grid%xmin = xmid - nsx * region%grid%dx
      region%grid%xmax = xmid + nsx * region%grid%dx
      DO i = 1, region%grid%nx
        region%grid%x( i) = region%grid%xmin + (i-1)*region%grid%dx
      END DO

      ! y
      region%grid%ymin = ymid - nsy * region%grid%dx
      region%grid%ymax = ymid + nsy * region%grid%dx
      DO j = 1, region%grid%ny
        region%grid%y( j) = region%grid%ymin + (j-1)*region%grid%dx
      END DO

    END IF ! IF (par%master) THEN
    CALL sync

    ! Lat,lon coordinates
    CALL allocate_shared_dp_2D( region%grid%ny, region%grid%nx, region%grid%lat, region%grid%wlat)
    CALL allocate_shared_dp_2D( region%grid%ny, region%grid%nx, region%grid%lon, region%grid%wlon)

    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
      CALL inverse_oblique_sg_projection( region%grid%x( i), region%grid%y( j), region%grid%lambda_M, region%grid%phi_M, region%grid%beta_stereo, region%grid%lon( j,i), region%grid%lat( j,i))
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_model_grid
  SUBROUTINE initialise_GIA_model_grid( region)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_GIA_model_grid'
    INTEGER                                            :: nsx, nsy, i, j
    REAL(dp)                                           :: xmin, xmax, ymin, ymax, xmid, ymid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no GIA is used, no need to even create a grid for it
    IF (C%choice_GIA_model == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Assign dummy values to suppress compiler warnings
    xmin = 0._dp
    xmax = 0._dp
    ymin = 0._dp
    ymax = 0._dp
    xmid = 0._dp
    ymid = 0._dp
    nsx  = 0
    nsy  = 0

    ! Allocate shared memory
    CALL allocate_shared_int_0D( region%grid_GIA%nx,           region%grid_GIA%wnx          )
    CALL allocate_shared_int_0D( region%grid_GIA%ny,           region%grid_GIA%wny          )
    CALL allocate_shared_dp_0D(  region%grid_GIA%dx,           region%grid_GIA%wdx          )
    CALL allocate_shared_dp_0D(  region%grid_GIA%xmin,         region%grid_GIA%wxmin        )
    CALL allocate_shared_dp_0D(  region%grid_GIA%xmax,         region%grid_GIA%wxmax        )
    CALL allocate_shared_dp_0D(  region%grid_GIA%ymin,         region%grid_GIA%wymin        )
    CALL allocate_shared_dp_0D(  region%grid_GIA%ymax,         region%grid_GIA%wymax        )
    CALL allocate_shared_dp_0D(  region%grid_GIA%lambda_M,     region%grid_GIA%wlambda_M    )
    CALL allocate_shared_dp_0D(  region%grid_GIA%phi_M,        region%grid_GIA%wphi_M       )
    CALL allocate_shared_dp_0D(  region%grid_GIA%beta_stereo,  region%grid_GIA%wbeta_stereo )

    ! Let the Master do the work
    IF (par%master) THEN

      ! Resolution, domain size, and projection parameters for this region are determined from the config
      IF     (region%name == 'NAM') THEN
        xmin                         = C%xmin_NAM
        xmax                         = C%xmax_NAM
        ymin                         = C%ymin_NAM
        ymax                         = C%ymax_NAM
        region%grid_GIA%dx           = C%dx_GIA
        region%grid_GIA%lambda_M     = C%lambda_M_NAM
        region%grid_GIA%phi_M        = C%phi_M_NAM
        region%grid_GIA%beta_stereo  = C%beta_stereo_NAM
      ELSEIF (region%name == 'EAS') THEN
        xmin                         = C%xmin_EAS
        xmax                         = C%xmax_EAS
        ymin                         = C%ymin_EAS
        ymax                         = C%ymax_EAS
        region%grid_GIA%dx           = C%dx_GIA
        region%grid_GIA%lambda_M     = C%lambda_M_EAS
        region%grid_GIA%phi_M        = C%phi_M_EAS
        region%grid_GIA%beta_stereo  = C%beta_stereo_EAS
      ELSEIF (region%name == 'GRL') THEN
        xmin                         = C%xmin_GRL
        xmax                         = C%xmax_GRL
        ymin                         = C%ymin_GRL
        ymax                         = C%ymax_GRL
        region%grid_GIA%dx           = C%dx_GIA
        region%grid_GIA%lambda_M     = C%lambda_M_GRL
        region%grid_GIA%phi_M        = C%phi_M_GRL
        region%grid_GIA%beta_stereo  = C%beta_stereo_GRL
      ELSEIF (region%name == 'ANT') THEN
        xmin                         = C%xmin_ANT
        xmax                         = C%xmax_ANT
        ymin                         = C%ymin_ANT
        ymax                         = C%ymax_ANT
        region%grid_GIA%dx           = C%dx_GIA
        region%grid_GIA%lambda_M     = C%lambda_M_ANT
        region%grid_GIA%phi_M        = C%phi_M_ANT
        region%grid_GIA%beta_stereo  = C%beta_stereo_ANT
      END IF

      ! Determine the number of grid cells we can fit in this domain
      xmid = (xmax + xmin) / 2._dp
      ymid = (ymax + ymin) / 2._dp
      nsx = FLOOR( (xmax - xmid) / region%grid_GIA%dx) - 1
      nsy = FLOOR( (ymax - ymid) / region%grid_GIA%dx) - 1

      region%grid_GIA%nx = 1 + 2*nsx
      region%grid_GIA%ny = 1 + 2*nsy

    END IF ! IF (par%master) THEN
    CALL sync

    ! Assign range to each processor
    CALL partition_list( region%grid_GIA%nx, par%i, par%n, region%grid_GIA%i1, region%grid_GIA%i2)
    CALL partition_list( region%grid_GIA%ny, par%i, par%n, region%grid_GIA%j1, region%grid_GIA%j2)

    ! Allocate shared memory for x and y
    CALL allocate_shared_dp_1D( region%grid_GIA%nx, region%grid_GIA%x, region%grid_GIA%wx)
    CALL allocate_shared_dp_1D( region%grid_GIA%ny, region%grid_GIA%y, region%grid_GIA%wy)

    ! Fill in x and y
    IF (par%master) THEN

      ! x
      region%grid_GIA%xmin = xmid - nsx * region%grid_GIA%dx
      region%grid_GIA%xmax = xmid + nsx * region%grid_GIA%dx
      DO i = 1, region%grid_GIA%nx
        region%grid_GIA%x( i) = region%grid_GIA%xmin + (i-1)*region%grid_GIA%dx
      END DO

      ! y
      region%grid_GIA%ymin = ymid - nsy * region%grid_GIA%dx
      region%grid_GIA%ymax = ymid + nsy * region%grid_GIA%dx
      DO j = 1, region%grid_GIA%ny
        region%grid_GIA%y( j) = region%grid_GIA%ymin + (j-1)*region%grid_GIA%dx
      END DO

    END IF ! IF (par%master) THEN
    CALL sync

    ! Lat,lon coordinates
    CALL allocate_shared_dp_2D( region%grid_GIA%ny, region%grid_GIA%nx, region%grid_GIA%lat, region%grid_GIA%wlat)
    CALL allocate_shared_dp_2D( region%grid_GIA%ny, region%grid_GIA%nx, region%grid_GIA%lon, region%grid_GIA%wlon)

    DO i = region%grid_GIA%i1, region%grid_GIA%i2
    DO j = 1, region%grid_GIA%ny
      CALL inverse_oblique_sg_projection( region%grid_GIA%x( i), region%grid_GIA%y( j), region%grid_GIA%lambda_M, &
        region%grid_GIA%phi_M, region%grid_GIA%beta_stereo, region%grid_GIA%lon( j,i), region%grid_GIA%lat( j,i))
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_GIA_model_grid

  ! Calculate this region's ice sheet's volume and area
  SUBROUTINE calculate_icesheet_volume_and_area( region)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)          :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calculate_icesheet_volume_and_area'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: ice_area, ice_volume, thickness_above_flotation, ice_volume_above_flotation

    ! Add routine to path
    CALL init_routine( routine_name)

    ice_area                   = 0._dp
    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp

    ! Calculate ice area and volume for processor domain
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny

      IF (region%ice%mask_ice_a( j,i) == 1) THEN
        ice_volume = ice_volume + (region%ice%Hi_a(j,i) * region%grid%dx * region%grid%dx * ice_density / (seawater_density * ocean_area))
        ice_area   = ice_area   + region%grid%dx * region%grid%dx * 1.0E-06_dp ! [km^3]

        ! Thickness above flotation
        thickness_above_flotation = MAX(0._dp, region%ice%Hi_a( j,i) - MAX(0._dp, (region%ice%SL_a( j,i) - region%ice%Hb_a( j,i) * (seawater_density / ice_density))))

        ice_volume_above_flotation = ice_volume_above_flotation + thickness_above_flotation * region%grid%dx * region%grid%dx * ice_density / (seawater_density * ocean_area)
      END IF

    END DO
    END DO
    CALL sync

    CALL MPI_REDUCE( ice_area,                   region%ice_area,                   1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume,                 region%ice_volume,                 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume_above_flotation, region%ice_volume_above_flotation, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    ! Calculate GMSL contribution
    IF (par%master) region%GMSL_contribution = -1._dp * (region%ice_volume_above_flotation - region%ice_volume_above_flotation_PD)
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_icesheet_volume_and_area
  SUBROUTINE calculate_PD_sealevel_contribution( region)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)          :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calculate_PD_sealevel_contribution'
    INTEGER                                            :: i, j
    REAL(dp)                                           :: ice_volume, thickness_above_flotation, ice_volume_above_flotation

    ! Add routine to path
    CALL init_routine( routine_name)

    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp

    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny

      ! Thickness above flotation
      IF (region%refgeo_PD%Hi( j,i) > 0._dp) THEN
        thickness_above_flotation = MAX(0._dp, region%refgeo_PD%Hi( j,i) - MAX(0._dp, (0._dp - region%refgeo_PD%Hb( j,i)) * (seawater_density / ice_density)))
      ELSE
        thickness_above_flotation = 0._dp
      END IF

      ! Ice volume (above flotation) in m.s.l.e
      ice_volume                 = ice_volume                 + region%refgeo_PD%Hi( j,i) * region%grid%dx * region%grid%dx * ice_density / (seawater_density * ocean_area)
      ice_volume_above_flotation = ice_volume_above_flotation + thickness_above_flotation * region%grid%dx * region%grid%dx * ice_density / (seawater_density * ocean_area)

    END DO
    END DO

    CALL MPI_REDUCE( ice_volume                , region%ice_volume_PD,                 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume_above_flotation, region%ice_volume_above_flotation_PD, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_PD_sealevel_contribution

END MODULE IMAU_ICE_main_model
