MODULE IMAU_ICE_main_model

  ! Contains all the routines for initialising and running the IMAU-ICE regional ice-sheet model.

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             allocate_shared_bool_0D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_model_region, type_ice_model, type_reference_geometry, &
                                             type_SMB_model, type_BMB_model, type_forcing_data, type_grid, &
                                             type_climate_matrix_global, type_ocean_matrix_global
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             inverse_oblique_sg_projection, surface_elevation
  USE parameters_module,               ONLY: seawater_density, ice_density, T0
  USE reference_fields_module,         ONLY: initialise_reference_geometries
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, initialise_debug_fields, create_debug_file, associate_debug_fields, &
                                             create_restart_file, create_help_fields_file, write_to_restart_file, write_to_help_fields_file, &
                                             create_regional_scalar_output_file
  USE forcing_module,                  ONLY: forcing, initialise_geothermal_heat_flux_regional
  USE general_ice_model_data_module,   ONLY: update_general_ice_model_data, initialise_basins, initialise_mask_noice
  USE ice_velocity_module,             ONLY: solve_DIVA
  USE ice_dynamics_module,             ONLY: initialise_ice_model,              run_ice_model
  USE thermodynamics_module,           ONLY: initialise_ice_temperature,        run_thermo_model, calc_ice_rheology
  USE ocean_module,                    ONLY: initialise_ocean_model_regional,   run_ocean_model
  USE climate_module,                  ONLY: initialise_climate_model_regional, run_climate_model
  USE SMB_module,                      ONLY: initialise_SMB_model,              run_SMB_model
  USE BMB_module,                      ONLY: initialise_BMB_model,              run_BMB_model
  USE isotopes_module,                 ONLY: initialise_isotopes_model,         run_isotopes_model
  USE bedrock_ELRA_module,             ONLY: initialise_ELRA_model,             run_ELRA_model
  USE SELEN_main_module,               ONLY: apply_SELEN_bed_geoid_deformation_rates
  USE calving_module,                  ONLY: apply_calving_law
  USE scalar_data_output_module,       ONLY: write_regional_scalar_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE run_model( region, climate_matrix_global, t_end)
    ! Run the model until t_end (usually a 100 years further)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_climate_matrix_global),    INTENT(INOUT) :: climate_matrix_global
    REAL(dp),                            INTENT(IN)    :: t_end
    
    ! Local variables:
    REAL(dp)                                      :: tstart, tstop, t1, t2
    INTEGER                                       :: it
    
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,'(A,A3,A,A13,A,F9.3,A,F9.3,A)') '  Running model region ', region%name, ' (', TRIM(region%long_name), & 
                                                            ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' kyr'
    
    ! Set the intermediary pointers in "debug" to this region's debug data fields
    CALL associate_debug_fields(  region)
               
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
    DO WHILE (region%time < t_end)
      it = it + 1
      
    ! GIA
    ! ===
    
      t1 = MPI_WTIME()
      IF     (C%choice_GIA_model == 'none') THEN
        ! Nothing to be done
      ELSEIF (C%choice_GIA_model == 'ELRA') THEN
        CALL run_ELRA_model( region)
      ELSEIF (C%choice_GIA_model == 'SELEN') THEN
        CALL apply_SELEN_bed_geoid_deformation_rates( region)
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR - choice_GIA_model "', C%choice_GIA_model, '" not implemented in run_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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
      
    ! Climate , SMB and BMB
    ! =====================
    
      t1 = MPI_WTIME()
            
      ! Run the climate model
      IF (region%do_climate) THEN
        CALL run_climate_model( region, climate_matrix_global, region%time)
      END IF
      
      ! Run the ocean model
      IF (region%do_ocean) THEN
        CALL run_ocean_model( region%grid, region%ice, region%ocean_matrix, region%climate_matrix, region%name, region%time)
      END IF     
    
      ! Run the SMB model
      IF (region%do_SMB) THEN
        CALL run_SMB_model( region%grid, region%ice, region%climate_matrix, region%time, region%SMB, region%mask_noice)
      END IF
    
      ! Run the BMB model
      IF (region%do_BMB) THEN
        CALL run_BMB_model( region%grid, region%ice, region%ocean_matrix%applied, region%BMB, region%name, region%time)
      END IF
      
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_climate = region%tcomp_climate + t2 - t1
      
    ! Thermodynamics
    ! ==============
    
      t1 = MPI_WTIME()
      CALL run_thermo_model( region%grid, region%ice, region%climate_matrix%applied, region%ocean_matrix%applied, region%SMB, region%time, do_solve_heat_equation = region%do_thermo)
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_thermo = region%tcomp_thermo + t2 - t1
      
    ! Isotopes
    ! ========
    
      CALL run_isotopes_model( region)
      
    ! Time step and output
    ! ====================
                            
      ! Write output
      IF (region%do_output) THEN
        IF (par%master) CALL write_to_restart_file(     region, forcing)
        IF (par%master) CALL write_to_help_fields_file( region)
        CALL sync
      END IF
      
      ! Update ice geometry and advance region time
      t1 = MPI_WTIME()
      region%ice%Hi_a( :,region%grid%i1:region%grid%i2) = region%ice%Hi_tplusdt_a( :,region%grid%i1:region%grid%i2)
      CALL sync
      
      ! Save the previous ice mask, for use in thermodynamics
      region%ice%mask_ice_a_prev( :,region%grid%i1:region%grid%i2) = region%ice%mask_ice_a( :,region%grid%i1:region%grid%i2)
      CALL sync
      
      ! Update masks, apply calving
      CALL update_general_ice_model_data( region%grid, region%ice)
      CALL apply_calving_law( region%grid, region%ice, region%refgeo_PD)
      CALL update_general_ice_model_data( region%grid, region%ice)
      IF (par%master) region%time = region%time + region%dt
      CALL sync
      t2 = MPI_WTIME()
      IF (par%master) region%tcomp_ice = region%tcomp_ice + t2 - t1
      
      ! DENK DROM
      !region%time = t_end
    
    END DO ! DO WHILE (region%time < t_end)
                            
  ! ===========================================
  ! ===== End of the main model time loop =====
  ! ===========================================
    
    ! Write to NetCDF output one last time at the end of the simulation
    IF (region%time == C%end_time_of_run) THEN
      IF (par%master) CALL write_to_restart_file(     region, forcing)
      IF (par%master) CALL write_to_help_fields_file( region)
    END IF 
    
    ! Determine total ice sheet area, volume, volume-above-flotation and GMSL contribution,
    ! used for writing to text output and in the inverse routine
    CALL calculate_icesheet_volume_and_area(region)
    
    tstop = MPI_WTIME()
    region%tcomp_total = tstop - tstart
    
    ! Write to text output
    CALL write_regional_scalar_data( region, region%time)
    
  END SUBROUTINE run_model
  
  ! Initialise the entire model region - read initial and PD data, initialise the ice dynamics, climate and SMB sub models
  SUBROUTINE initialise_model( region, name, climate_matrix_global, ocean_matrix_global)
  
    USE climate_module, ONLY: map_glob_to_grid_2D
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),          INTENT(INOUT)     :: region
    CHARACTER(LEN=3),                 INTENT(IN)        :: name
    TYPE(type_climate_matrix_global), INTENT(INOUT)     :: climate_matrix_global
    TYPE(type_ocean_matrix_global),   INTENT(INOUT)     :: ocean_matrix_global
    
    ! Local variables:
    CHARACTER(LEN=20)                                   :: short_filename
    INTEGER                                             :: n
    
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
    
    ! ===== Set up debug fields and output files ======
    ! =================================================
    
    ! Debug fields
    CALL initialise_debug_fields( region)
    IF (par%master .AND. C%do_write_debug_data) CALL create_debug_file( region)
    CALL associate_debug_fields(  region)
    
    ! Restart file
    short_filename = 'restart_NAM.nc'
    short_filename(9:11) = region%name
    DO n = 1, 256
      region%restart%netcdf%filename(n:n) = ' '
    END DO
    region%restart%netcdf%filename = TRIM(C%output_dir)//TRIM(short_filename)
    
    ! Help fields file
    short_filename = 'help_fields_NAM.nc'
    short_filename(13:15) = region%name
    DO n = 1, 256
      region%help_fields%filename(n:n) = ' '
    END DO
    region%help_fields%filename = TRIM(C%output_dir)//TRIM(short_filename)

    ! Let the Master create the (empty) NetCDF files
    IF (par%master) CALL create_restart_file(     region, forcing)
    IF (par%master) CALL create_help_fields_file( region)
    
    ! ===== Initialise initial, present-day, and GIA equilibrium reference geometries =====
    ! =====================================================================================
    
    CALL initialise_reference_geometries( region%grid, region%refgeo_init, region%refgeo_PD, region%refgeo_GIAeq, region%name)
    
    ! ===== Define mask where no ice is allowed to form (i.e. Greenland in NAM and EAS, Ellesmere Island in GRL)
    ! ==========================================================================================================
    
    CALL initialise_mask_noice( region)
  
    ! ===== The ice dynamics model
    ! ============================
    
    CALL initialise_ice_model( region%grid, region%ice, region%refgeo_init, region%mask_noice)
    
    ! ===== Define ice basins =====
    ! =============================
    
    ! Allocate shared memory
    CALL allocate_shared_int_2D( region%grid%ny, region%grid%nx, region%ice%basin_ID, region%ice%wbasin_ID)
    CALL allocate_shared_int_0D(                                 region%ice%nbasins,  region%ice%wnbasins )
    
    ! Define basins
    CALL initialise_basins( region%grid, region%ice%basin_ID, region%ice%nbasins, region%name)
        
    ! ===== The climate model =====
    ! =============================
    
    CALL initialise_climate_model_regional( region, climate_matrix_global)
    
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
    ELSEIF (C%choice_GIA_model == 'SELEN') THEN
      CALL initialise_GIA_model_grid( region)
    ELSE
      WRITE(0,*) '  ERROR - choice_GIA_model "', C%choice_GIA_model, '" not implemented in initialise_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! ===== The isotopes model =====
    ! ==============================
   
    CALL initialise_isotopes_model( region)
    
    ! ===== Geothermal heat flux =====
    ! ================================
    
    CALL initialise_geothermal_heat_flux_regional( region%grid, region%ice)
    
    ! ===== Initialise the ice temperature profile =====
    ! ==================================================
    
    ! Run the climate and SMB models once, to get the correct surface temperature+SMB fields for the ice temperature initialisation
    CALL run_climate_model( region, climate_matrix_global, C%start_time_of_run)
    CALL run_SMB_model( region%grid, region%ice, region%climate_matrix, C%start_time_of_run, region%SMB, region%mask_noice)
    
    ! Initialise the temperature field
    CALL initialise_ice_temperature( region%grid, region%ice, region%climate_matrix%applied, region%ocean_matrix%applied, region%SMB, region%name)
    
    ! Initialise the rheology
    CALL calc_ice_rheology( region%grid, region%ice, C%start_time_of_run)
    
    ! ===== Scalar output (regionally integrated ice volume, SMB components, etc.)
    ! ============================================================================
    
    ! Create the file
    CALL create_regional_scalar_output_file( region)
    
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
    
  END SUBROUTINE initialise_model
  SUBROUTINE allocate_region_timers_and_scalars( region)
    ! Allocate shared memory for this region's timers (used for the asynchronous coupling between the
    ! ice dynamics and the secondary model components), and for the scalars (integrated ice volume and
    ! area, SMB components, computation times, etc.)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    
    ! Timers and time steps
    ! =====================
    
    CALL allocate_shared_dp_0D(   region%time,             region%wtime            )
    CALL allocate_shared_dp_0D(   region%dt,               region%wdt              )
    CALL allocate_shared_dp_0D(   region%dt_prev,          region%wdt_prev         )
    CALL allocate_shared_dp_0D(   region%dt_crit_SIA,      region%wdt_crit_SIA     )
    CALL allocate_shared_dp_0D(   region%dt_crit_SSA,      region%wdt_crit_SSA     )
    CALL allocate_shared_dp_0D(   region%dt_crit_ice,      region%wdt_crit_ice     )
    CALL allocate_shared_dp_0D(   region%dt_crit_ice_prev, region%wdt_crit_ice_prev)
    
    CALL allocate_shared_dp_0D(   region%t_last_SIA,       region%wt_last_SIA      )
    CALL allocate_shared_dp_0D(   region%t_next_SIA,       region%wt_next_SIA      )
    CALL allocate_shared_bool_0D( region%do_SIA,           region%wdo_SIA          )
    
    CALL allocate_shared_dp_0D(   region%t_last_SSA,       region%wt_last_SSA      )
    CALL allocate_shared_dp_0D(   region%t_next_SSA,       region%wt_next_SSA      )
    CALL allocate_shared_bool_0D( region%do_SSA,           region%wdo_SSA          )
    
    CALL allocate_shared_dp_0D(   region%t_last_DIVA,      region%wt_last_DIVA     )
    CALL allocate_shared_dp_0D(   region%t_next_DIVA,      region%wt_next_DIVA     )
    CALL allocate_shared_bool_0D( region%do_DIVA,          region%wdo_DIVA         )
    
    CALL allocate_shared_dp_0D(   region%t_last_thermo,    region%wt_last_thermo   )
    CALL allocate_shared_dp_0D(   region%t_next_thermo,    region%wt_next_thermo   )
    CALL allocate_shared_bool_0D( region%do_thermo,        region%wdo_thermo       )
    
    CALL allocate_shared_dp_0D(   region%t_last_output,    region%wt_last_output   )
    CALL allocate_shared_dp_0D(   region%t_next_output,    region%wt_next_output   )
    CALL allocate_shared_bool_0D( region%do_output,        region%wdo_output       )
    
    CALL allocate_shared_dp_0D(   region%t_last_climate,   region%wt_last_climate  )
    CALL allocate_shared_dp_0D(   region%t_next_climate,   region%wt_next_climate  )
    CALL allocate_shared_bool_0D( region%do_climate,       region%wdo_climate      )
    
    CALL allocate_shared_dp_0D(   region%t_last_ocean,     region%wt_last_ocean    )
    CALL allocate_shared_dp_0D(   region%t_next_ocean,     region%wt_next_ocean    )
    CALL allocate_shared_bool_0D( region%do_ocean,         region%wdo_ocean        )
    
    CALL allocate_shared_dp_0D(   region%t_last_SMB,       region%wt_last_SMB      )
    CALL allocate_shared_dp_0D(   region%t_next_SMB,       region%wt_next_SMB      )
    CALL allocate_shared_bool_0D( region%do_SMB,           region%wdo_SMB          )
    
    CALL allocate_shared_dp_0D(   region%t_last_BMB,       region%wt_last_BMB      )
    CALL allocate_shared_dp_0D(   region%t_next_BMB,       region%wt_next_BMB      )
    CALL allocate_shared_bool_0D( region%do_BMB,           region%wdo_BMB          )
    
    CALL allocate_shared_dp_0D(   region%t_last_ELRA,      region%wt_last_ELRA     )
    CALL allocate_shared_dp_0D(   region%t_next_ELRA,      region%wt_next_ELRA     )
    CALL allocate_shared_bool_0D( region%do_ELRA,          region%wdo_ELRA         ) 
    
    IF (par%master) THEN
      region%time           = C%start_time_of_run
      region%dt             = C%dt_min
      region%dt_prev        = C%dt_min
      
      region%t_last_SIA     = C%start_time_of_run
      region%t_next_SIA     = C%start_time_of_run
      region%do_SIA         = .TRUE.
      
      region%t_last_SSA     = C%start_time_of_run
      region%t_next_SSA     = C%start_time_of_run
      region%do_SSA         = .TRUE.
      
      region%t_last_DIVA    = C%start_time_of_run
      region%t_next_DIVA    = C%start_time_of_run
      region%do_DIVA        = .TRUE.
      
      region%t_last_thermo  = C%start_time_of_run
      region%t_next_thermo  = C%start_time_of_run + C%dt_thermo
      region%do_thermo      = .FALSE.
      
      region%t_last_climate = C%start_time_of_run
      region%t_next_climate = C%start_time_of_run
      region%do_climate     = .TRUE.
      
      region%t_last_ocean   = C%start_time_of_run
      region%t_next_ocean   = C%start_time_of_run
      region%do_ocean       = .TRUE.
      
      region%t_last_SMB     = C%start_time_of_run
      region%t_next_SMB     = C%start_time_of_run
      region%do_SMB         = .TRUE.
      
      region%t_last_BMB     = C%start_time_of_run
      region%t_next_BMB     = C%start_time_of_run
      region%do_BMB         = .TRUE.
      
      region%t_last_ELRA    = C%start_time_of_run
      region%t_next_ELRA    = C%start_time_of_run
      region%do_ELRA        = .TRUE.
      
      region%t_last_output  = C%start_time_of_run
      region%t_next_output  = C%start_time_of_run
      region%do_output      = .TRUE.
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
    
  END SUBROUTINE allocate_region_timers_and_scalars
  SUBROUTINE initialise_model_grid( region)
    ! Initialise the square grid for this model region
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    
    ! Local variables:
    INTEGER                                            :: nsx, nsy, i, j
    REAL(dp)                                           :: xmin, xmax, ymin, ymax, xmid, ymid
    
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
    CALL allocate_shared_dp_0D(  region%grid%alpha_stereo, region%grid%walpha_stereo)
    
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
        region%grid%alpha_stereo = C%alpha_stereo_NAM
      ELSEIF (region%name == 'EAS') THEN
        xmin                     = C%xmin_EAS
        xmax                     = C%xmax_EAS
        ymin                     = C%ymin_EAS
        ymax                     = C%ymax_EAS
        region%grid%dx           = C%dx_EAS
        region%grid%lambda_M     = C%lambda_M_EAS
        region%grid%phi_M        = C%phi_M_EAS
        region%grid%alpha_stereo = C%alpha_stereo_EAS
      ELSEIF (region%name == 'GRL') THEN
        xmin                     = C%xmin_GRL
        xmax                     = C%xmax_GRL
        ymin                     = C%ymin_GRL
        ymax                     = C%ymax_GRL
        region%grid%dx           = C%dx_GRL
        region%grid%lambda_M     = C%lambda_M_GRL
        region%grid%phi_M        = C%phi_M_GRL
        region%grid%alpha_stereo = C%alpha_stereo_GRL
      ELSEIF (region%name == 'ANT') THEN
        xmin                     = C%xmin_ANT
        xmax                     = C%xmax_ANT
        ymin                     = C%ymin_ANT
        ymax                     = C%ymax_ANT
        region%grid%dx           = C%dx_ANT
        region%grid%lambda_M     = C%lambda_M_ANT
        region%grid%phi_M        = C%phi_M_ANT
        region%grid%alpha_stereo = C%alpha_stereo_ANT
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
      DO i = 1, region%grid%nx
        region%grid%x( i) = -nsx*region%grid%dx + (i-1)*region%grid%dx
      END DO
      DO j = 1, region%grid%ny
        region%grid%y( j) = -nsy*region%grid%dx + (j-1)*region%grid%dx
      END DO
      
      region%grid%xmin = MINVAL(region%grid%x)
      region%grid%xmax = MAXVAL(region%grid%x)
      region%grid%ymin = MINVAL(region%grid%y)
      region%grid%ymax = MAXVAL(region%grid%y)
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Lat,lon coordinates
    CALL allocate_shared_dp_2D( region%grid%ny, region%grid%nx, region%grid%lat, region%grid%wlat)
    CALL allocate_shared_dp_2D( region%grid%ny, region%grid%nx, region%grid%lon, region%grid%wlon)
    
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
      CALL inverse_oblique_sg_projection( region%grid%x( i), region%grid%y( j), region%grid%lambda_M, region%grid%phi_M, region%grid%alpha_stereo, region%grid%lon( j,i), region%grid%lat( j,i))
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE initialise_model_grid
  SUBROUTINE initialise_GIA_model_grid( region)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    
    ! Local variables:
    INTEGER                                            :: nsx, nsy, i, j
    REAL(dp)                                           :: xmin, xmax, ymin, ymax, xmid, ymid
    
    ! If no GIA is used, no need to even create a grid for it
    IF (C%choice_GIA_model == 'none') RETURN
    
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
    CALL allocate_shared_dp_0D(  region%grid_GIA%alpha_stereo, region%grid_GIA%walpha_stereo)
    
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
        region%grid_GIA%alpha_stereo = C%alpha_stereo_NAM
      ELSEIF (region%name == 'EAS') THEN
        xmin                         = C%xmin_EAS
        xmax                         = C%xmax_EAS
        ymin                         = C%ymin_EAS
        ymax                         = C%ymax_EAS
        region%grid_GIA%dx           = C%dx_GIA
        region%grid_GIA%lambda_M     = C%lambda_M_EAS
        region%grid_GIA%phi_M        = C%phi_M_EAS
        region%grid_GIA%alpha_stereo = C%alpha_stereo_EAS
      ELSEIF (region%name == 'GRL') THEN
        xmin                         = C%xmin_GRL
        xmax                         = C%xmax_GRL
        ymin                         = C%ymin_GRL
        ymax                         = C%ymax_GRL
        region%grid_GIA%dx           = C%dx_GIA
        region%grid_GIA%lambda_M     = C%lambda_M_GRL
        region%grid_GIA%phi_M        = C%phi_M_GRL
        region%grid_GIA%alpha_stereo = C%alpha_stereo_GRL
      ELSEIF (region%name == 'ANT') THEN
        xmin                         = C%xmin_ANT
        xmax                         = C%xmax_ANT
        ymin                         = C%ymin_ANT
        ymax                         = C%ymax_ANT
        region%grid_GIA%dx           = C%dx_GIA
        region%grid_GIA%lambda_M     = C%lambda_M_ANT
        region%grid_GIA%phi_M        = C%phi_M_ANT
        region%grid_GIA%alpha_stereo = C%alpha_stereo_ANT
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
      DO i = 1, region%grid_GIA%nx
        region%grid_GIA%x( i) = -nsx*region%grid_GIA%dx + (i-1)*region%grid_GIA%dx
      END DO
      DO j = 1, region%grid_GIA%ny
        region%grid_GIA%y( j) = -nsy*region%grid_GIA%dx + (j-1)*region%grid_GIA%dx
      END DO
      
      region%grid_GIA%xmin = MINVAL(region%grid_GIA%x)
      region%grid_GIA%xmax = MAXVAL(region%grid_GIA%x)
      region%grid_GIA%ymin = MINVAL(region%grid_GIA%y)
      region%grid_GIA%ymax = MAXVAL(region%grid_GIA%y)
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Lat,lon coordinates
    CALL allocate_shared_dp_2D( region%grid_GIA%ny, region%grid_GIA%nx, region%grid_GIA%lat, region%grid_GIA%wlat)
    CALL allocate_shared_dp_2D( region%grid_GIA%ny, region%grid_GIA%nx, region%grid_GIA%lon, region%grid_GIA%wlon)
    
    DO i = region%grid_GIA%i1, region%grid_GIA%i2
    DO j = 1, region%grid_GIA%ny
      CALL inverse_oblique_sg_projection( region%grid_GIA%x( i), region%grid_GIA%y( j), region%grid_GIA%lambda_M, region%grid_GIA%phi_M, region%grid_GIA%alpha_stereo, region%grid_GIA%lon( j,i), region%grid_GIA%lat( j,i))
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE initialise_GIA_model_grid
  
  ! Calculate this region's ice sheet's volume and area
  SUBROUTINE calculate_icesheet_volume_and_area( region)
    
    USE parameters_module,           ONLY: ocean_area, seawater_density, ice_density
  
    IMPLICIT NONE  
    
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    
    INTEGER                                       :: i,j
    REAL(dp)                                      :: ice_area, ice_volume, thickness_above_flotation, ice_volume_above_flotation
    
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
    
  END SUBROUTINE calculate_icesheet_volume_and_area
  SUBROUTINE calculate_PD_sealevel_contribution( region)
    
    USE parameters_module,           ONLY: ocean_area, seawater_density, ice_density
  
    IMPLICIT NONE  
    
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    
    INTEGER                                       :: i, j
    REAL(dp)                                      :: ice_volume, thickness_above_flotation, ice_volume_above_flotation
    
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
    
  END SUBROUTINE calculate_PD_sealevel_contribution

END MODULE IMAU_ICE_main_model
