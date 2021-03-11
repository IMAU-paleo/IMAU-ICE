MODULE IMAU_ICE_main_model

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parallel_module,                 ONLY: par, sync, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_model_region, type_ice_model, type_PD_data_fields, type_init_data_fields, &
                                             type_climate_model, type_climate_matrix, type_SMB_model, type_BMB_model, type_forcing_data, type_grid
  USE utilities_module,                ONLY: inverse_oblique_sg_projection
  USE parameters_module,               ONLY: seawater_density, ice_density, T0
  USE reference_fields_module,         ONLY: initialise_PD_data_fields, initialise_init_data_fields, map_PD_data_to_model_grid, map_init_data_to_model_grid
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, initialise_debug_fields, create_debug_file, associate_debug_fields, &
                                             create_restart_file, create_help_fields_file, write_to_restart_file, write_to_help_fields_file
  USE forcing_module,                  ONLY: forcing
  USE general_ice_model_data_module,   ONLY: update_general_ice_model_data
  USE ice_dynamics_module,             ONLY: initialise_ice_model, calculate_ice_thickness_change, solve_SIA, solve_SSA
  USE thermodynamics_module,           ONLY: initialise_ice_temperature, update_ice_temperature
  USE climate_module,                  ONLY: initialise_climate_model,  run_climate_model
  USE SMB_module,                      ONLY: initialise_SMB_model,      run_SMB_model
  USE BMB_module,                      ONLY: initialise_BMB_model,      run_BMB_model
  USE isotopes_module,                 ONLY: initialise_isotopes_model, run_isotopes_model
  USE bedrock_ELRA_module,             ONLY: initialise_ELRA_model,     run_ELRA_model

  IMPLICIT NONE

CONTAINS

  SUBROUTINE run_model( region, t_end)
    ! Run the model until t_end (usually a 100 years further)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    REAL(dp),                   INTENT(IN)        :: t_end
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    REAL(dp)                                      :: tstart, tstop, t1, t2
    
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,'(A,A3,A,A13,A,F9.3,A,F9.3,A)') '  Running model region ', region%name, ' (', TRIM(region%long_name), & 
                                                            ') from t = ', region%time/1000._dp, ' to t = ', t_end/1000._dp, ' kyr'
    
    ! Set the intermediary pointers in "debug" to this region's debug data fields
    CALL associate_debug_fields(  region)
               
    ! Computation time tracking
    region%tcomp_total          = 0._dp
    region%tcomp_SIA            = 0._dp
    region%tcomp_SSA            = 0._dp
    region%tcomp_thermo         = 0._dp
    region%tcomp_climate        = 0._dp
    
    tstart = MPI_WTIME()

    ! Write to text output at t=0
    IF (region%time == C%start_time_of_run) THEN
      CALL write_text_output(region)
    END IF
                            
  ! ====================================
  ! ===== The main model time loop =====
  ! ====================================
    
    DO WHILE (region%time < t_end)
    
      ! Update bedrock elevation
      IF (C%choice_GIA_model == 'ELRA') THEN
        CALL run_ELRA_model( region)
      ELSE
        WRITE(0,*) '  ERROR - choice_GIA_model "', C%choice_GIA_model, '" not implemented in run_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
          
      ! Update ice thickness based on ice velocity and mass balance
      t1 = MPI_WTIME()
      CALL calculate_ice_thickness_change( region%grid, region%ice, region%SMB, region%BMB, region%dt, region%mask_noice)
                        
      t2 = MPI_WTIME()
      region%tcomp_SIA = region%tcomp_SIA + t2 - t1
      
      ! Update general ice model data - Hs, masks, ice physical properties, on both Aa and Ac mesh
      t1 = MPI_WTIME()
      CALL update_general_ice_model_data( region%grid, region%ice, region%time)
      t2 = MPI_WTIME()
      region%tcomp_SIA = region%tcomp_SIA + t2 - t1
      
    ! Ice dynamics
    ! ============
      
      ! Solve the SIA
      IF (region%do_solve_SIA) THEN ! .OR. region%do_solve_SSA) THEN
        t1 = MPI_WTIME()
        CALL solve_SIA( region%grid, region%ice)
        region%t0_SIA = region%time
        t2 = MPI_WTIME()
        region%tcomp_SIA = region%tcomp_SIA + t2 - t1
      END IF
      
      ! Solve the SSA
      IF (region%do_solve_SSA) THEN ! .OR. region%do_solve_SIA) THEN
        t1 = MPI_WTIME()
        CALL solve_SSA( region%grid, region%ice)
        region%t0_SSA = region%time
        t2 = MPI_WTIME()
        region%tcomp_SSA = region%tcomp_SSA + t2 - t1
      END IF
      
      ! Add velocities together
      region%ice%U_vav_Aa(  :,region%grid%i1:region%grid%i2) = region%ice%U_vav_SIA_Aa(  :,region%grid%i1:region%grid%i2) + region%ice%U_vav_SSA_Aa(  :,region%grid%i1:region%grid%i2)
      region%ice%V_vav_Aa(  :,region%grid%i1:region%grid%i2) = region%ice%V_vav_SIA_Aa(  :,region%grid%i1:region%grid%i2) + region%ice%V_vav_SSA_Aa(  :,region%grid%i1:region%grid%i2)
      region%ice%U_vav_Acx( :,region%grid%i1:MIN(region%grid%nx-1,region%grid%i2)) = &
        region%ice%U_vav_SIA_Acx( :,region%grid%i1:MIN(region%grid%nx-1,region%grid%i2)) + region%ice%U_vav_SSA_Acx( :,region%grid%i1:MIN(region%grid%nx-1,region%grid%i2))
      region%ice%V_vav_Acy( :,region%grid%i1:MIN(region%grid%nx-1,region%grid%i2)) = &
        region%ice%V_vav_SIA_Acy( :,region%grid%i1:MIN(region%grid%nx-1,region%grid%i2)) + region%ice%V_vav_SSA_Acy( :,region%grid%i1:MIN(region%grid%nx-1,region%grid%i2))
      region%ice%U_base_Aa( :,region%grid%i1:region%grid%i2) = region%ice%U_vav_SSA_Aa(  :,region%grid%i1:region%grid%i2)
      region%ice%V_base_Aa( :,region%grid%i1:region%grid%i2) = region%ice%V_vav_SSA_Aa(  :,region%grid%i1:region%grid%i2)
      region%ice%U_surf_Aa( :,region%grid%i1:region%grid%i2) = region%ice%U_3D_Aa(     1,:,region%grid%i1:region%grid%i2)
      region%ice%V_surf_Aa( :,region%grid%i1:region%grid%i2) = region%ice%V_3D_Aa(     1,:,region%grid%i1:region%grid%i2)
      CALL sync
      
    ! Climate , SMB and BMB
    ! =====================
    
      t1 = MPI_WTIME()
            
      ! Run the climate model
      IF (region%do_climate) THEN
        CALL run_climate_model( region%grid, region%ice, region%SMB, region%climate, region%name, region%time)
        region%t0_climate = region%time
      END IF
    
      ! Run the SMB model
      IF (region%do_SMB) THEN
        CALL run_SMB_model( region%grid, region%ice, region%climate%applied, region%time, region%SMB, region%mask_noice)
        region%t0_SMB = region%time
      END IF
    
      ! Run the BMB model
      IF (region%do_BMB) THEN
        CALL run_BMB_model( region%grid, region%ice, region%climate%applied, region%BMB, region%name)
        region%t0_BMB = region%time
      END IF
      
      t2 = MPI_WTIME()
      region%tcomp_climate = region%tcomp_climate + t2 - t1
      
    ! Thermodynamics
    ! ==============
      
      ! Solve thermodynamics
      IF (region%do_thermodynamics) THEN
        t1 = MPI_WTIME()
        CALL update_ice_temperature( region%grid, region%ice, region%climate%applied, region%SMB)
        region%t0_thermo = region%time
        t2 = MPI_WTIME()
        region%tcomp_thermo = region%tcomp_thermo + t2 - t1
      END IF
      
    ! Isotopes
    ! ========
    
      CALL run_isotopes_model( region)
      
    ! Time step and output
    ! ====================
                            
      ! Write output
      IF (region%do_write_output) THEN
        IF (par%master) CALL write_to_restart_file(     region, forcing)
        IF (par%master) CALL write_to_help_fields_file( region)
        region%t0_output = region%time
      END IF
      
      ! Calculate time steps, determine actions and advance region time
      t1 = MPI_WTIME()
      CALL determine_timesteps_and_actions( region, t_end)
      t2 = MPI_WTIME()
      region%tcomp_SIA = region%tcomp_SIA + t2 - t1
      
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
    CALL write_text_output(region)
    
  END SUBROUTINE run_model
  
  ! Initialise the entire model region - read initial and PD data, initialise the ice dynamics, climate and SMB sub models
  SUBROUTINE initialise_model( region, name, matrix)
  
    USE climate_module, ONLY: map_glob_to_grid_2D
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    CHARACTER(LEN=3),                INTENT(IN)        :: name
    TYPE(type_climate_matrix),       INTENT(IN)        :: matrix
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    CHARACTER(LEN=20)                                  :: short_filename
    INTEGER                                            :: n
              
    ! Basic initialisation
    ! ====================
    
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
    
    ! Initialise last update times
    ! ============================
    
    region%time               = C%start_time_of_run
    
    region%t0_SIA             = C%start_time_of_run
    region%t0_SSA             = C%start_time_of_run
    region%t0_thermo          = C%start_time_of_run
    region%t0_climate         = C%start_time_of_run
    region%t0_SMB             = C%start_time_of_run
    region%t0_BMB             = C%start_time_of_run
    region%t0_ELRA            = C%start_time_of_run
    region%t0_output          = C%start_time_of_run
    
    region%t1_SIA             = C%start_time_of_run
    region%t1_SSA             = C%start_time_of_run
    region%t1_thermo          = C%start_time_of_run + C%dt_thermo
    region%t1_climate         = C%start_time_of_run
    region%t1_SMB             = C%start_time_of_run
    region%t1_BMB             = C%start_time_of_run
    region%t1_ELRA            = C%start_time_of_run
    region%t1_output          = C%start_time_of_run
    
    region%dt                 = 0._dp
    region%dt_prev            = 1000._dp
    
    region%dt_SIA             = 0._dp
    region%dt_SSA             = 0._dp
    
    region%do_solve_SIA       = .TRUE.
    region%do_solve_SSA       = .TRUE.
    region%do_thermodynamics  = .FALSE.
    region%do_climate         = .TRUE.
    region%do_SMB             = .TRUE.
    region%do_BMB             = .TRUE.
    region%do_ELRA            = .TRUE.
    region%do_write_output    = .TRUE.
    
    CALL allocate_shared_dp_0D( region%ice_area                     , region%wice_area                     )
    CALL allocate_shared_dp_0D( region%ice_volume                   , region%wice_volume                   )
    CALL allocate_shared_dp_0D( region%ice_volume_PD                , region%wice_volume_PD                )
    CALL allocate_shared_dp_0D( region%ice_volume_above_flotation   , region%wice_volume_above_flotation   )
    CALL allocate_shared_dp_0D( region%ice_volume_above_flotation_PD, region%wice_volume_above_flotation_PD)
    CALL allocate_shared_dp_0D( region%GMSL_contribution            , region%wGMSL_contribution            )
    CALL allocate_shared_dp_0D( region%mean_isotope_content         , region%wmean_isotope_content         )
    CALL allocate_shared_dp_0D( region%mean_isotope_content_PD      , region%wmean_isotope_content_PD      )
    CALL allocate_shared_dp_0D( region%d18O_contribution            , region%wd18O_contribution            )
    CALL allocate_shared_dp_0D( region%d18O_contribution_PD         , region%wd18O_contribution_PD         )
    
    ! ===== PD and init reference data fields =====
    ! =============================================
    
    CALL initialise_PD_data_fields(   region%PD,   region%name)
    CALL initialise_init_data_fields( region%init, region%name)
    
    ! ===== Initialise this region's grid =====
    ! =========================================
    
    CALL initialise_model_grid( region)
    
    ! ===== Initialise dummy fields for debugging
    ! ===========================================
    
    CALL initialise_debug_fields( region)

    ! ===== Map PD and init data to the model grid =====
    ! ==================================================
    
    CALL map_PD_data_to_model_grid(   region%grid, region%PD  )
    CALL map_init_data_to_model_grid( region%grid, region%init)
    
    CALL calculate_PD_sealevel_contribution(region)
    
    ! ===== Define mask where no ice is allowed to form (i.e. Greenland in NAM and EAS, Ellesmere Island in GRL)
    ! ==========================================================================================================
    
    CALL initialise_mask_noice( region)
       
    ! ===== Output files =====
    ! ========================

    ! Set output filenames for this region
    short_filename = 'restart_NAM.nc'
    short_filename(9:11) = region%name
    DO n = 1, 256
      region%restart%filename(n:n) = ' '
    END DO
    region%restart%filename = TRIM(C%output_dir)//TRIM(short_filename)

    short_filename = 'help_fields_NAM.nc'
    short_filename(13:15) = region%name
    DO n = 1, 256
      region%help_fields%filename(n:n) = ' '
    END DO
    region%help_fields%filename = TRIM(C%output_dir)//TRIM(short_filename)

    ! Let the Master create the (empty) NetCDF files
    IF (par%master)                             CALL create_restart_file(     region, forcing)
    IF (par%master)                             CALL create_help_fields_file( region)
    IF (par%master .AND. C%do_write_debug_data) CALL create_debug_file(       region)
    CALL associate_debug_fields(  region)
        
    ! ===== The climate model =====
    ! =============================    
    
    CALL initialise_climate_model( region%grid, region%climate, matrix, region%PD, region%name, region%mask_noice)    
    
    ! ===== The SMB model =====
    ! =========================    
    
    CALL initialise_SMB_model( region%grid, region%init, region%SMB, region%name)     
    
    ! ===== The BMB model =====
    ! =========================    
    
    CALL initialise_BMB_model( region%grid, region%BMB, region%name)       
  
    ! ===== The ice dynamics model
    ! ============================
    
    CALL initialise_ice_model( region%grid, region%ice, region%init)
    
    ! Geothermal heat flux
    IF     (C%choice_geothermal_heat_flux == 'constant') THEN
      region%ice%GHF_Aa( :,region%grid%i1:region%grid%i2) = C%constant_geothermal_heat_flux
    ELSEIF (C%choice_geothermal_heat_flux == 'spatial') THEN
      CALL map_glob_to_grid_2D( forcing%ghf_nlat, forcing%ghf_nlon, forcing%ghf_lat, forcing%ghf_lon, region%grid, forcing%ghf_ghf, region%ice%GHF_Aa)
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_geothermal_heat_flux "', TRIM(C%choice_geothermal_heat_flux), '" not implemented in initialise_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Run the climate, BMB and SMB models once, to get the correct surface temperature field for the ice temperature initialisation,
    ! and so that all the relevant fields already have sensible values in the first time frame of the output file.
    CALL update_general_ice_model_data( region%grid, region%ice, region%time)
    CALL run_climate_model(             region%grid, region%ice, region%SMB, region%climate, region%name, C%start_time_of_run)
    CALL run_SMB_model(                 region%grid, region%ice, region%climate%applied, C%start_time_of_run, region%SMB, region%mask_noice)
    CALL run_BMB_model(                 region%grid, region%ice, region%climate%applied, region%BMB, region%name)
    
    ! Initialise the ice temperature profile
    CALL initialise_ice_temperature( region%grid, region%ice, region%climate%applied, region%init)
    
    ! Calculate physical properties again, now with the initialised temperature profile, determine the masks and slopes
    CALL update_general_ice_model_data( region%grid, region%ice, C%start_time_of_run)
    
    ! Calculate ice sheet metadata (volume, area, GMSL contribution) for writing to the first line of the output file
    CALL calculate_icesheet_volume_and_area(region)
    
    ! ===== The isotopes model =====
    ! ==============================
   
    CALL initialise_isotopes_model( region)
    
    ! ===== The ELRA GIA model =====
    ! ==============================
    
    IF (C%choice_GIA_model == 'ELRA') THEN
      CALL initialise_ELRA_model( region%grid, region%ice, region%PD)
    ELSE
      WRITE(0,*) '  ERROR - choice_GIA_model "', C%choice_GIA_model, '" not implemented in initialise_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! ===== ASCII output (computation time tracking, general output)
    IF (par%master) CALL create_text_output_files(region)
    CALL sync
    
    IF (par%master) WRITE (0,*) ' Finished initialising model region ', region%name, '.'
    
  END SUBROUTINE initialise_model
  SUBROUTINE initialise_model_grid( region)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: nsx, nsy, i, j
    
    nsx = 0
    nsy = 0
    
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
  
      ! Resolution and projection parameters for this region are determined from the config
      IF     (region%name == 'NAM') THEN
        region%grid%dx           = C%dx_NAM
        region%grid%lambda_M     = C%lambda_M_NAM
        region%grid%phi_M        = C%phi_M_NAM
        region%grid%alpha_stereo = C%alpha_stereo_NAM
      ELSEIF (region%name == 'EAS') THEN
        region%grid%dx           = C%dx_EAS
        region%grid%lambda_M     = C%lambda_M_EAS
        region%grid%phi_M        = C%phi_M_EAS
        region%grid%alpha_stereo = C%alpha_stereo_EAS
      ELSEIF (region%name == 'GRL') THEN
        region%grid%dx           = C%dx_GRL
        region%grid%lambda_M     = C%lambda_M_GRL
        region%grid%phi_M        = C%phi_M_GRL
        region%grid%alpha_stereo = C%alpha_stereo_GRL
      ELSEIF (region%name == 'ANT') THEN
        region%grid%dx           = C%dx_ANT
        region%grid%lambda_M     = C%lambda_M_ANT
        region%grid%phi_M        = C%phi_M_ANT
        region%grid%alpha_stereo = C%alpha_stereo_ANT
      END IF
      
      ! Determine the number of grid cells we can fit in this domain (based on specified resolution, and the domain covered by the initial file)
      nsx = FLOOR( (1 / region%grid%dx) * (MAXVAL(region%init%x) + (region%init%x(2)-region%init%x(1))/2) - 1)
      nsy = FLOOR( (1 / region%grid%dx) * (MAXVAL(region%init%y) + (region%init%y(2)-region%init%y(1))/2) - 1)
      
      IF (C%do_benchmark_experiment .AND. C%choice_benchmark_experiment == 'SSA_icestream') nsx = 3
      
      region%grid%nx = 1 + 2*nsx
      region%grid%ny = 1 + 2*nsy
    
      ! If doing a restart from a previous run with the same resolution, just take that grid
      ! (due to rounding, off-by-one errors are nasty here)
      IF (C%is_restart) THEN
        IF (ABS( 1._dp - region%grid%dx / ABS(region%init%x(2) - region%init%x(1))) < 1E-4_dp) THEN
          region%grid%nx = region%init%nx
          region%grid%ny = region%init%ny
          nsx = INT((REAL(region%grid%nx,dp) - 1._dp) / 2._dp)
          nsy = INT((REAL(region%grid%ny,dp) - 1._dp) / 2._dp)
        ELSE
          WRITE(0,*) '  ERROR - config resolution of ', region%grid%dx, ' is different from restart file resolution of ', ABS(region%init%x(2) - region%init%x(1))
          CALl MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END IF
      
      ! If prescribed model resolution is the same as the initial file resolution, just use the initial file frid
      IF (ABS(1._dp - (region%grid%dx / (region%init%x(2) - region%init%x(1)))) < 1E-4_dp) THEN
        region%grid%nx = region%init%nx
        region%grid%ny = region%init%ny
      END IF
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    IF (par%master) WRITE (0,'(A,F5.2,A,I4,A,I4,A)') '   Initialising model grid at ', region%grid%dx/1000._dp, ' km resolution: [', region%grid%nx, ' x ', region%grid%ny, '] pixels'
    
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
  SUBROUTINE initialise_mask_noice( region)
    ! Mask a certain area where no ice is allowed to grow. This is used to "remove"
    ! Greenland from NAM and EAS, and Ellesmere Island from GRL.
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(2)                             :: pa, pb, pc, pd
    REAL(dp)                                           :: yl_ab, yl_bc, yl_cd
    
    ! Allocate shared memory
    CALL allocate_shared_int_2D( region%grid%ny, region%grid%nx, region%mask_noice, region%wmask_noice)
    
    region%mask_noice(    :,region%grid%i1:region%grid%i2) = 0
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_mask_remove_land!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    IF (region%name == 'NAM') THEN
      ! North America: remove Greenland
      
      pa = [ 490000._dp, 1530000._dp]
      pb = [2030000._dp,  570000._dp]
      
      DO i = region%grid%i1, region%grid%i2
        yl_ab = pa(2) + (region%grid%x(i) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
        DO j = 1, region%grid%ny
          IF (region%grid%y(j) > yl_ab .AND. region%grid%x(i) > pa(1) .AND. region%grid%y(j) > pb(2)) THEN
            region%mask_noice( j,i) = 1
          END IF
        END DO
      END DO
      CALL sync
    
    ELSEIF (region%name == 'EAS') THEN
      ! Eurasia: remove Greenland
    
      pa = [-2900000._dp, 1300000._dp]
      pb = [-1895000._dp,  900000._dp]
      pc = [ -835000._dp, 1135000._dp]
      pd = [ -400000._dp, 1855000._dp]
      
      DO i = region%grid%i1, region%grid%i2
        yl_ab = pa(2) + (region%grid%x(i) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
        yl_bc = pb(2) + (region%grid%x(i) - pb(1))*(pc(2)-pb(2))/(pc(1)-pb(1))
        yl_cd = pc(2) + (region%grid%x(i) - pc(1))*(pd(2)-pc(2))/(pd(1)-pc(1))
        DO j = 1, region%grid%ny
          IF ((region%grid%x(i) <  pa(1) .AND. region%grid%y(j) > pa(2)) .OR. &
              (region%grid%x(i) >= pa(1) .AND. region%grid%x(i) < pb(1) .AND. region%grid%y(j) > yl_ab) .OR. &
              (region%grid%x(i) >= pb(1) .AND. region%grid%x(i) < pc(1) .AND. region%grid%y(j) > yl_bc) .OR. &
              (region%grid%x(i) >= pc(1) .AND. region%grid%x(i) < pd(1) .AND. region%grid%y(j) > yl_cd)) THEN
            region%mask_noice( j,i) = 1
          END IF
        END DO
      END DO
      CALL sync
      
    ELSEIF (region%name == 'GRL') THEN
      ! Greenland: remove Ellesmere island
      
      pa = [-750000._dp,  900000._dp]
      pb = [-250000._dp, 1250000._dp]
      
      DO i = region%grid%i1, region%grid%i2
        yl_ab = pa(2) + (region%grid%x(i) - pa(1))*(pb(2)-pa(2))/(pb(1)-pa(1))
        DO j = 1, region%grid%ny
          IF (region%grid%y(j) > pa(2) .AND. region%grid%y(j) > yl_ab .AND. region%grid%x(i) < pb(1)) THEN
            region%mask_noice( j,i) = 1
          END IF
        END DO
      END DO
      CALL sync
      
    ELSEIF (region%name == 'ANT') THEN
      ! Antarctica: no changes needed
      
    END IF ! IF (region%name == 'NAM') THEN
    
  END SUBROUTINE initialise_mask_noice
  
  ! Calculate critical time steps from the SIA and SSA, and determine which actions should be taken
  SUBROUTINE determine_timesteps_and_actions( region, t_end)
    ! This subroutine calculates the SIA time step with the CFK-criterium
    ! Methods are similar to those in ANICE, which are based on PISM (see also Bueler  et al., 2007)
    
    ! It also determines whether the SIA velocities, SSA velocities or thermodynamics should be updated,
    ! or whether results should be written to the output file
    
    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end
    
    ! Local variables:
    INTEGER                                     :: i,j,k,cerr,ierr
    REAL(dp)                                    :: dt_D_2D,     dt_D_2D_min
    REAL(dp)                                    :: dt_V_3D_SIA, dt_V_3D_SIA_min
    REAL(dp)                                    :: dt_V_2D_SSA, dt_V_2D_SSA_min
    REAL(dp)                                    :: t_next_action
    
    REAL(dp), PARAMETER                         :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.
    
    ! Calculate the critical time steps resulting from the SIA and SSA velocities
    ! ===========================================================================

    dt_D_2D_min     = 1000._dp
    dt_V_3D_SIA_min = 1000._dp
    dt_V_2D_SSA_min = 1000._dp

    ! As in PISM we check the three 'domains' of diffusivity (SIA), the grounded 3-D SIA velocities
    ! and the shelf velocities (SSA). Equations are as in Bueler et al. (JoG, 2007)
    
    DO i = MAX(2,region%grid%i1), MIN(region%grid%nx-1,region%grid%i2)
    DO j = 2, region%grid%ny-1
    
      dt_D_2D = region%grid%dx * region%grid%dx / (-6._dp * MAX( region%ice%D_SIA_Acx( j,i) - 1E-09, region%ice%D_SIA_Acy( j,i) - 1E-09))
      dt_D_2D_min = MIN(dt_D_2D, dt_D_2D_min) 
      
     ! dt_V_2D_SSA = region%grid%dx / (ABS(region%ice%U_SSA_Aa( j,i)) + ABS(region%ice%V_SSA_Aa( j,i)))
     ! dt_V_2D_SSA_min = MIN(dt_V_2D_SSA, dt_V_2D_SSA_min)
      
      dt_V_2D_SSA = region%grid%dx / (ABS(region%ice%U_vav_SSA_Aa( j,i) + MAX(ABS(region%ice%U_vav_SIA_Acx( j,i-1)), ABS(region%ice%U_vav_SIA_Acx( j,i)))) + &
                                      ABS(region%ice%V_vav_SSA_Aa( j,i) + MAX(ABS(region%ice%V_vav_SIA_Acy( j-1,i)), ABS(region%ice%V_vav_SIA_Acy( j,i)))))
      dt_V_2D_SSA_min = MIN(dt_V_2D_SSA, dt_V_2D_SSA_min)
      
      DO k = 1, C%nZ
        dt_V_3D_SIA = SQRT(region%grid%dx**2 + region%grid%dx**2) / (ABS(region%ice%U_3D_Aa( k,j,i)) + ABS(region%ice%V_3D_Aa( k,j,i)))
        dt_V_3D_SIA_min = MIN(dt_V_3D_SIA, dt_V_3D_SIA_min)
      END DO
       
    END DO
    END DO
    
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dT_D_2D_min    , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_V_3D_SIA_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_V_2D_SSA_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
     
    ! Correction factor. Using the exact result of the CFK criterion result in some small not-quite-unstable oscillations
    ! in the solution; better to take a slightly larger time step to prevent this.
    dt_D_2D_min     = dt_D_2D_min     * dt_correction_factor
    dt_V_2D_SSA_min = dt_V_2D_SSA_min * dt_correction_factor
    dt_V_3D_SIA_min = dt_V_3D_SIA_min * dt_correction_factor
    
    ! Calculate ice dynamics critical time step - only for writing to screen, the actual time step can be smaller
    ! if we need to do thermodynamics or write to output...
    region%dt = MIN( MIN( MIN( dt_D_2D_min, dt_V_2D_SSA_min), dt_V_3D_SIA_min), C%dt_max)      
    IF (ABS(1._dp - region%dt / region%dt_prev) > 0.1_dp) THEN
      region%dt_prev = region%dt
      !IF (par%master) WRITE(0,'(A,F7.4,A)') '   Critical time step: ', region%dt, ' yr'
    END IF

    ! Determine which action should be taken next
    ! ===========================================
    
    region%dt_SIA = MIN(C%dt_max, MIN(dt_D_2D_min, dt_V_3D_SIA_min))
    region%dt_SSA = MIN(C%dt_max, dt_V_2D_SSA_min)
    
    region%t1_SIA      = region%t0_SIA      + region%dt_SIA
    region%t1_SSA      = region%t0_SSA      + region%dt_SSA
    region%t1_thermo   = region%t0_thermo   + C%dt_thermo
    region%t1_climate  = region%t0_climate  + C%dt_climate
    region%t1_SMB      = region%t0_SMB      + C%dt_SMB
    region%t1_BMB      = region%t0_BMB      + C%dt_BMB
    region%t1_ELRA     = region%t0_ELRA     + C%dt_bedrock_ELRA
    region%t1_output   = region%t0_output   + C%dt_output
    
    t_next_action = MINVAL( [region%t1_SIA, region%t1_SSA, region%t1_thermo, region%t1_climate, region%t1_SMB, region%t1_BMB, region%t1_ELRA, region%t1_output])
    
    region%dt = t_next_action - region%time
    
    region%do_solve_SIA       = .FALSE.
    region%do_solve_SSA       = .FALSE.
    region%do_thermodynamics  = .FALSE.
    region%do_climate         = .FALSE.
    region%do_SMB             = .FALSE.
    region%do_BMB             = .FALSE.
    region%do_ELRA            = .FALSE.
    region%do_write_output    = .FALSE.
    
    IF (t_next_action == region%t1_SIA     ) region%do_solve_SIA      = .TRUE.
    IF (t_next_action == region%t1_SSA     ) region%do_solve_SSA      = .TRUE.
    IF (t_next_action == region%t1_thermo  ) region%do_thermodynamics = .TRUE.
    IF (t_next_action == region%t1_climate ) region%do_climate        = .TRUE.
    IF (t_next_action == region%t1_SMB     ) region%do_SMB            = .TRUE.
    IF (t_next_action == region%t1_BMB     ) region%do_BMB            = .TRUE.
    IF (t_next_action == region%t1_ELRA    ) region%do_ELRA           = .TRUE.
    IF (t_next_action == region%t1_output  ) region%do_write_output   = .TRUE.
    
    ! If the next action occurs after the coupling interval, crop the time step and set all "do" logicals to TRUE
    ! (except for writing output, that one really should happen only at the prescribed interval)
    IF (t_next_action >= t_end) THEN
      region%dt = t_end - region%time
      region%do_solve_SIA       = .TRUE.
      region%do_solve_SSA       = .TRUE.
      region%do_thermodynamics  = .TRUE.   
      region%do_climate         = .TRUE.
      region%do_SMB             = .TRUE.
      region%do_BMB             = .TRUE.  
    END IF
            
    ! Advance region time
    region%time = region%time + region%dt
    
   ! WRITE(0,*) ' time  = ', region%time
   
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
              C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'Bueler' .OR. &
              C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        ! No changes here
      ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
        ! Stop the simulation here, as we're only interested in SSA velocities at t=0
        region%time = t_end
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in determine_timesteps_and_actions!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
        
  END SUBROUTINE determine_timesteps_and_actions
  
  ! Calculate this region's ice sheet's volume and area
  SUBROUTINE calculate_icesheet_volume_and_area( region)
    
    USE parameters_module,           ONLY: ocean_area, seawater_density, ice_density
  
    IMPLICIT NONE  
    
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    
    INTEGER                                       :: i,j, ierr
    REAL(dp)                                      :: ice_area, ice_volume, thickness_above_flotation, ice_volume_above_flotation
    
    ice_area                   = 0._dp
    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp
    
    ! Calculate ice area and volume for processor domain
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
    
      IF (region%ice%mask_ice_Aa( j,i) == 1) THEN
        ice_volume = ice_volume + (region%ice%Hi_Aa(j,i) * region%grid%dx * region%grid%dx * ice_density / (seawater_density * ocean_area))
        ice_area   = ice_area   + region%grid%dx * region%grid%dx * 1.0E-06_dp ! [km^3]
    
        ! Thickness above flotation
        thickness_above_flotation = MAX(0._dp, region%ice%Hi_Aa( j,i) - ((0._dp - region%ice%Hb_Aa( j,i)) * (seawater_density / ice_density)))
        
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
    
    INTEGER                                       :: i, j, ierr
    REAL(dp)                                      :: ice_volume, thickness_above_flotation, ice_volume_above_flotation
    
    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp
    
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
    
      ! Thickness above flotation
      IF (region%PD%Hi( j,i) > 0._dp) THEN
        thickness_above_flotation = MAX(0._dp, region%PD%Hi( j,i) - ((0._dp - region%PD%Hb( j,i)) * (seawater_density / ice_density)))
      ELSE
        thickness_above_flotation = 0._dp
      END IF
      
      ! Ice volume (above flotation) in m.s.l.e
      ice_volume                 = ice_volume                 + region%PD%Hi( j,i)        * region%grid%dx * region%grid%dx * ice_density / (seawater_density * ocean_area)
      ice_volume_above_flotation = ice_volume_above_flotation + thickness_above_flotation * region%grid%dx * region%grid%dx * ice_density / (seawater_density * ocean_area)
      
    END DO
    END DO
    
    CALL MPI_REDUCE( ice_volume                , region%ice_volume_PD,                 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume_above_flotation, region%ice_volume_above_flotation_PD, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
  END SUBROUTINE calculate_PD_sealevel_contribution
  
  ! Create and write to this region's text output files
  SUBROUTINE create_text_output_files( region)
    ! Creates the following text output files:
    !   time_log_REG.txt             - a log of how much computation time the different model parts take
    !   general_output_REG.txt       - some general info - ice sheet volume, average surface temperature, total mass balance, etc.
  
    IMPLICIT NONE  
    
    TYPE(type_model_region),    INTENT(IN)        :: region
    
    CHARACTER(LEN=256)                            :: filename
        
  ! Time log
  ! ========
    
    filename = TRIM(C%output_dir) // 'aa_time_log_' // region%name // '.txt'
    OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')
    
    WRITE(UNIT = 1337, FMT = '(A)') 'Time log for region ' // TRIM(region%long_name)
    WRITE(UNIT = 1337, FMT = '(A)') 'Computation time (in seconds) required by each model component'
    WRITE(UNIT = 1337, FMT = '(A)') ''
    WRITE(UNIT = 1337, FMT = '(A)') '     Time       total       SIA        SSA      thermo     climate'
    
    CLOSE(UNIT = 1337)
        
  ! General output
  ! ==============
    
    filename = TRIM(C%output_dir) // 'aa_general_output_' // region%name // '.txt'
    OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')
    
    WRITE(UNIT = 1337, FMT = '(A)') 'General output for region ' // TRIM(region%long_name)
    WRITE(UNIT = 1337, FMT = '(A)') ''
    WRITE(UNIT = 1337, FMT = '(A)') ' Columns in order:'
    WRITE(UNIT = 1337, FMT = '(A)') '   1)  Model time                  (years) '
    WRITE(UNIT = 1337, FMT = '(A)') '   2)  Ice volume                  (meter sea level equivalent)'
    WRITE(UNIT = 1337, FMT = '(A)') '   3)  Ice volume above flotation  (meter sea level equivalent)'
    WRITE(UNIT = 1337, FMT = '(A)') '   4)  Ice area                    (km^2)'
    WRITE(UNIT = 1337, FMT = '(A)') '   5)  Mean surface temperature    (Kelvin)'
    WRITE(UNIT = 1337, FMT = '(A)') '   6)  Total snowfall     over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '   7)  Total rainfall     over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '   8)  Total melt         over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '   9)  Total refreezing   over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '  10)  Total runoff       over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '  11)  Total SMB          over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '  12)  Total BMB          over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') '  13)  Total mass balance over ice (Gton/y)'
    WRITE(UNIT = 1337, FMT = '(A)') ''
    WRITE(UNIT = 1337, FMT = '(A)') '     Time     Ice  Ice-af     Ice-area     T2m       Snow       Rain       Melt   Refreeze     Runoff        SMB        BMB         MB'
    
    CLOSE(UNIT = 1337)
    
  ! Benchmark experiment output
  ! ===========================
    
    ! For many benchmark experiments, we're only interested in what happens at the ice divide (i.e. the central pixel),
    ! so instead of writing the entire NetCDF output every 100 years, write only some text for that pixel.
    
    IF (C%do_benchmark_experiment) THEN
    
      filename = TRIM(C%output_dir) // 'aa_benchmark_output.txt'
      OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')
      
      WRITE(UNIT = 1337, FMT = '(A)') 'Central pixel data, useful for the EISMINT (and other) benchmark experiments'
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A)') ' Columns in order:'
      WRITE(UNIT = 1337, FMT = '(A)') '   1)  Model time                  (years) '
      WRITE(UNIT = 1337, FMT = '(A)') '   2)  Ice thickness               (meter)'
      WRITE(UNIT = 1337, FMT = '(A)') '   3)  Basal temperature           (Kelvin)'
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A)') '     Time     Hi       Ti_basal'
      
      CLOSE(UNIT = 1337)
    
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  END SUBROUTINE create_text_output_files
  SUBROUTINE write_text_output( region)
    ! Write data to the following text output files:
    !   time_log_REG.txt             - a log of how much computation time the different model parts take
    !   general_output_REG.txt       - some general info - ice sheet volume, average surface temperature, total mass balance, etc.
    
    USE parameters_module,           ONLY: ocean_area, seawater_density, ice_density
  
    IMPLICIT NONE  
    
    TYPE(type_model_region),    INTENT(IN)        :: region
    
    CHARACTER(LEN=256)                            :: filename
    INTEGER                                       :: i,j,m,ierr
    REAL(dp)                                      :: T2m_mean
    REAL(dp)                                      :: total_snowfall
    REAL(dp)                                      :: total_rainfall
    REAL(dp)                                      :: total_melt
    REAL(dp)                                      :: total_refreezing
    REAL(dp)                                      :: total_runoff
    REAL(dp)                                      :: total_SMB
    REAL(dp)                                      :: total_BMB
    REAL(dp)                                      :: total_MB
        
  ! Time log
  ! ========
    
    IF (par%master) THEN
      filename = TRIM(C%output_dir) // 'aa_time_log_' // region%name // '.txt'
      OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')
      WRITE(UNIT = 1337, FMT = '(F10.1,F11.3,F11.3,F11.3,F11.3,F11.3,F11.3)') region%time, &
        region%tcomp_total, region%tcomp_SIA, region%tcomp_SSA, region%tcomp_thermo, region%tcomp_climate
      CLOSE(UNIT = 1337)
    END IF ! IF (par%master) THEN
    CALL sync
    
  ! General output
  ! ==============
    
    T2m_mean                   = 0._dp
    total_snowfall             = 0._dp
    total_rainfall             = 0._dp
    total_melt                 = 0._dp
    total_refreezing           = 0._dp
    total_runoff               = 0._dp
    total_SMB                  = 0._dp
    total_BMB                  = 0._dp
    total_MB                   = 0._dp
    
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
    
      IF (region%ice%Hi_Aa( j,i) > 0._dp) THEN
        
        total_BMB = total_BMB + (region%BMB%BMB( j,i) * region%grid%dx * region%grid%dx / 1E9_dp) / REAL(par%n,dp)
          
        DO m = 1, 12
          total_snowfall   = total_snowfall   + (region%SMB%Snowfall(   m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp) / REAL(par%n,dp)
          total_rainfall   = total_rainfall   + (region%SMB%Rainfall(   m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp) / REAL(par%n,dp)
          total_melt       = total_melt       + (region%SMB%Melt(       m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp) / REAL(par%n,dp)
          total_refreezing = total_refreezing + (region%SMB%Refreezing( m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp) / REAL(par%n,dp)
          total_runoff     = total_runoff     + (region%SMB%Runoff(     m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp) / REAL(par%n,dp)
          total_SMB        = total_SMB        + (region%SMB%SMB(        m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp) / REAL(par%n,dp)
        END DO
        
      END IF

      T2m_mean = T2m_mean + SUM(region%climate%applied%T2m( :,j,i)) / (12._dp * region%grid%nx * region%grid%ny) / REAL(par%n,dp)
      
    END DO
    END DO
    CALL sync
    
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, T2m_mean        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_snowfall  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_rainfall  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_melt      , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_refreezing, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_runoff    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_SMB       , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_BMB       , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_MB        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    total_MB = total_SMB + total_BMB
    
    IF (par%master) THEN
      filename = TRIM(C%output_dir) // 'aa_general_output_' // region%name // '.txt'
      OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')
      WRITE(UNIT = 1337, FMT = '(F10.1,2F8.2,F13.2,F8.2,8F11.2)') region%time, &
        region%ice_volume, region%ice_volume_above_flotation, region%ice_area, T2m_mean, &
        total_snowfall, total_rainfall, total_melt, total_refreezing, total_runoff, total_SMB, total_BMB, total_MB
      CLOSE(UNIT = 1337)
    END IF ! IF (par%master) THEN
    CALL sync
    
  ! Benchmark experiment output (ice thickness and temperature at the ice divide)
  ! ===========================
    
    i = CEILING(REAL(region%grid%nx,dp)/2._dp)
    j = CEILING(REAL(region%grid%ny,dp)/2._dp)
    
    IF (par%master .AND. C%do_benchmark_experiment) THEN
      filename = TRIM(C%output_dir) // 'aa_benchmark_output.txt'
      OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')
      WRITE(UNIT = 1337, FMT = '(F10.1,F12.3,F12.3)') region%time, region%ice%Hi_Aa( j,i), region%ice%Ti_Aa( C%nZ,j,i)
      CLOSE(UNIT = 1337)
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE write_text_output

END MODULE IMAU_ICE_main_model
