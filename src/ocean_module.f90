MODULE ocean_module

  ! Contains all the routines for calculating the climate forcing.

  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_climate_matrix, type_subocean_global, &
                                             type_climate_model, type_subclimate_region, type_init_data_fields
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, &
                                             inquire_PD_obs_data_file_ocean, read_PD_obs_data_file_ocean, &
                                             inquire_GCM_ocean_snapshot, read_GCM_ocean_snapshot
  USE forcing_module,                  ONLY: forcing
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             error_function, smooth_Gaussian_2D, smooth_Gaussian_3D, smooth_Shepard_2D, &
                                             map_glob_to_grid_2D, map_glob_to_grid_3D, &
                                             map_square_to_square_cons_2nd_order_2D, map_square_to_square_cons_2nd_order_3D, &
                                             remap_cons_2nd_order_1D, surface_elevation, extrapolate_Gaussian_floodfill

  IMPLICIT NONE
    
CONTAINS

  ! Run the climate model on a region grid
  SUBROUTINE run_ocean_model( grid, ice, climate, region_name, time)
    ! Run the regional climate model
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    REAL(dp),                            INTENT(IN)    :: time
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

  ! ================================================
  ! ===== Exceptions for benchmark experiments =====
  ! ================================================
  
    IF (C%do_benchmark_experiment) THEN
      IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
              C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'Bueler' .OR. &
              C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
              C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
              C%choice_benchmark_experiment == 'ISMIP_HOM_F' .OR. &
              C%choice_benchmark_experiment == 'MISMIPplus') THEN
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'MISOMIP1') THEN
        ! Set ocean temperature/salinity profiles according to the MISOMIP+ protocol
        CALL MISOMIP1_ocean_profiles( grid, climate%applied, time)
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_ocean_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  ! =======================================================
  ! ===== End of exceptions for benchmark experiments =====
  ! =======================================================
    
    
    IF (C%choice_ocean_temperature_model == 'fixed' .OR. &
        C%choice_ocean_temperature_model == 'scaled' .OR. &
        C%choice_ocean_temperature_model == 'schematic' .OR. &
        C%choice_ocean_temperature_model == 'WOA') THEN
      ! No need to update ocean temperature and salinity in these cases
      RETURN        
    ELSE IF (C%choice_ocean_temperature_model == 'matrix_warm_cold') THEN
      CALL run_ocean_model_matrix_warm_cold (grid, ice, climate, region_name)
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: ocean model"', TRIM(C%choice_ocean_temperature_model), '" not implemented in run_ocean_model!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr) 
    END IF
      
  END SUBROUTINE run_ocean_model
  
  ! Ocean matrix with warm and cold snapshots
  SUBROUTINE run_ocean_model_matrix_warm_cold ( grid, ice, climate, region_name )
  
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: CO2
    REAL(dp)                                           :: w_CO2, w_ice, w_tot
    REAL(dp)                                           :: w_cold, w_warm
    REAL(dp), PARAMETER                                :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]

      
    ! Find CO2 interpolation weight (use either prescribed or modelled CO2)
    ! =====================================================================
    
    IF (C%choice_forcing_method == 'CO2_direct') THEN
      CO2 = forcing%CO2_obs
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CO2 = forcing%CO2_mod
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      CO2 = 0._dp
      WRITE(0,*) '  ERROR - run_ocean_model_matrix_warm_cold must only be called with the correct forcing method, check your code!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSE
      CO2 = 0._dp
      WRITE(0,*) '  ERROR - choice_forcing_method "', C%choice_forcing_method, '" not implemented in run_ocean_model_matrix_warm_cold!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    w_CO2 = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (CO2 - C%matrix_low_CO2_level) / (C%matrix_high_CO2_level - C%matrix_low_CO2_level) ))
    
    ! Find ice interpolation weight 
    ! =============================
    
    ! First calculate the total ice volume term (second term in the equation)
    w_ice = 1._dp - MAX(-w_cutoff, MIN(1._dp + w_cutoff, (SUM(ice%Hs_a) - SUM(climate%GCM_warm%Hs)) / (SUM(climate%GCM_cold%Hs) - SUM(climate%GCM_warm%Hs)) ))
    
    ! Combine weigths CO2 and ice
    ! ===========================
    
    w_tot = 0._dp
    IF         (region_name == 'NAM') THEN
      w_tot  = (C%ocean_matrix_CO2vsice_NAM * w_CO2) + ((1._dp - C%ocean_matrix_CO2vsice_NAM) * w_ice) 
    ELSEIF     (region_name == 'EAS') THEN
      w_tot  = (C%ocean_matrix_CO2vsice_EAS * w_CO2) + ((1._dp - C%ocean_matrix_CO2vsice_EAS) * w_ice) 
    ELSEIF     (region_name == 'GRL') THEN
      w_tot  = (C%ocean_matrix_CO2vsice_GRL * w_CO2) + ((1._dp - C%ocean_matrix_CO2vsice_GRL) * w_ice) 
    ELSEIF     (region_name == 'ANT') THEN
      w_tot  = (C%ocean_matrix_CO2vsice_ANT * w_CO2) + ((1._dp - C%ocean_matrix_CO2vsice_ANT) * w_ice) 
    END IF
    w_warm = w_tot
    w_cold = 1._dp - w_warm
    
    ! Interpolate the GCM ocean snapshots
    ! =============================
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      climate%applied%T_ocean_corr_ext(:,j,i) = (w_tot * climate%GCM_warm%T_ocean_corr_ext( :,j,i)) + ((1._dp - w_tot) * climate%GCM_cold%T_ocean_corr_ext( :,j,i))
      climate%applied%S_ocean_corr_ext(:,j,i) = (w_tot * climate%GCM_warm%S_ocean_corr_ext( :,j,i)) + ((1._dp - w_tot) * climate%GCM_cold%S_ocean_corr_ext( :,j,i))      
    END DO
    END DO
    CALL sync
      
  END SUBROUTINE run_ocean_model_matrix_warm_cold

  ! Initialising the ocean matrix, containing all the global subclimates
  ! (PD observations and GCM snapshots) on their own lat-lon grids
  SUBROUTINE initialise_ocean_matrix( matrix)
    ! Allocate shared memory for the global ocean matrix
  
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_climate_matrix),      INTENT(INOUT) :: matrix
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_F' .OR. &
          C%choice_benchmark_experiment == 'MISMIPplus' .OR. &
          C%choice_benchmark_experiment == 'MISOMIP1') THEN
        RETURN
      ELSE 
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_ocean_matrix!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Initialising the ocean matrix...'
    
    ! The global WOA18 ocean
    CALL initialise_PD_obs_ocean_fields ( matrix%PD_obs_ocean, 'WOA18')
    
    ! The different GCM snapshots 
    IF ((C%choice_ocean_temperature_model == 'fixed') .OR. (C%choice_ocean_temperature_model == 'scaled') .OR. &
        (C%choice_ocean_temperature_model == 'schematic') .OR. (C%choice_ocean_temperature_model == 'WOA')) THEN
      ! These choices of forcing don't use any GCM (snapshot) ocean data
      RETURN
    ELSEIF (C%choice_ocean_temperature_model == 'matrix_warm_cold') THEN
      ! These two choices use the ocean matrix
            
      ! Initialise the GCM ocean snapshots
      CALL initialise_ocean_snapshot( matrix%GCM_PI_ocean,   name = 'ref_PI_ocean', nc_filename = C%filename_GCM_ocean_snapshot_PI,   CO2 = 280._dp,                 orbit_time =       0._dp)
      CALL initialise_ocean_snapshot( matrix%GCM_warm_ocean, name = 'warm_ocean',   nc_filename = C%filename_GCM_ocean_snapshot_warm, CO2 = C%matrix_high_CO2_level, orbit_time = C%matrix_warm_orbit_time)
      CALL initialise_ocean_snapshot( matrix%GCM_cold_ocean, name = 'cold_ocean',   nc_filename = C%filename_GCM_ocean_snapshot_cold, CO2 = C%matrix_low_CO2_level,  orbit_time = C%matrix_cold_orbit_time)
      
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_ocean_temperature_model "', TRIM(C%choice_ocean_temperature_model), '" not implemented in initialise_ocean_matrix!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE initialise_ocean_matrix
  
! == Initialising the global ocean matrix
  SUBROUTINE initialise_PD_obs_ocean_fields (PD_obs_ocean, name)
    ! Allocate shared memory for the global PD observed ocean data fields (stored in the climate matrix),
    ! read them from the specified NetCDF file (latter only done by master process). 
     
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_subocean_global),     INTENT(INOUT) :: PD_obs_ocean
    CHARACTER(LEN=*),               INTENT(IN)    :: name
    
    PD_obs_ocean%name            = name 
    PD_obs_ocean%netcdf%filename = C%filename_PD_obs_ocean
        
    ! Inquire if all required variables are present in the NetCDF file, and read the grid size.
    CALL allocate_shared_int_0D( PD_obs_ocean%nlon,     PD_obs_ocean%wnlon)
    CALL allocate_shared_int_0D( PD_obs_ocean%nlat,     PD_obs_ocean%wnlat)
    CALL allocate_shared_int_0D( PD_obs_ocean%nz_ocean, PD_obs_ocean%wnz_ocean)
    IF (par%master) CALL inquire_PD_obs_data_file_ocean(PD_obs_ocean)
    CALL sync
    
    ! Allocate memory  
    CALL allocate_shared_dp_1D( PD_obs_ocean%nlon,                                           PD_obs_ocean%lon,     PD_obs_ocean%wlon    )
    CALL allocate_shared_dp_1D(                    PD_obs_ocean%nlat,                        PD_obs_ocean%lat,     PD_obs_ocean%wlat    )
    CALL allocate_shared_dp_1D(                                       PD_obs_ocean%nz_ocean, PD_obs_ocean%z_ocean, PD_obs_ocean%wz_ocean)        
    
    CALL allocate_shared_dp_3D( PD_obs_ocean%nlon, PD_obs_ocean%nlat, PD_obs_ocean%nz_ocean, PD_obs_ocean%T_ocean, PD_obs_ocean%wT_ocean)
    CALL allocate_shared_dp_3D( PD_obs_ocean%nlon, PD_obs_ocean%nlat, PD_obs_ocean%nz_ocean, PD_obs_ocean%S_ocean, PD_obs_ocean%wS_ocean)
     
    ! Read data from the NetCDF file
    IF (par%master) WRITE(0,*) '   Reading PD observed ocean data from file ', TRIM(PD_obs_ocean%netcdf%filename), '...'
    IF (par%master) CALL read_PD_obs_data_file_ocean(PD_obs_ocean)
    CALL sync
      
    ! Determine process domains
    CALL partition_list( PD_obs_ocean%nlon, par%i, par%n, PD_obs_ocean%i1, PD_obs_ocean%i2)
    
    ! Map the data to the desired vertical grid
    CALL map_global_ocean_data_to_IMAUICE_vertical_grid( PD_obs_ocean)
      
  END SUBROUTINE initialise_PD_obs_ocean_fields  
  SUBROUTINE initialise_ocean_snapshot( snapshot, name, nc_filename, CO2, orbit_time)
    ! Allocate shared memory for the data fields of a GCM snapshot (stored in the climate matrix),
    ! read them from the specified NetCDF file (latter only done by master process).
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_subocean_global),     INTENT(INOUT) :: snapshot
    CHARACTER(LEN=*),               INTENT(IN)    :: name
    CHARACTER(LEN=*),               INTENT(IN)    :: nc_filename
    REAL(dp),                       INTENT(IN)    :: CO2
    REAL(dp),                       INTENT(IN)    :: orbit_time
    
    ! Metadata
    snapshot%name            = name 
    snapshot%netcdf%filename = nc_filename
    
    ! General forcing info
    CALL allocate_shared_dp_0D( snapshot%CO2,        snapshot%wCO2       )
    CALL allocate_shared_dp_0D( snapshot%orbit_time, snapshot%worbit_time)
    CALL allocate_shared_dp_0D( snapshot%orbit_ecc,  snapshot%worbit_ecc )
    CALL allocate_shared_dp_0D( snapshot%orbit_obl,  snapshot%worbit_obl )
    CALL allocate_shared_dp_0D( snapshot%orbit_pre,  snapshot%worbit_pre )
    
    snapshot%CO2        = CO2
    snapshot%orbit_time = orbit_time
    
    ! Inquire if all required variables are present in the NetCDF file, and read the grid size.
    CALL allocate_shared_int_0D( snapshot%nlon,     snapshot%wnlon)
    CALL allocate_shared_int_0D( snapshot%nlat,     snapshot%wnlat)
    CALL allocate_shared_int_0D( snapshot%nz_ocean, snapshot%wnz_ocean)
    IF (par%master) CALL inquire_GCM_ocean_snapshot( snapshot)
    CALL sync
    
    ! Allocate memory  
    CALL allocate_shared_dp_1D(                        snapshot%nlon,                    snapshot%lon,     snapshot%wlon    )
    CALL allocate_shared_dp_1D(                                           snapshot%nlat, snapshot%lat,     snapshot%wlat    )
    CALL allocate_shared_dp_1D( snapshot%nz_ocean,                                       snapshot%z_ocean, snapshot%wz_ocean)        
    !CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, snapshot%nz_ocean, snapshot%mask_ocean,      snapshot%wmask_ocean     )
    CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, snapshot%nz_ocean, snapshot%T_ocean,     snapshot%wT_ocean    )
    CALL allocate_shared_dp_3D( snapshot%nlon, snapshot%nlat, snapshot%nz_ocean, snapshot%S_ocean,     snapshot%wS_ocean    )
    
    ! Read data from the NetCDF file
    IF (par%master) WRITE(0,*) '   Reading GCM ocean snapshot ', TRIM(snapshot%name), ' from file ', TRIM(snapshot%netcdf%filename), '...'
    IF (par%master) CALL read_GCM_ocean_snapshot( snapshot)
    CALL sync
      
    ! Determine process domains
    CALL partition_list( snapshot%nlon, par%i, par%n, snapshot%i1, snapshot%i2)
    
    ! Map the data to the desired vertical grid
    CALL map_global_ocean_data_to_IMAUICE_vertical_grid( snapshot)
    
  END SUBROUTINE initialise_ocean_snapshot
  
! == Initialising the region-specific ocean data
  SUBROUTINE initialise_oceans_regional( grid, ice, climate, matrix, init)
    ! Allocate shared memory for the ocean data of the regional subclimates
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice  
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    TYPE(type_climate_matrix),           INTENT(INOUT) :: matrix
    TYPE(type_init_data_fields),         INTENT(IN)    :: init
        
    IF (par%master) WRITE (0,*) '  Initialising ocean model...'
    
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
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
          C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
        ! Nothing to be done here
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_oceans_regional!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Exception for schematic ocean temperature/salinity profile
    ! Set oceans to the ISOMIP+ "COLD" or "WARM" profile
    IF (C%choice_ocean_temperature_model == 'schematic') THEN
      IF (par%master) WRITE(*,*) '    Schematic"', TRIM(C%choice_schematic_ocean), '" ocean forcing used for BMB forcing!'
      IF     (C%choice_schematic_ocean == 'MISMIPplus_WARM') THEN
        CALL set_ocean_to_ISOMIPplus_WARM( grid, climate%PD_obs  )
        CALL set_ocean_to_ISOMIPplus_WARM( grid, climate%GCM_PI  )
        CALL set_ocean_to_ISOMIPplus_WARM( grid, climate%GCM_cold)
        CALL set_ocean_to_ISOMIPplus_WARM( grid, climate%GCM_warm)
        CALL set_ocean_to_ISOMIPplus_WARM( grid, climate%applied )
      ELSEIF (C%choice_schematic_ocean == 'MISMIPplus_COLD') THEN
        CALL set_ocean_to_ISOMIPplus_COLD( grid, climate%PD_obs  )
        CALL set_ocean_to_ISOMIPplus_COLD( grid, climate%GCM_PI  )
        CALL set_ocean_to_ISOMIPplus_COLD( grid, climate%GCM_cold)
        CALL set_ocean_to_ISOMIPplus_COLD( grid, climate%GCM_warm)
        CALL set_ocean_to_ISOMIPplus_COLD( grid, climate%applied )
      ELSEIF (C%choice_schematic_ocean == 'Reese2018') THEN
        CALL set_ocean_to_Reese2018( grid, ice, climate%PD_obs  )
        CALL set_ocean_to_Reese2018( grid, ice, climate%GCM_PI  )
        CALL set_ocean_to_Reese2018( grid, ice, climate%GCM_cold)
        CALL set_ocean_to_Reese2018( grid, ice, climate%GCM_warm)
        CALL set_ocean_to_Reese2018( grid, ice, climate%applied )
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_schematic_ocean "', TRIM(C%choice_schematic_ocean), '" not implemented in initialise_subclimate!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      RETURN  
    END IF
    
    ! Allocate memory for ocean data in all the regional subclimates
    CALL allocate_subclimate_regional_oceans( grid, climate%PD_obs )
    CALL allocate_subclimate_regional_oceans( grid, climate%applied)
    
    ! Map PD ocean data from the global lon/lat-grid to the regional x/y-grid
    CALL map_ocean_data_global_to_regional( grid, matrix%PD_obs_ocean, climate%PD_obs)
    
    ! Extend regional ocean data to cover the entire 3D domain
    CALL extend_regional_ocean_data_to_cover_domain( grid, ice, climate%PD_obs, init%Hi, init%Hb)
    
    IF (C%choice_ocean_temperature_model == 'WOA') THEN
      IF (par%master) WRITE(*,*) '    Constant present-day ocean forcing used for BMB forcing!'
    ELSE IF (C%choice_ocean_temperature_model == 'fixed' .OR. C%choice_ocean_temperature_model == 'scaled')  THEN
      IF (par%master) WRITE(*,*) '    Ocean forcing initialised, but not actually used for BMB forcing!'
    ELSE IF (C%choice_ocean_temperature_model == 'matrix_warm_cold') THEN 
    
      ! Allocate memory for regional ocean data     
      CALL allocate_subclimate_regional_oceans( grid, climate%GCM_PI  )
      CALL allocate_subclimate_regional_oceans( grid, climate%GCM_cold)
      CALL allocate_subclimate_regional_oceans( grid, climate%GCM_warm)
    
      ! Map ocean data from the global lon/lat-grid to the regional x/y-grid
      CALL map_ocean_data_global_to_regional( grid, matrix%GCM_PI_ocean,   climate%GCM_PI  )
      CALL map_ocean_data_global_to_regional( grid, matrix%GCM_cold_ocean, climate%GCM_cold)
      CALL map_ocean_data_global_to_regional( grid, matrix%GCM_warm_ocean, climate%GCM_warm)
    
      ! Extend regional ocean data to cover the entire 3D domain
      CALL extend_regional_ocean_data_to_cover_domain( grid, ice, climate%GCM_PI,   init%Hi, init%Hb)
      CALL extend_regional_ocean_data_to_cover_domain( grid, ice, climate%GCM_cold, init%Hi, init%Hb)
      CALL extend_regional_ocean_data_to_cover_domain( grid, ice, climate%GCM_warm, init%Hi, init%Hb)
    
      ! Correct regional ocean data for GCM bias
      CALL correct_GCM_bias_ocean( grid, climate, climate%GCM_warm)
      CALL correct_GCM_bias_ocean( grid, climate, climate%GCM_cold)
            
    ELSE  
      IF (par%master) WRITE(0,*) '  ERROR: choice_ocean_temperature_model "', TRIM(C%choice_ocean_temperature_model), '" not implemented in initialise_oceans_regional!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF !(C%choice_ocean_temperature_model == 'WOA')

    ! PI and PD do not have to be corrected 
    ! Initialise applied ocean forcing with present-day observations
    IF (par%master) THEN
      IF (C%choice_ocean_temperature_model == 'matrix_warm_cold') climate%GCM_PI%T_ocean_corr_ext  = climate%GCM_PI%T_ocean_ext
      IF (C%choice_ocean_temperature_model == 'matrix_warm_cold') climate%GCM_PI%S_ocean_corr_ext  = climate%GCM_PI%S_ocean_ext
      climate%PD_obs%T_ocean_corr_ext  = climate%PD_obs%T_ocean_ext
      climate%PD_obs%S_ocean_corr_ext  = climate%PD_obs%S_ocean_ext
    
      climate%applied%mask_ocean       = climate%PD_obs%mask_ocean
      climate%applied%T_ocean          = climate%PD_obs%T_ocean
      climate%applied%T_ocean_ext      = climate%PD_obs%T_ocean_ext
      climate%applied%T_ocean_corr_ext = climate%PD_obs%T_ocean_corr_ext
      climate%applied%S_ocean          = climate%PD_obs%S_ocean
      climate%applied%S_ocean_ext      = climate%PD_obs%S_ocean_ext
      climate%applied%S_ocean_corr_ext = climate%PD_obs%S_ocean_corr_ext
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE initialise_oceans_regional
  SUBROUTINE allocate_subclimate_regional_oceans( grid, climate)
    ! Allocate shared memory for the ocean data of a regional subclimate models
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    
    ! Set up the vertical coordinate
    CALL create_ocean_vertical_grid( climate%nz_ocean, climate%wnz_ocean, climate%z_ocean, climate%wz_ocean)
    
    ! Allocate shared memory
    CALL allocate_shared_int_3D( climate%nz_ocean, grid%ny, grid%nx, climate%mask_ocean,       climate%wmask_ocean      )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%T_ocean,          climate%wT_ocean         )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%S_ocean,          climate%wS_ocean         )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%T_ocean_ext,      climate%wT_ocean_ext     )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%S_ocean_ext,      climate%wS_ocean_ext     )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%T_ocean_corr_ext, climate%wT_ocean_corr_ext)
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%S_ocean_corr_ext, climate%wS_ocean_corr_ext)
  
  END SUBROUTINE allocate_subclimate_regional_oceans
  SUBROUTINE map_ocean_data_global_to_regional( grid, clim_glob, clim_reg)
    ! Map ocean data for a single subclimate from the global grid to the regional grid
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_subocean_global),          INTENT(IN)    :: clim_glob
    TYPE(type_subclimate_region),        INTENT(INOUT) :: clim_reg
    
    CALL map_glob_to_grid_3D ( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%T_ocean, clim_reg%T_ocean, clim_glob%nz_ocean)
    CALL map_glob_to_grid_3D ( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%S_ocean, clim_reg%S_ocean, clim_glob%nz_ocean)
  
  END SUBROUTINE map_ocean_data_global_to_regional
  SUBROUTINE extend_regional_ocean_data_to_cover_domain( grid, ice, ocean_data, Hi, Hb)
    ! Extend global ocean data over the whole grid, based on the procedure outlined in 
    ! Jourdain, N. C., Asay-Davis, X., Hattermann, T., Straneo, F., Seroussi, H., Little, C. M., & Nowicki, S. (2020). 
    ! A protocol for calculating basal melt rates in the ISMIP6 Antarctic ice sheet projections. The Cryosphere, 14(9), 3111-3134. 

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_subclimate_region),        INTENT(INOUT) :: ocean_data
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Hi, Hb   ! DENK DROM - this is the geometry for which the ocean data is valid, probably should be part of the ocean_data structure...
    
    ! Local variables:
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: NaN, Hs, z_bedrock, z_icebase, z
    INTEGER,  DIMENSION(:,:,:), POINTER                ::  mask_wetdry,  mask_hasdata
    INTEGER                                            :: wmask_wetdry, wmask_hasdata
    INTEGER                                            :: k1,k2,bi
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE            :: mask, mask_filled
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_T, d_S
    REAL(dp), PARAMETER                                :: sigma_Gaussian = 12000._dp
    
    ! Useful
    NaN = -1._dp
    NaN = SQRT( NaN)
    
    ! Initialise the extrapolated product with the provided ocean data
    ocean_data%T_ocean_ext( :,:,grid%i1:grid%i2) = ocean_data%T_ocean( :,:,grid%i1:grid%i2)
    ocean_data%S_ocean_ext( :,:,grid%i1:grid%i2) = ocean_data%S_ocean( :,:,grid%i1:grid%i2)
    CALL sync
    
    ! Define the two masks needed for the four extrapolation steps:
    ! 
    !  - mask_wetdry:
    !      1 = actually    wet (i.e. open ocean, sub-shelf cavity, above sea floor and beneath ice base)
    !      2 = potentially wet (i.e. grounded marine ice, above sea floor)
    !      3 =             dry (i.e. beneath bedrock                     )
    ! 
    !  - mask_hasdata:
    !      0 = has no data
    !      1 = has data provided
    !      2 = has data extrapolated
    
    CALL allocate_shared_int_3D( ocean_data%nz_ocean, grid%ny, grid%nx, mask_wetdry,  wmask_wetdry )
    CALL allocate_shared_int_3D( ocean_data%nz_ocean, grid%ny, grid%nx, mask_hasdata, wmask_hasdata)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      Hs = surface_elevation( Hi( j,i), Hb( j,i), 0._dp)
      z_bedrock = Hb( j,i)
      z_icebase = Hs - Hi( j,i)
      
      DO k = 1, ocean_data%nz_ocean
        
        z = -ocean_data%z_ocean( k)
        
        ! mask_wetdry
        IF (z < z_bedrock) THEN
          ! This 3D grid box is beneath the bedrock surface, so dry
          mask_wetdry( k,j,i) = 3
        ELSE
          ! This 3D grid box is above the bedrock surface so at least potentially wet
          IF (z < z_icebase) THEN
            ! This 3D grid box is above the bedrock surface and below the ice base, so it is actually wet
            mask_wetdry( k,j,i) = 1
          ELSE
            ! This 3D grid box is above the bedrock surface and above the ice base, so it is potentially wet (i.e. inside grounded marine ice)
            mask_wetdry( k,j,i) = 2
          END IF
        END IF
        
        ! mask_hasdata
        IF (ocean_data%T_ocean( k,j,i) /= ocean_data%T_ocean( k,j,i)) THEN
          ! This 3D grid box has no data (yet)
          mask_hasdata( k,j,i) = 0
        ELSE
          ! Data is already provided for this 3D grid box
          mask_hasdata( k,j,i) = 1
        END IF
        
      END DO
      
    END DO
    END DO
    CALL sync
        
!    ! DENK DROM
!    IF (par%master) THEN
!      debug%dp_3D_01( 1:10,:,:) = ocean_data%T_ocean_ext( 1:10,:,:)
!      CALL write_to_debug_file
!    END IF
!    CALL sync
    
  ! ================================================================
  ! ===== Step 1: horizontal extrapolation into shelf cavities =====
  ! ================================================================
  
    ! Here, we start with the ocean data as provided (i.e. only for open ocean), and
    ! perform a horizontal extrapolation (so for each vertical layer separately) into
    ! the shelf cavities. Only "actually wet" 3D grid cells are allowed to be filled,
    ! so the fill is limited by both bedrock sills and grounded ice.
    
    ! Allocate memory for the mask and data field of a single extrapolation step
    ALLOCATE( mask(        grid%ny, grid%nx))
    ALLOCATE( mask_filled( grid%ny, grid%nx))
    ALLOCATE( d_T(         grid%ny, grid%nx))
    ALLOCATE( d_S(         grid%ny, grid%nx))
    
    ! Parallelised by partitioning the vertical domain
    CALL partition_list( ocean_data%nz_ocean, par%i, par%n, k1, k2)
    
    DO k = k1, k2
    
      ! Extrapolate per basin
      DO bi = 1, ice%nbasins
        
        ! Define the mask and initial data fields for this particular flood-fill
        ! (i.e. this vertical layer and this basin)
        mask        = 0
        d_T         = NaN
        d_S         = NaN
        DO i = 1, grid%nx
        DO j = 1, grid%ny
          IF (ice%basin_ID( j,i) == bi) THEN
            IF (mask_hasdata( k,j,i) == 1) THEN
              ! This is where the source data comes from
              mask( j,i) = 2
              d_T(  j,i) = ocean_data%T_ocean_ext( k,j,i)
              d_S(  j,i) = ocean_data%S_ocean_ext( k,j,i)
            ELSEIF (mask_hasdata( k,j,i) == 0 .AND. mask_wetdry( k,j,i) == 1) THEN
              ! This is where we're supposed to fill it in
              mask( j,i) = 1
            END IF
          END IF
        END DO
        END DO
        
        ! Perform the flood-fill-based Gaussian extrapolation
        CALL extrapolate_Gaussian_floodfill( grid, mask, d_T, sigma_Gaussian, mask_filled)
        CALL extrapolate_Gaussian_floodfill( grid, mask, d_S, sigma_Gaussian, mask_filled)
        
        ! Copy extrapolated data to the data structure
        DO i = 1, grid%nx
        DO j = 1, grid%ny
          IF (mask_filled( j,i) == 1) THEN
            ocean_data%T_ocean_ext( k,j,i) = d_T( j,i)
            ocean_data%S_ocean_ext( k,j,i) = d_S( j,i)
            mask_hasdata(           k,j,i) = 2
          END IF
        END DO
        END DO
        
      END DO ! DO bi = 1, ice%nbasins
      
    END DO ! DO k = k1, k2
    CALL sync
    
    ! Clean up after yourself
    DEALLOCATE( mask       )
    DEALLOCATE( mask_filled)
    DEALLOCATE( d_T        )
    DEALLOCATE( d_S        )
        
!    ! DENK DROM
!    IF (par%master) THEN
!      debug%dp_3D_02( 1:10,:,:) = ocean_data%T_ocean_ext( 1:10,:,:)
!      CALL write_to_debug_file
!    END IF
!    CALL sync
    
  ! ===========================================================================
  ! ===== Step 2: vertical extrapolation into sill-blocked shelf cavities =====
  ! ===========================================================================
  
    ! Here, we start with the ocean data that has been horizontally extrapolated into
    ! the shelf cavities, allowing for bedrock topography to block the fill, so that
    ! for example the lower parts of the Filchner-Ronne and Ross cavities have not yet
    ! been filled. We now extrapolate the data vertically from the filled parts to
    ! fill those parts of the cavities. Barring any really weird geometry, the entire
    ! cavities will now be filled.
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! Move down through the vertical column
      DO k = 2, ocean_data%nz_ocean
        ! If this grid box is wet and has no data, but the one above it does,
        ! copy data from the one above it.
        IF (mask_wetdry( k,j,i) == 1 .AND. mask_hasdata( k,j,i) == 0) THEN
          ! This 3D grid box is wet but has no data
          IF (mask_hasdata( k-1,j,i) == 1 .OR. mask_hasdata( k-1,j,i) == 2) THEN
            ! The one above it has data; copy data
            mask_hasdata( k,j,i) = 2
            ocean_data%T_ocean_ext( k,j,i) = ocean_data%T_ocean_ext( k-1,j,i)
            ocean_data%S_ocean_ext( k,j,i) = ocean_data%S_ocean_ext( k-1,j,i)
          END IF
        END IF
      END DO
      
    END DO
    END DO
    CALL sync
        
!    ! DENK DROM
!    IF (par%master) THEN
!      debug%dp_3D_03( 1:10,:,:) = ocean_data%T_ocean_ext( 1:10,:,:)
!      CALL write_to_debug_file
!    END IF
!    CALL sync
    
  ! ===============================================================
  ! ===== Step 3: vertical extrapolation into ice and bedrock =====
  ! ===============================================================
  
    ! Extrapolate data vertically into 3D grid boxes that are occupied by ice
    ! or bedrock (since they might turn into ocean at some point during a simulation)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! Move down through the vertical column
      DO k = 2, ocean_data%nz_ocean
        ! If this grid box is wet and has no data, but the one above it does,
        ! copy data from the one above it.
        IF (mask_hasdata( k,j,i) == 0) THEN
          ! This 3D grid box is wet but has no data
          IF (mask_hasdata( k-1,j,i) == 1 .OR. mask_hasdata( k-1,j,i) == 2) THEN
            ! The one above it has data; copy data
            mask_hasdata( k,j,i) = 2
            ocean_data%T_ocean_ext( k,j,i) = ocean_data%T_ocean_ext( k-1,j,i)
            ocean_data%S_ocean_ext( k,j,i) = ocean_data%S_ocean_ext( k-1,j,i)
          END IF
        END IF
      END DO
      
      ! Move up through the vertical column
      DO k = ocean_data%nz_ocean-1, 1, -1
        ! If this grid box is wet and has no data, but the one above it does,
        ! copy data from the one above it.
        IF (mask_hasdata( k,j,i) == 0) THEN
          ! This 3D grid box is wet but has no data
          IF (mask_hasdata( k+1,j,i) == 1 .OR. mask_hasdata( k+1,j,i) == 2) THEN
            ! The one above it has data; copy data
            mask_hasdata( k,j,i) = 2
            ocean_data%T_ocean_ext( k,j,i) = ocean_data%T_ocean_ext( k+1,j,i)
            ocean_data%S_ocean_ext( k,j,i) = ocean_data%S_ocean_ext( k+1,j,i)
          END IF
        END IF
      END DO
      
    END DO
    END DO
    CALL sync
        
!    ! DENK DROM
!    IF (par%master) THEN
!      debug%dp_3D_04( 1:10,:,:) = ocean_data%T_ocean_ext( 1:10,:,:)
!      CALL write_to_debug_file
!    END IF
!    CALL sync
    
  ! =================================================================
  ! ===== Step 4: horizontal extrapolation into ice and bedrock =====
  ! =================================================================
  
    ! In the last step, extrapolate data horizontally into 3D
    ! grid boxes that are occupied by ice or bedrock
    
    ! Allocate memory for the mask and data field of a single extrapolation step
    ALLOCATE( mask(        grid%ny, grid%nx))
    ALLOCATE( mask_filled( grid%ny, grid%nx))
    ALLOCATE( d_T(         grid%ny, grid%nx))
    ALLOCATE( d_S(         grid%ny, grid%nx))
    
    ! Parallelised by partitioning the vertical domain
    CALL partition_list( ocean_data%nz_ocean, par%i, par%n, k1, k2)
    
    DO k = k1, k2
    
      ! Extrapolate per basin
      DO bi = 1, ice%nbasins
        
        ! Define the mask and initial data fields for this particular flood-fill
        ! (i.e. this vertical layer and this basin)
        mask        = 0
        d_T         = NaN
        d_S         = NaN
        DO i = 1, grid%nx
        DO j = 1, grid%ny
          IF (ice%basin_ID( j,i) == bi) THEN
            IF (mask_hasdata( k,j,i) == 1 .OR. mask_hasdata( k,j,i) == 2) THEN
              ! This is where the source data comes from
              mask( j,i) = 2
              d_T(  j,i) = ocean_data%T_ocean_ext( k,j,i)
              d_S(  j,i) = ocean_data%S_ocean_ext( k,j,i)
            ELSEIF (mask_hasdata( k,j,i) == 0) THEN
              ! This is where we're supposed to fill it in
              mask( j,i) = 1
            END IF
          END IF
        END DO
        END DO
        
        ! Perform the flood-fill-based Gaussian extrapolation
        CALL extrapolate_Gaussian_floodfill( grid, mask, d_T, sigma_Gaussian, mask_filled)
        CALL extrapolate_Gaussian_floodfill( grid, mask, d_S, sigma_Gaussian, mask_filled)
        
        ! Copy extrapolated data to the data structure
        DO i = 1, grid%nx
        DO j = 1, grid%ny
          IF (mask_filled( j,i) == 1) THEN
            ocean_data%T_ocean_ext( k,j,i) = d_T( j,i)
            ocean_data%S_ocean_ext( k,j,i) = d_S( j,i)
            mask_hasdata(           k,j,i) = 2
          END IF
        END DO
        END DO
        
      END DO ! DO bi = 1, ice%nbasins
      
    END DO ! DO k = k1, k2
    CALL sync
    
    ! Clean up after yourself
    DEALLOCATE( mask       )
    DEALLOCATE( mask_filled)
    DEALLOCATE( d_T        )
    DEALLOCATE( d_S        )
        
!    ! DENK DROM
!    IF (par%master) THEN
!      debug%dp_3D_05( 1:10,:,:) = ocean_data%T_ocean_ext( 1:10,:,:)
!      CALL write_to_debug_file
!    END IF
!    CALL sync
!    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    
    ! Clean up after yourself
    CALL deallocate_shared( wmask_wetdry )
    CALL deallocate_shared( wmask_hasdata)
    
  END SUBROUTINE extend_regional_ocean_data_to_cover_domain
  SUBROUTINE correct_GCM_bias_ocean( grid, climate, subclimate)
    ! Correct regional ocean data for GCM bias
    ! (must be done on regional grid, since the GCM grid and the World Ocean Atlas grid are generally not the same!)
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_model),            INTENT(IN)    :: climate
    TYPE(type_subclimate_region),        INTENT(INOUT) :: subclimate
    
    ! Local variables
    INTEGER                                            :: i,j,k
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO k = 1, subclimate%nz_ocean
      subclimate%T_ocean_corr_ext( k,j,i) = subclimate%T_ocean_ext( k,j,i) - (climate%GCM_PI%T_ocean_ext( k,j,i) - climate%PD_obs%T_ocean_ext( k,j,i))
      subclimate%S_ocean_corr_ext( k,j,i) = subclimate%S_ocean_ext( k,j,i) - (climate%GCM_PI%S_ocean_ext( k,j,i) - climate%PD_obs%S_ocean_ext( k,j,i))
    END DO
    END DO
    END DO
    CALL sync
  
  END SUBROUTINE correct_GCM_bias_ocean
  
! == Regrid 3-D ocean data fields in the vertical direction
  SUBROUTINE map_global_ocean_data_to_IMAUICE_vertical_grid( ocean_data)
    ! Map global 3-D ocean temperature and salinity from whatever vertical grid
    ! it has been provided on to the vertical grid used by IMAU-ICE
     
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_subocean_global),   INTENT(INOUT)   :: ocean_data
    
    ! Local variables:
    INTEGER,                    POINTER           :: nz_new
    REAL(dp), DIMENSION(:    ), POINTER           :: z_new
    REAL(dp), DIMENSION(:,:,:), POINTER           :: T_ocean_new, S_ocean_new
    INTEGER                                       :: wnz_new, wz_new, wT_ocean_new, wS_ocean_new
    INTEGER                                       :: i,j,k
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: z_mask_old, z_mask_new
    REAL(dp)                                      :: z_floor
    REAL(dp)                                      :: NaN
    
    ! A trick
    NaN = -1._dp
    NaN = SQRT( NaN)
    
    ! Set up the desired vertical grid
    CALL create_ocean_vertical_grid( nz_new, wnz_new, z_new, wz_new)
    
    ! Regrid 3-D ocean temperature and salinity
    CALL allocate_shared_dp_3D( ocean_data%nlon, ocean_data%nlat, nz_new, T_ocean_new, wT_ocean_new)
    CALL allocate_shared_dp_3D( ocean_data%nlon, ocean_data%nlat, nz_new, S_ocean_new, wS_ocean_new)
    
    ALLOCATE( z_mask_old( ocean_data%nz_ocean))
    ALLOCATE( z_mask_new( nz_new             ))
    
    DO i = ocean_data%i1, ocean_data%i2
    DO j = 1, ocean_data%nlat
      
      ! Determine local depth of the ocean floor, fill in both data masks
      IF (ocean_data%T_ocean( i,j,ocean_data%nz_ocean) == ocean_data%T_ocean( i,j,ocean_data%nz_ocean)) THEN
        ! Ocean floor lies below the vertical limit of the provided data
        z_mask_old = 1
        z_floor = ocean_data%z_ocean( ocean_data%nz_ocean) + (ocean_data%z_ocean( 2) - ocean_data%z_ocean( 1))
      ELSEIF (ocean_data%T_ocean( i,j,1) /= ocean_data%T_ocean( i,j,1)) THEN
        ! This grid cell isn't ocean at all
        z_mask_old = 0
        z_floor    = 0._dp
      ELSE
        z_mask_old = 1
        k = ocean_data%nz_ocean
        DO WHILE (ocean_data%T_ocean( i,j,k) /= ocean_data%T_ocean( i,j,k))
          z_mask_old( k) = 0
          z_floor = ocean_data%z_ocean( k)
          k = k - 1
        END DO
      END IF
      
      z_mask_new = 0
      DO k = 1, nz_new
        IF (z_new( k) < z_floor) z_mask_new = 1
      END DO
      
      ! Regrid vertical column
      CALL remap_cons_2nd_order_1D( ocean_data%z_ocean, z_mask_old, ocean_data%T_ocean( i,j,:), z_new, z_mask_new, T_ocean_new( i,j,:))
      CALL remap_cons_2nd_order_1D( ocean_data%z_ocean, z_mask_old, ocean_data%S_ocean( i,j,:), z_new, z_mask_new, S_ocean_new( i,j,:))
      
      ! Fill masked values with NaN
      DO k = 1, nz_new
        IF (z_mask_new( k) == 0) THEN
          T_ocean_new( i,j,k) = NaN
          S_ocean_new( i,j,k) = NaN
        END IF
      END DO
      
    END DO
    END DO
    CALL sync
    
    ! Deallocate old ocean data
    CALL deallocate_shared( ocean_data%wz_ocean )
    CALL deallocate_shared( ocean_data%wT_ocean )
    CALL deallocate_shared( ocean_data%wS_ocean )
    
    ! Allocate shared memory for regridded ocean data
    CALL allocate_shared_dp_1D(                                   nz_new, ocean_data%z_ocean, ocean_data%wz_ocean)  
    CALL allocate_shared_dp_3D( ocean_data%nlon, ocean_data%nlat, nz_new, ocean_data%T_ocean, ocean_data%wT_ocean)
    CALL allocate_shared_dp_3D( ocean_data%nlon, ocean_data%nlat, nz_new, ocean_data%S_ocean, ocean_data%wS_ocean)
    
    ! Copy regridded ocean data
    IF (par%master) THEN
      ocean_data%nz_ocean = nz_new
      ocean_data%z_ocean  = z_new
    END IF
    CALL sync
    ocean_data%T_ocean( ocean_data%i1:ocean_data%i2,:,:) = T_ocean_new( ocean_data%i1:ocean_data%i2,:,:)
    ocean_data%S_ocean( ocean_data%i1:ocean_data%i2,:,:) = S_ocean_new( ocean_data%i1:ocean_data%i2,:,:)
    
    ! Clean up after yourself
    CALL deallocate_shared( wnz_new)
    CALL deallocate_shared( wz_new )
    DEALLOCATE( z_mask_old)
    DEALLOCATE( z_mask_new)
    CALL deallocate_shared( wT_ocean_new)
    CALL deallocate_shared( wS_ocean_new)
    
  END SUBROUTINE map_global_ocean_data_to_IMAUICE_vertical_grid
  SUBROUTINE create_ocean_vertical_grid( nz, wnz, z, wz)
    ! Set up the vertical grid used for ocean data
     
    IMPLICIT NONE
      
    ! Input variables:
    INTEGER,                    POINTER, INTENT(OUT)   :: nz
    INTEGER,                             INTENT(OUT)   :: wnz
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   :: z
    INTEGER,                             INTENT(OUT)   :: wz
    
    IF (C%choice_ocean_vertical_grid == 'regular') THEN
      CALL create_ocean_vertical_grid_regular( nz, wnz, z, wz)
    ELSE
      IF (par%master) WRITE(0,*) '  create_ocean_vertical_grid - ERROR: choice_ocean_vertical_grid "', TRIM(C%choice_ocean_vertical_grid), '" not implemented!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE create_ocean_vertical_grid
  SUBROUTINE create_ocean_vertical_grid_regular( nz, wnz, z, wz)
    ! Set up the vertical grid used for ocean data - regular grid
     
    IMPLICIT NONE
      
    ! Input variables:
    INTEGER,                    POINTER, INTENT(OUT)   :: nz
    INTEGER,                             INTENT(OUT)   :: wnz
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   :: z
    INTEGER,                             INTENT(OUT)   :: wz
    
    ! Local variables:
    INTEGER                                       :: k
    
    ! Determine the number of vertical layers to be used
    CALL allocate_shared_int_0D( nz, wnz)
    IF (par%master) THEN
      nz = 1 + FLOOR( C%ocean_vertical_grid_max_depth / C%ocean_regular_grid_dz)
    END IF
    CALL sync
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( nz, z, wz)
    
    ! Fill in the values
    IF (par%master) THEN
      DO k = 1, nz
        z( k) = REAL(k-1,dp) * C%ocean_regular_grid_dz
      END DO
    END IF
    CALL sync
    
  END SUBROUTINE create_ocean_vertical_grid_regular
  
! == Some schematic ocean temperature/salinity profiles
  SUBROUTINE MISOMIP1_ocean_profiles( grid, climate, time)
    ! Set the ocean temperature and salinity according to the ISOMIP+ protocol
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    REAL(dp),                            INTENT(IN)    :: time
    
    IF     (C%MISOMIP1_scenario == 'IceOcean0') THEN
      ! Cold ocean always
      
      CALL set_ocean_to_ISOMIPplus_COLD( grid, climate)
      
    ELSEIF (C%MISOMIP1_scenario == 'IceOcean1ra' .OR. &
            C%MISOMIP1_scenario == 'IceOcean2ra') THEN
      ! Cold ocean during spin-up; warm ocean for 100 years, then cold ocean again
      
      IF (time < 0._dp) THEN
        CALL set_ocean_to_ISOMIPplus_COLD( grid, climate)
      ELSEIF (time < 100._dp) THEN
        CALL set_ocean_to_ISOMIPplus_WARM( grid, climate)
      ELSE
        CALL set_ocean_to_ISOMIPplus_COLD( grid, climate)
      END IF
      
    ELSEIF (C%MISOMIP1_scenario == 'IceOcean1rr' .OR. &
            C%MISOMIP1_scenario == 'IceOcean2rr') THEN
      ! Cold ocean during spin-up; warm ocean after t = 0
      
      IF (time < 0._dp) THEN
        CALL set_ocean_to_ISOMIPplus_COLD( grid, climate)
      ELSE
        CALL set_ocean_to_ISOMIPplus_WARM( grid, climate)
      END IF
      
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: MISOMIP1_scenario "', TRIM(C%MISOMIP1_scenario), '" not implemented in MISOMIP1_ocean_profiles!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
  
  END SUBROUTINE MISOMIP1_ocean_profiles
  SUBROUTINE set_ocean_to_ISOMIPplus_COLD( grid, climate)
    ! Set the ocean temperature and salinity to the ISOMIP+ "COLD" profile (Asay-Davis et al., 2016, Table 5)
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    
    ! Local variables
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: w
    REAL(dp), PARAMETER                                :: Tzero     = -1.9_dp    ! Sea surface temperature [degC] (originally T0, but that name is already taken...)
    REAL(dp), PARAMETER                                :: Tbot      = -1.9_dp    ! Sea floor   temperature [degC]
    REAL(dp), PARAMETER                                :: Szero     = 33.8_dp    ! Sea surface salinity    [PSU]
    REAL(dp), PARAMETER                                :: Sbot      = 34.55_dp   ! Sea floor   salinity    [PSU]
    REAL(dp), PARAMETER                                :: depth_max = 720._dp    ! Maximum depth for the profile (constant values below that)
    
    ! Allocate memory for the ocean data
    CALL allocate_subclimate_regional_oceans( grid, climate)
    
    ! Fill in the temperature and salinity profiles
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO k = 1, climate%nz_ocean
      
      ! Interpolation weight
      w = MIN(1._dp, MAX(0._dp, climate%z_ocean( k) / depth_max ))
      
      ! Temperature
      climate%T_ocean(          k,j,i) = w * Tbot + (1._dp - w) * Tzero
      climate%T_ocean_ext(      k,j,i) = w * Tbot + (1._dp - w) * Tzero
      climate%T_ocean_corr_ext( k,j,i) = w * Tbot + (1._dp - w) * Tzero
      
      ! Salinity
      climate%S_ocean(          k,j,i) = w * Sbot + (1._dp - w) * Szero
      climate%S_ocean_ext(      k,j,i) = w * Sbot + (1._dp - w) * Szero
      climate%S_ocean_corr_ext( k,j,i) = w * Sbot + (1._dp - w) * Szero
      
    END DO
    END DO
    END DO
    CALL sync
  
  END SUBROUTINE set_ocean_to_ISOMIPplus_COLD
  SUBROUTINE set_ocean_to_ISOMIPplus_WARM( grid, climate)
    ! Set the ocean temperature and salinity to the ISOMIP+ "WARM" profile (Asay-Davis et al., 2016, Table 6)
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    
    ! Local variables
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: w
    REAL(dp), PARAMETER                                :: Tzero     = -1.9_dp    ! Sea surface temperature [degC] (originally T0, but that name is already taken...)
    REAL(dp), PARAMETER                                :: Tbot      =  1.0_dp    ! Sea floor   temperature [degC]
    REAL(dp), PARAMETER                                :: Szero     = 33.8_dp    ! Sea surface salinity    [PSU]
    REAL(dp), PARAMETER                                :: Sbot      = 34.7_dp    ! Sea floor   salinity    [PSU]
    REAL(dp), PARAMETER                                :: depth_max = 720._dp    ! Maximum depth for the profile (constant values below that)
    
    ! Allocate memory for the ocean data
    CALL allocate_subclimate_regional_oceans( grid, climate)
    
    ! Fill in the temperature and salinity profiles
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    DO k = 1, climate%nz_ocean
      
      ! Interpolation weight
      w = MIN(1._dp, MAX(0._dp, climate%z_ocean( k) / depth_max ))
      
      ! Temperature
      climate%T_ocean(          k,j,i) = w * Tbot + (1._dp - w) * Tzero
      climate%T_ocean_ext(      k,j,i) = w * Tbot + (1._dp - w) * Tzero
      climate%T_ocean_corr_ext( k,j,i) = w * Tbot + (1._dp - w) * Tzero
      
      ! Salinity
      climate%S_ocean(          k,j,i) = w * Sbot + (1._dp - w) * Szero
      climate%S_ocean_ext(      k,j,i) = w * Sbot + (1._dp - w) * Szero
      climate%S_ocean_corr_ext( k,j,i) = w * Sbot + (1._dp - w) * Szero
      
    END DO
    END DO
    END DO
    CALL sync
  
  END SUBROUTINE set_ocean_to_ISOMIPplus_WARM
  SUBROUTINE set_ocean_to_Reese2018( grid, ice, climate)
    ! Set the ocean temperature and salinity to basin-dependent values
    ! provided by Reese et al. (2018) so that PICO gives realistic present-day melt rates
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    
    ! Local variables
    INTEGER                                            :: i,j
    REAL(dp)                                           :: T,S
    
    ! Safety
    IF (.NOT. (C%choice_basin_scheme_ANT == 'file' .AND. C%do_merge_basins_ANT)) THEN
      IF (par%master) THEN
        WRITE(0,*) ''
        WRITE(0,*) ' ===== '
        WRITE(0,*) 'set_ocean_to_Reese2018 - WARNING: This really only works when using the external Antarctic ice basins file "ant_full_drainagesystem_polygons.txt"'
        WRITE(0,*) '                                  This can be downloaded from https://earth.gsfc.nasa.gov/cryo/data/polar-altimetry/antarctic-and-greenland-drainage-systems'
        WRITE(0,*) '  and you will also need to set do_merge_basins_ANT_config = .TRUE.'
        WRITE(0,*) ' ===== '
        WRITE(0,*) ''
      END IF
    END IF
    
    ! Allocate memory for the ocean data
    CALL allocate_subclimate_regional_oceans( grid, climate)
    
    ! Fill in the temperature and salinity values
    T = 0._dp
    S = 0._dp
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      IF     (ice%basin_ID( j,i) == 1) THEN
        T = -1.76_dp
        S = 34.82_dp
      ELSEIF (ice%basin_ID( j,i) == 2) THEN
        T = -1.66_dp
        S = 34.70_dp
      ELSEIF (ice%basin_ID( j,i) == 3) THEN
        T = -1.65_dp
        S = 34.48_dp
      ELSEIF (ice%basin_ID( j,i) == 4) THEN
        T = -1.58_dp
        S = 34.49_dp
      ELSEIF (ice%basin_ID( j,i) == 5) THEN
        T = -1.51_dp
        S = 34.50_dp
      ELSEIF (ice%basin_ID( j,i) == 6) THEN
        T = -1.73_dp
        S = 34.70_dp
      ELSEIF (ice%basin_ID( j,i) == 7) THEN
        T = -1.68_dp
        S = 34.65_dp
      ELSEIF (ice%basin_ID( j,i) == 8) THEN
        T = -0.73_dp
        S = 34.73_dp
      ELSEIF (ice%basin_ID( j,i) == 9) THEN
        T = -1.61_dp
        S = 34.75_dp
      ELSEIF (ice%basin_ID( j,i) == 10) THEN
        T = -1.30_dp
        S = 34.84_dp
      ELSEIF (ice%basin_ID( j,i) == 11) THEN
        T = -1.58_dp
        S = 34.79_dp
      ELSEIF (ice%basin_ID( j,i) == 12) THEN
        T = -0.36_dp
        S = 34.58_dp
      ELSEIF (ice%basin_ID( j,i) == 13) THEN
        T =  0.80_dp
        S = 34.79_dp
      ELSEIF (ice%basin_ID( j,i) == 14) THEN
        T =  1.10_dp
        S = 34.85_dp
      ELSEIF (ice%basin_ID( j,i) == 15) THEN
        T =  0.23_dp
        S = 34.7_dp
      ELSEIF (ice%basin_ID( j,i) == 16) THEN
        T = -1.23_dp
        S = 34.67_dp
      ELSEIF (ice%basin_ID( j,i) == 17) THEN
        T = -1.80_dp
        S = 34.84_dp
      END IF
      
      ! Temperature
      climate%T_ocean(          :,j,i) = T
      climate%T_ocean_ext(      :,j,i) = T
      climate%T_ocean_corr_ext( :,j,i) = T
      
      ! Salinity
      climate%S_ocean(          :,j,i) = S
      climate%S_ocean_ext(      :,j,i) = S
      climate%S_ocean_corr_ext( :,j,i) = S
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE set_ocean_to_Reese2018

END MODULE ocean_module
