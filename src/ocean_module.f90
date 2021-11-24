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
                                             type_climate_model, type_subclimate_region, type_init_data_fields, &
                                             type_model_region, type_highres_ocean_data
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, &
                                             inquire_PD_obs_data_file_ocean, read_PD_obs_data_file_ocean, &
                                             inquire_GCM_ocean_snapshot, read_GCM_ocean_snapshot, &
                                             inquire_hires_geometry_file, read_hires_geometry_file, &
                                             create_extrapolated_ocean_file, inquire_extrapolated_ocean_file, &
                                             read_extrapolated_ocean_file
  USE forcing_module,                  ONLY: forcing
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             error_function, smooth_Gaussian_2D, smooth_Gaussian_3D, smooth_Shepard_2D, &
                                             inverse_oblique_sg_projection, map_glob_to_grid_2D, map_glob_to_grid_3D, &
                                             map_square_to_square_cons_2nd_order_2D, map_square_to_square_cons_2nd_order_3D, &
                                             remap_cons_2nd_order_1D, surface_elevation, extrapolate_Gaussian_floodfill, &
                                             transpose_dp_2D, transpose_dp_3D

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
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: a
    REAL(dp)                                           :: CO2
    REAL(dp)                                           :: w_CO2
    REAL(dp), DIMENSION(:,:  ), POINTER                :: w_ins, w_ins_smooth,  w_ice,  w_tot
    REAL(dp), DIMENSION(:,:,:), POINTER                :: w_tot_final
    INTEGER                                            :: ww_ins, ww_ins_smooth, ww_ice, ww_tot, ww_tot_final
    REAL(dp)                                           :: w_ins_av
    REAL(dp), PARAMETER                                :: w_cutoff = 0.25_dp        ! Crop weights to [-w_cutoff, 1 + w_cutoff]


    CALL allocate_shared_dp_2D(                          grid%ny, grid%nx, w_ins,        ww_ins         )
    CALL allocate_shared_dp_2D(                          grid%ny, grid%nx, w_ins_smooth, ww_ins_smooth  )
    CALL allocate_shared_dp_2D(                          grid%ny, grid%nx, w_ice,        ww_ice         )
    CALL allocate_shared_dp_2D(                          grid%ny, grid%nx, w_tot,        ww_tot         )
    CALL allocate_shared_dp_3D(climate%applied%nz_ocean, grid%ny, grid%nx, w_tot_final,  ww_tot_final   )
      
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
    !w_ice = 1._dp - MAX(-w_cutoff, MIN(1._dp + w_cutoff, (SUM(ice%Hs_a) - SUM(climate%GCM_warm%Hs)) / (SUM(climate%GCM_cold%Hs) - SUM(climate%GCM_warm%Hs)) ))
    
    ! Calculate weighting field
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      w_ins( j,i) = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (    climate%applied%I_abs(  j,i) -     climate%GCM_cold%I_abs( j,i)) / &  ! Berends et al., 2018 - Eq. 3
                                                           (    climate%GCM_warm%I_abs( j,i) -     climate%GCM_cold%I_abs( j,i)) ))
    END DO
    END DO
    CALL sync
    w_ins_av      = MAX( -w_cutoff, MIN( 1._dp + w_cutoff, (SUM(climate%applied%I_abs )      - SUM(climate%GCM_cold%I_abs)     ) / &
                                                           (SUM(climate%GCM_warm%I_abs)      - SUM(climate%GCM_cold%I_abs)     ) ))
   
    ! Smooth the weighting field
    w_ins_smooth( :,grid%i1:grid%i2) = w_ins( :,grid%i1:grid%i2)
    CALL smooth_Gaussian_2D( grid, w_ins_smooth, 200000._dp)
    
    ! Combine unsmoothed, smoothed, and regional average weighting fields (Berends et al., 2018, Eq. 4)
    w_ice( :,grid%i1:grid%i2) = (1._dp * w_ins_smooth( :,grid%i1:grid%i2) + 6._dp * w_ins_av) / 7._dp
    
    ! Combine weigths CO2 and ice
    ! ===========================
    
    IF         (region_name == 'NAM') THEN
      w_tot( :,grid%i1:grid%i2) = (C%ocean_matrix_CO2vsice_NAM * w_CO2) + ((1._dp - C%ocean_matrix_CO2vsice_NAM) * w_ice( :,grid%i1:grid%i2)) 
    ELSEIF     (region_name == 'EAS') THEN
      w_tot( :,grid%i1:grid%i2) = (C%ocean_matrix_CO2vsice_EAS * w_CO2) + ((1._dp - C%ocean_matrix_CO2vsice_EAS) * w_ice( :,grid%i1:grid%i2)) 
    ELSEIF     (region_name == 'GRL') THEN
      w_tot( :,grid%i1:grid%i2) = (C%ocean_matrix_CO2vsice_GRL * w_CO2) + ((1._dp - C%ocean_matrix_CO2vsice_GRL) * w_ice( :,grid%i1:grid%i2)) 
    ELSEIF     (region_name == 'ANT') THEN
      w_tot( :,grid%i1:grid%i2) = (C%ocean_matrix_CO2vsice_ANT * w_CO2) + ((1._dp - C%ocean_matrix_CO2vsice_ANT) * w_ice( :,grid%i1:grid%i2)) 
    END IF
    
    ! Update the history of the weighing fields
    ! =========================================
    
    ! 1st entry is the current value, 2nd is 1*dt_ocean ago, 3d is 2*dt_ocean ago, etc.
    climate%applied%w_tot_history( 2:climate%applied%nw_tot_history,:,grid%i1:grid%i2) = climate%applied%w_tot_history( 1:climate%applied%nw_tot_history-1,:,grid%i1:grid%i2)
    climate%applied%w_tot_history( 1, :,grid%i1:grid%i2) = w_tot( :,grid%i1:grid%i2)
    
    ! Interpolate the GCM ocean snapshots
    ! =============================
    
    DO k = 1,climate%applied%nz_ocean
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      a = ( ( climate%applied%z_ocean (k) / climate%applied%z_ocean (climate%applied%nz_ocean) ) * (climate%applied%nw_tot_history-1) ) + 1

      w_tot_final (k,j,i) = ( SUM(climate%applied%w_tot_history(1:FLOOR(a),j,i)) + ( (a - FLOOR(a)) * climate%applied%w_tot_history(CEILING(a),j,i) ) ) * (1._dp / a)
      
      climate%applied%T_ocean_corr_ext (k,j,i) = (           w_tot_final(k,j,i)  * climate%GCM_warm%T_ocean_corr_ext (k,j,i)  ) + &
                                                 (  (1._dp - w_tot_final(k,j,i) )* climate%GCM_cold%T_ocean_corr_ext (k,j,i)  )
      climate%applied%S_ocean_corr_ext (k,j,i) = (           w_tot_final(k,j,i)  * climate%GCM_warm%S_ocean_corr_ext (k,j,i)  ) + &
                                                 (  (1._dp - w_tot_final(k,j,i) )* climate%GCM_cold%S_ocean_corr_ext (k,j,i)  )     
    END DO
    END DO
    END DO
    CALL sync
    
    CALL deallocate_shared( ww_ins)
    CALL deallocate_shared( ww_ins_smooth)
    CALL deallocate_shared( ww_ice)
    CALL deallocate_shared( ww_tot)  
    CALL deallocate_shared( ww_tot_final)  
      
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
  
! == Initialising the region-specific ocean model
  SUBROUTINE initialise_ocean_model( region, matrix)
    ! Allocate shared memory for the ocean data of the regional subclimates
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_climate_matrix),           INTENT(IN)    :: matrix
        
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
    
    ! Exception for schematic ocean temperature/salinity profiles
    ! ===========================================================
    
    IF (C%choice_ocean_temperature_model == 'schematic') THEN
    
      IF (par%master) WRITE(*,*) '    Schematic"', TRIM(C%choice_schematic_ocean), '" ocean forcing used for BMB forcing!'
      
      IF     (C%choice_schematic_ocean == 'MISMIPplus_WARM') THEN
        CALL set_ocean_to_ISOMIPplus_WARM( region%grid, region%climate%PD_obs  )
        CALL set_ocean_to_ISOMIPplus_WARM( region%grid, region%climate%GCM_PI  )
        CALL set_ocean_to_ISOMIPplus_WARM( region%grid, region%climate%GCM_cold)
        CALL set_ocean_to_ISOMIPplus_WARM( region%grid, region%climate%GCM_warm)
        CALL set_ocean_to_ISOMIPplus_WARM( region%grid, region%climate%applied )
      ELSEIF (C%choice_schematic_ocean == 'MISMIPplus_COLD') THEN
        CALL set_ocean_to_ISOMIPplus_COLD( region%grid, region%climate%PD_obs  )
        CALL set_ocean_to_ISOMIPplus_COLD( region%grid, region%climate%GCM_PI  )
        CALL set_ocean_to_ISOMIPplus_COLD( region%grid, region%climate%GCM_cold)
        CALL set_ocean_to_ISOMIPplus_COLD( region%grid, region%climate%GCM_warm)
        CALL set_ocean_to_ISOMIPplus_COLD( region%grid, region%climate%applied )
      ELSEIF (C%choice_schematic_ocean == 'Reese2018') THEN
        CALL set_ocean_to_Reese2018(       region%grid, region%ice, region%climate%PD_obs  )
        CALL set_ocean_to_Reese2018(       region%grid, region%ice, region%climate%GCM_PI  )
        CALL set_ocean_to_Reese2018(       region%grid, region%ice, region%climate%GCM_cold)
        CALL set_ocean_to_Reese2018(       region%grid, region%ice, region%climate%GCM_warm)
        CALL set_ocean_to_Reese2018(       region%grid, region%ice, region%climate%applied )
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_schematic_ocean "', TRIM(C%choice_schematic_ocean), '" not implemented in initialise_subclimate!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
      RETURN  
      
    END IF ! IF (C%choice_ocean_temperature_model == 'schematic') THEN
    
    
      ! Allocate memory for regional ocean data 
      CALL allocate_subclimate_regional_oceans( region%grid, region%climate%PD_obs )
      CALL allocate_subclimate_regional_oceans( region%grid, region%climate%applied)
      
      ! Map ocean data from the global lon/lat-grid to the high-resolution regional x/y-grid,
      ! extrapolate mapped ocean data to cover the entire 3D domain, and finally
      ! map to the actual ice model resolution
      CALL get_extrapolated_ocean_data( region, matrix%PD_obs_ocean, region%climate%PD_obs, C%filename_PD_obs_ocean)

    IF (C%choice_ocean_temperature_model == 'WOA') THEN
      IF (par%master) WRITE(*,*) '    Constant present-day ocean forcing used for BMB forcing!'      
      ! PD_obs doesn't have a bias-corrected version
      region%climate%PD_obs%T_ocean_corr_ext  = region%climate%PD_obs%T_ocean_ext
      region%climate%PD_obs%S_ocean_corr_ext  = region%climate%PD_obs%S_ocean_ext
      
    ELSE IF (C%choice_ocean_temperature_model == 'fixed' .OR. C%choice_ocean_temperature_model == 'scaled')  THEN
      IF (par%master) WRITE(*,*) '    Ocean forcing initialised, but not actually used for BMB forcing!'
      
      ! Allocate memory for regional ocean data 
      CALL allocate_subclimate_regional_oceans( region%grid, region%climate%applied)
      
    ELSE IF (C%choice_ocean_temperature_model == 'matrix_warm_cold') THEN 
    
      CALL allocate_subclimate_regional_oceans( region%grid, region%climate%GCM_PI  )
      CALL allocate_subclimate_regional_oceans( region%grid, region%climate%GCM_cold)
      CALL allocate_subclimate_regional_oceans( region%grid, region%climate%GCM_warm)
      
      ! Map ocean data from the global lon/lat-grid to the high-resolution regional x/y-grid,
      ! extrapolate mapped ocean data to cover the entire 3D domain, and finally
      ! map to the actual ice model resolution
      CALL get_extrapolated_ocean_data( region, matrix%PD_obs_ocean,   region%climate%PD_obs,   C%filename_PD_obs_ocean           )
      CALL get_extrapolated_ocean_data( region, matrix%GCM_PI_ocean,   region%climate%GCM_PI,   C%filename_GCM_ocean_snapshot_PI  )
      CALL get_extrapolated_ocean_data( region, matrix%GCM_cold_ocean, region%climate%GCM_cold, C%filename_GCM_ocean_snapshot_cold)
      CALL get_extrapolated_ocean_data( region, matrix%GCM_warm_ocean, region%climate%GCM_warm, C%filename_GCM_ocean_snapshot_warm)
    
      ! Correct regional ocean data for GCM bias
      CALL correct_GCM_bias_ocean( region%grid, region%climate, region%climate%GCM_warm)
      CALL correct_GCM_bias_ocean( region%grid, region%climate, region%climate%GCM_cold)
      
      ! PI and PD_obs don't have a bias-corrected version
      region%climate%GCM_PI%T_ocean_corr_ext  = region%climate%GCM_PI%T_ocean_ext
      region%climate%GCM_PI%S_ocean_corr_ext  = region%climate%GCM_PI%S_ocean_ext
      region%climate%PD_obs%T_ocean_corr_ext  = region%climate%PD_obs%T_ocean_ext
      region%climate%PD_obs%S_ocean_corr_ext  = region%climate%PD_obs%S_ocean_ext
      
      ! Allocate memory for the weighing fields history, and initialise      
      CALL allocate_shared_int_0D ( region%climate%applied%nw_tot_history, region%climate%applied%wnw_tot_history)
      region%climate%applied%nw_tot_history = CEILING(C%w_tot_hist_averaging_window/C%dt_ocean)+1
      CALL allocate_shared_dp_3D  ( region%climate%applied%nw_tot_history, region%grid%ny, region%grid%nx, region%climate%applied%w_tot_history, region%climate%applied%ww_tot_history)
      region%climate%applied%w_tot_history = 0._dp ! Initiate at cold conditions
            
    ELSE  
      IF (par%master) WRITE(0,*) '  ERROR: choice_ocean_temperature_model "', TRIM(C%choice_ocean_temperature_model), '" not implemented in initialise_oceans_regional!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF !(C%choice_ocean_temperature_model == 'WOA')
 
    ! Initialise applied ocean forcing with present-day observations
    IF (par%master) THEN
      region%climate%applied%T_ocean          = region%climate%PD_obs%T_ocean
      region%climate%applied%T_ocean_ext      = region%climate%PD_obs%T_ocean_ext
      region%climate%applied%T_ocean_corr_ext = region%climate%PD_obs%T_ocean_corr_ext
      region%climate%applied%S_ocean          = region%climate%PD_obs%S_ocean
      region%climate%applied%S_ocean_ext      = region%climate%PD_obs%S_ocean_ext
      region%climate%applied%S_ocean_corr_ext = region%climate%PD_obs%S_ocean_corr_ext
    END IF ! IF (par%master) THEN
    CALL sync
  
  END SUBROUTINE initialise_ocean_model
  SUBROUTINE allocate_subclimate_regional_oceans( grid, climate)
    ! Allocate shared memory for the ocean data of a regional subclimate models
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    
    ! Set up the vertical coordinate
    CALL create_ocean_vertical_grid( climate%nz_ocean, climate%wnz_ocean, climate%z_ocean, climate%wz_ocean)
    
    ! Allocate shared memory
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%T_ocean,          climate%wT_ocean         )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%S_ocean,          climate%wS_ocean         )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%T_ocean_ext,      climate%wT_ocean_ext     )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%S_ocean_ext,      climate%wS_ocean_ext     )
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%T_ocean_corr_ext, climate%wT_ocean_corr_ext)
    CALL allocate_shared_dp_3D(  climate%nz_ocean, grid%ny, grid%nx, climate%S_ocean_corr_ext, climate%wS_ocean_corr_ext)
  
  END SUBROUTINE allocate_subclimate_regional_oceans
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

! == Extrapolating incomplete ocean data into the complete 3-D model domain
  SUBROUTINE get_extrapolated_ocean_data( region, ocean_glob, ocean_reg, filename_ocean_glob)
    ! Check if extrapolated ocean files for the current ice model setting exist. If so,
    ! read those. If not, perform the extrapolation and save the results to a new netCDF file.
        
    ! When creating a set of extrapolated files, a header file is created that describes
    ! the ice model settings for which those files were created. We check all existing header
    ! files, if any of them match the current settings, we read the extrapolated files listed
    ! there. If none of them match, we create a set of extrapolated files (and the accompanying
    ! header file) from scratch.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_subocean_global),          INTENT(IN)    :: ocean_glob
    TYPE(type_subclimate_region),        INTENT(INOUT) :: ocean_reg
    CHARACTER(LEN=256),                  INTENT(IN)    :: filename_ocean_glob
    
    ! Local variables:
    LOGICAL                                            :: foundmatch
    CHARACTER(LEN=256)                                 :: hires_ocean_foldername
    TYPE(type_highres_ocean_data)                      :: hires
    
  ! Initialise the high-resolution extrapolated ocean data
  ! ======================================================
    
    ! If a valid preprocessed file exists, read data from there. If not, perform
    ! the preprocessing and save the result to a file to save on future work
    
    ! First, check if any existing header matches the current ice model set-up.  
    CALL check_for_matching_ocean_header( region, ocean_reg, filename_ocean_glob, foundmatch, hires_ocean_foldername)
    
    IF (foundmatch) THEN
      IF (par%master) WRITE(0,*) '   Found valid extrapolated ocean data in folder "', TRIM( hires_ocean_foldername), '"'
      CALL get_hires_ocean_data_from_file( region, hires, hires_ocean_foldername)
    ELSE
      ! No header fitting the current ice model set-up was found. Create a new one describing
      ! the current set-up, and generate extrapolated ocean data files from scratch.
      IF (par%master) WRITE(0,*) '   Creating new extrapolated ocean data in folder "', TRIM( hires_ocean_foldername), '"'
      CALL map_and_extrapolate_hires_ocean_data( region, ocean_glob, hires)
      CALL write_hires_extrapolated_ocean_data_to_file( hires, filename_ocean_glob, hires_ocean_foldername)
    END IF ! IF (.NOT. foundmatch) THEN
    
  ! ===== Map extrapolated data from the high-resolution grid to the actual ice-model grid =====
  ! ============================================================================================
    
    IF (par%master) WRITE(0,*) '   Mapping high-resolution extrapolated ocean data to the ice-model grid...'
    CALL map_square_to_square_cons_2nd_order_3D( hires%grid%nx, hires%grid%ny, hires%grid%x, hires%grid%y, region%grid%nx, region%grid%ny, region%grid%x, region%grid%y, hires%T_ocean, ocean_reg%T_ocean_ext, hires%nz_ocean)
    CALL map_square_to_square_cons_2nd_order_3D( hires%grid%nx, hires%grid%ny, hires%grid%x, hires%grid%y, region%grid%nx, region%grid%ny, region%grid%x, region%grid%y, hires%S_ocean, ocean_reg%S_ocean_ext, hires%nz_ocean)
    
    ! Clean up after yourself
    CALL deallocate_shared( hires%grid%wnx          )
    CALL deallocate_shared( hires%grid%wny          )
    CALL deallocate_shared( hires%grid%wx           )
    CALL deallocate_shared( hires%grid%wy           )
    CALL deallocate_shared( hires%grid%wxmin        )
    CALL deallocate_shared( hires%grid%wxmax        )
    CALL deallocate_shared( hires%grid%wymin        )
    CALL deallocate_shared( hires%grid%wymax        )
    CALL deallocate_shared( hires%grid%wdx          )
    CALL deallocate_shared( hires%grid%wlambda_m    )
    CALL deallocate_shared( hires%grid%wphi_m       )
    CALL deallocate_shared( hires%grid%walpha_stereo)
    CALL deallocate_shared( hires%grid%wlat         )
    CALL deallocate_shared( hires%grid%wlon         )
    CALL deallocate_shared( hires%wnz_ocean         )
    CALL deallocate_shared( hires%wz_ocean          )
    CALL deallocate_shared( hires%wT_ocean          )
    CALL deallocate_shared( hires%wS_ocean          )
    
  END SUBROUTINE get_extrapolated_ocean_data
  SUBROUTINE get_hires_ocean_data_from_file( region, hires, hires_ocean_foldername)
    ! Read high-resolution extrapolated ocean data from an external file

    ! In/output variables:
    TYPE(type_model_region),             INTENT(IN)    :: region
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires
    CHARACTER(LEN=256),                  INTENT(IN)    :: hires_ocean_foldername
    
    ! Local variables
    INTEGER                                            :: i,j
    
    ! Check if the NetCDF file has all the required dimensions and variables
    CALL allocate_shared_int_0D( hires%grid%nx,  hires%grid%wnx)
    CALL allocate_shared_int_0D( hires%grid%ny,  hires%grid%wny)
    CALL allocate_shared_int_0D( hires%nz_ocean, hires%wnz_ocean)
    IF (par%master) THEN
      hires%netcdf%filename = TRIM( hires_ocean_foldername)//'/extrapolated_ocean_data.nc'
      CALL inquire_extrapolated_ocean_file( hires)
    END IF
    CALL sync
    
    ! Allocate shared memory for x,y and the actual data
    CALL allocate_shared_dp_1D( hires%grid%nx,                                hires%grid%x,  hires%grid%wx )
    CALL allocate_shared_dp_1D(                hires%grid%ny,                 hires%grid%y,  hires%grid%wy )
    CALL allocate_shared_dp_1D(                               hires%nz_ocean, hires%z_ocean, hires%wz_ocean)
    CALL allocate_shared_dp_3D( hires%grid%nx, hires%grid%ny, hires%nz_ocean, hires%T_ocean, hires%wT_ocean)
    CALL allocate_shared_dp_3D( hires%grid%nx, hires%grid%ny, hires%nz_ocean, hires%S_ocean, hires%wS_ocean)
    
    ! Read the data from the NetCDF file
    IF (par%master) THEN
      WRITE(0,*) '    Reading high-resolution extrapolated ocean data from file "', TRIM( hires%netcdf%filename), '"...'
      CALL read_extrapolated_ocean_file( hires)
    END IF
    CALL sync
    
    ! Transpose the data (since the file is [i,j] while the model is [j,i])
    CALL transpose_dp_3D( hires%T_ocean, hires%wT_ocean)
    CALL transpose_dp_3D( hires%S_ocean, hires%wS_ocean)
    
    ! Allocate shared memory for other grid parameters
    CALL allocate_shared_dp_0D(  hires%grid%dx,           hires%grid%wdx          )
    CALL allocate_shared_dp_0D(  hires%grid%xmin,         hires%grid%wxmin        )
    CALL allocate_shared_dp_0D(  hires%grid%xmax,         hires%grid%wxmax        )
    CALL allocate_shared_dp_0D(  hires%grid%ymin,         hires%grid%wymin        )
    CALL allocate_shared_dp_0D(  hires%grid%ymax,         hires%grid%wymax        )
    CALL allocate_shared_dp_0D(  hires%grid%lambda_M,     hires%grid%wlambda_M    )
    CALL allocate_shared_dp_0D(  hires%grid%phi_M,        hires%grid%wphi_M       )
    CALL allocate_shared_dp_0D(  hires%grid%alpha_stereo, hires%grid%walpha_stereo)
    
    ! Polar stereographic projection parameters and resolution
    IF (par%master) THEN
      ! Projection parameters are of course identical to those used for this ice model region
      hires%grid%lambda_M     = region%grid%lambda_M
      hires%grid%phi_M        = region%grid%phi_M
      hires%grid%alpha_stereo = region%grid%alpha_stereo
      ! But the resolution is different
      hires%grid%dx           = hires%grid%x( 2) - hires%grid%x( 1)
    END IF
    CALL sync
    
    ! Assign range to each processor
    CALL partition_list( hires%grid%nx, par%i, par%n, hires%grid%i1, hires%grid%i2)
    CALL partition_list( hires%grid%ny, par%i, par%n, hires%grid%j1, hires%grid%j2)
    
    ! Lat,lon coordinates
    CALL allocate_shared_dp_2D( hires%grid%ny, hires%grid%nx, hires%grid%lat, hires%grid%wlat)
    CALL allocate_shared_dp_2D( hires%grid%ny, hires%grid%nx, hires%grid%lon, hires%grid%wlon)
    
    DO i = hires%grid%i1, hires%grid%i2
    DO j = 1, hires%grid%ny
      CALL inverse_oblique_sg_projection( hires%grid%x( i), hires%grid%y( j), hires%grid%lambda_M, hires%grid%phi_M, hires%grid%alpha_stereo, hires%grid%lon( j,i), hires%grid%lat( j,i))
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE get_hires_ocean_data_from_file
  SUBROUTINE map_and_extrapolate_hires_ocean_data( region, ocean_glob, hires)
    ! Map ocean data from the global lon/lat-grid to the high-resolution regional x/y-grid,
    ! extrapolate mapped ocean data to cover the entire 3D domain, and finally
    ! map to the actual ice model resolution.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_model_region),             INTENT(IN)    :: region
    TYPE(type_subocean_global),          INTENT(IN)    :: ocean_glob
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  basin_ID_dp_lores,  basin_ID_dp_hires,  basin_ID_dp_hires_ext
    INTEGER                                            :: wbasin_ID_dp_lores, wbasin_ID_dp_hires, wbasin_ID_dp_hires_ext
    INTEGER                                            :: ii,jj,n
    LOGICAL                                            :: foundit
  
  ! ===== Read the high-resolution geometry data and generate the hi-res grid based on that=====
  ! ============================================================================================
    
    ! Determine which file to use for this region
    IF     (region%name == 'NAM') THEN
      hires%netcdf_geo%filename = C%ocean_extrap_hires_geo_filename_NAM
    ELSEIF (region%name == 'EAS') THEN
      hires%netcdf_geo%filename = C%ocean_extrap_hires_geo_filename_EAS
    ELSEIF (region%name == 'GRL') THEN
      hires%netcdf_geo%filename = C%ocean_extrap_hires_geo_filename_GRL
    ELSEIF (region%name == 'ANT') THEN
      hires%netcdf_geo%filename = C%ocean_extrap_hires_geo_filename_ANT
    END IF
    
    ! Check if the NetCDF file has all the required dimensions and variables
    CALL allocate_shared_int_0D( hires%grid%nx, hires%grid%wnx)
    CALL allocate_shared_int_0D( hires%grid%ny, hires%grid%wny)
    IF (par%master) THEN
      CALL inquire_hires_geometry_file( hires)
    END IF
    CALL sync
    
    ! Allocate shared memory for x,y and the actual data
    CALL allocate_shared_dp_1D( hires%grid%nx,                hires%grid%x, hires%grid%wx)
    CALL allocate_shared_dp_1D(                hires%grid%ny, hires%grid%y, hires%grid%wy)
    CALL allocate_shared_dp_2D( hires%grid%nx, hires%grid%ny, hires%Hi,     hires%wHi    )
    CALL allocate_shared_dp_2D( hires%grid%nx, hires%grid%ny, hires%Hb,     hires%wHb    )
    
    ! Read the data from the NetCDF file
    IF (par%master) THEN
      WRITE(0,*) '    Reading high-resolution geometry for ocean extrapolation from file "', TRIM( hires%netcdf_geo%filename), '"...'
      CALL read_hires_geometry_file( hires)
    END IF
    CALL sync
    
    ! Transpose the data (since the file is [i,j] while the model is [j,i])
    CALL transpose_dp_2D( hires%Hi, hires%wHi)
    CALL transpose_dp_2D( hires%Hb, hires%wHb)
    
    ! Allocate shared memory for other grid parameters
    CALL allocate_shared_dp_0D(  hires%grid%dx,           hires%grid%wdx          )
    CALL allocate_shared_dp_0D(  hires%grid%xmin,         hires%grid%wxmin        )
    CALL allocate_shared_dp_0D(  hires%grid%xmax,         hires%grid%wxmax        )
    CALL allocate_shared_dp_0D(  hires%grid%ymin,         hires%grid%wymin        )
    CALL allocate_shared_dp_0D(  hires%grid%ymax,         hires%grid%wymax        )
    CALL allocate_shared_dp_0D(  hires%grid%lambda_M,     hires%grid%wlambda_M    )
    CALL allocate_shared_dp_0D(  hires%grid%phi_M,        hires%grid%wphi_M       )
    CALL allocate_shared_dp_0D(  hires%grid%alpha_stereo, hires%grid%walpha_stereo)
    
    ! Polar stereographic projection parameters and resolution
    IF (par%master) THEN
      ! Projection parameters are of course identical to those used for this ice model region
      hires%grid%lambda_M     = region%grid%lambda_M
      hires%grid%phi_M        = region%grid%phi_M
      hires%grid%alpha_stereo = region%grid%alpha_stereo
      ! But the resolution is different
      hires%grid%dx           = hires%grid%x( 2) - hires%grid%x( 1)
      ! Check if this is the resolution we want
      IF (hires%grid%dx /= C%ocean_extrap_res) THEN
        WRITE(0,*) '  map_and_extrapolate_ocean_data - ERROR: high-resolution geometry file "', TRIM(hires%netcdf%filename), '" has a different resolution from C%ocean_extrap_res = ', C%ocean_extrap_res
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF
    CALL sync
    
    ! Assign range to each processor
    CALL partition_list( hires%grid%nx, par%i, par%n, hires%grid%i1, hires%grid%i2)
    CALL partition_list( hires%grid%ny, par%i, par%n, hires%grid%j1, hires%grid%j2)
    
    ! Lat,lon coordinates
    CALL allocate_shared_dp_2D( hires%grid%ny, hires%grid%nx, hires%grid%lat, hires%grid%wlat)
    CALL allocate_shared_dp_2D( hires%grid%ny, hires%grid%nx, hires%grid%lon, hires%grid%wlon)
    
    DO i = hires%grid%i1, hires%grid%i2
    DO j = 1, hires%grid%ny
      CALL inverse_oblique_sg_projection( hires%grid%x( i), hires%grid%y( j), hires%grid%lambda_M, hires%grid%phi_M, hires%grid%alpha_stereo, hires%grid%lon( j,i), hires%grid%lat( j,i))
    END DO
    END DO
    CALL sync
    
  ! ===== Map ocean data from the global lon/lat-grid to the high-resolution regional x/y-grid =====
  ! ================================================================================================
    
    IF (par%master) WRITE(0,'(A,F4.1,A)') '     Mapping ocean data from the global lat/lon-grid to the ', hires%grid%dx / 1000._dp, ' km regional x/y-grid...'
    
    ! Allocate shared memory for high-resolution ocean data
    CALL allocate_shared_dp_3D( ocean_glob%nz_ocean, hires%grid%ny, hires%grid%nx, hires%T_ocean, hires%wT_ocean)
    CALL allocate_shared_dp_3D( ocean_glob%nz_ocean, hires%grid%ny, hires%grid%nx, hires%S_ocean, hires%wS_ocean)
    
    ! Allocate memory for and copy vertical ocean grid
    CALL allocate_shared_int_0D(                      hires%nz_ocean, hires%wnz_ocean)
    CALL allocate_shared_dp_1D(  ocean_glob%nz_ocean, hires%z_ocean,  hires%wz_ocean )
    IF (par%master) THEN
      hires%nz_ocean = ocean_glob%nz_ocean
      hires%z_ocean  = ocean_glob%z_ocean
    END IF
    CALL sync
    
    ! Map the data from the global lon/lat-grid to the high-resolution regional x/y-grid
    CALL map_glob_to_grid_3D( ocean_glob%nlat, ocean_glob%nlon, ocean_glob%lat, ocean_glob%lon, hires%grid, ocean_glob%T_ocean, hires%T_ocean, hires%nz_ocean)
    CALL map_glob_to_grid_3D( ocean_glob%nlat, ocean_glob%nlon, ocean_glob%lat, ocean_glob%lon, hires%grid, ocean_glob%S_ocean, hires%S_ocean, hires%nz_ocean)
  
  ! ===== Perform the extrapolation on the high-resolution grid =====
  ! =================================================================
    
    IF (par%master) WRITE(0,'(A,F4.1,A)') '     Defining ice basins on the ', hires%grid%dx / 1000._dp, ' km regional x/y-grid...'
    
    ! Allocate shared memory for ice basins on the high-resolution grid
    CALL allocate_shared_int_2D( hires%grid%ny, hires%grid%nx, hires%basin_ID, hires%wbasin_ID)
    CALL allocate_shared_int_0D(                               hires%nbasins,  hires%wnbasins )
    
    IF (par%master) hires%nbasins = region%ice%nbasins
    CALL sync
    
    ! Instead of doing the "proper" basin definition on high resolution (which is insanely slow),
    ! just downscale the basin ID field from the ice model (using some tricks to get accurate values near the boundaries)
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( region%grid%ny, region%grid%nx, basin_ID_dp_lores,     wbasin_ID_dp_lores    )
    CALL allocate_shared_dp_2D( hires%grid%ny,  hires%grid%nx,  basin_ID_dp_hires,     wbasin_ID_dp_hires    )
    CALL allocate_shared_dp_2D( hires%grid%ny,  hires%grid%nx,  basin_ID_dp_hires_ext, wbasin_ID_dp_hires_ext)
    
    ! Convert basin ID field to double precision (for remapping)
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
      basin_ID_dp_lores( j,i) = REAL( region%ice%basin_ID( j,i), dp)
    END DO
    END DO
    CALL sync
    
    ! Map double-precision basin ID from ice-model grid to high-resolution grid
    CALL map_square_to_square_cons_2nd_order_2D( region%grid%nx, region%grid%ny, region%grid%x, region%grid%y, hires%grid%nx, hires%grid%ny, hires%grid%x, hires%grid%y, basin_ID_dp_lores, basin_ID_dp_hires)
    
    ! Remove all near-boundary cells
    DO i = hires%grid%i1, hires%grid%i2
    DO j = 1, hires%grid%ny
      IF (MODULO( basin_ID_dp_hires( j,i), 1._dp) > 0.01_dp) THEN
        basin_ID_dp_hires( j,i) = -1._dp
      END IF
    END DO
    END DO
    CALL sync
    
    ! For those, use extrapolation instead
    basin_ID_dp_hires_ext( :,hires%grid%i1:hires%grid%i2) = basin_ID_dp_hires( :,hires%grid%i1:hires%grid%i2)
    CALL sync
    
    DO i = hires%grid%i1, hires%grid%i2
    DO j = 1, hires%grid%ny
      IF (basin_ID_dp_hires_ext( j,i) == -1._dp) THEN
      
        n = 0
        foundit = .FALSE.
        DO WHILE (.NOT. foundit)
          
          n = n+1
          
          ! Take the value of the nearest non-boundary cell
          DO ii = MAX(1,i-n), MIN(hires%grid%nx,i+n)
          DO jj = MAX(1,j-n), MIN(hires%grid%ny,j+n)
            IF (basin_ID_dp_hires( jj,ii) > -1._dp) THEN
              basin_ID_dp_hires_ext( j,i) = basin_ID_dp_hires( jj,ii)
              foundit = .TRUE.
              EXIT
            END IF
          END DO
          IF (foundit) EXIT
          END DO
          
          ! Safety
          IF (n > MAX(hires%grid%nx, hires%grid%ny)) THEN
            WRITE(0,*) 'map_and_extrapolate_ocean_data - ERROR: basin ID downscaling got stuck!'
            CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
          END IF
          
        END DO ! DO WHILE (.NOT. foundit)
        
      END IF ! IF (basin_ID_dp_hires_ext( j,i) == -1._dp) THEN
    END DO
    END DO
    CALL sync
    
    ! Convert hi-resolution basin ID field back to integer precision
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
      hires%basin_ID( j,i) = NINT( basin_ID_dp_hires_ext( j,i))
    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wbasin_ID_dp_lores)
    CALL deallocate_shared( wbasin_ID_dp_hires)
    CALL deallocate_shared( wbasin_ID_dp_hires_ext)
    
    ! Perform the extrapolation on the high-resolution grid
    IF (par%master) WRITE(0,'(A,F4.1,A)') '     Performing ocean data extrapolation on the ', hires%grid%dx / 1000._dp, ' km regional x/y-grid...'
    CALL extend_regional_ocean_data_to_cover_domain( hires)
    
    ! Clean up fields that were needed only for the extrapolation
    CALL deallocate_shared( hires%wHi      )
    CALL deallocate_shared( hires%wHb      )
    CALL deallocate_shared( hires%wbasin_ID)
    CALL deallocate_shared( hires%wnbasins )
  
  END SUBROUTINE map_and_extrapolate_hires_ocean_data
  SUBROUTINE extend_regional_ocean_data_to_cover_domain( hires)
    ! Extend global ocean data over the whole grid, based on the procedure outlined in 
    ! Jourdain, N. C., Asay-Davis, X., Hattermann, T., Straneo, F., Seroussi, H., Little, C. M., & Nowicki, S. (2020). 
    ! A protocol for calculating basal melt rates in the ISMIP6 Antarctic ice sheet projections. The Cryosphere, 14(9), 3111-3134. 

    ! In/output variables:
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires
    
    ! Local variables:
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: NaN, Hs, z_bedrock, z_icebase, z
    INTEGER,  DIMENSION(:,:,:), POINTER                ::  mask_wetdry,  mask_hasdata
    INTEGER                                            :: wmask_wetdry, wmask_hasdata
    INTEGER                                            :: k1,k2,bi
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE            :: mask, mask_filled
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: d_T, d_S
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  T_ocean_ext,  S_ocean_ext
    INTEGER                                            :: wT_ocean_ext, wS_ocean_ext
    
    LOGICAL,  PARAMETER                                :: verbose = .FALSE.
    
    ! Useful
    NaN = -1._dp
    NaN = SQRT( NaN)
    
    ! Allocate shared memory
    CALL allocate_shared_dp_3D( hires%nz_ocean, hires%grid%ny, hires%grid%nx, T_ocean_ext, wT_ocean_ext)
    CALL allocate_shared_dp_3D( hires%nz_ocean, hires%grid%ny, hires%grid%nx, S_ocean_ext, wS_ocean_ext)
    
    ! Initialise the extrapolated product with the provided ocean data
    DO i = hires%grid%i1, hires%grid%i2
    DO j = 1, hires%grid%ny
    DO k = 1, hires%nz_ocean
      T_ocean_ext( k,j,i) = hires%T_ocean( k,j,i)
      S_ocean_ext( k,j,i) = hires%S_ocean( k,j,i)
    END DO
    END DO
    END DO
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
    
    CALL allocate_shared_int_3D( hires%nz_ocean, hires%grid%ny, hires%grid%nx, mask_wetdry,  wmask_wetdry )
    CALL allocate_shared_int_3D( hires%nz_ocean, hires%grid%ny, hires%grid%nx, mask_hasdata, wmask_hasdata)
    
    DO i = hires%grid%i1, hires%grid%i2
    DO j = 1, hires%grid%ny
      
      Hs = surface_elevation( hires%Hi( j,i), hires%Hb( j,i), 0._dp)
      z_bedrock = hires%Hb( j,i)
      z_icebase = Hs - hires%Hi( j,i)
      
      DO k = 1, hires%nz_ocean
        
        z = -hires%z_ocean( k)
        
        ! mask_wetdry
        IF (z < z_bedrock) THEN
          ! This 3D hires%grid box is beneath the bedrock surface, so dry
          mask_wetdry( k,j,i) = 3
        ELSE
          ! This 3D hires%grid box is above the bedrock surface so at least potentially wet
          IF (z < z_icebase) THEN
            ! This 3D hires%grid box is above the bedrock surface and below the ice base, so it is actually wet
            mask_wetdry( k,j,i) = 1
          ELSE
            ! This 3D hires%grid box is above the bedrock surface and above the ice base, so it is potentially wet (i.e. inside grounded marine ice)
            mask_wetdry( k,j,i) = 2
          END IF
        END IF
        
        ! mask_hasdata
        IF (hires%T_ocean( k,j,i) /= hires%T_ocean( k,j,i)) THEN
          ! This 3D hires%grid box has no data (yet)
          mask_hasdata( k,j,i) = 0
        ELSE
          ! Data is already provided for this 3D hires%grid box
          mask_hasdata( k,j,i) = 1
        END IF
        
      END DO
      
    END DO
    END DO
    CALL sync
    
  ! ================================================================
  ! ===== Step 1: horizontal extrapolation into shelf cavities =====
  ! ================================================================
    
    IF (par%master .AND. verbose) WRITE(0,*) '    extend_regional_ocean_data_to_cover_domain - step 1'
  
    ! Here, we start with the ocean data as provided (i.e. only for open ocean), and
    ! perform a horizontal extrapolation (so for each vertical layer separately) into
    ! the shelf cavities. Only "actually wet" 3D grid boxes are allowed to be filled,
    ! so the fill is limited by both bedrock sills and grounded ice.
    
    ! Allocate memory for the mask and data field of a single extrapolation step
    ALLOCATE( mask(        hires%grid%ny, hires%grid%nx))
    ALLOCATE( mask_filled( hires%grid%ny, hires%grid%nx))
    ALLOCATE( d_T(         hires%grid%ny, hires%grid%nx))
    ALLOCATE( d_S(         hires%grid%ny, hires%grid%nx))
    
    ! Parallelised by partitioning the vertical domain
    CALL partition_list( hires%nz_ocean, par%i, par%n, k1, k2)
    
    DO k = k1, k2
    
      ! Extrapolate per basin
      DO bi = 1, hires%nbasins
      
        IF (verbose) WRITE(0,'(A,I2,A,I3,A,I3,A,I3,A,I3)') '        process ', par%i, ': vertical layer ', k, '/', hires%nz_ocean, ', basin ', bi, '/', hires%nbasins
        
        ! Define the mask and initial data fields for this particular flood-fill
        ! (i.e. this vertical layer and this basin)
        mask        = 0
        d_T         = NaN
        d_S         = NaN
        DO i = 1, hires%grid%nx
        DO j = 1, hires%grid%ny
          IF (hires%basin_ID( j,i) == bi) THEN
            IF (mask_hasdata( k,j,i) == 1) THEN
              ! This is where the source data comes from
              mask( j,i) = 2
              d_T(  j,i) = T_ocean_ext( k,j,i)
              d_S(  j,i) = S_ocean_ext( k,j,i)
            ELSEIF (mask_hasdata( k,j,i) == 0 .AND. mask_wetdry( k,j,i) == 1) THEN
              ! This is where we're supposed to fill it in
              mask( j,i) = 1
            END IF
          END IF
        END DO
        END DO
        
        ! Perform the flood-fill-based Gaussian extrapolation
        CALL extrapolate_Gaussian_floodfill( hires%grid, mask, d_T, C%ocean_extrap_Gauss_sigma, mask_filled)
        CALL extrapolate_Gaussian_floodfill( hires%grid, mask, d_S, C%ocean_extrap_Gauss_sigma, mask_filled)
        
        ! Copy extrapolated data to the data structure
        DO i = 1, hires%grid%nx
        DO j = 1, hires%grid%ny
          IF (mask_filled( j,i) == 1) THEN
            T_ocean_ext(  k,j,i) = d_T( j,i)
            S_ocean_ext(  k,j,i) = d_S( j,i)
            mask_hasdata( k,j,i) = 2
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
    
  ! ===========================================================================
  ! ===== Step 2: vertical extrapolation into sill-blocked shelf cavities =====
  ! ===========================================================================
    
    IF (par%master .AND. verbose) WRITE(0,*) '    extend_regional_ocean_data_to_cover_domain - step 2'
  
    ! Here, we start with the ocean data that has been horizontally extrapolated into
    ! the shelf cavities, allowing for bedrock topography to block the fill, so that
    ! for example the lower parts of the Filchner-Ronne and Ross cavities have not yet
    ! been filled. We now extrapolate the data vertically from the filled parts to
    ! fill those parts of the cavities. Barring any really weird geometry, the entire
    ! cavities will now be filled.
    
    DO i = hires%grid%i1, hires%grid%i2
    DO j = 1, hires%grid%ny
      
      ! Move down through the vertical column
      DO k = 2, hires%nz_ocean
        ! If this grid box is wet and has no data, but the one above it does,
        ! copy data from the one above it.
        IF (mask_wetdry( k,j,i) == 1 .AND. mask_hasdata( k,j,i) == 0) THEN
          ! This 3D grid box is wet but has no data
          IF (mask_hasdata( k-1,j,i) == 1 .OR. mask_hasdata( k-1,j,i) == 2) THEN
            ! The one above it has data; copy data
            mask_hasdata( k,j,i) = 2
            T_ocean_ext( k,j,i) = T_ocean_ext( k-1,j,i)
            S_ocean_ext( k,j,i) = S_ocean_ext( k-1,j,i)
          END IF
        END IF
      END DO
      
    END DO
    END DO
    CALL sync
    
  ! ===============================================================
  ! ===== Step 3: vertical extrapolation into ice and bedrock =====
  ! ===============================================================
    
    IF (par%master .AND. verbose) WRITE(0,*) '    extend_regional_ocean_data_to_cover_domain - step 3'
  
    ! Extrapolate data vertically into 3D grid boxes that are occupied by ice
    ! or bedrock (since they might turn into ocean at some point during a simulation)
    
    DO i = hires%grid%i1, hires%grid%i2
    DO j = 1, hires%grid%ny
      
      ! Move down through the vertical column
      DO k = 2, hires%nz_ocean
        ! If this hires%grid box is wet and has no data, but the one above it does,
        ! copy data from the one above it.
        IF (mask_hasdata( k,j,i) == 0) THEN
          ! This 3D hires%grid box is wet but has no data
          IF (mask_hasdata( k-1,j,i) == 1 .OR. mask_hasdata( k-1,j,i) == 2) THEN
            ! The one above it has data; copy data
            mask_hasdata( k,j,i) = 2
            T_ocean_ext( k,j,i) = T_ocean_ext( k-1,j,i)
            S_ocean_ext( k,j,i) = S_ocean_ext( k-1,j,i)
          END IF
        END IF
      END DO
      
      ! Move up through the vertical column
      DO k = hires%nz_ocean-1, 1, -1
        ! If this hires%grid box is wet and has no data, but the one above it does,
        ! copy data from the one above it.
        IF (mask_hasdata( k,j,i) == 0) THEN
          ! This 3D hires%grid box is wet but has no data
          IF (mask_hasdata( k+1,j,i) == 1 .OR. mask_hasdata( k+1,j,i) == 2) THEN
            ! The one above it has data; copy data
            mask_hasdata( k,j,i) = 2
            T_ocean_ext(  k,j,i) = T_ocean_ext( k+1,j,i)
            S_ocean_ext(  k,j,i) = S_ocean_ext( k+1,j,i)
          END IF
        END IF
      END DO
      
    END DO
    END DO
    CALL sync
    
  ! =================================================================
  ! ===== Step 4: horizontal extrapolation into ice and bedrock =====
  ! =================================================================
    
    IF (par%master .AND. verbose) WRITE(0,*) '    extend_regional_ocean_data_to_cover_domain - step 4'
  
    ! In the last step, extrapolate data horizontally into 3D
    ! grid boxes that are occupied by ice or bedrock
    
    ! Allocate memory for the mask and data field of a single extrapolation step
    ALLOCATE( mask(        hires%grid%ny, hires%grid%nx))
    ALLOCATE( mask_filled( hires%grid%ny, hires%grid%nx))
    ALLOCATE( d_T(         hires%grid%ny, hires%grid%nx))
    ALLOCATE( d_S(         hires%grid%ny, hires%grid%nx))
    
    ! Parallelised by partitioning the vertical domain
    CALL partition_list( hires%nz_ocean, par%i, par%n, k1, k2)
    
    DO k = k1, k2
    
      ! Extrapolate per basin
      DO bi = 1, hires%nbasins
      
        IF (verbose) WRITE(0,'(A,I2,A,I3,A,I3,A,I3,A,I3)') '        process ', par%i, ': vertical layer ', k, '/', hires%nz_ocean, ', basin ', bi, '/', hires%nbasins
        
        ! Define the mask and initial data fields for this particular flood-fill
        ! (i.e. this vertical layer and this basin)
        mask        = 0
        d_T         = NaN
        d_S         = NaN
        DO i = 1, hires%grid%nx
        DO j = 1, hires%grid%ny
          IF (hires%basin_ID( j,i) == bi) THEN
            IF (mask_hasdata( k,j,i) == 1 .OR. mask_hasdata( k,j,i) == 2) THEN
              ! This is where the source data comes from
              mask( j,i) = 2
              d_T(  j,i) = T_ocean_ext( k,j,i)
              d_S(  j,i) = S_ocean_ext( k,j,i)
            ELSEIF (mask_hasdata( k,j,i) == 0) THEN
              ! This is where we're supposed to fill it in
              mask( j,i) = 1
            END IF
          END IF
        END DO
        END DO
        
        ! Perform the flood-fill-based Gaussian extrapolation
        CALL extrapolate_Gaussian_floodfill( hires%grid, mask, d_T, C%ocean_extrap_Gauss_sigma, mask_filled)
        CALL extrapolate_Gaussian_floodfill( hires%grid, mask, d_S, C%ocean_extrap_Gauss_sigma, mask_filled)
        
        ! Copy extrapolated data to the data structure
        DO i = 1, hires%grid%nx
        DO j = 1, hires%grid%ny
          IF (mask_filled( j,i) == 1) THEN
            T_ocean_ext(  k,j,i) = d_T( j,i)
            S_ocean_ext(  k,j,i) = d_S( j,i)
            mask_hasdata( k,j,i) = 2
          END IF
        END DO
        END DO
        
      END DO ! DO bi = 1, ice%nbasins
      
      ! One more pass without considering basins (sometimes, a handful of isolated "basin enclaves" can
      ! occur at high resolution, which will not be filled when using the basin-constrained flood-fill)
        
      ! Define the mask and initial data fields for this particular flood-fill
      ! (i.e. this vertical layer and this basin)
      mask        = 0
      d_T         = NaN
      d_S         = NaN
      DO i = 1, hires%grid%nx
      DO j = 1, hires%grid%ny
        IF (mask_hasdata( k,j,i) == 1 .OR. mask_hasdata( k,j,i) == 2) THEN
          ! This is where the source data comes from
          mask( j,i) = 2
          d_T(  j,i) = T_ocean_ext( k,j,i)
          d_S(  j,i) = S_ocean_ext( k,j,i)
        ELSEIF (mask_hasdata( k,j,i) == 0) THEN
          ! This is where we're supposed to fill it in
          mask( j,i) = 1
        END IF
      END DO
      END DO
      
      ! Perform the flood-fill-based Gaussian extrapolation
      CALL extrapolate_Gaussian_floodfill( hires%grid, mask, d_T, C%ocean_extrap_Gauss_sigma, mask_filled)
      CALL extrapolate_Gaussian_floodfill( hires%grid, mask, d_S, C%ocean_extrap_Gauss_sigma, mask_filled)
      
      ! Copy extrapolated data to the data structure
      DO i = 1, hires%grid%nx
      DO j = 1, hires%grid%ny
        IF (mask_filled( j,i) == 1) THEN
          T_ocean_ext(  k,j,i) = d_T( j,i)
          S_ocean_ext(  k,j,i) = d_S( j,i)
          mask_hasdata( k,j,i) = 2
        END IF
      END DO
      END DO
      
      ! Check if all pixels have now been filled
      DO i = 1, hires%grid%nx
      DO j = 1, hires%grid%ny
        IF (T_ocean_ext( k,j,i) /= T_ocean_ext( k,j,i)) THEN
          WRITE(0,*) '    extend_regional_ocean_data_to_cover_domain - ERROR: unfilled pixels remains at [i,j] = [', i,',',j,']'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END DO
      END DO
      
    END DO ! DO k = k1, k2
    CALL sync
    
    ! Copy data to the hires structure
    DO i = hires%grid%i1, hires%grid%i2
    DO j = 1, hires%grid%ny
    DO k = 1, hires%nz_ocean
      hires%T_ocean( k,j,i) = T_ocean_ext( k,j,i)
      hires%S_ocean( k,j,i) = S_ocean_ext( k,j,i)
    END DO
    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    DEALLOCATE( mask       )
    DEALLOCATE( mask_filled)
    DEALLOCATE( d_T        )
    DEALLOCATE( d_S        )
    CALL deallocate_shared( wmask_wetdry )
    CALL deallocate_shared( wmask_hasdata)
    CALL deallocate_shared( wT_ocean_ext )
    CALL deallocate_shared( wS_ocean_ext )
    
  END SUBROUTINE extend_regional_ocean_data_to_cover_domain
  SUBROUTINE write_hires_extrapolated_ocean_data_to_file( hires, filename_ocean_glob, hires_ocean_foldername)
    ! 1. Create a new folder inside the "extrapolated_ocean_files" folder
    ! 2. Create a header file inside this new folder listing the current model settings
    ! 3. Create a NetCDF file inside this new folder
    ! 4. Write the high-resolution extrapolated ocean data to this NetCDF file
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires
    CHARACTER(LEN=256),                  INTENT(IN)    :: filename_ocean_glob
    CHARACTER(LEN=256),                  INTENT(IN)    :: hires_ocean_foldername
    
    ! Local variables:
    CHARACTER(LEN=256)                                 :: hires_ocean_filename
    
    ! Create a new folder where the new extrapolated ocean file will be stored.
    IF (par%master) CALL system('mkdir ' // TRIM(hires_ocean_foldername))
    
    ! Create a header file describing the current ice-model set-up.
    CALL write_ocean_header( hires, filename_ocean_glob, hires_ocean_foldername)
    
    ! Create a NetCDF file and write data to it
    IF (par%master) hires_ocean_filename = TRIM(hires_ocean_foldername)//'/extrapolated_ocean_data.nc'
    IF (par%master) WRITE(0,*) '    Writing extrapolated ocean data to file "', TRIM(hires_ocean_filename), '"...'
    CALL create_extrapolated_ocean_file( hires, hires_ocean_filename)
    
  END SUBROUTINE write_hires_extrapolated_ocean_data_to_file
  SUBROUTINE check_for_matching_ocean_header( region, ocean_reg, filename_ocean_glob, foundmatch, hires_ocean_foldername)
    ! Inspect all the folder inside the "extrapolated_ocean_files" folder, read their header files,
    ! and see if any of them match the current model settings. If so, return the name of the folder
    ! where the matching header was found. If not, return the name of the folder where a new header
    ! should be created.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_model_region),             INTENT(IN)    :: region
    TYPE(type_subclimate_region),        INTENT(IN)    :: ocean_reg
    CHARACTER(LEN=256),                  INTENT(IN)    :: filename_ocean_glob
    LOGICAL,                             INTENT(OUT)   :: foundmatch
    CHARACTER(LEN=256),                  INTENT(OUT)   :: hires_ocean_foldername
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    CHARACTER(LEN=256)                                 :: header_filename
    LOGICAL                                            :: header_exists
    INTEGER                                            :: folder_i
    
    CHARACTER(LEN=256)                                 :: original_ocean_filename_read
    CHARACTER(LEN=256)                                 :: choice_ocean_vertical_grid_read
    INTEGER                                            :: nz_ocean_read
    REAL(dp)                                           :: ocean_vertical_grid_max_depth_read
    REAL(dp)                                           :: ocean_extrap_res_read
    REAL(dp)                                           :: ocean_extrap_Gauss_sigma_read
    REAL(dp)                                           :: lambda_M_read
    REAL(dp)                                           :: phi_M_read
    REAL(dp)                                           :: alpha_stereo_read
    
    IF (par%master) THEN
    
      folder_i = 1
      DO WHILE (folder_i < 1000)
      
        ! Generate a foldername to inspect
        IF     (folder_i < 10)   THEN
          WRITE( hires_ocean_foldername,'(A,A,A,A,I1)') TRIM(C%ocean_extrap_dir), '/', region%name, '_00', folder_i
        ELSEIF (folder_i < 100)  THEN
          WRITE( hires_ocean_foldername,'(A,A,A,A,I2)') TRIM(C%ocean_extrap_dir), '/', region%name, '_0',  folder_i
        ELSEIF (folder_i < 1000) THEN
          WRITE( hires_ocean_foldername,'(A,A,A,A,I3)') TRIM(C%ocean_extrap_dir), '/', region%name, '_',   folder_i
        ELSE
          IF (par%master) WRITE(0,*) ' ERROR: tried a thousand folders!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
        ! Check if a header in this folder exists. If not, then we've inspected all existing headers
        ! without finding the good one, so we must generate the extrapolated ocean files from scratch.
        header_filename = TRIM( hires_ocean_foldername)//'/'//'header.txt'
        INQUIRE( FILE = header_filename, EXIST = header_exists)
        
        IF (.NOT. header_exists) THEN
          ! No more headers exist to be inspected. Return foundmatch = .FALSE.
          
          foundmatch = .FALSE.
          EXIT
          
        ELSE
          ! If the header exists, read it and see if it fits the current ice-model set-up.
          
          CALL read_ocean_header( &
            header_filename,                    &
            original_ocean_filename_read,       &
            choice_ocean_vertical_grid_read,    &
            nz_ocean_read,                      &
            ocean_vertical_grid_max_depth_read, &
            ocean_extrap_res_read,              &
            ocean_extrap_Gauss_sigma_read,      &
            lambda_M_read,                      &
            phi_M_read,                         &
            alpha_stereo_read)
          
          IF ( TRIM(original_ocean_filename_read)    == TRIM(filename_ocean_glob)             .AND. &
               TRIM(choice_ocean_vertical_grid_read) == TRIM(choice_ocean_vertical_grid_read) .AND. &
               nz_ocean_read                         == ocean_reg%nz_ocean                    .AND. &
               ocean_vertical_grid_max_depth_read    == C%ocean_vertical_grid_max_depth       .AND. &
               ocean_extrap_res_read                 == C%ocean_extrap_res                    .AND. &
               ocean_extrap_Gauss_sigma_read         == C%ocean_extrap_Gauss_sigma            .AND. &
               lambda_M_read                         == region%grid%lambda_M                  .AND. &
               phi_M_read                            == region%grid%phi_M                     .AND. &
               alpha_stereo_read                     == region%grid%alpha_stereo) THEN
            ! This header matches the current model set-up!
            
            foundmatch = .TRUE.
            EXIT 
                     
          ELSE
            ! This header doesn't match the current model set-up. Try the next one.
            folder_i = folder_i + 1      
          END IF
          
        END IF
      
      END DO ! DO WHILE (header_i < 1000)
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    CALL MPI_BCAST( hires_ocean_foldername,                 256, MPI_CHAR,    0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( foundmatch,                 1,   MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( header_filename,            256, MPI_CHAR,    0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( folder_i,                   1,   MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
  END SUBROUTINE check_for_matching_ocean_header
  SUBROUTINE read_ocean_header( &
            header_filename,                    &
            original_ocean_filename,       &
            choice_ocean_vertical_grid,    &
            nz_ocean,                      &
            ocean_vertical_grid_max_depth, &
            ocean_extrap_res,              &
            ocean_extrap_Gauss_sigma,      &
            lambda_M,                      &
            phi_M,                         &
            alpha_stereo                   )
    ! Read a header file listing the model settings that were used to create a high-resolution extrapolated ocean data file
    
    IMPLICIT NONE
    
    ! In/output variables:
    CHARACTER(LEN=256),                  INTENT(IN)    :: header_filename
    CHARACTER(LEN=256),                  INTENT(OUT)   :: original_ocean_filename
    CHARACTER(LEN=256),                  INTENT(OUT)   :: choice_ocean_vertical_grid
    INTEGER,                             INTENT(OUT)   :: nz_ocean
    REAL(dp),                            INTENT(OUT)   :: ocean_vertical_grid_max_depth
    REAL(dp),                            INTENT(OUT)   :: ocean_extrap_res
    REAL(dp),                            INTENT(OUT)   :: ocean_extrap_Gauss_sigma
    REAL(dp),                            INTENT(OUT)   :: lambda_M
    REAL(dp),                            INTENT(OUT)   :: phi_M
    REAL(dp),                            INTENT(OUT)   :: alpha_stereo
    
    ! Local variables:
    INTEGER                                            :: ios, cerr, ierr
    
    ! The NAMELIST that's used to read the external header file.
    NAMELIST /HEADER/original_ocean_filename,             &
                     choice_ocean_vertical_grid,          &
                     nz_ocean,                            &
                     ocean_vertical_grid_max_depth,       &
                     ocean_extrap_res,                    &
                     ocean_extrap_Gauss_sigma,            &
                     lambda_M,                            &
                     phi_M,                               &
                     alpha_stereo
      
    OPEN( UNIT = 29, FILE = TRIM(header_filename), STATUS='OLD', ACTION='READ', iostat=ios)
    IF (ios /= 0) THEN
      WRITE(0,*) ' ERROR: could not open ""', TRIM(header_filename), '"'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! In the following statement the entire configuration file is read, using the namelist (NML=HEADER)
    READ(  UNIT = 29, NML = HEADER, IOSTAT = ios)
    CLOSE( UNIT = 29)

    IF (ios /= 0) THEN
      WRITE(0,*) ' ERROR: could not read ""', TRIM(header_filename), '"'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE read_ocean_header
  SUBROUTINE write_ocean_header( hires, filename_ocean_glob, hires_ocean_foldername)
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires
    CHARACTER(LEN=256),                  INTENT(IN)    :: filename_ocean_glob
    CHARACTER(LEN=256),                  INTENT(IN)    :: hires_ocean_foldername
    
    ! Local variables:
    INTEGER, DIMENSION(8)                              :: datevec
    CHARACTER(LEN=256)                                 :: header_filename
    
    ! Let the Master do the work
    IF (par%master) THEN
    
      header_filename = TRIM( hires_ocean_foldername)//'/header.txt'
      
      OPEN( UNIT = 1337, FILE = header_filename, STATUS = 'NEW')
      
      CALL date_and_time( VALUES = datevec)
      
      WRITE(UNIT = 1337, FMT = '(A)') '&HEADER'
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A,I4,A,I2,A,I2)') '! Icemodel-ocean header file, created on ', datevec(1), '-', datevec(2), '-', datevec(3)
      WRITE(UNIT = 1337, FMT = '(A)') '!'
      WRITE(UNIT = 1337, FMT = '(A)') '! This header describes the icemodel set-up that was used to created this'
      WRITE(UNIT = 1337, FMT = '(A)') '! extrapolated ocean data file. Since creating these is computationally intensive,'
      WRITE(UNIT = 1337, FMT = '(A)') '! reading them from files is preferred. These header files make ice model'
      WRITE(UNIT = 1337, FMT = '(A)') '! a bit more flexible when using different input files for the four model regions.'
      WRITE(UNIT = 1337, FMT = '(A)') ''
      
      WRITE(UNIT = 1337, FMT = '(A)')       '! The original global ocean file that was extrapolated'
      WRITE(UNIT = 1337, FMT = '(A,A,A)')   'original_ocean_filename       = ''', TRIM(filename_ocean_glob), ''''
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A)')       '! The vertical grid the ocean data was projected to'
      WRITE(UNIT = 1337, FMT = '(A,A,A)')   'choice_ocean_vertical_grid    = ''', TRIM(C%choice_ocean_vertical_grid), ''''
      WRITE(UNIT = 1337, FMT = '(A,I5)')    'nz_ocean                      = ', hires%nz_ocean
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'ocean_vertical_grid_max_depth = ', C%ocean_vertical_grid_max_depth
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A)')       '! Resolution and Gaussian smoothing radius used for the high-resolution extrapolation'
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'ocean_extrap_res              = ', C%ocean_extrap_res
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'ocean_extrap_Gauss_sigma      = ', C%ocean_extrap_Gauss_sigma
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A)')       '! Parameters of the high-resolution grid'
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'lambda_M                      = ', hires%grid%lambda_M
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'phi_M                         = ', hires%grid%phi_M
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'alpha_stereo                  = ', hires%grid%alpha_stereo
      
      WRITE(UNIT = 1337, FMT = '(A)') ''      
      WRITE(UNIT = 1337, FMT = '(A)') '/'
      
      CLOSE(UNIT = 1337)
    
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE write_ocean_header
  
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
