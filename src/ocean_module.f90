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
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_climate_matrix, type_subclimate_global, &
                                             type_climate_model, type_subclimate_region
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, &
                                             inquire_PD_obs_data_file_ocean, read_PD_obs_data_file_ocean
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             error_function, smooth_Gaussian_2D, smooth_Gaussian_3D, smooth_Shepard_2D, &
                                             extend_Gaussian_2D, extend_Gaussian_3D, check_for_NaN_dp_2D_return, &
                                             map_glob_to_grid_2D, map_glob_to_grid_3D, &
                                             map_square_to_square_cons_2nd_order_2D, map_square_to_square_cons_2nd_order_3D

  IMPLICIT NONE
    
CONTAINS
  
! == Initialising the region-specific ocean data
  SUBROUTINE initialise_oceans_regional( grid, ice, climate, matrix)
    ! Allocate shared memory for the ocean data of the regional subclimates
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice  
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
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
    
    ! Exception for schematic ocean temperature/salinity profile
    ! Set oceans to the ISOMIP+ "COLD" or "WARM" profile
    IF (C%choice_ocean_temperature_model == 'schematic') THEN
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
    CALL allocate_subclimate_regional_oceans( grid, climate%PD_obs,   matrix%PD_obs_ocean%z_ocean, matrix%PD_obs_ocean%nz_ocean)
    CALL allocate_subclimate_regional_oceans( grid, climate%applied,  matrix%PD_obs_ocean%z_ocean, matrix%PD_obs_ocean%nz_ocean)
    
    CALL map_ocean_data_global_to_regional( grid, matrix%PD_obs_ocean, climate%PD_obs  )    
    CALL extend_ocean_data_to_cover_grid ( grid, climate%PD_obs, ice )
    
    IF (C%choice_ocean_temperature_model == 'WOA') THEN
      IF (par%master) WRITE(*,*) '    Constant present-day ocean forcing used for BMB forcing!'
    ELSE IF (C%choice_ocean_temperature_model == 'PI' .OR. C%choice_ocean_temperature_model == 'scaled')  THEN
      IF (par%master) WRITE(*,*) '    Ocean forcing initialised, but not actually used for BMB forcing!'
    ELSE IF (C%choice_ocean_temperature_model == 'matrix_warm_cold') THEN      
      CALL allocate_subclimate_regional_oceans( grid, climate%GCM_PI,   matrix%PD_obs_ocean%z_ocean, matrix%PD_obs_ocean%nz_ocean)
      CALL allocate_subclimate_regional_oceans( grid, climate%GCM_cold, matrix%PD_obs_ocean%z_ocean, matrix%PD_obs_ocean%nz_ocean)
      CALL allocate_subclimate_regional_oceans( grid, climate%GCM_warm, matrix%PD_obs_ocean%z_ocean, matrix%PD_obs_ocean%nz_ocean)
    
      ! Map ocean data from the global subclimates to the regional subclimates
      ! TODO: z-coordinate interpolation    
      CALL map_ocean_data_global_to_regional( grid, matrix%GCM_PI,       climate%GCM_PI  )
      CALL map_ocean_data_global_to_regional( grid, matrix%GCM_cold,     climate%GCM_cold)
      CALL map_ocean_data_global_to_regional( grid, matrix%GCM_warm,     climate%GCM_warm)
    
      ! Correct regional ocean data for GCM bias
      CALL correct_GCM_bias_ocean( grid, climate, climate%GCM_PI)
      CALL correct_GCM_bias_ocean( grid, climate, climate%GCM_warm)
      CALL correct_GCM_bias_ocean( grid, climate, climate%GCM_cold)
      
    ELSE  
      IF (par%master) WRITE(0,*) '  ERROR: choice_ocean_temperature_model "', TRIM(C%choice_ocean_temperature_model), '" not implemented in initialise_oceans_regional!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF !(C%choice_ocean_temperature_model == 'WOA')
  
    ! Initialise applied ocean forcing with present-day observations
    IF (par%master) THEN
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
  SUBROUTINE allocate_subclimate_regional_oceans( grid, climate, z_ocean, nz_ocean)
    ! Allocate shared memory for the ocean data of a regional subclimate models
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: z_ocean
    INTEGER,                             INTENT(IN)    :: nz_ocean
    
    ! Set up the vertical coordinate
    CALL allocate_shared_int_0D( climate%nz_ocean, climate%wnz_ocean)
    IF (par%master) climate%nz_ocean = nz_ocean
    CALL sync
    CALL allocate_shared_dp_1D( climate%nz_ocean, climate%z_ocean, climate%wz_ocean)
    IF (par%master) climate%z_ocean = z_ocean
    CALL sync
    
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
    TYPE(type_subclimate_global),        INTENT(IN)    :: clim_glob
    TYPE(type_subclimate_region),        INTENT(INOUT) :: clim_reg
    
    ! TODO: instead of a taking the WOA vertical, remap to a predefined regular (e.g. 60-km) vertical coordinates 
    CALL map_glob_to_grid_3D ( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%T_ocean, clim_reg%T_ocean, clim_glob%nz_ocean)
    CALL map_glob_to_grid_3D ( clim_glob%nlat, clim_glob%nlon, clim_glob%lat, clim_glob%lon, grid, clim_glob%S_ocean, clim_reg%S_ocean, clim_glob%nz_ocean)    
    
    ! TODO: check if needed
    clim_reg%mask_ocean       = 0._dp   
  
  END SUBROUTINE map_ocean_data_global_to_regional
  SUBROUTINE extend_ocean_data_to_cover_grid (grid, clim_reg, ice)
    ! Extend the ocean data over the whole grid, based on the procedure outlined in 
    ! Jourdain, N. C., Asay-Davis, X., Hattermann, T., Straneo, F., Seroussi, H., Little, C. M., & Nowicki, S. (2020). 
    ! A protocol for calculating basal melt rates in the ISMIP6 Antarctic ice sheet projections. The Cryosphere, 14(9), 3111-3134. 

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_subclimate_region),        INTENT(INOUT) :: clim_reg
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    
    ! local variables
    INTEGER                                            :: l ! iterations
    LOGICAL                                            :: check = .TRUE.
    
    clim_reg%T_ocean_ext      = clim_reg%T_ocean
    clim_reg%S_ocean_ext      = clim_reg%S_ocean
    
    DO l=1,45
      WRITE(*,*) "l=",l
      CALL extend_Gaussian_3D ( grid, clim_reg%T_ocean_ext, 8000._dp, 12000._dp, 67, ice%basin_ID) ! TODO: 67 -> clim_reg%nz_ocean
      CALL extend_Gaussian_3D ( grid, clim_reg%S_ocean_ext, 8000._dp, 12000._dp, 67, ice%basin_ID) ! TODO: 67 -> clim_reg%nz_ocean
      
      CALL check_for_NaN_dp_2D_return ( clim_reg%T_ocean_ext(1,:,:) ,check ) 
      !WRITE(*,*) 'check=',check

    END DO
    
    !CALL extend_Gaussian_3D ( grid, clim_reg%T_ocean_ext, 1250000._dp, 2500000._dp, 67, ice%basin_ID)
    !CALL extend_Gaussian_3D ( grid, clim_reg%S_ocean_ext, 1250000._dp, 2500000._dp, 67, ice%basin_ID)
    
    ! Initialize the final bias-corrected value. Bias correction will follow later for the climates where this is desired.
    clim_reg%T_ocean_corr_ext      = clim_reg%T_ocean_ext
    clim_reg%S_ocean_corr_ext      = clim_reg%S_ocean_ext
    
  END SUBROUTINE extend_ocean_data_to_cover_grid
  SUBROUTINE correct_GCM_bias_ocean( grid, climate, subclimate)
    ! Correct regional ocean data for GCM bias
    ! (must be done on regional grid, since the GCM grid and the World Ocean Atlas grid are generally not the same!)
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_climate_model),            INTENT(IN)    :: climate
    TYPE(type_subclimate_region),        INTENT(INOUT) :: subclimate
    
    subclimate%T_ocean_corr_ext = subclimate%T_ocean_ext - (climate%GCM_PI%T_ocean_ext - climate%PD_obs%T_ocean_ext)
    subclimate%T_ocean_corr_ext = subclimate%S_ocean_ext - (climate%GCM_PI%S_ocean_ext - climate%PD_obs%S_ocean_ext)
       
  
  END SUBROUTINE correct_GCM_bias_ocean
  

  ! Allocate shared memory for the global PD observed ocean data fields (stored in the climate matrix),
  ! read them from the specified NetCDF file (latter only done by master process).  
  SUBROUTINE initialise_PD_obs_ocean_fields (PD_obs_ocean, name)
     
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_subclimate_global),   INTENT(INOUT) :: PD_obs_ocean
    CHARACTER(LEN=*),               INTENT(IN)    :: name
    
    PD_obs_ocean%name = name 
    PD_obs_ocean%netcdf%filename   = C%filename_PD_obs_ocean
        
    ! Inquire if all required variables are present in the NetCDF file, and read the grid size.
    CALL allocate_shared_int_0D( PD_obs_ocean%nlon,     PD_obs_ocean%wnlon)
    CALL allocate_shared_int_0D( PD_obs_ocean%nlat,     PD_obs_ocean%wnlat)
    CALL allocate_shared_int_0D( PD_obs_ocean%nz_ocean, PD_obs_ocean%wnz_ocean)
    IF (par%master) CALL inquire_PD_obs_data_file_ocean(PD_obs_ocean)
    CALL sync
    
    ! Allocate memory  
    CALL allocate_shared_dp_1D(                        PD_obs_ocean%nlon,                    PD_obs_ocean%lon,     PD_obs_ocean%wlon    )
    CALL allocate_shared_dp_1D(                                           PD_obs_ocean%nlat, PD_obs_ocean%lat,     PD_obs_ocean%wlat    )
    CALL allocate_shared_dp_1D( PD_obs_ocean%nz_ocean,                                       PD_obs_ocean%z_ocean, PD_obs_ocean%wz_ocean)        
    !CALL allocate_shared_dp_3D( PD_obs_ocean%nlon, PD_obs_ocean%nlat, PD_obs_ocean%nz_ocean, PD_obs_ocean%mask_ocean,      PD_obs_ocean%wPrecip     )
    CALL allocate_shared_dp_3D( PD_obs_ocean%nlon, PD_obs_ocean%nlat, PD_obs_ocean%nz_ocean, PD_obs_ocean%T_ocean,     PD_obs_ocean%wT_ocean    )
    CALL allocate_shared_dp_3D( PD_obs_ocean%nlon, PD_obs_ocean%nlat, PD_obs_ocean%nz_ocean, PD_obs_ocean%S_ocean,     PD_obs_ocean%wS_ocean    )
     
    ! Read data from the NetCDF file
    IF (par%master) WRITE(0,*) '   Reading PD observed ocean data from file ', TRIM(PD_obs_ocean%netcdf%filename), '...'
    IF (par%master) CALL read_PD_obs_data_file_ocean(PD_obs_ocean)
    CALL sync
      
    ! Determine process domains
    CALL partition_list( PD_obs_ocean%nlon, par%i, par%n, PD_obs_ocean%i1, PD_obs_ocean%i2)
      
  END SUBROUTINE initialise_PD_obs_ocean_fields
  
! == Some schematic ocean temperature/salinity profiles
  SUBROUTINE set_ocean_to_ISOMIPplus_COLD( grid, climate)
    ! Set the ocean temperature and salinity to the ISOMIP+ "COLD" profile (Asay-Davis et al., 2016, Table 5)
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    
    ! Local variables
    INTEGER                                            :: i,j,k
    INTEGER,  PARAMETER                                :: nz_ocean  = 30
    REAL(dp), DIMENSION(:    ), POINTER                ::  z_ocean
    INTEGER                                            :: wz_ocean
    REAL(dp)                                           :: w
    REAL(dp), PARAMETER                                :: Tzero     = -1.9_dp    ! Sea surface temperature [degC] (originally T0, but that name is already taken...)
    REAL(dp), PARAMETER                                :: Tbot      = -1.9_dp    ! Sea floor   temperature [degC]
    REAL(dp), PARAMETER                                :: Szero     = 33.8_dp    ! Sea surface salinity    [PSU]
    REAL(dp), PARAMETER                                :: Sbot      = 34.55_dp   ! Sea floor   salinity    [PSU]
    REAL(dp), PARAMETER                                :: depth_max = 720._dp    ! Maximum depth for the profile (constant values below that)
    
    ! Create vertical dimension
    CALL allocate_shared_dp_1D( nz_ocean, z_ocean, wz_ocean)
    IF (par%master) THEN
      DO k = 1, nz_ocean
        z_ocean( k) = depth_max * REAL(k-1,dp) / REAL(nz_ocean-1,dp)
      END DO
    END IF
    CALL sync
    
    ! Allocate memory for the ocean data
    CALL allocate_subclimate_regional_oceans( grid, climate, z_ocean, nz_ocean)
    
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
    
    ! Clean up after yourself
    CALL deallocate_shared( wz_ocean)
  
  END SUBROUTINE set_ocean_to_ISOMIPplus_COLD
  SUBROUTINE set_ocean_to_ISOMIPplus_WARM( grid, climate)
    ! Set the ocean temperature and salinity to the ISOMIP+ "WARM" profile (Asay-Davis et al., 2016, Table 6)
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid  
    TYPE(type_subclimate_region),        INTENT(INOUT) :: climate
    
    ! Local variables
    INTEGER                                            :: i,j,k
    INTEGER,  PARAMETER                                :: nz_ocean  = 30
    REAL(dp), DIMENSION(:    ), POINTER                ::  z_ocean
    INTEGER                                            :: wz_ocean
    REAL(dp)                                           :: w
    REAL(dp), PARAMETER                                :: Tzero     = -1.9_dp    ! Sea surface temperature [degC] (originally T0, but that name is already taken...)
    REAL(dp), PARAMETER                                :: Tbot      =  1.0_dp    ! Sea floor   temperature [degC]
    REAL(dp), PARAMETER                                :: Szero     = 33.8_dp    ! Sea surface salinity    [PSU]
    REAL(dp), PARAMETER                                :: Sbot      = 34.7_dp    ! Sea floor   salinity    [PSU]
    REAL(dp), PARAMETER                                :: depth_max = 720._dp    ! Maximum depth for the profile (constant values below that)
    
    ! Create vertical dimension
    CALL allocate_shared_dp_1D( nz_ocean, z_ocean, wz_ocean)
    IF (par%master) THEN
      DO k = 1, nz_ocean
        z_ocean( k) = depth_max * REAL(k-1,dp) / REAL(nz_ocean-1,dp)
      END DO
    END IF
    CALL sync
    
    ! Allocate memory for the ocean data
    CALL allocate_subclimate_regional_oceans( grid, climate, z_ocean, nz_ocean)
    
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
    
    ! Clean up after yourself
    CALL deallocate_shared( wz_ocean)
  
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
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: T,S
    INTEGER,  PARAMETER                                :: nz_ocean  = 2
    REAL(dp), DIMENSION(:    ), POINTER                ::  z_ocean
    INTEGER                                            :: wz_ocean
    REAL(dp), PARAMETER                                :: depth_max = 5000._dp    ! Maximum depth for the profile (constant values below that)
    
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
    
    ! Create vertical dimension
    CALL allocate_shared_dp_1D( nz_ocean, z_ocean, wz_ocean)
    IF (par%master) THEN
      DO k = 1, nz_ocean
        z_ocean( k) = depth_max * REAL(k-1,dp) / REAL(nz_ocean-1,dp)
      END DO
    END IF
    CALL sync
    
    ! Allocate memory for the ocean data
    CALL allocate_subclimate_regional_oceans( grid, climate, z_ocean, nz_ocean)
    
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
    
    ! Clean up after yourself
    CALL deallocate_shared( wz_ocean)
    
  END SUBROUTINE set_ocean_to_Reese2018

END MODULE ocean_module
