MODULE netcdf_output_module

! ===== Creating and writing to output files =====
! ================================================
!
! These routines create and write data to output files.

! ===== Preamble =====
! ====================

  ! Import basic functionality
! #include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, &
                                             allocate_shared_int_0D,   allocate_shared_dp_0D, &
                                             allocate_shared_int_1D,   allocate_shared_dp_1D, &
                                             allocate_shared_int_2D,   allocate_shared_dp_2D, &
                                             allocate_shared_int_3D,   allocate_shared_dp_3D, &
                                             allocate_shared_dp_4D,    deallocate_shared
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D

  ! Import specific functionality
  USE data_types_module,               ONLY: type_grid, type_grid_lonlat, type_model_region, type_ice_model, type_global_scalar_data, &
                                             type_highres_ocean_data, type_forcing_data, type_ocean_snapshot_regional
  USE netcdf,                          ONLY: NF90_NOERR, NF90_OPEN, NF90_CLOSE, NF90_NOWRITE, NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, &
                                             NF90_INQ_VARID, NF90_INQUIRE_VARIABLE, NF90_MAX_VAR_DIMS, NF90_GET_VAR, &
                                             NF90_CREATE, NF90_NOCLOBBER, NF90_NETCDF4, NF90_ENDDEF, NF90_REDEF, NF90_DEF_DIM, NF90_DEF_VAR, &
                                             NF90_PUT_ATT, NF90_WRITE, NF90_INT, NF90_INT64, NF90_FLOAT, NF90_DOUBLE, NF90_PUT_VAR, NF90_UNLIMITED, &
                                             NF90_INQUIRE_ATTRIBUTE
  USE utilities_module,                ONLY: deallocate_grid, deallocate_grid_lonlat, remap_zeta_grid_dp, &
                                             transpose_dp_2D, transpose_dp_3D, permute_2D_dp, permute_3D_dp,permute_3D_int
  USE netcdf_basic_module,             ONLY: nerr, field_name_options_x, field_name_options_y, field_name_options_zeta, field_name_options_z_ocean, &
                                             field_name_options_lon, field_name_options_lat, field_name_options_time, field_name_options_month, &
                                             field_name_options_Hi, field_name_options_Hb, field_name_options_Hs, field_name_options_dHb, &
                                             field_name_options_SL, field_name_options_Ti, &
                                             inquire_dim_multiple_options, inquire_var_multiple_options, &
                                             read_var_int_0D, read_var_int_1D, read_var_int_2D, read_var_int_3D, read_var_int_4D, &
                                             read_var_dp_0D , read_var_dp_1D , read_var_dp_2D , read_var_dp_3D , read_var_dp_4D, &
                                             check_x, check_y, check_lon, check_lat, check_zeta, check_z_ocean, find_timeframe, &
                                             create_new_netcdf_file_for_writing, create_dimension, create_variable, add_attribute_char, &
                                             write_var_int_0D, write_var_int_1D, write_var_int_2D, write_var_int_3D, write_var_int_4D, &
                                             write_var_dp_0D , write_var_dp_1D , write_var_dp_2D , write_var_dp_3D , write_var_dp_4D, &
                                             check_xy_grid_field_int_2D, check_xy_grid_field_dp_2D, check_xy_grid_field_dp_2D_monthly, check_xy_grid_field_dp_3D, &
                                             check_lonlat_grid_field_int_2D, check_lonlat_grid_field_dp_2D, check_lonlat_grid_field_dp_2D_monthly, check_lonlat_grid_field_dp_3D, &
                                             inquire_xy_grid, inquire_lonlat_grid, get_first_option_from_list, check_month, check_time, check_time_history, check_xy_grid_field_dp_3D_ocean

  IMPLICIT NONE

CONTAINS

  ! ===== Top-level functions =====
  ! ===============================

  SUBROUTINE write_to_regional_scalar_file( filename, region)
    ! Write model output to the regional scalar file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_regional_scalar_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write new time to file (thus extending the time dimension by one frame, making room for the new model data)
    CALL write_time_to_file( filename, region%time)

    ! Write model data
    ! ================
    CALL write_to_field_dp_0D( filename, 'ice_volume',    region%ice_volume)
    CALL write_to_field_dp_0D( filename, 'ice_volume_af', region%ice_volume_above_flotation)
    CALL write_to_field_dp_0D( filename, 'ice_area',      region%ice_area)
    CALL write_to_field_dp_0D( filename, 'T2m',           region%int_T2m)
    CALL write_to_field_dp_0D( filename, 'SMB',           region%int_SMB)
    CALL write_to_field_dp_0D( filename, 'BMB',           region%int_BMB)
    CALL write_to_field_dp_0D( filename, 'MB',            region%int_MB)

    ! Individual SMB components
    IF     (C%choice_SMB_model == 'uniform' .OR. &
            C%choice_SMB_model == 'idealised' .OR. &
            C%choice_SMB_model == 'direct_global' .OR. &
            C%choice_SMB_model == 'direct_regional' .OR. &
            C%choice_SMB_model == 'snapshot') THEN
      ! Do nothing
    ELSEIF (C%choice_SMB_model == 'IMAU-ITM' .OR. &
            C%choice_SMB_model == 'IMAU-ITM_wrongrefreezing') THEN
       CALL write_to_field_dp_0D( filename, 'snowfall',   region%int_snowfall)
       CALL write_to_field_dp_0D( filename, 'rainfall',   region%int_rainfall)
       CALL write_to_field_dp_0D( filename, 'melt',       region%int_melt)
       CALL write_to_field_dp_0D( filename, 'refreezing', region%int_refreezing)
       CALL write_to_field_dp_0D( filename, 'runoff',     region%int_runoff)
    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM(C%choice_SMB_model) // '"!')
    END IF

    ! Climate matrix
    IF (C%choice_climate_model == 'matrix') THEN
       CALL write_to_field_dp_0D( filename, 'w_EXT',   region%climate%matrix%w_EXT)
       CALL write_to_field_dp_0D( filename, 'w_tot_P', region%climate%matrix%w_tot_P)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_regional_scalar_file

  SUBROUTINE write_to_global_scalar_file( filename, global_data, time)
    ! Write model output to the global scalar file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_global_scalar_data),       INTENT(INOUT) :: global_data
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_global_scalar_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Write new time to file (thus extending the time dimension by one frame, making room for the new model data)
    CALL write_time_to_file( filename, time)

    ! Write global variables
    ! ================
    CALL write_to_field_dp_0D( filename, 'GMSL',     global_data%GMSL)
    CALL write_to_field_dp_0D( filename, 'GMSL_NAM', global_data%GMSL_NAM)
    CALL write_to_field_dp_0D( filename, 'GMSL_EAS', global_data%GMSL_EAS)
    CALL write_to_field_dp_0D( filename, 'GMSL_GRL', global_data%GMSL_GRL)
    CALL write_to_field_dp_0D( filename, 'GMSL_ANT', global_data%GMSL_ANT)

    ! CO2
    IF     (C%choice_forcing_method == 'none') THEN
    ELSEIF (C%choice_forcing_method == 'CO2_direct') THEN
      CALL write_to_field_dp_0D( filename, 'CO2_obs', global_data%CO2_obs)
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
!      CALL write_to_field_dp_0D( filename, 'CO2_obs', global_data%CO2_obs)
      CALL write_to_field_dp_0D( filename, 'CO2_mod', global_data%CO2_mod)
    ELSE
      CALL crash('unknown choice_forcing_method "' // TRIM(C%choice_forcing_method) // '"!')
    END IF

    ! d18O
    IF     (C%do_calculate_benthic_d18O) THEN
      CALL write_to_field_dp_0D( filename, 'dT_glob',  global_data%dT_glob)
      CALL write_to_field_dp_0D( filename, 'dT_dw',    global_data%dT_dw)
      CALL write_to_field_dp_0D( filename, 'd18O_mod', global_data%d18O_mod)
      CALL write_to_field_dp_0D( filename, 'd18O_ice', global_data%d18O_ice)
      CALL write_to_field_dp_0D( filename, 'd18O_Tdw', global_data%d18O_Tdw)
      CALL write_to_field_dp_0D( filename, 'd18O_NAM', global_data%d18O_NAM)
      CALL write_to_field_dp_0D( filename, 'd18O_EAS', global_data%d18O_EAS)
      CALL write_to_field_dp_0D( filename, 'd18O_GRL', global_data%d18O_GRL)
      CALL write_to_field_dp_0D( filename, 'd18O_ANT', global_data%d18O_ANT)
    END IF

    ! Computation time for different model components
    CALL write_to_field_dp_0D( filename, 'tcomp_total',   global_data%tcomp_total)
    CALL write_to_field_dp_0D( filename, 'tcomp_ice',     global_data%tcomp_ice)
    CALL write_to_field_dp_0D( filename, 'tcomp_thermo',  global_data%tcomp_thermo)
    CALL write_to_field_dp_0D( filename, 'tcomp_climate', global_data%tcomp_climate)
    CALL write_to_field_dp_0D( filename, 'tcomp_GIA',     global_data%tcomp_GIA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_global_scalar_file

    ! Write data to a grid output file
  SUBROUTINE write_to_field_dp_0D( filename, field_name_options, d)
    ! Write output data to a scalar field

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp),                            INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_dp_0D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    TYPE(type_grid)                                    :: grid
    REAL(dp), DIMENSION(:    ), POINTER                ::  d_grid_with_time
    INTEGER                                            :: wd_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( 1, d_grid_with_time, wd_grid_with_time)

    ! Copy data
    d_grid_with_time( 1) = d
    CALL sync

    ! Write data to the variable
    CALL write_var_dp_1D( filename, id_var, d_grid_with_time, start = (/ ti /), count = (/ 1 /) )

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_dp_0D

  SUBROUTINE create_regional_scalar_file(region_name, filename)
    ! Create an empty regional scalar file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    CHARACTER(LEN=*),                    INTENT(OUT)   :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_regional_scalar_file'
    LOGICAL                                            :: file_exists

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Generate filename
    filename = TRIM(C%output_dir) // '/scalar_output_' // TRIM(region_name) // '.nc'

    ! Create a new scalar file if none exists and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(filename))
    IF (par%master) THEN
      IF (file_exists) THEN
        CALL crash('file "' // TRIM( filename) // '" already exists!')
      END IF
    END IF

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Add time dimension
    CALL add_time_dimension_to_file(  filename)

    ! Create region variables
    ! ================
    CALL add_field_dp_0D( filename, 'ice_volume',    long_name='Ice volume', units='m.s.l.e')
    CALL add_field_dp_0D( filename, 'ice_volume_af', long_name='Ice volume above flotation', units='m.s.l.e')
    CALL add_field_dp_0D( filename, 'ice_area',      long_name='Ice volume', units='km^2')
    CALL add_field_dp_0D( filename, 'T2m',           long_name='Regionally averaged annual mean surface temperature', units='K')
    CALL add_field_dp_0D( filename, 'SMB',           long_name='Ice-sheet integrated surface mass balance', units='Gigaton yr^-1')
    CALL add_field_dp_0D( filename, 'BMB',           long_name='Ice-sheet integrated basal mass balance', units='Gigaton yr^-1')
    CALL add_field_dp_0D( filename, 'MB',            long_name='Ice-sheet integrated mass balance', units='Gigaton yr^-1')

    ! Individual SMB components
    IF     (C%choice_SMB_model == 'uniform' .OR. &
            C%choice_SMB_model == 'idealised' .OR. &
            C%choice_SMB_model == 'direct_global' .OR. &
            C%choice_SMB_model == 'direct_regional' .OR. &
            C%choice_SMB_model == 'snapshot') THEN
      ! Do nothing
    ELSEIF (C%choice_SMB_model == 'IMAU-ITM' .OR. &
            C%choice_SMB_model == 'IMAU-ITM_wrongrefreezing') THEN
       CALL add_field_dp_0D( filename, 'snowfall',   long_name='Ice-sheet integrated snowfall', units='Gigaton yr^-1')
       CALL add_field_dp_0D( filename, 'rainfall',   long_name='Ice-sheet integrated rainfall', units='Gigaton yr^-1')
       CALL add_field_dp_0D( filename, 'melt',       long_name='Ice-sheet integrated melt', units='Gigaton yr^-1')
       CALL add_field_dp_0D( filename, 'refreezing', long_name='Ice-sheet integrated refreezing', units='Gigaton yr^-1')
       CALL add_field_dp_0D( filename, 'runoff',     long_name='Ice-sheet integrated runoff', units='Gigaton yr^-1')
    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM(C%choice_SMB_model) // '"!')
    END IF

    ! Climate matrix
    IF (C%choice_climate_model == 'matrix') THEN
       CALL add_field_dp_0D( filename, 'w_EXT',   long_name='Weight factor from external forcing (CO2 + QTOA)', units='N/A')
       CALL add_field_dp_0D( filename, 'w_tot_P', long_name='Weigth factor for precipitation from total topography change', units='N/A')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_regional_scalar_file

  SUBROUTINE create_global_scalar_file( filename)
    ! Create an empty global scalar file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(OUT)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_global_scalar_file'
    LOGICAL                                            :: file_exists

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new scalar file if none exists and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    filename = TRIM(C%output_dir) // '/scalar_output_global.nc'
    INQUIRE(EXIST=file_exists, FILE = TRIM(filename))
    IF (par%master) THEN
      IF (file_exists) THEN
        CALL crash('file "' // TRIM( filename) // '" already exists!')
      END IF
    END IF

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Add time dimension
    CALL add_time_dimension_to_file(  filename)

    ! Create global variables
    ! ================
    CALL add_field_dp_0D( filename, 'GMSL',     long_name='Global mean sea level change', units='m')
    CALL add_field_dp_0D( filename, 'GMSL_NAM', long_name='Global mean sea level change from ice in North America', units='m')
    CALL add_field_dp_0D( filename, 'GMSL_EAS', long_name='Global mean sea level change from ice in Eurasia',       units='m')
    CALL add_field_dp_0D( filename, 'GMSL_GRL', long_name='Global mean sea level change from ice in Greenland',       units='m')
    CALL add_field_dp_0D( filename, 'GMSL_ANT', long_name='Global mean sea level change from ice in Antarctica',       units='m')

    ! CO2
    IF     (C%choice_forcing_method == 'none') THEN
    ELSEIF (C%choice_forcing_method == 'CO2_direct') THEN
      CALL add_field_dp_0D( filename, 'CO2_obs',  long_name='Observed atmospheric CO2 concentration', units='ppm')
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! CALL add_field_dp_0D( filename, 'CO2_obs', long_name='Observed atmospheric CO2 concentration', units='ppm')
      CALL add_field_dp_0D( filename, 'CO2_mod', long_name='Modelled atmospheric CO2 concentration', units='ppm')
    ELSE
      CALL crash('unknown choice_forcing_method "' // TRIM(C%choice_forcing_method) // '"!')
    END IF

    ! d18O
    IF     (C%do_calculate_benthic_d18O) THEN
      CALL add_field_dp_0D( filename, 'dT_glob',  long_name='Global annual mean surface temperature change', units='K')
      CALL add_field_dp_0D( filename, 'dT_dw',    long_name='Deep-water temperature change', units='K')
      CALL add_field_dp_0D( filename, 'd18O_mod', long_name='Modelled benthic d18O', units='per mil')
      CALL add_field_dp_0D( filename, 'd18O_ice', long_name='Modelled benthic d18O from global ice volume', units='per mil')
      CALL add_field_dp_0D( filename, 'd18O_Tdw', long_name='Modelled benthic d18O from deep-water temperature', units='per mil')
      CALL add_field_dp_0D( filename, 'd18O_NAM', long_name='Modelled benthic d18O from ice in North America', units='per mil')
      CALL add_field_dp_0D( filename, 'd18O_EAS', long_name='Modelled benthic d18O from ice in Eurasia', units='per mil')
      CALL add_field_dp_0D( filename, 'd18O_GRL', long_name='Modelled benthic d18O from ice in Greenland', units='per mil')
      CALL add_field_dp_0D( filename, 'd18O_ANT', long_name='Modelled benthic d18O from ice in Antarctica', units='per mil')
    END IF

    ! Computation time for different model components
    CALL add_field_dp_0D( filename, 'tcomp_total',   long_name='Total computation time', units='s')
    CALL add_field_dp_0D( filename, 'tcomp_ice',     long_name='Total computation time for ice dynamics', units='s')
    CALL add_field_dp_0D( filename, 'tcomp_thermo',  long_name='Total computation time for thermodynamics', units='s')
    CALL add_field_dp_0D( filename, 'tcomp_climate', long_name='Total computation time for climate+SMB+BMB', units='s')
    CALL add_field_dp_0D( filename, 'tcomp_GIA',     long_name='Total computation time for GIA', units='s')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_global_scalar_file

  SUBROUTINE add_field_dp_0D( filename, var_name, long_name, units)
    ! Add a 1-D variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_dp_0D'
    INTEGER                                            :: id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_time( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_dp_0D

  SUBROUTINE add_field_history_dp_1D( filename, var_name, ntime_history, long_name, units)
    ! Add a 1-D variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    INTEGER,                             INTENT(IN)    :: ntime_history
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_history_dp_1D'
    INTEGER                                            :: id_dim_time_history, id_dim_time, id_var
    CHARACTER(LEN=256)                                 :: var_name_time_history
    INTEGER                                            :: i
    REAL(dp), DIMENSION( :), POINTER                   :: time_history
    INTEGER                                            :: wtime_history

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Variable name for the time_history dimension
    var_name_time_history = 'time_' // TRIM( var_name)

    CALL allocate_shared_dp_1D(      ntime_history, time_history, wtime_history )

    DO i = 1,ntime_history
        time_history( i) = -i * C%dt_coupling
    END DO

    ! Add the time_history to file
    CALL add_time_history_dimension_to_file( filename, var_name_time_history, time_history)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_time_history( filename, var_name_time_history)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, var_name_time_history, id_dim_time_history)
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_time_history == -1) CALL crash('no time_history dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_time_history,  id_dim_time/), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Clean up after yourself
    CALL deallocate_shared(wtime_history)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_history_dp_1D

  SUBROUTINE create_restart_file_grid( filename, grid)
    ! Create an empty restart file with the specified grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN) :: filename
    TYPE(type_grid),                     INTENT(IN) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_restart_file_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Set up the x and y grid in figure
    CALL setup_xy_grid_in_netcdf_file( filename, grid)

    ! Add time, zeta, and month dimensions
    CALL add_time_dimension_to_file(  filename)
    CALL add_zeta_dimension_to_file(  filename)
    CALL add_month_dimension_to_file( filename)

    ! Add a time_history dimension for the inverse 18O simulation
    !IF (C%choice_forcing_method == 'd18O_inverse_CO2' .OR. &
    !   (C%choice_forcing_method == 'd18O_inverse_dT_glob')) THEN
    !    CALL add_time_history_dimension_to_file(  filename)
    !END IF

    ! Create variables
    ! ================

    ! Geometry
    CALL add_field_grid_dp_2D( filename, get_first_option_from_list( field_name_options_Hi ), long_name = 'Ice thickness'      , units = 'm')
    CALL add_field_grid_dp_2D( filename, get_first_option_from_list( field_name_options_Hb ), long_name = 'Bedrock elevation'  , units = 'm w.r.t. PD sea level')
    CALL add_field_grid_dp_2D( filename, get_first_option_from_list( field_name_options_Hs ), long_name = 'Surface elevation'  , units = 'm w.r.t. PD sea level')

    ! Thermodynamics
    CALL add_field_grid_dp_3D( filename, get_first_option_from_list( field_name_options_Ti ), long_name = 'Englacial temperature'    , units = 'K')

    ! GIA
    CALL add_field_grid_dp_2D( filename, get_first_option_from_list( field_name_options_SL ), long_name = 'Sea surface change' , units = 'm')
    CALL add_field_grid_dp_2D( filename, get_first_option_from_list( field_name_options_dHB), long_name = 'Bedrock deformation' , units = 'm')

    ! Velocities
    IF     (C%choice_ice_dynamics == 'SIA/SSA' .OR. C%choice_ice_dynamics == 'SSA') THEN
      CALL add_field_grid_dp_2D( filename, 'u_SSA_cx_a', long_name = 'SSA velocities in u direction' , units = 'm/yr')
      CALL add_field_grid_dp_2D( filename, 'v_SSA_cy_a', long_name = 'SSA velocities in v direction' , units = 'm/yr')
    ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
      CALL add_field_grid_dp_2D( filename, 'u_vav_cx_a', long_name = 'vav velocities in u direction' , units = 'm/yr')
      CALL add_field_grid_dp_2D( filename, 'v_vav_cy_a', long_name = 'vav velocities in v direction' , units = 'm/yr')
    ENDIF

    ! SMB
    IF     (C%choice_SMB_model == 'uniform') THEN
    ELSEIF (C%choice_SMB_model == 'idealised') THEN
    ELSEIF (C%choice_SMB_model == 'IMAU-ITM') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'FirnDepth', long_name = 'Monthly mean firn layer depth', units='m water equivalent')
      CALL add_field_grid_dp_2D( filename, 'MeltPreviousYear',  long_name = 'Melt during previous year' , units = 'mie')
    ELSEIF (C%choice_SMB_model == 'direct_global') THEN
    ELSEIF (C%choice_SMB_model == 'direct_regional') THEN
    ELSEIF (C%choice_SMB_model == 'snapshot') THEN
    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM(C%choice_SMB_model) // '"!')
    END IF

    ! Isotopes
    IF     (C%choice_ice_isotopes_model == 'none') THEN
    ELSEIF (C%choice_ice_isotopes_model == 'uniform') THEN
      CALL add_field_grid_dp_2D( filename, 'IsoIce', long_name='Vertically averaged 18O content', units='per mille')
    ELSEIF (C%choice_ice_isotopes_model == 'ANICE_legacy') THEN
      CALL add_field_grid_dp_2D( filename, 'IsoIce', long_name='Vertically averaged 18O content', units='per mille')
    ELSE
      CALL crash('unknown choice_ice_isotopes_model "' // TRIM(C%choice_ice_isotopes_model) // '"!')
    END IF

    ! Inverse routine data
    IF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
        CALL add_field_history_dp_1D( filename, 'CO2_inverse_history' ,  CEILING( C%CO2_inverse_averaging_window     / C%dt_coupling), long_name='inverse history of CO2' ,                    units='ppm')
        CALL add_field_history_dp_1D( filename, 'dT_glob_history' ,       CEILING( C%dT_deepwater_averaging_window / C%dt_coupling),   long_name= 'history of deep water temperature' , units='K')
    END IF

    IF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
        CALL add_field_history_dp_1D( filename, 'dT_glob_inverse_history' , CEILING( C%dT_glob_inverse_averaging_window / C%dt_coupling), long_name='inverse history of deep water temperature' , units='K')
    END IF


    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_grid

  SUBROUTINE create_help_fields_file_grid( filename, grid)
    ! Create an empty help fields file with the specified grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN) :: filename
    TYPE(type_grid),                     INTENT(IN) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_help_fields_file_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Set up the grids in this file
    CALL setup_xy_grid_in_netcdf_file( filename, grid)

    ! Add time, zeta, and month dimensions
    CALL add_time_dimension_to_file(  filename)
    CALL add_zeta_dimension_to_file(  filename)
    CALL add_month_dimension_to_file( filename)

    ! Create variables
    ! ================

    CALL create_help_field_grid( filename, C%help_field_01)
    CALL create_help_field_grid( filename, C%help_field_02)
    CALL create_help_field_grid( filename, C%help_field_03)
    CALL create_help_field_grid( filename, C%help_field_04)
    CALL create_help_field_grid( filename, C%help_field_05)
    CALL create_help_field_grid( filename, C%help_field_06)
    CALL create_help_field_grid( filename, C%help_field_07)
    CALL create_help_field_grid( filename, C%help_field_08)
    CALL create_help_field_grid( filename, C%help_field_09)
    CALL create_help_field_grid( filename, C%help_field_10)
    CALL create_help_field_grid( filename, C%help_field_11)
    CALL create_help_field_grid( filename, C%help_field_12)
    CALL create_help_field_grid( filename, C%help_field_13)
    CALL create_help_field_grid( filename, C%help_field_14)
    CALL create_help_field_grid( filename, C%help_field_15)
    CALL create_help_field_grid( filename, C%help_field_16)
    CALL create_help_field_grid( filename, C%help_field_17)
    CALL create_help_field_grid( filename, C%help_field_18)
    CALL create_help_field_grid( filename, C%help_field_19)
    CALL create_help_field_grid( filename, C%help_field_20)
    CALL create_help_field_grid( filename, C%help_field_21)
    CALL create_help_field_grid( filename, C%help_field_22)
    CALL create_help_field_grid( filename, C%help_field_23)
    CALL create_help_field_grid( filename, C%help_field_24)
    CALL create_help_field_grid( filename, C%help_field_25)
    CALL create_help_field_grid( filename, C%help_field_26)
    CALL create_help_field_grid( filename, C%help_field_27)
    CALL create_help_field_grid( filename, C%help_field_28)
    CALL create_help_field_grid( filename, C%help_field_29)
    CALL create_help_field_grid( filename, C%help_field_30)
    CALL create_help_field_grid( filename, C%help_field_31)
    CALL create_help_field_grid( filename, C%help_field_32)
    CALL create_help_field_grid( filename, C%help_field_33)
    CALL create_help_field_grid( filename, C%help_field_34)
    CALL create_help_field_grid( filename, C%help_field_35)
    CALL create_help_field_grid( filename, C%help_field_36)
    CALL create_help_field_grid( filename, C%help_field_37)
    CALL create_help_field_grid( filename, C%help_field_38)
    CALL create_help_field_grid( filename, C%help_field_39)
    CALL create_help_field_grid( filename, C%help_field_40)
    CALL create_help_field_grid( filename, C%help_field_41)
    CALL create_help_field_grid( filename, C%help_field_42)
    CALL create_help_field_grid( filename, C%help_field_43)
    CALL create_help_field_grid( filename, C%help_field_44)
    CALL create_help_field_grid( filename, C%help_field_45)
    CALL create_help_field_grid( filename, C%help_field_46)
    CALL create_help_field_grid( filename, C%help_field_47)
    CALL create_help_field_grid( filename, C%help_field_48)
    CALL create_help_field_grid( filename, C%help_field_49)
    CALL create_help_field_grid( filename, C%help_field_50)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_fields_file_grid

  SUBROUTINE create_help_field_grid( filename, field_name)
    ! Add a data field to the help_fields file

    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_help_field_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    ELSEIF (field_name == 'lon') THEN
      CALL warning('longitude is already by default written to output!')
    ELSEIF (field_name == 'lat') THEN
      CALL warning('latitude is already by default written to output!')

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL add_field_grid_dp_2D_notime( filename, 'GHF', long_name = 'Geothermal heat flux', units = 'J m^-2 yr^-1')

    ! Ice basins
    ! NOTE: needs new function that reads in integers instead of dp
    !ELSEIF (field_name == 'basin_ID') THEN
    !  CALL add_field_grid_dp_2D_notime( filename, 'basin_ID', long_name = 'Basin ID', units = '')

    ! Basal inversion target velocity
    ELSEIF (field_name == 'BIV_target_velocity') THEN
      CALL add_field_grid_dp_2D_notime( filename, 'BIV_target_velocity', long_name='Basal inversion target velocity', units='m yr^-1')

    ! Forcing climates
    ELSEIF (field_name == 'GCM_Warm_T2m') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Warm_T2m', long_name = 'Warm monthly mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'GCM_Warm_Precip') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Warm_Precip', long_name='Warm monthly total precipitation', units='mm')
    ELSEIF (field_name == 'GCM_Cold_T2m') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Cold_T2m', long_name='Cold monthly mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'GCM_Cold_Precip') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Cold_Precip', long_name = 'Cold monthly total precipitation', units='mm')
    ELSEIF (field_name == 'GCM_PI_T2m') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'PI_T2m', long_name= 'Ref PI monthly mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'GCM_PI_Precip') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'PI_Precip', long_name = 'Ref PI monthly total precipitation', units='mm')

    ! Climate matrix
    ELSEIF (field_name == 'w_ice_T') THEN
      CALL add_field_grid_dp_2D( filename, 'w_ice_T', long_name = 'Temperature weight factor based on regional albedo and insolation', units='N/A')
    ELSEIF (field_name == 'w_ins_T') THEN
      CALL add_field_grid_dp_2D( filename, 'w_ins_T', long_name = 'Temperature weight factor based on local absorbed insolation', units='N/A')
    ELSEIF (field_name == 'w_tot_T') THEN
      CALL add_field_grid_dp_2D( filename, 'w_tot_T', long_name = 'Temperature weight factor based on local absorbed insolation', units='N/A')
    ELSEIF (field_name == 'w_warm_P') THEN
      CALL add_field_grid_dp_2D( filename, 'w_warm_P', long_name = 'Precipitation weight factor contribution from GCM_warm snapshot', units='N/A')
    ELSEIF (field_name == 'w_cold_P') THEN
      CALL add_field_grid_dp_2D( filename, 'w_cold_P', long_name = 'Precipitation weight factor contribution from GCM_cold snapshot', units='N/A')

    ! Forcing oceans
    ELSEIF (field_name == 'GCM_Warm_T_ocean_3D') THEN
      CALL add_field_grid_dp_3D_notime( filename, 'Warm_T_ocean_3D',    long_name='Warm 3-D ocean temperature', units='K')
    ELSEIF (field_name == 'GCM_Warm_S_ocean_3D') THEN
      CALL add_field_grid_dp_3D_notime( filename, 'Warm_S_ocean_3D',    long_name='Warm 3-D ocean salinity', units='PSU')
    ELSEIF (field_name == 'GCM_Cold_T_ocean_3D') THEN
      CALL add_field_grid_dp_3D_notime( filename, 'Cold_T_ocean_3D',    long_name='Cold 3-D ocean temperature', units='K')
    ELSEIF (field_name == 'GCM_Cold_S_ocean_3D') THEN
      CALL add_field_grid_dp_3D_notime( filename, 'Cold_S_ocean_3D',    long_name='Cold 3-D ocean salinity', units='PSU')
    ELSEIF (field_name == 'GCM_PI_T_ocean_3D') THEN
      CALL add_field_grid_dp_3D_notime( filename, 'Ref_PI_T_ocean_3D',  long_name='Ref PI 3-D ocean temperature', units='K')
    ELSEIF (field_name == 'GCM_PI_S_ocean_3D') THEN
      CALL add_field_grid_dp_3D_notime( filename, 'Ref_PI_S_ocean_3D',  long_name='Ref PI 3-D ocean salinity', units='PSU')
    ELSEIF (field_name == 'PD_obs_T_ocean_3D') THEN
      CALL add_field_grid_dp_3D_notime( filename, 'Base_PD_T_ocean_3D', long_name='Base PD 3-D ocean temperature', units='K')
    ELSEIF (field_name == 'PD_obs_S_ocean_3D') THEN
      CALL add_field_grid_dp_3D_notime( filename, 'Base_PD_S_ocean_3D', long_name='Base PD 3-D ocean salinity', units='PSU')

    ! Fields with a time dimension
    ! ============================


    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL add_field_grid_dp_2D( filename, 'Hi', long_name = 'Ice thickness', units = 'm')
    ELSEIF (field_name == 'Hb') THEN
      CALL add_field_grid_dp_2D( filename, 'Hb', long_name = 'Bedrock elevation', units = 'm w.r.t. PD sea-level')
    ELSEIF (field_name == 'Hs') THEN
      CALL add_field_grid_dp_2D( filename, 'Hs', long_name = 'Surface elevation', units = 'm w.r.t. PD sea-level')
    ELSEIF (field_name == 'SL') THEN
      CALL add_field_grid_dp_2D( filename, 'SL', long_name = 'Geoid elevation', units = 'm w.r.t. PD sea-level')
    ELSEIF (field_name == 'dHs_dx') THEN
      CALL add_field_grid_dp_2D( filename, 'dHs_dx', long_name = 'Surface slope in x-direction', units = 'm/m')
    ELSEIF (field_name == 'dHs_dy') THEN
      CALL add_field_grid_dp_2D( filename, 'dHs_dy', long_name = 'Surface slope in y-direction', units = 'm/m')
    ELSEIF (field_name == 'dHi') THEN
      CALL add_field_grid_dp_2D( filename, 'dHi', long_name = 'Ice thickness difference w.r.t. PD', units = 'm')
    ELSEIF (field_name == 'dHs') THEN
      CALL add_field_grid_dp_2D( filename, 'dHs', long_name = 'Surface elevation difference w.r.t. PD', units = 'm')
    ELSEIF (field_name == 'dHi_dt') THEN
      CALL add_field_grid_dp_2D( filename, 'dHi_dt', long_name = 'Ice thickness rate of change', units = 'm/yr')
    ELSEIF (field_name == 'dHi_dt_target') THEN
      CALL add_field_grid_dp_2D( filename, 'dHi_dt_target', long_name = 'Target ice thickness rate of change', units = 'm/yr')
    ELSEIF (field_name == 'TAF') THEN
      CALL add_field_grid_dp_2D( filename, 'TAF', long_name = 'Ice thickness above floatation', units = 'm')
    ELSEIF (field_name == 'TAF_rel') THEN
      CALL add_field_grid_dp_2D( filename, 'TAF_rel', long_name = 'Ice thickness above floatation relative to present day', units = 'm')

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL add_field_grid_dp_3D( filename, 'Ti', long_name = 'Englacial temperature', units = 'K')
    ELSEIF (field_name == 'Cpi') THEN
      CALL add_field_grid_dp_3D( filename, 'Cpi', long_name = 'Ice heat capacity', units = 'J kg^-1 K^-1')
    ELSEIF (field_name == 'Ki') THEN
      CALL add_field_grid_dp_3D( filename, 'Ki', long_name = 'Ice thermal conductivity', units = 'J m^-1 K^-1 yr^-1')
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL add_field_grid_dp_2D( filename, 'Ti_basal', long_name = 'Ice basal temperature', units = 'K')
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL add_field_grid_dp_3D( filename, 'Ti_pmp', long_name = 'Englacial pressure melting point temperature', units = 'K')
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL add_field_grid_dp_3D( filename, 'A_flow_3D', long_name = 'Ice flow factor', units = 'Pa^-3 y^-1')
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL add_field_grid_dp_2D( filename, 'A_flow_vav', long_name = 'Vertically averaged ice flow factor', units = 'Pa^-3 y^-1')
    ELSEIF (field_name == 'Ti_base_rel') THEN
      CALL add_field_grid_dp_2D( filename, 'Ti_base_rel', long_name='Basal temperature relative to pressure melting point', units='K')

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL add_field_grid_dp_3D( filename, 'u_3D', long_name = '3D ice x-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'v_3D') THEN
      CALL add_field_grid_dp_3D( filename, 'v_3D', long_name = '3D ice y-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'w_3D') THEN
      CALL add_field_grid_dp_3D( filename, 'w_3D', long_name = '3D ice z-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'u_vav') THEN
      CALL add_field_grid_dp_2D( filename, 'u_vav', long_name = 'Vertically averaged ice x-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'v_vav') THEN
      CALL add_field_grid_dp_2D( filename, 'v_vav', long_name = 'Vertically averaged ice y-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL add_field_grid_dp_2D( filename, 'uabs_vav', long_name = 'Vertically averaged ice velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'u_surf') THEN
      CALL add_field_grid_dp_2D( filename, 'u_surf', long_name = 'Surface ice x-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'v_surf') THEN
      CALL add_field_grid_dp_2D( filename, 'v_surf', long_name = 'Surface ice y-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL add_field_grid_dp_2D( filename, 'uabs_surf', long_name = 'Surface ice velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'u_base') THEN
      CALL add_field_grid_dp_2D( filename, 'u_base', long_name = 'Basal ice x-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'v_base') THEN
      CALL add_field_grid_dp_2D( filename, 'v_base', long_name = 'Basal ice y-velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'uabs_base') THEN
      CALL add_field_grid_dp_2D( filename, 'uabs_base', long_name = 'Basal ice velocity', units = 'm yr^-1')
    ELSEIF (field_name == 'R_shear') THEN
      CALL add_field_grid_dp_2D( filename, 'R_shear', long_name='Shearing ratio; 1 = full shearing, 0 = full sliding')

    ! Climate
    ELSEIF (field_name == 'Q_TOA') THEN
      CALL add_field_grid_dp_3D( filename, 'Q_TOA', long_name='Monthly mean top-of-atmosphere insolation', units='J s^-1 m^-2')
    ELSEIF (field_name == 'T2m') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'T2m', long_name = 'Monthly mean 2-m air temperature', units = 'K')
    ELSEIF (field_name == 'T2m_year') THEN
      CALL add_field_grid_dp_2D( filename, 'T2m_year', long_name = 'Annual mean 2-m air temperature', units = 'K')
    ELSEIF (field_name == 'Precip') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Precip', long_name = 'Monthly total precipitation', units = 'mm')
    ELSEIF (field_name == 'Precip_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Precip_year', long_name = 'Annual total precipitation', units = 'mm')
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Wind_WE', long_name = 'Monthly mean zonal wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Wind_WE_year', long_name = 'Annual mean zonal wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Wind_SN', long_name = 'Monthly mean meridional wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Wind_SN_year', long_name = 'Annual mean meridional wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_LR') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Wind_LR', long_name = 'Monthly mean x-wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_LR_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Wind_LR_year', long_name = 'Annual mean x-wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_DU') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Wind_DU', long_name = 'Monthly mean y-wind', units = 'm s^-1')
    ELSEIF (field_name == 'Wind_DU_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Wind_DU_year', long_name = 'Annual mean y-wind', units = 'm s^-1')

    ! Surface mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'SMB', long_name = 'Monthly surface mass balance', units = 'm ice equivalent')
    ELSEIF (field_name == 'SMB_year') THEN
      CALL add_field_grid_dp_2D( filename, 'SMB_year', long_name = 'Annual surface mass balance', units = 'm ice equivalent')
    ELSEIF (field_name == 'Snowfall') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Snowfall', long_name = 'Monthly total snowfall', units = 'm water equivalent')
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Snowfall_year', long_name = 'Annual total snowfall', units = 'm water equivalent')
    ELSEIF (field_name == 'Rainfall') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Rainfall', long_name = 'Monthly total rainfall', units = 'm water equivalent')
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Rainfall_year', long_name = 'Annual total rainfall', units = 'm water equivalent')
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'AddedFirn', long_name = 'Monthly total added firn', units = 'm water equivalent')
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL add_field_grid_dp_2D( filename, 'AddedFirn_year', long_name = 'Annual total added firn', units = 'm water equivalent')
    ELSEIF (field_name == 'Refreezing') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Refreezing', long_name = 'Monthly total refreezing', units = 'm water equivalent')
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Refreezing_year', long_name = 'Annual total refreezing', units = 'm water equivalent')
    ELSEIF (field_name == 'Melt') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Melt', long_name = 'Monthly total melt', units = 'm water equivalent')
    ELSEIF (field_name == 'Melt_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Melt_year', long_name = 'Annual total melt', units = 'm water equivalent')
    ELSEIF (field_name == 'Runoff') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Runoff', long_name = 'Monthly total runoff', units = 'm water equivalent')
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Runoff_year', long_name = 'Annual total runoff', units = 'm water equivalent')
    ELSEIF (field_name == 'Albedo') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'Albedo', long_name = 'Monthly mean albedo')
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL add_field_grid_dp_2D( filename, 'Albedo_year', long_name = 'Annual mean albedo')
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL add_field_grid_dp_2D_monthly( filename, 'FirnDepth', long_name = 'Monthly mean firn layer depth', units='m water equivalent')
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL add_field_grid_dp_2D( filename, 'FirnDepth_year', long_name = 'Annual mean firn layer depth', units='m water equivalent')

    ! Basal mass balance and Oceans
    ELSEIF (field_name == 'BMB') THEN
      CALL add_field_grid_dp_2D( filename, 'BMB', long_name = 'Annual basal mass balance', units = 'm ice equivalent')
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL add_field_grid_dp_2D( filename, 'BMB_sheet', long_name = 'Annual basal mass balance for grounded ice', units = 'm ice equivalent')
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL add_field_grid_dp_2D( filename, 'BMB_shelf', long_name = 'Annual basal mass balance for floating ice', units = 'm ice equivalent')
    ELSEIF (field_name == 'T_ocean_base') THEN
      CALL add_field_grid_dp_2D( filename, 'T_ocean_base', long_name = 'Ocean temperature at shelf base', units = 'K')
    ELSEIF (field_name == 'T_ocean_3D') THEN
      CALL add_field_grid_dp_3D( filename, 'T_ocean_3D', long_name='3-D ocean temperature', units='K')
    ELSEIF (field_name == 'S_ocean_3D') THEN
      CALL add_field_grid_dp_3D( filename, 'S_ocean_3D', long_name='3-D ocean salinity', units='PSU')
    ELSEIF (field_name == 'dT_ocean') THEN
      CALL add_field_grid_dp_2D( filename, 'dT_ocean', long_name='Perturbation to ocean temperatures for inversion', units='K')
    ELSEIF (field_name == 'PICO_boxes') THEN
      CALL add_field_grid_int_2D( filename, 'PICO_boxes', long_name='PICO ocean boxes')

    ! Masks
    ELSEIF (field_name == 'mask') THEN
      CALL add_field_grid_int_2D( filename, 'mask', long_name = 'Mask')
    ELSEIF (field_name == 'mask_land') THEN
      CALL add_field_grid_int_2D( filename, 'mask_land', long_name = 'Land mask')
    ELSEIF (field_name == 'mask_ocean') THEN
      CALL add_field_grid_int_2D( filename, 'mask_ocean', long_name = 'Ocean mask')
    ELSEIF (field_name == 'mask_lake') THEN
      CALL add_field_grid_int_2D( filename, 'mask_lake', long_name = 'Lake mask')
    ELSEIF (field_name == 'mask_ice') THEN
      CALL add_field_grid_int_2D( filename, 'mask_ice', long_name = 'Ice mask')
    ELSEIF (field_name == 'mask_sheet') THEN
      CALL add_field_grid_int_2D( filename, 'mask_sheet', long_name = 'Sheet mask')
    ELSEIF (field_name == 'mask_shelf') THEN
      CALL add_field_grid_int_2D( filename, 'mask_shelf', long_name = 'Shelf mask')
    ELSEIF (field_name == 'mask_coast') THEN
      CALL add_field_grid_int_2D( filename, 'mask_coast', long_name = 'Coast mask')
    ELSEIF (field_name == 'mask_margin') THEN
      CALL add_field_grid_int_2D( filename, 'mask_margin', long_name = 'Margin mask')
    ELSEIF (field_name == 'mask_gl') THEN
      CALL add_field_grid_int_2D( filename, 'mask_gl', long_name = 'Grounding-line mask')
    ELSEIF (field_name == 'mask_cf') THEN
      CALL add_field_grid_int_2D( filename, 'mask_cf', long_name = 'Calving-front mask')

    ! Basal hydrology and roughness
    ELSEIF (field_name == 'pore_water_pressure') THEN
      CALL add_field_grid_dp_2D( filename, 'pore_water_pressure', long_name='pore water pressure', units='Pa')
    ELSEIF (field_name == 'effective_pressure' .OR. field_name == 'Neff') THEN
      CALL add_field_grid_dp_2D( filename, 'effective_pressure', long_name='effective basal pressure', units='Pa')
    ELSEIF (field_name == 'phi_fric') THEN
      CALL add_field_grid_dp_2D( filename, 'phi_fric', long_name = 'Till friction angle', units = 'degrees')
    ELSEIF (field_name == 'tau_yield') THEN
      CALL add_field_grid_dp_2D( filename, 'tau_yield', long_name = 'Basal yield stress', units = 'Pa')
    ELSEIF (field_name == 'alpha_sq') THEN
      CALL add_field_grid_dp_2D( filename, 'alpha_sq', long_name='Coulomb-law friction coefficient', units='unitless')
    ELSEIF (field_name == 'beta_sq') THEN
      CALL add_field_grid_dp_2D( filename, 'beta_sq', long_name='Power-law friction coefficient', units='Pa m^−1/3 yr^1/3')
    ELSEIF (field_name == 'beta') THEN
      CALL add_field_grid_dp_2D( filename, 'beta', long_name='Basal friction', units='Pa m-1 yr')
    ELSEIF (field_name == 'beta_eff') THEN
      CALL add_field_grid_dp_2D( filename, 'beta_eff', long_name='Effective friction', units='Pa m-1 yr')
    ELSEIF (field_name == 'f_grnd') THEN
      CALL add_field_grid_dp_2D( filename, 'f_grnd', long_name='Sub-grid grounded-area fraction', units='-')

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL add_field_grid_dp_2D( filename, 'iso_ice', long_name = 'Vertically averaged englacial d18O', units = 'per mille')
    ELSEIF (field_name == 'iso_surf') THEN
      CALL add_field_grid_dp_2D( filename, 'iso_surf', long_name = 'd18O of precipitation', units = 'per mille')

    ! GIA
    ELSEIF (field_name == 'dSL_dt') THEN
      CALL add_field_grid_dp_2D( filename, 'dSL_dt', long_name = 'Geoid deformation rate', units = 'm yr^-1')
    ELSEIF (field_name == 'dHb') THEN
      CALL add_field_grid_dp_2D( filename, 'dHb', long_name='Change in bedrock elevation w.r.t. PD', units='m')
    ELSEIF (field_name == 'surface_load_rel') THEN
      CALL add_field_grid_dp_2D( filename, 'surface_load_rel', long_name='surface_load_rel w.r.t. PD', units='m')
    ! Safety
    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_help_field_grid

  ! Write to restart and help fields files
  SUBROUTINE write_to_restart_file_grid( filename, region, forcing)
    ! Write model output to the grid restart file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_model_region),             INTENT(INOUT) :: region
    TYPE(type_forcing_data), OPTIONAL,   INTENT(IN)    :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_restart_file_grid'
    REAL(dp), DIMENSION( region%grid%ny, region%grid%nx) :: u_SSA_cx_a, v_SSA_cy_a, u_vav_cx_a, v_vav_cy_a

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN
      WRITE(0,'(A,F11.5,A)') '   t = ', region%time/1E3, ' kyr - writing to restart file...'
    END IF

    ! Write new time to file (thus extending the time dimension by one frame, making room for the new model data)
    CALL write_time_to_file( filename, region%time)

    ! Write model data
    ! ================

    ! Geometry
    CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, get_first_option_from_list(field_name_options_Hi)  , region%ice%Hi_a )
    CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, get_first_option_from_list(field_name_options_Hb)  , region%ice%Hb_a )
    CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, get_first_option_from_list(field_name_options_Hs)  , region%ice%Hs_a )

    ! Thermodynamics
    CALL write_to_field_multiple_options_grid_dp_3D( filename, region%grid, get_first_option_from_list(field_name_options_Ti)  , region%ice%Ti_a )

    ! GIA
    CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, get_first_option_from_list(field_name_options_SL)  , region%ice%SL_a )
    CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, get_first_option_from_list(field_name_options_dHb) , region%ice%dHb_a )

    ! Velocities
    IF     (C%choice_ice_dynamics == 'SIA/SSA' .OR. C%choice_ice_dynamics == 'SSA') THEN
      u_SSA_cx_a            = 0._dp
      v_SSA_cy_a            = 0._dp
      u_SSA_cx_a(:, 1:region%grid%nx-1) = region%ice%u_SSA_cx
      v_SSA_cy_a(1:region%grid%ny-1, :) = region%ice%v_SSA_cy
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'u_SSA_cx_a' , u_SSA_cx_a )
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'v_SSA_cy_a' , v_SSA_cy_a )
    ELSEIF     (C%choice_ice_dynamics == 'DIVA') THEN
      u_vav_cx_a            = 0._dp
      v_vav_cy_a            = 0._dp
      u_vav_cx_a(:, 1:region%grid%nx-1) = region%ice%u_vav_cx
      v_vav_cy_a(1:region%grid%ny-1, :) = region%ice%v_vav_cy
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'u_vav_cx_a' , region%ice%u_vav_cx )
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'v_vav_cy_a' , region%ice%v_vav_cy )
    END IF

    ! SMB
    IF     (C%choice_SMB_model == 'uniform') THEN
    ELSEIF (C%choice_SMB_model == 'idealised') THEN
    ELSEIF (C%choice_SMB_model == 'IMAU-ITM') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'FirnDepth',  region%SMB%FirnDepth)
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'MeltPreviousYear'  , region%SMB%MeltPreviousYear )
    ELSEIF (C%choice_SMB_model == 'direct_global') THEN
    ELSEIF (C%choice_SMB_model == 'direct_regional') THEN
    ELSEIF (C%choice_SMB_model == 'snapshot') THEN
    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM(C%choice_SMB_model) // '"!')
    END IF

    ! Isotopes
    IF     (C%choice_ice_isotopes_model == 'none') THEN
    ELSEIF (C%choice_ice_isotopes_model == 'uniform') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'IsoIce' , region%ice%IsoIce )
    ELSEIF (C%choice_ice_isotopes_model == 'ANICE_legacy') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'IsoIce' , region%ice%IsoIce )
    ELSE
      CALL crash('unknown choice_ice_isotopes_model "' // TRIM(C%choice_ice_isotopes_model) // '"!')
    END IF

    ! Inverse routine data
    IF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      IF (.NOT.(PRESENT( forcing))) CALL crash('write_to_restart_file_grid needs forcing field if d18O_inverse_CO2 is used')
      CALL write_to_field_history_dp_1D( filename, forcing%nCO2_inverse_history, 'CO2_inverse_history',  forcing%CO2_inverse_history)
      CALL write_to_field_history_dp_1D( filename, forcing%ndT_glob_history,     'dT_glob_history',      forcing%dT_glob_history)
    END IF

    IF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      IF (.NOT.(PRESENT( forcing))) CALL crash('write_to_restart_file_grid needs forcing field if d18O_inverse_dT_glob is used')
      CALL write_to_field_history_dp_1D( filename, forcing%ndT_glob_inverse_history, 'dT_glob_inverse_history' , forcing%dT_glob_inverse_history)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_grid

  SUBROUTINE write_to_help_fields_file_grid( filename, region)
    ! Write model output to the help fields file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_help_fields_file_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN
      WRITE(0,'(A,F11.5,A)') '   t = ', region%time/1E3, ' kyr - writing to help_fields file...'
    END IF

    ! Write new time to file (thus extending the time dimension by one frame, making room for the new model data)
    CALL write_time_to_file( filename, region%time)

    ! Write model data
    ! ================

    CALL write_help_field_grid( filename, region, C%help_field_01)
    CALL write_help_field_grid( filename, region, C%help_field_02)
    CALL write_help_field_grid( filename, region, C%help_field_03)
    CALL write_help_field_grid( filename, region, C%help_field_04)
    CALL write_help_field_grid( filename, region, C%help_field_05)
    CALL write_help_field_grid( filename, region, C%help_field_06)
    CALL write_help_field_grid( filename, region, C%help_field_07)
    CALL write_help_field_grid( filename, region, C%help_field_08)
    CALL write_help_field_grid( filename, region, C%help_field_09)
    CALL write_help_field_grid( filename, region, C%help_field_10)
    CALL write_help_field_grid( filename, region, C%help_field_11)
    CALL write_help_field_grid( filename, region, C%help_field_12)
    CALL write_help_field_grid( filename, region, C%help_field_13)
    CALL write_help_field_grid( filename, region, C%help_field_14)
    CALL write_help_field_grid( filename, region, C%help_field_15)
    CALL write_help_field_grid( filename, region, C%help_field_16)
    CALL write_help_field_grid( filename, region, C%help_field_17)
    CALL write_help_field_grid( filename, region, C%help_field_18)
    CALL write_help_field_grid( filename, region, C%help_field_19)
    CALL write_help_field_grid( filename, region, C%help_field_20)
    CALL write_help_field_grid( filename, region, C%help_field_21)
    CALL write_help_field_grid( filename, region, C%help_field_22)
    CALL write_help_field_grid( filename, region, C%help_field_23)
    CALL write_help_field_grid( filename, region, C%help_field_24)
    CALL write_help_field_grid( filename, region, C%help_field_25)
    CALL write_help_field_grid( filename, region, C%help_field_26)
    CALL write_help_field_grid( filename, region, C%help_field_27)
    CALL write_help_field_grid( filename, region, C%help_field_28)
    CALL write_help_field_grid( filename, region, C%help_field_29)
    CALL write_help_field_grid( filename, region, C%help_field_30)
    CALL write_help_field_grid( filename, region, C%help_field_31)
    CALL write_help_field_grid( filename, region, C%help_field_32)
    CALL write_help_field_grid( filename, region, C%help_field_33)
    CALL write_help_field_grid( filename, region, C%help_field_34)
    CALL write_help_field_grid( filename, region, C%help_field_35)
    CALL write_help_field_grid( filename, region, C%help_field_36)
    CALL write_help_field_grid( filename, region, C%help_field_37)
    CALL write_help_field_grid( filename, region, C%help_field_38)
    CALL write_help_field_grid( filename, region, C%help_field_39)
    CALL write_help_field_grid( filename, region, C%help_field_40)
    CALL write_help_field_grid( filename, region, C%help_field_41)
    CALL write_help_field_grid( filename, region, C%help_field_42)
    CALL write_help_field_grid( filename, region, C%help_field_43)
    CALL write_help_field_grid( filename, region, C%help_field_44)
    CALL write_help_field_grid( filename, region, C%help_field_45)
    CALL write_help_field_grid( filename, region, C%help_field_46)
    CALL write_help_field_grid( filename, region, C%help_field_47)
    CALL write_help_field_grid( filename, region, C%help_field_48)
    CALL write_help_field_grid( filename, region, C%help_field_49)
    CALL write_help_field_grid( filename, region, C%help_field_50)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_help_fields_file_grid

  SUBROUTINE write_help_field_grid( filename, region, field_name)
    ! Write a single data field to the help fields file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_model_region),             INTENT(INOUT) :: region
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_help_field_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (field_name == 'none') THEN
      CALL finalise_routine( routine_name)
      RETURN

    ! Fields with no time dimension
    ! =============================

    ! Lat/lon
    ELSEIF (field_name == 'lon') THEN
      CALL warning('longitude is already by default written to output!')
    ELSEIF (field_name == 'lat') THEN
      CALL warning('latitude is already by default written to output!')

    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_notime( filename, region%grid, 'GHF', region%ice%GHF_a)

    ! Ice basins
    ! NOTE: needs new function that reads in integers instead of dp
    !ELSEIF (field_name == 'basin_ID') THEN
    !  CALL write_to_field_multiple_options_grid_int_2D_notime( filename, region%grid, 'basin_ID', region%ice%basin_ID)

    ! Basal inversion target velocity
    ELSEIF (field_name == 'BIV_target_velocity') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_notime( filename, region%grid, 'BIV_target_velocity', region%ice%BIV_uabs_surf_target)

    ! Forcing climates
    ELSEIF (field_name == 'GCM_Warm_T2m') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'Warm_T2m', region%climate%matrix%GCM_warm%T2m)
    ELSEIF (field_name == 'GCM_Warm_Precip') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'Warm_Precip', region%climate%matrix%GCM_warm%Precip)
    ELSEIF (field_name == 'GCM_Cold_T2m') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'Cold_T2m', region%climate%matrix%GCM_cold%T2m)
    ELSEIF (field_name == 'GCM_Cold_Precip') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'Cold_Precip', region%climate%matrix%GCM_cold%Precip)
    ELSEIF (field_name == 'GCM_PI_T2m') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'PI_T2m', region%climate%matrix%GCM_PI%T2m)
    ELSEIF (field_name == 'GCM_PI_Precip') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'PI_Precip', region%climate%matrix%GCM_PI%Precip)

    ! Climate matrix
    ELSEIF (field_name == 'w_ice_T') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'w_ice_T', region%climate%matrix%w_ice_T)
    ELSEIF (field_name == 'w_ins_T') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'w_ins_T', region%climate%matrix%w_ins_T)
    ELSEIF (field_name == 'w_tot_T') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'w_tot_T', region%climate%matrix%w_tot_T)
    ELSEIF (field_name == 'w_warm_P') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'w_warm_P', region%climate%matrix%w_warm_P)
    ELSEIF (field_name == 'w_cold_P') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'w_cold_P', region%climate%matrix%w_cold_P)

    ! Forcing oceans
    ELSEIF (field_name == 'GCM_Warm_T_ocean_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D_notime( filename, region%grid, 'Warm_T_ocean_3D', region%ocean_matrix%GCM_Warm%T_ocean_ext)
    ELSEIF (field_name == 'GCM_Warm_S_ocean_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D_notime( filename, region%grid, 'Warm_S_ocean_3D', region%ocean_matrix%GCM_Warm%S_ocean_ext)
    ELSEIF (field_name == 'GCM_Cold_T_ocean_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D_notime( filename, region%grid, 'Cold_T_ocean_3D', region%ocean_matrix%GCM_Cold%T_ocean_ext)
    ELSEIF (field_name == 'GCM_Cold_S_ocean_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D_notime( filename, region%grid, 'Cold_S_ocean_3D', region%ocean_matrix%GCM_Cold%S_ocean_ext)
    ELSEIF (field_name == 'GCM_PI_T_ocean_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D_notime( filename, region%grid, 'Ref_PI_T_ocean_3D', region%ocean_matrix%GCM_PI%T_ocean_ext)
    ELSEIF (field_name == 'GCM_PI_S_ocean_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D_notime( filename, region%grid, 'Ref_PI_S_ocean_3D', region%ocean_matrix%GCM_PI%S_ocean_ext)
    ELSEIF (field_name == 'PD_obs_T_ocean_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D_notime( filename, region%grid, 'Base_PD_T_ocean_3D', region%ocean_matrix%PD_obs%T_ocean_ext)
    ELSEIF (field_name == 'PD_obs_S_ocean_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D_notime( filename, region%grid, 'Base_PD_S_ocean_3D', region%ocean_matrix%PD_obs%S_ocean_ext)

    ! Fields with a time dimension
    ! ============================

    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Hi', region%ice%Hi_a)
    ELSEIF (field_name == 'Hb') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Hb', region%ice%Hb_a)
    ELSEIF (field_name == 'Hs') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Hs', region%ice%Hs_a)
    ELSEIF (field_name == 'SL') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'SL', region%ice%SL_a)
    ELSEIF (field_name == 'dHs_dx') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'dHs_dx', region%ice%dHs_dx_a)
    ELSEIF (field_name == 'dHs_dy') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'dHs_dy', region%ice%dHs_dy_a)
    ELSEIF (field_name == 'dHi') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'dHi', region%ice%dHi_a)
    ELSEIF (field_name == 'dHs') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'dHs', region%ice%dHs_a)
    ELSEIF (field_name == 'dHi_dt') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'dHi_dt', region%ice%dHi_dt_a)
    ELSEIF (field_name == 'dHi_dt_target') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'dHi_dt_target', region%ice%dHi_dt_target)
    ELSEIF (field_name == 'TAF') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'TAF', region%ice%TAF_a)
    ELSEIF (field_name == 'TAF_rel') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'TAF_rel', region%ice%TAF_rel)

    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%grid, 'Ti', region%ice%Ti_a)
    ELSEIF (field_name == 'Cpi') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%grid, 'Cpi', region%ice%Cpi_a)
    ELSEIF (field_name == 'Ki') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%grid,  'Ki', region%ice%Ki_a)
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Ti_basal', region%ice%Ti_a( :,:,C%nz))
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%grid, 'Ti_pmp'  , region%ice%Ti_pmp_a)
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%grid, 'A_flow_3D', region%ice%A_flow_3D_a)
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'A_flow_vav', region%ice%A_flow_vav_a)
    ELSEIF (field_name == 'Ti_base_rel') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Ti_base_rel', (region%ice%Ti_a( C%nz,:,:) - region%ice%Ti_pmp_a( C%nz,:,:)))

    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%grid, 'u_3D', region%ice%u_3D_a)
    ELSEIF (field_name == 'v_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%grid, 'v_3D', region%ice%v_3D_a)
    ELSEIF (field_name == 'w_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%grid, 'w_3D', region%ice%w_3D_a)
    ELSEIF (field_name == 'u_vav') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'u_vav', region%ice%u_vav_a)
    ELSEIF (field_name == 'v_vav') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'v_vav', region%ice%v_vav_a)
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'uabs_vav' , region%ice%uabs_vav_a)
    ELSEIF (field_name == 'u_surf') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'u_surf'   , region%ice%u_surf_a)
    ELSEIF (field_name == 'v_surf') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'v_surf'   , region%ice%v_surf_a)
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'uabs_surf', region%ice%uabs_surf_a)
    ELSEIF (field_name == 'u_base') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'u_base'   , region%ice%u_base_a)
    ELSEIF (field_name == 'v_base') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'v_base'   , region%ice%v_base_a)
    ELSEIF (field_name == 'uabs_base') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'uabs_base', region%ice%uabs_base_a)
    ELSEIF (field_name == 'R_shear') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'R_shear', region%ice%R_shear)

    ! Climate
    ELSEIF (field_name == 'Q_TOA') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid,  'Q_TOA', region%climate%Q_TOA)
    ELSEIF (field_name == 'T2m') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid,  'T2m', region%climate%T2m)
    ELSEIF (field_name == 'T2m_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'T2m_year', SUM( region%climate%T2m,1) / 12._dp)
    ELSEIF (field_name == 'Precip') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid,  'Precip', region%climate%Precip)
    ELSEIF (field_name == 'Precip_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Precip_year', SUM( region%climate%Precip,1) )
    ELSEIF (field_name == 'Wind_LR') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'Wind_LR', region%climate%Wind_LR)
    ELSEIF (field_name == 'Wind_LR_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Wind_LR_year', SUM( region%climate%Wind_LR,1) / 12._dp)
    ELSEIF (field_name == 'Wind_DU') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'Wind_DU', region%climate%Wind_DU)
    ELSEIF (field_name == 'Wind_DU_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Wind_DU_year', SUM( region%climate%Wind_DU,1) / 12._dp)

    ! Surface mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'SMB', region%SMB%SMB)
    ELSEIF (field_name == 'SMB_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'SMB_year', region%SMB%SMB_year)
    ELSEIF (field_name == 'Snowfall') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid,  'Snowfall', region%SMB%Snowfall)
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Snowfall_year', SUM( region%SMB%Snowfall,1))
    ELSEIF (field_name == 'Rainfall') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'Rainfall', region%SMB%Rainfall)
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Rainfall_year', SUM( region%SMB%Rainfall,1))
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'AddedFirn', region%SMB%AddedFirn)
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'AddedFirn_year', SUM( region%SMB%AddedFirn,1))
    ELSEIF (field_name == 'Refreezing') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'Refreezing', region%SMB%Refreezing)
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Refreezing_year', SUM( region%SMB%Refreezing,1))
    ELSEIF (field_name == 'Melt') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'Melt', region%SMB%Melt)
    ELSEIF (field_name == 'Melt_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Melt_year', SUM( region%SMB%Melt,1))
    ELSEIF (field_name == 'Runoff') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'Runoff', region%SMB%Runoff)
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Runoff_year', SUM( region%SMB%Runoff,1))
    ELSEIF (field_name == 'Albedo') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'Albedo', region%SMB%Albedo)
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'Albedo_year', SUM( region%SMB%Albedo,1) / 12._dp)
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL write_to_field_multiple_options_grid_dp_2D_monthly( filename, region%grid, 'FirnDepth', region%SMB%FirnDepth)
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'FirnDepth_year', SUM( region%SMB%FirnDepth,1) / 12._dp)

    ! Basal mass balance
    ELSEIF (field_name == 'BMB') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'BMB', region%BMB%BMB)
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'BMB_sheet', region%BMB%BMB_sheet)
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'BMB_shelf', region%BMB%BMB_shelf)
    ELSEIF (field_name == 'T_ocean_base') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'T_ocean_base', region%BMB%T_ocean_base)
    ELSEIF (field_name == 'T_ocean_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%grid, 'T_ocean_3D', region%ocean_matrix%applied%T_ocean_corr_ext)
    ELSEIF (field_name == 'S_ocean_3D') THEN
      CALL write_to_field_multiple_options_grid_dp_3D( filename, region%grid, 'S_ocean_3D', region%ocean_matrix%applied%S_ocean_corr_ext)
    ELSEIF (field_name == 'dT_ocean') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'dT_ocean', region%ocean_matrix%applied%dT_ocean)
    ELSEIF (field_name == 'PICO_boxes') THEN
      CALL write_to_field_multiple_options_grid_int_2D( filename, region%grid, 'PICO_boxes', region%BMB%PICO_k)

    ! Masks
    ELSEIF (field_name == 'mask') THEN
      CALL write_to_field_multiple_options_grid_int_2D( filename, region%grid, 'mask', region%ice%mask_a)
    ELSEIF (field_name == 'mask_land') THEN
      CALL write_to_field_multiple_options_grid_int_2D( filename, region%grid, 'mask_land', region%ice%mask_land_a)
    ELSEIF (field_name == 'mask_ocean') THEN
      CALL write_to_field_multiple_options_grid_int_2D( filename, region%grid, 'mask_ocean', region%ice%mask_ocean_a)
    ELSEIF (field_name == 'mask_lake') THEN
      CALL write_to_field_multiple_options_grid_int_2D( filename, region%grid, 'mask_lake', region%ice%mask_lake_a)
    ELSEIF (field_name == 'mask_ice') THEN
      CALL write_to_field_multiple_options_grid_int_2D( filename, region%grid, 'mask_ice', region%ice%mask_ice_a)
    ELSEIF (field_name == 'mask_sheet') THEN
      CALL write_to_field_multiple_options_grid_int_2D( filename, region%grid, 'mask_sheet', region%ice%mask_sheet_a)
    ELSEIF (field_name == 'mask_shelf') THEN
      CALL write_to_field_multiple_options_grid_int_2D( filename, region%grid, 'mask_shelf', region%ice%mask_shelf_a)
    ELSEIF (field_name == 'mask_coast') THEN
      CALL write_to_field_multiple_options_grid_int_2D( filename, region%grid, 'mask_coast', region%ice%mask_coast_a)
    ELSEIF (field_name == 'mask_margin') THEN
      CALL write_to_field_multiple_options_grid_int_2D( filename, region%grid, 'mask_margin', region%ice%mask_margin_a)
    ELSEIF (field_name == 'mask_gl') THEN
      CALL write_to_field_multiple_options_grid_int_2D( filename, region%grid, 'mask_gl', region%ice%mask_gl_a)
    ELSEIF (field_name == 'mask_cf') THEN
      CALL write_to_field_multiple_options_grid_int_2D( filename, region%grid, 'mask_cf', region%ice%mask_cf_a)

    ! Basal conditions
    ELSEIF (field_name == 'pore_water_pressure') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'pore_water_pressure', region%ice%pore_water_pressure_a)
    ELSEIF (field_name == 'effective_pressure' .OR. field_name == 'Neff') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'effective_pressure', region%ice%Neff_a)
    ELSEIF (field_name == 'phi_fric') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'phi_fric', region%ice%phi_fric_a)
    ELSEIF (field_name == 'tau_yield') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'tau_yield', region%ice%tauc_a)
    ELSEIF (field_name == 'alpha_sq') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'alpha_sq', region%ice%alpha_sq_a)
    ELSEIF (field_name == 'beta_sq') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'beta_sq', region%ice%beta_sq_a)
    ELSEIF (field_name == 'beta') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'beta', region%ice%beta_a)
    ELSEIF (field_name == 'beta_eff') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'beta_eff', region%ice%beta_eff_a)
    ELSEIF (field_name == 'f_grnd') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'f_grnd', region%ice%f_grnd_a)

    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'iso_ice', region%ice%IsoIce)
    ELSEIF (field_name == 'iso_surf') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'iso_surf', region%ice%IsoSurf)

    ! GIA
    ELSEIF (field_name == 'dSL_dt') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'dSL_dt', region%ice%dSL_dt_a)
    ELSEIF (field_name == 'dHb') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'dHb', region%ice%dHb_a)
    ELSEIF (field_name == 'surface_load_rel') THEN
      CALL write_to_field_multiple_options_grid_dp_2D( filename, region%grid, 'surface_load_rel', region%ice%surface_load_rel)
    ! Safety
    ELSE
      CALL crash('unknown help field name "' // TRIM( field_name) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_help_field_grid

  ! Create/read an extrapolated ocean data file
  SUBROUTINE create_extrapolated_ocean_file(  hires, filename)
    ! Create a new folder extrapolated ocean data file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires
    CHARACTER(LEN=256),                  INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_extrapolated_ocean_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create a new file and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    CALL create_new_netcdf_file_for_writing( filename)

    ! Add ocean depth dimension
    CALL add_ocean_dimension_to_file(  filename)

    ! Set up the grids in this file
    CALL setup_xy_grid_in_netcdf_file( filename, hires%grid)

    ! Add T_ocean and S_ocean fields
    CALL add_field_grid_dp_3D_ocean_notime( filename, 'T_ocean', long_name = '3-D ocean temperature'  , units = 'K'  )
    CALL add_field_grid_dp_3D_ocean_notime( filename, 'S_ocean', long_name = '3-D ocean salinity'     , units = 'PSU')

    ! Write the T_ocean and S_ocean fields
    CALL write_to_field_multiple_options_grid_dp_3D_ocean_notime( filename, hires%grid, 'T_ocean', hires%T_ocean)
    CALL write_to_field_multiple_options_grid_dp_3D_ocean_notime( filename, hires%grid, 'S_ocean', hires%S_ocean)

    ! Finalise routine path
    CALL finalise_routine( routine_name)
  END SUBROUTINE create_extrapolated_ocean_file

  ! Inverted basal roughness
  SUBROUTINE create_BIV_bed_roughness_file( grid, ice)
    ! Create a new folder extrapolated ocean data file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_BIV_bed_roughness_file'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: id_var
    REAL(dp), DIMENSION(:,:), POINTER                  ::  d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, d_grid, wd_grid)

    ! Create a new file and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    filename = TRIM(C%output_dir) // TRIM(C%BIVgeo_filename_output)
    INQUIRE(EXIST=file_exists, FILE = TRIM( filename))
    IF (par%master) THEN
      IF (file_exists) THEN
        CALL crash('file "' // TRIM( filename) // '" already exists!')
      END IF
    END IF

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Writing inverted bed roughness to file "', TRIM( filename), '"...'

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Set up the grids in this file
    CALL setup_xy_grid_in_netcdf_file( filename, grid)

    ! Do it
    SELECT CASE(C%choice_sliding_law)

      CASE ('no_sliding')
        CALL crash('not defined for choice_sliding_law = "no_sliding"!')

      CASE ('Weertman')

        ! Copy data
        d_grid( :, grid%i1:grid%i2) = ice%beta_sq_a( :,grid%i1:grid%i2)
        CALL sync
        ! Transpose the output data
        CALL permute_2D_dp( d_grid, wd_grid, map = [2,1])
        ! Add field to output file
        CALL add_field_grid_dp_2D_notime( filename, 'beta_sq',  long_name = 'Power-law friction coefficient',   units = '[Pa m^−1/3 yr^1/3]')
        ! Write field to output file
        CALL inquire_var_multiple_options( filename, 'beta_sq', id_var)
        CALL write_var_dp_2D( filename, id_var, d_grid)

      CASE ('Coulomb','Coulomb_regularised','Zoet-Iverson')

        ! Copy data
        d_grid( :, grid%i1:grid%i2) = ice%phi_fric_a( :,grid%i1:grid%i2)
        CALL sync
        ! Transpose the output data
        CALL permute_2D_dp( d_grid, wd_grid, map = [2,1])
        ! Add field to output file
        CALL add_field_grid_dp_2D_notime( filename, 'phi_fric', long_name = 'Till friction angle', units = 'degrees')
        ! Write field to output file
        CALL inquire_var_multiple_options( filename, 'phi_fric', id_var)
        CALL write_var_dp_2D( filename, id_var, d_grid)

      CASE ('Tsai2015','Schoof2005')

        ! Copy data
        d_grid( :, grid%i1:grid%i2) = ice%alpha_sq_a( :,grid%i1:grid%i2)
        CALL sync
        ! Transpose the output data
        CALL permute_2D_dp( d_grid, wd_grid, map = [2,1])
        ! Add field to output file
        CALL add_field_grid_dp_2D_notime( filename, 'alpha_sq', long_name = 'Coulomb-law friction coefficient', units = '-')
        ! Write field to output file
        CALL inquire_var_multiple_options( filename, 'alpha_sq', id_var)
        CALL write_var_dp_2D( filename, id_var, d_grid)

        ! Copy data
        d_grid( :, grid%i1:grid%i2) = ice%beta_sq_a( :,grid%i1:grid%i2)
        CALL sync
        ! Transpose the output data
        CALL permute_2D_dp( d_grid, wd_grid, map = [2,1])
        ! Add field to output file
        CALL add_field_grid_dp_2D_notime( filename, 'beta_sq',  long_name = 'Power-law friction coefficient',   units = '[Pa m^−1/3 yr^1/3]')
        ! Write field to output file
        CALL inquire_var_multiple_options( filename, 'beta_sq', id_var)
        CALL write_var_dp_2D( filename, id_var, d_grid)

      CASE DEFAULT
        CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')

    END SELECT

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_BIV_bed_roughness_file

  ! Inverted ocean
  SUBROUTINE create_inverted_ocean_file( grid, ocean)
    ! Create a new folder extrapolated ocean data file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                    INTENT(IN)    :: grid
    TYPE(type_ocean_snapshot_regional), INTENT(INOUT) :: ocean
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'create_inverted_ocean_file'
    CHARACTER(LEN=256)                                :: filename
    LOGICAL                                           :: file_exists
    INTEGER                                           :: k, id_var
    REAL(dp), DIMENSION(:,:,:), POINTER               ::  d_grid
    INTEGER                                           :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( C%nz_ocean, grid%ny, grid%nx, d_grid, wd_grid)

    ! Add delta temperatures to baseline field
    DO k = 1, C%nz_ocean
      d_grid( k,:,grid%i1:grid%i2) = ocean%T_ocean_corr_ext( k,:,grid%i1:grid%i2) + ocean%dT_ocean( :,grid%i1:grid%i2)
    END DO
    CALL sync

    ! Create a new file and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    filename = TRIM(C%output_dir) // TRIM(C%inverted_ocean_filename_output)
    INQUIRE(EXIST=file_exists, FILE = TRIM( filename))
    IF (par%master) THEN
      IF (file_exists) THEN
        CALL crash('file "' // TRIM( filename) // '" already exists!')
      END IF
    END IF

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Writing inverted ocean to file "', TRIM( filename), '"...'

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Add ocean depth dimension
    CALL add_ocean_dimension_to_file(  filename)

    ! Set up the grids in this file
    CALL setup_xy_grid_in_netcdf_file( filename, grid)

    ! Add field to output file
    CALL add_field_grid_dp_3D_ocean_notime( filename, 'T_ocean', long_name = '3-D ocean temperature', units = 'degrees C')
    CALL add_field_grid_dp_3D_ocean_notime( filename, 'S_ocean', long_name = '3-D ocean salinity',    units = 'PSU')

    ! Write field to output file
    CALL inquire_var_multiple_options( filename, 'T_ocean', id_var)
    CALL inquire_var_multiple_options( filename, 'S_ocean', id_var)

    CALL write_to_field_multiple_options_grid_dp_3D_ocean_notime( filename, grid, 'T_ocean', d_grid)
    CALL write_to_field_multiple_options_grid_dp_3D_ocean_notime( filename, grid, 'S_ocean', ocean%S_ocean_corr_ext)

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_inverted_ocean_file

  ! ==== Write data to flexibly-defined fields =====
  ! ================================================

  ! Write data to a grid output file
  SUBROUTINE write_to_field_multiple_options_grid_dp_2D( filename, grid, field_name_options, d)
    ! Write a 2-D data field to a NetCDF file variable on an x/y-grid
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(:,:    ),        INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid_with_time
    INTEGER                                            :: wd_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_2D( filename, var_name, should_have_time = .TRUE.)

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%ny, grid%nx, 1, d_grid_with_time, wd_grid_with_time)

    ! Copy data
    d_grid_with_time( :, grid%i1:grid%i2, 1) = d( :,grid%i1:grid%i2)
    CALL sync

    ! Transpose the output data
    CALL permute_3D_dp( d_grid_with_time, wd_grid_with_time , map = [2,1,3])

    ! Write data to the variable
    CALL write_var_dp_3D( filename, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid%nx, grid%ny, 1 /) )

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_2D

  ! Write data to a grid output file
  SUBROUTINE write_to_field_history_dp_1D( filename, ntime_history, field_name_options, d)
    ! Write a 1-D data field to a NetCDF file variable on a time_history array
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: ntime_history
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    REAL(dp), DIMENSION(: ),             INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_history_dp_1D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:), POINTER                  ::  d_grid_with_time
    INTEGER                                            :: wd_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    ! CALL check_time_history( filename)

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( ntime_history, 1, d_grid_with_time, wd_grid_with_time)

    ! Copy data
    d_grid_with_time( :, 1) = d( :)
    CALL sync

    ! Write data to the variable
    CALL write_var_dp_2D( filename, id_var, d_grid_with_time, start = (/ 1, ti /), count = (/ ntime_history, 1 /) )

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_history_dp_1D


  SUBROUTINE write_to_field_multiple_options_grid_int_2D( filename, grid, field_name_options, d)
    ! Write a 2-D data field to a NetCDF file variable on an x/y-grid
    ! (Mind you, that's 2-D in the physical sense, so a 1-D array!)
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    INTEGER, DIMENSION(:,:    ),   INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_int_2D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER, DIMENSION(:,:,:), POINTER           ::  d_grid_with_time
    INTEGER                                            :: wd_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_int_2D( filename, var_name, should_have_time = .TRUE.)

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_int_3D( grid%ny, grid%nx, 1, d_grid_with_time, wd_grid_with_time)

    ! Copy data
    d_grid_with_time( :, grid%i1:grid%i2, 1) = d( :,grid%i1:grid%i2)
    CALL sync

    ! Transpose the output data
    CALL permute_3D_int( d_grid_with_time, wd_grid_with_time , map = [2,1,3])

    ! Write data to the variable
    CALL write_var_int_3D( filename, id_var, d_grid_with_time, start = (/ 1, 1, ti /), count = (/ grid%nx, grid%ny, 1 /) )

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_int_2D

  SUBROUTINE write_to_field_multiple_options_grid_dp_2D_monthly( filename, grid, field_name_options, d)
    ! Write a 2-D monthly data field to a NetCDF file variable on an x/y-grid
    ! (Mind you, that's 2-D monthly in the physical sense, so a 1-D array!)
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:  ),        INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_2D_monthly'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:,:,:), POINTER              :: d_grid_with_time
    REAL(dp), DIMENSION(:,:,:  ), POINTER              :: d_grid
    INTEGER                                            :: wd_grid, wd_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_2D_monthly( filename, var_name, should_have_time = .TRUE.)

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_dp_4D( grid%nx, grid%ny, 12, 1, d_grid_with_time, wd_grid_with_time)
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx,    d_grid, wd_grid)

    ! Swap the x and y grid
    d_grid( :, :, grid%i1:grid%i2) = d( :, :,grid%i1:grid%i2)
    CALL sync
    CALL permute_3D_dp( d_grid, wd_grid , map = [3,2,1])

    ! Copy data
    d_grid_with_time( grid%i1:grid%i2,:,:,1) = d_grid( grid%i1:grid%i2,:,:)
    CALL sync

    ! Write data to the variable
    CALL write_var_dp_4D( filename, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, 12, 1 /) )

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid_with_time)
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_2D_monthly

  SUBROUTINE write_to_field_multiple_options_grid_dp_3D( filename, grid, field_name_options, d)
    ! Write a 3-D data field to a NetCDF file variable on an x/y-grid
    !
    ! Write to the last time frame of the variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:  ),        INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_3D'
    INTEGER                                            :: id_var, id_dim_time, ti
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:,:,:), POINTER              :: d_grid_with_time
    REAL(dp), DIMENSION(:,:,:  ), POINTER              :: d_grid
    INTEGER                                            :: wd_grid, wd_grid_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_3D( filename, var_name, should_have_time = .TRUE.)

    ! Inquire length of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = ti)

    ! Allocate shared memory
    CALL allocate_shared_dp_4D( grid%nx, grid%ny, C%nz, 1, d_grid_with_time, wd_grid_with_time)
    CALL allocate_shared_dp_3D( C%nz, grid%ny, grid%nx,    d_grid, wd_grid)

    ! Swap the x and y grid
    d_grid( :, :, grid%i1:grid%i2) = d( :, :, grid%i1:grid%i2)
    CALL sync

    CALL permute_3D_dp( d_grid, wd_grid , map = [3,2,1])

    ! Copy data
    d_grid_with_time( grid%i1:grid%i2,:,:,1) = d_grid( grid%i1:grid%i2,:,:)
    CALL sync

    ! Write data to the variable
    CALL write_var_dp_4D( filename, id_var, d_grid_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, C%nz, 1 /) )

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid_with_time)
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_3D

  SUBROUTINE write_to_field_multiple_options_grid_dp_2D_notime( filename, grid, field_name_options, d)
    ! Write a 2-D data field defined on a grid to a NetCDF file variable on an x/y-grid
    !
    ! The variable in the NetCDF file has no time dimension.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:    ),        INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_2D_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:   ), POINTER               :: d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_2D( filename, var_name, should_have_time = .FALSE.)

    ! Transpose the output data
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, d_grid, wd_grid)

    d_grid( :,grid%i1:grid%i2) = d( :,grid%i1:grid%i2)

    CALL permute_2D_dp( d_grid, wd_grid , map = [2,1])

    ! Write data to the variable
    CALL write_var_dp_2D( filename, id_var, d_grid)

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_2D_notime

  SUBROUTINE write_to_field_multiple_options_grid_dp_2D_monthly_notime( filename, grid, field_name_options, d)
    ! Write a 2-D monthly data field defined on a grid to a NetCDF file variable on an x/y-grid
    !
    ! The variable in the NetCDF file has no time dimension.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:  ),        INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_2D_monthly_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:,:  ), POINTER              :: d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_2D_monthly( filename, var_name, should_have_time = .FALSE.)

    ! Transpose the output data
    CALL allocate_shared_dp_3D( 12, grid%ny, grid%nx, d_grid, wd_grid)

    d_grid( :, :, grid%i1:grid%i2) = d( :, :, grid%i1:grid%i2)

    CALL permute_3D_dp( d_grid, wd_grid , map = [3,2,1])

    ! Write data to the variable
    CALL write_var_dp_3D( filename, id_var, d_grid)

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_2D_monthly_notime

  SUBROUTINE write_to_field_multiple_options_grid_dp_3D_notime( filename, grid, field_name_options, d)
    ! Write a 3-D data field defined on a grid to a NetCDF file variable on an x/y-grid
    !
    ! The variable in the NetCDF file has no time dimension.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:  ),        INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_3D_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:,: ), POINTER               :: d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_3D( filename, var_name, should_have_time = .FALSE.)

    ! Transpose the output data
    CALL allocate_shared_dp_3D( C%nz, grid%ny, grid%nx, d_grid, wd_grid)

    d_grid( :, :,grid%i1:grid%i2) = d( :, :, grid%i1:grid%i2)
    CALL permute_3D_dp( d_grid, wd_grid , map = [3,2,1])

    ! Write data to the variable
    CALL write_var_dp_3D( filename, id_var, d_grid)

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_3D_notime

  SUBROUTINE write_to_field_multiple_options_grid_dp_3D_ocean_notime( filename, grid, field_name_options, d)
    ! Write a 3-D ocean data field defined on a grid to a NetCDF file variable on an x/y-grid
    !
    ! The variable in the NetCDF file has no time dimension.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:  ),        INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_field_multiple_options_grid_dp_3D_notime'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    REAL(dp), DIMENSION(:,:,: ), POINTER               :: d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire the variable
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('no variables for name options "' // TRIM( field_name_options) // '" were found in file "' // TRIM( filename) // '"!')

    ! Check if this variable has the correct type and dimensions
    CALL check_xy_grid_field_dp_3D_ocean( filename, var_name, should_have_time = .FALSE.)

    ! Transpose the output data
    CALL allocate_shared_dp_3D( C%nz_ocean, grid%ny, grid%nx, d_grid, wd_grid)

    d_grid( :, :,grid%i1:grid%i2) = d( :, :, grid%i1:grid%i2)
    CALL permute_3D_dp( d_grid, wd_grid , map = [3,2,1])

    ! Write data to the variable
    CALL write_var_dp_3D( filename, id_var, d_grid)

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_field_multiple_options_grid_dp_3D_ocean_notime

  ! Write new time value to file
  SUBROUTINE write_time_to_file( filename, time)
    ! Write new time value to file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_time_to_file'
    INTEGER                                            :: id_dim_time
    INTEGER                                            :: id_var_time
    INTEGER                                            :: nt

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check time dimension and variable
    CALL check_time( filename)

    ! Determine current length of time dimension in file
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = nt)

    ! Inquire variable id
    CALL inquire_var_multiple_options( filename, field_name_options_time, id_var_time)

    ! Write time
    nt = nt + 1
    CALL write_var_dp_1D( filename, id_var_time, (/ time /), start = (/ nt /), count = (/ 1 /) )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_time_to_file

  ! ===== Set up grid, time, and field variables in a NetCDF file =====
  ! ========================================================================

  ! Set up x/y-grid and gridded variables
  SUBROUTINE setup_xy_grid_in_netcdf_file( filename, grid)
    ! Set up a regular x/y-grid in an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid),                     INTENT(IN)    :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_xy_grid_in_netcdf_file'
    INTEGER                                            :: id_dim_x
    INTEGER                                            :: id_dim_y
    INTEGER                                            :: id_var_x
    INTEGER                                            :: id_var_y
    INTEGER                                            :: id_var_lon
    INTEGER                                            :: id_var_lat
    REAL(dp), DIMENSION(:,:), POINTER                  ::  d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create x/y dimensions
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_x), grid%nx, id_dim_x)
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_y), grid%ny, id_dim_y)

    ! Create and write x/y variables

    ! x
    CALL create_variable( filename, get_first_option_from_list( field_name_options_x), NF90_DOUBLE, (/ id_dim_x /), id_var_x)
    CALL add_attribute_char( filename, id_var_x, 'long_name', 'x-coordinate')
    CALL add_attribute_char( filename, id_var_x, 'units'    , 'm'           )
    CALL write_var_dp_1D( filename, id_var_x, grid%x)
    ! y
    CALL create_variable( filename, get_first_option_from_list( field_name_options_y), NF90_DOUBLE, (/ id_dim_y /), id_var_y)
    CALL add_attribute_char( filename, id_var_y, 'long_name', 'y-coordinate')
    CALL add_attribute_char( filename, id_var_y, 'units'    , 'm'           )
    CALL write_var_dp_1D( filename, id_var_y, grid%y)

    ! Create and write lon/lat variables

    ! lon
    CALL add_field_grid_dp_2D_notime( filename, get_first_option_from_list( field_name_options_lon), long_name = 'Longitude', units = 'degrees east')
    CALL inquire_var_multiple_options( filename, field_name_options_lon, id_var_lon)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, d_grid, wd_grid)
    d_grid( :, grid%i1:grid%i2) = grid%lon( :,grid%i1:grid%i2)
    CALL sync
    CALL permute_2D_dp( d_grid, wd_grid, map = [2,1])
    CALL write_var_dp_2D( filename, id_var_lon, d_grid)
    CALL deallocate_shared(wd_grid)

    ! lat
    CALL add_field_grid_dp_2D_notime( filename, get_first_option_from_list( field_name_options_lat), long_name = 'Latitude', units = 'degrees north')
    CALL inquire_var_multiple_options( filename, field_name_options_lat, id_var_lat)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, d_grid, wd_grid)
    d_grid( :, grid%i1:grid%i2) = grid%lat( :,grid%i1:grid%i2)
    CALL sync
    CALL permute_2D_dp( d_grid, wd_grid, map = [2,1])
    CALL write_var_dp_2D( filename, id_var_lat, d_grid)
    CALL deallocate_shared(wd_grid)

    ! Clean up after yourself
    CALL deallocate_shared( wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_xy_grid_in_netcdf_file

  SUBROUTINE add_field_grid_int_2D( filename, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_int_2D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(    filename)
    CALL check_y(    filename)
    CALL check_time( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multiple_options( filename, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_INT, (/ id_dim_x, id_dim_y, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_int_2D( filename, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_int_2D

  SUBROUTINE add_field_grid_dp_2D( filename, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_2D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(    filename)
    CALL check_y(    filename)
    CALL check_time( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multiple_options( filename, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_2D( filename, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_2D

  SUBROUTINE add_field_grid_dp_2D_monthly( filename, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_2D_monthly'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_month, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(     filename)
    CALL check_y(     filename)
    CALL check_month( filename)
    CALL check_time(  filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x    , id_dim_x    )
    CALL inquire_dim_multiple_options( filename, field_name_options_y    , id_dim_y    )
    CALL inquire_dim_multiple_options( filename, field_name_options_month, id_dim_month)
    CALL inquire_dim_multiple_options( filename, field_name_options_time , id_dim_time )

    ! Safety
    IF (id_dim_x     == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y     == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_month == -1) CALL crash('no month dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time  == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_month, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_2D_monthly( filename, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_2D_monthly

  SUBROUTINE add_field_grid_dp_3D( filename, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_3D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_zeta, id_dim_time, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(    filename)
    CALL check_y(    filename)
    CALL check_zeta( filename)
    CALL check_time( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multiple_options( filename, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multiple_options( filename, field_name_options_zeta, id_dim_zeta)
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_time == -1) CALL crash('no time dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_zeta, id_dim_time /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_3D( filename, var_name, should_have_time = .TRUE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_3D

  SUBROUTINE add_field_grid_int_2D_notime( filename, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_int_2D_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y dimensions and variables are there
    CALL check_x( filename)
    CALL check_y( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x, id_dim_x)
    CALL inquire_dim_multiple_options( filename, field_name_options_y, id_dim_y)

    ! Safety
    IF (id_dim_x == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_INT, (/ id_dim_x, id_dim_y /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_int_2D( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_int_2D_notime

  SUBROUTINE add_field_grid_dp_2D_notime( filename, var_name, long_name, units)
    ! Add a 2-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_2D_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x( filename)
    CALL check_y( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x, id_dim_x)
    CALL inquire_dim_multiple_options( filename, field_name_options_y, id_dim_y)

    ! Safety
    IF (id_dim_x == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_2D( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_2D_notime

  SUBROUTINE add_field_grid_dp_2D_monthly_notime( filename, var_name, long_name, units)
    ! Add a 2-D monthly variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_2D_monthly_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_month, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(     filename)
    CALL check_y(     filename)
    CALL check_month( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x    , id_dim_x    )
    CALL inquire_dim_multiple_options( filename, field_name_options_y    , id_dim_y    )
    CALL inquire_dim_multiple_options( filename, field_name_options_month, id_dim_month)

    ! Safety
    IF (id_dim_x     == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y     == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_month == -1) CALL crash('no month dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_month /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_2D_monthly( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_2D_monthly_notime

  SUBROUTINE add_field_grid_dp_3D_notime( filename, var_name, long_name, units)
    ! Add a 3-D variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_3D_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_zeta, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(    filename)
    CALL check_y(    filename)
    CALL check_zeta( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multiple_options( filename, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multiple_options( filename, field_name_options_zeta, id_dim_zeta)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_zeta /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_3D( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_3D_notime

  SUBROUTINE add_field_grid_dp_3D_ocean_notime( filename, var_name, long_name, units)
    ! Add a 3-D ocean variable to an existing NetCDF file with an x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: long_name
    CHARACTER(LEN=*),          OPTIONAL, INTENT(IN)    :: units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_field_grid_dp_3D_ocean_notime'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_zeta, id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if x,y, and time dimensions and variables are there
    CALL check_x(       filename)
    CALL check_y(       filename)
    CALL check_z_ocean( filename)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x      , id_dim_x   )
    CALL inquire_dim_multiple_options( filename, field_name_options_y      , id_dim_y   )
    CALL inquire_dim_multiple_options( filename, field_name_options_z_ocean, id_dim_zeta)

    ! Safety
    IF (id_dim_x    == -1) CALL crash('no x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_y    == -1) CALL crash('no y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (id_dim_zeta == -1) CALL crash('no zeta dimension could be found in file "' // TRIM( filename) // '"!')

    ! Create variable
    CALL create_variable( filename, var_name, NF90_DOUBLE, (/ id_dim_x, id_dim_y, id_dim_zeta /), id_var)

    ! Add attributes
    IF (PRESENT( long_name)) CALL add_attribute_char( filename, id_var, 'long_name', long_name)
    IF (PRESENT( units    )) CALL add_attribute_char( filename, id_var, 'units'    , units    )

    ! Final safety check
    CALL check_xy_grid_field_dp_3D_ocean( filename, var_name, should_have_time = .FALSE.)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_field_grid_dp_3D_ocean_notime

  ! Add extra dimensions
  SUBROUTINE add_time_dimension_to_file( filename)
    ! Add a time dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_time_dimension_to_file'
    INTEGER                                            :: id_dim_time
    INTEGER                                            :: id_var_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create time dimension
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_time), NF90_UNLIMITED, id_dim_time)

    ! Create time variable
    CALL create_variable(  filename, get_first_option_from_list( field_name_options_time), NF90_DOUBLE, (/ id_dim_time /), id_var_time)
    CALL add_attribute_char( filename, id_var_time, 'long_name', 'Time')
    CALL add_attribute_char( filename, id_var_time, 'units', 'years')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_time_dimension_to_file

  SUBROUTINE add_time_history_dimension_to_file( filename, var_name_time_history, time_history)
    ! Add a time_history dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=256),                  INTENT(IN)    :: var_name_time_history
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: time_history

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_time_history_dimension_to_file'
    INTEGER                                            :: id_dim_time_history
    INTEGER                                            :: id_var_time_history

    INTEGER                                            :: i
    INTEGER                                            :: ntime_history

    ! Number of time history
    ntime_history = SIZE(time_history,1)

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create time_history dimension
    CALL create_dimension( filename, var_name_time_history, ntime_history, id_dim_time_history)

    ! Create time_history variable
    CALL create_variable(    filename, var_name_time_history, NF90_DOUBLE, (/ id_dim_time_history /), id_var_time_history)
    CALL add_attribute_char( filename, id_var_time_history, 'long_name', 'time_history from current model time-step')
    CALL add_attribute_char( filename, id_var_time_history, 'units', 'years')

    ! Write time_history variable
    CALL write_var_dp_1D( filename, id_var_time_history, time_history)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_time_history_dimension_to_file

  SUBROUTINE add_month_dimension_to_file( filename)
    ! Add a month dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_month_dimension_to_file'
    INTEGER                                            :: id_dim_month
    INTEGER                                            :: id_var_month

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create month dimension
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_month), 12, id_dim_month)

    ! Create month variable
    CALL create_variable(  filename, get_first_option_from_list( field_name_options_month), NF90_INT, (/ id_dim_month /), id_var_month)
    CALL add_attribute_char( filename, id_var_month, 'long_name', 'Month')
    CALL add_attribute_char( filename, id_var_month, 'units', '1-12')
    CALL add_attribute_char( filename, id_var_month, 'description', '1 = Jan, 2 = Feb, ..., 12 = Dec')

    ! Write month variable
    CALL write_var_int_1D( filename, id_var_month, (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 /) )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_month_dimension_to_file

  SUBROUTINE add_zeta_dimension_to_file( filename)
    ! Add a zeta dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_zeta_dimension_to_file'
    INTEGER                                            :: id_dim_zeta
    INTEGER                                            :: id_var_zeta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create month dimension
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_zeta), C%nz, id_dim_zeta)

    ! Create month variable
    CALL create_variable(  filename, get_first_option_from_list( field_name_options_zeta), NF90_DOUBLE, (/ id_dim_zeta /), id_var_zeta)
    CALL add_attribute_char( filename, id_var_zeta, 'long_name', 'Scaled vertical coordinate')
    CALL add_attribute_char( filename, id_var_zeta, 'units', '0-1')
    CALL add_attribute_char( filename, id_var_zeta, 'transformation', 'zeta = (h - z) / H; zeta = 0 at the ice surface; zeta = 1 at the ice base')

    ! Write month variable
    CALL write_var_dp_1D( filename, id_var_zeta, C%zeta)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_zeta_dimension_to_file

  SUBROUTINE add_ocean_dimension_to_file( filename)
    ! Add a zeta dimension and variable to an existing NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_ocean_dimension_to_file'
    INTEGER                                            :: id_dim_zeta
    INTEGER                                            :: id_var_zeta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create ocean dimension
    CALL create_dimension( filename, get_first_option_from_list( field_name_options_z_ocean), C%nz_ocean, id_dim_zeta)

    ! Create ocean variable
    CALL create_variable(  filename, get_first_option_from_list( field_name_options_z_ocean), NF90_DOUBLE, (/ id_dim_zeta /), id_var_zeta)
    CALL add_attribute_char( filename, id_var_zeta, 'long_name', 'Depth in ocean')
    CALL add_attribute_char( filename, id_var_zeta, 'units', 'm')

    ! Write ocean variable
    CALL write_var_dp_1D( filename, id_var_zeta, C%z_ocean)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_ocean_dimension_to_file

END MODULE netcdf_output_module
