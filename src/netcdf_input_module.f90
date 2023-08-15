MODULE netcdf_input_module

! ===== Flexible reading of input files =====
! ===========================================
!
! These routines allow for flexible reading of input files. The three top-level functions
! allow you to read 2-D, 2-D monthly, and 3-D data fields from NetCDF files. The files can
! contain the data on a regular x/y-grid, a regular lon/lat-grid. The routines
! will automatically detect which one of these it is, and select the appropriate
! subroutine. The data will also be automatically mapped to the provided model grid.

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
                                             allocate_shared_dp_4D, deallocate_shared
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             map_square_to_square_cons_2nd_order_2D, map_square_to_square_cons_2nd_order_3D, &
                                             transpose_dp_2D, transpose_dp_3D
  ! Import specific functionality
  USE data_types_module,               ONLY: type_grid, type_grid_lonlat, type_model_region
  USE netcdf,                          ONLY: NF90_MAX_VAR_DIMS, nf90_inq_varid, nf90_noerr

  USE utilities_module,                ONLY: flip_1D_dp, flip_2D_x1_dp, flip_2D_x2_dp, flip_3D_x1_dp, flip_3D_x2_dp, flip_3D_x3_dp, &
                                             permute_2D_dp, permute_3D_dp, permute_2D_int, permute_3D_int, inverse_oblique_sg_projection, &
                                             deallocate_grid, deallocate_grid_lonlat, remap_zeta_grid_dp, map_glob_to_grid_2D, map_glob_to_grid_3D
  USE netcdf_basic_module,             ONLY: nerr, field_name_options_x, field_name_options_y, field_name_options_zeta, field_name_options_z_ocean, &
                                             field_name_options_lon, field_name_options_lat, field_name_options_time, field_name_options_month, &
                                             field_name_options_Hi, field_name_options_Hb, field_name_options_Hs, field_name_options_dHb, &
                                             field_name_options_SL, field_name_options_Ti, field_name_options_u_surf, &
                                             field_name_options_T_ocean, field_name_options_S_ocean, field_name_options_v_surf, parse_field_name_options, &
                                             inquire_dim_multiple_options, inquire_var_multiple_options, &
                                             read_var_int_0D, read_var_int_1D, read_var_int_2D, read_var_int_3D, read_var_int_4D, &
                                             read_var_dp_0D , read_var_dp_1D , read_var_dp_2D , read_var_dp_3D , read_var_dp_4D, &
                                             check_x, check_y, check_lon, check_lat, check_zeta, check_z_ocean, find_timeframe, &
                                             check_xy_grid_field_int_2D, check_xy_grid_field_dp_2D, check_xy_grid_field_dp_2D_monthly, check_xy_grid_field_dp_3D, &
                                             check_lonlat_grid_field_int_2D, check_lonlat_grid_field_dp_2D, check_lonlat_grid_field_dp_ocean_3D, &
                                             check_lonlat_grid_field_dp_2D_monthly, check_lonlat_grid_field_dp_3D, &
                                             inquire_xy_grid, inquire_lonlat_grid, check_xy_grid_field_dp_3D_ocean

  IMPLICIT NONE

CONTAINS

  ! ===== Top-level functions =====
  ! ===============================

  ! Read data
  SUBROUTINE read_field_from_file_2D(         filename, field_name_options, grid, d, region_name, time_to_read)
    ! Read a data field from a NetCDF file and map it to the grid
    !
    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_file_2D'
    LOGICAL                                            :: file_exists
    LOGICAL                                            :: has_xy_grid, has_lonlat_grid
    TYPE(type_grid)                                    :: grid_from_file
    TYPE(type_grid_lonlat)                             :: grid_lonlat_from_file
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_grid_from_file
    INTEGER                                            :: wd_grid_from_file
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_grid_lonlat_from_file
    INTEGER                                            :: wd_grid_lonlat_from_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Find out on what kind of grid the file is defined
    CALL inquire_xy_grid(     filename, has_xy_grid    )
    CALL inquire_lonlat_grid( filename, has_lonlat_grid)

    ! Files with more than one grid are not recognised
    IF (has_xy_grid     .AND. has_lonlat_grid) CALL crash('file "' // TRIM( filename) // '" contains both an x/y-grid and a lon/lat-grid!')

    ! Choose the appropriate subroutine
    IF (has_xy_grid) THEN
      ! Data is provided on an x/y-grid

      ! Read grid and gridded data
      CALL read_field_from_xy_file_2D( filename, field_name_options, region_name, grid_from_file, d_grid_from_file, wd_grid_from_file, time_to_read)

      ! Map (transposed) raw data to the model grid
      CALL map_square_to_square_cons_2nd_order_2D( grid_from_file%nx, grid_from_file%ny, grid_from_file%x, grid_from_file%y, grid%nx, grid%ny, grid%x, grid%y, d_grid_from_file,  d )

      ! Clean up after yourself
      CALL deallocate_grid(                             grid_from_file)
      CALL deallocate_shared(                        wd_grid_from_file)

    ELSEIF (has_lonlat_grid) THEN
      ! Data is provided on a lon/lat-grid

      ! Read grid and gridded data
      CALL read_field_from_lonlat_file_2D( filename, field_name_options, region_name, grid_lonlat_from_file, d_grid_lonlat_from_file, wd_grid_lonlat_from_file, time_to_read)

      ! Map data to the model grid
      CALL map_glob_to_grid_2D( grid_lonlat_from_file%nlat, grid_lonlat_from_file%nlon, grid_lonlat_from_file%lat, grid_lonlat_from_file%lon, grid, d_grid_lonlat_from_file  , d     )

      ! Clean up after yourself
      CALL deallocate_grid_lonlat(                     grid_lonlat_from_file)
      CALL deallocate_shared(                       wd_grid_lonlat_from_file)

    ELSE
      CALL crash('file "' // TRIM( filename) // '" does not contain a recognised x/y-grid or lon/lat-grid!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_file_2D

  SUBROUTINE read_field_from_file_2D_monthly( filename, field_name_options, grid, d, region_name, time_to_read)
    ! Read a data field from a NetCDF file, and map it to the model grid.
    !
    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_grid),                     INTENT(IN) :: grid
    REAL(dp), DIMENSION(:,:,:  ),        INTENT(OUT)   :: d
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_file_2D_monthly'
    LOGICAL                                            :: file_exists
    LOGICAL                                            :: has_xy_grid, has_lonlat_grid
    TYPE(type_grid)                                    :: grid_from_file
    TYPE(type_grid_lonlat)                             :: grid_lonlat_from_file
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid_from_file
    INTEGER                                            :: wd_grid_from_file
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid_lonlat_from_file
    INTEGER                                            :: wd_grid_lonlat_from_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Find out on what kind of grid the file is defined
    CALL inquire_xy_grid(     filename, has_xy_grid    )
    CALL inquire_lonlat_grid( filename, has_lonlat_grid)

    ! Files with more than one grid are not recognised
    IF (has_xy_grid     .AND. has_lonlat_grid) CALL crash('file "' // TRIM( filename) // '" contains both an x/y-grid and a lon/lat-grid!')

    ! Choose the appropriate subroutine
    IF (has_xy_grid) THEN
      ! Data is provided on an x/y-grid

      ! Read grid and gridded data
      CALL read_field_from_xy_file_2D_monthly( filename, field_name_options, region_name, grid_from_file, d_grid_from_file, wd_grid_from_file, time_to_read)

      ! Map (transposed) raw data to the model grid
      CALL map_square_to_square_cons_2nd_order_3D( grid_from_file%nx, grid_from_file%ny, grid_from_file%x, grid_from_file%y, grid%nx, grid%ny, grid%x, grid%y, d_grid_from_file,  d )

      ! Clean up after yourself
      CALL deallocate_grid(                             grid_from_file)
      CALL deallocate_shared(                        wd_grid_from_file)

    ELSEIF (has_lonlat_grid) THEN
      ! Data is provided on a lon/lat-grid

      ! Read grid and gridded data
      CALL read_field_from_lonlat_file_2D_monthly( filename, field_name_options, region_name, grid_lonlat_from_file, d_grid_lonlat_from_file, wd_grid_lonlat_from_file, time_to_read)

      ! Map data to the model grid
      CALL map_glob_to_grid_3D( grid_lonlat_from_file%nlat, grid_lonlat_from_file%nlon, grid_lonlat_from_file%lat, grid_lonlat_from_file%lon, grid, d_grid_lonlat_from_file     , d     )

      ! Clean up after yourself
      CALL deallocate_grid_lonlat(                     grid_lonlat_from_file)
      CALL deallocate_shared(                       wd_grid_lonlat_from_file)

    ELSE
      CALL crash('file "' // TRIM( filename) // '" does not contain a recognised x/y-grid, lon/lat-grid!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_file_2D_monthly

  SUBROUTINE read_field_from_file_3D(         filename, field_name_options, grid, d, region_name, time_to_read)
    ! Read a data field from a NetCDF file, and map it to the model grid.
    !
    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.
    !
    ! NOTE: only meant to be used for 3-D englacial data, not for 3-D ocean data!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:  ),        INTENT(OUT)   :: d
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_file_3D'
    LOGICAL                                            :: file_exists
    LOGICAL                                            :: has_xy_grid, has_lonlat_grid
    TYPE(type_grid)                                    :: grid_from_file
    TYPE(type_grid_lonlat)                             :: grid_lonlat_from_file
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid_from_file
    INTEGER                                            :: wd_grid_from_file
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid_lonlat_from_file
    INTEGER                                            :: wd_grid_lonlat_from_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Find out on what kind of grid the file is defined
    CALL inquire_xy_grid(     filename, has_xy_grid    )
    CALL inquire_lonlat_grid( filename, has_lonlat_grid)

    ! Files with more than one grid are not recognised
    IF (has_xy_grid     .AND. has_lonlat_grid) CALL crash('file "' // TRIM( filename) // '" contains both an x/y-grid and a lon/lat-grid!')

    ! Choose the appropriate subroutine
    IF (has_xy_grid) THEN
      ! Data is provided on an x/y-grid

      ! Read grid and gridded data
      CALL read_field_from_xy_file_3D( filename, field_name_options, region_name, grid_from_file, d_grid_from_file, wd_grid_from_file, time_to_read)

      ! Map (transposed) raw data to the model grid
      CALL map_square_to_square_cons_2nd_order_3D( grid_from_file%nx, grid_from_file%ny, grid_from_file%x, grid_from_file%y, grid%nx, grid%ny, grid%x, grid%y, d_grid_from_file,  d )

      ! Clean up after yourself
      CALL deallocate_grid(                             grid_from_file)
      CALL deallocate_shared(                        wd_grid_from_file)

    ELSEIF (has_lonlat_grid) THEN
      ! Data is provided on a lon/lat-grid

      ! Read grid and gridded data
      CALL read_field_from_lonlat_file_3D( filename, field_name_options, region_name, grid_lonlat_from_file, d_grid_lonlat_from_file, wd_grid_lonlat_from_file, time_to_read)

      ! Map data to the model grid
      CALL map_glob_to_grid_3D( grid_lonlat_from_file%nlat, grid_lonlat_from_file%nlon, grid_lonlat_from_file%lat, grid_lonlat_from_file%lon, grid, d_grid_lonlat_from_file     , d     )

      ! Clean up after yourself
      CALL deallocate_grid_lonlat(                     grid_lonlat_from_file)
      CALL deallocate_shared(                       wd_grid_lonlat_from_file)

    ELSE
      CALL crash('file "' // TRIM( filename) // '" does not contain a recognised x/y-grid or lon/lat-grid!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_file_3D

  SUBROUTINE read_field_from_file_ocean_3D(         filename, field_name_options, grid, d, region_name, time_to_read)
    ! Read a data field from a NetCDF file, and map it to the model grid.
    !
    ! Ultimate flexibility; the file can provide the data on a global lon/lat-grid,
    ! a regional x/y-grid - it matters not, all shall be fine.
    ! The order of dimensions ([x,y] or [y,x], [lon,lat] or [lat,lon]) and direction
    ! (increasing or decreasing) also does not matter any more.
    !
    ! NOTE: only meant to be used for 3-D englacial data, not for 3-D ocean data!

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:  ),        INTENT(OUT)   :: d
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_file_3D'
    LOGICAL                                            :: file_exists
    LOGICAL                                            :: has_xy_grid, has_lonlat_grid
    TYPE(type_grid)                                    :: grid_from_file
    TYPE(type_grid_lonlat)                             :: grid_lonlat_from_file
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid_from_file
    INTEGER                                            :: wd_grid_from_file
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_grid_lonlat_from_file
    INTEGER                                            :: wd_grid_lonlat_from_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Find out on what kind of grid the file is defined
    CALL inquire_xy_grid(     filename, has_xy_grid    )
    CALL inquire_lonlat_grid( filename, has_lonlat_grid)

    ! Files with more than one grid are not recognised
    IF (has_xy_grid     .AND. has_lonlat_grid) CALL crash('file "' // TRIM( filename) // '" contains both an x/y-grid and a lon/lat-grid!')

    ! Choose the appropriate subroutine
    IF (has_xy_grid) THEN
      ! Data is provided on an x/y-grid

      ! Read grid and gridded data
      CALL read_field_from_xy_file_ocean_3D( filename, field_name_options, region_name, grid_from_file, d_grid_from_file, wd_grid_from_file, time_to_read)

      ! Map (transposed) raw data to the model grid
      CALL map_square_to_square_cons_2nd_order_3D( grid_from_file%nx, grid_from_file%ny, grid_from_file%x, grid_from_file%y, grid%nx, grid%ny, grid%x, grid%y, d_grid_from_file,  d )

      ! Clean up after yourself
      CALL deallocate_grid(                             grid_from_file)
      CALL deallocate_shared(                        wd_grid_from_file)

    ELSEIF (has_lonlat_grid) THEN
      ! Data is provided on a lon/lat-grid

      ! Read grid and gridded data
      CALL read_field_from_lonlat_file_ocean_3D( filename, field_name_options, region_name, grid_lonlat_from_file, d_grid_lonlat_from_file, wd_grid_lonlat_from_file, time_to_read)

      ! Map data to the model grid
      CALL map_glob_to_grid_3D( grid_lonlat_from_file%nlat, grid_lonlat_from_file%nlon, grid_lonlat_from_file%lat, grid_lonlat_from_file%lon, grid, d_grid_lonlat_from_file     , d     )

      ! Clean up after yourself
      CALL deallocate_grid_lonlat(                     grid_lonlat_from_file)
      CALL deallocate_shared(                       wd_grid_lonlat_from_file)

    ELSE
      CALL crash('file "' // TRIM( filename) // '" does not contain a recognised x/y-grid or lon/lat-grid!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_file_ocean_3D

  SUBROUTINE read_field_from_xy_file_ocean_3D(             filename, field_name_options, region_name, grid, d_z_ocean, wd_z_ocean, time_to_read)
    ! Read a 3-D ocean data field from a NetCDF file on an x/y-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid),                     INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: d_z_ocean
    INTEGER,                             INTENT(OUT)   :: wd_z_ocean
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_xy_file_ocean_3D'
    INTEGER                                            :: nz_ocean
    REAL(dp), DIMENSION(:    ), POINTER                :: z_ocean
    INTEGER                                            :: wz_ocean
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, xdir, ydir
    REAL(dp), DIMENSION(:,:,:,:), POINTER              :: d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set up the grid from the file
    CALL setup_xy_grid_from_file( filename, grid, region_name)

    ! Set up the _ocean coordinate from the file
    CALL setup_z_ocean_from_file( filename, nz_ocean, z_ocean, wz_ocean)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_xy_grid_field_dp_3D_ocean( filename, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_xy_indexing( filename, var_name, indexing, xdir, ydir)

    IF     (indexing == 'xy') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nx, grid%ny, nz_ocean, d_z_ocean, wd_z_ocean)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, id_var, d_z_ocean)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nx, grid%ny, nz_ocean, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, nz_ocean, 1 /) )
        ! Copy to output memory
        d_z_ocean( grid%i1:grid%i2,:,:) = d_with_time( grid%i1:grid%i2,:,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'yx') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%ny, grid%nx, nz_ocean, d_z_ocean, wd_z_ocean)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, id_var, d_z_ocean)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%ny, grid%nx, nz_ocean, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%ny, grid%nx, nz_ocean, 1 /) )
        ! Copy to output memory
        d_z_ocean( :,grid%i1:grid%i2,:) = d_with_time( :,grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'xy') THEN
      ! No need to do anything
    ELSEIF (indexing == 'yx') THEN
      CALL permute_3D_dp( d_z_ocean, wd_z_ocean, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! xdir
    IF     (xdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (xdir == 'reverse') THEN
      CALL flip_1D_dp( grid%x)
      CALL flip_3D_x1_dp( d_z_ocean)
    ELSE
      CALL crash('unknown xdir = "' // TRIM( xdir) // '"!')
    END IF

    ! ydir
    IF     (ydir == 'normal') THEN
      ! No need to do anything
    ELSEIF (ydir == 'reverse') THEN
      CALL flip_1D_dp( grid%y)
      CALL flip_3D_x2_dp( d_z_ocean)
    ELSE
      CALL crash('unknown ydir = "' // TRIM( ydir) // '"!')
    END IF

    ! Transpose the input data
    CALL transpose_dp_3D( d_z_ocean, wd_z_ocean)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_xy_file_ocean_3D

  ! ===== Medium-level functions =====
  ! ==================================

  ! Read a field from an x/y-grid file
  SUBROUTINE read_field_from_xy_file_2D(             filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 2-D data field from a NetCDF file on an x/y-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid),                     INTENT(OUT)   :: grid
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_xy_file_2D'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, xdir, ydir
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set up the grid from the file
    CALL setup_xy_grid_from_file( filename, grid, region_name)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_xy_grid_field_dp_2D( filename, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_xy_indexing( filename, var_name, indexing, xdir, ydir)

    IF     (indexing == 'xy') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_2D( grid%nx, grid%ny, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_2D( filename, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_3D( grid%nx, grid%ny, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_3D( filename, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ grid%nx, grid%ny, 1 /) )
        ! Copy to output memory
        d( grid%i1:grid%i2,:) = d_with_time( grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'yx') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_2D( filename, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_3D( grid%ny, grid%nx, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_3D( filename, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ grid%ny, grid%nx, 1 /) )
        ! Copy to output memory
        d( :,grid%i1:grid%i2) = d_with_time( :,grid%i1:grid%i2,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'xy') THEN
      ! No need to do anything
    ELSEIF (indexing == 'yx') THEN
      CALL permute_2D_dp( d, wd, map = [2,1])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! xdir
    IF     (xdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (xdir == 'reverse') THEN
      CALL flip_1D_dp( grid%x)
      CALL flip_2D_x1_dp( d)
    ELSE
      CALL crash('unknown xdir = "' // TRIM( xdir) // '"!')
    END IF

    ! ydir
    IF     (ydir == 'normal') THEN
      ! No need to do anything
    ELSEIF (ydir == 'reverse') THEN
      CALL flip_1D_dp( grid%y)
      CALL flip_2D_x2_dp( d)
    ELSE
      CALL crash('unknown ydir = "' // TRIM( ydir) // '"!')
    END IF

    ! Transpose the input data
    CALL transpose_dp_2D( d, wd)

    ! Finalise routine path
    CALL finalise_routine( routine_name, 19)

  END SUBROUTINE read_field_from_xy_file_2D

  SUBROUTINE read_field_from_xy_file_2D_monthly(     filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 2-D monthly data field from a NetCDF file on an x/y-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid),                     INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_xy_file_2D_monthly'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, xdir, ydir
    REAL(dp), DIMENSION(:,:,:,:), POINTER              ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set up the grid from the file
    CALL setup_xy_grid_from_file( filename, grid, region_name)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_xy_grid_field_dp_2D_monthly( filename, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_xy_indexing( filename, var_name, indexing, xdir, ydir)

    IF     (indexing == 'xy') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nx, grid%ny, 12, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nx, grid%ny, 12, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, 12, 1 /) )
        ! Copy to output memory
        d( grid%i1:grid%i2,:,:) = d_with_time( grid%i1:grid%i2,:,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'yx') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%ny, grid%nx, 12, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%ny, grid%nx, 12, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%ny, grid%nx, 12, 1 /) )
        ! Copy to output memory
        d( :,grid%i1:grid%i2,:) = d_with_time( :,grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'xy') THEN
      ! No need to do anything
    ELSEIF (indexing == 'yx') THEN
      CALL permute_3D_dp( d, wd, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! xdir
    IF     (xdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (xdir == 'reverse') THEN
      CALL flip_1D_dp( grid%x)
      CALL flip_3D_x1_dp( d)
    ELSE
      CALL crash('unknown xdir = "' // TRIM( xdir) // '"!')
    END IF

    ! ydir
    IF     (ydir == 'normal') THEN
      ! No need to do anything
    ELSEIF (ydir == 'reverse') THEN
      CALL flip_1D_dp( grid%y)
      CALL flip_3D_x2_dp( d)
    ELSE
      CALL crash('unknown ydir = "' // TRIM( ydir) // '"!')
    END IF

    ! Transpose the input data
    CALL transpose_dp_3D( d, wd )

    ! Finalise routine path
    CALL finalise_routine( routine_name, 19)

  END SUBROUTINE read_field_from_xy_file_2D_monthly

  SUBROUTINE read_field_from_xy_file_3D(             filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 3-D data field from a NetCDF file on an x/y-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid),                     INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_xy_file_3D'
    INTEGER                                            :: nzeta_from_file
    REAL(dp), DIMENSION(:    ), POINTER                :: zeta_from_file
    INTEGER                                            :: wzeta_from_file
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, xdir, ydir
    REAL(dp), DIMENSION(:,:,:,:), POINTER              :: d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:,:,:), POINTER                :: d_zeta_from_file
    INTEGER                                            :: wd_zeta_from_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set up the grid from the file
    CALL setup_xy_grid_from_file( filename, grid, region_name)

    ! Set up the zeta coordinate from the file
    CALL setup_zeta_from_file( filename, nzeta_from_file, zeta_from_file, wzeta_from_file)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_xy_grid_field_dp_3D( filename, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_xy_indexing( filename, var_name, indexing, xdir, ydir)

    IF     (indexing == 'xy') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nx, grid%ny, nzeta_from_file, d_zeta_from_file, wd_zeta_from_file)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, id_var, d_zeta_from_file)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nx, grid%ny, nzeta_from_file, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nx, grid%ny, nzeta_from_file, 1 /) )
        ! Copy to output memory
        d_zeta_from_file( grid%i1:grid%i2,:,:) = d_with_time( grid%i1:grid%i2,:,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'yx') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%ny, grid%nx, nzeta_from_file, d_zeta_from_file, wd_zeta_from_file)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, id_var, d_zeta_from_file)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%ny, grid%nx, nzeta_from_file, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%ny, grid%nx, nzeta_from_file, 1 /) )
        ! Copy to output memory
        d_zeta_from_file( :,grid%i1:grid%i2,:) = d_with_time( :,grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'xy') THEN
      ! No need to do anything
    ELSEIF (indexing == 'yx') THEN
      CALL permute_3D_dp( d_zeta_from_file, wd_zeta_from_file, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! xdir
    IF     (xdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (xdir == 'reverse') THEN
      CALL flip_1D_dp( grid%x)
      CALL flip_3D_x1_dp( d_zeta_from_file)
    ELSE
      CALL crash('unknown xdir = "' // TRIM( xdir) // '"!')
    END IF

    ! ydir
    IF     (ydir == 'normal') THEN
      ! No need to do anything
    ELSEIF (ydir == 'reverse') THEN
      CALL flip_1D_dp( grid%y)
      CALL flip_3D_x2_dp( d_zeta_from_file)
    ELSE
      CALL crash('unknown ydir = "' // TRIM( ydir) // '"!')
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, C%nz, d, wd)

    ! Remap to the model vertical grid
    CALL remap_zeta_grid_dp( zeta_from_file, d_zeta_from_file, C%zeta, d)

    ! Transpose the input data
    CALL transpose_dp_3D( d, wd )

    ! Clean up after yourself
    CALL deallocate_shared(   wzeta_from_file)
    CALL deallocate_shared( wd_zeta_from_file)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_xy_file_3D

  ! Read a field from a lon/lat-grid file
  SUBROUTINE read_field_from_lonlat_file_2D(         filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 2-D data field from a NetCDF file on a lon/lat-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_lonlat_file_2D'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, londir, latdir
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set up the grid from the file
    CALL setup_lonlat_grid_from_file( filename, grid, region_name)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_lonlat_grid_field_dp_2D( filename, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_lonlat_indexing( filename, var_name, indexing, londir, latdir)

    IF     (indexing == 'lonlat') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_2D( grid%nlon, grid%nlat, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN

        CALL read_var_dp_2D( filename, id_var, d)

      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_3D( grid%nlon, grid%nlat, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_3D( filename, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ grid%nlon, grid%nlat, 1 /) )
        ! Copy to output memory
        d( grid%i1:grid%i2,:) = d_with_time( grid%i1:grid%i2,:,1)

        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'latlon') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_2D( grid%nlat, grid%nlon, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_2D( filename, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_3D( grid%nlat, grid%nlon, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_3D( filename, id_var, d_with_time, start = (/ 1, 1, ti /), count = (/ grid%nlat, grid%nlon, 1 /) )
        ! Copy to output memory
        d( :,grid%i1:grid%i2) = d_with_time( :,grid%i1:grid%i2,1)

        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF
    CALL sync

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'lonlat') THEN
      ! No need to do anything
    ELSEIF (indexing == 'latlon') THEN
      CALL permute_2D_dp( d, wd, map = [2,1])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF
    CALL sync

    ! londir
    IF     (londir == 'normal') THEN
      ! No need to do anything
    ELSEIF (londir == 'reverse') THEN
      CALL flip_1D_dp( grid%lon)
      CALL flip_2D_x1_dp( d)
    ELSE
      CALL crash('unknown londir = "' // TRIM( londir) // '"!')
    END IF
    CALL sync

    ! latdir
    IF     (latdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (latdir == 'reverse') THEN
      CALL flip_1D_dp( grid%lat)
      CALL flip_2D_x2_dp( d)
    ELSE
      CALL crash('unknown latdir = "' // TRIM( latdir) // '"!')
    END IF
    CALL sync

    ! Correct longitude shifts and range
    CALL correct_longitude_shifts_and_range_2D( filename, grid, d)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_lonlat_file_2D

  SUBROUTINE read_field_from_lonlat_file_2D_monthly( filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 2-D monthly data field from a NetCDF file on a lon/lat-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid_lonlat),              INTENT(INOUT)  :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_lonlat_file_2D'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, londir, latdir
    REAL(dp), DIMENSION(:,:,:,:), POINTER              ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set up the grid from the file
    CALL setup_lonlat_grid_from_file( filename, grid, region_name)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_lonlat_grid_field_dp_2D_monthly( filename, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_lonlat_indexing( filename, var_name, indexing, londir, latdir)

    IF     (indexing == 'lonlat') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nlon, grid%nlat, 12, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nlon, grid%nlat, 12, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlon, grid%nlat, 12, 1 /) )
        ! Copy to output memory
        d( grid%i1:grid%i2,:,:) = d_with_time( grid%i1:grid%i2,:,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'latlon') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nlat, grid%nlon, 12, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nlat, grid%nlon, 12, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlat, grid%nlon, 12, 1 /) )
        ! Copy to output memory
        d( :,grid%i1:grid%i2,:) = d_with_time( :,grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'lonlat') THEN
      ! No need to do anything
    ELSEIF (indexing == 'latlon') THEN
      CALL permute_3D_dp( d, wd, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! londir
    IF     (londir == 'normal') THEN
      ! No need to do anything
    ELSEIF (londir == 'reverse') THEN
      CALL flip_1D_dp( grid%lon)
      CALL flip_3D_x1_dp( d)
    ELSE
      CALL crash('unknown londir = "' // TRIM( londir) // '"!')
    END IF

    ! latdir
    IF     (latdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (latdir == 'reverse') THEN
      CALL flip_1D_dp( grid%lat)
      CALL flip_3D_x2_dp( d)
    ELSE
      CALL crash('unknown latdir = "' // TRIM( latdir) // '"!')
    END IF

    ! Correct longitude shifts and range
    CALL correct_longitude_shifts_and_range_3D( filename, grid, d)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_lonlat_file_2D_monthly

  SUBROUTINE read_field_from_lonlat_file_3D(         filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 3-D data field from a NetCDF file on a lon/lat-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_lonlat_file_3D'
    INTEGER                                            :: nzeta_from_file
    REAL(dp), DIMENSION(:    ), POINTER                ::  zeta_from_file
    INTEGER                                            :: wzeta_from_file
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, londir, latdir
    REAL(dp), DIMENSION(:,:,:,:), POINTER              ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  d_zeta_from_file
    INTEGER                                            :: wd_zeta_from_file

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set up the grid from the file
    CALL setup_lonlat_grid_from_file( filename, grid, region_name)

    ! Set up the zeta coordinate from the file
    CALL setup_zeta_from_file( filename, nzeta_from_file, zeta_from_file, wzeta_from_file)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_lonlat_grid_field_dp_3D( filename, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_lonlat_indexing( filename, var_name, indexing, londir, latdir)

    IF     (indexing == 'lonlat') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nlon, grid%nlat, nzeta_from_file, d_zeta_from_file, wd_zeta_from_file)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, id_var, d_zeta_from_file)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nlon, grid%nlat, nzeta_from_file, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlon, grid%nlat, nzeta_from_file, 1 /) )
        ! Copy to output memory
        d_zeta_from_file( grid%i1:grid%i2,:,:) = d_with_time( grid%i1:grid%i2,:,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'latlon') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nlat, grid%nlon, nzeta_from_file, d_zeta_from_file, wd_zeta_from_file)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, id_var, d_zeta_from_file)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nlat, grid%nlon, nzeta_from_file, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlat, grid%nlon, nzeta_from_file, 1 /) )
        ! Copy to output memory
        d_zeta_from_file( :,grid%i1:grid%i2,:) = d_with_time( :,grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'lonlat') THEN
      ! No need to do anything
    ELSEIF (indexing == 'latlon') THEN
      CALL permute_3D_dp( d_zeta_from_file, wd_zeta_from_file, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! londir
    IF     (londir == 'normal') THEN
      ! No need to do anything
    ELSEIF (londir == 'reverse') THEN
      CALL flip_1D_dp( grid%lon)
      CALL flip_3D_x1_dp( d_zeta_from_file)
    ELSE
      CALL crash('unknown londir = "' // TRIM( londir) // '"!')
    END IF

    ! latdir
    IF     (latdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (latdir == 'reverse') THEN
      CALL flip_1D_dp( grid%lat)
      CALL flip_3D_x2_dp( d_zeta_from_file)
    ELSE
      CALL crash('unknown latdir = "' // TRIM( latdir) // '"!')
    END IF

    ! Correct longitude shifts and range
    CALL correct_longitude_shifts_and_range_3D( filename, grid, d_zeta_from_file)

    ! Allocate shared memory
    CALL allocate_shared_dp_3D( grid%nlon, grid%nlat, C%nz, d, wd)

    ! Remap to the model vertical grid
    CALL remap_zeta_grid_dp( zeta_from_file, d_zeta_from_file, C%zeta, d)

    ! Clean up after yourself
    CALL deallocate_shared(   wzeta_from_file)
    CALL deallocate_shared( wd_zeta_from_file)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_lonlat_file_3D

  SUBROUTINE read_field_from_lonlat_file_ocean_3D(         filename, field_name_options, region_name, grid, d, wd, time_to_read)
    ! Read a 3-D data field from a NetCDF file on a lon/lat-grid,
    ! and return both the grid and the data.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: d
    INTEGER,                             INTENT(OUT)   :: wd
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_field_from_lonlat_file_ocean_3D'
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    CHARACTER(LEN=256)                                 :: indexing, londir, latdir
    REAL(dp), DIMENSION(:,:,:,:), POINTER              ::  d_with_time
    INTEGER                                            :: wd_with_time
    INTEGER                                            :: ti
    INTEGER                                            :: nz_ocean
    REAL(dp), DIMENSION(:    ),   POINTER              ::  z_ocean
    INTEGER                                            :: wz_ocean

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set up the grid from the file
    CALL setup_lonlat_grid_from_file( filename, grid, region_name)

    ! Set up the ocean coordinate from the file
    CALL setup_z_ocean_from_file( filename, nz_ocean, z_ocean, wz_ocean)

    ! Look for the specified variable in the file
    CALL inquire_var_multiple_options( filename, field_name_options, id_var, var_name = var_name)
    IF (id_var == -1) CALL crash('couldnt find any of the options "' // TRIM( field_name_options) // '" in file "' // TRIM( filename)  // '"!')

    ! Check if the variable has the required dimensions
    CALL check_lonlat_grid_field_dp_ocean_3D( filename, var_name, should_have_time = PRESENT( time_to_read))

    ! Determine indexing and dimension directions
    CALL determine_lonlat_indexing( filename, var_name, indexing, londir, latdir)

    IF     (indexing == 'lonlat') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nlon, grid%nlat, nz_ocean, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nlon, grid%nlat, nz_ocean, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlon, grid%nlat, nz_ocean, 1 /) )
        ! Copy to output memory
        d( grid%i1:grid%i2,:,:) = d_with_time( grid%i1:grid%i2,:,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSEIF (indexing == 'latlon') THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_3D( grid%nlat, grid%nlon, nz_ocean, d, wd)

      ! Read data from file
      IF (.NOT. PRESENT( time_to_read)) THEN
        CALL read_var_dp_3D( filename, id_var, d)
      ELSE
        ! Allocate shared memory
        CALL allocate_shared_dp_4D( grid%nlat, grid%nlon, nz_ocean, 1, d_with_time, wd_with_time)
        ! Find out which timeframe to read
        CALL find_timeframe( filename, time_to_read, ti)
        ! Read data
        CALL read_var_dp_4D( filename, id_var, d_with_time, start = (/ 1, 1, 1, ti /), count = (/ grid%nlat, grid%nlon, nz_ocean, 1 /) )
        ! Copy to output memory
        d( :,grid%i1:grid%i2,:) = d_with_time( :,grid%i1:grid%i2,:,1)
        ! Clean up after yourself
        CALL deallocate_shared( wd_with_time)
      END IF

    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! Perform necessary corrections to the gridded data

    ! Indexing
    IF     (indexing == 'lonlat') THEN
      ! No need to do anything
    ELSEIF (indexing == 'latlon') THEN
      CALL permute_3D_dp( d, wd, map = [2,1,3])
    ELSE
      CALL crash('unknown indexing = "' // TRIM( indexing) // '"!')
    END IF

    ! londir
    IF     (londir == 'normal') THEN
      ! No need to do anything
    ELSEIF (londir == 'reverse') THEN
      CALL flip_1D_dp( grid%lon)
      CALL flip_3D_x1_dp( d)
    ELSE
      CALL crash('unknown londir = "' // TRIM( londir) // '"!')
    END IF

    ! latdir
    IF     (latdir == 'normal') THEN
      ! No need to do anything
    ELSEIF (latdir == 'reverse') THEN
      CALL flip_1D_dp( grid%lat)
      CALL flip_3D_x2_dp( d)
    ELSE
      CALL crash('unknown latdir = "' // TRIM( latdir) // '"!')
    END IF

    ! Correct longitude shifts and range
    CALL correct_longitude_shifts_and_range_3D( filename, grid, d)

    CALL deallocate_shared( wz_ocean)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_field_from_lonlat_file_ocean_3D

  ! ===== Set up grids from a NetCDF file =====
  ! ================================================

  SUBROUTINE setup_xy_grid_from_file(     filename, grid, region_name)
    ! Set up an x/y-grid from a NetCDF file
    !
    ! Assumes no memory has yet been allocated for the grid at all

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid),                     INTENT(INOUT) :: grid
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_xy_grid_from_file'
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp
    INTEGER                                            :: id_dim_x, id_dim_y
    INTEGER                                            :: id_var_x, id_var_y
    INTEGER                                            :: i,j,n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check grid dimensions and variables for validity
    CALL check_x( filename)
    CALL check_y( filename)

    ! Allocate memory for the grid size
    CALL allocate_shared_int_0D( grid%nx, grid%wnx)
    CALL allocate_shared_int_0D( grid%ny, grid%wny)
    CALL allocate_shared_int_0D( grid%n , grid%wn )

    ! Inquire x and y dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x, id_dim_x, dim_length = grid%nx)
    CALL inquire_dim_multiple_options( filename, field_name_options_y, id_dim_y, dim_length = grid%ny)
    grid%n = grid%nx * grid%ny

    ! Allocate memory for x and y
    CALL allocate_shared_dp_1D( grid%nx, grid%x, grid%wx)
    CALL allocate_shared_dp_1D( grid%ny, grid%y, grid%wy)

    ! Inquire x and y variables
    CALL inquire_var_multiple_options( filename, field_name_options_x, id_var_x)
    CALL inquire_var_multiple_options( filename, field_name_options_y, id_var_y)

    ! Read x and y
    CALL read_var_dp_1D(  filename, id_var_x, grid%x)
    CALL read_var_dp_1D(  filename, id_var_y, grid%y)

    ! Allocate memory for, and calculate, some secondary grid properties
    CALL allocate_shared_dp_0D(                    grid%xmin       , grid%wxmin       )
    CALL allocate_shared_dp_0D(                    grid%xmax       , grid%wxmax       )
    CALL allocate_shared_dp_0D(                    grid%ymin       , grid%wymin       )
    CALL allocate_shared_dp_0D(                    grid%ymax       , grid%wymax       )
    CALL allocate_shared_dp_0D(                    grid%dx         , grid%wdx         )
    CALL allocate_shared_dp_0D(                    grid%tol_dist   , grid%wtol_dist   )
    CALL allocate_shared_int_2D( grid%nx, grid%ny, grid%ij2n       , grid%wij2n       )
    CALL allocate_shared_int_2D( grid%n , 2,       grid%n2ij       , grid%wn2ij       )
    CALL allocate_shared_dp_0D(                    grid%lambda_m   , grid%wlambda_m   )
    CALL allocate_shared_dp_0D(                    grid%phi_m      , grid%wphi_m      )
    CALL allocate_shared_dp_0D(                    grid%beta_stereo, grid%wbeta_stereo)
    CALL allocate_shared_dp_2D(  grid%nx, grid%ny, grid%lon        , grid%wlon        )
    CALL allocate_shared_dp_2D(  grid%nx, grid%ny, grid%lat        , grid%wlat        )

    ! Calculate secondary grid data
    IF (par%master) THEN

      ! Resolution
      grid%dx   = ABS( grid%x( 2) - grid%x( 1))

      ! Safety
      DO i = 1, grid%nx-1
        IF (1._dp - ABS(grid%x( i+1) - grid%x( i)) / grid%dx > 1E-6_dp) CALL crash('file "' // TRIM( filename) // '" has an irregular x-dimension!')
      END DO
      DO j = 1, grid%ny-1
        IF (1._dp - ABS(grid%y( j+1) - grid%y( j)) / grid%dx > 1E-6_dp) CALL crash('file "' // TRIM( filename) // '" has an irregular y-dimension!')
      END DO

      ! Domain size
      grid%xmin = MINVAL( grid%x)
      grid%xmax = MAXVAL( grid%x)
      grid%ymin = MINVAL( grid%y)
      grid%ymax = MAXVAL( grid%y)

      ! Tolerance; points lying within this distance of each other are treated as identical
      grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp

      ! Conversion tables for grid-form vs. vector-form data
      n = 0
      DO i = 1, grid%nx
        IF (MOD(i,2) == 1) THEN
          DO j = 1, grid%ny
            n = n+1
            grid%ij2n( i,j) = n
            grid%n2ij( n,:) = [i,j]
          END DO
        ELSE
          DO j = grid%ny, 1, -1
            n = n+1
            grid%ij2n( i,j) = n
            grid%n2ij( n,:) = [i,j]
          END DO
        END IF
      END DO

    END IF ! IF (par%master) THEN
    CALL sync

    ! Set up parallelisation domains
    CALL partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)

    ! Projection parameters for this region
    IF     (region_name == 'NAM') THEN
      grid%lambda_M     = C%lambda_M_NAM
      grid%phi_M        = C%phi_M_NAM
      grid%beta_stereo  = C%beta_stereo_NAM
    ELSEIF (region_name == 'EAS') THEN
      grid%lambda_M     = C%lambda_M_EAS
      grid%phi_M        = C%phi_M_EAS
      grid%beta_stereo  = C%beta_stereo_EAS
    ELSEIF (region_name == 'GRL') THEN
      grid%lambda_M     = C%lambda_M_GRL
      grid%phi_M        = C%phi_M_GRL
      grid%beta_stereo  = C%beta_stereo_GRL
    ELSEIF (region_name == 'ANT') THEN
      grid%lambda_M     = C%lambda_M_ANT
      grid%phi_M        = C%phi_M_ANT
      grid%beta_stereo  = C%beta_stereo_ANT
    END IF

    ! Lon/lat coordinates
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      CALL inverse_oblique_sg_projection( grid%x( i), grid%y( j), grid%lambda_M, grid%phi_M, grid%beta_stereo, grid%lon( i,j), grid%lat( i,j))
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name, 18)

  END SUBROUTINE setup_xy_grid_from_file

  SUBROUTINE setup_lonlat_grid_from_file( filename, grid, region_name)
    ! Set up a lon/lat-grid from a NetCDF file
    !
    ! Assumes no memory has yet been allocated for the grid at all

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_lonlat_grid_from_file'
    INTEGER                                            :: id_dim_lon, id_dim_lat
    INTEGER                                            :: id_var_lon, id_var_lat
    INTEGER                                            :: i,j
    REAL(dp)                                           :: dlon

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check grid dimensions and variables for validity
    CALL check_lon( filename)
    CALL check_lat( filename)

    ! Allocate memory for the grid size
    CALL allocate_shared_int_0D( grid%nlon, grid%wnlon)
    CALL allocate_shared_int_0D( grid%nlat, grid%wnlat)

    ! Inquire lon and lat dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_lon, id_dim_lon, dim_length = grid%nlon)
    CALL inquire_dim_multiple_options( filename, field_name_options_lat, id_dim_lat, dim_length = grid%nlat)

    ! Allocate memory for lon and lat
    CALL allocate_shared_dp_1D( grid%nlon, grid%lon, grid%wlon)
    CALL allocate_shared_dp_1D( grid%nlat, grid%lat, grid%wlat)

    ! Inquire lon and lat variables
    CALL inquire_var_multiple_options( filename, field_name_options_lon, id_var_lon)
    CALL inquire_var_multiple_options( filename, field_name_options_lat, id_var_lat)

    ! Read lon and lat
    CALL read_var_dp_1D(  filename, id_var_lon, grid%lon)
    CALL read_var_dp_1D(  filename, id_var_lat, grid%lat)

    ! Allocate memory for, and calculate, some secondary grid properties
    CALL allocate_shared_dp_0D( grid%lonmin     , grid%wlonmin     )
    CALL allocate_shared_dp_0D( grid%lonmax     , grid%wlonmax     )
    CALL allocate_shared_dp_0D( grid%latmin     , grid%wlatmin     )
    CALL allocate_shared_dp_0D( grid%latmax     , grid%wlatmax     )
    CALL allocate_shared_dp_0D( grid%dlon       , grid%wdlon       )
    CALL allocate_shared_dp_0D( grid%dlat       , grid%wdlat       )
    CALL allocate_shared_dp_0D( grid%lambda_m   , grid%wlambda_m   )
    CALL allocate_shared_dp_0D( grid%phi_m      , grid%wphi_m      )
    CALL allocate_shared_dp_0D( grid%beta_stereo, grid%wbeta_stereo)

    ! Calculate secondary grid data
    IF (par%master) THEN

      ! Resolution
      grid%dlon = ABS( grid%lon( 2) - grid%lon( 1))
      grid%dlat = ABS( grid%lat( 2) - grid%lat( 1))

      ! Safety
      DO i = 1, grid%nlon - 1
        ! Check for regularity in longitude, but allow for 360-degree jumps
        dlon = MIN( ABS( grid%lon( i+1) - grid%lon( i)), ABS( grid%lon( i+1) + 360._dp - grid%lon( i)))
        IF (ABS( 1._dp - dlon / grid%dlon) > 1E-6_dp) CALL crash('file "' // TRIM( filename) // '" has an irregular longitude dimension!')
      END DO
      DO j = 1, grid%nlat - 1
        IF (ABS( 1._dp - ABS( grid%lat( j+1) - grid%lat( j)) / grid%dlat) > 1E-6_dp) &
          CALL crash('file "' // TRIM( filename) // '" has an irregular latitude dimension!')
      END DO

      ! Domain size
      grid%lonmin = MINVAL( grid%lon)
      grid%lonmax = MAXVAL( grid%lon)
      grid%latmin = MINVAL( grid%lat)
      grid%latmax = MAXVAL( grid%lat)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Projection parameters for this region
    IF     (region_name == 'NAM') THEN
      grid%lambda_M     = C%lambda_M_NAM
      grid%phi_M        = C%phi_M_NAM
      grid%beta_stereo  = C%beta_stereo_NAM
    ELSEIF (region_name == 'EAS') THEN
      grid%lambda_M     = C%lambda_M_EAS
      grid%phi_M        = C%phi_M_EAS
      grid%beta_stereo  = C%beta_stereo_EAS
    ELSEIF (region_name == 'GRL') THEN
      grid%lambda_M     = C%lambda_M_GRL
      grid%phi_M        = C%phi_M_GRL
      grid%beta_stereo  = C%beta_stereo_GRL
    ELSEIF (region_name == 'ANT') THEN
      grid%lambda_M     = C%lambda_M_ANT
      grid%phi_M        = C%phi_M_ANT
      grid%beta_stereo  = C%beta_stereo_ANT
    END IF

    ! Set up parallelisation domains
    CALL partition_list( grid%nlon, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%nlat, par%i, par%n, grid%j1, grid%j2)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_lonlat_grid_from_file

  SUBROUTINE setup_zeta_from_file( filename, nzeta, zeta, wzeta)
    ! Set up a zeta coordinate from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(OUT)   :: nzeta
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   ::  zeta
    INTEGER,                             INTENT(OUT)   :: wzeta

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_zeta_from_file'
    INTEGER                                            :: id_dim_zeta, id_var_zeta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check zeta dimension and variable for validity
    CALL check_zeta( filename)

    ! Inquire zeta dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_zeta, id_dim_zeta, dim_length = nzeta)

    ! Inquire zeta variable
    CALL inquire_var_multiple_options( filename, field_name_options_zeta, id_var_zeta)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( nzeta, zeta, wzeta)

    ! Read zeta from file
    CALL read_var_dp_1D( filename, id_var_zeta, zeta)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_zeta_from_file

  SUBROUTINE setup_z_ocean_from_file( filename, nz_ocean, z_ocean, wz_ocean)
    ! Set up a z_ocean coordinate from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(OUT)   :: nz_ocean
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   ::  z_ocean
    INTEGER,                             INTENT(OUT)   :: wz_ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_z_ocean_from_file'
    INTEGER                                            :: id_dim_z_ocean, id_var_z_ocean

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check z_ocean dimension and variable for validity
    CALL check_z_ocean( filename)

    ! Inquire z_ocean dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_z_ocean, id_dim_z_ocean, dim_length = nz_ocean)

    ! Inquire z_ocean variable
    CALL inquire_var_multiple_options( filename, field_name_options_z_ocean, id_var_z_ocean)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( nz_ocean, z_ocean, wz_ocean)

    ! Read z_ocean from file
    CALL read_var_dp_1D( filename, id_var_z_ocean, z_ocean)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_z_ocean_from_file

  ! ===== Determine indexing and dimension directions =====
  ! =======================================================

  SUBROUTINE determine_xy_indexing( filename, var_name, indexing, xdir, ydir)
    ! Determine the indexing and dimension directions of a variable in an x/y-grid file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=256),                  INTENT(OUT)   :: indexing
    CHARACTER(LEN=256),                  INTENT(OUT)   :: xdir
    CHARACTER(LEN=256),                  INTENT(OUT)   :: ydir

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_xy_indexing'
    INTEGER                                            :: id_dim_x, id_dim_y
    INTEGER                                            :: nx, ny
    REAL(dp), DIMENSION(:    ), POINTER                :: x, y
    INTEGER                                            :: wx, wy
    INTEGER                                            :: id_var_x, id_var_y, id_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if the x and y dimensions and variables of this file are valid
    CALL check_x( filename)
    CALL check_y( filename)

    ! Inquire x and y dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x, id_dim_x, dim_length = nx)
    CALL inquire_dim_multiple_options( filename, field_name_options_y, id_dim_y, dim_length = ny)

    ! Allocate memory for x and y
    CALL allocate_shared_dp_1D( nx, x, wx)
    CALL allocate_shared_dp_1D( ny, y, wy)

    ! Inquire x and y variables
    CALL inquire_var_multiple_options( filename, field_name_options_x, id_var_x)
    CALL inquire_var_multiple_options( filename, field_name_options_y, id_var_y)

    ! Read x and y
    CALL read_var_dp_1D(  filename, id_var_x, x)
    CALL read_var_dp_1D(  filename, id_var_y, y)

    ! Determine directions of x and y
    IF (x( 2) > x( 1)) THEN
      xdir = 'normal'
    ELSE
      xdir = 'reverse'
    END IF
    IF (y( 2) > y( 1)) THEN
      ydir = 'normal'
    ELSE
      ydir = 'reverse'
    END IF

    ! Inquire dimensions of the specified field variable
    CALL inquire_var_multiple_options( filename, var_name, id_var, dims_of_var = dims_of_var)

    ! Determine indexing
    IF     (dims_of_var( 1) == id_dim_x .AND. dims_of_var( 2) == id_dim_y) THEN
      indexing = 'xy'
    ELSEIF (dims_of_var( 1) == id_dim_y .AND. dims_of_var( 2) == id_dim_x) THEN
      indexing = 'yx'
    ELSE
      CALL crash('x and y are not the first two dimensions of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wx)
    CALL deallocate_shared( wy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_xy_indexing

  SUBROUTINE determine_lonlat_indexing( filename, var_name, indexing, londir, latdir)
    ! Determine the indexing and dimension directions of a variable in a lon/lat-grid file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    CHARACTER(LEN=256),                  INTENT(OUT)   :: indexing
    CHARACTER(LEN=256),                  INTENT(OUT)   :: londir
    CHARACTER(LEN=256),                  INTENT(OUT)   :: latdir

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_lonlat_indexing'
    INTEGER                                            :: id_dim_lon, id_dim_lat
    INTEGER                                            :: nlon, nlat
    REAL(dp), DIMENSION(:    ), POINTER                :: lon, lat
    INTEGER                                            :: wlon, wlat
    INTEGER                                            :: id_var_lon, id_var_lat, id_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if the lon and lat dimensions and variables of this file are valid
    CALL check_lon( filename)
    CALL check_lat( filename)

    ! Inquire lon and lat dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_lon, id_dim_lon, dim_length = nlon)
    CALL inquire_dim_multiple_options( filename, field_name_options_lat, id_dim_lat, dim_length = nlat)

    ! Allocate memory for lon and lat
    CALL allocate_shared_dp_1D( nlon, lon, wlon)
    CALL allocate_shared_dp_1D( nlat, lat, wlat)

    ! Inquire lon and lon variables
    CALL inquire_var_multiple_options( filename, field_name_options_lon, id_var_lon)
    CALL inquire_var_multiple_options( filename, field_name_options_lat, id_var_lat)

    ! Read lon and lat
    CALL read_var_dp_1D(  filename, id_var_lon, lon)
    CALL read_var_dp_1D(  filename, id_var_lat, lat)

    ! Determine directions of x and y
    IF (lon( 2) > lon( 1)) THEN
      londir = 'normal'
    ELSE
      londir = 'reverse'
    END IF
    IF (lat( 2) > lat( 1)) THEN
      latdir = 'normal'
    ELSE
      latdir = 'reverse'
    END IF

    ! Inquire dimensions of the specified field variable
    CALL inquire_var_multiple_options( filename, var_name, id_var, dims_of_var = dims_of_var)

    ! Determine indexing
    IF     (dims_of_var( 1) == id_dim_lon .AND. dims_of_var( 2) == id_dim_lat) THEN
      indexing = 'lonlat'
    ELSEIF (dims_of_var( 1) == id_dim_lat .AND. dims_of_var( 2) == id_dim_lon) THEN
      indexing = 'latlon'
    ELSE
      CALL crash('longitude and latitude are not the first two dimensions of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wlon)
    CALL deallocate_shared( wlat)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_lonlat_indexing

  ! ===== Longitude corrections to a lon/lat grid + data field =====
  ! ================================================================

  SUBROUTINE correct_longitude_shifts_and_range_2D( filename, grid_lonlat, d_grid)
    ! Make sure longitude is bounded between 0 and 360, and increases monotonically

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid_lonlat
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'correct_longitude_shifts_and_range_2D'
    INTEGER                                            :: i,j
    LOGICAL                                            :: is_correct
    INTEGER                                            :: n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the grid is already correct, do nothing
    is_correct = .TRUE.
    IF (ANY( grid_lonlat%lon < 0._dp) .OR. ANY( grid_lonlat%lon > 360._dp)) is_correct = .FALSE.
    DO i = 1, grid_lonlat%nlon - 1
      IF (grid_lonlat%lon( i+1) <= grid_lonlat%lon( i)) is_correct = .FALSE.
    END DO
    CALL sync

    IF (is_correct) THEN
       CALL finalise_routine( routine_name)
       RETURN
    END IF

    ! Limit values to [0,360]
    IF (par%master) THEN
      DO i = 1, grid_lonlat%nlon
        IF     (grid_lonlat%lon( i) <   0._dp) THEN
          grid_lonlat%lon( i) = grid_lonlat%lon( i) + 360._dp
        ELSEIF (grid_lonlat%lon( i) > 360._dp) THEN
          grid_lonlat%lon( i) = grid_lonlat%lon( i) - 360._dp
        END IF
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Fix shifts
    n = 0
    DO i = 1, grid_lonlat%nlon-1
      IF (grid_lonlat%lon( i) > grid_lonlat%lon( i+1)) THEN
        n = i
        EXIT
      END IF
    END DO

    IF (n > 0) THEN

      ! Fix lon
      IF (par%master) THEN
        grid_lonlat%lon = [grid_lonlat%lon( n+1:grid_lonlat%nlon), grid_lonlat%lon( 1:n)]
      END IF
      CALL sync

      ! Fix data field
      DO j = grid_lonlat%j1, grid_lonlat%j2
        d_grid( :,j) = [d_grid( n+1:grid_lonlat%nlon,j), d_grid( 1:n,j)]
      END DO
      CALL sync

    END IF ! IF (n > 0) THEN

    ! The grid should now be correct
    is_correct = .TRUE.
    IF (ANY( grid_lonlat%lon < 0._dp) .OR. ANY( grid_lonlat%lon > 360._dp)) is_correct = .FALSE.
    DO i = 1, grid_lonlat%nlon - 1
      IF (grid_lonlat%lon( i+1) <= grid_lonlat%lon( i)) is_correct = .FALSE.
    END DO
    IF (.NOT. is_correct) CALL crash('something is seriously wrong with the longitude of file "' // TRIM( filename) // '"!')

    PRINT*,('This should not appear for now')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE correct_longitude_shifts_and_range_2D

  SUBROUTINE correct_longitude_shifts_and_range_3D( filename, grid_lonlat, d_grid)
    ! Make sure longitude is bounded between 0 and 360, and increases monotonically

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    TYPE(type_grid_lonlat),              INTENT(INOUT) :: grid_lonlat
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'correct_longitude_shifts_and_range_3D'
    INTEGER                                            :: i,j,k
    LOGICAL                                            :: is_correct
    INTEGER                                            :: n

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If the grid is already correct, do nothing
    is_correct = .TRUE.
    IF (ANY( grid_lonlat%lon < 0._dp) .OR. ANY( grid_lonlat%lon > 360._dp)) is_correct = .FALSE.
    DO i = 1, grid_lonlat%nlon - 1
      IF (grid_lonlat%lon( i+1) <= grid_lonlat%lon( i)) is_correct = .FALSE.
    END DO
    CALL sync
    IF (is_correct) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Limit values to [0,360]
    IF (par%master) THEN
      DO i = 1, grid_lonlat%nlon
        IF     (grid_lonlat%lon( i) <   0._dp) THEN
          grid_lonlat%lon( i) = grid_lonlat%lon( i) + 360._dp
        ELSEIF (grid_lonlat%lon( i) > 360._dp) THEN
          grid_lonlat%lon( i) = grid_lonlat%lon( i) - 360._dp
        END IF
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Fix shifts
    n = 0
    DO i = 1, grid_lonlat%nlon-1
      IF (grid_lonlat%lon( i) > grid_lonlat%lon( i+1)) THEN
        n = i
        EXIT
      END IF
    END DO

    IF (n > 0) THEN

      ! Fix lon
      IF (par%master) THEN
        grid_lonlat%lon = [grid_lonlat%lon( n+1:grid_lonlat%nlon), grid_lonlat%lon( 1:n)]
      END IF
      CALL sync

      ! Fix data field
      DO j = grid_lonlat%j1, grid_lonlat%j2
      DO k = 1, SIZE( d_grid,3)
        d_grid( :,j,k) = [d_grid( n+1:grid_lonlat%nlon,j,k), d_grid( 1:n,j,k)]
      END DO
      END DO
      CALL sync

    END IF ! IF (n > 0) THEN

    ! The grid should now be correct
    is_correct = .TRUE.
    IF (ANY( grid_lonlat%lon < 0._dp) .OR. ANY( grid_lonlat%lon > 360._dp)) is_correct = .FALSE.
    DO i = 1, grid_lonlat%nlon - 1
      IF (grid_lonlat%lon( i+1) <= grid_lonlat%lon( i)) is_correct = .FALSE.
    END DO
    IF (.NOT. is_correct) CALL crash('something is seriously wrong with the longitude of file "' // TRIM( filename) // '"!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE correct_longitude_shifts_and_range_3D

END MODULE netcdf_input_module
