MODULE netcdf_basic_module

  ! Basic NetCDF routines
  ! =====================
  !
  ! These routines handle some very basic stuff; opening and closing
  ! NetCDF files in different modes, inquiring dimensions and variables,
  ! creating dimensions, variables and attributes;.
  !
  ! Also contains the routines for setting up an x/y-grid, lon/lat-grid,
  ! from a NetCDF file.

! ===== Preamble =====
! ====================

  ! Import basic functionality
! include <petsc/finclude/petscksp.h> 
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, &
                                             allocate_shared_int_0D,   allocate_shared_dp_0D, &
                                             allocate_shared_int_1D,   allocate_shared_dp_1D, &
                                             allocate_shared_int_2D,   allocate_shared_dp_2D, &
                                             allocate_shared_int_3D,   allocate_shared_dp_3D, &
                                             deallocate_shared
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D

  ! Import specific functionality
  USE data_types_module,               ONLY: type_grid, type_grid_lonlat
  USE netcdf,                          ONLY: NF90_NOERR, NF90_OPEN, NF90_CLOSE, NF90_NOWRITE, NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, &
                                             NF90_INQ_VARID, NF90_INQUIRE_VARIABLE, NF90_MAX_VAR_DIMS, NF90_GET_VAR, &
                                             NF90_CREATE, NF90_NOCLOBBER, NF90_NETCDF4, NF90_ENDDEF, NF90_REDEF, NF90_DEF_DIM, NF90_DEF_VAR, &
                                             NF90_PUT_ATT, NF90_WRITE, NF90_INT, NF90_INT64, NF90_FLOAT, NF90_DOUBLE, NF90_PUT_VAR, NF90_UNLIMITED, &
                                             NF90_INQUIRE_ATTRIBUTE
  USE utilities_module,                ONLY: inverse_oblique_sg_projection

  IMPLICIT NONE

  ! NetCDF error code
  INTEGER                 :: nerr

  ! Possible names for different dimensions and variables
  ! =====================================================

  ! Different options for the name of a dimension or variable can now be tried.
  ! They are separated by a double vertical bar ||

  ! Dimensions
  CHARACTER(LEN=256), PARAMETER :: field_name_options_x              = 'x||X||x1||X1||nx||NX||x-coordinate||X-coordinate||easting||Easting'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_y              = 'y||Y||y1||Y1||ny||NY||y-coordinate||Y-coordinate||northing||Northing'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_zeta           = 'zeta||Zeta'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_z_ocean        = 'depth||Depth||z_ocean'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_lon            = 'lon||Lon||long||Long||longitude||Longitude'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_lat            = 'lat||Lat||latitude||Latitude'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_time           = 'time||Time||t||nt'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_month          = 'month||Month'
  
  ! Variables
  CHARACTER(LEN=256), PARAMETER :: field_name_options_Hi             = 'Hi||thickness||lithk'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_Hb             = 'Hb||bed||topg'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_Hs             = 'Hs||surface||orog'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_dHb            = 'dHb'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_SL             = 'SL'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_Ti             = 'Ti'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_v_surf         = 'v_surf||VX'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_u_surf         = 'u_surf||UX'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_T_ocean        = 'T_ocean||t_an'
  CHARACTER(LEN=256), PARAMETER :: field_name_options_S_ocean        = 'S_ocean||s_an'

CONTAINS

  ! Check if a file contains all the variables and dimensions describing
  ! an x/y-grid, a lon/lat-grid
  SUBROUTINE inquire_xy_grid( filename, has_xy_grid)
    ! Inquire if a NetCDF file contains all the dimensions and variables
    ! describing a regular x/y-grid.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_xy_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_xy_grid'
    INTEGER                                            :: id_dim_x, id_dim_y
    INTEGER                                            :: id_var_x, id_var_y

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Look for x and y dimensions and variables
    CALL inquire_dim_multiple_options( filename, field_name_options_x, id_dim_x)
    CALL inquire_dim_multiple_options( filename, field_name_options_y, id_dim_y)
    CALL inquire_var_multiple_options( filename, field_name_options_x, id_var_x)
    CALL inquire_var_multiple_options( filename, field_name_options_y, id_var_y)

    ! Check if everything is there
    has_xy_grid = .TRUE.

    IF (id_dim_x              == -1) has_xy_grid = .FALSE.
    IF (id_dim_y              == -1) has_xy_grid = .FALSE.
    IF (id_var_x              == -1) has_xy_grid = .FALSE.
    IF (id_var_y              == -1) has_xy_grid = .FALSE.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_xy_grid

  SUBROUTINE inquire_lonlat_grid( filename, has_lonlat_grid)
    ! Inquire if a NetCDF file contains all the dimensions and variables
    ! describing a regular lon/lat-grid.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_lonlat_grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_lonlat_grid'
    INTEGER                                            :: id_dim_lon, id_dim_lat
    INTEGER                                            :: id_var_lon, id_var_lat

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Look for x and y dimensions and variables
    CALL inquire_dim_multiple_options( filename, field_name_options_lon, id_dim_lon)
    CALL inquire_dim_multiple_options( filename, field_name_options_lat, id_dim_lat)
    CALL inquire_var_multiple_options( filename, field_name_options_lon, id_var_lon)
    CALL inquire_var_multiple_options( filename, field_name_options_lat, id_var_lat)

    ! Check if everything is there
    has_lonlat_grid = .TRUE.

    IF (id_dim_lon            == -1) has_lonlat_grid = .FALSE.
    IF (id_dim_lat            == -1) has_lonlat_grid = .FALSE.
    IF (id_var_lon            == -1) has_lonlat_grid = .FALSE.
    IF (id_var_lat            == -1) has_lonlat_grid = .FALSE.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_lonlat_grid

  ! Inquire if a file contains the variable and dimension for
  ! zeta, z_ocean, time, or months
  SUBROUTINE inquire_zeta( filename, has_zeta)
    ! Inquire if a NetCDF file contains a zeta dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_zeta

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_zeta'
    INTEGER                                            :: id_dim_zeta, id_var_zeta

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Look for zeta dimension and variable
    CALL inquire_dim_multiple_options( filename, field_name_options_zeta, id_dim_zeta)
    CALL inquire_var_multiple_options( filename, field_name_options_zeta, id_var_zeta)

    ! Check if everything is there
    has_zeta = .TRUE.

    IF (id_dim_zeta           == -1) has_zeta = .FALSE.
    IF (id_var_zeta           == -1) has_zeta = .FALSE.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_zeta

  SUBROUTINE inquire_z_ocean( filename, has_z_ocean)
    ! Inquire if a NetCDF file contains a z_ocean dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_z_ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_z_ocean'
    INTEGER                                            :: id_dim_z_ocean, id_var_z_ocean

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Look for z_ocean dimension and variable
    CALL inquire_dim_multiple_options( filename, field_name_options_z_ocean, id_dim_z_ocean)
    CALL inquire_var_multiple_options( filename, field_name_options_z_ocean, id_var_z_ocean)

    ! Check if everything is there
    has_z_ocean = .TRUE.

    IF (id_dim_z_ocean           == -1) has_z_ocean = .FALSE.
    IF (id_var_z_ocean           == -1) has_z_ocean = .FALSE.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_z_ocean

  SUBROUTINE inquire_month( filename, has_month)
    ! Inquire if a NetCDF file contains a month dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_month

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_month'
    INTEGER                                            :: id_dim_month, id_var_month

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Look for month dimension and variable
    CALL inquire_dim_multiple_options( filename, field_name_options_month, id_dim_month)
    CALL inquire_var_multiple_options( filename, field_name_options_month, id_var_month)

    ! Check if everything is there
    has_month = .TRUE.

    IF (id_dim_month           == -1) has_month = .FALSE.
    IF (id_var_month           == -1) has_month = .FALSE.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_month

  SUBROUTINE inquire_time( filename, has_time)
    ! Inquire if a NetCDF file contains a time dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    LOGICAL,                             INTENT(OUT)   :: has_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_time'
    INTEGER                                            :: id_dim_time, id_var_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Look for time dimension and variable
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)
    CALL inquire_var_multiple_options( filename, field_name_options_time, id_var_time)

    ! Check if everything is there
    has_time = .TRUE.

    IF (id_dim_time           == -1) has_time = .FALSE.
    IF (id_var_time           == -1) has_time = .FALSE.

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_time

  SUBROUTINE find_timeframe( filename, time, ti)
    ! Find the timeframe in the file that is closest to the desired time.
    ! If the file has no time dimension or variable, throw an error.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    REAL(dp),                            INTENT(IN)    :: time
    INTEGER,                             INTENT(OUT)   :: ti

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_timeframe'
    INTEGER                                            :: nt, id_dim_time, id_var_time
    INTEGER                                            :: var_type
    REAL(dp), DIMENSION(:    ), POINTER                ::  time_from_file
    INTEGER                                            :: wtime_from_file
    INTEGER                                            :: tii
    REAL(dp)                                           :: dt_min
    INTEGER,  DIMENSION(:    ), POINTER                ::  time_from_file_int
    INTEGER                                            :: wtime_from_file_int

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file contains a valid time dimension and variable
    CALL check_time( filename)

    ! Inquire size of time dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time, dim_length = nt)

    ! Inquire time variable ID
    CALL inquire_var_multiple_options( filename, field_name_options_time, id_var_time, var_type = var_type)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( nt, time_from_file, wtime_from_file)

    ! Read time from file
    IF (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE) THEN
      CALL read_var_dp_1D( filename, id_var_time, time_from_file)
    ELSEIF (var_type == NF90_INT .OR. var_type == NF90_INT64) THEN
      CALL allocate_shared_int_1D( nt, time_from_file_int, wtime_from_file_int)
      CALL read_var_int_1D( filename, id_var_time, time_from_file_int)
      DO tii = 1, nt
        time_from_file( tii) = REAL( time_from_file_int( tii),dp)
      END DO
      CALL deallocate_shared( wtime_from_file_int)
    ELSE
      CALL crash('whaa!')
    END IF

    ! Find timeframe closest to desired time
    IF (time_from_file( 1) > time) THEN
      ! Desired time beyond lower limit
      CALL warning('desired timeframe at t = {dp_01} before start of file time for file "' // TRIM( filename) // '"; reading data from t = {dp_02} instead!', dp_01 = time, dp_02 = time_from_file( 1))
      ti = 1
    ELSEIF (time_from_file( nt) < time) THEN
      ! Desired time beyond upper limit
      CALL warning('desired timeframe at t = {dp_01} after end of file time for file "' // TRIM( filename) // '"; reading data from t = {dp_02} instead!', dp_01 = time, dp_02 = time_from_file( nt))
      ti = nt
    ELSE
      ! Desired time is within the file time
      dt_min = HUGE( 1._dp)
      DO tii = 1, nt
        IF (ABS( time_from_file( tii) - time) < dt_min) THEN
          ti = tii
          dt_min = ABS( time_from_file( tii) - time)
        END IF
      END DO
      IF (dt_min > 0._dp) THEN
        CALL warning('desired timeframe at t = {dp_01} not present in file "' // TRIM( filename) // '"; reading data from closest match at t = {dp_02} instead!', dp_01 = time, dp_02 = time_from_file( ti))
      END IF
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wtime_from_file)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE find_timeframe

  ! ===== Safety checks on variables and dimensions =====
  ! =====================================================

  ! x/y-grid dimensions
  SUBROUTINE check_x( filename)
    ! Check if this file contains a valid x dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_x'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), POINTER                :: x
    INTEGER                                            :: wx
    REAL(dp)                                           :: dx, dxp
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_x, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid x dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('dimension "' // TRIM( dim_name) // '" in file "' // TRIM( filename) // '" is unlimited!')
    IF (n < 1) CALL crash('dimension "' // TRIM( dim_name) // '" in file "' // TRIM( filename) // '" has length {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multiple_options( filename, field_name_options_x, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid x variable could be found in file "' // TRIM( filename) // '"!')

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('variable "' // TRIM( var_name) // &
      '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check variable dimension
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have ' // TRIM( dim_name) // ' as a dimension!')

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( n, x, wx)

    ! Read variable
    CALL read_var_dp_1D( filename, id_var, x)

    ! Check validity
    CALL check_for_NaN_dp_1D( x, 'x')

    ! Check grid spacing
    dx = x( 2) - x( 1)
    DO i = 2, n
      dxp = x( i) - x( i-1)
      IF (ABS( 1._dp - dxp / dx) > 1E-5_dp) CALL crash('x coordinate in file "' // TRIM( filename) // '" is irregular!')
    END DO

    ! Clean up after yourself
    CALL deallocate_shared( wx)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_x

  SUBROUTINE check_y( filename)
    ! Check if this file contains a valid y dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_y'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), POINTER                :: y
    INTEGER                                            :: wy
    REAL(dp)                                           :: dy, dyp
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_y, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid y dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('dimension "' // TRIM( dim_name) // '" in file "' // TRIM( filename) // '" is unlimited!')
    IF (n < 1) CALL crash('dimension "' // TRIM( dim_name) // '" in file "' // TRIM( filename) // '" has length {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multiple_options( filename, field_name_options_y, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid y variable could be found in file "' // TRIM( filename) // '"!')

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('variable "' // TRIM( var_name) // &
      '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check variable dimension
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have ' // TRIM( dim_name) // ' as a dimension!')

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( n, y, wy)

    ! Read variable
    CALL read_var_dp_1D( filename, id_var, y)

    ! Check validity
    CALL check_for_NaN_dp_1D( y, 'y')

    ! Check grid spacing
    dy = y( 2) - y( 1)
    DO i = 2, n
      dyp = y( i) - y( i-1)
      IF (ABS( 1._dp - dyp / dy) > 1E-5_dp) CALL crash('y coordinate in file "' // TRIM( filename) // '" is irregular!')
    END DO

    ! Clean up after yourself
    CALL deallocate_shared( wy)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_y

  ! lon/lat-grid dimensions
  SUBROUTINE check_lon( filename)
    ! Check if this file contains a valid longitude dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lon'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), POINTER                :: lon
    INTEGER                                            :: wlon
    REAL(dp)                                           :: dlon, dlonp
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_lon, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid longitude dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('longitude dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (n < 1) CALL crash('longitude dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multiple_options( filename, field_name_options_lon, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid longitude variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('longitude variable in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 1) CALL crash('longitude variable in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('longitude variable in file "' // TRIM( filename) // '" does not have longitude as a dimension!')

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( n, lon, wlon)

    ! Read variable
    CALL read_var_dp_1D( filename, id_var, lon)

    ! Check validity
    CALL check_for_NaN_dp_1D( lon, 'lon')

    ! Check grid spacing
    dlon = MINVAL([ ABS(lon( 2) - lon( 1)), ABS( lon( 2) + 360._dp - lon( 1)), ABS( lon( 2) - 360._dp - lon( 1)) ])
    DO i = 2, n
      dlonp = MINVAL([ ABS(lon( i) - lon( i-1)), ABS( lon( i) + 360._dp - lon( i-1)), ABS( lon( i) - 360._dp - lon( i-1)) ])
      IF (ABS( 1._dp - dlonp / dlon) > 1E-5_dp) CALL crash('lon coordinate in file "' // TRIM( filename) // '" is irregular!')
    END DO

    ! Clean up after yourself
    CALL deallocate_shared( wlon)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lon

  SUBROUTINE check_lat( filename)
    ! Check if this file contains a valid latitude dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lat'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), POINTER                :: lat
    INTEGER                                            :: wlat
    REAL(dp)                                           :: dlat, dlatp
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_lat, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid latitude dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('latitude dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (n < 1) CALL crash('latitude dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multiple_options( filename, field_name_options_lat, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid latitude variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('latitude variable in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 1) CALL crash('latitude variable in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('latitude variable in file "' // TRIM( filename) // '" does not have latitude as a dimension!')

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( n, lat, wlat)

    ! Read variable
    CALL read_var_dp_1D( filename, id_var, lat)

    ! Check validity
    CALL check_for_NaN_dp_1D( lat, 'lat')

    ! Check grid spacing
    dlat = lat( 2) - lat( 1)
    DO i = 2, n
      dlatp = lat( i) - lat( i-1)
      IF (ABS( 1._dp - dlatp / dlat) > 1E-5_dp) CALL crash('latitude coordinate in file "' // TRIM( filename) // '" is irregular!')
    END DO

    ! Clean up after yourself
    CALL deallocate_shared( wlat)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lat
  
  ! Zeta, z_ocean, month, time dimensions
  SUBROUTINE check_zeta( filename)
    ! Check if this file contains a valid zeta dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_zeta'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), POINTER                :: zeta
    INTEGER                                            :: wzeta
    INTEGER                                            :: k

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_zeta, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid zeta dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('zeta dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (n < 1) CALL crash('zeta dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multiple_options( filename, field_name_options_zeta, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid zeta variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('zeta variable in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 1) CALL crash('zeta variable in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('zeta variable in file "' // TRIM( filename) // '" does not have zeta as a dimension!')

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( n, zeta, wzeta)

    ! Read variable
    CALL read_var_dp_1D( filename, id_var, zeta)

    ! Check validity
    CALL check_for_NaN_dp_1D( zeta, 'zeta')

    IF (zeta( 1) /= 0._dp) CALL crash('zeta in file "' // TRIM( filename) // '" does not start at zero!')
    IF (zeta( n) /= 1._dp) CALL crash('zeta in file "' // TRIM( filename) // '" does not end at one!')

    DO k = 2, n
      IF (zeta( k) <= zeta( k-1)) CALL crash('zeta in file "' // TRIM( filename) // '" does not increase monotonously!')
    END DO

    ! Clean up after yourself
    CALL deallocate_shared( wzeta)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_zeta

  SUBROUTINE check_z_ocean( filename)
    ! Check if this file contains a valid z_ocean dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_z_ocean'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), POINTER                :: z_ocean
    INTEGER                                            :: wz_ocean

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_z_ocean, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid z_ocean dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('z_ocean dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (n < 1) CALL crash('z_ocean dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multiple_options( filename, field_name_options_z_ocean, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid z_ocean variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) CALL crash('z_ocean variable in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    IF (ndims_of_var /= 1) CALL crash('z_ocean variable in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('z_ocean variable in file "' // TRIM( filename) // '" does not have z_ocean as a dimension!')

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( n, z_ocean, wz_ocean)

    ! Read variable
    CALL read_var_dp_1D( filename, id_var, z_ocean)

    ! Check validity
    CALL check_for_NaN_dp_1D( z_ocean, 'z_ocean')

    ! Clean up after yourself
    CALL deallocate_shared( wz_ocean)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_z_ocean

  SUBROUTINE check_month( filename)
    ! Check if this file contains a valid month dimension (we don't really care about the variable)

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_month'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_month, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid month dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n == NF90_UNLIMITED) CALL crash('month dimension in file "' // TRIM( filename) // '" is unlimited!')
    IF (n /= 12) CALL crash('month dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_month

  SUBROUTINE check_time( filename)
    ! Check if this file contains a valid time dimension and variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_time'
    INTEGER                                            :: id_dim
    INTEGER                                            :: n
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: id_var
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    REAL(dp), DIMENSION(:    ), POINTER                :: time
    INTEGER                                            :: wtime
    INTEGER,  DIMENSION(:    ), POINTER                :: time_int
    INTEGER                                            :: wtime_int
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire dimension
    CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim, dim_length = n, dim_name = dim_name)

    ! Safety checks on dimension
    IF (id_dim == -1) CALL crash('no valid time dimension could be found in file "' // TRIM( filename) // '"!')
    IF (n < 0) CALL crash('time dimension in file "' // TRIM( filename) // '" has length n = {int_01}!', int_01  = n)

    ! Inquire variable
    CALL inquire_var_multiple_options( filename, field_name_options_time, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    IF (id_var == -1) CALL crash('no valid time variable could be found in file "' // TRIM( filename) // '"!')
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE .OR. var_type == NF90_INT .OR. var_type == NF90_INT64)) &
      CALL crash('time variable in file "' // TRIM( filename) // '" is not of type NF90_FLOAT, NF90_DOUBLE, or NF90_INT!')
    IF (ndims_of_var /= 1) CALL crash('time variable in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    IF (dims_of_var( 1) /= id_dim) CALL crash('time variable in file "' // TRIM( filename) // '" does not have time as a dimension!')

    ! For new output files, time is still empty. If it's not, check if entries are valid
    IF (n > 0) THEN

      ! Allocate shared memory
      CALL allocate_shared_dp_1D( n, time, wtime)

      ! Read variable
      IF (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE) THEN
        CALL read_var_dp_1D( filename, id_var, time)
      ELSEIF (var_type == NF90_INT .OR. var_type == NF90_INT64) THEN
        CALL allocate_shared_int_1D( n, time_int, wtime_int)
        CALL read_var_int_1D( filename, id_var, time_int)
        DO ti = 1, n
          time( ti) = REAL( time_int( ti),dp)
        END DO
        CALL deallocate_shared( wtime_int)
      ELSE
        CALL crash('whaa!')
      END IF

      ! Check validity
      CALL check_for_NaN_dp_1D( time, 'time')

      ! Clean up after yourself
      CALL deallocate_shared( wtime)

    END IF ! IF (n > 0) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_time

  ! x/y-grid field variables
  SUBROUTINE check_xy_grid_field_int_2D(            filename, var_name, should_have_time)
    ! Check if this file contains a 2-D x/y-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_xy_grid_field_int_2D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid x and y dimensions and variables
    CALL check_x( filename)
    CALL check_y( filename)

    ! Inquire x,y dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x, id_dim_x)
    CALL inquire_dim_multiple_options( filename, field_name_options_y, id_dim_y)

    ! Inquire variable
    CALL inquire_var( filename, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. var_type == NF90_INT) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    END IF

    ! Check x,y dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_x)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have x as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_y)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have y as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 2) THEN
          ! The variable only has x,y as dimensions.
        ELSE
          IF (ndims_of_var == 3) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has three dimensions, but the third one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have x,y as dimensions
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename)

        ! Inquire the time dimension
        CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have x,y as dimensions

        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_xy_grid_field_int_2D

  SUBROUTINE check_xy_grid_field_dp_2D(             filename, var_name, should_have_time)
    ! Check if this file contains a 2-D x/y-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_xy_grid_field_dp_2D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid x and y dimensions and variables
    CALL check_x( filename)
    CALL check_y( filename)

    ! Inquire x,y dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x, id_dim_x)
    CALL inquire_dim_multiple_options( filename, field_name_options_y, id_dim_y)

    ! Inquire variable
    CALL inquire_var( filename, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check x,y dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_x)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have x as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_y)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have y as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 2) THEN
          ! The variable only has x,y as dimensions.
        ELSE
          IF (ndims_of_var == 3) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has three dimensions, but the third one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have x,y as dimensions
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename)

        ! Inquire the time dimension
        CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have x,y as dimensions

        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_xy_grid_field_dp_2D

  SUBROUTINE check_xy_grid_field_dp_2D_monthly(     filename, var_name, should_have_time)
    ! Check if this file contains a 2-D monthly x/y-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_xy_grid_field_dp_2D_monthly'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_month, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid x and y dimensions and variables
    CALL check_x(     filename)
    CALL check_y(     filename)
    CALL check_month( filename)

    ! Inquire x,y dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x    , id_dim_x    )
    CALL inquire_dim_multiple_options( filename, field_name_options_y    , id_dim_y    )
    CALL inquire_dim_multiple_options( filename, field_name_options_month, id_dim_month)

    ! Inquire variable
    CALL inquire_var( filename, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check x,y dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_x    )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have x as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_y    )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have y as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_month)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have month as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 3) THEN
          ! The variable only has x,y,m as dimensions.
        ELSE
          IF (ndims_of_var == 4) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has four dimensions, but the fourth one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have x,y,m as dimensions
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename)

        ! Inquire the time dimension
        CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have x,y,m as dimensions

        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_xy_grid_field_dp_2D_monthly

  SUBROUTINE check_xy_grid_field_dp_3D(             filename, var_name, should_have_time)
    ! Check if this file contains a 3-D x/y-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_xy_grid_field_dp_3D'
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_zeta, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid x and y dimensions and variables
    CALL check_x(    filename)
    CALL check_y(    filename)
    CALL check_zeta( filename)

    ! Inquire x,y dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_x   , id_dim_x   )
    CALL inquire_dim_multiple_options( filename, field_name_options_y   , id_dim_y   )
    CALL inquire_dim_multiple_options( filename, field_name_options_zeta, id_dim_zeta)

    ! Inquire variable
    CALL inquire_var( filename, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check x,y dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_x   )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have x as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_y   )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have y as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_zeta)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have zeta as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 3) THEN
          ! The variable only has x,y,zeta as dimensions.
        ELSE
          IF (ndims_of_var == 4) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has four dimensions, but the fourth one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have x,y,zeta as dimensions
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename)

        ! Inquire the time dimension
        CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have x,y,zeta as dimensions

        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_xy_grid_field_dp_3D

  ! lon/lat-grid field variables
  SUBROUTINE check_lonlat_grid_field_int_2D(        filename, var_name, should_have_time)
    ! Check if this file contains a 2-D lon/lat-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lonlat_grid_field_int_2D'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid lon and lat dimensions and variables
    CALL check_lon( filename)
    CALL check_lat( filename)

    ! Inquire lon,lat dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_lon, id_dim_lon)
    CALL inquire_dim_multiple_options( filename, field_name_options_lat, id_dim_lat)

    ! Inquire variable
    CALL inquire_var( filename, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. var_type == NF90_INT) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')
    END IF

    ! Check lon,lat dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_lon)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have longitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_lat)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have latitude as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 2) THEN
          ! The variable only has lon,lat as dimensions.
        ELSE
          IF (ndims_of_var == 3) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has three dimensions, but the third one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have lon,lat as dimensions
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename)

        ! Inquire the time dimension
        CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have lon,lat as dimensions

        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lonlat_grid_field_int_2D

  SUBROUTINE check_lonlat_grid_field_dp_2D(         filename, var_name, should_have_time)
    ! Check if this file contains a 2-D lon/lat-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lonlat_grid_field_dp_2D'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid lon and lat dimensions and variables
    CALL check_lon( filename)
    CALL check_lat( filename)

    ! Inquire lon,lat dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_lon, id_dim_lon)
    CALL inquire_dim_multiple_options( filename, field_name_options_lat, id_dim_lat)

    ! Inquire variable
    CALL inquire_var( filename, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check lon,lat dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_lon)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have longitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_lat)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have latitude as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 2) THEN
          ! The variable only has lon,lat as dimensions.
        ELSE
          IF (ndims_of_var == 3) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has three dimensions, but the third one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have lon,lat as dimensions
        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename)

        ! Inquire the time dimension
        CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have lon,lat as dimensions

        IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lonlat_grid_field_dp_2D

  SUBROUTINE check_lonlat_grid_field_dp_2D_monthly( filename, var_name, should_have_time)
    ! Check if this file contains a 2-D monthly lon/lat-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lonlat_grid_field_dp_2D_monthly'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_month, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid lon and lat dimensions and variables
    CALL check_lon(   filename)
    CALL check_lat(   filename)
    CALL check_month( filename)

    ! Inquire lon,lat dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_lon  , id_dim_lon  )
    CALL inquire_dim_multiple_options( filename, field_name_options_lat  , id_dim_lat  )
    CALL inquire_dim_multiple_options( filename, field_name_options_month, id_dim_month)

    ! Inquire variable
    CALL inquire_var( filename, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check lon,lat dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_lon  )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have longitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_lat  )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have latitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_month)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have month as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 3) THEN
          ! The variable only has lon,lat,m as dimensions.
        ELSE
          IF (ndims_of_var == 4) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has four dimensions, but the fourth one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have lon,lat,m as dimensions
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename)

        ! Inquire the time dimension
        CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have lon,lat,m as dimensions

        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lonlat_grid_field_dp_2D_monthly

  SUBROUTINE check_lonlat_grid_field_dp_3D(         filename, var_name, should_have_time)
    ! Check if this file contains a 3-D lon/lat-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lonlat_grid_field_dp_3D'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_zeta, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid lon and lat dimensions and variables
    CALL check_lon(  filename)
    CALL check_lat(  filename)
    CALL check_zeta( filename)

    ! Inquire lon,lat dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_lon , id_dim_lon )
    CALL inquire_dim_multiple_options( filename, field_name_options_lat , id_dim_lat )
    CALL inquire_dim_multiple_options( filename, field_name_options_zeta, id_dim_zeta)

    ! Inquire variable
    CALL inquire_var( filename, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check lon,lat dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_lon )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have longitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_lat )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have latitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_zeta)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have zeta as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 3) THEN
          ! The variable only has lon,lat,zeta as dimensions.
        ELSE
          IF (ndims_of_var == 4) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has four dimensions, but the fourth one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have lon,lat,zeta as dimensions
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename)

        ! Inquire the time dimension
        CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have lon,lat,zeta as dimensions

        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lonlat_grid_field_dp_3D

  SUBROUTINE check_lonlat_grid_field_dp_ocean_3D(         filename, var_name, should_have_time)
    ! Check if this file contains a 3-D lon/lat-grid variable by this name

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: should_have_time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_lonlat_grid_field_dp_ocean_3D'
    INTEGER                                            :: id_dim_lon, id_dim_lat, id_dim_z_ocean, id_dim_time, id_var
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS)             :: dims_of_var
    LOGICAL                                            :: file_has_time

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if the file has valid lon and lat dimensions and variables
    CALL check_lon(     filename)
    CALL check_lat(     filename)
    CALL check_z_ocean( filename)

    ! Inquire lon,lat dimensions
    CALL inquire_dim_multiple_options( filename, field_name_options_lon , id_dim_lon )
    CALL inquire_dim_multiple_options( filename, field_name_options_lat , id_dim_lat )
    CALL inquire_dim_multiple_options( filename, field_name_options_z_ocean, id_dim_z_ocean)

    ! Inquire variable
    CALL inquire_var( filename, var_name, id_var)
    IF (id_var == -1) CALL crash('variable "' // TRIM( var_name) // '" could not be found in file "' // TRIM( filename) // '"!')

    ! Inquire variable info
    CALL inquire_var_info( filename, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) THEN
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')
    END IF

    ! Check lon,lat dimensions
    IF (.NOT. ANY( dims_of_var == id_dim_lon )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have longitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_lat )) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have latitude as a dimension!')
    IF (.NOT. ANY( dims_of_var == id_dim_z_ocean)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have z_ocean as a dimension!')

    IF (.NOT. PRESENT( should_have_time)) THEN
      ! This variable is allowed to either have or not have a time dimension

      ! Check if the file contains a time dimension
      CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)
      IF (id_dim_time == -1) THEN
        file_has_time = .FALSE.
      ELSE
        file_has_time = .TRUE.
      END IF

      IF (file_has_time) THEN
        ! Check if the variable has time as a dimension
        IF (ndims_of_var == 3) THEN
          ! The variable only has lon,lat,z_ocean as dimensions.
        ELSE
          IF (ndims_of_var == 4) THEN
            IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' &
              // TRIM( filename) // '" has four dimensions, but the fourth one is not time!')
          ELSE
            CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
          END IF
        END IF
      ELSE ! IF (file_has_time) THEN
        ! The file does not have a time dimension; the variable should only have lon,lat,z_ocean as dimensions
        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
      END IF ! IF (file_has_time) THEN

    ELSE ! IF (.NOT. PRESENT( should_have_time)) THEN
      IF (should_have_time) THEN
        ! This variable should have a time dimension

        ! Check if the file has a valid time dimension
        CALL check_time( filename)

        ! Inquire the time dimension
        CALL inquire_dim_multiple_options( filename, field_name_options_time, id_dim_time)

        ! Check if the variable has time as a dimension
        IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
        IF (.NOT. ANY( dims_of_var == id_dim_time)) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" does not have time as a dimension!')

      ELSE ! IF (should_have_time) THEN
        ! This variable should not have a time dimension; the variable should only have lon,lat,z_ocean as dimensions

        IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

      END IF ! IF (should_have_time) THEN
    END IF ! IF (.NOT. PRESENT( should_have_time)) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_lonlat_grid_field_dp_ocean_3D


  ! ===== Flexible looking for dimensions and variables =====
  ! =========================================================

  ! Look for dimensions
  SUBROUTINE inquire_dim_multiple_options( filename, dim_name_options, id_dim, dim_length, dim_name)
    ! Inquire if this file contains a dimension by name of dim_name.
    ! If so, return its length and identifier. If not, return -1 for both.
    !
    ! Supports providing multiple options for the dimension name, separated by two
    ! vertical bars || e.g. if we're looking for an X-dimension, we could do something like:
    !
    ! CALL inquire_dim_multiple_options( ncid, dim_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', dim_length, id_dim)
    !
    ! IF more than one match is found, crash.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: dim_name_options
    INTEGER,                             INTENT(OUT)   :: id_dim

    INTEGER,                                INTENT(OUT), OPTIONAL :: dim_length
    CHARACTER(LEN=256),                     INTENT(OUT), OPTIONAL :: dim_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_dim_multiple_options'
    CHARACTER(LEN=256)                                 :: dim_name_options_parsed
    CHARACTER(LEN=256)                                 :: dim_name_options_redux
    INTEGER                                            :: i, n_matches
    INTEGER                                            :: dim_length_try, dim_length_match
    INTEGER                                            :: id_dim_try, id_dim_match
    CHARACTER(LEN=256)                                 :: dim_name_try, dim_name_match

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Parse field name options
    CALL parse_field_name_options( dim_name_options, dim_name_options_parsed)

    ! Try all options provided in dim_name_options

    dim_name_options_redux = TRIM( dim_name_options_parsed)
    n_matches = 0

    DO WHILE (.TRUE.)

      i = INDEX( dim_name_options_redux, '||')

      IF (i > 0) THEN
        ! More than one option is left over; take the last one

        dim_name_try = dim_name_options_redux( 1:i-1)
        dim_name_options_redux = dim_name_options_redux( i+2:LEN_TRIM( dim_name_options_redux))

      ELSE
        ! Only one option is left over

        dim_name_try = dim_name_options_redux
        dim_name_options_redux( 1:LEN( dim_name_options_redux)) = ''

      END IF

      ! Try the selected name option
      CALL inquire_dim( filename, dim_name_try, dim_length_try, id_dim_try)

      IF (id_dim_try == -1) THEN
        ! No dimension by this name was found; try the next option
      ELSE
        ! A dimension by this name was found; hurray!
        n_matches  = n_matches + 1
        dim_length_match = dim_length_try
        id_dim_match     = id_dim_try
        dim_name_match   = dim_name_try
      END IF

      ! If the list of options is now empty, exit
      IF (LEN_TRIM( dim_name_options_redux) == 0) EXIT

    END DO

    IF (n_matches == 0) THEN
      ! None of the optional dimension names were found in the NetCDF file
      dim_length_match = -1
      id_dim_match     = -1
    ELSEIF (n_matches > 1) THEN
      ! More than one match was found
      CALL crash('more than one of the provided dimension names were found in file "' // TRIM( filename) // '"!')
    ELSE
      ! We found exactly one match; hurray!
    END IF

    ! Copy to output arguments
    id_dim = id_dim_match
    IF (PRESENT( dim_name  )) dim_name   = dim_name_match
    IF (PRESENT( dim_length)) dim_length = dim_length_match

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_dim_multiple_options

  ! Look for variables
  SUBROUTINE inquire_var_multiple_options( filename, var_name_options, id_var, var_name, var_type, ndims_of_var, dims_of_var)
    ! Inquire if this file contains a variable by name of var_name.
    ! If so, return its identifier. If not, return -1.
    !
    ! Supports providing multiple options for the variable name, separated by two
    ! vertical bars || e.g. if we're looking for an X-variable, we could do something like:
    !
    ! CALL inquire_var_multiple_options( ncid, var_name_options = 'x||X||x1||X1||x-coordinate||X-coordinate||easting', id_var)
    !
    ! IF more than one match is found, crash.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name_options
    INTEGER,                             INTENT(OUT)   :: id_var

    CHARACTER(LEN=256),                     INTENT(OUT), OPTIONAL :: var_name
    INTEGER,                                INTENT(OUT), OPTIONAL :: var_type
    INTEGER,                                INTENT(OUT), OPTIONAL :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS), INTENT(OUT), OPTIONAL :: dims_of_var

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_var_multiple_options'
    CHARACTER(LEN=256)                                 :: var_name_options_parsed
    CHARACTER(LEN=256)                                 :: var_name_options_redux
    INTEGER                                            :: i, n_matches, id_var_try
    CHARACTER(LEN=256)                                 :: var_name_try, var_name_match

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Parse field name options
    CALL parse_field_name_options( var_name_options, var_name_options_parsed)

    ! Try all options provided in var_name_options
    var_name_options_redux = TRIM( var_name_options_parsed)
    n_matches = 0

    DO WHILE (.TRUE.)

      i = INDEX( var_name_options_redux, '||')

      IF (i > 0) THEN
        ! More than one option is left over; take the last one

        var_name_try = var_name_options_redux( 1:i-1)
        var_name_options_redux = var_name_options_redux( i+2:LEN_TRIM( var_name_options_redux))

      ELSE
        ! Only one option is left over

        var_name_try = TRIM( var_name_options_redux)
        var_name_options_redux( 1:LEN( var_name_options_redux)) = ''

      END IF

      ! Try the selected name option
      CALL inquire_var( filename, var_name_try, id_var_try)

      IF (id_var_try == -1) THEN
        ! No variable by this name was found; try the next option
      ELSE
        ! A variable by this name was found; hurray!
        n_matches      = n_matches + 1
        id_var         = id_var_try
        var_name_match = var_name_try
      END IF

      ! If the list of options is now empty, exit
      IF (LEN_TRIM( var_name_options_redux) == 0) EXIT

    END DO

    IF (n_matches == 0) THEN
      ! None of the optional variable names were found in the NetCDF file
      id_var     = -1
    ELSEIF (n_matches > 1) THEN
      ! More than one match was found
      CALL crash('more than one of the provided variable names were found in file "' // TRIM( filename) // '"!')
    ELSE
      ! We found exactly one match. Inquire additional info on this variable.
      CALL inquire_var_info( filename, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    END IF

    ! Copy to output arguments
    IF (PRESENT( var_name)) var_name = var_name_match

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_var_multiple_options

  ! ===== Parse flexible dimension/variable names =====
  ! ===================================================

  SUBROUTINE parse_field_name_options( field_name_options, field_name_options_parsed)
    ! Check if a default set of field name options should be used.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=256),                  INTENT(OUT)   :: field_name_options_parsed

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'parse_field_name_options'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    field_name_options_parsed = field_name_options

    IF (INDEX( field_name_options,'default_options_') > 0) THEN
      ! Use one of the default options

      ! Dimensions
      IF     (field_name_options == 'default_options_x') THEN
        field_name_options_parsed = field_name_options_x
      ELSEIF (field_name_options == 'default_options_y') THEN
        field_name_options_parsed = field_name_options_y
      ELSEIF (field_name_options == 'default_options_zeta') THEN
        field_name_options_parsed = field_name_options_zeta
      ELSEIF (field_name_options == 'default_options_z_ocean') THEN
        field_name_options_parsed = field_name_options_z_ocean
      ELSEIF (field_name_options == 'default_options_lon') THEN
        field_name_options_parsed = field_name_options_lon
      ELSEIF (field_name_options == 'default_options_lat') THEN
        field_name_options_parsed = field_name_options_lat
      ELSEIF (field_name_options == 'default_options_time') THEN
        field_name_options_parsed = field_name_options_time

      ! Variables
      ELSEIF (field_name_options == 'default_options_Hi') THEN
        field_name_options_parsed = field_name_options_Hi
      ELSEIF (field_name_options == 'default_options_Hb') THEN
        field_name_options_parsed = field_name_options_Hb
      ELSEIF (field_name_options == 'default_options_Hs') THEN
        field_name_options_parsed = field_name_options_Hs
      ELSEIF (field_name_options == 'default_options_u_surf') THEN
        field_name_options_parsed = field_name_options_u_surf
      ELSEIF (field_name_options == 'default_options_v_surf') THEN
        field_name_options_parsed = field_name_options_v_surf
      ELSEIF (field_name_options == 'default_options_T_ocean') THEN
        field_name_options_parsed = field_name_options_T_ocean
      ELSEIF (field_name_options == 'default_options_S_ocean') THEN
        field_name_options_parsed = field_name_options_S_ocean
      ! Unrecognised default options
      ELSE
        CALL crash('unregocnised default field name option "' // TRIM( field_name_options) // '"')
      END IF

    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE parse_field_name_options

  FUNCTION get_first_option_from_list( field_name_options) RESULT( field_name)
    ! Get the first option from a list of field name options

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: field_name_options
    CHARACTER(LEN=256)                                 :: field_name

    ! Local variables:
    INTEGER                                            :: i

    field_name( 1:256) = ' '

    i = INDEX( field_name_options,'||')

    IF (i > 0) THEN
      field_name = field_name_options( 1:i-1)
    ELSE
      field_name = TRIM( field_name_options)
    END IF

  END FUNCTION get_first_option_from_list

  ! ===== Read data from variables =====
  ! ====================================

  SUBROUTINE read_var_int_0D(  filename, id_var, d)
    ! Read data from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,                             INTENT(OUT)   :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_int_0D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_int_0D

  SUBROUTINE read_var_int_1D(  filename, id_var, d, start, count)
    ! Read data from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:    ),          INTENT(OUT)   :: d
    INTEGER,  DIMENSION(1    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_int_1D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 1)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT .OR. var_type == NF90_INT64)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_int_1D

  SUBROUTINE read_var_int_2D(  filename, id_var, d, start, count)
    ! Read data from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:,:  ),          INTENT(OUT)   :: d
    INTEGER,  DIMENSION(2    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_int_2D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 2)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1), SIZE( d,2) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_int_2D

  SUBROUTINE read_var_int_3D(  filename, id_var, d, start, count)
    ! Read data from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:,:,:),          INTENT(OUT)   :: d
    INTEGER,  DIMENSION(3    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_int_3D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 3)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1), SIZE( d,2), SIZE( d,3) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_int_3D

  SUBROUTINE read_var_int_4D(  filename, id_var, d, start, count)
    ! Read data from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:,:,:,:),        INTENT(OUT)   :: d
    INTEGER,  DIMENSION(4    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_int_4D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 4)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1, 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1), SIZE( d,2), SIZE( d,3), SIZE( d,4) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_int_4D

  SUBROUTINE read_var_dp_0D(  filename, id_var, d)
    ! Read data from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp),                            INTENT(OUT)   :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_dp_0D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_dp_0D

  SUBROUTINE read_var_dp_1D(  filename, id_var, d, start, count)
    ! Read data from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d
    INTEGER,  DIMENSION(1    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_dp_1D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 1)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_dp_1D

  SUBROUTINE read_var_dp_2D(  filename, id_var, d, start, count)
    ! Read data from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d
    INTEGER,  DIMENSION(2    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_dp_2D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 2)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1), SIZE( d,2) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_dp_2D

  SUBROUTINE read_var_dp_3D(  filename, id_var, d, start, count)
    ! Read data from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:,:,:),          INTENT(OUT)   :: d
    INTEGER,  DIMENSION(3    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_dp_3D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 3)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1), SIZE( d,2), SIZE( d,3) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '": start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_dp_3D

  SUBROUTINE read_var_dp_4D(  filename, id_var, d, start, count)
    ! Read data from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:,:,:,:),        INTENT(OUT)   :: d
    INTEGER,  DIMENSION(4    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_var_dp_4D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 4)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1, 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1), SIZE( d,2), SIZE( d,3), SIZE( d,4) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    ! Read the data
    IF (par%master) THEN
      nerr = NF90_GET_VAR( ncid, id_var, d, start, count)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_GET_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_var_dp_4D

  ! ===== Write data to variables =====
  ! ===================================

  SUBROUTINE write_var_int_0D(  filename, id_var, d)
    ! Write data to a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,                             INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_int_0D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Open the file
    CALL open_existing_netcdf_file_in_data_mode( filename, ncid)

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_int_0D

  SUBROUTINE write_var_int_1D(  filename, id_var, d, start, count)
    ! Write data to a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: d
    INTEGER,  DIMENSION(1    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_int_1D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 1)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_in_data_mode( filename, ncid)

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_int_1D

  SUBROUTINE write_var_int_2D(  filename, id_var, d, start, count)
    ! Write data to a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: d
    INTEGER,  DIMENSION(2    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_int_2D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 2)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1), SIZE( d,2) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_in_data_mode( filename, ncid)

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_int_2D

  SUBROUTINE write_var_int_3D(  filename, id_var, d, start, count)
    ! Write data to a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:,:,:),          INTENT(IN)    :: d
    INTEGER,  DIMENSION(3    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_int_3D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 3)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1), SIZE( d,2), SIZE( d,3) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_in_data_mode( filename, ncid)

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_int_3D

  SUBROUTINE write_var_int_4D(  filename, id_var, d, start, count)
    ! Write data to a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    INTEGER,  DIMENSION(:,:,:,:),        INTENT(IN)    :: d
    INTEGER,  DIMENSION(4    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_int_4D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 4)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_INT)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_INT!')

    ! Check number of dimensions
    IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1, 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1), SIZE( d,2), SIZE( d,3), SIZE( d,4) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_in_data_mode( filename, ncid)

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_int_4D

  SUBROUTINE write_var_dp_0D(  filename, id_var, d)
    ! Write data to a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp),                            INTENT(IN)    :: d

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_dp_0D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 0) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Open the file
    CALL open_existing_netcdf_file_in_data_mode( filename, ncid)

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_dp_0D

  SUBROUTINE write_var_dp_1D(  filename, id_var, d, start, count)
    ! Write data to a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d
    INTEGER,  DIMENSION(1    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_dp_1D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 1)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 1) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (var_name /= 'time') THEN
        ! Exception for time, because there the dimension is usually unlimited
        IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
          TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)
      END IF

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_in_data_mode( filename, ncid)

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_dp_1D

  SUBROUTINE write_var_dp_2D(  filename, id_var, d, start, count)
    ! Write data to a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    INTEGER,  DIMENSION(2    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_dp_2D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 2)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 2) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1), SIZE( d,2) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_in_data_mode( filename, ncid)

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_dp_2D

  SUBROUTINE write_var_dp_3D(  filename, id_var, d, start, count)
    ! Write data to a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:,:,:),          INTENT(IN)    :: d
    INTEGER,  DIMENSION(3    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_dp_3D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 3)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 3) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1), SIZE( d,2), SIZE( d,3) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_in_data_mode( filename, ncid)

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_dp_3D

  SUBROUTINE write_var_dp_4D(  filename, id_var, d, start, count)
    ! Write data to a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    REAL(dp), DIMENSION(:,:,:,:),        INTENT(IN)    :: d
    INTEGER,  DIMENSION(4    ), OPTIONAL,INTENT(IN)    :: start, count

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_var_dp_4D'
    INTEGER                                            :: ncid
    CHARACTER(LEN=256)                                 :: var_name
    INTEGER                                            :: var_type
    INTEGER                                            :: ndims_of_var
    INTEGER,  DIMENSION( NF90_MAX_VAR_DIMS)            :: dims_of_var
    INTEGER                                            :: di
    CHARACTER(LEN=256)                                 :: dim_name
    INTEGER                                            :: dim_length
    INTEGER,  DIMENSION( 4)                            :: start_applied, count_applied

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Inquire some info on this variable
    CALL inquire_var_info( filename, id_var, var_name = var_name, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Check variable type
    IF (.NOT. (var_type == NF90_FLOAT .OR. var_type == NF90_DOUBLE)) &
      CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" is not of type NF90_FLOAT or NF90_DOUBLE!')

    ! Check number of dimensions
    IF (ndims_of_var /= 4) CALL crash('variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)

    ! Set start and count
    IF (PRESENT( start)) THEN
      start_applied = start
    ELSE
      start_applied = (/ 1, 1, 1, 1 /)
    END IF
    IF (ANY( start_applied == 0)) CALL crash('start must be positive!')

    IF (PRESENT( count)) THEN
      count_applied = count
    ELSE
      count_applied = (/ SIZE( d,1), SIZE( d,2), SIZE( d,3), SIZE( d,4) /)
    END IF
    IF (ANY( count_applied == 0)) CALL crash('count must be positive!')

    ! Check sizes of dimensions
    DO di = 1, ndims_of_var

      ! Check size of this dimension in the file
      CALL inquire_dim_info( filename, dims_of_var( di), dim_name = dim_name, dim_length = dim_length)

      ! Check if the combination of dimension size, start, and count, matches the size of d
      IF (count_applied( di) /= SIZE( d,di)) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // &
        '": count({int_01}) = {int_02}, but SIZE(d,{int_03}) = {int_04}!', int_01 = di, int_02 = count_applied( di), int_03 = di, int_04 = SIZE( d,di))

      ! Check if this dimension is large enough to read this amount of data
      IF (start_applied( di) + count_applied( di) - 1 > dim_length) CALL crash('error for dimension "' // TRIM( dim_name) // '" of variable "' // TRIM( var_name) // '" in file "' // &
        TRIM( filename) // '"start + count - 1 = {int_01}, but dim_length = {int_02}!', int_01 = start_applied( di) + count_applied( di) - 1, int_02 = dim_length)

    END DO

    ! Open the file
    CALL open_existing_netcdf_file_in_data_mode( filename, ncid)

    ! Write the data
    IF (par%master) THEN
      nerr = NF90_PUT_VAR( ncid, id_var, d, start_applied, count_applied)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_VAR failed for variable "' // TRIM( var_name) // '" in file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_var_dp_4D

  ! ===== Basic NetCDF wrapper functions =====
  ! ==========================================

  ! Open and close a NetCDF file
  SUBROUTINE open_existing_netcdf_file_for_reading( filename, ncid)
    ! Open the NetCDF file in the specified location for reading only,
    ! and return its identifier.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(OUT)   :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'open_netcdf_file_for_reading'
    LOGICAL                                            :: file_exists

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Open the NetCDF file with read-only access
    IF (par%master) THEN
      nerr = NF90_OPEN( TRIM( filename), NF90_NOWRITE, ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_OPEN failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN

    CALL MPI_BCAST( ncid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE open_existing_netcdf_file_for_reading

  SUBROUTINE open_existing_netcdf_file_in_define_mode( filename, ncid)
    ! Open an existing NetCDF file in define mode
    ! In define mode, new dimensions, variables, or attributes can be created,
    ! but no data can be written to existing variables.
    ! When opening an existing NetCDF file, it is by default in data mode.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(OUT)   :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'open_existing_netcdf_file_in_define_mode'
    LOGICAL                                            :: file_exists

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Open the NetCDF file with read+write access
    IF (par%master) THEN
      nerr = NF90_OPEN( TRIM( filename), NF90_WRITE, ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_OPEN failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN

    CALL MPI_BCAST( ncid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Put the NetCDF file in define mode
    IF (par%master) THEN
      nerr = NF90_REDEF( ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_REDEF failed for file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE open_existing_netcdf_file_in_define_mode

  SUBROUTINE open_existing_netcdf_file_in_data_mode( filename, ncid)
    ! Open an existing NetCDF file in data mode
    ! In data mode, no new dimensions, variables, or attributes can be created,
    ! but data can be written to existing variables.
    ! When opening an existing NetCDF file, it is by default in data mode.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(OUT)   :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'open_existing_netcdf_file_in_data_mode'
    LOGICAL                                            :: file_exists

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if this file actually exists
    INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
    IF (.NOT. file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" not found!')
    END IF

    ! Open the NetCDF file with read+write access
    IF (par%master) THEN
      nerr = NF90_OPEN( TRIM( filename), NF90_WRITE, ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_OPEN failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN

    CALL MPI_BCAST( ncid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE open_existing_netcdf_file_in_data_mode

  SUBROUTINE close_netcdf_file( ncid)
    ! Close an opened NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,                             INTENT(OUT)   :: ncid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'close_netcdf_file'

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Close netCDF file:
    IF (par%master) THEN
      nerr = NF90_CLOSE( ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_CLOSE failed!')
    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE close_netcdf_file

  ! Inquire dimensions and variables
  SUBROUTINE inquire_dim( filename, dim_name, dim_length, id_dim)
    ! Inquire if this file contains a dimension by name of dim_name.
    ! If so, return its length and identifier; if not, return -1 for both.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: dim_name
    INTEGER,                             INTENT(OUT)   :: dim_length
    INTEGER,                             INTENT(OUT)   :: id_dim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_dim'
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    IF (par%master) THEN

      ! Check if a dimension of this name exists in the file
      nerr = NF90_INQ_DIMID( ncid, dim_name, id_dim)

      IF (nerr /= NF90_NOERR) THEN
        ! If a dimension by this name does not exist, return -1 for the length and ID
        id_dim     = -1
        dim_length = -1
      ELSE
        ! If a dimension by this name exists, find its length
        nerr = NF90_INQUIRE_DIMENSION( ncid, id_dim, len = dim_length)
        IF (nerr /= NF90_NOERR) CALL crash('NF90_INQUIRE_DIMENSION failed for file "' // TRIM( filename) // '"!')
      END IF

    END IF ! IF (par%master) THEN
    CALL sync

    CALL MPI_BCAST( id_dim    , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( dim_length, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_dim

  SUBROUTINE inquire_dim_info( filename, id_dim, dim_name, dim_length)
    ! Inquire some info of a dimension

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_dim

    CHARACTER(LEN=256),                     INTENT(OUT), OPTIONAL :: dim_name
    INTEGER,                                INTENT(OUT), OPTIONAL :: dim_length

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_dim_info'
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    IF (par%master) THEN
      ! Inquire some info on this variable
      nerr = NF90_INQUIRE_DIMENSION( ncid, id_dim, name = dim_name, len = dim_length)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_INQUIRE_DIMENSION failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN

    IF (PRESENT( dim_name  )) CALL MPI_BCAST( dim_name  , 256, MPI_CHAR   , 0, MPI_COMM_WORLD, ierr)
    IF (PRESENT( dim_length)) CALL MPI_BCAST( dim_length, 1  , MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_dim_info

  SUBROUTINE inquire_var( filename, var_name, id_var)
    ! Inquire if this file contains a variable by name of var_name.
    ! If so, return its identifier. If not, return -1 for the identifier.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    INTEGER,                             INTENT(OUT)   :: id_var

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_var'
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    IF (par%master) THEN

      ! Check if a variable of this name exists in the file
      nerr = NF90_INQ_VARID( ncid, var_name, id_var)

      IF (nerr /= NF90_NOERR) THEN
        ! If a variable by this name does not exist, return -1 for the ID
        id_var = -1
      END IF

    END IF ! IF (par%master) THEN
    CALL sync

    CALL MPI_BCAST( id_var, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_var

  SUBROUTINE inquire_var_info( filename, id_var, var_name, var_type, ndims_of_var, dims_of_var)
    ! Inquire some info of a variable

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var

    CHARACTER(LEN=256),                     INTENT(OUT), OPTIONAL :: var_name
    INTEGER,                                INTENT(OUT), OPTIONAL :: var_type
    INTEGER,                                INTENT(OUT), OPTIONAL :: ndims_of_var
    INTEGER, DIMENSION( NF90_MAX_VAR_DIMS), INTENT(OUT), OPTIONAL :: dims_of_var

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_var_info'
    INTEGER                                            :: ncid
    
    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Open the file
    CALL open_existing_netcdf_file_for_reading( filename, ncid)

    IF (par%master) THEN
      ! Inquire some info on this variable
      nerr = NF90_INQUIRE_VARIABLE( ncid, id_var, name = var_name, xtype = var_type, ndims = ndims_of_var, dimids = dims_of_var)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_INQUIRE_VARIABLE failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN

    IF (PRESENT( var_name    )) CALL MPI_BCAST( var_name    , 256              , MPI_CHAR   , 0, MPI_COMM_WORLD, ierr)
    IF (PRESENT( var_type    )) CALL MPI_BCAST( var_type    , 1                , MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    IF (PRESENT( ndims_of_var)) CALL MPI_BCAST( ndims_of_var, 1                , MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    IF (PRESENT(  dims_of_var)) CALL MPI_BCAST( dims_of_var , NF90_MAX_VAR_DIMS, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_var_info

  ! Create new NetCDF file
  SUBROUTINE create_new_netcdf_file_for_writing( filename)
    ! Create a new NetCDF file in the specified location for writing.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_new_netcdf_file_for_writing'
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Check if this file already exists
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL crash('file "' // TRIM( filename) // '" already exists!')
      END IF
    END IF

    ! Create the NetCDF file
    IF (par%master) THEN
      nerr = NF90_CREATE( filename, IOR( NF90_NOCLOBBER, NF90_NETCDF4), ncid)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_CREATE failed for file "' // TRIM( filename) // '"!')
    END IF ! IF (par%master) THEN
    CALL MPI_BCAST( ncid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_new_netcdf_file_for_writing

  ! Create dimensions, variables, and attributes
  SUBROUTINE create_dimension( filename, dim_name, dim_length, id_dim)
    ! Create a new dimension in a NetCDF file.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: dim_name
    INTEGER,                             INTENT(IN)    :: dim_length
    INTEGER,                             INTENT(OUT)   :: id_dim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_dimension'
    INTEGER                                            :: dim_length_present, ncid

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety: check if a dimension by this name is already present in this file
    CALL inquire_dim( filename, dim_name, dim_length_present, id_dim)
    IF (id_dim /= -1) CALL crash('file "' // TRIM( filename) // '" already contains dimension "' // TRIM( dim_name) // '"!')

    ! Open the NetCDF file in define mode
    CALL open_existing_netcdf_file_in_define_mode( filename, ncid)

    ! Add the dimension
    IF (par%master) THEN
      nerr = NF90_DEF_DIM( ncid, dim_name, dim_length, id_dim)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_DEF_DIM failed for file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    CALL MPI_BCAST( id_dim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_dimension

  SUBROUTINE create_variable( filename, var_name, var_type, dim_ids, id_var)
    ! Create a new variable in a NetCDF file.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    CHARACTER(LEN=*),                    INTENT(IN)    :: var_name
    INTEGER,                             INTENT(IN)    :: var_type
    INTEGER, DIMENSION(:),               INTENT(IN)    :: dim_ids
    INTEGER,                             INTENT(OUT)   :: id_var

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_variable'
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Safety: check if a variable by this name is already present in this file
    CALL inquire_var( filename, var_name, id_var)
    IF (id_var /= -1) CALL crash('file "' // TRIM( filename) // '" already contains variable "' // TRIM( var_name) // '"!')

    ! Open the NetCDF file in define mode
    CALL open_existing_netcdf_file_in_define_mode( filename, ncid)

    ! Add the variable
    IF (par%master) THEN
      nerr = NF90_DEF_VAR( ncid, name = var_name, xtype = var_type, dimids = dim_ids, varid = id_var)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_DEF_VAR failed for file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    CALL MPI_BCAST( id_var, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_variable

  SUBROUTINE add_attribute_int( filename, id_var, att_name, att_val)
    ! Add an integer-valued attributes to a variable.
    ! Assume the file is in data mode; put it in define mode,
    ! add the attribute, and put it back in data mode.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    CHARACTER(LEN=*),                    INTENT(IN)    :: att_name
    INTEGER,                             INTENT(IN)    :: att_val

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_attribute_int'
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Open the NetCDF file in define mode
    CALL open_existing_netcdf_file_in_define_mode( filename, ncid)

    ! Add the attribute
    IF (par%master) THEN
      nerr = NF90_PUT_ATT( ncid, id_var, att_name, att_val)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_ATT failed for file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_attribute_int

  SUBROUTINE add_attribute_dp( filename, id_var, att_name, att_val)
    ! Add a double-precision-valued attributes to a variable.
    ! Assume the file is in data mode; put it in define mode,
    ! add the attribute, and put it back in data mode.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    CHARACTER(LEN=*),                    INTENT(IN)    :: att_name
    REAL(dp),                            INTENT(IN)    :: att_val

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_attribute_dp'
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Open the NetCDF file in define mode
    CALL open_existing_netcdf_file_in_define_mode( filename, ncid)

    ! Add the attribute
    IF (par%master) THEN
      nerr = NF90_PUT_ATT( ncid, id_var, att_name, att_val)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_ATT failed for file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_attribute_dp

  SUBROUTINE add_attribute_char( filename, id_var, att_name, att_val)
    ! Add a character-valued attributes to a variable.
    ! Assume the file is in data mode; put it in define mode,
    ! add the attribute, and put it back in data mode.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename
    INTEGER,                             INTENT(IN)    :: id_var
    CHARACTER(LEN=*),                    INTENT(IN)    :: att_name
    CHARACTER(LEN=*),                    INTENT(IN)    :: att_val

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'add_attribute_char'
    INTEGER                                            :: ncid

    ! Add routine to path
    CALL init_routine( routine_name, do_track_resource_use = .FALSE.)

    ! Open the NetCDF file in define mode
    CALL open_existing_netcdf_file_in_define_mode( filename, ncid)

    ! Add the attribute
    IF (par%master) THEN
      nerr = NF90_PUT_ATT( ncid, id_var, att_name, att_val)
      IF (nerr /= NF90_NOERR) CALL crash('NF90_PUT_ATT failed for file "' // TRIM( filename) // '"!')
    END IF
    CALL sync

    ! Close the NetCDF file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE add_attribute_char

END MODULE netcdf_basic_module
