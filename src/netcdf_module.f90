MODULE netcdf_module

  ! Contains routines for creating, reading, and writing all the NetCDF files
  ! involved in the ice model.

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE netcdf,                          ONLY: nf90_max_var_dims, nf90_create, nf90_close, nf90_clobber, nf90_share, nf90_unlimited , &
                                             nf90_enddef, nf90_put_var, nf90_sync, nf90_def_var, nf90_int, nf90_put_att, nf90_def_dim, &
                                             nf90_open, nf90_write, nf90_inq_dimid, nf90_inquire_dimension, nf90_inquire, nf90_double, &
                                             nf90_inq_varid, nf90_inquire_variable, nf90_get_var, nf90_noerr, nf90_strerror, nf90_float
  USE data_types_module,               ONLY: type_grid, type_grid_lonlat, type_ice_model, type_model_region, type_reference_geometry, &
                                             type_restart_data, type_forcing_data, type_climate_snapshot, type_ocean_snapshot_global, &
                                             type_SELEN_global, type_global_scalar_data, type_highres_ocean_data, &
                                             type_BIV_target_velocity, type_BIV_bed_roughness, type_sparse_matrix_CSR
  USE data_types_netcdf_module

  IMPLICIT NONE


CONTAINS


  SUBROUTINE write_to_resource_tracking_file( netcdf, time)
    ! Write to the resource tracking output file

    USE configuration_module, ONLY: resource_tracker

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_resource_tracker), INTENT(INOUT) :: netcdf
    REAL(dp),                           INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'write_to_resource_tracking_file'
    INTEGER                                           :: i,n
    INTEGER,  DIMENSION(1024)                         :: path_int_enc

    IF (.NOT. par%master) RETURN

    ! Open the file for writing
    CALL open_netcdf_file( netcdf%filename, netcdf%ncid)

    ! Time
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_time, time, start = (/netcdf%ti/)))

    ! Actual variables
    ! ================

    n = SIZE( resource_tracker)

    DO i = 1, n

      ! Subroutine name
      CALL encode_subroutine_path_as_integer( resource_tracker( i)%routine_path, path_int_enc)
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_names( i), path_int_enc ))

      ! Computation time
      CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_vars( i), resource_tracker( i)%tcomp, start = (/ netcdf%ti /) ))

    END DO

    ! Close the file
    CALL close_netcdf_file( netcdf%ncid)

    ! Increase time frame counter
    netcdf%ti = netcdf%ti + 1

  END SUBROUTINE write_to_resource_tracking_file
  SUBROUTINE encode_subroutine_path_as_integer( subroutine_path, path_int_enc)
    ! Encode the current subroutine path as an integer array so it can be saved as a NetCDF variable
    !
    ! Use the simplest possible encoding:
    !
    !  ' ' = -1 (empty character)
    !
    !    0 = 0
    !    1 = 1
    !    ...
    !    9 = 9
    !
    !    a = 10
    !    b = 11
    !    c = 12
    !    ...
    !    z = 36
    !
    !    A = 37
    !    B = 38
    !    C = 39
    !    ...
    !    Z = 62
    !
    !    _ = 63 (underscore)
    !    / = 64 (forward slash)
    !    ( = 65 (left  bracket)
    !    ) = 66 (right bracket)

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=1024),                INTENT(IN)    :: subroutine_path
    INTEGER,  DIMENSION(1024),          INTENT(OUT)   :: path_int_enc

    ! Local variables:
    INTEGER                                           :: i

    path_int_enc = 0

    DO i = 1, 1024

      SELECT CASE ( subroutine_path( i:i))
      CASE( ' ')
        path_int_enc( i) = -1
      CASE( '0')
        path_int_enc( i) = 0
      CASE( '1')
        path_int_enc( i) = 1
      CASE( '2')
        path_int_enc( i) = 2
      CASE( '3')
        path_int_enc( i) = 3
      CASE( '4')
        path_int_enc( i) = 4
      CASE( '5')
        path_int_enc( i) = 5
      CASE( '6')
        path_int_enc( i) = 6
      CASE( '7')
        path_int_enc( i) = 7
      CASE( '8')
        path_int_enc( i) = 8
      CASE( '9')
        path_int_enc( i) = 9
      CASE( 'a')
        path_int_enc( i) = 11
      CASE( 'b')
        path_int_enc( i) = 12
      CASE( 'c')
        path_int_enc( i) = 13
      CASE( 'd')
        path_int_enc( i) = 14
      CASE( 'e')
        path_int_enc( i) = 15
      CASE( 'f')
        path_int_enc( i) = 16
      CASE( 'g')
        path_int_enc( i) = 17
      CASE( 'h')
        path_int_enc( i) = 18
      CASE( 'i')
        path_int_enc( i) = 19
      CASE( 'j')
        path_int_enc( i) = 20
      CASE( 'k')
        path_int_enc( i) = 21
      CASE( 'l')
        path_int_enc( i) = 22
      CASE( 'm')
        path_int_enc( i) = 23
      CASE( 'n')
        path_int_enc( i) = 24
      CASE( 'o')
        path_int_enc( i) = 25
      CASE( 'p')
        path_int_enc( i) = 26
      CASE( 'q')
        path_int_enc( i) = 27
      CASE( 'r')
        path_int_enc( i) = 28
      CASE( 's')
        path_int_enc( i) = 29
      CASE( 't')
        path_int_enc( i) = 30
      CASE( 'u')
        path_int_enc( i) = 31
      CASE( 'v')
        path_int_enc( i) = 32
      CASE( 'w')
        path_int_enc( i) = 33
      CASE( 'x')
        path_int_enc( i) = 34
      CASE( 'y')
        path_int_enc( i) = 35
      CASE( 'z')
        path_int_enc( i) = 36
      CASE( 'A')
        path_int_enc( i) = 37
      CASE( 'B')
        path_int_enc( i) = 38
      CASE( 'C')
        path_int_enc( i) = 39
      CASE( 'D')
        path_int_enc( i) = 40
      CASE( 'E')
        path_int_enc( i) = 41
      CASE( 'F')
        path_int_enc( i) = 42
      CASE( 'G')
        path_int_enc( i) = 43
      CASE( 'H')
        path_int_enc( i) = 44
      CASE( 'I')
        path_int_enc( i) = 45
      CASE( 'J')
        path_int_enc( i) = 46
      CASE( 'K')
        path_int_enc( i) = 47
      CASE( 'L')
        path_int_enc( i) = 48
      CASE( 'M')
        path_int_enc( i) = 49
      CASE( 'N')
        path_int_enc( i) = 50
      CASE( 'O')
        path_int_enc( i) = 51
      CASE( 'P')
        path_int_enc( i) = 52
      CASE( 'Q')
        path_int_enc( i) = 53
      CASE( 'R')
        path_int_enc( i) = 54
      CASE( 'S')
        path_int_enc( i) = 55
      CASE( 'T')
        path_int_enc( i) = 56
      CASE( 'U')
        path_int_enc( i) = 57
      CASE( 'V')
        path_int_enc( i) = 58
      CASE( 'W')
        path_int_enc( i) = 59
      CASE( 'X')
        path_int_enc( i) = 60
      CASE( 'Y')
        path_int_enc( i) = 61
      CASE( 'Z')
        path_int_enc( i) = 62
      CASE( '_')
        path_int_enc( i) = 63
      CASE( '/')
        path_int_enc( i) = 64
      CASE( '(')
        path_int_enc( i) = 65
      CASE( ')')
        path_int_enc( i) = 66
      CASE DEFAULT
        CALL crash('unknown character in routine_path "' // TRIM( subroutine_path) // '"!')
      END SELECT

    END DO

  END SUBROUTINE encode_subroutine_path_as_integer
  SUBROUTINE write_to_SELEN_output_file( SELEN, time)
    ! Write the current model state to the existing output file

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_to_SELEN_output_file'
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN
    REAL(dp),                       INTENT(IN)    :: time

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the file for writing
    CALL open_netcdf_file( SELEN%output%filename, SELEN%output%ncid)

    ! Time
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_time, time, start = (/ SELEN%output%ti /)))

    ! Model data
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_Hi,             SELEN%Hi_glob,     start = (/ 1, SELEN%output%ti/) ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_Hi_rel,         SELEN%Hi_rel_glob, start = (/ 1, SELEN%output%ti/) ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_U,              SELEN%U_glob,      start = (/ 1, SELEN%output%ti/) ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_N,              SELEN%N_glob,      start = (/ 1, SELEN%output%ti/) ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_ocean_function, SELEN%of_glob,     start = (/ 1, SELEN%output%ti/) ))

    ! Close the file
    CALL close_netcdf_file(SELEN%output%ncid)

    ! Increase time frame counter
    SELEN%output%ti = SELEN%output%ti + 1

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_SELEN_output_file

  SUBROUTINE create_resource_tracking_file( netcdf)
    ! Create the resource tracking output file

    USE configuration_module, ONLY: resource_tracker

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_netcdf_resource_tracker), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'create_resource_tracking_file'
    LOGICAL                                           :: file_exists
    INTEGER                                           :: t,nl
    INTEGER                                           :: i,n
    CHARACTER(LEN=256)                                :: var_name, long_name

    IF (.NOT. par%master) RETURN

    ! Set time frame index to 1
    netcdf%ti = 1

    ! Create a new file if none exists and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    netcdf%filename = TRIM(C%output_dir) // '/resource_tracking.nc'
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( netcdf%filename) // '" already exists!')
    END IF

    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
    CALL handle_error( nf90_create( netcdf%filename, IOR( nf90_clobber, nf90_share), netcdf%ncid))

    ! Define dimensions:
    CALL create_dim( netcdf%ncid, netcdf%name_dim_time       , nf90_unlimited, netcdf%id_dim_time       )
    CALL create_dim( netcdf%ncid, netcdf%name_dim_name_length, 1024          , netcdf%id_dim_name_length)

    ! Placeholders for the dimension ID's, for shorter code
    t  = netcdf%id_dim_time
    nl = netcdf%id_dim_name_length

    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.

    ! Dimension variables: time
    CALL create_double_var( netcdf%ncid, netcdf%name_var_time, [t], netcdf%id_var_time, long_name='Time', units='years'   )

    ! Actual variables
    ! ================

    n = SIZE( resource_tracker)

    ALLOCATE( netcdf%id_vars(      n))
    ALLOCATE( netcdf%id_var_names( n))

    DO i = 1, n

      ! Subroutine name
      ! ===============

      ! Generate variable name (name_00001, name_00002, etc.)
      var_name(  1:256) = ' '
      long_name( 1:256) = ' '
      IF     (i < 10) THEN
        WRITE( var_name ,'(A,I1)') 'name_0000', i
      ELSEIF (i < 100) THEN
        WRITE( var_name,'(A,I2)') 'name_000', i
      ELSEIF (i < 1000) THEN
        WRITE( var_name,'(A,I3)') 'name_00', i
      ELSEIF (i < 10000) THEN
        WRITE( var_name,'(A,I4)') 'name_0', i
      ELSEIF (i < 100000) THEN
        WRITE( var_name,'(A,I5)') 'name_', i
      END IF

      WRITE( long_name,'(A,I1)') 'Full name of subroutine #', i

      ! Create the variable in the NetCDF file
      CALL create_int_var( netcdf%ncid, var_name, [nl], netcdf%id_var_names( i),  long_name = long_name)

      ! Computation time
      ! ================

      ! Generate variable name (tcomp_00001, tcomp_00002, etc.)
      var_name(  1:256) = ' '
      long_name( 1:256) = ' '
      IF     (i < 10) THEN
        WRITE( var_name ,'(A,I1)') 'tcomp_0000', i
      ELSEIF (i < 100) THEN
        WRITE( var_name,'(A,I2)') 'tcomp_000', i
      ELSEIF (i < 1000) THEN
        WRITE( var_name,'(A,I3)') 'tcomp_00', i
      ELSEIF (i < 10000) THEN
        WRITE( var_name,'(A,I4)') 'tcomp_0', i
      ELSEIF (i < 100000) THEN
        WRITE( var_name,'(A,I5)') 'tcomp_', i
      END IF

      WRITE( long_name,'(A,I1)') 'Computation time for subroutine #', i

      ! Create the variable in the NetCDF file
      CALL create_double_var( netcdf%ncid, var_name, [t], netcdf%id_vars( i),  long_name = long_name, units = 's', missing_value = 0._dp)

    END DO

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( netcdf%ncid))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( netcdf%ncid))

    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)

  END SUBROUTINE create_resource_tracking_file
  SUBROUTINE create_SELEN_output_file( SELEN)
    ! Create a new NetCDF output file for SELEN (on the irregular global mesh)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_SELEN_output_file'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, three, time

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set time frame index to 1
    SELEN%output%ti = 1

    ! Set output filename
    SELEN%output%filename = TRIM(C%output_dir) // 'SELEN_output.nc'

    ! Create a new restart file if none exists and, to prevent loss of data,
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(SELEN%output%filename))
    IF(file_exists) THEN
      CALL crash('file "' // TRIM( SELEN%output%filename) // '" already exists!')
    END IF

    ! Create netCDF file
    CALL handle_error(nf90_create(SELEN%output%filename,IOR(nf90_clobber,nf90_share),SELEN%output%ncid))

    ! Mesh data
    ! =========

    ! Define dimensions
    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_vi,           SELEN%mesh%nV,           SELEN%output%id_dim_vi          ) ! Vertex indices
    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_ti,           SELEN%mesh%nTri,         SELEN%output%id_dim_ti          ) ! Triangle indices
    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_ci,           SELEN%mesh%nC_mem,       SELEN%output%id_dim_ci          ) ! Connection indices
    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_three,        3,                       SELEN%output%id_dim_three       ) ! 3 (each vertex has three coordinates, each triangle has three vertices)

    ! Placeholders for the dimension ID's, for shorter code
    vi        = SELEN%output%id_dim_vi
    ti        = SELEN%output%id_dim_ti
    ci        = SELEN%output%id_dim_ci
    three     = SELEN%output%id_dim_three

    ! Define variables
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_V,                [vi,  three], SELEN%output%id_var_V,                long_name='Vertex coordinates', units='m')
    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_Tri,              [ti,  three], SELEN%output%id_var_Tri,              long_name='Vertex indices')
    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_nC,               [vi        ], SELEN%output%id_var_nC,               long_name='Number of connected vertices')
    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_C,                [vi,  ci   ], SELEN%output%id_var_C,                long_name='Indices of connected vertices')
    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_niTri,            [vi        ], SELEN%output%id_var_niTri,            long_name='Number of inverse triangles')
    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_iTri,             [vi,  ci   ], SELEN%output%id_var_iTri,             long_name='Indices of inverse triangles')

    ! Model output
    ! ============

    ! Define dimensions
    CALL create_dim( SELEN%output%ncid, SELEN%output%name_dim_time,  nf90_unlimited, SELEN%output%id_dim_time ) ! Time frames

    ! Placeholders for the dimension ID's, for shorter code
    time  = SELEN%output%id_dim_time

    ! Define dimension variables
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_time,  [time  ], SELEN%output%id_var_time,  long_name='Time', units='years')

    ! Define model data variables
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_lat,              [vi             ], SELEN%output%id_var_lat,              long_name='Latitude', units='degrees north')
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_lon,              [vi             ], SELEN%output%id_var_lon,              long_name='Longtitude', units='degrees east')
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_Hi,               [vi,        time], SELEN%output%id_var_Hi,               long_name='Surface load', units='mie')
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_Hi_rel,           [vi,        time], SELEN%output%id_var_Hi_rel,           long_name='Relative surface load', units='mie')
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_U,                [vi,        time], SELEN%output%id_var_U,                long_name='Land surface change', units='m')
    CALL create_double_var( SELEN%output%ncid, SELEN%output%name_var_N,                [vi,        time], SELEN%output%id_var_N,                long_name='Sea surface change', units='m')
    CALL create_int_var(    SELEN%output%ncid, SELEN%output%name_var_ocean_function,   [vi,        time], SELEN%output%id_var_ocean_function,   long_name='Ocean function (1 = ocean)')

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( SELEN%output%ncid))

    ! Write mesh data
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_V,               SELEN%mesh%V             ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_Tri,             SELEN%mesh%Tri           ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_nC,              SELEN%mesh%nC            ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_C,               SELEN%mesh%C             ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_niTri,           SELEN%mesh%niTri         ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_iTri,            SELEN%mesh%iTri          ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_lat,             SELEN%mesh%lat           ))
    CALL handle_error( nf90_put_var( SELEN%output%ncid, SELEN%output%id_var_lon,             SELEN%mesh%lon           ))

    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( SELEN%output%ncid))

    ! Close the file
    CALL close_netcdf_file(SELEN%output%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_SELEN_output_file

!  SUBROUTINE inquire_restart_file_isotopes(    restart)
!    ! Check if the right dimensions and variables are present in the file.
!
!    IMPLICIT NONE
!
!    ! Input variables:
!    TYPE(type_restart_data), INTENT(INOUT)        :: restart
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_restart_file_isotopes'
!    INTEGER                                       :: x, y, t
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    ! Open the netcdf file
!    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)
!
!    ! Inquire dimensions id's. Check that all required dimensions exist, return their lengths.
!    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_x,     restart%nx, restart%netcdf%id_dim_x    )
!    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_y,     restart%ny, restart%netcdf%id_dim_y    )
!    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_time,  restart%nt, restart%netcdf%id_dim_time )
!
!    ! Abbreviations for shorter code
!    x = restart%netcdf%id_dim_x
!    y = restart%netcdf%id_dim_y
!    t = restart%netcdf%id_dim_time
!
!    ! Inquire variable ID's; make sure that each variable has the correct dimensions.
!
!    ! Dimensions
!    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_x,                (/ x             /), restart%netcdf%id_var_x   )
!    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_y,                (/    y          /), restart%netcdf%id_var_y   )
!    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_time,             (/             t /), restart%netcdf%id_var_time)
!
!    ! Data
!    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_IsoIce,           (/ x, y,       t /), restart%netcdf%id_var_IsoIce )
!
!    ! Close the netcdf file
!    CALL close_netcdf_file( restart%netcdf%ncid)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE inquire_restart_file_isotopes
!
!  SUBROUTINE read_restart_file_SMB(            restart, time_to_restart_from)
!    ! Read the restart netcdf file
!
!    IMPLICIT NONE
!
!    ! Input variables:
!    TYPE(type_restart_data),        INTENT(INOUT) :: restart
!    REAL(dp),                       INTENT(IN)    :: time_to_restart_from
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_restart_file_SMB'
!    INTEGER                                       :: ti, ti_min
!    REAL(dp)                                      :: dt, dt_min
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    ! Open the netcdf file
!    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)
!
!    ! Read x,y
!    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_x, restart%x, start=(/1/) ))
!    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_y, restart%y, start=(/1/) ))
!
!    ! Read time, determine which time frame to read
!    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_time, restart%time, start=(/1/) ))
!
!    IF (time_to_restart_from < MINVAL(restart%time) .OR. time_to_restart_from > MAXVAL(restart%time)) THEN
!      CALL crash('time_to_restart_from = {dp_01} outside range of restart file!', dp_01 = time_to_restart_from)
!    END IF
!
!    ti_min = 0
!    dt_min = 1E8_dp
!    DO ti = 1, restart%nt
!      dt = ABS(restart%time( ti) - time_to_restart_from)
!      IF (dt < dt_min) THEN
!        ti_min = ti
!        dt_min = dt
!      END IF
!    END DO
!    ti = ti_min
!
!    IF (dt_min > 0._dp) THEN
!      CALL warning('no exact match for time_to_restart_from = {dp_01} yr in restart file! Reading closest match = {dp_02} yr instead.', &
!        dp_01 = time_to_restart_from, dp_02 = restart%time( ti))
!    END IF
!    IF (time_to_restart_from /= C%start_time_of_run) THEN
!      CALL warning('starting run at t = {dp_01} yr with restart data at t = {dp_02} yr!', &
!        dp_01 = C%start_time_of_run, dp_02 = time_to_restart_from)
!    END IF
!
!    ! Read the data
!    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_FirnDepth,        restart%FirnDepth,        start = (/ 1, 1, 1, ti /), count = (/ restart%nx, restart%ny, 12,         1 /) ))
!    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_MeltPreviousYear, restart%MeltPreviousYear, start = (/ 1, 1,    ti /), count = (/ restart%nx, restart%ny,             1 /) ))
!
!    ! Close the netcdf file
!    CALL close_netcdf_file( restart%netcdf%ncid)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE read_restart_file_SMB
!  SUBROUTINE read_restart_file_isotopes(       restart, time_to_restart_from)
!    ! Read the restart netcdf file
!
!    IMPLICIT NONE
!
!    ! Input variables:
!    TYPE(type_restart_data),        INTENT(INOUT) :: restart
!    REAL(dp),                       INTENT(IN)    :: time_to_restart_from
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_restart_file_isotopes'
!    INTEGER                                       :: ti, ti_min
!    REAL(dp)                                      :: dt, dt_min
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    ! Open the netcdf file
!    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)
!
!    ! Read x,y
!    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_x, restart%x, start=(/1/) ))
!    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_y, restart%y, start=(/1/) ))
!
!    ! Read time, determine which time frame to read
!    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_time, restart%time, start=(/1/) ))
!
!    IF (time_to_restart_from < MINVAL(restart%time) .OR. time_to_restart_from > MAXVAL(restart%time)) THEN
!      CALL crash('time_to_restart_from = {dp_01} outside range of restart file!', dp_01 = time_to_restart_from)
!    END IF
!
!    ti_min = 0
!    dt_min = 1E8_dp
!    DO ti = 1, restart%nt
!      dt = ABS(restart%time( ti) - time_to_restart_from)
!      IF (dt < dt_min) THEN
!        ti_min = ti
!        dt_min = dt
!      END IF
!    END DO
!    ti = ti_min
!
!    IF (dt_min > 0._dp) THEN
!      CALL warning('no exact match for time_to_restart_from = {dp_01} yr in restart file! Reading closest match = {dp_02} yr instead.', &
!        dp_01 = time_to_restart_from, dp_02 = restart%time( ti))
!    END IF
!    IF (time_to_restart_from /= C%start_time_of_run) THEN
!      CALL warning('starting run at t = {dp_01} yr with restart data at t = {dp_02} yr!', &
!        dp_01 = C%start_time_of_run, dp_02 = time_to_restart_from)
!    END IF
!
!    ! Read the data
!    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_IsoIce,           restart%IsoIce,           start = (/ 1, 1,    ti /), count = (/ restart%nx, restart%ny,             1 /) ))
!
!    ! Close the netcdf file
!    CALL close_netcdf_file( restart%netcdf%ncid)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE read_restart_file_isotopes

!  SUBROUTINE read_inverse_routine_history_dT_glob(         forcing, filename)
!    ! Read the inverse routine history from the specified NetCDF file
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
!    CHARACTER(LEN=256),             INTENT(IN)    :: filename
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_inverse_routine_history_dT_glob'
!    CHARACTER(LEN=256) :: dummy_char
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    dummy_char = forcing%netcdf_ins%filename
!    dummy_char = filename
!
!    CALL crash('need to fix the inverse routine stuff to cope with restarting!')
!
!    ! Local variables:
!    INTEGER                                       :: ncid, id_dim_nH, id_dim_time, nt, id_var_H, id_var_time
!    INTEGER                                       :: int_dummy
!    INTEGER                                       :: ti, ti_min
!    REAL(dp)                                      :: dt, dt_min
!    REAL(dp), DIMENSION(:), ALLOCATABLE           :: time
!
!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    WRITE(0,*) ' Reading inverse routine dT_glob_history from file "', TRIM(filename), '"...'
!
!    ! Check if this file contains the required history
!
!    ! Open the netcdf file
!    CALL open_netcdf_file( filename, ncid)
!
!    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
!    CALL inquire_dim( ncid, 'ndT_glob_history', int_dummy, id_dim_nH)
!    IF (int_dummy /= forcing%ndT_glob_history) THEN
!      WRITE(0,*) '  ERROR - ndT_glob_history in file "', filename, '" doesnt match forcing configuration!'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!    CALL inquire_dim( ncid, 'time', nt, id_dim_time)
!
!    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
!    CALL inquire_double_var( ncid, 'dT_glob_history',  (/ id_dim_nH, id_dim_time/), id_var_H)
!    CALL inquire_double_var( ncid, 'time',             (/            id_dim_time/), id_var_time)
!
!    ! Read time, determine which time frame to read
!    ALLOCATE( time( nt))
!
!    CALL handle_error(nf90_get_var( ncid, id_var_time, time))
!
!    IF (C%time_to_restart_from < MINVAL(time) .OR. C%time_to_restart_from > MAXVAL(time)) THEN
!      WRITE(0,*) '  ERROR - time_to_restart_from ', C%time_to_restart_from, ' outside range of restart file!'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!
!    ti_min = 0
!    dt_min = 1E8_dp
!    DO ti = 1, nt
!      dt = ABS(time( ti) - C%time_to_restart_from)
!      IF (dt < dt_min) THEN
!        ti_min = ti
!        dt_min = dt
!      END IF
!    END DO
!    ti = ti_min
!
!    IF (dt_min > 0._dp) THEN
!      WRITE(0,*) '  WARNING - no exact match for time_to_restart_from ', C%time_to_restart_from, ' in restart file! Reading closest match ', time( ti), ' instead.'
!    END IF
!    IF (C%time_to_restart_from /= C%start_time_of_run) THEN
!      WRITE(0,*) '  WARNING - starting run at t = ', C%start_time_of_run, ' with restart data at t = ', C%time_to_restart_from
!    END IF
!
!    DEALLOCATE( time)
!
!    ! Read inverse routine history
!    CALL handle_error(nf90_get_var( ncid, id_var_H, forcing%dT_glob_history, start=(/1,ti/) ))
!
!    ! Close the netcdf file
!    CALL close_netcdf_file( ncid)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE read_inverse_routine_history_dT_glob
!  SUBROUTINE read_inverse_routine_history_dT_glob_inverse( forcing, filename)
!    ! Read the inverse routine history from the specified NetCDF file
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
!    CHARACTER(LEN=256),             INTENT(IN)    :: filename
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_inverse_routine_history_dT_glob_inverse'
!    CHARACTER(LEN=256) :: dummy_char
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    dummy_char = forcing%netcdf_ins%filename
!    dummy_char = filename
!
!    CALL crash('need to fix the inverse routine stuff to cope with restarting!')

!    ! Local variables:
!    INTEGER                                       :: ncid, id_dim_nH, id_dim_time, nt, id_var_H, id_var_time
!    INTEGER                                       :: int_dummy
!    INTEGER                                       :: ti, ti_min
!    REAL(dp)                                      :: dt, dt_min
!    REAL(dp), DIMENSION(:), ALLOCATABLE           :: time
!
!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    WRITE(0,*) ' Reading inverse routine dT_glob_inverse_history from file "', TRIM(filename), '"...'
!
!    ! Check if this file contains the required history
!
!    ! Open the netcdf file
!    CALL open_netcdf_file( filename, ncid)
!
!    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
!    CALL inquire_dim( ncid, 'ndT_glob_inverse_history', int_dummy, id_dim_nH)
!    IF (int_dummy /= forcing%ndT_glob_inverse_history) THEN
!      WRITE(0,*) '  ERROR - ndT_glob_inverse_history in file "', filename, '" doesnt match forcing configuration!'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!    CALL inquire_dim( ncid, 'time', nt, id_dim_time)
!
!    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
!    CALL inquire_double_var( ncid, 'dT_glob_inverse_history',  (/ id_dim_nH, id_dim_time/), id_var_H)
!    CALL inquire_double_var( ncid, 'time',                     (/            id_dim_time/), id_var_time)
!
!    ! Read time, determine which time frame to read
!    ALLOCATE( time( nt))
!
!    CALL handle_error(nf90_get_var( ncid, id_var_time, time))
!
!    IF (C%time_to_restart_from < MINVAL(time) .OR. C%time_to_restart_from > MAXVAL(time)) THEN
!      WRITE(0,*) '  ERROR - time_to_restart_from ', C%time_to_restart_from, ' outside range of restart file!'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!
!    ti_min = 0
!    dt_min = 1E8_dp
!    DO ti = 1, nt
!      dt = ABS(time( ti) - C%time_to_restart_from)
!      IF (dt < dt_min) THEN
!        ti_min = ti
!        dt_min = dt
!      END IF
!    END DO
!    ti = ti_min
!
!    IF (dt_min > 0._dp) THEN
!      WRITE(0,*) '  WARNING - no exact match for time_to_restart_from ', C%time_to_restart_from, ' in restart file! Reading closest match ', time( ti), ' instead.'
!    END IF
!    IF (C%time_to_restart_from /= C%start_time_of_run) THEN
!      WRITE(0,*) '  WARNING - starting run at t = ', C%start_time_of_run, ' with restart data at t = ', C%time_to_restart_from
!    END IF
!
!    DEALLOCATE( time)
!
!    ! Read inverse routine history
!    CALL handle_error(nf90_get_var( ncid, id_var_H, forcing%dT_glob_inverse_history, start=(/1,ti/) ))
!
!    ! Close the netcdf file
!    CALL close_netcdf_file( ncid)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE read_inverse_routine_history_dT_glob_inverse
!  SUBROUTINE read_inverse_routine_history_CO2_inverse(     forcing, filename)
!    ! Read the inverse routine history from the specified NetCDF file
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
!    CHARACTER(LEN=256),             INTENT(IN)    :: filename
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_inverse_routine_history_CO2_inverse'
!    CHARACTER(LEN=256) :: dummy_char
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    dummy_char = forcing%netcdf_ins%filename
!    dummy_char = filename
!
!    CALL crash('need to fix the inverse routine stuff to cope with restarting!')
!
!    ! Local variables:
!    INTEGER                                       :: ncid, id_dim_nH, id_dim_time, nt, id_var_H, id_var_time
!    INTEGER                                       :: int_dummy
!    INTEGER                                       :: ti, ti_min
!    REAL(dp)                                      :: dt, dt_min
!    REAL(dp), DIMENSION(:), ALLOCATABLE           :: time
!
!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    WRITE(0,*) ' Reading inverse routine CO2_inverse_history from file "', TRIM(filename), '"...'
!
!    ! Check if this file contains the required history
!
!    ! Open the netcdf file
!    CALL open_netcdf_file( filename, ncid)
!
!    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
!    CALL inquire_dim( ncid, 'nCO2_inverse_history', int_dummy, id_dim_nH)
!    IF (int_dummy /= forcing%nCO2_inverse_history) THEN
!      WRITE(0,*) '  ERROR - nCO2_inverse_history in file "', filename, '" doesnt match forcing configuration!'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!    CALL inquire_dim( ncid, 'time', nt, id_dim_time)
!
!    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
!    CALL inquire_double_var( ncid, 'CO2_inverse_history',  (/ id_dim_nH, id_dim_time/), id_var_H)
!    CALL inquire_double_var( ncid, 'time',                 (/            id_dim_time/), id_var_time)
!
!    ! Read time, determine which time frame to read
!    ALLOCATE( time( nt))
!
!    CALL handle_error(nf90_get_var( ncid, id_var_time, time))
!
!    IF (C%time_to_restart_from < MINVAL(time) .OR. C%time_to_restart_from > MAXVAL(time)) THEN
!      WRITE(0,*) '  ERROR - time_to_restart_from ', C%time_to_restart_from, ' outside range of restart file!'
!      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!    END IF
!
!    ti_min = 0
!    dt_min = 1E8_dp
!    DO ti = 1, nt
!      dt = ABS(time( ti) - C%time_to_restart_from)
!      IF (dt < dt_min) THEN
!        ti_min = ti
!        dt_min = dt
!      END IF
!    END DO
!    ti = ti_min
!
!    IF (dt_min > 0._dp) THEN
!      WRITE(0,*) '  WARNING - no exact match for time_to_restart_from ', C%time_to_restart_from, ' in restart file! Reading closest match ', time( ti), ' instead.'
!    END IF
!    IF (C%time_to_restart_from /= C%start_time_of_run) THEN
!      WRITE(0,*) '  WARNING - starting run at t = ', C%start_time_of_run, ' with restart data at t = ', C%time_to_restart_from
!    END IF
!
!    DEALLOCATE( time)
!
!    ! Read inverse routine history
!    CALL handle_error(nf90_get_var( ncid, id_var_H, forcing%CO2_inverse_history, start=(/1,ti/) ))
!
!    ! Close the netcdf file
!    CALL close_netcdf_file( ncid)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE read_inverse_routine_history_CO2_inverse

! Read all kinds of input files
! =============================

  SUBROUTINE setup_grid_from_file( filename, grid)
    ! Set up a regional x/y-grid from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                 INTENT(IN)    :: filename
    TYPE(type_grid),                    INTENT(INOUT) :: grid

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'setup_grid_from_file'
    INTEGER                                           :: ncid, id_dim_x, id_dim_y, id_var_x, id_var_y

!    ! Add routine to path
!    CALL init_routine( routine_name)

    ! Allocate memory for nx/ny
    CALL allocate_shared_int_0D( grid%nx, grid%wnx)
    CALL allocate_shared_int_0D( grid%ny, grid%wny)

    ! Inquire grid size from file
    IF (par%master) THEN

      CALL open_netcdf_file( filename, ncid)

      CALL inquire_dim( ncid, 'x', grid%nx, id_dim_x)
      CALL inquire_dim( ncid, 'y', grid%ny, id_dim_y)

      CALL close_netcdf_file( ncid)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Allocate memory for x/y
    CALL allocate_shared_dp_1D( grid%nx, grid%x, grid%wx)
    CALL allocate_shared_dp_1D( grid%ny, grid%y, grid%wy)
    CALL allocate_shared_dp_0D( grid%dx, grid%wdx)

    ! Read grid from file
    IF (par%master) THEN

      CALL open_netcdf_file( filename, ncid)

      CALL inquire_single_or_double_var( ncid, 'x', [id_dim_x], id_var_x)
      CALL inquire_single_or_double_var( ncid, 'y', [id_dim_y], id_var_y)

      ! Read x,y variables
      CALL handle_error( nf90_get_var( ncid, id_var_x, grid%x ))
      CALL handle_error( nf90_get_var( ncid, id_var_y, grid%y ))

      grid%dx = grid%x( 2) - grid%x( 1)

      CALL close_netcdf_file( ncid)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Assign range to each processor
    CALL partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)

!    ! Finalise routine path
!    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_grid_from_file
  SUBROUTINE setup_lonlat_grid_from_file( filename, grid)
    ! Set up a global lon/lat-grid from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                 INTENT(IN)    :: filename
    TYPE(type_grid_lonlat),             INTENT(INOUT) :: grid

    ! Local variables:
    INTEGER                                           :: ncid, id_dim_lon, id_dim_lat, id_var_lon, id_var_lat

    ! Allocate memory for nlon/nlat
    CALL allocate_shared_int_0D( grid%nlon, grid%wnlon)
    CALL allocate_shared_int_0D( grid%nlat, grid%wnlat)

    ! Inquire grid size from file
    IF (par%master) THEN

      CALL open_netcdf_file( filename, ncid)

      ! Check different possible dimension names
      CALL inquire_dim( ncid, 'lon', grid%nlon, id_dim_lon)
      CALL inquire_dim( ncid, 'lat', grid%nlat, id_dim_lat)

      CALL close_netcdf_file( ncid)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Allocate memory for lon/lat
    CALL allocate_shared_dp_1D( grid%nlon, grid%lon, grid%wlon)
    CALL allocate_shared_dp_1D( grid%nlat, grid%lat, grid%wlat)

    ! Read grid from file
    IF (par%master) THEN

      CALL open_netcdf_file( filename, ncid)

      CALL inquire_single_or_double_var( ncid, 'lon', [id_dim_lon], id_var_lon)
      CALL inquire_single_or_double_var( ncid, 'lat', [id_dim_lat], id_var_lat)

      ! Read lon,lat variables
      CALL handle_error( nf90_get_var( ncid, id_var_lon, grid%lon))
      CALL handle_error( nf90_get_var( ncid, id_var_lat, grid%lat))

      CALL close_netcdf_file( ncid)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Assign range to each processor
    CALL partition_list( grid%nlon, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%nlat, par%i, par%n, grid%j1, grid%j2)

  END SUBROUTINE setup_lonlat_grid_from_file


!  ! Direct regional SMB forcing
!  SUBROUTINE inquire_direct_regional_SMB_forcing_file( clim)
!
!    IMPLICIT NONE
!
!    ! Output variable
!    TYPE(type_direct_SMB_forcing_regional), INTENT(INOUT) :: clim
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_direct_regional_SMB_forcing_file'
!    INTEGER                                :: x,y,t
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    ! Safety
!    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
!      CALL crash('should only be called when choice_SMB_model = "direct_regional"!')
!    END IF
!
!    ! Open the netcdf file
!    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
!
!    ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
!    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_x,     clim%nx_raw, clim%netcdf%id_dim_x    )
!    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_y,     clim%ny_raw, clim%netcdf%id_dim_y    )
!    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_time,  clim%nyears, clim%netcdf%id_dim_time )
!
!    ! Abbreviate dimension ID's for more readable code
!    x = clim%netcdf%id_dim_x
!    y = clim%netcdf%id_dim_y
!    t = clim%netcdf%id_dim_time
!
!    ! Inquire variable id's. Make sure that each variable has the correct dimensions.
!    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_x,        (/ x       /), clim%netcdf%id_var_x       )
!    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_y,        (/    y    /), clim%netcdf%id_var_y       )
!    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_time,     (/       t /), clim%netcdf%id_var_time    )
!
!    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m_year, (/ x, y, t /), clim%netcdf%id_var_T2m_year)
!    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_SMB_year, (/ x, y, t /), clim%netcdf%id_var_SMB_year)
!
!    ! Close the netcdf file
!    CALL close_netcdf_file( clim%netcdf%ncid)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE inquire_direct_regional_SMB_forcing_file
!  SUBROUTINE read_direct_regional_SMB_file_timeframes( clim, ti0, ti1)
!
!    IMPLICIT NONE
!
!    ! In/output variables:
!    TYPE(type_direct_SMB_forcing_regional), INTENT(INOUT) :: clim
!    INTEGER,                        INTENT(IN)    :: ti0, ti1
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_regional_SMB_file_timeframes'
!    INTEGER                                       :: i,j
!    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE     :: T2m_temp0, T2m_temp1, SMB_temp0, SMB_temp1
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    ! Safety
!    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
!      CALL crash('should only be called when choice_SMB_model = "direct_regional"!')
!    END IF
!
!    ! Temporary memory to store the data read from the netCDF file
!    ALLOCATE( T2m_temp0( clim%nx_raw, clim%ny_raw, 1))
!    ALLOCATE( T2m_temp1( clim%nx_raw, clim%ny_raw, 1))
!    ALLOCATE( SMB_temp0( clim%nx_raw, clim%ny_raw, 1))
!    ALLOCATE( SMB_temp1( clim%nx_raw, clim%ny_raw, 1))
!
!    ! Open netcdf file
!    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
!
!    ! Read the data
!    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m_year, T2m_temp0, start = (/ 1, 1, ti0 /), count = (/ clim%nx_raw, clim%ny_raw, 1 /), stride = (/ 1, 1, 1 /) ))
!    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m_year, T2m_temp1, start = (/ 1, 1, ti1 /), count = (/ clim%nx_raw, clim%nx_raw, 1 /), stride = (/ 1, 1, 1 /) ))
!    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_SMB_year, SMB_temp0, start = (/ 1, 1, ti0 /), count = (/ clim%nx_raw, clim%ny_raw, 1 /), stride = (/ 1, 1, 1 /) ))
!    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_SMB_year, SMB_temp1, start = (/ 1, 1, ti1 /), count = (/ clim%nx_raw, clim%nx_raw, 1 /), stride = (/ 1, 1, 1 /) ))
!
!     ! Close netcdf file
!    CALL close_netcdf_file( clim%netcdf%ncid)
!
!    ! Store the data in the shared memory structure
!    DO i = 1, clim%nx_raw
!    DO j = 1, clim%ny_raw
!      clim%T2m_year0_raw( j,i) = T2m_temp0( i,j,1)
!      clim%T2m_year1_raw( j,i) = T2m_temp1( i,j,1)
!      clim%SMB_year0_raw( j,i) = SMB_temp0( i,j,1)
!      clim%SMB_year1_raw( j,i) = SMB_temp1( i,j,1)
!    END DO
!    END DO
!
!    ! Clean up after yourself
!    DEALLOCATE( T2m_temp0)
!    DEALLOCATE( T2m_temp1)
!    DEALLOCATE( SMB_temp0)
!    DEALLOCATE( SMB_temp1)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE read_direct_regional_SMB_file_timeframes
!  SUBROUTINE read_direct_regional_SMB_file_time_xy( clim)
!
!    IMPLICIT NONE
!
!    ! Output variable
!    TYPE(type_direct_SMB_forcing_regional), INTENT(INOUT) :: clim
!
!    ! Local variables:
!    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_regional_SMB_file_time_xy'
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    IF (.NOT. par%master) THEN
!      CALL finalise_routine( routine_name)
!      RETURN
!    END IF
!
!    ! Safety
!    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
!      CALL crash('should only be called when choice_SMB_model = "direct_regional"!')
!    END IF
!
!    ! Open the netcdf file
!    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
!
!    ! Read the data
!    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_time, clim%time,  start = (/ 1 /) ))
!    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_x,    clim%x_raw, start = (/ 1 /) ))
!    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_y,    clim%y_raw, start = (/ 1 /) ))
!
!    ! Close the netcdf file
!    CALL close_netcdf_file( clim%netcdf%ncid)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
!
!  END SUBROUTINE read_direct_regional_SMB_file_time_xy

  ! Global topography for SELEN
  SUBROUTINE inquire_SELEN_global_topo_file( SELEN)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_SELEN_global_topo_file'
    INTEGER                                       :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file( SELEN%netcdf_topo%filename, SELEN%netcdf_topo%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_dim_vi,    SELEN%mesh%nV,     SELEN%netcdf_topo%id_dim_vi)
    CALL inquire_dim( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_dim_ti,    SELEN%mesh%nTri,   SELEN%netcdf_topo%id_dim_ti)
    CALL inquire_dim( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_dim_ci,    SELEN%mesh%nC_mem, SELEN%netcdf_topo%id_dim_ci)
    CALL inquire_dim( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_dim_three, int_dummy,         SELEN%netcdf_topo%id_dim_three)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_V,     (/ SELEN%netcdf_topo%id_dim_vi, SELEN%netcdf_topo%id_dim_three /),  SELEN%netcdf_topo%id_var_V       )
    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_Tri,   (/ SELEN%netcdf_topo%id_dim_ti, SELEN%netcdf_topo%id_dim_three /),  SELEN%netcdf_topo%id_var_Tri     )
    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_nC,    (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_nC      )
    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_C,     (/ SELEN%netcdf_topo%id_dim_vi, SELEN%netcdf_topo%id_dim_ci    /),  SELEN%netcdf_topo%id_var_C       )
    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_niTri, (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_niTri   )
    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_iTri,  (/ SELEN%netcdf_topo%id_dim_vi, SELEN%netcdf_topo%id_dim_ci    /),  SELEN%netcdf_topo%id_var_iTri    )
    CALL inquire_double_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_lat,   (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_lat     )
    CALL inquire_double_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_lon,   (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_lon     )
    CALL inquire_double_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_Hb,    (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_Hb      )
    CALL inquire_int_var(    SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%name_var_ianc,  (/ SELEN%netcdf_topo%id_dim_vi                                 /),  SELEN%netcdf_topo%id_var_ianc    )

    ! Close the netcdf file
    CALL close_netcdf_file(SELEN%netcdf_topo%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_SELEN_global_topo_file
  SUBROUTINE read_SELEN_global_topo_file( SELEN)
    ! Read the init netcdf file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_SELEN_global_topo_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Open the netcdf file
    CALL open_netcdf_file(SELEN%netcdf_topo%filename, SELEN%netcdf_topo%ncid)

    ! Read the data
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_V,     SELEN%mesh%V,     start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_Tri,   SELEN%mesh%Tri,   start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_nC,    SELEN%mesh%nC,    start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_C,     SELEN%mesh%C,     start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_niTri, SELEN%mesh%niTri, start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_iTri,  SELEN%mesh%iTri,  start = (/ 1, 1 /) ))

    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_lat,   SELEN%mesh%lat,   start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_lon,   SELEN%mesh%lon,   start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_Hb,    SELEN%topo_ref,   start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( SELEN%netcdf_topo%ncid, SELEN%netcdf_topo%id_var_ianc,  SELEN%mesh%ianc,  start = (/ 1    /) ))

    ! Close the netcdf file
    CALL close_netcdf_file(SELEN%netcdf_topo%ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_SELEN_global_topo_file

  ! Write CSR matrix to file
  SUBROUTINE write_CSR_matrix_to_NetCDF( A_CSR, filename)
    ! Write a CSR matrix to a NetCDF file

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: A_CSR
    CHARACTER(LEN=*),                    INTENT(IN)    :: filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_CSR_matrix_to_NetCDF'
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_m, id_dim_mp1, id_dim_n, id_dim_nnz
    INTEGER                                            :: id_var_ptr, id_var_index, id_var_val

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN

      ! Safety
      INQUIRE(EXIST=file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL crash('file "' // TRIM( filename) // '" already exists!')
      END IF

      WRITE(0,*) '   NOTE: writing CSR matrix to file "', TRIM(filename), '"'

      ! Create netCDF file
      CALL handle_error( nf90_create( TRIM(C%output_dir)//TRIM(filename), IOR(nf90_clobber,nf90_share), ncid))

      ! Define dimensions:
      CALL create_dim( ncid, 'm',      A_CSR%m      , id_dim_m  )
      CALL create_dim( ncid, 'mplus1', A_CSR%m+1    , id_dim_mp1)
      CALL create_dim( ncid, 'n',      A_CSR%n      , id_dim_n  )
      CALL create_dim( ncid, 'nnz',    A_CSR%nnz_max, id_dim_nnz)

      ! Define variables:
      ! The order of the CALL statements for the different variables determines their
      ! order of appearence in the netcdf file.

      CALL create_int_var(    ncid, 'ptr',   [id_dim_mp1], id_var_ptr  , long_name = 'ptr'  )
      CALL create_int_var(    ncid, 'index', [id_dim_nnz], id_var_index, long_name = 'index')
      CALL create_double_var( ncid, 'val',   [id_dim_nnz], id_var_val  , long_name = 'val'  )

      ! Leave definition mode
      CALL handle_error( nf90_enddef( ncid))

      ! Write data
      CALL handle_error( nf90_put_var( ncid, id_var_ptr  , A_CSR%a_ptr  ))
      CALL handle_error( nf90_put_var( ncid, id_var_index, A_CSR%a_index))
      CALL handle_error( nf90_put_var( ncid, id_var_val  , A_CSR%a_val  ))

      ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
      CALL handle_error(nf90_sync( ncid))

      ! Close the file
      CALL close_netcdf_file( ncid)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_CSR_matrix_to_NetCDF

! Some general useful stuff
! =========================

  ! "Flip" data fields from internal [x,y] = [j,i] to output [x,y] = [i,j] and write them to output
  SUBROUTINE write_data_to_file_int_2D( ncid, nx, ny,     var_id, d, start_vec)
    ! Take a 2D data fields, flip it back from [z,y,x] to [x,y,z], and write it to an output NetCDF file

    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid, var_id
    INTEGER,  DIMENSION(:    ), OPTIONAL, INTENT(IN)    :: start_vec
    INTEGER,  DIMENSION(ny,nx), INTENT(IN)    :: d
    INTEGER,                    INTENT(IN)    :: nx, ny

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_data_to_file_int_2D'
    INTEGER                                   :: i,j
    REAL(dp), DIMENSION(nx,ny)                :: d_flip

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = 1, nx
    DO j = 1, ny
      d_flip( i,j) = d( j,i)
    END DO
    END DO

    CALL handle_error( nf90_put_var( ncid, var_id, d_flip, start=start_vec))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_data_to_file_int_2D
  SUBROUTINE write_data_to_file_dp_2D(  ncid, nx, ny,     var_id, d, start_vec)
    ! Take a 2D data fields, flip it back from [z,y,x] to [x,y,z], and write it to an output NetCDF file

    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid, var_id
    INTEGER,  DIMENSION(:    ), OPTIONAL, INTENT(IN)    :: start_vec
    REAL(dp), DIMENSION(ny,nx), INTENT(IN)    :: d
    INTEGER,                    INTENT(IN)    :: nx, ny

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_data_to_file_dp_2D'
    INTEGER                                   :: i,j
    REAL(dp), DIMENSION(nx,ny)                :: d_flip

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = 1, nx
    DO j = 1, ny
      d_flip( i,j) = d( j,i)
    END DO
    END DO

    CALL handle_error( nf90_put_var( ncid, var_id, d_flip, start=start_vec))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_data_to_file_dp_2D
  SUBROUTINE write_data_to_file_dp_3D(  ncid, nx, ny, nz, var_id, d, start_vec)
    ! Take a 2D data fields, flip it back from [z,y,x] to [x,y,z], and write it to an output NetCDF file

    IMPLICIT NONE

    ! Input variables:
    INTEGER,                       INTENT(IN)    :: ncid, var_id
    INTEGER,  DIMENSION(:       ), INTENT(IN)    :: start_vec
    REAL(dp), DIMENSION(nz,ny,nx), INTENT(IN)    :: d
    INTEGER,                       INTENT(IN)    :: nx, ny, nz

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_data_to_file_dp_3D'
    INTEGER                                      :: i,j,k
    REAL(dp), DIMENSION(nx,ny,nz)                :: d_flip

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = 1, nx
    DO j = 1, ny
    DO k = 1, nz
      d_flip( i,j,k) = d( k,j,i)
    END DO
    END DO
    END DO

    CALL handle_error( nf90_put_var( ncid, var_id, d_flip, start=start_vec))

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_data_to_file_dp_3D

  ! Determine if a NetCDF file contains data on a global lon/lat-grid or a regional x/y-grid
  SUBROUTINE determine_file_grid_type( filename, file_grid_type)
    ! Determine if a NetCDF file contains data on a global lon/lat-grid or a regional x/y-grid

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                 INTENT(IN)    :: filename
    CHARACTER(LEN=256),                 INTENT(OUT)   :: file_grid_type

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'determine_file_grid_type'
    INTEGER                                           :: ncid
    LOGICAL                                           :: hasx, hasy, haslon, haslat

    ! Add routine to path
    CALL init_routine( routine_name)

    file_grid_type = ''

    IF (par%master) THEN

      ! Open the NetCDF file
      CALL open_netcdf_file( filename, ncid)

      ! See which dimensions the file has
      CALL inquire_if_dim_exists_nocrash( ncid, 'x'  , hasx  )
      CALL inquire_if_dim_exists_nocrash( ncid, 'x'  , hasy  )
      CALL inquire_if_dim_exists_nocrash( ncid, 'lon', haslon)
      CALL inquire_if_dim_exists_nocrash( ncid, 'lat', haslat)

      IF (hasx .AND. hasy) THEN
        file_grid_type = 'x/y'
      ELSEIF (haslon .AND. haslat) THEN
        file_grid_type = 'lon/lat'
      ELSE
        CALL crash('couldnt find x/y or lon/lat dimensions in file "' // TRIM( filename) // '"!')
      END IF

      ! Close the NetCDF file
      CALL close_netcdf_file( ncid)

    END IF ! IF (par%master) THEN
    CALL MPI_BCAST( file_grid_type, 256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_file_grid_type

  ! Basic NetCDF wrapper functions
  SUBROUTINE open_netcdf_file( filename, ncid)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: filename
    INTEGER,          INTENT(OUT) :: ncid

    ! Open netCDF file:
    CALL handle_error(nf90_open(filename, IOR(nf90_write,nf90_share), ncid))

  END SUBROUTINE open_netcdf_file
  SUBROUTINE close_netcdf_file( ncid)
    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: ncid

    ! Close netCDF file:
    CALL handle_error(nf90_close(ncid))

  END SUBROUTINE close_netcdf_file
  SUBROUTINE create_dim( ncid, dim_name, length, id_dim)
    ! Subroutine for creating netCDF dimensions more convenient:
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN) :: ncid
    CHARACTER(LEN=*),           INTENT(IN) :: dim_name
    INTEGER,                    INTENT(IN) :: length

    ! Output variables:
    INTEGER, INTENT(OUT)               :: id_dim

    CALL handle_error(nf90_def_dim(ncid,dim_name,length,id_dim))

  END SUBROUTINE create_dim
  SUBROUTINE create_int_var(    ncid, var_name, id_dims, id_var, long_name, units, missing_value)
    ! Subroutine for creating netCDF variables of type nf90_int more convenient:

    ! Input variables:
    INTEGER,                      INTENT(IN)  :: ncid
    CHARACTER(LEN=*),             INTENT(IN)  :: var_name
    INTEGER, DIMENSION(:),        INTENT(IN)  :: id_dims
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: long_name
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: units
    REAL(dp),           OPTIONAL, INTENT(IN)  :: missing_value

    ! Output variables:
    INTEGER,                      INTENT(OUT) :: id_var

    CALL handle_error(nf90_def_var(ncid,var_name,nf90_int,id_dims,id_var))
    IF(PRESENT(long_name))     CALL handle_error(nf90_put_att(ncid,id_var,'long_name',long_name))
    IF(PRESENT(units))         CALL handle_error(nf90_put_att(ncid,id_var,'units',units))
    IF(PRESENT(missing_value)) CALL handle_error(nf90_put_att(ncid,id_var,'missing_value',missing_value))

  END SUBROUTINE create_int_var
  SUBROUTINE create_single_var( ncid, var_name, id_dims, id_var, long_name, units, missing_value)
    ! Subroutine for creating netCDF variables of type nf90_FLOAT more convenient:

    ! Input variables:
    INTEGER,                      INTENT(IN)  :: ncid
    CHARACTER(LEN=*),             INTENT(IN)  :: var_name
    INTEGER, DIMENSION(:),        INTENT(IN)  :: id_dims
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: long_name
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: units
    REAL(dp),           OPTIONAL, INTENT(IN)  :: missing_value

    ! Output variables:
    INTEGER,                      INTENT(OUT) :: id_var

    CALL handle_error( nf90_def_var( ncid, var_name, nf90_float, id_dims, id_var))
    IF(PRESENT(long_name))     CALL handle_error( nf90_put_att( ncid, id_var ,'long_name'    ,long_name    ))
    IF(PRESENT(units))         CALL handle_error( nf90_put_att( ncid, id_var, 'units'        ,units        ))
    IF(PRESENT(missing_value)) CALL handle_error( nf90_put_att( ncid, id_var, 'missing_value',missing_value))

  END SUBROUTINE create_single_var
  SUBROUTINE create_double_var( ncid, var_name, id_dims, id_var, long_name, units, missing_value)
    ! Subroutine for creating netCDF variables of type nf90_DOUBLE more convenient:

    ! Input variables:
    INTEGER,                      INTENT(IN)  :: ncid
    CHARACTER(LEN=*),             INTENT(IN)  :: var_name
    INTEGER, DIMENSION(:),        INTENT(IN)  :: id_dims
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: long_name
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: units
    REAL(dp),           OPTIONAL, INTENT(IN)  :: missing_value

    ! Output variables:
    INTEGER,                      INTENT(OUT) :: id_var

    CALL handle_error(nf90_def_var(ncid,var_name,nf90_double,id_dims,id_var))
    IF(PRESENT(long_name))     CALL handle_error(nf90_put_att(ncid,id_var,'long_name',long_name))
    IF(PRESENT(units))         CALL handle_error(nf90_put_att(ncid,id_var,'units',units))
    IF(PRESENT(missing_value)) CALL handle_error(nf90_put_att(ncid,id_var,'missing_value',missing_value))

  END SUBROUTINE create_double_var
  SUBROUTINE inquire_dim( ncid, dim_name, dim_length, id_dim)
    ! Inquire the id of a dimension and return its length.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)  :: ncid
    CHARACTER(LEN=*),           INTENT(IN)  :: dim_name

    ! Output variables:
    INTEGER,                    INTENT(OUT) :: dim_length
    INTEGER,                    INTENT(OUT) :: id_dim

    CALL handle_error(nf90_inq_dimid(ncid,dim_name,id_dim))
    CALL handle_error(nf90_inquire_dimension(ncid, id_dim, len=dim_length))

  END SUBROUTINE inquire_dim
  SUBROUTINE inquire_int_var( ncid, var_name, id_dims, id_var)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_int.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: var_name
    INTEGER, DIMENSION(:),      INTENT(IN)    :: id_dims

    ! Output variables:
    INTEGER,                INTENT(OUT)   :: id_var

    ! Local variables:
    INTEGER                               :: xtype, ndims
    INTEGER, DIMENSION(nf90_max_var_dims) :: actual_id_dims

    CALL handle_error(nf90_inq_varid(ncid, var_name, id_var))
    CALL handle_error(nf90_inquire_variable(ncid, id_var, xtype=xtype,ndims=ndims,dimids=actual_id_dims))
    IF (xtype /= nf90_int) THEN
      CALL crash('Actual type of variable "' // TRIM( var_name) // '" is not nf90_int!')
    END IF
    IF (ndims /= SIZE( id_dims)) THEN
      CALL crash('Actual number of dimensions = {int_01} of variable "' // TRIM( var_name) // '" does not match required number of dimensions = {int_02}', &
        int_01 = ndims, int_02 = SIZE( id_dims))
    END IF
    IF (ANY( actual_id_dims( 1:ndims) /= id_dims)) THEN
      CALL crash('Actual dimensions of variable "' // TRIM( var_name) // '" does not match required dimensions!')
    END IF

  END SUBROUTINE inquire_int_var
  SUBROUTINE inquire_single_var( ncid, var_name, id_dims, id_var)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_DOUBLE.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: var_name
    INTEGER, DIMENSION(:),      INTENT(IN)    :: id_dims

    ! Output variables:
    INTEGER,                INTENT(OUT)   :: id_var

    ! Local variables:
    INTEGER                               :: xtype, ndims
    INTEGER, DIMENSION(nf90_max_var_dims) :: actual_id_dims

    CALL handle_error(nf90_inq_varid(ncid, var_name, id_var))
    CALL handle_error(nf90_inquire_variable(ncid, id_var, xtype=xtype,ndims=ndims,dimids=actual_id_dims))
    IF (xtype /= nf90_float) THEN
      CALL crash('Actual type of variable "' // TRIM( var_name) // '" is not nf90_float!')
    END IF
    IF (ndims /= SIZE( id_dims)) THEN
      CALL crash('Actual number of dimensions = {int_01} of variable "' // TRIM( var_name) // '" does not match required number of dimensions = {int_02}', &
        int_01 = ndims, int_02 = SIZE( id_dims))
    END IF
    IF (ANY( actual_id_dims( 1:ndims) /= id_dims)) THEN
      CALL crash('Actual dimensions of variable "' // TRIM( var_name) // '" does not match required dimensions!')
    END IF

  END SUBROUTINE inquire_single_var
  SUBROUTINE inquire_double_var( ncid, var_name, id_dims, id_var)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_DOUBLE.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: var_name
    INTEGER, DIMENSION(:),      INTENT(IN)    :: id_dims

    ! Output variables:
    INTEGER,                INTENT(OUT)   :: id_var

    ! Local variables:
    INTEGER                               :: xtype, ndims
    INTEGER, DIMENSION(nf90_max_var_dims) :: actual_id_dims

    CALL handle_error(nf90_inq_varid(ncid, var_name, id_var))
    CALL handle_error(nf90_inquire_variable(ncid, id_var, xtype=xtype,ndims=ndims,dimids=actual_id_dims))
    IF(xtype /= nf90_double) THEN
      CALL crash('Actual type of variable "' // TRIM( var_name) // '" is not nf90_double!')
    END IF
    IF (ndims /= SIZE( id_dims)) THEN
      CALL crash('Actual number of dimensions = {int_01} of variable "' // TRIM( var_name) // '" does not match required number of dimensions = {int_02}', &
        int_01 = ndims, int_02 = SIZE( id_dims))
    END IF
    IF (ANY( actual_id_dims( 1:ndims) /= id_dims)) THEN
      CALL crash('Actual dimensions of variable "' // TRIM( var_name) // '" does not match required dimensions!')
    END IF

  END SUBROUTINE inquire_double_var
  SUBROUTINE inquire_single_or_double_var( ncid, var_name, id_dims, id_var)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_FLOAT or nf90_DOUBLE.
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)    :: ncid
    CHARACTER(LEN=*),           INTENT(IN)    :: var_name
    INTEGER, DIMENSION(:),      INTENT(IN)    :: id_dims

    ! Output variables:
    INTEGER,                INTENT(OUT)   :: id_var

    ! Local variables:
    INTEGER                               :: xtype, ndims
    INTEGER, DIMENSION(nf90_max_var_dims) :: actual_id_dims

    CALL handle_error(nf90_inq_varid(ncid, var_name, id_var))
    CALL handle_error(nf90_inquire_variable(ncid, id_var, xtype=xtype,ndims=ndims,dimids=actual_id_dims))
    IF (xtype /= nf90_double .AND. xtype /= nf90_float) THEN
      CALL crash('Actual type of variable "' // TRIM( var_name) // '" is neither nf90_float nor nf90_double!')
    END IF
    IF (ndims /= SIZE( id_dims)) THEN
      CALL crash('Actual number of dimensions = {int_01} of variable "' // TRIM( var_name) // '" does not match required number of dimensions = {int_02}', &
        int_01 = ndims, int_02 = SIZE( id_dims))
    END IF
    IF (ANY( actual_id_dims( 1:ndims) /= id_dims)) THEN
      CALL crash('Actual dimensions of variable "' // TRIM( var_name) // '" does not match required dimensions!')
    END IF

  END SUBROUTINE inquire_single_or_double_var
  SUBROUTINE handle_error( stat, message)
    USE netcdf, ONLY: nf90_noerr, nf90_strerror
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN) :: stat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: message

    IF (stat /= nf90_noerr) THEN
      IF (PRESENT( message)) THEN
        CALL crash( message)
      ELSE
        CALL crash( 'netcdf error')
      END IF
    END IF

  END SUBROUTINE handle_error

  SUBROUTINE inquire_if_dim_exists_nocrash( ncid, dim_name, dim_exists)
    ! Inquire if a dimension called "dim_name" exists in the file identified by ncid
    IMPLICIT NONE

    ! Input variables:
    INTEGER,                    INTENT(IN)  :: ncid
    CHARACTER(LEN=*),           INTENT(IN)  :: dim_name

    ! Output variables:
    LOGICAL,                    INTENT(OUT) :: dim_exists

    ! Local variables:
    INTEGER                                 :: ncerr, id_dim

    ncerr = nf90_inq_dimid( ncid, dim_name, id_dim)

    dim_exists = ncerr == nf90_noerr

  END SUBROUTINE inquire_if_dim_exists_nocrash
  SUBROUTINE inquire_if_var_exists_nocrash( ncid, var_name, var_exists)
    ! Inquire the id of a variable and check that the dimensions of the variable match the dimensions given by the user and
    ! that the variable is of type nf90_FLOAT or nf90_DOUBLE.

    ! Input variables:
    INTEGER,                    INTENT(IN)  :: ncid
    CHARACTER(LEN=*),           INTENT(IN)  :: var_name

    ! Output variables:
    LOGICAL,                    INTENT(OUT) :: var_exists

    ! Local variables:
    INTEGER                                 :: ncerr, id_var

    ncerr = nf90_inq_varid( ncid, var_name, id_var)

    var_exists = ncerr == nf90_noerr

  END SUBROUTINE inquire_if_var_exists_nocrash

END MODULE netcdf_module
