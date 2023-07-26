MODULE netcdf_extra_module

! This module contains routines for handling NetCDF files
! that don't conform to the typical x/y-grid, lon/lat-grid, or mesh,
! 2-D, 2-D monthly, or 3-D data fields.
!
! For example, the Laskar insolation solution.
!
! ===== Preamble =====
! ====================

  ! Import basic functionality
! #include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE netcdf_basic_module,             ONLY: nerr, field_name_options_x, field_name_options_y, field_name_options_zeta, field_name_options_z_ocean, &
                                             field_name_options_lon, field_name_options_lat, field_name_options_time, field_name_options_month, &
                                             field_name_options_Hi, field_name_options_Hb, field_name_options_Hs, field_name_options_dHb, &
                                             field_name_options_SL, field_name_options_Ti, read_var_dp_1D, &
                                             read_var_dp_3D, inquire_dim_multiple_options, inquire_var_multiple_options
  USE netcdf,                          ONLY: NF90_NOERR, NF90_OPEN, NF90_CLOSE, NF90_NOWRITE, NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, &
                                             NF90_INQ_VARID, NF90_INQUIRE_VARIABLE, NF90_MAX_VAR_DIMS, NF90_GET_VAR, &
                                             NF90_CREATE, NF90_NOCLOBBER, NF90_NETCDF4, NF90_ENDDEF, NF90_REDEF, NF90_DEF_DIM, NF90_DEF_VAR, &
                                             NF90_PUT_ATT, NF90_WRITE, NF90_INT, NF90_INT64, NF90_FLOAT, NF90_DOUBLE, NF90_PUT_VAR, NF90_UNLIMITED, &
                                             NF90_INQUIRE_ATTRIBUTE
  USE data_types_module,               ONLY: type_grid, type_grid_lonlat, type_ice_model, type_model_region, type_reference_geometry, &
                                             type_restart_data, type_forcing_data, type_climate_snapshot, type_ocean_snapshot_global, &
                                             type_SELEN_global, type_global_scalar_data, type_highres_ocean_data, &
                                             type_BIV_target_velocity, type_BIV_bed_roughness, type_sparse_matrix_CSR
  USE utilities_module,                ONLY: inverse_oblique_sg_projection
  IMPLICIT NONE

CONTAINS

  ! Insolation solution (e.g. Laskar 2004)
  SUBROUTINE inquire_insolation_data_file( forcing)

    IMPLICIT NONE

    ! Output variable
    TYPE(type_forcing_data),             INTENT(INOUT) :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inquire_insolation_data_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire dimensions
    CALL inquire_dim_multiple_options( forcing%netcdf_ins%filename, field_name_options_time , &
      forcing%netcdf_ins%id_dim_time , dim_length = forcing%ins_nyears)
    CALL inquire_dim_multiple_options( forcing%netcdf_ins%filename, field_name_options_month, &
      forcing%netcdf_ins%id_dim_month)
    CALL inquire_dim_multiple_options( forcing%netcdf_ins%filename, field_name_options_lat  , &
      forcing%netcdf_ins%id_dim_lat  , dim_length = forcing%ins_nlat)

    ! Inquire variables
    CALL inquire_var_multiple_options( forcing%netcdf_ins%filename, field_name_options_time , &
      forcing%netcdf_ins%id_var_time )
    CALL inquire_var_multiple_options( forcing%netcdf_ins%filename, field_name_options_month, &
      forcing%netcdf_ins%id_var_month)
    CALL inquire_var_multiple_options( forcing%netcdf_ins%filename, field_name_options_lat  , &
      forcing%netcdf_ins%id_var_lat  )
    CALL inquire_var_multiple_options( forcing%netcdf_ins%filename, 'Q_TOA', &
      forcing%netcdf_ins%id_var_Q_TOA)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_insolation_data_file

  SUBROUTINE read_insolation_data_file_timeframes( forcing, ti0, ti1, ins_Q_TOA0, ins_Q_TOA1)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_forcing_data),             INTENT(INOUT) :: forcing
    INTEGER,                             INTENT(IN)    :: ti0, ti1
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: ins_Q_TOA0, ins_Q_TOA1

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_insolation_data_file_timeframes'
    INTEGER                                            :: mi, li
    REAL(dp), DIMENSION(:,:,:), POINTER                ::  Q0_with_time,  Q1_with_time
    INTEGER                                            :: wQ0_with_time, wQ1_with_time

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Temporary memory to store the data read from the netCDF file
    CALL allocate_shared_dp_3D( 1, 12, forcing%ins_nlat, Q0_with_time, wQ0_with_time)
    CALL allocate_shared_dp_3D( 1, 12, forcing%ins_nlat, Q1_with_time, wQ1_with_time)
     
    CALL read_var_dp_3D( forcing%netcdf_ins%filename, forcing%netcdf_ins%id_var_Q_TOA, Q0_with_time, &
      start = (/ ti0, 1, 1 /), count = (/ 1, 12, forcing%ins_nlat /) )
    
    CALL read_var_dp_3D( forcing%netcdf_ins%filename, forcing%netcdf_ins%id_var_Q_TOA, Q1_with_time, &
      start = (/ ti1, 1, 1 /), count = (/ 1, 12, forcing%ins_nlat /) )

    ! Remove the time dimension
    IF (par%master) THEN
      DO mi = 1, 12
      DO li = 1, forcing%ins_nlat
        ins_Q_TOA0( li,mi) = Q0_with_time( 1,mi,li)
        ins_Q_TOA1( li,mi) = Q1_with_time( 1,mi,li)
      END DO
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Clean up temporary memory
    CALL deallocate_shared( wQ0_with_time)
    CALL deallocate_shared( wQ1_with_time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_insolation_data_file_timeframes

  SUBROUTINE read_insolation_data_file_time_lat( forcing)
    ! Read only the time and latitude variables from the insolation file,
    !so that later on we know which time frames to read and how to map them.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_forcing_data),             INTENT(INOUT) :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_insolation_data_file_time_lat'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Read the data
    CALL read_var_dp_1D( forcing%netcdf_ins%filename, forcing%netcdf_ins%id_var_time, forcing%ins_time)
    CALL read_var_dp_1D( forcing%netcdf_ins%filename, forcing%netcdf_ins%id_var_lat , forcing%ins_lat )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_insolation_data_file_time_lat

END MODULE netcdf_extra_module
