MODULE netcdf_debug_module

! ===== Creating and writing to debug files =====
! ===============================================
!
! These routines create and write data to the NetCDF debug files.
! At pretty much any point in the model code, any data field
! on the grid (both A-, B-, and C-grids) can be immediately
! written to the NetCDF debug file for inspection. This is by
! far the most useful developer's tool available in IMAU_ICE.
!
! Example: say we want to see what's going on with the field ice%Hi_a:
!
! CALL save_variable_as_netcdf_dp_2D(ice%Hi_a,'ice%Hi_a') 
!
! And voila, there you go.
! NOTE: the routines in netcdf_debug_module are like all NetCDF routines
!       shared; if called by only a single process, the program will freeze.
!
! To use the debug fields, add these lines to the module list: 
!   USE netcdf_debug_module,             ONLY: save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, save_variable_as_netcdf_int_3D, &
!                                              save_variable_as_netcdf_dp_1D,  save_variable_as_netcdf_dp_2D,  save_variable_as_netcdf_dp_3D
!
! =====================================
!
! Additionally: single variables can now also be immediately saved as
!               NetCDF files. They will have no dimensions at all, so beware!
!
! ===== Preamble =====
! ====================

  ! Import basic functionality
!#include <petsc/finclude/petscksp.h>
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
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             transpose_dp_2D, transpose_dp_3D, transpose_int_2D, transpose_int_3D
  ! Import specific functionality
  USE data_types_module,               ONLY: type_model_region, type_grid
  USE netcdf,                          ONLY: NF90_INT, NF90_DOUBLE
  USE netcdf_basic_module,             ONLY: create_new_netcdf_file_for_writing, create_dimension, create_variable, &
                                             write_var_int_1D, write_var_dp_1D, &
                                             write_var_int_2D, write_var_dp_2D, &
                                             write_var_int_3D, write_var_dp_3D
  USE netcdf_output_module,            ONLY: setup_xy_grid_in_netcdf_file, add_zeta_dimension_to_file, &
                                             add_field_grid_int_2D_notime, add_field_grid_dp_2D_notime, &
                                             add_field_grid_dp_2D_monthly_notime, add_field_grid_dp_3D_notime


  IMPLICIT NONE

CONTAINS

! ===== Single-variable NetCDF files =====
! ========================================

  SUBROUTINE save_variable_as_netcdf_int_1D( d, field_name, region_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    INTEGER,  DIMENSION(:    ),      INTENT(IN)        :: d
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name
    CHARACTER(LEN=*), INTENT(IN),      OPTIONAL        :: region_name    
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_int_1D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: id_dim_n1
    INTEGER                                            :: id_var

     ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine if region_name needs to be left empty
    ! IF (PRESENT( region_name)) THEN
    !     filename = TRIM( C%output_dir) // TRIM( field_name) // TRIM(region_name) // '.nc'
    ! ELSE
        filename = TRIM( C%output_dir) // TRIM( field_name) // '.nc'
    ! END IF 

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Create dimensions
    CALL create_dimension( filename, 'n1', SIZE( d,1), id_dim_n1)

    ! Create variable
    CALL create_variable( filename, field_name, NF90_INT, (/ id_dim_n1 /), id_var)

    ! Write data
    CALL write_var_int_1D(  filename, id_var, d)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_int_1D

  SUBROUTINE save_variable_as_netcdf_int_2D( d, field_name, region_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    INTEGER,  DIMENSION(:,:  ),      INTENT(IN)        :: d
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name
    CHARACTER(LEN=*), INTENT(IN),      OPTIONAL        :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_int_2D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: id_dim_n1, id_dim_n2
    INTEGER                                            :: id_var
    INTEGER, DIMENSION(:,:  ),          POINTER        :: d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine if region_name needs to be left empty
    !IF (PRESENT( region_name)) THEN
    !    filename = TRIM( C%output_dir) // TRIM( field_name) // TRIM(region_name) // '.nc'
    !ELSE
        filename = TRIM( C%output_dir) // TRIM( field_name) // '.nc'
    !END IF

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Transpose so the orientation matches the help fields
    CALL allocate_shared_int_2D( SIZE( d,1), SIZE( d,2), d_grid, wd_grid)
    
    IF (par%master) THEN           
         d_grid = d
    END IF 

    CALL transpose_int_2D( d_grid, wd_grid )

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Create dimensions
    CALL create_dimension( filename, 'n1', SIZE( d_grid,1), id_dim_n1)
    CALL create_dimension( filename, 'n2', SIZE( d_grid,2), id_dim_n2)

    ! Create variable
    CALL create_variable( filename, field_name, NF90_INT, (/ id_dim_n1, id_dim_n2 /), id_var)

    ! Write data
    CALL write_var_int_2D(  filename, id_var, d_grid)

    ! Clean up after yourself
    CALL deallocate_shared(wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_int_2D

  SUBROUTINE save_variable_as_netcdf_int_3D( d, field_name, region_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    INTEGER,  DIMENSION(:,:,:),      INTENT(IN)        :: d
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name
    CHARACTER(LEN=*), INTENT(IN),      OPTIONAL        :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_int_3D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: id_dim_n1, id_dim_n2, id_dim_n3
    INTEGER                                            :: id_var
    INTEGER, DIMENSION(:,:,:  ),        POINTER        :: d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine if region_name needs to be left empty
    !IF (PRESENT( region_name)) THEN
    !    filename = TRIM( C%output_dir) // TRIM( field_name) // TRIM(region_name) // '.nc'
    !ELSE
        filename = TRIM( C%output_dir) // TRIM( field_name) // '.nc'
    !END IF

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Transpose so the orientation matches the help fields
    CALL allocate_shared_int_3D( SIZE( d,1), SIZE( d,2), SIZE( d,3), d_grid, wd_grid)

    IF (par%master) THEN
         d_grid = d
    END IF

    CALL transpose_int_3D( d_grid, wd_grid )

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Create dimensions
    CALL create_dimension( filename, 'n1', SIZE( d_grid,1), id_dim_n1)
    CALL create_dimension( filename, 'n2', SIZE( d_grid,2), id_dim_n2)
    CALL create_dimension( filename, 'n3', SIZE( d_grid,3), id_dim_n3)

    ! Create variable
    CALL create_variable( filename, field_name, NF90_INT, (/ id_dim_n1, id_dim_n2, id_dim_n3 /), id_var)

    ! Write data
    CALL write_var_int_3D(  filename, id_var, d_grid)
    
    ! Clean up after yourself
    CALL deallocate_shared(wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_int_3D

  SUBROUTINE save_variable_as_netcdf_dp_1D( d, field_name, region_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(:    ),      INTENT(IN)        :: d
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name
    CHARACTER(LEN=*), INTENT(IN),      OPTIONAL        :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_dp_1D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: id_dim_n1
    INTEGER                                            :: id_var

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine if region_name needs to be left empty
    ! IF (PRESENT( region_name)) THEN
    !     filename = TRIM( C%output_dir) // TRIM( field_name) // TRIM(region_name) // '.nc'
    ! ELSE
        filename = TRIM( C%output_dir) // TRIM( field_name) // '.nc'
    ! END IF

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Create dimensions
    CALL create_dimension( filename, 'n1', SIZE( d,1), id_dim_n1)

    ! Create variable
    CALL create_variable( filename, field_name, NF90_DOUBLE, (/ id_dim_n1 /), id_var)

    ! Write data
    CALL write_var_dp_1D(  filename, id_var, d)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_dp_1D

  SUBROUTINE save_variable_as_netcdf_dp_2D( d, field_name, region_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(:,:  ),      INTENT(IN)        :: d
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name
    CHARACTER(LEN=*),   OPTIONAL,    INTENT(IN)        :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_dp_2D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: id_dim_n1, id_dim_n2
    INTEGER                                            :: id_var
    REAL(dp), DIMENSION(:,:  ),         POINTER        :: d_grid
    INTEGER                                            :: wd_grid

    ! Add routine to path
    CALL init_routine( routine_name)
    CALL sync

    ! Determine if region_name needs to be left empty
    !IF (PRESENT( region_name)) THEN
    !  filename = TRIM( C%output_dir) // TRIM( field_name) // TRIM(region_name) // '.nc'
    !  PRINT*,('Wrong')
    !ELSE
      filename = TRIM( C%output_dir) // TRIM( field_name) // '.nc'
    !  PRINT*,('Right')
    !END IF

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Transpose so the orientation matches the help fields
    CALL allocate_shared_dp_2D( SIZE( d,1), SIZE( d,2), d_grid, wd_grid)

    IF (par%master) THEN
         d_grid = d
    END IF
    CALL sync 

    CALL transpose_dp_2D( d_grid, wd_grid )

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)
    
    ! Create dimensions
    CALL create_dimension( filename, 'n1', SIZE( d_grid,1), id_dim_n1)
    CALL create_dimension( filename, 'n2', SIZE( d_grid,2), id_dim_n2)

    ! Create variable
    CALL create_variable( filename, field_name, NF90_DOUBLE, (/ id_dim_n1, id_dim_n2 /), id_var)

    ! Write data
    CALL write_var_dp_2D(  filename, id_var, d_grid)

    ! Clean up after yourself
    CALL deallocate_shared(wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_dp_2D

  SUBROUTINE save_variable_as_netcdf_dp_3D( d, field_name, region_name)
    ! Save a single variable to a NetCDF file

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(:,:,:),      INTENT(IN)        :: d
    CHARACTER(LEN=*),                INTENT(IN)        :: field_name
    CHARACTER(LEN=*), INTENT(IN),      OPTIONAL        :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'save_variable_as_netcdf_dp_3D'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: id_dim_n1, id_dim_n2, id_dim_n3
    INTEGER                                            :: id_var
    REAL(dp), DIMENSION(:,:,:  ),       POINTER        :: d_grid
    INTEGER                                            :: wd_grid
    
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine if region_name needs to be left empty FJESSE: Succes!
    !IF (PRESENT( region_name)) THEN
    !    filename = TRIM( C%output_dir) // TRIM( field_name) // TRIM(region_name) // '.nc'
    !ELSE
        filename = TRIM( C%output_dir) // TRIM( field_name) // '.nc'
    !END IF

    ! Delete existing file
    IF (par%master) THEN
      INQUIRE( EXIST = file_exists, FILE = TRIM( filename))
      IF (file_exists) THEN
        CALL system('rm -f ' // filename)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync

    ! Transpose so the orientation matches the help fields
    CALL allocate_shared_dp_3D( SIZE( d,1), SIZE( d,2), SIZE( d,3), d_grid, wd_grid)

    IF (par%master) THEN
         d_grid = d
    END IF

    CALL transpose_dp_3D( d_grid, wd_grid )

    ! Create a new NetCDF file
    CALL create_new_netcdf_file_for_writing( filename)

    ! Create dimensions
    CALL create_dimension( filename, 'n1', SIZE( d_grid,1), id_dim_n1)
    CALL create_dimension( filename, 'n2', SIZE( d_grid,2), id_dim_n2)
    CALL create_dimension( filename, 'n3', SIZE( d_grid,3), id_dim_n3)

    ! Create variable
    CALL create_variable( filename, field_name, NF90_DOUBLE, (/ id_dim_n1, id_dim_n2, id_dim_n3 /), id_var)

    ! Write data
    CALL write_var_dp_3D(  filename, id_var, d_grid)

    ! Clean up after yourself
    CALL deallocate_shared(wd_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE save_variable_as_netcdf_dp_3D

  END MODULE netcdf_debug_module
