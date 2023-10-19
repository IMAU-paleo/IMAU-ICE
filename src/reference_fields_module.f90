MODULE reference_fields_module

  ! Contains the routines for setting up the "PD", "init", and "topo" reference data fields.

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE data_types_module,               ONLY: type_grid, type_reference_geometry, type_restart_data
  USE parameters_module,               ONLY: seawater_density, ice_density, sec_per_year, pi
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             map_square_to_square_cons_2nd_order_2D, map_square_to_square_cons_2nd_order_3D, &
                                             transpose_dp_2D, transpose_dp_3D, remove_Lake_Vostok, surface_elevation, &
                                             smooth_Gaussian_2D, is_floating, deallocate_grid
  USE netcdf_input_module,             ONLY: read_field_from_xy_file_2D, read_field_from_file_2D
  USE netcdf_basic_module,             ONLY: field_name_options_Hi, field_name_options_Hb, field_name_options_Hs
  USE netcdf_debug_module,             ONLY: save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, save_variable_as_netcdf_int_3D, &
                                             save_variable_as_netcdf_dp_1D,  save_variable_as_netcdf_dp_2D,  save_variable_as_netcdf_dp_3D

  IMPLICIT NONE

CONTAINS

  ! Initialise all three reference geometries
  SUBROUTINE initialise_reference_geometries( grid, refgeo_init, refgeo_PD, refgeo_GIAeq, region_name)
    ! Initialise all three reference geometries

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo_init
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo_PD
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo_GIAeq
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometries'
    CHARACTER(LEN=256)                            :: choice_refgeo_init, choice_refgeo_PD, choice_refgeo_GIAeq
    CHARACTER(LEN=256)                            :: filename_refgeo_init, filename_refgeo_PD, filename_refgeo_GIAeq
    REAL(dp)                                      :: time_to_restart_from

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine parameters for this region
    IF     (region_name == 'NAM') THEN
      choice_refgeo_init    = C%choice_refgeo_init_NAM
      choice_refgeo_PD      = C%choice_refgeo_PD_NAM
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_NAM
      filename_refgeo_init  = C%filename_refgeo_init_NAM
      filename_refgeo_PD    = C%filename_refgeo_PD_NAM
      filename_refgeo_GIAeq = C%filename_refgeo_GIAeq_NAM
      time_to_restart_from  = C%time_to_restart_from_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_refgeo_init    = C%choice_refgeo_init_EAS
      choice_refgeo_PD      = C%choice_refgeo_PD_EAS
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_EAS
      filename_refgeo_init  = C%filename_refgeo_init_EAS
      filename_refgeo_PD    = C%filename_refgeo_PD_EAS
      filename_refgeo_GIAeq = C%filename_refgeo_GIAeq_EAS
      time_to_restart_from  = C%time_to_restart_from_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_refgeo_init    = C%choice_refgeo_init_GRL
      choice_refgeo_PD      = C%choice_refgeo_PD_GRL
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_GRL
      filename_refgeo_init  = C%filename_refgeo_init_GRL
      filename_refgeo_PD    = C%filename_refgeo_PD_GRL
      filename_refgeo_GIAeq = C%filename_refgeo_GIAeq_GRL
      time_to_restart_from  = C%time_to_restart_from_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_refgeo_init    = C%choice_refgeo_init_ANT
      choice_refgeo_PD      = C%choice_refgeo_PD_ANT
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_ANT
      filename_refgeo_init  = C%filename_refgeo_init_ANT
      filename_refgeo_PD    = C%filename_refgeo_PD_ANT
      filename_refgeo_GIAeq = C%filename_refgeo_GIAeq_ANT
      time_to_restart_from  = C%time_to_restart_from_ANT
    END IF

    ! Initial ice-sheet geometry
    ! ==========================

    IF     (choice_refgeo_init == 'idealised') THEN
      IF (par%master) WRITE(0,*) '  Initialising initial         reference geometry from idealised case "', TRIM(C%choice_refgeo_init_idealised), '"...'
      CALL initialise_reference_geometry_idealised( grid, refgeo_init, C%choice_refgeo_init_idealised)
    ELSEIF (choice_refgeo_init == 'realistic') THEN
      IF (par%master) WRITE(0,*) '  Initialising initial         reference geometry from file ', TRIM( filename_refgeo_init), '...'
      CALL initialise_reference_geometry_from_file( grid, refgeo_init, filename_refgeo_init, region_name)
    ELSEIF (choice_refgeo_init == 'restart') THEN
      IF (par%master) WRITE(0,*) '  Initialising initial         reference geometry from restart file ', TRIM( filename_refgeo_init), '...'
      CALL initialise_reference_geometry_from_restart_file( grid, refgeo_init, filename_refgeo_init, region_name, time_to_restart_from)
    ELSE
      CALL crash('unknown choice_refgeo_init "' // TRIM( choice_refgeo_init) // '"!')
    END IF

    ! Present-day ice-sheet geometry
    ! ==============================

    IF     (choice_refgeo_PD == 'idealised') THEN
      IF (par%master) WRITE(0,*) '  Initialising present-day     reference geometry from idealised case "', TRIM(C%choice_refgeo_PD_idealised), '"...'
      CALL initialise_reference_geometry_idealised( grid, refgeo_PD, C%choice_refgeo_PD_idealised)
    ELSEIF (choice_refgeo_PD == 'realistic') THEN
      IF (par%master) WRITE(0,*) '  Initialising present-day     reference geometry from file ', TRIM( filename_refgeo_PD), '...'
      CALL initialise_reference_geometry_from_file( grid, refgeo_PD, filename_refgeo_PD, region_name)
    ELSE
      CALL crash('unknown choice_refgeo_PD "' // TRIM( choice_refgeo_PD) // '"!')
    END IF

    ! GIA equilibrium ice-sheet geometry
    ! ==================================

    IF     (choice_refgeo_GIAeq == 'idealised') THEN
      IF (par%master) WRITE(0,*) '  Initialising GIA equilibrium reference geometry from idealised case "', TRIM(C%choice_refgeo_GIAeq_idealised), '"...'
      CALL initialise_reference_geometry_idealised( grid, refgeo_GIAeq, C%choice_refgeo_GIAeq_idealised)
    ELSEIF (choice_refgeo_GIAeq == 'realistic') THEN
      IF (par%master) WRITE(0,*) '  Initialising GIA equilibrium reference geometry from file ', TRIM( filename_refgeo_GIAeq), '...'
      CALL initialise_reference_geometry_from_file( grid, refgeo_GIAeq, filename_refgeo_GIAeq, region_name)
    ELSE
      CALL crash('unknown choice_refgeo_GIAeq "' // TRIM( choice_refgeo_GIAeq) // '"!')
    END IF

    ! Smooth input geometry (bed and ice)
    IF (C%do_smooth_geometry) THEN
      CALL smooth_model_geometry( grid, refgeo_PD%Hi,    refgeo_PD%Hb,    refgeo_PD%Hs   )
      CALL smooth_model_geometry( grid, refgeo_init%Hi,  refgeo_init%Hb,  refgeo_init%Hs )
      CALL smooth_model_geometry( grid, refgeo_GIAeq%Hi, refgeo_GIAeq%Hb, refgeo_GIAeq%Hs)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometries

  ! Initialise a reference geometry with data from a (timeless) NetCDF file (e.g. BedMachine)
  SUBROUTINE initialise_reference_geometry_from_file( grid, refgeo, filename_refgeo, region_name)
    ! Initialise a reference geometry with data from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    CHARACTER(LEN=256),             INTENT(IN)    :: filename_refgeo
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_from_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Read data from input file
    CALL read_reference_geometry_file( refgeo, filename_refgeo, region_name)
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( refgeo%Hi_raw, 'refgeo%Hi_raw')
    CALL check_for_NaN_dp_2D( refgeo%Hb_raw, 'refgeo%Hb_raw')
    CALL check_for_NaN_dp_2D( refgeo%Hs_raw, 'refgeo%Hs_raw')

    ! Remove Lake Vostok from Antarctica (because it's annoying)
    IF (region_name == 'ANT'.AND. C%remove_Lake_Vostok) THEN
      CALL remove_Lake_Vostok( refgeo%grid%x, refgeo%grid%y, refgeo%Hi_raw, refgeo%Hb_raw, refgeo%Hs_raw)
    END IF

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, refgeo%Hi, refgeo%wHi)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, refgeo%Hb, refgeo%wHb)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, refgeo%Hs, refgeo%wHs)

    ! Map (transposed) raw data to the model grid
    CALL map_square_to_square_cons_2nd_order_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%grid%x, refgeo%grid%y, grid%nx, grid%ny, grid%x, grid%y, refgeo%Hi_raw, refgeo%Hi)
    CALL map_square_to_square_cons_2nd_order_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%grid%x, refgeo%grid%y, grid%nx, grid%ny, grid%x, grid%y, refgeo%Hb_raw, refgeo%Hb)
    CALL map_square_to_square_cons_2nd_order_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%grid%x, refgeo%grid%y, grid%nx, grid%ny, grid%x, grid%y, refgeo%Hs_raw, refgeo%Hs)

    ! Deallocate raw data
    CALL deallocate_shared( refgeo%wHi_raw)
    CALL deallocate_shared( refgeo%wHb_raw)
    CALL deallocate_shared( refgeo%wHs_raw)
    CALL deallocate_grid(   refgeo%grid   )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_from_file

  ! Initialise a reference geometry with data from a previous simulation's restart file
  SUBROUTINE initialise_reference_geometry_from_restart_file( grid, refgeo, filename_refgeo, region_name, time_to_restart_from)
    ! Initialise a reference geometry with data from a previous simulation's restart file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    CHARACTER(LEN=256),             INTENT(IN)    :: filename_refgeo
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    REAL(dp),                       INTENT(IN)    :: time_to_restart_from

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_from_restart_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, refgeo%Hi, refgeo%wHi)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, refgeo%Hb, refgeo%wHb)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, refgeo%Hs, refgeo%wHs)

    ! Read data from input file
    CALL read_restart_file_geometry( filename_refgeo, grid, refgeo, region_name, time_to_restart_from)
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_2D( refgeo%Hi, 'refgeo%Hi')
    CALL check_for_NaN_dp_2D( refgeo%Hb, 'refgeo%Hb')
    CALL check_for_NaN_dp_2D( refgeo%Hs, 'refgeo%Hs')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_from_restart_file
  SUBROUTINE adapt_initial_geometry_from_restart_file( grid, refgeo_PD, refgeo_init, filename_refgeo_init, region_name, time_to_restart_from)
    ! This subroutine that is not being used anymore.

    ! Restarting a run can mean the initial bedrock is deformed, which should be accounted for.
    ! Also, the current model resolution might be higher than that which was used to generate
    ! the restart file. Both fo these problems are solved by adding the restart dHb to the PD Hb.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(IN)    :: refgeo_PD
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo_init
    CHARACTER(LEN=256),             INTENT(IN)    :: filename_refgeo_init
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    REAL(dp),                       INTENT(IN)    :: time_to_restart_from

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'adapt_initial_geometry_from_restart_file'
    INTEGER                                       :: i,j
    TYPE(type_restart_data)                       :: restart
    REAL(dp), DIMENSION(:,:  ), POINTER           ::  dHb,  SL
    INTEGER                                       :: wdHb, wSL
    REAL(dp)                                      :: Hs, Hs_max_float

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate the additional input fields
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, SL , wSL )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dHb, wdHb)

    ! Read data from input file
    CALL read_field_from_file_2D(    filename_refgeo_init, 'SL',  grid,  SL,   region_name, time_to_restart_from)
    CALL read_field_from_file_2D(    filename_refgeo_init, 'dHb', grid,  dHb,  region_name, time_to_restart_from)
    CALL sync

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Assume the mapped surface elevation is correct
      Hs = surface_elevation( refgeo_init%Hi( j,i), refgeo_init%Hb( j,i), SL( j,i))

      ! Define bedrock as PD bedrock + restarted deformation
      refgeo_init%Hb( j,i) = refgeo_PD%Hb( j,i) + dHb( j,i)

      ! Surface elevation cannot be below bedrock
      Hs = MAX( Hs, refgeo_init%Hb( j,i))

      ! Define ice thickness
      Hs_max_float = refgeo_init%Hb( j,i) + MAX( 0._dp, (SL( j,i) - refgeo_init%Hb( j,i)) * seawater_density / ice_density)
      IF (Hs > Hs_max_float) THEN
        ! Ice here must be grounded
        refgeo_init%Hi( j,i) = Hs - refgeo_init%Hb( j,i)
      ELSE
        ! Ice here must be floating
        refgeo_init%Hi( j,i) = MIN( Hs - refgeo_init%Hb( j,i), MAX( 0._dp, (Hs - SL( j,i))) / (1._dp - ice_density / seawater_density))
      END IF

      ! Recalculate surface elevation
      refgeo_init%Hs( j,i) = surface_elevation( refgeo_init%Hi( j,i), refgeo_init%Hb( j,i), SL( j,i))

    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wSL    )
    CALL deallocate_shared( wdHb   )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE adapt_initial_geometry_from_restart_file

  ! Initialise a reference geometry according to an idealised world
  SUBROUTINE initialise_reference_geometry_idealised( grid, refgeo, choice_refgeo_idealised)
    ! Initialise a reference geometry according to an idealised world

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    CHARACTER(LEN=256),             INTENT(IN)    :: choice_refgeo_idealised

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, refgeo%Hi, refgeo%wHi)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, refgeo%Hb, refgeo%wHb)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, refgeo%Hs, refgeo%wHs)

    IF     (choice_refgeo_idealised == 'flatearth') THEN
      ! Simply a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments
      CALL initialise_reference_geometry_idealised_flatearth( grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'Halfar') THEN
      ! The Halfar dome solution at t = 0
      CALL initialise_reference_geometry_idealised_Halfar( grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'Bueler') THEN
      ! The Bueler dome solution at t = 0
      CALL initialise_reference_geometry_idealised_Bueler( grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'SSA_icestream') THEN
      ! The SSA_icestream infinite slab on a flat slope
      CALL initialise_reference_geometry_idealised_SSA_icestream( grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'MISMIP_mod') THEN
      ! The MISMIP_mod cone-shaped island
      CALL initialise_reference_geometry_idealised_MISMIP_mod( grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_A') THEN
      ! The ISMIP-HOM A bumpy slope
      CALL initialise_reference_geometry_idealised_ISMIP_HOM_A( grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_B') THEN
      ! The ISMIP-HOM B bumpy slope
      CALL initialise_reference_geometry_idealised_ISMIP_HOM_B( grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_C' .OR. &
            choice_refgeo_idealised == 'ISMIP_HOM_D') THEN
      ! The ISMIP-HOM C/D bumpy slope
      CALL initialise_reference_geometry_idealised_ISMIP_HOM_CD( grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_E') THEN
      ! The ISMIP-HOM E Glacier d'Arolla geometry
      CALL initialise_reference_geometry_idealised_ISMIP_HOM_E( grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_F') THEN
      ! The ISMIP-HOM A bumpy slope
      CALL initialise_reference_geometry_idealised_ISMIP_HOM_F( grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'MISMIP+') THEN
      ! The MISMIP+ fjord geometry
      CALL initialise_reference_geometry_idealised_MISMIPplus( grid, refgeo)
    ELSE
      CALL crash('unknown choice_refgeo_idealised "' // TRIM( choice_refgeo_idealised) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised
  SUBROUTINE initialise_reference_geometry_idealised_flatearth(     grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! Simply a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_flatearth'
    INTEGER                                       :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hi( j,i) = 0._dp
      refgeo%Hb( j,i) = 0._dp
      refgeo%Hs( j,i) = 0._dp
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised_flatearth
  SUBROUTINE initialise_reference_geometry_idealised_Halfar(        grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The Halfar dome solution at t = 0

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_Halfar'
    INTEGER                                       :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hi( j,i) = Halfar_solution( grid%x( i), grid%y( j), C%start_time_of_run)
      refgeo%Hb( j,i) = 0._dp
      refgeo%Hs( j,i) = refgeo%Hi( j,i)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised_Halfar
  SUBROUTINE initialise_reference_geometry_idealised_Bueler(        grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The Bueler dome solution at t = 0

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_Bueler'
    INTEGER                                       :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hi( j,i) = Bueler_solution( grid%x( i), grid%y( j), C%start_time_of_run)
      refgeo%Hb( j,i) = 0._dp
      refgeo%Hs( j,i) = refgeo%Hi( j,i)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised_Bueler
  SUBROUTINE initialise_reference_geometry_idealised_SSA_icestream( grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The SSA_icestream infinite slab on a flat slope

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_SSA_icestream'
    INTEGER                                       :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hi( j,i) = C%SSA_icestream_H
      refgeo%Hb( j,i) = -C%SSA_icestream_tantheta * grid%x( i)
      refgeo%Hs( j,i) = refgeo%Hb( j,i) + refgeo%Hi( j,i)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised_SSA_icestream
  SUBROUTINE initialise_reference_geometry_idealised_MISMIP_mod(    grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The MISMIP_mod cone-shaped island

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_MISMIP_mod'
    INTEGER                                       :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hi( j,i) = 100._dp
      refgeo%Hb( j,i) = 720._dp - 778.5_dp * SQRT( grid%x(i)**2 + grid%y(j)**2)/ 750000._dp
      refgeo%Hs( j,i) = refgeo%Hb( j,i) + refgeo%Hi( j,i)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised_MISMIP_mod
  SUBROUTINE initialise_reference_geometry_idealised_ISMIP_HOM_A(   grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM A bumpy slope

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_ISMIP_HOM_A'
    INTEGER                                       :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hs( j,i) = 2000._dp - grid%x( i) * TAN( 0.5_dp * pi / 180._dp)
      refgeo%Hb( j,i) = refgeo%Hs( j,i) - 1000._dp + 500._dp * SIN( grid%x( i) * 2._dp * pi / C%ISMIP_HOM_L) * SIN( grid%y( j) * 2._dp * pi / C%ISMIP_HOM_L)
      refgeo%Hi( j,i) = refgeo%Hs( j,i) - refgeo%Hb( j,i)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised_ISMIP_HOM_A
  SUBROUTINE initialise_reference_geometry_idealised_ISMIP_HOM_B(   grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM B bumpy slope

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_ISMIP_HOM_B'
    INTEGER                                       :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hs( j,i) = 2000._dp - grid%x( i) * TAN( 0.5_dp * pi / 180._dp)
      refgeo%Hb( j,i) = refgeo%Hs( j,i) - 1000._dp + 500._dp * SIN( grid%x(i) * 2._dp * pi / C%ISMIP_HOM_L)
      refgeo%Hi( j,i) = refgeo%Hs( j,i) - refgeo%Hb( j,i)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised_ISMIP_HOM_B
  SUBROUTINE initialise_reference_geometry_idealised_ISMIP_HOM_CD(  grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM C/D bumpy slope

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_ISMIP_HOM_CD'
    INTEGER                                       :: i,j

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
        refgeo%Hs( j,i) = 2000._dp - grid%x( i) * TAN( 0.1_dp * pi / 180._dp)
        refgeo%Hb( j,i) = refgeo%Hs( j,i) - 1000._dp
        refgeo%Hi( j,i) = refgeo%Hs( j,i) - refgeo%Hb( j,i)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised_ISMIP_HOM_CD
  SUBROUTINE initialise_reference_geometry_idealised_ISMIP_HOM_E(   grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM E Glacier d'Arolla geometry

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_ISMIP_HOM_E'
    INTEGER                                       :: i,j
    REAL(dp)                                      :: x,Hs,Hb
    INTEGER                                       :: ios,slides

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Read data from external file
    IF (par%master) THEN

      OPEN( UNIT = 1337, FILE=C%ISMIP_HOM_E_Arolla_filename, ACTION = 'READ')
      DO i = 1, 51
        READ( UNIT = 1337, FMT=*, IOSTAT=ios) x, Hb, Hs, slides
        DO j = 1, grid%ny
          refgeo%Hb( j,i) = Hb
          refgeo%Hi( j,i) = Hs - Hb
          refgeo%Hs( j,i) = Hs
        END DO
        IF (ios /= 0) THEN
          CALL crash('length of text file "' // TRIM(C%ISMIP_HOM_E_Arolla_filename) // '" should be 51 lines!')
        END IF
      END DO
      CLOSE( UNIT  = 1337)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised_ISMIP_HOM_E
  SUBROUTINE initialise_reference_geometry_idealised_ISMIP_HOM_F(   grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM A bumpy slope

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_ISMIP_HOM_F'
    INTEGER                                       :: i,j

    REAL(dp), PARAMETER                           :: H0    = 1000._dp
    REAL(dp), PARAMETER                           :: a0    = 100._dp
    REAL(dp), PARAMETER                           :: sigma = 10000._dp

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hs( j,i) = 5000._dp - grid%x( i) * TAN( 3._dp * pi / 180._dp)
      refgeo%Hb( j,i) = refgeo%Hs( j,i) - H0 + a0 * EXP( -(grid%x( i)**2 + grid%y( j)**2) / sigma**2)
      refgeo%Hi( j,i) = refgeo%Hs( j,i) - refgeo%Hb( j,i)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised_ISMIP_HOM_F
  SUBROUTINE initialise_reference_geometry_idealised_MISMIPplus(    grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The MISMIpplus fjord geometry

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_MISMIPplus'
    INTEGER                                       :: i,j
    REAL(dp)                                      :: x,y,xtilde,Bx,By
    REAL(dp), PARAMETER                           :: B0     = -150._dp
    REAL(dp), PARAMETER                           :: B2     = -728.8_dp
    REAL(dp), PARAMETER                           :: B4     = 343.91_dp
    REAL(dp), PARAMETER                           :: B6     = -50.57_dp
    REAL(dp), PARAMETER                           :: xbar   = 300000._dp
    REAL(dp), PARAMETER                           :: fc     = 4000._dp
    REAL(dp), PARAMETER                           :: dc     = 500._dp
    REAL(dp), PARAMETER                           :: wc     = 24000._dp
    REAL(dp), PARAMETER                           :: zbdeep = -720._dp

    ! Add routine to path
    CALL init_routine( routine_name)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      x = grid%x( i) + 400000._dp
      y = -40000._dp +  80000._dp * REAL(j-1,dp) / REAL(grid%ny-1,dp)
      xtilde = x / xbar
      Bx = B0 + (B2 * xtilde**2._dp) + (B4 * xtilde**4._dp) + (B6 * xtilde**6._dp)
      By = (dc / (1 + EXP(-2._dp*(y - wc)/fc))) + &
           (dc / (1 + EXP( 2._dp*(y + wc)/fc)))
      refgeo%Hi( j,i) = 100._dp
      refgeo%Hb( j,i) = MAX( Bx + By, zbdeep)
      refgeo%Hs( j,i) = surface_elevation( refgeo%Hi( j,i), refgeo%Hb( j,i), 0._dp)
    END DO
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_reference_geometry_idealised_MISMIPplus

  ! Apply some light smoothing to the initial geometry to improve numerical stability
  SUBROUTINE smooth_model_geometry( grid, Hi, Hb, Hs)
    ! Apply some light smoothing to the initial geometry to improve numerical stability

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: Hi, Hb, Hs

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_model_geometry'
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  Hb_old,  dHb
    INTEGER                                            :: wHb_old, wdHb
    REAL(dp)                                           :: r_smooth

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Smooth with a 2-D Gaussian filter with a standard deviation of 1/2 grid cell
    r_smooth = grid%dx * C%r_smooth_geometry

    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, Hb_old, wHb_old)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dHb,    wdHb   )

    ! Calculate original surface elevation
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      Hs( j,i) = surface_elevation( Hi( j,i), Hb( j,i), 0._dp)
    END DO
    END DO
    CALL sync

    ! Store the unsmoothed bed topography so we can determine the smoothing anomaly later
    Hb_old( :,grid%i1:grid%i2) = Hb( :,grid%i1:grid%i2)
    CALL sync

    ! Apply smoothing to the bed topography
    CALL smooth_Gaussian_2D( grid, Hb, r_smooth)

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Calculate the smoothing anomaly
      dHb( j,i) = Hb( j,i) - Hb_old( j,i)

      IF (.NOT. is_floating( Hi( j,i), Hb( j,i), 0._dp) .AND. Hi( j,i) > 0._dp) THEN

        ! Correct the ice thickness so the ice surface remains unchanged (only relevant for floating ice)
        Hi( j,i) = Hi( j,i) - dHb( j,i)

        ! Don't allow negative ice thickness
        Hi( j,i) = MAX(0._dp, Hi( j,i))

        ! Correct the surface elevation for this if necessary
        Hs( j,i) = Hi( j,i) + MAX( 0._dp - ice_density / seawater_density * Hi( j,i), Hb( j,i))

      END IF

    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wHb_old)
    CALL deallocate_shared( wdHb   )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE smooth_model_geometry

  ! Analytical solutions used to initialise some benchmark experiments
  FUNCTION Halfar_solution( x, y, t) RESULT(H)
    ! Describes an ice-sheet at time t (in years) conforming to the Halfar similarity function
    ! with dome thickness H0 and margin radius R0 at t0. Used to initialise the model
    ! for the Halfar solution test run

    IMPLICIT NONE

    ! Input variables
    REAL(dp), INTENT(IN) :: x  ! x coordinate [m]
    REAL(dp), INTENT(IN) :: y  ! y coordinate [m]
    REAL(dp), INTENT(IN) :: t  ! Time from t0 [years]

    ! Result
    REAL(dp)             :: H  ! Ice thickness at [x,y] at t=0 [m]

    ! Local variables
    REAL(dp) :: A_flow, rho, g, Gamma, t0, r, f1, f2, f3, tp

    REAL(dp), PARAMETER :: H0 = 5000._dp   ! Ice dome thickness at t=0 [m]
    REAL(dp), PARAMETER :: R0 = 300000._dp ! Ice margin radius  at t=0 [m]

    A_flow  = 1E-16_dp
    rho     = 910._dp
    g       = 9.81_dp

    Gamma = (2._dp / 5._dp) * (A_flow / sec_per_year) * (rho * g)**3._dp
    t0 = 1._dp / (18._dp * Gamma) * (7._dp/4._dp)**3._dp * (R0**4._dp)/(H0**7._dp)

    tp = (t * sec_per_year) + t0

    r = SQRT(x**2._dp + y**2._dp)

    f1 = (t0/tp)**(1._dp/9._dp)
    f2 = (t0/tp)**(1._dp/18._dp)
    f3 = (r/R0)

    H = H0 * f1 * MAX(0._dp, (1._dp - (f2*f3)**(4._dp/3._dp)))**(3._dp/7._dp)

  END FUNCTION Halfar_solution
  FUNCTION Bueler_solution( x, y, t) RESULT(H)
    ! Describes an ice-sheet at time t (in years) conforming to the Bueler solution
    ! with dome thickness H0 and margin radius R0 at t0, with a surface mass balance
    ! determined by lambda. Used to intialise the model for the Bueler solution test run

    IMPLICIT NONE

    ! Input variables
    REAL(dp), INTENT(IN) :: x       ! x coordinate [m]
    REAL(dp), INTENT(IN) :: y       ! y coordinate [m]
    REAL(dp), INTENT(IN) :: t       ! Time from t0 [years]

    ! Result
    REAL(dp)             :: H  ! Ice thickness at [x,y] at t=0 [m]

    ! Local variables
    REAL(dp) :: A_flow, rho, g, n, alpha, beta, Gamma, f1, f2, t0, tp, f3, f4

    REAL(dp), PARAMETER :: H0     = 3000._dp    ! Ice dome thickness at t=0 [m]
    REAL(dp), PARAMETER :: R0     = 500000._dp  ! Ice margin radius  at t=0 [m]
    REAL(dp), PARAMETER :: lambda = 5.0_dp      ! Mass balance parameter

    A_flow  = 1E-16_dp
    rho     = 910._dp
    g       = 9.81_dp
    n       = 3._dp

    alpha = (2._dp - (n+1._dp)*lambda) / ((5._dp*n)+3._dp)
    beta  = (1._dp + ((2._dp*n)+1._dp)*lambda) / ((5._dp*n)+3._dp)
    Gamma = 2._dp/5._dp * (A_flow/sec_per_year) * (rho * g)**n

    f1 = ((2._dp*n)+1)/(n+1._dp)
    f2 = (R0**(n+1._dp))/(H0**((2._dp*n)+1._dp))
    t0 = (beta / Gamma) * (f1**n) * f2

    !tp = (t * sec_per_year) + t0; % Acutal equation needs t in seconds from zero , but we want to supply t in years from t0
    tp = t * sec_per_year

    f1 = (tp / t0)**(-alpha)
    f2 = (tp / t0)**(-beta)
    f3 = SQRT( (x**2._dp) + (y**2._dp) )/R0
    f4 = MAX(0._dp, 1._dp - (f2*f3)**((n+1._dp)/n))
    H = H0 * f1 * f4**(n/((2._dp*n)+1._dp))

    !M = (lambda / tp) * H * sec_per_year

  END FUNCTION Bueler_solution

  SUBROUTINE read_reference_geometry_file( refgeo, filename_refgeo, region_name)
    ! Read reference geometry dat from a NetCDF file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_reference_geometry), INTENT(INOUT) :: refgeo
    CHARACTER(LEN=256),            INTENT(IN)    :: filename_refgeo
    CHARACTER(LEN=3),              INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'read_reference_geometry_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Read the data
    CALL read_field_from_xy_file_2D( filename_refgeo, field_name_options_Hi, region_name, refgeo%grid, refgeo%Hi_raw, refgeo%wHi_raw)
    CALL deallocate_grid(refgeo%grid) ! Deallocate the grid, as it is automatically being allocated. Could use a clean-up
    CALL read_field_from_xy_file_2D( filename_refgeo, field_name_options_Hb, region_name, refgeo%grid, refgeo%Hb_raw, refgeo%wHb_raw)
    CALL deallocate_grid(refgeo%grid) ! Deallocate the grid, as it is automatically being allocated. Could use a clean-up
    CALL read_field_from_xy_file_2D( filename_refgeo, field_name_options_Hs, region_name, refgeo%grid, refgeo%Hs_raw, refgeo%wHs_raw)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_reference_geometry_file

  SUBROUTINE read_restart_file_geometry( filename_refgeo, grid, refgeo, region_name, time_to_read)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),               INTENT(IN)    :: grid
    TYPE(type_reference_geometry), INTENT(INOUT) :: refgeo
    CHARACTER(LEN=256),            INTENT(IN)    :: filename_refgeo
    CHARACTER(LEN=3),              INTENT(IN)    :: region_name
    REAL(dp),                      INTENT(IN)    :: time_to_read

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'read_restart_file_geometry'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Read the data
    CALL read_field_from_file_2D(   filename_refgeo, field_name_options_Hi, grid,  refgeo%Hi,  region_name, time_to_read)
    CALL read_field_from_file_2D(   filename_refgeo, field_name_options_Hb, grid,  refgeo%Hb,  region_name, time_to_read)
    CALL read_field_from_file_2D(   filename_refgeo, field_name_options_Hs, grid,  refgeo%Hs,  region_name, time_to_read)

    ! NOTE: SL and dHB were loaded in IMAU-ICE 2.1, however were deallocated before used. I have therefore
    ! commented them out for now. If these fields are requiered they should be added to the refgeo type.
    ! CALL read_field_from_file_2D(         filename, 'SL', grid,  refgeo%SL,  region_name, time_to_read)
    ! CALL read_field_from_file_2D(         filename, 'dHb', grid, refgeo%dHb, region_name, time_to_read)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_restart_file_geometry

END MODULE reference_fields_module
