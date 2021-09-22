MODULE reference_fields_module

  ! Contains the routines for setting up the "PD" and "init" reference data fields.

  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE data_types_module,               ONLY: type_grid, type_init_data_fields, type_PD_data_fields
  USE parameters_module,               ONLY: seawater_density, ice_density, sec_per_year, pi
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, inquire_PD_data_file, inquire_init_data_file, &
                                             inquire_restart_file, read_PD_data_file, read_init_data_file, read_restart_file
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             map_square_to_square_cons_2nd_order_2D, map_square_to_square_cons_2nd_order_3D, &
                                             transpose_dp_2D, transpose_dp_3D, remove_Lake_Vostok, surface_elevation

  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE initialise_PD_data_fields( PD, region_name)
    ! Allocate memory for the reference data fields, read them from the specified NetCDF file (latter only done by master process).
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_PD_data_fields),      INTENT(INOUT) :: PD
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    
    IF (C%do_benchmark_experiment) THEN
      CALL initialise_PD_data_fields_schematic_benchmarks( PD)
      RETURN
    END IF
    
    IF      (region_name == 'NAM') THEN
      PD%netcdf%filename   = C%filename_PD_NAM
    ELSE IF (region_name == 'EAS') THEN
      PD%netcdf%filename   = C%filename_PD_EAS
    ELSE IF (region_name == 'GRL') THEN
      PD%netcdf%filename   = C%filename_PD_GRL
    ELSE IF (region_name == 'ANT') THEN
      PD%netcdf%filename   = C%filename_PD_ANT
    END IF
    
    ! Inquire if all the required fields are present in the specified NetCDF file,
    ! and determine the dimensions of the memory to be allocated.
    CALL allocate_shared_int_0D( PD%nx, PD%wnx)
    CALL allocate_shared_int_0D( PD%ny, PD%wny)
    IF (par%master) CALL inquire_PD_data_file(PD)
    CALL sync
    
    ! Allocate memory
    CALL allocate_shared_dp_1D( PD%nx,        PD%x,      PD%wx     )
    CALL allocate_shared_dp_1D(        PD%ny, PD%y,      PD%wy     )
    CALL allocate_shared_dp_2D( PD%nx, PD%ny, PD%Hi_raw, PD%wHi_raw)
    CALL allocate_shared_dp_2D( PD%nx, PD%ny, PD%Hb_raw, PD%wHb_raw)
    CALL allocate_shared_dp_2D( PD%nx, PD%ny, PD%Hs_raw, PD%wHs_raw)
  
    ! Read data from input file
    IF (par%master) WRITE(0,*) '  Reading PD      data from file ', TRIM(PD%netcdf%filename), '...'
    IF (par%master) CALL read_PD_data_file( PD)
    
    ! Safety
    CALL check_for_NaN_dp_2D( PD%Hi_raw, 'PD%Hi_raw', 'initialise_PD_data_fields')
    CALL check_for_NaN_dp_2D( PD%Hb_raw, 'PD%Hb_raw', 'initialise_PD_data_fields')
    CALL check_for_NaN_dp_2D( PD%Hs_raw, 'PD%Hs_raw', 'initialise_PD_data_fields')
    
    ! Since we want data represented as [j,i] internally, transpose the data we just read.
    CALL transpose_dp_2D( PD%Hi_raw, PD%wHi_raw)
    CALL transpose_dp_2D( PD%Hb_raw, PD%wHb_raw)
    CALL transpose_dp_2D( PD%Hs_raw, PD%wHs_raw)
    
    ! Remove Lake Vostok from Antarctica (because it's annoying)
    IF (region_name == 'ANT') CALL remove_Lake_Vostok( PD%x, PD%y, PD%Hi_raw, PD%Hb_raw, PD%Hs_raw)
    
  END SUBROUTINE initialise_PD_data_fields
  SUBROUTINE initialise_PD_data_fields_schematic_benchmarks( PD)
    ! Allocate memory for the reference data fields, initialise them
    ! for schematic benchmark experiments (Halfar dome, EISMINT, MISMIP, etc.)
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_PD_data_fields),      INTENT(INOUT) :: PD
    
    ! Local variables:
    INTEGER                                       :: i,j
    
    REAL(dp), PARAMETER                           :: EISMINT_xmin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_xmax =  750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymax =  750000._dp
    
    CALL allocate_shared_int_0D( PD%nx, PD%wnx)
    CALL allocate_shared_int_0D( PD%ny, PD%wny)
    
    IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
        C%choice_benchmark_experiment == 'Halfar' .OR. &
        C%choice_benchmark_experiment == 'Bueler') THEN
      PD%nx = 51
      PD%ny = 51
    ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
      PD%nx =  41
      PD%ny = 241
    ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
      PD%nx = 751
      PD%ny = 751
    ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
      PD%nx = CEILING( C%ISMIP_HOM_L / C%dx_ANT)
      PD%ny = CEILING( C%ISMIP_HOM_L / C%dx_ANT)
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_data_fields_schematic_benchmarks!'
      STOP
    END IF
    
    ! Allocate memory - PD
    CALL allocate_shared_dp_1D(        PD%nx, PD%x,      PD%wx     )
    CALL allocate_shared_dp_1D( PD%ny,        PD%y,      PD%wy     )
    CALL allocate_shared_dp_2D( PD%ny, PD%nx, PD%Hi_raw, PD%wHi_raw)
    CALL allocate_shared_dp_2D( PD%ny, PD%nx, PD%Hb_raw, PD%wHb_raw)
    CALL allocate_shared_dp_2D( PD%ny, PD%nx, PD%Hs_raw, PD%wHs_raw)
  
    IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
        C%choice_benchmark_experiment == 'Halfar' .OR. &
        C%choice_benchmark_experiment == 'Bueler') THEN
  
      ! Simple square grid with no data
      IF (par%master) THEN
        DO i = 1, PD%nx
          PD%x(i) = EISMINT_xmin + (EISMINT_xmax - EISMINT_xmin) * (i-1) / (PD%nx-1)
        END DO
        DO j = 1, PD%ny
          PD%y(j) = EISMINT_ymin + (EISMINT_ymax - EISMINT_ymin) * (j-1) / (PD%ny-1)
        END DO
        
        ! Note: data set to zero for now, filled in after mapping to model grid to circumvent mapping errors
        PD%Hi_raw = 0._dp
        PD%Hb_raw = 0._dp
        PD%Hs_raw = 0._dp
        
      END IF ! IF (par%master) THEN
      CALL sync
      
    ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
  
      ! Simple square grid with no data
      IF (par%master) THEN
        DO i = 1, PD%nx
          PD%x(i) =  -20000._dp +  40000._dp * (i-1) / (PD%nx-1)
        END DO
        DO j = 1, PD%ny
          PD%y(j) = -120000._dp + 240000._dp * (j-1) / (PD%ny-1)
        END DO
        
        ! Note: data set to zero for now, filled in after mapping to model grid to circumvent mapping errors
        PD%Hi_raw = 0._dp
        PD%Hb_raw = 0._dp
        PD%Hs_raw = 0._dp
        
      END IF ! IF (par%master) THEN
      CALL sync
    
    ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
    
      IF (par%master) THEN
        DO i = 1, PD%nx
          PD%x(i) = -1500000._dp + 3000000._dp * (i-1) / (PD%nx-1)
        END DO
        DO j = 1, PD%ny
          PD%y(j) = -1500000._dp + 3000000._dp * (j-1) / (PD%ny-1)
        END DO
        
        ! Note: data set to zero for now, filled in after mapping to model grid to circumvent mapping errors
        PD%Hi_raw = 0._dp
        PD%Hb_raw = 0._dp
        PD%Hs_raw = 0._dp
        
      END IF ! IF (par%master) THEN
      CALL sync
    
    ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_E' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
    
      IF (par%master) THEN
        DO i = 1, PD%nx
          PD%x(i) = -0.5_dp * C%ISMIP_HOM_L + C%ISMIP_HOM_L * (i-1) / (PD%nx-1)
        END DO
        DO j = 1, PD%ny
          PD%y(j) = -0.5_dp * C%ISMIP_HOM_L + C%ISMIP_HOM_L * (j-1) / (PD%ny-1)
        END DO
        
        ! Note: data set to zero for now, filled in after mapping to model grid to circumvent mapping errors
        PD%Hi_raw = 0._dp
        PD%Hb_raw = 0._dp
        PD%Hs_raw = 0._dp
        
      END IF ! IF (par%master) THEN
      CALL sync
      
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_data_fields_schematic_benchmarks!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE initialise_PD_data_fields_schematic_benchmarks
  SUBROUTINE initialise_init_data_fields( init, region_name)
    ! Allocate memory for the reference data fields, read them from the specified NetCDF file (latter only done by master process).
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    
    IF (C%do_benchmark_experiment) THEN
      CALL initialise_init_data_fields_schematic_benchmarks( init)
      RETURN
    END IF
    
    ! Select the file to read from
    IF      (region_name == 'NAM') THEN
      init%netcdf%filename   = C%filename_init_NAM
    ELSE IF (region_name == 'EAS') THEN
      init%netcdf%filename   = C%filename_init_EAS
    ELSE IF (region_name == 'GRL') THEN
      init%netcdf%filename   = C%filename_init_GRL
    ELSE IF (region_name == 'ANT') THEN
      init%netcdf%filename   = C%filename_init_ANT
    END IF
    
    CALL allocate_shared_int_0D( init%nx, init%wnx)
    CALL allocate_shared_int_0D( init%ny, init%wny)
    CALL allocate_shared_int_0D( init%nz, init%wnz)
    CALL allocate_shared_int_0D( init%nt, init%wnt)
  
    ! Inquire if all the required fields are present in the specified NetCDF file,
    ! and determine the dimensions of the memory to be allocated.
    IF (.NOT. C%is_restart) THEN
      ! Use an externally generated initial file, containing only Hi,Hb,Hs without zeta or time dimensions
      IF (par%master) CALL inquire_init_data_file( init)
    ELSE
      ! Use the restart file of an earlier IMAU-ICE run, which also contains Ti, FirnDepth and MeltPreviousYear (and zeta, time and month dimensions)
      IF (par%master) CALL inquire_restart_file(   init)
    END IF
    CALL sync
    
    ! Read data from input file
    IF (.NOT. C%is_restart) THEN
      ! Use an externally generated initial file, containing only Hi,Hb,Hs without zeta or time dimensions
  
      ! Allocate memory - init
      CALL allocate_shared_dp_1D( init%nx,          init%x,         init%wx     )
      CALL allocate_shared_dp_1D(          init%ny, init%y,         init%wy     )
      CALL allocate_shared_dp_2D( init%nx, init%ny, init%Hi_raw,    init%wHi_raw)
      CALL allocate_shared_dp_2D( init%nx, init%ny, init%Hb_raw,    init%wHb_raw)
      CALL allocate_shared_dp_2D( init%nx, init%ny, init%Hs_raw,    init%wHs_raw)
    
      ! Read data from input file
      IF (par%master) WRITE(0,*) '  Reading init    data from file ', TRIM(init%netcdf%filename), '...'
      IF (par%master) CALL read_init_data_file( init)
    
      ! Safety
      CALL check_for_NaN_dp_2D( init%Hi_raw, 'init%Hi_raw', 'initialise_init_data_fields')
      CALL check_for_NaN_dp_2D( init%Hb_raw, 'init%Hb_raw', 'initialise_init_data_fields')
      CALL check_for_NaN_dp_2D( init%Hs_raw, 'init%Hs_raw', 'initialise_init_data_fields')
      
      ! Since we want data represented as [j,i] internally, transpose the data we just read.
      CALL transpose_dp_2D( init%Hi_raw, init%wHi_raw)
      CALL transpose_dp_2D( init%Hb_raw, init%wHb_raw)
      CALL transpose_dp_2D( init%Hs_raw, init%wHs_raw)
    
      ! Remove Lake Vostok from Antarctica (because it's annoying)
      IF (region_name == 'ANT') CALL remove_Lake_Vostok( init%x, init%y, init%Hi_raw, init%Hb_raw, init%Hs_raw)
    
    ELSE ! IF (.NOT. C%is_restart) THEN
      ! Use the restart file of an earlier IMAU-ICE run, which also contains Ti, FirnDepth and MeltPreviousYear (and zeta, time and month dimensions)
  
      ! Allocate memory - init
      CALL allocate_shared_dp_1D( init%nx,                   init%x,                    init%wx                   )
      CALL allocate_shared_dp_1D(          init%ny,          init%y,                    init%wy                   )
      CALL allocate_shared_dp_1D(                   init%nz, init%zeta,                 init%wzeta                )
      CALL allocate_shared_dp_1D( init%nt,                   init%time,                 init%wtime                )
      
      CALL allocate_shared_dp_2D( init%nx, init%ny,          init%Hi_raw,               init%wHi_raw              )
      CALL allocate_shared_dp_2D( init%nx, init%ny,          init%Hb_raw,               init%wHb_raw              )
      CALL allocate_shared_dp_2D( init%nx, init%ny,          init%SL_raw,               init%wSL_raw              )
      CALL allocate_shared_dp_2D( init%nx, init%ny,          init%dHb_raw,              init%wdHb_raw             )
      CALL allocate_shared_dp_3D( init%nx, init%ny, init%nz, init%Ti_raw,               init%wTi_raw              )
      CALL allocate_shared_dp_3D( init%nx, init%ny, 12,      init%FirnDepth_raw,        init%wFirnDepth_raw       )
      CALL allocate_shared_dp_2D( init%nx, init%ny,          init%MeltPreviousYear_raw, init%wMeltPreviousYear_raw)
    
      ! Read data from input file
      IF (par%master) WRITE(0,*) '  Reading restart data from file ', TRIM(init%netcdf%filename), '...'
      IF (par%master) CALL read_restart_file( init)
      
      ! Since we want data represented as [j,i] internally, transpose the data we just read.
      CALL transpose_dp_2D( init%Hi_raw,               init%wHi_raw              )
      CALL transpose_dp_2D( init%Hb_raw,               init%wHb_raw              )
      CALL transpose_dp_2D( init%SL_raw,               init%wSL_raw              )
      CALL transpose_dp_2D( init%dHb_raw,              init%wdHb_raw             )
      CALL transpose_dp_3D( init%Ti_raw,               init%wTi_raw              )
      CALL transpose_dp_3D( init%FirnDepth_raw,        init%wFirnDepth_raw       )
      CALL transpose_dp_2D( init%MeltPreviousYear_raw, init%wMeltPreviousYear_raw)
    
    END IF ! IF (.NOT. C%is_restart) THEN
    
  END SUBROUTINE initialise_init_data_fields
  SUBROUTINE initialise_init_data_fields_schematic_benchmarks( init)
    ! Allocate memory for the reference data fields, initialise them
    ! for schematic benchmark experiments (Halfar dome, EISMINT, MISMIP, etc.)
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    
    ! Local variables:
    INTEGER                                       :: i,j
    
    REAL(dp), PARAMETER                           :: EISMINT_xmin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_xmax =  750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymax =  750000._dp
    
    ! For the benchmark experiments, use dummy input data.
    ! For realistic experiments, read the provided input file.
    
    CALL allocate_shared_int_0D( init%nx, init%wnx)
    CALL allocate_shared_int_0D( init%ny, init%wny)
    
    IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
        C%choice_benchmark_experiment == 'Halfar' .OR. &
        C%choice_benchmark_experiment == 'Bueler') THEN
      init%nx = 51
      init%ny = 51
    ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
      init%nx =  41
      init%ny = 241
    ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
      init%nx = 51
      init%ny = 51
    ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
      init%nx = CEILING( C%ISMIP_HOM_L / C%dx_ANT)
      init%ny = CEILING( C%ISMIP_HOM_L / C%dx_ANT)
    ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_E') THEN
      init%nx = 51
      init%ny = 251
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_init_data_fields_schematic_benchmarks!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate memory - init
    CALL allocate_shared_dp_1D(          init%nx, init%x,      init%wx     )
    CALL allocate_shared_dp_1D( init%ny,          init%y,      init%wy     )
    CALL allocate_shared_dp_2D( init%ny, init%nx, init%Hi_raw, init%wHi_raw)
    CALL allocate_shared_dp_2D( init%ny, init%nx, init%Hb_raw, init%wHb_raw)
    CALL allocate_shared_dp_2D( init%ny, init%nx, init%Hs_raw, init%wHs_raw)
  
    IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
        C%choice_benchmark_experiment == 'Halfar' .OR. &
        C%choice_benchmark_experiment == 'Bueler') THEN
  
      ! Simple square grid with no data
      IF (par%master) THEN
        DO i = 1, init%nx
          init%x(i) = EISMINT_xmin + (EISMINT_xmax - EISMINT_xmin) * (i-1) / (init%nx-1)
        END DO
        DO j = 1, init%ny
          init%y(j) = EISMINT_ymin + (EISMINT_ymax - EISMINT_ymin) * (j-1) / (init%ny-1)
        END DO
        
        ! Note: data set to zero for now, filled in after mapping to model grid to circumvent mapping errors
        init%Hi_raw = 0._dp
        init%Hb_raw = 0._dp
        init%Hs_raw = 0._dp
        
      END IF ! IF (par%master) THEN
      CALL sync
      
    ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
  
      ! Simple square grid with no data
      IF (par%master) THEN
        DO i = 1, init%nx
          init%x(i) =  -20000._dp +  40000._dp * (i-1) / (init%nx-1)
        END DO
        DO j = 1, init%ny
          init%y(j) = -120000._dp + 240000._dp * (j-1) / (init%ny-1)
        END DO
        
        ! Note: data set to zero for now, filled in after mapping to model grid to circumvent mapping errors
        init%Hi_raw = 0._dp
        init%Hb_raw = 0._dp
        init%Hs_raw = 0._dp
        
      END IF ! IF (par%master) THEN
      CALL sync
    
    ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
    
      IF (par%master) THEN
        DO i = 1, init%nx
        DO j = 1, init%ny
          init%x(i) = -1500000._dp + 3000000._dp * (i-1) / (init%nx-1)
          init%y(j) = -1500000._dp + 3000000._dp * (j-1) / (init%ny-1)
        END DO
        END DO
        
        ! Note: data set to zero for now, filled in after mapping to model grid to circumvent mapping errors
        init%Hi_raw = 0._dp
        init%Hb_raw = 0._dp
        init%Hs_raw = 0._dp
        
      END IF ! IF (par%master) THEN
      CALL sync
    
    ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_A' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_B' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_D' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
    
      IF (par%master) THEN
        DO i = 1, init%nx
          init%x(i) = -0.5_dp * C%ISMIP_HOM_L + C%ISMIP_HOM_L * (i-1) / (init%nx-1)
        END DO
        DO j = 1, init%ny
          init%y(j) = -0.5_dp * C%ISMIP_HOM_L + C%ISMIP_HOM_L * (j-1) / (init%ny-1)
        END DO
        
        ! Note: data set to zero for now, filled in after mapping to model grid to circumvent mapping errors
        init%Hi_raw = 0._dp
        init%Hb_raw = 0._dp
        init%Hs_raw = 0._dp
        
      END IF ! IF (par%master) THEN
      CALL sync
    
    ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_E') THEN
    
      IF (par%master) THEN
        DO i = 1, init%nx
          init%x(i) = 5000._dp * (i-1) / (init%nx-1)
        END DO
        DO j = 1, init%ny
          init%y(j) = 25000._dp * (j-1) / (init%ny-1)
        END DO
        
        ! Note: data set to zero for now, filled in after mapping to model grid to circumvent mapping errors
        init%Hi_raw = 0._dp
        init%Hb_raw = 0._dp
        init%Hs_raw = 0._dp
        
      END IF ! IF (par%master) THEN
      CALL sync
      
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_init_data_fields_schematic_benchmarks!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE initialise_init_data_fields_schematic_benchmarks
  
  SUBROUTINE map_PD_data_to_model_grid(   grid, PD)
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_PD_data_fields),      INTENT(INOUT) :: PD
    
    IF (par%master) WRITE(0,*) '  Mapping PD      data to model grid...'
    
    ! Map the PD data from the provided grid to the model grid
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, PD%Hi, PD%wHi)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, PD%Hb, PD%wHb)
    
    CALL map_square_to_square_cons_2nd_order_2D( PD%nx, PD%ny, PD%x, PD%y, grid%nx, grid%ny, grid%x, grid%y, PD%Hi_raw, PD%Hi)
    CALL map_square_to_square_cons_2nd_order_2D( PD%nx, PD%ny, PD%x, PD%y, grid%nx, grid%ny, grid%x, grid%y, PD%Hb_raw, PD%Hb)
    
  END SUBROUTINE map_PD_data_to_model_grid
  SUBROUTINE map_init_data_to_model_grid( grid, init)
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    
    ! Local variables:
    INTEGER                                       :: i,j
    
    IF (C%do_benchmark_experiment) THEN
      ! No need to map, just initialise the schematic world directly on the model grid
      CALL initialise_initial_model_state_schematic_benchmarks( grid, init)
      RETURN
    END IF
    
    IF (par%master) WRITE(0,*) '  Mapping init    data to model grid...'
    
    ! Map the init data from the provided grid to the model grid
    IF (.NOT. C%is_restart) THEN
      ! External initial file (Bedmachine) has only geometry, no ice temperature or SMB history
    
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, init%Hi, init%wHi)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, init%Hb, init%wHb)
      
      CALL map_square_to_square_cons_2nd_order_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%Hi_raw, init%Hi)
      CALL map_square_to_square_cons_2nd_order_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%Hb_raw, init%Hb)
    
    ELSE
      ! IMAU-ICE restart files have ice geometry, as well as temperature and SMB history
    
      CALL allocate_shared_dp_2D(       grid%ny, grid%nx, init%Hi,               init%wHi              )
      CALL allocate_shared_dp_2D(       grid%ny, grid%nx, init%Hb,               init%wHb              )
      CALL allocate_shared_dp_2D(       grid%ny, grid%nx, init%SL,               init%wSL              )
      CALL allocate_shared_dp_2D(       grid%ny, grid%nx, init%dHb,              init%wdHb             )
      CALL allocate_shared_dp_3D( C%nz, grid%ny, grid%nx, init%Ti,               init%wTi              )
      CALL allocate_shared_dp_3D( 12,   grid%ny, grid%nx, init%FirnDepth,        init%wFirnDepth       )
      CALL allocate_shared_dp_2D(       grid%ny, grid%nx, init%MeltPreviousYear, init%wMeltPreviousYear)
    
      IF (ABS( 1._dp - grid%dx / ABS(init%x(2) - init%x(1))) < 1E-4_dp) THEN
        ! If the model grid is identical to the restart file grid, no remapping is needed
        
        init%Hi(                 :,grid%i1:grid%i2) = init%Hi_raw(                 :,grid%i1:grid%i2)
        init%Hb(                 :,grid%i1:grid%i2) = init%Hb_raw(                 :,grid%i1:grid%i2)
        init%SL(                 :,grid%i1:grid%i2) = init%SL_raw(                 :,grid%i1:grid%i2)
        init%dHb(                :,grid%i1:grid%i2) = init%dHb_raw(                :,grid%i1:grid%i2)
        init%Ti(               :,:,grid%i1:grid%i2) = init%Ti_raw(               :,:,grid%i1:grid%i2)
        init%FirnDepth(        :,:,grid%i1:grid%i2) = init%FirnDepth_raw(        :,:,grid%i1:grid%i2)
        init%MeltPreviousYear(   :,grid%i1:grid%i2) = init%MeltPreviousYear_raw(   :,grid%i1:grid%i2)
        CALL sync
        
      ELSE
        ! If the model grid is different from the restart file grid, remap the data
      
        CALL map_square_to_square_cons_2nd_order_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%Hi_raw,               init%Hi                    )
        CALL map_square_to_square_cons_2nd_order_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%Hb_raw,               init%Hb                    )
        CALL map_square_to_square_cons_2nd_order_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%SL_raw,               init%SL                    )
        CALL map_square_to_square_cons_2nd_order_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%dHb_raw,              init%dHb                   )
        CALL map_square_to_square_cons_2nd_order_3D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%Ti_raw,               init%Ti              , C%nz)
        CALL map_square_to_square_cons_2nd_order_3D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%FirnDepth_raw,        init%FirnDepth       , 12  )
        CALL map_square_to_square_cons_2nd_order_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%MeltPreviousYear_raw, init%MeltPreviousYear      )
        
        ! Some small manual corrections
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          
          IF (init%Hi( j,i) < 0._dp) THEN
            ! Don't allow negative ice thickness
            init%Hi(                 j,i) = 0._dp
            init%Ti(               :,j,i) = init%Ti( 1,j,i)
            init%FirnDepth(        :,j,i) = 0._dp
            init%MeltPreviousYear(   j,i) = 0._dp
          END IF
          
        END DO
        END DO
        CALL sync
        
      END IF ! IF (ABS( 1._dp - grid%dx / ABS(init%x(2) - init%x(1))) < 1E-4_dp) THEN
      
    END IF ! IF (.NOT. C%is_restart) THEN
    
  END SUBROUTINE map_init_data_to_model_grid
  SUBROUTINE initialise_initial_model_state_schematic_benchmarks( grid, init)
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    
    ! Local variables:
    INTEGER                                       :: i,j
    REAL(dp)                                      :: x,Hs,Hb
    INTEGER                                       :: ios,slides
    REAL(dp)                                      :: H0,a0,sigma
    
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, init%Hi, init%wHi)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, init%Hb, init%wHb)
    
    IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
        C%choice_benchmark_experiment == 'EISMINT_6') THEN
      ! These all start ice-free
      
    ELSEIF (C%choice_benchmark_experiment == 'Halfar') THEN
    
      ! Start with the Halfar solution for the given parameters
      init%Hb(:,grid%i1:grid%i2) = 0._dp
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        init%Hb( j,i) = 0._dp
        init%Hi( j,i) = Halfar_solution( grid%x(i), grid%y(j), C%start_time_of_run)
      END DO
      END DO
    
    ELSEIF (C%choice_benchmark_experiment == 'Bueler') THEN
    
      ! Start with the Bueler solution for the given parameters
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        init%Hb( j,i) = 0._dp
        init%Hi( j,i) = Bueler_solution( grid%x(i), grid%y(j), C%start_time_of_run)
      END DO
      END DO
    
    ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
    
      ! Start with the set-up of the Schoof2006 SSA solution (according to Bueler&Brown 2009)
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        init%Hb( j,i) = -0.001_dp * grid%x( i)
        init%Hi( j,i) = 2000._dp
      END DO
      END DO
    
    ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        init%Hi( j,i) = 100._dp
        init%Hb( j,i) = 720._dp - 778.5_dp * SQRT( grid%x(i)**2 + grid%y(j)**2)/ 750000._dp
      END DO
      END DO
      
    ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_A') THEN
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        Hs = 2000._dp - grid%x( i) * TAN( 0.5_dp * pi / 180._dp)
        init%Hb( j,i) = Hs - 1000._dp + 500._dp * SIN( grid%x(i) * 2._dp * pi / C%ISMIP_HOM_L) * SIN( grid%y(j) * 2._dp * pi / C%ISMIP_HOM_L)
        init%Hi( j,i) = Hs - init%Hb( j,i)
      END DO
      END DO
      
    ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_B') THEN
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        Hs = 2000._dp - grid%x( i) * TAN( 0.5_dp * pi / 180._dp)
        init%Hb( j,i) = Hs - 1000._dp + 500._dp * SIN( grid%x(i) * 2._dp * pi / C%ISMIP_HOM_L)
        init%Hi( j,i) = Hs - init%Hb( j,i)
      END DO
      END DO
      
    ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_C' .OR. &
            C%choice_benchmark_experiment == 'ISMIP_HOM_D') THEN
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        Hs = 2000._dp - grid%x( i) * TAN( 0.1_dp * pi / 180._dp)
        init%Hb( j,i) = Hs - 1000._dp
        init%Hi( j,i) = Hs - init%Hb( j,i)
      END DO
      END DO
      
    ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_E') THEN
      ! Read data from external file
      
      IF (par%master) THEN
        
        OPEN(   UNIT = 1337, FILE=C%ISMIP_HOM_E_Arolla_filename, ACTION='READ')
        DO i = 1, 51
          READ( UNIT = 1337, FMT=*, IOSTAT=ios) x, Hb, Hs, slides
          DO j = 1, init%ny
            init%Hb( j,i) = Hb
            init%Hi( j,i) = Hs - Hb
          END DO
          IF (ios /= 0) THEN
            WRITE(0,*) ' initialise_initial_model_state_schematic_benchmarks - ERROR: length of text file "', TRIM(C%ISMIP_HOM_E_Arolla_filename), '" should be 51 lines!'
            CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
          END IF
        END DO
        CLOSE( UNIT  = 1337)
        
      END IF
      
    ELSEIF (C%choice_benchmark_experiment == 'ISMIP_HOM_F') THEN
      
      H0    = 1000._dp
      a0    = 100._dp
      sigma = 10000._dp
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        Hs = 5000._dp - grid%x( i) * TAN( 3._dp * pi / 180._dp)
        init%Hb( j,i) = Hs - H0 + a0 * EXP( -(init%x( i)**2 + init%y( j)**2) / sigma**2)
        init%Hi( j,i) = Hs - init%Hb( j,i)
      END DO
      END DO
      
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_initial_model_state_schematic_benchmarks!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    CALL sync
    
  END SUBROUTINE initialise_initial_model_state_schematic_benchmarks
  
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

END MODULE reference_fields_module
