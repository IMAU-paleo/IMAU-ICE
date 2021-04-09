MODULE reference_fields_module

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parallel_module,                 ONLY: par, sync, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_init_data_fields, type_PD_data_fields
  USE parameters_module,               ONLY: seawater_density, ice_density, sec_per_year
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, inquire_PD_data_file, inquire_init_data_file, &
                                             inquire_restart_file, read_PD_data_file, read_init_data_file, read_restart_file

  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE initialise_PD_data_fields(   PD,   region_name)
    ! Allocate memory for the reference data fields, read them from the specified NetCDF file (latter only done by master process).
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_PD_data_fields),      INTENT(INOUT) :: PD
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER           ::  Hi_raw_temp,  Hb_raw_temp,  Hs_raw_temp
    INTEGER                                       :: wHi_raw_temp, wHb_raw_temp, wHs_raw_temp
    
    REAL(dp), PARAMETER                           :: EISMINT_xmin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_xmax =  750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymax =  750000._dp
    
    IF      (region_name == 'NAM') THEN
      PD%netcdf%filename   = C%filename_PD_NAM
    ELSE IF (region_name == 'EAS') THEN
      PD%netcdf%filename   = C%filename_PD_EAS
    ELSE IF (region_name == 'GRL') THEN
      PD%netcdf%filename   = C%filename_PD_GRL
    ELSE IF (region_name == 'ANT') THEN
        PD%netcdf%filename   = C%filename_PD_ANT
    END IF
    
    ! For the benchmark experiments, use dummy input data.
    ! For realistic experiments, read the provided input file.
    
    CALL allocate_shared_int_0D( PD%nx, PD%wnx)
    CALL allocate_shared_int_0D( PD%ny, PD%wny)
    
    IF (C%do_benchmark_experiment) THEN
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
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_data_fields!'
        STOP
      END IF
      
    ELSE ! IF (C%do_benchmark_experiment) THEN
    
      ! Read data from input file
      IF (par%master) CALL inquire_PD_data_file(PD)
      
    END IF ! IF (C%do_benchmark_experiment) THEN
    CALL sync
    
    ! Read PD data
    IF (C%do_benchmark_experiment) THEN
    
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
          
          PD%Hi_raw = 0._dp
          PD%Hb_raw = 0._dp
          PD%Hs_raw = 0._dp
        END IF ! IF (par%master) THEN
        CALL sync
      
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
      
        IF (par%master) THEN
          DO i = 1, PD%nx
          DO j = 1, PD%ny
            PD%x(i) = -1500000._dp + 3000000._dp * (i-1) / (PD%nx-1)
            PD%y(j) = -1500000._dp + 3000000._dp * (j-1) / (PD%ny-1)
            PD%Hb_raw( i,j) = 720._dp - 778.5_dp * SQRT( PD%x(i)**2 + PD%y(j)**2)/ 750000._dp
          END DO
          END DO
          PD%Hi_raw = 0._dp
          PD%Hs_raw = MAX(0._dp, PD%Hb_raw)
        END IF ! IF (par%master) THEN
        CALL sync
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_data_fields!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSE ! IF (C%do_benchmark_experiment) THEN
    
      ! Allocate memory - PD
      CALL allocate_shared_dp_1D( PD%nx,        PD%x,      PD%wx     )
      CALL allocate_shared_dp_1D(        PD%ny, PD%y,      PD%wy     )
      CALL allocate_shared_dp_2D( PD%nx, PD%ny, PD%Hi_raw, PD%wHi_raw)
      CALL allocate_shared_dp_2D( PD%nx, PD%ny, PD%Hb_raw, PD%wHb_raw)
      CALL allocate_shared_dp_2D( PD%nx, PD%ny, PD%Hs_raw, PD%wHs_raw)
    
      ! Read data from input file
      IF (par%master) WRITE(0,*) '  Reading PD      data from file ', TRIM(PD%netcdf%filename), '...'
      IF (par%master) CALL read_PD_data_file( PD)
      
      ! Since we want data represented as [j,i] internally, transpose the data we just read.
      CALL allocate_shared_dp_2D( PD%nx, PD%ny, Hi_raw_temp, wHi_raw_temp)
      CALL allocate_shared_dp_2D( PD%nx, PD%ny, Hb_raw_temp, wHb_raw_temp)
      CALL allocate_shared_dp_2D( PD%nx, PD%ny, Hs_raw_temp, wHs_raw_temp)
      
      IF (par%master) THEN
        Hi_raw_temp = PD%Hi_raw
        Hb_raw_temp = PD%Hb_raw
        Hs_raw_temp = PD%Hs_raw
      END IF
      CALL sync
      
      CALL deallocate_shared( PD%wHi_raw)
      CALL deallocate_shared( PD%wHb_raw)
      CALL deallocate_shared( PD%wHs_raw)
      
      CALL allocate_shared_dp_2D( PD%ny, PD%nx, PD%Hi_raw, PD%wHi_raw)
      CALL allocate_shared_dp_2D( PD%ny, PD%nx, PD%Hb_raw, PD%wHb_raw)
      CALL allocate_shared_dp_2D( PD%ny, PD%nx, PD%Hs_raw, PD%wHs_raw)
      
      IF (par%master) THEN
        DO i = 1, PD%nx
        DO j = 1, PD%ny
          PD%Hi_raw( j,i) = Hi_raw_temp( i,j)
          PD%Hb_raw( j,i) = Hb_raw_temp( i,j)
          PD%Hs_raw( j,i) = Hs_raw_temp( i,j)
        END DO
        END DO
      END IF
      CALL sync
      
      CALL deallocate_shared( wHi_raw_temp)
      CALL deallocate_shared( wHb_raw_temp)
      CALL deallocate_shared( wHs_raw_temp)
      
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  END SUBROUTINE initialise_PD_data_fields
  SUBROUTINE initialise_topo_data_fields(   topo,  PD,   region_name)
    ! Allocate memory for the reference data fields, read them from the specified NetCDF file (latter only done by master process).
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_PD_data_fields),      INTENT(INOUT) :: topo
    TYPE(type_PD_data_fields),      INTENT(IN)    :: PD
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
  
    ! Local variables:
    INTEGER                                       :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER           ::  Hi_raw_temp,  Hb_raw_temp,  Hs_raw_temp
    INTEGER                                       :: wHi_raw_temp, wHb_raw_temp, wHs_raw_temp
  
    IF      (region_name == 'NAM') THEN
      topo%netcdf%filename   = C%filename_topo_NAM
    ELSE IF (region_name == 'EAS') THEN
      topo%netcdf%filename   = C%filename_topo_EAS
    ELSE IF (region_name == 'GRL') THEN
      topo%netcdf%filename   = C%filename_topo_GRL
    ELSE IF (region_name == 'ANT') THEN
      topo%netcdf%filename   = C%filename_topo_ANT
    END IF
  
    CALL allocate_shared_int_0D( topo%nx, topo%wnx)
    CALL allocate_shared_int_0D( topo%ny, topo%wny)
    
    IF (C%do_benchmark_experiment) THEN
      ! Just use the same field as PD
      IF (par%master) THEN
        topo%nx = PD%nx
        topo%ny = PD%ny
      END IF
    ELSE
      ! Read data from input file
      IF (par%master) CALL inquire_PD_data_file(topo)
    END IF ! (C%do_benchmark_experiment)
    CALL sync

    ! Read paleo data

    ! Allocate memory - PD
    CALL allocate_shared_dp_1D( topo%nx,               topo%x,      topo%wx)
    CALL allocate_shared_dp_1D(               topo%ny, topo%y,      topo%wy)
    CALL allocate_shared_dp_2D( topo%nx, topo%ny, topo%Hi_raw, topo%wHi_raw)
    CALL allocate_shared_dp_2D( topo%nx, topo%ny, topo%Hb_raw, topo%wHb_raw)
    CALL allocate_shared_dp_2D( topo%nx, topo%ny, topo%Hs_raw, topo%wHs_raw)

    IF (C%do_benchmark_experiment) THEN
  
      IF (par%master) THEN
        ! Just use the same field as PD
        topo%x      = PD%x
        topo%y      = PD%y
        topo%Hi_raw = PD%Hi_raw
        topo%Hs_raw = PD%Hs_raw
        topo%Hb_raw = PD%Hb_raw
      END IF ! IF (par%master) THEN
      CALL sync

    ELSE
      ! Read data from input file
      IF (par%master) WRITE(0,*) '  Reading topo    data from file ', TRIM(topo%netcdf%filename), '...'
      IF (par%master) CALL read_PD_data_file( topo)
    
      ! Since we want data represented as [j,i] internally, transpose the data we just read.
      CALL allocate_shared_dp_2D( topo%nx, topo%ny, Hi_raw_temp, wHi_raw_temp)
      CALL allocate_shared_dp_2D( topo%nx, topo%ny, Hb_raw_temp, wHb_raw_temp)
      CALL allocate_shared_dp_2D( topo%nx, topo%ny, Hs_raw_temp, wHs_raw_temp)
    
      IF (par%master) THEN
        Hi_raw_temp = topo%Hi_raw
        Hb_raw_temp = topo%Hb_raw
        Hs_raw_temp = topo%Hs_raw
      END IF
      CALL sync
    
      CALL deallocate_shared( topo%wHi_raw)
      CALL deallocate_shared( topo%wHb_raw)
      CALL deallocate_shared( topo%wHs_raw)
    
      CALL allocate_shared_dp_2D( topo%ny, topo%nx, topo%Hi_raw, topo%wHi_raw)
      CALL allocate_shared_dp_2D( topo%ny, topo%nx, topo%Hb_raw, topo%wHb_raw)
      CALL allocate_shared_dp_2D( topo%ny, topo%nx, topo%Hs_raw, topo%wHs_raw)
    
      IF (par%master) THEN
        DO i = 1, topo%nx
        DO j = 1, topo%ny
          topo%Hi_raw( j,i) = Hi_raw_temp( i,j)
          topo%Hb_raw( j,i) = Hb_raw_temp( i,j)
          topo%Hs_raw( j,i) = Hs_raw_temp( i,j)
        END DO
        END DO
      END IF
      CALL sync
    
      CALL deallocate_shared( wHi_raw_temp)
      CALL deallocate_shared( wHb_raw_temp)
      CALL deallocate_shared( wHs_raw_temp)

    END IF !(C%do_benchmark_experiment)
  
  END SUBROUTINE initialise_topo_data_fields
  SUBROUTINE initialise_init_data_fields( init, region_name)
    ! Allocate memory for the reference data fields, read them from the specified NetCDF file (latter only done by master process).
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER           :: Hi_raw_temp
    REAL(dp), DIMENSION(:,:  ), POINTER           :: Hb_raw_temp
    REAL(dp), DIMENSION(:,:  ), POINTER           :: Hs_raw_temp
    REAL(dp), DIMENSION(:,:,:), POINTER           :: Ti_raw_temp
    REAL(dp), DIMENSION(:,:  ), POINTER           :: U_SSA_raw_temp
    REAL(dp), DIMENSION(:,:  ), POINTER           :: V_SSA_raw_temp
    REAL(dp), DIMENSION(:,:,:), POINTER           :: FirnDepth_raw_temp
    REAL(dp), DIMENSION(:,:  ), POINTER           :: MeltPreviousYear_raw_temp
    INTEGER                                       :: wHi_raw_temp, wHb_raw_temp, wHs_raw_temp, wTi_raw_temp, wU_SSA_raw_temp, wV_SSA_raw_temp
    INTEGER                                       :: wFirnDepth_raw_temp, wMeltPreviousYear_raw_temp
    
    REAL(dp), PARAMETER                           :: EISMINT_xmin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_xmax =  750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymax =  750000._dp
    
    ! For the benchmark experiments, use dummy input data.
    ! For realistic experiments, read the provided input file.
    
    CALL allocate_shared_int_0D( init%nx, init%wnx)
    CALL allocate_shared_int_0D( init%ny, init%wny)
    CALL allocate_shared_int_0D( init%nz, init%wnz)
    CALL allocate_shared_int_0D( init%nt, init%wnt)
    
    IF (C%do_benchmark_experiment) THEN
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
        init%nx = 751
        init%ny = 751
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_init_data_fields!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSE ! IF (C%do_benchmark_experiment) THEN
    
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
    
      ! Read data from input file
      IF (par%master) THEN
        IF (.NOT. C%is_restart) THEN
          ! Use an externally generated initial file, containing only Hi,Hb,Hs without zeta or time dimensions
          CALL inquire_init_data_file( init)
        ELSE
          ! Use the restart file of an earlier IMAU-ICE run, which also contains Ti, U_SSA, V_SSA, FirnDepth and MeltPreviousYear (and zeta, time and month dimensions)
          CALL inquire_restart_file(   init)
        END IF
      END IF
      CALL sync
      
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Read init data
    IF (C%do_benchmark_experiment) THEN
    
      ! Allocate memory - init
      CALL allocate_shared_dp_1D(        init%nx, init%x,      init%wx     )
      CALL allocate_shared_dp_1D( init%ny,        init%y,      init%wy     )
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
            init%Hb_raw( i,j) = 720._dp - 778.5_dp * SQRT( init%x(i)**2 + init%y(j)**2)/ 750000._dp
          END DO
          END DO
          init%Hi_raw = 0._dp
          init%Hs_raw = MAX(0._dp, init%Hb_raw)
        END IF ! IF (par%master) THEN
        CALL sync
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_init_data_fields!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSE ! IF (C%do_benchmark_experiment) THEN
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
        
        ! Since we want data represented as [j,i] internally, transpose the data we just read.
        CALL allocate_shared_dp_2D( init%nx, init%ny, Hi_raw_temp, wHi_raw_temp)
        CALL allocate_shared_dp_2D( init%nx, init%ny, Hb_raw_temp, wHb_raw_temp)
        CALL allocate_shared_dp_2D( init%nx, init%ny, Hs_raw_temp, wHs_raw_temp)
        
        IF (par%master) THEN
          Hi_raw_temp = init%Hi_raw
          Hb_raw_temp = init%Hb_raw
          Hs_raw_temp = init%Hs_raw
        END IF
        CALL sync
        
        CALL deallocate_shared( init%wHi_raw)
        CALL deallocate_shared( init%wHb_raw)
        CALL deallocate_shared( init%wHs_raw)
        
        CALL allocate_shared_dp_2D( init%ny, init%nx, init%Hi_raw, init%wHi_raw)
        CALL allocate_shared_dp_2D( init%ny, init%nx, init%Hb_raw, init%wHb_raw)
        CALL allocate_shared_dp_2D( init%ny, init%nx, init%Hs_raw, init%wHs_raw)
        
        IF (par%master) THEN
          DO i = 1, init%nx
          DO j = 1, init%ny
            init%Hi_raw( j,i) = Hi_raw_temp( i,j)
            init%Hb_raw( j,i) = Hb_raw_temp( i,j)
            init%Hs_raw( j,i) = Hs_raw_temp( i,j)
          END DO
          END DO
        END IF
        CALL sync
        
        CALL deallocate_shared( wHi_raw_temp)
        CALL deallocate_shared( wHb_raw_temp)
        CALL deallocate_shared( wHs_raw_temp)
      
      ELSE ! IF (.NOT. C%is_restart) THEN
        ! Use the restart file of an earlier IMAU-ICE run, which also contains Ti, U_SSA, V_SSA, FirnDepth and MeltPreviousYear (and zeta, time and month dimensions)
    
        ! Allocate memory - init
        CALL allocate_shared_dp_1D( init%nx,                   init%x,                    init%wx                   )
        CALL allocate_shared_dp_1D(          init%ny,          init%y,                    init%wy                   )
        CALL allocate_shared_dp_1D(                   init%nz, init%zeta,                 init%wzeta                )
        CALL allocate_shared_dp_1D( init%nt,                   init%time,                 init%wtime                )
        CALL allocate_shared_dp_2D( init%nx, init%ny,          init%Hi_raw,               init%wHi_raw              )
        CALL allocate_shared_dp_2D( init%nx, init%ny,          init%Hb_raw,               init%wHb_raw              )
        CALL allocate_shared_dp_2D( init%nx, init%ny,          init%Hs_raw,               init%wHs_raw              )
        CALL allocate_shared_dp_3D( init%nx, init%ny, init%nz, init%Ti_raw,               init%wTi_raw              )
        CALL allocate_shared_dp_2D( init%nx, init%ny,          init%U_SSA_raw,            init%wU_SSA_raw           )
        CALL allocate_shared_dp_2D( init%nx, init%ny,          init%V_SSA_raw,            init%wV_SSA_raw           )
        CALL allocate_shared_dp_3D( init%nx, init%ny, 12,      init%FirnDepth_raw,        init%wFirnDepth_raw       )
        CALL allocate_shared_dp_2D( init%nx, init%ny,          init%MeltPreviousYear_raw, init%wMeltPreviousYear_raw)
      
        ! Read data from input file
        IF (par%master) WRITE(0,*) '  Reading restart data from file ', TRIM(init%netcdf%filename), '...'
        IF (par%master) CALL read_restart_file( init)
        
        ! Since we want data represented as [j,i] internally, transpose the data we just read.
        CALL allocate_shared_dp_2D( init%nx, init%ny,          Hi_raw_temp,               wHi_raw_temp              )
        CALL allocate_shared_dp_2D( init%nx, init%ny,          Hb_raw_temp,               wHb_raw_temp              )
        CALL allocate_shared_dp_2D( init%nx, init%ny,          Hs_raw_temp,               wHs_raw_temp              )
        CALL allocate_shared_dp_3D( init%nx, init%ny, init%nz, Ti_raw_temp,               wTi_raw_temp              )
        CALL allocate_shared_dp_2D( init%nx, init%ny,          U_SSA_raw_temp,            wU_SSA_raw_temp           )
        CALL allocate_shared_dp_2D( init%nx, init%ny,          V_SSA_raw_temp,            wV_SSA_raw_temp           )
        CALL allocate_shared_dp_3D( init%nx, init%ny, 12,      FirnDepth_raw_temp,        wFirnDepth_raw_temp       )
        CALL allocate_shared_dp_2D( init%nx, init%ny,          MeltPreviousYear_raw_temp, wMeltPreviousYear_raw_temp)
        
        IF (par%master) THEN
          Hi_raw_temp               = init%Hi_raw
          Hb_raw_temp               = init%Hb_raw
          Hs_raw_temp               = init%Hs_raw
          Ti_raw_temp               = init%Ti_raw
          U_SSA_raw_temp            = init%U_SSA_raw
          V_SSA_raw_temp            = init%V_SSA_raw
          FirnDepth_raw_temp        = init%FirnDepth_raw
          MeltPreviousYear_raw_temp = init%MeltPreviousYear_raw
        END IF
        CALL sync
        
        CALL deallocate_shared( init%wHi_raw)
        CALL deallocate_shared( init%wHb_raw)
        CALL deallocate_shared( init%wHs_raw)
        CALL deallocate_shared( init%wTi_raw)
        CALL deallocate_shared( init%wU_SSA_raw)
        CALL deallocate_shared( init%wV_SSA_raw)
        CALL deallocate_shared( init%wFirnDepth_raw)
        CALL deallocate_shared( init%wMeltPreviousYear_raw)
        
        CALL allocate_shared_dp_2D(          init%ny, init%nx, init%Hi_raw,               init%wHi_raw              )
        CALL allocate_shared_dp_2D(          init%ny, init%nx, init%Hb_raw,               init%wHb_raw              )
        CALL allocate_shared_dp_2D(          init%ny, init%nx, init%Hs_raw,               init%wHs_raw              )
        CALL allocate_shared_dp_3D( init%nz, init%ny, init%nx, init%Ti_raw,               init%wTi_raw              )
        CALL allocate_shared_dp_2D(          init%ny, init%nx, init%U_SSA_raw,            init%wU_SSA_raw           )
        CALL allocate_shared_dp_2D(          init%ny, init%nx, init%V_SSA_raw,            init%wV_SSA_raw           )
        CALL allocate_shared_dp_3D( 12,      init%ny, init%nx, init%FirnDepth_raw,        init%wFirnDepth_raw       )
        CALL allocate_shared_dp_2D(          init%ny, init%nx, init%MeltPreviousYear_raw, init%wMeltPreviousYear_raw)
        
        IF (par%master) THEN
          DO i = 1, init%nx
          DO j = 1, init%ny
            init%Hi_raw(                 j,i) = Hi_raw_temp(               i,j  )
            init%Hb_raw(                 j,i) = Hb_raw_temp(               i,j  )
            init%Hs_raw(                 j,i) = Hs_raw_temp(               i,j  )
            init%Ti_raw(               :,j,i) = Ti_raw_temp(               i,j,:)
            init%U_SSA_raw(              j,i) = U_SSA_raw_temp(            i,j  )
            init%V_SSA_raw(              j,i) = V_SSA_raw_temp(            i,j  )
            init%FirnDepth_raw(        :,j,i) = FirnDepth_raw_temp(        i,j,:)
            init%MeltPreviousYear_raw(   j,i) = MeltPreviousYear_raw_temp( i,j  )
          END DO
          END DO
        END IF
        CALL sync
        
        CALL deallocate_shared( wHi_raw_temp              )
        CALL deallocate_shared( wHb_raw_temp              )
        CALL deallocate_shared( wHs_raw_temp              )
        CALL deallocate_shared( wTi_raw_temp              )
        CALL deallocate_shared( wU_SSA_raw_temp           )
        CALL deallocate_shared( wV_SSA_raw_temp           )
        CALL deallocate_shared( wFirnDepth_raw_temp       )
        CALL deallocate_shared( wMeltPreviousYear_raw_temp)
      
      END IF ! IF (.NOT. C%is_restart) THEN
      
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  END SUBROUTINE initialise_init_data_fields
  
  SUBROUTINE map_PD_data_to_model_grid(   grid, PD)
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_PD_data_fields),      INTENT(INOUT) :: PD
    
    IF (par%master) WRITE(0,*) '  Mapping PD      data to model grid...'
    
    ! Map the PD data from the provided grid to the model grid
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, PD%Hi, PD%wHi)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, PD%Hb, PD%wHb)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, PD%Hs, PD%wHs)
    
    CALL map_square_to_square_2D( PD%nx, PD%ny, PD%x, PD%y, grid%nx, grid%ny, grid%x, grid%y, PD%Hi_raw, PD%Hi)
    CALL map_square_to_square_2D( PD%nx, PD%ny, PD%x, PD%y, grid%nx, grid%ny, grid%x, grid%y, PD%Hb_raw, PD%Hb)
    CALL map_square_to_square_2D( PD%nx, PD%ny, PD%x, PD%y, grid%nx, grid%ny, grid%x, grid%y, PD%Hs_raw, PD%Hs)
    
    ! Local variables:
  END SUBROUTINE map_PD_data_to_model_grid

  SUBROUTINE map_topo_data_to_model_grid(   grid, topo)
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_PD_data_fields),      INTENT(INOUT) :: topo
  
    IF (par%master) WRITE(0,*) '  Mapping topo    data to model grid...'
  
    ! Map the PD data from the provided grid to the model grid
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, topo%Hi, topo%wHi)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, topo%Hb, topo%wHb)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, topo%Hs, topo%wHs)
  
    CALL map_square_to_square_2D( topo%nx, topo%ny, topo%x, topo%y, grid%nx, grid%ny, grid%x, grid%y, topo%Hi_raw, topo%Hi)
    CALL map_square_to_square_2D( topo%nx, topo%ny, topo%x, topo%y, grid%nx, grid%ny, grid%x, grid%y, topo%Hb_raw, topo%Hb)
    CALL map_square_to_square_2D( topo%nx, topo%ny, topo%x, topo%y, grid%nx, grid%ny, grid%x, grid%y, topo%Hs_raw, topo%Hs)
  
  END SUBROUTINE map_topo_data_to_model_grid

  SUBROUTINE map_init_data_to_model_grid( grid, init)
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    
    ! Local variables:
    INTEGER                                       :: i,j,cerr,ierr
    
    IF (par%master) WRITE(0,*) '  Mapping init    data to model grid...'
    
    ! Map the init data from the provided grid to the model grid
    IF (.NOT. C%is_restart) THEN
    
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, init%Hi, init%wHi)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, init%Hb, init%wHb)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, init%Hs, init%wHs)
      
      CALL map_square_to_square_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%Hi_raw, init%Hi)
      CALL map_square_to_square_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%Hb_raw, init%Hb)
      CALL map_square_to_square_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%Hs_raw, init%Hs)
    
    ELSE
    
      CALL allocate_shared_dp_2D(       grid%ny, grid%nx, init%Hi,               init%wHi              )
      CALL allocate_shared_dp_2D(       grid%ny, grid%nx, init%Hb,               init%wHb              )
      CALL allocate_shared_dp_2D(       grid%ny, grid%nx, init%Hs,               init%wHs              )
      CALL allocate_shared_dp_3D( C%nz, grid%ny, grid%nx, init%Ti,               init%wTi              )
      CALL allocate_shared_dp_2D(       grid%ny, grid%nx, init%U_SSA,            init%wU_SSA           )
      CALL allocate_shared_dp_2D(       grid%ny, grid%nx, init%V_SSA,            init%wV_SSA           )
      CALL allocate_shared_dp_3D( 12,   grid%ny, grid%nx, init%FirnDepth,        init%wFirnDepth       )
      CALL allocate_shared_dp_2D(       grid%ny, grid%nx, init%MeltPreviousYear, init%wMeltPreviousYear)
      
      CALL map_square_to_square_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%Hi_raw,               init%Hi              )
      CALL map_square_to_square_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%Hb_raw,               init%Hb              )
      CALL map_square_to_square_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%Hs_raw,               init%Hs              )
      CALL map_square_to_square_3D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%Ti_raw,               init%Ti              , C%nZ)
      CALL map_square_to_square_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%U_SSA_raw,            init%U_SSA           )
      CALL map_square_to_square_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%V_SSA_raw,            init%V_SSA           )
      CALL map_square_to_square_3D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%FirnDepth_raw,        init%FirnDepth       , 12)
      CALL map_square_to_square_2D( init%nx, init%ny, init%x, init%y, grid%nx, grid%ny, grid%x, grid%y, init%MeltPreviousYear_raw, init%MeltPreviousYear)
      
    END IF
    
    ! Make sure that Hi, Hb and Hs match
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (init%Hi( j,i) > 0._dp .AND. init%Hi( j,i) > -init%Hb( j,i) * seawater_density / ice_density) THEN
        init%Hi( j,i) = init%Hs( j,i) - init%Hb( j,i)
      END IF
    END DO
    END DO
    CALL sync
    
    ! For the Halfar and Bueler experiments, make sure we start with the EXACT solution
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6') THEN
        ! No changes here
      ELSEIF (C%choice_benchmark_experiment == 'Halfar') THEN
      
        ! Start with the Halfar solution for the given parameters
        init%Hb(:,grid%i1:grid%i2) = 0._dp
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          init%Hi(i,j) = Halfar_solution( grid%x(i), grid%y(j), C%start_time_of_run)
        END DO
        END DO
        init%Hs = init%Hi
      
      ELSEIF (C%choice_benchmark_experiment == 'Bueler') THEN
      
        ! Start with the Bueler solution for the given parameters
        init%Hb(:,grid%i1:grid%i2) = 0._dp
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          init%Hi(i,j) = Bueler_solution( grid%x(i), grid%y(j), C%start_time_of_run)
        END DO
        END DO
        init%Hs = init%Hi
      
      ELSEIF (C%choice_benchmark_experiment == 'SSA_icestream') THEN
      
        ! Start with the set-up of the Schoof2006 SSA solution (according to Bueler&Brown 2009)
        init%Hi(:,grid%i1:grid%i2) = 2000._dp
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          init%Hb( j,i) = -0.001_dp * grid%x( i)
        END DO
        END DO
        init%Hs = init%Hb + init%Hi
      
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
      
        init%Hi(:,grid%i1:grid%i2) = 10._dp
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          init%Hb( j,i) = 720._dp - 778.5_dp * SQRT( grid%x(i)**2 + grid%y(j)**2)/ 750000._dp
          init%Hs( j,i) = init%Hi( j,i) + MAX( -ice_density / seawater_density * init%Hi( j,i), init%Hb( j,i))
        END DO
        END DO
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in map_init_data_to_model_grid!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  END SUBROUTINE map_init_data_to_model_grid
  
  SUBROUTINE map_square_to_square_2D( nx_src, ny_src, x_src, y_src, nx_dst, ny_dst, x_dst, y_dst, d_src, d_dst)
    ! Map data from one square grid to another (e.g. PD ice thickness from the square grid in the input file to the model square grid)
    
    IMPLICIT NONE
  
    ! Input and output variables
    INTEGER,                            INTENT(IN)    :: nx_src
    INTEGER,                            INTENT(IN)    :: ny_src
    REAL(dp), DIMENSION(nx_src),        INTENT(IN)    :: x_src
    REAL(dp), DIMENSION(ny_src),        INTENT(IN)    :: y_src
    INTEGER,                            INTENT(IN)    :: nx_dst
    INTEGER,                            INTENT(IN)    :: ny_dst
    REAL(dp), DIMENSION(nx_dst),        INTENT(IN)    :: x_dst
    REAL(dp), DIMENSION(ny_dst),        INTENT(IN)    :: y_dst
    REAL(dp), DIMENSION(ny_src,nx_src), INTENT(IN)    :: d_src
    REAL(dp), DIMENSION(ny_dst,nx_dst), INTENT(OUT)   :: d_dst
    
    ! Local variables
    INTEGER                                       :: i, j, i_src, j_src, i1, i2
    REAL(dp)                                      :: dx_src, dy_src, dx_dst, dy_dst, xcmin, xcmax, ycmin, ycmax
    INTEGER,  DIMENSION(nx_dst,2)                 :: ir_src
    INTEGER,  DIMENSION(ny_dst,2)                 :: jr_src
    REAL(dp)                                      :: xomin, xomax, yomin, yomax
    REAL(dp)                                      :: A_overlap, A_dst, Aint
    
    ! Find grid spacings
    dx_src = x_src(2) - x_src(1)
    dy_src = y_src(2) - y_src(1)
    dx_dst = x_dst(2) - x_dst(1)
    dy_dst = y_dst(2) - y_dst(1)
    A_dst = dx_dst * dy_dst
    
    ! If the grids are equal, the solution is trivial; just copy the data
    IF (dx_src == dx_dst .AND. dy_src == dy_dst .AND. nx_src == nx_dst .AND. ny_src == ny_dst) THEN
      CALL partition_list( nx_dst, par%i, par%n, i1, i2)
      d_dst( :,i1:i2) = d_src( :,i1:i2)
      CALL sync
      RETURN
    END IF
    
    ! Find overlaps between grids
    DO i = 1, nx_dst
      ! Dst cell i overlaps with src cells ir_src( i,1) to ir_src( i,2)
      xcmin = x_dst( i) - dx_dst/2._dp
      xcmax = x_dst( i) + dx_dst/2._dp
      ir_src( i,:) = MAX( 1, MIN( nx_src, [CEILING(-1.5_dp + FLOOR(nx_src/2._dp) + xcmin / dx_src), &
                                           CEILING( 1.5_dp + FLOOR(nx_src/2._dp) + xcmax / dx_src)] ))
    END DO ! DO i = 1, nx_dst
    DO j = 1, ny_dst
      ! Dst cell j overlaps with src cells jr_src( j,1) to ir_src( j,2)
      ycmin = y_dst( j) - dy_dst/2._dp
      ycmax = y_dst( j) + dy_dst/2._dp
      jr_src( j,:) = MAX( 1, MIN( ny_src, [CEILING(-1.5_dp + FLOOR(ny_src/2._dp) + ycmin / dy_src), &
                                           CEILING( 1.5_dp + FLOOR(ny_src/2._dp) + ycmax / dy_src)] ))
    END DO ! DO j = 1, ny_dst
    
    ! Find parallelisation domains
    CALL partition_list( nx_dst, par%i, par%n, i1, i2)
    
    DO i = i1, i2
    DO j = 1, ny_dst
    
      Aint = 0._dp
      
      DO i_src = ir_src( i,1), ir_src( i,2)
      DO j_src = jr_src( j,1), jr_src( j,2)
        
        xomin = MAX( x_dst( i) - dx_dst/2._dp, x_src( i_src) - dx_src/2._dp)
        xomax = MIN( x_dst( i) + dx_dst/2._dp, x_src( i_src) + dx_src/2._dp)
        yomin = MAX( y_dst( j) - dy_dst/2._dp, y_src( j_src) - dy_src/2._dp)
        yomax = MIN( y_dst( j) + dy_dst/2._dp, y_src( j_src) + dy_src/2._dp)
        
        IF (xomax <= xomin .OR. yomax <= yomin) CYCLE
        
        A_overlap = (xomax - xomin) * (yomax - yomin)
        Aint = Aint + A_overlap
        
        d_dst( j,i) = d_dst( j,i) + (A_overlap / A_dst) * d_src( j_src,i_src)
        
      END DO ! DO j_src = jr_src( j,1), jr_src( j,2)
      END DO ! DO i_src = ir_src( i,1), ir_src( i,2)
      
    END DO ! DO j = 1, ny_dst
    END DO ! DO i = 1, nx_dst
    CALL sync
  
  END SUBROUTINE map_square_to_square_2D
  SUBROUTINE map_square_to_square_3D( nx_src, ny_src, x_src, y_src, nx_dst, ny_dst, x_dst, y_dst, d_src, d_dst, nz)
    ! Map data from one square grid to another (e.g. PD ice thickness from the square grid in the input file to the model square grid)
    
    IMPLICIT NONE
  
    ! Input and output variables
    INTEGER,                            INTENT(IN)    :: nx_src
    INTEGER,                            INTENT(IN)    :: ny_src
    REAL(dp), DIMENSION(nx_src),        INTENT(IN)    :: x_src
    REAL(dp), DIMENSION(ny_src),        INTENT(IN)    :: y_src
    INTEGER,                            INTENT(IN)    :: nx_dst
    INTEGER,                            INTENT(IN)    :: ny_dst
    REAL(dp), DIMENSION(nx_dst),        INTENT(IN)    :: x_dst
    REAL(dp), DIMENSION(ny_dst),        INTENT(IN)    :: y_dst
    REAL(dp), DIMENSION(nz,ny_src,nx_src), INTENT(IN)    :: d_src
    REAL(dp), DIMENSION(nz,ny_dst,nx_dst), INTENT(OUT)   :: d_dst
    INTEGER,                            INTENT(IN)    :: nz
    
    ! Local variables
    INTEGER                                       :: i, j, i_src, j_src, i1, i2
    REAL(dp)                                      :: dx_src, dy_src, dx_dst, dy_dst, xcmin, xcmax, ycmin, ycmax
    INTEGER,  DIMENSION(nx_dst,2)                 :: ir_src
    INTEGER,  DIMENSION(ny_dst,2)                 :: jr_src
    REAL(dp)                                      :: xomin, xomax, yomin, yomax
    REAL(dp)                                      :: A_overlap, A_dst, Aint
    
    ! Find grid spacings
    dx_src = x_src(2) - x_src(1)
    dy_src = y_src(2) - y_src(1)
    dx_dst = x_dst(2) - x_dst(1)
    dy_dst = y_dst(2) - y_dst(1)
    A_dst = dx_dst * dy_dst
    
    ! If the grids are equal, the solution is trivial; just copy the data
    IF (dx_src == dx_dst .AND. dy_src == dy_dst .AND. nx_src == nx_dst .AND. ny_src == ny_dst) THEN
      CALL partition_list( nx_dst, par%i, par%n, i1, i2)
      d_dst( :,:,i1:i2) = d_src( :,:,i1:i2)
      CALL sync
      RETURN
    END IF
    
    ! Find overlaps between grids
    DO i = 1, nx_dst
      ! Dst cell i overlaps with src cells ir_src( i,1) to ir_src( i,2)
      xcmin = x_dst( i) - dx_dst/2._dp
      xcmax = x_dst( i) + dx_dst/2._dp
      ir_src( i,:) = MAX( 1, MIN( nx_src, [CEILING(-1.5_dp + FLOOR(nx_src/2._dp) + xcmin / dx_src), &
                                           CEILING( 1.5_dp + FLOOR(nx_src/2._dp) + xcmax / dx_src)] ))
    END DO ! DO i = 1, nx_dst
    DO j = 1, ny_dst
      ! Dst cell j overlaps with src cells jr_src( j,1) to ir_src( j,2)
      ycmin = y_dst( j) - dy_dst/2._dp
      ycmax = y_dst( j) + dy_dst/2._dp
      jr_src( j,:) = MAX( 1, MIN( ny_src, [CEILING(-1.5_dp + FLOOR(ny_src/2._dp) + ycmin / dy_src), &
                                           CEILING( 1.5_dp + FLOOR(ny_src/2._dp) + ycmax / dy_src)] ))
    END DO ! DO j = 1, ny_dst
    
    ! Find parallelisation domains
    CALL partition_list( nx_dst, par%i, par%n, i1, i2)
    
    DO i = i1, i2
    DO j = 1, ny_dst
    
      Aint = 0._dp
      
      DO i_src = ir_src( i,1), ir_src( i,2)
      DO j_src = jr_src( j,1), jr_src( j,2)
        
        xomin = MAX( x_dst( i) - dx_dst/2._dp, x_src( i_src) - dx_src/2._dp)
        xomax = MIN( x_dst( i) + dx_dst/2._dp, x_src( i_src) + dx_src/2._dp)
        yomin = MAX( y_dst( j) - dy_dst/2._dp, y_src( j_src) - dy_src/2._dp)
        yomax = MIN( y_dst( j) + dy_dst/2._dp, y_src( j_src) + dy_src/2._dp)
        
        IF (xomax <= xomin .OR. yomax <= yomin) CYCLE
        
        A_overlap = (xomax - xomin) * (yomax - yomin)
        Aint = Aint + A_overlap
        
        d_dst( :,j,i) = d_dst( :,j,i) + (A_overlap / A_dst) * d_src( :,j_src,i_src)
        
      END DO ! DO j_src = jr_src( j,1), jr_src( j,2)
      END DO ! DO i_src = ir_src( i,1), ir_src( i,2)
      
    END DO ! DO j = 1, ny_dst
    END DO ! DO i = 1, nx_dst
    CALL sync
  
  END SUBROUTINE map_square_to_square_3D
  
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
