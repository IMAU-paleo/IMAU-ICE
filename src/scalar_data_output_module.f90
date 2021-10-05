MODULE scalar_data_output_module

  ! Contains all the routines for writing to the scalar output NetCDF files.

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_global_scalar_data, type_model_region, type_forcing_data
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, create_global_scalar_output_file, write_to_global_scalar_output_file, &
                                             write_to_regional_scalar_output_file
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D

  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE write_regional_scalar_data( region, time)
    ! Write some regionally integrated scalar values to the NetCDF output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables:
    INTEGER                                       :: i,j,m
    REAL(dp)                                      :: T2m_mean
    REAL(dp)                                      :: total_snowfall
    REAL(dp)                                      :: total_rainfall
    REAL(dp)                                      :: total_melt
    REAL(dp)                                      :: total_refreezing
    REAL(dp)                                      :: total_runoff
    REAL(dp)                                      :: total_SMB
    REAL(dp)                                      :: total_BMB
    REAL(dp)                                      :: total_MB
    
  ! ================================================
  ! ===== Exceptions for benchmark experiments =====
  ! ================================================
  
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
        ! Do nothing; global scalar output not needed in these experiments
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in gather_global_scalar_data!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  ! =======================================================
  ! ===== End of exceptions for benchmark experiments =====
  ! =======================================================
  
  
      
    ! Regionally integrated output
    ! ============================
      
    T2m_mean                   = 0._dp
    total_snowfall             = 0._dp
    total_rainfall             = 0._dp
    total_melt                 = 0._dp
    total_refreezing           = 0._dp
    total_runoff               = 0._dp
    total_SMB                  = 0._dp
    total_BMB                  = 0._dp
    
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
    
      IF (region%ice%Hi_a( j,i) > 0._dp) THEN
        
        total_BMB = total_BMB + (region%BMB%BMB( j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
          
        DO m = 1, 12
          total_snowfall   = total_snowfall   + (region%SMB%Snowfall(   m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
          total_rainfall   = total_rainfall   + (region%SMB%Rainfall(   m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
          total_melt       = total_melt       + (region%SMB%Melt(       m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
          total_refreezing = total_refreezing + (region%SMB%Refreezing( m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
          total_runoff     = total_runoff     + (region%SMB%Runoff(     m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
          total_SMB        = total_SMB        + (region%SMB%SMB(        m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
        END DO
        
      END IF

      T2m_mean = T2m_mean + SUM(region%climate%applied%T2m( :,j,i)) / (12._dp * region%grid%nx * region%grid%ny)
      
    END DO
    END DO
    CALL sync
    
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, T2m_mean        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_snowfall  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_rainfall  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_melt      , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_refreezing, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_runoff    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_SMB       , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_BMB       , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    total_MB = total_SMB + total_BMB
    
    IF (par%master) THEN
      region%int_T2m        = T2m_mean
      region%int_snowfall   = total_snowfall
      region%int_rainfall   = total_rainfall
      region%int_melt       = total_melt
      region%int_refreezing = total_refreezing
      region%int_runoff     = total_runoff
      region%int_SMB        = total_SMB
      region%int_BMB        = total_BMB
      region%int_MB         = total_MB
    END IF
    CALL sync
    
    ! Write to NetCDF file
    CALL write_to_regional_scalar_output_file( region, time)
    
  END SUBROUTINE write_regional_scalar_data
  SUBROUTINE write_global_scalar_data( global_data, NAM, EAS, GRL, ANT, forcing, time)
    ! Collect some global scalar data values and write to the NetCDF output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_global_scalar_data),       INTENT(INOUT) :: global_data
    TYPE(type_model_region),             INTENT(IN)    :: NAM, EAS, GRL, ANT
    TYPE(type_forcing_data),             INTENT(IN)    :: forcing
    REAL(dp),                            INTENT(IN)    :: time
    
  ! ================================================
  ! ===== Exceptions for benchmark experiments =====
  ! ================================================
  
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
        ! Do nothing; global scalar output not needed in these experiments
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in gather_global_scalar_data!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  ! =======================================================
  ! ===== End of exceptions for benchmark experiments =====
  ! =======================================================
  
    IF (par%master) THEN
    
      ! Sea level: variables have already been set in the IMAU_ICE_program time loop
      
      ! CO2
      IF     (C%choice_forcing_method == 'CO2_direct') THEN
        global_data%CO2_obs        = forcing%CO2_obs
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
        global_data%CO2_obs        = forcing%CO2_obs
        global_data%CO2_mod        = forcing%CO2_mod
      ELSEIF (C%choice_forcing_method == 'climate_direct') THEN
      ELSEIF (C%choice_forcing_method == 'SMB_direct') THEN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in write_global_scalar_data!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
      ! d18O
      IF     (C%choice_forcing_method == 'CO2_direct') THEN
        global_data%d18O_mod       = forcing%d18O_mod
        global_data%d18O_ice       = forcing%d18O_from_ice_volume_mod
        global_data%d18O_Tdw       = forcing%d18O_from_temperature_mod
        global_data%d18O_NAM       = forcing%d18O_NAM
        global_data%d18O_EAS       = forcing%d18O_EAS
        global_data%d18O_GRL       = forcing%d18O_GRL
        global_data%d18O_ANT       = forcing%d18O_ANT
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
        global_data%d18O_obs       = forcing%d18O_obs
        global_data%d18O_mod       = forcing%d18O_mod
        global_data%d18O_ice       = forcing%d18O_from_ice_volume_mod
        global_data%d18O_Tdw       = forcing%d18O_from_temperature_mod
        global_data%d18O_NAM       = forcing%d18O_NAM
        global_data%d18O_EAS       = forcing%d18O_EAS
        global_data%d18O_GRL       = forcing%d18O_GRL
        global_data%d18O_ANT       = forcing%d18O_ANT
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
        global_data%d18O_obs       = forcing%d18O_obs
        global_data%d18O_mod       = forcing%d18O_mod
        global_data%d18O_ice       = forcing%d18O_from_ice_volume_mod
        global_data%d18O_Tdw       = forcing%d18O_from_temperature_mod
        global_data%d18O_NAM       = forcing%d18O_NAM
        global_data%d18O_EAS       = forcing%d18O_EAS
        global_data%d18O_GRL       = forcing%d18O_GRL
        global_data%d18O_ANT       = forcing%d18O_ANT
      ELSEIF (C%choice_forcing_method == 'climate_direct') THEN
        global_data%d18O_mod       = forcing%d18O_mod
        global_data%d18O_ice       = forcing%d18O_from_ice_volume_mod
        global_data%d18O_Tdw       = forcing%d18O_from_temperature_mod
        global_data%d18O_NAM       = forcing%d18O_NAM
        global_data%d18O_EAS       = forcing%d18O_EAS
        global_data%d18O_GRL       = forcing%d18O_GRL
        global_data%d18O_ANT       = forcing%d18O_ANT
      ELSEIF (C%choice_forcing_method == 'SMB_direct') THEN
        global_data%d18O_mod       = forcing%d18O_mod
        global_data%d18O_ice       = forcing%d18O_from_ice_volume_mod
        global_data%d18O_Tdw       = forcing%d18O_from_temperature_mod
        global_data%d18O_NAM       = forcing%d18O_NAM
        global_data%d18O_EAS       = forcing%d18O_EAS
        global_data%d18O_GRL       = forcing%d18O_GRL
        global_data%d18O_ANT       = forcing%d18O_ANT
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in write_global_scalar_data!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
      ! Temperature (surface and deep-water)
      global_data%dT_glob        = forcing%dT_glob
      global_data%dT_dw          = forcing%dT_deepwater
      
      ! Computation times
      global_data%tcomp_total    = 0._dp
      global_data%tcomp_ice      = 0._dp
      global_data%tcomp_thermo   = 0._dp
      global_data%tcomp_climate  = 0._dp
      global_data%tcomp_GIA      = 0._dp
      
      IF (C%do_NAM) THEN
        global_data%tcomp_total   = global_data%tcomp_total   + NAM%tcomp_total
        global_data%tcomp_ice     = global_data%tcomp_ice     + NAM%tcomp_ice
        global_data%tcomp_thermo  = global_data%tcomp_thermo  + NAM%tcomp_thermo
        global_data%tcomp_climate = global_data%tcomp_climate + NAM%tcomp_climate
        global_data%tcomp_GIA     = global_data%tcomp_GIA     + NAM%tcomp_GIA
      END IF
      IF (C%do_EAS) THEN
        global_data%tcomp_total   = global_data%tcomp_total   + EAS%tcomp_total
        global_data%tcomp_ice     = global_data%tcomp_ice     + EAS%tcomp_ice
        global_data%tcomp_thermo  = global_data%tcomp_thermo  + EAS%tcomp_thermo
        global_data%tcomp_climate = global_data%tcomp_climate + EAS%tcomp_climate
        global_data%tcomp_GIA     = global_data%tcomp_GIA     + EAS%tcomp_GIA
      END IF
      IF (C%do_GRL) THEN
        global_data%tcomp_total   = global_data%tcomp_total   + GRL%tcomp_total
        global_data%tcomp_ice     = global_data%tcomp_ice     + GRL%tcomp_ice
        global_data%tcomp_thermo  = global_data%tcomp_thermo  + GRL%tcomp_thermo
        global_data%tcomp_climate = global_data%tcomp_climate + GRL%tcomp_climate
        global_data%tcomp_GIA     = global_data%tcomp_GIA     + GRL%tcomp_GIA
      END IF
      IF (C%do_ANT) THEN
        global_data%tcomp_total   = global_data%tcomp_total   + ANT%tcomp_total
        global_data%tcomp_ice     = global_data%tcomp_ice     + ANT%tcomp_ice
        global_data%tcomp_thermo  = global_data%tcomp_thermo  + ANT%tcomp_thermo
        global_data%tcomp_climate = global_data%tcomp_climate + ANT%tcomp_climate
        global_data%tcomp_GIA     = global_data%tcomp_GIA     + ANT%tcomp_GIA
      END IF
      
      ! Write to output file
      CALL write_to_global_scalar_output_file( global_data, time)
      
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE write_global_scalar_data
  
  SUBROUTINE initialise_global_scalar_data(   global_data)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_global_scalar_data),       INTENT(INOUT) :: global_data
    
  ! ================================================
  ! ===== Exceptions for benchmark experiments =====
  ! ================================================
  
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
        ! Do nothing; global scalar output not needed in these experiments
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in gather_global_scalar_data!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
  ! =======================================================
  ! ===== End of exceptions for benchmark experiments =====
  ! =======================================================
    
    ! Allocate shared memory
    
    ! Sea level
    CALL allocate_shared_dp_0D( global_data%GMSL          , global_data%wGMSL          )
    CALL allocate_shared_dp_0D( global_data%GMSL_NAM      , global_data%wGMSL_NAM      )
    CALL allocate_shared_dp_0D( global_data%GMSL_EAS      , global_data%wGMSL_EAS      )
    CALL allocate_shared_dp_0D( global_data%GMSL_GRL      , global_data%wGMSL_GRL      )
    CALL allocate_shared_dp_0D( global_data%GMSL_ANT      , global_data%wGMSL_ANT      )
    
    ! CO2
    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      CALL allocate_shared_dp_0D( global_data%CO2_obs       , global_data%wCO2_obs       )
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CALL allocate_shared_dp_0D( global_data%CO2_obs       , global_data%wCO2_obs       )
      CALL allocate_shared_dp_0D( global_data%CO2_mod       , global_data%wCO2_mod       )
    ELSEIF (C%choice_forcing_method == 'climate_direct') THEN
    ELSEIF (C%choice_forcing_method == 'SMB_direct') THEN
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in initialise_global_scalar_data!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! d18O
    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      CALL allocate_shared_dp_0D( global_data%d18O_mod      , global_data%wd18O_mod      )
      CALL allocate_shared_dp_0D( global_data%d18O_ice      , global_data%wd18O_ice      )
      CALL allocate_shared_dp_0D( global_data%d18O_Tdw      , global_data%wd18O_Tdw      )
      CALL allocate_shared_dp_0D( global_data%d18O_NAM      , global_data%wd18O_NAM      )
      CALL allocate_shared_dp_0D( global_data%d18O_EAS      , global_data%wd18O_EAS      )
      CALL allocate_shared_dp_0D( global_data%d18O_GRL      , global_data%wd18O_GRL      )
      CALL allocate_shared_dp_0D( global_data%d18O_ANT      , global_data%wd18O_ANT      )
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      CALL allocate_shared_dp_0D( global_data%d18O_obs      , global_data%wd18O_obs      )
      CALL allocate_shared_dp_0D( global_data%d18O_mod      , global_data%wd18O_mod      )
      CALL allocate_shared_dp_0D( global_data%d18O_ice      , global_data%wd18O_ice      )
      CALL allocate_shared_dp_0D( global_data%d18O_Tdw      , global_data%wd18O_Tdw      )
      CALL allocate_shared_dp_0D( global_data%d18O_NAM      , global_data%wd18O_NAM      )
      CALL allocate_shared_dp_0D( global_data%d18O_EAS      , global_data%wd18O_EAS      )
      CALL allocate_shared_dp_0D( global_data%d18O_GRL      , global_data%wd18O_GRL      )
      CALL allocate_shared_dp_0D( global_data%d18O_ANT      , global_data%wd18O_ANT      )
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CALL allocate_shared_dp_0D( global_data%d18O_obs      , global_data%wd18O_obs      )
      CALL allocate_shared_dp_0D( global_data%d18O_mod      , global_data%wd18O_mod      )
      CALL allocate_shared_dp_0D( global_data%d18O_ice      , global_data%wd18O_ice      )
      CALL allocate_shared_dp_0D( global_data%d18O_Tdw      , global_data%wd18O_Tdw      )
      CALL allocate_shared_dp_0D( global_data%d18O_NAM      , global_data%wd18O_NAM      )
      CALL allocate_shared_dp_0D( global_data%d18O_EAS      , global_data%wd18O_EAS      )
      CALL allocate_shared_dp_0D( global_data%d18O_GRL      , global_data%wd18O_GRL      )
      CALL allocate_shared_dp_0D( global_data%d18O_ANT      , global_data%wd18O_ANT      )
    ELSEIF (C%choice_forcing_method == 'climate_direct') THEN
      CALL allocate_shared_dp_0D( global_data%d18O_mod      , global_data%wd18O_mod      )
      CALL allocate_shared_dp_0D( global_data%d18O_ice      , global_data%wd18O_ice      )
      CALL allocate_shared_dp_0D( global_data%d18O_Tdw      , global_data%wd18O_Tdw      )
      CALL allocate_shared_dp_0D( global_data%d18O_NAM      , global_data%wd18O_NAM      )
      CALL allocate_shared_dp_0D( global_data%d18O_EAS      , global_data%wd18O_EAS      )
      CALL allocate_shared_dp_0D( global_data%d18O_GRL      , global_data%wd18O_GRL      )
      CALL allocate_shared_dp_0D( global_data%d18O_ANT      , global_data%wd18O_ANT      )
    ELSEIF (C%choice_forcing_method == 'SMB_direct') THEN
      CALL allocate_shared_dp_0D( global_data%d18O_mod      , global_data%wd18O_mod      )
      CALL allocate_shared_dp_0D( global_data%d18O_ice      , global_data%wd18O_ice      )
      CALL allocate_shared_dp_0D( global_data%d18O_Tdw      , global_data%wd18O_Tdw      )
      CALL allocate_shared_dp_0D( global_data%d18O_NAM      , global_data%wd18O_NAM      )
      CALL allocate_shared_dp_0D( global_data%d18O_EAS      , global_data%wd18O_EAS      )
      CALL allocate_shared_dp_0D( global_data%d18O_GRL      , global_data%wd18O_GRL      )
      CALL allocate_shared_dp_0D( global_data%d18O_ANT      , global_data%wd18O_ANT      )
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in initialise_global_scalar_data!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Temperature (surface and deep-water)
    CALL allocate_shared_dp_0D( global_data%dT_glob       , global_data%wdT_glob       )
    CALL allocate_shared_dp_0D( global_data%dT_dw         , global_data%wdT_dw         )
    
    ! Computation times
    CALL allocate_shared_dp_0D( global_data%tcomp_total   , global_data%wtcomp_total   )
    CALL allocate_shared_dp_0D( global_data%tcomp_ice     , global_data%wtcomp_ice     )
    CALL allocate_shared_dp_0D( global_data%tcomp_thermo  , global_data%wtcomp_thermo  )
    CALL allocate_shared_dp_0D( global_data%tcomp_climate , global_data%wtcomp_climate )
    CALL allocate_shared_dp_0D( global_data%tcomp_GIA     , global_data%wtcomp_GIA     )
    
    ! Create the netcdf file
    CALL create_global_scalar_output_file( global_data%netcdf)
    
  END SUBROUTINE initialise_global_scalar_data
  
END MODULE scalar_data_output_module
