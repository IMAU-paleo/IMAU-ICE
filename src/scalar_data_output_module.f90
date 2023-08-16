MODULE scalar_data_output_module

  ! Contains all the routines for writing to the scalar output NetCDF files.

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_global_scalar_data, type_model_region, type_forcing_data
  USE netcdf_output_module,            ONLY: create_global_scalar_file, write_to_global_scalar_file, write_to_regional_scalar_file
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  USE netcdf_debug_module,             ONLY: save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, save_variable_as_netcdf_int_3D, &
                                             save_variable_as_netcdf_dp_1D,  save_variable_as_netcdf_dp_2D,  save_variable_as_netcdf_dp_3D

  IMPLICIT NONE

CONTAINS

  SUBROUTINE write_regional_scalar_data( region, time)
    ! Write some regionally integrated scalar values to the NetCDF output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_regional_scalar_data'
    INTEGER                                            :: i,j,m
    REAL(dp)                                           :: T2m_mean
    REAL(dp)                                           :: total_snowfall
    REAL(dp)                                           :: total_rainfall
    REAL(dp)                                           :: total_melt
    REAL(dp)                                           :: total_refreezing
    REAL(dp)                                           :: total_runoff
    REAL(dp)                                           :: total_SMB
    REAL(dp)                                           :: total_BMB
    REAL(dp)                                           :: total_MB

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (.NOT. C%do_write_regional_scalar_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Region-wide annual mean surface temperature
    IF (C%choice_climate_model == 'none') THEN
      ! In this case, no surface temperature is calculated at all
    ELSE

      T2m_mean = 0._dp

      DO i = region%grid%i1, region%grid%i2
      DO j = 1, region%grid%ny
        T2m_mean = T2m_mean + SUM( region%climate%T2m( :,j,i)) / (12._dp * region%grid%nx * region%grid%ny)
      END DO
      END DO
      CALL sync

      CALL MPI_ALLREDUCE( MPI_IN_PLACE, T2m_mean , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      IF (par%master) THEN
        region%int_T2m = T2m_mean
      END IF
      CALL sync

    END IF ! IF (C%choice_climate_model == 'none') THEN

    ! Ice-sheet-integrated surface/basal mass balance
    total_SMB = 0._dp
    total_BMB = 0._dp

    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
      IF (region%ice%mask_ice_a( j,i) == 1) THEN
        total_SMB = total_SMB + (region%SMB%SMB_year( j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
        total_BMB = total_BMB + (region%BMB%BMB(      j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
      END IF
    END DO
    END DO
    CALL sync

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, T2m_mean , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_SMB, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_BMB, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    total_MB = total_SMB + total_BMB

    IF (par%master) THEN
      region%int_SMB = total_SMB
      region%int_BMB = total_BMB
      region%int_MB  = total_MB
    END IF
    CALL sync

    ! Individual SMB components
    IF     (C%choice_SMB_model == 'uniform' .OR. &
            C%choice_SMB_model == 'idealised' .OR. &
            C%choice_SMB_model == 'direct_global' .OR. &
            C%choice_SMB_model == 'direct_regional' .OR. &
            C%choice_SMB_model == 'snapshot' .OR. &
            C%choice_SMB_model == 'ISMIP_style') THEN
      ! Do nothing
    ELSEIF (C%choice_SMB_model == 'IMAU-ITM' .OR. &
            C%choice_SMB_model == 'IMAU-ITM_wrongrefreezing') THEN

      total_snowfall             = 0._dp
      total_rainfall             = 0._dp
      total_melt                 = 0._dp
      total_refreezing           = 0._dp
      total_runoff               = 0._dp

      DO i = region%grid%i1, region%grid%i2
      DO j = 1, region%grid%ny

        IF (region%ice%Hi_a( j,i) > 0._dp) THEN

          DO m = 1, 12
            total_snowfall   = total_snowfall   + (region%SMB%Snowfall(   m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
            total_rainfall   = total_rainfall   + (region%SMB%Rainfall(   m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
            total_melt       = total_melt       + (region%SMB%Melt(       m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
            total_refreezing = total_refreezing + (region%SMB%Refreezing( m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
            total_runoff     = total_runoff     + (region%SMB%Runoff(     m,j,i) * region%grid%dx * region%grid%dx / 1E9_dp)
          END DO

        END IF

      END DO
      END DO
      CALL sync

      CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_snowfall  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_rainfall  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_melt      , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_refreezing, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_runoff    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      IF (par%master) THEN
        region%int_snowfall   = total_snowfall
        region%int_rainfall   = total_rainfall
        region%int_melt       = total_melt
        region%int_refreezing = total_refreezing
        region%int_runoff     = total_runoff
      END IF
      CALL sync

    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM( C%choice_SMB_model) // '"!')
    END IF

    ! Write to NetCDF file
    CALL write_to_regional_scalar_file( region%scalar_filename, region)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_regional_scalar_data
  SUBROUTINE write_global_scalar_data( global_data, NAM, EAS, GRL, ANT, forcing, time)
    ! Collect some global scalar data values and write to the NetCDF output file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_global_scalar_data),       INTENT(INOUT) :: global_data
    TYPE(type_model_region),             INTENT(IN)    :: NAM, EAS, GRL, ANT
    TYPE(type_forcing_data),             INTENT(IN)    :: forcing
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_global_scalar_data'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN

      ! Sea level: variables have already been set in the IMAU_ICE_program time loop

      ! CO2
      IF     (C%choice_forcing_method == 'none') THEN
      ELSEIF (C%choice_forcing_method == 'CO2_direct') THEN
        global_data%CO2_obs        = forcing%CO2_obs
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
        global_data%CO2_mod        = forcing%CO2_mod
      ELSE
        CALL crash('unknown choice_forcing_method "' // TRIM( C%choice_forcing_method) // '"!')
      END IF

      ! d18O
      IF     (C%do_calculate_benthic_d18O) THEN
        global_data%dT_glob        = forcing%dT_glob
        global_data%dT_dw          = forcing%dT_deepwater
        global_data%d18O_obs       = forcing%d18O_obs
        global_data%d18O_mod       = forcing%d18O_mod
        global_data%d18O_ice       = forcing%d18O_from_ice_volume_mod
        global_data%d18O_Tdw       = forcing%d18O_from_temperature_mod
        global_data%d18O_NAM       = forcing%d18O_NAM
        global_data%d18O_EAS       = forcing%d18O_EAS
        global_data%d18O_GRL       = forcing%d18O_GRL
        global_data%d18O_ANT       = forcing%d18O_ANT
      END IF

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

    END IF ! IF (par%master) THEN
    CALL sync

    ! Write to output file
    CALL write_to_global_scalar_file( global_data%filename, global_data, time)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_global_scalar_data

  SUBROUTINE initialise_global_scalar_data(   global_data)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_global_scalar_data),       INTENT(INOUT) :: global_data

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_global_scalar_data'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory

    ! Sea level
    CALL allocate_shared_dp_0D( global_data%GMSL          , global_data%wGMSL          )
    CALL allocate_shared_dp_0D( global_data%GMSL_NAM      , global_data%wGMSL_NAM      )
    CALL allocate_shared_dp_0D( global_data%GMSL_EAS      , global_data%wGMSL_EAS      )
    CALL allocate_shared_dp_0D( global_data%GMSL_GRL      , global_data%wGMSL_GRL      )
    CALL allocate_shared_dp_0D( global_data%GMSL_ANT      , global_data%wGMSL_ANT      )

    ! CO2
    IF     (C%choice_forcing_method == 'none') THEN
    ELSEIF (C%choice_forcing_method == 'CO2_direct') THEN
      CALL allocate_shared_dp_0D( global_data%CO2_obs       , global_data%wCO2_obs       )
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CALL allocate_shared_dp_0D( global_data%CO2_mod       , global_data%wCO2_mod       )
    ELSE
      CALL crash('unknown choice_forcing_method "' // TRIM( C%choice_forcing_method) // '"!')
    END IF

    ! d18O
    IF     (C%do_calculate_benthic_d18O) THEN
      ! Temperature (surface and deep-water)
      CALL allocate_shared_dp_0D( global_data%dT_glob       , global_data%wdT_glob       )
      CALL allocate_shared_dp_0D( global_data%dT_dw         , global_data%wdT_dw         )
    
      ! d18O      
      CALL allocate_shared_dp_0D( global_data%d18O_obs      , global_data%wd18O_obs      )
      CALL allocate_shared_dp_0D( global_data%d18O_mod      , global_data%wd18O_mod      )
      CALL allocate_shared_dp_0D( global_data%d18O_ice      , global_data%wd18O_ice      )
      CALL allocate_shared_dp_0D( global_data%d18O_Tdw      , global_data%wd18O_Tdw      )
      CALL allocate_shared_dp_0D( global_data%d18O_NAM      , global_data%wd18O_NAM      )
      CALL allocate_shared_dp_0D( global_data%d18O_EAS      , global_data%wd18O_EAS      )
      CALL allocate_shared_dp_0D( global_data%d18O_GRL      , global_data%wd18O_GRL      )
      CALL allocate_shared_dp_0D( global_data%d18O_ANT      , global_data%wd18O_ANT      )
    END IF

    ! Computation times
    CALL allocate_shared_dp_0D( global_data%tcomp_total   , global_data%wtcomp_total   )
    CALL allocate_shared_dp_0D( global_data%tcomp_ice     , global_data%wtcomp_ice     )
    CALL allocate_shared_dp_0D( global_data%tcomp_thermo  , global_data%wtcomp_thermo  )
    CALL allocate_shared_dp_0D( global_data%tcomp_climate , global_data%wtcomp_climate )
    CALL allocate_shared_dp_0D( global_data%tcomp_GIA     , global_data%wtcomp_GIA     )

    ! Create the netcdf file
    CALL create_global_scalar_file( global_data%filename)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_global_scalar_data

END MODULE scalar_data_output_module
