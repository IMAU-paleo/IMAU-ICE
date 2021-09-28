MODULE netcdf_module

  ! Contains routines for creating, reading, and writing all the NetCDF files
  ! involved in the ice model.

  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D
  USE netcdf,                          ONLY: nf90_max_var_dims, nf90_create, nf90_close, nf90_clobber, nf90_share, nf90_unlimited , &
                                             nf90_enddef, nf90_put_var, nf90_sync, nf90_def_var, nf90_int, nf90_put_att, nf90_def_dim, &
                                             nf90_open, nf90_write, nf90_inq_dimid, nf90_inquire_dimension, nf90_inquire, nf90_double, &
                                             nf90_inq_varid, nf90_inquire_variable, nf90_get_var, nf90_noerr, nf90_strerror, nf90_float
  USE data_types_module,               ONLY: type_model_region, type_PD_data_fields, type_forcing_data, type_ICE5G_timeframe, &
                                             type_init_data_fields, type_subclimate_global, type_debug_fields, type_SELEN_global, &
                                             type_netcdf_scalars_global, type_netcdf_scalars_regional, type_global_scalar_data
  IMPLICIT NONE
  
  TYPE(type_debug_fields) :: debug_NAM, debug_EAS, debug_GRL, debug_ANT, debug

CONTAINS

  ! Write to output NetCDF files
  SUBROUTINE write_to_restart_file( region, forcing)
    ! Write the current model state to the existing output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_forcing_data),        INTENT(IN)    :: forcing
    
    ! Local variables:
    INTEGER                                       :: ncid, nx, ny, nz, ti
    
    IF (.NOT. par%master) RETURN
    
    WRITE(0,'(A,F8.2,A)') '   t = ', region%time/1E3, ' kyr - writing output...'
    
    ! Open the file for writing
    CALL open_netcdf_file( region%restart%filename, region%restart%ncid)
        
    ! Time
    CALL handle_error( nf90_put_var( region%restart%ncid, region%restart%id_var_time,             region%time,                    start=(/       region%restart%ti/)))
    
    ! Placeholders for the dimension ID's, for shorter code
    ncid = region%restart%ncid
    nx   = region%grid%nx
    ny   = region%grid%ny
    nz   = C%nZ
    ti   = region%restart%ti
    
    CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%id_var_Hi,               region%ice%Hi_a,             (/1, 1,    ti/))
    CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%id_var_Hb,               region%ice%Hb_a,             (/1, 1,    ti/))
    CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%id_var_dHb,              region%ice%dHb_a,            (/1, 1,    ti/))
    CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%id_var_SL,               region%ice%SL_a,             (/1, 1,    ti/))
    CALL write_data_to_file_dp_3D( ncid, nx, ny, nz, region%restart%id_var_Ti,               region%ice%Ti_a,             (/1, 1, 1, ti/))
    CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, region%restart%id_var_FirnDepth,        region%SMB%FirnDepth,        (/1, 1, 1, ti/))
    CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%id_var_MeltPreviousYear, region%SMB%MeltPreviousYear, (/1, 1,    ti/))
    
    ! Inverse routine data
    
    ! Exceptions for benchmark experiments
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
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in write_to_restart_file!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSE
    
      IF     (C%choice_forcing_method == 'CO2_direct' .OR. C%choice_forcing_method == 'SMB_direct' .OR. C%choice_forcing_method == 'climate_direct') THEN
        ! No inverse routine used in these forcing methods
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
        ! Need to write dT_glob_history and dT_glob_inverse_history
        
        CALL handle_error( nf90_put_var( ncid, region%restart%id_var_dT_glob_history,         forcing%dT_glob_history,         start=(/1, ti/) ))
        CALL handle_error( nf90_put_var( ncid, region%restart%id_var_dT_glob_inverse_history, forcing%dT_glob_inverse_history, start=(/1, ti/) ))
        
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
        ! Need to write dT_glob_history and CO2_inverse_history
        
        CALL handle_error( nf90_put_var( ncid, region%restart%id_var_dT_glob_history,         forcing%dT_glob_history,         start=(/1, ti/) ))
        CALL handle_error( nf90_put_var( ncid, region%restart%id_var_CO2_inverse_history,     forcing%CO2_inverse_history,     start=(/1, ti/) ))
        
      ELSE
        WRITE(0,*) '  ERROR - choice_forcing_method "', C%choice_forcing_method, '" not implemented in write_to_restart_file!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Close the file
    CALL close_netcdf_file(region%restart%ncid)
    
    ! Increase time frame counter
    region%restart%ti = region%restart%ti + 1
        
  END SUBROUTINE write_to_restart_file
  SUBROUTINE write_to_help_fields_file( region)
    ! Write the current model state to the existing output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region), INTENT(INOUT) :: region
    
    IF (.NOT. par%master) RETURN
    
    ! Open the file for writing
    CALL open_netcdf_file( region%help_fields%filename, region%help_fields%ncid)
        
    ! Time
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_time,             region%time,                    start=(/       region%help_fields%ti/)))
    
    ! Data
    CALL write_help_field( region, region%help_fields%id_help_field_01, C%help_field_01)
    CALL write_help_field( region, region%help_fields%id_help_field_02, C%help_field_02)
    CALL write_help_field( region, region%help_fields%id_help_field_03, C%help_field_03)
    CALL write_help_field( region, region%help_fields%id_help_field_04, C%help_field_04)
    CALL write_help_field( region, region%help_fields%id_help_field_05, C%help_field_05)
    CALL write_help_field( region, region%help_fields%id_help_field_06, C%help_field_06)
    CALL write_help_field( region, region%help_fields%id_help_field_07, C%help_field_07)
    CALL write_help_field( region, region%help_fields%id_help_field_08, C%help_field_08)
    CALL write_help_field( region, region%help_fields%id_help_field_09, C%help_field_09)
    CALL write_help_field( region, region%help_fields%id_help_field_10, C%help_field_10)
    CALL write_help_field( region, region%help_fields%id_help_field_11, C%help_field_11)
    CALL write_help_field( region, region%help_fields%id_help_field_12, C%help_field_12)
    CALL write_help_field( region, region%help_fields%id_help_field_13, C%help_field_13)
    CALL write_help_field( region, region%help_fields%id_help_field_14, C%help_field_14)
    CALL write_help_field( region, region%help_fields%id_help_field_15, C%help_field_15)
    CALL write_help_field( region, region%help_fields%id_help_field_16, C%help_field_16)
    CALL write_help_field( region, region%help_fields%id_help_field_17, C%help_field_17)
    CALL write_help_field( region, region%help_fields%id_help_field_18, C%help_field_18)
    CALL write_help_field( region, region%help_fields%id_help_field_19, C%help_field_19)
    CALL write_help_field( region, region%help_fields%id_help_field_20, C%help_field_20)
    CALL write_help_field( region, region%help_fields%id_help_field_21, C%help_field_21)
    CALL write_help_field( region, region%help_fields%id_help_field_22, C%help_field_22)
    CALL write_help_field( region, region%help_fields%id_help_field_23, C%help_field_23)
    CALL write_help_field( region, region%help_fields%id_help_field_24, C%help_field_24)
    CALL write_help_field( region, region%help_fields%id_help_field_25, C%help_field_25)
    CALL write_help_field( region, region%help_fields%id_help_field_26, C%help_field_26)
    CALL write_help_field( region, region%help_fields%id_help_field_27, C%help_field_27)
    CALL write_help_field( region, region%help_fields%id_help_field_28, C%help_field_28)
    CALL write_help_field( region, region%help_fields%id_help_field_29, C%help_field_29)
    CALL write_help_field( region, region%help_fields%id_help_field_30, C%help_field_30)
    CALL write_help_field( region, region%help_fields%id_help_field_31, C%help_field_31)
    CALL write_help_field( region, region%help_fields%id_help_field_32, C%help_field_32)
    CALL write_help_field( region, region%help_fields%id_help_field_33, C%help_field_33)
    CALL write_help_field( region, region%help_fields%id_help_field_34, C%help_field_34)
    CALL write_help_field( region, region%help_fields%id_help_field_35, C%help_field_35)
    CALL write_help_field( region, region%help_fields%id_help_field_36, C%help_field_36)
    CALL write_help_field( region, region%help_fields%id_help_field_37, C%help_field_37)
    CALL write_help_field( region, region%help_fields%id_help_field_38, C%help_field_38)
    CALL write_help_field( region, region%help_fields%id_help_field_39, C%help_field_39)
    CALL write_help_field( region, region%help_fields%id_help_field_40, C%help_field_40)
    CALL write_help_field( region, region%help_fields%id_help_field_41, C%help_field_41)
    CALL write_help_field( region, region%help_fields%id_help_field_42, C%help_field_42)
    CALL write_help_field( region, region%help_fields%id_help_field_43, C%help_field_43)
    CALL write_help_field( region, region%help_fields%id_help_field_44, C%help_field_44)
    CALL write_help_field( region, region%help_fields%id_help_field_45, C%help_field_45)
    CALL write_help_field( region, region%help_fields%id_help_field_46, C%help_field_46)
    CALL write_help_field( region, region%help_fields%id_help_field_47, C%help_field_47)
    CALL write_help_field( region, region%help_fields%id_help_field_48, C%help_field_48)
    CALL write_help_field( region, region%help_fields%id_help_field_49, C%help_field_49)
    CALL write_help_field( region, region%help_fields%id_help_field_50, C%help_field_50)

    ! Close the file
    CALL close_netcdf_file(region%help_fields%ncid)
    
    ! Increase time frame counter
    region%help_fields%ti = region%help_fields%ti + 1
        
  END SUBROUTINE write_to_help_fields_file
  SUBROUTINE write_help_field( region, id_var, field_name)
    ! Write a single data field to the help_fields file
    ! (NOTE: if you want to add an extra option, be sure to also add it to create_help_field!)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(IN)    :: region
    INTEGER,                        INTENT(IN)    :: id_var
    CHARACTER(LEN=*),               INTENT(IN)    :: field_name
    
    ! Local variables:
    INTEGER                                       :: ncid, nx, ny, ti, nz, i, j
    REAL(dp), DIMENSION( region%grid%ny, region%grid%nx) :: d_temp
    
    IF (.NOT. par%master) RETURN
    
    ! Placeholders for the dimension ID's, for shorter code
    ncid = region%help_fields%ncid
    nx   = region%grid%nx
    ny   = region%grid%ny
    nz   = C%nZ
    ti   = region%help_fields%ti
    
    IF (field_name == 'none') THEN
      RETURN
      
    ! Fields with no time dimension
    ! =============================
      
    ! Lat/lon
    ELSEIF (field_name == 'lat') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%grid%lat,           (/1, 1        /))
    ELSEIF (field_name == 'lon') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%grid%lon,           (/1, 1        /))
    
    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%GHF_a,          (/1, 1        /))

    ! Forcing climates
    ELSEIF (field_name == 'GCM_Warm_T2m') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate%GCM_warm%T2m,    (/1, 1, 1/))
    ELSEIF (field_name == 'GCM_Warm_Precip') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate%GCM_warm%Precip, (/1, 1, 1/))
    ELSEIF (field_name == 'GCM_Cold_T2m') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate%GCM_cold%T2m,    (/1, 1, 1/))
    ELSEIF (field_name == 'GCM_Cold_Precip') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate%GCM_cold%Precip, (/1, 1, 1/))
    ELSEIF (field_name == 'GCM_PI_T2m') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate%GCM_PI%T2m,      (/1, 1, 1/))
    ELSEIF (field_name == 'GCM_PI_Precip') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate%GCM_PI%Precip,   (/1, 1, 1 /))
    ELSEIF (field_name == 'PD_obs_T2m') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate%PD_obs%T2m,      (/1, 1, 1/))
    ELSEIF (field_name == 'PD_obs_Precip') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate%PD_obs%Precip,   (/1, 1, 1/))
      
    ! Fields with a time dimension
    ! ============================
      
    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%Hi_a,           (/1, 1,    ti /))
    ELSEIF (field_name == 'Hb') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%Hb_a,           (/1, 1,    ti /))
    ELSEIF (field_name == 'Hs') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%Hs_a,           (/1, 1,    ti /))
    ELSEIF (field_name == 'SL') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%SL_a,           (/1, 1,    ti /))
    ELSEIF (field_name == 'dHs_dx') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%dHs_dx_a,       (/1, 1,    ti /))
    ELSEIF (field_name == 'dHs_dy') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%dHs_dy_a,       (/1, 1,    ti /))
      
    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nz, id_var,               region%ice%Ti_a,           (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Cpi') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nz, id_var,               region%ice%Cpi_a,          (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Ki') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nz, id_var,               region%ice%Ki_a,           (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Ti_basal') THEN
      d_temp = region%ice%Ti_a( C%nZ,:,:)
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                    (/1, 1,    ti /))
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nz, id_var,               region%ice%Ti_pmp_a,       (/1, 1, 1, ti /))
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nz, id_var,               region%ice%A_flow_3D_a,    (/1, 1, 1, ti /))
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nz, id_var,               region%ice%A_flow_vav_a,   (/1, 1, 1, ti /))
      
    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nz, id_var,               region%ice%u_3D_a,         (/1, 1, 1, ti /))
    ELSEIF (field_name == 'v_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nz, id_var,               region%ice%v_3D_a,         (/1, 1, 1, ti /))
    ELSEIF (field_name == 'w_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nz, id_var,               region%ice%w_3D_a,         (/1, 1, 1, ti /))
    ELSEIF (field_name == 'u_vav') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%u_vav_a,        (/1, 1,    ti /))
    ELSEIF (field_name == 'v_vav') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%v_vav_a,        (/1, 1,    ti /))
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%uabs_vav_a,     (/1, 1,    ti /))
    ELSEIF (field_name == 'u_surf') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%u_surf_a,       (/1, 1,    ti /))
    ELSEIF (field_name == 'v_surf') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%v_surf_a,       (/1, 1,    ti /))
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%uabs_surf_a,    (/1, 1,    ti /))
    ELSEIF (field_name == 'u_base') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%u_base_a,       (/1, 1,    ti /))
    ELSEIF (field_name == 'v_base') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%v_base_a,       (/1, 1,    ti /))
    ELSEIF (field_name == 'uabs_base') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%uabs_base_a,    (/1, 1,    ti /))
      
    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate%applied%T2m, (/1, 1, 1, ti /))
    ELSEIF (field_name == 'T2m_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%climate%applied%T2m( :,j,i)) / 12._dp
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))
    ELSEIF (field_name == 'Precip') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate%applied%Precip, (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Precip_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%climate%applied%Precip( :,j,i))
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate%applied%Wind_WE, (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Wind_WE_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%climate%applied%Wind_WE( :,j,i)) / 12._dp
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate%applied%Wind_SN, (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Wind_SN_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%climate%applied%Wind_SN( :,j,i)) / 12._dp
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))    
      
    ! Mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%SMB%SMB,       (/1, 1, 1, ti /))
    ELSEIF (field_name == 'SMB_year') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%SMB%SMB_year,  (/1, 1,    ti /))
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%BMB%BMB_sheet, (/1, 1,   ti /))
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%BMB%BMB_shelf, (/1, 1,   ti /))
    ELSEIF (field_name == 'BMB') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%BMB%BMB,       (/1, 1,   ti /))
    ELSEIF (field_name == 'Snowfall') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%SMB%Snowfall,  (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Snowfall_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%SMB%Snowfall( :,j,i))
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))
    ELSEIF (field_name == 'Rainfall') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%SMB%Rainfall,  (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Rainfall_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%SMB%Rainfall( :,j,i))
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%SMB%AddedFirn,  (/1, 1, 1, ti /))
    ELSEIF (field_name == 'AddedFirn_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%SMB%AddedFirn( :,j,i))
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))
    ELSEIF (field_name == 'Refreezing') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%SMB%Refreezing,       (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%SMB%Refreezing_year,  (/1, 1,    ti /))
    ELSEIF (field_name == 'Runoff') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%SMB%Runoff,       (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Runoff_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%SMB%Runoff( :,j,i))
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))
    ELSEIF (field_name == 'Albedo') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%SMB%Albedo,  (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Albedo_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%SMB%Albedo( :,j,i)) / 12._dp
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%SMB%FirnDepth,       (/1, 1, 1, ti /))
    ELSEIF (field_name == 'FirnDepth_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%SMB%FirnDepth( :,j,i)) / 12._dp
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))
      
    ! Masks
    ELSEIF (field_name == 'mask') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_a,                     (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_land') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_land_a,                (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_ocean') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_ocean_a,               (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_lake') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_lake_a,                (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_ice') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_ice_a,                 (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_sheet') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_sheet_a,               (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_shelf') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_shelf_a,               (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_coast') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_coast_a,               (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_margin') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_margin_a,              (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_gl') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_gl_a,                  (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_cf') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_cf_a,                  (/1, 1,    ti /))
      
    ! Basal conditions
    ELSEIF (field_name == 'A_slid') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%A_slid_a,         (/1, 1,    ti /))
    ELSEIF (field_name == 'phi_fric') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%phi_fric_a,       (/1, 1,    ti /))
    ELSEIF (field_name == 'tau_yield') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%tauc_a,           (/1, 1,    ti /))
      
    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%IsoIce,            (/1, 1,    ti /))
    ELSEIF (field_name == 'iso_surf') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%IsoSurf,           (/1, 1,    ti /))
    
    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%dHb_a, (/1, 1,    ti /))
    
    ELSE
      WRITE(0,*) ' ERROR: help field "', TRIM(field_name), '" not implemented in write_help_field!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE write_help_field
  SUBROUTINE write_to_debug_file
    ! Write the current set of debug data fields to the debug NetCDF file
   
    IMPLICIT NONE
    
    ! Local variables:
    INTEGER                                :: ncid, nx, ny
    
    IF (.NOT. par%master) RETURN
    
    IF (.NOT. C%do_write_debug_data) RETURN
    
    nx = debug%nx
    ny = debug%ny
    
    ! Open the file for writing
    CALL open_netcdf_file( debug%netcdf%filename, ncid)
    
    CALL write_data_to_file_int_2D( ncid, nx, ny,       debug%netcdf%id_var_int_2D_01,        debug%int_2D_01,         (/1, 1    /))
    CALL write_data_to_file_int_2D( ncid, nx, ny,       debug%netcdf%id_var_int_2D_02,        debug%int_2D_02,         (/1, 1    /))
    CALL write_data_to_file_int_2D( ncid, nx, ny,       debug%netcdf%id_var_int_2D_03,        debug%int_2D_03,         (/1, 1    /))
    CALL write_data_to_file_int_2D( ncid, nx, ny,       debug%netcdf%id_var_int_2D_04,        debug%int_2D_04,         (/1, 1    /))
    CALL write_data_to_file_int_2D( ncid, nx, ny,       debug%netcdf%id_var_int_2D_05,        debug%int_2D_05,         (/1, 1    /))
    CALL write_data_to_file_int_2D( ncid, nx, ny,       debug%netcdf%id_var_int_2D_06,        debug%int_2D_06,         (/1, 1    /))
    CALL write_data_to_file_int_2D( ncid, nx, ny,       debug%netcdf%id_var_int_2D_07,        debug%int_2D_07,         (/1, 1    /))
    CALL write_data_to_file_int_2D( ncid, nx, ny,       debug%netcdf%id_var_int_2D_08,        debug%int_2D_08,         (/1, 1    /))
    CALL write_data_to_file_int_2D( ncid, nx, ny,       debug%netcdf%id_var_int_2D_09,        debug%int_2D_09,         (/1, 1    /))
    CALL write_data_to_file_int_2D( ncid, nx, ny,       debug%netcdf%id_var_int_2D_10,        debug%int_2D_10,         (/1, 1    /))

    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_01,         debug%dp_2D_01,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_02,         debug%dp_2D_02,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_03,         debug%dp_2D_03,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_04,         debug%dp_2D_04,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_05,         debug%dp_2D_05,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_06,         debug%dp_2D_06,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_07,         debug%dp_2D_07,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_08,         debug%dp_2D_08,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_09,         debug%dp_2D_09,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_10,         debug%dp_2D_10,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_11,         debug%dp_2D_11,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_12,         debug%dp_2D_12,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_13,         debug%dp_2D_13,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_14,         debug%dp_2D_14,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_15,         debug%dp_2D_15,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_16,         debug%dp_2D_16,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_17,         debug%dp_2D_17,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_18,         debug%dp_2D_18,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_19,         debug%dp_2D_19,          (/1, 1    /))
    CALL write_data_to_file_dp_2D(  ncid, nx, ny,       debug%netcdf%id_var_dp_2D_20,         debug%dp_2D_20,          (/1, 1    /))
    
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, C%nZ, debug%netcdf%id_var_dp_3D_01,         debug%dp_3D_01,          (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, C%nZ, debug%netcdf%id_var_dp_3D_02,         debug%dp_3D_02,          (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, C%nZ, debug%netcdf%id_var_dp_3D_03,         debug%dp_3D_03,          (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, C%nZ, debug%netcdf%id_var_dp_3D_04,         debug%dp_3D_04,          (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, C%nZ, debug%netcdf%id_var_dp_3D_05,         debug%dp_3D_05,          (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, C%nZ, debug%netcdf%id_var_dp_3D_06,         debug%dp_3D_06,          (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, C%nZ, debug%netcdf%id_var_dp_3D_07,         debug%dp_3D_07,          (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, C%nZ, debug%netcdf%id_var_dp_3D_08,         debug%dp_3D_08,          (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, C%nZ, debug%netcdf%id_var_dp_3D_09,         debug%dp_3D_09,          (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, C%nZ, debug%netcdf%id_var_dp_3D_10,         debug%dp_3D_10,          (/1, 1, 1 /))
    
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, 12,   debug%netcdf%id_var_dp_2D_monthly_01, debug%dp_2D_monthly_01,  (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, 12,   debug%netcdf%id_var_dp_2D_monthly_02, debug%dp_2D_monthly_02,  (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, 12,   debug%netcdf%id_var_dp_2D_monthly_03, debug%dp_2D_monthly_03,  (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, 12,   debug%netcdf%id_var_dp_2D_monthly_04, debug%dp_2D_monthly_04,  (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, 12,   debug%netcdf%id_var_dp_2D_monthly_05, debug%dp_2D_monthly_05,  (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, 12,   debug%netcdf%id_var_dp_2D_monthly_06, debug%dp_2D_monthly_06,  (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, 12,   debug%netcdf%id_var_dp_2D_monthly_07, debug%dp_2D_monthly_07,  (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, 12,   debug%netcdf%id_var_dp_2D_monthly_08, debug%dp_2D_monthly_08,  (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, 12,   debug%netcdf%id_var_dp_2D_monthly_09, debug%dp_2D_monthly_09,  (/1, 1, 1 /))
    CALL write_data_to_file_dp_3D(  ncid, nx, ny, 12,   debug%netcdf%id_var_dp_2D_monthly_10, debug%dp_2D_monthly_10,  (/1, 1, 1 /))

    ! Close the file
    CALL close_netcdf_file( ncid)
        
  END SUBROUTINE write_to_debug_file
  SUBROUTINE write_to_SELEN_output_file( SELEN, time)
    ! Write the current model state to the existing output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN
    REAL(dp),                       INTENT(IN)    :: time
    
    IF (.NOT. par%master) RETURN
    
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
        
  END SUBROUTINE write_to_SELEN_output_file
  SUBROUTINE write_to_global_scalar_output_file( global_data, time)
    ! Write data to the global scalar output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_global_scalar_data),   INTENT(INOUT) :: global_data
    REAL(dp),                        INTENT(IN)    :: time
    
    IF (.NOT. par%master) RETURN
    
    ! Open the file for writing
    CALL open_netcdf_file( global_data%netcdf%filename, global_data%netcdf%ncid)
        
    ! Time
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_time,             time,                       start = (/global_data%netcdf%ti/)))
    
    ! Sea level
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_GMSL,             global_data%GMSL,           start = (/global_data%netcdf%ti/)))
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_GMSL_NAM,         global_data%GMSL_NAM,       start = (/global_data%netcdf%ti/)))
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_GMSL_EAS,         global_data%GMSL_EAS,       start = (/global_data%netcdf%ti/)))
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_GMSL_GRL,         global_data%GMSL_GRL,       start = (/global_data%netcdf%ti/)))
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_GMSL_ANT,         global_data%GMSL_ANT,       start = (/global_data%netcdf%ti/)))
    
    ! CO2
    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_CO2_obs,          global_data%CO2_obs,        start = (/global_data%netcdf%ti/)))
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_CO2_obs,          global_data%CO2_obs,        start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_CO2_mod,          global_data%CO2_mod,        start = (/global_data%netcdf%ti/)))
    ELSEIF (C%choice_forcing_method == 'climate_direct') THEN
    ELSEIF (C%choice_forcing_method == 'SMB_direct') THEN
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in write_to_global_scalar_output_file!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! d18O
    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_mod,         global_data%d18O_mod,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ice,         global_data%d18O_ice,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_Tdw,         global_data%d18O_Tdw,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_NAM,         global_data%d18O_NAM,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_EAS,         global_data%d18O_EAS,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_GRL,         global_data%d18O_GRL,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ANT,         global_data%d18O_ANT,       start = (/global_data%netcdf%ti/)))
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_obs,         global_data%d18O_obs,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_mod,         global_data%d18O_mod,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ice,         global_data%d18O_ice,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_Tdw,         global_data%d18O_Tdw,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_NAM,         global_data%d18O_NAM,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_EAS,         global_data%d18O_EAS,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_GRL,         global_data%d18O_GRL,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ANT,         global_data%d18O_ANT,       start = (/global_data%netcdf%ti/)))
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_obs,         global_data%d18O_obs,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_mod,         global_data%d18O_mod,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ice,         global_data%d18O_ice,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_Tdw,         global_data%d18O_Tdw,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_NAM,         global_data%d18O_NAM,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_EAS,         global_data%d18O_EAS,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_GRL,         global_data%d18O_GRL,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ANT,         global_data%d18O_ANT,       start = (/global_data%netcdf%ti/)))
    ELSEIF (C%choice_forcing_method == 'climate_direct') THEN
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_mod,         global_data%d18O_mod,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ice,         global_data%d18O_ice,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_Tdw,         global_data%d18O_Tdw,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_NAM,         global_data%d18O_NAM,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_EAS,         global_data%d18O_EAS,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_GRL,         global_data%d18O_GRL,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ANT,         global_data%d18O_ANT,       start = (/global_data%netcdf%ti/)))
    ELSEIF (C%choice_forcing_method == 'SMB_direct') THEN
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_mod,         global_data%d18O_mod,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ice,         global_data%d18O_ice,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_Tdw,         global_data%d18O_Tdw,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_NAM,         global_data%d18O_NAM,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_EAS,         global_data%d18O_EAS,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_GRL,         global_data%d18O_GRL,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ANT,         global_data%d18O_ANT,       start = (/global_data%netcdf%ti/)))
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in write_to_global_scalar_output_file!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Temperature (surface and deep-water)
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_dT_glob,          global_data%dT_glob,        start = (/global_data%netcdf%ti/)))
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_dT_dw,            global_data%dT_dw,          start = (/global_data%netcdf%ti/)))
    
    ! Computation time for different model components
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_tcomp_total,      global_data%tcomp_total,    start = (/global_data%netcdf%ti/)))
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_tcomp_ice,        global_data%tcomp_ice,      start = (/global_data%netcdf%ti/)))
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_tcomp_thermo,     global_data%tcomp_thermo,   start = (/global_data%netcdf%ti/)))
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_tcomp_climate,    global_data%tcomp_climate,  start = (/global_data%netcdf%ti/)))
    CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_tcomp_GIA  ,      global_data%tcomp_GIA,      start = (/global_data%netcdf%ti/)))
    
    ! Close the file
    CALL close_netcdf_file(global_data%netcdf%ncid)
    
    ! Increase time frame counter
    global_data%netcdf%ti = global_data%netcdf%ti + 1
        
  END SUBROUTINE write_to_global_scalar_output_file
  SUBROUTINE write_to_regional_scalar_output_file( region, time)
    ! Write data to the regional scalar output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),         INTENT(INOUT) :: region
    REAL(dp),                        INTENT(IN)    :: time
    
    IF (.NOT. par%master) RETURN
    
    ! Open the file for writing
    CALL open_netcdf_file( region%scalars%filename, region%scalars%ncid)
        
    ! Time
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_time,             time,                       start = (/region%scalars%ti/)))
    
    ! Variables
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_ice_volume,    region%ice_volume,                 start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_ice_volume_af, region%ice_volume_above_flotation, start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_ice_area,      region%ice_area,                   start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_T2m,           region%int_T2m,                    start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_snowfall,      region%int_snowfall,               start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_rainfall,      region%int_rainfall,               start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_melt,          region%int_melt,                   start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_refreezing,    region%int_refreezing,             start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_runoff,        region%int_runoff,                 start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_SMB,           region%int_SMB,                    start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_BMB,           region%int_BMB,                    start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_MB,            region%int_MB,                     start = (/region%scalars%ti/)))
    
    ! Close the file
    CALL close_netcdf_file(region%scalars%ncid)
    
    ! Increase time frame counter
    region%scalars%ti = region%scalars%ti + 1
        
  END SUBROUTINE write_to_regional_scalar_output_file
  
  ! Create output netcdf files
  SUBROUTINE create_restart_file( region, forcing)
    ! Create a new restart file, containing the key model output required to start another run.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_forcing_data),        INTENT(IN)    :: forcing

    ! Local variables:
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x, y, z, m, t
    
    IF (.NOT. par%master) RETURN
    
    ! Set time frame index to 1
    region%restart%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(region%restart%filename))
    IF (file_exists) THEN
      WRITE(0,*) '  create_restart_file - ERROR: ', TRIM(region%restart%filename), ' already exists!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( region%restart%filename)
    CALL handle_error(nf90_create(region%restart%filename,IOR(nf90_clobber,nf90_share),region%restart%ncid))
        
    ! Define dimensions:
    CALL create_dim( region%restart%ncid, region%restart%name_dim_x,         region%grid%nx,     region%restart%id_dim_x         )
    CALL create_dim( region%restart%ncid, region%restart%name_dim_y,         region%grid%ny,     region%restart%id_dim_y         )
    CALL create_dim( region%restart%ncid, region%restart%name_dim_zeta,      C%nZ,               region%restart%id_dim_zeta      ) ! Scaled vertical coordinate
    CALL create_dim( region%restart%ncid, region%restart%name_dim_month,     12,                 region%restart%id_dim_month     ) ! Months (for monthly data)
    CALL create_dim( region%restart%ncid, region%restart%name_dim_time,      nf90_unlimited,     region%restart%id_dim_time      ) ! Time frames
    
    ! Placeholders for the dimension ID's, for shorter code
    x = region%restart%id_dim_x
    y = region%restart%id_dim_y
    z = region%restart%id_dim_zeta
    m = region%restart%id_dim_month
    t = region%restart%id_dim_time
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.
    
    ! Dimension variables: zeta, month, time
    CALL create_double_var( region%restart%ncid, region%restart%name_var_x,                [x            ], region%restart%id_var_x,                long_name='X-coordinate', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_y,                [   y         ], region%restart%id_var_y,                long_name='Y-coordinate', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_zeta,             [      z      ], region%restart%id_var_zeta,             long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_month,            [         m   ], region%restart%id_var_month,            long_name='Month', units='1-12'    )
    CALL create_double_var( region%restart%ncid, region%restart%name_var_time,             [            t], region%restart%id_var_time,             long_name='Time', units='years'   )
    
    ! Ice model data
    CALL create_double_var( region%restart%ncid, region%restart%name_var_Hi,               [x, y,       t], region%restart%id_var_Hi,               long_name='Ice thickness', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_Hb,               [x, y,       t], region%restart%id_var_Hb,               long_name='Bedrock elevation', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_SL,               [x, y,       t], region%restart%id_var_SL,               long_name='Sea-level change', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_dHb,              [x, y,       t], region%restart%id_var_dHb,              long_name='Bedrock deformation', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_Ti,               [x, y, z,    t], region%restart%id_var_Ti,               long_name='Ice temperature', units='K')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_FirnDepth,        [x, y,    m, t], region%restart%id_var_FirnDepth,        long_name='Firn depth', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_MeltPreviousYear, [x, y,       t], region%restart%id_var_MeltPreviousYear, long_name='Melt during previous year', units='mie')
    
    ! Inverse routine data
    IF (C%do_benchmark_experiment) THEN
      ! Exceptions for benchmark experiments
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
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in create_restart_file!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSE
    
      IF     (C%choice_forcing_method == 'CO2_direct' .OR. C%choice_forcing_method == 'SMB_direct' .OR. C%choice_forcing_method == 'climate_direct') THEN
        ! No inverse routine used in these forcing methods
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
        ! Need to write dT_glob_history and dT_glob_inverse_history
        
        CALL create_dim( region%restart%ncid, region%restart%name_dim_ndT_glob_history,         forcing%ndT_glob_history,         region%restart%id_dim_ndT_glob_history        )
        CALL create_dim( region%restart%ncid, region%restart%name_dim_ndT_glob_inverse_history, forcing%ndT_glob_inverse_history, region%restart%id_dim_ndT_glob_inverse_history)
        
        CALL create_double_var( region%restart%ncid, region%restart%name_var_dT_glob_history,         [region%restart%id_dim_ndT_glob_history,         t], region%restart%id_var_dT_glob_history,         long_name='dT_glob history')
        CALL create_double_var( region%restart%ncid, region%restart%name_var_dT_glob_inverse_history, [region%restart%id_dim_ndT_glob_inverse_history, t], region%restart%id_var_dT_glob_inverse_history, long_name='dT_glob_inverse history')
        
      ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
        ! Need to write dT_glob_history and CO2_inverse_history
        
        CALL create_dim( region%restart%ncid, region%restart%name_dim_ndT_glob_history,         forcing%ndT_glob_history,         region%restart%id_dim_ndT_glob_history    )
        CALL create_dim( region%restart%ncid, region%restart%name_dim_nCO2_inverse_history,     forcing%nCO2_inverse_history,     region%restart%id_dim_nCO2_inverse_history)
        
        CALL create_double_var( region%restart%ncid, region%restart%name_var_dT_glob_history,         [region%restart%id_dim_ndT_glob_history,         t], region%restart%id_var_dT_glob_history,     long_name='dT_glob history')
        CALL create_double_var( region%restart%ncid, region%restart%name_var_CO2_inverse_history,     [region%restart%id_dim_nCO2_inverse_history,     t], region%restart%id_var_CO2_inverse_history, long_name='CO2_inverse history')
        
      ELSE
        WRITE(0,*) '  create_restart_file - ERROR - choice_forcing_method "', C%choice_forcing_method, '" not implemented in create_restart_file!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Leave definition mode:
    CALL handle_error(nf90_enddef( region%restart%ncid))
    
    ! Write the x, y, zeta and months variable data
    CALL handle_error( nf90_put_var( region%restart%ncid, region%restart%id_var_x,        region%grid%x                            ))
    CALL handle_error( nf90_put_var( region%restart%ncid, region%restart%id_var_y,        region%grid%y                            ))
    CALL handle_error( nf90_put_var( region%restart%ncid, region%restart%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( region%restart%ncid, region%restart%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( region%restart%ncid))
    
    ! Close the file
    CALL close_netcdf_file(region%restart%ncid)
    
  END SUBROUTINE create_restart_file
  SUBROUTINE create_help_fields_file( region)
    ! Create a new help fields file, containing secondary model output (not needed for a restart, but interesting to look at)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region

    ! Local variables:
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x, y, z, m, t
    
    IF (.NOT. par%master) RETURN
    
    ! Set time frame index to 1
    region%help_fields%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(region%help_fields%filename))
    IF (file_exists) THEN
      WRITE(0,*) '  create_help_fields_file - ERROR: ', TRIM(region%help_fields%filename), ' already exists!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( region%help_fields%filename)
    CALL handle_error(nf90_create(region%help_fields%filename,IOR(nf90_clobber,nf90_share),region%help_fields%ncid))
        
    ! Define dimensions:
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_x,     region%grid%nx, region%help_fields%id_dim_x    )
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_y,     region%grid%ny, region%help_fields%id_dim_y    )
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_zeta,  C%nZ,           region%help_fields%id_dim_zeta )
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_month, 12,             region%help_fields%id_dim_month)
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_time,  nf90_unlimited, region%help_fields%id_dim_time )
    
    ! Placeholders for the dimension ID's, for shorter code
    x = region%help_fields%id_dim_x
    y = region%help_fields%id_dim_y
    z = region%help_fields%id_dim_zeta
    m = region%help_fields%id_dim_month
    t = region%help_fields%id_dim_time
    
    ! Dimension variables: zeta, month, time
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_x,     [region%help_fields%id_dim_x    ], region%help_fields%id_var_x,     long_name='X-coordinate', units='m')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_y,     [region%help_fields%id_dim_y    ], region%help_fields%id_var_y,     long_name='Y-coordinate', units='m')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_zeta,  [region%help_fields%id_dim_zeta ], region%help_fields%id_var_zeta,  long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_month, [region%help_fields%id_dim_month], region%help_fields%id_var_month, long_name='Month', units='1-12'    )
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_time,  [region%help_fields%id_dim_time ], region%help_fields%id_var_time,  long_name='Time', units='years'   )
    
    ! Define data variables
    CALL create_help_field( region, region%help_fields%id_help_field_01, C%help_field_01)
    CALL create_help_field( region, region%help_fields%id_help_field_02, C%help_field_02)
    CALL create_help_field( region, region%help_fields%id_help_field_03, C%help_field_03)
    CALL create_help_field( region, region%help_fields%id_help_field_04, C%help_field_04)
    CALL create_help_field( region, region%help_fields%id_help_field_05, C%help_field_05)
    CALL create_help_field( region, region%help_fields%id_help_field_06, C%help_field_06)
    CALL create_help_field( region, region%help_fields%id_help_field_07, C%help_field_07)
    CALL create_help_field( region, region%help_fields%id_help_field_08, C%help_field_08)
    CALL create_help_field( region, region%help_fields%id_help_field_09, C%help_field_09)
    CALL create_help_field( region, region%help_fields%id_help_field_10, C%help_field_10)
    CALL create_help_field( region, region%help_fields%id_help_field_11, C%help_field_11)
    CALL create_help_field( region, region%help_fields%id_help_field_12, C%help_field_12)
    CALL create_help_field( region, region%help_fields%id_help_field_13, C%help_field_13)
    CALL create_help_field( region, region%help_fields%id_help_field_14, C%help_field_14)
    CALL create_help_field( region, region%help_fields%id_help_field_15, C%help_field_15)
    CALL create_help_field( region, region%help_fields%id_help_field_16, C%help_field_16)
    CALL create_help_field( region, region%help_fields%id_help_field_17, C%help_field_17)
    CALL create_help_field( region, region%help_fields%id_help_field_18, C%help_field_18)
    CALL create_help_field( region, region%help_fields%id_help_field_19, C%help_field_19)
    CALL create_help_field( region, region%help_fields%id_help_field_20, C%help_field_20)
    CALL create_help_field( region, region%help_fields%id_help_field_21, C%help_field_21)
    CALL create_help_field( region, region%help_fields%id_help_field_22, C%help_field_22)
    CALL create_help_field( region, region%help_fields%id_help_field_23, C%help_field_23)
    CALL create_help_field( region, region%help_fields%id_help_field_24, C%help_field_24)
    CALL create_help_field( region, region%help_fields%id_help_field_25, C%help_field_25)
    CALL create_help_field( region, region%help_fields%id_help_field_26, C%help_field_26)
    CALL create_help_field( region, region%help_fields%id_help_field_27, C%help_field_27)
    CALL create_help_field( region, region%help_fields%id_help_field_28, C%help_field_28)
    CALL create_help_field( region, region%help_fields%id_help_field_29, C%help_field_29)
    CALL create_help_field( region, region%help_fields%id_help_field_30, C%help_field_30)
    CALL create_help_field( region, region%help_fields%id_help_field_31, C%help_field_31)
    CALL create_help_field( region, region%help_fields%id_help_field_32, C%help_field_32)
    CALL create_help_field( region, region%help_fields%id_help_field_33, C%help_field_33)
    CALL create_help_field( region, region%help_fields%id_help_field_34, C%help_field_34)
    CALL create_help_field( region, region%help_fields%id_help_field_35, C%help_field_35)
    CALL create_help_field( region, region%help_fields%id_help_field_36, C%help_field_36)
    CALL create_help_field( region, region%help_fields%id_help_field_37, C%help_field_37)
    CALL create_help_field( region, region%help_fields%id_help_field_38, C%help_field_38)
    CALL create_help_field( region, region%help_fields%id_help_field_39, C%help_field_39)
    CALL create_help_field( region, region%help_fields%id_help_field_40, C%help_field_40)
    CALL create_help_field( region, region%help_fields%id_help_field_41, C%help_field_41)
    CALL create_help_field( region, region%help_fields%id_help_field_42, C%help_field_42)
    CALL create_help_field( region, region%help_fields%id_help_field_43, C%help_field_43)
    CALL create_help_field( region, region%help_fields%id_help_field_44, C%help_field_44)
    CALL create_help_field( region, region%help_fields%id_help_field_45, C%help_field_45)
    CALL create_help_field( region, region%help_fields%id_help_field_46, C%help_field_46)
    CALL create_help_field( region, region%help_fields%id_help_field_47, C%help_field_47)
    CALL create_help_field( region, region%help_fields%id_help_field_48, C%help_field_48)
    CALL create_help_field( region, region%help_fields%id_help_field_49, C%help_field_49)
    CALL create_help_field( region, region%help_fields%id_help_field_50, C%help_field_50)
    
    ! Leave definition mode:
    CALL handle_error(nf90_enddef( region%help_fields%ncid))
    
    ! Write the x, y, zeta, months, and lat/lon variable data
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_x,        region%grid%x                            ))
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_y,        region%grid%y                            ))
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( region%help_fields%ncid))
    
    ! Close the file
    CALL close_netcdf_file(region%help_fields%ncid)
    
  END SUBROUTINE create_help_fields_file
  SUBROUTINE create_help_field( region, id_var, field_name)
    ! Add a data field to the help_fields file
    ! (NOTE: if you want to add an extra option, be sure to also add it to write_help_field!)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    INTEGER,                        INTENT(INOUT) :: id_var
    CHARACTER(LEN=*),               INTENT(IN)    :: field_name
    
    ! Local variables:
    INTEGER                                       :: x, y, z, m, t
    
    IF (.NOT. par%master) RETURN
    
    ! Placeholders for the dimension ID's, for shorter code
    x = region%help_fields%id_dim_x
    y = region%help_fields%id_dim_y
    z = region%help_fields%id_dim_zeta
    m = region%help_fields%id_dim_month
    t = region%help_fields%id_dim_time
    
    IF (field_name == 'none') THEN
      RETURN
      
    ! Fields with no time dimension
    ! =============================
      
    ! Lat/lon
    ELSEIF (field_name == 'lat') THEN
      CALL create_double_var( region%help_fields%ncid, 'lat',                      [x, y      ], id_var, long_name='Latitude',  units='degrees north')
    ELSEIF (field_name == 'lon') THEN
      CALL create_double_var( region%help_fields%ncid, 'lon',                      [x, y      ], id_var, long_name='Longitude', units='degrees east')
      
    ! Geothermal heat flux
    ELSEIF (field_name == 'GHF') THEN
      CALL create_double_var( region%help_fields%ncid, 'GHF',                      [x, y      ], id_var, long_name='Geothermal heat flux', units='J m^-2 yr^-1')
 
    ! Forcing climates
    ELSEIF (field_name == 'GCM_Warm_T2m') THEN
      CALL create_double_var( region%help_fields%ncid, 'Warm_T2m',                 [x, y, m], id_var, long_name='Warm monthly mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'GCM_Warm_Precip') THEN
      CALL create_double_var( region%help_fields%ncid, 'Warm_Precip',              [x, y, m], id_var, long_name='Warm monthly total precipitation', units='mm')
    ELSEIF (field_name == 'GCM_Cold_T2m') THEN
      CALL create_double_var( region%help_fields%ncid, 'Cold_T2m',                 [x, y, m], id_var, long_name='Cold monthly mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'GCM_Cold_Precip') THEN
      CALL create_double_var( region%help_fields%ncid, 'Cold_Precip',              [x, y, m], id_var, long_name='Cold monthly total precipitation', units='mm')
    ELSEIF (field_name == 'GCM_PI_T2m') THEN
      CALL create_double_var( region%help_fields%ncid, 'Ref_PI_T2m',               [x, y, m], id_var, long_name='Ref PI monthly mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'GCM_PI_Precip') THEN
      CALL create_double_var( region%help_fields%ncid, 'Ref_PI_Precip',            [x, y, m], id_var, long_name='Ref PI monthly total precipitation', units='mm')
    ELSEIF (field_name == 'PD_obs_T2m') THEN
      CALL create_double_var( region%help_fields%ncid, 'Base_PD_T2m',              [x, y, m], id_var, long_name='Base PD monthly mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'PD_obs_Precip') THEN
      CALL create_double_var( region%help_fields%ncid, 'Base_PD_Precip',           [x, y, m], id_var, long_name='Base PD monthly total precipitation', units='mm')
     
    ! Fields with a time dimension
    ! ============================
      
    ! Geometry
    ELSEIF (field_name == 'Hi') THEN
      CALL create_double_var( region%help_fields%ncid, 'Hi',                       [x, y,    t], id_var, long_name='Ice thickness', units='m')
    ELSEIF (field_name == 'Hb') THEN
      CALL create_double_var( region%help_fields%ncid, 'Hb',                       [x, y,    t], id_var, long_name='Bedrock elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'Hs') THEN
      CALL create_double_var( region%help_fields%ncid, 'Hs',                       [x, y,    t], id_var, long_name='Surface elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'SL') THEN
      CALL create_double_var( region%help_fields%ncid, 'SL',                       [x, y,    t], id_var, long_name='Geoid elevation', units='m w.r.t PD sealevel')
    ELSEIF (field_name == 'dHs_dx') THEN
      CALL create_double_var( region%help_fields%ncid, 'dHs_dx',                   [x, y,    t], id_var, long_name='Surface slope in x-direction', units='m/m')
    ELSEIF (field_name == 'dHs_dy') THEN
      CALL create_double_var( region%help_fields%ncid, 'dHs_dy',                   [x, y,    t], id_var, long_name='Surface slope in y-direction', units='m/m')
      
    ! Thermal properties
    ELSEIF (field_name == 'Ti') THEN
      CALL create_double_var( region%help_fields%ncid, 'Ti',                       [x, y, z, t], id_var, long_name='Englacial temperature', units='K')
    ELSEIF (field_name == 'Cpi') THEN
      CALL create_double_var( region%help_fields%ncid, 'Cpi',                      [x, y, z, t], id_var, long_name='Ice heat capacity', units='J kg^-1 K^-1')
    ELSEIF (field_name == 'Ki') THEN
      CALL create_double_var( region%help_fields%ncid, 'Ki',                       [x, y, z, t], id_var, long_name='Ice thermal conductivity', units='J m^-1 K^-1 yr^-1')
    ELSEIF (field_name == 'Ti_basal') THEN
      CALL create_double_var( region%help_fields%ncid, 'Ti_basal',                 [x, y,    t], id_var, long_name='Ice basal temperature', units='K')
    ELSEIF (field_name == 'Ti_pmp') THEN
      CALL create_double_var( region%help_fields%ncid, 'Ti_pmp',                   [x, y, z, t], id_var, long_name='Ice pressure melting point temperature', units='K')
    ELSEIF (field_name == 'A_flow') THEN
      CALL create_double_var( region%help_fields%ncid, 'A_flow_3D',                [x, y, z, t], id_var, long_name='Ice flow factor', units='Pa^-3 y^-1')
    ELSEIF (field_name == 'A_flow_mean') THEN
      CALL create_double_var( region%help_fields%ncid, 'A_flow_vav',               [x, y,    t], id_var, long_name='Vertically averaged ice flow factor', units='Pa^-3 y^-1')
      
    ! Velocity fields
    ELSEIF (field_name == 'u_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'u_3D',                     [x, y, z, t], id_var, long_name='3D ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'v_3D',                     [x, y, z, t], id_var, long_name='3D ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'w_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'w_3D',                     [x, y, z, t], id_var, long_name='3D ice z-velocity', units='m/yr')
    ELSEIF (field_name == 'u_vav') THEN
      CALL create_double_var( region%help_fields%ncid, 'u_vav',                    [x, y,    t], id_var, long_name='Vertically averaged ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_vav') THEN
      CALL create_double_var( region%help_fields%ncid, 'v_vav',                    [x, y,    t], id_var, long_name='Vertically averaged ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_vav') THEN
      CALL create_double_var( region%help_fields%ncid, 'uabs_vav',                 [x, y,    t], id_var, long_name='Vertically averaged ice velocity', units='m/yr')
    ELSEIF (field_name == 'u_surf') THEN
      CALL create_double_var( region%help_fields%ncid, 'u_surf',                   [x, y,    t], id_var, long_name='Surface ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_surf') THEN
      CALL create_double_var( region%help_fields%ncid, 'v_surf',                   [x, y,    t], id_var, long_name='Surface ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_surf') THEN
      CALL create_double_var( region%help_fields%ncid, 'uabs_surf',                [x, y,    t], id_var, long_name='Surface ice velocity', units='m/yr')
    ELSEIF (field_name == 'u_base') THEN
      CALL create_double_var( region%help_fields%ncid, 'u_base',                   [x, y,    t], id_var, long_name='Basal ice x-velocity', units='m/yr')
    ELSEIF (field_name == 'v_base') THEN
      CALL create_double_var( region%help_fields%ncid, 'v_base',                   [x, y,    t], id_var, long_name='Basal ice y-velocity', units='m/yr')
    ELSEIF (field_name == 'uabs_base') THEN
      CALL create_double_var( region%help_fields%ncid, 'uabs_base',                [x, y,    t], id_var, long_name='Basal ice velocity', units='m/yr')
      
    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL create_double_var( region%help_fields%ncid, 'T2m',                      [x, y, m, t], id_var, long_name='Monthly mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'T2m_year') THEN
      CALL create_double_var( region%help_fields%ncid, 'T2m_year',                 [x, y,    t], id_var, long_name='Annual mean 2-m air temperature', units='K')
    ELSEIF (field_name == 'Precip') THEN
      CALL create_double_var( region%help_fields%ncid, 'Precip',                   [x, y, m, t], id_var, long_name='Monthly total precipitation', units='mm')
    ELSEIF (field_name == 'Precip_year') THEN
      CALL create_double_var( region%help_fields%ncid, 'Precip_year',              [x, y,    t], id_var, long_name='Annual total precipitation', units='mm')
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL create_double_var( region%help_fields%ncid, 'Wind_WE',                  [x, y, m, t], id_var, long_name='Monthly mean zonal wind', units='m/s')
    ELSEIF (field_name == 'Wind_WE_year') THEN
      CALL create_double_var( region%help_fields%ncid, 'Wind_WE_year',             [x, y,    t], id_var, long_name='Annual mean zonal wind', units='m/s')
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL create_double_var( region%help_fields%ncid, 'Wind_SN',                  [x, y, m, t], id_var, long_name='Monthly mean meridional wind', units='m/s')
    ELSEIF (field_name == 'Wind_SN_year') THEN
      CALL create_double_var( region%help_fields%ncid, 'Wind_SN_year',             [x, y,    t], id_var, long_name='Annual mean meridional wind', units='m/s')
      
    ! Mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL create_double_var( region%help_fields%ncid, 'SMB',                      [x, y, m, t], id_var, long_name='Monthly surface mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'SMB_year') THEN
      CALL create_double_var( region%help_fields%ncid, 'SMB_year',                 [x, y,    t], id_var, long_name='Annual surface mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL create_double_var( region%help_fields%ncid, 'BMB_sheet',                [x, y,    t], id_var, long_name='Annual basal mass balance for grounded ice', units='m ice equivalent')
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL create_double_var( region%help_fields%ncid, 'BMB_shelf',                [x, y,    t], id_var, long_name='Annual basal mass balance for floating ice', units='m ice equivalent')
    ELSEIF (field_name == 'BMB') THEN
      CALL create_double_var( region%help_fields%ncid, 'BMB',                      [x, y,    t], id_var, long_name='Annual basal mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'Snowfall') THEN
      CALL create_double_var( region%help_fields%ncid, 'Snowfall',                 [x, y, m, t], id_var, long_name='Monthly total snowfall', units='m water equivalent')
    ELSEIF (field_name == 'Snowfall_year') THEN
      CALL create_double_var( region%help_fields%ncid, 'Snowfall_year',            [x, y,    t], id_var, long_name='Annual total snowfall', units='m water equivalent')
    ELSEIF (field_name == 'Rainfall') THEN
      CALL create_double_var( region%help_fields%ncid, 'Rainfall',                 [x, y, m, t], id_var, long_name='Monthly total rainfall', units='m water equivalent')
    ELSEIF (field_name == 'Rainfall_year') THEN
      CALL create_double_var( region%help_fields%ncid, 'Rainfall_year',            [x, y,    t], id_var, long_name='Annual total rainfall', units='m water equivalent')
    ELSEIF (field_name == 'AddedFirn') THEN
      CALL create_double_var( region%help_fields%ncid, 'AddedFirn',                [x, y, m, t], id_var, long_name='Monthly total added firn', units='m water equivalent')
    ELSEIF (field_name == 'AddedFirn_year') THEN
      CALL create_double_var( region%help_fields%ncid, 'AddedFirn_year',           [x, y,    t], id_var, long_name='Annual total added firn', units='m water equivalent')
    ELSEIF (field_name == 'Refreezing') THEN
      CALL create_double_var( region%help_fields%ncid, 'Refreezing',               [x, y, m, t], id_var, long_name='Monthly total refreezing', units='m water equivalent')
    ELSEIF (field_name == 'Refreezing_year') THEN
      CALL create_double_var( region%help_fields%ncid, 'Refreezing_year',          [x, y,    t], id_var, long_name='Annual total refreezing', units='m water equivalent')
    ELSEIF (field_name == 'Runoff') THEN
      CALL create_double_var( region%help_fields%ncid, 'Runoff',                   [x, y, m, t], id_var, long_name='Monthly total runoff', units='m water equivalent')
    ELSEIF (field_name == 'Runoff_year') THEN
      CALL create_double_var( region%help_fields%ncid, 'Runoff_year',              [x, y,    t], id_var, long_name='Annual total runoff', units='m water equivalent')
    ELSEIF (field_name == 'Albedo') THEN
      CALL create_double_var( region%help_fields%ncid, 'Albedo',                   [x, y, m, t], id_var, long_name='Monthly mean albedo')
    ELSEIF (field_name == 'Albedo_year') THEN
      CALL create_double_var( region%help_fields%ncid, 'Albedo_year',              [x, y,    t], id_var, long_name='Annual mean albedo')
    ELSEIF (field_name == 'FirnDepth') THEN
      CALL create_double_var( region%help_fields%ncid, 'FirnDepth',                [x, y, m, t], id_var, long_name='Monthly mean firn layer depth', units='m water equivalent')
    ELSEIF (field_name == 'FirnDepth_year') THEN
      CALL create_double_var( region%help_fields%ncid, 'FirnDepth_year',           [x, y,    t], id_var, long_name='Annual mean firn layer depth', units='m water equivalent')
      
    ! Masks
    ELSEIF (field_name == 'mask') THEN
      CALL create_int_var(    region%help_fields%ncid, 'mask',                     [x, y,    t], id_var, long_name='mask')
    ELSEIF (field_name == 'mask_land') THEN
      CALL create_int_var(    region%help_fields%ncid, 'mask_land',                [x, y,    t], id_var, long_name='land mask')
    ELSEIF (field_name == 'mask_ocean') THEN
      CALL create_int_var(    region%help_fields%ncid, 'mask_ocean',               [x, y,    t], id_var, long_name='ocean mask')
    ELSEIF (field_name == 'mask_lake') THEN
      CALL create_int_var(    region%help_fields%ncid, 'mask_lake',                [x, y,    t], id_var, long_name='lake mask')
    ELSEIF (field_name == 'mask_ice') THEN
      CALL create_int_var(    region%help_fields%ncid, 'mask_ice',                 [x, y,    t], id_var, long_name='ice mask')
    ELSEIF (field_name == 'mask_sheet') THEN
      CALL create_int_var(    region%help_fields%ncid, 'mask_sheet',               [x, y,    t], id_var, long_name='sheet mask')
    ELSEIF (field_name == 'mask_shelf') THEN
      CALL create_int_var(    region%help_fields%ncid, 'mask_shelf',               [x, y,    t], id_var, long_name='shelf mask')
    ELSEIF (field_name == 'mask_coast') THEN
      CALL create_int_var(    region%help_fields%ncid, 'mask_coast',               [x, y,    t], id_var, long_name='coast mask')
    ELSEIF (field_name == 'mask_margin') THEN
      CALL create_int_var(    region%help_fields%ncid, 'mask_margin',              [x, y,    t], id_var, long_name='margin mask')
    ELSEIF (field_name == 'mask_gl') THEN
      CALL create_int_var(    region%help_fields%ncid, 'mask_gl',                  [x, y,    t], id_var, long_name='grounding-line mask')
    ELSEIF (field_name == 'mask_cf') THEN
      CALL create_int_var(    region%help_fields%ncid, 'mask_cf',                  [x, y,    t], id_var, long_name='calving-front mask')
      
    ! Basal conditions
    ELSEIF (field_name == 'A_slid') THEN
      CALL create_double_var( region%help_fields%ncid, 'A_slid',                   [x, y,    t], id_var, long_name='Basal sliding coefficient', units='Pa^-m yr^-1')
    ELSEIF (field_name == 'phi_fric') THEN
      CALL create_double_var( region%help_fields%ncid, 'phi_fric',                 [x, y,    t], id_var, long_name='till friction angle', units='degrees')
    ELSEIF (field_name == 'tau_yield') THEN
      CALL create_double_var( region%help_fields%ncid, 'tau_yield',                [x, y,    t], id_var, long_name='basal yield stress', units='Pa')
      
    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL create_double_var( region%help_fields%ncid, 'iso_ice',                  [x, y,    t], id_var, long_name='Vertically averaged ice d18O', units='per mille')
    ELSEIF (field_name == 'iso_surf') THEN
      CALL create_double_var( region%help_fields%ncid, 'iso_surf',                 [x, y,    t], id_var, long_name='d18O of precipitation', units='per mille')
    
    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      CALL create_double_var( region%help_fields%ncid, 'dHb',                      [x, y,    t], id_var, long_name='Change in bedrock elevation w.r.t. PD', units='m')
      
      
    ELSE
      WRITE(0,*) ' ERROR: help field "', TRIM(field_name), '" not implemented in create_help_field!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE create_help_field
  SUBROUTINE create_debug_file( region)
    ! Create the debug NetCDF file; a lot of data fields but no time dimension.
    
    USE data_types_netcdf_module, ONLY: type_netcdf_debug
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region

    ! Local variables:
    TYPE(type_netcdf_debug)                       :: debug_temp
    CHARACTER(LEN=20)                             :: short_filename
    INTEGER                                       :: n
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x, y, z, m, ncid
    
    IF (.NOT. par%master) RETURN
    
    IF (.NOT. C%do_write_debug_data) RETURN

    ! Determine debug NetCDF filename for this model region
    short_filename = 'debug_NAM.nc'
    short_filename(7:9) = region%name
    DO n = 1, 256
      debug_temp%filename(n:n) = ' '
    END DO
    debug_temp%filename = TRIM(C%output_dir)//TRIM(short_filename)

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(debug_temp%filename))
    IF (file_exists) THEN
      WRITE(0,*) '  create_debug_file - ERROR: ', TRIM(debug_temp%filename), ' already exists!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( debug_temp%filename)
    CALL handle_error(nf90_create(debug_temp%filename,IOR(nf90_clobber,nf90_share),ncid))
    
    debug_temp%ncid = ncid
        
    ! Define dimensions:
    CALL create_dim( ncid, debug_temp%name_dim_x,         region%grid%nx,     debug_temp%id_dim_x         )
    CALL create_dim( ncid, debug_temp%name_dim_y,         region%grid%ny,     debug_temp%id_dim_y         )
    CALL create_dim( ncid, debug_temp%name_dim_zeta,      C%nZ,               debug_temp%id_dim_zeta      ) ! Scaled vertical coordinate
    CALL create_dim( ncid, debug_temp%name_dim_month,     12,                 debug_temp%id_dim_month     ) ! Months (for monthly data)
    
    ! Placeholders for the dimension ID's, for shorter code
    x = debug_temp%id_dim_x
    y = debug_temp%id_dim_y
    z = debug_temp%id_dim_zeta
    m = debug_temp%id_dim_month
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.
    
    ! Dimension variables: zeta, month, time
    CALL create_double_var( ncid, debug_temp%name_var_x,                [x         ], debug_temp%id_var_x,                long_name='X-coordinate', units='m')
    CALL create_double_var( ncid, debug_temp%name_var_y,                [   y      ], debug_temp%id_var_y,                long_name='Y-coordinate', units='m')
    CALL create_double_var( ncid, debug_temp%name_var_zeta,             [      z   ], debug_temp%id_var_zeta,             long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( ncid, debug_temp%name_var_month,            [         m], debug_temp%id_var_month,            long_name='Month', units='1-12'    )
    
    ! Dummy variables (useful for debugging)
    CALL create_int_var(    ncid, debug_temp%name_var_int_2D_01,        [x, y      ], debug_temp%id_var_int_2D_01,        long_name='2D int variable 01')
    CALL create_int_var(    ncid, debug_temp%name_var_int_2D_02,        [x, y      ], debug_temp%id_var_int_2D_02,        long_name='2D int variable 02')
    CALL create_int_var(    ncid, debug_temp%name_var_int_2D_03,        [x, y      ], debug_temp%id_var_int_2D_03,        long_name='2D int variable 03')
    CALL create_int_var(    ncid, debug_temp%name_var_int_2D_04,        [x, y      ], debug_temp%id_var_int_2D_04,        long_name='2D int variable 04')
    CALL create_int_var(    ncid, debug_temp%name_var_int_2D_05,        [x, y      ], debug_temp%id_var_int_2D_05,        long_name='2D int variable 05')
    CALL create_int_var(    ncid, debug_temp%name_var_int_2D_06,        [x, y      ], debug_temp%id_var_int_2D_06,        long_name='2D int variable 06')
    CALL create_int_var(    ncid, debug_temp%name_var_int_2D_07,        [x, y      ], debug_temp%id_var_int_2D_07,        long_name='2D int variable 07')
    CALL create_int_var(    ncid, debug_temp%name_var_int_2D_08,        [x, y      ], debug_temp%id_var_int_2D_08,        long_name='2D int variable 08')
    CALL create_int_var(    ncid, debug_temp%name_var_int_2D_09,        [x, y      ], debug_temp%id_var_int_2D_09,        long_name='2D int variable 09')
    CALL create_int_var(    ncid, debug_temp%name_var_int_2D_10,        [x, y      ], debug_temp%id_var_int_2D_10,        long_name='2D int variable 10')
    
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_01,         [x, y      ], debug_temp%id_var_dp_2D_01,         long_name='2D dp variable 01')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_02,         [x, y      ], debug_temp%id_var_dp_2D_02,         long_name='2D dp variable 02')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_03,         [x, y      ], debug_temp%id_var_dp_2D_03,         long_name='2D dp variable 03')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_04,         [x, y      ], debug_temp%id_var_dp_2D_04,         long_name='2D dp variable 04')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_05,         [x, y      ], debug_temp%id_var_dp_2D_05,         long_name='2D dp variable 05')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_06,         [x, y      ], debug_temp%id_var_dp_2D_06,         long_name='2D dp variable 06')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_07,         [x, y      ], debug_temp%id_var_dp_2D_07,         long_name='2D dp variable 07')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_08,         [x, y      ], debug_temp%id_var_dp_2D_08,         long_name='2D dp variable 08')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_09,         [x, y      ], debug_temp%id_var_dp_2D_09,         long_name='2D dp variable 09')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_10,         [x, y      ], debug_temp%id_var_dp_2D_10,         long_name='2D dp variable 10')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_11,         [x, y      ], debug_temp%id_var_dp_2D_11,         long_name='2D dp variable 11')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_12,         [x, y      ], debug_temp%id_var_dp_2D_12,         long_name='2D dp variable 12')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_13,         [x, y      ], debug_temp%id_var_dp_2D_13,         long_name='2D dp variable 13')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_14,         [x, y      ], debug_temp%id_var_dp_2D_14,         long_name='2D dp variable 14')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_15,         [x, y      ], debug_temp%id_var_dp_2D_15,         long_name='2D dp variable 15')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_16,         [x, y      ], debug_temp%id_var_dp_2D_16,         long_name='2D dp variable 16')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_17,         [x, y      ], debug_temp%id_var_dp_2D_17,         long_name='2D dp variable 17')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_18,         [x, y      ], debug_temp%id_var_dp_2D_18,         long_name='2D dp variable 18')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_19,         [x, y      ], debug_temp%id_var_dp_2D_19,         long_name='2D dp variable 19')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_20,         [x, y      ], debug_temp%id_var_dp_2D_20,         long_name='2D dp variable 20')
    
    CALL create_double_var( ncid, debug_temp%name_var_dp_3D_01,         [x, y, z   ], debug_temp%id_var_dp_3D_01,         long_name='3D dp variable 01')
    CALL create_double_var( ncid, debug_temp%name_var_dp_3D_02,         [x, y, z   ], debug_temp%id_var_dp_3D_02,         long_name='3D dp variable 02')
    CALL create_double_var( ncid, debug_temp%name_var_dp_3D_03,         [x, y, z   ], debug_temp%id_var_dp_3D_03,         long_name='3D dp variable 03')
    CALL create_double_var( ncid, debug_temp%name_var_dp_3D_04,         [x, y, z   ], debug_temp%id_var_dp_3D_04,         long_name='3D dp variable 04')
    CALL create_double_var( ncid, debug_temp%name_var_dp_3D_05,         [x, y, z   ], debug_temp%id_var_dp_3D_05,         long_name='3D dp variable 05')
    CALL create_double_var( ncid, debug_temp%name_var_dp_3D_06,         [x, y, z   ], debug_temp%id_var_dp_3D_06,         long_name='3D dp variable 06')
    CALL create_double_var( ncid, debug_temp%name_var_dp_3D_07,         [x, y, z   ], debug_temp%id_var_dp_3D_07,         long_name='3D dp variable 07')
    CALL create_double_var( ncid, debug_temp%name_var_dp_3D_08,         [x, y, z   ], debug_temp%id_var_dp_3D_08,         long_name='3D dp variable 08')
    CALL create_double_var( ncid, debug_temp%name_var_dp_3D_09,         [x, y, z   ], debug_temp%id_var_dp_3D_09,         long_name='3D dp variable 09')
    CALL create_double_var( ncid, debug_temp%name_var_dp_3D_10,         [x, y, z   ], debug_temp%id_var_dp_3D_10,         long_name='3D dp variable 10')
    
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_monthly_01, [x, y,    m], debug_temp%id_var_dp_2D_monthly_01, long_name='2D_monthly dp variable 01')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_monthly_02, [x, y,    m], debug_temp%id_var_dp_2D_monthly_02, long_name='2D_monthly dp variable 02')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_monthly_03, [x, y,    m], debug_temp%id_var_dp_2D_monthly_03, long_name='2D_monthly dp variable 03')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_monthly_04, [x, y,    m], debug_temp%id_var_dp_2D_monthly_04, long_name='2D_monthly dp variable 04')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_monthly_05, [x, y,    m], debug_temp%id_var_dp_2D_monthly_05, long_name='2D_monthly dp variable 05')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_monthly_06, [x, y,    m], debug_temp%id_var_dp_2D_monthly_06, long_name='2D_monthly dp variable 06')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_monthly_07, [x, y,    m], debug_temp%id_var_dp_2D_monthly_07, long_name='2D_monthly dp variable 07')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_monthly_08, [x, y,    m], debug_temp%id_var_dp_2D_monthly_08, long_name='2D_monthly dp variable 08')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_monthly_09, [x, y,    m], debug_temp%id_var_dp_2D_monthly_09, long_name='2D_monthly dp variable 09')
    CALL create_double_var( ncid, debug_temp%name_var_dp_2D_monthly_10, [x, y,    m], debug_temp%id_var_dp_2D_monthly_10, long_name='2D_monthly dp variable 10')
    
    ! Leave definition mode:
    CALL handle_error(nf90_enddef( ncid))
    
    ! Write the x, y, zeta, months, and lat/lon variable data
    CALL handle_error( nf90_put_var( ncid, debug_temp%id_var_x,        region%grid%x                            ))
    CALL handle_error( nf90_put_var( ncid, debug_temp%id_var_y,        region%grid%y                            ))
    CALL handle_error( nf90_put_var( ncid, debug_temp%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( ncid))
    
    ! Close the file
    CALL close_netcdf_file(ncid)
    
    ! Copy data to the relevant debug structure
    IF     (region%name == 'NAM') THEN
      debug_NAM%netcdf = debug_temp
    ELSEIF (region%name == 'EAS') THEN
      debug_EAS%netcdf = debug_temp
    ELSEIF (region%name == 'GRL') THEN
      debug_GRL%netcdf = debug_temp
    ELSEIF (region%name == 'ANT') THEN
      debug_ANT%netcdf = debug_temp
    END IF
    
  END SUBROUTINE create_debug_file
  SUBROUTINE create_SELEN_output_file( SELEN)
    ! Create a new NetCDF output file for SELEN (on the irregular global mesh)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN

    ! Local variables:
    LOGICAL                                       :: file_exists
    INTEGER                                       :: vi, ti, ci, three, time
    
    IF (.NOT. par%master) RETURN
    
    ! Set time frame index to 1
    SELEN%output%ti = 1

    ! Set output filename
    SELEN%output%filename = TRIM(C%output_dir)//'SELEN_output.nc'

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(SELEN%output%filename))
    IF(file_exists) THEN
      WRITE(0,*) 'ERROR: ', TRIM(SELEN%output%filename), ' already exists!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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
    
  END SUBROUTINE create_SELEN_output_file
  SUBROUTINE create_global_scalar_output_file( netcdf)
    ! Create a new global scalar output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_netcdf_scalars_global), INTENT(INOUT) :: netcdf

    ! Local variables:
    LOGICAL                                         :: file_exists
    INTEGER                                         :: t
    
    IF (.NOT. par%master) RETURN
    
    ! Set time frame index to 1
    netcdf%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    netcdf%filename = TRIM(C%output_dir)//'/scalar_output_global.nc'
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      WRITE(0,*) '  create_global_scalar_output_file - ERROR: ', TRIM(netcdf%filename), ' already exists!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( netcdf%filename)
    CALL handle_error(nf90_create(netcdf%filename,IOR(nf90_clobber,nf90_share),netcdf%ncid))
        
    ! Define dimensions:
    CALL create_dim( netcdf%ncid, netcdf%name_dim_time,      nf90_unlimited,     netcdf%id_dim_time      ) ! Time frames
    
    ! Placeholders for the dimension ID's, for shorter code
    t = netcdf%id_dim_time
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.
    
    ! Dimension variables: zeta, month, time
    CALL create_double_var( netcdf%ncid, netcdf%name_var_time,          [t], netcdf%id_var_time,          long_name='Time', units='years'   )
    
    ! Sea level
    CALL create_double_var( netcdf%ncid, netcdf%name_var_GMSL,          [t], netcdf%id_var_GMSL,          long_name='Global mean sea level change', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_GMSL_NAM,      [t], netcdf%id_var_GMSL_NAM,      long_name='Global mean sea level change from ice in North America', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_GMSL_EAS,      [t], netcdf%id_var_GMSL_EAS,      long_name='Global mean sea level change from ice in Eurasia',       units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_GMSL_GRL,      [t], netcdf%id_var_GMSL_GRL,      long_name='Global mean sea level change from ice in Greenland',       units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_GMSL_ANT,      [t], netcdf%id_var_GMSL_ANT,      long_name='Global mean sea level change from ice in Antarctica',       units='m')
    
    ! CO2
    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_CO2_obs,       [t], netcdf%id_var_CO2_obs,       long_name='Observed atmospheric CO2 concentration', units='ppm')
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_CO2_obs,       [t], netcdf%id_var_CO2_obs,       long_name='Observed atmospheric CO2 concentration', units='ppm')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_CO2_mod,       [t], netcdf%id_var_CO2_mod,       long_name='Modelled atmospheric CO2 concentration', units='ppm')
    ELSEIF (C%choice_forcing_method == 'climate_direct') THEN
    ELSEIF (C%choice_forcing_method == 'SMB_direct') THEN
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in create_global_scalar_output_file!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! d18O
    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_mod,      [t], netcdf%id_var_d18O_mod,      long_name='Modelled benthic d18O', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_ice,      [t], netcdf%id_var_d18O_ice,      long_name='Modelled benthic d18O from global ice volume', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_Tdw,      [t], netcdf%id_var_d18O_Tdw,      long_name='Modelled benthic d18O from deep-water temperature', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_NAM,      [t], netcdf%id_var_d18O_NAM,      long_name='Modelled benthic d18O from ice in North America', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_EAS,      [t], netcdf%id_var_d18O_EAS,      long_name='Modelled benthic d18O from ice in Eurasia', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_GRL,      [t], netcdf%id_var_d18O_GRL,      long_name='Modelled benthic d18O from ice in Greenland', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_ANT,      [t], netcdf%id_var_d18O_ANT,      long_name='Modelled benthic d18O from ice in Antarctica', units='per mil')
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_obs,      [t], netcdf%id_var_d18O_obs,      long_name='Observed benthic d18O', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_mod,      [t], netcdf%id_var_d18O_mod,      long_name='Modelled benthic d18O', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_ice,      [t], netcdf%id_var_d18O_ice,      long_name='Modelled benthic d18O from global ice volume', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_Tdw,      [t], netcdf%id_var_d18O_Tdw,      long_name='Modelled benthic d18O from deep-water temperature', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_NAM,      [t], netcdf%id_var_d18O_NAM,      long_name='Modelled benthic d18O from ice in North America', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_EAS,      [t], netcdf%id_var_d18O_EAS,      long_name='Modelled benthic d18O from ice in Eurasia', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_GRL,      [t], netcdf%id_var_d18O_GRL,      long_name='Modelled benthic d18O from ice in Greenland', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_ANT,      [t], netcdf%id_var_d18O_ANT,      long_name='Modelled benthic d18O from ice in Antarctica', units='per mil')
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_obs,      [t], netcdf%id_var_d18O_obs,      long_name='Observed benthic d18O', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_mod,      [t], netcdf%id_var_d18O_mod,      long_name='Modelled benthic d18O', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_ice,      [t], netcdf%id_var_d18O_ice,      long_name='Modelled benthic d18O from global ice volume', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_Tdw,      [t], netcdf%id_var_d18O_Tdw,      long_name='Modelled benthic d18O from deep-water temperature', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_NAM,      [t], netcdf%id_var_d18O_NAM,      long_name='Modelled benthic d18O from ice in North America', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_EAS,      [t], netcdf%id_var_d18O_EAS,      long_name='Modelled benthic d18O from ice in Eurasia', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_GRL,      [t], netcdf%id_var_d18O_GRL,      long_name='Modelled benthic d18O from ice in Greenland', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_ANT,      [t], netcdf%id_var_d18O_ANT,      long_name='Modelled benthic d18O from ice in Antarctica', units='per mil')
    ELSEIF (C%choice_forcing_method == 'climate_direct') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_mod,      [t], netcdf%id_var_d18O_mod,      long_name='Modelled benthic d18O', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_ice,      [t], netcdf%id_var_d18O_ice,      long_name='Modelled benthic d18O from global ice volume', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_Tdw,      [t], netcdf%id_var_d18O_Tdw,      long_name='Modelled benthic d18O from deep-water temperature', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_NAM,      [t], netcdf%id_var_d18O_NAM,      long_name='Modelled benthic d18O from ice in North America', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_EAS,      [t], netcdf%id_var_d18O_EAS,      long_name='Modelled benthic d18O from ice in Eurasia', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_GRL,      [t], netcdf%id_var_d18O_GRL,      long_name='Modelled benthic d18O from ice in Greenland', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_ANT,      [t], netcdf%id_var_d18O_ANT,      long_name='Modelled benthic d18O from ice in Antarctica', units='per mil')
    ELSEIF (C%choice_forcing_method == 'SMB_direct') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_mod,      [t], netcdf%id_var_d18O_mod,      long_name='Modelled benthic d18O', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_ice,      [t], netcdf%id_var_d18O_ice,      long_name='Modelled benthic d18O from global ice volume', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_Tdw,      [t], netcdf%id_var_d18O_Tdw,      long_name='Modelled benthic d18O from deep-water temperature', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_NAM,      [t], netcdf%id_var_d18O_NAM,      long_name='Modelled benthic d18O from ice in North America', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_EAS,      [t], netcdf%id_var_d18O_EAS,      long_name='Modelled benthic d18O from ice in Eurasia', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_GRL,      [t], netcdf%id_var_d18O_GRL,      long_name='Modelled benthic d18O from ice in Greenland', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_ANT,      [t], netcdf%id_var_d18O_ANT,      long_name='Modelled benthic d18O from ice in Antarctica', units='per mil')
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in create_global_scalar_output_file!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Temperature (surface and deep-water)
    CALL create_double_var( netcdf%ncid, netcdf%name_var_dT_glob,       [t], netcdf%id_var_dT_glob,       long_name='Global annual mean surface temperature change', units='K')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_dT_dw,         [t], netcdf%id_var_dT_dw,         long_name='Deep-water temperature change', units='K')
    
    ! Computation time for different model components
    CALL create_double_var( netcdf%ncid, netcdf%name_var_tcomp_total,   [t], netcdf%id_var_tcomp_total,   long_name='Total computation time', units='s')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_tcomp_ice,     [t], netcdf%id_var_tcomp_ice,     long_name='Total computation time for ice dynamics', units='s')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_tcomp_thermo,  [t], netcdf%id_var_tcomp_thermo,  long_name='Total computation time for thermodynamics', units='s')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_tcomp_climate, [t], netcdf%id_var_tcomp_climate, long_name='Total computation time for climate+SMB+BMB', units='s')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_tcomp_GIA,     [t], netcdf%id_var_tcomp_GIA,     long_name='Total computation time for GIA', units='s')
    
    ! Leave definition mode:
    CALL handle_error(nf90_enddef( netcdf%ncid))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( netcdf%ncid))
    
    ! Close the file
    CALL close_netcdf_file(netcdf%ncid)
    
  END SUBROUTINE create_global_scalar_output_file
  SUBROUTINE create_regional_scalar_output_file( region)
    ! Create a new regional scalar output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),          INTENT(INOUT) :: region

    ! Local variables:
    LOGICAL                                         :: file_exists
    INTEGER                                         :: t
    
    IF (.NOT. par%master) RETURN
    
    ! Set time frame index to 1
    region%scalars%ti = 1
    
    ! Generate filename
    region%scalars%filename = TRIM(C%output_dir)//'/scalar_output_'//TRIM(regioN%name)//'.nc'

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(region%scalars%filename))
    IF (file_exists) THEN
      WRITE(0,*) '  create_global_scalar_output_file - ERROR: ', TRIM(region%scalars%filename), ' already exists!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Create netCDF file
    CALL handle_error(nf90_create(region%scalars%filename,IOR(nf90_clobber,nf90_share),region%scalars%ncid))
        
    ! Define dimensions:
    CALL create_dim( region%scalars%ncid, region%scalars%name_dim_time,      nf90_unlimited,     region%scalars%id_dim_time      ) ! Time frames
    
    ! Placeholders for the dimension ID's, for shorter code
    t = region%scalars%id_dim_time
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.
    
    ! Dimension variables: zeta, month, time
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_time,          [t], region%scalars%id_var_time,          long_name='Time', units='years'   )
    
    ! Variables
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_ice_volume,    [t], region%scalars%id_var_ice_volume,    long_name='Ice volume', units='m.s.l.e')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_ice_volume_af, [t], region%scalars%id_var_ice_volume_af, long_name='Ice volume above flotation', units='m.s.l.e')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_ice_area,      [t], region%scalars%id_var_ice_area,      long_name='Ice volume', units='km^2')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_T2m,           [t], region%scalars%id_var_T2m,           long_name='Regionally averaged annual mean surface temperature', units='K')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_snowfall,      [t], region%scalars%id_var_snowfall,      long_name='Ice-sheet integrated snowfall', units='Gigaton yr^-1')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_rainfall,      [t], region%scalars%id_var_rainfall,      long_name='Ice-sheet integrated rainfall', units='Gigaton yr^-1')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_melt,          [t], region%scalars%id_var_melt,          long_name='Ice-sheet integrated melt', units='Gigaton yr^-1')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_refreezing,    [t], region%scalars%id_var_refreezing,    long_name='Ice-sheet integrated refreezing', units='Gigaton yr^-1')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_runoff,        [t], region%scalars%id_var_runoff,        long_name='Ice-sheet integrated runoff', units='Gigaton yr^-1')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_SMB,           [t], region%scalars%id_var_SMB,           long_name='Ice-sheet integrated surface mass balance', units='Gigaton yr^-1')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_BMB,           [t], region%scalars%id_var_BMB,           long_name='Ice-sheet integrated basal mass balance', units='Gigaton yr^-1')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_MB,            [t], region%scalars%id_var_MB,            long_name='Ice-sheet integrated mass balance', units='Gigaton yr^-1')
    
    ! Leave definition mode:
    CALL handle_error(nf90_enddef( region%scalars%ncid))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( region%scalars%ncid))
    
    ! Close the file
    CALL close_netcdf_file( region%scalars%ncid)
    
  END SUBROUTINE create_regional_scalar_output_file
  
  ! Manage memory for the debug data fields
  SUBROUTINE associate_debug_fields( region)
    ! Since the dimensions vary, each region needs its own set of debug fields. However, if
    ! we make them part of the "region" TYPE, they need to be passed to every subroutine as an
    ! argument before they can be used, which is a lot of hassle. So instead they are saved as
    ! global variables of this module, where they can be accessed from anywhere. This is done
    ! via the "intermediary" set of pointers, which are bound to the region-specific debug structure
    ! with this here subroutine.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(IN)    :: region
    
    ! Copy the netcdf ID's
    IF (par%master) THEN
      IF     (region%name == 'NAM') THEN
        debug%netcdf = debug_NAM%netcdf
      ELSEIF (region%name == 'EAS') THEN
        debug%netcdf = debug_EAS%netcdf
      ELSEIF (region%name == 'GRL') THEN
        debug%netcdf = debug_GRL%netcdf
      ELSEIF (region%name == 'ANT') THEN
        debug%netcdf = debug_ANT%netcdf
      END IF
    END IF
    CALL sync
    
    ! If necessary (i.e. every time except the first ever time this subroutine is called), de-associate the intermediary pointers first.
    IF (ASSOCIATED(debug%dp_2D_01)) THEN
    
      NULLIFY( debug%nx)
      NULLIFY( debug%ny)
    
      NULLIFY( debug%int_2D_01)
      NULLIFY( debug%int_2D_02)
      NULLIFY( debug%int_2D_03)
      NULLIFY( debug%int_2D_04)
      NULLIFY( debug%int_2D_05)
      NULLIFY( debug%int_2D_06)
      NULLIFY( debug%int_2D_07)
      NULLIFY( debug%int_2D_08)
      NULLIFY( debug%int_2D_09)
      NULLIFY( debug%int_2D_10)
      
      NULLIFY( debug%dp_2D_01)
      NULLIFY( debug%dp_2D_02)
      NULLIFY( debug%dp_2D_03)
      NULLIFY( debug%dp_2D_04)
      NULLIFY( debug%dp_2D_05)
      NULLIFY( debug%dp_2D_06)
      NULLIFY( debug%dp_2D_07)
      NULLIFY( debug%dp_2D_08)
      NULLIFY( debug%dp_2D_09)
      NULLIFY( debug%dp_2D_10)
      NULLIFY( debug%dp_2D_11)
      NULLIFY( debug%dp_2D_12)
      NULLIFY( debug%dp_2D_13)
      NULLIFY( debug%dp_2D_14)
      NULLIFY( debug%dp_2D_15)
      NULLIFY( debug%dp_2D_16)
      NULLIFY( debug%dp_2D_17)
      NULLIFY( debug%dp_2D_18)
      NULLIFY( debug%dp_2D_19)
      NULLIFY( debug%dp_2D_20)
      
      NULLIFY( debug%dp_3D_01)
      NULLIFY( debug%dp_3D_02)
      NULLIFY( debug%dp_3D_03)
      NULLIFY( debug%dp_3D_04)
      NULLIFY( debug%dp_3D_05)
      NULLIFY( debug%dp_3D_06)
      NULLIFY( debug%dp_3D_07)
      NULLIFY( debug%dp_3D_08)
      NULLIFY( debug%dp_3D_09)
      NULLIFY( debug%dp_3D_10)
      
      NULLIFY( debug%dp_2D_monthly_01)
      NULLIFY( debug%dp_2D_monthly_02)
      NULLIFY( debug%dp_2D_monthly_03)
      NULLIFY( debug%dp_2D_monthly_04)
      NULLIFY( debug%dp_2D_monthly_05)
      NULLIFY( debug%dp_2D_monthly_06)
      NULLIFY( debug%dp_2D_monthly_07)
      NULLIFY( debug%dp_2D_monthly_08)
      NULLIFY( debug%dp_2D_monthly_09)
      NULLIFY( debug%dp_2D_monthly_10)
      
    END IF
    
    ! Bind to the actual memory for this region
    IF (region%name == 'NAM') THEN
    
      debug%nx => debug_NAM%nx
      debug%ny => debug_NAM%ny
    
      debug%int_2D_01 => debug_NAM%int_2D_01    
      debug%int_2D_02 => debug_NAM%int_2D_02    
      debug%int_2D_03 => debug_NAM%int_2D_03    
      debug%int_2D_04 => debug_NAM%int_2D_04    
      debug%int_2D_05 => debug_NAM%int_2D_05    
      debug%int_2D_06 => debug_NAM%int_2D_06    
      debug%int_2D_07 => debug_NAM%int_2D_07    
      debug%int_2D_08 => debug_NAM%int_2D_08    
      debug%int_2D_09 => debug_NAM%int_2D_09    
      debug%int_2D_10 => debug_NAM%int_2D_10
    
      debug%dp_2D_01 => debug_NAM%dp_2D_01    
      debug%dp_2D_02 => debug_NAM%dp_2D_02    
      debug%dp_2D_03 => debug_NAM%dp_2D_03    
      debug%dp_2D_04 => debug_NAM%dp_2D_04    
      debug%dp_2D_05 => debug_NAM%dp_2D_05    
      debug%dp_2D_06 => debug_NAM%dp_2D_06    
      debug%dp_2D_07 => debug_NAM%dp_2D_07    
      debug%dp_2D_08 => debug_NAM%dp_2D_08    
      debug%dp_2D_09 => debug_NAM%dp_2D_09    
      debug%dp_2D_10 => debug_NAM%dp_2D_10
      debug%dp_2D_11 => debug_NAM%dp_2D_11    
      debug%dp_2D_12 => debug_NAM%dp_2D_12    
      debug%dp_2D_13 => debug_NAM%dp_2D_13    
      debug%dp_2D_14 => debug_NAM%dp_2D_14    
      debug%dp_2D_15 => debug_NAM%dp_2D_15    
      debug%dp_2D_16 => debug_NAM%dp_2D_16    
      debug%dp_2D_17 => debug_NAM%dp_2D_17    
      debug%dp_2D_18 => debug_NAM%dp_2D_18    
      debug%dp_2D_19 => debug_NAM%dp_2D_19    
      debug%dp_2D_20 => debug_NAM%dp_2D_20
      
      debug%dp_3D_01 => debug_NAM%dp_3D_01    
      debug%dp_3D_02 => debug_NAM%dp_3D_02    
      debug%dp_3D_03 => debug_NAM%dp_3D_03    
      debug%dp_3D_04 => debug_NAM%dp_3D_04    
      debug%dp_3D_05 => debug_NAM%dp_3D_05    
      debug%dp_3D_06 => debug_NAM%dp_3D_06    
      debug%dp_3D_07 => debug_NAM%dp_3D_07    
      debug%dp_3D_08 => debug_NAM%dp_3D_08    
      debug%dp_3D_09 => debug_NAM%dp_3D_09    
      debug%dp_3D_10 => debug_NAM%dp_3D_10
      
      debug%dp_2D_monthly_01 => debug_NAM%dp_2D_monthly_01    
      debug%dp_2D_monthly_02 => debug_NAM%dp_2D_monthly_02    
      debug%dp_2D_monthly_03 => debug_NAM%dp_2D_monthly_03    
      debug%dp_2D_monthly_04 => debug_NAM%dp_2D_monthly_04    
      debug%dp_2D_monthly_05 => debug_NAM%dp_2D_monthly_05    
      debug%dp_2D_monthly_06 => debug_NAM%dp_2D_monthly_06    
      debug%dp_2D_monthly_07 => debug_NAM%dp_2D_monthly_07    
      debug%dp_2D_monthly_08 => debug_NAM%dp_2D_monthly_08    
      debug%dp_2D_monthly_09 => debug_NAM%dp_2D_monthly_09    
      debug%dp_2D_monthly_10 => debug_NAM%dp_2D_monthly_10
      
    ELSEIF (region%name == 'EAS') THEN
    
      debug%nx => debug_EAS%nx
      debug%ny => debug_EAS%ny
    
      debug%int_2D_01 => debug_EAS%int_2D_01    
      debug%int_2D_02 => debug_EAS%int_2D_02    
      debug%int_2D_03 => debug_EAS%int_2D_03    
      debug%int_2D_04 => debug_EAS%int_2D_04    
      debug%int_2D_05 => debug_EAS%int_2D_05    
      debug%int_2D_06 => debug_EAS%int_2D_06    
      debug%int_2D_07 => debug_EAS%int_2D_07    
      debug%int_2D_08 => debug_EAS%int_2D_08    
      debug%int_2D_09 => debug_EAS%int_2D_09    
      debug%int_2D_10 => debug_EAS%int_2D_10
    
      debug%dp_2D_01 => debug_EAS%dp_2D_01    
      debug%dp_2D_02 => debug_EAS%dp_2D_02    
      debug%dp_2D_03 => debug_EAS%dp_2D_03    
      debug%dp_2D_04 => debug_EAS%dp_2D_04    
      debug%dp_2D_05 => debug_EAS%dp_2D_05    
      debug%dp_2D_06 => debug_EAS%dp_2D_06    
      debug%dp_2D_07 => debug_EAS%dp_2D_07    
      debug%dp_2D_08 => debug_EAS%dp_2D_08    
      debug%dp_2D_09 => debug_EAS%dp_2D_09    
      debug%dp_2D_10 => debug_EAS%dp_2D_10
      debug%dp_2D_11 => debug_EAS%dp_2D_11    
      debug%dp_2D_12 => debug_EAS%dp_2D_12    
      debug%dp_2D_13 => debug_EAS%dp_2D_13    
      debug%dp_2D_14 => debug_EAS%dp_2D_14    
      debug%dp_2D_15 => debug_EAS%dp_2D_15    
      debug%dp_2D_16 => debug_EAS%dp_2D_16    
      debug%dp_2D_17 => debug_EAS%dp_2D_17    
      debug%dp_2D_18 => debug_EAS%dp_2D_18    
      debug%dp_2D_19 => debug_EAS%dp_2D_19    
      debug%dp_2D_20 => debug_EAS%dp_2D_20
      
      debug%dp_3D_01 => debug_EAS%dp_3D_01    
      debug%dp_3D_02 => debug_EAS%dp_3D_02    
      debug%dp_3D_03 => debug_EAS%dp_3D_03    
      debug%dp_3D_04 => debug_EAS%dp_3D_04    
      debug%dp_3D_05 => debug_EAS%dp_3D_05    
      debug%dp_3D_06 => debug_EAS%dp_3D_06    
      debug%dp_3D_07 => debug_EAS%dp_3D_07    
      debug%dp_3D_08 => debug_EAS%dp_3D_08    
      debug%dp_3D_09 => debug_EAS%dp_3D_09    
      debug%dp_3D_10 => debug_EAS%dp_3D_10
      
      debug%dp_2D_monthly_01 => debug_EAS%dp_2D_monthly_01    
      debug%dp_2D_monthly_02 => debug_EAS%dp_2D_monthly_02    
      debug%dp_2D_monthly_03 => debug_EAS%dp_2D_monthly_03    
      debug%dp_2D_monthly_04 => debug_EAS%dp_2D_monthly_04    
      debug%dp_2D_monthly_05 => debug_EAS%dp_2D_monthly_05    
      debug%dp_2D_monthly_06 => debug_EAS%dp_2D_monthly_06    
      debug%dp_2D_monthly_07 => debug_EAS%dp_2D_monthly_07    
      debug%dp_2D_monthly_08 => debug_EAS%dp_2D_monthly_08    
      debug%dp_2D_monthly_09 => debug_EAS%dp_2D_monthly_09    
      debug%dp_2D_monthly_10 => debug_EAS%dp_2D_monthly_10
      
    ELSEIF (region%name == 'GRL') THEN
    
      debug%nx => debug_GRL%nx
      debug%ny => debug_GRL%ny
    
      debug%int_2D_01 => debug_GRL%int_2D_01    
      debug%int_2D_02 => debug_GRL%int_2D_02    
      debug%int_2D_03 => debug_GRL%int_2D_03    
      debug%int_2D_04 => debug_GRL%int_2D_04    
      debug%int_2D_05 => debug_GRL%int_2D_05    
      debug%int_2D_06 => debug_GRL%int_2D_06    
      debug%int_2D_07 => debug_GRL%int_2D_07    
      debug%int_2D_08 => debug_GRL%int_2D_08    
      debug%int_2D_09 => debug_GRL%int_2D_09    
      debug%int_2D_10 => debug_GRL%int_2D_10
    
      debug%dp_2D_01 => debug_GRL%dp_2D_01    
      debug%dp_2D_02 => debug_GRL%dp_2D_02    
      debug%dp_2D_03 => debug_GRL%dp_2D_03    
      debug%dp_2D_04 => debug_GRL%dp_2D_04    
      debug%dp_2D_05 => debug_GRL%dp_2D_05    
      debug%dp_2D_06 => debug_GRL%dp_2D_06    
      debug%dp_2D_07 => debug_GRL%dp_2D_07    
      debug%dp_2D_08 => debug_GRL%dp_2D_08    
      debug%dp_2D_09 => debug_GRL%dp_2D_09    
      debug%dp_2D_10 => debug_GRL%dp_2D_10
      debug%dp_2D_11 => debug_GRL%dp_2D_11    
      debug%dp_2D_12 => debug_GRL%dp_2D_12    
      debug%dp_2D_13 => debug_GRL%dp_2D_13    
      debug%dp_2D_14 => debug_GRL%dp_2D_14    
      debug%dp_2D_15 => debug_GRL%dp_2D_15    
      debug%dp_2D_16 => debug_GRL%dp_2D_16    
      debug%dp_2D_17 => debug_GRL%dp_2D_17    
      debug%dp_2D_18 => debug_GRL%dp_2D_18    
      debug%dp_2D_19 => debug_GRL%dp_2D_19    
      debug%dp_2D_20 => debug_GRL%dp_2D_20
      
      debug%dp_3D_01 => debug_GRL%dp_3D_01    
      debug%dp_3D_02 => debug_GRL%dp_3D_02    
      debug%dp_3D_03 => debug_GRL%dp_3D_03    
      debug%dp_3D_04 => debug_GRL%dp_3D_04    
      debug%dp_3D_05 => debug_GRL%dp_3D_05    
      debug%dp_3D_06 => debug_GRL%dp_3D_06    
      debug%dp_3D_07 => debug_GRL%dp_3D_07    
      debug%dp_3D_08 => debug_GRL%dp_3D_08    
      debug%dp_3D_09 => debug_GRL%dp_3D_09    
      debug%dp_3D_10 => debug_GRL%dp_3D_10
      
      debug%dp_2D_monthly_01 => debug_GRL%dp_2D_monthly_01    
      debug%dp_2D_monthly_02 => debug_GRL%dp_2D_monthly_02    
      debug%dp_2D_monthly_03 => debug_GRL%dp_2D_monthly_03    
      debug%dp_2D_monthly_04 => debug_GRL%dp_2D_monthly_04    
      debug%dp_2D_monthly_05 => debug_GRL%dp_2D_monthly_05    
      debug%dp_2D_monthly_06 => debug_GRL%dp_2D_monthly_06    
      debug%dp_2D_monthly_07 => debug_GRL%dp_2D_monthly_07    
      debug%dp_2D_monthly_08 => debug_GRL%dp_2D_monthly_08    
      debug%dp_2D_monthly_09 => debug_GRL%dp_2D_monthly_09    
      debug%dp_2D_monthly_10 => debug_GRL%dp_2D_monthly_10
      
    ELSEIF (region%name == 'ANT') THEN
    
      debug%nx => debug_ANT%nx
      debug%ny => debug_ANT%ny
    
      debug%int_2D_01 => debug_ANT%int_2D_01    
      debug%int_2D_02 => debug_ANT%int_2D_02    
      debug%int_2D_03 => debug_ANT%int_2D_03    
      debug%int_2D_04 => debug_ANT%int_2D_04    
      debug%int_2D_05 => debug_ANT%int_2D_05    
      debug%int_2D_06 => debug_ANT%int_2D_06    
      debug%int_2D_07 => debug_ANT%int_2D_07    
      debug%int_2D_08 => debug_ANT%int_2D_08    
      debug%int_2D_09 => debug_ANT%int_2D_09    
      debug%int_2D_10 => debug_ANT%int_2D_10
    
      debug%dp_2D_01 => debug_ANT%dp_2D_01    
      debug%dp_2D_02 => debug_ANT%dp_2D_02    
      debug%dp_2D_03 => debug_ANT%dp_2D_03    
      debug%dp_2D_04 => debug_ANT%dp_2D_04    
      debug%dp_2D_05 => debug_ANT%dp_2D_05    
      debug%dp_2D_06 => debug_ANT%dp_2D_06    
      debug%dp_2D_07 => debug_ANT%dp_2D_07    
      debug%dp_2D_08 => debug_ANT%dp_2D_08    
      debug%dp_2D_09 => debug_ANT%dp_2D_09    
      debug%dp_2D_10 => debug_ANT%dp_2D_10
      debug%dp_2D_11 => debug_ANT%dp_2D_11    
      debug%dp_2D_12 => debug_ANT%dp_2D_12    
      debug%dp_2D_13 => debug_ANT%dp_2D_13    
      debug%dp_2D_14 => debug_ANT%dp_2D_14    
      debug%dp_2D_15 => debug_ANT%dp_2D_15    
      debug%dp_2D_16 => debug_ANT%dp_2D_16    
      debug%dp_2D_17 => debug_ANT%dp_2D_17    
      debug%dp_2D_18 => debug_ANT%dp_2D_18    
      debug%dp_2D_19 => debug_ANT%dp_2D_19    
      debug%dp_2D_20 => debug_ANT%dp_2D_20
      
      debug%dp_3D_01 => debug_ANT%dp_3D_01    
      debug%dp_3D_02 => debug_ANT%dp_3D_02    
      debug%dp_3D_03 => debug_ANT%dp_3D_03    
      debug%dp_3D_04 => debug_ANT%dp_3D_04    
      debug%dp_3D_05 => debug_ANT%dp_3D_05    
      debug%dp_3D_06 => debug_ANT%dp_3D_06    
      debug%dp_3D_07 => debug_ANT%dp_3D_07    
      debug%dp_3D_08 => debug_ANT%dp_3D_08    
      debug%dp_3D_09 => debug_ANT%dp_3D_09    
      debug%dp_3D_10 => debug_ANT%dp_3D_10
      
      debug%dp_2D_monthly_01 => debug_ANT%dp_2D_monthly_01    
      debug%dp_2D_monthly_02 => debug_ANT%dp_2D_monthly_02    
      debug%dp_2D_monthly_03 => debug_ANT%dp_2D_monthly_03    
      debug%dp_2D_monthly_04 => debug_ANT%dp_2D_monthly_04    
      debug%dp_2D_monthly_05 => debug_ANT%dp_2D_monthly_05    
      debug%dp_2D_monthly_06 => debug_ANT%dp_2D_monthly_06    
      debug%dp_2D_monthly_07 => debug_ANT%dp_2D_monthly_07    
      debug%dp_2D_monthly_08 => debug_ANT%dp_2D_monthly_08    
      debug%dp_2D_monthly_09 => debug_ANT%dp_2D_monthly_09    
      debug%dp_2D_monthly_10 => debug_ANT%dp_2D_monthly_10
      
    END IF
    
  END SUBROUTINE associate_debug_fields
  SUBROUTINE initialise_debug_fields( region)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    
    IF     (region%name == 'NAM') THEN
      CALL initialise_debug_fields_region( debug_NAM, region%grid%nx, region%grid%ny)
    ELSEIF (region%name == 'EAS') THEN
      CALL initialise_debug_fields_region( debug_EAS, region%grid%nx, region%grid%ny)
    ELSEIF (region%name == 'GRL') THEN
      CALL initialise_debug_fields_region( debug_GRL, region%grid%nx, region%grid%ny)
    ELSEIF (region%name == 'ANT') THEN
      CALL initialise_debug_fields_region( debug_ANT, region%grid%nx, region%grid%ny)
    END IF
    
  END SUBROUTINE initialise_debug_fields
  SUBROUTINE initialise_debug_fields_region( debug, nx, ny)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_debug_fields),         INTENT(INOUT)     :: debug
    INTEGER,                         INTENT(IN)        :: nx, ny
    
    ! Grid size
    CALL allocate_shared_int_0D( debug%nx, debug%wnx)
    CALL allocate_shared_int_0D( debug%ny, debug%wny)
    
    IF (par%master) THEN
      debug%nx = nx
      debug%ny = ny
    END IF
    CALL sync
    
    ! Data
    CALL allocate_shared_int_2D(      ny, nx, debug%int_2D_01,        debug%wint_2D_01       )
    CALL allocate_shared_int_2D(      ny, nx, debug%int_2D_02,        debug%wint_2D_02       )
    CALL allocate_shared_int_2D(      ny, nx, debug%int_2D_03,        debug%wint_2D_03       )
    CALL allocate_shared_int_2D(      ny, nx, debug%int_2D_04,        debug%wint_2D_04       )
    CALL allocate_shared_int_2D(      ny, nx, debug%int_2D_05,        debug%wint_2D_05       )
    CALL allocate_shared_int_2D(      ny, nx, debug%int_2D_06,        debug%wint_2D_06       )
    CALL allocate_shared_int_2D(      ny, nx, debug%int_2D_07,        debug%wint_2D_07       )
    CALL allocate_shared_int_2D(      ny, nx, debug%int_2D_08,        debug%wint_2D_08       )
    CALL allocate_shared_int_2D(      ny, nx, debug%int_2D_09,        debug%wint_2D_09       )
    CALL allocate_shared_int_2D(      ny, nx, debug%int_2D_10,        debug%wint_2D_10       )
    
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_01,         debug%wdp_2D_01        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_02,         debug%wdp_2D_02        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_03,         debug%wdp_2D_03        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_04,         debug%wdp_2D_04        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_05,         debug%wdp_2D_05        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_06,         debug%wdp_2D_06        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_07,         debug%wdp_2D_07        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_08,         debug%wdp_2D_08        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_09,         debug%wdp_2D_09        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_10,         debug%wdp_2D_10        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_11,         debug%wdp_2D_11        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_12,         debug%wdp_2D_12        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_13,         debug%wdp_2D_13        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_14,         debug%wdp_2D_14        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_15,         debug%wdp_2D_15        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_16,         debug%wdp_2D_16        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_17,         debug%wdp_2D_17        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_18,         debug%wdp_2D_18        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_19,         debug%wdp_2D_19        )
    CALL allocate_shared_dp_2D(       ny, nx, debug%dp_2D_20,         debug%wdp_2D_20        )
    
    CALL allocate_shared_dp_3D( C%nZ, ny, nx, debug%dp_3D_01,         debug%wdp_3D_01        )
    CALL allocate_shared_dp_3D( C%nZ, ny, nx, debug%dp_3D_02,         debug%wdp_3D_02        )
    CALL allocate_shared_dp_3D( C%nZ, ny, nx, debug%dp_3D_03,         debug%wdp_3D_03        )
    CALL allocate_shared_dp_3D( C%nZ, ny, nx, debug%dp_3D_04,         debug%wdp_3D_04        )
    CALL allocate_shared_dp_3D( C%nZ, ny, nx, debug%dp_3D_05,         debug%wdp_3D_05        )
    CALL allocate_shared_dp_3D( C%nZ, ny, nx, debug%dp_3D_06,         debug%wdp_3D_06        )
    CALL allocate_shared_dp_3D( C%nZ, ny, nx, debug%dp_3D_07,         debug%wdp_3D_07        )
    CALL allocate_shared_dp_3D( C%nZ, ny, nx, debug%dp_3D_08,         debug%wdp_3D_08        )
    CALL allocate_shared_dp_3D( C%nZ, ny, nx, debug%dp_3D_09,         debug%wdp_3D_09        )
    CALL allocate_shared_dp_3D( C%nZ, ny, nx, debug%dp_3D_10,         debug%wdp_3D_10        )
    
    CALL allocate_shared_dp_3D( 12,   ny, nx, debug%dp_2D_monthly_01, debug%wdp_2D_monthly_01)
    CALL allocate_shared_dp_3D( 12,   ny, nx, debug%dp_2D_monthly_02, debug%wdp_2D_monthly_02)
    CALL allocate_shared_dp_3D( 12,   ny, nx, debug%dp_2D_monthly_03, debug%wdp_2D_monthly_03)
    CALL allocate_shared_dp_3D( 12,   ny, nx, debug%dp_2D_monthly_04, debug%wdp_2D_monthly_04)
    CALL allocate_shared_dp_3D( 12,   ny, nx, debug%dp_2D_monthly_05, debug%wdp_2D_monthly_05)
    CALL allocate_shared_dp_3D( 12,   ny, nx, debug%dp_2D_monthly_06, debug%wdp_2D_monthly_06)
    CALL allocate_shared_dp_3D( 12,   ny, nx, debug%dp_2D_monthly_07, debug%wdp_2D_monthly_07)
    CALL allocate_shared_dp_3D( 12,   ny, nx, debug%dp_2D_monthly_08, debug%wdp_2D_monthly_08)
    CALL allocate_shared_dp_3D( 12,   ny, nx, debug%dp_2D_monthly_09, debug%wdp_2D_monthly_09)
    CALL allocate_shared_dp_3D( 12,   ny, nx, debug%dp_2D_monthly_10, debug%wdp_2D_monthly_10)
    
  END SUBROUTINE initialise_debug_fields_region
  
! Read data to restart a run
! ==========================
  
  SUBROUTINE inquire_restart_file( init) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields), INTENT(INOUT) :: init
 
    ! Local variables:
    INTEGER                               :: int_dummy
    
    IF (.NOT. par%master) RETURN
        
    ! Open the netcdf file
    CALL open_netcdf_file(init%netcdf%filename, init%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist, return their lengths.
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_x,     init%nx,   init%netcdf%id_dim_x    )
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_y,     init%ny,   init%netcdf%id_dim_y    )
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_zeta,  init%nz,   init%netcdf%id_dim_zeta )
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_time,  init%nt,   init%netcdf%id_dim_time )
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_month, int_dummy, init%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_x,                (/ init%netcdf%id_dim_x                                                                          /), init%netcdf%id_var_x   )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_y,                (/                       init%netcdf%id_dim_y                                                    /), init%netcdf%id_var_y   )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_zeta,             (/                                             init%netcdf%id_dim_zeta                           /), init%netcdf%id_var_zeta)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_time,             (/                                                                       init%netcdf%id_dim_time /), init%netcdf%id_var_time)
    
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Hi,               (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y                          , init%netcdf%id_dim_time /), init%netcdf%id_var_Hi  )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Hb,               (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y                          , init%netcdf%id_dim_time /), init%netcdf%id_var_Hb  )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_SL,               (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y                          , init%netcdf%id_dim_time /), init%netcdf%id_var_SL  )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_dHb,              (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y                          , init%netcdf%id_dim_time /), init%netcdf%id_var_dHb )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Ti,               (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y, init%netcdf%id_dim_zeta , init%netcdf%id_dim_time /), init%netcdf%id_var_Ti  )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_FirnDepth,        (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y, init%netcdf%id_dim_month, init%netcdf%id_dim_time /), init%netcdf%id_var_FirnDepth)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_MeltPreviousYear, (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y                          , init%netcdf%id_dim_time /), init%netcdf%id_var_MeltPreviousYear)
        
    ! Close the netcdf file
    CALL close_netcdf_file(init%netcdf%ncid)
    
  END SUBROUTINE inquire_restart_file
  SUBROUTINE read_restart_file(    init)
    ! Read the init netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields), INTENT(INOUT) :: init

    ! Local variables:
    INTEGER                                       :: k,ti,ti_min
    REAL(dp)                                      :: dt, dt_min
    
    IF (.NOT. par%master) RETURN
    
    ! Open the netcdf file
    CALL open_netcdf_file(init%netcdf%filename, init%netcdf%ncid)
    
    ! Read zeta, check if it matches the config zeta
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_zeta, init%zeta, start=(/1/) ))
    IF (init%nz /= C%nz) THEN
      WRITE(0,*) '  read_restart_file - ERROR: vertical coordinate zeta in restart file doesnt match zeta in config!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSE
      DO k = 1, C%nz
        IF (ABS(C%zeta(k) - init%zeta(k)) > 0.0001_dp) THEN
      WRITE(0,*) '  read_restart_file - ERROR: vertical coordinate zeta in restart file doesnt match zeta in config!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END DO
    END IF
    
    ! Read x,y
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_x, init%x, start=(/1/) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_y, init%y, start=(/1/) ))
    
    ! Read time, determine which time frame to read
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_time, init%time, start=(/1/) ))
    
    IF (C%time_to_restart_from < MINVAL(init%time) .OR. C%time_to_restart_from > MAXVAL(init%time)) THEN
      WRITE(0,*) '  ERROR - time_to_restart_from ', C%time_to_restart_from, ' outside range of restart file!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF    
    
    ti_min = 0
    dt_min = 1E8_dp
    DO ti = 1, init%nt
      dt = ABS(init%time( ti) - C%time_to_restart_from)
      IF (dt < dt_min) THEN
        ti_min = ti
        dt_min = dt
      END IF
    END DO
    ti = ti_min
    
    IF (dt_min > 0._dp) THEN
      WRITE(0,*) ' ======== '
      WRITE(0,*) '  WARNING - no exact match for time_to_restart_from ', C%time_to_restart_from, ' in restart file! Reading closest match ', init%time( ti), ' instead.'
      WRITE(0,*) ' ======== '
    END IF
    IF (C%time_to_restart_from /= C%start_time_of_run) THEN
      WRITE(0,*) ' ======== '
      WRITE(0,*) '  WARNING - starting run at t = ', C%start_time_of_run, ' with restart data at t = ', C%time_to_restart_from
      WRITE(0,*) ' ======== '
    END IF
    
    ! Read the data
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hi,               init%Hi_raw,               start = (/ 1, 1,    ti /), count = (/ init%nx, init%ny,          1 /) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hb,               init%Hb_raw,               start = (/ 1, 1,    ti /), count = (/ init%nx, init%ny,          1 /) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_SL,               init%SL_raw,               start = (/ 1, 1,    ti /), count = (/ init%nx, init%ny,          1 /) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_dHb,              init%dHb_raw,              start = (/ 1, 1,    ti /), count = (/ init%nx, init%ny,          1 /) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Ti,               init%Ti_raw,               start = (/ 1, 1, 1, ti /), count = (/ init%nx, init%ny, init%nz, 1 /) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_FirnDepth,        init%FirnDepth_raw,        start = (/ 1, 1, 1, ti /), count = (/ init%nx, init%ny, 12,      1 /) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_MeltPreviousYear, init%MeltPreviousYear_raw, start = (/ 1, 1,    ti /), count = (/ init%nx, init%ny,          1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file(init%netcdf%ncid)
    
  END SUBROUTINE read_restart_file
  
  SUBROUTINE read_inverse_routine_history_dT_glob(         forcing, filename)
    ! Read the inverse routine history from the specified NetCDF file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    CHARACTER(LEN=256),             INTENT(IN)    :: filename
    
    ! Local variables:
    INTEGER                                       :: ncid, id_dim_nH, id_dim_time, nt, id_var_H, id_var_time
    INTEGER                                       :: int_dummy
    INTEGER                                       :: ti, ti_min
    REAL(dp)                                      :: dt, dt_min
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: time
    
    IF (.NOT. par%master) RETURN
    
    WRITE(0,*) ' Reading inverse routine dT_glob_history from file "', TRIM(filename), '"...'
    
    ! Check if this file contains the required history
        
    ! Open the netcdf file
    CALL open_netcdf_file( filename, ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( ncid, 'ndT_glob_history', int_dummy, id_dim_nH)
    IF (int_dummy /= forcing%ndT_glob_history) THEN
      WRITE(0,*) '  ERROR - ndT_glob_history in file "', filename, '" doesnt match forcing configuration!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    CALL inquire_dim( ncid, 'time', nt, id_dim_time)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( ncid, 'dT_glob_history',  (/ id_dim_nH, id_dim_time/), id_var_H)
    CALL inquire_double_var( ncid, 'time',             (/            id_dim_time/), id_var_time)
    
    ! Read time, determine which time frame to read
    ALLOCATE( time( nt))
    
    CALL handle_error(nf90_get_var( ncid, id_var_time, time))
    
    IF (C%time_to_restart_from < MINVAL(time) .OR. C%time_to_restart_from > MAXVAL(time)) THEN
      WRITE(0,*) '  ERROR - time_to_restart_from ', C%time_to_restart_from, ' outside range of restart file!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF    
    
    ti_min = 0
    dt_min = 1E8_dp
    DO ti = 1, nt
      dt = ABS(time( ti) - C%time_to_restart_from)
      IF (dt < dt_min) THEN
        ti_min = ti
        dt_min = dt
      END IF
    END DO
    ti = ti_min
    
    IF (dt_min > 0._dp) THEN
      WRITE(0,*) '  WARNING - no exact match for time_to_restart_from ', C%time_to_restart_from, ' in restart file! Reading closest match ', time( ti), ' instead.'
    END IF
    IF (C%time_to_restart_from /= C%start_time_of_run) THEN
      WRITE(0,*) '  WARNING - starting run at t = ', C%start_time_of_run, ' with restart data at t = ', C%time_to_restart_from
    END IF
    
    DEALLOCATE( time)
    
    ! Read inverse routine history
    CALL handle_error(nf90_get_var( ncid, id_var_H, forcing%dT_glob_history, start=(/1,ti/) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( ncid)
    
  END SUBROUTINE read_inverse_routine_history_dT_glob
  SUBROUTINE read_inverse_routine_history_dT_glob_inverse( forcing, filename)
    ! Read the inverse routine history from the specified NetCDF file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    CHARACTER(LEN=256),             INTENT(IN)    :: filename
    
    ! Local variables:
    INTEGER                                       :: ncid, id_dim_nH, id_dim_time, nt, id_var_H, id_var_time
    INTEGER                                       :: int_dummy
    INTEGER                                       :: ti, ti_min
    REAL(dp)                                      :: dt, dt_min
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: time
    
    IF (.NOT. par%master) RETURN
    
    WRITE(0,*) ' Reading inverse routine dT_glob_inverse_history from file "', TRIM(filename), '"...'
    
    ! Check if this file contains the required history
        
    ! Open the netcdf file
    CALL open_netcdf_file( filename, ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( ncid, 'ndT_glob_inverse_history', int_dummy, id_dim_nH)
    IF (int_dummy /= forcing%ndT_glob_inverse_history) THEN
      WRITE(0,*) '  ERROR - ndT_glob_inverse_history in file "', filename, '" doesnt match forcing configuration!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    CALL inquire_dim( ncid, 'time', nt, id_dim_time)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( ncid, 'dT_glob_inverse_history',  (/ id_dim_nH, id_dim_time/), id_var_H)
    CALL inquire_double_var( ncid, 'time',                     (/            id_dim_time/), id_var_time)
    
    ! Read time, determine which time frame to read
    ALLOCATE( time( nt))
    
    CALL handle_error(nf90_get_var( ncid, id_var_time, time))
    
    IF (C%time_to_restart_from < MINVAL(time) .OR. C%time_to_restart_from > MAXVAL(time)) THEN
      WRITE(0,*) '  ERROR - time_to_restart_from ', C%time_to_restart_from, ' outside range of restart file!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF    
    
    ti_min = 0
    dt_min = 1E8_dp
    DO ti = 1, nt
      dt = ABS(time( ti) - C%time_to_restart_from)
      IF (dt < dt_min) THEN
        ti_min = ti
        dt_min = dt
      END IF
    END DO
    ti = ti_min
    
    IF (dt_min > 0._dp) THEN
      WRITE(0,*) '  WARNING - no exact match for time_to_restart_from ', C%time_to_restart_from, ' in restart file! Reading closest match ', time( ti), ' instead.'
    END IF
    IF (C%time_to_restart_from /= C%start_time_of_run) THEN
      WRITE(0,*) '  WARNING - starting run at t = ', C%start_time_of_run, ' with restart data at t = ', C%time_to_restart_from
    END IF
    
    DEALLOCATE( time)
    
    ! Read inverse routine history
    CALL handle_error(nf90_get_var( ncid, id_var_H, forcing%dT_glob_inverse_history, start=(/1,ti/) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( ncid)
    
  END SUBROUTINE read_inverse_routine_history_dT_glob_inverse
  SUBROUTINE read_inverse_routine_history_CO2_inverse(     forcing, filename)
    ! Read the inverse routine history from the specified NetCDF file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    CHARACTER(LEN=256),             INTENT(IN)    :: filename
    
    ! Local variables:
    INTEGER                                       :: ncid, id_dim_nH, id_dim_time, nt, id_var_H, id_var_time
    INTEGER                                       :: int_dummy
    INTEGER                                       :: ti, ti_min
    REAL(dp)                                      :: dt, dt_min
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: time
    
    IF (.NOT. par%master) RETURN
    
    WRITE(0,*) ' Reading inverse routine CO2_inverse_history from file "', TRIM(filename), '"...'
    
    ! Check if this file contains the required history
        
    ! Open the netcdf file
    CALL open_netcdf_file( filename, ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( ncid, 'nCO2_inverse_history', int_dummy, id_dim_nH)
    IF (int_dummy /= forcing%nCO2_inverse_history) THEN
      WRITE(0,*) '  ERROR - nCO2_inverse_history in file "', filename, '" doesnt match forcing configuration!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    CALL inquire_dim( ncid, 'time', nt, id_dim_time)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( ncid, 'CO2_inverse_history',  (/ id_dim_nH, id_dim_time/), id_var_H)
    CALL inquire_double_var( ncid, 'time',                 (/            id_dim_time/), id_var_time)
    
    ! Read time, determine which time frame to read
    ALLOCATE( time( nt))
    
    CALL handle_error(nf90_get_var( ncid, id_var_time, time))
    
    IF (C%time_to_restart_from < MINVAL(time) .OR. C%time_to_restart_from > MAXVAL(time)) THEN
      WRITE(0,*) '  ERROR - time_to_restart_from ', C%time_to_restart_from, ' outside range of restart file!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF    
    
    ti_min = 0
    dt_min = 1E8_dp
    DO ti = 1, nt
      dt = ABS(time( ti) - C%time_to_restart_from)
      IF (dt < dt_min) THEN
        ti_min = ti
        dt_min = dt
      END IF
    END DO
    ti = ti_min
    
    IF (dt_min > 0._dp) THEN
      WRITE(0,*) '  WARNING - no exact match for time_to_restart_from ', C%time_to_restart_from, ' in restart file! Reading closest match ', time( ti), ' instead.'
    END IF
    IF (C%time_to_restart_from /= C%start_time_of_run) THEN
      WRITE(0,*) '  WARNING - starting run at t = ', C%start_time_of_run, ' with restart data at t = ', C%time_to_restart_from
    END IF
    
    DEALLOCATE( time)
    
    ! Read inverse routine history
    CALL handle_error(nf90_get_var( ncid, id_var_H, forcing%CO2_inverse_history, start=(/1,ti/) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( ncid)
    
  END SUBROUTINE read_inverse_routine_history_CO2_inverse
  
! Read all kinds of input files
! =============================
  
  ! Present-day observed geometry (ice thickness & bed topography)
  SUBROUTINE inquire_PD_data_file( PD) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_PD_data_fields), INTENT(INOUT) :: PD
    
    IF (.NOT. par%master) RETURN
        
    ! Open the netcdf file
    CALL open_netcdf_file(PD%netcdf%filename, PD%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( PD%netcdf%ncid, PD%netcdf%name_dim_x, PD%nx, PD%netcdf%id_dim_x)
    CALL inquire_dim( PD%netcdf%ncid, PD%netcdf%name_dim_y, PD%ny, PD%netcdf%id_dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_var_x,  (/ PD%netcdf%id_dim_x                     /), PD%netcdf%id_var_x )
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_var_y,  (/ PD%netcdf%id_dim_y                     /), PD%netcdf%id_var_y )
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_var_Hi, (/ PD%netcdf%id_dim_x, PD%netcdf%id_dim_y /), PD%netcdf%id_var_Hi)
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_var_Hb, (/ PD%netcdf%id_dim_x, PD%netcdf%id_dim_y /), PD%netcdf%id_var_Hb)
    CALL inquire_double_var( PD%netcdf%ncid, PD%netcdf%name_var_Hs, (/ PD%netcdf%id_dim_x, PD%netcdf%id_dim_y /), PD%netcdf%id_var_Hs)
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD%netcdf%ncid)
    
  END SUBROUTINE inquire_PD_data_file
  SUBROUTINE read_PD_data_file(    PD)
    ! Read the PD netcdf file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_PD_data_fields), INTENT(INOUT) :: PD
    
    IF (.NOT. par%master) RETURN
    
    ! Open the netcdf file
    CALL open_netcdf_file(PD%netcdf%filename, PD%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_x,      PD%x,      start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_y,      PD%y,      start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_Hi,     PD%Hi_raw, start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_Hb,     PD%Hb_raw, start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_Hs,     PD%Hs_raw, start = (/ 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD%netcdf%ncid)
    
  END SUBROUTINE read_PD_data_file 
  
  ! Initial observed geometry (ice thickness % bed topography)
  SUBROUTINE inquire_init_data_file( init) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields), INTENT(INOUT) :: init
    
    IF (.NOT. par%master) RETURN
        
    ! Open the netcdf file
    CALL open_netcdf_file(init%netcdf%filename, init%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_x, init%nx, init%netcdf%id_dim_x)
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_y, init%ny, init%netcdf%id_dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_x,  (/ init%netcdf%id_dim_x                       /), init%netcdf%id_var_x)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_y,  (/ init%netcdf%id_dim_y                       /), init%netcdf%id_var_y)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Hi, (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y /), init%netcdf%id_var_Hi)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Hb, (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y /), init%netcdf%id_var_Hb)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Hs, (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y /), init%netcdf%id_var_Hs)
        
    ! Close the netcdf file
    CALL close_netcdf_file(init%netcdf%ncid)
    
  END SUBROUTINE inquire_init_data_file
  SUBROUTINE read_init_data_file(    init)
    ! Read the init netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields), INTENT(INOUT) :: init
    
    IF (.NOT. par%master) RETURN
    
    ! Open the netcdf file
    CALL open_netcdf_file(init%netcdf%filename, init%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_x,      init%x,      start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_y,      init%y,      start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hi,     init%Hi_raw, start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hb,     init%Hb_raw, start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hs,     init%Hs_raw, start = (/ 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file(init%netcdf%ncid)
    
  END SUBROUTINE read_init_data_file
  
  ! Present-day observed global climate (e.g. ERA-40)
  SUBROUTINE inquire_PD_obs_data_file( PD_obs) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: PD_obs
 
    ! Local variables:
    INTEGER                               :: int_dummy
    
    IF (.NOT. par%master) RETURN
        
    ! Open the netcdf file
    CALL open_netcdf_file(PD_obs%netcdf%filename, PD_obs%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( PD_obs%netcdf%ncid, PD_obs%netcdf%name_dim_lat,     PD_obs%nlat,  PD_obs%netcdf%id_dim_lat)
    CALL inquire_dim( PD_obs%netcdf%ncid, PD_obs%netcdf%name_dim_lon,     PD_obs%nlon,  PD_obs%netcdf%id_dim_lon)
    CALL inquire_dim( PD_obs%netcdf%ncid, PD_obs%netcdf%name_dim_month,   int_dummy,    PD_obs%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_lat,      (/ PD_obs%netcdf%id_dim_lat                                                       /),  PD_obs%netcdf%id_var_lat)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_lon,      (/ PD_obs%netcdf%id_dim_lon                                                       /),  PD_obs%netcdf%id_var_lon)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_T2m,      (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_T2m)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Precip,   (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_Precip)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Hs,       (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat                             /),  PD_obs%netcdf%id_var_Hs)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Wind_WE,  (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_Wind_WE)
    CALL inquire_double_var( PD_obs%netcdf%ncid, PD_obs%netcdf%name_var_Wind_SN,  (/ PD_obs%netcdf%id_dim_lon, PD_obs%netcdf%id_dim_lat, PD_obs%netcdf%id_dim_month /),  PD_obs%netcdf%id_var_Wind_SN)
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD_obs%netcdf%ncid)
    
  END SUBROUTINE inquire_PD_obs_data_file
  SUBROUTINE read_PD_obs_data_file(    PD_obs)
    ! Read the PD_obs0 netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: PD_obs
    
    IF (.NOT. par%master) RETURN
    
    ! Open the netcdf file
    CALL open_netcdf_file(PD_obs%netcdf%filename, PD_obs%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_lon,     PD_obs%lon,     start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_lat,     PD_obs%lat,     start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_T2m,     PD_obs%T2m,     start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Precip,  PD_obs%Precip,  start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Hs,      PD_obs%Hs_ref,  start = (/ 1, 1    /) ))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Wind_WE, PD_obs%Wind_WE, start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Wind_SN, PD_obs%Wind_SN, start = (/ 1, 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD_obs%netcdf%ncid)
    
  END SUBROUTINE read_PD_obs_data_file
  
  ! GCM global climate (climate matrix snapshots)
  SUBROUTINE inquire_GCM_snapshot( snapshot) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: snapshot
 
    ! Local variables:
    INTEGER                                     :: int_dummy
    
    IF (.NOT. par%master) RETURN
        
    ! Open the netcdf file
    CALL open_netcdf_file( snapshot%netcdf%filename, snapshot%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( snapshot%netcdf%ncid, snapshot%netcdf%name_dim_lat,     snapshot%nlat,  snapshot%netcdf%id_dim_lat)
    CALL inquire_dim( snapshot%netcdf%ncid, snapshot%netcdf%name_dim_lon,     snapshot%nlon,  snapshot%netcdf%id_dim_lon)
    CALL inquire_dim( snapshot%netcdf%ncid, snapshot%netcdf%name_dim_month,   int_dummy,      snapshot%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_lat,      (/ snapshot%netcdf%id_dim_lat                                                           /),  snapshot%netcdf%id_var_lat)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_lon,      (/ snapshot%netcdf%id_dim_lon                                                           /),  snapshot%netcdf%id_var_lon)
    !CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Hi,       (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat                               /),  snapshot%netcdf%id_var_Hi)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Hs,       (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat                               /),  snapshot%netcdf%id_var_Hs)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_T2m,      (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat, snapshot%netcdf%id_dim_month /),  snapshot%netcdf%id_var_T2m)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Precip,   (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat, snapshot%netcdf%id_dim_month /),  snapshot%netcdf%id_var_Precip)
    !CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Wind_WE,  (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat, snapshot%netcdf%id_dim_month /),  snapshot%netcdf%id_var_Wind_WE)
    !CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Wind_SN,  (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat, snapshot%netcdf%id_dim_month /),  snapshot%netcdf%id_var_Wind_SN)
        
    ! Close the netcdf file
    CALL close_netcdf_file(snapshot%netcdf%ncid)
    
  END SUBROUTINE inquire_GCM_snapshot
  SUBROUTINE read_GCM_snapshot(    snapshot)
    ! Read the PD_obs0 netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: snapshot
    
    IF (.NOT. par%master) RETURN
    
    ! Open the netcdf file
    CALL open_netcdf_file(snapshot%netcdf%filename, snapshot%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_lon,     snapshot%lon,     start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_lat,     snapshot%lat,     start = (/ 1       /) ))
   !CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Hi,      snapshot%Hi,      start = (/ 1, 1    /) ))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Hs,      snapshot%Hs_ref,  start = (/ 1, 1    /) ))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_T2m,     snapshot%T2m,     start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Precip,  snapshot%Precip,  start = (/ 1, 1, 1 /) ))
   !CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Wind_WE, snapshot%Wind_WE, start = (/ 1, 1, 1 /) ))
   !CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Wind_SN, snapshot%Wind_SN, start = (/ 1, 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file(snapshot%netcdf%ncid)
    
  END SUBROUTINE read_GCM_snapshot
  
  ! ICE5G ice geometry (needed for GCM snapshots)
  SUBROUTINE inquire_ICE5G_data( ICE5G)
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_ICE5G_timeframe), INTENT(INOUT) :: ICE5G
    
    IF (.NOT. par%master) RETURN
        
    ! Open the netcdf file
    CALL open_netcdf_file( ICE5G%netcdf%filename, ICE5G%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( ICE5G%netcdf%ncid, ICE5G%netcdf%name_dim_lat, ICE5G%nlat, ICE5G%netcdf%id_dim_lat)
    CALL inquire_dim( ICE5G%netcdf%ncid, ICE5G%netcdf%name_dim_lon, ICE5G%nlon, ICE5G%netcdf%id_dim_lon)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_single_var( ICE5G%netcdf%ncid, ICE5G%netcdf%name_var_lat,      (/ ICE5G%netcdf%id_dim_lat                         /),  ICE5G%netcdf%id_var_lat     )
    CALL inquire_single_var( ICE5G%netcdf%ncid, ICE5G%netcdf%name_var_lon,      (/ ICE5G%netcdf%id_dim_lon                         /),  ICE5G%netcdf%id_var_lon     )
    CALL inquire_single_var( ICE5G%netcdf%ncid, ICE5G%netcdf%name_var_Hi,       (/ ICE5G%netcdf%id_dim_lon, ICE5G%netcdf%id_dim_lat/),  ICE5G%netcdf%id_var_Hi      )
    CALL inquire_single_var( ICE5G%netcdf%ncid, ICE5G%netcdf%name_var_Hb,       (/ ICE5G%netcdf%id_dim_lon, ICE5G%netcdf%id_dim_lat/),  ICE5G%netcdf%id_var_Hb      )
    CALL inquire_single_var( ICE5G%netcdf%ncid, ICE5G%netcdf%name_var_mask_ice, (/ ICE5G%netcdf%id_dim_lon, ICE5G%netcdf%id_dim_lat/),  ICE5G%netcdf%id_var_mask_ice)
        
    ! Close the netcdf file
    CALL close_netcdf_file(ICE5G%netcdf%ncid)
    
  END SUBROUTINE inquire_ICE5G_data
  SUBROUTINE read_ICE5G_data(    ICE5G)
    ! Read the PD_obs0 netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_ICE5G_timeframe), INTENT(INOUT) :: ICE5G
    
    IF (.NOT. par%master) RETURN
    
    ! Open the netcdf file
    CALL open_netcdf_file( ICE5G%netcdf%filename, ICE5G%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( ICE5G%netcdf%ncid, ICE5G%netcdf%id_var_lon,      ICE5G%lon,      start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( ICE5G%netcdf%ncid, ICE5G%netcdf%id_var_lat,      ICE5G%lat,      start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( ICE5G%netcdf%ncid, ICE5G%netcdf%id_var_Hi,       ICE5G%Hi,       start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( ICE5G%netcdf%ncid, ICE5G%netcdf%id_var_Hb,       ICE5G%Hb,       start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( ICE5G%netcdf%ncid, ICE5G%netcdf%id_var_mask_ice, ICE5G%mask_ice, start = (/ 1, 1 /) ))
    
    ! For some reason, "orog" in ICE5G is Hb+Hi (so like Hs without any oceans)
    ICE5G%Hb = ICE5G%Hb - ICE5G%Hi
        
    ! Close the netcdf file
    CALL close_netcdf_file(ICE5G%netcdf%ncid)
    
  END SUBROUTINE read_ICE5G_data
  
  ! Insolation solution (e.g. Laskar 2004)
  SUBROUTINE inquire_insolation_data_file( forcing)
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
 
    ! Local variables:  
    INTEGER                                :: int_dummy
    
    IF (.NOT. par%master) RETURN
            
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_ins%filename, forcing%netcdf_ins%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_dim_time,     forcing%ins_nyears,        forcing%netcdf_ins%id_dim_time)
    CALL inquire_dim( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_dim_month,    int_dummy,                 forcing%netcdf_ins%id_dim_month)  
    CALL inquire_dim( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_dim_lat,      forcing%ins_nlat,          forcing%netcdf_ins%id_dim_lat)
    
    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_time,  (/ forcing%netcdf_ins%id_dim_time                                                                 /), forcing%netcdf_ins%id_var_time)
    CALL inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_month, (/ forcing%netcdf_ins%id_dim_month                                                                /), forcing%netcdf_ins%id_var_month)
    CALL inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_lat,   (/ forcing%netcdf_ins%id_dim_lat                                                                  /), forcing%netcdf_ins%id_var_lat)
    CALL inquire_double_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%name_var_Q_TOA, (/ forcing%netcdf_ins%id_dim_time, forcing%netcdf_ins%id_dim_month, forcing%netcdf_ins%id_dim_lat /), forcing%netcdf_ins%id_var_Q_TOA)
        
    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf_ins%ncid)
    
  END SUBROUTINE inquire_insolation_data_file
  SUBROUTINE read_insolation_data_file( forcing, ti0, ti1, ins_Q_TOA0, ins_Q_TOA1) 
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    INTEGER,                        INTENT(IN)    :: ti0, ti1
    REAL(dp), DIMENSION(:,:),       INTENT(OUT)   :: ins_Q_TOA0, ins_Q_TOA1
    
    ! Local variables:
    INTEGER                                       :: mi, li
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Q_temp0, Q_temp1
    
    IF (.NOT. par%master) RETURN
    
    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE( Q_temp0(1, 12, forcing%ins_nlat))
    ALLOCATE( Q_temp1(1, 12, forcing%ins_nlat))
        
    ! Read data
    CALL open_netcdf_file(forcing%netcdf_ins%filename, forcing%netcdf_ins%ncid)
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_Q_TOA, Q_temp0, start = (/ ti0, 1, 1 /), count = (/ 1, 12, forcing%ins_nlat /), stride = (/ 1, 1, 1 /) ))    
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_Q_TOA, Q_temp1, start = (/ ti1, 1, 1 /), count = (/ 1, 12, forcing%ins_nlat /), stride = (/ 1, 1, 1 /) ))
    CALL close_netcdf_file(forcing%netcdf_ins%ncid) 
    
    ! Store the data in the shared memory structure
    DO mi = 1, 12
    DO li = 1, forcing%ins_nlat    
      ins_Q_TOA0( li,mi) = Q_temp0( 1,mi,li)
      ins_Q_TOA1( li,mi) = Q_temp1( 1,mi,li)
    END DO
    END DO
        
    ! Clean up temporary memory
    DEALLOCATE(Q_temp0)
    DEALLOCATE(Q_temp1)
   
  END SUBROUTINE read_insolation_data_file
  SUBROUTINE read_insolation_data_file_time_lat( forcing) 
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
    
    IF (.NOT. par%master) RETURN
    
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_ins%filename, forcing%netcdf_ins%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_time,    forcing%ins_time,    start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_lat,     forcing%ins_lat,     start = (/ 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf_ins%ncid)
    
  END SUBROUTINE read_insolation_data_file_time_lat
  
  ! Geothermal heat flux
  SUBROUTINE inquire_geothermal_heat_flux_file( forcing)
    IMPLICIT NONE

    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
    
    IF (.NOT. par%master) RETURN

    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_ghf%filename, forcing%netcdf_ghf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_dim_lon,      forcing%ghf_nlon,                 forcing%netcdf_ghf%id_dim_lon)
    CALL inquire_dim( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_dim_lat,      forcing%ghf_nlat,                 forcing%netcdf_ghf%id_dim_lat)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_var_lon, (/ forcing%netcdf_ghf%id_dim_lon                                /), forcing%netcdf_ghf%id_var_lon)
    CALL inquire_double_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_var_lat, (/ forcing%netcdf_ghf%id_dim_lat                                /), forcing%netcdf_ghf%id_var_lat)
    CALL inquire_double_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_var_ghf, (/ forcing%netcdf_ghf%id_dim_lon, forcing%netcdf_ghf%id_dim_lat /), forcing%netcdf_ghf%id_var_ghf)

    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf_ghf%ncid)

  END SUBROUTINE inquire_geothermal_heat_flux_file
  SUBROUTINE read_geothermal_heat_flux_file( forcing)
  
    USE parameters_module, ONLY: sec_per_year
  
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    
    IF (.NOT. par%master) RETURN

    ! Read data
    CALL open_netcdf_file(forcing%netcdf_ghf%filename, forcing%netcdf_ghf%ncid)
    CALL handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_lon, forcing%ghf_lon, start=(/1   /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_lat, forcing%ghf_lat, start=(/1   /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_ghf, forcing%ghf_ghf, start=(/1, 1/), count=(/forcing%ghf_nlon, forcing%ghf_nlat/), stride=(/1, 1/) ))
    CALL close_netcdf_file(forcing%netcdf_ghf%ncid)

    ! Convert from W m-2 (J m-2 s-1) to J m-2 yr-1
    forcing%ghf_ghf = forcing%ghf_ghf * sec_per_year

  END SUBROUTINE read_geothermal_heat_flux_file

  ! Climate forcing
  SUBROUTINE inquire_climate_SMB_forcing_data_file( forcing)
  
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
 
    ! Local variables:
    INTEGER                                :: lon,lat,x,y,t,m
    INTEGER                                :: int_dummy
    
    IF (.NOT. par%master) RETURN
            
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_clim%filename, forcing%netcdf_clim%ncid)
    
    IF     (C%choice_forcing_method == 'climate_direct') THEN
      ! Monthly climate fields are prescribed as forcing
      
      IF     (C%domain_climate_forcing == 'global') THEN
        ! Forcing is provided on a global lon/lat-grid
        
        ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_lon,   forcing%clim_nlon,   forcing%netcdf_clim%id_dim_lon  )
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_lat,   forcing%clim_nlat,   forcing%netcdf_clim%id_dim_lat  )
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_month, int_dummy,           forcing%netcdf_clim%id_dim_month)
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_time,  forcing%clim_nyears, forcing%netcdf_clim%id_dim_time )
        
        ! Abbreviate dimension ID's for more readable code
        lon = forcing%netcdf_clim%id_dim_lon
        lat = forcing%netcdf_clim%id_dim_lat
        t   = forcing%netcdf_clim%id_dim_time
        m   = forcing%netcdf_clim%id_dim_month
        
        ! Inquire variable id's. Make sure that each variable has the correct dimensions.
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_lon,    (/ lon            /), forcing%netcdf_clim%id_var_lon   )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_lat,    (/      lat       /), forcing%netcdf_clim%id_var_lat   )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_month,  (/           m    /), forcing%netcdf_clim%id_var_month )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_time,   (/              t /), forcing%netcdf_clim%id_var_time  )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_T2m,    (/ lon, lat, m, t /), forcing%netcdf_clim%id_var_T2m   )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_Precip, (/ lon, lat, m, t /), forcing%netcdf_clim%id_var_Precip)

      ELSEIF (C%domain_climate_forcing == 'regional') THEN
        ! Forcing is provided on a regional x/y-grid
        
        ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_x,     forcing%clim_nx,     forcing%netcdf_clim%id_dim_x    )
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_y,     forcing%clim_ny,     forcing%netcdf_clim%id_dim_y    )
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_month, int_dummy,           forcing%netcdf_clim%id_dim_month)
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_time,  forcing%clim_nyears, forcing%netcdf_clim%id_dim_time )
        
        ! Abbreviate dimension ID's for more readable code
        x   = forcing%netcdf_clim%id_dim_x
        y   = forcing%netcdf_clim%id_dim_y
        t   = forcing%netcdf_clim%id_dim_time
        m   = forcing%netcdf_clim%id_dim_month
        
        ! Inquire variable id's. Make sure that each variable has the correct dimensions.
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_x,      (/ x              /), forcing%netcdf_clim%id_var_x     )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_y,      (/      y         /), forcing%netcdf_clim%id_var_y     )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_month,  (/           m    /), forcing%netcdf_clim%id_var_month )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_time,   (/              t /), forcing%netcdf_clim%id_var_time  )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_T2m,    (/ x,   y,   m, t /), forcing%netcdf_clim%id_var_T2m   )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_Precip, (/ x,   y,   m, t /), forcing%netcdf_clim%id_var_Precip)
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: domain_climate_forcing "', TRIM(C%choice_benchmark_experiment), '" not implemented in inquire_climate_SMB_forcing_data_file!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSEIF (C%choice_forcing_method == 'SMB_direct') THEN
      ! Yearly SMB and surface temperature are prescribed as forcing
      
      IF     (C%domain_climate_forcing == 'global') THEN
        ! Forcing is provided on a global lon/lat-grid
        
        ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_lon,   forcing%clim_nlon,   forcing%netcdf_clim%id_dim_lon  )
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_lat,   forcing%clim_nlat,   forcing%netcdf_clim%id_dim_lat  )
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_time,  forcing%clim_nyears, forcing%netcdf_clim%id_dim_time )
        
        ! Abbreviate dimension ID's for more readable code
        lon = forcing%netcdf_clim%id_dim_lon
        lat = forcing%netcdf_clim%id_dim_lat
        t   = forcing%netcdf_clim%id_dim_time
        
        ! Inquire variable id's. Make sure that each variable has the correct dimensions.
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_lon,    (/ lon         /), forcing%netcdf_clim%id_var_lon   )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_lat,    (/      lat    /), forcing%netcdf_clim%id_var_lat   )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_time,   (/           t /), forcing%netcdf_clim%id_var_time  )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_SMB,    (/ lon, lat, t /), forcing%netcdf_clim%id_var_SMB   )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_T2m,    (/ lon, lat, t /), forcing%netcdf_clim%id_var_T2m   )

      ELSEIF (C%domain_climate_forcing == 'regional') THEN
        ! Forcing is provided on a regional x/y-grid
        
        ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_x,     forcing%clim_nx,     forcing%netcdf_clim%id_dim_x    )
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_y,     forcing%clim_ny,     forcing%netcdf_clim%id_dim_y    )
        CALL inquire_dim( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_dim_time,  forcing%clim_nyears, forcing%netcdf_clim%id_dim_time )
        
        ! Abbreviate dimension ID's for more readable code
        x   = forcing%netcdf_clim%id_dim_x
        y   = forcing%netcdf_clim%id_dim_y
        t   = forcing%netcdf_clim%id_dim_time
        
        ! Inquire variable id's. Make sure that each variable has the correct dimensions.
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_x,      (/ x           /), forcing%netcdf_clim%id_var_x     )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_y,      (/      y      /), forcing%netcdf_clim%id_var_y     )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_time,   (/           t /), forcing%netcdf_clim%id_var_time  )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_SMB,    (/ x,   y,   t /), forcing%netcdf_clim%id_var_SMB   )
        CALL inquire_double_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%name_var_T2m,    (/ x,   y,   t /), forcing%netcdf_clim%id_var_T2m   )
        
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: domain_climate_forcing "', TRIM(C%choice_benchmark_experiment), '" not implemented in inquire_climate_SMB_forcing_data_file!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in inquire_climate_SMB_forcing_data_file!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
        
    ! Close the netcdf file
    CALL close_netcdf_file( forcing%netcdf_clim%ncid)
    
  END SUBROUTINE inquire_climate_SMB_forcing_data_file
  SUBROUTINE read_climate_forcing_data_file_SMB( forcing, ti0, ti1)
  
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    INTEGER,                        INTENT(IN)    :: ti0, ti1

    ! Local variables:
    INTEGER                                       :: loni,lati,i,j
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: SMB_temp0, SMB_temp1, T2m_year_temp0, T2m_year_temp1
    
    IF (.NOT. par%master) RETURN


    IF (C%domain_climate_forcing == 'global') THEN
      ! Forcing is provided on a global lon/lat-grid
      
      ! Temporary memory to store the data read from the netCDF file
      ALLOCATE( SMB_temp0(      forcing%clim_nlon, forcing%clim_nlat, 1))
      ALLOCATE( SMB_temp1(      forcing%clim_nlon, forcing%clim_nlat, 1))
      ALLOCATE( T2m_year_temp0( forcing%clim_nlon, forcing%clim_nlat, 1))
      ALLOCATE( T2m_year_temp1( forcing%clim_nlon, forcing%clim_nlat, 1))
        
      ! Read data
      CALL open_netcdf_file(forcing%netcdf_clim%filename, forcing%netcdf_clim%ncid)
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_SMB,    SMB_temp0,      start = (/ 1, 1, ti0 /), count = (/ forcing%clim_nlon, forcing%clim_nlat, 1 /), stride = (/ 1, 1, 1 /) ))
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_SMB,    SMB_temp1,      start = (/ 1, 1, ti1 /), count = (/ forcing%clim_nlon, forcing%clim_nlat, 1 /), stride = (/ 1, 1, 1 /) ))
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_T2m,    T2m_year_temp0, start = (/ 1, 1, ti0 /), count = (/ forcing%clim_nlon, forcing%clim_nlat, 1 /), stride = (/ 1, 1, 1 /) ))
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_T2m,    T2m_year_temp1, start = (/ 1, 1, ti1 /), count = (/ forcing%clim_nlon, forcing%clim_nlat, 1 /), stride = (/ 1, 1, 1 /) ))
      CALL close_netcdf_file(forcing%netcdf_clim%ncid)

      ! Store the data in the shared memory structure
      DO loni = 1, forcing%clim_nlon
      DO lati = 1, forcing%clim_nlat 
        forcing%clim_SMB0(      loni,lati) =      SMB_temp0( loni,lati,1)
        forcing%clim_SMB1(      loni,lati) =      SMB_temp1( loni,lati,1)
        forcing%clim_T2m_year0( loni,lati) = T2m_year_temp0( loni,lati,1)
        forcing%clim_T2m_year1( loni,lati) = T2m_year_temp1( loni,lati,1)
      END DO
      END DO

    ELSEIF (C%domain_climate_forcing == 'regional') THEN
      ! Forcing is provided on a regional x/y-grid 
      
      ! Temporary memory to store the data read from the netCDF file
      ALLOCATE( SMB_temp0(      forcing%clim_nx, forcing%clim_ny, 1))
      ALLOCATE( SMB_temp1(      forcing%clim_nx, forcing%clim_ny, 1))
      ALLOCATE( T2m_year_temp0( forcing%clim_nx, forcing%clim_ny, 1))
      ALLOCATE( T2m_year_temp1( forcing%clim_nx, forcing%clim_ny, 1))
        
      ! Read data
      CALL open_netcdf_file(forcing%netcdf_clim%filename, forcing%netcdf_clim%ncid)
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_SMB,    SMB_temp0,      start = (/ 1, 1, ti0 /), count = (/ forcing%clim_nx, forcing%clim_ny, 1 /), stride = (/ 1, 1, 1 /) ))
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_SMB,    SMB_temp1,      start = (/ 1, 1, ti1 /), count = (/ forcing%clim_nx, forcing%clim_ny, 1 /), stride = (/ 1, 1, 1 /) ))
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_T2m,    T2m_year_temp0, start = (/ 1, 1, ti0 /), count = (/ forcing%clim_nx, forcing%clim_ny, 1 /), stride = (/ 1, 1, 1 /) ))
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_T2m,    T2m_year_temp1, start = (/ 1, 1, ti1 /), count = (/ forcing%clim_nx, forcing%clim_ny, 1 /), stride = (/ 1, 1, 1 /) ))
      CALL close_netcdf_file(forcing%netcdf_clim%ncid)

      ! Store the data in the shared memory structure (and convert from provided [x,y] to model [y,x] order)
      DO i = 1, forcing%clim_nx
      DO j = 1, forcing%clim_ny 
        forcing%clim_SMB0(      j,i) =      SMB_temp0( i,j,1)
        forcing%clim_SMB1(      j,i) =      SMB_temp1( i,j,1)
        forcing%clim_T2m_year0( j,i) = T2m_year_temp0( i,j,1)
        forcing%clim_T2m_year1( j,i) = T2m_year_temp1( i,j,1)
      END DO
      END DO

    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: domain_climate_forcing "', TRIM(C%choice_benchmark_experiment), '" not implemented in read_climate_forcing_data_file_SMB!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
      
    ! Clean up after yourself
    DEALLOCATE( SMB_temp0)
    DEALLOCATE( SMB_temp1)
    DEALLOCATE( T2m_year_temp0)
    DEALLOCATE( T2m_year_temp1)
      
  END SUBROUTINE read_climate_forcing_data_file_SMB
  SUBROUTINE read_climate_forcing_data_file_climate( forcing, ti0, ti1)
  
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    INTEGER,                        INTENT(IN)    :: ti0, ti1

    ! Local variables:
    INTEGER                                       :: loni,lati,i,j,m
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: T2m_temp0, T2m_temp1, Precip_temp0, Precip_temp1
    
    IF (.NOT. par%master) RETURN


    IF (C%domain_climate_forcing == 'global') THEN
      ! Forcing is provided on a global lon/lat-grid
      
      ! Temporary memory to store the data read from the netCDF file
      ALLOCATE(    T2m_temp0( forcing%clim_nlon, forcing%clim_nlat, 12, 1))
      ALLOCATE(    T2m_temp1( forcing%clim_nlon, forcing%clim_nlat, 12, 1))
      ALLOCATE( Precip_temp0( forcing%clim_nlon, forcing%clim_nlat, 12, 1))
      ALLOCATE( Precip_temp1( forcing%clim_nlon, forcing%clim_nlat, 12, 1))
        
      ! Read data
      CALL open_netcdf_file(forcing%netcdf_clim%filename, forcing%netcdf_clim%ncid)
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_T2m,    T2m_temp0,    start = (/ 1, 1, 1, ti0 /), count = (/ forcing%clim_nlon, forcing%clim_nlat, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_T2m,    T2m_temp1,    start = (/ 1, 1, 1, ti1 /), count = (/ forcing%clim_nlon, forcing%clim_nlat, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_Precip, Precip_temp0, start = (/ 1, 1, 1, ti0 /), count = (/ forcing%clim_nlon, forcing%clim_nlat, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_Precip, Precip_temp1, start = (/ 1, 1, 1, ti1 /), count = (/ forcing%clim_nlon, forcing%clim_nlat, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
      CALL close_netcdf_file(forcing%netcdf_clim%ncid)

      ! Store the data in the shared memory structure
      DO m    = 1, 12
      DO loni = 1, forcing%clim_nlon
      DO lati = 1, forcing%clim_nlat 
        forcing%clim_T2m0(    loni,lati,m) =    T2m_temp0( loni,lati,m,1)
        forcing%clim_T2m1(    loni,lati,m) =    T2m_temp1( loni,lati,m,1)
        forcing%clim_Precip0( loni,lati,m) = Precip_temp0( loni,lati,m,1)
        forcing%clim_Precip1( loni,lati,m) = Precip_temp1( loni,lati,m,1)
      END DO
      END DO
      END DO

    ELSEIF (C%domain_climate_forcing == 'regional') THEN
      ! Forcing is provided on a regional x/y-grid 
      
      ! Temporary memory to store the data read from the netCDF file
      ALLOCATE(    T2m_temp0( forcing%clim_nx,forcing%clim_ny, 12, 1))
      ALLOCATE(    T2m_temp1( forcing%clim_nx,forcing%clim_ny, 12, 1))
      ALLOCATE( Precip_temp0( forcing%clim_nx,forcing%clim_ny, 12, 1))
      ALLOCATE( Precip_temp1( forcing%clim_nx,forcing%clim_ny, 12, 1))
        
      ! Read data
      CALL open_netcdf_file(forcing%netcdf_clim%filename, forcing%netcdf_clim%ncid)
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_T2m,    T2m_temp0,    start = (/ 1, 1, 1, ti0 /), count = (/ forcing%clim_nx, forcing%clim_ny, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_T2m,    T2m_temp1,    start = (/ 1, 1, 1, ti1 /), count = (/ forcing%clim_nx, forcing%clim_ny, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_Precip, Precip_temp0, start = (/ 1, 1, 1, ti0 /), count = (/ forcing%clim_nx, forcing%clim_ny, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
      CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_Precip, Precip_temp1, start = (/ 1, 1, 1, ti1 /), count = (/ forcing%clim_nx, forcing%clim_ny, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
      CALL close_netcdf_file(forcing%netcdf_clim%ncid)

      ! Store the data in the shared memory structure (and convert from provided [x,y,m] to model [m,y,x] order)
      DO m = 1, 12
      DO i = 1, forcing%clim_nx
      DO j = 1, forcing%clim_ny 
        forcing%clim_T2m0(    m,j,i) =    T2m_temp0( i,j,m,1)
        forcing%clim_T2m1(    m,j,i) =    T2m_temp1( i,j,m,1)
        forcing%clim_Precip0( m,j,i) = Precip_temp0( i,j,m,1)
       forcing% clim_Precip1( m,j,i) = Precip_temp1( i,j,m,1)
      END DO
      END DO
      END DO

    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: domain_climate_forcing "', TRIM(C%choice_benchmark_experiment), '" not implemented in read_climate_forcing_data_file_climate!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Clean up after yourself
    DEALLOCATE(    T2m_temp0)
    DEALLOCATE(    T2m_temp1)
    DEALLOCATE( Precip_temp0)
    DEALLOCATE( Precip_temp1)
      
  END SUBROUTINE read_climate_forcing_data_file_climate
  SUBROUTINE read_climate_forcing_data_file_time_latlon( forcing)
  
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
    
    IF (.NOT. par%master) RETURN
    
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_clim%filename, forcing%netcdf_clim%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_time,    forcing%clim_time,    start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_lat,     forcing%clim_lat,     start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_lon,     forcing%clim_lon,     start = (/ 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( forcing%netcdf_clim%ncid)
    
  END SUBROUTINE read_climate_forcing_data_file_time_latlon
  SUBROUTINE read_climate_forcing_data_file_time_xy( forcing)
  
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
    
    IF (.NOT. par%master) RETURN
    
    ! Open the netcdf file
    CALL open_netcdf_file( forcing%netcdf_clim%filename, forcing%netcdf_clim%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_time,    forcing%clim_time,    start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_x,       forcing%clim_x,       start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_clim%ncid, forcing%netcdf_clim%id_var_y,       forcing%clim_y,       start = (/ 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( forcing%netcdf_clim%ncid)
    
  END SUBROUTINE read_climate_forcing_data_file_time_xy
  
  ! Global topography for SELEN
  SUBROUTINE inquire_SELEN_global_topo_file( SELEN)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN
    
    ! Local variables:
    INTEGER                                       :: int_dummy
    
    IF (.NOT. par%master) RETURN
    
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
    
  END SUBROUTINE inquire_SELEN_global_topo_file
  SUBROUTINE read_SELEN_global_topo_file( SELEN)
    ! Read the init netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_SELEN_global),        INTENT(INOUT) :: SELEN
    
    IF (.NOT. par%master) RETURN
    
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
    
  END SUBROUTINE read_SELEN_global_topo_file
  
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
    INTEGER                                   :: i,j
    REAL(dp), DIMENSION(nx,ny)                :: d_flip
    
    DO i = 1, nx
    DO j = 1, ny
      d_flip( i,j) = d( j,i)
    END DO
    END DO
    
    CALL handle_error( nf90_put_var( ncid, var_id, d_flip, start=start_vec))
    
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
    INTEGER                                   :: i,j
    REAL(dp), DIMENSION(nx,ny)                :: d_flip
    
    DO i = 1, nx
    DO j = 1, ny
      d_flip( i,j) = d( j,i)
    END DO
    END DO
    
    CALL handle_error( nf90_put_var( ncid, var_id, d_flip, start=start_vec))
    
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
    INTEGER                                      :: i,j,k
    REAL(dp), DIMENSION(nx,ny,nz)                :: d_flip
    
    DO i = 1, nx
    DO j = 1, ny
    DO k = 1, nz
      d_flip( i,j,k) = d( k,j,i)
    END DO
    END DO
    END DO
    
    CALL handle_error( nf90_put_var( ncid, var_id, d_flip, start=start_vec))
    
  END SUBROUTINE write_data_to_file_dp_3D

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
  SUBROUTINE create_int_var( ncid, var_name, id_dims, id_var, long_name, units, missing_value)
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
    IF(xtype /= nf90_int) THEN
     WRITE(0,'(3A)') 'ERROR: Actual type of variable "',var_name,'" is not nf90_int.'
     STOP
    END IF
    IF(ndims /= SIZE(id_dims)) THEN
     WRITE(0,'(A,I5,3A,I5,A)') 'ERROR: Actual number of dimensions(', &
            ndims,') of variable "',var_name,'": does not match required number of dimensions (',SIZE(id_dims),').'
     STOP
    END IF
    IF(ANY(actual_id_dims(1:ndims) /= id_dims)) THEN
     WRITE(0,'(3A)') 'ERROR: Actual dimensions of variable "',var_name,'" does not match required dimensions.'
     STOP
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
    IF(xtype /= nf90_float) THEN
     WRITE(0,'(3A)') 'ERROR: Actual type of variable "',var_name,'" is not nf90_float.'
     STOP
    END IF
    IF(ndims /= SIZE(id_dims)) THEN
     WRITE(0,'(A,I5,3A,I5,A)') 'ERROR: Actual number of dimensions(', &
            ndims,') of variable "',var_name,'": does not match required number of dimensions (',SIZE(id_dims),').'
     STOP
    END IF
    IF(ANY(actual_id_dims(1:ndims) /= id_dims)) THEN
     WRITE(0,'(3A)') 'ERROR: Actual dimensions of variable "',var_name,'" does not match required dimensions.'
     STOP
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
     WRITE(0,'(3A)') 'ERROR: Actual type of variable "',var_name,'" is not nf90_DOUBLE.'
     STOP
    END IF
    IF(ndims /= SIZE(id_dims)) THEN
     WRITE(0,'(A,I5,3A,I5,A)') 'ERROR: Actual number of dimensions(', &
            ndims,') of variable "',var_name,'": does not match required number of dimensions (',SIZE(id_dims),').'
     STOP
    END IF
    IF(ANY(actual_id_dims(1:ndims) /= id_dims)) THEN
     WRITE(0,'(3A)') 'ERROR: Actual dimensions of variable "',var_name,'" does not match required dimensions.'
     STOP
    END IF
    
  END SUBROUTINE inquire_double_var
  SUBROUTINE handle_error( stat, message)
    USE netcdf, ONLY: nf90_noerr, nf90_strerror
    IMPLICIT NONE
    
    ! Input variables:
    INTEGER,                    INTENT(IN) :: stat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: message

    IF(stat /= nf90_noerr) THEN
     IF(PRESENT(message)) THEN
      WRITE(0,'(A,A,A,A)') 'ERROR: ', TRIM(nf90_strerror(stat)), ' concerning: ', message
     ELSE
      WRITE(0,'(A,A)')     'ERROR: ', TRIM(nf90_strerror(stat))
     END IF
     STOP
    END IF
    
  END SUBROUTINE handle_error
    
END MODULE netcdf_module
