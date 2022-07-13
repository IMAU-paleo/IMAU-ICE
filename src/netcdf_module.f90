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
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_model_region, type_reference_geometry, type_restart_data, &
                                             type_forcing_data, type_climate_snapshot_global, type_ocean_snapshot_global, type_debug_fields, &
                                             type_SELEN_global, type_global_scalar_data, type_highres_ocean_data, &
                                             type_direct_climate_forcing_global, type_direct_climate_forcing_regional, &
                                             type_direct_SMB_forcing_global, type_direct_SMB_forcing_regional, type_BIV_target_velocity, &
                                             type_BIV_bed_roughness
  USE data_types_netcdf_module
  
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
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_to_restart_file'
    INTEGER                                       :: ncid, nx, ny, nz, ti
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    WRITE(0,'(A,F9.3,A)') '   t = ', region%time/1E3, ' kyr - writing output...'
    
    ! Open the file for writing
    CALL open_netcdf_file( region%restart%netcdf%filename, region%restart%netcdf%ncid)
        
    ! Time
    CALL handle_error( nf90_put_var( region%restart%netcdf%ncid, region%restart%netcdf%id_var_time, region%time, start = (/ region%restart%netcdf%ti /)))
    
    ! Placeholders for the dimension ID's, for shorter code
    ncid = region%restart%netcdf%ncid
    nx   = region%grid%nx
    ny   = region%grid%ny
    nz   = C%nZ
    ti   = region%restart%netcdf%ti
    
    ! Geometry
    CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%netcdf%id_var_Hi,               region%ice%Hi_a,             (/1, 1,    ti/))
    CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%netcdf%id_var_Hb,               region%ice%Hb_a,             (/1, 1,    ti/))
    CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%netcdf%id_var_Hs,               region%ice%Hs_a,             (/1, 1,    ti/))
    CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%netcdf%id_var_SL,               region%ice%SL_a,             (/1, 1,    ti/))
    CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%netcdf%id_var_dHb,              region%ice%dHb_a,            (/1, 1,    ti/))
    
    ! Temperature
    CALL write_data_to_file_dp_3D( ncid, nx, ny, nz, region%restart%netcdf%id_var_Ti,               region%ice%Ti_a,             (/1, 1, 1, ti/))
    
    ! SMB
    IF     (C%choice_SMB_model == 'uniform') THEN
    ELSEIF (C%choice_SMB_model == 'idealised') THEN
    ELSEIF (C%choice_SMB_model == 'IMAU-ITM') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, region%restart%netcdf%id_var_FirnDepth,        region%SMB%FirnDepth,        (/1, 1, 1, ti/))
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%netcdf%id_var_MeltPreviousYear, region%SMB%MeltPreviousYear, (/1, 1,    ti/))
    ELSEIF (C%choice_SMB_model == 'direct_global') THEN
    ELSEIF (C%choice_SMB_model == 'direct_regional') THEN
    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM(C%choice_SMB_model) // '"!')
    END IF
    
    ! Isotopes
    IF     (C%choice_ice_isotopes_model == 'none') THEN
    ELSEIF (C%choice_ice_isotopes_model == 'uniform') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%netcdf%id_var_IsoIce,           region%ice%IsoIce,           (/1, 1,    ti/))
    ELSEIF (C%choice_ice_isotopes_model == 'ANICE_legacy') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     region%restart%netcdf%id_var_IsoIce,           region%ice%IsoIce,           (/1, 1,    ti/))
    ELSE
      CALL crash('unknown choice_ice_isotopes_model "' // TRIM(C%choice_ice_isotopes_model) // '"!')
    END IF
    
    ! Inverse routine
    IF     (C%choice_forcing_method == 'none' .OR. &
            C%choice_forcing_method == 'CO2_direct') THEN
      ! No inverse routine used in these forcing methods
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      ! Need to write dT_glob_history and dT_glob_inverse_history
      CALL handle_error( nf90_put_var( ncid, region%restart%netcdf%id_var_dT_glob_history,         forcing%dT_glob_history,         start = (/1, ti/) ))
      CALL handle_error( nf90_put_var( ncid, region%restart%netcdf%id_var_dT_glob_inverse_history, forcing%dT_glob_inverse_history, start = (/1, ti/) ))
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! Need to write dT_glob_history and CO2_inverse_history
      CALL handle_error( nf90_put_var( ncid, region%restart%netcdf%id_var_dT_glob_history,         forcing%dT_glob_history,         start = (/1, ti/) ))
      CALL handle_error( nf90_put_var( ncid, region%restart%netcdf%id_var_CO2_inverse_history,     forcing%CO2_inverse_history,     start = (/1, ti/) ))
    ELSE
      CALL crash('unknown choice_forcing_method "' // TRIM(C%choice_forcing_method) // '"!')
    END IF
    
    ! Close the file
    CALL close_netcdf_file(region%restart%netcdf%ncid)
    
    ! Increase time frame counter
    region%restart%netcdf%ti = region%restart%netcdf%ti + 1
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
        
  END SUBROUTINE write_to_restart_file
  SUBROUTINE write_to_help_fields_file( region)
    ! Write the current model state to the existing output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region), INTENT(INOUT) :: region
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_to_help_fields_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
        
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
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_help_field'
    INTEGER                                       :: ncid, nx, ny, ti, nz, nzo, i, j
    REAL(dp), DIMENSION( region%grid%ny, region%grid%nx) :: d_temp
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Placeholders for the dimension ID's, for shorter code
    ncid = region%help_fields%ncid
    nx   = region%grid%nx
    ny   = region%grid%ny
    nz   = C%nz
    ti   = region%help_fields%ti
    nzo  = C%nz_ocean
    
    IF (field_name == 'none') THEN
      
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
    
    ! Ice basins
    ELSEIF (field_name == 'basin_ID') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,    id_var,               region%ice%basin_ID,       (/1, 1        /))
      
    ! Basal inversion target velocity
    ELSEIF (field_name == 'BIV_target_velocity') THEN
      CALL write_data_to_file_dp_2D( ncid,  nx, ny,    id_var,               region%ice%BIV_uabs_surf_target, (/1, 1        /))

    ! Forcing climates
    ELSEIF (field_name == 'GCM_Warm_T2m') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate_matrix%GCM_warm%T2m,    (/1, 1, 1/))
    ELSEIF (field_name == 'GCM_Warm_Precip') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate_matrix%GCM_warm%Precip, (/1, 1, 1/))
    ELSEIF (field_name == 'GCM_Cold_T2m') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate_matrix%GCM_cold%T2m,    (/1, 1, 1/))
    ELSEIF (field_name == 'GCM_Cold_Precip') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate_matrix%GCM_cold%Precip, (/1, 1, 1/))
    ELSEIF (field_name == 'GCM_PI_T2m') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate_matrix%GCM_PI%T2m,      (/1, 1, 1/))
    ELSEIF (field_name == 'GCM_PI_Precip') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate_matrix%GCM_PI%Precip,   (/1, 1, 1 /))
    ELSEIF (field_name == 'PD_obs_T2m') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate_matrix%PD_obs%T2m,      (/1, 1, 1/))
    ELSEIF (field_name == 'PD_obs_Precip') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate_matrix%PD_obs%Precip,   (/1, 1, 1/))
    
    ! Forcing ocean data  
    ELSEIF (field_name == 'GCM_Warm_T_ocean_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nzo, id_var,              region%ocean_matrix%GCM_Warm%T_ocean_ext, (/1, 1, 1 /))
    ELSEIF (field_name == 'GCM_Warm_S_ocean_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nzo, id_var,              region%ocean_matrix%GCM_Warm%S_ocean_ext, (/1, 1, 1 /))
    ELSEIF (field_name == 'GCM_Cold_T_ocean_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nzo, id_var,              region%ocean_matrix%GCM_Cold%T_ocean_ext, (/1, 1, 1 /))
    ELSEIF (field_name == 'GCM_Cold_S_ocean_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nzo, id_var,              region%ocean_matrix%GCM_Cold%S_ocean_ext, (/1, 1, 1 /))
    ELSEIF (field_name == 'GCM_PI_T_ocean_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nzo, id_var,              region%ocean_matrix%GCM_PI%T_ocean_ext,   (/1, 1, 1 /))
    ELSEIF (field_name == 'GCM_PI_S_ocean_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nzo, id_var,              region%ocean_matrix%GCM_PI%S_ocean_ext,   (/1, 1, 1 /))
    ELSEIF (field_name == 'PD_obs_T_ocean_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nzo, id_var,              region%ocean_matrix%PD_obs%T_ocean_ext,   (/1, 1, 1 /))
    ELSEIF (field_name == 'PD_obs_S_ocean_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nzo, id_var,              region%ocean_matrix%PD_obs%S_ocean_ext,   (/1, 1, 1 /)) 
      
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
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%A_flow_vav_a,   (/1, 1,    ti /))
    ELSEIF (field_name == 'Ti_base_rel') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               (region%ice%Ti_a( C%nz,:,:) - region%ice%Ti_pmp_a( C%nz,:,:)),           (/1, 1, ti /))
      
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
    ELSEIF (field_name == 'R_shear') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%R_shear,        (/1, 1,    ti /))
      
    ! Climate
    ELSEIF (field_name == 'T2m') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate_matrix%applied%T2m, (/1, 1, 1, ti /))
    ELSEIF (field_name == 'T2m_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%climate_matrix%applied%T2m( :,j,i)) / 12._dp
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))
    ELSEIF (field_name == 'Precip') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate_matrix%applied%Precip, (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Precip_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%climate_matrix%applied%Precip( :,j,i))
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))
    ELSEIF (field_name == 'Wind_WE') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate_matrix%applied%Wind_WE, (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Wind_WE_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%climate_matrix%applied%Wind_WE( :,j,i)) / 12._dp
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))
    ELSEIF (field_name == 'Wind_SN') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%climate_matrix%applied%Wind_SN, (/1, 1, 1, ti /))
    ELSEIF (field_name == 'Wind_SN_year') THEN
      DO i = 1, region%grid%nx
      DO j = 1, region%grid%ny
        d_temp( j,i) = SUM( region%climate_matrix%applied%Wind_SN( :,j,i)) / 12._dp
      END DO
      END DO
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               d_temp,                     (/1, 1,    ti /))    
      
    ! Mass balance
    ELSEIF (field_name == 'SMB') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, 12, id_var,               region%SMB%SMB,       (/1, 1, 1, ti /))
    ELSEIF (field_name == 'SMB_year') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%SMB%SMB_year,  (/1, 1,    ti /))
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
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_a,         (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_land') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_land_a,    (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_ocean') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_ocean_a,   (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_lake') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_lake_a,    (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_ice') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_ice_a,     (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_sheet') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_sheet_a,   (/1, 1,    ti /))
    ELSEIF (field_name == 'mask_shelf') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,               region%ice%mask_shelf_a,   (/1, 1,    ti /))
      
    ! Basal conditions
    ELSEIF (field_name == 'pore_water_pressure') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%pore_water_pressure_a, (/1, 1,    ti /))
    ELSEIF (field_name == 'effective_pressure' .OR. field_name == 'Neff') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%Neff_a,           (/1, 1,    ti /))
    ELSEIF (field_name == 'phi_fric') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%phi_fric_a,       (/1, 1,    ti /))
    ELSEIF (field_name == 'tau_yield') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%tauc_a,           (/1, 1,    ti /))
    ELSEIF (field_name == 'alpha_sq') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%alpha_sq_a,       (/1, 1,    ti /))
    ELSEIF (field_name == 'beta_sq') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%beta_sq_a,        (/1, 1,    ti /))
      
    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%IsoIce,            (/1, 1,    ti /))
    ELSEIF (field_name == 'iso_surf') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%IsoSurf,           (/1, 1,    ti /))
    
    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%ice%dHb_a, (/1, 1,    ti /))
      
    ! Oceans and basal melt
    ELSEIF (field_name == 'BMB') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%BMB%BMB,       (/1, 1,   ti /))
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%BMB%BMB_sheet, (/1, 1,   ti /))
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL write_data_to_file_dp_2D( ncid, nx, ny,     id_var,               region%BMB%BMB_shelf, (/1, 1,   ti /))
    ELSEIF (field_name == 'T_ocean_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nzo, id_var,              region%ocean_matrix%applied%T_ocean_corr_ext,  (/1, 1, 1, ti /))
    ELSEIF (field_name == 'S_ocean_3D') THEN
      CALL write_data_to_file_dp_3D( ncid, nx, ny, nzo, id_var,              region%ocean_matrix%applied%S_ocean_corr_ext,  (/1, 1, 1, ti /))   
    ELSEIF (field_name == 'PICO_boxes') THEN
      CALL write_data_to_file_int_2D( ncid, nx, ny,     id_var,              region%BMB%PICO_k, (/ 1, 1,   ti /))
      
    
    ELSE
      CALL crash('unknown help field "' // TRIM(field_name) // '"!')
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE write_help_field
  SUBROUTINE write_to_debug_file
    ! Write the current set of debug data fields to the debug NetCDF file
   
    IMPLICIT NONE
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_to_debug_file'
    INTEGER                                       :: ncid, nx, ny
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    IF (.NOT. C%do_write_debug_data) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
        
  END SUBROUTINE write_to_debug_file
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
  SUBROUTINE write_to_global_scalar_output_file( global_data, time)
    ! Write data to the global scalar output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_global_scalar_data),   INTENT(INOUT) :: global_data
    REAL(dp),                        INTENT(IN)    :: time
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_to_global_scalar_output_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    IF (.NOT. C%do_write_global_scalar_output) RETURN
    
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
    IF     (C%choice_forcing_method == 'none') THEN
    ELSEIF (C%choice_forcing_method == 'CO2_direct') THEN
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_CO2_obs,          global_data%CO2_obs,        start = (/global_data%netcdf%ti/)))
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_CO2_obs,          global_data%CO2_obs,        start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_CO2_mod,          global_data%CO2_mod,        start = (/global_data%netcdf%ti/)))
    ELSE
      CALL crash('unknown choice_forcing_method "' // TRIM(C%choice_forcing_method) // '"!')
    END IF
    
    ! d18O
    IF     (C%do_calculate_benthic_d18O) THEN
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_dT_glob,          global_data%dT_glob,        start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_dT_dw,            global_data%dT_dw,          start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_mod,         global_data%d18O_mod,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ice,         global_data%d18O_ice,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_Tdw,         global_data%d18O_Tdw,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_NAM,         global_data%d18O_NAM,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_EAS,         global_data%d18O_EAS,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_GRL,         global_data%d18O_GRL,       start = (/global_data%netcdf%ti/)))
      CALL handle_error( nf90_put_var( global_data%netcdf%ncid, global_data%netcdf%id_var_d18O_ANT,         global_data%d18O_ANT,       start = (/global_data%netcdf%ti/)))
    END IF
    
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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
        
  END SUBROUTINE write_to_global_scalar_output_file
  SUBROUTINE write_to_regional_scalar_output_file( region, time)
    ! Write data to the regional scalar output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),         INTENT(INOUT) :: region
    REAL(dp),                        INTENT(IN)    :: time
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'write_to_regional_scalar_output_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    IF (.NOT. C%do_write_regional_scalar_output) RETURN
    
    ! Open the file for writing
    CALL open_netcdf_file( region%scalars%filename, region%scalars%ncid)
        
    ! Time
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_time,             time,                       start = (/region%scalars%ti/)))
    
    ! Variables
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_ice_volume,    region%ice_volume,                 start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_ice_volume_af, region%ice_volume_above_flotation, start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_ice_area,      region%ice_area,                   start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_T2m,           region%int_T2m,                    start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_SMB,           region%int_SMB,                    start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_BMB,           region%int_BMB,                    start = (/region%scalars%ti/)))
    CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_MB,            region%int_MB,                     start = (/region%scalars%ti/)))
    
    ! Individual SMB components
    IF     (C%choice_SMB_model == 'uniform' .OR. &
            C%choice_SMB_model == 'idealised' .OR. &
            C%choice_SMB_model == 'direct_global' .OR. &
            C%choice_SMB_model == 'direct_regional') THEN
      ! Do nothing
    ELSEIF (C%choice_SMB_model == 'IMAU-ITM' .OR. &
            C%choice_SMB_model == 'IMAU-ITM_wrongrefreezing') THEN
      CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_snowfall,      region%int_snowfall,               start = (/region%scalars%ti/)))
      CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_rainfall,      region%int_rainfall,               start = (/region%scalars%ti/)))
      CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_melt,          region%int_melt,                   start = (/region%scalars%ti/)))
      CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_refreezing,    region%int_refreezing,             start = (/region%scalars%ti/)))
      CALL handle_error( nf90_put_var( region%scalars%ncid, region%scalars%id_var_runoff,        region%int_runoff,                 start = (/region%scalars%ti/)))
    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM(C%choice_SMB_model) // '"!')
    END IF
    
    ! Close the file
    CALL close_netcdf_file(region%scalars%ncid)
    
    ! Increase time frame counter
    region%scalars%ti = region%scalars%ti + 1
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
        
  END SUBROUTINE write_to_regional_scalar_output_file
  
  ! Create output netcdf files
  SUBROUTINE create_restart_file( region, forcing)
    ! Create a new restart file, containing the key model output required to start another run.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_forcing_data),        INTENT(IN)    :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_restart_file'
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x, y, z, m, t
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Set time frame index to 1
    region%restart%netcdf%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM( region%restart%netcdf%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( region%restart%netcdf%filename) // '" already exists!')
    END IF
    
    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( region%restart%netcdf%filename)
    CALL handle_error(nf90_create( region%restart%netcdf%filename,IOR(nf90_clobber,nf90_share), region%restart%netcdf%ncid))
        
    ! Define dimensions:
    CALL create_dim( region%restart%netcdf%ncid, region%restart%netcdf%name_dim_x,         region%grid%nx,     region%restart%netcdf%id_dim_x         )
    CALL create_dim( region%restart%netcdf%ncid, region%restart%netcdf%name_dim_y,         region%grid%ny,     region%restart%netcdf%id_dim_y         )
    CALL create_dim( region%restart%netcdf%ncid, region%restart%netcdf%name_dim_zeta,      C%nZ,               region%restart%netcdf%id_dim_zeta      ) ! Scaled vertical coordinate
    CALL create_dim( region%restart%netcdf%ncid, region%restart%netcdf%name_dim_month,     12,                 region%restart%netcdf%id_dim_month     ) ! Months (for monthly data)
    CALL create_dim( region%restart%netcdf%ncid, region%restart%netcdf%name_dim_time,      nf90_unlimited,     region%restart%netcdf%id_dim_time      ) ! Time frames
    
    ! Placeholders for the dimension ID's, for shorter code
    x = region%restart%netcdf%id_dim_x
    y = region%restart%netcdf%id_dim_y
    z = region%restart%netcdf%id_dim_zeta
    m = region%restart%netcdf%id_dim_month
    t = region%restart%netcdf%id_dim_time
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.
    
    ! Dimension variables: zeta, month, time
    CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_x,                [x            ], region%restart%netcdf%id_var_x,                long_name='X-coordinate', units='m')
    CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_y,                [   y         ], region%restart%netcdf%id_var_y,                long_name='Y-coordinate', units='m')
    CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_zeta,             [      z      ], region%restart%netcdf%id_var_zeta,             long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_month,            [         m   ], region%restart%netcdf%id_var_month,            long_name='Month', units='1-12'    )
    CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_time,             [            t], region%restart%netcdf%id_var_time,             long_name='Time', units='years'   )
    
  ! ==== Create fields for the different model components =====
  ! ===========================================================
    
    ! Geometry
    CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_Hi,               [x, y,       t], region%restart%netcdf%id_var_Hi,               long_name='Ice thickness', units='m')
    CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_Hb,               [x, y,       t], region%restart%netcdf%id_var_Hb,               long_name='Bedrock elevation', units='m')
    CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_Hs,               [x, y,       t], region%restart%netcdf%id_var_Hs,               long_name='Surface elevation', units='m')
    CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_SL,               [x, y,       t], region%restart%netcdf%id_var_SL,               long_name='Sea surface change',  units='m')
    CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_dHb,              [x, y,       t], region%restart%netcdf%id_var_dHb,              long_name='Bedrock deformation', units='m')
    
    ! Temperature
    CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_Ti,               [x, y, z,    t], region%restart%netcdf%id_var_Ti,               long_name='Ice temperature', units='K')
    
    ! SMB
    IF     (C%choice_SMB_model == 'uniform') THEN
    ELSEIF (C%choice_SMB_model == 'idealised') THEN
    ELSEIF (C%choice_SMB_model == 'IMAU-ITM') THEN
      CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_FirnDepth,        [x, y,    m, t], region%restart%netcdf%id_var_FirnDepth,        long_name='Firn depth', units='m')
      CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_MeltPreviousYear, [x, y,       t], region%restart%netcdf%id_var_MeltPreviousYear, long_name='Melt during previous year', units='mie')
    ELSEIF (C%choice_SMB_model == 'direct_global') THEN
    ELSEIF (C%choice_SMB_model == 'direct_regional') THEN
    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM(C%choice_SMB_model) // '"!')
    END IF
    
    ! Isotopes
    IF     (C%choice_ice_isotopes_model == 'none') THEN
    ELSEIF (C%choice_ice_isotopes_model == 'uniform') THEN
      CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_IsoIce,           [x, y,       t], region%restart%netcdf%id_var_IsoIce,           long_name='Vertically averaged 18O content', units='per mille')
    ELSEIF (C%choice_ice_isotopes_model == 'ANICE_legacy') THEN
      CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_IsoIce,           [x, y,       t], region%restart%netcdf%id_var_IsoIce,           long_name='Vertically averaged 18O content', units='per mille')
    ELSE
      CALL crash('unknown choice_ice_isotopes_model "' // TRIM(C%choice_ice_isotopes_model) // '"!')
    END IF
    
    ! Inverse routine data
    IF     (C%choice_forcing_method == 'none' .OR. &
            C%choice_forcing_method == 'CO2_direct') THEN
      ! No inverse routine used in these forcing methods
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      ! Need to write dT_glob_history and dT_glob_inverse_history
      
      CALL create_dim( region%restart%netcdf%ncid, region%restart%netcdf%name_dim_ndT_glob_history,         forcing%ndT_glob_history,         region%restart%netcdf%id_dim_ndT_glob_history        )
      CALL create_dim( region%restart%netcdf%ncid, region%restart%netcdf%name_dim_ndT_glob_inverse_history, forcing%ndT_glob_inverse_history, region%restart%netcdf%id_dim_ndT_glob_inverse_history)
      
      CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_dT_glob_history,         [region%restart%netcdf%id_dim_ndT_glob_history,         t], region%restart%netcdf%id_var_dT_glob_history,         long_name='dT_glob history')
      CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_dT_glob_inverse_history, [region%restart%netcdf%id_dim_ndT_glob_inverse_history, t], region%restart%netcdf%id_var_dT_glob_inverse_history, long_name='dT_glob_inverse history')
      
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      ! Need to write dT_glob_history and CO2_inverse_history
      
      CALL create_dim( region%restart%netcdf%ncid, region%restart%netcdf%name_dim_ndT_glob_history,         forcing%ndT_glob_history,         region%restart%netcdf%id_dim_ndT_glob_history    )
      CALL create_dim( region%restart%netcdf%ncid, region%restart%netcdf%name_dim_nCO2_inverse_history,     forcing%nCO2_inverse_history,     region%restart%netcdf%id_dim_nCO2_inverse_history)
      
      CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_dT_glob_history,         [region%restart%netcdf%id_dim_ndT_glob_history,         t], region%restart%netcdf%id_var_dT_glob_history,     long_name='dT_glob history')
      CALL create_double_var( region%restart%netcdf%ncid, region%restart%netcdf%name_var_CO2_inverse_history,     [region%restart%netcdf%id_dim_nCO2_inverse_history,     t], region%restart%netcdf%id_var_CO2_inverse_history, long_name='CO2_inverse history')
      
    ELSE
      CALL crash('unknown choice_forcing_method "' // TRIM(C%choice_forcing_method) // '"!')
    END IF
    
  ! ===== End of fields definition =====
  ! ====================================
    
    ! Leave definition mode:
    CALL handle_error(nf90_enddef( region%restart%netcdf%ncid))
    
    ! Write the x, y, zeta and months variable data
    CALL handle_error( nf90_put_var( region%restart%netcdf%ncid, region%restart%netcdf%id_var_x,        region%grid%x                            ))
    CALL handle_error( nf90_put_var( region%restart%netcdf%ncid, region%restart%netcdf%id_var_y,        region%grid%y                            ))
    CALL handle_error( nf90_put_var( region%restart%netcdf%ncid, region%restart%netcdf%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( region%restart%netcdf%ncid, region%restart%netcdf%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( region%restart%netcdf%ncid))
    
    ! Close the file
    CALL close_netcdf_file( region%restart%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE create_restart_file
  SUBROUTINE create_help_fields_file( region)
    ! Create a new help fields file, containing secondary model output (not needed for a restart, but interesting to look at)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_help_fields_file'
    LOGICAL                                       :: file_exists
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Set time frame index to 1
    region%help_fields%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM( region%help_fields%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( region%help_fields%filename) // '" already exists!')
    END IF
    
    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( region%help_fields%filename)
    CALL handle_error( nf90_create( region%help_fields%filename,IOR(nf90_clobber,nf90_share),region%help_fields%ncid))
        
    ! Define dimensions:
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_x,       region%grid%nx, region%help_fields%id_dim_x      )
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_y,       region%grid%ny, region%help_fields%id_dim_y      )
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_zeta,    C%nZ,           region%help_fields%id_dim_zeta   )
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_month,   12,             region%help_fields%id_dim_month  )
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_time,    nf90_unlimited, region%help_fields%id_dim_time   )
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_z_ocean, C%nz_ocean,     region%help_fields%id_dim_z_ocean)
    
    ! Dimension variables: zeta, month, time
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_x,       [region%help_fields%id_dim_x      ], region%help_fields%id_var_x,       long_name='X-coordinate', units='m')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_y,       [region%help_fields%id_dim_y      ], region%help_fields%id_var_y,       long_name='Y-coordinate', units='m')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_zeta,    [region%help_fields%id_dim_zeta   ], region%help_fields%id_var_zeta,    long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_month,   [region%help_fields%id_dim_month  ], region%help_fields%id_var_month,   long_name='Month', units='1-12'    )
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_time,    [region%help_fields%id_dim_time   ], region%help_fields%id_var_time,    long_name='Time', units='years'   )
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_z_ocean, [region%help_fields%id_dim_z_ocean], region%help_fields%id_var_z_ocean, long_name='Depth in ocean', units='m')
    
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
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_x,        region%grid%x                             ))
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_y,        region%grid%y                             ))
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_zeta,     C%zeta                                    ))
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/) ))
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_z_ocean,  C%z_ocean                                 ))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( region%help_fields%ncid))
    
    ! Close the file
    CALL close_netcdf_file( region%help_fields%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
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
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_help_field'
    INTEGER                                       :: x, y, z, m, t, zo
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Placeholders for the dimension ID's, for shorter code
    x  = region%help_fields%id_dim_x
    y  = region%help_fields%id_dim_y
    z  = region%help_fields%id_dim_zeta
    m  = region%help_fields%id_dim_month
    t  = region%help_fields%id_dim_time
    zo = region%help_fields%id_dim_z_ocean
    
    IF (field_name == 'none') THEN
      
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
      
    ! Ice basins
    ELSEIF (field_name == 'basin_ID') THEN
      CALL create_int_var(    region%help_fields%ncid, 'basin_ID',                 [x, y      ], id_var, long_name='Basin ID')
      
    ! Basal inversion target velocity
    ELSEIF (field_name == 'BIV_target_velocity') THEN
      CALL create_double_var( region%help_fields%ncid, 'BIV_target_velocity',      [x, y      ], id_var, long_name='Basal inversion target velocity', units='m yr^-1')
 
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
    
    ! Forcing oceans
    ELSEIF (field_name == 'GCM_Warm_T_ocean_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'Warm_T_ocean_3D',          [x, y, zo], id_var, long_name='Warm 3-D ocean temperature', units='K')
    ELSEIF (field_name == 'GCM_Warm_S_ocean_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'Warm_S_ocean_3D',          [x, y, zo], id_var, long_name='Warm 3-D ocean salinity', units='PSU')
    ELSEIF (field_name == 'GCM_Cold_T_ocean_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'Cold_T_ocean_3D',          [x, y, zo], id_var, long_name='Cold 3-D ocean temperature', units='K')
    ELSEIF (field_name == 'GCM_Cold_S_ocean_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'Cold_S_ocean_3D',          [x, y, zo], id_var, long_name='Cold 3-D ocean salinity', units='PSU')
    ELSEIF (field_name == 'GCM_PI_T_ocean_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'Ref_PI_T_ocean_3D',        [x, y, zo], id_var, long_name='Ref PI 3-D ocean temperature', units='K')
    ELSEIF (field_name == 'GCM_PI_S_ocean_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'Ref_PI_S_ocean_3D',        [x, y, zo], id_var, long_name='Ref PI 3-D ocean salinity', units='PSU')
    ELSEIF (field_name == 'PD_obs_T_ocean_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'Base_PD_T_ocean_3D',       [x, y, zo], id_var, long_name='Base PD 3-D ocean temperature', units='K')
    ELSEIF (field_name == 'PD_obs_S_ocean_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'Base_PD_S_ocean_3D',       [x, y, zo], id_var, long_name='Base PD 3-D ocean salinity', units='PSU')
      
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
    ELSEIF (field_name == 'A_flow_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'A_flow_3D',                [x, y, z, t], id_var, long_name='Ice flow factor', units='Pa^-3 y^-1')
    ELSEIF (field_name == 'A_flow_vav') THEN
      CALL create_double_var( region%help_fields%ncid, 'A_flow_vav',               [x, y,    t], id_var, long_name='Vertically averaged ice flow factor', units='Pa^-3 y^-1')
    ELSEIF (field_name == 'Ti_base_rel') THEN
      CALL create_double_var( region%help_fields%ncid, 'Ti_base_rel',              [x, y,    t], id_var, long_name='Basal temperature relative to pressure melting point', units='K')
      
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
    ELSEIF (field_name == 'R_shear') THEN
      CALL create_double_var( region%help_fields%ncid, 'R_shear',                  [x, y,    t], id_var, long_name='Shearing ratio; 1 = full shearing, 0 = full sliding')
      
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
      
    ! Basal hydrology and roughness
    ELSEIF (field_name == 'pore_water_pressure') THEN
      CALL create_double_var( region%help_fields%ncid, 'pore_water_pressure',      [x, y,    t], id_var, long_name='pore water pressure', units='Pa')
    ELSEIF (field_name == 'effective_pressure' .OR. field_name == 'Neff') THEN
      CALL create_double_var( region%help_fields%ncid, 'effective_pressure',       [x, y,    t], id_var, long_name='effective basal pressure', units='Pa')
    ELSEIF (field_name == 'phi_fric') THEN
      CALL create_double_var( region%help_fields%ncid, 'phi_fric',                 [x, y,    t], id_var, long_name='till friction angle', units='degrees')
    ELSEIF (field_name == 'tau_yield') THEN
      CALL create_double_var( region%help_fields%ncid, 'tau_yield',                [x, y,    t], id_var, long_name='basal yield stress', units='Pa')
    ELSEIF (field_name == 'alpha_sq') THEN
      CALL create_double_var( region%help_fields%ncid, 'alpha_sq',                 [x, y,    t], id_var, long_name='Coulomb-law friction coefficient', units='unitless')
    ELSEIF (field_name == 'beta_sq') THEN
      CALL create_double_var( region%help_fields%ncid, 'beta_sq',                  [x, y,    t], id_var, long_name='Power-law friction coefficient', units='Pa m^−1/3 yr^1/3')
      
    ! Isotopes
    ELSEIF (field_name == 'iso_ice') THEN
      CALL create_double_var( region%help_fields%ncid, 'iso_ice',                  [x, y,    t], id_var, long_name='Vertically averaged ice d18O', units='per mille')
    ELSEIF (field_name == 'iso_surf') THEN
      CALL create_double_var( region%help_fields%ncid, 'iso_surf',                 [x, y,    t], id_var, long_name='d18O of precipitation', units='per mille')
    
    ! GIA
    ELSEIF (field_name == 'dHb') THEN
      CALL create_double_var( region%help_fields%ncid, 'dHb',                      [x, y,    t], id_var, long_name='Change in bedrock elevation w.r.t. PD', units='m')
    
    ! Oceans and basal melt
    ELSEIF (field_name == 'BMB') THEN
      CALL create_double_var( region%help_fields%ncid, 'BMB',                      [x, y,    t], id_var, long_name='Annual basal mass balance', units='m ice equivalent')
    ELSEIF (field_name == 'BMB_sheet') THEN
      CALL create_double_var( region%help_fields%ncid, 'BMB_sheet',                [x, y,    t], id_var, long_name='Annual basal mass balance for grounded ice', units='m ice equivalent')
    ELSEIF (field_name == 'BMB_shelf') THEN
      CALL create_double_var( region%help_fields%ncid, 'BMB_shelf',                [x, y,    t], id_var, long_name='Annual basal mass balance for floating ice', units='m ice equivalent')
    ELSEIF (field_name == 'T_ocean_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'T_ocean_3D',               [x, y, zo, t], id_var, long_name='3-D ocean temperature', units='K')
    ELSEIF (field_name == 'S_ocean_3D') THEN
      CALL create_double_var( region%help_fields%ncid, 'S_ocean_3D',               [x, y, zo, t], id_var, long_name='3-D ocean salinity', units='PSU')
    ELSEIF (field_name == 'PICO_boxes') THEN
      CALL create_int_var(    region%help_fields%ncid, 'PICO_boxes',               [x, y,    t], id_Var, long_name='PICO ocean boxes')
      
    ELSE
      CALL crash('unknown field name "' // TRIM(field_name) // '"!')
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE create_help_field
  SUBROUTINE create_debug_file( region)
    ! Create the debug NetCDF file; a lot of data fields but no time dimension.
    
    USE data_types_netcdf_module, ONLY: type_netcdf_debug
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_debug_file'
    TYPE(type_netcdf_debug)                       :: debug_temp
    CHARACTER(LEN=20)                             :: short_filename
    INTEGER                                       :: n
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x, y, z, m, ncid
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    IF (.NOT. C%do_write_debug_data) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Determine debug NetCDF filename for this model region
    short_filename = 'debug_NAM.nc'
    short_filename(7:9) = region%name
    DO n = 1, 256
      debug_temp%filename(n:n) = ' '
    END DO
    debug_temp%filename = TRIM(C%output_dir) // TRIM(short_filename)

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(debug_temp%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( debug_temp%filename) // '" already exists!')
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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE create_debug_file
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
  SUBROUTINE create_global_scalar_output_file( netcdf)
    ! Create a new global scalar output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_netcdf_scalars_global), INTENT(INOUT) :: netcdf

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_global_scalar_output_file'
    LOGICAL                                         :: file_exists
    INTEGER                                         :: t
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master .OR. .NOT. C%do_write_global_scalar_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Set time frame index to 1
    netcdf%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    netcdf%filename = TRIM(C%output_dir) // '/scalar_output_global.nc'
    INQUIRE(EXIST=file_exists, FILE = TRIM(netcdf%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( netcdf%filename) // '" already exists!')
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
    IF     (C%choice_forcing_method == 'none') THEN
    ELSEIF (C%choice_forcing_method == 'CO2_direct') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_CO2_obs,       [t], netcdf%id_var_CO2_obs,       long_name='Observed atmospheric CO2 concentration', units='ppm')
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_CO2_obs,       [t], netcdf%id_var_CO2_obs,       long_name='Observed atmospheric CO2 concentration', units='ppm')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_CO2_mod,       [t], netcdf%id_var_CO2_mod,       long_name='Modelled atmospheric CO2 concentration', units='ppm')
    ELSE
      CALL crash('unknown choice_forcing_method "' // TRIM(C%choice_forcing_method) // '"!')
    END IF
    
    ! d18O
    IF     (C%do_calculate_benthic_d18O) THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_dT_glob,       [t], netcdf%id_var_dT_glob,       long_name='Global annual mean surface temperature change', units='K')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_dT_dw,         [t], netcdf%id_var_dT_dw,         long_name='Deep-water temperature change', units='K')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_mod,      [t], netcdf%id_var_d18O_mod,      long_name='Modelled benthic d18O', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_ice,      [t], netcdf%id_var_d18O_ice,      long_name='Modelled benthic d18O from global ice volume', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_Tdw,      [t], netcdf%id_var_d18O_Tdw,      long_name='Modelled benthic d18O from deep-water temperature', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_NAM,      [t], netcdf%id_var_d18O_NAM,      long_name='Modelled benthic d18O from ice in North America', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_EAS,      [t], netcdf%id_var_d18O_EAS,      long_name='Modelled benthic d18O from ice in Eurasia', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_GRL,      [t], netcdf%id_var_d18O_GRL,      long_name='Modelled benthic d18O from ice in Greenland', units='per mil')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_d18O_ANT,      [t], netcdf%id_var_d18O_ANT,      long_name='Modelled benthic d18O from ice in Antarctica', units='per mil')
    END IF
    
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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE create_global_scalar_output_file
  SUBROUTINE create_regional_scalar_output_file( region)
    ! Create a new regional scalar output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),          INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_regional_scalar_output_file'
    LOGICAL                                         :: file_exists
    INTEGER                                         :: t
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master .OR. .NOT. C%do_write_regional_scalar_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Set time frame index to 1
    region%scalars%ti = 1
    
    ! Generate filename
    region%scalars%filename = TRIM(C%output_dir) // '/scalar_output_' // TRIM(regioN%name) // '.nc'

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(region%scalars%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( region%scalars%filename) // '" already exists!')
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
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_SMB,           [t], region%scalars%id_var_SMB,           long_name='Ice-sheet integrated surface mass balance', units='Gigaton yr^-1')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_BMB,           [t], region%scalars%id_var_BMB,           long_name='Ice-sheet integrated basal mass balance', units='Gigaton yr^-1')
    CALL create_double_var( region%scalars%ncid, region%scalars%name_var_MB,            [t], region%scalars%id_var_MB,            long_name='Ice-sheet integrated mass balance', units='Gigaton yr^-1')
    
    ! Individual SMB components
    IF     (C%choice_SMB_model == 'uniform' .OR. &
            C%choice_SMB_model == 'idealised' .OR. &
            C%choice_SMB_model == 'direct_global' .OR. &
            C%choice_SMB_model == 'direct_regional') THEN
      ! Do nothing
    ELSEIF (C%choice_SMB_model == 'IMAU-ITM' .OR. &
            C%choice_SMB_model == 'IMAU-ITM_wrongrefreezing') THEN
      CALL create_double_var( region%scalars%ncid, region%scalars%name_var_snowfall,      [t], region%scalars%id_var_snowfall,      long_name='Ice-sheet integrated snowfall', units='Gigaton yr^-1')
      CALL create_double_var( region%scalars%ncid, region%scalars%name_var_rainfall,      [t], region%scalars%id_var_rainfall,      long_name='Ice-sheet integrated rainfall', units='Gigaton yr^-1')
      CALL create_double_var( region%scalars%ncid, region%scalars%name_var_melt,          [t], region%scalars%id_var_melt,          long_name='Ice-sheet integrated melt', units='Gigaton yr^-1')
      CALL create_double_var( region%scalars%ncid, region%scalars%name_var_refreezing,    [t], region%scalars%id_var_refreezing,    long_name='Ice-sheet integrated refreezing', units='Gigaton yr^-1')
      CALL create_double_var( region%scalars%ncid, region%scalars%name_var_runoff,        [t], region%scalars%id_var_runoff,        long_name='Ice-sheet integrated runoff', units='Gigaton yr^-1')
    ELSE
      CALL crash('unknown choice_SMB_model "' // TRIM(C%choice_SMB_model) // '"!')
    END IF
    
    ! Leave definition mode:
    CALL handle_error(nf90_enddef( region%scalars%ncid))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( region%scalars%ncid))
    
    ! Close the file
    CALL close_netcdf_file( region%scalars%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
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
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'associate_debug_fields'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE associate_debug_fields
  SUBROUTINE initialise_debug_fields( region)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_debug_fields'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF     (region%name == 'NAM') THEN
      CALL initialise_debug_fields_region( debug_NAM, region%grid%nx, region%grid%ny)
    ELSEIF (region%name == 'EAS') THEN
      CALL initialise_debug_fields_region( debug_EAS, region%grid%nx, region%grid%ny)
    ELSEIF (region%name == 'GRL') THEN
      CALL initialise_debug_fields_region( debug_GRL, region%grid%nx, region%grid%ny)
    ELSEIF (region%name == 'ANT') THEN
      CALL initialise_debug_fields_region( debug_ANT, region%grid%nx, region%grid%ny)
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE initialise_debug_fields
  SUBROUTINE initialise_debug_fields_region( debug, nx, ny)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_debug_fields),         INTENT(INOUT)     :: debug
    INTEGER,                         INTENT(IN)        :: nx, ny
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_debug_fields_region'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE initialise_debug_fields_region
  
! Read data to restart a run
! ==========================
  
  SUBROUTINE inquire_restart_file_geometry(    restart) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_restart_data),        INTENT(INOUT) :: restart
 
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_restart_file_geometry'
    INTEGER                                       :: x, y, t
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
        
    ! Open the netcdf file
    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist, return their lengths.
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_x,     restart%nx, restart%netcdf%id_dim_x    )
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_y,     restart%ny, restart%netcdf%id_dim_y    )
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_time,  restart%nt, restart%netcdf%id_dim_time )
    
    ! Abbreviations for shorter code
    x = restart%netcdf%id_dim_x
    y = restart%netcdf%id_dim_y
    t = restart%netcdf%id_dim_time

    ! Inquire variable ID's; make sure that each variable has the correct dimensions.
    
    ! Dimensions
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_x,                (/ x             /), restart%netcdf%id_var_x   )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_y,                (/    y          /), restart%netcdf%id_var_y   )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_time,             (/             t /), restart%netcdf%id_var_time)
    
    ! Data
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_Hi,               (/ x, y   ,    t /), restart%netcdf%id_var_Hi  )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_Hb,               (/ x, y   ,    t /), restart%netcdf%id_var_Hb  )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_Hs,               (/ x, y   ,    t /), restart%netcdf%id_var_Hs  )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_SL,               (/ x, y,       t /), restart%netcdf%id_var_SL  )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_dHb,              (/ x, y,       t /), restart%netcdf%id_var_dHb )

    ! Close the netcdf file
    CALL close_netcdf_file( restart%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_restart_file_geometry
  SUBROUTINE inquire_restart_file_temperature( restart) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_restart_data), INTENT(INOUT)        :: restart
 
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_restart_file_temperature'
    INTEGER                                       :: x, y, z, t
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
        
    ! Open the netcdf file
    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist, return their lengths.
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_x,     restart%nx, restart%netcdf%id_dim_x    )
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_y,     restart%ny, restart%netcdf%id_dim_y    )
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_zeta,  restart%nz, restart%netcdf%id_dim_zeta )
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_time,  restart%nt, restart%netcdf%id_dim_time )
    
    ! Abbreviations for shorter code
    x = restart%netcdf%id_dim_x
    y = restart%netcdf%id_dim_y
    z = restart%netcdf%id_dim_zeta
    t = restart%netcdf%id_dim_time

    ! Inquire variable ID's; make sure that each variable has the correct dimensions.
    
    ! Dimensions
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_x,                (/ x             /), restart%netcdf%id_var_x   )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_y,                (/    y          /), restart%netcdf%id_var_y   )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_zeta,             (/       z       /), restart%netcdf%id_var_zeta)
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_time,             (/             t /), restart%netcdf%id_var_time)
    
    ! Data
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_Ti,               (/ x, y, z,    t /), restart%netcdf%id_var_Ti  )

    ! Close the netcdf file
    CALL close_netcdf_file( restart%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_restart_file_temperature
  SUBROUTINE inquire_restart_file_SMB(         restart) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_restart_data), INTENT(INOUT)        :: restart
 
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_restart_file_SMB'
    INTEGER                                       :: x, y, m, t, int_dummy
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
        
    ! Open the netcdf file
    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist, return their lengths.
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_x,     restart%nx, restart%netcdf%id_dim_x    )
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_y,     restart%ny, restart%netcdf%id_dim_y    )
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_time,  restart%nt, restart%netcdf%id_dim_time )
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_month, int_dummy,  restart%netcdf%id_dim_month)
    
    ! Abbreviations for shorter code
    x = restart%netcdf%id_dim_x
    y = restart%netcdf%id_dim_y
    m = restart%netcdf%id_dim_month
    t = restart%netcdf%id_dim_time

    ! Inquire variable ID's; make sure that each variable has the correct dimensions.
    
    ! Dimensions
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_x,                (/ x             /), restart%netcdf%id_var_x   )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_y,                (/    y          /), restart%netcdf%id_var_y   )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_time,             (/             t /), restart%netcdf%id_var_time)
    
    ! Data
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_FirnDepth,        (/ x, y,    m, t /), restart%netcdf%id_var_FirnDepth)
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_MeltPreviousYear, (/ x, y,       t /), restart%netcdf%id_var_MeltPreviousYear)

    ! Close the netcdf file
    CALL close_netcdf_file( restart%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_restart_file_SMB
  SUBROUTINE inquire_restart_file_isotopes(    restart) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_restart_data), INTENT(INOUT)        :: restart
 
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_restart_file_isotopes'
    INTEGER                                       :: x, y, t
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
        
    ! Open the netcdf file
    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist, return their lengths.
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_x,     restart%nx, restart%netcdf%id_dim_x    )
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_y,     restart%ny, restart%netcdf%id_dim_y    )
    CALL inquire_dim( restart%netcdf%ncid, restart%netcdf%name_dim_time,  restart%nt, restart%netcdf%id_dim_time )
    
    ! Abbreviations for shorter code
    x = restart%netcdf%id_dim_x
    y = restart%netcdf%id_dim_y
    t = restart%netcdf%id_dim_time

    ! Inquire variable ID's; make sure that each variable has the correct dimensions.
    
    ! Dimensions
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_x,                (/ x             /), restart%netcdf%id_var_x   )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_y,                (/    y          /), restart%netcdf%id_var_y   )
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_time,             (/             t /), restart%netcdf%id_var_time)
    
    ! Data
    CALL inquire_double_var( restart%netcdf%ncid, restart%netcdf%name_var_IsoIce,           (/ x, y,       t /), restart%netcdf%id_var_IsoIce )

    ! Close the netcdf file
    CALL close_netcdf_file( restart%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_restart_file_isotopes
  
  SUBROUTINE read_restart_file_geometry(       restart, time_to_restart_from)
    ! Read the restart netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_restart_data),        INTENT(INOUT) :: restart
    REAL(dp),                       INTENT(IN)    :: time_to_restart_from

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_restart_file_geometry'
    INTEGER                                       :: ti, ti_min
    REAL(dp)                                      :: dt, dt_min
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)
    
    ! Read x,y
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_x, restart%x, start=(/1/) ))
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_y, restart%y, start=(/1/) ))
    
    ! Read time, determine which time frame to read
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_time, restart%time, start=(/1/) ))
    
    IF (time_to_restart_from < MINVAL(restart%time) .OR. time_to_restart_from > MAXVAL(restart%time)) THEN
      CALL crash('time_to_restart_from = {dp_01} outside range of restart file!', dp_01 = time_to_restart_from)
    END IF    
    
    ti_min = 0
    dt_min = 1E8_dp
    DO ti = 1, restart%nt
      dt = ABS(restart%time( ti) - time_to_restart_from)
      IF (dt < dt_min) THEN
        ti_min = ti
        dt_min = dt
      END IF
    END DO
    ti = ti_min
    
    IF (dt_min > 0._dp) THEN
      CALL warning('no exact match for time_to_restart_from = {dp_01} yr in restart file! Reading closest match = {dp_02} yr instead.', &
        dp_01 = time_to_restart_from, dp_02 = restart%time( ti))
    END IF
    IF (time_to_restart_from /= C%start_time_of_run) THEN
      CALL warning('starting run at t = {dp_01} yr with restart data at t = {dp_02} yr!', &
        dp_01 = C%start_time_of_run, dp_02 = time_to_restart_from)
    END IF
    
    ! Read the data
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_Hi,               restart%Hi,               start = (/ 1, 1,    ti /), count = (/ restart%nx, restart%ny,             1 /) ))
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_Hb,               restart%Hb,               start = (/ 1, 1,    ti /), count = (/ restart%nx, restart%ny,             1 /) ))
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_Hs,               restart%Hs,               start = (/ 1, 1,    ti /), count = (/ restart%nx, restart%ny,             1 /) ))
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_SL,               restart%SL,               start = (/ 1, 1,    ti /), count = (/ restart%nx, restart%ny,             1 /) ))
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_dHb,              restart%dHb,              start = (/ 1, 1,    ti /), count = (/ restart%nx, restart%ny,             1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( restart%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_restart_file_geometry
  SUBROUTINE read_restart_file_temperature(    restart, time_to_restart_from)
    ! Read the restart netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_restart_data),        INTENT(INOUT) :: restart
    REAL(dp),                       INTENT(IN)    :: time_to_restart_from

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_restart_file_temperature'
    INTEGER                                       :: ti, ti_min, k
    REAL(dp)                                      :: dt, dt_min
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)
    
    ! Read zeta, check if it matches the config zeta
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_zeta, restart%zeta, start=(/1/) ))
    IF (restart%nz /= C%nz) THEN
      CALL crash('vertical coordinate zeta in restart file doesnt match zeta in config!')
    ELSE
      DO k = 1, C%nz
        IF (ABS(C%zeta(k) - restart%zeta(k)) > 0.0001_dp) THEN
          CALL crash('vertical coordinate zeta in restart file doesnt match zeta in config!')
        END IF
      END DO
    END IF
    
    ! Read x,y
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_x, restart%x, start=(/1/) ))
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_y, restart%y, start=(/1/) ))
    
    ! Read time, determine which time frame to read
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_time, restart%time, start=(/1/) ))
    
    IF (time_to_restart_from < MINVAL(restart%time) .OR. time_to_restart_from > MAXVAL(restart%time)) THEN
      CALL crash('time_to_restart_from = {dp_01} outside range of restart file!', dp_01 = time_to_restart_from)
    END IF    
    
    ti_min = 0
    dt_min = 1E8_dp
    DO ti = 1, restart%nt
      dt = ABS(restart%time( ti) - time_to_restart_from)
      IF (dt < dt_min) THEN
        ti_min = ti
        dt_min = dt
      END IF
    END DO
    ti = ti_min
    
    IF (dt_min > 0._dp) THEN
      CALL warning('no exact match for time_to_restart_from = {dp_01} yr in restart file! Reading closest match = {dp_02} yr instead.', &
        dp_01 = time_to_restart_from, dp_02 = restart%time( ti))
    END IF
    IF (time_to_restart_from /= C%start_time_of_run) THEN
      CALL warning('starting run at t = {dp_01} yr with restart data at t = {dp_02} yr!', &
        dp_01 = C%start_time_of_run, dp_02 = time_to_restart_from)
    END IF
        
    ! Read the data
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_Ti,               restart%Ti,               start = (/ 1, 1, 1, ti /), count = (/ restart%nx, restart%ny, restart%nz, 1 /) ))
    
    ! Close the netcdf file
    CALL close_netcdf_file( restart%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_restart_file_temperature
  SUBROUTINE read_restart_file_SMB(            restart, time_to_restart_from)
    ! Read the restart netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_restart_data),        INTENT(INOUT) :: restart
    REAL(dp),                       INTENT(IN)    :: time_to_restart_from

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_restart_file_SMB'
    INTEGER                                       :: ti, ti_min
    REAL(dp)                                      :: dt, dt_min
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)
    
    ! Read x,y
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_x, restart%x, start=(/1/) ))
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_y, restart%y, start=(/1/) ))
    
    ! Read time, determine which time frame to read
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_time, restart%time, start=(/1/) ))
    
    IF (time_to_restart_from < MINVAL(restart%time) .OR. time_to_restart_from > MAXVAL(restart%time)) THEN
      CALL crash('time_to_restart_from = {dp_01} outside range of restart file!', dp_01 = time_to_restart_from)
    END IF    
    
    ti_min = 0
    dt_min = 1E8_dp
    DO ti = 1, restart%nt
      dt = ABS(restart%time( ti) - time_to_restart_from)
      IF (dt < dt_min) THEN
        ti_min = ti
        dt_min = dt
      END IF
    END DO
    ti = ti_min
    
    IF (dt_min > 0._dp) THEN
      CALL warning('no exact match for time_to_restart_from = {dp_01} yr in restart file! Reading closest match = {dp_02} yr instead.', &
        dp_01 = time_to_restart_from, dp_02 = restart%time( ti))
    END IF
    IF (time_to_restart_from /= C%start_time_of_run) THEN
      CALL warning('starting run at t = {dp_01} yr with restart data at t = {dp_02} yr!', &
        dp_01 = C%start_time_of_run, dp_02 = time_to_restart_from)
    END IF
        
    ! Read the data
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_FirnDepth,        restart%FirnDepth,        start = (/ 1, 1, 1, ti /), count = (/ restart%nx, restart%ny, 12,         1 /) ))
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_MeltPreviousYear, restart%MeltPreviousYear, start = (/ 1, 1,    ti /), count = (/ restart%nx, restart%ny,             1 /) ))
    
    ! Close the netcdf file
    CALL close_netcdf_file( restart%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_restart_file_SMB
  SUBROUTINE read_restart_file_isotopes(       restart, time_to_restart_from)
    ! Read the restart netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_restart_data),        INTENT(INOUT) :: restart
    REAL(dp),                       INTENT(IN)    :: time_to_restart_from

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_restart_file_isotopes'
    INTEGER                                       :: ti, ti_min
    REAL(dp)                                      :: dt, dt_min
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( restart%netcdf%filename, restart%netcdf%ncid)
    
    ! Read x,y
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_x, restart%x, start=(/1/) ))
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_y, restart%y, start=(/1/) ))
    
    ! Read time, determine which time frame to read
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_time, restart%time, start=(/1/) ))
    
    IF (time_to_restart_from < MINVAL(restart%time) .OR. time_to_restart_from > MAXVAL(restart%time)) THEN
      CALL crash('time_to_restart_from = {dp_01} outside range of restart file!', dp_01 = time_to_restart_from)
    END IF    
    
    ti_min = 0
    dt_min = 1E8_dp
    DO ti = 1, restart%nt
      dt = ABS(restart%time( ti) - time_to_restart_from)
      IF (dt < dt_min) THEN
        ti_min = ti
        dt_min = dt
      END IF
    END DO
    ti = ti_min
    
    IF (dt_min > 0._dp) THEN
      CALL warning('no exact match for time_to_restart_from = {dp_01} yr in restart file! Reading closest match = {dp_02} yr instead.', &
        dp_01 = time_to_restart_from, dp_02 = restart%time( ti))
    END IF
    IF (time_to_restart_from /= C%start_time_of_run) THEN
      CALL warning('starting run at t = {dp_01} yr with restart data at t = {dp_02} yr!', &
        dp_01 = C%start_time_of_run, dp_02 = time_to_restart_from)
    END IF
        
    ! Read the data
    CALL handle_error(nf90_get_var( restart%netcdf%ncid, restart%netcdf%id_var_IsoIce,           restart%IsoIce,           start = (/ 1, 1,    ti /), count = (/ restart%nx, restart%ny,             1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( restart%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_restart_file_isotopes
  
  SUBROUTINE read_inverse_routine_history_dT_glob(         forcing, filename)
    ! Read the inverse routine history from the specified NetCDF file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    CHARACTER(LEN=256),             INTENT(IN)    :: filename
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_inverse_routine_history_dT_glob'
    CHARACTER(LEN=256) :: dummy_char
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    dummy_char = forcing%netcdf_ins%filename
    dummy_char = filename
    
    CALL crash('need to fix the inverse routine stuff to cope with restarting!')
    
    ! Local variables:
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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_inverse_routine_history_dT_glob
  SUBROUTINE read_inverse_routine_history_dT_glob_inverse( forcing, filename)
    ! Read the inverse routine history from the specified NetCDF file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    CHARACTER(LEN=256),             INTENT(IN)    :: filename
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_inverse_routine_history_dT_glob_inverse'
    CHARACTER(LEN=256) :: dummy_char
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    dummy_char = forcing%netcdf_ins%filename
    dummy_char = filename
    
    CALL crash('need to fix the inverse routine stuff to cope with restarting!')
    
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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_inverse_routine_history_dT_glob_inverse
  SUBROUTINE read_inverse_routine_history_CO2_inverse(     forcing, filename)
    ! Read the inverse routine history from the specified NetCDF file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    CHARACTER(LEN=256),             INTENT(IN)    :: filename
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_inverse_routine_history_CO2_inverse'
    CHARACTER(LEN=256) :: dummy_char
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    dummy_char = forcing%netcdf_ins%filename
    dummy_char = filename
    
    CALL crash('need to fix the inverse routine stuff to cope with restarting!')
    
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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_inverse_routine_history_CO2_inverse
  
! Read all kinds of input files
! =============================
  
  ! Reference ice-sheet geometry (ice thickness, bed topography, and surface elevation)
  SUBROUTINE inquire_reference_geometry_file( refgeo) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_reference_geometry), INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_reference_geometry_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
      
    ! Open the netcdf file
    CALL open_netcdf_file( refgeo%netcdf%filename, refgeo%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( refgeo%netcdf%ncid, refgeo%netcdf%name_dim_x, refgeo%nx_raw, refgeo%netcdf%id_dim_x)
    CALL inquire_dim( refgeo%netcdf%ncid, refgeo%netcdf%name_dim_y, refgeo%ny_Raw, refgeo%netcdf%id_dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_x,  (/ refgeo%netcdf%id_dim_x                         /), refgeo%netcdf%id_var_x )
    CALL inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_y,  (/ refgeo%netcdf%id_dim_y                         /), refgeo%netcdf%id_var_y )
    CALL inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_Hi, (/ refgeo%netcdf%id_dim_x, refgeo%netcdf%id_dim_y /), refgeo%netcdf%id_var_Hi)
    CALL inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_Hb, (/ refgeo%netcdf%id_dim_x, refgeo%netcdf%id_dim_y /), refgeo%netcdf%id_var_Hb)
    CALL inquire_double_var( refgeo%netcdf%ncid, refgeo%netcdf%name_var_Hs, (/ refgeo%netcdf%id_dim_x, refgeo%netcdf%id_dim_y /), refgeo%netcdf%id_var_Hs)
        
    ! Close the netcdf file
    CALL close_netcdf_file( refgeo%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_reference_geometry_file
  SUBROUTINE read_reference_geometry_file( refgeo)
    ! Read reference geometry dat from a NetCDF file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_reference_geometry), INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_reference_geometry_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( refgeo%netcdf%filename, refgeo%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_x,  refgeo%x_raw,  start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_y,  refgeo%y_raw,  start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_Hi, refgeo%Hi_raw, start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_Hb, refgeo%Hb_raw, start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( refgeo%netcdf%ncid, refgeo%netcdf%id_var_Hs, refgeo%Hs_raw, start = (/ 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( refgeo%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_reference_geometry_file
  
  ! Present-day observed global climate (e.g. ERA-40)
  SUBROUTINE inquire_PD_obs_global_climate_file( clim) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_climate_snapshot_global), INTENT(INOUT) :: clim
 
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_PD_obs_global_climate_file'
    INTEGER                               :: int_dummy
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
        
    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lat,     clim%nlat,  clim%netcdf%id_dim_lat)
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lon,     clim%nlon,  clim%netcdf%id_dim_lon)
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_month,   int_dummy,  clim%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lat,      (/ clim%netcdf%id_dim_lat                                                   /), clim%netcdf%id_var_lat)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lon,      (/ clim%netcdf%id_dim_lon                                                   /), clim%netcdf%id_var_lon)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Hs,       (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat                           /), clim%netcdf%id_var_Hs)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m,      (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /), clim%netcdf%id_var_T2m)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Precip,   (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /), clim%netcdf%id_var_Precip)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Wind_WE,  (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /), clim%netcdf%id_var_Wind_WE)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Wind_SN,  (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /), clim%netcdf%id_var_Wind_SN)
        
    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_PD_obs_global_climate_file
  SUBROUTINE read_PD_obs_global_climate_file(    clim)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_climate_snapshot_global), INTENT(INOUT) :: clim
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_PD_obs_global_climate_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lon,     clim%lon,     start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lat,     clim%lat,     start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Hs,      clim%Hs,      start = (/ 1, 1    /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,     clim%T2m,     start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip,  clim%Precip,  start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Wind_WE, clim%Wind_WE, start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Wind_SN, clim%Wind_SN, start = (/ 1, 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_PD_obs_global_climate_file
  
  ! Present-day observed global ocean (e.g. WOA18)
  SUBROUTINE inquire_PD_obs_global_ocean_file( ocn) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_ocean_snapshot_global), INTENT(INOUT) :: ocn
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_PD_obs_global_ocean_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
        
    ! Open the netcdf file
    CALL open_netcdf_file( ocn%netcdf%filename, ocn%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_lat,     ocn%nlat,         ocn%netcdf%id_dim_lat    )
    CALL inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_lon,     ocn%nlon,         ocn%netcdf%id_dim_lon    )
    CALL inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_z_ocean, ocn%nz_ocean_raw, ocn%netcdf%id_dim_z_ocean)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_lat,     (/ ocn%netcdf%id_dim_lat     /), ocn%netcdf%id_var_lat    )
    CALL inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_lon,     (/ ocn%netcdf%id_dim_lon     /), ocn%netcdf%id_var_lon    )
    CALL inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_z_ocean, (/ ocn%netcdf%id_dim_z_ocean /), ocn%netcdf%id_var_z_ocean)

    CALL inquire_double_var( ocn%netcdf%ncid, TRIM(C%name_ocean_temperature), (/ ocn%netcdf%id_dim_lon, ocn%netcdf%id_dim_lat, ocn%netcdf%id_dim_z_ocean /),  ocn%netcdf%id_var_T_ocean)
    CALL inquire_double_var( ocn%netcdf%ncid, TRIM(C%name_ocean_salinity)   , (/ ocn%netcdf%id_dim_lon, ocn%netcdf%id_dim_lat, ocn%netcdf%id_dim_z_ocean /),  ocn%netcdf%id_var_S_ocean)
        
    ! Close the netcdf file
    CALL close_netcdf_file( ocn%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_PD_obs_global_ocean_file
  SUBROUTINE read_PD_obs_global_ocean_file(    ocn)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_ocean_snapshot_global), INTENT(INOUT) :: ocn
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_PD_obs_global_ocean_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( ocn%netcdf%filename, ocn%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_lon,     ocn%lon,         start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_lat,     ocn%lat,         start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_z_ocean, ocn%z_ocean_raw, start = (/ 1       /) ))

    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_T_ocean, ocn%T_ocean_raw, start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_S_ocean, ocn%S_ocean_raw, start = (/ 1, 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( ocn%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_PD_obs_global_ocean_file
  
  ! GCM global climate (climate matrix snapshots)
  SUBROUTINE inquire_GCM_global_climate_file( clim) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_climate_snapshot_global), INTENT(INOUT) :: clim
 
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_GCM_global_climate_file'
    INTEGER                                     :: int_dummy
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
        
    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lat,     clim%nlat,  clim%netcdf%id_dim_lat)
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lon,     clim%nlon,  clim%netcdf%id_dim_lon)
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_month,   int_dummy,  clim%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lat,      (/ clim%netcdf%id_dim_lat                                                   /),  clim%netcdf%id_var_lat)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lon,      (/ clim%netcdf%id_dim_lon                                                   /),  clim%netcdf%id_var_lon)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Hs,       (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat                           /),  clim%netcdf%id_var_Hs)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m,      (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /),  clim%netcdf%id_var_T2m)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Precip,   (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /),  clim%netcdf%id_var_Precip)
    !CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Wind_WE,  (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /),  clim%netcdf%id_var_Wind_WE)
    !CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Wind_SN,  (/ clim%netcdf%id_dim_lon, clim%netcdf%id_dim_lat, clim%netcdf%id_dim_month /),  clim%netcdf%id_var_Wind_SN)
        
    ! Close the netcdf file
    CALL close_netcdf_file(clim%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_GCM_global_climate_file
  SUBROUTINE read_GCM_global_climate_file(    clim)
  
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_climate_snapshot_global), INTENT(INOUT) :: clim
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_GCM_global_climate_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file(clim%netcdf%filename, clim%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lon,     clim%lon,     start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lat,     clim%lat,     start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Hs,      clim%Hs,      start = (/ 1, 1    /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,     clim%T2m,     start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip,  clim%Precip,  start = (/ 1, 1, 1 /) ))
    !CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Wind_WE, clim%Wind_WE, start = (/ 1, 1, 1 /) ))
    !CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Wind_SN, clim%Wind_SN, start = (/ 1, 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file(clim%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_GCM_global_climate_file
  
  ! GCM global ocean (ocean matrix snapshots)
  SUBROUTINE inquire_GCM_global_ocean_file( ocn) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_ocean_snapshot_global), INTENT(INOUT) :: ocn
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_GCM_global_ocean_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
        
    ! Open the netcdf file
    CALL open_netcdf_file( ocn%netcdf%filename, ocn%netcdf%ncid)
 
     ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_lat,     ocn%nlat,         ocn%netcdf%id_dim_lat    )
    CALL inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_lon,     ocn%nlon,         ocn%netcdf%id_dim_lon    )
    CALL inquire_dim( ocn%netcdf%ncid, ocn%netcdf%name_dim_z_ocean, ocn%nz_ocean_raw, ocn%netcdf%id_dim_z_ocean)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_lat,     (/ ocn%netcdf%id_dim_lat     /), ocn%netcdf%id_var_lat    )
    CALL inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_lon,     (/ ocn%netcdf%id_dim_lon     /), ocn%netcdf%id_var_lon    )
    CALL inquire_double_var( ocn%netcdf%ncid, ocn%netcdf%name_var_z_ocean, (/ ocn%netcdf%id_dim_z_ocean /), ocn%netcdf%id_var_z_ocean)

    CALL inquire_double_var( ocn%netcdf%ncid, TRIM(C%name_ocean_temperature), (/ ocn%netcdf%id_dim_lon, ocn%netcdf%id_dim_lat, ocn%netcdf%id_dim_z_ocean /), ocn%netcdf%id_var_T_ocean)
    CALL inquire_double_var( ocn%netcdf%ncid, TRIM(C%name_ocean_salinity)   , (/ ocn%netcdf%id_dim_lon, ocn%netcdf%id_dim_lat, ocn%netcdf%id_dim_z_ocean /), ocn%netcdf%id_var_S_ocean)
   
    ! Close the netcdf file
    CALL close_netcdf_file( ocn%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_GCM_global_ocean_file
  SUBROUTINE read_GCM_global_ocean_file(    ocn)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_ocean_snapshot_global), INTENT(INOUT) :: ocn
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_GCM_global_ocean_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( ocn%netcdf%filename, ocn%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_lon,     ocn%lon,         start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_lat,     ocn%lat,         start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_z_ocean, ocn%z_ocean_raw, start = (/ 1       /) ))

    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_T_ocean, ocn%T_ocean_raw, start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( ocn%netcdf%ncid, ocn%netcdf%id_var_S_ocean, ocn%S_ocean_raw, start = (/ 1, 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( ocn%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_GCM_global_ocean_file  
  
  ! High-resolution geometry used for extrapolating ocean data
  SUBROUTINE inquire_hires_geometry_file( hires)
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_highres_ocean_data), INTENT(INOUT) :: hires
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_hires_geometry_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
        
    ! Open the netcdf file
    CALL open_netcdf_file( hires%netcdf_geo%filename, hires%netcdf_geo%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( hires%netcdf_geo%ncid, hires%netcdf_geo%name_dim_x, hires%grid%nx, hires%netcdf_geo%id_dim_x)
    CALL inquire_dim( hires%netcdf_geo%ncid, hires%netcdf_geo%name_dim_y, hires%grid%ny, hires%netcdf_geo%id_dim_y)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( hires%netcdf_geo%ncid, hires%netcdf_geo%name_var_x,  (/ hires%netcdf_geo%id_dim_x                        /), hires%netcdf_geo%id_var_x )
    CALL inquire_double_var( hires%netcdf_geo%ncid, hires%netcdf_geo%name_var_y,  (/ hires%netcdf_geo%id_dim_y                        /), hires%netcdf_geo%id_var_y )
    CALL inquire_double_var( hires%netcdf_geo%ncid, hires%netcdf_geo%name_var_Hi, (/ hires%netcdf_geo%id_dim_x, hires%netcdf_geo%id_dim_y /), hires%netcdf_geo%id_var_Hi)
    CALL inquire_double_var( hires%netcdf_geo%ncid, hires%netcdf_geo%name_var_Hb, (/ hires%netcdf_geo%id_dim_x, hires%netcdf_geo%id_dim_y /), hires%netcdf_geo%id_var_Hb)
        
    ! Close the netcdf file
    CALL close_netcdf_file( hires%netcdf_geo%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_hires_geometry_file
  SUBROUTINE read_hires_geometry_file(    hires)
    ! Read the high-resolution geometry netcdf file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_highres_ocean_data), INTENT(INOUT) :: hires
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_hires_geometry_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( hires%netcdf_geo%filename, hires%netcdf_geo%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( hires%netcdf_geo%ncid, hires%netcdf_geo%id_var_x,      hires%grid%x, start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( hires%netcdf_geo%ncid, hires%netcdf_geo%id_var_y,      hires%grid%y, start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( hires%netcdf_geo%ncid, hires%netcdf_geo%id_var_Hi,     hires%Hi,     start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( hires%netcdf_geo%ncid, hires%netcdf_geo%id_var_Hb,     hires%Hb,     start = (/ 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( hires%netcdf_geo%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_hires_geometry_file
  
  ! Create/read an extrapolated ocean data file
  SUBROUTINE create_extrapolated_ocean_file(  hires, hires_ocean_filename)
    ! Create a new folder extrapolated ocean data file
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_highres_ocean_data),       INTENT(INOUT) :: hires
    CHARACTER(LEN=256),                  INTENT(IN)    :: hires_ocean_filename

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_extrapolated_ocean_file'
    LOGICAL                                            :: file_exists
    INTEGER                                            :: x, y, z
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Create a new file and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    
    hires%netcdf%filename = hires_ocean_filename
    INQUIRE(EXIST=file_exists, FILE = TRIM( hires%netcdf%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( hires%netcdf%filename) // '" already exists!')
    END IF
    
    ! Create hires%netcdf file
    !WRITE(0,*) ' Creating new hires%netcdf file at ', TRIM( hires%netcdf%filename)
    CALL handle_error(nf90_create( hires%netcdf%filename, IOR(nf90_clobber,nf90_share), hires%netcdf%ncid))
        
    ! Define dimensions:
    CALL create_dim( hires%netcdf%ncid, hires%netcdf%name_dim_x,       hires%grid%nx, hires%netcdf%id_dim_x      )
    CALL create_dim( hires%netcdf%ncid, hires%netcdf%name_dim_y,       hires%grid%ny, hires%netcdf%id_dim_y      )
    CALL create_dim( hires%netcdf%ncid, hires%netcdf%name_dim_z_ocean, C%nz_ocean,    hires%netcdf%id_dim_z_ocean)
    
    ! Placeholders for the dimension ID's, for shorter code
    x = hires%netcdf%id_dim_x
    y = hires%netcdf%id_dim_y
    z = hires%netcdf%id_dim_z_ocean
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the hires%netcdf file.
    
    ! Dimension variables
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_x,       [x      ], hires%netcdf%id_var_x,       long_name='X-coordinate', units='m')
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_y,       [   y   ], hires%netcdf%id_var_y,       long_name='Y-coordinate', units='m')
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_z_ocean, [      z], hires%netcdf%id_var_z_ocean, long_name='Depth in ocean', units='m')
    
    ! Extrapolated ocean data
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_T_ocean, [x, y, z], hires%netcdf%id_var_T_ocean, long_name='3-D ocean temperature', units='K')
    CALL create_double_var( hires%netcdf%ncid, hires%netcdf%name_var_S_ocean, [x, y, z], hires%netcdf%id_var_S_ocean, long_name='3-D ocean salinity', units='PSU')

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( hires%netcdf%ncid))
    
    ! Write the data
    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_x,        hires%grid%x ))
    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_y,        hires%grid%y ))
    CALL handle_error( nf90_put_var( hires%netcdf%ncid, hires%netcdf%id_var_z_ocean,  C%z_ocean    ))
    
    CALL write_data_to_file_dp_3D(  hires%netcdf%ncid, hires%grid%nx, hires%grid%ny, C%nz_ocean, hires%netcdf%id_var_T_ocean, hires%T_ocean, (/ 1,1,1 /) )
    CALL write_data_to_file_dp_3D(  hires%netcdf%ncid, hires%grid%nx, hires%grid%ny, C%nz_ocean, hires%netcdf%id_var_S_ocean, hires%S_ocean, (/ 1,1,1 /) )
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( hires%netcdf%ncid))
    
    ! Close the file
    CALL close_netcdf_file( hires%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE create_extrapolated_ocean_file
  SUBROUTINE inquire_extrapolated_ocean_file( hires) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_highres_ocean_data), INTENT(INOUT) :: hires
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_extrapolated_ocean_file'
    INTEGER                                      :: int_dummy
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
        
    ! Open the netcdf file
    CALL open_netcdf_file( hires%netcdf%filename, hires%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist, return their lengths.
    CALL inquire_dim( hires%netcdf%ncid, hires%netcdf%name_dim_x,       hires%grid%nx,   hires%netcdf%id_dim_x      )
    CALL inquire_dim( hires%netcdf%ncid, hires%netcdf%name_dim_y,       hires%grid%ny,   hires%netcdf%id_dim_y      )
    CALL inquire_dim( hires%netcdf%ncid, hires%netcdf%name_dim_z_ocean, int_dummy,  hires%netcdf%id_dim_z_ocean)
    
    ! Safety
    IF (int_dummy /= C%nz_ocean) THEN
      CALL crash('nz_ocean in file "' // TRIM( hires%netcdf%filename) // '" doesnt match ice model settings!')
    END IF

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_x,       (/ hires%netcdf%id_dim_x                                                     /), hires%netcdf%id_var_x      )
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_y,       (/                        hires%netcdf%id_dim_y                              /), hires%netcdf%id_var_y      )
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_z_ocean, (/                                               hires%netcdf%id_dim_z_ocean /), hires%netcdf%id_var_z_ocean)
    
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_T_ocean, (/ hires%netcdf%id_dim_x, hires%netcdf%id_dim_y, hires%netcdf%id_dim_z_ocean /), hires%netcdf%id_var_T_ocean)
    CALL inquire_double_var( hires%netcdf%ncid, hires%netcdf%name_var_S_ocean, (/ hires%netcdf%id_dim_x, hires%netcdf%id_dim_y, hires%netcdf%id_dim_z_ocean /), hires%netcdf%id_var_S_ocean)
        
    ! Close the netcdf file
    CALL close_netcdf_file( hires%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_extrapolated_ocean_file
  SUBROUTINE read_extrapolated_ocean_file(    hires)
    ! Read the extrapolated ocean data netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_highres_ocean_data), INTENT(INOUT) :: hires
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_extrapolated_ocean_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( hires%netcdf%filename, hires%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( hires%netcdf%ncid, hires%netcdf%id_var_x,       hires%grid%x,  start = (/ 1       /) ))
    CALL handle_error(nf90_get_var( hires%netcdf%ncid, hires%netcdf%id_var_y,       hires%grid%y,  start = (/ 1       /) ))
    
    CALL handle_error(nf90_get_var( hires%netcdf%ncid, hires%netcdf%id_var_T_ocean, hires%T_ocean, start = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( hires%netcdf%ncid, hires%netcdf%id_var_S_ocean, hires%S_ocean, start = (/ 1, 1, 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( hires%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_extrapolated_ocean_file
  
  ! Insolation solution (e.g. Laskar 2004)
  SUBROUTINE inquire_insolation_file( forcing)
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
 
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_insolation_file'
    INTEGER                                :: int_dummy
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
            
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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_insolation_file
  SUBROUTINE read_insolation_file_timeframes( forcing, ti0, ti1) 
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    INTEGER,                        INTENT(IN)    :: ti0, ti1
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_insolation_file_timeframes'
    INTEGER                                       :: mi, li
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Q_temp0, Q_temp1
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
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
      forcing%ins_Q_TOA0( li,mi) = Q_temp0( 1,mi,li)
      forcing%ins_Q_TOA1( li,mi) = Q_temp1( 1,mi,li)
    END DO
    END DO
        
    ! Clean up temporary memory
    DEALLOCATE(Q_temp0)
    DEALLOCATE(Q_temp1)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
   
  END SUBROUTINE read_insolation_file_timeframes
  SUBROUTINE read_insolation_file_time_lat( forcing) 
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_insolation_file_time_lat'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_ins%filename, forcing%netcdf_ins%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_time,    forcing%ins_time,    start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ins%ncid, forcing%netcdf_ins%id_var_lat,     forcing%ins_lat,     start = (/ 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf_ins%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_insolation_file_time_lat
  
  ! Geothermal heat flux
  SUBROUTINE inquire_geothermal_heat_flux_file( forcing)
    IMPLICIT NONE

    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_geothermal_heat_flux_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

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
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE inquire_geothermal_heat_flux_file
  SUBROUTINE read_geothermal_heat_flux_file( forcing)
  
    USE parameters_module, ONLY: sec_per_year
  
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_geothermal_heat_flux_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Read data
    CALL open_netcdf_file(forcing%netcdf_ghf%filename, forcing%netcdf_ghf%ncid)
    CALL handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_lon, forcing%ghf_lon, start=(/1   /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_lat, forcing%ghf_lat, start=(/1   /) ))
    CALL handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_ghf, forcing%ghf_ghf, start=(/1, 1/), count=(/forcing%ghf_nlon, forcing%ghf_nlat/), stride=(/1, 1/) ))
    CALL close_netcdf_file(forcing%netcdf_ghf%ncid)

    ! Convert from W m-2 (J m-2 s-1) to J m-2 yr-1
    forcing%ghf_ghf = forcing%ghf_ghf * sec_per_year
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_geothermal_heat_flux_file

  ! Direct global climate forcing
  SUBROUTINE inquire_direct_global_climate_forcing_file( clim)
  
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_direct_climate_forcing_global), INTENT(INOUT) :: clim
 
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_direct_global_climate_forcing_file'
    INTEGER                                :: lon,lat,t,m
    INTEGER                                :: int_dummy
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_global') THEN
      CALL crash('should only be called when choice_climate_model = "direct_global"!')
    END IF
            
    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
        
    ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lon,   clim%nlon,   clim%netcdf%id_dim_lon  )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lat,   clim%nlat,   clim%netcdf%id_dim_lat  )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_month, int_dummy,   clim%netcdf%id_dim_month)
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_time,  clim%nyears, clim%netcdf%id_dim_time )
    
    ! Abbreviate dimension ID's for more readable code
    lon = clim%netcdf%id_dim_lon
    lat = clim%netcdf%id_dim_lat
    t   = clim%netcdf%id_dim_time
    m   = clim%netcdf%id_dim_month
    
    ! Inquire variable id's. Make sure that each variable has the correct dimensions.
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lon,    (/ lon            /), clim%netcdf%id_var_lon   )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lat,    (/      lat       /), clim%netcdf%id_var_lat   )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_month,  (/           m    /), clim%netcdf%id_var_month )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_time,   (/              t /), clim%netcdf%id_var_time  )
    
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m,    (/ lon, lat, m, t /), clim%netcdf%id_var_T2m   )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Precip, (/ lon, lat, m, t /), clim%netcdf%id_var_Precip)
        
    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_direct_global_climate_forcing_file
  SUBROUTINE read_direct_global_climate_file_timeframes( clim, ti0, ti1)
  
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_direct_climate_forcing_global), INTENT(INOUT) :: clim
    INTEGER,                        INTENT(IN)    :: ti0, ti1

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_global_climate_file_timeframes'
    INTEGER                                       :: loni,lati,m
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: T2m_temp0, T2m_temp1, Precip_temp0, Precip_temp1
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_global') THEN
      CALL crash('should only be called when choice_climate_model = "direct_global"!')
    END IF
      
    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE(    T2m_temp0( clim%nlon, clim%nlat, 12, 1))
    ALLOCATE(    T2m_temp1( clim%nlon, clim%nlat, 12, 1))
    ALLOCATE( Precip_temp0( clim%nlon, clim%nlat, 12, 1))
    ALLOCATE( Precip_temp1( clim%nlon, clim%nlat, 12, 1))
      
    ! Open netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,    T2m_temp0,    start = (/ 1, 1, 1, ti0 /), count = (/ clim%nlon, clim%nlat, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,    T2m_temp1,    start = (/ 1, 1, 1, ti1 /), count = (/ clim%nlon, clim%nlat, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip, Precip_temp0, start = (/ 1, 1, 1, ti0 /), count = (/ clim%nlon, clim%nlat, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip, Precip_temp1, start = (/ 1, 1, 1, ti1 /), count = (/ clim%nlon, clim%nlat, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    
     ! Close netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Store the data in the shared memory structure
    DO m    = 1, 12
    DO loni = 1, clim%nlon
    DO lati = 1, clim%nlat 
      clim%T2m0(    loni,lati,m) =    T2m_temp0( loni,lati,m,1)
      clim%T2m1(    loni,lati,m) =    T2m_temp1( loni,lati,m,1)
      clim%Precip0( loni,lati,m) = Precip_temp0( loni,lati,m,1)
      clim%Precip1( loni,lati,m) = Precip_temp1( loni,lati,m,1)
    END DO
    END DO
    END DO

    ! Clean up after yourself
    DEALLOCATE(    T2m_temp0)
    DEALLOCATE(    T2m_temp1)
    DEALLOCATE( Precip_temp0)
    DEALLOCATE( Precip_temp1)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
      
  END SUBROUTINE read_direct_global_climate_file_timeframes
  SUBROUTINE read_direct_global_climate_file_time_latlon( clim)
  
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_direct_climate_forcing_global), INTENT(INOUT) :: clim
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_global_climate_file_time_latlon'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_global') THEN
      CALL crash('should only be called when choice_climate_model = "direct_global"!')
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_time, clim%time, start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lat,  clim%lat,  start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lon,  clim%lon,  start = (/ 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_direct_global_climate_file_time_latlon
  
  ! Direct regional climate forcing
  SUBROUTINE inquire_direct_regional_climate_forcing_file( clim)
  
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_direct_climate_forcing_regional), INTENT(INOUT) :: clim
 
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_direct_regional_climate_forcing_file'
    INTEGER                                :: x,y,t,m
    INTEGER                                :: int_dummy
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_regional') THEN
      CALL crash('should only be called when choice_climate_model = "direct_regional"!')
    END IF
            
    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
        
    ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_x,     clim%nx_raw, clim%netcdf%id_dim_x    )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_y,     clim%ny_raw, clim%netcdf%id_dim_y    )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_month, int_dummy,   clim%netcdf%id_dim_month)
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_time,  clim%nyears, clim%netcdf%id_dim_time )
    
    ! Abbreviate dimension ID's for more readable code
    x = clim%netcdf%id_dim_x
    y = clim%netcdf%id_dim_y
    t = clim%netcdf%id_dim_time
    m = clim%netcdf%id_dim_month
    
    ! Inquire variable id's. Make sure that each variable has the correct dimensions.
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_x,      (/ x          /), clim%netcdf%id_var_x     )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_y,      (/    y       /), clim%netcdf%id_var_y     )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_month,  (/       m    /), clim%netcdf%id_var_month )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_time,   (/          t /), clim%netcdf%id_var_time  )
    
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m,    (/ x, y, m, t /), clim%netcdf%id_var_T2m   )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_Precip, (/ x, y, m, t /), clim%netcdf%id_var_Precip)
        
    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_direct_regional_climate_forcing_file
  SUBROUTINE read_direct_regional_climate_file_timeframes( clim, ti0, ti1)
  
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_direct_climate_forcing_regional), INTENT(INOUT) :: clim
    INTEGER,                        INTENT(IN)    :: ti0, ti1

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_regional_climate_file_timeframes'
    INTEGER                                       :: i,j,m
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE     :: T2m_temp0, T2m_temp1, Precip_temp0, Precip_temp1
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_regional') THEN
      CALL crash('should only be called when choice_climate_model = "direct_regional"!')
    END IF
      
    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE(    T2m_temp0( clim%nx_raw, clim%ny_raw, 12, 1))
    ALLOCATE(    T2m_temp1( clim%nx_raw, clim%ny_raw, 12, 1))
    ALLOCATE( Precip_temp0( clim%nx_raw, clim%ny_raw, 12, 1))
    ALLOCATE( Precip_temp1( clim%nx_raw, clim%ny_raw, 12, 1))
      
    ! Open netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,    T2m_temp0,    start = (/ 1, 1, 1, ti0 /), count = (/ clim%nx_raw, clim%ny_raw, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m,    T2m_temp1,    start = (/ 1, 1, 1, ti1 /), count = (/ clim%nx_raw, clim%nx_raw, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip, Precip_temp0, start = (/ 1, 1, 1, ti0 /), count = (/ clim%nx_raw, clim%nx_raw, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_Precip, Precip_temp1, start = (/ 1, 1, 1, ti1 /), count = (/ clim%nx_raw, clim%nx_raw, 12, 1 /), stride = (/ 1, 1, 1, 1 /) ))
    
     ! Close netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Store the data in the shared memory structure
    DO m = 1, 12
    DO i = 1, clim%nx_raw
    DO j = 1, clim%ny_raw
      clim%T2m0_raw(    m,j,i) =    T2m_temp0( i,j,m,1)
      clim%T2m1_raw(    m,j,i) =    T2m_temp1( i,j,m,1)
      clim%Precip0_raw( m,j,i) = Precip_temp0( i,j,m,1)
      clim%Precip1_raw( m,j,i) = Precip_temp1( i,j,m,1)
    END DO
    END DO
    END DO

    ! Clean up after yourself
    DEALLOCATE(    T2m_temp0)
    DEALLOCATE(    T2m_temp1)
    DEALLOCATE( Precip_temp0)
    DEALLOCATE( Precip_temp1)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
      
  END SUBROUTINE read_direct_regional_climate_file_timeframes
  SUBROUTINE read_direct_regional_climate_file_time_xy( clim)
  
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_direct_climate_forcing_regional), INTENT(INOUT) :: clim
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_regional_climate_file_time_xy'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Safety
    IF (.NOT. C%choice_climate_model == 'direct_regional') THEN
      CALL crash('should only be called when choice_climate_model = "direct_regional"!')
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_time, clim%time,  start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_x,    clim%x_raw, start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_y,    clim%y_raw, start = (/ 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_direct_regional_climate_file_time_xy

  ! Direct global SMB forcing
  SUBROUTINE inquire_direct_global_SMB_forcing_file( clim)
  
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_direct_SMB_forcing_global), INTENT(INOUT) :: clim
 
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_direct_global_SMB_forcing_file'
    INTEGER                                :: lon,lat,t
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_global') THEN
      CALL crash('should only be called when choice_SMB_model = "direct_global"!')
    END IF
            
    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
        
    ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lon,   clim%nlon,   clim%netcdf%id_dim_lon  )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_lat,   clim%nlat,   clim%netcdf%id_dim_lat  )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_time,  clim%nyears, clim%netcdf%id_dim_time )
    
    ! Abbreviate dimension ID's for more readable code
    lon = clim%netcdf%id_dim_lon
    lat = clim%netcdf%id_dim_lat
    t   = clim%netcdf%id_dim_time
    
    ! Inquire variable id's. Make sure that each variable has the correct dimensions.
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lon,      (/ lon         /), clim%netcdf%id_var_lon     )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_lat,      (/      lat    /), clim%netcdf%id_var_lat     )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_time,     (/           t /), clim%netcdf%id_var_time    )
    
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m_year, (/ lon, lat, t /), clim%netcdf%id_var_T2m_year)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_SMB_year, (/ lon, lat, t /), clim%netcdf%id_var_SMB_year)
        
    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_direct_global_SMB_forcing_file
  SUBROUTINE read_direct_global_SMB_file_timeframes( clim, ti0, ti1)
  
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_direct_SMB_forcing_global), INTENT(INOUT) :: clim
    INTEGER,                        INTENT(IN)    :: ti0, ti1

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_global_SMB_file_timeframes'
    INTEGER                                       :: loni,lati
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE     :: T2m_temp0, T2m_temp1, SMB_temp0, SMB_temp1
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_global') THEN
      CALL crash('should only be called when choice_SMB_model = "direct_global"!')
    END IF
      
    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE( T2m_temp0( clim%nlon, clim%nlat, 1))
    ALLOCATE( T2m_temp1( clim%nlon, clim%nlat, 1))
    ALLOCATE( SMB_temp0( clim%nlon, clim%nlat, 1))
    ALLOCATE( SMB_temp1( clim%nlon, clim%nlat, 1))
      
    ! Open netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m_year, T2m_temp0, start = (/ 1, 1, ti0 /), count = (/ clim%nlon, clim%nlat, 1 /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m_year, T2m_temp1, start = (/ 1, 1, ti1 /), count = (/ clim%nlon, clim%nlat, 1 /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_SMB_year, SMB_temp0, start = (/ 1, 1, ti0 /), count = (/ clim%nlon, clim%nlat, 1 /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_SMB_year, SMB_temp1, start = (/ 1, 1, ti1 /), count = (/ clim%nlon, clim%nlat, 1 /), stride = (/ 1, 1, 1 /) ))
    
     ! Close netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Store the data in the shared memory structure
    DO loni = 1, clim%nlon
    DO lati = 1, clim%nlat 
      clim%T2m_year0( loni,lati) = T2m_temp0( loni,lati,1)
      clim%T2m_year1( loni,lati) = T2m_temp1( loni,lati,1)
      clim%SMB_year0( loni,lati) = SMB_temp0( loni,lati,1)
      clim%SMB_year1( loni,lati) = SMB_temp1( loni,lati,1)
    END DO
    END DO

    ! Clean up after yourself
    DEALLOCATE( T2m_temp0)
    DEALLOCATE( T2m_temp1)
    DEALLOCATE( SMB_temp0)
    DEALLOCATE( SMB_temp1)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
      
  END SUBROUTINE read_direct_global_SMB_file_timeframes
  SUBROUTINE read_direct_global_SMB_file_time_latlon( clim)
  
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_direct_SMB_forcing_global), INTENT(INOUT) :: clim
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_global_SMB_file_time_latlon'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_global') THEN
      CALL crash('should only be called when choice_SMB_model = "direct_global"!')
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_time, clim%time, start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lat,  clim%lat,  start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_lon,  clim%lon,  start = (/ 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_direct_global_SMB_file_time_latlon
  
  ! Direct regional SMB forcing
  SUBROUTINE inquire_direct_regional_SMB_forcing_file( clim)
  
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_direct_SMB_forcing_regional), INTENT(INOUT) :: clim
 
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_direct_regional_SMB_forcing_file'
    INTEGER                                :: x,y,t
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
      CALL crash('should only be called when choice_SMB_model = "direct_regional"!')
    END IF
            
    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
        
    ! Inquire dimensions id's. Check that all required dimensions exist, and return their lengths.
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_x,     clim%nx_raw, clim%netcdf%id_dim_x    )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_y,     clim%ny_raw, clim%netcdf%id_dim_y    )
    CALL inquire_dim( clim%netcdf%ncid, clim%netcdf%name_dim_time,  clim%nyears, clim%netcdf%id_dim_time )
    
    ! Abbreviate dimension ID's for more readable code
    x = clim%netcdf%id_dim_x
    y = clim%netcdf%id_dim_y
    t = clim%netcdf%id_dim_time
    
    ! Inquire variable id's. Make sure that each variable has the correct dimensions.
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_x,        (/ x       /), clim%netcdf%id_var_x       )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_y,        (/    y    /), clim%netcdf%id_var_y       )
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_time,     (/       t /), clim%netcdf%id_var_time    )
    
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_T2m_year, (/ x, y, t /), clim%netcdf%id_var_T2m_year)
    CALL inquire_double_var( clim%netcdf%ncid, clim%netcdf%name_var_SMB_year, (/ x, y, t /), clim%netcdf%id_var_SMB_year)
        
    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_direct_regional_SMB_forcing_file
  SUBROUTINE read_direct_regional_SMB_file_timeframes( clim, ti0, ti1)
  
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_direct_SMB_forcing_regional), INTENT(INOUT) :: clim
    INTEGER,                        INTENT(IN)    :: ti0, ti1

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_regional_SMB_file_timeframes'
    INTEGER                                       :: i,j
    REAL(dp), DIMENSION(:,:,:  ), ALLOCATABLE     :: T2m_temp0, T2m_temp1, SMB_temp0, SMB_temp1
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
      CALL crash('should only be called when choice_SMB_model = "direct_regional"!')
    END IF
      
    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE( T2m_temp0( clim%nx_raw, clim%ny_raw, 1))
    ALLOCATE( T2m_temp1( clim%nx_raw, clim%ny_raw, 1))
    ALLOCATE( SMB_temp0( clim%nx_raw, clim%ny_raw, 1))
    ALLOCATE( SMB_temp1( clim%nx_raw, clim%ny_raw, 1))
      
    ! Open netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m_year, T2m_temp0, start = (/ 1, 1, ti0 /), count = (/ clim%nx_raw, clim%ny_raw, 1 /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_T2m_year, T2m_temp1, start = (/ 1, 1, ti1 /), count = (/ clim%nx_raw, clim%nx_raw, 1 /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_SMB_year, SMB_temp0, start = (/ 1, 1, ti0 /), count = (/ clim%nx_raw, clim%ny_raw, 1 /), stride = (/ 1, 1, 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_SMB_year, SMB_temp1, start = (/ 1, 1, ti1 /), count = (/ clim%nx_raw, clim%nx_raw, 1 /), stride = (/ 1, 1, 1 /) ))
    
     ! Close netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)

    ! Store the data in the shared memory structure
    DO i = 1, clim%nx_raw
    DO j = 1, clim%ny_raw
      clim%T2m_year0_raw( j,i) = T2m_temp0( i,j,1)
      clim%T2m_year1_raw( j,i) = T2m_temp1( i,j,1)
      clim%SMB_year0_raw( j,i) = SMB_temp0( i,j,1)
      clim%SMB_year1_raw( j,i) = SMB_temp1( i,j,1)
    END DO
    END DO

    ! Clean up after yourself
    DEALLOCATE( T2m_temp0)
    DEALLOCATE( T2m_temp1)
    DEALLOCATE( SMB_temp0)
    DEALLOCATE( SMB_temp1)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
      
  END SUBROUTINE read_direct_regional_SMB_file_timeframes
  SUBROUTINE read_direct_regional_SMB_file_time_xy( clim)
  
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_direct_SMB_forcing_regional), INTENT(INOUT) :: clim
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_direct_regional_SMB_file_time_xy'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Safety
    IF (.NOT. C%choice_SMB_model == 'direct_regional') THEN
      CALL crash('should only be called when choice_SMB_model = "direct_regional"!')
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( clim%netcdf%filename, clim%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_time, clim%time,  start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_x,    clim%x_raw, start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( clim%netcdf%ncid, clim%netcdf%id_var_y,    clim%y_raw, start = (/ 1 /) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file( clim%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_direct_regional_SMB_file_time_xy
  
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
  
  ! Inverted basal roughness
  SUBROUTINE create_BIV_bed_roughness_file( grid, ice)
    ! Create a new folder extrapolated ocean data file
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_BIV_bed_roughness_file'
    TYPE(type_netcdf_BIV_bed_roughness)                :: netcdf
    LOGICAL                                            :: file_exists
    INTEGER                                            :: x, y
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Create a new file and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    
    netcdf%filename = TRIM(C%output_dir) // TRIM(C%BIVgeo_filename_output)
    INQUIRE(EXIST=file_exists, FILE = TRIM( netcdf%filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( netcdf%filename) // '" already exists!')
    END IF
    
    ! Create netcdf file
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Writing inverted bed roughness to file "', TRIM( netcdf%filename), '"...'
    CALL handle_error(nf90_create( netcdf%filename, IOR(nf90_clobber,nf90_share), netcdf%ncid))
        
    ! Define dimensions:
    CALL create_dim( netcdf%ncid, netcdf%name_dim_x, grid%nx, netcdf%id_dim_x)
    CALL create_dim( netcdf%ncid, netcdf%name_dim_y, grid%ny, netcdf%id_dim_y)
    
    ! Placeholders for the dimension ID's, for shorter code
    x = netcdf%id_dim_x
    y = netcdf%id_dim_y
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.
    
    ! Dimension variables
    CALL create_double_var( netcdf%ncid, netcdf%name_var_x,       [x      ], netcdf%id_var_x,       long_name='X-coordinate', units='m')
    CALL create_double_var( netcdf%ncid, netcdf%name_var_y,       [   y   ], netcdf%id_var_y,       long_name='Y-coordinate', units='m')
    
    ! Bed roughness
    IF     (C%choice_sliding_law == 'no_sliding') THEN
      CALL crash('not defined for choice_sliding_law = "no_sliding"!')
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_beta_sq,  [x, y], netcdf%id_var_beta_sq,  long_name='Power-law friction coefficient', units='[Pa m^−1/3 yr^1/3]')
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_phi_fric, [x, y], netcdf%id_var_phi_fric, long_name='Till friction angle', units='degrees')
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_beta_sq,  [x, y], netcdf%id_var_beta_sq,  long_name='Coulomb-law friction coefficient', units='unitless')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_beta_sq,  [x, y], netcdf%id_var_beta_sq,  long_name='Power-law friction coefficient', units='[Pa m^−1/3 yr^1/3]')
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_beta_sq,  [x, y], netcdf%id_var_beta_sq,  long_name='Coulomb-law friction coefficient', units='unitless')
      CALL create_double_var( netcdf%ncid, netcdf%name_var_beta_sq,  [x, y], netcdf%id_var_beta_sq,  long_name='Power-law friction coefficient', units='[Pa m^−1/3 yr^1/3]')
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      CALL create_double_var( netcdf%ncid, netcdf%name_var_phi_fric, [x, y], netcdf%id_var_phi_fric, long_name='Till friction angle', units='degrees')
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! Leave definition mode:
    CALL handle_error(nf90_enddef( netcdf%ncid))
    
    ! Write the data
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_x,        grid%x ))
    CALL handle_error( nf90_put_var( netcdf%ncid, netcdf%id_var_y,        grid%y ))
    
    IF     (C%choice_sliding_law == 'no_sliding') THEN
      CALL crash('not defined for choice_sliding_law = "no_sliding"!')
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      CALL write_data_to_file_dp_2D(  netcdf%ncid, grid%nx, grid%ny, netcdf%id_var_beta_sq,  ice%beta_sq_a,  (/ 1,1 /) )
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      CALL write_data_to_file_dp_2D(  netcdf%ncid, grid%nx, grid%ny, netcdf%id_var_phi_fric, ice%phi_fric_a, (/ 1,1 /) )
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      CALL write_data_to_file_dp_2D(  netcdf%ncid, grid%nx, grid%ny, netcdf%id_var_alpha_sq, ice%alpha_sq_a, (/ 1,1 /) )
      CALL write_data_to_file_dp_2D(  netcdf%ncid, grid%nx, grid%ny, netcdf%id_var_beta_sq,  ice%beta_sq_a,  (/ 1,1 /) )
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      CALL write_data_to_file_dp_2D(  netcdf%ncid, grid%nx, grid%ny, netcdf%id_var_alpha_sq, ice%alpha_sq_a, (/ 1,1 /) )
      CALL write_data_to_file_dp_2D(  netcdf%ncid, grid%nx, grid%ny, netcdf%id_var_beta_sq,  ice%beta_sq_a,  (/ 1,1 /) )
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      CALL write_data_to_file_dp_2D(  netcdf%ncid, grid%nx, grid%ny, netcdf%id_var_phi_fric, ice%phi_fric_a, (/ 1,1 /) )
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( netcdf%ncid))
    
    ! Close the file
    CALL close_netcdf_file( netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE create_BIV_bed_roughness_file
  SUBROUTINE inquire_BIV_bed_roughness_file( BIV)
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_BIV_bed_roughness),        INTENT(INOUT) :: BIV
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_BIV_bed_roughness_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
        
    ! Open the netcdf file
    CALL open_netcdf_file( BIV%netcdf%filename, BIV%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist, return their lengths.
    CALL inquire_dim( BIV%netcdf%ncid, BIV%netcdf%name_dim_x,       BIV%nx,   BIV%netcdf%id_dim_x      )
    CALL inquire_dim( BIV%netcdf%ncid, BIV%netcdf%name_dim_y,       BIV%ny,   BIV%netcdf%id_dim_y      )

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( BIV%netcdf%ncid, BIV%netcdf%name_var_x,       (/ BIV%netcdf%id_dim_x                        /), BIV%netcdf%id_var_x      )
    CALL inquire_double_var( BIV%netcdf%ncid, BIV%netcdf%name_var_y,       (/                        BIV%netcdf%id_dim_y /), BIV%netcdf%id_var_y      )
    
    IF     (C%choice_sliding_law == 'no_sliding') THEN
      CALL crash('not defined for choice_sliding_law = "no_sliding"!')
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      CALL inquire_double_var( BIV%netcdf%ncid, BIV%netcdf%name_var_beta_sq,  (/ BIV%netcdf%id_dim_x, BIV%netcdf%id_dim_y /), BIV%netcdf%id_var_beta_sq )
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      CALL inquire_double_var( BIV%netcdf%ncid, BIV%netcdf%name_var_phi_fric, (/ BIV%netcdf%id_dim_x, BIV%netcdf%id_dim_y /), BIV%netcdf%id_var_phi_fric)
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      CALL inquire_double_var( BIV%netcdf%ncid, BIV%netcdf%name_var_alpha_sq, (/ BIV%netcdf%id_dim_x, BIV%netcdf%id_dim_y /), BIV%netcdf%id_var_alpha_sq)
      CALL inquire_double_var( BIV%netcdf%ncid, BIV%netcdf%name_var_beta_sq,  (/ BIV%netcdf%id_dim_x, BIV%netcdf%id_dim_y /), BIV%netcdf%id_var_beta_sq )
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      CALL inquire_double_var( BIV%netcdf%ncid, BIV%netcdf%name_var_alpha_sq, (/ BIV%netcdf%id_dim_x, BIV%netcdf%id_dim_y /), BIV%netcdf%id_var_alpha_sq)
      CALL inquire_double_var( BIV%netcdf%ncid, BIV%netcdf%name_var_beta_sq,  (/ BIV%netcdf%id_dim_x, BIV%netcdf%id_dim_y /), BIV%netcdf%id_var_beta_sq )
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      CALL inquire_double_var( BIV%netcdf%ncid, BIV%netcdf%name_var_phi_fric, (/ BIV%netcdf%id_dim_x, BIV%netcdf%id_dim_y /), BIV%netcdf%id_var_phi_fric)
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF
        
    ! Close the netcdf file
    CALL close_netcdf_file( BIV%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_BIV_bed_roughness_file
  SUBROUTINE read_BIV_bed_roughness_file(    BIV)
    ! Read the extrapolated ocean data netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_BIV_bed_roughness),        INTENT(INOUT) :: BIV
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_BIV_bed_roughness_file'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( BIV%netcdf%filename, BIV%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( BIV%netcdf%ncid, BIV%netcdf%id_var_x, BIV%x,  start = (/ 1 /) ))
    CALL handle_error(nf90_get_var( BIV%netcdf%ncid, BIV%netcdf%id_var_y, BIV%y,  start = (/ 1 /) ))
    
    IF     (C%choice_sliding_law == 'no_sliding') THEN
      CALL crash('not defined for choice_sliding_law = "no_sliding"!')
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      CALL handle_error(nf90_get_var( BIV%netcdf%ncid, BIV%netcdf%id_var_beta_sq,  BIV%beta_sq,  start = (/ 1, 1 /) ))
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      CALL handle_error(nf90_get_var( BIV%netcdf%ncid, BIV%netcdf%id_var_phi_fric, BIV%phi_fric, start = (/ 1, 1 /) ))
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      CALL handle_error(nf90_get_var( BIV%netcdf%ncid, BIV%netcdf%id_var_alpha_sq, BIV%alpha_sq, start = (/ 1, 1 /) ))
      CALL handle_error(nf90_get_var( BIV%netcdf%ncid, BIV%netcdf%id_var_beta_sq,  BIV%beta_sq,  start = (/ 1, 1 /) ))
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      CALL handle_error(nf90_get_var( BIV%netcdf%ncid, BIV%netcdf%id_var_alpha_sq, BIV%alpha_sq, start = (/ 1, 1 /) ))
      CALL handle_error(nf90_get_var( BIV%netcdf%ncid, BIV%netcdf%id_var_beta_sq,  BIV%beta_sq,  start = (/ 1, 1 /) ))
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      CALL handle_error(nf90_get_var( BIV%netcdf%ncid, BIV%netcdf%id_var_phi_fric, BIV%phi_fric, start = (/ 1, 1 /) ))
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF
        
    ! Close the netcdf file
    CALL close_netcdf_file( BIV%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_BIV_bed_roughness_file
  
  ! Target velocity fields for basal inversion
  SUBROUTINE inquire_BIV_target_velocity( BIV_target)
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_BIV_target_velocity), INTENT(INOUT) :: BIV_target
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'inquire_BIV_target_velocity'
    INTEGER                                       :: i
    LOGICAL                                       :: is_Rignot2011
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
        
    ! Open the netcdf file
    CALL open_netcdf_file( BIV_target%netcdf%filename, BIV_target%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( BIV_target%netcdf%ncid, BIV_target%netcdf%name_dim_x, BIV_target%nx, BIV_target%netcdf%id_dim_x)
    CALL inquire_dim( BIV_target%netcdf%ncid, BIV_target%netcdf%name_dim_y, BIV_target%ny, BIV_target%netcdf%id_dim_y)
    
    ! Exception for the Rignot et al. (2011) velocity file (downloadable from https: // nsidc.org/data/NSIDC-0484/versions/2)
    is_Rignot2011 = .FALSE.
    DO i = 1, 256-36
      IF (BIV_target%netcdf%filename(i:i+33) == 'antarctica_ice_velocity_450m_v2.nc') THEN
        is_Rignot2011 = .TRUE.
      END IF
    END DO
    
    IF (is_Rignot2011) THEN
      ! Exception for the Rignot et al. (2011) velocity file (downloadable from https: // nsidc.org/data/NSIDC-0484/versions/2)

      BIV_target%netcdf%name_var_u_surf = 'VX'
      BIV_target%netcdf%name_var_v_surf = 'VY'
      
      ! Inquire variable id's. Make sure that each variable has the correct dimensions:
      CALL inquire_double_var( BIV_target%netcdf%ncid, BIV_target%netcdf%name_var_x,      (/ BIV_target%netcdf%id_dim_x                             /), BIV_target%netcdf%id_var_x     )
      CALL inquire_double_var( BIV_target%netcdf%ncid, BIV_target%netcdf%name_var_y,      (/                             BIV_target%netcdf%id_dim_y /), BIV_target%netcdf%id_var_y     )
      
      CALL inquire_single_var( BIV_target%netcdf%ncid, BIV_target%netcdf%name_var_u_surf, (/ BIV_target%netcdf%id_dim_x, BIV_target%netcdf%id_dim_y /), BIV_target%netcdf%id_var_u_surf)
      CALL inquire_single_var( BIV_target%netcdf%ncid, BIV_target%netcdf%name_var_v_surf, (/ BIV_target%netcdf%id_dim_x, BIV_target%netcdf%id_dim_y /), BIV_target%netcdf%id_var_v_surf)
      
    ELSE ! IF (is_Rignot2011) THEN
      ! Generic velocity file
      
      ! Inquire variable id's. Make sure that each variable has the correct dimensions:
      CALL inquire_double_var( BIV_target%netcdf%ncid, BIV_target%netcdf%name_var_x,      (/ BIV_target%netcdf%id_dim_x                             /), BIV_target%netcdf%id_var_x     )
      CALL inquire_double_var( BIV_target%netcdf%ncid, BIV_target%netcdf%name_var_y,      (/                             BIV_target%netcdf%id_dim_y /), BIV_target%netcdf%id_var_y     )
      
      CALL inquire_double_var( BIV_target%netcdf%ncid, BIV_target%netcdf%name_var_u_surf, (/ BIV_target%netcdf%id_dim_x, BIV_target%netcdf%id_dim_y /), BIV_target%netcdf%id_var_u_surf)
      CALL inquire_double_var( BIV_target%netcdf%ncid, BIV_target%netcdf%name_var_v_surf, (/ BIV_target%netcdf%id_dim_x, BIV_target%netcdf%id_dim_y /), BIV_target%netcdf%id_var_v_surf)
      
    END IF ! IF (is_Rignot2011) THEN
        
    ! Close the netcdf file
    CALL close_netcdf_file( BIV_target%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE inquire_BIV_target_velocity
  SUBROUTINE read_BIV_target_velocity(    BIV_target)
    ! Read the target velocity NetCDF file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_BIV_target_velocity), INTENT(INOUT) :: BIV_target
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'read_BIV_target_velocity'
    INTEGER                                       :: i,j
    LOGICAL                                       :: is_Rignot2011
    REAL(dp), DIMENSION(:    ), ALLOCATABLE       :: y_temp
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: u_temp, v_temp
    REAL(dp)                                      :: NaN
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Open the netcdf file
    CALL open_netcdf_file( BIV_target%netcdf%filename, BIV_target%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( BIV_target%netcdf%ncid, BIV_target%netcdf%id_var_x,      BIV_target%x,      start = (/ 1    /) ))
    CALL handle_error(nf90_get_var( BIV_target%netcdf%ncid, BIV_target%netcdf%id_var_y,      BIV_target%y,      start = (/ 1    /) ))
    
    CALL handle_error(nf90_get_var( BIV_target%netcdf%ncid, BIV_target%netcdf%id_var_u_surf, BIV_target%u_surf, start = (/ 1, 1 /) ))
    CALL handle_error(nf90_get_var( BIV_target%netcdf%ncid, BIV_target%netcdf%id_var_v_surf, BIV_target%v_surf, start = (/ 1, 1 /) ))
    
    ! Exception: for some reason, the Rignot 2011 data has the y-axis reversed...
    is_Rignot2011 = .FALSE.
    DO i = 1, 256-36
      IF (BIV_target%netcdf%filename(i:i+33) == 'antarctica_ice_velocity_450m_v2.nc') THEN
        is_Rignot2011 = .TRUE.
      END IF
    END DO
    IF (is_Rignot2011) THEN
      
      ! Allocate temporary memory for storing the upside-down data
      ALLOCATE( y_temp( BIV_target%ny))
      ALLOCATE( u_temp( BIV_target%nx, BIV_target%ny))
      ALLOCATE( v_temp( BIV_target%nx, BIV_target%ny))
      
      ! Copy the upside-down data to temporary memory
      y_temp = BIV_target%y
      u_temp = BIV_target%u_surf
      v_temp = BIV_target%v_surf
      
      ! Flip the data
      DO j = 1, BIV_target%ny
        BIV_target%y(        j) = y_temp(   BIV_target%ny + 1 - j)
        BIV_target%u_surf( :,j) = u_temp( :,BIV_target%ny + 1 - j)
        BIV_target%v_surf( :,j) = v_temp( :,BIV_target%ny + 1 - j)
      END DO
      
      ! Deallocate temporary memory
      DEALLOCATE( y_temp)
      DEALLOCATE( u_temp)
      DEALLOCATE( v_temp)
      
      ! Set missing values to NaN
      NaN = 0._dp
      NaN = 0._dp / NaN
      DO i = 1, BIV_target%nx
      DO j = 1, BIV_target%ny
        IF (BIV_target%u_surf( i,j) == 0._dp .AND. BIV_target%v_surf( i,j) == 0._dp) THEN
          BIV_target%u_surf( i,j) = NaN
          BIV_target%v_surf( i,j) = NaN
        END IF
      END DO
      END DO
      
    END IF ! IF (is_Rignot2011) THEN
        
    ! Close the netcdf file
    CALL close_netcdf_file( BIV_target%netcdf%ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE read_BIV_target_velocity
  
  ! ISMIP6 output
  SUBROUTINE create_ISMIP6_output_files( region)
    ! Create all the ISMIP6 output files
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_model_region), INTENT(IN)    :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_ISMIP6_output_files'
    CHARACTER(LEN=256)                                 :: icesheet_code, foldername
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Code for the ice sheet name in the ISMIP6 output file names
    IF     (region%name == 'NAM') THEN
      icesheet_code = 'NAIS'
    ELSEIF (region%name == 'EAS') THEN
      icesheet_code = 'EUIS'
    ELSEIF (region%name == 'GRL') THEN
      icesheet_code = 'GIS'
    ELSEIF (region%name == 'ANT') THEN
      icesheet_code = 'AIS'
    ELSEIF (region%name == 'PAT') THEN
      icesheet_code = 'PIS'
    ELSE
      icesheet_code = 'beep'
      CALL crash('unknown region "' // TRIM( region%name) // '"!')
    END IF
    
    ! Create a subdirectory within the output directory
    foldername = TRIM( C%output_dir) // TRIM(                 icesheet_code  ) // '_' // &
                                        TRIM( C%ISMIP6_output_group_code     ) // '_' // &
                                        TRIM( C%ISMIP6_output_model_code     ) // '_' // &
                                        TRIM( C%ISMIP6_output_experiment_code)
    CALL system('mkdir ' // foldername)

    ! Create all the ISMIP6 output files
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'lithk'                , 'land_ice_thickness'                                               , 'm'             )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'orog'                 , 'surface_altitude'                                                 , 'm'             )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'topg'                 , 'bedrock_altitude'                                                 , 'm'             )
    CALL create_ISMIP6_output_file_field_notime( foldername, icesheet_code, region%grid, 'hfgeoubed'            , 'upward_geothermal_heat_flux_at_ground_level'                      , 'W m-2'         )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'acabf'                , 'land_ice_surface_specific_mass_balance_flux'                      , 'kg m-2 s-1'    )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'libmassbfgr'          , 'land_ice_basal_specific_mass_balance_flux'                        , 'kg m-2 s-1'    )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'libmassbffl'          , 'land_ice_basal_specific_mass_balance_flux'                        , 'kg m-2 s-1'    )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'dlithkdt'             , 'tendency_of_land_ice_thickness'                                   , 'm s-1'         )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'xvelsurf'             , 'land_ice_surface_x_velocity'                                      , 'm s-1'         )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'yvelsurf'             , 'land_ice_surface_y_velocity'                                      , 'm s-1'         )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'zvelsurf'             , 'land_ice_surface_upward_velocity'                                 , 'm s-1'         )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'xvelbase'             , 'land_ice_basal_x_velocity'                                        , 'm s-1'         )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'yvelbase'             , 'land_ice_basal_y_velocity'                                        , 'm s-1'         )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'zvelbase'             , 'land_ice_basal_upward_velocity'                                   , 'm s-1'         )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'xvelmean'             , 'land_ice_vertical_mean_x_velocity'                                , 'm s-1'         )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'yvelmean'             , 'land_ice_vertical_mean_y_velocity'                                , 'm s-1'         )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'litemptop'            , 'temperature_at_top_of_ice_sheet_model'                            , 'K'             )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'litempbotgr'          , 'temperature_at_base_of_ice_sheet_model'                           , 'K'             )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'litempbotfl'          , 'temperature_at_base_of_ice_sheet_model'                           , 'K'             )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'strbasemag'           , 'land_ice_basal_drag'                                              , 'Pa'            )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'licalvf'              , 'land_ice_specific_mass_flux_due_to_calving'                       , 'kg m-2 s-1'    )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'lifmassbf'            , 'land_ice_specific_mass_flux_due_to_calving_and_ice_front_melting' , 'kg m-2 s-1'    )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'sftgif'               , 'land_ice_area_fraction'                                           , '1'             )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'sftgrf'               , 'grounded_ice_sheet_area_fraction'                                 , '1'             )
    CALL create_ISMIP6_output_file_field(        foldername, icesheet_code, region%grid, 'sftflf'               , 'floating_ice_shelf_area_fraction'                                 , '1'             )
    CALL create_ISMIP6_output_file_scalar(       foldername, icesheet_code,              'lim'                  , 'land_ice_mass'                                                    , 'kg'            )
    CALL create_ISMIP6_output_file_scalar(       foldername, icesheet_code,              'limnsw'               , 'land_ice_mass_not_displacing_sea_water'                           , 'kg'            )
    CALL create_ISMIP6_output_file_scalar(       foldername, icesheet_code,              'iareagr'              , 'grounded_ice_sheet_area'                                          , 'm2'            )
    CALL create_ISMIP6_output_file_scalar(       foldername, icesheet_code,              'iareafl'              , 'floating_ice_sheet_area'                                          , 'm2'            )
    CALL create_ISMIP6_output_file_scalar(       foldername, icesheet_code,              'tendacabf'            , 'tendency_of_land_ice_mass_due_to_surface_mass_balance'            , 'kg s-1'        )
    CALL create_ISMIP6_output_file_scalar(       foldername, icesheet_code,              'tendlibmassbf'        , 'tendency_of_land_ice_mass_due_to_basal_mass_balance'              , 'kg s-1'        )
    CALL create_ISMIP6_output_file_scalar(       foldername, icesheet_code,              'tendlibmassbffl'      , 'tendency_of_land_ice_mass_due_to_basal_mass_balance'              , 'kg s-1'        )
    CALL create_ISMIP6_output_file_scalar(       foldername, icesheet_code,              'tendlicalvf'          , 'tendency_of_land_ice_mass_due_to_calving'                         , 'kg s-1'        )
    CALL create_ISMIP6_output_file_scalar(       foldername, icesheet_code,              'tendlifmassbf'        , 'tendency_of_land_ice_mass_due_to_calving_and_ice_front_melting'   , 'kg s-1'        )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE create_ISMIP6_output_files
  SUBROUTINE create_ISMIP6_output_file_scalar( foldername, icesheet_code, variable_name, standard_name, units)
    ! Create a single ISMIP6 output file for a scalar
    
    IMPLICIT NONE
    
    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: foldername, icesheet_code, variable_name, standard_name, units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_ISMIP6_output_file_scalar'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_t
    INTEGER                                            :: id_var_t, id_var_scalar
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Filename
    filename = TRIM(foldername) // '/' // TRIM( variable_name                  ) // '_' // &
                                          TRIM(                 icesheet_code  ) // '_' // &
                                          TRIM( C%ISMIP6_output_group_code     ) // '_' // &
                                          TRIM( C%ISMIP6_output_model_code     ) // '_' // &
                                          TRIM( C%ISMIP6_output_experiment_code) // '.nc'
    
    ! Check if a file by this name already exists; if so, crash the model.
    INQUIRE(EXIST=file_exists, FILE = TRIM( filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" already exists!')
    END IF
    
    ! Create netcdf file
    CALL handle_error( nf90_create( filename, IOR( nf90_clobber,nf90_share), ncid))
        
    ! Define dimensions:
    CALL create_dim( ncid, 'time',      nf90_unlimited,     id_dim_t      )
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.
    
    ! time variable (needs some attributes that are not in the standard subroutine)
    CALL handle_error( nf90_def_var( ncid, 'time', nf90_float, [id_dim_t], id_var_t))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'standard_name', 'time'))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'long_name'    , 'time'))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'units'        , 'days since ' // TRIM( C%ISMIP6_output_basetime)))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'calendar'     , '360_day'))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'axis'         , 'T'))
    
    ! Field variable (needs some attributes that are not in the standard subroutine)
    CALL handle_error( nf90_def_var( ncid, variable_name, nf90_float, [id_dim_t], id_var_scalar))
    CALL handle_error( nf90_put_att( ncid, id_var_scalar, 'standard_name', standard_name))
    CALL handle_error( nf90_put_att( ncid, id_var_scalar, 'units'        , units))
    CALL handle_error( nf90_put_att( ncid, id_var_scalar, 'missing_value', 1.E20))
    
    ! Leave definition mode:
    CALL handle_error( nf90_enddef( ncid))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( ncid))
    
    ! Close the file
    CALL close_netcdf_file( ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE create_ISMIP6_output_file_scalar
  SUBROUTINE create_ISMIP6_output_file_field( foldername, icesheet_code, grid, variable_name, standard_name, units)
    ! Create a single ISMIP6 output file for an [x,y,t] data field
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: foldername, icesheet_code, variable_name, standard_name, units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_ISMIP6_output_file_field'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_t
    INTEGER                                            :: id_var_x, id_var_y, id_var_t, id_var_field
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Filename
    filename = TRIM(foldername) // '/' // TRIM( variable_name                  ) // '_' // &
                                          TRIM(                 icesheet_code  ) // '_' // &
                                          TRIM( C%ISMIP6_output_group_code     ) // '_' // &
                                          TRIM( C%ISMIP6_output_model_code     ) // '_' // &
                                          TRIM( C%ISMIP6_output_experiment_code) // '.nc'
    
    ! Check if a file by this name already exists; if so, crash the model.
    INQUIRE(EXIST=file_exists, FILE = TRIM( filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" already exists!')
    END IF
    
    ! Create netcdf file
    CALL handle_error( nf90_create( filename, IOR( nf90_clobber,nf90_share), ncid))
        
    ! Define dimensions:
    CALL create_dim( ncid, 'x',         grid%nx,            id_dim_x      )
    CALL create_dim( ncid, 'y',         grid%ny,            id_dim_y      )
    CALL create_dim( ncid, 'time',      nf90_unlimited,     id_dim_t      )
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.
    
    ! x,y variables
    CALL create_single_var( ncid, 'x', [id_dim_x], id_var_x, long_name = 'x-coordinate', units = 'm')
    CALL create_single_var( ncid, 'y', [id_dim_y], id_var_y, long_name = 'y-coordinate', units = 'm')
    
    ! time variable (needs some attributes that are not in the standard subroutine)
    CALL handle_error( nf90_def_var( ncid, 'time', nf90_float, [id_dim_t], id_var_t))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'standard_name', 'time'))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'long_name'    , 'time'))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'units'        , 'days since ' // TRIM( C%ISMIP6_output_basetime)))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'calendar'     , '360_day'))
    CALL handle_error( nf90_put_att( ncid, id_var_t, 'axis'         , 'T'))
    
    ! Field variable (needs some attributes that are not in the standard subroutine)
    CALL handle_error( nf90_def_var( ncid, variable_name, nf90_float, [id_dim_x, id_dim_y, id_dim_t], id_var_field))
    CALL handle_error( nf90_put_att( ncid, id_var_field, 'standard_name', standard_name))
    CALL handle_error( nf90_put_att( ncid, id_var_field, 'units'        , units))
    CALL handle_error( nf90_put_att( ncid, id_var_field, 'missing_value', 1.E20))
    
    ! Leave definition mode:
    CALL handle_error( nf90_enddef( ncid))
    
    ! Write the x, y variable data
    CALL handle_error( nf90_put_var( ncid, id_var_x, grid%x))
    CALL handle_error( nf90_put_var( ncid, id_var_y, grid%y))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( ncid))
    
    ! Close the file
    CALL close_netcdf_file( ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE create_ISMIP6_output_file_field
  SUBROUTINE create_ISMIP6_output_file_field_notime( foldername, icesheet_code, grid, variable_name, standard_name, units)
    ! Create a single ISMIP6 output file for an [x,y] data field
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    CHARACTER(LEN=*),                    INTENT(IN)    :: foldername, icesheet_code, variable_name, standard_name, units

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_ISMIP6_output_file_field_notime'
    CHARACTER(LEN=256)                                 :: filename
    LOGICAL                                            :: file_exists
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_x, id_dim_y
    INTEGER                                            :: id_var_x, id_var_y, id_var_field
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Filename
    filename = TRIM(foldername) // '/' // TRIM( variable_name                  ) // '_' // &
                                          TRIM(                 icesheet_code  ) // '_' // &
                                          TRIM( C%ISMIP6_output_group_code     ) // '_' // &
                                          TRIM( C%ISMIP6_output_model_code     ) // '_' // &
                                          TRIM( C%ISMIP6_output_experiment_code) // '.nc'
    
    ! Check if a file by this name already exists; if so, crash the model.
    INQUIRE(EXIST=file_exists, FILE = TRIM( filename))
    IF (file_exists) THEN
      CALL crash('file "' // TRIM( filename) // '" already exists!')
    END IF
    
    ! Create netcdf file
    CALL handle_error( nf90_create( filename, IOR( nf90_clobber,nf90_share), ncid))
        
    ! Define dimensions:
    CALL create_dim( ncid, 'x',         grid%nx,            id_dim_x      )
    CALL create_dim( ncid, 'y',         grid%ny,            id_dim_y      )
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.
    
    ! x,y variables
    CALL create_single_var( ncid, 'x', [id_dim_x], id_var_x, long_name = 'x-coordinate', units = 'm')
    CALL create_single_var( ncid, 'y', [id_dim_y], id_var_y, long_name = 'y-coordinate', units = 'm')
    
    ! Field variable (needs some attributes that are not in the standard subroutine)
    CALL handle_error( nf90_def_var( ncid, variable_name, nf90_float, [id_dim_x, id_dim_y], id_var_field))
    CALL handle_error( nf90_put_att( ncid, id_var_field, 'standard_name', standard_name))
    CALL handle_error( nf90_put_att( ncid, id_var_field, 'units'        , units))
    CALL handle_error( nf90_put_att( ncid, id_var_field, 'missing_value', 1.E20))
    
    ! Leave definition mode:
    CALL handle_error( nf90_enddef( ncid))
    
    ! Write the x, y variable data
    CALL handle_error( nf90_put_var( ncid, id_var_x, grid%x))
    CALL handle_error( nf90_put_var( ncid, id_var_y, grid%y))
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( ncid))
    
    ! Close the file
    CALL close_netcdf_file( ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE create_ISMIP6_output_file_field_notime
  SUBROUTINE write_to_ISMIP6_output_files( region)
    ! Write to all the ISMIP6 output files
  
    USE parameters_module, ONLY: ice_density, sec_per_year
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_model_region), INTENT(IN)    :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                        :: routine_name = 'write_to_ISMIP6_output_files'
    INTEGER                                              :: i,j
    CHARACTER(LEN=256)                                   :: icesheet_code, foldername
    REAL(dp), PARAMETER                                  :: missing_value = 1.E20
    REAL(dp), DIMENSION( region%grid%ny, region%grid%nx) :: Ti_base_gr
    REAL(dp), DIMENSION( region%grid%ny, region%grid%nx) :: Ti_base_fl
    REAL(dp), DIMENSION( region%grid%ny, region%grid%nx) :: basal_drag
    REAL(dp), DIMENSION( region%grid%ny, region%grid%nx) :: calving_flux
    REAL(dp), DIMENSION( region%grid%ny, region%grid%nx) :: calving_and_front_melt_flux
    REAL(dp), DIMENSION( region%grid%ny, region%grid%nx) :: land_ice_area_fraction
    REAL(dp), DIMENSION( region%grid%ny, region%grid%nx) :: grounded_ice_sheet_area_fraction
    REAL(dp), DIMENSION( region%grid%ny, region%grid%nx) :: floating_ice_shelf_area_fraction
    REAL(dp)                                             :: land_ice_mass
    REAL(dp)                                             :: mass_above_floatation
    REAL(dp)                                             :: grounded_ice_sheet_area
    REAL(dp)                                             :: floating_ice_sheet_area
    REAL(dp)                                             :: total_SMB
    REAL(dp)                                             :: total_BMB
    REAL(dp)                                             :: total_BMB_shelf
    REAL(dp)                                             :: total_calving_flux
    REAL(dp)                                             :: total_calving_and_front_melt_flux
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Code for the ice sheet name in the ISMIP6 output file names
    IF     (region%name == 'NAM') THEN
      icesheet_code = 'NAIS'
    ELSEIF (region%name == 'EAS') THEN
      icesheet_code = 'EUIS'
    ELSEIF (region%name == 'GRL') THEN
      icesheet_code = 'GIS'
    ELSEIF (region%name == 'ANT') THEN
      icesheet_code = 'AIS'
    ELSEIF (region%name == 'PAT') THEN
      icesheet_code = 'PIS'
    ELSE
      icesheet_code = 'beep'
      CALL crash('unknown region "' // TRIM( region%name) // '"!')
    END IF
    
    ! The folder where the ISMIp6 output files are located
    foldername = TRIM( C%output_dir) // TRIM(                 icesheet_code  ) // '_' // &
                                        TRIM( C%ISMIP6_output_group_code     ) // '_' // &
                                        TRIM( C%ISMIP6_output_model_code     ) // '_' // &
                                        TRIM( C%ISMIP6_output_experiment_code)
                                            
    ! Calculate some quantities that are not natively in the ice model
    
    ! 2-D fields
    Ti_base_gr                        = missing_value
    Ti_base_fl                        = missing_value
    basal_drag                        = missing_value
    calving_flux                      = missing_value
    calving_and_front_melt_flux       = missing_value
    land_ice_area_fraction            = missing_value
    grounded_ice_sheet_area_fraction  = missing_value
    floating_ice_shelf_area_fraction  = missing_value
    
    ! Scalars (integrated values)
    land_ice_mass                     = 0._dp
    mass_above_floatation             = 0._dp
    grounded_ice_sheet_area           = 0._dp
    floating_ice_sheet_area           = 0._dp
    total_SMB                         = 0._dp
    total_BMB                         = 0._dp
    total_BMB_shelf                   = 0._dp
    total_calving_flux                = 0._dp
    total_calving_and_front_melt_flux = 0._dp
    
    DO i = 1, region%grid%nx
    DO j = 1, region%grid%ny
      
      ! Ice base temperature separate for sheet and shelf
      ! =================================================
      
      IF (region%ice%mask_sheet_a( j,i) == 1) THEN
        Ti_base_gr( j,i) = region%ice%Ti_a( C%nz,j,i)
      END IF
      
      IF (region%ice%mask_shelf_a( j,i) == 1) THEN
        Ti_base_fl( j,i) = region%ice%Ti_a( C%nz,j,i)
      END IF
      
      ! Basal drag
      ! ==========
      
      IF (region%ice%mask_ice_a( j,i) == 1 .AND. region%ice%f_grnd_a( j,i) > 0._dp) THEN
        basal_drag( j,i) = region%ice%uabs_base_a( j,i) * region%ice%beta_a( j,i) * region%ice%f_grnd_a( j,i)**2
      END IF
      
      ! Calving and front melting fluxes
      ! ================================
      
      calving_flux(                j,i) = 0._dp ! FIXME
      calving_and_front_melt_flux( j,i) = 0._dp ! FIXME
      
      ! Ice fractions
      ! =============
      
      IF (region%ice%mask_cf_a( j,i) == 0) THEN
        land_ice_area_fraction( j,i) = REAL( region%ice%mask_ice_a( j,i), dp)
      ELSE
        land_ice_area_fraction( j,i) = region%ice%float_margin_frac_a( j,i)
      END IF
      
      IF (region%ice%mask_ice_a( j,i) == 1) THEN
        grounded_ice_sheet_area_fraction( j,i) = region%ice%f_grnd_a( j,i)
      ELSE
        grounded_ice_sheet_area_fraction( j,i) = 0._dp
      END IF
      
      floating_ice_shelf_area_fraction( j,i) = REAL( region%ice%mask_ice_a( j,i), dp) * MAX( (1._dp - region%ice%f_grnd_a( j,i)), region%ice%float_margin_frac_a( j,i))
      
      ! Integrated values
      ! =================
      
      land_ice_mass                     = land_ice_mass                     + (region%ice%Hi_a(                  j,i) * region%grid%dx**2 * ice_density)    ! kg
      mass_above_floatation             = mass_above_floatation             + (region%ice%TAF_a(                 j,i) * region%grid%dx**2 * ice_density)    ! kg
      grounded_ice_sheet_area           = grounded_ice_sheet_area           + (grounded_ice_sheet_area_fraction( j,i) * region%grid%dx**2)                  ! m2
      floating_ice_sheet_area           = floating_ice_sheet_area           + (floating_ice_shelf_area_fraction( j,i) * region%grid%dx**2)                  ! m2
      total_SMB                         = total_SMB                         + (land_ice_area_fraction( j,i) * region%SMB%SMB_year(  j,i) * region%grid%dx**2 * ice_density / sec_per_year) ! kg s-1
      total_BMB                         = total_BMB                         + (land_ice_area_fraction( j,i) * region%BMB%BMB(       j,i) * region%grid%dx**2 * ice_density / sec_per_year) ! kg s-1
      total_BMB_shelf                   = total_BMB_shelf                   + (land_ice_area_fraction( j,i) * region%BMB%BMB_shelf( j,i) * region%grid%dx**2 * ice_density / sec_per_year) ! kg s-1
      total_calving_flux                = total_calving_flux                + (calving_flux( j,i)                                        * region%grid%dx**2 * ice_density / sec_per_year) ! kg s-1
      total_calving_and_front_melt_flux = total_calving_and_front_melt_flux + (calving_and_front_melt_flux( j,i)                         * region%grid%dx**2 * ice_density / sec_per_year) ! kg s-1
      
    END DO
    END DO
    
    ! Write to all the ISMIP6 output files
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%Hi_a                    , 'lithk'                    )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%Hs_a                    , 'orog'                     )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%Hb_a                    , 'topg'                     )
    CALL write_to_ISMIP6_output_file_field_notime(  foldername, icesheet_code,              region%ice%GHF_a                   , 'hfgeoubed'                )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%SMB%SMB_year                , 'acabf'                    )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%BMB%BMB_sheet               , 'libmassbfgr'              )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%BMB%BMB_shelf               , 'libmassbffl'              )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%dHs_dt_a                , 'dlithkdt'                 )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%u_surf_a                , 'xvelsurf'                 )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%v_surf_a                , 'yvelsurf'                 )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%w_surf_a                , 'zvelsurf'                 )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%u_base_a                , 'xvelbase'                 )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%v_base_a                , 'yvelbase'                 )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%w_base_a                , 'zvelbase'                 )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%u_vav_a                 , 'xvelmean'                 )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%v_vav_a                 , 'yvelmean'                 )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, region%ice%Ti_a( 1   ,:,:)         , 'litemptop'                )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, Ti_base_gr                         , 'litempbotgr'              )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, Ti_base_fl                         , 'litempbotfl'              )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, basal_drag                         , 'strbasemag'               )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, calving_flux                       , 'licalvf'                  )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, calving_and_front_melt_flux        , 'lifmassbf'                )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, land_ice_area_fraction             , 'sftgif'                   )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, grounded_ice_sheet_area_fraction   , 'sftgrf'                   )
    CALL write_to_ISMIP6_output_file_field(         foldername, icesheet_code, region%time, floating_ice_shelf_area_fraction   , 'sftflf'                   )
    CALL write_to_ISMIP6_output_file_scalar(        foldername, icesheet_code, region%time, land_ice_mass                      , 'lim'                      )
    CALL write_to_ISMIP6_output_file_scalar(        foldername, icesheet_code, region%time, mass_above_floatation              , 'limnsw'                   )
    CALL write_to_ISMIP6_output_file_scalar(        foldername, icesheet_code, region%time, grounded_ice_sheet_area            , 'iareagr'                  )
    CALL write_to_ISMIP6_output_file_scalar(        foldername, icesheet_code, region%time, floating_ice_sheet_area            , 'iareafl'                  )
    CALL write_to_ISMIP6_output_file_scalar(        foldername, icesheet_code, region%time, total_SMB                          , 'tendacabf'                )
    CALL write_to_ISMIP6_output_file_scalar(        foldername, icesheet_code, region%time, total_BMB                          , 'tendlibmassbf'            )
    CALL write_to_ISMIP6_output_file_scalar(        foldername, icesheet_code, region%time, total_BMB_shelf                    , 'tendlibmassbffl'          )
    CALL write_to_ISMIP6_output_file_scalar(        foldername, icesheet_code, region%time, total_calving_flux                 , 'tendlicalvf'              )
    CALL write_to_ISMIP6_output_file_scalar(        foldername, icesheet_code, region%time, total_calving_and_front_melt_flux  , 'tendlifmassbf'            )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE write_to_ISMIP6_output_files
  SUBROUTINE write_to_ISMIP6_output_file_scalar( foldername, icesheet_code, time, d, variable_name)
    ! Write a single scalar to the corresponding ISMIP6 output file
   
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time
    REAL(dp),                            INTENT(IN)    :: d
    CHARACTER(LEN=*),                    INTENT(IN)    :: foldername, icesheet_code, variable_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_ISMIP6_output_file_scalar'
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_t
    INTEGER                                            :: id_var_t, id_var
    INTEGER                                            :: nt
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Filename
    filename = TRIM(foldername) // '/' // TRIM( variable_name                  ) // '_' // &
                                          TRIM(                 icesheet_code  ) // '_' // &
                                          TRIM( C%ISMIP6_output_group_code     ) // '_' // &
                                          TRIM( C%ISMIP6_output_model_code     ) // '_' // &
                                          TRIM( C%ISMIP6_output_experiment_code) // '.nc'
        
    ! Open the netcdf file
    CALL open_netcdf_file( filename, ncid)
 
    ! Inquire for dimension IDs
    CALL inquire_dim( ncid, 'time', nt, id_dim_t)
    
    ! Inquire for variable IDs
    CALL inquire_single_var( ncid, 'time'       , (/ id_dim_t /), id_var_t)
    CALL inquire_single_var( ncid, variable_name, (/ id_dim_t /), id_var  )
        
    ! Write time
    CALL handle_error( nf90_put_var( ncid, id_var_t, days_since_ISMIP6_basetime( time), start = (/ nt+1 /)))
      
    ! Write data to the NetCDF file
    CALL handle_error( nf90_put_var( ncid, id_var, d, start = (/ nt+1 /)))
    
    ! Close the file
    CALL close_netcdf_file( ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE write_to_ISMIP6_output_file_scalar
  SUBROUTINE write_to_ISMIP6_output_file_field( foldername, icesheet_code, time, d, variable_name)
    ! Write a single [x,y,t] data field to the corresponding ISMIP6 output file
   
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    CHARACTER(LEN=*),                    INTENT(IN)    :: foldername, icesheet_code, variable_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_ISMIP6_output_file_field'
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_x, id_dim_y, id_dim_t
    INTEGER                                            :: id_var_t, id_var
    INTEGER                                            :: nx, ny, nt
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Filename
    filename = TRIM(foldername) // '/' // TRIM( variable_name                  ) // '_' // &
                                          TRIM(                 icesheet_code  ) // '_' // &
                                          TRIM( C%ISMIP6_output_group_code     ) // '_' // &
                                          TRIM( C%ISMIP6_output_model_code     ) // '_' // &
                                          TRIM( C%ISMIP6_output_experiment_code) // '.nc'
        
    ! Open the netcdf file
    CALL open_netcdf_file( filename, ncid)
 
    ! Inquire for dimension IDs
    CALL inquire_dim( ncid, 'x'   , nx, id_dim_x)
    CALL inquire_dim( ncid, 'y'   , ny, id_dim_y)
    CALL inquire_dim( ncid, 'time', nt, id_dim_t)
    
    ! Inquire for variable IDs
    CALL inquire_single_var( ncid, 'time'       , (/                     id_dim_t /), id_var_t)
    CALL inquire_single_var( ncid, variable_name, (/ id_dim_x, id_dim_y, id_dim_t /), id_var  )
        
    ! Write time
    CALL handle_error( nf90_put_var( ncid, id_var_t, days_since_ISMIP6_basetime( time), start = (/ nt+1 /)))
      
    ! Write data to the NetCDF file
    CALL write_data_to_file_dp_2D( ncid, nx, ny, id_var, ISMIP6_unit_conversion_field( d, variable_name), (/ 1, 1, nt+1 /) )
    
    ! Close the file
    CALL close_netcdf_file( ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE write_to_ISMIP6_output_file_field
  SUBROUTINE write_to_ISMIP6_output_file_field_notime( foldername, icesheet_code, d, variable_name)
    ! Write a single [x,y,t] data field to the corresponding ISMIP6 output file
   
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    CHARACTER(LEN=*),                    INTENT(IN)    :: foldername, icesheet_code, variable_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_to_ISMIP6_output_file_field_notime'
    CHARACTER(LEN=256)                                 :: filename
    INTEGER                                            :: ncid
    INTEGER                                            :: id_dim_x, id_dim_y
    INTEGER                                            :: id_var
    INTEGER                                            :: nx, ny
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (.NOT. par%master) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF
    
    ! Filename
    filename = TRIM(foldername) // '/' // TRIM( variable_name                  ) // '_' // &
                                          TRIM(                 icesheet_code  ) // '_' // &
                                          TRIM( C%ISMIP6_output_group_code     ) // '_' // &
                                          TRIM( C%ISMIP6_output_model_code     ) // '_' // &
                                          TRIM( C%ISMIP6_output_experiment_code) // '.nc'
        
    ! Open the netcdf file
    CALL open_netcdf_file( filename, ncid)
 
    ! Inquire for dimension IDs
    CALL inquire_dim( ncid, 'x'   , nx, id_dim_x)
    CALL inquire_dim( ncid, 'y'   , ny, id_dim_y)
    
    ! Inquire for variable IDs
    CALL inquire_single_var( ncid, variable_name, (/ id_dim_x, id_dim_y /), id_var  )
      
    ! Write data to the NetCDF file
    CALL write_data_to_file_dp_2D( ncid, nx, ny, id_var, ISMIP6_unit_conversion_field( d, variable_name), (/ 1, 1 /) )
    
    ! Close the file
    CALL close_netcdf_file( ncid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE write_to_ISMIP6_output_file_field_notime
  FUNCTION ISMIP6_unit_conversion_field( d, variable_name) RESULT( d_conv)
    ! Convert data fields from IMAU-ICE units to SI units
  
    USE parameters_module, ONLY: sec_per_year, ice_density
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    CHARACTER(LEN=*),                    INTENT(IN)    :: variable_name
    REAL(dp), DIMENSION(SIZE(d,1),SIZE(d,2))           :: d_conv
    
    IF     (variable_name == 'lithk') THEN
      ! land_ice_thickness
      ! Ice model units: m
      ! SI units: m
      d_conv = d
    ELSEIF (variable_name == 'orog') THEN
      ! surface_altitude
      ! Ice model units: m
      ! SI units: m
      d_conv = d
    ELSEIF (variable_name == 'topg') THEN
      ! bedrock_altitude
      ! Ice model units: m
      ! SI units: m
      d_conv = d
    ELSEIF (variable_name == 'hfgeoubed') THEN
      ! upward_geothermal_heat_flux_at_ground_level
      ! Ice model units: J m-2 yr-1
      ! SI units: J m-2 s-1
      d_conv = d / sec_per_year
    ELSEIF (variable_name == 'acabf') THEN
      ! land_ice_surface_specific_mass_balance_flux
      ! Ice model units: m.i.e./yr
      ! SI units: kg m-2 s-1
      d_conv = d * ice_density / sec_per_year
    ELSEIF (variable_name == 'libmassbfgr') THEN
      ! land_ice_basal_specific_mass_balance_flux
      ! Ice model units: m.i.e./yr
      ! SI units: kg m-2 s-1
      d_conv = d * ice_density / sec_per_year
    ELSEIF (variable_name == 'libmassbffl') THEN
      ! land_ice_basal_specific_mass_balance_flux
      ! Ice model units: m.i.e./yr
      ! SI units: kg m-2 s-1
      d_conv = d * ice_density / sec_per_year
    ELSEIF (variable_name == 'dlithkdt') THEN
      ! tendency_of_land_ice_thickness
      ! Ice model units: m/yr
      ! SI units: m s-1
      d_conv = d / sec_per_year
    ELSEIF (variable_name == 'xvelsurf' .OR. &
            variable_name == 'yvelsurf' .OR. &
            variable_name == 'zvelsurf' .OR. &
            variable_name == 'xvelbase' .OR. &
            variable_name == 'yvelbase' .OR. &
            variable_name == 'zvelbase' .OR. &
            variable_name == 'xvelmean' .OR. &
            variable_name == 'yvelmean') THEN
      ! different ice velocities
      ! Ice model units: m/yr
      ! SI units: m s-1
      d_conv = d / sec_per_year
    ELSEIF (variable_name == 'litemptop') THEN
      ! temperature_at_top_of_ice_sheet_model
      ! Ice model units: K
      ! SI units: K
      d_conv = d
    ELSEIF (variable_name == 'litempbotgr') THEN
      ! temperature_at_base_of_ice_sheet_model
      ! Ice model units: K
      ! SI units: K
      d_conv = d
    ELSEIF (variable_name == 'litempbotfl') THEN
      ! temperature_at_base_of_ice_sheet_model
      ! Ice model units: K
      ! SI units: K
      d_conv = d
    ELSEIF (variable_name == 'strbasemag') THEN
      ! land_ice_basal_drag
      ! Ice model units: Pa
      ! SI units: Pa
      d_conv = d
    ELSEIF (variable_name == 'licalvf' .OR. &
            variable_name == 'lifmassbf') THEN
      ! land_ice_specific_mass_flux_due_to_calving (possibly also front melting)
      ! Ice model units: m/yr
      ! SI units: kg m-2 s-1
      d_conv = d * ice_density / sec_per_year
    ELSEIF (variable_name == 'sftgif' .OR. &
            variable_name == 'sftgrf' .OR. &
            variable_name == 'sftflf') THEN
      ! ice fractions
      ! Ice model units:
      ! SI units:
      d_conv = d
    ELSE
      ! Unknown variable name
      CALL crash('ISMIP6_unit_conversion: unknown variable name "' // TRIM( variable_name) // '"!')
    END IF
    
  END FUNCTION ISMIP6_unit_conversion_field
  FUNCTION days_since_ISMIP6_basetime( time) RESULT( ndays)
    ! Calculate the number of days since ISMIP6 basetime
    !
    ! Assume basetime equals t = 0
  
    USE parameters_module, ONLY: sec_per_year
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time
    REAL(dp)                                           :: ndays
    
    ndays = time * 360._dp
    
  END FUNCTION days_since_ISMIP6_basetime
  
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
    
END MODULE netcdf_module
