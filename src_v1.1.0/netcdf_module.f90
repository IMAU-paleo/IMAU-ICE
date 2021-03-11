MODULE netcdf_module
  ! This module defines methods to read and write the current state of the model from or to a NetCDF file.

  USE mpi
  USE parallel_module,                 ONLY: par, sync
  USE netcdf,                          ONLY: nf90_max_var_dims, nf90_create, nf90_close, nf90_clobber, nf90_share, nf90_unlimited , &
                                             nf90_enddef, nf90_put_var, nf90_sync, nf90_def_var, nf90_int, nf90_put_att, nf90_def_dim, &
                                             nf90_open, nf90_write, nf90_inq_dimid, nf90_inquire_dimension, nf90_inquire, nf90_double, &
                                             nf90_inq_varid, nf90_inquire_variable, nf90_get_var, nf90_noerr, nf90_strerror
  USE configuration_module,            ONLY: dp, C
  USE data_types_module,               ONLY: type_model_region, type_PD_data_fields, type_forcing_data, &
                                             type_init_data_fields, type_subclimate_global, type_netcdf_output
  USE parameters_module,               ONLY: sec_per_year
  
  IMPLICIT NONE

CONTAINS

  ! Create and write to an output NetCDF file
  SUBROUTINE write_to_restart_file(     region, forcing)
    ! Write the current model state to the existing output file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_forcing_data),        INTENT(IN)    :: forcing
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: ncid, nx, ny, nz, ti
    
    WRITE(0,'(A,F8.2,A)') '   t = ', region%time/1E3, ' kyr - writing output...'
    
    ! Open the file for writing
    CALL open_netcdf_file( region%restart%filename, region%restart%ncid)
        
    ! Time
    CALL handle_error( nf90_put_var( region%restart%ncid, region%restart%id_var_time,             region%time,                    start=(/       region%restart%ti/)))
    
    ! Key output
    ! ==========
    
    ! Placeholders for the dimension ID's, for shorter code
    ncid = region%restart%ncid
    nx   = region%grid%nx
    ny   = region%grid%ny
    nz   = C%nZ
    ti   = region%restart%ti
    
    CALL write_data_to_file_2D( ncid, nx, ny,     region%restart%id_var_Hi,               region%ice%Hi_Aa,            (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%restart%id_var_Hb,               region%ice%Hb_Aa,            (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%restart%id_var_Hs,               region%ice%Hs_Aa,            (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%restart%id_var_U_SIA,            region%ice%U_vav_SIA_Aa,     (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%restart%id_var_V_SIA,            region%ice%V_vav_SIA_Aa,     (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%restart%id_var_U_SSA,            region%ice%U_vav_SSA_Aa,     (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%restart%id_var_V_SSA,            region%ice%V_vav_SSA_Aa,     (/1, 1,    ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nz, region%restart%id_var_Ti,               region%ice%Ti_Aa,            (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, 12, region%restart%id_var_FirnDepth,        region%SMB%FirnDepth,        (/1, 1, 1, ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%restart%id_var_MeltPreviousYear, region%SMB%MeltPreviousYear, (/1, 1,    ti/))
    
    ! Inverse routine data
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod') THEN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in write_to_restart_file!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSE
    
      IF     (C%choice_forcing_method == 'CO2_direct') THEN
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
    
    ! Local variables:
    INTEGER                                :: ncid, nx, ny, ti, nz
    
    ! Open the file for writing
    CALL open_netcdf_file( region%help_fields%filename, region%help_fields%ncid)
        
    ! Time
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_time,             region%time,                    start=(/       region%help_fields%ti/)))
    
    ! Key output
    ! ==========
    
    ! Placeholders for the dimension ID's, for shorter code
    ncid = region%help_fields%ncid
    nx   = region%grid%nx
    ny   = region%grid%ny
    nz   = C%nZ
    ti   = region%help_fields%ti
    
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_mask,        REAL(region%ice%mask_Aa,dp),          (/1, 1,    ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, 12, region%help_fields%id_var_T2m,              region%climate%applied%T2m,      (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, 12, region%help_fields%id_var_Precip,           region%climate%applied%Precip,   (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, 12, region%help_fields%id_var_Albedo,           region%SMB%Albedo,               (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, 12, region%help_fields%id_var_SMB,              region%SMB%SMB,                  (/1, 1, 1, ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_BMB,              region%BMB%BMB,                  (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_dHs_dx,           region%ice%dHs_dx_Aa,            (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_dHs_dy,           region%ice%dHs_dy_Aa,            (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_D_SIA,            region%ice%D_SIA_Aa,             (/1, 1,    ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_U_3D,             region%ice%U_3D_Aa,              (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_V_3D,             region%ice%V_3D_Aa,              (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_W_3D,             region%ice%W_3D_Aa,              (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_A_flow,           region%ice%A_flow_Aa,            (/1, 1, 1, ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_A_flow_mean,      region%ice%A_flow_mean_Aa,       (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_Fr_Aa,            region%ice%Fr_Aa,                (/1, 1,    ti/))
    
    IF (C%do_write_dummy_output) THEN
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_dummy2D_01,       region%ice%dummy2D_01,           (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_dummy2D_02,       region%ice%dummy2D_02,           (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_dummy2D_03,       region%ice%dummy2D_03,           (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_dummy2D_04,       region%ice%dummy2D_04,           (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_dummy2D_05,       region%ice%dummy2D_05,           (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_dummy2D_06,       region%ice%dummy2D_06,           (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_dummy2D_07,       region%ice%dummy2D_07,           (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_dummy2D_08,       region%ice%dummy2D_08,           (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_dummy2D_09,       region%ice%dummy2D_09,           (/1, 1,    ti/))
    CALL write_data_to_file_2D( ncid, nx, ny,     region%help_fields%id_var_dummy2D_10,       region%ice%dummy2D_10,           (/1, 1,    ti/))
    
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_dummy3D_01,       region%ice%dummy3D_01,           (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_dummy3D_02,       region%ice%dummy3D_02,           (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_dummy3D_03,       region%ice%dummy3D_03,           (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_dummy3D_04,       region%ice%dummy3D_04,           (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_dummy3D_05,       region%ice%dummy3D_05,           (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_dummy3D_06,       region%ice%dummy3D_06,           (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_dummy3D_07,       region%ice%dummy3D_07,           (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_dummy3D_08,       region%ice%dummy3D_08,           (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_dummy3D_09,       region%ice%dummy3D_09,           (/1, 1, 1, ti/))
    CALL write_data_to_file_3D( ncid, nx, ny, nZ, region%help_fields%id_var_dummy3D_10,       region%ice%dummy3D_10,           (/1, 1, 1, ti/))
    END IF ! IF (C%do_write_dummy_output) THEN

    ! Close the file
    CALL close_netcdf_file(region%help_fields%ncid)
    
    ! Increase time frame counter
    region%help_fields%ti = region%help_fields%ti + 1
        
  END SUBROUTINE write_to_help_fields_file
  
  SUBROUTINE create_restart_file(       region, forcing)
    ! Create a new restart file, containing the key model output required to start another run.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_forcing_data),        INTENT(IN)    :: forcing

    ! Local variables:
    INTEGER                                       :: cerr, ierr
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x, y, z, m, t
    
    ! Set time frame index to 1
    region%restart%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(region%restart%filename))
    IF(file_exists) THEN
      WRITE(0,'(5A)') 'ERROR: ', TRIM(region%restart%filename), ' already exists!'
      STOP
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

    ! Primary output - everything that's needed to restart a new run, and to plot/analyse data
    ! ========================================================================================
    
    ! Dimension variables: zeta, month, time
    CALL create_double_var( region%restart%ncid, region%restart%name_var_x,                [x            ], region%restart%id_var_x,                long_name='X-coordinate', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_y,                [   y         ], region%restart%id_var_y,                long_name='Y-coordinate', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_zeta,             [      z      ], region%restart%id_var_zeta,             long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_month,            [         m   ], region%restart%id_var_month,            long_name='Month', units='1-12'    )
    CALL create_double_var( region%restart%ncid, region%restart%name_var_time,             [            t], region%restart%id_var_time,             long_name='Time', units='years'   )
    
    ! Ice model data
    CALL create_double_var( region%restart%ncid, region%restart%name_var_Hi,               [x, y,       t], region%restart%id_var_Hi,               long_name='Ice Thickness', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_Hb,               [x, y,       t], region%restart%id_var_Hb,               long_name='Bedrock Height', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_Hs,               [x, y,       t], region%restart%id_var_Hs,               long_name='Surface Height', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_U_SIA,            [x, y,       t], region%restart%id_var_U_SIA,            long_name='SIA ice x-velocity', units='m/yr')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_V_SIA,            [x, y,       t], region%restart%id_var_V_SIA,            long_name='SIA ice y-velocity', units='m/yr') 
    CALL create_double_var( region%restart%ncid, region%restart%name_var_U_SSA,            [x, y,       t], region%restart%id_var_U_SSA,            long_name='SSA ice x-velocity', units='m/yr')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_V_SSA,            [x, y,       t], region%restart%id_var_V_SSA,            long_name='SSA ice y-velocity', units='m/yr')              
    CALL create_double_var( region%restart%ncid, region%restart%name_var_Ti,               [x, y, z,    t], region%restart%id_var_Ti,               long_name='Ice temperature', units='K')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_FirnDepth,        [x, y,    m, t], region%restart%id_var_FirnDepth,        long_name='Firn depth', units='m')
    CALL create_double_var( region%restart%ncid, region%restart%name_var_MeltPreviousYear, [x, y,       t], region%restart%id_var_MeltPreviousYear, long_name='Melt during previous year', units='mie')
    
    ! Inverse routine data
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod') THEN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in create_restart_file!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    ELSE
    
      IF     (C%choice_forcing_method == 'CO2_direct') THEN
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
        WRITE(0,*) '  ERROR - choice_forcing_method "', C%choice_forcing_method, '" not implemented in create_restart_file!'
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
  SUBROUTINE create_help_fields_file(   region)
    ! Create a new help fields file, containing secondary model output (not needed for a restart, but interesting to look at)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region),        INTENT(INOUT) :: region

    ! Local variables:
    LOGICAL                                       :: file_exists
    INTEGER                                       :: x, y, z, m, t
    
    ! Set time frame index to 1
    region%help_fields%ti = 1

    ! Create a new restart file if none exists and, to prevent loss of data, 
    ! stop with an error message if one already exists (not when differences are considered):
    INQUIRE(EXIST=file_exists, FILE = TRIM(region%help_fields%filename))
    IF(file_exists) THEN
      WRITE(0,'(5A)') 'ERROR: ', TRIM(region%help_fields%filename), ' already exists!'
      STOP
    END IF
    
    ! Create netCDF file
    !WRITE(0,*) ' Creating new NetCDF output file at ', TRIM( region%help_fields%filename)
    CALL handle_error(nf90_create(region%help_fields%filename,IOR(nf90_clobber,nf90_share),region%help_fields%ncid))
        
    ! Define dimensions:
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_x,         region%grid%nx,     region%help_fields%id_dim_x         )
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_y,         region%grid%ny,     region%help_fields%id_dim_y         )
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_zeta,      C%nZ,               region%help_fields%id_dim_zeta      ) ! Scaled vertical coordinate
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_month,     12,                 region%help_fields%id_dim_month     ) ! Months (for monthly data)
    CALL create_dim( region%help_fields%ncid, region%help_fields%name_dim_time,      nf90_unlimited,     region%help_fields%id_dim_time      ) ! Time frames
    
    ! Placeholders for the dimension ID's, for shorter code
    x = region%help_fields%id_dim_x
    y = region%help_fields%id_dim_y
    z = region%help_fields%id_dim_zeta
    m = region%help_fields%id_dim_month
    t = region%help_fields%id_dim_time
    
    ! Define variables:
    ! The order of the CALL statements for the different variables determines their
    ! order of appearence in the netcdf file.    

    ! Primary output - everything that's needed to restart a new run, and to plot/analyse data
    ! ========================================================================================
    
    ! Dimension variables: zeta, month, time
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_x,      [x            ], region%help_fields%id_var_x,      long_name='X-coordinate', units='m')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_y,      [   y         ], region%help_fields%id_var_y,      long_name='Y-coordinate', units='m')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_zeta,   [      z      ], region%help_fields%id_var_zeta,   long_name='Vertical scaled coordinate', units='unitless')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_month,  [         m   ], region%help_fields%id_var_month,  long_name='Month', units='1-12'    )
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_time,   [            t], region%help_fields%id_var_time,   long_name='Time', units='years'   )
    
    ! Ice model data
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_lat,         [x, y         ], region%help_fields%id_var_lat,         long_name='Latitude',  units='degrees north')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_lon,         [x, y         ], region%help_fields%id_var_lon,         long_name='Longitude', units='degrees east' )
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_mask,        [x, y,       t], region%help_fields%id_var_mask,        long_name='mask', units='unitless')   
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_T2m,         [x, y,    m, t], region%help_fields%id_var_T2m,         long_name='Monthly mean 2-, air temperature', units='K')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_Precip,      [x, y,    m, t], region%help_fields%id_var_Precip,      long_name='Monthly total precipitation', units='m')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_Albedo,      [x, y,    m, t], region%help_fields%id_var_Albedo,      long_name='Surface albedo', units='unitless')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_SMB,         [x, y,    m, t], region%help_fields%id_var_SMB,         long_name='Surface mass balance', units='mie/yr')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_BMB,         [x, y,       t], region%help_fields%id_var_BMB,         long_name='Basal mass balance', units='mie/yr')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dHs_dx,      [x, y,       t], region%help_fields%id_var_dHs_dx,      long_name='Surface slope in x direction', units='m/m')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dHs_dy,      [x, y,       t], region%help_fields%id_var_dHs_dy,      long_name='Surface slope in y direction', units='m/m') 
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_D_SIA,       [x, y,       t], region%help_fields%id_var_D_SIA,       long_name='SIA ice diffusivity', units='m/yr')  
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_U_3D,        [x, y, z,    t], region%help_fields%id_var_U_3D,        long_name='3D SIA x ice velocity (m/yr)') 
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_V_3D,        [x, y, z,    t], region%help_fields%id_var_V_3D,        long_name='3D SIA y ice velocity (m/yr)') 
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_W_3D,        [x, y, z,    t], region%help_fields%id_var_W_3D,        long_name='3D SIA z ice velocity (m/yr)')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_A_flow,      [x, y, z,    t], region%help_fields%id_var_A_flow,      long_name='3D ice flow factor')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_A_flow_mean, [x, y,       t], region%help_fields%id_var_A_flow_mean, long_name='Vertically averaged ice flow factor')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_Fr_Aa,       [x, y,       t], region%help_fields%id_var_Fr_Aa, long_name='Geothermal heat flux')
    
    ! Dummy variables (useful for debugging)
    IF (C%do_write_dummy_output) THEN
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy2D_01,   [x, y,       t], region%help_fields%id_var_dummy2D_01, long_name='2D dummy variable 01')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy2D_02,   [x, y,       t], region%help_fields%id_var_dummy2D_02, long_name='2D dummy variable 02')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy2D_03,   [x, y,       t], region%help_fields%id_var_dummy2D_03, long_name='2D dummy variable 03')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy2D_04,   [x, y,       t], region%help_fields%id_var_dummy2D_04, long_name='2D dummy variable 04')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy2D_05,   [x, y,       t], region%help_fields%id_var_dummy2D_05, long_name='2D dummy variable 05')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy2D_06,   [x, y,       t], region%help_fields%id_var_dummy2D_06, long_name='2D dummy variable 06')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy2D_07,   [x, y,       t], region%help_fields%id_var_dummy2D_07, long_name='2D dummy variable 07')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy2D_08,   [x, y,       t], region%help_fields%id_var_dummy2D_08, long_name='2D dummy variable 08')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy2D_09,   [x, y,       t], region%help_fields%id_var_dummy2D_09, long_name='2D dummy variable 09')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy2D_10,   [x, y,       t], region%help_fields%id_var_dummy2D_10, long_name='2D dummy variable 10')
    
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy3D_01,   [x, y, z,    t], region%help_fields%id_var_dummy3D_01, long_name='3D dummy variable 01')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy3D_02,   [x, y, z,    t], region%help_fields%id_var_dummy3D_02, long_name='3D dummy variable 02')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy3D_03,   [x, y, z,    t], region%help_fields%id_var_dummy3D_03, long_name='3D dummy variable 03')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy3D_04,   [x, y, z,    t], region%help_fields%id_var_dummy3D_04, long_name='3D dummy variable 04')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy3D_05,   [x, y, z,    t], region%help_fields%id_var_dummy3D_05, long_name='3D dummy variable 05')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy3D_06,   [x, y, z,    t], region%help_fields%id_var_dummy3D_06, long_name='3D dummy variable 06')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy3D_07,   [x, y, z,    t], region%help_fields%id_var_dummy3D_07, long_name='3D dummy variable 07')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy3D_08,   [x, y, z,    t], region%help_fields%id_var_dummy3D_08, long_name='3D dummy variable 08')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy3D_09,   [x, y, z,    t], region%help_fields%id_var_dummy3D_09, long_name='3D dummy variable 09')
    CALL create_double_var( region%help_fields%ncid, region%help_fields%name_var_dummy3D_10,   [x, y, z,    t], region%help_fields%id_var_dummy3D_10, long_name='3D dummy variable 10')
    END IF ! IF (C%do_write_dummy_output) THEN
    
    ! Leave definition mode:
    CALL handle_error(nf90_enddef( region%help_fields%ncid))
    
    ! Write the x, y, zeta, months, and lat/lon variable data
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_x,        region%grid%x                            ))
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_y,        region%grid%y                            ))
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_zeta,     C%zeta                                   ))
    CALL handle_error( nf90_put_var( region%help_fields%ncid, region%help_fields%id_var_month,    (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)))
    CALL write_data_to_file_2D( region%help_fields%ncid, region%grid%nx, region%grid%ny, region%help_fields%id_var_lat, region%grid%lat)
    CALL write_data_to_file_2D( region%help_fields%ncid, region%grid%nx, region%grid%ny, region%help_fields%id_var_lon, region%grid%lon)
        
    ! Synchronize with disk (otherwise it doesn't seem to work on a MAC)
    CALL handle_error(nf90_sync( region%help_fields%ncid))
    
    ! Close the file
    CALL close_netcdf_file(region%help_fields%ncid)
    
  END SUBROUTINE create_help_fields_file
  
  SUBROUTINE write_data_to_file_2D( ncid, nx, ny,     var_id, d, start_vec)
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
    
  END SUBROUTINE write_data_to_file_2D
  SUBROUTINE write_data_to_file_3D( ncid, nx, ny, nz, var_id, d, start_vec)
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
    
  END SUBROUTINE write_data_to_file_3D
  
  SUBROUTINE read_inverse_routine_history_dT_glob(         forcing, filename)
    ! Read the inverse routine history from the specified NetCDF file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    CHARACTER(LEN=256),             INTENT(IN)    :: filename
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: ncid, id_dim_nH, id_dim_time, nt, id_var_H, id_var_time
    INTEGER                                       :: int_dummy
    INTEGER                                       :: ti, ti_min
    REAL(dp)                                      :: dt, dt_min
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: time
    
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
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: ncid, id_dim_nH, id_dim_time, nt, id_var_H, id_var_time
    INTEGER                                       :: int_dummy
    INTEGER                                       :: ti, ti_min
    REAL(dp)                                      :: dt, dt_min
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: time
    
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
    INTEGER                                       :: cerr, ierr
    INTEGER                                       :: ncid, id_dim_nH, id_dim_time, nt, id_var_H, id_var_time
    INTEGER                                       :: int_dummy
    INTEGER                                       :: ti, ti_min
    REAL(dp)                                      :: dt, dt_min
    REAL(dp), DIMENSION(:), ALLOCATABLE           :: time
    
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
  
  SUBROUTINE read_PD_data_file( PD)
    ! Read the PD netcdf file
   
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_PD_data_fields), INTENT(INOUT) :: PD
    
    ! Open the netcdf file
    CALL open_netcdf_file(PD%netcdf%filename, PD%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_x,      PD%x,      start=(/1   /)))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_y,      PD%y,      start=(/1   /)))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_Hi,     PD%Hi_raw, start=(/1, 1/)))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_Hb,     PD%Hb_raw, start=(/1, 1/)))
    CALL handle_error(nf90_get_var( PD%netcdf%ncid, PD%netcdf%id_var_Hs,     PD%Hs_raw, start=(/1, 1/)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD%netcdf%ncid)
    
  END SUBROUTINE read_PD_data_file 
  SUBROUTINE read_init_data_file( init)
    ! Read the init netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields), INTENT(INOUT) :: init
    
    ! Open the netcdf file
    CALL open_netcdf_file(init%netcdf%filename, init%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_x,      init%x,      start=(/1   /)))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_y,      init%y,      start=(/1   /)))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hi,     init%Hi_raw, start=(/1, 1/)))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hb,     init%Hb_raw, start=(/1, 1/)))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hs,     init%Hs_raw, start=(/1, 1/)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(init%netcdf%ncid)
    
  END SUBROUTINE read_init_data_file
  SUBROUTINE read_restart_file( init)
    ! Read the init netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields), INTENT(INOUT) :: init

    ! Local variables:
    INTEGER                                       :: cerr,ierr
    INTEGER                                       :: k,ti,ti_min
    REAL(dp)                                      :: dt, dt_min
    
    ! Open the netcdf file
    CALL open_netcdf_file(init%netcdf%filename, init%netcdf%ncid)
    
    ! Read zeta, check if it matches the config zeta
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_zeta, init%zeta, start=(/1/) ))
    IF (init%nz /= C%nz) THEN
      WRITE(0,*) '  WARNING - vertical coordinate zeta in restart file doesnt match zeta in config!'
    ELSE
      DO k = 1, C%nz
        IF (ABS(C%zeta(k) - init%zeta(k)) > 0.0001_dp) THEN
          WRITE(0,*) '  WARNING - vertical coordinate zeta in restart file doesnt match zeta in config!'
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
      WRITE(0,*) '  WARNING - no exact match for time_to_restart_from ', C%time_to_restart_from, ' in restart file! Reading closest match ', init%time( ti), ' instead.'
    END IF
    IF (C%time_to_restart_from /= C%start_time_of_run) THEN
      WRITE(0,*) '  WARNING - starting run at t = ', C%start_time_of_run, ' with restart data at t = ', C%time_to_restart_from
    END IF
    
    ! Read the data
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hi,               init%Hi_raw,               start=(/1, 1,    ti/), count=(/init%nx, init%ny,          1/) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hb,               init%Hb_raw,               start=(/1, 1,    ti/), count=(/init%nx, init%ny,          1/) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Hs,               init%Hs_raw,               start=(/1, 1,    ti/), count=(/init%nx, init%ny,          1/) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_Ti,               init%Ti_raw,               start=(/1, 1, 1, ti/), count=(/init%nx, init%ny, init%nz, 1/) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_U_SSA,            init%U_SSA_raw,            start=(/1, 1,    ti/), count=(/init%nx, init%ny,          1/) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_V_SSA,            init%V_SSA_raw,            start=(/1, 1,    ti/), count=(/init%nx, init%ny,          1/) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_FirnDepth,        init%FirnDepth_raw,        start=(/1, 1, 1, ti/), count=(/init%nx, init%ny, 12,      1/) ))
    CALL handle_error(nf90_get_var( init%netcdf%ncid, init%netcdf%id_var_MeltPreviousYear, init%MeltPreviousYear_raw, start=(/1, 1,    ti/), count=(/init%nx, init%ny,          1/) ))
        
    ! Close the netcdf file
    CALL close_netcdf_file(init%netcdf%ncid)
    
  END SUBROUTINE read_restart_file
  SUBROUTINE read_PD_obs_data_file( PD_obs)
    ! Read the PD_obs0 netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: PD_obs
    
    ! Open the netcdf file
    CALL open_netcdf_file(PD_obs%netcdf%filename, PD_obs%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_lon,     PD_obs%lon,     start=(/1      /)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_lat,     PD_obs%lat,     start=(/1      /)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_T2m,     PD_obs%T2m,     start=(/1, 1, 1/)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Precip,  PD_obs%Precip,  start=(/1, 1, 1/)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Hs,      PD_obs%Hs,      start=(/1, 1   /)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Wind_WE, PD_obs%Wind_WE, start=(/1, 1, 1/)))
    CALL handle_error(nf90_get_var( PD_obs%netcdf%ncid, PD_obs%netcdf%id_var_Wind_SN, PD_obs%Wind_SN, start=(/1, 1, 1/)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(PD_obs%netcdf%ncid)
    
  END SUBROUTINE read_PD_obs_data_file
  SUBROUTINE read_GCM_snapshot( snapshot)
    ! Read the PD_obs0 netcdf file
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: snapshot
    
    ! Open the netcdf file
    CALL open_netcdf_file(snapshot%netcdf%filename, snapshot%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_lon,     snapshot%lon,     start=(/1      /)))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_lat,     snapshot%lat,     start=(/1      /)))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Hi,      snapshot%Hi,      start=(/1, 1   /)))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Hs,      snapshot%Hs,      start=(/1, 1   /)))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_T2m,     snapshot%T2m,     start=(/1, 1, 1/)))
    CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Precip,  snapshot%Precip,  start=(/1, 1, 1/)))
    !CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Wind_WE, snapshot%Wind_WE, start=(/1, 1, 1/)))
    !CALL handle_error(nf90_get_var( snapshot%netcdf%ncid, snapshot%netcdf%id_var_Wind_SN, snapshot%Wind_SN, start=(/1, 1, 1/)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(snapshot%netcdf%ncid)
    
  END SUBROUTINE read_GCM_snapshot
  SUBROUTINE read_insolation_data_file( forcing, ti0, ti1, ins_Q_TOA0, ins_Q_TOA1) 
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    INTEGER,                        INTENT(IN)    :: ti0, ti1
    REAL(dp), DIMENSION(:,:),       INTENT(OUT)   :: ins_Q_TOA0, ins_Q_TOA1
    
    ! Local variables:
    INTEGER                                       :: mi, li
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: Q_temp0, Q_temp1
    
    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE( Q_temp0(1, 12, forcing%ins_nlat))
    ALLOCATE( Q_temp1(1, 12, forcing%ins_nlat))
        
    ! Read data
    CALL open_netcdf_file(forcing%netcdf%filename, forcing%netcdf%ncid)
    CALL handle_error(nf90_get_var( forcing%netcdf%ncid, forcing%netcdf%id_var_Q_TOA, Q_temp0, start=(/ti0, 1, 1/), count=(/ 1, 12, forcing%ins_nlat/), stride=(/1, 1, 1/) ))    
    CALL handle_error(nf90_get_var( forcing%netcdf%ncid, forcing%netcdf%id_var_Q_TOA, Q_temp1, start=(/ti1, 1, 1/), count=(/ 1, 12, forcing%ins_nlat/), stride=(/1, 1, 1/) ))
    CALL close_netcdf_file(forcing%netcdf%ncid) 
    
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
    
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf%filename, forcing%netcdf%ncid)
    
    ! Read the data
    CALL handle_error(nf90_get_var( forcing%netcdf%ncid, forcing%netcdf%id_var_time,    forcing%ins_time,    start=(/1      /)))
    CALL handle_error(nf90_get_var( forcing%netcdf%ncid, forcing%netcdf%id_var_lat,     forcing%ins_lat,     start=(/1      /)))
        
    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf%ncid)
    
  END SUBROUTINE read_insolation_data_file_time_lat

  SUBROUTINE inquire_PD_data_file( PD) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_PD_data_fields), INTENT(INOUT) :: PD
        
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
  SUBROUTINE inquire_init_data_file( init) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields), INTENT(INOUT) :: init
        
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
  SUBROUTINE inquire_restart_file( init) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_init_data_fields), INTENT(INOUT) :: init
 
    ! Local variables:
    INTEGER                               :: int_dummy
        
    ! Open the netcdf file
    CALL open_netcdf_file(init%netcdf%filename, init%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist, return their lengths.
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_x,    init%nx, init%netcdf%id_dim_x)
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_y,    init%ny, init%netcdf%id_dim_y)
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_zeta, init%nz, init%netcdf%id_dim_zeta)
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_time, init%nt, init%netcdf%id_dim_time)
    CALL inquire_dim( init%netcdf%ncid, init%netcdf%name_dim_month, int_dummy, init%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_x,                (/ init%netcdf%id_dim_x                                                                          /), init%netcdf%id_var_x   )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_y,                (/                       init%netcdf%id_dim_y                                                    /), init%netcdf%id_var_y   )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_zeta,             (/                                             init%netcdf%id_dim_zeta                           /), init%netcdf%id_var_zeta)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_time,             (/                                                                       init%netcdf%id_dim_time /), init%netcdf%id_var_time)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Hi,               (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y                          , init%netcdf%id_dim_time /), init%netcdf%id_var_Hi  )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Hb,               (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y                          , init%netcdf%id_dim_time /), init%netcdf%id_var_Hb  )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Hs,               (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y                          , init%netcdf%id_dim_time /), init%netcdf%id_var_Hs  )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_Ti,               (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y, init%netcdf%id_dim_zeta , init%netcdf%id_dim_time /), init%netcdf%id_var_Ti  )
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_U_SSA,            (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y                          , init%netcdf%id_dim_time /), init%netcdf%id_var_U_SSA)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_V_SSA,            (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y                          , init%netcdf%id_dim_time /), init%netcdf%id_var_V_SSA)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_FirnDepth,        (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y, init%netcdf%id_dim_month, init%netcdf%id_dim_time /), init%netcdf%id_var_FirnDepth)
    CALL inquire_double_var( init%netcdf%ncid, init%netcdf%name_var_MeltPreviousYear, (/ init%netcdf%id_dim_x, init%netcdf%id_dim_y                          , init%netcdf%id_dim_time /), init%netcdf%id_var_MeltPreviousYear)
        
    ! Close the netcdf file
    CALL close_netcdf_file(init%netcdf%ncid)
    
  END SUBROUTINE inquire_restart_file
  SUBROUTINE inquire_PD_obs_data_file( PD_obs) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: PD_obs
 
    ! Local variables:
    INTEGER                               :: int_dummy
        
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
  SUBROUTINE inquire_GCM_snapshot( snapshot) 
    ! Check if the right dimensions and variables are present in the file.
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_subclimate_global), INTENT(INOUT) :: snapshot
 
    ! Local variables:
    INTEGER                               :: int_dummy
        
    ! Open the netcdf file
    CALL open_netcdf_file( snapshot%netcdf%filename, snapshot%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( snapshot%netcdf%ncid, snapshot%netcdf%name_dim_lat,     snapshot%nlat,  snapshot%netcdf%id_dim_lat)
    CALL inquire_dim( snapshot%netcdf%ncid, snapshot%netcdf%name_dim_lon,     snapshot%nlon,  snapshot%netcdf%id_dim_lon)
    CALL inquire_dim( snapshot%netcdf%ncid, snapshot%netcdf%name_dim_month,   int_dummy,      snapshot%netcdf%id_dim_month)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_lat,      (/ snapshot%netcdf%id_dim_lat                                                           /),  snapshot%netcdf%id_var_lat)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_lon,      (/ snapshot%netcdf%id_dim_lon                                                           /),  snapshot%netcdf%id_var_lon)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Hi,       (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat                               /),  snapshot%netcdf%id_var_Hi)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Hs,       (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat                               /),  snapshot%netcdf%id_var_Hs)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_T2m,      (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat, snapshot%netcdf%id_dim_month /),  snapshot%netcdf%id_var_T2m)
    CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Precip,   (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat, snapshot%netcdf%id_dim_month /),  snapshot%netcdf%id_var_Precip)
    !CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Wind_WE,  (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat, snapshot%netcdf%id_dim_month /),  snapshot%netcdf%id_var_Wind_WE)
    !CALL inquire_double_var( snapshot%netcdf%ncid, snapshot%netcdf%name_var_Wind_SN,  (/ snapshot%netcdf%id_dim_lon, snapshot%netcdf%id_dim_lat, snapshot%netcdf%id_dim_month /),  snapshot%netcdf%id_var_Wind_SN)
        
    ! Close the netcdf file
    CALL close_netcdf_file(snapshot%netcdf%ncid)
    
  END SUBROUTINE inquire_GCM_snapshot 
  SUBROUTINE inquire_insolation_data_file( forcing)
    IMPLICIT NONE
    
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
 
    ! Local variables:  
    INTEGER                                :: int_dummy
            
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf%filename, forcing%netcdf%ncid)
    
    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( forcing%netcdf%ncid, forcing%netcdf%name_dim_time,     forcing%ins_nyears,        forcing%netcdf%id_dim_time)
    CALL inquire_dim( forcing%netcdf%ncid, forcing%netcdf%name_dim_month,    int_dummy,                 forcing%netcdf%id_dim_month)  
    CALL inquire_dim( forcing%netcdf%ncid, forcing%netcdf%name_dim_lat,      forcing%ins_nlat,          forcing%netcdf%id_dim_lat)
    
    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( forcing%netcdf%ncid, forcing%netcdf%name_var_time,  (/ forcing%netcdf%id_dim_time                                                         /), forcing%netcdf%id_var_time)
    CALL inquire_double_var( forcing%netcdf%ncid, forcing%netcdf%name_var_month, (/ forcing%netcdf%id_dim_month                                                        /), forcing%netcdf%id_var_month)
    CALL inquire_double_var( forcing%netcdf%ncid, forcing%netcdf%name_var_lat,   (/ forcing%netcdf%id_dim_lat                                                          /), forcing%netcdf%id_var_lat)
    CALL inquire_double_var( forcing%netcdf%ncid, forcing%netcdf%name_var_Q_TOA, (/ forcing%netcdf%id_dim_time, forcing%netcdf%id_dim_month, forcing%netcdf%id_dim_lat /), forcing%netcdf%id_var_Q_TOA)
        
    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf%ncid)
    
  END SUBROUTINE inquire_insolation_data_file

  SUBROUTINE inquire_geothermal_heat_flux_file( forcing)
    IMPLICIT NONE
  
    ! Output variable
    TYPE(type_forcing_data), INTENT(INOUT) :: forcing
          
    ! Open the netcdf file
    CALL open_netcdf_file(forcing%netcdf_ghf%filename, forcing%netcdf_ghf%ncid)

    ! Inquire dimensions id's. Check that all required dimensions exist return their lengths.
    CALL inquire_dim( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_dim_lon,      forcing%ghf_nlon,                 forcing%netcdf_ghf%id_dim_lon)
    CALL inquire_dim( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_dim_lat,      forcing%ghf_nlat,                 forcing%netcdf_ghf%id_dim_lat)

    ! Inquire variable id's. Make sure that each variable has the correct dimensions:
    CALL inquire_double_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_var_lon,   (/ forcing%netcdf_ghf%id_dim_lon                                                          /), forcing%netcdf_ghf%id_var_lon)
    CALL inquire_double_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_var_lat,   (/ forcing%netcdf_ghf%id_dim_lat                                                          /), forcing%netcdf_ghf%id_var_lat)
    CALL inquire_double_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%name_var_ghf, (/ forcing%netcdf_ghf%id_dim_lon, forcing%netcdf_ghf%id_dim_lat /), forcing%netcdf_ghf%id_var_ghf)
      
    ! Close the netcdf file
    CALL close_netcdf_file(forcing%netcdf_ghf%ncid)
  
  END SUBROUTINE inquire_geothermal_heat_flux_file

  SUBROUTINE read_geothermal_heat_flux_file( forcing, ghf_ghf)
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_forcing_data),        INTENT(INOUT) :: forcing
    REAL(dp), DIMENSION(:,:),       INTENT(OUT)   :: ghf_ghf
  
    ! Local variables:
    INTEGER                                       :: mi, li
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: GHF
  
    ! Temporary memory to store the data read from the netCDF file
    ALLOCATE( GHF(forcing%ghf_nlon, forcing%ghf_nlat))
      
    ! Read data
    CALL open_netcdf_file(forcing%netcdf_ghf%filename, forcing%netcdf_ghf%ncid)
    CALL handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_lon,     forcing%ghf_lon,    start=(/1      /)))
    CALL handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_lat,     forcing%ghf_lat,    start=(/1      /)))
    CALL handle_error(nf90_get_var( forcing%netcdf_ghf%ncid, forcing%netcdf_ghf%id_var_ghf, GHF, start=(/1, 1/), count=(/forcing%ghf_nlon, forcing%ghf_nlat/), stride=(/1, 1/) ))
    CALL close_netcdf_file(forcing%netcdf_ghf%ncid)
  
    ! Store the data in the shared memory structure
    DO mi = 1, forcing%ghf_nlat
    DO li = 1, forcing%ghf_nlon
      ghf_ghf( li,mi) = GHF( li,mi) * sec_per_year ! Convert from W m-2 (J m-2 s-1) to J m-2 yr-1
    END DO
    END DO
      
    ! Clean up temporary memory
    DEALLOCATE(GHF)
 
  END SUBROUTINE read_geothermal_heat_flux_file

  SUBROUTINE get_output_filenames( region)
   
    IMPLICIT NONE
    
    ! Input variables:
    TYPE(type_model_region), INTENT(INOUT) :: region

    CHARACTER(LEN=20)          :: short_filename
    INTEGER                    :: n

    short_filename = 'restart_NAM.nc'
    short_filename(9:11) = region%name
    DO n = 1, 256
      region%restart%filename(n:n) = ' '
    END DO
    region%restart%filename = TRIM(C%output_dir)//TRIM(short_filename)

    short_filename = 'help_fields_NAM.nc'
    short_filename(13:15) = region%name
    DO n = 1, 256
      region%help_fields%filename(n:n) = ' '
    END DO
    region%help_fields%filename = TRIM(C%output_dir)//TRIM(short_filename)

  END SUBROUTINE get_output_filenames

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
