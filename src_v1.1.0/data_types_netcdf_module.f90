MODULE data_types_netcdf_module
  ! Contains the TYPES for different NetCDf files read and written by UFEMISM.

  USE configuration_module,        ONLY: dp, C

  IMPLICIT NONE
    
  TYPE type_netcdf_output
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! Index of time frame to be written to
    INTEGER :: ti
    
  ! ID's for dimensions:
  ! ===================
    
    ! Dimensions
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_zeta
    INTEGER :: id_dim_time
    INTEGER :: id_dim_month
  
    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    
  ! Restart file data - everything that's needed to restart a new run
  ! =================================================================
  
    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_time
    INTEGER :: id_var_zeta
    INTEGER :: id_var_month
    
    INTEGER :: id_var_Hi
    INTEGER :: id_var_Hb
    INTEGER :: id_var_Hs
    INTEGER :: id_var_U_SIA
    INTEGER :: id_var_V_SIA
    INTEGER :: id_var_U_SSA
    INTEGER :: id_var_V_SSA
    INTEGER :: id_var_Ti
    INTEGER :: id_var_FirnDepth
    INTEGER :: id_var_MeltPreviousYear
    
    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '
  
    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_var_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   '
    CHARACTER(LEN=256) :: name_var_U_SIA                 = 'U_SIA                '
    CHARACTER(LEN=256) :: name_var_V_SIA                 = 'V_SIA                '
    CHARACTER(LEN=256) :: name_var_U_SSA                 = 'U_SSA                '
    CHARACTER(LEN=256) :: name_var_V_SSA                 = 'V_SSA                '
    CHARACTER(LEN=256) :: name_var_Ti                    = 'Ti                   '
    CHARACTER(LEN=256) :: name_var_FirnDepth             = 'FirnDepth            '
    CHARACTER(LEN=256) :: name_var_MeltPreviousYear      = 'MeltPreviousYear     '
    
  ! Restart data for the inverse routine
  ! ====================================
  
    INTEGER :: id_dim_ndT_glob_history
    INTEGER :: id_dim_ndT_glob_inverse_history
    INTEGER :: id_dim_nCO2_inverse_history
    
    CHARACTER(LEN=256) :: name_dim_ndT_glob_history         = 'ndT_glob_history        '
    CHARACTER(LEN=256) :: name_dim_ndT_glob_inverse_history = 'ndT_glob_inverse_history'
    CHARACTER(LEN=256) :: name_dim_nCO2_inverse_history     = 'nCO2_inverse_history    '
  
    INTEGER :: id_var_dT_glob_history
    INTEGER :: id_var_dT_glob_inverse_history
    INTEGER :: id_var_CO2_inverse_history
    
    CHARACTER(LEN=256) :: name_var_dT_glob_history          = 'dT_glob_history        '
    CHARACTER(LEN=256) :: name_var_dT_glob_inverse_history  = 'dT_glob_inverse_history'
    CHARACTER(LEN=256) :: name_var_CO2_inverse_history      = 'CO2_inverse_history    '
    
  ! Help fields file data - everything that's not in the restart file, but is still interesting
  ! ===========================================================================================
    
    INTEGER :: id_var_lat
    INTEGER :: id_var_lon
    INTEGER :: id_var_U_3D
    INTEGER :: id_var_V_3D
    INTEGER :: id_var_W_3D
    INTEGER :: id_var_mask
    INTEGER :: id_var_T2m
    INTEGER :: id_var_Precip
    INTEGER :: id_var_Albedo
    INTEGER :: id_var_SMB
    INTEGER :: id_var_BMB
    INTEGER :: id_var_dHs_dx
    INTEGER :: id_var_dHs_dy 
    INTEGER :: id_var_D_SIA
    INTEGER :: id_var_Wind_WE
    INTEGER :: id_var_Wind_SN
    INTEGER :: id_var_Snowfall
    INTEGER :: id_var_Melt
    INTEGER :: id_var_Refreezing
    INTEGER :: id_var_Runoff
    INTEGER :: id_var_A_flow
    INTEGER :: id_var_A_flow_mean
    
    INTEGER :: id_var_dummy2D_01
    INTEGER :: id_var_dummy2D_02
    INTEGER :: id_var_dummy2D_03
    INTEGER :: id_var_dummy2D_04
    INTEGER :: id_var_dummy2D_05
    INTEGER :: id_var_dummy2D_06
    INTEGER :: id_var_dummy2D_07
    INTEGER :: id_var_dummy2D_08
    INTEGER :: id_var_dummy2D_09
    INTEGER :: id_var_dummy2D_10
    
    INTEGER :: id_var_dummy3D_01
    INTEGER :: id_var_dummy3D_02
    INTEGER :: id_var_dummy3D_03
    INTEGER :: id_var_dummy3D_04
    INTEGER :: id_var_dummy3D_05
    INTEGER :: id_var_dummy3D_06
    INTEGER :: id_var_dummy3D_07
    INTEGER :: id_var_dummy3D_08
    INTEGER :: id_var_dummy3D_09
    INTEGER :: id_var_dummy3D_10
  
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_U_3D                  = 'U_3D                 '
    CHARACTER(LEN=256) :: name_var_V_3D                  = 'V_3D                 '
    CHARACTER(LEN=256) :: name_var_W_3D                  = 'W_3D                 ' 
    CHARACTER(LEN=256) :: name_var_mask                  = 'mask                 '
    CHARACTER(LEN=256) :: name_var_T2m                   = 'T2m                  '
    CHARACTER(LEN=256) :: name_var_Precip                = 'Precip               '
    CHARACTER(LEN=256) :: name_var_Albedo                = 'Albedo               '
    CHARACTER(LEN=256) :: name_var_SMB                   = 'SMB                  '
    CHARACTER(LEN=256) :: name_var_BMB                   = 'BMB                  '
    CHARACTER(LEN=256) :: name_var_dHs_dx                = 'dHs_dx               '
    CHARACTER(LEN=256) :: name_var_dHs_dy                = 'dHs_dy               '
    CHARACTER(LEN=256) :: name_var_D_SIA                 = 'D_SIA                '
    CHARACTER(LEN=256) :: name_var_Wind_WE               = 'Wind_WE              '
    CHARACTER(LEN=256) :: name_var_Wind_SN               = 'Wind_SN              '
    CHARACTER(LEN=256) :: name_var_Snowfall              = 'Snowfall             '
    CHARACTER(LEN=256) :: name_var_Melt                  = 'Melt                 '
    CHARACTER(LEN=256) :: name_var_Refreezing            = 'Refreezing           '
    CHARACTER(LEN=256) :: name_var_Runoff                = 'Runoff               '
    CHARACTER(LEN=256) :: name_var_A_flow                = 'A_flow               '
    CHARACTER(LEN=256) :: name_var_A_flow_mean           = 'A_flow_mean          '
    
    CHARACTER(LEN=256) :: name_var_dummy2D_01            = 'dummy2D_01           '
    CHARACTER(LEN=256) :: name_var_dummy2D_02            = 'dummy2D_02           '
    CHARACTER(LEN=256) :: name_var_dummy2D_03            = 'dummy2D_03           '
    CHARACTER(LEN=256) :: name_var_dummy2D_04            = 'dummy2D_04           '
    CHARACTER(LEN=256) :: name_var_dummy2D_05            = 'dummy2D_05           '
    CHARACTER(LEN=256) :: name_var_dummy2D_06            = 'dummy2D_06           '
    CHARACTER(LEN=256) :: name_var_dummy2D_07            = 'dummy2D_07           '
    CHARACTER(LEN=256) :: name_var_dummy2D_08            = 'dummy2D_08           '
    CHARACTER(LEN=256) :: name_var_dummy2D_09            = 'dummy2D_09           '
    CHARACTER(LEN=256) :: name_var_dummy2D_10            = 'dummy2D_10           '
    
    CHARACTER(LEN=256) :: name_var_dummy3D_01            = 'dummy3D_01           '
    CHARACTER(LEN=256) :: name_var_dummy3D_02            = 'dummy3D_02           '
    CHARACTER(LEN=256) :: name_var_dummy3D_03            = 'dummy3D_03           '
    CHARACTER(LEN=256) :: name_var_dummy3D_04            = 'dummy3D_04           '
    CHARACTER(LEN=256) :: name_var_dummy3D_05            = 'dummy3D_05           '
    CHARACTER(LEN=256) :: name_var_dummy3D_06            = 'dummy3D_06           '
    CHARACTER(LEN=256) :: name_var_dummy3D_07            = 'dummy3D_07           '
    CHARACTER(LEN=256) :: name_var_dummy3D_08            = 'dummy3D_08           '
    CHARACTER(LEN=256) :: name_var_dummy3D_09            = 'dummy3D_09           '
    CHARACTER(LEN=256) :: name_var_dummy3D_10            = 'dummy3D_10           '
        
  END TYPE type_netcdf_output
    
  TYPE type_netcdf_help_fields
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! Index of time frame to be written to
    INTEGER :: ti
    
  ! ID's for dimensions:
  ! ===================
    
    ! Dimensions
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_zeta
    INTEGER :: id_dim_time
    INTEGER :: id_dim_month
  
    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    
  END TYPE type_netcdf_help_fields
    
  TYPE type_netcdf_PD_data
    ! For reading an input file describing a present-day model region, on a Cartesian grid
  
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! ID's for variables:
    ! ===================
    
    ! Dimensions
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_month
    
    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    
    ! Variables:
    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_Hi
    INTEGER :: id_var_Hb
    INTEGER :: id_var_Hs
    
    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_var_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   ' 
        
  END TYPE type_netcdf_PD_data
    
  TYPE type_netcdf_init_data
    ! For reading an input file describing the initial state of a model region, on a Cartesian grid
    
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! ID's for variables:
    ! ===================
    
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_zeta
    INTEGER :: id_dim_time
    INTEGER :: id_dim_month
    
    ! Variable names
    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    
    ! Key output
    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_zeta
    INTEGER :: id_var_time
    INTEGER :: id_var_month
    INTEGER :: id_var_Hi
    INTEGER :: id_var_Hb
    INTEGER :: id_var_Hs
    INTEGER :: id_var_Ti
    INTEGER :: id_var_U_SSA
    INTEGER :: id_var_V_SSA
    INTEGER :: id_var_FirnDepth
    INTEGER :: id_var_MeltPreviousYear
    
    ! Variable names
    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '
    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_var_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   '
    CHARACTER(LEN=256) :: name_var_Ti                    = 'Ti                   '
    CHARACTER(LEN=256) :: name_var_U_SSA                 = 'U_SSA                '
    CHARACTER(LEN=256) :: name_var_V_SSA                 = 'V_SSA                '
    CHARACTER(LEN=256) :: name_var_FirnDepth             = 'FirnDepth            '
    CHARACTER(LEN=256) :: name_var_MeltPreviousYear      = 'MeltPreviousYear     '
        
  END TYPE type_netcdf_init_data
    
  TYPE type_netcdf_climate_data
    ! For reading an input file containing either a GCM snapshot or a PD observations data set (e.g. ERA-40),
    ! describing the global climate with monthly fields on a lat/lon grid
  
    ! Integers describing open ports to different variables in an opened NetCDF file.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! ID's for variables:
    ! ===================
    
    ! Dimensions
    INTEGER :: id_dim_lon
    INTEGER :: id_dim_lat
    INTEGER :: id_dim_month
    
    CHARACTER(LEN=256) :: name_dim_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    
    ! Variables
    INTEGER :: id_var_lon
    INTEGER :: id_var_lat
    INTEGER :: id_var_Hi
    INTEGER :: id_var_Hs
    INTEGER :: id_var_T2m
    INTEGER :: id_var_Precip
    INTEGER :: id_var_Wind_WE
    INTEGER :: id_var_Wind_SN
    
    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   '
    CHARACTER(LEN=256) :: name_var_T2m                   = 'T2m                  '
    CHARACTER(LEN=256) :: name_var_Precip                = 'Precip               '
    CHARACTER(LEN=256) :: name_var_Wind_WE               = 'Wind_WE              '
    CHARACTER(LEN=256) :: name_var_Wind_SN               = 'Wind_SN              '
        
  END TYPE type_netcdf_climate_data
    
  TYPE type_netcdf_insolation
    ! For reading an input file containing an insolation history reconstruction (e.g. Lasker et al., 2004),
    ! describing top-of-the-atmosphere insolation for every month of the year at a latitudinal grid.
  
    ! Integers describing open ports to different variables in an opened NetCDF file.
    
    CHARACTER(LEN=256) :: filename
    
    ! ID for NetCDF file:
    INTEGER :: ncid
    
    ! ID's for variables:
    ! ===================
    
    ! Dimensions
    INTEGER :: id_dim_time
    INTEGER :: id_dim_month
    INTEGER :: id_dim_lat
    
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'lat                  '
    
    ! Variables
    INTEGER :: id_var_time
    INTEGER :: id_var_month
    INTEGER :: id_var_lat
    INTEGER :: id_var_Q_TOA
    
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_Q_TOA                 = 'Q_TOA                '
        
  END TYPE type_netcdf_insolation
  
CONTAINS

END MODULE data_types_netcdf_module
