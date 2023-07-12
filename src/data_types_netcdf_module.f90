MODULE data_types_netcdf_module

  ! Contains the TYPEs for different NetCDf files read and written by IMAU-ICE.

  USE configuration_module,        ONLY: dp, C

  IMPLICIT NONE

  TYPE type_netcdf_resource_tracker
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! Index of time frame to be written to
    INTEGER :: ti

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_time
    INTEGER :: id_dim_name_length

    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_name_length           = 'name_length          '

    INTEGER :: id_var_time

    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '

  ! Variables
  ! =========

    INTEGER, DIMENSION(:), ALLOCATABLE :: id_vars
    INTEGER, DIMENSION(:), ALLOCATABLE :: id_var_names

  END TYPE type_netcdf_resource_tracker

  TYPE type_netcdf_reference_geometry
    ! For reading an input file describing a reference ice-sheet geometry on a Cartesian grid

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! ID's for variables:
    ! ===================

    ! Dimensions
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '

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

  END TYPE type_netcdf_reference_geometry

  TYPE type_netcdf_climate_snapshot
    ! For reading an input file containing monthly climate data on a lon/lat-grid

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_lon
    INTEGER :: id_dim_lat
    INTEGER :: id_dim_month

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '

  ! Variables
  ! =========

    ! Dimensions
    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_lon
    INTEGER :: id_var_lat
    INTEGER :: id_var_month

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '

    ! Field variables
    INTEGER :: id_var_Hs
    INTEGER :: id_var_T2m
    INTEGER :: id_var_Precip
    INTEGER :: id_var_Wind_WE
    INTEGER :: id_var_Wind_SN
    INTEGER :: id_var_Wind_LR
    INTEGER :: id_var_Wind_DU

    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   '
    CHARACTER(LEN=256) :: name_var_T2m                   = 'T2m                  '
    CHARACTER(LEN=256) :: name_var_Precip                = 'Precip               '
    CHARACTER(LEN=256) :: name_var_Wind_WE               = 'Wind_WE              '
    CHARACTER(LEN=256) :: name_var_Wind_SN               = 'Wind_SN              '
    CHARACTER(LEN=256) :: name_var_Wind_LR               = 'Wind_LR              '
    CHARACTER(LEN=256) :: name_var_Wind_DU               = 'Wind_DU              '

  END TYPE type_netcdf_climate_snapshot

  TYPE type_netcdf_global_ocean_data
    ! For reading an input file containing ocean data on a lon/lat-grid

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_lon
    INTEGER :: id_dim_lat
    INTEGER :: id_dim_z_ocean

    CHARACTER(LEN=256) :: name_dim_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_dim_z_ocean               = 'depth                '

  ! Variables
  ! =========

    ! Dimensions
    INTEGER :: id_var_lon
    INTEGER :: id_var_lat
    INTEGER :: id_var_z_ocean

    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_z_ocean               = 'depth                '

    ! Field variables
    INTEGER :: id_var_T_ocean
    INTEGER :: id_var_S_ocean

    ! Variables names for ocean temperature and salinity are set through the config!

  END TYPE type_netcdf_global_ocean_data


  TYPE type_netcdf_insolation
    ! For reading an input file containing an insolation history reconstruction (e.g. Lasker et al., 2004),
    ! describing top-of-the-atmosphere insolation for every month of the year at a latitudinal grid.

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

  TYPE type_netcdf_regional_SMB_data
    ! For reading an input file containing monthly SMB and 2-m air temperature on an x/y-grid

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_time
    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_x                     = 'NX                   '
    CHARACTER(LEN=256) :: name_dim_y                     = 'NY                   '

  ! Variables
  ! ==========

    ! Dimensions
    INTEGER :: id_var_time
    INTEGER :: id_var_x
    INTEGER :: id_var_y

    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '

    ! Field variables
    INTEGER :: id_var_T2m_year
    INTEGER :: id_var_SMB_year

    CHARACTER(LEN=256) :: name_var_T2m_year              = 'T2m                  '
    CHARACTER(LEN=256) :: name_var_SMB_year              = 'SMB                  '

  END TYPE type_netcdf_regional_SMB_data

  TYPE type_netcdf_SELEN_global_topo
    ! A NETCDF file containing global topography data for SELEN on an irregular global mesh

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_vi
    INTEGER :: id_dim_ti
    INTEGER :: id_dim_ci
    INTEGER :: id_dim_three

    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '

    INTEGER :: id_var_V
    INTEGER :: id_var_Tri
    INTEGER :: id_var_nC
    INTEGER :: id_var_C
    INTEGER :: id_var_niTri
    INTEGER :: id_var_iTri

    CHARACTER(LEN=256) :: name_var_V                     = 'V                    '
    CHARACTER(LEN=256) :: name_var_Tri                   = 'Tri                  '
    CHARACTER(LEN=256) :: name_var_nC                    = 'nC                   '
    CHARACTER(LEN=256) :: name_var_C                     = 'C                    '
    CHARACTER(LEN=256) :: name_var_niTri                 = 'niTri                '
    CHARACTER(LEN=256) :: name_var_iTri                  = 'iTri                 '

  ! Variables
  ! =========

    INTEGER :: id_var_lat
    INTEGER :: id_var_lon
    INTEGER :: id_var_Hb
    INTEGER :: id_var_ianc

    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_var_ianc                  = 'ianc                 '

  END TYPE type_netcdf_SELEN_global_topo

  TYPE type_netcdf_SELEN_output
    ! A NETCDF file containing output of SELEN (ice loading, bed topography and geoid perturbation)
    ! on the global SELEN mesh

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! Index of time frame to be written to
    INTEGER :: ti

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_vi
    INTEGER :: id_dim_ti
    INTEGER :: id_dim_ci
    INTEGER :: id_dim_three
    INTEGER :: id_dim_time

    CHARACTER(LEN=256) :: name_dim_vi                    = 'vi                   '
    CHARACTER(LEN=256) :: name_dim_ti                    = 'ti                   '
    CHARACTER(LEN=256) :: name_dim_ci                    = 'ci                   '
    CHARACTER(LEN=256) :: name_dim_three                 = 'three                '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '

    INTEGER :: id_var_V
    INTEGER :: id_var_Tri
    INTEGER :: id_var_nC
    INTEGER :: id_var_C
    INTEGER :: id_var_niTri
    INTEGER :: id_var_iTri
    INTEGER :: id_var_time

    CHARACTER(LEN=256) :: name_var_V                     = 'V                    '
    CHARACTER(LEN=256) :: name_var_Tri                   = 'Tri                  '
    CHARACTER(LEN=256) :: name_var_nC                    = 'nC                   '
    CHARACTER(LEN=256) :: name_var_C                     = 'C                    '
    CHARACTER(LEN=256) :: name_var_niTri                 = 'niTri                '
    CHARACTER(LEN=256) :: name_var_iTri                  = 'iTri                 '
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '

  ! Variables
  ! =========

    INTEGER :: id_var_lat
    INTEGER :: id_var_lon
    INTEGER :: id_var_Hi
    INTEGER :: id_var_Hi_rel
    INTEGER :: id_var_N
    INTEGER :: id_var_U
    INTEGER :: id_var_ocean_function

    CHARACTER(LEN=256) :: name_var_lat                   = 'lat                  '
    CHARACTER(LEN=256) :: name_var_lon                   = 'lon                  '
    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_var_Hi_rel                = 'Hi_rel               '
    CHARACTER(LEN=256) :: name_var_N                     = 'N                    '
    CHARACTER(LEN=256) :: name_var_U                     = 'U                    '
    CHARACTER(LEN=256) :: name_var_ocean_function        = 'ocean_function       '

  END TYPE type_netcdf_SELEN_output

  TYPE type_netcdf_extrapolated_ocean_data
    ! A NetCDF file containing extrapolated ocean data on a high-resolution grid

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! Index of time frame to be written to
    INTEGER :: ti

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_z_ocean

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_z_ocean               = 'z_ocean              '

    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_z_ocean

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_z_ocean               = 'z_ocean              '

  ! Variables
  ! =========

    INTEGER :: id_var_T_ocean
    INTEGER :: id_var_S_ocean

    CHARACTER(LEN=256) :: name_var_T_ocean               = 'T_ocean              '
    CHARACTER(LEN=256) :: name_var_S_ocean               = 'S_ocean              '

  END TYPE type_netcdf_extrapolated_ocean_data


CONTAINS

END MODULE data_types_netcdf_module
