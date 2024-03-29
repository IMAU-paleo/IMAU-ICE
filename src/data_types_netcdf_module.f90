MODULE data_types_netcdf_module

  ! Contains the TYPEs for different NetCDf files read and written by IMAU-ICE.

  USE configuration_module,        ONLY: dp, C

  IMPLICIT NONE

  TYPE type_netcdf_restart
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

    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_time
    INTEGER :: id_var_zeta
    INTEGER :: id_var_month

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '

  ! Restart file data - everything that's needed to restart a new run
  ! =================================================================

    ! Ice dynamics
    INTEGER :: id_var_Hi
    INTEGER :: id_var_Hb
    INTEGER :: id_var_Hs

    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '
    CHARACTER(LEN=256) :: name_var_Hb                    = 'Hb                   '
    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   '

    ! Thermodynamics
    INTEGER :: id_var_Ti

    CHARACTER(LEN=256) :: name_var_Ti                    = 'Ti                   '

    ! GIA
    INTEGER :: id_var_SL
    INTEGER :: id_var_dHb

    CHARACTER(LEN=256) :: name_var_SL                    = 'SL                   '
    CHARACTER(LEN=256) :: name_var_dHb                   = 'dHb                  '

    ! SMB
    INTEGER :: id_var_FirnDepth
    INTEGER :: id_var_MeltPreviousYear

    CHARACTER(LEN=256) :: name_var_FirnDepth             = 'FirnDepth            '
    CHARACTER(LEN=256) :: name_var_MeltPreviousYear      = 'MeltPreviousYear     '

    ! Isotopes
    INTEGER :: id_var_IsoIce

    CHARACTER(LEN=256) :: name_var_IsoIce                = 'IsoIce               '

    ! Inverse routine
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

  END TYPE type_netcdf_restart

  TYPE type_netcdf_help_fields
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! Index of time frame to be written to
    INTEGER :: ti

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_zeta
    INTEGER :: id_dim_time
    INTEGER :: id_dim_month
    INTEGER :: id_dim_z_ocean

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_dim_month                 = 'month                '
    CHARACTER(LEN=256) :: name_dim_z_ocean               = 'z_ocean              '

    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_time
    INTEGER :: id_var_zeta
    INTEGER :: id_var_month
    INTEGER :: id_var_z_ocean

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '
    CHARACTER(LEN=256) :: name_var_z_ocean               = 'z_ocean              '

  ! Variables
  ! =========

    INTEGER :: id_help_field_01
    INTEGER :: id_help_field_02
    INTEGER :: id_help_field_03
    INTEGER :: id_help_field_04
    INTEGER :: id_help_field_05
    INTEGER :: id_help_field_06
    INTEGER :: id_help_field_07
    INTEGER :: id_help_field_08
    INTEGER :: id_help_field_09
    INTEGER :: id_help_field_10
    INTEGER :: id_help_field_11
    INTEGER :: id_help_field_12
    INTEGER :: id_help_field_13
    INTEGER :: id_help_field_14
    INTEGER :: id_help_field_15
    INTEGER :: id_help_field_16
    INTEGER :: id_help_field_17
    INTEGER :: id_help_field_18
    INTEGER :: id_help_field_19
    INTEGER :: id_help_field_20
    INTEGER :: id_help_field_21
    INTEGER :: id_help_field_22
    INTEGER :: id_help_field_23
    INTEGER :: id_help_field_24
    INTEGER :: id_help_field_25
    INTEGER :: id_help_field_26
    INTEGER :: id_help_field_27
    INTEGER :: id_help_field_28
    INTEGER :: id_help_field_29
    INTEGER :: id_help_field_30
    INTEGER :: id_help_field_31
    INTEGER :: id_help_field_32
    INTEGER :: id_help_field_33
    INTEGER :: id_help_field_34
    INTEGER :: id_help_field_35
    INTEGER :: id_help_field_36
    INTEGER :: id_help_field_37
    INTEGER :: id_help_field_38
    INTEGER :: id_help_field_39
    INTEGER :: id_help_field_40
    INTEGER :: id_help_field_41
    INTEGER :: id_help_field_42
    INTEGER :: id_help_field_43
    INTEGER :: id_help_field_44
    INTEGER :: id_help_field_45
    INTEGER :: id_help_field_46
    INTEGER :: id_help_field_47
    INTEGER :: id_help_field_48
    INTEGER :: id_help_field_49
    INTEGER :: id_help_field_50

  END TYPE type_netcdf_help_fields

  TYPE type_netcdf_scalars_global
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

    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '

    INTEGER :: id_var_time

    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '

  ! Variables
  ! =========

    ! Global model stuff: sealevel, CO2, d18O components, etc.
    INTEGER :: id_var_GMSL
    INTEGER :: id_var_GMSL_NAM
    INTEGER :: id_var_GMSL_EAS
    INTEGER :: id_var_GMSL_GRL
    INTEGER :: id_var_GMSL_ANT
    INTEGER :: id_var_CO2_obs
    INTEGER :: id_var_CO2_mod
    INTEGER :: id_var_d18O_obs
    INTEGER :: id_var_d18O_mod
    INTEGER :: id_var_d18O_ice
    INTEGER :: id_var_d18O_Tdw
    INTEGER :: id_var_d18O_NAM
    INTEGER :: id_var_d18O_EAS
    INTEGER :: id_var_d18O_GRL
    INTEGER :: id_var_d18O_ANT
    INTEGER :: id_var_dT_glob
    INTEGER :: id_var_dT_dw

    CHARACTER(LEN=256) :: name_var_GMSL                  = 'GMSL                 '
    CHARACTER(LEN=256) :: name_var_GMSL_NAM              = 'GMSL_NAM             '
    CHARACTER(LEN=256) :: name_var_GMSL_EAS              = 'GMSL_EAS             '
    CHARACTER(LEN=256) :: name_var_GMSL_GRL              = 'GMSL_GRL             '
    CHARACTER(LEN=256) :: name_var_GMSL_ANT              = 'GMSL_ANT             '
    CHARACTER(LEN=256) :: name_var_CO2_obs               = 'CO2_obs              '
    CHARACTER(LEN=256) :: name_var_CO2_mod               = 'CO2_mod              '
    CHARACTER(LEN=256) :: name_var_d18O_obs              = 'd18O_obs             '
    CHARACTER(LEN=256) :: name_var_d18O_mod              = 'd18O_mod             '
    CHARACTER(LEN=256) :: name_var_d18O_ice              = 'd18O_ice             '
    CHARACTER(LEN=256) :: name_var_d18O_Tdw              = 'd18O_Tdw             '
    CHARACTER(LEN=256) :: name_var_d18O_NAM              = 'd18O_NAM             '
    CHARACTER(LEN=256) :: name_var_d18O_EAS              = 'd18O_EAS             '
    CHARACTER(LEN=256) :: name_var_d18O_GRL              = 'd18O_GRL             '
    CHARACTER(LEN=256) :: name_var_d18O_ANT              = 'd18O_ANT             '
    CHARACTER(LEN=256) :: name_var_dT_glob               = 'dT_glob              '
    CHARACTER(LEN=256) :: name_var_dT_dw                 = 'dT_dw                '

    ! Computation time for different model components
    INTEGER :: id_var_tcomp_total
    INTEGER :: id_var_tcomp_ice
    INTEGER :: id_var_tcomp_thermo
    INTEGER :: id_var_tcomp_climate
    INTEGER :: id_var_tcomp_GIA

    CHARACTER(LEN=256) :: name_var_tcomp_total           = 'tcomp_total          '
    CHARACTER(LEN=256) :: name_var_tcomp_ice             = 'tcomp_ice            '
    CHARACTER(LEN=256) :: name_var_tcomp_thermo          = 'tcomp_thermo         '
    CHARACTER(LEN=256) :: name_var_tcomp_climate         = 'tcomp_climate        '
    CHARACTER(LEN=256) :: name_var_tcomp_GIA             = 'tcomp_GIA            '

  END TYPE type_netcdf_scalars_global

  TYPE type_netcdf_scalars_regional
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

    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '

    INTEGER :: id_var_time

    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '

  ! Variables
  ! =========

    ! Regionally integrated stuff: ice volume, SMB components, etc.
    INTEGER :: id_var_ice_volume
    INTEGER :: id_var_ice_volume_af
    INTEGER :: id_var_ice_area
    INTEGER :: id_var_T2m
    INTEGER :: id_var_snowfall
    INTEGER :: id_var_rainfall
    INTEGER :: id_var_melt
    INTEGER :: id_var_refreezing
    INTEGER :: id_var_runoff
    INTEGER :: id_var_SMB
    INTEGER :: id_var_BMB
    INTEGER :: id_var_MB

    CHARACTER(LEN=256) :: name_var_ice_volume        = 'ice_volume           '
    CHARACTER(LEN=256) :: name_var_ice_volume_af     = 'ice_volume_af        '
    CHARACTER(LEN=256) :: name_var_ice_area          = 'ice_area             '
    CHARACTER(LEN=256) :: name_var_T2m               = 'T2m                  '
    CHARACTER(LEN=256) :: name_var_snowfall          = 'snowfall             '
    CHARACTER(LEN=256) :: name_var_rainfall          = 'rainfall             '
    CHARACTER(LEN=256) :: name_var_melt              = 'melt                 '
    CHARACTER(LEN=256) :: name_var_refreezing        = 'refreezing           '
    CHARACTER(LEN=256) :: name_var_runoff            = 'runoff               '
    CHARACTER(LEN=256) :: name_var_SMB               = 'SMB                  '
    CHARACTER(LEN=256) :: name_var_BMB               = 'BMB                  '
    CHARACTER(LEN=256) :: name_var_MB                = 'MB                   '

  END TYPE type_netcdf_scalars_regional

  TYPE type_netcdf_debug
    ! Integers describing open ports to different variables in an opened NetCDF file,
    ! plus character strings describing the names of those variables.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

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

    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_time
    INTEGER :: id_var_zeta
    INTEGER :: id_var_month

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '
    CHARACTER(LEN=256) :: name_var_zeta                  = 'zeta                 '
    CHARACTER(LEN=256) :: name_var_month                 = 'month                '

  ! Variables
  ! =========

    INTEGER :: id_var_int_2D_01
    INTEGER :: id_var_int_2D_02
    INTEGER :: id_var_int_2D_03
    INTEGER :: id_var_int_2D_04
    INTEGER :: id_var_int_2D_05
    INTEGER :: id_var_int_2D_06
    INTEGER :: id_var_int_2D_07
    INTEGER :: id_var_int_2D_08
    INTEGER :: id_var_int_2D_09
    INTEGER :: id_var_int_2D_10

    CHARACTER(LEN=256) :: name_var_int_2D_01        = 'int_2D_01      '
    CHARACTER(LEN=256) :: name_var_int_2D_02        = 'int_2D_02      '
    CHARACTER(LEN=256) :: name_var_int_2D_03        = 'int_2D_03      '
    CHARACTER(LEN=256) :: name_var_int_2D_04        = 'int_2D_04      '
    CHARACTER(LEN=256) :: name_var_int_2D_05        = 'int_2D_05      '
    CHARACTER(LEN=256) :: name_var_int_2D_06        = 'int_2D_06      '
    CHARACTER(LEN=256) :: name_var_int_2D_07        = 'int_2D_07      '
    CHARACTER(LEN=256) :: name_var_int_2D_08        = 'int_2D_08      '
    CHARACTER(LEN=256) :: name_var_int_2D_09        = 'int_2D_09      '
    CHARACTER(LEN=256) :: name_var_int_2D_10        = 'int_2D_10      '

    INTEGER :: id_var_dp_2D_01
    INTEGER :: id_var_dp_2D_02
    INTEGER :: id_var_dp_2D_03
    INTEGER :: id_var_dp_2D_04
    INTEGER :: id_var_dp_2D_05
    INTEGER :: id_var_dp_2D_06
    INTEGER :: id_var_dp_2D_07
    INTEGER :: id_var_dp_2D_08
    INTEGER :: id_var_dp_2D_09
    INTEGER :: id_var_dp_2D_10
    INTEGER :: id_var_dp_2D_11
    INTEGER :: id_var_dp_2D_12
    INTEGER :: id_var_dp_2D_13
    INTEGER :: id_var_dp_2D_14
    INTEGER :: id_var_dp_2D_15
    INTEGER :: id_var_dp_2D_16
    INTEGER :: id_var_dp_2D_17
    INTEGER :: id_var_dp_2D_18
    INTEGER :: id_var_dp_2D_19
    INTEGER :: id_var_dp_2D_20

    CHARACTER(LEN=256) :: name_var_dp_2D_01         = 'dp_2D_01       '
    CHARACTER(LEN=256) :: name_var_dp_2D_02         = 'dp_2D_02       '
    CHARACTER(LEN=256) :: name_var_dp_2D_03         = 'dp_2D_03       '
    CHARACTER(LEN=256) :: name_var_dp_2D_04         = 'dp_2D_04       '
    CHARACTER(LEN=256) :: name_var_dp_2D_05         = 'dp_2D_05       '
    CHARACTER(LEN=256) :: name_var_dp_2D_06         = 'dp_2D_06       '
    CHARACTER(LEN=256) :: name_var_dp_2D_07         = 'dp_2D_07       '
    CHARACTER(LEN=256) :: name_var_dp_2D_08         = 'dp_2D_08       '
    CHARACTER(LEN=256) :: name_var_dp_2D_09         = 'dp_2D_09       '
    CHARACTER(LEN=256) :: name_var_dp_2D_10         = 'dp_2D_10       '
    CHARACTER(LEN=256) :: name_var_dp_2D_11         = 'dp_2D_11       '
    CHARACTER(LEN=256) :: name_var_dp_2D_12         = 'dp_2D_12       '
    CHARACTER(LEN=256) :: name_var_dp_2D_13         = 'dp_2D_13       '
    CHARACTER(LEN=256) :: name_var_dp_2D_14         = 'dp_2D_14       '
    CHARACTER(LEN=256) :: name_var_dp_2D_15         = 'dp_2D_15       '
    CHARACTER(LEN=256) :: name_var_dp_2D_16         = 'dp_2D_16       '
    CHARACTER(LEN=256) :: name_var_dp_2D_17         = 'dp_2D_17       '
    CHARACTER(LEN=256) :: name_var_dp_2D_18         = 'dp_2D_18       '
    CHARACTER(LEN=256) :: name_var_dp_2D_19         = 'dp_2D_19       '
    CHARACTER(LEN=256) :: name_var_dp_2D_20         = 'dp_2D_20       '

    INTEGER :: id_var_dp_3D_01
    INTEGER :: id_var_dp_3D_02
    INTEGER :: id_var_dp_3D_03
    INTEGER :: id_var_dp_3D_04
    INTEGER :: id_var_dp_3D_05
    INTEGER :: id_var_dp_3D_06
    INTEGER :: id_var_dp_3D_07
    INTEGER :: id_var_dp_3D_08
    INTEGER :: id_var_dp_3D_09
    INTEGER :: id_var_dp_3D_10

    CHARACTER(LEN=256) :: name_var_dp_3D_01         = 'dp_3D_01      '
    CHARACTER(LEN=256) :: name_var_dp_3D_02         = 'dp_3D_02      '
    CHARACTER(LEN=256) :: name_var_dp_3D_03         = 'dp_3D_03      '
    CHARACTER(LEN=256) :: name_var_dp_3D_04         = 'dp_3D_04      '
    CHARACTER(LEN=256) :: name_var_dp_3D_05         = 'dp_3D_05      '
    CHARACTER(LEN=256) :: name_var_dp_3D_06         = 'dp_3D_06      '
    CHARACTER(LEN=256) :: name_var_dp_3D_07         = 'dp_3D_07      '
    CHARACTER(LEN=256) :: name_var_dp_3D_08         = 'dp_3D_08      '
    CHARACTER(LEN=256) :: name_var_dp_3D_09         = 'dp_3D_09      '
    CHARACTER(LEN=256) :: name_var_dp_3D_10         = 'dp_3D_10      '

    INTEGER :: id_var_dp_2D_monthly_01
    INTEGER :: id_var_dp_2D_monthly_02
    INTEGER :: id_var_dp_2D_monthly_03
    INTEGER :: id_var_dp_2D_monthly_04
    INTEGER :: id_var_dp_2D_monthly_05
    INTEGER :: id_var_dp_2D_monthly_06
    INTEGER :: id_var_dp_2D_monthly_07
    INTEGER :: id_var_dp_2D_monthly_08
    INTEGER :: id_var_dp_2D_monthly_09
    INTEGER :: id_var_dp_2D_monthly_10

    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_01 = 'dp_2D_monthly_01'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_02 = 'dp_2D_monthly_02'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_03 = 'dp_2D_monthly_03'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_04 = 'dp_2D_monthly_04'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_05 = 'dp_2D_monthly_05'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_06 = 'dp_2D_monthly_06'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_07 = 'dp_2D_monthly_07'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_08 = 'dp_2D_monthly_08'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_09 = 'dp_2D_monthly_09'
    CHARACTER(LEN=256) :: name_var_dp_2D_monthly_10 = 'dp_2D_monthly_10'

  END TYPE type_netcdf_debug

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

  TYPE type_netcdf_geothermal_heat_flux
    ! For reading an input file containing geothermal heat flux (e.g. Shapiro and Ritzwoller, 2004),
    ! describing geothermal heat flux at a lon-lat grid.

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

    ! ID's for variables:
    ! ===================

    ! Dimensions
    INTEGER :: id_dim_time
    INTEGER :: id_dim_month
    INTEGER :: id_dim_lon
    INTEGER :: id_dim_lat

    CHARACTER(LEN=256) :: name_dim_lon                   = 'Longitude            '
    CHARACTER(LEN=256) :: name_dim_lat                   = 'Latitude             '

    ! Variables
    INTEGER :: id_var_lon
    INTEGER :: id_var_lat
    INTEGER :: id_var_ghf

    CHARACTER(LEN=256) :: name_var_lon                   = 'Longitude            '
    CHARACTER(LEN=256) :: name_var_lat                   = 'Latitude             '
    CHARACTER(LEN=256) :: name_var_ghf                   = 'hflux                '

  END TYPE type_netcdf_geothermal_heat_flux

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

  TYPE type_netcdf_BIV_bed_roughness
    ! A NetCDF file containing bed roughness resulting from a basal inversion routine

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '

  ! Variables
  ! =========

    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_phi_fric
    INTEGER :: id_var_alpha_sq
    INTEGER :: id_var_beta_sq

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_phi_fric              = 'phi_fric             '
    CHARACTER(LEN=256) :: name_var_alpha_sq              = 'alpha_sq             '
    CHARACTER(LEN=256) :: name_var_beta_sq               = 'beta_sq              '

  END TYPE type_netcdf_BIV_bed_roughness

  TYPE type_netcdf_BIV_target_velocity
    ! A NetCDF file containing surface [u,v]-velocity fields to be used as
    ! the target in a velocity-based basal inversion routine

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '

  ! Variables
  ! =========

    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_u_surf
    INTEGER :: id_var_v_surf

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_u_surf                = 'u_surf               '
    CHARACTER(LEN=256) :: name_var_v_surf                = 'v_surf               '

  END TYPE type_netcdf_BIV_target_velocity

  TYPE type_netcdf_ISMIP_style_forcing
    ! NetCDF files containing the different ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing fields

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_time

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '

  ! Variables
  ! =========

    ! Dimensions
    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_time

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '

    ! Field variables
    INTEGER :: id_var_aSMB
    INTEGER :: id_var_dSMBdz
    INTEGER :: id_var_aST
    INTEGER :: id_var_dSTdz

    CHARACTER(LEN=256) :: name_var_aSMB                  = 'aSMB                 '
    CHARACTER(LEN=256) :: name_var_dSMBdz                = 'dSMBdz               '
    CHARACTER(LEN=256) :: name_var_aST                   = 'aST                  '
    CHARACTER(LEN=256) :: name_var_dSTdz                 = 'dSTdz                '

  END TYPE type_netcdf_ISMIP_style_forcing

  TYPE type_netcdf_ISMIP_style_baseline
    ! NetCDF file containing the baseline climate and orography for the ISMIP-style climate forcing

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '

  ! Variables
  ! =========

    ! Dimensions
    INTEGER :: id_var_x
    INTEGER :: id_var_y

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '

    ! Field variables
    INTEGER :: id_var_Hs
    INTEGER :: id_var_SMB
    INTEGER :: id_var_ST

    CHARACTER(LEN=256) :: name_var_Hs                    = 'Hs                   '
    CHARACTER(LEN=256) :: name_var_SMB                   = 'SMB                  '
    CHARACTER(LEN=256) :: name_var_ST                    = 'ST                   '

  END TYPE type_netcdf_ISMIP_style_baseline

  TYPE type_netcdf_prescribed_retreat_mask
    ! NetCDF files containing a prescribed retreat mask

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y
    INTEGER :: id_dim_time

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_dim_time                  = 'time                 '

  ! Variables
  ! =========

    ! Dimensions
    INTEGER :: id_var_x
    INTEGER :: id_var_y
    INTEGER :: id_var_time

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '
    CHARACTER(LEN=256) :: name_var_time                  = 'time                 '

    CHARACTER(LEN=256) :: time_units  ! ISMIP uses "days since XXX", paleo stuff just uses "years"

    ! Field variables
    INTEGER :: id_var_ice_fraction

  END TYPE type_netcdf_prescribed_retreat_mask

  TYPE type_netcdf_prescribed_retreat_mask_refice
    ! NetCDF files containing the reference ice thickness for a prescribed retreat mask

    CHARACTER(LEN=256) :: filename

    ! ID for NetCDF file:
    INTEGER :: ncid

  ! Dimensions
  ! ==========

    INTEGER :: id_dim_x
    INTEGER :: id_dim_y

    CHARACTER(LEN=256) :: name_dim_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_dim_y                     = 'y                    '

  ! Variables
  ! =========

    ! Dimensions
    INTEGER :: id_var_x
    INTEGER :: id_var_y

    CHARACTER(LEN=256) :: name_var_x                     = 'x                    '
    CHARACTER(LEN=256) :: name_var_y                     = 'y                    '

    ! Field variables
    INTEGER :: id_var_Hi

    CHARACTER(LEN=256) :: name_var_Hi                    = 'Hi                   '

  END TYPE type_netcdf_prescribed_retreat_mask_refice

CONTAINS

END MODULE data_types_netcdf_module
