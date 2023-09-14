MODULE configuration_module

  ! The way it's done right now:
  ! Each config variable has two versions: one with the "_config" extension, which is
  ! an actual variable in this module only, and one without the extension, which is
  ! a field in the "C" type. The "_config" variables are used to create a NAMELIST,
  ! which makes reading an external config file really easy - anything in the file that
  ! matches a variable in the namelist overwrites the default value. After that's done,
  ! the fields in the "C" type are replaced with the values of the "_config" variables,
  ! which now have either the default values, or those specified in the external config
  ! file.
  !
  ! While this is certainly very convenient when running the model (allowing for very
  ! short config files, where only the non-default variables are listed), it does make
  ! adding new config parameters a bit tedious - you have to add the "_config" variable,
  ! add it as a field in the "C" type, add it to the namelist, and let the "C" type field
  ! be overwritten in the end.
  !
  ! Some day I'll figure out a more elegant solution for this...

  USE mpi

  IMPLICIT NONE

  INTEGER, PARAMETER  :: dp  = KIND(1.0D0)  ! Kind of double precision numbers. Reals should be declared as: REAL(dp) :: example

  ! === Error messaging / debugging / profiling system ===

  CHARACTER(LEN=1024) :: routine_path
  INTEGER             :: n_MPI_windows

  TYPE subroutine_resource_tracker
    ! Track the resource use (computation time, memory) of a single subroutine
    CHARACTER(LEN = 1024) :: routine_path
    REAL(dp)              :: tstart, tcomp
    INTEGER               :: n_MPI_windows_init, n_MPI_windows_final
  END TYPE subroutine_resource_tracker

  TYPE( subroutine_resource_tracker), DIMENSION(:), ALLOCATABLE :: resource_tracker


  ! ===================================================================================
  ! "_config  variables, which will be collected into a NAMELIST, and possibly replaced
  ! by the values in the external config file. Remember the "_config" extension!
  ! ===================================================================================

  ! Time steps and range
  ! ====================

  REAL(dp)            :: start_time_of_run_config                    = 0.0_dp                           ! Start time (in years) of the simulations
  REAL(dp)            :: end_time_of_run_config                      = 50000.0_dp                       ! End   time (in years) of the simulations
  REAL(dp)            :: dt_coupling_config                          = 100._dp                          ! Interval of coupling (in years) between the four ice-sheets
  REAL(dp)            :: dt_max_config                               = 10.0_dp                          ! Maximum time step (in years) of the ice model
  REAL(dp)            :: dt_min_config                               = 0.01_dp                          ! Smallest allowed time step [yr]
  REAL(dp)            :: dt_startup_phase_config                     = 10._dp                           ! Length of time window (in years) after start_time and before end_time when dt = dt_min, to ensure smooth restarts
  REAL(dp)            :: dt_thermo_config                            = 10.0_dp                          ! Time step (in years) for updating thermodynamics
  REAL(dp)            :: dt_climate_config                           = 10._dp                           ! Time step (in years) for updating the climate
  REAL(dp)            :: dt_ocean_config                             = 10._dp                           ! Time step (in years) for updating the ocean
  REAL(dp)            :: dt_SMB_config                               = 10._dp                           ! Time step (in years) for updating the SMB
  REAL(dp)            :: dt_BMB_config                               = 10._dp                           ! Time step (in years) for updating the BMB
  REAL(dp)            :: dt_bedrock_ELRA_config                      = 100._dp                          ! Time step (in years) for updating the bedrock deformation rate with the ELRA model
  REAL(dp)            :: dt_SELEN_config                             = 1000._dp                         ! Time step (in years) for calling SELEN
  REAL(dp)            :: dt_output_config                            = 5000.0_dp                        ! Time step (in years) for writing help output
  REAL(dp)            :: dt_output_restart_config                    = 10000.0_dp                       ! Time step (in years) for writing restart output

  ! Which ice sheets do we simulate?
  ! ================================

  LOGICAL             :: do_NAM_config                               = .FALSE.                          ! North America
  LOGICAL             :: do_EAS_config                               = .FALSE.                          ! Eurasia
  LOGICAL             :: do_GRL_config                               = .FALSE.                          ! Greenland
  LOGICAL             :: do_ANT_config                               = .TRUE.                           ! Antarctica

  ! Benchmark experiments
  ! =====================

  ! SSA_icestream (see Schoof 2006, and also Bueler and Brown 2009)
  REAL(dp)            :: SSA_icestream_A_config                      = 1.0E-17_dp                       ! Glen's flow law factor (Pa yr^-1)
  REAL(dp)            :: SSA_icestream_L_config                      = 100000._dp                       ! Ice-stream half-width (m)
  REAL(dp)            :: SSA_icestream_m_config                      = 0.5_dp                           ! Determines the "smoothness" of the transition between the low-friction ice stream centre and the high-friction margin areas; low value (~0.5) = gradual, high value ( ~10) = sudden
  REAL(dp)            :: SSA_icestream_tantheta_config               = 0.0003_dp                        ! Slope
  REAL(dp)            :: SSA_icestream_H_config                      = 2000._dp                         ! Ice thickness (m)

  ! ISMIP-HOM (see Pattyn et al. 2008)
  REAL(dp)            :: ISMIP_HOM_L_config                          = 160000.0                         ! Domain size of the ISMIP-HOM benchmarks
  CHARACTER(LEN=256)  :: ISMIP_HOM_E_Arolla_filename_config          = 'arolla100.dat'                  ! Path to the Haut Glacier d'Arolla input file

  ! MISMIP+ (see Asay-Davis et al., 2016)
  LOGICAL             :: MISMIPplus_do_tune_A_for_GL_config          = .FALSE.                          ! Whether or not the flow factor A should be tuned for the GL position
  REAL(dp)            :: MISMIPplus_xGL_target_config                = 450000._dp                       ! Mid-channel GL position to tune the flow factor A for
  REAL(dp)            :: MISMIPplus_A_flow_initial_config            = 2.0E-17_dp                       ! Initial flow factor before tuning (or throughout the run when tuning is not used)
  CHARACTER(LEN=256)  :: MISMIPplus_scenario_config                  = ''                               ! Choose between the five MISMIP+  scenarios from Cornford   et al. (2020): ice0, ice1ra, ice1rr, ice2ra, ice2rr

  ! MISOMIP1 (see Asay-Davis et al., 2016)
  CHARACTER(LEN=256)  :: MISOMIP1_scenario_config                    = ''                               ! Choose between the four MISOMIP+ scenarios from Asay-Davis et al. (2016): IceOcean1ra, IceOcean1rr, IceOcean2ra, IceOcean2rr

  ! Whether or not to let IMAU_ICE dynamically create its own output folder.
  ! This works fine locally, on LISA its better to use a fixed folder name.
  ! =======================================================================

  LOGICAL             :: create_procedural_output_dir_config            = .TRUE.                           ! Automatically create an output directory with a procedural name (e.g. results_20210720_001/)
  CHARACTER(LEN=256)  :: fixed_output_dir_config                        = 'results_IMAU_ICE'               ! If not, create a directory with this name instead (stops the program if this directory already exists)
  CHARACTER(LEN=256)  :: fixed_output_dir_suffix_config                 = ''                               ! Suffix to put after the fixed output directory name, useful when doing ensemble runs with the template+variation set-up
  LOGICAL             :: do_write_regional_scalar_output_config         = .TRUE.
  LOGICAL             :: do_write_global_scalar_output_config           = .TRUE.
  LOGICAL             :: do_write_regional_scalar_every_timestep_config = .FALSE.

  ! Debugging
  ! =========

  LOGICAL             :: do_check_for_NaN_config                     = .FALSE.                          ! Whether or not fields should be checked for NaN values
  LOGICAL             :: do_time_display_config                      = .FALSE.                          ! Print current model time to screen

  ! ISMIP-style output
  ! ==================

  LOGICAL             :: do_write_ISMIP_output_config                = .FALSE.                          ! Whether or not to create a set of ISMIP output files
  CHARACTER(LEN=256)  :: ISMIP_output_group_code_config              = 'IMAU'                           ! Code for the group      name in the ISMIP output file names
  CHARACTER(LEN=256)  :: ISMIP_output_model_code_config              = 'IMAUICE'                        ! Code for the model      name in the ISMIP output file names
  CHARACTER(LEN=256)  :: ISMIP_output_experiment_code_config         = 'test'                           ! Code for the experiment name in the ISMIP output file names
  CHARACTER(LEN=256)  :: ISMIP_output_basetime_config                = 'YYYY-MM-DD'                     ! Basetime for the ISMIP output files (e.g. '1900-01-01')

  ! Grids
  ! =====

  ! North America
  REAL(dp)            :: lambda_M_NAM_config                         = 265._dp                          ! Longitude of the pole of the stereographic projection for the North America domain [degrees east]
  REAL(dp)            :: phi_M_NAM_config                            = 62._dp                           ! Latitude  of the pole of the stereographic projection for the North America domain [degrees north]
  REAL(dp)            :: beta_stereo_NAM_config                      = 71._dp                           ! Standard parallel     of the stereographic projection for the North America domain [degrees]
  REAL(dp)            :: xmin_NAM_config                             = -3600000._dp                     ! Western  boundary        of the North America domain [m]
  REAL(dp)            :: xmax_NAM_config                             =  3600000._dp                     ! Eastern  boundary     of the North America domain [m]
  REAL(dp)            :: ymin_NAM_config                             = -2400000._dp                     ! Southern boundary     of the North America domain [m]
  REAL(dp)            :: ymax_NAM_config                             =  2400000._dp                     ! Northern boundary     of the North America domain [m]
  REAL(dp)            :: dx_NAM_config                               = 40000._dp                        ! Horizontal resolution of the North America domain [m]

  ! Eurasia
  REAL(dp)            :: lambda_M_EAS_config                         = 40._dp                           ! Longitude of the pole of the stereographic projection for the Eurasia domain [degrees east]
  REAL(dp)            :: phi_M_EAS_config                            = 70._dp                           ! Latitude  of the pole of the stereographic projection for the Eurasia domain [degrees north]
  REAL(dp)            :: beta_stereo_EAS_config                      = 71._dp                           ! Standard parallel     of the stereographic projection for the Eurasia domain [degrees]
  REAL(dp)            :: xmin_EAS_config                             = -3400000._dp                     ! Western  boundary     of the Eurasia domain [m]
  REAL(dp)            :: xmax_EAS_config                             =  3400000._dp                     ! Eastern  boundary     of the Eurasia domain [m]
  REAL(dp)            :: ymin_EAS_config                             = -2080000._dp                     ! Southern boundary     of the Eurasia domain [m]
  REAL(dp)            :: ymax_EAS_config                             =  2080000._dp                     ! Northern boundary     of the Eurasia domain [m]
  REAL(dp)            :: dx_EAS_config                               = 40000._dp                        ! Horizontal resolution of the Eurasia domain [m]

  ! Greenland
  REAL(dp)            :: lambda_M_GRL_config                         = -45._dp                          ! Longitude of the pole of the stereographic projection for the Greenland domain [degrees east]
  REAL(dp)            :: phi_M_GRL_config                            = 90._dp                           ! Latitude  of the pole of the stereographic projection for the Greenland domain [degrees north]
  REAL(dp)            :: beta_stereo_GRL_config                      = 70._dp                           ! Standard parallel     of the stereographic projection for the Greenland domain [degrees]
  REAL(dp)            :: xmin_GRL_config                             =  -720000._dp                     ! Western  boundary     of the Greenland domain [m]
  REAL(dp)            :: xmax_GRL_config                             =   960000._dp                     ! Eastern  boundary     of the Greenland domain [m]
  REAL(dp)            :: ymin_GRL_config                             = -3450000._dp                     ! Southern boundary     of the Greenland domain [m]
  REAL(dp)            :: ymax_GRL_config                             =  -570000._dp                     ! Northern boundary     of the Greenland domain [m]
  REAL(dp)            :: dx_GRL_config                               = 20000._dp                        ! Horizontal resolution of the Greenland domain [m]

  ! Antarctica
  REAL(dp)            :: lambda_M_ANT_config                         = 0._dp                            ! Longitude of the pole of the stereographic projection for the Antarctica domain [degrees east]
  REAL(dp)            :: phi_M_ANT_config                            = -90._dp                          ! Latitude  of the pole of the stereographic projection for the Antarctica domain [degrees north]
  REAL(dp)            :: beta_stereo_ANT_config                      = 71._dp                           ! Standard parallel     of the stereographic projection for the Antarctica domain [degrees]
  REAL(dp)            :: xmin_ANT_config                             = -3300000._dp                     ! Western  boundary     of the Antarctica domain [m]
  REAL(dp)            :: xmax_ANT_config                             =  3300000._dp                     ! Eastern  boundary     of the Antarctica domain [m]
  REAL(dp)            :: ymin_ANT_config                             = -3300000._dp                     ! Southern boundary     of the Antarctica domain [m]
  REAL(dp)            :: ymax_ANT_config                             =  3300000._dp                     ! Northern boundary     of the Antarctica domain [m]
  REAL(dp)            :: dx_ANT_config                               = 40000._dp                        ! Horizontal resolution of the Antarctica domain [m]

  ! The scaled vertical coordinate zeta, used mainly in thermodynamics
  ! ==================================================================

  INTEGER                        :: nz_config                        = 15
  REAL(dp), DIMENSION(210), SAVE :: zeta_config                      = &
   (/0.00_dp, 0.10_dp, 0.20_dp, 0.30_dp, 0.40_dp, 0.50_dp, 0.60_dp, 0.70_dp, 0.80_dp, 0.90_dp, 0.925_dp, 0.95_dp, 0.975_dp, 0.99_dp, 1.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp, &
     0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp,  0.00_dp, 0.00_dp  /)

  ! Reference geometries (initial, present-day, and GIA equilibrium)
  ! ================================================================

  ! Initial geometry
  CHARACTER(LEN=256)  :: choice_refgeo_init_NAM_config               = 'realistic'                      ! Choice of initial geometry for North America; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_init_EAS_config               = 'realistic'                      ! Choice of initial geometry for Eurasia      ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_init_GRL_config               = 'realistic'                      ! Choice of initial geometry for Greenland    ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_init_ANT_config               = 'realistic'                      ! Choice of initial geometry for Antarctica   ; can be "idealised", "realistic", or "restart"
  REAL(dp)            :: time_to_restart_from_NAM_config             = 0._dp                            ! Can be different from C%start_time_of_run, though this will issue a warning
  REAL(dp)            :: time_to_restart_from_EAS_config             = 0._dp
  REAL(dp)            :: time_to_restart_from_GRL_config             = 0._dp
  REAL(dp)            :: time_to_restart_from_ANT_config             = 0._dp
  CHARACTER(LEN=256)  :: choice_refgeo_init_idealised_config         = 'flatearth'                      ! Choice of schematic initial geometry; see "generate_idealised_geometry" in reference_fields_module for options
  CHARACTER(LEN=256)  :: filename_refgeo_init_NAM_config             = '/Users/berends/Documents/Datasets/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_init_EAS_config             = '/Users/berends/Documents/Datasets/ETOPO1/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_init_GRL_config             = '/Users/berends/Documents/Datasets/Bedmachine_Greenland_v4/BedMachineGreenland-2021-04-20.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_init_ANT_config             = '/Users/berends/Documents/Datasets/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'

  ! Present-day geometry
  CHARACTER(LEN=256)  :: choice_refgeo_PD_NAM_config                 = 'realistic'                      ! Choice of present-day geometry for North America; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_PD_EAS_config                 = 'realistic'                      ! Choice of present-day geometry for Eurasia      ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_PD_GRL_config                 = 'realistic'                      ! Choice of present-day geometry for Greenland    ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_PD_ANT_config                 = 'realistic'                      ! Choice of present-day geometry for Antarctica   ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_PD_idealised_config           = 'flatearth'                      ! Choice of schematic present-day geometry; see "generate_idealised_geometry" in reference_fields_module for options
  CHARACTER(LEN=256)  :: filename_refgeo_PD_NAM_config               = '/Users/berends/Documents/Datasets/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_PD_EAS_config               = '/Users/berends/Documents/Datasets/ETOPO1/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_PD_GRL_config               = '/Users/berends/Documents/Datasets/Bedmachine_Greenland_v4/BedMachineGreenland-2021-04-20.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_PD_ANT_config               = '/Users/berends/Documents/Datasets/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'

  ! GIA equilibrium geometry
  CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_NAM_config              = 'realistic'                      ! Choice of GIA equilibrium geometry for North America; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_EAS_config              = 'realistic'                      ! Choice of GIA equilibrium geometry for Eurasia      ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_GRL_config              = 'realistic'                      ! Choice of GIA equilibrium geometry for Greenland    ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_ANT_config              = 'realistic'                      ! Choice of GIA equilibrium geometry for Antarctica   ; can be "idealised", "realistic", or "restart"
  CHARACTER(LEN=256)  :: choice_refgeo_GIAeq_idealised_config        = 'flatearth'                      ! Choice of schematic GIA equilibrium geometry; see "generate_idealised_geometry" in reference_fields_module for options
  CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_NAM_config            = '/Users/berends/Documents/Datasets/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_EAS_config            = '/Users/berends/Documents/Datasets/ETOPO1/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_GRL_config            = '/Users/berends/Documents/Datasets/Bedmachine_Greenland_v4/BedMachineGreenland-2021-04-20.nc'
  CHARACTER(LEN=256)  :: filename_refgeo_GIAeq_ANT_config            = '/Users/berends/Documents/Datasets/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'

  LOGICAL             :: remove_Lake_Vostok_config                   = .TRUE.

  ! Global forcing (insolation, CO2, d18O, geothermal heat flux)
  ! ============================================================

  ! Possible choice_forcing_method options:
  ! 'none'                 : No global forcing used at all; climate or SMB are fully parameterised or directly prescribed
  ! 'd18O_inverse_dT_glob' : Use the inverse routine with the specified d18O record to calculate a global temperature offset (e.g. de Boer et al., 2013)
  ! 'CO2_direct'           : Use the specified CO2 record to force the climate matrix (e.g. Berends et al., 2018)
  ! 'd18O_inverse_CO2'     : Use the inverse routine with the specified d18O record to calculate CO2 and then force the climate matrix (e.g. Berends et al., 2019)
  CHARACTER(LEN=256)  :: choice_forcing_method_config                = 'd18O_inverse_dT_glob'

  ! Insolation forcing (NetCDF)
  CHARACTER(LEN=256)  :: choice_insolation_forcing_config            = 'realistic'                      ! Choice of insolation forcing: "none", "static", "realistic"
  REAL(dp)            :: static_insolation_time_config               = 0._dp                            ! Keep insolation values fixed to this time when choice_insolation_forcing = 'static'
  CHARACTER(LEN=256)  :: filename_insolation_config                  = '/Users/berends/Documents/Datasets/Insolation/Laskar_etal_2004_insolation.nc'

  ! CO2 record (ASCII text file in kyr, so the number of rows needs to be specified)
  CHARACTER(LEN=256)  :: filename_CO2_record_config                  = '/Users/berends/Documents/Datasets/CO2/EPICA_CO2_Bereiter_2015_100yr.dat'
  INTEGER             :: CO2_record_length_config                    = 8001

  ! d18O record (ASCII text file, so the number of rows needs to be specified)
  CHARACTER(LEN=256)  :: filename_d18O_record_config                 = '/Users/berends/Documents/Datasets/d18O/Ahn2017_d18O.dat'
  INTEGER             :: d18O_record_length_config                   = 2051

  ! Parameters for combining insolation and CO2 in the matrix method
  LOGICAL             :: do_combine_CO2_and_insolation_config        = .FALSE.                          ! Combine the effect of insolation and CO2 to calculate the external forcing in the matrix method
  REAL(dp)            :: insolation_weigth_mean_NAM_config           = 440                              ! (W/m2) insolation at which insolation does not alter the external forcing (w_INS = 0)
  REAL(dp)            :: insolation_weigth_amplitude_NAM_config      = 70                               ! (W/m2) the amplitude at which insolation affects the interpolation weight.
  REAL(dp)            :: insolation_weigth_mean_EAS_config           = 440                              ! (W/m2) insolation at which insolation does not alter the external forcing (w_INS = 0)
  REAL(dp)            :: insolation_weigth_amplitude_EAS_config      = 70                               ! (W/m2) the amplitude at which insolation affects the interpolation weight.
  REAL(dp)            :: insolation_weigth_mean_GRL_config           = 440                              ! (W/m2) insolation at which insolation does not alter the external forcing (w_INS = 0)
  REAL(dp)            :: insolation_weigth_amplitude_GRL_config      = 70                               ! (W/m2) the amplitude at which insolation affects the interpolation weight.
  REAL(dp)            :: insolation_weigth_mean_ANT_config           = 440                              ! (W/m2) insolation at which insolation does not alter the external forcing (w_INS = 0)
  REAL(dp)            :: insolation_weigth_amplitude_ANT_config      = 70                               ! (W/m2) the amplitude at which insolation affects the interpolation weight.

  ! Determine if the GCM wind is used in the climate matrix method
  LOGICAL             :: do_climate_matrix_wind_config               = .TRUE.                            ! If TRUE, use the wind from the climate forcing. If FALSE, use the wind from the reference climate.

  ! Geothermal heat flux
  CHARACTER(LEN=256)  :: choice_geothermal_heat_flux_config          = 'spatial'                        ! Choice of geothermal heat flux; can be 'constant' or 'spatial'
  REAL(dp)            :: constant_geothermal_heat_flux_config        = 1.72E06_dp                       ! Geothermal Heat flux [J m^-2 yr^-1] Sclater et al. (1980)
  CHARACTER(LEN=256)  :: filename_geothermal_heat_flux_config        = '/Users/berends/Documents/Datasets/GHF/geothermal_heatflux_ShapiroRitzwoller2004_global_1x1_deg.nc'

  ! Parameters for calculating modelled benthic d18O
  LOGICAL             :: do_calculate_benthic_d18O_config            = .FALSE.                           ! Whether or not to calculate modelled benthic d18O (set to .FALSE. for e.g. idealised-geometry experiments, future projections)
  REAL(dp)            :: dT_deepwater_averaging_window_config        = 3000                             ! Time window (in yr) over which global mean temperature anomaly is averaged to find the deep-water temperature anomaly
  REAL(dp)            :: dT_deepwater_dT_surf_ratio_config           = 0.25_dp                          ! Ratio between global mean surface temperature change and deep-water temperature change
  REAL(dp)            :: d18O_dT_deepwater_ratio_config              = -0.28_dp                         ! Ratio between deep-water temperature change and benthic d18O change

  ! Parameters for the inverse routine
  REAL(dp)            :: dT_glob_inverse_averaging_window_config     = 2000._dp                         ! Time window (in yr) over which global mean temperature anomaly is averaged before changing it with the inverse routine
  REAL(dp)            :: inverse_d18O_to_dT_glob_scaling_config      = 20._dp                           ! Scaling factor between modelled d18O anomaly and prescribed temperature anomaly change (value from de Boer et al., 2013)
  REAL(dp)            :: CO2_inverse_averaging_window_config         = 2000._dp                         ! Time window (in yr) over which CO2                             is averaged before changing it with the inverse routine
  REAL(dp)            :: inverse_d18O_to_CO2_scaling_config          = 68._dp                           ! Scaling factor between modelled d18O anomaly and modelled CO2 change (value from Berends et al., 2019)
  REAL(dp)            :: inverse_d18O_to_CO2_initial_CO2_config      = 280._dp                          ! CO2 value at the start of the simulation when using the inverse method to calculate CO2

  ! Initial file for the inverse d18O
  CHARACTER(LEN=256)  :: choice_d18O_inverse_init_NAM_config       =  'init'
  CHARACTER(LEN=256)  :: choice_d18O_inverse_init_EAS_config       =  'init'
  CHARACTER(LEN=256)  :: choice_d18O_inverse_init_GRL_config       =  'init'
  CHARACTER(LEN=256)  :: choice_d18O_inverse_init_ANT_config       =  'init'

  CHARACTER(LEN=256)  :: filename_d18O_inverse_init_NAM_config       = ''
  CHARACTER(LEN=256)  :: filename_d18O_inverse_init_EAS_config       = ''
  CHARACTER(LEN=256)  :: filename_d18O_inverse_init_GRL_config       = ''
  CHARACTER(LEN=256)  :: filename_d18O_inverse_init_ANT_config       = ''

  ! Ice dynamics - velocity
  ! =======================

  CHARACTER(LEN=256)  :: choice_ice_dynamics_config                  = 'DIVA'                           ! Choice of ice-dynamica approximation: "none" (= fixed geometry), "SIA", "SSA", "SIA/SSA", "DIVA"
  REAL(dp)            :: n_flow_config                               = 3.0_dp                           ! Exponent in Glen's flow law
  REAL(dp)            :: m_enh_sheet_config                          = 1.0_dp                           ! Ice flow enhancement factor for grounded ice
  REAL(dp)            :: m_enh_shelf_config                          = 1.0_dp                           ! Ice flow enhancement factor for floating ice
  CHARACTER(LEN=256)  :: choice_ice_margin_config                    = 'infinite_slab'                  ! Choice of ice margin boundary conditions: "BC", "infinite_slab"
  LOGICAL             :: include_SSADIVA_crossterms_config           = .TRUE.                           ! Whether or not to include the "cross-terms" of the SSA/DIVA
  LOGICAL             :: do_GL_subgrid_friction_config               = .TRUE.                           ! Whether or not to scale basal friction with the sub-grid grounded fraction (needed to get proper GL migration; only turn this off for showing the effect on the MISMIP_mod results!)
  LOGICAL             :: do_smooth_geometry_config                   = .FALSE.                          ! Whether or not to smooth the model geometry (bedrock + initial ice thickness)
  REAL(dp)            :: r_smooth_geometry_config                    = 0.5_dp                           ! Geometry smoothing radius (in number of grid cells)

  ! Some parameters for numerically solving the SSA/DIVA
  REAL(dp)            :: DIVA_visc_it_norm_dUV_tol_config            = 1E-2_dp                          ! Successive solutions of UV in the effective viscosity iteration must not differ by more than this amount (on average)
  INTEGER             :: DIVA_visc_it_nit_config                     = 50                               ! Maximum number of effective viscosity iterations
  REAL(dp)            :: DIVA_visc_it_relax_config                   = 0.4_dp                           ! Relaxation parameter for subsequent viscosity iterations (for improved stability)
  REAL(dp)            :: DIVA_epsilon_sq_0_config                    = 1E-15_dp                         ! Normalisation term so that zero velocity gives non-zero viscosity
  REAL(dp)            :: DIVA_visc_eff_min_config                    = 1E3_dp                           ! Minimum value for effective viscosity
  REAL(dp)            :: DIVA_beta_max_config                        = 1E20_dp                          ! Maximum value for effective basal friction
  REAL(dp)            :: DIVA_beta_min_thin_ice_config               = 0.0_dp                           ! Minimum value for effective basal friction over thin-ice areas
  REAL(dp)            :: DIVA_beta_min_thin_ice_H0_config            = 1E20_dp                          ! Threshold thinness below which a minimum beta limit starts to apply
  REAL(dp)            :: DIVA_vel_max_config                         = 5000._dp                         ! DIVA velocities are limited to this value
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_u_west_config              = 'infinite'                       ! Boundary conditions for the ice velocity field at the domain boundary in the DIVA
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_u_east_config              = 'infinite'                       ! Allowed choices: "infinite", "periodic", "zero"
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_u_south_config             = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_u_north_config             = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_v_west_config              = 'infinite'                       ! Boundary conditions for the ice velocity field at the domain boundary in the DIVA
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_v_east_config              = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_v_south_config             = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_v_north_config             = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_choice_matrix_solver_config            = 'PETSc'                          ! Choice of matrix solver for the ice velocity equations: "SOR", "PETSc"
  INTEGER             :: DIVA_SOR_nit_config                         = 10000                            ! DIVA SOR   solver - maximum number of iterations
  REAL(dp)            :: DIVA_SOR_tol_config                         = 2.5_dp                           ! DIVA SOR   solver - stop criterion, absolute difference
  REAL(dp)            :: DIVA_SOR_omega_config                       = 1.3_dp                           ! DIVA SOR   solver - over-relaxation parameter
  REAL(dp)            :: DIVA_PETSc_rtol_config                      = 0.01_dp                          ! DIVA PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
  REAL(dp)            :: DIVA_PETSc_abstol_config                    = 2.5_dp                           ! DIVA PETSc solver - stop criterion, absolute difference
  LOGICAL             :: do_read_velocities_from_restart_config      = .FALSE.                          ! Whether or not to read velocities from restart file to initialise model

  ! Ice dynamics - time integration
  ! ===============================

  CHARACTER(LEN=256)  :: choice_timestepping_config                  = 'pc'                             ! Choice of timestepping method: "direct", "pc" (NOTE: 'direct' does not work with DIVA ice dynamcis!)
  CHARACTER(LEN=256)  :: choice_ice_integration_method_config        = 'explicit'                       ! Choice of ice thickness integration scheme: "none" (i.e. unchanging geometry), "explicit", "semi-implicit"
  CHARACTER(LEN=256)  :: dHi_choice_matrix_solver_config             = 'SOR'                            ! Choice of matrix solver for the semi-implicit ice thickness equation: "SOR", "PETSc"
  INTEGER             :: dHi_SOR_nit_config                          = 3000                             ! dHi SOR   solver - maximum number of iterations
  REAL(dp)            :: dHi_SOR_tol_config                          = 2.5_dp                           ! dHi SOR   solver - stop criterion, absolute difference
  REAL(dp)            :: dHi_SOR_omega_config                        = 1.3_dp                           ! dHi SOR   solver - over-relaxation parameter
  REAL(dp)            :: dHi_PETSc_rtol_config                       = 0.001_dp                         ! dHi PETSc solver - stop criterion, relative difference (iteration stops if rtol OR abstol is reached)
  REAL(dp)            :: dHi_PETSc_abstol_config                     = 0.001_dp                         ! dHi PETSc solver - stop criterion, absolute difference

  ! Predictor-corrector ice-thickness update
  REAL(dp)            :: pc_epsilon_config                           = 0.5_dp                           ! Target truncation error in dHi_dt [m/yr] (epsilon in Robinson et al., 2020, Eq. 33)
  REAL(dp)            :: pc_k_I_config                               = 0.2_dp                           ! Exponent k_I in  Robinson et al., 2020, Eq. 33
  REAL(dp)            :: pc_k_p_config                               = 0.2_dp                           ! Exponent k_p in  Robinson et al., 2020, Eq. 33
  REAL(dp)            :: pc_eta_min_config                           = 1E-8_dp                          ! Normalisation term in estimation of the truncation error (Robinson et al., Eq. 32)
  INTEGER             :: pc_max_timestep_iterations_config           = 5                                ! Maximum number of iterations of each time step
  REAL(dp)            :: pc_redo_tol_config                          = 10._dp                           ! Maximum allowed truncation error (any higher and the timestep is decreased)

  ! Ice thickness boundary conditions
  CHARACTER(LEN=256)  :: ice_thickness_west_BC_config                = 'zero'                           ! Choice of boundary conditions for ice thickness at the domain boundary: "infinite", "periodic", "zero", "fixed"
  CHARACTER(LEN=256)  :: ice_thickness_east_BC_config                = 'zero'
  CHARACTER(LEN=256)  :: ice_thickness_south_BC_config               = 'zero'
  CHARACTER(LEN=256)  :: ice_thickness_north_BC_config               = 'zero'
  CHARACTER(LEN=256)  :: choice_mask_noice_NAM_config                = 'NAM_remove_GRL'                 ! Choice of mask_noice configuration
  CHARACTER(LEN=256)  :: choice_mask_noice_EAS_config                = 'EAS_remove_GRL'
  CHARACTER(LEN=256)  :: choice_mask_noice_GRL_config                = 'GRL_remove_Ellesmere'
  CHARACTER(LEN=256)  :: choice_mask_noice_ANT_config                = 'none'                           ! For Antarctica, additional choices are included for certain idealised-geometry experiments: "MISMIP_mod", "MISMIP+"

  ! Fixed/delayed ice thickness evolution
  REAL(dp)            :: fixed_sheet_geometry_config                 = 0.0_dp                           ! Fix (1), release (0), or delay grounded ice geometry evolution
  REAL(dp)            :: fixed_shelf_geometry_config                 = 0.0_dp                           ! Fix (1), release (0), or delay floating ice geometry evolution
  REAL(dp)            :: fixed_grounding_line_g_config               = 0.0_dp                           ! Fix (1), release (0), or delay GL geometry evolution (grounded side)
  REAL(dp)            :: fixed_grounding_line_f_config               = 0.0_dp                           ! Fix (1), release (0), or delay GL geometry evolution (floating side)
  REAL(dp)            :: fixed_decay_t_start_config                  = 0.0_dp                           ! Start time of linear transition between fixed/delayed and free evolution
  REAL(dp)            :: fixed_decay_t_end_config                    = 0.0_dp                           ! End   time of linear transition between fixed/delayed and free evolution
  REAL(dp)            :: relax_thick_t_start_config                  = -9.9e99_dp                       ! Start time of ice thickness relaxation: no fixiness/MB; refgeo_PD set to relaxed geometry
  REAL(dp)            :: relax_thick_t_end_config                    = -8.8e88_dp                       ! End   time of ice thickness relaxation: no fixiness/MB; refgeo_PD set to relaxed geometry

  ! Target dHi_dt during model spinup
  LOGICAL             :: do_target_dhdt_config                       = .FALSE.                          ! Whether or not to perform spinup using a target dHi_dt field
  CHARACTER(LEN=256)  :: target_dhdt_filename_config                 = 'dhdt_target.nc'                 ! NetCDF file containing target dHi_dt for model spinup
  REAL(dp)            :: target_dhdt_t_end_config                    = -9.9e99_dp                       ! Time at which the target dHi_dt field is released into the wild

  ! Ice dynamics - basal conditions and sliding
  ! ===========================================

  ! Sliding laws
  CHARACTER(LEN=256)  :: choice_sliding_law_config                   = 'Coulomb_regularised'            ! Choice of sliding law: "no_sliding", "idealised", "Coulomb", "Coulomb_regularised", "Weertman", "Tsai2015", "Schoof2005", "Zoet-Iverson"
  CHARACTER(LEN=256)  :: choice_idealised_sliding_law_config         = ''                               ! "ISMIP_HOM_C", "ISMIP_HOM_D", "ISMIP_HOM_E", "ISMIP_HOM_F"
  REAL(dp)            :: slid_delta_v_config                         = 1.0E-3_dp                        ! Normalisation parameter to prevent errors when velocity is zero
  REAL(dp)            :: slid_Weertman_m_config                      = 3._dp                            ! Exponent in Weertman sliding law
  REAL(dp)            :: slid_Coulomb_reg_q_plastic_config           = 0.3_dp                           ! Scaling exponent   in regularised Coulomb sliding law
  REAL(dp)            :: slid_Coulomb_reg_u_threshold_config         = 100._dp                          ! Threshold velocity in regularised Coulomb sliding law
  REAL(dp)            :: slid_ZI_ut_config                           = 200._dp                          ! (uniform) transition velocity used in the Zoet-Iverson sliding law [m/yr]
  REAL(dp)            :: slid_ZI_p_config                            = 5._dp                            ! Velocity exponent             used in the Zoet-Iverson sliding law
  LOGICAL             :: do_slid_ZI_no_angle_config                  = .FALSE.                          ! If .TRUE., use phi_fric as a fraction (i.e./e.g. [0 1]) instead of an angle. Adjust its limits accordingly elsewhere.
  LOGICAL             :: include_basal_freezing_config               = .TRUE.                           ! If .TRUE., no basal sliding is allowed when the basal temperature is more than [deltaT_basal_freezing] below the pressure melting point
  REAL(dp)            :: deltaT_basal_freezing_config                = 2._dp                            ! See above.
  REAL(dp)            :: subgrid_friction_exponent_config            = 2._dp                            ! Exponent to which f_grnd should be raised before being used to scale beta

  ! Basal hydrology
  CHARACTER(LEN=256)  :: choice_basal_hydrology_config               = 'Martin2011'                     ! Choice of basal conditions: "saturated", "Martin2011"
  REAL(dp)            :: Martin2011_hydro_N_lim_config               = 0.96_dp                          ! Martin et al. (2011) basal hydrology model: limit pore water pressure w.r.t. overburden [1=100% allowed; 0=no hydrology; 0.96 used in ref. paper]
  REAL(dp)            :: Martin2011_hydro_Hb_min_config              = 0._dp                            ! Martin et al. (2011) basal hydrology model: low-end  Hb  value of bedrock-dependent pore-water pressure
  REAL(dp)            :: Martin2011_hydro_Hb_max_config              = 1000._dp                         ! Martin et al. (2011) basal hydrology model: high-end Hb  value of bedrock-dependent pore-water pressure

  ! Basal roughness / friction
  CHARACTER(LEN=256)  :: choice_basal_roughness_config               = 'parameterised'                  ! "uniform"", parameterised", "prescribed"
  REAL(dp)            :: uniform_Weertman_beta_sq_config             = 1.0E4_dp                         ! Uniform value for beta_sq  in Weertman sliding law
  REAL(dp)            :: uniform_Coulomb_phi_fric_config             = 15._dp                           ! Uniform value for phi_fric in (regularised) Coulomb sliding law
  REAL(dp)            :: uniform_Tsai2015_alpha_sq_config            = 0.5_dp                           ! Uniform value for alpha_sq in Tsai2015 sliding law
  REAL(dp)            :: uniform_Tsai2015_beta_sq_config             = 1.0E4_dp                         ! Uniform value for beta_sq  in Tsai2015 sliding law
  REAL(dp)            :: uniform_Schoof2005_alpha_sq_config          = 0.5_dp                           ! Uniform value for alpha_sq in Schoof2005 sliding law
  REAL(dp)            :: uniform_Schoof2005_beta_sq_config           = 1.0E4_dp                         ! Uniform value for beta_sq  in Schoof2005 sliding law
  CHARACTER(LEN=256)  :: choice_param_basal_roughness_config         = 'Martin2011'                     ! "Martin2011", "SSA_icestream", "MISMIP+"
  REAL(dp)            :: Martin2011till_phi_Hb_min_config            = -1000._dp                        ! Martin et al. (2011) bed roughness model: low-end  Hb  value of bedrock-dependent till friction angle
  REAL(dp)            :: Martin2011till_phi_Hb_max_config            = 0._dp                            ! Martin et al. (2011) bed roughness model: high-end Hb  value of bedrock-dependent till friction angle
  REAL(dp)            :: Martin2011till_phi_min_config               = 5._dp                            ! Martin et al. (2011) bed roughness model: low-end  phi value of bedrock-dependent till friction angle
  REAL(dp)            :: Martin2011till_phi_max_config               = 20._dp                           ! Martin et al. (2011) bed roughness model: high-end phi value of bedrock-dependent till friction angle
  CHARACTER(LEN=256)  :: basal_roughness_filename_config             = ''                               ! NetCDF file containing a basal roughness field for the chosen sliding law
  LOGICAL             :: do_smooth_phi_restart_config                = .FALSE.                          ! Whether or not to smooth the prescribed bed roughness once (crucial for downscaling runs)
  REAL(dp)            :: r_smooth_phi_restart_config                 = 0.5_dp                           ! Prescribed bed roughness smoothing radius (in number of grid cells)

  ! Basal inversion
  LOGICAL             :: do_BIVgeo_config                            = .FALSE.                          ! Whether or not to perform a geometry-based basal inversion (following Pollard & DeConto, 2012)
  REAL(dp)            :: BIVgeo_t_start_config                       = -9.9E9_dp                        ! Minimum model time when the inversion is allowed
  REAL(dp)            :: BIVgeo_t_end_config                         = +9.9E9_dp                        ! Maximum model time when the inversion is allowed
  CHARACTER(LEN=256)  :: choice_BIVgeo_method_config                 = 'Berends2022'                    ! Choice of geometry-based inversion method: "PDC2012", "Lipscomb2021", "CISM+", "Berends2022", "Bernales2017"
  REAL(dp)            :: BIVgeo_dt_config                            = 5._dp                            ! Time step      for bed roughness updates in the PDC2012 geometry-based basal inversion method [yr]
  REAL(dp)            :: BIVgeo_PDC2012_hinv_config                  = 500._dp                          ! Scaling factor for bed roughness updates in the PDC2012 geometry-based basal inversion method [m]
  REAL(dp)            :: BIVgeo_Lipscomb2021_tauc_config             = 500._dp                          ! Timescale       in the Lipscomb2021 geometry-based basal inversion method [yr]
  REAL(dp)            :: BIVgeo_Lipscomb2021_H0_config               = 100._dp                          ! Thickness scale in the Lipscomb2021 geometry-based basal inversion method [m]
  REAL(dp)            :: BIVgeo_CISMplus_wH_config                   = 1._dp                            ! Weighting factor for ice thickness in the CISM+ geometry/velocity-based basal inversion method
  REAL(dp)            :: BIVgeo_CISMplus_wu_config                   = 1._dp                            ! Weighting factor for velocity      in the CISM+ geometry/velocity-based basal inversion method
  REAL(dp)            :: BIVgeo_CISMplus_tauc_config                 = 500._dp                          ! Timescale       in the CISM+ geometry/velocity-based basal inversion method [yr]
  REAL(dp)            :: BIVgeo_CISMplus_H0_config                   = 100._dp                          ! Thickness scale in the CISM+ geometry/velocity-based basal inversion method [m]
  REAL(dp)            :: BIVgeo_CISMplus_u0_config                   = 10._dp                           ! Velocity  scale in the CISM+ geometry/velocity-based basal inversion method [m/yr]
  REAL(dp)            :: BIVgeo_Berends2022_tauc_config              = 10._dp                           ! Timescale       in the Berends2022 geometry/velocity-based basal inversion method [yr]
  REAL(dp)            :: BIVgeo_Berends2022_H0_config                = 100._dp                          ! First  thickness scale in the Berends2022 geometry/velocity-based basal inversion method [m]
  REAL(dp)            :: BIVgeo_Berends2022_u0_config                = 250._dp                          ! First  velocity  scale in the Berends2022 geometry/velocity-based basal inversion method [m/yr]
  REAL(dp)            :: BIVgeo_Berends2022_Hi_scale_config          = 300._dp                          ! Second thickness scale in the Berends2022 geometry/velocity-based basal inversion method [m]
  REAL(dp)            :: BIVgeo_Berends2022_u_scale_config           = 3000._dp                         ! Second velocity  scale in the Berends2022 geometry/velocity-based basal inversion method [m/yr]
  REAL(dp)            :: BIVgeo_Berends2022_phimin_config            = 0.1_dp                           ! Smallest allowed value for the inverted till friction angle phi
  REAL(dp)            :: BIVgeo_Berends2022_phimax_config            = 30._dp                           ! Largest  allowed value for the inverted till friction angle phi
  CHARACTER(LEN=256)  :: BIVgeo_target_velocity_filename_config      = ''                               ! NetCDF file where the target velocities are read in the CISM+ and Berends2022 geometry/velocity-based basal inversion methods
  CHARACTER(LEN=256)  :: BIVgeo_filename_output_config               = 'bed_roughness_inv.nc'           ! NetCDF file where the final inverted basal roughness will be saved
  REAL(dp)            :: BIVgeo_Bernales2017_hinv_config             = 100._dp                          ! Scaling factor for bed roughness updates in the Bernales (2017) geometry-based basal inversion method [m]
  REAL(dp)            :: BIVgeo_Bernales2017_tol_diff_config         = 100._dp                          ! Minimum ice thickness difference [m] that triggers inversion (.OR. &)
  REAL(dp)            :: BIVgeo_Bernales2017_tol_frac_config         = 1.0_dp                           ! Minimum ratio between ice thickness difference and reference value that triggers inversion
  REAL(dp)            :: BIVgeo_Pien2023_H0_config                   = 200.0_dp                         ! Ice thickness scale for regularisation in Pien's friction inversion method
  REAL(dp)            :: BIVgeo_Pien2023_tau_config                  = 200.0_dp                         ! Time scale for regularisation in Pien's friction inversion method
  REAL(dp)            :: BIVgeo_Pien2023_r_config                    = 0.0_dp                           ! Relaxation scale for regularisation in Pien's friction inversion method
  REAL(dp)            :: BIVgeo_Pien2023_lowerHb_config              = -1000.0_dp                       ! Lower bedrock limit for target relaxation friction in Pien's friction inversion method
  REAL(dp)            :: BIVgeo_Pien2023_upperHb_config              = +1000.0_dp                       ! Upper bedrock limit for target relaxation friction in Pien's friction inversion method
  REAL(dp)            :: BIVgeo_Pien2023_min_config                  = 0.1_dp                           ! Lower till friction angle limit in Pien's friction inversion method
  REAL(dp)            :: BIVgeo_Pien2023_max_config                  = 30.0_dp                          ! Upper till friction angle limit in Pien's friction inversion method

  ! Ice dynamics - calving
  ! ======================

  CHARACTER(LEN=256)  :: choice_calving_law_config                   = 'threshold_thickness'            ! Choice of calving law: "none", "threshold_thickness"
  REAL(dp)            :: calving_threshold_thickness_config          = 200._dp                          ! Threshold ice thickness in the "threshold_thickness" calving law (200m taken from ANICE)
  LOGICAL             :: do_remove_shelves_config                    = .FALSE.                          ! If set to TRUE, all floating ice is always instantly removed (used in the ABUMIP-ABUK experiment)
  LOGICAL             :: remove_shelves_larger_than_PD_config        = .FALSE.                          ! If set to TRUE, all floating ice beyond the present-day calving front is removed (used for some Antarctic spin-ups)
  LOGICAL             :: continental_shelf_calving_config            = .FALSE.                          ! If set to TRUE, all ice beyond the continental shelf edge (set by a maximum depth) is removed
  REAL(dp)            :: continental_shelf_min_height_config         = -2000._dp                        ! Maximum depth of the continental shelf

  ! Thermodynamics and rheology
  ! ===========================

  CHARACTER(LEN=256)  :: choice_initial_ice_temperature_config       = 'Robin'                          ! Choice of initial ice temperature profile: "uniform", "linear", "Robin", "restart"
  REAL(dp)            :: uniform_ice_temperature_config              = 270._dp                          ! Uniform ice temperature (applied when choice_initial_ice_temperature_config = "uniform")
  CHARACTER(LEN=256)  :: choice_thermo_model_config                  = '3D_heat_equation'               ! Choice of thermodynamical model: "none", "3D_heat_equation"
  CHARACTER(LEN=256)  :: choice_ice_rheology_config                  = 'Huybrechts1992'                 ! Choice of ice rheology model: "uniform", "Huybrechts1992", "MISMIP_mod"
  REAL(dp)            :: uniform_flow_factor_config                  = 1E-16_dp                         ! Uniform ice flow factor (applied when choice_ice_rheology_model_config = "uniform")
  CHARACTER(LEN=256)  :: choice_ice_heat_capacity_config             = 'Pounder1965'                    ! Choice of ice heat capacity model: "uniform", "Pounder1965"
  REAL(dp)            :: uniform_ice_heat_capacity_config            = 2009._dp                         ! Uniform ice heat capacity (applied when choice_ice_heat_capacity_config = "uniform")
  CHARACTER(LEN=256)  :: choice_ice_thermal_conductivity_config      = 'Ritz1987'                       ! Choice of ice heat capacity model: "uniform", "Ritz1987"
  REAL(dp)            :: uniform_ice_thermal_conductivity_config     = 6.626958E7_dp                    ! Uniform ice thermal conductivity (applied when choice_ice_thermal_conductivity_config = "uniform")

  ! Climate
  ! =======

  CHARACTER(LEN=256)  :: choice_climate_model_config                 = 'matrix_warm_cold'               ! Choice of climate model: "none", "idealised", "PD_obs", "PD_dTglob", "matrix_warm_cold", "direct_global", "direct_regional"
  CHARACTER(LEN=256)  :: choice_idealised_climate_config             = 'EISMINT1_A'

  ! Folder with NetCDF files containing direct climate forcing
  CHARACTER(LEN=256)  :: direct_climate_foldername_NAM_config        = ''
  CHARACTER(LEN=256)  :: direct_climate_foldername_EAS_config        = ''
  CHARACTER(LEN=256)  :: direct_climate_foldername_GRL_config        = ''
  CHARACTER(LEN=256)  :: direct_climate_foldername_ANT_config        = ''

  CHARACTER(LEN=256)  :: direct_climate_basefilename_NAM_config      = ''
  CHARACTER(LEN=256)  :: direct_climate_basefilename_EAS_config      = ''
  CHARACTER(LEN=256)  :: direct_climate_basefilename_GRL_config      = ''
  CHARACTER(LEN=256)  :: direct_climate_basefilename_ANT_config      = ''

  ! Whether or not to apply geometry corrections to direct temperature/precipitation
  LOGICAL             :: do_direct_climate_geo_corr_config           = .TRUE.

  ! NetCDF file containing the present-day observed climate (e.g. ERA40)
  CHARACTER(LEN=256)  :: filename_PD_obs_climate_config              = '/Users/berends/Documents/Datasets/ERA40/ERA40_climate_global.nc'

  ! GCM snapshots in the matrix_warm_cold option
  CHARACTER(LEN=256)  :: filename_climate_snapshot_PI_config         = '/Users/berends/Documents/Datasets/GCM_snapshots/Singarayer_Valdes_2010_PI_Control.nc'
  CHARACTER(LEN=256)  :: filename_climate_snapshot_warm_config       = '/Users/berends/Documents/Datasets/GCM_snapshots/Singarayer_Valdes_2010_PI_Control.nc'
  CHARACTER(LEN=256)  :: filename_climate_snapshot_cold_config       = '/Users/berends/Documents/Datasets/GCM_snapshots/Singarayer_Valdes_2010_LGM.nc'

  ! Ice and ocean mask from GCM snapshots
  CHARACTER(LEN=256)  :: reference_mask_method_config                = 'estimate'                        ! 'estimate' (based on GCM topography and temperature) or 'file' from a selected file
  
  CHARACTER(LEN=256)  :: filename_snapshot_mask_PD_obs_config        = ''   
  CHARACTER(LEN=256)  :: filename_snapshot_mask_GCM_PI_config        = ''
  CHARACTER(LEN=256)  :: filename_snapshot_mask_GCM_warm_config      = ''
  CHARACTER(LEN=256)  :: filename_snapshot_mask_GCM_cold_config      = ''
  
  ! Lapse rate
  REAL(dp)            :: constant_lapserate_config                   = 0.008_dp                         ! Constant atmospheric lapse rate [K m^-1]

  ! Scaling factor for CO2 vs ice weights
  REAL(dp)            :: climate_matrix_CO2vsice_NAM_config          = 0.5_dp                           ! Weight factor for the influence of CO2 vs ice cover on temperature
  REAL(dp)            :: climate_matrix_CO2vsice_EAS_config          = 0.5_dp                           ! Can be set separately for different regions
  REAL(dp)            :: climate_matrix_CO2vsice_GRL_config          = 0.75_dp                          ! Default values are from Berends et al, 2018
  REAL(dp)            :: climate_matrix_CO2vsice_ANT_config          = 0.75_dp                          ! 1.0_dp equals glacial index method

  ! Orbit time and CO2 concentration of the warm and cold snapshots
  REAL(dp)            :: matrix_high_CO2_level_config                = 280._dp                          ! CO2 level  pertaining to the warm climate (PI  level default)
  REAL(dp)            :: matrix_low_CO2_level_config                 = 190._dp                          ! CO2 level  pertaining to the cold climate (LGM level default)
  REAL(dp)            :: matrix_warm_orbit_time_config               = 0._dp                            ! Orbit time pertaining to the warm climate (PI default)
  REAL(dp)            :: matrix_cold_orbit_time_config               = -21000._dp                       ! Orbit time pertaining to the cold climate (LGM default)

  ! Whether or not to apply a bias correction to the GCM snapshots
  LOGICAL             :: climate_matrix_biascorrect_warm_config      = .TRUE.                           ! Whether or not to apply a bias correction (modelled vs observed PI climate) to the "warm" GCM snapshot
  LOGICAL             :: climate_matrix_biascorrect_cold_config      = .FALSE.                          ! Whether or not to apply a bias correction (modelled vs observed PI climate) to the "cold" GCM snapshot

  LOGICAL             :: switch_glacial_index_precip_config          = .FALSE.                          ! If a glacial index is used for the precipitation forcing, it will only depend on CO2

  ! Ocean
  ! =====

  CHARACTER(LEN=256)  :: choice_ocean_model_config                   = 'matrix_warm_cold'               ! Choice of ocean model: "none", "idealised", "uniform_warm_cold", "PD_obs", "matrix_warm_cold"
  CHARACTER(LEN=256)  :: choice_idealised_ocean_config               = 'MISMIP+_warm'                   ! Choice of idealised ocean: 'MISMIP+_warm', 'MISMIP+_cold', 'MISOMIP1', 'Reese2018_ANT'

  ! Delta ocean temperature inversion
  LOGICAL             :: do_ocean_temperature_inversion_config       = .FALSE.                          ! Whether or not to apply the inversion
  REAL(dp)            :: ocean_temperature_inv_t_start_config        = -9.9E9_dp                        ! Minimum model time when the inversion is allowed
  REAL(dp)            :: ocean_temperature_inv_t_end_config          = +9.9E9_dp                        ! Maximum model time when the inversion is allowed
  REAL(dp)            :: ocean_temperature_inv_dT0_config            = 1.0_dp                           ! Adjustment scale for regularisation in Pien's delta inversion method
  REAL(dp)            :: ocean_temperature_inv_H0_config             = 200.0_dp                         ! Ice thickness scale for regularisation in Pien's delta inversion method
  REAL(dp)            :: ocean_temperature_inv_tau_config            = 200.0_dp                         ! Time scale for regularisation in Pien's delta inversion method
  REAL(dp)            :: ocean_temperature_inv_min_config            = -3.0_dp                          ! Lower delta temperature limit in Pien's delta inversion method
  REAL(dp)            :: ocean_temperature_inv_max_config            = +3.0_dp                          ! Upper delta temperature limit in Pien's delta inversion method
  CHARACTER(LEN=256)  :: inverted_ocean_filename_output_config       = 'ocean_inv.nc'                   ! NetCDF file where the final inverted ocean state will be saved

  ! NetCDF file containing the present-day observed ocean (WOA18) (NetCDF)
  CHARACTER(LEN=256)  :: filename_PD_obs_ocean_config                = '/Users/berends/Documents/Datasets/WOA/woa18_decav_ts00_04_remapcon_r360x180_NaN.nc'
  CHARACTER(LEN=256)  :: name_ocean_temperature_obs_config           = 't_an' ! E.g. objectively analysed mean (t_an) or statistical mean (t_mn)
  CHARACTER(LEN=256)  :: name_ocean_salinity_obs_config              = 's_an' ! E.g. objectively analysed mean (s_an) or statistical mean (s_mn)

  ! GCM snapshots in the matrix_warm_cold option
  CHARACTER(LEN=256)  :: filename_GCM_ocean_snapshot_PI_config       = '/Users/berends/Documents/Datasets/COSMOS_ocean_examples/COSMOS_PI_oceanTS_prep.nc'
  CHARACTER(LEN=256)  :: filename_GCM_ocean_snapshot_warm_config     = '/Users/berends/Documents/Datasets/COSMOS_ocean_examples/COSMOS_PI_oceanTS_prep.nc'
  CHARACTER(LEN=256)  :: filename_GCM_ocean_snapshot_cold_config     = '/Users/berends/Documents/Datasets/COSMOS_ocean_examples/COSMOS_LGM_oceanTS_prep.nc'
  CHARACTER(LEN=256)  :: name_ocean_temperature_GCM_config           = 't_an'
  CHARACTER(LEN=256)  :: name_ocean_salinity_GCM_config              = 's_an'

  ! Uniform ocean temperature values used when choice_ocean_model = "uniform_warm_cold"
  REAL(dp)            :: ocean_temperature_PD_config                 = 271.46_dp                        ! present day temperature of the ocean beneath the shelves [K; -1.7 Celsius]
  REAL(dp)            :: ocean_temperature_cold_config               = 268.16_dp                        ! cold period temperature of the ocean beneath the shelves [K; -5.0 Celcius]
  REAL(dp)            :: ocean_temperature_warm_config               = 275.16_dp                        ! warm period temperature of the ocean beneath the shelves [K;  2.0 Celcius]

  ! Parameters used when choice_idealised_ocean = "matrix_warm_cold"
  CHARACTER(LEN=256)  :: choice_ocean_vertical_grid_config           = 'regular'                        ! Choice of vertical grid to be used for ocean data
  REAL(dp)            :: ocean_vertical_grid_max_depth_config        = 1500._dp                         ! Maximum depth           to be used for ocean data
  REAL(dp)            :: ocean_regular_grid_dz_config                = 150._dp                          ! Vertical grid spacing   to be used for ocean data when choice_ocean_vertical_grid_config = 'regular'
  LOGICAL             :: do_ocean_extrap_config                      = .FALSE.                          ! Create extrapolated ocean file
  CHARACTER(LEN=256)  :: ocean_extrap_dir_config                     = 'extrapolated_ocean_files'       ! Directory where extrapolated ocean files are stored
  REAL(dp)            :: ocean_extrap_res_config                     = 5000._dp                         ! High resolution at which the ocean data extrapolation should be performed
  REAL(dp)            :: ocean_extrap_Gauss_sigma_config             = 8000._dp                         ! 1-sigma of the Gaussian smoothing operation used to extrapolate the ocean data
  CHARACTER(LEN=256)  :: ocean_extrap_hires_geo_filename_NAM_config  = ''                               ! Path to a NetCDF file containing (present-day) geometry at high resolution, used for ocean data extrapolation
  CHARACTER(LEN=256)  :: ocean_extrap_hires_geo_filename_EAS_config  = ''
  CHARACTER(LEN=256)  :: ocean_extrap_hires_geo_filename_GRL_config  = ''
  CHARACTER(LEN=256)  :: ocean_extrap_hires_geo_filename_ANT_config  = ''
  REAL(dp)            :: ocean_w_tot_hist_averaging_window_config    = 1500._dp                         ! Time window (in yr) over which the weighing fields for sea-water temperature at maximum depth are averaged

  ! Scaling factor for CO2 vs ice weights
  REAL(dp)            :: ocean_matrix_CO2vsice_NAM_config            = 0.5_dp                           ! Weight factor for the influence of CO2 vs ice cover on ocean T and S
  REAL(dp)            :: ocean_matrix_CO2vsice_EAS_config            = 0.5_dp                           ! Can be set separately for different regions
  REAL(dp)            :: ocean_matrix_CO2vsice_GRL_config            = 0.75_dp
  REAL(dp)            :: ocean_matrix_CO2vsice_ANT_config            = 0.75_dp

  ! Basin-dependent linear temperature profiles (used when choice_idealised_ocean = "linear_per_basin")
  REAL(dp), DIMENSION(100) :: ocean_T_surf_per_basin_config = &
      (/ -1.69_dp, -1.59_dp, -1.51_dp, -1.41_dp, -1.30_dp, -1.38_dp, -1.13_dp, -1.50_dp, -1.40_dp, -1.10_dp, &
         -1.58_dp, -1.16_dp, -1.27_dp, -0.80_dp, -0.80_dp, -1.30_dp, -1.37_dp,  0._dp  ,  0._dp  ,  0._dp  , &
          0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  , &
          0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  , &
          0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  , &
          0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  , &
          0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  , &
          0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  , &
          0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  , &
          0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  ,  0._dp  /)
  REAL(dp), DIMENSION(100) :: ocean_dT_dz_per_basin_config = &
      (/ -0.000419_dp, -0.000293_dp, -0.000602_dp, -0.001165_dp, -0.000600_dp, -0.000618_dp, -0.002153_dp, -0.000400_dp, -0.000400_dp, -0.000400_dp, &
         -0.000487_dp, -0.001189_dp, -0.000238_dp,  0.004000_dp, -0.000400_dp, -0.000400_dp, -0.001392_dp,  0._dp      ,  0._dp      ,  0._dp      , &
          0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      , &
          0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      , &
          0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      , &
          0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      , &
          0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      , &
          0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      , &
          0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      , &
          0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      ,  0._dp      /)
  LOGICAL             :: do_invert_linear_per_basin_config           = .FALSE.                          ! Whether or not to invert for ocean temperature based on thinning rates

  ! Apply anomalies to baseline (used when choice_idealised_ocean = "anomalies")
  CHARACTER(LEN=256)  :: ocean_filename_baseline_config        = ''                              ! NetCDF file containing the baseline ocean for the ISMIP-style ocean
  CHARACTER(LEN=256)  :: ocean_foldername_aTO_config           = ''                              ! Folder containing the single timeframe NetCDF files of the ocean temperature anomaly
  CHARACTER(LEN=256)  :: ocean_basefilename_aTO_config         = ''                              ! Filename without the year (e.g. if the actual file is "aTO-1950.nc",   then this variable should be "aTO-"
  CHARACTER(LEN=256)  :: ocean_foldername_aSO_config           = ''                              ! Folder containing the single timeframe NetCDF files of the ocean salinity anomaly
  CHARACTER(LEN=256)  :: ocean_basefilename_aSO_config         = ''                              ! Filename without the year (e.g. if the actual file is "aSO-1950.nc", then this variable should be "aSO-"


  ! Surface mass balance
  ! ====================

  CHARACTER(LEN=256)  :: choice_SMB_model_config                     = 'IMAU-ITM'                       ! Choice of SMB model: "uniform", "idealised", "IMAU-ITM", "direct_global", "direct_regional", "snapshot"
  CHARACTER(LEN=256)  :: choice_idealised_SMB_config                 = 'EISMINT1_A'
  REAL(dp)            :: SMB_uniform_config                          = 0._dp                            ! Uniform SMB, applied when choice_SMB_model = "uniform" [mie/yr]

  ! NetCDF file containing direct global/regional climate forcing
  CHARACTER(LEN=256)  :: filename_direct_global_SMB_config           = ''
  CHARACTER(LEN=256)  :: filename_direct_regional_SMB_NAM_config     = ''
  CHARACTER(LEN=256)  :: filename_direct_regional_SMB_EAS_config     = ''
  CHARACTER(LEN=256)  :: filename_direct_regional_SMB_GRL_config     = ''
  CHARACTER(LEN=256)  :: filename_direct_regional_SMB_ANT_config     = ''

  ! NetCDF file containing a direct time-less SMB snapshot
  CHARACTER(LEN=256)  :: filename_SMB_snapshot_NAM_config            = ''
  CHARACTER(LEN=256)  :: filename_SMB_snapshot_EAS_config            = ''
  CHARACTER(LEN=256)  :: filename_SMB_snapshot_GRL_config            = ''
  CHARACTER(LEN=256)  :: filename_SMB_snapshot_ANT_config            = ''

  ! Tuning parameters for the IMAU-ITM SMB model
  CHARACTER(LEN=256)  :: SMB_IMAUITM_choice_init_firn_NAM_config     = 'uniform'                        ! How to initialise the firn layer in the IMAU-ITM SMB model: "uniform", "restart"
  CHARACTER(LEN=256)  :: SMB_IMAUITM_choice_init_firn_EAS_config     = 'uniform'
  CHARACTER(LEN=256)  :: SMB_IMAUITM_choice_init_firn_GRL_config     = 'uniform'
  CHARACTER(LEN=256)  :: SMB_IMAUITM_choice_init_firn_ANT_config     = 'uniform'
  REAL(dp)            :: SMB_IMAUITM_initial_firn_thickness_config   = 1._dp                            ! Initial firn thickness of the IMAU-ITEM SMB model [m] (used when SMB_IMAUITM_choice_init_firn = "uniform")
  REAL(dp)            :: SMB_IMAUITM_C_abl_constant_NAM_config       = -49._dp                          ! 34._dp    (commented values are old ANICE defaults, but since refreezing was not calculated right
  REAL(dp)            :: SMB_IMAUITM_C_abl_constant_EAS_config       = -49._dp                          !            and this has since been fixed, these values will still not give the same results as
  REAL(dp)            :: SMB_IMAUITM_C_abl_constant_GRL_config       = -49._dp                          !            they used to in ANICE.)
  REAL(dp)            :: SMB_IMAUITM_C_abl_constant_ANT_config       = -49._dp
  REAL(dp)            :: SMB_IMAUITM_C_abl_Ts_NAM_config             = 10._dp                           ! 10._dp
  REAL(dp)            :: SMB_IMAUITM_C_abl_Ts_EAS_config             = 10._dp
  REAL(dp)            :: SMB_IMAUITM_C_abl_Ts_GRL_config             = 10._dp
  REAL(dp)            :: SMB_IMAUITM_C_abl_Ts_ANT_config             = 10._dp
  REAL(dp)            :: SMB_IMAUITM_C_abl_Q_NAM_config              = 0.0227_dp                        ! 0.513_dp
  REAL(dp)            :: SMB_IMAUITM_C_abl_Q_EAS_config              = 0.0227_dp
  REAL(dp)            :: SMB_IMAUITM_C_abl_Q_GRL_config              = 0.0227_dp
  REAL(dp)            :: SMB_IMAUITM_C_abl_Q_ANT_config              = 0.0227_dp
  REAL(dp)            :: SMB_IMAUITM_C_refr_NAM_config               = 0.051_dp                         ! 0.012_dp
  REAL(dp)            :: SMB_IMAUITM_C_refr_EAS_config               = 0.051_dp
  REAL(dp)            :: SMB_IMAUITM_C_refr_GRL_config               = 0.051_dp
  REAL(dp)            :: SMB_IMAUITM_C_refr_ANT_config               = 0.051_dp

  ! ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing
  ! ==============================================================

  CHARACTER(LEN=256)  :: ISMIP_forcing_filename_baseline_config      = ''                              ! NetCDF file containing the baseline climate for the ISMIP-style forcing
  CHARACTER(LEN=256)  :: ISMIP_forcing_foldername_aSMB_config        = ''                              ! Folder containing the single-year NetCDF files of the SMB anomaly
  CHARACTER(LEN=256)  :: ISMIP_forcing_basefilename_aSMB_config      = ''                              ! Filename without the year (e.g. if the actual file is "aSMB_MARv3.12-yearly-CESM2-ssp585-1950.nc",   then this variable should be "aSMB_MARv3.12-yearly-CESM2-ssp585-"
  CHARACTER(LEN=256)  :: ISMIP_forcing_foldername_dSMBdz_config      = ''                              ! Folder containing the single-year NetCDF files of the SMB lapse rate
  CHARACTER(LEN=256)  :: ISMIP_forcing_basefilename_dSMBdz_config    = ''                              ! Filename without the year (e.g. if the actual file is "dSMBdz_MARv3.12-yearly-CESM2-ssp585-1950.nc", then this variable should be "dSMBdz_MARv3.12-yearly-CESM2-ssp585-"
  CHARACTER(LEN=256)  :: ISMIP_forcing_foldername_aST_config         = ''                              ! Folder containing the single-year NetCDF files of the temperature anomaly
  CHARACTER(LEN=256)  :: ISMIP_forcing_basefilename_aST_config       = ''                              ! Filename without the year (e.g. if the actual file is "aST_MARv3.12-yearly-CESM2-ssp585-1950.nc",    then this variable should be "aST_MARv3.12-yearly-CESM2-ssp585-"
  CHARACTER(LEN=256)  :: ISMIP_forcing_foldername_dSTdz_config       = ''                              ! Folder containing the single-year NetCDF files of the temperature lapse rate
  CHARACTER(LEN=256)  :: ISMIP_forcing_basefilename_dSTdz_config     = ''                              ! Filename without the year (e.g. if the actual file is "dSTdz_MARv3.12-yearly-CESM2-ssp585-1950.nc",  then this variable should be "dSTdz_MARv3.12-yearly-CESM2-ssp585-"

  ! Basal mass balance
  ! ==================

  CHARACTER(LEN=256)  :: choice_BMB_shelf_model_config               = 'ANICE_legacy'                   ! Choice of shelf BMB: "uniform", "idealised", "ANICE_legacy", "Favier2019_lin", "Favier2019_quad", "Favier2019_Mplus", "Lazeroms2018_plume", "PICO", "PICOP"
  CHARACTER(LEN=256)  :: choice_idealised_BMB_shelf_config           = 'MISMIP+'
  CHARACTER(LEN=256)  :: choice_BMB_sheet_model_config               = 'uniform'                        ! Choice of sheet BMB: "uniform"
  REAL(dp)            :: BMB_shelf_uniform_config                    = 0._dp                            ! Uniform shelf BMB, applied when choice_BMB_shelf_model = "uniform" [mie/yr]
  REAL(dp)            :: BMB_sheet_uniform_config                    = 0._dp                            ! Uniform sheet BMB, applied when choice_BMB_sheet_model = "uniform" [mie/yr]
  CHARACTER(LEN=256)  :: choice_BMB_subgrid_config                   = 'FCMP'                           ! Choice of sub-grid BMB scheme: "FCMP", "PMP", "NMP" (following Leguy et al., 2021)
  LOGICAL             :: do_asynchronous_BMB_config                  = .FALSE.                          ! Whether or not to run the BMB asynchronously from the ice dynamics (if so, run it at dt_BMB; if not, run it in every ice dynamics time step)
  REAL(dp)            :: BMB_max_config                              = 50._dp                           ! Maximum amount of allowed basal melt     [mie/yr]
  REAL(dp)            :: BMB_min_config                              = 5._dp                            ! Maximum amount of allowed basal freezing [mie/yr]

  CHARACTER(LEN=256)  :: choice_basin_scheme_NAM_config              = 'none'                           ! Choice of basin ID scheme; can be 'none' or 'file'
  CHARACTER(LEN=256)  :: choice_basin_scheme_EAS_config              = 'none'
  CHARACTER(LEN=256)  :: choice_basin_scheme_GRL_config              = 'none'
  CHARACTER(LEN=256)  :: choice_basin_scheme_ANT_config              = 'none'
  CHARACTER(LEN=256)  :: filename_basins_NAM_config                  = ''                               ! Path to a text file containing polygons of drainage basins
  CHARACTER(LEN=256)  :: filename_basins_EAS_config                  = ''
  CHARACTER(LEN=256)  :: filename_basins_GRL_config                  = ''
  CHARACTER(LEN=256)  :: filename_basins_ANT_config                  = ''
  LOGICAL             :: do_merge_basins_ANT_config                  = .TRUE.                           ! Whether or not to merge some of the Antarctic basins
  LOGICAL             :: do_merge_basins_GRL_config                  = .TRUE.                           ! Whether or not to merge some of the Greenland basins

  CHARACTER(LEN=256)       ::  choice_BMB_shelf_amplification_config        = 'basin'                   ! Choice of method to determine BMB amplification factors: "uniform", "basin"
  INTEGER                  ::  basin_BMB_amplification_n_ANT_config         = 17                        ! Number of basins used for ANT
  REAL(dp), DIMENSION(17)  ::  basin_BMB_amplification_factor_ANT_config    = &                         ! BMB amplification factor for each basin for ANT
    (/ 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, &
       1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp /)
  INTEGER                  ::  basin_BMB_amplification_n_GRL_config         = 8                         ! Number of basins used for GRL
  REAL(dp), DIMENSION(8)   ::  basin_BMB_amplification_factor_GRL_config    = &                         ! BMB amplification factor for each basin for GRL
    (/ 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp /)

  ! Parameters for the three simple melt parameterisations from Favier et al. (2019)
  REAL(dp)            :: BMB_Favier2019_lin_GammaT_config            = 3.3314E-05_dp  ! 2.03E-5_dp      ! Heat exchange velocity [m s^-1]
  REAL(dp)            :: BMB_Favier2019_quad_GammaT_config           = 111.6E-5_dp    ! 99.32E-5_dp     ! Commented values are from Favier et al. (2019), Table 3
  REAL(dp)            :: BMB_Favier2019_Mplus_GammaT_config          = 108.6E-5_dp    ! 132.9E-5_dp     ! Actual value are re-tuned for IMAU-ICE, following the same approach (see Asay-Davis et al., 2016, ISOMIP+)

  ! Parameters for the Lazeroms et al. (2018) plume-parameterisation BMB model
  REAL(dp)            :: BMB_Lazeroms2018_GammaT_config              = 3.7506E-04_dp  ! 1.1E-3_dp       ! Thermal exchange velocity; tuned following ISOMIP+ protocol (Asay-Davis et al., 2016, Sect. 3.2.1), commented value from Lazeroms et al. (2018)
  CHARACTER(LEN=256)  :: BMB_Lazeroms2018_find_GL_scheme_config      = 'GL_average_Tijn'                 ! How to determine the GL origin of a plume: "GL_average", "along_ice_flow"

  ! Parameters for the PICO BMB model
  INTEGER             :: BMB_PICO_nboxes_config                      = 5                                ! Number of sub-shelf ocean boxes used by PICO
  REAL(dp)            :: BMB_PICO_GammaTstar_config                  = 3.6131E-05_dp  ! 2.0E-5_dp       ! Effective turbulent temperature exchange velocity [m s^-1]; tuned following ISOMIP+ protocol (Asay-Davis et al., 2016, Sect. 3.2.1), commented value from Reese et al. (2018)

  ! File path for the LADDIE model output
  CHARACTER(LEN=256)  :: filename_BMB_LADDIE_config                  = ''                               ! Path to a netcdf file containing melt pattern computed by LADDIE

  ! Parameters for the ANICE_legacy sub-shelf melt model
  REAL(dp)            :: T_ocean_mean_PD_NAM_config                  = -1.7_dp                          ! Present day temperature of the ocean beneath the shelves [Celcius]
  REAL(dp)            :: T_ocean_mean_PD_EAS_config                  = -1.7_dp
  REAL(dp)            :: T_ocean_mean_PD_GRL_config                  =  2.0_dp
  REAL(dp)            :: T_ocean_mean_PD_ANT_config                  = -1.7_dp
  REAL(dp)            :: T_ocean_mean_cold_NAM_config                = -5.0_dp                          ! Cold period temperature of the ocean beneath the shelves [Celcius]
  REAL(dp)            :: T_ocean_mean_cold_EAS_config                = -5.0_dp
  REAL(dp)            :: T_ocean_mean_cold_GRL_config                =  0.0_dp
  REAL(dp)            :: T_ocean_mean_cold_ANT_config                = -5.0_dp
  REAL(dp)            :: T_ocean_mean_warm_NAM_config                =  2.0_dp                          ! Warm period temperature of the ocean beneath the shelves [Celcius]
  REAL(dp)            :: T_ocean_mean_warm_EAS_config                =  2.0_dp
  REAL(dp)            :: T_ocean_mean_warm_GRL_config                =  4.0_dp
  REAL(dp)            :: T_ocean_mean_warm_ANT_config                =  2.0_dp

  REAL(dp)            :: BMB_deepocean_PD_NAM_config                 =  -5._dp                          ! Present-day sub-shelf melt rate for deep-ocean areas [m/year]
  REAL(dp)            :: BMB_deepocean_PD_EAS_config                 =  -5._dp
  REAL(dp)            :: BMB_deepocean_PD_GRL_config                 =  -5._dp
  REAL(dp)            :: BMB_deepocean_PD_ANT_config                 =  -5._dp
  REAL(dp)            :: BMB_deepocean_cold_NAM_config               =  -2._dp                          ! Cold period sub-shelf melt rate for deep-ocean areas [m/year]
  REAL(dp)            :: BMB_deepocean_cold_EAS_config               =  -2._dp
  REAL(dp)            :: BMB_deepocean_cold_GRL_config               =  -2._dp
  REAL(dp)            :: BMB_deepocean_cold_ANT_config               =  -2._dp
  REAL(dp)            :: BMB_deepocean_warm_NAM_config               = -10._dp                          ! Warm period sub-shelf melt rate for deep-ocean areas [m/year]
  REAL(dp)            :: BMB_deepocean_warm_EAS_config               = -10._dp
  REAL(dp)            :: BMB_deepocean_warm_GRL_config               = -10._dp
  REAL(dp)            :: BMB_deepocean_warm_ANT_config               = -10._dp

  REAL(dp)            :: BMB_shelf_exposed_PD_NAM_config             =  -3._dp                          ! Present-day sub-shelf melt rate for exposed areas    [m/year]
  REAL(dp)            :: BMB_shelf_exposed_PD_EAS_config             =  -3._dp
  REAL(dp)            :: BMB_shelf_exposed_PD_GRL_config             =  -3._dp
  REAL(dp)            :: BMB_shelf_exposed_PD_ANT_config             =  -3._dp
  REAL(dp)            :: BMB_shelf_exposed_cold_NAM_config           =  -0._dp                          ! Cold period sub-shelf melt rate for exposed areas    [m/year]
  REAL(dp)            :: BMB_shelf_exposed_cold_EAS_config           =  -0._dp
  REAL(dp)            :: BMB_shelf_exposed_cold_GRL_config           =  -0._dp
  REAL(dp)            :: BMB_shelf_exposed_cold_ANT_config           =  -0._dp
  REAL(dp)            :: BMB_shelf_exposed_warm_NAM_config           =  -6._dp                          ! Warm period sub-shelf melt rate for exposed areas    [m/year]
  REAL(dp)            :: BMB_shelf_exposed_warm_EAS_config           =  -6._dp
  REAL(dp)            :: BMB_shelf_exposed_warm_GRL_config           =  -6._dp
  REAL(dp)            :: BMB_shelf_exposed_warm_ANT_config           =  -6._dp

  REAL(dp)            :: subshelf_melt_factor_NAM_config             = 0.005_dp                         ! Overall tuning factor for sub-shelf melt rate
  REAL(dp)            :: subshelf_melt_factor_EAS_config             = 0.005_dp
  REAL(dp)            :: subshelf_melt_factor_GRL_config             = 0.005_dp
  REAL(dp)            :: subshelf_melt_factor_ANT_config             = 0.005_dp

  REAL(dp)            :: deep_ocean_threshold_depth_NAM_config       = 1200._dp                         ! Threshold water depth for "deep ocean" (as opposed to continental shelf);
  REAL(dp)            :: deep_ocean_threshold_depth_EAS_config       = 800._dp                          ! this mostly prevents ice shelves from growing beyond the continental shelf
  REAL(dp)            :: deep_ocean_threshold_depth_GRL_config       = 800._dp                          ! Different depths for different regions is a bit ad hoc, but in reality
  REAL(dp)            :: deep_ocean_threshold_depth_ANT_config       = 1800._dp                         ! the different surface ocean temperatures probably result in the same effect...

  ! Englacial isotope tracing
  ! =========================

  CHARACTER(LEN=256)  :: choice_ice_isotopes_model_config             = 'ANICE_legacy'                  ! Choice of englacial isotopes model: "none", "uniform", "ANICE_legacy"
  REAL(dp)            :: uniform_ice_d18O_config                      = 0._dp                           ! Uniform englacial d18O-value (used when choice_ice_isotope_model_config = "uniform")

  ! Sea level and GIA
  ! =================

  LOGICAL             :: do_ocean_floodfill_config                   = .TRUE.                           ! Use a flood-fill to determine the ocean mask, so that (pro-/sub-glacial) lakes dont exist
  CHARACTER(LEN=256)  :: choice_sealevel_model_config                = 'eustatic'                       ! Can be "fixed", "prescribed", "eustatic", or "SELEN"
  REAL(dp)            :: fixed_sealevel_config                       = 0._dp                            ! Height of fixed sealevel w.r.t. PD
  CHARACTER(LEN=256)  :: filename_sealevel_record_config             = 'name_of_file.dat'               ! Name of a file containing a sealevel record (in years)
  INTEGER             :: sealevel_record_length_config               = 1

  CHARACTER(LEN=256)  :: choice_GIA_model_config                     = 'ELRA'                           ! Can be "none", "ELRA", or "SELEN"
  REAL(dp)            :: dx_GIA_config                               = 100000._dp                       ! Horizontal resolution of the square grid used for the GIA model
  REAL(dp)            :: ELRA_lithosphere_flex_rigidity_config       = 1.0E+25_dp                       ! Lithospheric flexural rigidity [kg m^2 s^-2]
  REAL(dp)            :: ELRA_bedrock_relaxation_time_config         = 3000.0_dp                        ! Relaxation time for bedrock adjustment [yr]
  REAL(dp)            :: ELRA_mantle_density_config                  = 3300.0_dp                        ! Mantle density [kg m^-3]

  ! SELEN
  ! =====

  LOGICAL             :: SELEN_run_at_t_start_config                  = .FALSE.                         ! Whether or not to run SELEN in the first coupling loop (needed for some benchmark experiments)
  INTEGER             :: SELEN_n_TDOF_iterations_config               = 1                               ! Number of Time-Dependent Ocean Function iterations
  INTEGER             :: SELEN_n_recursion_iterations_config          = 1                               ! Number of recursion iterations
  LOGICAL             :: SELEN_use_rotational_feedback_config         = .FALSE.                         ! If TRUE, rotational feedback is included
  INTEGER             :: SELEN_n_harmonics_config                     = 128                             ! Maximum number of harmonic degrees
  LOGICAL             :: SELEN_display_progress_config                = .FALSE.                         ! Whether or not to display the progress of the big loops to the screen (doesn't work on Cartesius!)

  CHARACTER(LEN=256)  :: SELEN_dir_config                             = 'SELEN_files'                   ! Directory where SELEN initial files and spherical harmonics are stored
  CHARACTER(LEN=256)  :: SELEN_global_topo_filename_config            = 'SELEN_global_topography.nc'    ! Filename for the SELEN global topography file (located in SELEN_dir)
  CHARACTER(LEN=256)  :: SELEN_TABOO_init_filename_config             = 'SELEN_TABOO_initial_file.dat'  ! Filename for the TABOO initial file           (idem                )
  CHARACTER(LEN=256)  :: SELEN_LMJ_VALUES_filename_config             = 'SELEN_lmj_values.bin'          ! Filename for the LJ and MJ values file        (idem                )

  INTEGER                  :: SELEN_irreg_time_n_config               = 15                              ! Number of entries in the irregular moving time window
  REAL(dp), DIMENSION(50)  :: SELEN_irreg_time_window_config          = &                               ! Values of entries in the irregular moving time window
   (/20._dp, 20._dp, 20._dp, 5._dp, 5._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp  /)

  REAL(dp)            :: SELEN_lith_thickness_config                  = 100._dp                         ! Thickness of the elastic lithosphere [km]
  INTEGER             :: SELEN_visc_n_config                          = 3                               ! Number      of viscous asthenosphere layers
  REAL(dp), DIMENSION(3) :: SELEN_visc_prof_config                    = (/ 3._dp, 0.6_dp, 0.3_dp /)     ! Viscosities of viscous asthenosphere layers [?]

  ! Settings for the TABOO Earth deformation model
  INTEGER             :: SELEN_TABOO_CDE_config                       = 0                               ! code of the model (see taboo for explanation)
  INTEGER             :: SELEN_TABOO_TLOVE_config                     = 1                               ! Tidal love numbers yes/no
  INTEGER             :: SELEN_TABOO_DEG1_config                      = 1                               ! Tidal love numbers degree
  REAL(dp)            :: SELEN_TABOO_RCMB_config                      = 3480._dp                        ! Radius of CMB (km)

  ! Which data fields will be written to the help_fields output file
  ! ================================================================

  CHARACTER(LEN=256)  :: help_field_01_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_02_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_03_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_04_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_05_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_06_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_07_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_08_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_09_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_10_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_11_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_12_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_13_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_14_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_15_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_16_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_17_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_18_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_19_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_20_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_21_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_22_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_23_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_24_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_25_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_26_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_27_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_28_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_29_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_30_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_31_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_32_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_33_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_34_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_35_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_36_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_37_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_38_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_39_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_40_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_41_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_42_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_43_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_44_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_45_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_46_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_47_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_48_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_49_config                        = 'none'
  CHARACTER(LEN=256)  :: help_field_50_config                        = 'none'


  ! ==========================================================================
  ! The "C" type, which contains all the config parameters as fields.
  ! These will all be overwritten with the values of the "_config" variables,
  ! which are either the default values specified above, are the values
  ! specified from the external config file.
  ! ==========================================================================

  TYPE constants_type

    ! Time steps and range
    ! =====================

    REAL(dp)                            :: start_time_of_run
    REAL(dp)                            :: end_time_of_run
    REAL(dp)                            :: dt_coupling
    REAL(dp)                            :: dt_max
    REAL(dp)                            :: dt_min
    REAL(dp)                            :: dt_startup_phase
    REAL(dp)                            :: dt_thermo
    REAL(dp)                            :: dt_climate
    REAL(dp)                            :: dt_ocean
    REAL(dp)                            :: dt_SMB
    REAL(dp)                            :: dt_BMB
    REAL(dp)                            :: dt_bedrock_ELRA
    REAL(dp)                            :: dt_SELEN
    REAL(dp)                            :: dt_output
    REAL(dp)                            :: dt_output_restart

    ! Which ice sheets do we simulate?
    ! ================================

    LOGICAL                             :: do_NAM
    LOGICAL                             :: do_EAS
    LOGICAL                             :: do_GRL
    LOGICAL                             :: do_ANT

    ! Benchmark experiments
    ! =====================

    ! SSA_icestream (see Schoof 2006, and also Bueler and Brown 2009)
    REAL(dp)                            :: SSA_icestream_A
    REAL(dp)                            :: SSA_icestream_L
    REAL(dp)                            :: SSA_icestream_m
    REAL(dp)                            :: SSA_icestream_tantheta
    REAL(dp)                            :: SSA_icestream_H

    ! ISMIP-HOM (see Pattyn et al. 2008)
    REAL(dp)                            :: ISMIP_HOM_L
    CHARACTER(LEN=256)                  :: ISMIP_HOM_E_Arolla_filename

    ! MISMIP+ (see Asay-Davis et al., 2016)
    LOGICAL                             :: MISMIPplus_do_tune_A_for_GL
    REAL(dp)                            :: MISMIPplus_xGL_target
    REAL(dp)                            :: MISMIPplus_A_flow_initial
    CHARACTER(LEN=256)                  :: MISMIPplus_scenario

    ! MISOMIP1 (see Asay-Davis et al., 2016)
    CHARACTER(LEN=256)                  :: MISOMIP1_scenario

    ! Whether or not to let IMAU_ICE dynamically create its own output folder.
    ! This works fine locally, on LISA its better to use a fixed folder name.
    ! =======================================================================

    LOGICAL                             :: create_procedural_output_dir
    CHARACTER(LEN=256)                  :: fixed_output_dir
    CHARACTER(LEN=256)                  :: fixed_output_dir_suffix
    LOGICAL                             :: do_write_regional_scalar_output
    LOGICAL                             :: do_write_global_scalar_output
    LOGICAL                             :: do_write_regional_scalar_every_timestep

    ! Debugging
    ! =========

    LOGICAL                             :: do_check_for_NaN
    LOGICAL                             :: do_time_display

    ! ISMIP output
    ! =============

    LOGICAL                             :: do_write_ISMIP_output
    CHARACTER(LEN=256)                  :: ISMIP_output_group_code
    CHARACTER(LEN=256)                  :: ISMIP_output_model_code
    CHARACTER(LEN=256)                  :: ISMIP_output_experiment_code
    CHARACTER(LEN=256)                  :: ISMIP_output_basetime

    ! Grids
    ! =====

    ! North America
    REAL(dp)                            :: lambda_M_NAM
    REAL(dp)                            :: phi_M_NAM
    REAL(dp)                            :: beta_stereo_NAM
    REAL(dp)                            :: xmin_NAM
    REAL(dp)                            :: xmax_NAM
    REAL(dp)                            :: ymin_NAM
    REAL(dp)                            :: ymax_NAM
    REAL(dp)                            :: dx_NAM

    ! Eurasia
    REAL(dp)                            :: lambda_M_EAS
    REAL(dp)                            :: phi_M_EAS
    REAL(dp)                            :: beta_stereo_EAS
    REAL(dp)                            :: xmin_EAS
    REAL(dp)                            :: xmax_EAS
    REAL(dp)                            :: ymin_EAS
    REAL(dp)                            :: ymax_EAS
    REAL(dp)                            :: dx_EAS

    ! Greenland
    REAL(dp)                            :: lambda_M_GRL
    REAL(dp)                            :: phi_M_GRL
    REAL(dp)                            :: beta_stereo_GRL
    REAL(dp)                            :: xmin_GRL
    REAL(dp)                            :: xmax_GRL
    REAL(dp)                            :: ymin_GRL
    REAL(dp)                            :: ymax_GRL
    REAL(dp)                            :: dx_GRL

    ! Antarctica
    REAL(dp)                            :: lambda_M_ANT
    REAL(dp)                            :: phi_M_ANT
    REAL(dp)                            :: beta_stereo_ANT
    REAL(dp)                            :: xmin_ANT
    REAL(dp)                            :: xmax_ANT
    REAL(dp)                            :: ymin_ANT
    REAL(dp)                            :: ymax_ANT
    REAL(dp)                            :: dx_ANT

    ! Scaled vertical coordinate zeta
    ! ===============================

    INTEGER                             :: nz
    REAL(dp), DIMENSION(:), ALLOCATABLE :: zeta

    ! Reference geometries (initial, present-day, and GIA equilibrium)
    ! ================================================================

    ! Initial geometry
    CHARACTER(LEN=256)                  :: choice_refgeo_init_NAM
    CHARACTER(LEN=256)                  :: choice_refgeo_init_EAS
    CHARACTER(LEN=256)                  :: choice_refgeo_init_GRL
    CHARACTER(LEN=256)                  :: choice_refgeo_init_ANT
    REAL(dp)                            :: time_to_restart_from_NAM
    REAL(dp)                            :: time_to_restart_from_EAS
    REAL(dp)                            :: time_to_restart_from_GRL
    REAL(dp)                            :: time_to_restart_from_ANT
    CHARACTER(LEN=256)                  :: choice_refgeo_init_idealised
    CHARACTER(LEN=256)                  :: filename_refgeo_init_NAM
    CHARACTER(LEN=256)                  :: filename_refgeo_init_EAS
    CHARACTER(LEN=256)                  :: filename_refgeo_init_GRL
    CHARACTER(LEN=256)                  :: filename_refgeo_init_ANT

    ! Present-day geometry
    CHARACTER(LEN=256)                  :: choice_refgeo_PD_NAM
    CHARACTER(LEN=256)                  :: choice_refgeo_PD_EAS
    CHARACTER(LEN=256)                  :: choice_refgeo_PD_GRL
    CHARACTER(LEN=256)                  :: choice_refgeo_PD_ANT
    CHARACTER(LEN=256)                  :: choice_refgeo_PD_idealised
    CHARACTER(LEN=256)                  :: filename_refgeo_PD_NAM
    CHARACTER(LEN=256)                  :: filename_refgeo_PD_EAS
    CHARACTER(LEN=256)                  :: filename_refgeo_PD_GRL
    CHARACTER(LEN=256)                  :: filename_refgeo_PD_ANT

    ! GIA equilibrium geometry
    CHARACTER(LEN=256)                  :: choice_refgeo_GIAeq_NAM
    CHARACTER(LEN=256)                  :: choice_refgeo_GIAeq_EAS
    CHARACTER(LEN=256)                  :: choice_refgeo_GIAeq_GRL
    CHARACTER(LEN=256)                  :: choice_refgeo_GIAeq_ANT
    CHARACTER(LEN=256)                  :: choice_refgeo_GIAeq_idealised
    CHARACTER(LEN=256)                  :: filename_refgeo_GIAeq_NAM
    CHARACTER(LEN=256)                  :: filename_refgeo_GIAeq_EAS
    CHARACTER(LEN=256)                  :: filename_refgeo_GIAeq_GRL
    CHARACTER(LEN=256)                  :: filename_refgeo_GIAeq_ANT

    LOGICAL                             :: remove_Lake_Vostok

    ! Global forcing (insolation, CO2, d18O, geothermal heat flux)
    ! ============================================================

    CHARACTER(LEN=256)                  :: choice_forcing_method

    ! Insolation forcing (NetCDF)
    CHARACTER(LEN=256)                  :: choice_insolation_forcing
    REAL(dp)                            :: static_insolation_time
    CHARACTER(LEN=256)                  :: filename_insolation

    ! CO2 record (ASCII text file, so the number of rows needs to be specified)
    CHARACTER(LEN=256)                  :: filename_CO2_record
    INTEGER                             :: CO2_record_length

    ! d18O record (ASCII text file, so the number of rows needs to be specified)
    CHARACTER(LEN=256)                  :: filename_d18O_record
    INTEGER                             :: d18O_record_length

    LOGICAL                             :: do_combine_CO2_and_insolation
    REAL(dp)                            :: insolation_weigth_mean_NAM         
    REAL(dp)                            :: insolation_weigth_amplitude_NAM    
    REAL(dp)                            :: insolation_weigth_mean_EAS         
    REAL(dp)                            :: insolation_weigth_amplitude_EAS    
    REAL(dp)                            :: insolation_weigth_mean_GRL         
    REAL(dp)                            :: insolation_weigth_amplitude_GRL    
    REAL(dp)                            :: insolation_weigth_mean_ANT         
    REAL(dp)                            :: insolation_weigth_amplitude_ANT

    ! Climate matrix wind
    LOGICAL                             :: do_climate_matrix_wind

    ! Geothermal heat flux
    CHARACTER(LEN=256)                  :: choice_geothermal_heat_flux
    REAL(dp)                            :: constant_geothermal_heat_flux
    CHARACTER(LEN=256)                  :: filename_geothermal_heat_flux

    ! Parameters for calculating modelled benthic d18O
    LOGICAL                             :: do_calculate_benthic_d18O
    REAL(dp)                            :: dT_deepwater_averaging_window
    REAL(dp)                            :: dT_deepwater_dT_surf_ratio
    REAL(dp)                            :: d18O_dT_deepwater_ratio

    ! Parameters for the inverse routine
    REAL(dp)                            :: dT_glob_inverse_averaging_window
    REAL(dp)                            :: inverse_d18O_to_dT_glob_scaling
    REAL(dp)                            :: CO2_inverse_averaging_window
    REAL(dp)                            :: inverse_d18O_to_CO2_scaling
    REAL(dp)                            :: inverse_d18O_to_CO2_initial_CO2

    ! Files and choices for initializing the inverse d18O routines
    CHARACTER(LEN=256)                  :: choice_d18O_inverse_init_NAM
    CHARACTER(LEN=256)                  :: choice_d18O_inverse_init_EAS
    CHARACTER(LEN=256)                  :: choice_d18O_inverse_init_GRL
    CHARACTER(LEN=256)                  :: choice_d18O_inverse_init_ANT

    CHARACTER(LEN=256)                  :: filename_d18O_inverse_init_NAM
    CHARACTER(LEN=256)                  :: filename_d18O_inverse_init_EAS
    CHARACTER(LEN=256)                  :: filename_d18O_inverse_init_GRL
    CHARACTER(LEN=256)                  :: filename_d18O_inverse_init_ANT

    ! Ice dynamics - velocity
    ! =======================

    CHARACTER(LEN=256)                  :: choice_ice_dynamics
    REAL(dp)                            :: n_flow
    REAL(dp)                            :: m_enh_sheet
    REAL(dp)                            :: m_enh_shelf
    CHARACTER(LEN=256)                  :: choice_ice_margin
    LOGICAL                             :: include_SSADIVA_crossterms
    LOGICAL                             :: do_GL_subgrid_friction
    LOGICAL                             :: do_smooth_geometry
    REAL(dp)                            :: r_smooth_geometry

    ! Some parameters for numerically solving the SSA/DIVA
    REAL(dp)                            :: DIVA_visc_it_norm_dUV_tol
    INTEGER                             :: DIVA_visc_it_nit
    REAL(dp)                            :: DIVA_visc_it_relax
    REAL(dp)                            :: DIVA_epsilon_sq_0
    REAL(dp)                            :: DIVA_visc_eff_min
    REAL(dp)                            :: DIVA_beta_max
    REAL(dp)                            :: DIVA_beta_min_thin_ice
    REAL(dp)                            :: DIVA_beta_min_thin_ice_H0
    REAL(dp)                            :: DIVA_vel_max
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_u_west
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_u_east
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_u_south
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_u_north
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_v_west
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_v_east
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_v_south
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_v_north
    CHARACTER(LEN=256)                  :: DIVA_choice_matrix_solver
    INTEGER                             :: DIVA_SOR_nit
    REAL(dp)                            :: DIVA_SOR_tol
    REAL(dp)                            :: DIVA_SOR_omega
    REAL(dp)                            :: DIVA_PETSc_rtol
    REAL(dp)                            :: DIVA_PETSc_abstol
    LOGICAL                             :: do_read_velocities_from_restart
    ! Ice dynamics - time integration
    ! ===============================

    CHARACTER(LEN=256)                  :: choice_timestepping
    CHARACTER(LEN=256)                  :: choice_ice_integration_method
    CHARACTER(LEN=256)                  :: dHi_choice_matrix_solver
    INTEGER                             :: dHi_SOR_nit
    REAL(dp)                            :: dHi_SOR_tol
    REAL(dp)                            :: dHi_SOR_omega
    REAL(dp)                            :: dHi_PETSc_rtol
    REAL(dp)                            :: dHi_PETSc_abstol

    ! Predictor-corrector ice-thickness update
    REAL(dp)                            :: pc_epsilon
    REAL(dp)                            :: pc_k_I
    REAL(dp)                            :: pc_k_p
    REAL(dp)                            :: pc_eta_min
    INTEGER                             :: pc_max_timestep_iterations
    REAL(dp)                            :: pc_redo_tol

    ! Ice thickness boundary conditions
    CHARACTER(LEN=256)                  :: ice_thickness_west_BC
    CHARACTER(LEN=256)                  :: ice_thickness_east_BC
    CHARACTER(LEN=256)                  :: ice_thickness_south_BC
    CHARACTER(LEN=256)                  :: ice_thickness_north_BC
    CHARACTER(LEN=256)                  :: choice_mask_noice_NAM
    CHARACTER(LEN=256)                  :: choice_mask_noice_EAS
    CHARACTER(LEN=256)                  :: choice_mask_noice_GRL
    CHARACTER(LEN=256)                  :: choice_mask_noice_ANT

    ! Fixed/delayed ice thickness evolution
    REAL(dp)                            :: fixed_sheet_geometry
    REAL(dp)                            :: fixed_shelf_geometry
    REAL(dp)                            :: fixed_grounding_line_g
    REAL(dp)                            :: fixed_grounding_line_f
    REAL(dp)                            :: fixed_decay_t_start
    REAL(dp)                            :: fixed_decay_t_end
    REAL(dp)                            :: relax_thick_t_start
    REAL(dp)                            :: relax_thick_t_end

    ! Target dHi_dt during model spinup
    LOGICAL                             :: do_target_dhdt
    CHARACTER(LEN=256)                  :: target_dhdt_filename
    REAL(dp)                            :: target_dhdt_t_end

    ! Ice dynamics - basal conditions and sliding
    ! ===========================================

    ! Sliding laws
    CHARACTER(LEN=256)                  :: choice_sliding_law
    CHARACTER(LEN=256)                  :: choice_idealised_sliding_law
    REAL(dp)                            :: slid_delta_v
    REAL(dp)                            :: slid_Weertman_m
    REAL(dp)                            :: slid_Coulomb_reg_q_plastic
    REAL(dp)                            :: slid_Coulomb_reg_u_threshold
    REAL(dp)                            :: slid_ZI_ut
    REAL(dp)                            :: slid_ZI_p
    LOGICAL                             :: do_slid_ZI_no_angle
    LOGICAL                             :: include_basal_freezing
    REAL(dp)                            :: deltaT_basal_freezing
    REAL(dp)                            :: subgrid_friction_exponent

    ! Basal hydrology
    CHARACTER(LEN=256)                  :: choice_basal_hydrology
    REAL(dp)                            :: Martin2011_hydro_N_lim
    REAL(dp)                            :: Martin2011_hydro_Hb_min
    REAL(dp)                            :: Martin2011_hydro_Hb_max

    ! Basal roughness / friction
    CHARACTER(LEN=256)                  :: choice_basal_roughness
    REAL(dp)                            :: uniform_Weertman_beta_sq
    REAL(dp)                            :: uniform_Coulomb_phi_fric
    REAL(dp)                            :: uniform_Tsai2015_alpha_sq
    REAL(dp)                            :: uniform_Tsai2015_beta_sq
    REAL(dp)                            :: uniform_Schoof2005_alpha_sq
    REAL(dp)                            :: uniform_Schoof2005_beta_sq
    CHARACTER(LEN=256)                  :: choice_param_basal_roughness
    REAL(dp)                            :: Martin2011till_phi_Hb_min
    REAL(dp)                            :: Martin2011till_phi_Hb_max
    REAL(dp)                            :: Martin2011till_phi_min
    REAL(dp)                            :: Martin2011till_phi_max
    CHARACTER(LEN=256)                  :: basal_roughness_filename
    LOGICAL                             :: do_smooth_phi_restart
    REAL(dp)                            :: r_smooth_phi_restart

    ! Basal inversion
    LOGICAL                             :: do_BIVgeo
    REAL(dp)                            :: BIVgeo_t_start
    REAL(dp)                            :: BIVgeo_t_end
    CHARACTER(LEN=256)                  :: choice_BIVgeo_method
    REAL(dp)                            :: BIVgeo_dt
    REAL(dp)                            :: BIVgeo_PDC2012_hinv
    REAL(dp)                            :: BIVgeo_Lipscomb2021_tauc
    REAL(dp)                            :: BIVgeo_Lipscomb2021_H0
    REAL(dp)                            :: BIVgeo_CISMplus_wH
    REAL(dp)                            :: BIVgeo_CISMplus_wu
    REAL(dp)                            :: BIVgeo_CISMplus_tauc
    REAL(dp)                            :: BIVgeo_CISMplus_H0
    REAL(dp)                            :: BIVgeo_CISMplus_u0
    REAL(dp)                            :: BIVgeo_Berends2022_tauc
    REAL(dp)                            :: BIVgeo_Berends2022_H0
    REAL(dp)                            :: BIVgeo_Berends2022_u0
    REAL(dp)                            :: BIVgeo_Berends2022_Hi_scale
    REAL(dp)                            :: BIVgeo_Berends2022_u_scale
    REAL(dp)                            :: BIVgeo_Berends2022_phimin
    REAL(dp)                            :: BIVgeo_Berends2022_phimax
    CHARACTER(LEN=256)                  :: BIVgeo_target_velocity_filename
    CHARACTER(LEN=256)                  :: BIVgeo_filename_output
    REAL(dp)                            :: BIVgeo_Bernales2017_hinv
    REAL(dp)                            :: BIVgeo_Bernales2017_tol_diff
    REAL(dp)                            :: BIVgeo_Bernales2017_tol_frac
    REAL(dp)                            :: BIVgeo_Pien2023_H0
    REAL(dp)                            :: BIVgeo_Pien2023_tau
    REAL(dp)                            :: BIVgeo_Pien2023_r
    REAL(dp)                            :: BIVgeo_Pien2023_lowerHb
    REAL(dp)                            :: BIVgeo_Pien2023_upperHb
    REAL(dp)                            :: BIVgeo_Pien2023_min
    REAL(dp)                            :: BIVgeo_Pien2023_max

    ! Ice dynamics - calving
    ! ======================

    CHARACTER(LEN=256)                  :: choice_calving_law
    REAL(dp)                            :: calving_threshold_thickness
    LOGICAL                             :: do_remove_shelves
    LOGICAL                             :: remove_shelves_larger_than_PD
    LOGICAL                             :: continental_shelf_calving
    REAL(dp)                            :: continental_shelf_min_height

    ! Thermodynamics and rheology
    ! ===========================

    CHARACTER(LEN=256)                  :: choice_initial_ice_temperature
    REAL(dp)                            :: uniform_ice_temperature
    CHARACTER(LEN=256)                  :: choice_thermo_model
    CHARACTER(LEN=256)                  :: choice_ice_rheology
    REAL(dp)                            :: uniform_flow_factor
    CHARACTER(LEN=256)                  :: choice_ice_heat_capacity
    REAL(dp)                            :: uniform_ice_heat_capacity
    CHARACTER(LEN=256)                  :: choice_ice_thermal_conductivity
    REAL(dp)                            :: uniform_ice_thermal_conductivity

    ! Climate
    ! =======

    CHARACTER(LEN=256)                  :: choice_climate_model
    CHARACTER(LEN=256)                  :: choice_idealised_climate

    ! Folder with NetCDF files containing direct climate forcing
    CHARACTER(LEN=256)                  :: direct_climate_foldername_NAM
    CHARACTER(LEN=256)                  :: direct_climate_foldername_EAS
    CHARACTER(LEN=256)                  :: direct_climate_foldername_GRL
    CHARACTER(LEN=256)                  :: direct_climate_foldername_ANT

    CHARACTER(LEN=256)                  :: direct_climate_basefilename_NAM
    CHARACTER(LEN=256)                  :: direct_climate_basefilename_EAS
    CHARACTER(LEN=256)                  :: direct_climate_basefilename_GRL
    CHARACTER(LEN=256)                  :: direct_climate_basefilename_ANT

    ! Whether or not to apply geometry corrections to direct temperature/precipitation
    LOGICAL                             :: do_direct_climate_geo_corr

    ! NetCDF file containing the present-day observed climate (e.g. ERA40)
    CHARACTER(LEN=256)                  :: filename_PD_obs_climate

    ! GCM snapshots in the matrix_warm_cold option
    CHARACTER(LEN=256)                  :: filename_climate_snapshot_PI
    CHARACTER(LEN=256)                  :: filename_climate_snapshot_warm
    CHARACTER(LEN=256)                  :: filename_climate_snapshot_cold

    ! Ice and ocean mask from GCM snapshots
    CHARACTER(LEN=256)                  :: reference_mask_method
    
    CHARACTER(LEN=256)                  :: filename_snapshot_mask_PD_obs 
    CHARACTER(LEN=256)                  :: filename_snapshot_mask_GCM_PI
    CHARACTER(LEN=256)                  :: filename_snapshot_mask_GCM_warm
    CHARACTER(LEN=256)                  :: filename_snapshot_mask_GCM_cold  
  
    ! Lapse rate
    REAL(dp)                            :: constant_lapserate

    ! Scaling factor for CO2 vs ice weights
    REAL(dp)                            :: climate_matrix_CO2vsice_NAM
    REAL(dp)                            :: climate_matrix_CO2vsice_EAS
    REAL(dp)                            :: climate_matrix_CO2vsice_GRL
    REAL(dp)                            :: climate_matrix_CO2vsice_ANT

    ! Orbit time and CO2 concentration of the warm and cold snapshots
    REAL(dp)                            :: matrix_high_CO2_level
    REAL(dp)                            :: matrix_low_CO2_level
    REAL(dp)                            :: matrix_warm_orbit_time
    REAL(dp)                            :: matrix_cold_orbit_time

    ! Whether or not to apply a bias correction to the GCM snapshots
    LOGICAL                             :: climate_matrix_biascorrect_warm
    LOGICAL                             :: climate_matrix_biascorrect_cold

    LOGICAL                             :: switch_glacial_index_precip

    ! Ocean
    ! =====

    CHARACTER(LEN=256)                  :: choice_ocean_model
    CHARACTER(LEN=256)                  :: choice_idealised_ocean
    LOGICAL                             :: do_ocean_temperature_inversion
    REAL(dp)                            :: ocean_temperature_inv_t_start
    REAL(dp)                            :: ocean_temperature_inv_t_end
    REAL(dp)                            :: ocean_temperature_inv_dT0
    REAL(dp)                            :: ocean_temperature_inv_H0
    REAL(dp)                            :: ocean_temperature_inv_tau
    REAL(dp)                            :: ocean_temperature_inv_min
    REAL(dp)                            :: ocean_temperature_inv_max
    CHARACTER(LEN=256)                  :: inverted_ocean_filename_output

    ! NetCDF file containing the present-day observed ocean (WOA18) (NetCDF)
    CHARACTER(LEN=256)                  :: filename_PD_obs_ocean
    CHARACTER(LEN=256)                  :: name_ocean_temperature_obs
    CHARACTER(LEN=256)                  :: name_ocean_salinity_obs

    ! GCM snapshots in the matrix_warm_cold option
    CHARACTER(LEN=256)                  :: filename_GCM_ocean_snapshot_PI
    CHARACTER(LEN=256)                  :: filename_GCM_ocean_snapshot_warm
    CHARACTER(LEN=256)                  :: filename_GCM_ocean_snapshot_cold
    CHARACTER(LEN=256)                  :: name_ocean_temperature_GCM
    CHARACTER(LEN=256)                  :: name_ocean_salinity_GCM

    ! Uniform ocean temperature values used when choice_ocean_model = "uniform_warm_cold"
    REAL(dp)                            :: ocean_temperature_PD
    REAL(dp)                            :: ocean_temperature_cold
    REAL(dp)                            :: ocean_temperature_warm

    ! Parameters used when choice_ocean_model = "matrix_warm_cold"
    CHARACTER(LEN=256)                  :: choice_ocean_vertical_grid
    REAL(dp)                            :: ocean_vertical_grid_max_depth
    REAL(dp)                            :: ocean_regular_grid_dz
    INTEGER                             :: nz_ocean ! NOTE: nz_ocean and z_ocean cannot be set through the config file, but are filled in by the "initialise_ocean_vertical_grid" in the ocean_module!
    REAL(dp), DIMENSION(:), ALLOCATABLE :: z_ocean
    LOGICAL                             :: do_ocean_extrap
    CHARACTER(LEN=256)                  :: ocean_extrap_dir
    REAL(dp)                            :: ocean_extrap_res
    REAL(dp)                            :: ocean_extrap_Gauss_sigma
    CHARACTER(LEN=256)                  :: ocean_extrap_hires_geo_filename_NAM
    CHARACTER(LEN=256)                  :: ocean_extrap_hires_geo_filename_EAS
    CHARACTER(LEN=256)                  :: ocean_extrap_hires_geo_filename_GRL
    CHARACTER(LEN=256)                  :: ocean_extrap_hires_geo_filename_ANT
    REAL(dp)                            :: ocean_w_tot_hist_averaging_window

    ! Scaling factor for CO2 vs ice weights
    REAL(dp)                            :: ocean_matrix_CO2vsice_NAM
    REAL(dp)                            :: ocean_matrix_CO2vsice_EAS
    REAL(dp)                            :: ocean_matrix_CO2vsice_GRL
    REAL(dp)                            :: ocean_matrix_CO2vsice_ANT

    ! Basin-dependent linear temperature profiles (used when choice_idealised_ocean = "linear_per_basin")
    REAL(dp), DIMENSION(100)            :: ocean_T_surf_per_basin
    REAL(dp), DIMENSION(100)            :: ocean_dT_dz_per_basin
    LOGICAL                             :: do_invert_linear_per_basin

    ! Apply anomalies to baseline (used when choice_idealised_ocean = "anomalies")
    CHARACTER(LEN=256)                  :: ocean_filename_baseline
    CHARACTER(LEN=256)                  :: ocean_foldername_aTO
    CHARACTER(LEN=256)                  :: ocean_basefilename_aTO
    CHARACTER(LEN=256)                  :: ocean_foldername_aSO
    CHARACTER(LEN=256)                  :: ocean_basefilename_aSO

    ! Surface mass balance
    ! ====================

    CHARACTER(LEN=256)                  :: choice_SMB_model
    CHARACTER(LEN=256)                  :: choice_idealised_SMB
    REAL(dp)                            :: SMB_uniform

    ! NetCDF file containing direct global/regional SMB forcing
    CHARACTER(LEN=256)                  :: filename_direct_global_SMB
    CHARACTER(LEN=256)                  :: filename_direct_regional_SMB_NAM
    CHARACTER(LEN=256)                  :: filename_direct_regional_SMB_EAS
    CHARACTER(LEN=256)                  :: filename_direct_regional_SMB_GRL
    CHARACTER(LEN=256)                  :: filename_direct_regional_SMB_ANT

    ! NetCDF file containing a direct time-less SMB snapshot
    CHARACTER(LEN=256)                  :: filename_SMB_snapshot_NAM
    CHARACTER(LEN=256)                  :: filename_SMB_snapshot_EAS
    CHARACTER(LEN=256)                  :: filename_SMB_snapshot_GRL
    CHARACTER(LEN=256)                  :: filename_SMB_snapshot_ANT

    ! Tuning parameters for the IMAU-ITM SMB model
    CHARACTER(LEN=256)                  :: SMB_IMAUITM_choice_init_firn_NAM
    CHARACTER(LEN=256)                  :: SMB_IMAUITM_choice_init_firn_EAS
    CHARACTER(LEN=256)                  :: SMB_IMAUITM_choice_init_firn_GRL
    CHARACTER(LEN=256)                  :: SMB_IMAUITM_choice_init_firn_ANT
    REAL(dp)                            :: SMB_IMAUITM_initial_firn_thickness
    REAL(dp)                            :: SMB_IMAUITM_C_abl_constant_NAM
    REAL(dp)                            :: SMB_IMAUITM_C_abl_constant_EAS
    REAL(dp)                            :: SMB_IMAUITM_C_abl_constant_GRL
    REAL(dp)                            :: SMB_IMAUITM_C_abl_constant_ANT
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Ts_NAM
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Ts_EAS
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Ts_GRL
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Ts_ANT
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Q_NAM
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Q_EAS
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Q_GRL
    REAL(dp)                            :: SMB_IMAUITM_C_abl_Q_ANT
    REAL(dp)                            :: SMB_IMAUITM_C_refr_NAM
    REAL(dp)                            :: SMB_IMAUITM_C_refr_EAS
    REAL(dp)                            :: SMB_IMAUITM_C_refr_GRL
    REAL(dp)                            :: SMB_IMAUITM_C_refr_ANT

    ! ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing
    ! ==============================================================

    CHARACTER(LEN=256)                  :: ISMIP_forcing_filename_baseline
    CHARACTER(LEN=256)                  :: ISMIP_forcing_foldername_aSMB
    CHARACTER(LEN=256)                  :: ISMIP_forcing_basefilename_aSMB
    CHARACTER(LEN=256)                  :: ISMIP_forcing_foldername_dSMBdz
    CHARACTER(LEN=256)                  :: ISMIP_forcing_basefilename_dSMBdz
    CHARACTER(LEN=256)                  :: ISMIP_forcing_foldername_aST
    CHARACTER(LEN=256)                  :: ISMIP_forcing_basefilename_aST
    CHARACTER(LEN=256)                  :: ISMIP_forcing_foldername_dSTdz
    CHARACTER(LEN=256)                  :: ISMIP_forcing_basefilename_dSTdz

    ! Basal mass balance - sub-shelf melt
    ! ===================================

    CHARACTER(LEN=256)                  :: choice_BMB_shelf_model
    CHARACTER(LEN=256)                  :: choice_idealised_BMB_shelf
    CHARACTER(LEN=256)                  :: choice_BMB_sheet_model
    REAL(dp)                            :: BMB_shelf_uniform
    REAL(dp)                            :: BMB_sheet_uniform
    CHARACTER(LEN=256)                  :: choice_BMB_subgrid
    LOGICAL                             :: do_asynchronous_BMB
    REAL(dp)                            :: BMB_max
    REAL(dp)                            :: BMB_min

    CHARACTER(LEN=256)                  :: choice_basin_scheme_NAM
    CHARACTER(LEN=256)                  :: choice_basin_scheme_EAS
    CHARACTER(LEN=256)                  :: choice_basin_scheme_GRL
    CHARACTER(LEN=256)                  :: choice_basin_scheme_ANT
    CHARACTER(LEN=256)                  :: filename_basins_NAM
    CHARACTER(LEN=256)                  :: filename_basins_EAS
    CHARACTER(LEN=256)                  :: filename_basins_GRL
    CHARACTER(LEN=256)                  :: filename_basins_ANT
    LOGICAL                             :: do_merge_basins_ANT
    LOGICAL                             :: do_merge_basins_GRL

    CHARACTER(LEN=256)                  :: choice_BMB_shelf_amplification
    INTEGER                             :: basin_BMB_amplification_n_ANT
    REAL(dp), DIMENSION(:), ALLOCATABLE :: basin_BMB_amplification_factor_ANT
    INTEGER                             :: basin_BMB_amplification_n_GRL
    REAL(dp), DIMENSION(:), ALLOCATABLE :: basin_BMB_amplification_factor_GRL

    ! Parameters for the three simple melt parameterisations from Favier et al. (2019)
    REAL(dp)                            :: BMB_Favier2019_lin_GammaT
    REAL(dp)                            :: BMB_Favier2019_quad_GammaT
    REAL(dp)                            :: BMB_Favier2019_Mplus_GammaT

    ! Parameters for the Lazeroms et al. (2018) plume-parameterisation BMB model
    REAL(dp)                            :: BMB_Lazeroms2018_GammaT
    CHARACTER(LEN=256)                  :: BMB_Lazeroms2018_find_GL_scheme

    ! Parameters for the PICO BMB model
    INTEGER                             :: BMB_PICO_nboxes
    REAL(dp)                            :: BMB_PICO_GammaTstar

    ! Parameters for the LADDIE model
    CHARACTER(LEN=256)                  :: filename_BMB_LADDIE

    ! Parameters for the ANICE_legacy sub-shelf melt model
    REAL(dp)                            :: T_ocean_mean_PD_NAM
    REAL(dp)                            :: T_ocean_mean_PD_EAS
    REAL(dp)                            :: T_ocean_mean_PD_GRL
    REAL(dp)                            :: T_ocean_mean_PD_ANT
    REAL(dp)                            :: T_ocean_mean_cold_NAM
    REAL(dp)                            :: T_ocean_mean_cold_EAS
    REAL(dp)                            :: T_ocean_mean_cold_GRL
    REAL(dp)                            :: T_ocean_mean_cold_ANT
    REAL(dp)                            :: T_ocean_mean_warm_NAM
    REAL(dp)                            :: T_ocean_mean_warm_EAS
    REAL(dp)                            :: T_ocean_mean_warm_GRL
    REAL(dp)                            :: T_ocean_mean_warm_ANT

    REAL(dp)                            :: BMB_deepocean_PD_NAM
    REAL(dp)                            :: BMB_deepocean_PD_EAS
    REAL(dp)                            :: BMB_deepocean_PD_GRL
    REAL(dp)                            :: BMB_deepocean_PD_ANT
    REAL(dp)                            :: BMB_deepocean_cold_NAM
    REAL(dp)                            :: BMB_deepocean_cold_EAS
    REAL(dp)                            :: BMB_deepocean_cold_GRL
    REAL(dp)                            :: BMB_deepocean_cold_ANT
    REAL(dp)                            :: BMB_deepocean_warm_NAM
    REAL(dp)                            :: BMB_deepocean_warm_EAS
    REAL(dp)                            :: BMB_deepocean_warm_GRL
    REAL(dp)                            :: BMB_deepocean_warm_ANT

    REAL(dp)                            :: BMB_shelf_exposed_PD_NAM
    REAL(dp)                            :: BMB_shelf_exposed_PD_EAS
    REAL(dp)                            :: BMB_shelf_exposed_PD_GRL
    REAL(dp)                            :: BMB_shelf_exposed_PD_ANT
    REAL(dp)                            :: BMB_shelf_exposed_cold_NAM
    REAL(dp)                            :: BMB_shelf_exposed_cold_EAS
    REAL(dp)                            :: BMB_shelf_exposed_cold_GRL
    REAL(dp)                            :: BMB_shelf_exposed_cold_ANT
    REAL(dp)                            :: BMB_shelf_exposed_warm_NAM
    REAL(dp)                            :: BMB_shelf_exposed_warm_EAS
    REAL(dp)                            :: BMB_shelf_exposed_warm_GRL
    REAL(dp)                            :: BMB_shelf_exposed_warm_ANT

    REAL(dp)                            :: subshelf_melt_factor_NAM
    REAL(dp)                            :: subshelf_melt_factor_EAS
    REAL(dp)                            :: subshelf_melt_factor_GRL
    REAL(dp)                            :: subshelf_melt_factor_ANT

    REAL(dp)                            :: deep_ocean_threshold_depth_NAM
    REAL(dp)                            :: deep_ocean_threshold_depth_EAS
    REAL(dp)                            :: deep_ocean_threshold_depth_GRL
    REAL(dp)                            :: deep_ocean_threshold_depth_ANT

    ! Englacial isotope tracing
    ! ========================

    CHARACTER(LEN=256)                  :: choice_ice_isotopes_model
    REAL(dp)                            :: uniform_ice_d18O

    ! Sea level and GIA
    ! =================

    LOGICAL                             :: do_ocean_floodfill
    CHARACTER(LEN=256)                  :: choice_sealevel_model
    REAL(dp)                            :: fixed_sealevel
    CHARACTER(LEN=256)                  :: filename_sealevel_record
    INTEGER                             :: sealevel_record_length

    CHARACTER(LEN=256)                  :: choice_GIA_model
    REAL(dp)                            :: dx_GIA
    REAL(dp)                            :: ELRA_lithosphere_flex_rigidity
    REAL(dp)                            :: ELRA_bedrock_relaxation_time
    REAL(dp)                            :: ELRA_mantle_density

    ! SELEN
    ! =====

    LOGICAL                             :: SELEN_run_at_t_start
    INTEGER                             :: SELEN_n_TDOF_iterations
    INTEGER                             :: SELEN_n_recursion_iterations
    LOGICAL                             :: SELEN_use_rotational_feedback
    INTEGER                             :: SELEN_n_harmonics
    LOGICAL                             :: SELEN_display_progress

    CHARACTER(LEN=256)                  :: SELEN_dir
    CHARACTER(LEN=256)                  :: SELEN_global_topo_filename
    CHARACTER(LEN=256)                  :: SELEN_TABOO_init_filename
    CHARACTER(LEN=256)                  :: SELEN_LMJ_VALUES_filename

    INTEGER                             :: SELEN_irreg_time_n
    REAL(dp), DIMENSION(:), ALLOCATABLE :: SELEN_irreg_time_window

    REAL(dp)                            :: SELEN_lith_thickness
    INTEGER                             :: SELEN_visc_n
    REAL(dp), DIMENSION(:), ALLOCATABLE :: SELEN_visc_prof

    INTEGER                             :: SELEN_TABOO_CDE
    INTEGER                             :: SELEN_TABOO_TLOVE
    INTEGER                             :: SELEN_TABOO_DEG1
    REAL(dp)                            :: SELEN_TABOO_RCMB

    ! Some derived values
    INTEGER                             :: SELEN_i1, SELEN_i2           ! Parallelisation of loops over global grid pixels
    INTEGER                             :: SELEN_j1, SELEN_j2           ! Parallelisation of loops over harmonic degrees
    REAL(dp)                            :: SELEN_alfa
    INTEGER                             :: SELEN_jmax
    INTEGER                             :: SELEN_reg_time_n

    ! Which data fields will be written to the help_fields output file
    ! ================================================================

    CHARACTER(LEN=256)                  :: help_field_01
    CHARACTER(LEN=256)                  :: help_field_02
    CHARACTER(LEN=256)                  :: help_field_03
    CHARACTER(LEN=256)                  :: help_field_04
    CHARACTER(LEN=256)                  :: help_field_05
    CHARACTER(LEN=256)                  :: help_field_06
    CHARACTER(LEN=256)                  :: help_field_07
    CHARACTER(LEN=256)                  :: help_field_08
    CHARACTER(LEN=256)                  :: help_field_09
    CHARACTER(LEN=256)                  :: help_field_10
    CHARACTER(LEN=256)                  :: help_field_11
    CHARACTER(LEN=256)                  :: help_field_12
    CHARACTER(LEN=256)                  :: help_field_13
    CHARACTER(LEN=256)                  :: help_field_14
    CHARACTER(LEN=256)                  :: help_field_15
    CHARACTER(LEN=256)                  :: help_field_16
    CHARACTER(LEN=256)                  :: help_field_17
    CHARACTER(LEN=256)                  :: help_field_18
    CHARACTER(LEN=256)                  :: help_field_19
    CHARACTER(LEN=256)                  :: help_field_20
    CHARACTER(LEN=256)                  :: help_field_21
    CHARACTER(LEN=256)                  :: help_field_22
    CHARACTER(LEN=256)                  :: help_field_23
    CHARACTER(LEN=256)                  :: help_field_24
    CHARACTER(LEN=256)                  :: help_field_25
    CHARACTER(LEN=256)                  :: help_field_26
    CHARACTER(LEN=256)                  :: help_field_27
    CHARACTER(LEN=256)                  :: help_field_28
    CHARACTER(LEN=256)                  :: help_field_29
    CHARACTER(LEN=256)                  :: help_field_30
    CHARACTER(LEN=256)                  :: help_field_31
    CHARACTER(LEN=256)                  :: help_field_32
    CHARACTER(LEN=256)                  :: help_field_33
    CHARACTER(LEN=256)                  :: help_field_34
    CHARACTER(LEN=256)                  :: help_field_35
    CHARACTER(LEN=256)                  :: help_field_36
    CHARACTER(LEN=256)                  :: help_field_37
    CHARACTER(LEN=256)                  :: help_field_38
    CHARACTER(LEN=256)                  :: help_field_39
    CHARACTER(LEN=256)                  :: help_field_40
    CHARACTER(LEN=256)                  :: help_field_41
    CHARACTER(LEN=256)                  :: help_field_42
    CHARACTER(LEN=256)                  :: help_field_43
    CHARACTER(LEN=256)                  :: help_field_44
    CHARACTER(LEN=256)                  :: help_field_45
    CHARACTER(LEN=256)                  :: help_field_46
    CHARACTER(LEN=256)                  :: help_field_47
    CHARACTER(LEN=256)                  :: help_field_48
    CHARACTER(LEN=256)                  :: help_field_49
    CHARACTER(LEN=256)                  :: help_field_50

    ! Values to be filled into the total mask (used only for diagnostic output)
    ! ==========================================================================

    INTEGER                             :: type_land
    INTEGER                             :: type_ocean
    INTEGER                             :: type_lake
    INTEGER                             :: type_sheet
    INTEGER                             :: type_shelf
    INTEGER                             :: type_coast
    INTEGER                             :: type_margin
    INTEGER                             :: type_groundingline
    INTEGER                             :: type_calvingfront

    ! The output directory
    ! ====================

    CHARACTER(LEN=256)                  :: output_dir

  END TYPE constants_type




  ! Since some of the TABOO routines have variables named C (thanks, Giorgio...),
  ! we cannot use the regular config structure there. Collect the required config
  ! parameters into a smaller separate structure called C_TABOO
  TYPE constants_type_TABOO

    INTEGER                             :: IMODE            ! SELEN integration mode
    INTEGER                             :: NV               ! Number of viscoelastic layers
    REAL(dp), DIMENSION(:), ALLOCATABLE :: VSC              ! Viscosity profile
    INTEGER                             :: CDE              ! Code of the model (see taboo for explanation)
    INTEGER                             :: TLOVE            ! Tidal love numbers yes/no
    INTEGER                             :: DEG1             ! Tidal love numbers degree
    REAL(dp)                            :: LTH              ! Lithospheric thickness [km]
    REAL(dp)                            :: RCMB             ! Radius of CMB (km)

  END TYPE constants_type_TABOO


  ! ===============================================
  ! "C" is an instance of the "constants_type" type
  ! ===============================================

  TYPE(constants_type      ), SAVE :: C
  TYPE(constants_type_TABOO), SAVE :: C_TABOO

CONTAINS

  SUBROUTINE initialise_model_configuration( version_number)
    ! Initialise the C (configuration) structure from one or two external config text files,
    ! set up the output directory (either procedurally from the current date, or directly
    ! from the config-specified folder name), and copy the config file(s) there.

    ! In/output variables:
    CHARACTER(LEN=256),                  INTENT(IN)    :: version_number

    ! Local variables:
    INTEGER                                            :: ierr, process_rank, number_of_processes, p
    LOGICAL                                            :: master
    CHARACTER(LEN=256)                                 :: config_filename, template_filename, variation_filename, config_mode
    INTEGER                                            :: i,n
    CHARACTER(LEN=20)                                  :: output_dir_procedural
    LOGICAL                                            :: ex

    ! Get rank of current process and total number of processes
    ! (needed because the configuration_module cannot access the par structure)
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, process_rank, ierr)
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, number_of_processes, ierr)
    master = (process_rank == 0)

  ! ===== Set up the config structure =====
  ! =======================================

    ! The name(s) of the config file(s) are provided as input arguments when calling the IMAU_ICE_program
    ! executable. After calling MPI_INIT, only the master process "sees" these arguments, so they need to be
    ! broadcast to the other processes.

    IF (master) THEN

      config_filename       = ''
      template_filename     = ''
      variation_filename    = ''
      config_mode           = ''

      IF     (iargc() == 0) THEN
        CALL crash('IMAU-ICE v' // TRIM( version_number) // ' needs at least one config file to run!')
      ELSEIF (iargc() == 1) THEN
        ! Run the model with a single config file

        CALL getarg( 1, config_filename)
        config_mode = 'single_config'

      ELSEIF (iargc() == 2) THEN
        ! Run the model with two config files (template+variation)

        CALL getarg( 1, template_filename )
        CALL getarg( 2, variation_filename)
        config_mode = 'template+variation'

      ELSE
        CALL crash('IMAU-ICE v' // TRIM( version_number) // ' can take either one or two config files to run!')
      END IF

    END IF ! IF (master) THEN

    CALL MPI_BCAST( config_filename,    256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( template_filename,  256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( variation_filename, 256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( config_mode,        256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Let each of the processors read the config file in turns so there's no access conflicts
    IF (config_mode == 'single_config') THEN
      ! Read only a single config file

      DO p = 0, number_of_processes-1
        IF (p == process_rank) THEN

          ! Read the external file, use a Fortran NAMELIST to overwrite the default
          ! values of the XXX_config variables
          CALL read_config_file( config_filename)

          ! Copy values from the XXX_config variables to the C structure
          CALL copy_variables_to_struct

        END IF
        CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
      END DO

    ELSEIF (config_mode == 'template+variation') THEN
      ! Read two config file consecutively: one "template" and one "variation"

      DO p = 0, number_of_processes-1
        IF (p == process_rank) THEN

          ! Read the external file, use a Fortran NAMELIST to overwrite the default
          ! values of the XXX_config variables

          ! First the template, then the variation
          CALL read_config_file( template_filename)
          CALL read_config_file( variation_filename)

          ! Copy values from the XXX_config variables to the C structure
          CALL copy_variables_to_struct

        END IF
        CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
      END DO

    ELSE ! IF (config_mode == 'single_config') THEN
      CALL crash(' unknown config_mode "' // TRIM( config_mode) // '"!')
    END IF ! IF (config_mode == 'single_config') THEN

  ! ===== Set up the output directory =====
  ! =======================================

    ! First get the name of the output directory (either procedural, or provided in the config file)

    DO n = 1, 256
      C%output_dir(n:n) = ' '
    END DO

    IF (C%create_procedural_output_dir) THEN
      ! Automatically create an output directory with a procedural name (e.g. results_20210720_001/)

      IF (master) THEN
        CALL get_procedural_output_dir_name( output_dir_procedural)
        C%output_dir(1:21) = TRIM(output_dir_procedural) // '/'
      END IF
      CALL MPI_BCAST( C%output_dir, 256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ELSE
      ! Use the provided name (return an error if this directory already exists)

      C%output_dir = TRIM(C%fixed_output_dir) // TRIM(C%fixed_output_dir_suffix) // '/'

      INQUIRE( FILE = TRIM(C%output_dir)//'/.', EXIST=ex)
      IF (ex) THEN
        CALL crash(' fixed_output_dir_config "' // TRIM( C%output_dir) // '" already exists!')
      END IF

    END IF

    ! Create the directory
    IF (master) THEN
      CALL system('mkdir ' // TRIM(C%output_dir))
      WRITE(0,*) ''
      WRITE(0,*) ' Output directory: ', TRIM(C%output_dir)
      WRITE(0,*) ''
    END IF
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)

    ! Copy the config file to the output directory
    IF (master) THEN
      IF     (config_mode == 'single_config') THEN
        CALL system('cp ' // config_filename    // ' ' // TRIM(C%output_dir))
      ELSEIF (config_mode == 'template+variation') THEN
        CALL system('cp ' // template_filename  // ' ' // TRIM(C%output_dir))
        CALL system('cp ' // variation_filename // ' ' // TRIM(C%output_dir))
      ELSE
        CALL crash(' unknown config_mode "' // TRIM( config_mode) // '"!')
      END IF ! IF (config_mode == 'single_config') THEN
    END IF ! IF (master) THEN
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)

    ! Set up the subroutine resource tracker
    ! ======================================

    ! Allocate space to track up to 1,000 subroutines. That should be enough for a while...
    n = 7000
    ALLOCATE( resource_tracker( n))

    ! Initialise values
    DO i = 1, n
      resource_tracker( i)%routine_path = 'subroutine_placeholder'
      resource_tracker( i)%tstart       = 0._dp
      resource_tracker( i)%tcomp       = 0._dp
    END DO

  END SUBROUTINE initialise_model_configuration

  SUBROUTINE read_config_file( config_filename)
    ! Use a NAMELIST containing all the "_config" variables to read
    ! an external config file, and overwrite the default values of
    ! the specified variables with the values from the file.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),INTENT(IN) :: config_filename
    ! Local variables:
    CHARACTER(LEN=256)            :: namelist_filename
    INTEGER, PARAMETER            :: config_unit   = 1337
    INTEGER, PARAMETER            :: namelist_unit = 1338
    INTEGER                       :: ios, ierr, cerr

    ! The NAMELIST that's used to read the external config file.

    NAMELIST /CONFIG/start_time_of_run_config,                        &
                     end_time_of_run_config,                          &
                     dt_coupling_config,                              &
                     dt_max_config,                                   &
                     dt_min_config,                                   &
                     dt_startup_phase_config,                         &
                     dt_thermo_config,                                &
                     dt_climate_config,                               &
                     dt_ocean_config,                                 &
                     dt_SMB_config,                                   &
                     dt_BMB_config,                                   &
                     dt_bedrock_ELRA_config,                          &
                     dt_SELEN_config,                                 &
                     dt_output_config,                                &
                     dt_output_restart_config,                        &
                     do_NAM_config,                                   &
                     do_EAS_config,                                   &
                     do_GRL_config,                                   &
                     do_ANT_config,                                   &
                     SSA_icestream_A_config,                          &
                     SSA_icestream_L_config,                          &
                     SSA_icestream_m_config,                          &
                     SSA_icestream_tantheta_config,                   &
                     SSA_icestream_H_config,                          &
                     ISMIP_HOM_L_config,                              &
                     ISMIP_HOM_E_Arolla_filename_config,              &
                     MISMIPplus_do_tune_A_for_GL_config,              &
                     MISMIPplus_xGL_target_config,                    &
                     MISMIPplus_A_flow_initial_config,                &
                     MISMIPplus_scenario_config,                      &
                     MISOMIP1_scenario_config,                        &
                     create_procedural_output_dir_config,             &
                     fixed_output_dir_config,                         &
                     fixed_output_dir_suffix_config,                  &
                     do_write_regional_scalar_output_config,          &
                     do_write_global_scalar_output_config,            &
                     do_write_regional_scalar_every_timestep_config,  &
                     do_check_for_NaN_config,                         &
                     do_time_display_config,                          &
                     do_write_ISMIP_output_config,                    &
                     ISMIP_output_group_code_config,                  &
                     ISMIP_output_model_code_config,                  &
                     ISMIP_output_experiment_code_config,             &
                     ISMIP_output_basetime_config,                    &
                     lambda_M_NAM_config,                             &
                     phi_M_NAM_config,                                &
                     beta_stereo_NAM_config,                          &
                     xmin_NAM_config,                                 &
                     xmax_NAM_config,                                 &
                     ymin_NAM_config,                                 &
                     ymax_NAM_config,                                 &
                     dx_NAM_config,                                   &
                     lambda_M_EAS_config,                             &
                     phi_M_EAS_config,                                &
                     beta_stereo_EAS_config,                          &
                     xmin_EAS_config,                                 &
                     xmax_EAS_config,                                 &
                     ymin_EAS_config,                                 &
                     ymax_EAS_config,                                 &
                     dx_EAS_config,                                   &
                     lambda_M_GRL_config,                             &
                     phi_M_GRL_config,                                &
                     beta_stereo_GRL_config,                          &
                     xmin_GRL_config,                                 &
                     xmax_GRL_config,                                 &
                     ymin_GRL_config,                                 &
                     ymax_GRL_config,                                 &
                     dx_GRL_config,                                   &
                     lambda_M_ANT_config,                             &
                     phi_M_ANT_config,                                &
                     beta_stereo_ANT_config,                          &
                     xmin_ANT_config,                                 &
                     xmax_ANT_config,                                 &
                     ymin_ANT_config,                                 &
                     ymax_ANT_config,                                 &
                     dx_ANT_config,                                   &
                     nz_config,                                       &
                     zeta_config,                                     &
                     choice_refgeo_init_NAM_config,                   &
                     choice_refgeo_init_EAS_config,                   &
                     choice_refgeo_init_GRL_config,                   &
                     choice_refgeo_init_ANT_config,                   &
                     time_to_restart_from_NAM_config,                 &
                     time_to_restart_from_EAS_config,                 &
                     time_to_restart_from_GRL_config,                 &
                     time_to_restart_from_ANT_config,                 &
                     choice_refgeo_init_idealised_config,             &
                     filename_refgeo_init_NAM_config,                 &
                     filename_refgeo_init_EAS_config,                 &
                     filename_refgeo_init_GRL_config,                 &
                     filename_refgeo_init_ANT_config,                 &
                     choice_refgeo_PD_NAM_config,                     &
                     choice_refgeo_PD_EAS_config,                     &
                     choice_refgeo_PD_GRL_config,                     &
                     choice_refgeo_PD_ANT_config,                     &
                     choice_refgeo_PD_idealised_config,               &
                     filename_refgeo_PD_NAM_config,                   &
                     filename_refgeo_PD_EAS_config,                   &
                     filename_refgeo_PD_GRL_config,                   &
                     filename_refgeo_PD_ANT_config,                   &
                     choice_refgeo_GIAeq_NAM_config,                  &
                     choice_refgeo_GIAeq_EAS_config,                  &
                     choice_refgeo_GIAeq_GRL_config,                  &
                     choice_refgeo_GIAeq_ANT_config,                  &
                     choice_refgeo_GIAeq_idealised_config,            &
                     filename_refgeo_GIAeq_NAM_config,                &
                     filename_refgeo_GIAeq_EAS_config,                &
                     filename_refgeo_GIAeq_GRL_config,                &
                     filename_refgeo_GIAeq_ANT_config,                &
                     remove_Lake_Vostok_config,                       &
                     choice_forcing_method_config,                    &
                     choice_insolation_forcing_config,                &
                     static_insolation_time_config,                   &
                     filename_insolation_config,                      &
                     filename_CO2_record_config,                      &
                     CO2_record_length_config,                        &
                     filename_d18O_record_config,                     &
                     d18O_record_length_config,                       &
                     do_combine_CO2_and_insolation_config,            &
                     insolation_weigth_mean_NAM_config,               &
                     insolation_weigth_amplitude_NAM_config,          &
                     insolation_weigth_mean_EAS_config,               &
                     insolation_weigth_amplitude_EAS_config,          &
                     insolation_weigth_mean_GRL_config ,              &
                     insolation_weigth_amplitude_GRL_config,          &
                     insolation_weigth_mean_ANT_config,               &
                     insolation_weigth_amplitude_ANT_config,          &
                     do_climate_matrix_wind_config,                   &
                     choice_geothermal_heat_flux_config,              &
                     constant_geothermal_heat_flux_config,            &
                     filename_geothermal_heat_flux_config,            &
                     do_calculate_benthic_d18O_config,                &
                     dT_deepwater_averaging_window_config,            &
                     dT_deepwater_dT_surf_ratio_config,               &
                     d18O_dT_deepwater_ratio_config,                  &
                     dT_glob_inverse_averaging_window_config,         &
                     inverse_d18O_to_dT_glob_scaling_config,          &
                     CO2_inverse_averaging_window_config,             &
                     inverse_d18O_to_CO2_scaling_config,              &
                     inverse_d18O_to_CO2_initial_CO2_config,          &
                     choice_d18O_inverse_init_NAM_config,             &
                     choice_d18O_inverse_init_EAS_config,             &
                     choice_d18O_inverse_init_GRL_config,             &
                     choice_d18O_inverse_init_ANT_config,             &
                     filename_d18O_inverse_init_NAM_config,           &
                     filename_d18O_inverse_init_EAS_config,           &
                     filename_d18O_inverse_init_GRL_config,           &
                     filename_d18O_inverse_init_ANT_config,           &
                     choice_ice_dynamics_config,                      &
                     n_flow_config,                                   &
                     m_enh_sheet_config,                              &
                     m_enh_shelf_config,                              &
                     choice_ice_margin_config,                        &
                     include_SSADIVA_crossterms_config,               &
                     do_GL_subgrid_friction_config,                   &
                     do_smooth_geometry_config,                       &
                     r_smooth_geometry_config,                        &
                     DIVA_visc_it_norm_dUV_tol_config,                &
                     DIVA_visc_it_nit_config,                         &
                     DIVA_visc_it_relax_config,                       &
                     DIVA_epsilon_sq_0_config,                        &
                     DIVA_visc_eff_min_config,                        &
                     DIVA_beta_max_config,                            &
                     DIVA_beta_min_thin_ice_config,                   &
                     DIVA_beta_min_thin_ice_H0_config,                &
                     DIVA_vel_max_config,                             &
                     DIVA_boundary_BC_u_west_config,                  &
                     DIVA_boundary_BC_u_east_config,                  &
                     DIVA_boundary_BC_u_south_config,                 &
                     DIVA_boundary_BC_u_north_config,                 &
                     DIVA_boundary_BC_v_west_config,                  &
                     DIVA_boundary_BC_v_east_config,                  &
                     DIVA_boundary_BC_v_south_config,                 &
                     DIVA_boundary_BC_v_north_config,                 &
                     DIVA_choice_matrix_solver_config,                &
                     DIVA_SOR_nit_config,                             &
                     DIVA_SOR_tol_config,                             &
                     DIVA_SOR_omega_config,                           &
                     DIVA_PETSc_rtol_config,                          &
                     DIVA_PETSc_abstol_config,                        &
                     do_read_velocities_from_restart_config,          &
                     choice_timestepping_config,                      &
                     choice_ice_integration_method_config,            &
                     dHi_choice_matrix_solver_config,                 &
                     dHi_SOR_nit_config,                              &
                     dHi_SOR_tol_config,                              &
                     dHi_SOR_omega_config,                            &
                     dHi_PETSc_rtol_config,                           &
                     dHi_PETSc_abstol_config,                         &
                     pc_epsilon_config,                               &
                     pc_k_I_config,                                   &
                     pc_k_p_config,                                   &
                     pc_eta_min_config,                               &
                     pc_max_timestep_iterations_config,               &
                     pc_redo_tol_config,                              &
                     ice_thickness_west_BC_config,                    &
                     ice_thickness_east_BC_config,                    &
                     ice_thickness_south_BC_config,                   &
                     ice_thickness_north_BC_config,                   &
                     choice_mask_noice_NAM_config,                    &
                     choice_mask_noice_EAS_config,                    &
                     choice_mask_noice_GRL_config,                    &
                     choice_mask_noice_ANT_config,                    &
                     fixed_sheet_geometry_config,                     &
                     fixed_shelf_geometry_config,                     &
                     fixed_grounding_line_g_config,                   &
                     fixed_grounding_line_f_config,                   &
                     fixed_decay_t_start_config,                      &
                     fixed_decay_t_end_config,                        &
                     relax_thick_t_start_config,                      &
                     relax_thick_t_end_config,                        &
                     do_target_dhdt_config,                           &
                     target_dhdt_filename_config,                     &
                     target_dhdt_t_end_config,                        &
                     choice_sliding_law_config,                       &
                     choice_idealised_sliding_law_config,             &
                     slid_delta_v_config,                             &
                     slid_Weertman_m_config,                          &
                     slid_Coulomb_reg_q_plastic_config,               &
                     slid_Coulomb_reg_u_threshold_config,             &
                     slid_ZI_ut_config,                               &
                     slid_ZI_p_config,                                &
                     do_slid_ZI_no_angle_config,                      &
                     include_basal_freezing_config,                   &
                     deltaT_basal_freezing_config,                    &
                     subgrid_friction_exponent_config,                &
                     choice_basal_hydrology_config,                   &
                     Martin2011_hydro_N_lim_config,                   &
                     Martin2011_hydro_Hb_min_config,                  &
                     Martin2011_hydro_Hb_max_config,                  &
                     choice_basal_roughness_config,                   &
                     uniform_Weertman_beta_sq_config,                 &
                     uniform_Coulomb_phi_fric_config,                 &
                     uniform_Tsai2015_alpha_sq_config,                &
                     uniform_Tsai2015_beta_sq_config,                 &
                     uniform_Schoof2005_alpha_sq_config,              &
                     uniform_Schoof2005_beta_sq_config,               &
                     choice_param_basal_roughness_config,             &
                     Martin2011till_phi_Hb_min_config,                &
                     Martin2011till_phi_Hb_max_config,                &
                     Martin2011till_phi_min_config,                   &
                     Martin2011till_phi_max_config,                   &
                     basal_roughness_filename_config,                 &
                     do_smooth_phi_restart_config,                    &
                     r_smooth_phi_restart_config,                     &
                     do_BIVgeo_config,                                &
                     BIVgeo_t_start_config,                           &
                     BIVgeo_t_end_config,                             &
                     choice_BIVgeo_method_config,                     &
                     BIVgeo_dt_config,                                &
                     BIVgeo_PDC2012_hinv_config,                      &
                     BIVgeo_Lipscomb2021_tauc_config,                 &
                     BIVgeo_Lipscomb2021_H0_config,                   &
                     BIVgeo_CISMplus_wH_config,                       &
                     BIVgeo_CISMplus_wu_config,                       &
                     BIVgeo_CISMplus_tauc_config,                     &
                     BIVgeo_CISMplus_H0_config,                       &
                     BIVgeo_CISMplus_u0_config,                       &
                     BIVgeo_Berends2022_tauc_config,                  &
                     BIVgeo_Berends2022_H0_config,                    &
                     BIVgeo_Berends2022_u0_config,                    &
                     BIVgeo_Berends2022_Hi_scale_config,              &
                     BIVgeo_Berends2022_u_scale_config,               &
                     BIVgeo_Berends2022_phimin_config,                &
                     BIVgeo_Berends2022_phimax_config,                &
                     BIVgeo_target_velocity_filename_config,          &
                     BIVgeo_filename_output_config,                   &
                     BIVgeo_Bernales2017_hinv_config,                 &
                     BIVgeo_Bernales2017_tol_diff_config,             &
                     BIVgeo_Bernales2017_tol_frac_config,             &
                     BIVgeo_Pien2023_H0_config,                       &
                     BIVgeo_Pien2023_tau_config,                      &
                     BIVgeo_Pien2023_r_config,                        &
                     BIVgeo_Pien2023_lowerHb_config,                  &
                     BIVgeo_Pien2023_upperHb_config,                  &
                     BIVgeo_Pien2023_min_config,                      &
                     BIVgeo_Pien2023_max_config,                      &
                     choice_calving_law_config,                       &
                     calving_threshold_thickness_config,              &
                     do_remove_shelves_config,                        &
                     remove_shelves_larger_than_PD_config,            &
                     continental_shelf_calving_config,                &
                     continental_shelf_min_height_config,             &
                     choice_initial_ice_temperature_config,           &
                     uniform_ice_temperature_config,                  &
                     choice_thermo_model_config,                      &
                     choice_ice_rheology_config,                      &
                     uniform_flow_factor_config,                      &
                     choice_ice_heat_capacity_config,                 &
                     uniform_ice_heat_capacity_config,                &
                     choice_ice_thermal_conductivity_config,          &
                     uniform_ice_thermal_conductivity_config,         &
                     choice_climate_model_config,                     &
                     choice_idealised_climate_config,                 &
                     direct_climate_foldername_NAM_config,            &
                     direct_climate_foldername_EAS_config,            &
                     direct_climate_foldername_GRL_config,            &
                     direct_climate_foldername_ANT_config,            &
                     direct_climate_basefilename_NAM_config,          &
                     direct_climate_basefilename_EAS_config,          &
                     direct_climate_basefilename_GRL_config,          &
                     direct_climate_basefilename_ANT_config,          &
                     do_direct_climate_geo_corr_config,               &
                     filename_PD_obs_climate_config,                  &
                     filename_climate_snapshot_PI_config,             &
                     filename_climate_snapshot_warm_config,           &
                     filename_climate_snapshot_cold_config,           &
                     reference_mask_method_config,                    &      
                     filename_snapshot_mask_PD_obs_config,            &
                     filename_snapshot_mask_GCM_PI_config,            &
                     filename_snapshot_mask_GCM_warm_config,          &
                     filename_snapshot_mask_GCM_cold_config,          &
                     constant_lapserate_config,                       &
                     climate_matrix_CO2vsice_NAM_config,              &
                     climate_matrix_CO2vsice_EAS_config,              &
                     climate_matrix_CO2vsice_GRL_config,              &
                     climate_matrix_CO2vsice_ANT_config,              &
                     matrix_high_CO2_level_config,                    &
                     matrix_low_CO2_level_config,                     &
                     matrix_warm_orbit_time_config,                   &
                     matrix_cold_orbit_time_config,                   &
                     climate_matrix_biascorrect_warm_config,          &
                     climate_matrix_biascorrect_cold_config,          &
                     switch_glacial_index_precip_config,              &
                     choice_ocean_model_config,                       &
                     choice_idealised_ocean_config,                   &
                     do_ocean_temperature_inversion_config,           &
                     ocean_temperature_inv_t_start_config,            &
                     ocean_temperature_inv_t_end_config,              &
                     ocean_temperature_inv_dT0_config,                &
                     ocean_temperature_inv_H0_config,                 &
                     ocean_temperature_inv_tau_config,                &
                     ocean_temperature_inv_min_config,                &
                     ocean_temperature_inv_max_config,                &
                     inverted_ocean_filename_output_config,           &
                     filename_PD_obs_ocean_config,                    &
                     name_ocean_temperature_obs_config,               &
                     name_ocean_salinity_obs_config,                  &
                     filename_GCM_ocean_snapshot_PI_config,           &
                     filename_GCM_ocean_snapshot_warm_config,         &
                     filename_GCM_ocean_snapshot_cold_config,         &
                     name_ocean_temperature_GCM_config,               &
                     name_ocean_salinity_GCM_config,                  &
                     ocean_temperature_PD_config,                     &
                     ocean_temperature_cold_config,                   &
                     ocean_temperature_warm_config,                   &
                     choice_ocean_vertical_grid_config,               &
                     ocean_vertical_grid_max_depth_config,            &
                     ocean_regular_grid_dz_config,                    &
                     do_ocean_extrap_config,                          &
                     ocean_extrap_dir_config,                         &
                     ocean_extrap_res_config,                         &
                     ocean_extrap_Gauss_sigma_config,                 &
                     ocean_extrap_hires_geo_filename_NAM_config,      &
                     ocean_extrap_hires_geo_filename_EAS_config,      &
                     ocean_extrap_hires_geo_filename_GRL_config,      &
                     ocean_extrap_hires_geo_filename_ANT_config,      &
                     ocean_w_tot_hist_averaging_window_config,        &
                     ocean_matrix_CO2vsice_NAM_config,                &
                     ocean_matrix_CO2vsice_EAS_config,                &
                     ocean_matrix_CO2vsice_GRL_config,                &
                     ocean_matrix_CO2vsice_ANT_config,                &
                     ocean_T_surf_per_basin_config,                   &
                     ocean_dT_dz_per_basin_config,                    &
                     ocean_filename_baseline_config,                  &
                     ocean_foldername_aTO_config,                     &
                     ocean_basefilename_aTO_config,                   &
                     ocean_foldername_aSO_config,                     &
                     ocean_basefilename_aSO_config,                   &
                     do_invert_linear_per_basin_config,               &
                     choice_SMB_model_config,                         &
                     choice_idealised_SMB_config,                     &
                     SMB_uniform_config,                              &
                     filename_direct_global_SMB_config,               &
                     filename_direct_regional_SMB_NAM_config,         &
                     filename_direct_regional_SMB_EAS_config,         &
                     filename_direct_regional_SMB_GRL_config,         &
                     filename_direct_regional_SMB_ANT_config,         &
                     filename_SMB_snapshot_NAM_config,                &
                     filename_SMB_snapshot_EAS_config,                &
                     filename_SMB_snapshot_GRL_config,                &
                     filename_SMB_snapshot_ANT_config,                &
                     SMB_IMAUITM_choice_init_firn_NAM_config,         &
                     SMB_IMAUITM_choice_init_firn_EAS_config,         &
                     SMB_IMAUITM_choice_init_firn_GRL_config,         &
                     SMB_IMAUITM_choice_init_firn_ANT_config,         &
                     SMB_IMAUITM_initial_firn_thickness_config,       &
                     SMB_IMAUITM_C_abl_constant_NAM_config,           &
                     SMB_IMAUITM_C_abl_constant_EAS_config,           &
                     SMB_IMAUITM_C_abl_constant_GRL_config,           &
                     SMB_IMAUITM_C_abl_constant_ANT_config,           &
                     SMB_IMAUITM_C_abl_Ts_NAM_config,                 &
                     SMB_IMAUITM_C_abl_Ts_EAS_config,                 &
                     SMB_IMAUITM_C_abl_Ts_GRL_config,                 &
                     SMB_IMAUITM_C_abl_Ts_ANT_config,                 &
                     SMB_IMAUITM_C_abl_Q_NAM_config,                  &
                     SMB_IMAUITM_C_abl_Q_EAS_config,                  &
                     SMB_IMAUITM_C_abl_Q_GRL_config,                  &
                     SMB_IMAUITM_C_abl_Q_ANT_config,                  &
                     SMB_IMAUITM_C_refr_NAM_config,                   &
                     SMB_IMAUITM_C_refr_EAS_config,                   &
                     SMB_IMAUITM_C_refr_GRL_config,                   &
                     SMB_IMAUITM_C_refr_ANT_config,                   &
                     ISMIP_forcing_filename_baseline_config,          &
                     ISMIP_forcing_foldername_aSMB_config,            &
                     ISMIP_forcing_basefilename_aSMB_config,          &
                     ISMIP_forcing_foldername_dSMBdz_config,          &
                     ISMIP_forcing_basefilename_dSMBdz_config,        &
                     ISMIP_forcing_foldername_aST_config,             &
                     ISMIP_forcing_basefilename_aST_config,           &
                     ISMIP_forcing_foldername_dSTdz_config,           &
                     ISMIP_forcing_basefilename_dSTdz_config,         &
                     choice_BMB_shelf_model_config,                   &
                     choice_idealised_BMB_shelf_config,               &
                     choice_BMB_sheet_model_config,                   &
                     BMB_shelf_uniform_config,                        &
                     BMB_sheet_uniform_config,                        &
                     choice_BMB_subgrid_config,                       &
                     do_asynchronous_BMB_config,                      &
                     BMB_max_config,                                  &
                     BMB_min_config,                                  &
                     choice_basin_scheme_NAM_config,                  &
                     choice_basin_scheme_EAS_config,                  &
                     choice_basin_scheme_GRL_config,                  &
                     choice_basin_scheme_ANT_config,                  &
                     filename_basins_NAM_config,                      &
                     filename_basins_EAS_config,                      &
                     filename_basins_GRL_config,                      &
                     filename_basins_ANT_config,                      &
                     do_merge_basins_ANT_config,                      &
                     do_merge_basins_GRL_config,                      &
                     choice_BMB_shelf_amplification_config,           &
                     basin_BMB_amplification_n_ANT_config,            &
                     basin_BMB_amplification_factor_ANT_config,       &
                     basin_BMB_amplification_n_GRL_config,            &
                     basin_BMB_amplification_factor_GRL_config,       &
                     BMB_Favier2019_lin_GammaT_config,                &
                     BMB_Favier2019_quad_GammaT_config,               &
                     BMB_Favier2019_Mplus_GammaT_config,              &
                     BMB_Lazeroms2018_GammaT_config,                  &
                     BMB_Lazeroms2018_find_GL_scheme_config,          &
                     BMB_PICO_nboxes_config,                          &
                     BMB_PICO_GammaTstar_config,                      &
                     filename_BMB_LADDIE_config,                      &
                     T_ocean_mean_PD_NAM_config,                      &
                     T_ocean_mean_PD_EAS_config,                      &
                     T_ocean_mean_PD_GRL_config,                      &
                     T_ocean_mean_PD_ANT_config,                      &
                     T_ocean_mean_cold_NAM_config,                    &
                     T_ocean_mean_cold_EAS_config,                    &
                     T_ocean_mean_cold_GRL_config,                    &
                     T_ocean_mean_cold_ANT_config,                    &
                     T_ocean_mean_warm_NAM_config,                    &
                     T_ocean_mean_warm_EAS_config,                    &
                     T_ocean_mean_warm_GRL_config,                    &
                     T_ocean_mean_warm_ANT_config,                    &
                     BMB_deepocean_PD_NAM_config,                     &
                     BMB_deepocean_PD_EAS_config,                     &
                     BMB_deepocean_PD_GRL_config,                     &
                     BMB_deepocean_PD_ANT_config,                     &
                     BMB_deepocean_cold_NAM_config,                   &
                     BMB_deepocean_cold_EAS_config,                   &
                     BMB_deepocean_cold_GRL_config,                   &
                     BMB_deepocean_cold_ANT_config,                   &
                     BMB_deepocean_warm_NAM_config,                   &
                     BMB_deepocean_warm_EAS_config,                   &
                     BMB_deepocean_warm_GRL_config,                   &
                     BMB_deepocean_warm_ANT_config,                   &
                     BMB_shelf_exposed_PD_NAM_config,                 &
                     BMB_shelf_exposed_PD_EAS_config,                 &
                     BMB_shelf_exposed_PD_GRL_config,                 &
                     BMB_shelf_exposed_PD_ANT_config,                 &
                     BMB_shelf_exposed_cold_NAM_config,               &
                     BMB_shelf_exposed_cold_EAS_config,               &
                     BMB_shelf_exposed_cold_GRL_config,               &
                     BMB_shelf_exposed_cold_ANT_config,               &
                     BMB_shelf_exposed_warm_NAM_config,               &
                     BMB_shelf_exposed_warm_EAS_config,               &
                     BMB_shelf_exposed_warm_GRL_config,               &
                     BMB_shelf_exposed_warm_ANT_config,               &
                     subshelf_melt_factor_NAM_config,                 &
                     subshelf_melt_factor_EAS_config,                 &
                     subshelf_melt_factor_GRL_config,                 &
                     subshelf_melt_factor_ANT_config,                 &
                     deep_ocean_threshold_depth_NAM_config,           &
                     deep_ocean_threshold_depth_EAS_config,           &
                     deep_ocean_threshold_depth_GRL_config,           &
                     deep_ocean_threshold_depth_ANT_config,           &
                     choice_ice_isotopes_model_config,                &
                     uniform_ice_d18O_config,                         &
                     do_ocean_floodfill_config,                       &
                     choice_sealevel_model_config,                    &
                     fixed_sealevel_config,                           &
                     filename_sealevel_record_config,                 &
                     sealevel_record_length_config,                   &
                     choice_GIA_model_config,                         &
                     dx_GIA_config,                                   &
                     ELRA_lithosphere_flex_rigidity_config,           &
                     ELRA_bedrock_relaxation_time_config,             &
                     ELRA_mantle_density_config,                      &
                     SELEN_run_at_t_start_config,                     &
                     SELEN_n_TDOF_iterations_config,                  &
                     SELEN_n_recursion_iterations_config,             &
                     SELEN_use_rotational_feedback_config,            &
                     SELEN_n_harmonics_config,                        &
                     SELEN_display_progress_config,                   &
                     SELEN_dir_config,                                &
                     SELEN_global_topo_filename_config,               &
                     SELEN_TABOO_init_filename_config,                &
                     SELEN_LMJ_VALUES_filename_config,                &
                     SELEN_irreg_time_n_config,                       &
                     SELEN_irreg_time_window_config,                  &
                     SELEN_lith_thickness_config,                     &
                     SELEN_visc_n_config,                             &
                     SELEN_visc_prof_config,                          &
                     SELEN_TABOO_CDE_config,                          &
                     SELEN_TABOO_TLOVE_config,                        &
                     SELEN_TABOO_DEG1_config,                         &
                     SELEN_TABOO_RCMB_config,                         &
                     help_field_01_config,                            &
                     help_field_02_config,                            &
                     help_field_03_config,                            &
                     help_field_04_config,                            &
                     help_field_05_config,                            &
                     help_field_06_config,                            &
                     help_field_07_config,                            &
                     help_field_08_config,                            &
                     help_field_09_config,                            &
                     help_field_10_config,                            &
                     help_field_11_config,                            &
                     help_field_12_config,                            &
                     help_field_13_config,                            &
                     help_field_14_config,                            &
                     help_field_15_config,                            &
                     help_field_16_config,                            &
                     help_field_17_config,                            &
                     help_field_18_config,                            &
                     help_field_19_config,                            &
                     help_field_20_config,                            &
                     help_field_21_config,                            &
                     help_field_22_config,                            &
                     help_field_23_config,                            &
                     help_field_24_config,                            &
                     help_field_25_config,                            &
                     help_field_26_config,                            &
                     help_field_27_config,                            &
                     help_field_28_config,                            &
                     help_field_29_config,                            &
                     help_field_30_config,                            &
                     help_field_31_config,                            &
                     help_field_32_config,                            &
                     help_field_33_config,                            &
                     help_field_34_config,                            &
                     help_field_35_config,                            &
                     help_field_36_config,                            &
                     help_field_37_config,                            &
                     help_field_38_config,                            &
                     help_field_39_config,                            &
                     help_field_40_config,                            &
                     help_field_41_config,                            &
                     help_field_42_config,                            &
                     help_field_43_config,                            &
                     help_field_44_config,                            &
                     help_field_45_config,                            &
                     help_field_46_config,                            &
                     help_field_47_config,                            &
                     help_field_48_config,                            &
                     help_field_49_config,                            &
                     help_field_50_config

    IF (config_filename == '') RETURN

    ! Write the CONFIG namelist to a temporary file
    namelist_filename = TRIM( config_filename) // '_namelist_temp.txt'
    OPEN(  UNIT = namelist_unit, FILE = TRIM( namelist_filename))
    WRITE( UNIT = namelist_unit, NML  = CONFIG)
    CLOSE( UNIT = namelist_unit)

    ! Check the config file for validity
    CALL check_config_file_validity( config_filename, namelist_filename)

    ! Delete the temporary CONFIG namelist file
    CALL system('rm -f ' // TRIM( namelist_filename))

    ! Open the config file
    OPEN(  UNIT = config_unit, FILE = TRIM( config_filename), STATUS = 'OLD', ACTION = 'READ', IOSTAT = ios)
    IF (ios /= 0) THEN
      WRITE(0,'(A,A,A)') colour_string('ERROR: config file "' // TRIM( config_filename),'red') // '" not found!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Read the config file using the CONFIG namelist
    READ(  UNIT = config_unit, NML = CONFIG, IOSTAT = ios)
    IF (ios /= 0) THEN
      WRITE(0,'(A,A,A)') colour_string('ERROR: error while reading config file "' // TRIM( config_filename),'red') // '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Close the config file
    CLOSE( UNIT = config_unit)

  END SUBROUTINE read_config_file

  SUBROUTINE check_config_file_validity( config_filename, namelist_filename)
    ! Check if the provided config file is valid
    !
    ! Do this by reading one line at a time of the config file, determining the name of the variable
    ! declared in that line, and checking if that variable also exists in the namelist file
    !
    ! Assumes that the CONFIG namelist has already been written to the specified file.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),INTENT(IN) :: config_filename, namelist_filename

    ! Local variables:
    INTEGER, PARAMETER            :: config_unit   = 1337
    INTEGER, PARAMETER            :: namelist_unit = 1338
    INTEGER                       :: ios, ierr, cerr
    LOGICAL                       :: found_end_of_file_config, found_end_of_file_namelist
    CHARACTER(256)                :: single_line_config      , single_line_namelist
    INTEGER                       :: line_counter_config     , line_counter_namelist
    LOGICAL                       :: found_match, found_mismatch

    ! Open the config and namelist files
    OPEN( UNIT = config_unit, FILE = config_filename, IOSTAT = ios)
    IF (ios /= 0) THEN
      WRITE(0,'(A)') colour_string('ERROR','red') // ': config file "' // TRIM( config_filename) // '" not found!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Read one line at a time of the config file, determine the name of the variable
    ! declared in that line, and check if that variable also exists in the namelist file

    found_end_of_file_config = .FALSE.
    line_counter_config      = 0
    found_mismatch           = .FALSE.

    DO WHILE (.NOT. found_end_of_file_config)

      line_counter_config = line_counter_config + 1

      ! Read a single line from the config file
      READ( UNIT = config_unit, FMT = '(A)', IOSTAT = ios) single_line_config

      ! If we've reached the end of the file before finding the terminating forward slash, this config file is not valid.
      IF (ios < 0) THEN
        WRITE(0,'(A)') colour_string('ERROR','red') // ': config file "' // TRIM( config_filename) // '" is not terminated with a forward slash!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

      ! Remove all leading spaces
      CALL remove_leading_spaces( single_line_config)

      ! The variable name is the part of the string left of the first (, =, or space.
      single_line_config = single_line_config( 1:SCAN( single_line_config, '( =')-1)

      ! Get config variable in all caps for case-insensitive comparison
      CALL capitalise_string( single_line_config)

      ! The forward slash at the end terminates the config file
      IF (single_line_config == '/') THEN
        found_end_of_file_config = .TRUE.
      END IF

      ! Disregard empty lines, commented lines, and the header line
      IF (single_line_config == '' .OR. single_line_config == '&CONFIG' .OR. single_line_config( 1:1) == '!') THEN
        CYCLE
      END IF

      ! Open the namelist file
      OPEN( UNIT = namelist_unit, FILE = namelist_filename)
      IF (ios /= 0) THEN
        WRITE(0,'(A)') colour_string('ERROR','red') // ': namelist file "' // TRIM( namelist_filename) // '" not found!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

      ! Read all variables from the namelist file and check if any of them match the current config variable

      found_end_of_file_namelist = .FALSE.
      line_counter_namelist      = 0
      found_match                = .FALSE.

      DO WHILE ((.NOT. found_end_of_file_namelist) .AND. (.NOT. found_match))

        line_counter_namelist = line_counter_namelist + 1

        ! Read a single line from the namelist file
        READ( UNIT = namelist_unit, FMT = '(A)', IOSTAT = ios) single_line_namelist

        ! If we've reached the end of the file before finding the terminating forward slash, this namelist file is not valid.
        IF (ios < 0) THEN
          WRITE(0,'(A)') colour_string('ERROR','red') // ': namelist file "' // TRIM( namelist_filename) // '" is not terminated with a forward slash!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF

        ! Remove all leading spaces
        CALL remove_leading_spaces( single_line_namelist)

        ! The variable name is the part of the string left of the first (, =, or space.
        single_line_namelist = single_line_namelist( 1:SCAN( single_line_namelist, '( =')-1)

        ! Get namelist variable in all caps for case-insensitive comparison
        CALL capitalise_string( single_line_namelist)

        ! The forward slash at the end terminates the config file
        IF (single_line_namelist == '/') THEN
          found_end_of_file_namelist = .TRUE.
        END IF

        ! Disregard empty lines, commented lines, and the header line
        IF (single_line_namelist == '' .OR. single_line_namelist == '&CONFIG' .OR. single_line_namelist( 1:1) == '!') THEN
          CYCLE
        END IF

        ! Check if this namelist variable matches the config variable
        IF (single_line_namelist == single_line_config) THEN
          found_match = .TRUE.
        END IF

      END DO ! DO WHILE ((.NOT. found_end_of_file_namelist) .AND. (.NOT. found_match))

      ! If no matching variable was found in the namelist file, print an error
      IF (.NOT. found_match) THEN
        WRITE(0,'(A,I4)') colour_string('ERROR','red') // ': invalid config variable "' // TRIM( single_line_config) // &
          '" in file "' // TRIM( config_filename) // '", line ', line_counter_config
        found_mismatch = .TRUE.
      END IF

      ! Close the namelist file
      CLOSE( UNIT = namelist_unit)

    END DO ! DO WHILE (.NOT. found_end_of_file_config)

    ! Close the config file
    CLOSE( UNIT = config_unit)

    ! If an invalid config variable was found, crash.
    IF (found_mismatch) CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

  END SUBROUTINE check_config_file_validity

  SUBROUTINE copy_variables_to_struct
    ! Overwrite the values in the fields of the "C" type with the values
    ! of the "_config" variables, some which by now have had their default
    ! values overwritten by the values specified in the external config file.

    IMPLICIT NONE

   ! Time steps and range
   !=====================

    C%start_time_of_run                        = start_time_of_run_config
    C%end_time_of_run                          = end_time_of_run_config
    C%dt_coupling                              = dt_coupling_config
    C%dt_max                                   = dt_max_config
    C%dt_min                                   = dt_min_config
    C%dt_startup_phase                         = dt_startup_phase_config
    C%dt_thermo                                = dt_thermo_config
    C%dt_climate                               = dt_climate_config
    C%dt_ocean                                 = dt_ocean_config
    C%dt_SMB                                   = dt_SMB_config
    C%dt_BMB                                   = dt_BMB_config
    C%dt_bedrock_ELRA                          = dt_bedrock_ELRA_config
    C%dt_SELEN                                 = dt_SELEN_config
    C%dt_output                                = dt_output_config
    C%dt_output_restart                        = dt_output_restart_config

    ! Which ice sheets do we simulate?
    ! ================================

    C%do_NAM                                   = do_NAM_config
    C%do_EAS                                   = do_EAS_config
    C%do_GRL                                   = do_GRL_config
    C%do_ANT                                   = do_ANT_config

    ! Benchmark experiments
    ! =====================

    ! SSA_icestream (see Schoof 2006, and also Bueler and Brown 2009)
    C%SSA_icestream_A                          = SSA_icestream_A_config
    C%SSA_icestream_L                          = SSA_icestream_L_config
    C%SSA_icestream_m                          = SSA_icestream_m_config
    C%SSA_icestream_tantheta                   = SSA_icestream_tantheta_config
    C%SSA_icestream_H                          = SSA_icestream_H_config

    ! ISMIP-HOM (see Pattyn et al. 2008)
    C%ISMIP_HOM_L                              = ISMIP_HOM_L_config
    C%ISMIP_HOM_E_Arolla_filename              = ISMIP_HOM_E_Arolla_filename_config

    ! MISMIP+ (see Asay-Davis et al., 2016)
    C%MISMIPplus_do_tune_A_for_GL              = MISMIPplus_do_tune_A_for_GL_config
    C%MISMIPplus_xGL_target                    = MISMIPplus_xGL_target_config
    C%MISMIPplus_A_flow_initial                = MISMIPplus_A_flow_initial_config
    C%MISMIPplus_scenario                      = MISMIPplus_scenario_config

    ! MISOMIP1 (see Asay-Davis et al., 2016)
    C%MISOMIP1_scenario                        = MISOMIP1_scenario_config

    ! Whether or not to let IMAU_ICE dynamically create its own output folder
    ! =======================================================================

    C%create_procedural_output_dir             = create_procedural_output_dir_config
    C%fixed_output_dir                         = fixed_output_dir_config
    C%fixed_output_dir_suffix                  = fixed_output_dir_suffix_config
    C%do_write_regional_scalar_output          = do_write_regional_scalar_output_config
    C%do_write_global_scalar_output            = do_write_global_scalar_output_config
    C%do_write_regional_scalar_every_timestep  = do_write_regional_scalar_every_timestep_config
    ! Debugging
    ! =========

    C%do_check_for_NaN                         = do_check_for_NaN_config
    C%do_time_display                          = do_time_display_config

    ! ISMIP-style output
    ! ==================

    C%do_write_ISMIP_output                    = do_write_ISMIP_output_config
    C%ISMIP_output_group_code                  = ISMIP_output_group_code_config
    C%ISMIP_output_model_code                  = ISMIP_output_model_code_config
    C%ISMIP_output_experiment_code             = ISMIP_output_experiment_code_config
    C%ISMIP_output_basetime                    = ISMIP_output_basetime_config

    ! Grids
    ! =====

    ! North America
    C%lambda_M_NAM                             = lambda_M_NAM_config
    C%phi_M_NAM                                = phi_M_NAM_config
    C%beta_stereo_NAM                          = beta_stereo_NAM_config
    C%xmin_NAM                                 = xmin_NAM_config
    C%xmax_NAM                                 = xmax_NAM_config
    C%ymin_NAM                                 = ymin_NAM_config
    C%ymax_NAM                                 = ymax_NAM_config
    C%dx_NAM                                   = dx_NAM_config

    ! Eurasia
    C%lambda_M_EAS                             = lambda_M_EAS_config
    C%phi_M_EAS                                = phi_M_EAS_config
    C%beta_stereo_EAS                          = beta_stereo_EAS_config
    C%xmin_EAS                                 = xmin_EAS_config
    C%xmax_EAS                                 = xmax_EAS_config
    C%ymin_EAS                                 = ymin_EAS_config
    C%ymax_EAS                                 = ymax_EAS_config
    C%dx_EAS                                   = dx_EAS_config

    ! Greenland
    C%lambda_M_GRL                             = lambda_M_GRL_config
    C%phi_M_GRL                                = phi_M_GRL_config
    C%beta_stereo_GRL                          = beta_stereo_GRL_config
    C%xmin_GRL                                 = xmin_GRL_config
    C%xmax_GRL                                 = xmax_GRL_config
    C%ymin_GRL                                 = ymin_GRL_config
    C%ymax_GRL                                 = ymax_GRL_config
    C%dx_GRL                                   = dx_GRL_config

    ! Antarctica
    C%lambda_M_ANT                             = lambda_M_ANT_config
    C%phi_M_ANT                                = phi_M_ANT_config
    C%beta_stereo_ANT                          = beta_stereo_ANT_config
    C%xmin_ANT                                 = xmin_ANT_config
    C%xmax_ANT                                 = xmax_ANT_config
    C%ymin_ANT                                 = ymin_ANT_config
    C%ymax_ANT                                 = ymax_ANT_config
    C%dx_ANT                                   = dx_ANT_config

    ! Scaled vertical coordinate zeta
    ! ===============================

    C%nz                                       = nz_config
    ALLOCATE( C%zeta( C%nz))
    C%zeta                                     = zeta_config( 1:C%nz)

    ! Reference geometries (initial, present-day, and GIA equilibrium)
    ! ================================================================

    ! Initial geometry
    C%choice_refgeo_init_NAM                   = choice_refgeo_init_NAM_config
    C%choice_refgeo_init_EAS                   = choice_refgeo_init_EAS_config
    C%choice_refgeo_init_GRL                   = choice_refgeo_init_GRL_config
    C%choice_refgeo_init_ANT                   = choice_refgeo_init_ANT_config
    C%time_to_restart_from_NAM                 = time_to_restart_from_NAM_config
    C%time_to_restart_from_EAS                 = time_to_restart_from_EAS_config
    C%time_to_restart_from_GRL                 = time_to_restart_from_GRL_config
    C%time_to_restart_from_ANT                 = time_to_restart_from_ANT_config
    C%choice_refgeo_init_idealised             = choice_refgeo_init_idealised_config
    C%filename_refgeo_init_NAM                 = filename_refgeo_init_NAM_config
    C%filename_refgeo_init_EAS                 = filename_refgeo_init_EAS_config
    C%filename_refgeo_init_GRL                 = filename_refgeo_init_GRL_config
    C%filename_refgeo_init_ANT                 = filename_refgeo_init_ANT_config

    ! Present-day geometry
    C%choice_refgeo_PD_NAM                     = choice_refgeo_PD_NAM_config
    C%choice_refgeo_PD_EAS                     = choice_refgeo_PD_EAS_config
    C%choice_refgeo_PD_GRL                     = choice_refgeo_PD_GRL_config
    C%choice_refgeo_PD_ANT                     = choice_refgeo_PD_ANT_config
    C%choice_refgeo_PD_idealised               = choice_refgeo_PD_idealised_config
    C%filename_refgeo_PD_NAM                   = filename_refgeo_PD_NAM_config
    C%filename_refgeo_PD_EAS                   = filename_refgeo_PD_EAS_config
    C%filename_refgeo_PD_GRL                   = filename_refgeo_PD_GRL_config
    C%filename_refgeo_PD_ANT                   = filename_refgeo_PD_ANT_config

    ! GIA equilibrium geometry
    C%choice_refgeo_GIAeq_NAM                  = choice_refgeo_GIAeq_NAM_config
    C%choice_refgeo_GIAeq_EAS                  = choice_refgeo_GIAeq_EAS_config
    C%choice_refgeo_GIAeq_GRL                  = choice_refgeo_GIAeq_GRL_config
    C%choice_refgeo_GIAeq_ANT                  = choice_refgeo_GIAeq_ANT_config
    C%choice_refgeo_GIAeq_idealised            = choice_refgeo_GIAeq_idealised_config
    C%filename_refgeo_GIAeq_NAM                = filename_refgeo_GIAeq_NAM_config
    C%filename_refgeo_GIAeq_EAS                = filename_refgeo_GIAeq_EAS_config
    C%filename_refgeo_GIAeq_GRL                = filename_refgeo_GIAeq_GRL_config
    C%filename_refgeo_GIAeq_ANT                = filename_refgeo_GIAeq_ANT_config

    C%remove_Lake_Vostok                       = remove_Lake_Vostok_config

    ! Global forcing (insolation, CO2, d18O, geothermal heat flux)
    ! ============================================================

    C%choice_forcing_method                    = choice_forcing_method_config

    ! Insolation forcing (NetCDF)
    C%choice_insolation_forcing                = choice_insolation_forcing_config
    C%static_insolation_time                   = static_insolation_time_config
    C%filename_insolation                      = filename_insolation_config

    ! CO2 record (ASCII text file, so the number of rows needs to be specified)
    C%filename_CO2_record                      = filename_CO2_record_config
    C%CO2_record_length                        = CO2_record_length_config

    ! d18O record (ASCII text file, so the number of rows needs to be specified)
    C%filename_d18O_record                     = filename_d18O_record_config
    C%d18O_record_length                       = d18O_record_length_config

    ! Parameters for combining insolation and CO2 in the matrix method
    C%do_combine_CO2_and_insolation            = do_combine_CO2_and_insolation_config
    
    C%insolation_weigth_mean_NAM               = insolation_weigth_mean_NAM_config
    C%insolation_weigth_amplitude_NAM          = insolation_weigth_amplitude_NAM_config
    C%insolation_weigth_mean_EAS               = insolation_weigth_mean_EAS_config
    C%insolation_weigth_amplitude_EAS          = insolation_weigth_amplitude_EAS_config
    C%insolation_weigth_mean_GRL               = insolation_weigth_mean_GRL_config
    C%insolation_weigth_amplitude_GRL          = insolation_weigth_amplitude_GRL_config
    C%insolation_weigth_mean_ANT               = insolation_weigth_mean_ANT_config
    C%insolation_weigth_amplitude_ANT          = insolation_weigth_amplitude_ANT_config

    ! Climate matrix wind
    C%do_climate_matrix_wind                   = do_climate_matrix_wind_config

    ! Geothermal heat flux
    C%choice_geothermal_heat_flux              = choice_geothermal_heat_flux_config
    C%constant_geothermal_heat_flux            = constant_geothermal_heat_flux_config
    C%filename_geothermal_heat_flux            = filename_geothermal_heat_flux_config

    ! Parameters for calculating modelled benthic d18O
    C%do_calculate_benthic_d18O                = do_calculate_benthic_d18O_config
    C%dT_deepwater_averaging_window            = dT_deepwater_averaging_window_config
    C%dT_deepwater_dT_surf_ratio               = dT_deepwater_dT_surf_ratio_config
    C%d18O_dT_deepwater_ratio                  = d18O_dT_deepwater_ratio_config

    ! Parameters for the inverse routine
    C%dT_glob_inverse_averaging_window         = dT_glob_inverse_averaging_window_config
    C%inverse_d18O_to_dT_glob_scaling          = inverse_d18O_to_dT_glob_scaling_config
    C%CO2_inverse_averaging_window             = CO2_inverse_averaging_window_config
    C%inverse_d18O_to_CO2_scaling              = inverse_d18O_to_CO2_scaling_config
    C%inverse_d18O_to_CO2_initial_CO2          = inverse_d18O_to_CO2_initial_CO2_config

    ! Intitial inverse d18O history (from restart file netcdf)
    C%choice_d18O_inverse_init_NAM             = choice_d18O_inverse_init_NAM_config
    C%choice_d18O_inverse_init_EAS             = choice_d18O_inverse_init_EAS_config
    C%choice_d18O_inverse_init_GRL             = choice_d18O_inverse_init_GRL_config
    C%choice_d18O_inverse_init_ANT             = choice_d18O_inverse_init_ANT_config

    C%filename_d18O_inverse_init_NAM           = filename_d18O_inverse_init_NAM_config
    C%filename_d18O_inverse_init_EAS           = filename_d18O_inverse_init_EAS_config
    C%filename_d18O_inverse_init_GRL           = filename_d18O_inverse_init_GRL_config
    C%filename_d18O_inverse_init_ANT           = filename_d18O_inverse_init_ANT_config

    ! Ice dynamics - velocity
    ! =======================

    C%choice_ice_dynamics                      = choice_ice_dynamics_config
    C%n_flow                                   = n_flow_config
    C%m_enh_sheet                              = m_enh_sheet_config
    C%m_enh_shelf                              = m_enh_shelf_config
    C%choice_ice_margin                        = choice_ice_margin_config
    C%include_SSADIVA_crossterms               = include_SSADIVA_crossterms_config
    C%do_GL_subgrid_friction                   = do_GL_subgrid_friction_config
    C%do_smooth_geometry                       = do_smooth_geometry_config
    C%r_smooth_geometry                        = r_smooth_geometry_config

    ! Some parameters for numerically solving the SSA/DIVA
    C%DIVA_visc_it_norm_dUV_tol                = DIVA_visc_it_norm_dUV_tol_config
    C%DIVA_visc_it_nit                         = DIVA_visc_it_nit_config
    C%DIVA_visc_it_relax                       = DIVA_visc_it_relax_config
    C%DIVA_epsilon_sq_0                        = DIVA_epsilon_sq_0_config
    C%DIVA_visc_eff_min                        = DIVA_visc_eff_min_config
    C%DIVA_beta_max                            = DIVA_beta_max_config
    C%DIVA_beta_min_thin_ice                   = DIVA_beta_min_thin_ice_config
    C%DIVA_beta_min_thin_ice_H0                = DIVA_beta_min_thin_ice_H0_config
    C%DIVA_vel_max                             = DIVA_vel_max_config
    C%DIVA_boundary_BC_u_west                  = DIVA_boundary_BC_u_west_config
    C%DIVA_boundary_BC_u_east                  = DIVA_boundary_BC_u_east_config
    C%DIVA_boundary_BC_u_south                 = DIVA_boundary_BC_u_south_config
    C%DIVA_boundary_BC_u_north                 = DIVA_boundary_BC_u_north_config
    C%DIVA_boundary_BC_v_west                  = DIVA_boundary_BC_v_west_config
    C%DIVA_boundary_BC_v_east                  = DIVA_boundary_BC_v_east_config
    C%DIVA_boundary_BC_v_south                 = DIVA_boundary_BC_v_south_config
    C%DIVA_boundary_BC_v_north                 = DIVA_boundary_BC_v_north_config
    C%DIVA_choice_matrix_solver                = DIVA_choice_matrix_solver_config
    C%DIVA_SOR_nit                             = DIVA_SOR_nit_config
    C%DIVA_SOR_tol                             = DIVA_SOR_tol_config
    C%DIVA_SOR_omega                           = DIVA_SOR_omega_config
    C%DIVA_PETSc_rtol                          = DIVA_PETSc_rtol_config
    C%DIVA_PETSc_abstol                        = DIVA_PETSc_abstol_config
    C%do_read_velocities_from_restart          = do_read_velocities_from_restart_config

    ! Ice dynamics - time integration
    ! ===============================

    C%choice_timestepping                      = choice_timestepping_config
    C%choice_ice_integration_method            = choice_ice_integration_method_config
    C%dHi_choice_matrix_solver                 = dHi_choice_matrix_solver_config
    C%dHi_SOR_nit                              = dHi_SOR_nit_config
    C%dHi_SOR_tol                              = dHi_SOR_tol_config
    C%dHi_SOR_omega                            = dHi_SOR_omega_config
    C%dHi_PETSc_rtol                           = dHi_PETSc_rtol_config
    C%dHi_PETSc_abstol                         = dHi_PETSc_abstol_config

    ! Predictor-corrector ice-thickness update
    C%pc_epsilon                               = pc_epsilon_config
    C%pc_k_I                                   = pc_k_I_config
    C%pc_k_p                                   = pc_k_p_config
    C%pc_eta_min                               = pc_eta_min_config
    C%pc_max_timestep_iterations               = pc_max_timestep_iterations_config
    C%pc_redo_tol                              = pc_redo_tol_config

    ! Ice thickness boundary conditions
    C%ice_thickness_west_BC                    = ice_thickness_west_BC_config
    C%ice_thickness_east_BC                    = ice_thickness_east_BC_config
    C%ice_thickness_south_BC                   = ice_thickness_south_BC_config
    C%ice_thickness_north_BC                   = ice_thickness_north_BC_config
    C%choice_mask_noice_NAM                    = choice_mask_noice_NAM_config
    C%choice_mask_noice_EAS                    = choice_mask_noice_EAS_config
    C%choice_mask_noice_GRL                    = choice_mask_noice_GRL_config
    C%choice_mask_noice_ANT                    = choice_mask_noice_ANT_config

    ! Fixed/delayed ice thickness evolution
    C%fixed_sheet_geometry                     = fixed_sheet_geometry_config
    C%fixed_shelf_geometry                     = fixed_shelf_geometry_config
    C%fixed_grounding_line_g                   = fixed_grounding_line_g_config
    C%fixed_grounding_line_f                   = fixed_grounding_line_f_config
    C%fixed_decay_t_start                      = fixed_decay_t_start_config
    C%fixed_decay_t_end                        = fixed_decay_t_end_config
    C%relax_thick_t_start                      = relax_thick_t_start_config
    C%relax_thick_t_end                        = relax_thick_t_end_config

    ! Target dHi_dt during model spinup
    C%do_target_dhdt                           = do_target_dhdt_config
    C%target_dhdt_filename                     = target_dhdt_filename_config
    C%target_dhdt_t_end                        = target_dhdt_t_end_config

    ! Ice dynamics - basal conditions and sliding
    ! ===========================================

    ! Sliding laws
    C%choice_sliding_law                       = choice_sliding_law_config
    C%choice_idealised_sliding_law             = choice_idealised_sliding_law_config
    C%slid_delta_v                             = slid_delta_v_config
    C%slid_Weertman_m                          = slid_Weertman_m_config
    C%slid_Coulomb_reg_q_plastic               = slid_Coulomb_reg_q_plastic_config
    C%slid_Coulomb_reg_u_threshold             = slid_Coulomb_reg_u_threshold_config
    C%slid_ZI_ut                               = slid_ZI_ut_config
    C%slid_ZI_p                                = slid_ZI_p_config
    C%do_slid_ZI_no_angle                      = do_slid_ZI_no_angle_config
    C%include_basal_freezing                   = include_basal_freezing_config
    C%deltaT_basal_freezing                    = deltaT_basal_freezing_config
    C%subgrid_friction_exponent                = subgrid_friction_exponent_config

    ! Basal hydrology
    C%choice_basal_hydrology                   = choice_basal_hydrology_config
    C%Martin2011_hydro_N_lim                   = Martin2011_hydro_N_lim_config
    C%Martin2011_hydro_Hb_min                  = Martin2011_hydro_Hb_min_config
    C%Martin2011_hydro_Hb_max                  = Martin2011_hydro_Hb_max_config

    ! Basal roughness / friction
    C%choice_basal_roughness                   = choice_basal_roughness_config
    C%uniform_Weertman_beta_sq                 = uniform_Weertman_beta_sq_config
    C%uniform_Coulomb_phi_fric                 = uniform_Coulomb_phi_fric_config
    C%uniform_Tsai2015_alpha_sq                = uniform_Tsai2015_alpha_sq_config
    C%uniform_Tsai2015_beta_sq                 = uniform_Tsai2015_beta_sq_config
    C%uniform_Schoof2005_alpha_sq              = uniform_Schoof2005_alpha_sq_config
    C%uniform_Schoof2005_beta_sq               = uniform_Schoof2005_beta_sq_config
    C%choice_param_basal_roughness             = choice_param_basal_roughness_config
    C%Martin2011till_phi_Hb_min                = Martin2011till_phi_Hb_min_config
    C%Martin2011till_phi_Hb_max                = Martin2011till_phi_Hb_max_config
    C%Martin2011till_phi_min                   = Martin2011till_phi_min_config
    C%Martin2011till_phi_max                   = Martin2011till_phi_max_config
    C%basal_roughness_filename                 = basal_roughness_filename_config
    C%do_smooth_phi_restart                    = do_smooth_phi_restart_config
    C%r_smooth_phi_restart                     = r_smooth_phi_restart_config

    ! Basal inversion
    C%do_BIVgeo                                = do_BIVgeo_config
    C%BIVgeo_t_start                           = BIVgeo_t_start_config
    C%BIVgeo_t_end                             = BIVgeo_t_end_config
    C%choice_BIVgeo_method                     = choice_BIVgeo_method_config
    C%BIVgeo_dt                                = BIVgeo_dt_config
    C%BIVgeo_PDC2012_hinv                      = BIVgeo_PDC2012_hinv_config
    C%BIVgeo_Lipscomb2021_tauc                 = BIVgeo_Lipscomb2021_tauc_config
    C%BIVgeo_Lipscomb2021_H0                   = BIVgeo_Lipscomb2021_H0_config
    C%BIVgeo_CISMplus_wH                       = BIVgeo_CISMplus_wH_config
    C%BIVgeo_CISMplus_wu                       = BIVgeo_CISMplus_wu_config
    C%BIVgeo_CISMplus_tauc                     = BIVgeo_CISMplus_tauc_config
    C%BIVgeo_CISMplus_H0                       = BIVgeo_CISMplus_H0_config
    C%BIVgeo_CISMplus_u0                       = BIVgeo_CISMplus_u0_config
    C%BIVgeo_Berends2022_tauc                  = BIVgeo_Berends2022_tauc_config
    C%BIVgeo_Berends2022_H0                    = BIVgeo_Berends2022_H0_config
    C%BIVgeo_Berends2022_u0                    = BIVgeo_Berends2022_u0_config
    C%BIVgeo_Berends2022_Hi_scale              = BIVgeo_Berends2022_Hi_scale_config
    C%BIVgeo_Berends2022_u_scale               = BIVgeo_Berends2022_u_scale_config
    C%BIVgeo_Berends2022_phimin                = BIVgeo_Berends2022_phimin_config
    C%BIVgeo_Berends2022_phimax                = BIVgeo_Berends2022_phimax_config
    C%BIVgeo_target_velocity_filename          = BIVgeo_target_velocity_filename_config
    C%BIVgeo_filename_output                   = BIVgeo_filename_output_config
    C%BIVgeo_Bernales2017_hinv                 = BIVgeo_Bernales2017_hinv_config
    C%BIVgeo_Bernales2017_tol_diff             = BIVgeo_Bernales2017_tol_diff_config
    C%BIVgeo_Bernales2017_tol_frac             = BIVgeo_Bernales2017_tol_frac_config
    C%BIVgeo_Pien2023_H0                       = BIVgeo_Pien2023_H0_config
    C%BIVgeo_Pien2023_tau                      = BIVgeo_Pien2023_tau_config
    C%BIVgeo_Pien2023_r                        = BIVgeo_Pien2023_r_config
    C%BIVgeo_Pien2023_lowerHb                  = BIVgeo_Pien2023_lowerHb_config
    C%BIVgeo_Pien2023_upperHb                  = BIVgeo_Pien2023_upperHb_config
    C%BIVgeo_Pien2023_min                      = BIVgeo_Pien2023_min_config
    C%BIVgeo_Pien2023_max                      = BIVgeo_Pien2023_max_config

    ! Ice dynamics - calving
    ! ======================

    C%choice_calving_law                       = choice_calving_law_config
    C%calving_threshold_thickness              = calving_threshold_thickness_config
    C%do_remove_shelves                        = do_remove_shelves_config
    C%remove_shelves_larger_than_PD            = remove_shelves_larger_than_PD_config
    C%continental_shelf_calving                = continental_shelf_calving_config
    C%continental_shelf_min_height             = continental_shelf_min_height_config

    ! Thermodynamics and rheology
    ! ===========================

    C%choice_initial_ice_temperature           = choice_initial_ice_temperature_config
    C%uniform_ice_temperature                  = uniform_ice_temperature_config
    C%choice_thermo_model                      = choice_thermo_model_config
    C%choice_ice_rheology                      = choice_ice_rheology_config
    C%uniform_flow_factor                      = uniform_flow_factor_config
    C%choice_ice_heat_capacity                 = choice_ice_heat_capacity_config
    C%uniform_ice_heat_capacity                = uniform_ice_heat_capacity_config
    C%choice_ice_thermal_conductivity          = choice_ice_thermal_conductivity_config
    C%uniform_ice_thermal_conductivity         = uniform_ice_thermal_conductivity_config

    ! Climate
    ! =======

    C%choice_climate_model                     = choice_climate_model_config
    C%choice_idealised_climate                 = choice_idealised_climate_config

    ! Folder with NetCDF files containing direct climate forcing
    C%direct_climate_foldername_NAM            = direct_climate_foldername_NAM_config
    C%direct_climate_foldername_EAS            = direct_climate_foldername_EAS_config
    C%direct_climate_foldername_GRL            = direct_climate_foldername_GRL_config
    C%direct_climate_foldername_ANT            = direct_climate_foldername_ANT_config

    C%direct_climate_basefilename_NAM          = direct_climate_basefilename_NAM_config
    C%direct_climate_basefilename_EAS          = direct_climate_basefilename_EAS_config
    C%direct_climate_basefilename_GRL          = direct_climate_basefilename_GRL_config
    C%direct_climate_basefilename_ANT          = direct_climate_basefilename_ANT_config

    ! Whether or not to apply geometry corrections to direct temperature/precipitation
    C%do_direct_climate_geo_corr               = do_direct_climate_geo_corr_config

    ! NetCDF file containing the present-day observed climate (e.g. ERA40)
    C%filename_PD_obs_climate                  = filename_PD_obs_climate_config

    ! GCM snapshots in the matrix_warm_cold option
    C%filename_climate_snapshot_PI             = filename_climate_snapshot_PI_config
    C%filename_climate_snapshot_warm           = filename_climate_snapshot_warm_config
    C%filename_climate_snapshot_cold           = filename_climate_snapshot_cold_config
    
    ! GCM snapshots in the matrix_warm_cold option
    C%reference_mask_method                    = reference_mask_method_config 
    
    C%filename_snapshot_mask_PD_obs            = filename_snapshot_mask_PD_obs_config
    C%filename_snapshot_mask_GCM_PI            = filename_snapshot_mask_GCM_PI_config
    C%filename_snapshot_mask_GCM_warm          = filename_snapshot_mask_GCM_warm_config
    C%filename_snapshot_mask_GCM_cold          = filename_snapshot_mask_GCM_cold_config

    ! Constant lapse rate
    C%constant_lapserate                       = constant_lapserate_config

    ! Scaling factor for CO2 vs ice weights
    C%climate_matrix_CO2vsice_NAM              = climate_matrix_CO2vsice_NAM_config
    C%climate_matrix_CO2vsice_EAS              = climate_matrix_CO2vsice_EAS_config
    C%climate_matrix_CO2vsice_GRL              = climate_matrix_CO2vsice_GRL_config
    C%climate_matrix_CO2vsice_ANT              = climate_matrix_CO2vsice_ANT_config

    ! Orbit time and CO2 concentration of the warm and cold snapshots
    C%matrix_high_CO2_level                    = matrix_high_CO2_level_config
    C%matrix_low_CO2_level                     = matrix_low_CO2_level_config
    C%matrix_warm_orbit_time                   = matrix_warm_orbit_time_config
    C%matrix_cold_orbit_time                   = matrix_cold_orbit_time_config

    ! Whether or not to apply a bias correction to the GCM snapshots
    C%climate_matrix_biascorrect_warm          = climate_matrix_biascorrect_warm_config
    C%climate_matrix_biascorrect_cold          = climate_matrix_biascorrect_cold_config

    C%switch_glacial_index_precip              = switch_glacial_index_precip_config

    ! Ocean
    ! =====

    C%choice_ocean_model                       = choice_ocean_model_config
    C%choice_idealised_ocean                   = choice_idealised_ocean_config
    C%do_ocean_temperature_inversion           = do_ocean_temperature_inversion_config
    C%ocean_temperature_inv_t_start            = ocean_temperature_inv_t_start_config
    C%ocean_temperature_inv_t_end              = ocean_temperature_inv_t_end_config
    C%ocean_temperature_inv_dT0                = ocean_temperature_inv_dT0_config
    C%ocean_temperature_inv_H0                 = ocean_temperature_inv_H0_config
    C%ocean_temperature_inv_tau                = ocean_temperature_inv_tau_config
    C%ocean_temperature_inv_min                = ocean_temperature_inv_min_config
    C%ocean_temperature_inv_max                = ocean_temperature_inv_max_config
    C%inverted_ocean_filename_output           = inverted_ocean_filename_output_config

    ! NetCDF file containing the present-day observed ocean (WOA18) (NetCDF)
    C%filename_PD_obs_ocean                    = filename_PD_obs_ocean_config
    C%name_ocean_temperature_obs               = name_ocean_temperature_obs_config
    C%name_ocean_salinity_obs                  = name_ocean_salinity_obs_config

    ! GCM snapshots in the matrix_warm_cold option
    C%filename_GCM_ocean_snapshot_PI           = filename_GCM_ocean_snapshot_PI_config
    C%filename_GCM_ocean_snapshot_warm         = filename_GCM_ocean_snapshot_warm_config
    C%filename_GCM_ocean_snapshot_cold         = filename_GCM_ocean_snapshot_cold_config
    C%name_ocean_temperature_GCM               = name_ocean_temperature_GCM_config
    C%name_ocean_salinity_GCM                  = name_ocean_salinity_GCM_config

    ! Uniform ocean temperature values used when choice_ocean_model = "uniform_warm_cold"
    C%ocean_temperature_PD                     = ocean_temperature_PD_config
    C%ocean_temperature_cold                   = ocean_temperature_cold_config
    C%ocean_temperature_warm                   = ocean_temperature_warm_config

    ! Parameters used when choice_ocean_model = "matrix_warm_cold"
    C%choice_ocean_vertical_grid               = choice_ocean_vertical_grid_config
    C%ocean_vertical_grid_max_depth            = ocean_vertical_grid_max_depth_config
    C%ocean_regular_grid_dz                    = ocean_regular_grid_dz_config
    C%do_ocean_extrap                          = do_ocean_extrap_config
    C%ocean_extrap_dir                         = ocean_extrap_dir_config
    C%ocean_extrap_res                         = ocean_extrap_res_config
    C%ocean_extrap_Gauss_sigma                 = ocean_extrap_Gauss_sigma_config
    C%ocean_extrap_hires_geo_filename_NAM      = ocean_extrap_hires_geo_filename_NAM_config
    C%ocean_extrap_hires_geo_filename_EAS      = ocean_extrap_hires_geo_filename_EAS_config
    C%ocean_extrap_hires_geo_filename_GRL      = ocean_extrap_hires_geo_filename_GRL_config
    C%ocean_extrap_hires_geo_filename_ANT      = ocean_extrap_hires_geo_filename_ANT_config
    C%ocean_w_tot_hist_averaging_window        = ocean_w_tot_hist_averaging_window_config

    ! Scaling factor for CO2 vs ice weights
    C%ocean_matrix_CO2vsice_NAM                = ocean_matrix_CO2vsice_NAM_config
    C%ocean_matrix_CO2vsice_EAS                = ocean_matrix_CO2vsice_EAS_config
    C%ocean_matrix_CO2vsice_GRL                = ocean_matrix_CO2vsice_GRL_config
    C%ocean_matrix_CO2vsice_ANT                = ocean_matrix_CO2vsice_ANT_config

    ! Basin-dependent linear temperature profiles (used when choice_idealised_ocean = "linear_per_basin")
    C%ocean_T_surf_per_basin                   = ocean_T_surf_per_basin_config
    C%ocean_dT_dz_per_basin                    = ocean_dT_dz_per_basin_config
    C%do_invert_linear_per_basin               = do_invert_linear_per_basin_config

    ! Apply anomalies to baseline (used when choice_idealised_ocean = "anomalies")
    C%ocean_filename_baseline                  = ocean_filename_baseline_config
    C%ocean_foldername_aTO                     = ocean_foldername_aTO_config
    C%ocean_basefilename_aTO                   = ocean_basefilename_aTO_config
    C%ocean_foldername_aSO                     = ocean_foldername_aSO_config
    C%ocean_basefilename_aSO                   = ocean_basefilename_aSO_config

    ! Surface mass balance
    ! ====================

    C%choice_SMB_model                         = choice_SMB_model_config
    C%choice_idealised_SMB                     = choice_idealised_SMB_config
    C%SMB_uniform                              = SMB_uniform_config

    ! NetCDF file containing direct global/regional SMB forcing
    C%filename_direct_global_SMB               = filename_direct_global_SMB_config
    C%filename_direct_regional_SMB_NAM         = filename_direct_regional_SMB_NAM_config
    C%filename_direct_regional_SMB_EAS         = filename_direct_regional_SMB_EAS_config
    C%filename_direct_regional_SMB_GRL         = filename_direct_regional_SMB_GRL_config
    C%filename_direct_regional_SMB_ANT         = filename_direct_regional_SMB_ANT_config

    ! NetCDF file containing a direct time-less SMB snapshot
    C%filename_SMB_snapshot_NAM                = filename_SMB_snapshot_NAM_config
    C%filename_SMB_snapshot_EAS                = filename_SMB_snapshot_EAS_config
    C%filename_SMB_snapshot_GRL                = filename_SMB_snapshot_GRL_config
    C%filename_SMB_snapshot_ANT                = filename_SMB_snapshot_ANT_config

    ! Tuning parameters for the IMAU-ITM SMB model
    C%SMB_IMAUITM_choice_init_firn_NAM         = SMB_IMAUITM_choice_init_firn_NAM_config
    C%SMB_IMAUITM_choice_init_firn_EAS         = SMB_IMAUITM_choice_init_firn_EAS_config
    C%SMB_IMAUITM_choice_init_firn_GRL         = SMB_IMAUITM_choice_init_firn_GRL_config
    C%SMB_IMAUITM_choice_init_firn_ANT         = SMB_IMAUITM_choice_init_firn_ANT_config
    C%SMB_IMAUITM_initial_firn_thickness       = SMB_IMAUITM_initial_firn_thickness_config
    C%SMB_IMAUITM_C_abl_constant_NAM           = SMB_IMAUITM_C_abl_constant_NAM_config
    C%SMB_IMAUITM_C_abl_constant_EAS           = SMB_IMAUITM_C_abl_constant_EAS_config
    C%SMB_IMAUITM_C_abl_constant_GRL           = SMB_IMAUITM_C_abl_constant_GRL_config
    C%SMB_IMAUITM_C_abl_constant_ANT           = SMB_IMAUITM_C_abl_constant_ANT_config
    C%SMB_IMAUITM_C_abl_Ts_NAM                 = SMB_IMAUITM_C_abl_Ts_NAM_config
    C%SMB_IMAUITM_C_abl_Ts_EAS                 = SMB_IMAUITM_C_abl_Ts_EAS_config
    C%SMB_IMAUITM_C_abl_Ts_GRL                 = SMB_IMAUITM_C_abl_Ts_GRL_config
    C%SMB_IMAUITM_C_abl_Ts_ANT                 = SMB_IMAUITM_C_abl_Ts_ANT_config
    C%SMB_IMAUITM_C_abl_Q_NAM                  = SMB_IMAUITM_C_abl_Q_NAM_config
    C%SMB_IMAUITM_C_abl_Q_EAS                  = SMB_IMAUITM_C_abl_Q_EAS_config
    C%SMB_IMAUITM_C_abl_Q_GRL                  = SMB_IMAUITM_C_abl_Q_GRL_config
    C%SMB_IMAUITM_C_abl_Q_ANT                  = SMB_IMAUITM_C_abl_Q_ANT_config
    C%SMB_IMAUITM_C_refr_NAM                   = SMB_IMAUITM_C_refr_NAM_config
    C%SMB_IMAUITM_C_refr_EAS                   = SMB_IMAUITM_C_refr_EAS_config
    C%SMB_IMAUITM_C_refr_GRL                   = SMB_IMAUITM_C_refr_GRL_config
    C%SMB_IMAUITM_C_refr_ANT                   = SMB_IMAUITM_C_refr_ANT_config

    ! ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing
    ! ==============================================================

    C%ISMIP_forcing_filename_baseline          = ISMIP_forcing_filename_baseline_config
    C%ISMIP_forcing_foldername_aSMB            = ISMIP_forcing_foldername_aSMB_config
    C%ISMIP_forcing_basefilename_aSMB          = ISMIP_forcing_basefilename_aSMB_config
    C%ISMIP_forcing_foldername_dSMBdz          = ISMIP_forcing_foldername_dSMBdz_config
    C%ISMIP_forcing_basefilename_dSMBdz        = ISMIP_forcing_basefilename_dSMBdz_config
    C%ISMIP_forcing_foldername_aST             = ISMIP_forcing_foldername_aST_config
    C%ISMIP_forcing_basefilename_aST           = ISMIP_forcing_basefilename_aST_config
    C%ISMIP_forcing_foldername_dSTdz           = ISMIP_forcing_foldername_dSTdz_config
    C%ISMIP_forcing_basefilename_dSTdz         = ISMIP_forcing_basefilename_dSTdz_config

    ! Basal mass balance - sub-shelf melt
    ! ===================================

    C%choice_BMB_shelf_model                   = choice_BMB_shelf_model_config
    C%choice_idealised_BMB_shelf               = choice_idealised_BMB_shelf_config
    C%choice_BMB_sheet_model                   = choice_BMB_sheet_model_config
    C%BMB_shelf_uniform                        = BMB_shelf_uniform_config
    C%BMB_sheet_uniform                        = BMB_sheet_uniform_config
    C%choice_BMB_subgrid                       = choice_BMB_subgrid_config
    C%do_asynchronous_BMB                      = do_asynchronous_BMB_config
    C%BMB_max                                  = BMB_max_config
    C%BMB_min                                  = BMB_min_config

    C%choice_basin_scheme_NAM                  = choice_basin_scheme_NAM_config
    C%choice_basin_scheme_EAS                  = choice_basin_scheme_EAS_config
    C%choice_basin_scheme_GRL                  = choice_basin_scheme_GRL_config
    C%choice_basin_scheme_ANT                  = choice_basin_scheme_ANT_config
    C%filename_basins_NAM                      = filename_basins_NAM_config
    C%filename_basins_EAS                      = filename_basins_EAS_config
    C%filename_basins_GRL                      = filename_basins_GRL_config
    C%filename_basins_ANT                      = filename_basins_ANT_config
    C%do_merge_basins_ANT                      = do_merge_basins_ANT_config
    C%do_merge_basins_GRL                      = do_merge_basins_GRL_config

    C%choice_BMB_shelf_amplification           = choice_BMB_shelf_amplification_config
    C%basin_BMB_amplification_n_ANT            = basin_BMB_amplification_n_ANT_config
    ALLOCATE( C%basin_BMB_amplification_factor_ANT( C%basin_BMB_amplification_n_ANT))
    C%basin_BMB_amplification_factor_ANT       = basin_BMB_amplification_factor_ANT_config( 1:C%basin_BMB_amplification_n_ANT)
    C%basin_BMB_amplification_n_GRL            = basin_BMB_amplification_n_GRL_config
    ALLOCATE( C%basin_BMB_amplification_factor_GRL( C%basin_BMB_amplification_n_GRL))
    C%basin_BMB_amplification_factor_GRL       = basin_BMB_amplification_factor_GRL_config( 1:C%basin_BMB_amplification_n_GRL)

    ! Parameters for the three simple melt parameterisations from Favier et al. (2019)
    C%BMB_Favier2019_lin_GammaT                = BMB_Favier2019_lin_GammaT_config
    C%BMB_Favier2019_quad_GammaT               = BMB_Favier2019_quad_GammaT_config
    C%BMB_Favier2019_Mplus_GammaT              = BMB_Favier2019_Mplus_GammaT_config

    ! Parameters for the Lazeroms et al. (2018) plume-parameterisation BMB model
    C%BMB_Lazeroms2018_GammaT                  = BMB_Lazeroms2018_GammaT_config
    C%BMB_Lazeroms2018_find_GL_scheme          = BMB_Lazeroms2018_find_GL_scheme_config

    ! Parameters for the PICO BMB model
    C%BMB_PICO_nboxes                          = BMB_PICO_nboxes_config
    C%BMB_PICO_GammaTstar                      = BMB_PICO_GammaTstar_config

    ! Parameters for the LADDIE model
    C%filename_BMB_LADDIE                      = filename_BMB_LADDIE_config

    ! Parameters for the ANICE_legacy sub-shelf melt model
    C%T_ocean_mean_PD_NAM                      = T_ocean_mean_PD_NAM_config
    C%T_ocean_mean_PD_EAS                      = T_ocean_mean_PD_EAS_config
    C%T_ocean_mean_PD_GRL                      = T_ocean_mean_PD_GRL_config
    C%T_ocean_mean_PD_ANT                      = T_ocean_mean_PD_ANT_config
    C%T_ocean_mean_cold_NAM                    = T_ocean_mean_cold_NAM_config
    C%T_ocean_mean_cold_EAS                    = T_ocean_mean_cold_EAS_config
    C%T_ocean_mean_cold_GRL                    = T_ocean_mean_cold_GRL_config
    C%T_ocean_mean_cold_ANT                    = T_ocean_mean_cold_ANT_config
    C%T_ocean_mean_warm_NAM                    = T_ocean_mean_warm_NAM_config
    C%T_ocean_mean_warm_EAS                    = T_ocean_mean_warm_EAS_config
    C%T_ocean_mean_warm_GRL                    = T_ocean_mean_warm_GRL_config
    C%T_ocean_mean_warm_ANT                    = T_ocean_mean_warm_ANT_config

    C%BMB_deepocean_PD_NAM                     = BMB_deepocean_PD_NAM_config
    C%BMB_deepocean_PD_EAS                     = BMB_deepocean_PD_EAS_config
    C%BMB_deepocean_PD_GRL                     = BMB_deepocean_PD_GRL_config
    C%BMB_deepocean_PD_ANT                     = BMB_deepocean_PD_ANT_config
    C%BMB_deepocean_cold_NAM                   = BMB_deepocean_cold_NAM_config
    C%BMB_deepocean_cold_EAS                   = BMB_deepocean_cold_EAS_config
    C%BMB_deepocean_cold_GRL                   = BMB_deepocean_cold_GRL_config
    C%BMB_deepocean_cold_ANT                   = BMB_deepocean_cold_ANT_config
    C%BMB_deepocean_warm_NAM                   = BMB_deepocean_warm_NAM_config
    C%BMB_deepocean_warm_EAS                   = BMB_deepocean_warm_EAS_config
    C%BMB_deepocean_warm_GRL                   = BMB_deepocean_warm_GRL_config
    C%BMB_deepocean_warm_ANT                   = BMB_deepocean_warm_ANT_config

    C%BMB_shelf_exposed_PD_NAM                 = BMB_shelf_exposed_PD_NAM_config
    C%BMB_shelf_exposed_PD_EAS                 = BMB_shelf_exposed_PD_EAS_config
    C%BMB_shelf_exposed_PD_GRL                 = BMB_shelf_exposed_PD_GRL_config
    C%BMB_shelf_exposed_PD_ANT                 = BMB_shelf_exposed_PD_ANT_config
    C%BMB_shelf_exposed_cold_NAM               = BMB_shelf_exposed_cold_NAM_config
    C%BMB_shelf_exposed_cold_EAS               = BMB_shelf_exposed_cold_EAS_config
    C%BMB_shelf_exposed_cold_GRL               = BMB_shelf_exposed_cold_GRL_config
    C%BMB_shelf_exposed_warm_NAM               = BMB_shelf_exposed_warm_NAM_config
    C%BMB_shelf_exposed_cold_ANT               = BMB_shelf_exposed_cold_ANT_config
    C%BMB_shelf_exposed_warm_EAS               = BMB_shelf_exposed_warm_EAS_config
    C%BMB_shelf_exposed_warm_GRL               = BMB_shelf_exposed_warm_GRL_config
    C%BMB_shelf_exposed_warm_ANT               = BMB_shelf_exposed_warm_ANT_config

    C%subshelf_melt_factor_NAM                 = subshelf_melt_factor_NAM_config
    C%subshelf_melt_factor_EAS                 = subshelf_melt_factor_EAS_config
    C%subshelf_melt_factor_GRL                 = subshelf_melt_factor_GRL_config
    C%subshelf_melt_factor_ANT                 = subshelf_melt_factor_ANT_config

    C%deep_ocean_threshold_depth_NAM           = deep_ocean_threshold_depth_NAM_config
    C%deep_ocean_threshold_depth_EAS           = deep_ocean_threshold_depth_EAS_config
    C%deep_ocean_threshold_depth_GRL           = deep_ocean_threshold_depth_GRL_config
    C%deep_ocean_threshold_depth_ANT           = deep_ocean_threshold_depth_ANT_config

    ! Englacial isotope tracing
    ! ========================

    C%choice_ice_isotopes_model                = choice_ice_isotopes_model_config
    C%uniform_ice_d18O                         = uniform_ice_d18O_config

    ! Sea level and GIA
    ! =================

    C%do_ocean_floodfill                       = do_ocean_floodfill_config
    C%choice_sealevel_model                    = choice_sealevel_model_config
    C%fixed_sealevel                           = fixed_sealevel_config
    C%filename_sealevel_record                 = filename_sealevel_record_config
    C%sealevel_record_length                   = sealevel_record_length_config

    C%choice_GIA_model                         = choice_GIA_model_config
    C%dx_GIA                                   = dx_GIA_config
    C%ELRA_lithosphere_flex_rigidity           = ELRA_lithosphere_flex_rigidity_config
    C%ELRA_bedrock_relaxation_time             = ELRA_bedrock_relaxation_time_config
    C%ELRA_mantle_density                      = ELRA_mantle_density_config

    ! SELEN
    ! =====

    C%SELEN_run_at_t_start                     = SELEN_run_at_t_start_config
    C%SELEN_n_TDOF_iterations                  = SELEN_n_TDOF_iterations_config
    C%SELEN_n_recursion_iterations             = SELEN_n_recursion_iterations_config
    C%SELEN_use_rotational_feedback            = SELEN_use_rotational_feedback_config
    C%SELEN_n_harmonics                        = SELEN_n_harmonics_config
    C%SELEN_display_progress                   = SELEN_display_progress_config

    C%SELEN_dir                                = SELEN_dir_config
    C%SELEN_global_topo_filename               = SELEN_global_topo_filename_config
    C%SELEN_TABOO_init_filename                = SELEN_TABOO_init_filename_config
    C%SELEN_LMJ_VALUES_filename                = SELEN_LMJ_VALUES_filename_config

    C%SELEN_irreg_time_n                       = SELEN_irreg_time_n_config
    ALLOCATE( C%SELEN_irreg_time_window( C%SELEN_irreg_time_n))
    C%SELEN_irreg_time_window                  = SELEN_irreg_time_window_config( 1:C%SELEN_irreg_time_n)

    C%SELEN_lith_thickness                     = SELEN_lith_thickness_config
    C%SELEN_visc_n                             = SELEN_visc_n_config
    ALLOCATE( C%SELEN_visc_prof( C%SELEN_visc_n))
    C%SELEN_visc_prof      = SELEN_visc_prof_config( 1:C%SELEN_visc_n)

    C%SELEN_TABOO_CDE                          = SELEN_TABOO_CDE_config
    C%SELEN_TABOO_TLOVE                        = SELEN_TABOO_TLOVE_config
    C%SELEN_TABOO_DEG1                         = SELEN_TABOO_DEG1_config
    C%SELEN_TABOO_RCMB                         = SELEN_TABOO_RCMB_config

    ! Fill in some derived values
    C%SELEN_jmax       = (C%SELEN_n_harmonics + 1) * (C%SELEN_n_harmonics + 2) / 2
    C%SELEN_reg_time_n = MAX(1, INT(SUM(C%SELEN_irreg_time_window( 1:C%SELEN_irreg_time_n)))) * INT(1000. / C%dt_SELEN)

    CALL initialize_TABOO_config

    ! Which data fields will be written to the help_fields output file
    ! ================================================================

    C%help_field_01                            = help_field_01_config
    C%help_field_02                            = help_field_02_config
    C%help_field_03                            = help_field_03_config
    C%help_field_04                            = help_field_04_config
    C%help_field_05                            = help_field_05_config
    C%help_field_06                            = help_field_06_config
    C%help_field_07                            = help_field_07_config
    C%help_field_08                            = help_field_08_config
    C%help_field_09                            = help_field_09_config
    C%help_field_10                            = help_field_10_config
    C%help_field_11                            = help_field_11_config
    C%help_field_12                            = help_field_12_config
    C%help_field_13                            = help_field_13_config
    C%help_field_14                            = help_field_14_config
    C%help_field_15                            = help_field_15_config
    C%help_field_16                            = help_field_16_config
    C%help_field_17                            = help_field_17_config
    C%help_field_18                            = help_field_18_config
    C%help_field_19                            = help_field_19_config
    C%help_field_20                            = help_field_20_config
    C%help_field_21                            = help_field_21_config
    C%help_field_22                            = help_field_22_config
    C%help_field_23                            = help_field_23_config
    C%help_field_24                            = help_field_24_config
    C%help_field_25                            = help_field_25_config
    C%help_field_26                            = help_field_26_config
    C%help_field_27                            = help_field_27_config
    C%help_field_28                            = help_field_28_config
    C%help_field_29                            = help_field_29_config
    C%help_field_30                            = help_field_30_config
    C%help_field_31                            = help_field_31_config
    C%help_field_32                            = help_field_32_config
    C%help_field_33                            = help_field_33_config
    C%help_field_34                            = help_field_34_config
    C%help_field_35                            = help_field_35_config
    C%help_field_36                            = help_field_36_config
    C%help_field_37                            = help_field_37_config
    C%help_field_38                            = help_field_38_config
    C%help_field_39                            = help_field_39_config
    C%help_field_40                            = help_field_40_config
    C%help_field_41                            = help_field_41_config
    C%help_field_42                            = help_field_42_config
    C%help_field_43                            = help_field_43_config
    C%help_field_44                            = help_field_44_config
    C%help_field_45                            = help_field_45_config
    C%help_field_46                            = help_field_46_config
    C%help_field_47                            = help_field_47_config
    C%help_field_48                            = help_field_48_config
    C%help_field_49                            = help_field_49_config
    C%help_field_50                            = help_field_50_config

    ! Values to be filled into the total mask (used only for diagnostic output)
    ! ==========================================================================

    C%type_land                                = 0
    C%type_ocean                               = 1
    C%type_lake                                = 2
    C%type_sheet                               = 3
    C%type_shelf                               = 4
    C%type_coast                               = 5
    C%type_margin                              = 6
    C%type_groundingline                       = 7
    C%type_calvingfront                        = 8

  END SUBROUTINE copy_variables_to_struct

  SUBROUTINE initialize_TABOO_config

    ! Taboo settings
    C_TABOO%IMODE            = 1
    C_TABOO%NV               = C%SELEN_visc_n
    ALLOCATE( C_TABOO%VSC( C%SELEN_visc_n))
    C_TABOO%VSC              = C%SELEN_visc_prof
    C_TABOO%CDE              = C%SELEN_TABOO_CDE
    C_TABOO%TLOVE            = C%SELEN_TABOO_TLOVE
    C_TABOO%DEG1             = C%SELEN_TABOO_DEG1
    C_TABOO%LTH              = C%SELEN_lith_thickness
    C_TABOO%RCMB             = C%SELEN_TABOO_RCMB

  END SUBROUTINE initialize_TABOO_config

  SUBROUTINE get_procedural_output_dir_name( output_dir)
    ! Generate a procedural output directory for the current date (e.g. results_20210721_001)
    ! Keep increasing the counter at the end until a directory is available.

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(20),                       INTENT(INOUT) :: output_dir

    ! Local variables:
    INTEGER,  DIMENSION(8)                             :: values
    LOGICAL                                            :: ex

    CALL date_and_time(VALUES=values)

    ! Get proper year (assume we're still in the 21st century...)
    output_dir(1:10) = 'results_20'
    SELECT CASE( FLOOR(REAL(values(1))/10._dp)-200)
    CASE(0)
      output_dir(11:11) = '0'
    CASE(1)
      output_dir(11:11) = '1'
    CASE(2)
      output_dir(11:11) = '2'
    CASE(3)
      output_dir(11:11) = '3'
    CASE(4)
      output_dir(11:11) = '4'
    CASE(5)
      output_dir(11:11) = '5'
    CASE(6)
      output_dir(11:11) = '6'
    CASE(7)
      output_dir(11:11) = '7'
    CASE(8)
      output_dir(11:11) = '8'
    CASE(9)
      output_dir(11:11) = '9'
    CASE DEFAULT
      WRITE(0,*) 'get_procedural_output_dir: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( MOD(values(1),10))
    CASE(0)
      output_dir(12:12) = '0'
    CASE(1)
      output_dir(12:12) = '1'
    CASE(2)
      output_dir(12:12) = '2'
    CASE(3)
      output_dir(12:12) = '3'
    CASE(4)
      output_dir(12:12) = '4'
    CASE(5)
      output_dir(12:12) = '5'
    CASE(6)
      output_dir(12:12) = '6'
    CASE(7)
      output_dir(12:12) = '7'
    CASE(8)
      output_dir(12:12) = '8'
    CASE(9)
      output_dir(12:12) = '9'
    CASE DEFAULT
      WRITE(0,*) 'get_procedural_output_dir: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( values(2))
    CASE(1)
      output_dir(13:14) = '01'
    CASE(2)
      output_dir(13:14) = '02'
    CASE(3)
      output_dir(13:14) = '03'
    CASE(4)
      output_dir(13:14) = '04'
    CASE(5)
      output_dir(13:14) = '05'
    CASE(6)
      output_dir(13:14) = '06'
    CASE(7)
      output_dir(13:14) = '07'
    CASE(8)
      output_dir(13:14) = '08'
    CASE(9)
      output_dir(13:14) = '09'
    CASE(10)
      output_dir(13:14) = '10'
    CASE(11)
      output_dir(13:14) = '11'
    CASE(12)
      output_dir(13:14) = '12'
    CASE DEFAULT
      WRITE(0,*) 'get_procedural_output_dir: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( FLOOR(REAL(values(3))/10._dp))
    CASE(0)
      output_dir(15:15) = '0'
    CASE(1)
      output_dir(15:15) = '1'
    CASE(2)
      output_dir(15:15) = '2'
    CASE(3)
      output_dir(15:15) = '3'
    CASE DEFAULT
      WRITE(0,*) 'get_procedural_output_dir: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( MOD(values(3),10))
    CASE(0)
      output_dir(16:16) = '0'
    CASE(1)
      output_dir(16:16) = '1'
    CASE(2)
      output_dir(16:16) = '2'
    CASE(3)
      output_dir(16:16) = '3'
    CASE(4)
      output_dir(16:16) = '4'
    CASE(5)
      output_dir(16:16) = '5'
    CASE(6)
      output_dir(16:16) = '6'
    CASE(7)
      output_dir(16:16) = '7'
    CASE(8)
      output_dir(16:16) = '8'
    CASE(9)
      output_dir(16:16) = '9'
    CASE DEFAULT
      WRITE(0,*) 'get_procedural_output_dir: ERROR retrieving date and time!'
    END SELECT

    output_dir(17:20) = '_001'

    INQUIRE( FILE = TRIM(output_dir)//'/.', EXIST=ex )

    DO WHILE (ex)

     IF      (output_dir(20:20) == '0') THEN
       output_dir(20:20) = '1'
     ELSE IF (output_dir(20:20) == '1') THEN
       output_dir(20:20) = '2'
     ELSE IF (output_dir(20:20) == '2') THEN
       output_dir(20:20) = '3'
     ELSE IF (output_dir(20:20) == '3') THEN
       output_dir(20:20) = '4'
     ELSE IF (output_dir(20:20) == '4') THEN
       output_dir(20:20) = '5'
     ELSE IF (output_dir(20:20) == '5') THEN
       output_dir(20:20) = '6'
     ELSE IF (output_dir(20:20) == '6') THEN
       output_dir(20:20) = '7'
     ELSE IF (output_dir(20:20) == '7') THEN
       output_dir(20:20) = '8'
     ELSE IF (output_dir(20:20) == '8') THEN
       output_dir(20:20) = '9'
     ELSE IF (output_dir(20:20) == '9') THEN
       output_dir(20:20) = '0'

       IF      (output_dir(19:19) == '0') THEN
         output_dir(19:19) = '1'
       ELSE IF (output_dir(19:19) == '1') THEN
         output_dir(19:19) = '2'
       ELSE IF (output_dir(19:19) == '2') THEN
         output_dir(19:19) = '3'
       ELSE IF (output_dir(19:19) == '3') THEN
         output_dir(19:19) = '4'
       ELSE IF (output_dir(19:19) == '4') THEN
         output_dir(19:19) = '5'
       ELSE IF (output_dir(19:19) == '5') THEN
         output_dir(19:19) = '6'
       ELSE IF (output_dir(19:19) == '6') THEN
         output_dir(19:19) = '7'
       ELSE IF (output_dir(19:19) == '7') THEN
         output_dir(19:19) = '8'
       ELSE IF (output_dir(19:19) == '8') THEN
         output_dir(19:19) = '9'
       ELSE IF (output_dir(19:19) == '9') THEN
         output_dir(19:19) = '0'

         IF      (output_dir(18:18) == '0') THEN
           output_dir(18:18) = '1'
         ELSE IF (output_dir(18:18) == '1') THEN
           output_dir(18:18) = '2'
         ELSE IF (output_dir(18:18) == '2') THEN
           output_dir(18:18) = '3'
         ELSE IF (output_dir(18:18) == '3') THEN
           output_dir(18:18) = '4'
         ELSE IF (output_dir(18:18) == '4') THEN
           output_dir(18:18) = '5'
         ELSE IF (output_dir(18:18) == '5') THEN
           output_dir(18:18) = '6'
         ELSE IF (output_dir(18:18) == '6') THEN
           output_dir(18:18) = '7'
         ELSE IF (output_dir(18:18) == '7') THEN
           output_dir(18:18) = '8'
         ELSE IF (output_dir(18:18) == '8') THEN
           output_dir(18:18) = '9'
         ELSE IF (output_dir(18:18) == '9') THEN
           output_dir(18:18) = '0'
         END IF

       END IF

     END IF

     INQUIRE( FILE=TRIM(output_dir)//'/.', EXIST=ex )

    END DO

  END SUBROUTINE get_procedural_output_dir_name

  SUBROUTINE write_total_model_time_to_screen( tstart, tstop)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: tstart, tstop

    ! Local variables
    REAL(dp)                                           :: dt
    INTEGER                                            :: nr, ns, nm, nh, nd

    dt = tstop - tstart

    ns = CEILING(dt)

    nr = MOD(ns, 60*60*24)
    nd = (ns - nr) / (60*60*24)
    ns = ns - (nd*60*60*24)

    nr = MOD(ns, 60*60)
    nh = (ns - nr) / (60*60)
    ns = ns - (nh*60*60)

    nr = MOD(ns, 60)
    nm = (ns - nr) / (60)
    ns = ns - (nm*60)

    WRITE(0,'(A)') ''
    WRITE(0,'(A)') ' ================================================================================'
    WRITE(0,'(A,I2,A,I2,A,I2,A,I2,A)') ' ===== Simulation finished in ', nd, ' days, ', nh, ' hours, ', nm, ' minutes and ', ns, ' seconds! ====='
    WRITE(0,'(A)') ' ================================================================================'
    WRITE(0,'(A)') ''

  END SUBROUTINE write_total_model_time_to_screen

  ! Routines for the extended error messaging / debugging system
  SUBROUTINE init_routine( routine_name, do_track_resource_use)
    ! Initialise an IMAU-ICE subroutine; update the routine path

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                  INTENT(IN)    :: routine_name
    LOGICAL,                   OPTIONAL, INTENT(IN)    :: do_track_resource_use

    ! Local variables:
    INTEGER                                            :: len_path_tot, len_path_used, len_name
    INTEGER                                            :: ierr, cerr
    INTEGER                                            :: i
    LOGICAL                                            :: do_track_resource_use_loc

    ! Check if routine_name has enough memory
    len_path_tot  = LEN(      routine_path)
    len_path_used = LEN_TRIM( routine_path)
    len_name      = LEN_TRIM( routine_name)

    IF (len_path_used + 1 + len_name > len_path_tot) THEN
      WRITE(0,*) 'init_routine - ERROR: routine_path = "', TRIM( routine_path), '", no more space to append routine_name = "', TRIM( routine_name), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Append this routine to the routine path
    routine_path = TRIM( routine_path) // '/' // TRIM( routine_name)

    ! Check if resource use for this subroutine should be tracked
    ! (especially for the NetCDF routines we don't want to do this, as there are
    ! a great many of them and the resource tracker output file will become annoyingly big)

    IF (PRESENT( do_track_resource_use)) THEN
      do_track_resource_use_loc = do_track_resource_use
    ELSE
      do_track_resource_use_loc = .TRUE.
    END IF

    IF (do_track_resource_use_loc) THEN

      ! Initialise the computation time tracker
      CALL find_subroutine_in_resource_tracker( i)
      resource_tracker( i)%tstart = MPI_WTIME()

      ! Check maximum MPI window at the start of the routine
      resource_tracker( i)%n_MPI_windows_init = n_MPI_windows

    ELSE

      routine_path = TRIM( routine_path) // '_NOTRACK'

    END IF


  END SUBROUTINE init_routine

  SUBROUTINE finalise_routine( routine_name, n_extra_windows_expected)
    ! Finalise; remove the current routine name from the routine path

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                  INTENT(IN)    :: routine_name
    INTEGER,  INTENT(IN), OPTIONAL                     :: n_extra_windows_expected

    ! Local variables:
    LOGICAL                                            :: do_track_resource_use
    INTEGER                                            :: len_path_tot, i, ii
    INTEGER                                            :: ierr, cerr
    REAL(dp)                                           :: dt
    INTEGER                                            :: n_extra_windows_expected_loc, n_extra_windows_found

    ! Check if resource use should be tracked for this subroutine
    i = INDEX( routine_path, '_NOTRACK')
    IF ( i == 0) THEN
      do_track_resource_use = .TRUE.
    ELSE
      do_track_resource_use = .FALSE.
    END IF

    IF (do_track_resource_use) THEN
      ! Add computation time to the resource tracker
      CALL find_subroutine_in_resource_tracker( i)
      dt = MPI_WTIME() - resource_tracker( i)%tstart
      resource_tracker( i)%tcomp = resource_tracker( i)%tcomp + dt

      ! Check maximum MPI window at the end of the routine
      resource_tracker( i)%n_MPI_windows_final = n_MPI_windows

      ! If it is larger than expected, warn that there might be a memory leak
      n_extra_windows_expected_loc = 0
      IF (PRESENT( n_extra_windows_expected)) n_extra_windows_expected_loc = n_extra_windows_expected
      n_extra_windows_found = resource_tracker( i)%n_MPI_windows_final - resource_tracker( i)%n_MPI_windows_init

      ii = INDEX( routine_path, 'IMAU_ICE_program/initialise_')
      IF (ii == 0 .AND. n_extra_windows_found > n_extra_windows_expected_loc) THEN
        ! This subroutine has more memory allocated at the start than at the beginning.
        CALL warning('more memory was allocated and not freed than expected; possible memory leak! (expected {int_01} extra windows, found {int_02})', &
          int_01 = n_extra_windows_expected_loc, int_02 = n_extra_windows_found)
      END IF

      ! Find where in the string exactly the current routine name is located
      i = INDEX( routine_path, routine_name)

      IF (i == 0) THEN
        WRITE(0,*) 'finalise_routine - ERROR: routine_name = "', TRIM( routine_name), '" not found in routine_path = "', TRIM( routine_path), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

      ! Remove the current routine name from the routine path
      len_path_tot = LEN( routine_path)
      routine_path( i-1:len_path_tot) = ' '

    ELSE ! IF (do_track_resource_use) THEN
      ! Resource use for this subroutine should not be tracked

      ! Find where in the string exactly the current routine name is located
      i = INDEX( routine_path, TRIM( routine_name) // '_NOTRACK')

      IF (i == 0) THEN
        WRITE(0,*) 'finalise_routine - ERROR: routine_name = "', TRIM( routine_name), '" not found in routine_path = "', TRIM( routine_path), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF

      ! Remove the current routine name from the routine path
      len_path_tot = LEN( routine_path)
      routine_path( i-1:len_path_tot) = ' '

    END IF ! IF (do_track_resource_use) THEN

  END SUBROUTINE finalise_routine
  SUBROUTINE find_subroutine_in_resource_tracker( i)
    ! Find the current subroutine in the resource tracker. If it's not there yet, add it.

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,                             INTENT(OUT)   :: i

    ! Local variables:
    INTEGER                                            :: n

    n = SIZE( resource_tracker)

    DO i = 1, n
      IF     (resource_tracker( i)%routine_path == routine_path) THEN
        ! The current subroutine is listed at this position in the resource tracker
        RETURN
      ELSEIF (resource_tracker( i)%routine_path == 'subroutine_placeholder') THEN
        ! We've checked all listed subroutines and haven't found the current one; add it
        resource_tracker( i)%routine_path = routine_path
        RETURN
      END IF
    END DO

    ! If we've reached this point, then the resource tracker is overflowing
    CALL crash('Resource tracker overflows! Allocate more memory for it in initialise_model_configuration.')

  END SUBROUTINE find_subroutine_in_resource_tracker
  SUBROUTINE reset_computation_times
    ! Reset the computation times in the resource tracker

    IMPLICIT NONE

    ! Local variables:
    INTEGER                                            :: i,n

    n = SIZE( resource_tracker)

    DO i = 1, n
      resource_tracker( i)%tstart = 0._dp
      resource_tracker( i)%tcomp  = 0._dp
    END DO

  END SUBROUTINE reset_computation_times
  SUBROUTINE crash( err_msg, int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
                              dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)
    ! Crash the model, write the error message to the screen

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: err_msg
    INTEGER,  INTENT(IN), OPTIONAL                     :: int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10
    REAL(dp), INTENT(IN), OPTIONAL                     ::  dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10

    ! Local variables:
    CHARACTER(LEN=1024)                                :: err_msg_loc
    INTEGER                                            :: ierr, cerr, pari, parn, nc
    CHARACTER(LEN=9)                                   :: fmt
    CHARACTER(LEN=:), ALLOCATABLE                      :: process_str

    ! Get local, edit-able copy of error message string
    err_msg_loc = err_msg

    ! Get rank of current process and total number of processes
    ! (needed because the configuration_module cannot access the par structure)
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, pari, ierr)
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, parn, ierr)

    ! Set the process string (e.g. "05/16")
    IF     (parn < 10) THEN
      nc = 1
    ELSEIF (parn < 100) THEN
      nc = 2
    ELSEIF (parn < 1000) THEN
      nc = 3
    ELSEIF (parn < 10000) THEN
      nc = 4
    ELSE
      nc = 5
    END IF

    WRITE( fmt,'(A,I1,A,I1,A)') '(I', nc, ',A,I', nc, ')'
    ALLOCATE(CHARACTER(2*nc+1) :: process_str)
    WRITE( process_str,fmt) pari, '/', parn

    ! Insert numbers into string if needed
    IF (PRESENT( int_01)) CALL insert_val_into_string_int( err_msg_loc, '{int_01}', int_01)
    IF (PRESENT( int_02)) CALL insert_val_into_string_int( err_msg_loc, '{int_02}', int_02)
    IF (PRESENT( int_03)) CALL insert_val_into_string_int( err_msg_loc, '{int_03}', int_03)
    IF (PRESENT( int_04)) CALL insert_val_into_string_int( err_msg_loc, '{int_04}', int_04)
    IF (PRESENT( int_05)) CALL insert_val_into_string_int( err_msg_loc, '{int_05}', int_05)
    IF (PRESENT( int_06)) CALL insert_val_into_string_int( err_msg_loc, '{int_06}', int_06)
    IF (PRESENT( int_07)) CALL insert_val_into_string_int( err_msg_loc, '{int_07}', int_07)
    IF (PRESENT( int_08)) CALL insert_val_into_string_int( err_msg_loc, '{int_08}', int_08)
    IF (PRESENT( int_09)) CALL insert_val_into_string_int( err_msg_loc, '{int_09}', int_09)
    IF (PRESENT( int_10)) CALL insert_val_into_string_int( err_msg_loc, '{int_10}', int_10)

    IF (PRESENT( dp_01 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_01}' , dp_01 )
    IF (PRESENT( dp_02 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_02}' , dp_02 )
    IF (PRESENT( dp_03 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_03}' , dp_03 )
    IF (PRESENT( dp_04 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_04}' , dp_04 )
    IF (PRESENT( dp_05 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_05}' , dp_05 )
    IF (PRESENT( dp_06 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_06}' , dp_06 )
    IF (PRESENT( dp_07 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_07}' , dp_07 )
    IF (PRESENT( dp_08 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_08}' , dp_08 )
    IF (PRESENT( dp_09 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_09}' , dp_09 )
    IF (PRESENT( dp_10 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_10}' , dp_10 )

    ! Write the error to the screen
    WRITE(0,'(A,A,A,A,A,A)') colour_string('ERROR: ' // TRIM( err_msg_loc),'red') // ' in ' // colour_string( TRIM(routine_path),'light blue') // &
      ' on process ', colour_string( process_str,'light blue'), ' (0 = master)'

    ! Stop the program
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

  END SUBROUTINE crash
  SUBROUTINE warning( err_msg, int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
                                dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)
    ! Write the warning message to the screen, but don't crash the model

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: err_msg
    INTEGER,  INTENT(IN), OPTIONAL                     :: int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10
    REAL(dp), INTENT(IN), OPTIONAL                     ::  dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10

    ! Local variables:
    CHARACTER(LEN=1024)                                :: err_msg_loc
    INTEGER                                            :: ierr, pari, parn, nc
    CHARACTER(LEN=9)                                   :: fmt
    CHARACTER(LEN=:), ALLOCATABLE                      :: process_str

    ! Get local, edit-able copy of error message string
    err_msg_loc = err_msg

    ! Get rank of current process and total number of processes
    ! (needed because the configuration_module cannot access the par structure)
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, pari, ierr)
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, parn, ierr)

    ! Set the process string (e.g. "05/16")
    IF     (parn < 10) THEN
      nc = 1
    ELSEIF (parn < 100) THEN
      nc = 2
    ELSEIF (parn < 1000) THEN
      nc = 3
    ELSEIF (parn < 10000) THEN
      nc = 4
    ELSE
      nc = 5
    END IF

    WRITE( fmt,'(A,I1,A,I1,A)') '(I', nc, ',A,I', nc, ')'
    ALLOCATE(CHARACTER(2*nc+1) :: process_str)
    WRITE( process_str,fmt) pari, '/', parn

    ! Insert numbers into string if needed
    IF (PRESENT( int_01)) CALL insert_val_into_string_int( err_msg_loc, '{int_01}', int_01)
    IF (PRESENT( int_02)) CALL insert_val_into_string_int( err_msg_loc, '{int_02}', int_02)
    IF (PRESENT( int_03)) CALL insert_val_into_string_int( err_msg_loc, '{int_03}', int_03)
    IF (PRESENT( int_04)) CALL insert_val_into_string_int( err_msg_loc, '{int_04}', int_04)
    IF (PRESENT( int_05)) CALL insert_val_into_string_int( err_msg_loc, '{int_05}', int_05)
    IF (PRESENT( int_06)) CALL insert_val_into_string_int( err_msg_loc, '{int_06}', int_06)
    IF (PRESENT( int_07)) CALL insert_val_into_string_int( err_msg_loc, '{int_07}', int_07)
    IF (PRESENT( int_08)) CALL insert_val_into_string_int( err_msg_loc, '{int_08}', int_08)
    IF (PRESENT( int_09)) CALL insert_val_into_string_int( err_msg_loc, '{int_09}', int_09)
    IF (PRESENT( int_10)) CALL insert_val_into_string_int( err_msg_loc, '{int_10}', int_10)

    IF (PRESENT( dp_01 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_01}' , dp_01 )
    IF (PRESENT( dp_02 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_02}' , dp_02 )
    IF (PRESENT( dp_03 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_03}' , dp_03 )
    IF (PRESENT( dp_04 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_04}' , dp_04 )
    IF (PRESENT( dp_05 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_05}' , dp_05 )
    IF (PRESENT( dp_06 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_06}' , dp_06 )
    IF (PRESENT( dp_07 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_07}' , dp_07 )
    IF (PRESENT( dp_08 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_08}' , dp_08 )
    IF (PRESENT( dp_09 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_09}' , dp_09 )
    IF (PRESENT( dp_10 )) CALL insert_val_into_string_dp(  err_msg_loc, '{dp_10}' , dp_10 )

    ! Write the error to the screen
    WRITE(0,'(A,A,A,A,A,A)') colour_string('WARNING: ' // TRIM( err_msg_loc),'yellow') // ' in ' // colour_string( TRIM(routine_path),'light blue') // &
      ' on process ', colour_string( process_str,'light blue'), ' (0 = master)'

    ! Clean up after yourself
    DEALLOCATE( process_str)

  END SUBROUTINE warning
  FUNCTION colour_string( str, col) RESULT( str_col)
    ! Add colour to a string for writing to the terminal

    IMPLICIT NONE

    ! Input variables:
    CHARACTER(LEN=*),                    INTENT(IN)    :: str, col

    ! Result variables:
    CHARACTER(LEN=:), ALLOCATABLE                      :: str_col

    ALLOCATE(CHARACTER(LEN(str)+9) :: str_col)       ! The +9 is just enough to store the color characters

    ! The 91m gives red, 0m sets the default back
    ! Available colors: 90:gray, 91:red, 92:green, 93:yellow, 94:blue, 95:pink, 96:light blue
    IF     (col == 'gray') THEN
      str_col = achar(27)//'[90m'//str//achar(27)//'[0m'
    ELSEIF (col == 'red') THEN
      str_col = achar(27)//'[91m'//str//achar(27)//'[0m'
    ELSEIF (col == 'green') THEN
      str_col = achar(27)//'[92m'//str//achar(27)//'[0m'
    ELSEIF (col == 'yellow') THEN
      str_col = achar(27)//'[93m'//str//achar(27)//'[0m'
    ELSEIF (col == 'blue') THEN
      str_col = achar(27)//'[94m'//str//achar(27)//'[0m'
    ELSEIF (col == 'pink') THEN
      str_col = achar(27)//'[95m'//str//achar(27)//'[0m'
    ELSEIF (col == 'light blue') THEN
      str_col = achar(27)//'[96m'//str//achar(27)//'[0m'
    ELSE
      WRITE(0,*) ''
    END IF

  END FUNCTION colour_string
  SUBROUTINE insert_val_into_string_int( str, marker, val)
    ! Replace marker in str with val (where val is an integer)
    !
    ! Example: str    = 'Johnny has {int_01} apples.'
    !          marker = '{int_01}'
    !          val    = 5
    !
    ! This returns: str = 'Johnny has 5 apples'

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(INOUT) :: str
    CHARACTER(LEN=*),                    INTENT(IN)    :: marker
    INTEGER,                             INTENT(IN)    :: val

    ! Local variables:
    INTEGER                                            :: ci
    INTEGER                                            :: nc
    CHARACTER(LEN=4)                                   :: fmt
    CHARACTER(LEN=:), ALLOCATABLE                      :: val_str
    INTEGER                                            :: len_str, len_marker

    ! Find position ci in str where i_str occurs
    ci = INDEX( str, marker)

    ! Safety
    IF (ci == 0) CALL crash('insert_val_into_string_int: couldnt find marker "' // TRIM( marker) // '" in string "' // TRIM( str) // '"!')

    ! Write val to a string
    IF     (ABS( val) < 10) THEN
      nc = 1
    ELSEIF (ABS( val) < 100) THEN
      nc = 2
    ELSEIF (ABS( val) < 1000) THEN
      nc = 3
    ELSEIF (ABS( val) < 10000) THEN
      nc = 4
    ELSEIF (ABS( val) < 100000) THEN
      nc = 5
    ELSEIF (ABS( val) < 1000000) THEN
      nc = 6
    ELSEIF (ABS( val) < 10000000) THEN
      nc = 7
    ELSEIF (ABS( val) < 100000000) THEN
      nc = 8
    ELSE
      nc = 9
    END IF
    ! Add room for a minus sign if needed
    IF (val < 0) nc = nc + 1

    WRITE( fmt,'(A,I1,A)') '(I', nc, ')'
    ALLOCATE(CHARACTER(nc) :: val_str)
    WRITE( val_str,fmt) val

    ! Find total string length right now
    len_str    = LEN( str)
    len_marker = LEN( marker)

    ! Insert the integer string into the string
    str = str(1:ci-1) // val_str // str(ci+len_marker:len_str)

    ! Clean up after yourself
    DEALLOCATE( val_str)

  END SUBROUTINE insert_val_into_string_int
  SUBROUTINE insert_val_into_string_dp( str, marker, val)
    ! Replace marker in str with val (where val is a double-precision number)
    !
    ! Example: str    = 'Johnny weighs {dp_01} kg.'
    !          marker = '{dp_01}'
    !          val    = 57.098
    !
    ! This returns: str = 'Johnny weighs 57.098 kg'

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=*),                    INTENT(INOUT) :: str
    CHARACTER(LEN=*),                    INTENT(IN)    :: marker
    REAL(dp),                            INTENT(IN)    :: val

    ! Local variables:
    INTEGER                                            :: ci
    CHARACTER(LEN=11)                                  :: val_str
    INTEGER                                            :: len_str, len_marker

    ! Find position ci in str where i_str occurs
    ci = INDEX( str, marker)

    ! Safety
    IF (ci == 0) CALL crash('insert_val_into_string_dp: couldnt find marker "' // TRIM( marker) // '" in string "' // TRIM( str) // '"!')

    ! Write val to a string
    WRITE( val_str,'(E11.5)') val

    ! Find total string length right now
    len_str    = LEN( str)
    len_marker = LEN( marker)

    ! Insert the integer string into the string
    str = str(1:ci-1) // val_str // str(ci+len_marker:len_str)

  END SUBROUTINE insert_val_into_string_dp

  SUBROUTINE capitalise_string( str)
    ! Changes a string which contains lower case letters to a string with upper case letters
    ! Useful for case-insensitive string comparison

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                   INTENT(INOUT) :: str

    ! Local variables:
    INTEGER                                             :: i, index_cap

    CHARACTER(26), PARAMETER                            :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(26), PARAMETER                            :: low = 'abcdefghijklmnopqrstuvwxyz'

    DO i = 1, LEN_TRIM( str)
      index_cap = INDEX( low, str( i:i))
      IF (index_cap > 0) str( i:i) = cap( index_cap:index_cap)
    END DO

  END SUBROUTINE capitalise_string
  SUBROUTINE remove_leading_spaces( str)
    ! Changes a string which contains lower case letters to a string with upper case letters
    ! Useful for case-insensitive string comparison

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                   INTENT(INOUT) :: str

    ! Local variables:
    INTEGER                                             :: lstr

    DO WHILE (str( 1:1) == ' ' .AND. LEN_TRIM( str) > 0)
      lstr = LEN_TRIM( str)
      str( 1:lstr-1) = str( 2:lstr)
      str( lstr:lstr) = ' '
    END DO

  END SUBROUTINE remove_leading_spaces

END MODULE configuration_module
