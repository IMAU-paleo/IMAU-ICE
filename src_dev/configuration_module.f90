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

  IMPLICIT NONE
  
  INTEGER, PARAMETER  :: dp  = KIND(1.0D0)  ! Kind of double precision numbers. Reals should be declared as: REAL(dp) :: example

  ! ===================================================================================
  ! "_config  variables, which will be collected into a NAMELIST, and possibly replaced
  ! by the values in the external config file. Remember the "_config" extension!
  ! ===================================================================================

  ! Time steps and range
  ! ====================
  
  REAL(dp)            :: start_time_of_run_config                = 0.0_dp                           ! Start time (in years) of the simulations
  REAL(dp)            :: end_time_of_run_config                  = 50000.0_dp                       ! End   time (in years) of the simulations
  REAL(dp)            :: dt_coupling_config                      = 100._dp                          ! Interval of coupling (in years) between the four ice-sheets  
  REAL(dp)            :: dt_max_config                           = 10.0_dp                          ! Maximum time step (in years) of the ice model
  REAL(dp)            :: dt_thermo_config                        = 10.0_dp                          ! Time step (in years) for updating thermodynamics
  REAL(dp)            :: dt_climate_config                       = 10._dp                           ! Time step (in years) for updating the climate
  REAL(dp)            :: dt_SMB_config                           = 10._dp                           ! Time step (in years) for updating the SMB
  REAL(dp)            :: dt_BMB_config                           = 10._dp                           ! Time step (in years) for updating the BMB
  REAL(dp)            :: dt_bedrock_ELRA_config                  = 100._dp                          ! Time step (in years) for updating the bedrock deformation rate with the ELRA model
  REAL(dp)            :: dt_SELEN_config                         = 1000._dp                         ! Time step (in years) for calling SELEN
  REAL(dp)            :: dt_output_config                        = 5000.0_dp                        ! Time step (in years) for writing output
  
  ! Which ice sheets do we simulate?
  ! ================================
  
  LOGICAL             :: do_NAM_config                           = .FALSE.                          ! North America
  LOGICAL             :: do_EAS_config                           = .FALSE.                          ! Eurasia
  LOGICAL             :: do_GRL_config                           = .FALSE.                          ! Greenland
  LOGICAL             :: do_ANT_config                           = .TRUE.                           ! Antarctica
  
  ! Benchmark experiments
  ! =====================
  
  LOGICAL             :: do_benchmark_experiment_config          = .TRUE.
  CHARACTER(LEN=256)  :: choice_benchmark_experiment_config      = 'EISMINT_I'
  REAL(dp)            :: SSA_icestream_m_config                  = 1                                ! Values tested by Schoof are 1, 10, and 20
  REAL(dp)            :: ISMIP_HOM_L_config                      = 160000.0                         ! Domain size of the ISMIP-HOM benchmarks
  CHARACTER(LEN=256)  :: ISMIP_HOM_E_Arolla_filename_config      = 'arolla100.dat'                  ! Path to the Haut Glacier d'Arolla input file
  
  ! Whether or not to let IMAU_ICE dynamically create its own output folder.
  ! This works fine locally, on LISA its better to use a fixed folder name.
  ! =======================================================================
  
  LOGICAL             :: create_new_output_dir_config            = .TRUE.
  CHARACTER(LEN=256)  :: output_dir_config                       = 'results_IMAU_ICE'
  LOGICAL             :: do_write_debug_data_config              = .FALSE.
  
  ! Whether or not the simulation is a restart of a previous simulation
  ! ===================================================================
  
  LOGICAL             :: is_restart_config                       = .FALSE.
  REAL(dp)            :: time_to_restart_from_config             = 0._dp                            ! Can be different from C%start_time_of_run, though this will issue a warning
  
  ! Horizontal grid spacing and size for the four regions
  ! =====================================================
  
  REAL(dp)            :: dx_NAM_config                           = 40000._dp
  REAL(dp)            :: dx_EAS_config                           = 40000._dp
  REAL(dp)            :: dx_GRL_config                           = 20000._dp
  REAL(dp)            :: dx_ANT_config                           = 40000._dp

  ! The scaled vertical coordinate zeta, used mainly in thermodynamics
  ! ==================================================================
  
  INTEGER                        :: nz_config                    = 15
  REAL(dp), DIMENSION(210), SAVE :: zeta_config                  = &
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

  ! Input data file paths
  ! =====================
  
  ! Initial model state (NetCDF)
  CHARACTER(LEN=256)  :: filename_init_NAM_config                = 'Datasets/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_init_EAS_config                = 'Datasets/ETOPO1/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_init_GRL_config                = 'Datasets/Bedmachine_Greenland/Greenland_BedMachine_5km_noEllesmere.nc'
  CHARACTER(LEN=256)  :: filename_init_ANT_config                = 'Datasets/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'
  
  ! PD reference data (NetCDF)
  CHARACTER(LEN=256)  :: filename_PD_NAM_config                  = 'Datasets/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_PD_EAS_config                  = 'Datasets/ETOPO1/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256)  :: filename_PD_GRL_config                  = 'Datasets/Bedmachine_Greenland/Greenland_BedMachine_5km_noEllesmere.nc'
  CHARACTER(LEN=256)  :: filename_PD_ANT_config                  = 'Datasets/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc' 
   
  ! Insolation forcing (NetCDF) (Laskar et al., 2004)
  CHARACTER(LEN=256)  :: filename_insolation_config              = '/Datasets/Insolation_laskar/Insolation_Laskar_etal_2004.nc'
  
  ! CO2 record (ASCII text file, so the number of rows needs to be specified)
  CHARACTER(LEN=256)  :: filename_CO2_record_config              = 'Datasets/CO2/EPICA_CO2_Bereiter_2015_100yr.dat'
  INTEGER             :: CO2_record_length_config                = 8001
  
  ! d18O record (ASCII text file, so the number of rows needs to be specified)
  CHARACTER(LEN=256)  :: filename_d18O_record_config             = 'Datasets/d18O/Ahn2017_d18O.dat'
  INTEGER             :: d18O_record_length_config               = 2051
  
  ! Ice5G ice-sheet geometry
  CHARACTER(LEN=256)  :: filename_ICE5G_PD_config                = 'Datasets/ICE5G/ice5g_v1.2_00.0k_1deg.nc'
  CHARACTER(LEN=256)  :: filename_ICE5G_LGM_config               = 'Datasets/ICE5G/ice5g_v1.2_21.0k_1deg.nc'

  ! Ice dynamics - velocity
  ! =======================
  
  CHARACTER(LEN=256)  :: choice_ice_dynamics_config              = 'DIVA'                           ! Can be 'SIA', 'SSA', 'SIA/SSA', and 'DIVA'
  REAL(dp)            :: n_flow_config                           = 3.0_dp                           ! Exponent in Glen's flow law
  REAL(dp)            :: m_enh_sheet_config                      = 1.0_dp                           ! Ice flow enhancement factor for grounded ice
  REAL(dp)            :: m_enh_shelf_config                      = 1.0_dp                           ! Ice flow enhancement factor for floating ice
  LOGICAL             :: use_analytical_GL_flux_config           = .FALSE.                          ! Whether or not the analytical grounding line flux solution is used
  
  ! Some parameters for numerically solving the SSA/DIVA
  CHARACTER(LEN=256)  :: DIVA_choice_ice_margin_config           = 'infinite_slab'                   ! Choice of ice margin boundary conditions; can be "BC" or "infinite_slab"
  REAL(dp)            :: DIVA_norm_dUV_tol_config                = 1E-2_dp                          ! Successive solutions of UV in the effective viscosity iteration must not differ by more than this amount (on average)
  INTEGER             :: DIVA_max_outer_loops_config             = 50                               ! Maximum number of effective viscosity iterations
  REAL(dp)            :: DIVA_max_residual_UV_config             = 2.5_dp                           ! The maximum residual in U and V in the SOR solver of the linearised DIVA must drop below this value [m yr^-1]
  REAL(dp)            :: DIVA_SOR_omega_config                   = 1.0_dp                           ! The over-relaxation parameter in the SOR solver of the linearised DIVA
  INTEGER             :: DIVA_max_inner_loops_config             = 3000                             ! Maximum number of iterations in the SOR solver of the linearised DIVA
  REAL(dp)            :: DIVA_beta_max_config                    = 1E20_dp                          ! The DIVA is not solved (i.e. velocities are assumed to be zero) for beta values larger than this
  REAL(dp)            :: DIVA_err_lim_config                     = 1E-5_dp                          ! The DIVA is not refined (i.e. velocities are no longer updated with SOR) wherever successive velocity iterations change the velocity by less than this amount
  REAL(dp)            :: DIVA_vel_max_config                     = 5000._dp                         ! DIVA velocities are limited to this value (u,v evaluated separately)
  REAL(dp)            :: DIVA_vel_min_config                     = 1E-5_dp                          ! DIVA velocities below this value are set to zero (u,v evaluated separately)
  REAL(dp)            :: DIVA_visc_iter_relax_config             = 0.7_dp                           ! Relaxation parameter for subsequent viscosity iterations (for improved stability)
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_u_west_config          = 'infinite'                       ! Boundary conditions for the ice velocity field at the domain boundary in the DIVA
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_u_east_config          = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_u_south_config         = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_u_north_config         = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_v_west_config          = 'infinite'                       ! Boundary conditions for the ice velocity field at the domain boundary in the DIVA
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_v_east_config          = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_v_south_config         = 'infinite'
  CHARACTER(LEN=256)  :: DIVA_boundary_BC_v_north_config         = 'infinite'

  ! Ice dynamics - time integration
  ! ===============================
  
  CHARACTER(LEN=256)  :: choice_timestepping_config              = 'pc'                             ! Can be 'direct' or 'pc' (NOTE: 'direct' does not work with DIVA ice dynamcis!)
  CHARACTER(LEN=256)  :: choice_ice_integration_method_config    = 'explicit'                       ! Can be "explicit" or "implicit"
  REAL(dp)            :: dHi_SOR_omega_config                    = 1.2_dp                           ! The over-relaxation parameter    in the SOR solver of the mass conservation matrix equation
  REAL(dp)            :: dHi_SOR_max_res_config                  = 10._dp                           ! The stop criterion               in the SOR solver of the mass conservation matrix equation
  INTEGER             :: dHi_SOR_max_loops_config                = 10000                            ! The maximum number of iterations in the SOR solver of the mass conservation matrix equation
  
  ! Predictor-corrector ice-thickness update
  REAL(dp)            :: pc_epsilon_config                       = 3._dp                            ! Target truncation error in dHi_dt [m/yr] (epsilon in Robinson et al., 2020, Eq. 33)
  REAL(dp)            :: pc_k_I_config                           = 0.2_dp                           ! Exponent k_I in  Robinson et al., 2020, Eq. 33
  REAL(dp)            :: pc_k_p_config                           = 0.2_dp                           ! Exponent k_p in  Robinson et al., 2020, Eq. 33
  REAL(dp)            :: pc_eta_min_config                       = 1E-8_dp                          ! Normalisation term in estimation of the truncation error (Robinson et al., Eq. 32)
  INTEGER             :: pc_max_timestep_iterations_config       = 5                                ! Maximum number of iterations of each time step
  REAL(dp)            :: pc_redo_tol_config                      = 10._dp                           ! Maximum allowed truncation error (any higher and the timestep is decreased)
  REAL(dp)            :: dt_min_config                           = 0.01_dp                          ! Smallest allowed time step [yr]

  ! Ice dynamics - sliding
  ! ======================
  LOGICAL             :: no_sliding_config                       = .FALSE.                          ! If set to TRUE, no sliding is allowed anywhere (used for some schematic experiments)
  CHARACTER(LEN=256)  :: choice_sliding_law_config               = 'Coulomb_regularised'            ! Choice of sliding law (currently implemented: "Coulomb", "Coulomb_regularised")
  REAL(dp)            :: C_sliding_config                        = 1.0E7_dp                         ! Factor   in Weertman sliding law
  REAL(dp)            :: m_sliding_config                        = 1._dp/3._dp                      ! Exponent in Weertman sliding law

  ! Ice dynamics - calving
  ! ======================
  
  CHARACTER(LEN=256)  :: choice_calving_law_config               = 'threshold_thickness'            ! Choice of calving law (currently only "threshold_thickness" is implemented)
  REAL(dp)            :: calving_threshold_thickness_config      = 200._dp                          ! Threshold ice thickness in the "threshold_thickness" calving law (200m taken from ANICE)

  ! Thermodynamics
  ! ==============
  
  CHARACTER(LEN=256)  :: choice_geothermal_heat_flux_config      = 'constant'                       ! Choice of geothermal heat flux; can be 'constant' or 'spatial'
  REAL(dp)            :: constant_geothermal_heat_flux_config    = 1.72E06_dp                       ! Geothermal Heat flux [J m^-2 yr^-1] Sclater et al. (1980)
  CHARACTER(LEN=256)  :: filename_geothermal_heat_flux_config    = 'Datasets/GHF/geothermal_heatflux_ShapiroRitzwoller2004_global_1x1_deg.nc'
  
  ! Climate matrix
  ! ==============
  
  ! Present-day observed climate (ERA40) (NetCDF)
  CHARACTER(LEN=256)  :: filename_PD_obs_climate_config          = 'Datasets/ERA40/ERA40_climate_global.nc'
  
  ! GCM snapshots
  CHARACTER(LEN=256)  :: choice_climate_matrix_config            = 'PI_LGM'                         ! 'PI_LGM' uses 2 snapshots
  CHARACTER(LEN=256)  :: filename_GCM_snapshot_PI_config         = 'Datasets/GCM_snapshots/Singarayer_Valdes_2010_PI_Control.nc'
  CHARACTER(LEN=256)  :: filename_GCM_snapshot_LGM_config        = 'Datasets/GCM_snapshots/Singarayer_Valdes_2010_LGM.nc'
  
  ! Ocean temperature (used for both thermodynamics and basal melt)
  CHARACTER(LEN=256)  :: choice_ocean_temperature_model_config   = 'scaled'                         ! Can be "fixed" (use PD value) or "scaled" (scale between "PD", "warm", and "cold" values based on forcing (prescribed or inverse-modelled))
  REAL(dp)            :: ocean_temperature_PD_config             = 271.46_dp                        ! present day temperature of the ocean beneath the shelves [K; -1.7 Celsius]
  REAL(dp)            :: ocean_temperature_cold_config           = 268.16_dp                        ! cold period temperature of the ocean beneath the shelves [K; -5.0 Celcius]
  REAL(dp)            :: ocean_temperature_warm_config           = 275.16_dp                        ! warm period temperature of the ocean beneath the shelves [K;  2.0 Celcius]
  
  REAL(dp)            :: constant_lapserate_config               = 0.008_dp                         ! Constant atmospheric lapse rate [K m^-1]
  
  ! Forcing
  ! =======
  
  ! The choice of forcing:
  ! 'd18O_inverse_dT_glob' : Use the inverse routine with the specified d18O record to calculate a global temperature offset (e.g. de Boer et al., 2013)
  ! 'CO2_direct'           : Use the specified CO2 record to force the climate matrix (e.g. Berends et al., 2018)
  ! 'd18O_inverse_CO2'     : Use the inverse routine with the specified d18O record to calculate CO2 and then force the climate matrix (e.g. Berends et al., 2019)
  CHARACTER(LEN=256)  :: choice_forcing_method_config            = 'd18O_inverse_dT_glob'
  
  REAL(dp)            :: dT_deepwater_averaging_window_config    = 3000                             ! Time window (in yr) over which global mean temperature anomaly is averaged to find the deep-water temperature anomaly
  REAL(dp)            :: dT_deepwater_dT_surf_ratio_config       = 0.25_dp                          ! Ratio between global mean surface temperature change and deep-water temperature change
  REAL(dp)            :: d18O_dT_deepwater_ratio_config          = -0.28_dp                         ! Ratio between deep-water temperature change and benthic d18O change
  
  REAL(dp)            :: dT_glob_inverse_averaging_window_config = 2000._dp                         ! Time window (in yr) over which global mean temperature anomaly is averaged before changing it with the inverse routine
  REAL(dp)            :: inverse_d18O_to_dT_glob_scaling_config  = 20._dp                           ! Scaling factor between modelled d18O anomaly and prescribed temperature anomaly change (value from de Boer et al., 2013)
  REAL(dp)            :: CO2_inverse_averaging_window_config     = 2000._dp                         ! Time window (in yr) over which CO2                             is averaged before changing it with the inverse routine
  REAL(dp)            :: inverse_d18O_to_CO2_scaling_config      = 68._dp                           ! Scaling factor between modelled d18O anomaly and modelled CO2 change (value from Berends et al., 2019)
  REAL(dp)            :: inverse_d18O_to_CO2_initial_CO2_config  = 280._dp                          ! CO2 value at the start of the simulation when using the inverse method to calculate CO2
  
  ! SMB tuning
  ! ==========
  
  REAL(dp)            :: C_abl_constant_NAM_config               = -49._dp                          ! 34._dp    (commented values are old ANICE defaults, but since refreezing was not calculated right
  REAL(dp)            :: C_abl_constant_EAS_config               = -49._dp                          !            and this has since been fixed, these values will still not give the same results as
  REAL(dp)            :: C_abl_constant_GRL_config               = -49._dp                          !            they used to in ANICE.)
  REAL(dp)            :: C_abl_constant_ANT_config               = -49._dp
  REAL(dp)            :: C_abl_Ts_NAM_config                     = 10._dp                           ! 10._dp
  REAL(dp)            :: C_abl_Ts_EAS_config                     = 10._dp
  REAL(dp)            :: C_abl_Ts_GRL_config                     = 10._dp
  REAL(dp)            :: C_abl_Ts_ANT_config                     = 10._dp
  REAL(dp)            :: C_abl_Q_NAM_config                      = 0.0227_dp                        ! 0.513_dp
  REAL(dp)            :: C_abl_Q_EAS_config                      = 0.0227_dp
  REAL(dp)            :: C_abl_Q_GRL_config                      = 0.0227_dp
  REAL(dp)            :: C_abl_Q_ANT_config                      = 0.0227_dp
  REAL(dp)            :: C_refr_NAM_config                       = 0.051_dp                         ! 0.012_dp
  REAL(dp)            :: C_refr_EAS_config                       = 0.051_dp 
  REAL(dp)            :: C_refr_GRL_config                       = 0.051_dp 
  REAL(dp)            :: C_refr_ANT_config                       = 0.051_dp
  
  ! Sub-shelf melt parameterisation
  ! ===============================
  
  ! Mean ocean temperatures and deep ocean & exposed shelf melt rates were tuned by Bas de Boer in 2011.
  REAL(dp)            :: T_ocean_mean_PD_NAM_config              = -1.7_dp                          ! Present day temperature of the ocean beneath the shelves [Celcius]
  REAL(dp)            :: T_ocean_mean_PD_EAS_config              = -1.7_dp
  REAL(dp)            :: T_ocean_mean_PD_GRL_config              =  2.0_dp
  REAL(dp)            :: T_ocean_mean_PD_ANT_config              = -1.7_dp
  REAL(dp)            :: T_ocean_mean_cold_NAM_config            = -5.0_dp                          ! Cold period temperature of the ocean beneath the shelves [Celcius]
  REAL(dp)            :: T_ocean_mean_cold_EAS_config            = -5.0_dp
  REAL(dp)            :: T_ocean_mean_cold_GRL_config            =  0.0_dp
  REAL(dp)            :: T_ocean_mean_cold_ANT_config            = -5.0_dp
  REAL(dp)            :: T_ocean_mean_warm_NAM_config            =  2.0_dp                          ! Warm period temperature of the ocean beneath the shelves [Celcius]
  REAL(dp)            :: T_ocean_mean_warm_EAS_config            =  2.0_dp
  REAL(dp)            :: T_ocean_mean_warm_GRL_config            =  4.0_dp
  REAL(dp)            :: T_ocean_mean_warm_ANT_config            =  2.0_dp
            
  REAL(dp)            :: BMB_deepocean_PD_NAM_config             =  -5._dp                          ! Present-day sub-shelf melt rate for deep-ocean areas [m/year]
  REAL(dp)            :: BMB_deepocean_PD_EAS_config             =  -5._dp
  REAL(dp)            :: BMB_deepocean_PD_GRL_config             =  -5._dp
  REAL(dp)            :: BMB_deepocean_PD_ANT_config             =  -5._dp
  REAL(dp)            :: BMB_deepocean_cold_NAM_config           =  -2._dp                          ! Cold period sub-shelf melt rate for deep-ocean areas [m/year]
  REAL(dp)            :: BMB_deepocean_cold_EAS_config           =  -2._dp
  REAL(dp)            :: BMB_deepocean_cold_GRL_config           =  -2._dp
  REAL(dp)            :: BMB_deepocean_cold_ANT_config           =  -2._dp
  REAL(dp)            :: BMB_deepocean_warm_NAM_config           = -10._dp                          ! Warm period sub-shelf melt rate for deep-ocean areas [m/year]    
  REAL(dp)            :: BMB_deepocean_warm_EAS_config           = -10._dp
  REAL(dp)            :: BMB_deepocean_warm_GRL_config           = -10._dp
  REAL(dp)            :: BMB_deepocean_warm_ANT_config           = -10._dp

  REAL(dp)            :: BMB_shelf_exposed_PD_NAM_config         =  -3._dp                          ! Present-day sub-shelf melt rate for exposed areas    [m/year]
  REAL(dp)            :: BMB_shelf_exposed_PD_EAS_config         =  -3._dp
  REAL(dp)            :: BMB_shelf_exposed_PD_GRL_config         =  -3._dp
  REAL(dp)            :: BMB_shelf_exposed_PD_ANT_config         =  -3._dp
  REAL(dp)            :: BMB_shelf_exposed_cold_NAM_config       =  -0._dp                          ! Cold period sub-shelf melt rate for exposed areas    [m/year]
  REAL(dp)            :: BMB_shelf_exposed_cold_EAS_config       =  -0._dp
  REAL(dp)            :: BMB_shelf_exposed_cold_GRL_config       =  -0._dp
  REAL(dp)            :: BMB_shelf_exposed_cold_ANT_config       =  -0._dp
  REAL(dp)            :: BMB_shelf_exposed_warm_NAM_config       =  -6._dp                          ! Warm period sub-shelf melt rate for exposed areas    [m/year]
  REAL(dp)            :: BMB_shelf_exposed_warm_EAS_config       =  -6._dp
  REAL(dp)            :: BMB_shelf_exposed_warm_GRL_config       =  -6._dp
  REAL(dp)            :: BMB_shelf_exposed_warm_ANT_config       =  -6._dp
    
  REAL(dp)            :: subshelf_melt_factor_NAM_config         = 0.005_dp                         ! Overall tuning factor for sub-shelf melt rate
  REAL(dp)            :: subshelf_melt_factor_EAS_config         = 0.005_dp
  REAL(dp)            :: subshelf_melt_factor_GRL_config         = 0.005_dp
  REAL(dp)            :: subshelf_melt_factor_ANT_config         = 0.005_dp
  
  REAL(dp)            :: deep_ocean_threshold_depth_NAM_config   = 1200._dp                         ! Threshold water depth for "deep ocean" (as opposed to continental shelf);
  REAL(dp)            :: deep_ocean_threshold_depth_EAS_config   = 800._dp                          ! this mostly prevents ice shelves from growing beyond the continental shelf
  REAL(dp)            :: deep_ocean_threshold_depth_GRL_config   = 800._dp                          ! Different depths for different regions is a bit ad hoc, but in reality
  REAL(dp)            :: deep_ocean_threshold_depth_ANT_config   = 1800._dp                         ! the different surface ocean temperatures probably result in the same effect...
    
  ! Sea level and GIA
  ! =================
  
  LOGICAL             :: do_ocean_floodfill_config               = .TRUE.                           ! Use a flood-fill to determine the ocean mask, so that (pro-/sub-glacial) lakes dont exist
  CHARACTER(LEN=256)  :: choice_sealevel_model_config            = 'eustatic'                       ! Can be "fixed", "prescribed", "eustatic", or "SELEN"
  REAL(dp)            :: fixed_sealevel_config                   = 0._dp                            ! Height of fixed sealevel w.r.t. PD
  CHARACTER(LEN=256)  :: filename_sealevel_record_config         = 'name_of_file.dat'               ! Filename of a file containing a sealevel record
  INTEGER             :: sealevel_record_length_config           = 1
  
  CHARACTER(LEN=256)  :: choice_GIA_model_config                 = 'ELRA'                           ! Can be "none", "ELRA", or "SELEN"
  REAL(dp)            :: dx_GIA_config                           = 100000._dp                       ! Horizontal resolution of the square grid used for the GIA model
  REAL(dp)            :: ELRA_lithosphere_flex_rigidity_config   = 1.0E+25_dp                       ! Lithospheric flexural rigidity [kg m^2 s^-2]
  REAL(dp)            :: ELRA_bedrock_relaxation_time_config     = 3000.0_dp                        ! Relaxation time for bedrock adjustment [yr]
  REAL(dp)            :: ELRA_mantle_density_config              = 3300.0_dp                        ! Mantle density [kg m^-3]

  ! SELEN
  ! =====
  
  LOGICAL             :: SELEN_run_at_t_start_config              = .FALSE.                         ! Whether or not to run SELEN in the first coupling loop (needed for some benchmark experiments)
  INTEGER             :: SELEN_n_TDOF_iterations_config           = 1                               ! Number of Time-Dependent Ocean Function iterations
  INTEGER             :: SELEN_n_recursion_iterations_config      = 1                               ! Number of recursion iterations
  LOGICAL             :: SELEN_use_rotational_feedback_config     = .FALSE.                         ! If TRUE, rotational feedback is included
  INTEGER             :: SELEN_n_harmonics_config                 = 128                             ! Maximum number of harmonic degrees
  LOGICAL             :: SELEN_display_progress_config            = .FALSE.                         ! Whether or not to display the progress of the big loops to the screen (doesn't work on Cartesius!)
  
  CHARACTER(LEN=256)  :: SELEN_dir_config                         = 'SELEN_files'                   ! Directory where SELEN initial files and spherical harmonics are stored
  CHARACTER(LEN=256)  :: SELEN_global_topo_filename_config        = 'SELEN_global_topography.nc'    ! Filename for the SELEN global topography file (located in SELEN_dir)
  CHARACTER(LEN=256)  :: SELEN_TABOO_init_filename_config         = 'SELEN_TABOO_initial_file.dat'  ! Filename for the TABOO initial file           (idem                )
  CHARACTER(LEN=256)  :: SELEN_LMJ_VALUES_filename_config         = 'SELEN_lmj_values.bin'          ! Filename for the LJ and MJ values file        (idem                )
  
  INTEGER                  :: SELEN_irreg_time_n_config           = 15                              ! Number of entries in the irregular moving time window
  REAL(dp), DIMENSION(50)  :: SELEN_irreg_time_window_config      = &                               ! Values of entries in the irregular moving time window
   (/20._dp, 20._dp, 20._dp, 5._dp, 5._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, 1._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, &
      0._dp,  0._dp,  0._dp, 0._dp, 0._dp  /)
      
  REAL(dp)            :: SELEN_lith_thickness_config              = 100._dp                         ! Thickness of the elastic lithosphere [km]
  INTEGER             :: SELEN_visc_n_config                      = 3                               ! Number      of viscous asthenosphere layers
  REAL(dp), DIMENSION(3) :: SELEN_visc_prof_config                = (/ 3._dp, 0.6_dp, 0.3_dp /)     ! Viscosities of viscous asthenosphere layers [?]
    
  ! Settings for the TABOO Earth deformation model
  INTEGER             :: SELEN_TABOO_CDE_config                   = 0                               ! code of the model (see taboo for explanation)
  INTEGER             :: SELEN_TABOO_TLOVE_config                 = 1                               ! Tidal love numbers yes/no
  INTEGER             :: SELEN_TABOO_DEG1_config                  = 1                               ! Tidal love numbers degree
  REAL(dp)            :: SELEN_TABOO_RCMB_config                  = 3480._dp                        ! Radius of CMB (km)
  
  ! Which data fields will be written to the help_fields output file
  ! ================================================================
  
  CHARACTER(LEN=256)  :: help_field_01_config                     = 'lat'
  CHARACTER(LEN=256)  :: help_field_02_config                     = 'lon'
  CHARACTER(LEN=256)  :: help_field_03_config                     = 'uabs_surf'
  CHARACTER(LEN=256)  :: help_field_04_config                     = 'uabs_base'
  CHARACTER(LEN=256)  :: help_field_05_config                     = 'uabs_vav'
  CHARACTER(LEN=256)  :: help_field_06_config                     = 'u_3D'
  CHARACTER(LEN=256)  :: help_field_07_config                     = 'v_3D'
  CHARACTER(LEN=256)  :: help_field_08_config                     = 'w_3D'
  CHARACTER(LEN=256)  :: help_field_09_config                     = 'mask'
  CHARACTER(LEN=256)  :: help_field_10_config                     = 'dHs_dx'
  CHARACTER(LEN=256)  :: help_field_11_config                     = 'dHs_dy'
  CHARACTER(LEN=256)  :: help_field_12_config                     = 'T2m_year'
  CHARACTER(LEN=256)  :: help_field_13_config                     = 'Precip_year'
  CHARACTER(LEN=256)  :: help_field_14_config                     = 'Albedo_year'
  CHARACTER(LEN=256)  :: help_field_15_config                     = 'SMB_year'
  CHARACTER(LEN=256)  :: help_field_16_config                     = 'BMB'
  CHARACTER(LEN=256)  :: help_field_17_config                     = 'T2m'
  CHARACTER(LEN=256)  :: help_field_18_config                     = 'Precip'
  CHARACTER(LEN=256)  :: help_field_19_config                     = 'Albedo'
  CHARACTER(LEN=256)  :: help_field_20_config                     = 'SMB'
  CHARACTER(LEN=256)  :: help_field_21_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_22_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_23_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_24_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_25_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_26_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_27_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_28_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_29_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_30_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_31_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_32_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_33_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_34_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_35_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_36_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_37_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_38_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_39_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_40_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_41_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_42_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_43_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_44_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_45_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_46_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_47_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_48_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_49_config                     = 'none'
  CHARACTER(LEN=256)  :: help_field_50_config                     = 'none'


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
    REAL(dp)                            :: dt_thermo
    REAL(dp)                            :: dt_climate
    REAL(dp)                            :: dt_SMB
    REAL(dp)                            :: dt_BMB
    REAL(dp)                            :: dt_bedrock_ELRA
    REAL(dp)                            :: dt_SELEN
    REAL(dp)                            :: dt_output
    
    ! Which ice sheets do we simulate?
    ! ================================
    
    LOGICAL                             :: do_NAM
    LOGICAL                             :: do_EAS
    LOGICAL                             :: do_GRL
    LOGICAL                             :: do_ANT
    
    ! Benchmark experiments
    ! =====================
    
    LOGICAL                             :: do_benchmark_experiment
    CHARACTER(LEN=256)                  :: choice_benchmark_experiment
    REAL(dp)                            :: SSA_icestream_m
    REAL(dp)                            :: ISMIP_HOM_L
    CHARACTER(LEN=256)                  :: ISMIP_HOM_E_Arolla_filename

    ! Whether or not to let IMAU_ICE dynamically create its own output folder
   ! =======================================================================
   
    LOGICAL                             :: create_new_output_dir
    CHARACTER(LEN=256)                  :: output_dir
    LOGICAL                             :: do_write_debug_data
  
    ! Whether or not the simulation is a restart of a previous simulation
    ! ===================================================================
   
    LOGICAL                             :: is_restart
    REAL(dp)                            :: time_to_restart_from
  
    ! Horizontal grid spacing and size for the four regions
    ! =====================================================
    
    REAL(dp)                            :: dx_NAM
    REAL(dp)                            :: dx_EAS
    REAL(dp)                            :: dx_GRL
    REAL(dp)                            :: dx_ANT

    ! Scaled vertical coordinate zeta  
    ! ===============================
     
    INTEGER                             :: nz
    REAL(dp), DIMENSION(:), ALLOCATABLE :: zeta

    ! Input data file paths
    ! =====================
                            
    CHARACTER(LEN=256)                  :: filename_init_NAM
    CHARACTER(LEN=256)                  :: filename_init_EAS
    CHARACTER(LEN=256)                  :: filename_init_GRL
    CHARACTER(LEN=256)                  :: filename_init_ANT
    
    CHARACTER(LEN=256)                  :: filename_PD_NAM
    CHARACTER(LEN=256)                  :: filename_PD_EAS
    CHARACTER(LEN=256)                  :: filename_PD_GRL
    CHARACTER(LEN=256)                  :: filename_PD_ANT
    
    CHARACTER(LEN=256)                  :: filename_insolation
    
    CHARACTER(LEN=256)                  :: filename_CO2_record
    INTEGER                             :: CO2_record_length
    CHARACTER(LEN=256)                  :: filename_d18O_record
    INTEGER                             :: d18O_record_length
  
    CHARACTER(LEN=256)                  :: filename_ICE5G_PD
    CHARACTER(LEN=256)                  :: filename_ICE5G_LGM

    ! Ice dynamics - velocity
    ! =======================
    
    CHARACTER(LEN=256)                  :: choice_ice_dynamics
    REAL(dp)                            :: n_flow
    REAL(dp)                            :: m_enh_sheet
    REAL(dp)                            :: m_enh_shelf
    LOGICAL                             :: use_analytical_GL_flux
    
    ! Some parameters for numerically solving the SSA/DIVA
    CHARACTER(LEN=256)                  :: DIVA_choice_ice_margin
    REAL(dp)                            :: DIVA_norm_dUV_tol
    INTEGER                             :: DIVA_max_outer_loops
    REAL(dp)                            :: DIVA_max_residual_UV
    REAL(dp)                            :: DIVA_SOR_omega
    INTEGER                             :: DIVA_max_inner_loops
    REAL(dp)                            :: DIVA_beta_max
    REAL(dp)                            :: DIVA_err_lim
    REAL(dp)                            :: DIVA_vel_max
    REAL(dp)                            :: DIVA_vel_min
    REAL(dp)                            :: DIVA_visc_iter_relax
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_u_west
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_u_east
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_u_south
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_u_north
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_v_west
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_v_east
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_v_south
    CHARACTER(LEN=256)                  :: DIVA_boundary_BC_v_north
  
    ! Ice dynamics - time integration
    ! ===============================
    
    CHARACTER(LEN=256)                  :: choice_timestepping
    CHARACTER(LEN=256)                  :: choice_ice_integration_method
    REAL(dp)                            :: dHi_SOR_omega
    REAL(dp)                            :: dHi_SOR_max_res
    INTEGER                             :: dHi_SOR_max_loops
    
    ! Predictor-corrector ice-thickness update
    REAL(dp)                            :: pc_epsilon
    REAL(dp)                            :: pc_k_I
    REAL(dp)                            :: pc_k_p
    REAL(dp)                            :: pc_eta_min
    INTEGER                             :: pc_max_timestep_iterations
    REAL(dp)                            :: pc_redo_tol
    REAL(dp)                            :: dt_min
  
    ! Ice dynamics - sliding
    ! ======================
    LOGICAL                             :: no_sliding
    CHARACTER(LEN=256)                  :: choice_sliding_law
    REAL(dp)                            :: C_sliding
    REAL(dp)                            :: m_sliding
  
    ! Ice dynamics - calving
    ! ======================
    
    CHARACTER(LEN=256)                  :: choice_calving_law
    REAL(dp)                            :: calving_threshold_thickness
  
    ! Thermodynamics
    ! ==============
    
    CHARACTER(LEN=256)                  :: choice_geothermal_heat_flux
    REAL(dp)                            :: constant_geothermal_heat_flux
    CHARACTER(LEN=256)                  :: filename_geothermal_heat_flux
    
    ! Climate matrix
    ! ==============
    
    CHARACTER(LEN=256)                  :: filename_PD_obs_climate
    CHARACTER(LEN=256)                  :: choice_climate_matrix
    CHARACTER(LEN=256)                  :: filename_GCM_snapshot_PI
    CHARACTER(LEN=256)                  :: filename_GCM_snapshot_LGM
    
    CHARACTER(LEN=256)                  :: choice_ocean_temperature_model
    REAL(dp)                            :: ocean_temperature_PD
    REAL(dp)                            :: ocean_temperature_cold
    REAL(dp)                            :: ocean_temperature_warm
    
    REAL(dp)                            :: constant_lapserate
    
    ! Forcing
    ! =======
    
    CHARACTER(LEN=256)                  :: choice_forcing_method
    
    REAL(dp)                            :: dT_deepwater_averaging_window
    REAL(dp)                            :: dT_deepwater_dT_surf_ratio
    REAL(dp)                            :: d18O_dT_deepwater_ratio
    
    REAL(dp)                            :: dT_glob_inverse_averaging_window
    REAL(dp)                            :: inverse_d18O_to_dT_glob_scaling
    REAL(dp)                            :: CO2_inverse_averaging_window
    REAL(dp)                            :: inverse_d18O_to_CO2_scaling
    REAL(dp)                            :: inverse_d18O_to_CO2_initial_CO2
    
    ! SMB melt tuning
    ! ===============
    
    REAL(dp)                            :: C_abl_constant_NAM
    REAL(dp)                            :: C_abl_constant_EAS
    REAL(dp)                            :: C_abl_constant_GRL
    REAL(dp)                            :: C_abl_constant_ANT
    REAL(dp)                            :: C_abl_Ts_NAM
    REAL(dp)                            :: C_abl_Ts_EAS
    REAL(dp)                            :: C_abl_Ts_GRL
    REAL(dp)                            :: C_abl_Ts_ANT
    REAL(dp)                            :: C_abl_Q_NAM
    REAL(dp)                            :: C_abl_Q_EAS
    REAL(dp)                            :: C_abl_Q_GRL
    REAL(dp)                            :: C_abl_Q_ANT
    REAL(dp)                            :: C_refr_NAM
    REAL(dp)                            :: C_refr_EAS
    REAL(dp)                            :: C_refr_GRL
    REAL(dp)                            :: C_refr_ANT
  
    ! Sub-shelf melt parameterisation
    ! ===============================
    
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
    
    ! Parameters of the polar stereographic projections of the four model regions
    ! (These have to match the values used to create the input files!)
    ! ===========================================================================     
                          
    REAL(dp)                            :: lambda_M_NAM                           
    REAL(dp)                            :: lambda_M_EAS                           
    REAL(dp)                            :: lambda_M_GRL                           
    REAL(dp)                            :: lambda_M_ANT
    REAL(dp)                            :: phi_M_NAM
    REAL(dp)                            :: phi_M_EAS
    REAL(dp)                            :: phi_M_GRL
    REAL(dp)                            :: phi_M_ANT
    REAL(dp)                            :: alpha_stereo_NAM
    REAL(dp)                            :: alpha_stereo_EAS
    REAL(dp)                            :: alpha_stereo_GRL
    REAL(dp)                            :: alpha_stereo_ANT
    
    
    ! Some useful constants for the scaled vertical coordinate transformation
    ! =======================================================================

    REAL(dp), DIMENSION(:), ALLOCATABLE :: a_k
    REAL(dp), DIMENSION(:), ALLOCATABLE :: b_k
    REAL(dp), DIMENSION(:), ALLOCATABLE :: c_k
    REAL(dp), DIMENSION(:), ALLOCATABLE :: d_k
    REAL(dp), DIMENSION(:), ALLOCATABLE :: a_zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: b_zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: c_zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: a_zetazeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: b_zetazeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: c_zetazeta

    REAL(dp), DIMENSION(:), ALLOCATABLE :: z_zeta_minus
    REAL(dp), DIMENSION(:), ALLOCATABLE :: a_zeta_minus
    REAL(dp), DIMENSION(:), ALLOCATABLE :: b_zeta_minus

    REAL(dp), DIMENSION(:), ALLOCATABLE :: b_zeta_plus
    REAL(dp), DIMENSION(:), ALLOCATABLE :: c_zeta_plus
    REAL(dp), DIMENSION(:), ALLOCATABLE :: d_zeta_plus

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

  SUBROUTINE read_config_file( config_filename)
    ! Use a NAMELIST containing all the "_config" variables to read
    ! an external config file, and overwrite the default values of
    ! the specified variables with the values from the file.
    
    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256),INTENT(IN) :: config_filename
    
    INTEGER, PARAMETER :: config_unit = 28 ! Unit number which is used for the configuration file.
    INTEGER            :: ios
    
    ! The NAMELIST that's used to read the external config file.
    NAMELIST /CONFIG/start_time_of_run_config,                   &                
                     end_time_of_run_config,                     &
                     dt_coupling_config,                         &
                     dt_max_config,                              &
                     dt_thermo_config,                           &
                     dt_climate_config,                          &
                     dt_SMB_config,                              &
                     dt_BMB_config,                              &
                     dt_bedrock_ELRA_config,                     &
                     dt_SELEN_config,                            &
                     dt_output_config,                           &
                     do_NAM_config,                              &
                     do_EAS_config,                              &
                     do_GRL_config,                              &
                     do_ANT_config,                              &
                     do_benchmark_experiment_config,             &
                     choice_benchmark_experiment_config,         &
                     SSA_icestream_m_config,                     &
                     ISMIP_HOM_L_config,                         &
                     ISMIP_HOM_E_Arolla_filename_config,         &
                     create_new_output_dir_config,               &
                     output_dir_config,                          &
                     do_write_debug_data_config,                 &
                     is_restart_config,                          &
                     time_to_restart_from_config,                &
                     dx_NAM_config,                              &
                     dx_EAS_config,                              &
                     dx_GRL_config,                              &
                     dx_ANT_config,                              &
                     nz_config,                                  &
                     zeta_config,                                &
                     filename_init_NAM_config,                   &
                     filename_init_EAS_config,                   &
                     filename_init_GRL_config,                   &
                     filename_init_ANT_config,                   &
                     filename_PD_NAM_config,                     &
                     filename_PD_EAS_config,                     &
                     filename_PD_GRL_config,                     &
                     filename_PD_ANT_config,                     &
                     filename_insolation_config,                 &
                     filename_CO2_record_config,                 &
                     CO2_record_length_config,                   &
                     filename_d18O_record_config,                &
                     d18O_record_length_config,                  &
                     filename_ICE5G_PD_config,                   &
                     filename_ICE5G_LGM_config,                  &
                     choice_ice_dynamics_config,                 &
                     n_flow_config,                              &
                     m_enh_sheet_config,                         &
                     m_enh_shelf_config,                         &
                     use_analytical_GL_flux_config,              &
                     DIVA_choice_ice_margin_config,              &
                     DIVA_norm_dUV_tol_config,                   &
                     DIVA_max_outer_loops_config,                &
                     DIVA_max_residual_UV_config,                &
                     DIVA_SOR_omega_config,                      &
                     DIVA_max_inner_loops_config,                &
                     DIVA_beta_max_config,                       &
                     DIVA_err_lim_config,                        &
                     DIVA_vel_max_config,                        &
                     DIVA_vel_min_config,                        &
                     DIVA_visc_iter_relax_config,                &
                     DIVA_boundary_BC_u_west_config,             &
                     DIVA_boundary_BC_u_east_config,             &
                     DIVA_boundary_BC_u_south_config,            &
                     DIVA_boundary_BC_u_north_config,            &
                     DIVA_boundary_BC_v_west_config,             &
                     DIVA_boundary_BC_v_east_config,             &
                     DIVA_boundary_BC_v_south_config,            &
                     DIVA_boundary_BC_v_north_config,            &
                     choice_timestepping_config,                 &
                     choice_ice_integration_method_config,       &
                     dHi_SOR_omega_config,                       &
                     dHi_SOR_max_res_config,                     &
                     dHi_SOR_max_loops_config,                   &
                     pc_epsilon_config,                          &
                     pc_k_I_config,                              &
                     pc_k_p_config,                              &
                     pc_eta_min_config,                          &
                     pc_max_timestep_iterations_config,          &
                     pc_redo_tol_config,                         &
                     dt_min_config,                              &
                     no_sliding_config,                          &
                     choice_sliding_law_config,                  &
                     C_sliding_config,                           &
                     m_sliding_config,                           &
                     choice_calving_law_config,                  &
                     calving_threshold_thickness_config,         &
                     choice_geothermal_heat_flux_config,         &
                     constant_geothermal_heat_flux_config,       &
                     filename_geothermal_heat_flux_config,       &
                     filename_PD_obs_climate_config,             &
                     choice_climate_matrix_config,               &
                     filename_GCM_snapshot_PI_config,            &
                     filename_GCM_snapshot_LGM_config,           &
                     choice_ocean_temperature_model_config,      &
                     ocean_temperature_PD_config,                &
                     ocean_temperature_cold_config,              &
                     ocean_temperature_warm_config,              &
                     constant_lapserate_config,                  &
                     choice_forcing_method_config,               &
                     dT_deepwater_averaging_window_config,       &
                     dT_deepwater_dT_surf_ratio_config,          &
                     d18O_dT_deepwater_ratio_config,             &
                     dT_glob_inverse_averaging_window_config,    &
                     inverse_d18O_to_dT_glob_scaling_config,     &
                     CO2_inverse_averaging_window_config,        &
                     inverse_d18O_to_CO2_scaling_config,         &
                     inverse_d18O_to_CO2_initial_CO2_config,     &
                     C_abl_constant_NAM_config,                  &
                     C_abl_constant_EAS_config,                  &
                     C_abl_constant_GRL_config,                  &
                     C_abl_constant_ANT_config,                  &
                     C_abl_Ts_NAM_config,                        &
                     C_abl_Ts_EAS_config,                        &
                     C_abl_Ts_GRL_config,                        &
                     C_abl_Ts_ANT_config,                        &
                     C_abl_Q_NAM_config,                         &
                     C_abl_Q_EAS_config,                         &
                     C_abl_Q_GRL_config,                         &
                     C_abl_Q_ANT_config,                         &
                     C_refr_NAM_config,                          &
                     C_refr_EAS_config,                          &
                     C_refr_GRL_config,                          &
                     C_refr_ANT_config,                          &
                     T_ocean_mean_PD_NAM_config,                 &
                     T_ocean_mean_PD_EAS_config,                 &
                     T_ocean_mean_PD_GRL_config,                 &
                     T_ocean_mean_PD_ANT_config,                 &
                     T_ocean_mean_cold_NAM_config,               &
                     T_ocean_mean_cold_EAS_config,               &
                     T_ocean_mean_cold_GRL_config,               &
                     T_ocean_mean_cold_ANT_config,               &
                     T_ocean_mean_warm_NAM_config,               &
                     T_ocean_mean_warm_EAS_config,               &
                     T_ocean_mean_warm_GRL_config,               &
                     T_ocean_mean_warm_ANT_config,               &
                     BMB_deepocean_PD_NAM_config,                &
                     BMB_deepocean_PD_EAS_config,                &
                     BMB_deepocean_PD_GRL_config,                &
                     BMB_deepocean_PD_ANT_config,                &
                     BMB_deepocean_cold_NAM_config,              &
                     BMB_deepocean_cold_EAS_config,              &
                     BMB_deepocean_cold_GRL_config,              &
                     BMB_deepocean_cold_ANT_config,              &
                     BMB_deepocean_warm_NAM_config,              &
                     BMB_deepocean_warm_EAS_config,              &
                     BMB_deepocean_warm_GRL_config,              &
                     BMB_deepocean_warm_ANT_config,              &
                     BMB_shelf_exposed_PD_NAM_config,            &
                     BMB_shelf_exposed_PD_EAS_config,            &
                     BMB_shelf_exposed_PD_GRL_config,            &
                     BMB_shelf_exposed_PD_ANT_config,            &
                     BMB_shelf_exposed_cold_NAM_config,          &
                     BMB_shelf_exposed_cold_EAS_config,          &
                     BMB_shelf_exposed_cold_GRL_config,          &
                     BMB_shelf_exposed_cold_ANT_config,          &
                     BMB_shelf_exposed_warm_NAM_config,          &
                     BMB_shelf_exposed_warm_EAS_config,          &
                     BMB_shelf_exposed_warm_GRL_config,          &
                     BMB_shelf_exposed_warm_ANT_config,          &
                     subshelf_melt_factor_NAM_config,            &
                     subshelf_melt_factor_EAS_config,            &
                     subshelf_melt_factor_GRL_config,            &
                     subshelf_melt_factor_ANT_config,            &
                     deep_ocean_threshold_depth_NAM_config,      &
                     deep_ocean_threshold_depth_EAS_config,      &
                     deep_ocean_threshold_depth_GRL_config,      &
                     deep_ocean_threshold_depth_ANT_config,      &
                     do_ocean_floodfill_config,                  &
                     choice_sealevel_model_config,               &
                     fixed_sealevel_config,                      &
                     filename_sealevel_record_config,            &
                     sealevel_record_length_config,              &
                     choice_GIA_model_config,                    &
                     dx_GIA_config,                              &
                     ELRA_lithosphere_flex_rigidity_config,      &
                     ELRA_bedrock_relaxation_time_config,        &
                     ELRA_mantle_density_config,                 &
                     SELEN_run_at_t_start_config,                &
                     SELEN_n_TDOF_iterations_config,             &
                     SELEN_n_recursion_iterations_config,        &
                     SELEN_use_rotational_feedback_config,       &
                     SELEN_n_harmonics_config,                   &
                     SELEN_display_progress_config,              &
                     SELEN_dir_config,                           &
                     SELEN_global_topo_filename_config,          &
                     SELEN_TABOO_init_filename_config,           &
                     SELEN_LMJ_VALUES_filename_config,           &
                     SELEN_irreg_time_n_config,                  &
                     SELEN_irreg_time_window_config,             &
                     SELEN_lith_thickness_config,                &
                     SELEN_visc_n_config,                        &
                     SELEN_visc_prof_config,                     &
                     SELEN_TABOO_CDE_config,                     & 
                     SELEN_TABOO_TLOVE_config,                   & 
                     SELEN_TABOO_DEG1_config,                    & 
                     SELEN_TABOO_RCMB_config,                    &
                     help_field_01_config,                       &
                     help_field_02_config,                       &
                     help_field_03_config,                       &
                     help_field_04_config,                       &
                     help_field_05_config,                       &
                     help_field_06_config,                       &
                     help_field_07_config,                       &
                     help_field_08_config,                       &
                     help_field_09_config,                       &
                     help_field_10_config,                       &
                     help_field_11_config,                       &
                     help_field_12_config,                       &
                     help_field_13_config,                       &
                     help_field_14_config,                       &
                     help_field_15_config,                       &
                     help_field_16_config,                       &
                     help_field_17_config,                       &
                     help_field_18_config,                       &
                     help_field_19_config,                       &
                     help_field_20_config,                       &
                     help_field_21_config,                       &
                     help_field_22_config,                       &
                     help_field_23_config,                       &
                     help_field_24_config,                       &
                     help_field_25_config,                       &
                     help_field_26_config,                       &
                     help_field_27_config,                       &
                     help_field_28_config,                       &
                     help_field_29_config,                       &
                     help_field_30_config,                       &
                     help_field_31_config,                       &
                     help_field_32_config,                       &
                     help_field_33_config,                       &
                     help_field_34_config,                       &
                     help_field_35_config,                       &
                     help_field_36_config,                       &
                     help_field_37_config,                       &
                     help_field_38_config,                       &
                     help_field_39_config,                       &
                     help_field_40_config,                       &
                     help_field_41_config,                       &
                     help_field_42_config,                       &
                     help_field_43_config,                       &
                     help_field_44_config,                       &
                     help_field_45_config,                       &
                     help_field_46_config,                       &
                     help_field_47_config,                       &
                     help_field_48_config,                       &
                     help_field_49_config,                       &
                     help_field_50_config
                      
     IF (config_filename == '') RETURN
      
     OPEN(UNIT=config_unit, FILE=TRIM(config_filename), STATUS='OLD', ACTION='READ', iostat=ios)
     IF(ios /= 0) THEN
       WRITE(UNIT=*, FMT='(/3A/)') ' ERROR: Could not open the configuration file: ', TRIM(config_filename)
       STOP
     END IF

     ! In the following statement the entire configuration file is read, using the namelist (NML=CONFIG)
     READ(UNIT=config_unit, NML=CONFIG, IOSTAT=ios)
     CLOSE(UNIT=config_unit)

     IF(ios /= 0) THEN
       WRITE(UNIT=*, FMT='(/3A)') ' ERROR while reading configuration file: ', TRIM(config_filename)
       STOP
     END IF     

  END SUBROUTINE read_config_file

  SUBROUTINE initialise_config
    ! Overwrite the values in the fields of the "C" type with the values
    ! of the "_config" variables, some which by now have had their default
    ! values overwritten by the values specified in the external config file.
    
    IMPLICIT NONE
    
   ! Time steps and range
   !=====================
   
    C%start_time_of_run                   = start_time_of_run_config
    C%end_time_of_run                     = end_time_of_run_config
    C%dt_coupling                         = dt_coupling_config
    C%dt_max                              = dt_max_config
    C%dt_thermo                           = dt_thermo_config
    C%dt_climate                          = dt_climate_config
    C%dt_SMB                              = dt_SMB_config
    C%dt_BMB                              = dt_BMB_config
    C%dt_bedrock_ELRA                     = dt_bedrock_ELRA_config
    C%dt_SELEN                            = dt_SELEN_config
    C%dt_output                           = dt_output_config
    
    ! Which ice sheets do we simulate?
    ! ================================
    
    C%do_NAM                              = do_NAM_config
    C%do_EAS                              = do_EAS_config
    C%do_GRL                              = do_GRL_config
    C%do_ANT                              = do_ANT_config
    
    ! Benchmark experiments
    ! =====================
    
    C%do_benchmark_experiment             = do_benchmark_experiment_config
    C%choice_benchmark_experiment         = choice_benchmark_experiment_config
    C%SSA_icestream_m                     = SSA_icestream_m_config
    C%ISMIP_HOM_L                         = ISMIP_HOM_L_config
    C%ISMIP_HOM_E_Arolla_filename         = ISMIP_HOM_E_Arolla_filename_config
    
    ! Whether or not to let IMAU_ICE dynamically create its own output folder
   ! =======================================================================
   
    C%create_new_output_dir               = create_new_output_dir_config
    C%output_dir                          = output_dir_config
    C%do_write_debug_data                 = do_write_debug_data_config
  
    ! Whether or not the simulation is a restart of a previous simulation
    ! ===================================================================
    
    C%is_restart                          = is_restart_config
    C%time_to_restart_from                = time_to_restart_from_config
  
    ! Horizontal grid spacing and size for the four regions
    ! =====================================================
    
    C%dx_NAM                              = dx_NAM_config
    C%dx_EAS                              = dx_EAS_config
    C%dx_GRL                              = dx_GRL_config
    C%dx_ANT                              = dx_ANT_config

    ! Scaled vertical coordinate zeta  
    ! ===============================
    
    C%nz                                  = nz_config
    ALLOCATE( C%zeta( C%nz))
    C%zeta                                = zeta_config( 1:C%nz)

    ! Input data file paths
    ! =====================
    
    C%filename_init_NAM                   = filename_init_NAM_config
    C%filename_init_EAS                   = filename_init_EAS_config
    C%filename_init_GRL                   = filename_init_GRL_config
    C%filename_init_ANT                   = filename_init_ANT_config
    
    C%filename_PD_NAM                     = filename_PD_NAM_config
    C%filename_PD_EAS                     = filename_PD_EAS_config
    C%filename_PD_GRL                     = filename_PD_GRL_config
    C%filename_PD_ANT                     = filename_PD_ANT_config
    
    C%filename_insolation                 = filename_insolation_config
    
    C%filename_CO2_record                 = filename_CO2_record_config
    C%CO2_record_length                   = CO2_record_length_config
    C%filename_d18O_record                = filename_d18O_record_config
    C%d18O_record_length                  = d18O_record_length_config
    
    C%filename_ICE5G_PD                   = filename_ICE5G_PD_config
    C%filename_ICE5G_LGM                  = filename_ICE5G_LGM_config

    ! Ice dynamics - velocity
    ! =======================
    
    C%choice_ice_dynamics                 = choice_ice_dynamics_config
    C%n_flow                              = n_flow_config
    C%m_enh_sheet                         = m_enh_sheet_config
    C%m_enh_shelf                         = m_enh_shelf_config
    C%use_analytical_GL_flux              = use_analytical_GL_flux_config
    
    ! Some parameters for numerically solving the SSA/DIVA
    C%DIVA_choice_ice_margin              = DIVA_choice_ice_margin_config
    C%DIVA_norm_dUV_tol                   = DIVA_norm_dUV_tol_config
    C%DIVA_max_outer_loops                = DIVA_max_outer_loops_config
    C%DIVA_max_residual_UV                = DIVA_max_residual_UV_config
    C%DIVA_SOR_omega                      = DIVA_SOR_omega_config
    C%DIVA_max_inner_loops                = DIVA_max_inner_loops_config
    C%DIVA_beta_max                       = DIVA_beta_max_config
    C%DIVA_err_lim                        = DIVA_err_lim_config
    C%DIVA_vel_max                        = DIVA_vel_max_config
    C%DIVA_vel_min                        = DIVA_vel_min_config
    C%DIVA_visc_iter_relax                = DIVA_visc_iter_relax_config
    C%DIVA_boundary_BC_u_west             = DIVA_boundary_BC_u_west_config
    C%DIVA_boundary_BC_u_east             = DIVA_boundary_BC_u_east_config
    C%DIVA_boundary_BC_u_south            = DIVA_boundary_BC_u_south_config
    C%DIVA_boundary_BC_u_north            = DIVA_boundary_BC_u_north_config
    C%DIVA_boundary_BC_v_west             = DIVA_boundary_BC_v_west_config
    C%DIVA_boundary_BC_v_east             = DIVA_boundary_BC_v_east_config
    C%DIVA_boundary_BC_v_south            = DIVA_boundary_BC_v_south_config
    C%DIVA_boundary_BC_v_north            = DIVA_boundary_BC_v_north_config
  
    ! Ice dynamics - time integration
    ! ===============================
    
    C%choice_timestepping                 = choice_timestepping_config
    C%choice_ice_integration_method       = choice_ice_integration_method_config
    C%dHi_SOR_omega                       = dHi_SOR_omega_config
    C%dHi_SOR_max_res                     = dHi_SOR_max_res_config
    C%dHi_SOR_max_loops                   = dHi_SOR_max_loops_config
    
    ! Predictor-corrector ice-thickness update
    C%pc_epsilon                          = pc_epsilon_config
    C%pc_k_I                              = pc_k_I_config
    C%pc_k_p                              = pc_k_p_config
    C%pc_eta_min                          = pc_eta_min_config
    C%pc_max_timestep_iterations          = pc_max_timestep_iterations_config
    C%pc_redo_tol                         = pc_redo_tol_config
    C%dt_min                              = dt_min_config
  
    ! Ice dynamics - sliding
    ! ======================
    
    C%no_sliding                          = no_sliding_config
    C%choice_sliding_law                  = choice_sliding_law_config
    C%C_sliding                           = C_sliding_config
    C%m_sliding                           = m_sliding_config
  
    ! Ice dynamics - calving
    ! ======================
    
    C%choice_calving_law                  = choice_calving_law_config
    C%calving_threshold_thickness         = calving_threshold_thickness_config
  
    ! Thermodynamics
    ! ==============
    
    C%choice_geothermal_heat_flux         = choice_geothermal_heat_flux_config
    C%constant_geothermal_heat_flux       = constant_geothermal_heat_flux_config
    C%filename_geothermal_heat_flux       = filename_geothermal_heat_flux_config
    
    ! Climate matrix
    ! ==============
    
    C%filename_PD_obs_climate             = filename_PD_obs_climate_config
    C%choice_climate_matrix               = choice_climate_matrix_config
    C%filename_GCM_snapshot_PI            = filename_GCM_snapshot_PI_config
    C%filename_GCM_snapshot_LGM           = filename_GCM_snapshot_LGM_config
    
    C%choice_ocean_temperature_model      = choice_ocean_temperature_model_config
    C%ocean_temperature_PD                = ocean_temperature_PD_config
    C%ocean_temperature_cold              = ocean_temperature_cold_config
    C%ocean_temperature_warm              = ocean_temperature_warm_config
    
    C%constant_lapserate                  = constant_lapserate_config
    
    ! Forcing
    ! =======
    
    C%choice_forcing_method               = choice_forcing_method_config
    
    C%dT_deepwater_averaging_window       = dT_deepwater_averaging_window_config
    C%dT_deepwater_dT_surf_ratio          = dT_deepwater_dT_surf_ratio_config
    C%d18O_dT_deepwater_ratio             = d18O_dT_deepwater_ratio_config
    
    C%dT_glob_inverse_averaging_window    = dT_glob_inverse_averaging_window_config
    C%inverse_d18O_to_dT_glob_scaling     = inverse_d18O_to_dT_glob_scaling_config
    C%CO2_inverse_averaging_window        = CO2_inverse_averaging_window_config
    C%inverse_d18O_to_CO2_scaling         = inverse_d18O_to_CO2_scaling_config
    C%inverse_d18O_to_CO2_initial_CO2     = inverse_d18O_to_CO2_initial_CO2_config
    
    ! SMB melt tuning
    ! ===============
    
    C%C_abl_constant_NAM                  = C_abl_constant_NAM_config
    C%C_abl_constant_EAS                  = C_abl_constant_EAS_config
    C%C_abl_constant_GRL                  = C_abl_constant_GRL_config
    C%C_abl_constant_ANT                  = C_abl_constant_ANT_config
    C%C_abl_Ts_NAM                        = C_abl_Ts_NAM_config
    C%C_abl_Ts_EAS                        = C_abl_Ts_EAS_config
    C%C_abl_Ts_GRL                        = C_abl_Ts_GRL_config
    C%C_abl_Ts_ANT                        = C_abl_Ts_ANT_config
    C%C_abl_Q_NAM                         = C_abl_Q_NAM_config
    C%C_abl_Q_EAS                         = C_abl_Q_EAS_config
    C%C_abl_Q_GRL                         = C_abl_Q_GRL_config
    C%C_abl_Q_ANT                         = C_abl_Q_ANT_config
    C%C_refr_NAM                          = C_refr_NAM_config
    C%C_refr_EAS                          = C_refr_EAS_config
    C%C_refr_GRL                          = C_refr_GRL_config
    C%C_refr_ANT                          = C_refr_ANT_config
    
    ! Sub-shelf melt parameterisation
    ! ===============================
    
    C%T_ocean_mean_PD_NAM                 = T_ocean_mean_PD_NAM_config
    C%T_ocean_mean_PD_EAS                 = T_ocean_mean_PD_EAS_config
    C%T_ocean_mean_PD_GRL                 = T_ocean_mean_PD_GRL_config
    C%T_ocean_mean_PD_ANT                 = T_ocean_mean_PD_ANT_config
    C%T_ocean_mean_cold_NAM               = T_ocean_mean_cold_NAM_config
    C%T_ocean_mean_cold_EAS               = T_ocean_mean_cold_EAS_config
    C%T_ocean_mean_cold_GRL               = T_ocean_mean_cold_GRL_config
    C%T_ocean_mean_cold_ANT               = T_ocean_mean_cold_ANT_config
    C%T_ocean_mean_warm_NAM               = T_ocean_mean_warm_NAM_config
    C%T_ocean_mean_warm_EAS               = T_ocean_mean_warm_EAS_config
    C%T_ocean_mean_warm_GRL               = T_ocean_mean_warm_GRL_config
    C%T_ocean_mean_warm_ANT               = T_ocean_mean_warm_ANT_config
    
    C%BMB_deepocean_PD_NAM                = BMB_deepocean_PD_NAM_config
    C%BMB_deepocean_PD_EAS                = BMB_deepocean_PD_EAS_config
    C%BMB_deepocean_PD_GRL                = BMB_deepocean_PD_GRL_config
    C%BMB_deepocean_PD_ANT                = BMB_deepocean_PD_ANT_config
    C%BMB_deepocean_cold_NAM              = BMB_deepocean_cold_NAM_config
    C%BMB_deepocean_cold_EAS              = BMB_deepocean_cold_EAS_config
    C%BMB_deepocean_cold_GRL              = BMB_deepocean_cold_GRL_config
    C%BMB_deepocean_cold_ANT              = BMB_deepocean_cold_ANT_config
    C%BMB_deepocean_warm_NAM              = BMB_deepocean_warm_NAM_config
    C%BMB_deepocean_warm_EAS              = BMB_deepocean_warm_EAS_config
    C%BMB_deepocean_warm_GRL              = BMB_deepocean_warm_GRL_config
    C%BMB_deepocean_warm_ANT              = BMB_deepocean_warm_ANT_config
    
    C%BMB_shelf_exposed_PD_NAM            = BMB_shelf_exposed_PD_NAM_config
    C%BMB_shelf_exposed_PD_EAS            = BMB_shelf_exposed_PD_EAS_config
    C%BMB_shelf_exposed_PD_GRL            = BMB_shelf_exposed_PD_GRL_config
    C%BMB_shelf_exposed_PD_ANT            = BMB_shelf_exposed_PD_ANT_config
    C%BMB_shelf_exposed_cold_NAM          = BMB_shelf_exposed_cold_NAM_config
    C%BMB_shelf_exposed_cold_EAS          = BMB_shelf_exposed_cold_EAS_config
    C%BMB_shelf_exposed_cold_GRL          = BMB_shelf_exposed_cold_GRL_config
    C%BMB_shelf_exposed_warm_NAM          = BMB_shelf_exposed_warm_NAM_config
    C%BMB_shelf_exposed_cold_ANT          = BMB_shelf_exposed_cold_ANT_config
    C%BMB_shelf_exposed_warm_EAS          = BMB_shelf_exposed_warm_EAS_config
    C%BMB_shelf_exposed_warm_GRL          = BMB_shelf_exposed_warm_GRL_config
    C%BMB_shelf_exposed_warm_ANT          = BMB_shelf_exposed_warm_ANT_config
    
    C%subshelf_melt_factor_NAM            = subshelf_melt_factor_NAM_config
    C%subshelf_melt_factor_EAS            = subshelf_melt_factor_EAS_config
    C%subshelf_melt_factor_GRL            = subshelf_melt_factor_GRL_config
    C%subshelf_melt_factor_ANT            = subshelf_melt_factor_ANT_config
    
    C%deep_ocean_threshold_depth_NAM      = deep_ocean_threshold_depth_NAM_config
    C%deep_ocean_threshold_depth_EAS      = deep_ocean_threshold_depth_EAS_config
    C%deep_ocean_threshold_depth_GRL      = deep_ocean_threshold_depth_GRL_config
    C%deep_ocean_threshold_depth_ANT      = deep_ocean_threshold_depth_ANT_config
  
    ! Sea level and GIA
    ! =================
    
    C%do_ocean_floodfill                  = do_ocean_floodfill_config
    C%choice_sealevel_model               = choice_sealevel_model_config
    C%fixed_sealevel                      = fixed_sealevel_config
    C%filename_sealevel_record            = filename_sealevel_record_config
    C%sealevel_record_length              = sealevel_record_length_config
  
    C%choice_GIA_model                    = choice_GIA_model_config
    C%dx_GIA                              = dx_GIA_config
    C%ELRA_lithosphere_flex_rigidity      = ELRA_lithosphere_flex_rigidity_config
    C%ELRA_bedrock_relaxation_time        = ELRA_bedrock_relaxation_time_config
    C%ELRA_mantle_density                 = ELRA_mantle_density_config

    ! SELEN
    ! =====
    
    C%SELEN_run_at_t_start                = SELEN_run_at_t_start_config
    C%SELEN_n_TDOF_iterations             = SELEN_n_TDOF_iterations_config
    C%SELEN_n_recursion_iterations        = SELEN_n_recursion_iterations_config
    C%SELEN_use_rotational_feedback       = SELEN_use_rotational_feedback_config
    C%SELEN_n_harmonics                   = SELEN_n_harmonics_config
    C%SELEN_display_progress              = SELEN_display_progress_config
    
    C%SELEN_dir                           = SELEN_dir_config
    C%SELEN_global_topo_filename          = SELEN_global_topo_filename_config
    C%SELEN_TABOO_init_filename           = SELEN_TABOO_init_filename_config
    C%SELEN_LMJ_VALUES_filename           = SELEN_LMJ_VALUES_filename_config
    
    C%SELEN_irreg_time_n                  = SELEN_irreg_time_n_config
    ALLOCATE( C%SELEN_irreg_time_window( C%SELEN_irreg_time_n))
    C%SELEN_irreg_time_window             = SELEN_irreg_time_window_config( 1:C%SELEN_irreg_time_n)
    
    C%SELEN_lith_thickness                = SELEN_lith_thickness_config
    C%SELEN_visc_n                        = SELEN_visc_n_config
    ALLOCATE( C%SELEN_visc_prof( C%SELEN_visc_n))
    C%SELEN_visc_prof = SELEN_visc_prof_config( 1:C%SELEN_visc_n)
    
    C%SELEN_TABOO_CDE                     = SELEN_TABOO_CDE_config
    C%SELEN_TABOO_TLOVE                   = SELEN_TABOO_TLOVE_config
    C%SELEN_TABOO_DEG1                    = SELEN_TABOO_DEG1_config
    C%SELEN_TABOO_RCMB                    = SELEN_TABOO_RCMB_config
    
    ! Fill in some derived values
    C%SELEN_jmax        = (C%SELEN_n_harmonics + 1) * (C%SELEN_n_harmonics + 2) / 2
    C%SELEN_reg_time_n  = MAX(1, INT(SUM(C%SELEN_irreg_time_window( 1:C%SELEN_irreg_time_n)))) * INT(1000. / C%dt_SELEN)
    
    CALL initialize_TABOO_config
  
    ! Which data fields will be written to the help_fields output file
    ! ================================================================
    
    C%help_field_01                       = help_field_01_config
    C%help_field_02                       = help_field_02_config
    C%help_field_03                       = help_field_03_config
    C%help_field_04                       = help_field_04_config
    C%help_field_05                       = help_field_05_config
    C%help_field_06                       = help_field_06_config
    C%help_field_07                       = help_field_07_config
    C%help_field_08                       = help_field_08_config
    C%help_field_09                       = help_field_09_config
    C%help_field_10                       = help_field_10_config
    C%help_field_11                       = help_field_11_config
    C%help_field_12                       = help_field_12_config
    C%help_field_13                       = help_field_13_config
    C%help_field_14                       = help_field_14_config
    C%help_field_15                       = help_field_15_config
    C%help_field_16                       = help_field_16_config
    C%help_field_17                       = help_field_17_config
    C%help_field_18                       = help_field_18_config
    C%help_field_19                       = help_field_19_config
    C%help_field_20                       = help_field_20_config
    C%help_field_21                       = help_field_21_config
    C%help_field_22                       = help_field_22_config
    C%help_field_23                       = help_field_23_config
    C%help_field_24                       = help_field_24_config
    C%help_field_25                       = help_field_25_config
    C%help_field_26                       = help_field_26_config
    C%help_field_27                       = help_field_27_config
    C%help_field_28                       = help_field_28_config
    C%help_field_29                       = help_field_29_config
    C%help_field_30                       = help_field_30_config
    C%help_field_31                       = help_field_31_config
    C%help_field_32                       = help_field_32_config
    C%help_field_33                       = help_field_33_config
    C%help_field_34                       = help_field_34_config
    C%help_field_35                       = help_field_35_config
    C%help_field_36                       = help_field_36_config
    C%help_field_37                       = help_field_37_config
    C%help_field_38                       = help_field_38_config
    C%help_field_39                       = help_field_39_config
    C%help_field_40                       = help_field_40_config
    C%help_field_41                       = help_field_41_config
    C%help_field_42                       = help_field_42_config
    C%help_field_43                       = help_field_43_config
    C%help_field_44                       = help_field_44_config
    C%help_field_45                       = help_field_45_config
    C%help_field_46                       = help_field_46_config
    C%help_field_47                       = help_field_47_config
    C%help_field_48                       = help_field_48_config
    C%help_field_49                       = help_field_49_config
    C%help_field_50                       = help_field_50_config
    
    ! Values to be filled into the total mask (used only for diagnostic output)
    ! ==========================================================================

    C%type_land                           = 0
    C%type_ocean                          = 1
    C%type_lake                           = 2
    C%type_sheet                          = 3
    C%type_shelf                          = 4
    C%type_coast                          = 5
    C%type_margin                         = 6
    C%type_groundingline                  = 7
    C%type_calvingfront                   = 8
    
    ! Parameters of the polar stereographic projections of the four model regions
    ! (These have to match the values used to create the input files!)
    ! ===========================================================================  
  
    C%lambda_M_NAM                        = 265._dp
    C%lambda_M_EAS                        = 40._dp
    C%lambda_M_GRL                        = 320._dp
    C%lambda_M_ANT                        = 0._dp
    C%phi_M_NAM                           = 62._dp
    C%phi_M_EAS                           = 70._dp
    C%phi_M_GRL                           = 72._dp
    C%phi_M_ANT                           = -90._dp
    C%alpha_stereo_NAM                    = 165.0923_dp
    C%alpha_stereo_EAS                    = 165.04_dp
    C%alpha_stereo_GRL                    = 164.85_dp
    C%alpha_stereo_ANT                    = 165.0263_dp

  END SUBROUTINE initialise_config
  
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

  SUBROUTINE create_output_dir

    IMPLICIT NONE

    CHARACTER(20)              :: output_folder_name

    INTEGER,    DIMENSION(8)   :: values
    LOGICAL                    :: ex

    CALL date_and_time(VALUES=values)

    ! Get proper year (assume we're still in the 21st century...)
    output_folder_name(1:10) = 'results_20'
    SELECT CASE( FLOOR(REAL(values(1))/10._dp)-200)
     CASE(0)
     output_folder_name(11:11) = '0'
     CASE(1)
     output_folder_name(11:11) = '1'
     CASE(2)
     output_folder_name(11:11) = '2'
     CASE(3)
     output_folder_name(11:11) = '3'
     CASE(4)
     output_folder_name(11:11) = '4'
     CASE(5)
     output_folder_name(11:11) = '5'
     CASE(6)
     output_folder_name(11:11) = '6'
     CASE(7)
     output_folder_name(11:11) = '7'
     CASE(8)
     output_folder_name(11:11) = '8'
     CASE(9)
     output_folder_name(11:11) = '9'
     CASE DEFAULT
     WRITE(0,*) 'make_output_folder: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( MOD(values(1),10))
     CASE(0)
     output_folder_name(12:12) = '0'
     CASE(1)
     output_folder_name(12:12) = '1'
     CASE(2)
     output_folder_name(12:12) = '2'
     CASE(3)
     output_folder_name(12:12) = '3'
     CASE(4)
     output_folder_name(12:12) = '4'
     CASE(5)
     output_folder_name(12:12) = '5'
     CASE(6)
     output_folder_name(12:12) = '6'
     CASE(7)
     output_folder_name(12:12) = '7'
     CASE(8)
     output_folder_name(12:12) = '8'
     CASE(9)
     output_folder_name(12:12) = '9'
     CASE DEFAULT
     WRITE(0,*) 'make_output_folder: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( values(2))
     CASE(1)
     output_folder_name(13:14) = '01'
     CASE(2)
     output_folder_name(13:14) = '02'
     CASE(3)
     output_folder_name(13:14) = '03'
     CASE(4)
     output_folder_name(13:14) = '04'
     CASE(5)
     output_folder_name(13:14) = '05'
     CASE(6)
     output_folder_name(13:14) = '06'
     CASE(7)
     output_folder_name(13:14) = '07'
     CASE(8)
     output_folder_name(13:14) = '08'
     CASE(9)
     output_folder_name(13:14) = '09'
     CASE(10)
     output_folder_name(13:14) = '10'
     CASE(11)
     output_folder_name(13:14) = '11'
     CASE(12)
     output_folder_name(13:14) = '12'
     CASE DEFAULT
     WRITE(0,*) 'make_output_folder: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( FLOOR(REAL(values(3))/10._dp))
     CASE(0)
     output_folder_name(15:15) = '0'
     CASE(1)
     output_folder_name(15:15) = '1'
     CASE(2)
     output_folder_name(15:15) = '2'
     CASE(3)
     output_folder_name(15:15) = '3'
     CASE DEFAULT
     WRITE(0,*) 'make_output_folder: ERROR retrieving date and time!'
    END SELECT

    SELECT CASE( MOD(values(3),10))
     CASE(0)
     output_folder_name(16:16) = '0'
     CASE(1)
     output_folder_name(16:16) = '1'
     CASE(2)
     output_folder_name(16:16) = '2'
     CASE(3)
     output_folder_name(16:16) = '3'
     CASE(4)
     output_folder_name(16:16) = '4'
     CASE(5)
     output_folder_name(16:16) = '5'
     CASE(6)
     output_folder_name(16:16) = '6'
     CASE(7)
     output_folder_name(16:16) = '7'
     CASE(8)
     output_folder_name(16:16) = '8'
     CASE(9)
     output_folder_name(16:16) = '9'
     CASE DEFAULT
     WRITE(0,*) 'make_output_folder: ERROR retrieving date and time!'
    END SELECT

    output_folder_name(17:20) = '_001'

    INQUIRE( FILE=TRIM(output_folder_name)//'/.', EXIST=ex )

    DO WHILE (ex)

     IF      (output_folder_name(20:20) == '0') THEN
      output_folder_name(20:20) = '1'
     ELSE IF (output_folder_name(20:20) == '1') THEN
      output_folder_name(20:20) = '2'
     ELSE IF (output_folder_name(20:20) == '2') THEN
      output_folder_name(20:20) = '3'
     ELSE IF (output_folder_name(20:20) == '3') THEN
      output_folder_name(20:20) = '4'
     ELSE IF (output_folder_name(20:20) == '4') THEN
      output_folder_name(20:20) = '5'
     ELSE IF (output_folder_name(20:20) == '5') THEN
      output_folder_name(20:20) = '6'
     ELSE IF (output_folder_name(20:20) == '6') THEN
      output_folder_name(20:20) = '7'
     ELSE IF (output_folder_name(20:20) == '7') THEN
      output_folder_name(20:20) = '8'
     ELSE IF (output_folder_name(20:20) == '8') THEN
      output_folder_name(20:20) = '9'
     ELSE IF (output_folder_name(20:20) == '9') THEN
      output_folder_name(20:20) = '0'

      IF      (output_folder_name(19:19) == '0') THEN
       output_folder_name(19:19) = '1'
      ELSE IF (output_folder_name(19:19) == '1') THEN
       output_folder_name(19:19) = '2'
      ELSE IF (output_folder_name(19:19) == '2') THEN
       output_folder_name(19:19) = '3'
      ELSE IF (output_folder_name(19:19) == '3') THEN
       output_folder_name(19:19) = '4'
      ELSE IF (output_folder_name(19:19) == '4') THEN
       output_folder_name(19:19) = '5'
      ELSE IF (output_folder_name(19:19) == '5') THEN
       output_folder_name(19:19) = '6'
      ELSE IF (output_folder_name(19:19) == '6') THEN
       output_folder_name(19:19) = '7'
      ELSE IF (output_folder_name(19:19) == '7') THEN
       output_folder_name(19:19) = '8'
      ELSE IF (output_folder_name(19:19) == '8') THEN
       output_folder_name(19:19) = '9'
      ELSE IF (output_folder_name(19:19) == '9') THEN
       output_folder_name(19:19) = '0'

       IF      (output_folder_name(18:18) == '0') THEN
        output_folder_name(18:18) = '1'
       ELSE IF (output_folder_name(18:18) == '1') THEN
        output_folder_name(18:18) = '2'
       ELSE IF (output_folder_name(18:18) == '2') THEN
        output_folder_name(18:18) = '3'
       ELSE IF (output_folder_name(18:18) == '3') THEN
        output_folder_name(18:18) = '4'
       ELSE IF (output_folder_name(18:18) == '4') THEN
        output_folder_name(18:18) = '5'
       ELSE IF (output_folder_name(18:18) == '5') THEN
        output_folder_name(18:18) = '6'
       ELSE IF (output_folder_name(18:18) == '6') THEN
        output_folder_name(18:18) = '7'
       ELSE IF (output_folder_name(18:18) == '7') THEN
        output_folder_name(18:18) = '8'
       ELSE IF (output_folder_name(18:18) == '8') THEN
        output_folder_name(18:18) = '9'
       ELSE IF (output_folder_name(18:18) == '9') THEN
        output_folder_name(18:18) = '0'
       END IF

      END IF

     END IF

     INQUIRE( FILE=TRIM(output_folder_name)//'/.', EXIST=ex )

    END DO

    IF (C%create_new_output_dir) THEN
      C%output_dir = TRIM(output_folder_name)
    END IF
    
    C%output_dir = TRIM(C%output_dir(1:255)) // '/'
    CALL system('mkdir ' // TRIM(C%output_dir))

  END SUBROUTINE create_output_dir
  
  SUBROUTINE initialise_zeta_discretisation
  
    IMPLICIT NONE

    ! Local variables:
    INTEGER :: k

    ALLOCATE(C%a_k(2:C%NZ  ))
    ALLOCATE(C%b_k(1:C%NZ-1))
    ALLOCATE(C%c_k(3:C%NZ  ))
    ALLOCATE(C%d_k(1:C%NZ-2))
    ALLOCATE(C%a_zeta(2:C%NZ-1))
    ALLOCATE(C%b_zeta(2:C%NZ-1))
    ALLOCATE(C%c_zeta(2:C%NZ-1))
    ALLOCATE(C%a_zetazeta(2:C%NZ-1))
    ALLOCATE(C%b_zetazeta(2:C%NZ-1))
    ALLOCATE(C%c_zetazeta(2:C%NZ-1))

    ALLOCATE(C%z_zeta_minus(3:C%NZ))
    ALLOCATE(C%a_zeta_minus(3:C%NZ))
    ALLOCATE(C%b_zeta_minus(3:C%NZ))
    ALLOCATE(C%b_zeta_plus(1:C%NZ-2))
    ALLOCATE(C%c_zeta_plus(1:C%NZ-2))
    ALLOCATE(C%d_zeta_plus(1:C%NZ-2))

    DO k = 2, C%NZ
     C%a_k(k) = C%zeta(k)   - C%zeta(k-1)
    END DO
    DO k = 1, C%NZ-1
     C%b_k(k) = C%zeta(k+1) - C%zeta(k)
    END DO
    DO k = 3, C%NZ
     C%c_k(k) = C%zeta(k)   - C%zeta(k-2)
    END DO
    DO k = 1, C%NZ-2
     C%d_k(k) = C%zeta(k+2) - C%zeta(k)
    END DO

    DO k = 2, C%NZ-1
     C%a_zeta(k)     =           - C%b_k(k)  / (C%a_k(k) * (C%a_k(k) + C%b_k(k)))
     C%b_zeta(k)     = (C%b_k(k) - C%a_k(k)) / (C%a_k(k) *  C%b_k(k)            )
     C%c_zeta(k)     =             C%a_k(k)  / (C%b_k(k) * (C%a_k(k) + C%b_k(k)))
     C%a_zetazeta(k) =             2.0_dp    / (C%a_k(k) * (C%a_k(k) + C%b_k(k)))
     C%b_zetazeta(k) =           - 2.0_dp    / (C%a_k(k) *  C%b_k(k)            )
     C%c_zetazeta(k) =             2.0_dp    / (C%b_k(k) * (C%a_k(k) + C%b_k(k)))
    END DO

    ! Not all of these are in use:
    DO k = 1, C%NZ-2
      C%b_zeta_plus(k) = -(C%b_k(k) + C%d_k(k)) / (C%b_k(k) *  C%d_k(k)            )
      C%c_zeta_plus(k) =              C%d_k(k)  / (C%b_k(k) * (C%d_k(k) - C%b_k(k)))
      C%d_zeta_plus(k) =              C%b_k(k)  / (C%d_k(k) * (C%b_k(k) - C%d_k(k)))
    END DO

    ! Not all of these are in use:
    DO k = 3, C%NZ
      C%z_zeta_minus(k) =             C%a_k(k)  / (C%c_k(k) * (C%c_k(k) - C%a_k(k)))
      C%a_zeta_minus(k) =             C%c_k(k)  / (C%a_k(k) * (C%a_k(k) - C%c_k(k)))
      C%b_zeta_minus(k) = (C%a_k(k) + C%c_k(k)) / (C%a_k(k) *  C%c_k(k)            )
    END DO
    
  END SUBROUTINE initialise_zeta_discretisation

END MODULE configuration_module