MODULE configuration_module
  
  ! The way it's done right now:
  ! Each config variable has two versions: one with the "_config" extension, which is
  ! an actual variable in this module only, and one without the extension, which is
  ! a field in the "C" type. THe "_config" variables are used to create a NAMELIST,
  ! which makes reading an external config file really easy - anything in the file that
  ! matches a variable in the namelist overwrites the default value. After that's done,
  ! the fields in the "C" type are replaced with the values of the "_config" variables,
  ! which now have either the default values, or those specified in the external config
  ! file.
  !
  ! While this is certainly very convenient when running the model, it does make adding
  ! new config parameters a bit tedious - you have to add the "_config" variable, add it
  ! as a field in the "C" type, add it to the namelist, and let the "C" type field be
  ! overwritten in the end.
  !
  ! Some day I'll figure out a more elegant solution for this...

  IMPLICIT NONE
  
  INTEGER, PARAMETER                :: dp  = KIND(1.0D0)  ! Kind of double precision numbers. Reals should be declared as: REAL(dp) :: example

  ! ===================================================================================
  ! "_config  variables, which will be collected into a NAMELIST, and possibly replaced
  ! by the values in the external config file. Remember the "_config" extension!
  ! ===================================================================================

  ! Time steps and range
  ! ====================
  
  REAL(dp) :: start_time_of_run_config   = 0.0_dp       ! Start time (in years) of the simulations
  REAL(dp) :: end_time_of_run_config     = 50000.0_dp   ! End   time (in years) of the simulations
  REAL(dp) :: dt_coupling_config         = 100._dp      ! Interval of coupling (in years) between the four ice-sheets  
  REAL(dp) :: dt_max_config              = 10.0_dp      ! Maximum time step (in years) of the ice model
  REAL(dp) :: dt_thermo_config           = 10.0_dp      ! Time step (in years) for updating thermodynamics
  REAL(dp) :: dt_climate_config          = 10._dp       ! Time step (in years) for updating the climate
  REAL(dp) :: dt_SMB_config              = 10._dp       ! Time step (in years) for updating the SMB
  REAL(dp) :: dt_BMB_config              = 10._dp       ! Time step (in years) for updating the BMB
  REAL(dp) :: dt_bedrock_ELRA_config     = 100._dp      ! Time step (in years) for updating the bedrock deformation rate with the ELRA model
  REAL(dp) :: dt_output_config           = 5000.0_dp    ! Time step (in years) for writing output
  
  ! Which ice sheets do we simulate?
  ! ================================
  
  LOGICAL :: do_NAM_config               = .FALSE.      ! North America
  LOGICAL :: do_EAS_config               = .FALSE.      ! Eurasia
  LOGICAL :: do_GRL_config               = .FALSE.      ! Greenland
  LOGICAL :: do_ANT_config               = .TRUE.       ! Antarctica
  
  ! Benchmark experiments
  ! =====================
  
  LOGICAL            :: do_benchmark_experiment_config          = .TRUE.
  CHARACTER(LEN=256) :: choice_benchmark_experiment_config      = 'EISMINT_I'
  REAL(dp)           :: m_SSA_icestream_config                  = 1   ! Values tested by Schoof are 1, 10, and 20
  
  ! Whether or not to let IMAU_ICE dynamically create its own output folder.
  ! This works fine locally, on LISA its better to use a fixed folder name.
  ! =======================================================================
  
  LOGICAL            :: create_new_output_dir_config = .TRUE.
  CHARACTER(LEN=256) :: output_dir_config            = 'results_IMAU_ICE'
  LOGICAL            :: do_write_dummy_output_config = .FALSE.
  
  ! Whether or not the simulation is a restart of a previous simulation
  ! ===================================================================
  
  LOGICAL            :: is_restart_config           = .FALSE.
  REAL(dp)           :: time_to_restart_from_config = 0._dp     ! Can be different from C%start_time_of_run, though this will issue a warning
  
  ! Horizontal grid spacing and size for the four regions
  ! =====================================================
  
  REAL(dp)           :: dx_NAM_config   = 40000._dp
  REAL(dp)           :: dx_EAS_config   = 40000._dp
  REAL(dp)           :: dx_GRL_config   = 20000._dp
  REAL(dp)           :: dx_ANT_config   = 40000._dp

  ! The scaled vertical coordinate zeta, used mainly in thermodynamics
  ! ==================================================================
  
  INTEGER                        :: nz_config   = 15
  REAL(dp), DIMENSION(210), SAVE :: zeta_config = &
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
  CHARACTER(LEN=256) :: filename_init_NAM_config = '/Users/berends/Documents/Datasets/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256) :: filename_init_EAS_config = '/Users/berends/Documents/Datasets/ETOPO1/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256) :: filename_init_GRL_config = '/Users/berends/Documents/Datasets/Bedmachine_Greenland/Greenland_BedMachine_5km_noEllesmere.nc'
  CHARACTER(LEN=256) :: filename_init_ANT_config = '/Users/berends/Documents/Datasets/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc'
  
  ! PD reference data (NetCDF)
  CHARACTER(LEN=256) :: filename_PD_NAM_config   = '/Users/berends/Documents/Datasets/ETOPO1/NorthAmerica_ETOPO1_5km.nc'
  CHARACTER(LEN=256) :: filename_PD_EAS_config   = '/Users/berends/Documents/Datasets/ETOPO1/Eurasia_ETOPO1_5km.nc'
  CHARACTER(LEN=256) :: filename_PD_GRL_config   = '/Users/berends/Documents/Datasets/Bedmachine_Greenland/Greenland_BedMachine_5km_noEllesmere.nc'
  CHARACTER(LEN=256) :: filename_PD_ANT_config   = '/Users/berends/Documents/Datasets/Bedmachine_Antarctica/Bedmachine_v1_Antarctica_5km.nc' 
   
  ! Insolation forcing (NetCDF) (Laskar et al., 2004)
  CHARACTER(LEN=256) :: filename_insolation_config = '/Users/berends/Documents/Datasets/Insolation_laskar/Insolation_Laskar_etal_2004.nc'
  
  ! CO2 record (ASCII text file, so the number of rows needs to be specified)
  CHARACTER(LEN=256) :: filename_CO2_record_config = '/Users/berends/Documents/Datasets/CO2/EPICA_CO2_Bereiter_2015_100yr.dat'
  INTEGER            :: CO2_record_length_config   = 8001
  
  ! d18O record (ASCII text file, so the number of rows needs to be specified)
  CHARACTER(LEN=256) :: filename_d18O_record_config = '/Users/berends/Documents/Datasets/d18O/Ahn2017_d18O.dat'
  INTEGER            :: d18O_record_length_config   = 2051

  ! Ice dynamics and thermodynamics
  ! ===============================
  REAL(dp)           :: m_enh_sia_config                        = 1.0_dp                 ! Ice flow enhancement factor in the SIA (5.0 in ANICE)
  REAL(dp)           :: m_enh_ssa_config                        = 1.0_dp                 ! Ice flow enhancement factor in the SSA (0.7 in ANICE)
  CHARACTER(LEN=256) :: choice_sliding_law_config               = 'Coulomb_regularised'  ! Choice of sliding law (currently only "Coulomb_regularised" is implemented)
  REAL(dp)           :: C_sliding_config                        = 1.0E7_dp               ! Factor   in Weertman sliding law
  REAL(dp)           :: m_sliding_config                        = 1._dp/3._dp            ! Exponent in Weertman sliding law
  LOGICAL            :: use_analytical_GL_flux_config           = .FALSE.                ! Whether or not the analytical grounding line flux solution is used
  REAL(dp)           :: geothermal_heat_flux_config             = 1.72E06_dp             ! Geothermal Heat flux [J m^-2 yr^-1] Sclater et al. (1980)
  CHARACTER(LEN=256) :: choice_calving_law_config               = 'threshold_thickness'  ! Choice of calving law (currently only "threshold_thickness" is implemented)
  REAL(dp)           :: calving_threshold_thickness_config      = 200._dp                ! Threshold ice thickness in the "threshold_thickness" calving law (200m taken from ANICE)
  
  ! Some parameters for numerically solving the SSA
  REAL(dp)           :: SSA_RN_tol_config                       = 1E-5_dp                ! Successive solutions of the effective viscosity iteration must not differ by more than this amount (norm(N-Nprev) / norm(N) < RN_tol)
  REAL(dp)           :: SSA_norm_dUV_tol_config                 = 0.5_dp                 ! Successive solutions of UV in the effective viscosity iteration must not differ by more than this amount (on average)
  INTEGER            :: SSA_max_outer_loops_config              = 50                     ! Maximum number of effective viscosity iterations
  REAL(dp)           :: SSA_max_grad_N_config                   = 315.5694_dp            ! = 1e-5 * sec_per_year according to Frank Pattyn. If set to zero, the solver becomes identical to the old ANICE solver.
  REAL(dp)           :: SSA_max_residual_UV_config              = 2.5_dp                 ! The maximum residual in U and V in the SOR solver of the linearised SSA must drop below this value [m yr^-1]
  REAL(dp)           :: SSA_SOR_omega_config                    = 1.3_dp                 ! The over-relaxation parameter in the SOR solver of the linearised SSA
  INTEGER            :: SSA_max_inner_loops_config              = 10000                  ! Maximum number of iterations in the SOR solver of the linearised SSA
  
  ! Sea level and GIA
  ! =================
  
  LOGICAL            :: do_ocean_floodfill_config               = .TRUE.                 ! Use a flood-fill to determine the ocean mask, so that (pro-/sub-glacial) lakes dont exist
  CHARACTER(LEN=256) :: choice_sealevel_model_config            = 'eustatic'             ! Can be "fixed", "prescribed", "eustatic", or "SELEN"
  REAL(dp)           :: fixed_sealevel_config                   = 0._dp
  CHARACTER(LEN=256) :: filename_sealevel_record_config         = 'name_of_file.dat'
  INTEGER            :: sealevel_record_length_config           = 1
  
  CHARACTER(LEN=256) :: choice_GIA_model_config                 = 'ELRA'                 ! Can be "none", "ELRA", or "SELEN"
  REAL(dp)           :: ELRA_lithosphere_flex_rigidity_config   = 1.0E+25_dp             ! Lithospheric flexural rigidity [kg m^2 s^-2]
  REAL(dp)           :: ELRA_bedrock_relaxation_time_config     = 3000.0_dp              ! Relaxation time for bedrock adjustment [yr]
  REAL(dp)           :: mantle_density_config                   = 3300.0_dp              ! Mantle density [kg m^-3]
  
  ! Climate matrix
  ! ==============
  
  ! Present-day observed climate (ERA40) (NetCDF)
  CHARACTER(LEN=256) :: filename_PD_obs_climate_config          = '/Users/berends/Documents/Datasets/ERA40/ERA40_climate_global.nc'
  
  ! GCM snapshots
  CHARACTER(LEN=256) :: choice_climate_matrix_config            = 'PI_LGM'               ! 'PI_LGM' uses 2 snapshots
  CHARACTER(LEN=256) :: filename_GCM_snapshot_PI_config         = 'GCM_snapshots/PI_Valdes.nc'
  CHARACTER(LEN=256) :: filename_GCM_snapshot_LGM_config        = 'GCM_snapshots/LGM_Valdes.nc'
  
  ! Ocean temperature (used for both thermodynamics and basal melt)
  CHARACTER(LEN=256) :: choice_ocean_temperature_model_config   = 'scaled'               ! Can be "fixed" (use PD value) or "scaled" (scale between "PD", "warm", and "cold" values based on forcing (prescribed or inverse-modelled))
  REAL(dp)           :: ocean_temperature_PD_config             = 271.46_dp              ! present day temperature of the ocean beneath the shelves [K; -1.7 Celsius]
  REAL(dp)           :: ocean_temperature_cold_config           = 268.16_dp              ! cold period temperature of the ocean beneath the shelves [K; -5.0 Celcius]
  REAL(dp)           :: ocean_temperature_warm_config           = 275.16_dp              ! warm period temperature of the ocean beneath the shelves [K;  2.0 Celcius]
  
  REAL(dp)           :: constant_lapserate_config               = -0.008_dp              ! Constant atmospheric lapse rate [K m^-1]
  
  ! Forcing
  ! =======
  
  ! The choice of forcing:
  ! 'd18O_inverse_dT_glob' : Use the inverse routine with the specified d18O record to calculate a global temperature offset (e.g. de Boer et al., 2013)
  ! 'CO2_direct'           : Use the specified CO2 record to force the climate matrix (e.g. Berends et al., 2018)
  ! 'd18O_inverse_CO2'     : Use the inverse routine with the specified d18O record to calculate CO2 and then force the climate matrix (e.g. Berends et al., 2019)
  CHARACTER(LEN=256) :: choice_forcing_method_config            = 'd18O_inverse_dT_glob'
  
  REAL(dp)           :: dT_deepwater_averaging_window_config    = 3000      ! Time window (in yr) over which global mean temperature anomaly is averaged to find the deep-water temperature anomaly
  REAL(dp)           :: dT_deepwater_dT_surf_ratio_config       = 0.25_dp   ! Ratio between global mean surface temperature change and deep-water temperature change
  REAL(dp)           :: d18O_dT_deepwater_ratio_config          = -0.28_dp  ! Ratio between deep-water temperature change and benthic d18O change
  
  REAL(dp)           :: dT_glob_inverse_averaging_window_config = 2000._dp      ! Time window (in yr) over which global mean temperature anomaly is averaged before changing it with the inverse routine
  REAL(dp)           :: inverse_d18O_to_dT_glob_scaling_config  = 20._dp        ! Scaling factor between modelled d18O anomaly and prescribed temperature anomaly change (value from de Boer et al., 2013)
  REAL(dp)           :: CO2_inverse_averaging_window_config     = 2000._dp      ! Time window (in yr) over which CO2                             is averaged before changing it with the inverse routine
  REAL(dp)           :: inverse_d18O_to_CO2_scaling_config      = 68._dp        ! Scaling factor between modelled d18O anomaly and modelled CO2 change (value from Berends et al., 2019)
  REAL(dp)           :: inverse_d18O_to_CO2_initial_CO2_config  = 280._dp       ! CO2 value at the start of the simulation when using the inverse method to calculate CO2
  
  ! SMB tuning
  ! ==========
  
  REAL(dp)           :: C_abl_constant_NAM_config = -49._dp       ! 34._dp    (commented values are old ANICE defaults, but since refreezing was not calculated right
  REAL(dp)           :: C_abl_constant_EAS_config = -49._dp       !            and this has since been fixed, these values will still not give the same results as
  REAL(dp)           :: C_abl_constant_GRL_config = -49._dp       !            they used to in ANICE.)
  REAL(dp)           :: C_abl_constant_ANT_config = -49._dp
  REAL(dp)           :: C_abl_Ts_NAM_config       = 10._dp        ! 10._dp
  REAL(dp)           :: C_abl_Ts_EAS_config       = 10._dp
  REAL(dp)           :: C_abl_Ts_GRL_config       = 10._dp
  REAL(dp)           :: C_abl_Ts_ANT_config       = 10._dp
  REAL(dp)           :: C_abl_Q_NAM_config        = 0.0227_dp     ! 0.513_dp
  REAL(dp)           :: C_abl_Q_EAS_config        = 0.0227_dp
  REAL(dp)           :: C_abl_Q_GRL_config        = 0.0227_dp
  REAL(dp)           :: C_abl_Q_ANT_config        = 0.0227_dp
  REAL(dp)           :: C_refr_NAM_config         = 0.051_dp      ! 0.012_dp
  REAL(dp)           :: C_refr_EAS_config         = 0.051_dp 
  REAL(dp)           :: C_refr_GRL_config         = 0.051_dp 
  REAL(dp)           :: C_refr_ANT_config         = 0.051_dp
  
  ! Sub-shelf melt parameterisation
  ! ===============================
            
  REAL(dp)           :: BMB_deepocean_PD_config       =  -5._dp       ! Present-day value for deep-ocean areas [m/year]
  REAL(dp)           :: BMB_deepocean_cold_config     =  -2._dp       ! Cold period value for deep-ocean areas [m/year]
  REAL(dp)           :: BMB_deepocean_warm_config     = -10._dp       ! Warm period value for deep-ocean areas [m/year]    

  REAL(dp)           :: BMB_shelf_exposed_PD_config   =  -3._dp       ! Present-day value for exposed areas    [m/year]
  REAL(dp)           :: BMB_shelf_exposed_cold_config =  -0._dp       ! Cold period value for exposed areas    [m/year]
  REAL(dp)           :: BMB_shelf_exposed_warm_config =  -6._dp       ! Warm period value for exposed areas    [m/year]
    
  REAL(dp)           :: subshelf_melt_factor_config = 0.005_dp       ! Overall tuning factor for sub-shelf melt


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
    REAL(dp)                            :: m_SSA_icestream

    ! Whether or not to let IMAU_ICE dynamically create its own output folder
   ! =======================================================================
   
    LOGICAL                  :: create_new_output_dir
    CHARACTER(LEN=256)       :: output_dir
    LOGICAL                  :: do_write_dummy_output
  
    ! Whether or not the simulation is a restart of a previous simulation
    ! ===================================================================
   
    LOGICAL                  :: is_restart
    REAL(dp)                 :: time_to_restart_from
  
    ! Horizontal grid spacing and size for the four regions
    ! =====================================================
    
    REAL(dp)           :: dx_NAM
    REAL(dp)           :: dx_EAS
    REAL(dp)           :: dx_GRL
    REAL(dp)           :: dx_ANT

    ! Scaled vertical coordinate zeta  
    ! ===============================
     
    INTEGER                             :: nz       ! Number of grid points in vertical direction for thermodynamics in ice sheet
    REAL(dp), DIMENSION(:), ALLOCATABLE :: zeta

    ! Input data file paths
    ! =====================
                            
    CHARACTER(LEN=256)       :: filename_init_NAM
    CHARACTER(LEN=256)       :: filename_init_EAS
    CHARACTER(LEN=256)       :: filename_init_GRL
    CHARACTER(LEN=256)       :: filename_init_ANT
    
    CHARACTER(LEN=256)       :: filename_PD_NAM
    CHARACTER(LEN=256)       :: filename_PD_EAS
    CHARACTER(LEN=256)       :: filename_PD_GRL
    CHARACTER(LEN=256)       :: filename_PD_ANT
    
    CHARACTER(LEN=256)       :: filename_insolation
    
    CHARACTER(LEN=256)       :: filename_CO2_record
    INTEGER                  :: CO2_record_length
    CHARACTER(LEN=256)       :: filename_d18O_record
    INTEGER                  :: d18O_record_length

    ! Ice dynamics & thermodynamics
    ! =============================
    
    REAL(dp)                 :: m_enh_sia
    REAL(dp)                 :: m_enh_ssa
    CHARACTER(LEN=256)       :: choice_sliding_law
    REAL(dp)                 :: C_sliding
    REAL(dp)                 :: m_sliding
    LOGICAL                  :: use_analytical_GL_flux
    REAL(dp)                 :: geothermal_heat_flux
    CHARACTER(LEN=256)       :: choice_calving_law
    REAL(dp)                 :: calving_threshold_thickness
  
    ! Some parameters for numerically solving the SSA
    REAL(dp)                 :: SSA_RN_tol
    REAL(dp)                 :: SSA_norm_dUV_tol
    INTEGER                  :: SSA_max_outer_loops
    REAL(dp)                 :: SSA_max_grad_N
    REAL(dp)                 :: SSA_max_residual_UV
    REAL(dp)                 :: SSA_SOR_omega
    INTEGER                  :: SSA_max_inner_loops
  
    ! Sea level and GIA
    ! =================
    
    LOGICAL                  :: do_ocean_floodfill
    CHARACTER(LEN=256)       :: choice_sealevel_model
    REAL(dp)                 :: fixed_sealevel
    CHARACTER(LEN=256)       :: filename_sealevel_record
    INTEGER                  :: sealevel_record_length
  
    CHARACTER(LEN=256)       :: choice_GIA_model
    REAL(dp)                 :: ELRA_lithosphere_flex_rigidity
    REAL(dp)                 :: ELRA_bedrock_relaxation_time
    REAL(dp)                 :: mantle_density
    
    ! Climate matrix
    ! ==============
    
    CHARACTER(LEN=256)       :: filename_PD_obs_climate
    CHARACTER(LEN=256)       :: choice_climate_matrix
    CHARACTER(LEN=256)       :: filename_GCM_snapshot_PI
    CHARACTER(LEN=256)       :: filename_GCM_snapshot_LGM
    
    CHARACTER(LEN=256)       :: choice_ocean_temperature_model
    REAL(dp)                 :: ocean_temperature_PD
    REAL(dp)                 :: ocean_temperature_cold
    REAL(dp)                 :: ocean_temperature_warm
    
    REAL(dp)                 :: constant_lapserate
    
    ! Forcing
    ! =======
    
    CHARACTER(LEN=256)       :: choice_forcing_method
    
    REAL(dp)                 :: dT_deepwater_averaging_window
    REAL(dp)                 :: dT_deepwater_dT_surf_ratio
    REAL(dp)                 :: d18O_dT_deepwater_ratio
    
    REAL(dp)                 :: dT_glob_inverse_averaging_window
    REAL(dp)                 :: inverse_d18O_to_dT_glob_scaling
    REAL(dp)                 :: CO2_inverse_averaging_window
    REAL(dp)                 :: inverse_d18O_to_CO2_scaling
    REAL(dp)                 :: inverse_d18O_to_CO2_initial_CO2
    
    ! SMB melt tuning
    ! ===============
    
    REAL(dp)                 :: C_abl_constant_NAM
    REAL(dp)                 :: C_abl_constant_EAS
    REAL(dp)                 :: C_abl_constant_GRL
    REAL(dp)                 :: C_abl_constant_ANT
    REAL(dp)                 :: C_abl_Ts_NAM
    REAL(dp)                 :: C_abl_Ts_EAS
    REAL(dp)                 :: C_abl_Ts_GRL
    REAL(dp)                 :: C_abl_Ts_ANT
    REAL(dp)                 :: C_abl_Q_NAM
    REAL(dp)                 :: C_abl_Q_EAS
    REAL(dp)                 :: C_abl_Q_GRL
    REAL(dp)                 :: C_abl_Q_ANT
    REAL(dp)                 :: C_refr_NAM
    REAL(dp)                 :: C_refr_EAS
    REAL(dp)                 :: C_refr_GRL
    REAL(dp)                 :: C_refr_ANT
  
    ! Sub-shelf melt parameterisation
    ! ===============================
              
    REAL(dp)                 :: BMB_deepocean_PD
    REAL(dp)                 :: BMB_deepocean_cold
    REAL(dp)                 :: BMB_deepocean_warm 
  
    REAL(dp)                 :: BMB_shelf_exposed_PD
    REAL(dp)                 :: BMB_shelf_exposed_cold
    REAL(dp)                 :: BMB_shelf_exposed_warm
      
    REAL(dp)                 :: subshelf_melt_factor
    
    ! Values to be filled into the total mask (used only for diagnostic output)
    ! ==========================================================================
    
    INTEGER                  :: type_land
    INTEGER                  :: type_ocean
    INTEGER                  :: type_lake
    INTEGER                  :: type_sheet
    INTEGER                  :: type_shelf
    INTEGER                  :: type_coast
    INTEGER                  :: type_margin
    INTEGER                  :: type_groundingline
    INTEGER                  :: type_calvingfront
    
    ! Parameters of the polar stereographic projections of the four model regions
    ! (These have to match the values used to create the input files!)
    ! ===========================================================================     
                          
    REAL(dp)                 :: lambda_M_NAM                           
    REAL(dp)                 :: lambda_M_EAS                           
    REAL(dp)                 :: lambda_M_GRL                           
    REAL(dp)                 :: lambda_M_ANT
    REAL(dp)                 :: phi_M_NAM
    REAL(dp)                 :: phi_M_EAS
    REAL(dp)                 :: phi_M_GRL
    REAL(dp)                 :: phi_M_ANT
    REAL(dp)                 :: alpha_stereo_NAM
    REAL(dp)                 :: alpha_stereo_EAS
    REAL(dp)                 :: alpha_stereo_GRL
    REAL(dp)                 :: alpha_stereo_ANT
    
    
    ! Some useful constants for the scaled vertical coordinate transformation
    ! =======================================================================

    REAL(dp), DIMENSION(:), ALLOCATABLE        :: a_k
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: b_k
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: c_k
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: d_k
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: a_zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: b_zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: c_zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: a_zetazeta
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: b_zetazeta
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: c_zetazeta

    REAL(dp), DIMENSION(:), ALLOCATABLE        :: z_zeta_minus
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: a_zeta_minus
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: b_zeta_minus

    REAL(dp), DIMENSION(:), ALLOCATABLE        :: b_zeta_plus
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: c_zeta_plus
    REAL(dp), DIMENSION(:), ALLOCATABLE        :: d_zeta_plus

  END TYPE constants_type


  ! ===============================================
  ! "C" is an instance of the "constants_type" type
  ! ===============================================
  
  TYPE(constants_type), SAVE :: C


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
                     dt_output_config,                           &
                     do_NAM_config,                              &
                     do_EAS_config,                              &
                     do_GRL_config,                              &
                     do_ANT_config,                              &
                     do_benchmark_experiment_config,             &
                     choice_benchmark_experiment_config,         &
                     m_SSA_icestream_config,                     &
                     create_new_output_dir_config,               &
                     output_dir_config,                          &
                     do_write_dummy_output_config,               &
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
                     m_enh_sia_config,                           &
                     m_enh_ssa_config,                           &
                     choice_sliding_law_config,                  &
                     C_sliding_config,                           &
                     m_sliding_config,                           &
                     use_analytical_GL_flux_config,              &
                     geothermal_heat_flux_config,                &
                     choice_calving_law_config,                  &
                     calving_threshold_thickness_config,         &
                     SSA_RN_tol_config,                          &
                     SSA_norm_dUV_tol_config,                    &
                     SSA_max_outer_loops_config,                 &
                     SSA_max_grad_N_config,                      &
                     SSA_max_residual_UV_config,                 &
                     SSA_SOR_omega_config,                       &
                     SSA_max_inner_loops_config,                 &
                     do_ocean_floodfill_config,                  &
                     choice_sealevel_model_config,               &
                     fixed_sealevel_config,                      &
                     filename_sealevel_record_config,            &
                     sealevel_record_length_config,              &
                     choice_GIA_model_config,                    &
                     ELRA_lithosphere_flex_rigidity_config,      &
                     ELRA_bedrock_relaxation_time_config,        &
                     mantle_density_config,                      &
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
                     BMB_deepocean_PD_config,                    &
                     BMB_deepocean_cold_config,                  &
                     BMB_deepocean_warm_config,                  &
                     BMB_shelf_exposed_PD_config,                &
                     BMB_shelf_exposed_cold_config,              &
                     BMB_shelf_exposed_warm_config,              &
                     subshelf_melt_factor_config
                      
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
    C%m_SSA_icestream                     = m_SSA_icestream_config
    
    ! A quick check to see if the specified benchmark experiment
    ! has actually been implemented
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
        C%do_NAM = .FALSE.
        C%do_EAS = .FALSE.
        C%do_GRL = .FALSE.
        C%do_ANT = .TRUE.
      ELSE
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialize_main_constants!'
        STOP
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Whether or not to let IMAU_ICE dynamically create its own output folder
   ! =======================================================================
   
    C%create_new_output_dir               = create_new_output_dir_config
    C%output_dir                          = output_dir_config
    C%do_write_dummy_output               = do_write_dummy_output_config
  
    ! Whether or not the simulation is a restart of a previous simulation
    ! ===================================================================
    
    C%is_restart                          = is_restart_config
    C%time_to_restart_from                = time_to_restart_from_config
  
    ! Horizontal grid spacing and size for the four regions
    ! =====================================================
    
    C%dx_NAM     = dx_NAM_config
    C%dx_EAS     = dx_EAS_config
    C%dx_GRL     = dx_GRL_config
    C%dx_ANT     = dx_ANT_config

    ! Scaled vertical coordinate zeta  
    ! ===============================
    
    C%nz     = nz_config
    ALLOCATE( C%zeta( C%nz))
    C%zeta   = zeta_config( 1:C%nz)

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

    ! Ice dynamics & thermodynamics
    ! =============================
    
    C%m_enh_sia                           = m_enh_sia_config
    C%m_enh_ssa                           = m_enh_ssa_config    
    C%choice_sliding_law                  = choice_sliding_law_config
    C%C_sliding                           = C_sliding_config
    C%m_sliding                           = m_sliding_config
    C%use_analytical_GL_flux              = use_analytical_GL_flux_config
    C%geothermal_heat_flux                = geothermal_heat_flux_config
    C%choice_calving_law                  = choice_calving_law_config
    C%calving_threshold_thickness         = calving_threshold_thickness_config
  
    ! Some parameters for numerically solving the SSA
    C%SSA_RN_tol                          = SSA_RN_tol_config
    C%SSA_norm_dUV_tol                    = SSA_norm_dUV_tol_config
    C%SSA_max_outer_loops                 = SSA_max_outer_loops_config
    C%SSA_max_grad_N                      = SSA_max_grad_N_config
    C%SSA_max_residual_UV                 = SSA_max_residual_UV_config
    C%SSA_SOR_omega                       = SSA_SOR_omega_config
    C%SSA_max_inner_loops                 = SSA_max_inner_loops_config
  
    ! Sea level and GIA
    ! =================
    
    C%do_ocean_floodfill                  = do_ocean_floodfill_config
    C%choice_sealevel_model               = choice_sealevel_model_config
    C%fixed_sealevel                      = fixed_sealevel_config
    C%filename_sealevel_record            = filename_sealevel_record_config
    C%sealevel_record_length              = sealevel_record_length_config
  
    C%choice_GIA_model                    = choice_GIA_model_config
    C%ELRA_lithosphere_flex_rigidity      = ELRA_lithosphere_flex_rigidity_config
    C%ELRA_bedrock_relaxation_time        = ELRA_bedrock_relaxation_time_config
    C%mantle_density                      = mantle_density_config
    
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
    
    C%BMB_deepocean_PD                    = BMB_deepocean_PD_config
    C%BMB_deepocean_cold                  = BMB_deepocean_cold_config
    C%BMB_deepocean_warm                  = BMB_deepocean_warm_config
    C%BMB_shelf_exposed_PD                = BMB_shelf_exposed_PD_config
    C%BMB_shelf_exposed_cold              = BMB_shelf_exposed_cold_config
    C%BMB_shelf_exposed_warm              = BMB_shelf_exposed_warm_config
    C%subshelf_melt_factor                = subshelf_melt_factor_config
    
    ! Values to be filled into the total mask (used only for diagnostic output)
    ! ==========================================================================

    C%type_land           = 0
    C%type_ocean          = 1
    C%type_lake           = 2
    C%type_sheet          = 3
    C%type_shelf          = 4
    C%type_coast          = 5
    C%type_margin         = 6
    C%type_groundingline  = 7
    C%type_calvingfront   = 8
    
    ! Parameters of the polar stereographic projections of the four model regions
    ! (These have to match the values used to create the input files!)
    ! ===========================================================================  
  
    C%lambda_M_NAM     = 265._dp
    C%lambda_M_EAS     = 40._dp
    C%lambda_M_GRL     = 320._dp
    C%lambda_M_ANT     = 0._dp
    C%phi_M_NAM        = 62._dp
    C%phi_M_EAS        = 70._dp
    C%phi_M_GRL        = 72._dp
    C%phi_M_ANT        = -90._dp
    C%alpha_stereo_NAM = 165.0923_dp
    C%alpha_stereo_EAS = 165.04_dp
    C%alpha_stereo_GRL = 164.85_dp
    C%alpha_stereo_ANT = 165.0263_dp

  END SUBROUTINE initialise_config

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
    ! Calculation of the discretization coefficients (on an Arakawa A-grid) which are used in the
    ! thermodynamics_module and in the velocity_full_stokes_module. See table 4.
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