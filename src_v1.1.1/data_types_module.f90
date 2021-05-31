MODULE data_types_module
  ! Contains all the different types for storing data. Put all together in a separate module so that
  ! all subroutines can use all types without interdependency conflicts, and also to make the 
  ! modules with the actual physics code more readable.
  ! If only Types could be collapsed in BBEdit...

  USE configuration_module,        ONLY: dp, C
  USE data_types_netcdf_module,    ONLY: type_netcdf_climate_data, type_netcdf_PD_data, type_netcdf_init_data, &
                                         type_netcdf_insolation, type_netcdf_restart, type_netcdf_help_fields, &
                                         type_netcdf_ICE5G_data, type_netcdf_debug, type_netcdf_geothermal_heat_flux

  IMPLICIT NONE
  
  TYPE type_grid
    ! The regular square grid
    
    INTEGER,                    POINTER     :: nx, ny
    REAL(dp),                   POINTER     :: dx
    REAL(dp), DIMENSION(:    ), POINTER     :: x, y
    REAL(dp),                   POINTER     :: xmin, xmax, ymin, ymax
    INTEGER :: wnx, wny, wdx, wx, wy, wxmin, wxmax, wymin, wymax
    INTEGER                                 :: i1, i2, j1, j2 ! Parallelisation by domain decomposition
    
    REAL(dp),                   POINTER     :: lambda_M
    REAL(dp),                   POINTER     :: phi_M
    REAL(dp),                   POINTER     :: alpha_stereo
    REAL(dp), DIMENSION(:,:  ), POINTER     :: lat, lon
    INTEGER :: wlambda_M, wphi_M, walpha_stereo, wlat, wlon
  
  END TYPE type_grid
  
  TYPE type_ice_model
    ! The ice dynamics sub-model data structure.
    
    ! Basic data - ice thickness, bedrock & surface elevation, sea level (geoid elevation), englacial temperature
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_Aa                  ! Ice thickness [m]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_Acx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_Acy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_Ab
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb_Aa                  ! Bedrock elevation [m w.r.t. PD sea level]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb_Acx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb_Acy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb_Ab
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_Aa                  ! Surface elevation [m w.r.t. PD sea level]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_Acx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_Acy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_Ab
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SL_Aa                  ! Sea level (geoid elevation) [m w.r.t. PD sea level]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SL_Acx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SL_Acy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SL_Ab
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti_Aa                  ! Englacial temperature [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti_Acx
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti_Acy
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti_Ab
    INTEGER :: wHi_Aa, wHi_Acx, wHi_Acy, wHi_Ab
    INTEGER :: wHb_Aa, wHb_Acx, wHb_Acy, wHb_Ab
    INTEGER :: wHs_Aa, wHs_Acx, wHs_Acy, wHs_Ab
    INTEGER :: wSL_Aa, wSL_Acx, wSL_Acy, wSL_Ab
    INTEGER :: wTi_Aa, wTi_Acx, wTi_Acy, wTi_Ab
    
    ! Ice velocities
    REAL(dp), DIMENSION(:,:  ), POINTER     :: U_vav_SIA_Aa           ! Vertically averaged ice velocity resulting from the SIA [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V_vav_SIA_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: U_vav_SIA_Acx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V_vav_SIA_Acy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: U_vav_SSA_Aa           ! Vertically averaged ice velocity resulting from the SSA [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V_vav_SSA_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: U_vav_SSA_Acx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V_vav_SSA_Acy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: U_vav_Aa               ! Vertically averaged ice velocity [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V_vav_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: U_vav_Acx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V_vav_Acy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: U_surf_Aa              ! Ice velocity at the surface [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V_surf_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: U_base_Aa              ! Ice velocity at the base [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V_base_Aa
    REAL(dp), DIMENSION(:,:,:), POINTER     :: U_3D_Aa                ! 3-D ice velocity [m yr^-1]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: V_3D_Aa
    REAL(dp), DIMENSION(:,:,:), POINTER     :: W_3D_Aa
    INTEGER :: wU_vav_SIA_Aa, wV_vav_SIA_Aa, wU_vav_SIA_Acx, wV_vav_SIA_Acy
    INTEGER :: wU_vav_SSA_Aa, wV_vav_SSA_Aa, wU_vav_SSA_Acx, wV_vav_SSA_Acy
    INTEGER :: wU_vav_Aa, wV_vav_Aa, wU_vav_Acx, wV_vav_Acy
    INTEGER :: wU_surf_Aa, wV_surf_Aa, wU_base_Aa, wV_base_Aa, wU_3D_Aa, wV_3D_Aa, wW_3D_Aa
    
    ! Different masks
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_land_Aa
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_ocean_Aa
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_ocean_Acx
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_ocean_Acy
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_lake_Aa
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_ice_Aa
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_sheet_Aa
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_sheet_Acx
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_sheet_Acy
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_shelf_Aa
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_shelf_Acx
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_shelf_Acy
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_coast_Aa
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_coast_Acx
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_coast_Acy
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_margin_Aa
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_margin_Acx
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_margin_Acy
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_gl_Aa
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_gl_Acx
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_gl_Acy
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_cf_Aa
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_cf_Acx
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_cf_Acy
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_Aa
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_Acx
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_Acy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: grounded_fraction
    INTEGER :: wmask_land_Aa, wmask_ocean_Aa, wmask_ocean_Acx, wmask_ocean_Acy, wmask_lake_Aa, wmask_ice_Aa
    INTEGER :: wmask_sheet_Aa, wmask_sheet_Acx, wmask_sheet_Acy, wmask_shelf_Aa, wmask_shelf_Acx, wmask_shelf_Acy
    INTEGER :: wmask_coast_Aa, wmask_coast_Acx, wmask_coast_Acy, wmask_margin_Aa, wmask_margin_Acx, wmask_margin_Acy
    INTEGER :: wmask_gl_Aa, wmask_gl_Acx, wmask_gl_Acy, wmask_cf_Aa, wmask_cf_Acx, wmask_cf_Acy, wmask_Aa, wmask_Acx, wmask_Acy
    INTEGER :: wgrounded_fraction
    
    ! Ice physical properties
    REAL(dp), DIMENSION(:,:,:), POINTER     :: A_flow_Aa      ! Flow parameter [Pa^-3 y^-1]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: A_flow_Acx
    REAL(dp), DIMENSION(:,:,:), POINTER     :: A_flow_Acy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: A_flow_mean_Aa ! Vertically averaged flow parameter [Pa^-3 y^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: A_flow_mean_Acx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: A_flow_mean_Acy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: A_flow_mean_Ab
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti_pmp_Aa      ! The pressure melting point temperature [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Cpi_Aa         ! Specific heat capacity of ice [J kg^-1 K^-1].
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ki_Aa          ! Conductivity of ice [J m^-1 K^-1 yr^-1].
    INTEGER :: wA_flow_Aa, wA_flow_Acx, wA_flow_Acy, wA_flow_mean_Aa, wA_flow_mean_Acx, wA_flow_mean_Acy, wA_flow_mean_Ab
    INTEGER :: wTi_pmp_Aa, wCpi_Aa, wKi_Aa
    
    ! Spatial and temporal derivatives
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHi_dt_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHi_dx_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHi_dy_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_dt_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_dx_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_dy_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dt_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dx_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dy_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dx_Acx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dy_Acx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dx_Acy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dy_Acy
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dTi_dx_Aa
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dTi_dy_Aa
    INTEGER :: wdHi_dt_Aa, wdHi_dx_Aa, wdHi_dy_Aa, wdHb_dt_Aa, wdHb_dx_Aa, wdHb_dy_Aa
    INTEGER :: wdHs_dt_Aa, wdHs_dx_Aa, wdHs_dy_Aa, wdHs_dx_Acx, wdHs_dy_Acx, wdHs_dx_Acy, wdHs_dy_Acy
    INTEGER :: wdTi_dx_Aa, wdTi_dy_Aa
    
    ! Zeta derivatives
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dzeta_dt_Aa
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dzeta_dx_Aa
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dzeta_dy_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dzeta_dz_Aa
    INTEGER :: wdzeta_dt_Aa, wdzeta_dx_Aa, wdzeta_dy_Aa, wdzeta_dz_Aa
        
    ! Ice dynamics - SIA
    REAL(dp), DIMENSION(:,:,:), POINTER     :: D_SIA_3D_Acx
    REAL(dp), DIMENSION(:,:,:), POINTER     :: D_SIA_3D_Acy
    REAL(dp), DIMENSION(:,:,:), POINTER     :: D_SIA_3D_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: D_SIA_Acx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: D_SIA_Acy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: D_SIA_Aa
    INTEGER :: wD_SIA_3D_Acx, wD_SIA_3D_Acy, wD_SIA_3D_Aa, wD_SIA_Acx, wD_SIA_Acy, wD_SIA_Aa
    
    ! Ice dynamics - SSA
    REAL(dp), DIMENSION(:,:  ), POINTER     :: phi_fric_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: tau_c_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Ux_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Uy_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Vx_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Vy_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: eta_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: N_Aa, N_Aa_prev
    REAL(dp), DIMENSION(:,:  ), POINTER     :: S_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dx_shelf_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dy_shelf_Aa
    REAL(dp), DIMENSION(:,:  ), POINTER     :: RHSx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: RHSy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: LHSx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: LHSy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: resU
    REAL(dp), DIMENSION(:,:  ), POINTER     :: resV
    REAL(dp), DIMENSION(:,:  ), POINTER     :: eu_ij
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ev_ij
    INTEGER :: wphi_fric_Aa, wtau_c_Aa, wUx_Aa, wUy_Aa, wVx_Aa, wVy_Aa, weta_Aa, wN_Aa, wN_Aa_prev, wS_Aa
    INTEGER :: wdHs_dx_shelf_Aa, wdHs_dy_shelf_Aa, wRHSx, wRHSy, wLHSx, wLHSy, wresU, wresV, weu_ij, wev_ij
    
    ! Ice dynamics - ice thickness calculation
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Qx_Acx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Qy_Acy
    INTEGER :: wQx_Acx, wQy_Acy
    
    ! Thermodynamics
    REAL(dp), DIMENSION(:,:  ), POINTER     :: frictional_heating_Aa   ! Friction heating due to basal sliding
    REAL(dp), DIMENSION(:,:  ), POINTER     :: GHF_Aa                  ! Geothermal heat flux
    INTEGER :: wfrictional_heating_Aa, wGHF_Aa
    
    ! Isotope content
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_Aa_prev
    REAL(dp), DIMENSION(:,:  ), POINTER     :: IsoRef
    REAL(dp), DIMENSION(:,:  ), POINTER     :: IsoSurf
    REAL(dp), DIMENSION(:,:  ), POINTER     :: IsoIce
    REAL(dp), DIMENSION(:,:  ), POINTER     :: MB_iso
    INTEGER :: wHi_Aa_prev, wIsoRef, wIsoSurf, wIsoIce, wMB_iso
    
    ! ELRA GIA model
    INTEGER,                    POINTER     :: flex_prof_rad
    REAL(dp), DIMENSION(:,:  ), POINTER     :: flex_prof
    REAL(dp), DIMENSION(:,:  ), POINTER     :: surface_load_PD
    REAL(dp), DIMENSION(:,:  ), POINTER     :: surface_load
    REAL(dp), DIMENSION(:,:  ), POINTER     :: surface_load_rel
    REAL(dp), DIMENSION(:,:  ), POINTER     :: surface_load_rel_ext
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_eq
    INTEGER :: wflex_prof_rad, wflex_prof, wsurface_load_PD, wsurface_load, wsurface_load_rel, wsurface_load_rel_ext, wdHb_eq
    
  END TYPE type_ice_model
  
  TYPE type_debug_fields
    ! Dummy variables for debugging
    
    ! NetCDF debug file
    TYPE(type_netcdf_debug)                 :: netcdf
    
    ! Grid size
    INTEGER,                    POINTER     :: nx, ny
    INTEGER :: wnx, wny
    
    ! Data
    INTEGER,  DIMENSION(:,:  ), POINTER     :: int_2D_01
    INTEGER,  DIMENSION(:,:  ), POINTER     :: int_2D_02
    INTEGER,  DIMENSION(:,:  ), POINTER     :: int_2D_03
    INTEGER,  DIMENSION(:,:  ), POINTER     :: int_2D_04
    INTEGER,  DIMENSION(:,:  ), POINTER     :: int_2D_05
    INTEGER,  DIMENSION(:,:  ), POINTER     :: int_2D_06
    INTEGER,  DIMENSION(:,:  ), POINTER     :: int_2D_07
    INTEGER,  DIMENSION(:,:  ), POINTER     :: int_2D_08
    INTEGER,  DIMENSION(:,:  ), POINTER     :: int_2D_09
    INTEGER,  DIMENSION(:,:  ), POINTER     :: int_2D_10
    INTEGER :: wint_2D_01, wint_2D_02, wint_2D_03, wint_2D_04, wint_2D_05, wint_2D_06, wint_2D_07, wint_2D_08, wint_2D_09, wint_2D_10
    
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_01
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_02
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_03
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_04
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_05
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_06
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_07
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_08
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_09
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_10
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_11
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_12
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_13
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_14
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_15
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_16
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_17
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_18
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_19
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dp_2D_20
    INTEGER :: wdp_2D_01, wdp_2D_02, wdp_2D_03, wdp_2D_04, wdp_2D_05, wdp_2D_06, wdp_2D_07, wdp_2D_08, wdp_2D_09, wdp_2D_10
    INTEGER :: wdp_2D_11, wdp_2D_12, wdp_2D_13, wdp_2D_14, wdp_2D_15, wdp_2D_16, wdp_2D_17, wdp_2D_18, wdp_2D_19, wdp_2D_20
    
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_3D_01
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_3D_02
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_3D_03
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_3D_04
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_3D_05
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_3D_06
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_3D_07
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_3D_08
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_3D_09
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_3D_10
    INTEGER :: wdp_3D_01, wdp_3D_02, wdp_3D_03, wdp_3D_04, wdp_3D_05, wdp_3D_06, wdp_3D_07, wdp_3D_08, wdp_3D_09, wdp_3D_10
    
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_2D_monthly_01
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_2D_monthly_02
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_2D_monthly_03
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_2D_monthly_04
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_2D_monthly_05
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_2D_monthly_06
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_2D_monthly_07
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_2D_monthly_08
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_2D_monthly_09
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dp_2D_monthly_10
    INTEGER :: wdp_2D_monthly_01, wdp_2D_monthly_02, wdp_2D_monthly_03, wdp_2D_monthly_04, wdp_2D_monthly_05, wdp_2D_monthly_06, wdp_2D_monthly_07, wdp_2D_monthly_08, wdp_2D_monthly_09, wdp_2D_monthly_10
    
  END TYPE type_debug_fields
  
  TYPE type_subclimate_global
    ! Global climate data, either from present-day observations (ERA40) or from a GCM snapshot.
    
    CHARACTER(LEN=256)                      :: name                          ! 'ERA40', 'HadCM3_PI', etc.
    
    ! NetCDF file containing the data
    TYPE(type_netcdf_climate_data)          :: netcdf
    
    ! Grid
    INTEGER,                    POINTER     :: nlat, nlon
    REAL(dp), DIMENSION(:    ), POINTER     :: lat
    REAL(dp), DIMENSION(:    ), POINTER     :: lon
    INTEGER :: wnlat, wnlon, wlat, wlon
        
    ! General forcing info (not relevant for PD observations)
    REAL(dp),                   POINTER     :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    REAL(dp),                   POINTER     :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp),                   POINTER     :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    REAL(dp),                   POINTER     :: orbit_obl
    REAL(dp),                   POINTER     :: orbit_pre
    INTEGER :: wCO2, worbit_time, worbit_ecc, worbit_obl, worbit_pre ! MPI windows to all these memory spaces
    
    ! Actual GCM data
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T2m                           ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Precip                        ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_ref                        ! Orography     that was used to force the GCM (m w.r.t. PD sea-level)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_mask                       !  LBS: Mask of the ice sheet, to calculate absorbed insolation 
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_WE                       ! Monthly mean west-east wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_SN                       ! Monthly mean south_north wind speed (m/s)
    INTEGER :: wT2m, wPrecip, wHs_ref, wHi_mask, wWind_WE, wWind_SN ! MPI windows to all these memory spaces
    
    ! Paralelisation
    INTEGER                                 :: i1, i2                        ! Grid domain (:,i1:i2) of each process
  
  END TYPE type_subclimate_global
  
  TYPE type_ICE5G_timeframe
    ! Global ICE5G data, used to create reference orography for GCM snapshots so they can be downscaled to the model orography
    
    ! Time of the reconstruction
    REAL(dp)                                :: time        ! 0 for PD, -21000 for LGM
    
    ! NetCDF file containing the data
    TYPE(type_netcdf_ICE5G_data)            :: netcdf
    
    ! Grid
    INTEGER,                    POINTER     :: nlat, nlon
    REAL(dp), DIMENSION(:    ), POINTER     :: lat
    REAL(dp), DIMENSION(:    ), POINTER     :: lon
    INTEGER :: wnlat, wnlon, wlat, wlon
    
    ! Data
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb
    REAL(dp), DIMENSION(:,:  ), POINTER     :: mask_ice
    INTEGER :: wHi, wHb, wmask_ice
    
    ! Paralelisation
    INTEGER                                 :: i1, i2                        ! Grid domain (:,i1:i2) of each process
    
  END TYPE type_ICE5G_timeframe
  
  TYPE type_climate_matrix
    ! The climate matrix data structure. Contains all the different global GCM snapshots.
    
    ! The present-day observed climate (ERA40)
    TYPE(type_subclimate_global)            :: PD_obs
    TYPE(type_subclimate_global)            :: PD_mod
    
    ! The GCM snapshots.
    TYPE(type_subclimate_global)            :: GCM_PI
    TYPE(type_subclimate_global)            :: GCM_LGM
    
    ! ICE5G data
    TYPE(type_ICE5G_timeframe)              :: ICE5G_PD
    TYPE(type_ICE5G_timeframe)              :: ICE5G_LGM
  
  END TYPE type_climate_matrix
  
  TYPE type_subclimate_region
    ! Data from PD observations or a GCM snapshot, projected from their initial global grid onto the model region grid.
    
    CHARACTER(LEN=256)                      :: name                          ! 'ERA40', 'HadCM3_PI', etc.
        
    ! General forcing info (not relevant for PD observations)
    REAL(dp),                   POINTER     :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    REAL(dp),                   POINTER     :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp),                   POINTER     :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    REAL(dp),                   POINTER     :: orbit_obl
    REAL(dp),                   POINTER     :: orbit_pre
    REAL(dp),                   POINTER     :: sealevel
    INTEGER :: wCO2, worbit_time, worbit_ecc, worbit_obl, worbit_pre, wsealevel
    
    ! Actual observed / GCM data (read from external NetCDF file)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T2m                           ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Precip                        ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_ref                        ! Orography      that was used to force the GCM (m w.r.t. PD sea-level) (coarse GCM resolution)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi                            ! Ice thickness  that was used to force the GCM (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_mask                       ! LBS: Use a mask, rather than the actual height ! 
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb                            ! Bed topography that was used to force the GCM (m w.r.t. PD sea-level)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs                            ! Orography      that was used to force the GCM (m w.r.t. PD sea-level) (fine   ISM resolution)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_WE                       ! Monthly mean west-east   wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_SN                       ! Monthly mean south-north wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_LR                       ! Monthly mean wind speed in the x-direction (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_DU                       ! Monthly mean wind speed in the y-direction (m/s)
    INTEGER :: wT2m, wPrecip, wHs_ref, wHi, wHi_mask, wHb, wHs, wWind_WE, wWind_SN, wWind_LR, wWind_DU ! MPI windows to all these memory spaces
    
    ! Spatially variable lapse rate for GCM snapshots (see Berends et al., 2018)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: lambda
    INTEGER :: wlambda
    
    ! Reference absorbed insolation (for GCM snapshots), or insolation at model time for the applied climate
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Q_TOA                         ! Monthly mean insolation at the top of the atmosphere (W/m2) (taken from the prescribed insolation solution at orbit_time)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Albedo                        ! Monthly mean surface albedo (calculated using our own SMB scheme for consistency)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: I_abs                         ! Total yearly absorbed insolation, used in the climate matrix for interpolation
    REAL(dp),                   POINTER     :: Q_TOA_jun_65N, Q_TOA_jan_80S
    INTEGER :: wQ_TOA, wAlbedo, wI_abs, wQ_TOA_jun_65N, wQ_TOA_jan_80S
    
    ! Ocean temperature (parameterised right now)
    REAL(dp),                   POINTER     :: T_ocean_mean
    INTEGER :: wT_ocean_mean
  
  END TYPE type_subclimate_region
  
  TYPE type_climate_model
    ! All the relevant climate data fields (PD observations, GCM snapshots, and final, applied climate) on the model region grid

    TYPE(type_subclimate_region)            :: PD_obs         ! PD observations (e.g. ERA40)
    TYPE(type_subclimate_region)            :: PD_mod         ! Modeled PD (PI) climate, for calculating the bias
    TYPE(type_subclimate_region)            :: GCM_PI         ! Pre-industrial  (e.g. HadCM3, Singarayer & Valdes, 2010)
    TYPE(type_subclimate_region)            :: GCM_LGM        ! LGM             (e.g. HadCM3, Singarayer & Valdes, 2010)
    TYPE(type_subclimate_region)            :: applied        ! Final applied climate
    
    ! GCM bias
    REAL(dp), DIMENSION(:,:,:), POINTER     :: GCM_bias_T2m
    REAL(dp), DIMENSION(:,:,:), POINTER     :: GCM_bias_Precip
    INTEGER :: wGCM_bias_T2m, wGCM_bias_Precip
          
  END TYPE type_climate_model
  
  TYPE type_SMB_model
    ! The different SMB components, calculated from the prescribed climate
    
    ! Tuning parameters (different for each region, set from config)
    REAL(dp),                   POINTER     :: C_abl_constant
    REAL(dp),                   POINTER     :: C_abl_Ts
    REAL(dp),                   POINTER     :: C_abl_Q
    REAL(dp),                   POINTER     :: C_refr
    INTEGER :: wC_abl_constant, wC_abl_Ts, wC_abl_Q, wC_refr
    
    ! Data fields
    REAL(dp), DIMENSION(:,:  ), POINTER     :: AlbedoSurf                    ! Surface albedo underneath the snow layer (water, rock or ice)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: MeltPreviousYear              ! Total melt that occurred during the previous year (m)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: FirnDepth                     ! Depth of the firn layer (m)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Rainfall                      ! Monthly rainfall (m)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Snowfall                      ! Monthly snowfall (m)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: AddedFirn                     ! Monthly added firn (m)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Melt                          ! Monthly melt (m)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Refreezing                    ! Monthly refreezing (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Refreezing_year               ! Yearly  refreezing (m)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Runoff                        ! Monthly runoff (m)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Albedo                        ! Monthly albedo
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Albedo_year                   ! Yearly albedo
    REAL(dp), DIMENSION(:,:,:), POINTER     :: SMB                           ! Monthly SMB (m)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SMB_year                      ! Yearly  SMB (m)
    INTEGER :: wAlbedoSUrf, wMeltPreviousYear, wFirnDepth, wRainfall, wSnowfall, wAddedFirn, wMelt
    INTEGER :: wRefreezing, wRefreezing_year, wRunoff, wAlbedo, wAlbedo_year, wSMB, wSMB_year
  
  END TYPE type_SMB_model
  
  TYPE type_BMB_model
    ! The different BMB components
    
    ! Tuning parameters (different for each region, set from config)
    REAL(dp),                   POINTER     :: T_ocean_mean_PD
    REAL(dp),                   POINTER     :: T_ocean_mean_cold
    REAL(dp),                   POINTER     :: T_ocean_mean_warm
    REAL(dp),                   POINTER     :: BMB_deepocean_PD
    REAL(dp),                   POINTER     :: BMB_deepocean_cold
    REAL(dp),                   POINTER     :: BMB_deepocean_warm
    REAL(dp),                   POINTER     :: BMB_shelf_exposed_PD
    REAL(dp),                   POINTER     :: BMB_shelf_exposed_cold
    REAL(dp),                   POINTER     :: BMB_shelf_exposed_warm
    REAL(dp),                   POINTER     :: subshelf_melt_factor
    REAL(dp),                   POINTER     :: deep_ocean_threshold_depth
    INTEGER :: wT_ocean_mean_PD, wT_ocean_mean_cold, wT_ocean_mean_warm
    INTEGER :: wBMB_deepocean_PD, wBMB_deepocean_cold, wBMB_deepocean_warm
    INTEGER :: wBMB_shelf_exposed_PD, wBMB_shelf_exposed_cold, wBMB_shelf_exposed_warm
    INTEGER :: wsubshelf_melt_factor, wdeep_ocean_threshold_depth
    
    ! Data fields
    REAL(dp), DIMENSION(:,:  ), POINTER     :: BMB                           ! The basal mass balance (same as SMB: negative means melt)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: BMB_sheet                     ! The basal mass balance underneath the land-based ice sheet
    REAL(dp), DIMENSION(:,:  ), POINTER     :: BMB_shelf                     ! The basal mass balance underneath the floating   ice shelf
    REAL(dp), DIMENSION(:,:  ), POINTER     :: sub_angle                     ! "subtended angle"      for the sub-shelf melt parameterisation
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dist_open                     ! distance to open ocean for the sub-shelf melt parameterisation
    INTEGER :: wBMB, wBMB_sheet, wBMB_shelf, wsub_angle, wdist_open
  
  END TYPE type_BMB_model
  
  TYPE type_PD_data_fields
    ! Data structure containing data fields describing the present-day world.
    
    ! NetCDF file containing the data
    TYPE(type_netcdf_PD_data)               :: netcdf
    
    ! The grid on which the input file is defined (can be different from the model grid)
    INTEGER,                    POINTER     :: nx, ny
    REAL(dp), DIMENSION(:    ), POINTER     :: x, y 
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_raw
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb_raw
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_raw
    INTEGER :: wnx, wny, wx, wy, wHi_raw, wHb_raw, wHs_raw
    
    ! The data mapped to the model grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs
    INTEGER :: wHi, wHb, wHs
              
  END TYPE type_PD_data_fields
  
  TYPE type_init_data_fields
    ! Data structure containing data fields describing model state at the start of the simulation.
    
    ! NetCDF file containing the data
    TYPE(type_netcdf_init_data)             :: netcdf
    
    ! The grid on which the input file is defined (can be different from the model grid)
    INTEGER,                    POINTER     :: nx, ny, nz, nt
    REAL(dp), DIMENSION(:    ), POINTER     :: x, y, zeta, time
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_raw
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb_raw
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_raw
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti_raw
    REAL(dp), DIMENSION(:,:  ), POINTER     :: U_SSA_raw
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V_SSA_raw
    REAL(dp), DIMENSION(:,:,:), POINTER     :: FirnDepth_raw
    REAL(dp), DIMENSION(:,:  ), POINTER     :: MeltPreviousYear_raw
    INTEGER :: wnx, wny, wnz, wnt, wx, wy, wzeta, wtime, wHi_raw, wHb_raw, wHs_raw, wTi_raw, wU_SSA_raw, wV_SSA_raw, wFirnDepth_raw, wMeltPreviousYear_raw
    
    ! The data mapped to the model grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti
    REAL(dp), DIMENSION(:,:  ), POINTER     :: U_SSA
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V_SSA
    REAL(dp), DIMENSION(:,:,:), POINTER     :: FirnDepth
    REAL(dp), DIMENSION(:,:  ), POINTER     :: MeltPreviousYear
    INTEGER :: wHi, wHb, wHs, wTi, wU_SSA, wV_SSA, wFirnDepth, wMeltPreviousYear
              
  END TYPE type_init_data_fields
  
  TYPE type_forcing_data
    ! Data structure containing model forcing data - CO2 record, d18O record, (global) insolation record
    
    ! Data for the inverse routine
    REAL(dp),                   POINTER     :: dT_glob                                   ! Modelled global mean annual surface temperature anomaly w.r.t. PD
    REAL(dp), DIMENSION(:    ), POINTER     :: dT_glob_history                           ! Time window (length set in config) listing previous dT_glob values
    INTEGER,                    POINTER     :: ndT_glob_history                          ! Number of entries (= length of time window / dt_coupling)
    REAL(dp),                   POINTER     :: dT_deepwater                              ! Modelled deep-water temperature anomaly (= window averaged dT_glob * scaling factor), scaling factor set in config
    INTEGER :: wdT_glob, wdT_glob_history, wndT_glob_history, wdT_deepwater
    
    REAL(dp),                   POINTER     :: d18O_NAM, d18O_EAS, d18O_GRL, d18O_ANT    ! Modelled benthic d18O contributions from ice volume in the model regions
    REAL(dp),                   POINTER     :: d18O_from_ice_volume_mod                  ! Modelled benthic d18O contribution from global ice volume
    REAL(dp),                   POINTER     :: d18O_from_temperature_mod                 ! Modelled benthic d18O contribution from deep-water temperature change
    REAL(dp),                   POINTER     :: d18O_mod                                  ! Modelled benthic d18O
    INTEGER :: wd18O_NAM, wd18O_EAS, wd18O_GRL, wd18O_ANT, wd18O_from_ice_volume_mod, wd18O_from_temperature_mod, wd18O_mod
    
    REAL(dp),                   POINTER     :: dT_glob_inverse                           ! Global mean annual surface temperature anomaly resulting from the inverse method
    REAL(dp), DIMENSION(:    ), POINTER     :: dT_glob_inverse_history                   ! Time window (length set in config) listing previous dT_glob values
    INTEGER,                    POINTER     :: ndT_glob_inverse_history                  ! Number of entries (= length of time window / dt_coupling)
    REAL(dp),                   POINTER     :: CO2_inverse                               ! CO2 resulting from the inverse method
    REAL(dp), DIMENSION(:    ), POINTER     :: CO2_inverse_history                       ! Time window (length set in config) listing previous CO2_inverse values
    INTEGER,                    POINTER     :: nCO2_inverse_history                      ! Number of entries (= length of time window / dt_coupling)
    REAL(dp),                   POINTER     :: CO2_mod                                   ! Either equal to CO2_obs or to CO2_inverse, for easier writing to output.
    INTEGER :: wdT_glob_inverse, wdT_glob_inverse_history, wndT_glob_inverse_history, wCO2_inverse, wCO2_inverse_history, wnCO2_inverse_history, wCO2_mod
    
    ! External forcing: CO2 record
    REAL(dp), DIMENSION(:    ), POINTER     :: CO2_time
    REAL(dp), DIMENSION(:    ), POINTER     :: CO2_record
    REAL(dp),                   POINTER     :: CO2_obs
    INTEGER :: wCO2_time, wCO2_record, wCO2_obs
    
    ! External forcing: d18O record
    REAL(dp), DIMENSION(:    ), POINTER     :: d18O_time
    REAL(dp), DIMENSION(:    ), POINTER     :: d18O_record
    REAL(dp),                   POINTER     :: d18O_obs
    REAL(dp),                   POINTER     :: d18O_obs_PD
    INTEGER :: wd18O_time, wd18O_record, wd18O_obs, wd18O_obs_PD
    
    ! External forcing: insolation
    TYPE(type_netcdf_insolation)            :: netcdf_ins
    INTEGER,                    POINTER     :: ins_nyears
    INTEGER,                    POINTER     :: ins_nlat
    REAL(dp), DIMENSION(:    ), POINTER     :: ins_time
    REAL(dp), DIMENSION(:    ), POINTER     :: ins_lat
    REAL(dp),                   POINTER     :: ins_t0, ins_t1
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ins_Q_TOA0, ins_Q_TOA1
    INTEGER :: wins_nyears, wins_nlat, wins_time, wins_lat, wins_t0, wins_t1, wins_Q_TOA0, wins_Q_TOA1
    
    ! External forcing: geothermal heat flux
    TYPE(type_netcdf_geothermal_heat_flux)  :: netcdf_ghf
    INTEGER,                    POINTER     :: ghf_nlon
    INTEGER,                    POINTER     :: ghf_nlat
    REAL(dp), DIMENSION(:    ), POINTER     :: ghf_lon
    REAL(dp), DIMENSION(:    ), POINTER     :: ghf_lat
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ghf_ghf
    INTEGER :: wghf_nlon, wghf_nlat, wghf_lon, wghf_lat, wghf_ghf
    
  END TYPE type_forcing_data
  
  TYPE type_model_region
    ! Contains all the different data structures, organised by sub-model (ice, climate)
    
    ! Metadata
    CHARACTER(LEN=3)                        :: name                  ! NAM, EAS, GRL, ANT
    CHARACTER(LEN=256)                      :: long_name             ! North America, Eurasia, Greenland, Antarctica
    
    ! Timers and switches for determining which modules need to be called at what points in time during the simulation
    REAL(dp)                                :: time                   ! The current time of this particular region.
    REAL(dp)                                :: dt, dt_prev            ! The (dynamic) time step    
    REAL(dp)                                :: dt_SIA                 ! The critical timestep following from the SIA ice velocities
    REAL(dp)                                :: dt_SSA                 ! The critical timestep following from the SSA ice velocities    
    REAL(dp)                                :: t0_SIA,     t1_SIA     ! Time of last and next update of SIA velocity
    REAL(dp)                                :: t0_SSA,     t1_SSA     !                                 SSA velocity
    REAL(dp)                                :: t0_thermo,  t1_Thermo  !                                 thermodynamics
    REAL(dp)                                :: t0_output,  t1_output  !                                 output file
    REAL(dp)                                :: t0_climate, t1_climate !                                 climate
    REAL(dp)                                :: t0_SMB,     t1_SMB     !                                 SMB
    REAL(dp)                                :: t0_BMB,     t1_BMB     !                                 BMB
    REAL(dp)                                :: t0_ELRA,    t1_ELRA    !                                 ELRA bedrock deformation rates
    
    LOGICAL                                 :: do_solve_SIA
    LOGICAL                                 :: do_solve_SSA
    LOGICAL                                 :: do_thermodynamics
    LOGICAL                                 :: do_climate
    LOGICAL                                 :: do_SMB
    LOGICAL                                 :: do_BMB
    LOGICAL                                 :: do_write_output
    LOGICAL                                 :: do_ELRA
    
    ! The region's ice sheet's volume and volume above flotation (in mSLE, so the second one is the ice sheets GMSL contribution)
    REAL(dp), POINTER                       :: ice_area
    REAL(dp), POINTER                       :: ice_volume
    REAL(dp), POINTER                       :: ice_volume_PD
    REAL(dp), POINTER                       :: ice_volume_above_flotation
    REAL(dp), POINTER                       :: ice_volume_above_flotation_PD
    REAL(dp), POINTER                       :: GMSL_contribution
    REAL(dp), POINTER                       :: mean_isotope_content
    REAL(dp), POINTER                       :: mean_isotope_content_PD
    REAL(dp), POINTER                       :: d18O_contribution
    REAL(dp), POINTER                       :: d18O_contribution_PD
    INTEGER :: wice_area, wice_volume, wice_volume_PD, wice_volume_above_flotation, wice_volume_above_flotation_PD, wGMSL_contribution
    INTEGER :: wmean_isotope_content, wmean_isotope_content_PD, wd18O_contribution, wd18O_contribution_PD
        
    ! Reference data fields
    TYPE(type_PD_data_fields)               :: PD               ! The present-day data fields for this model region, on a high-res Cartesian grid
    TYPE(type_PD_data_fields)               :: topo             ! The topo data fields for this model region, on a high-res Cartesian grid
    TYPE(type_init_data_fields)             :: init             ! The initial     data fields for this model region, on a high-res Cartesian grid
    
    ! Mask where ice is not allowed to form (so Greenland is not included in NAM and EAS, and Ellesmere is not included in GRL)
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_noice
    INTEGER                                 :: wmask_noice
        
    ! Sub-models
    TYPE(type_grid)                         :: grid             ! This region's x,y grid
    TYPE(type_ice_model)                    :: ice              ! All the ice model data for this model region
    TYPE(type_climate_model)                :: climate          ! All the climate data for this model region
    TYPE(type_SMB_model)                    :: SMB              ! The different SMB components for this model region
    TYPE(type_BMB_model)                    :: BMB              ! The different BMB components for this model region
    
    ! Output netcdf files
    TYPE(type_netcdf_restart)               :: restart
    TYPE(type_netcdf_help_fields)           :: help_fields
    
    ! Computation times
    REAL(dp)                                :: tcomp_total
    REAL(dp)                                :: tcomp_SIA
    REAL(dp)                                :: tcomp_SSA
    REAL(dp)                                :: tcomp_thermo
    REAL(dp)                                :: tcomp_climate
    
  END TYPE type_model_region
  
CONTAINS

END MODULE data_types_module
