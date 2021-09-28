MODULE data_types_module

  ! Contains all the different TYPEs for storing data. Put all together in a separate module so that
  ! all subroutines can use all types without interdependency conflicts, and also to make the 
  ! modules with the actual physics code more readable.
  ! If only Types could be collapsed in BBEdit...

  USE configuration_module,        ONLY: dp, C
  USE data_types_netcdf_module,    ONLY: type_netcdf_climate_data, type_netcdf_PD_data, type_netcdf_init_data, &
                                         type_netcdf_insolation, type_netcdf_restart, type_netcdf_help_fields, &
                                         type_netcdf_ICE5G_data, type_netcdf_debug, type_netcdf_geothermal_heat_flux, &
                                         type_netcdf_SELEN_output, type_netcdf_SELEN_global_topo, type_netcdf_climate_forcing, &
                                         type_netcdf_scalars_global, type_netcdf_scalars_regional

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
  
  TYPE type_sparse_matrix_CSR
    ! A matrix equation Ax=b, represented in the Compressed Sparse Row (CSR) format
    
    INTEGER,                    POINTER     :: m,n             ! A = [m-by-n]
    INTEGER,                    POINTER     :: nnz_per_row_max ! Maximum number of non-zero entries per row in A (determines how much memory is allocated)
    INTEGER,                    POINTER     :: nnz_max         ! Maximum number of non-zero entries         in A (determines how much memory is allocated)
    INTEGER,  DIMENSION(:    ), POINTER     :: A_ptr
    INTEGER,  DIMENSION(:    ), POINTER     :: A_index
    REAL(dp), DIMENSION(:    ), POINTER     :: A_val
    REAL(dp), DIMENSION(:    ), POINTER     :: x
    REAL(dp), DIMENSION(:    ), POINTER     :: b
    INTEGER :: wm, wn, wnnz_per_row_max, wnnz_max, wA_ptr, wA_index, wA_val, wx, wb
    
  END TYPE type_sparse_matrix_CSR
  
  TYPE type_ice_model
    ! The ice dynamics sub-model data structure.
    
    ! Basic data - ice thickness, bedrock & surface elevation, sea level (geoid elevation), englacial temperature
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_a                  ! Ice thickness [m]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_tplusdt_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb_a                  ! Bedrock elevation [m w.r.t. PD sea level]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_a                  ! Surface elevation [m w.r.t. PD sea level]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SL_a                  ! Sea level (geoid elevation) [m w.r.t. PD sea level]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SL_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SL_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SL_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: TAF_a                 ! Thickness above flotation [m]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: TAF_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: TAF_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: TAF_b
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti_a                  ! Englacial temperature [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti_cx
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti_cy
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti_b
    INTEGER :: wHi_a,  wHi_cx,  wHi_cy,  wHi_b, wHi_tplusdt_a
    INTEGER :: wHb_a,  wHb_cx,  wHb_cy,  wHb_b
    INTEGER :: wHs_a,  wHs_cx,  wHs_cy,  wHs_b
    INTEGER :: wSL_a,  wSL_cx,  wSL_cy,  wSL_b
    INTEGER :: wTAF_a, wTAF_cx, wTAF_cy, wTAF_b
    INTEGER :: wTi_a,  wTi_cx,  wTi_cy,  wTi_b
    
    ! Ice velocities
    REAL(dp), DIMENSION(:,:,:), POINTER     :: u_3D_a                ! 3-D ice velocity [m yr^-1]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: v_3D_a
    REAL(dp), DIMENSION(:,:,:), POINTER     :: u_3D_cx
    REAL(dp), DIMENSION(:,:,:), POINTER     :: v_3D_cy
    REAL(dp), DIMENSION(:,:,:), POINTER     :: w_3D_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_vav_a               ! Vertically averaged ice velocity [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_vav_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_vav_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_vav_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: uabs_vav_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_surf_a              ! Ice velocity at the surface [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_surf_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_surf_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_surf_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: uabs_surf_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_base_a              ! Ice velocity at the base [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_base_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_base_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_base_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: uabs_base_a
    REAL(dp), DIMENSION(:,:,:), POINTER     :: u_3D_SIA_cx           ! Separate fields for the SIA/SSA components, required for the old hybrid method
    REAL(dp), DIMENSION(:,:,:), POINTER     :: v_3D_SIA_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_SSA_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_SSA_cy
    INTEGER :: wu_3D_a, wv_3D_a, wu_3D_cx, wv_3D_cy, ww_3D_a
    INTEGER :: wu_vav_a,  wv_vav_a,  wu_vav_cx,  wv_vav_cy,  wuabs_vav_a
    INTEGER :: wu_surf_a, wv_surf_a, wu_surf_cx, wv_surf_cy, wuabs_surf_a
    INTEGER :: wu_base_a, wv_base_a, wu_base_cx, wv_base_cy, wuabs_base_a
    INTEGER :: wu_3D_SIA_cx, wv_3D_SIA_cy, wu_SSA_cx, wv_SSA_cy
    
    ! Different masks
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_land_a           ! Land touching air or ice
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_ocean_a          ! Land covered by ocean (possibly covered by floating ice)
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_lake_a           ! Land covered by lake  (possibly covered by floating ice)
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_ice_a            ! Any      ice
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_sheet_a          ! Grounded ice
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_shelf_a          ! Floating ice
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_coast_a          ! On A-grid: land bordering ocean
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_coast_cx         ! On C-grid: in between A-land and A-ocean
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_coast_cy
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_margin_a         ! On A-grid: ice  bordering non-ice
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_margin_cx        ! On C-grid: in between A-ice and A-non-ice
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_margin_cy
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_gl_a             ! On A-grid: sheet bordering shelf
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_gl_cx            ! On C-grid: in between A-sheet and A-shelf
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_gl_cy
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_cf_a             ! On A-grid: shelf bordering ocean
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_cf_cx            ! On C-grid: in between A-shelf and A-ocean
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_cf_cy
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_a                ! Multi-info mask, only used for writing to output
    REAL(dp), DIMENSION(:,:  ), POINTER     :: f_grnd_a              ! Grounded fraction (used to determine basal friction in DIVA)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: f_grnd_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: f_grnd_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: f_grnd_b
    INTEGER :: wmask_land_a, wmask_ocean_a, wmask_lake_a, wmask_ice_a, wmask_sheet_a, wmask_shelf_a
    INTEGER :: wmask_coast_a, wmask_coast_cx, wmask_coast_cy, wmask_margin_a, wmask_margin_cx, wmask_margin_cy
    INTEGER :: wmask_gl_a, wmask_gl_cx, wmask_gl_cy, wmask_cf_a, wmask_cf_cx, wmask_cf_cy, wmask_a
    INTEGER :: wf_grnd_a, wf_grnd_cx, wf_grnd_cy, wf_grnd_b
    
    ! Ice physical properties
    REAL(dp), DIMENSION(:,:,:), POINTER     :: A_flow_3D_a           ! Flow parameter [Pa^-3 y^-1]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: A_flow_3D_cx
    REAL(dp), DIMENSION(:,:,:), POINTER     :: A_flow_3D_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: A_flow_vav_a          ! Vertically averaged flow parameter [Pa^-3 y^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: A_flow_vav_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: A_flow_vav_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: A_flow_vav_b
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti_pmp_a              ! The pressure melting point temperature [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Cpi_a                 ! Specific heat capacity of ice [J kg^-1 K^-1].
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ki_a                  ! Conductivity of ice [J m^-1 K^-1 yr^-1].
    INTEGER :: wA_flow_3D_a, wA_flow_3D_cx, wA_flow_3D_cy, wA_flow_vav_a, wA_flow_vav_cx, wA_flow_vav_cy, wA_flow_vav_b
    INTEGER :: wTi_pmp_a, wCpi_a, wKi_a
    
    ! Spatial and temporal derivatives
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHi_dt_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHi_dx_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHi_dx_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHi_dy_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHi_dy_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_dt_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_dx_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_dy_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dt_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dx_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dy_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dx_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dy_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dx_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHs_dy_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dSL_dt_a
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dTi_dx_a
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dTi_dy_a
    INTEGER :: wdHi_dt_a, wdHi_dx_a, wdHi_dy_a, wdHi_dx_cx, wdHi_dy_cy, wdHb_a, wdHb_dt_a, wdHb_dx_a, wdHb_dy_a
    INTEGER :: wdHs_dt_a, wdHs_dx_a, wdHs_dy_a, wdHs_dx_cx, wdHs_dy_cx, wdHs_dx_cy, wdHs_dy_cy, wdSL_dt_a
    INTEGER :: wdTi_dx_a, wdTi_dy_a
    
    ! Zeta derivatives
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dzeta_dt_a
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dzeta_dx_a
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dzeta_dy_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dzeta_dz_a
    INTEGER :: wdzeta_dt_a, wdzeta_dx_a, wdzeta_dy_a, wdzeta_dz_a
    
    ! Ice dynamics - basal conditions and driving stress
    REAL(dp), DIMENSION(:,:  ), POINTER     :: taudx_cx              ! Driving stress taud in the x-direction (on the cx-grid)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: taudy_cy              !      "    "      "     "   y-direction ( "  "  cy-grid)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: phi_fric_a            ! Till friction angle (degrees)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: tauc_a                ! Till yield stress tauc   (used when choice_sliding_law = 'Coloumb' or 'Coulomb_regularised')
    REAL(dp), DIMENSION(:,:  ), POINTER     :: A_slid_a              ! Sliding factor           (used when choice_sliding_law = 'Weertman')
    INTEGER :: wtaudx_cx, wtaudy_cy, wphi_fric_a, wtauc_a, wA_slid_a
    
    ! Ice dynamics - physical terms in the SSA/DIVA
    REAL(dp), DIMENSION(:,:  ), POINTER     :: du_dx_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: du_dy_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dv_dx_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dv_dy_b
    REAL(dp), DIMENSION(:,:,:), POINTER     :: du_dz_3D_cx
    REAL(dp), DIMENSION(:,:,:), POINTER     :: dv_dz_3D_cy
    REAL(dp), DIMENSION(:,:,:), POINTER     :: visc_eff_3D_a
    REAL(dp), DIMENSION(:,:,:), POINTER     :: visc_eff_3D_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: visc_eff_int_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: visc_eff_int_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: N_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: N_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: N_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: N_b
    REAL(dp), DIMENSION(:,:  ), POINTER     :: beta_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: F2_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: beta_eff_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: beta_eff_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: beta_eff_cy
    REAL(dp), DIMENSION(:,:,:), POINTER     :: F1_3D_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: taub_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: taub_cy
    INTEGER :: wdu_dx_b, wdu_dy_b, wdv_dx_b, wdv_dy_b, wdu_dz_3D_cx, wdv_dz_3D_cy
    INTEGER :: wvisc_eff_3D_a, wvisc_eff_3D_b, wvisc_eff_int_a, wvisc_eff_int_b, wN_a, wN_cx, wN_cy, wN_b
    INTEGER :: wbeta_a, wF2_a, wbeta_eff_a, wbeta_eff_cx, wbeta_eff_cy, wF1_3D_a, wtaub_cx, wtaub_cy
    
    ! Ice dynamics - additional solver fields for the SSA/DIVA
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_cx_prev
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_cy_prev
    REAL(dp), DIMENSION(:,:  ), POINTER     :: DIVA_err_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: DIVA_err_cy
    INTEGER,  DIMENSION(:,:  ), POINTER     :: DIVA_mask_combi
    INTEGER,  DIMENSION(:,:  ), POINTER     :: DIVA_isfront_inner
    INTEGER,  DIMENSION(:,:  ), POINTER     :: DIVA_isfront_outer
    INTEGER,  DIMENSION(:,:  ), POINTER     :: DIVA_m_ij2n_u
    INTEGER,  DIMENSION(:,:  ), POINTER     :: DIVA_m_ij2n_v
    INTEGER,  DIMENSION(:,:  ), POINTER     :: DIVA_m_n2ij_uv
    TYPE(type_sparse_matrix_CSR)            :: DIVA_m
    INTEGER :: wu_cx_prev, wv_cy_prev, wDIVA_err_cx, wDIVA_err_cy
    INTEGER :: wDIVA_mask_combi, wDIVA_isfront_inner, wDIVA_isfront_outer
    INTEGER :: wDIVA_m_ij2n_u, wDIVA_m_ij2n_v, wDIVA_m_n2ij_uv
    
    ! Ice dynamics - ice fluxes
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Qx_cx         ! Ice flux in the x-direction on the cx-grid, used for explicit ice thickness update
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Qy_cy         ! Ice flux in the y-direction on the cy-grid, used for explicit ice thickness update
    TYPE(type_sparse_matrix_CSR)            :: dHi_m         ! Sparse matrix for mass conservation,        used for implicit ice thickness update
    INTEGER,  DIMENSION(:,:  ), POINTER     :: dHi_ij2n      ! Matrix-vector translation table,            used for implicit ice thickness update
    INTEGER :: wQx_cx, wQy_cy, wdHi_ij2n
    
    ! Ice dynamics - calving
    REAL(dp), DIMENSION(:,:  ), POINTER     :: float_margin_frac_a   ! Ice-covered fraction for calving front pixels
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_actual_cf_a        ! "actual" ice thickness at calving front pixels (= Hi of thinnest non-calving-front neighbour)
    INTEGER :: wfloat_margin_frac_a, wHi_actual_cf_a
    
    ! Ice dynamics - predictor/corrector ice thickness update
    REAL(dp),                   POINTER     :: pc_zeta
    REAL(dp), DIMENSION(:,:  ), POINTER     :: pc_tau
    REAL(dp), DIMENSION(:,:  ), POINTER     :: pc_fcb
    REAL(dp),                   POINTER     :: pc_eta
    REAL(dp),                   POINTER     :: pc_eta_prev
    REAL(dp),                   POINTER     :: pc_beta1
    REAL(dp),                   POINTER     :: pc_beta2
    REAL(dp),                   POINTER     :: pc_beta3
    REAL(dp),                   POINTER     :: pc_beta4
    REAL(dp), DIMENSION(:,:  ), POINTER     :: pc_f1
    REAL(dp), DIMENSION(:,:  ), POINTER     :: pc_f2
    REAL(dp), DIMENSION(:,:  ), POINTER     :: pc_f3
    REAL(dp), DIMENSION(:,:  ), POINTER     :: pc_f4
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_old
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_pred
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_corr
    INTEGER :: wpc_zeta, wpc_tau, wpc_fcb, wpc_eta, wpc_eta_prev, wpc_beta1, wpc_beta2, wpc_beta3, wpc_beta4
    INTEGER :: wpc_f1, wpc_f2, wpc_f3, wpc_f4, wHi_old, wHi_pred, wHi_corr
    
    ! Thermodynamics
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_ice_a_prev        ! Ice mask from previous time step
    REAL(dp), DIMENSION(:,:  ), POINTER     :: frictional_heating_a   ! Friction heating due to basal sliding
    REAL(dp), DIMENSION(:,:  ), POINTER     :: GHF_a                  ! Geothermal heat flux
    INTEGER :: wmask_ice_a_prev, wfrictional_heating_a, wGHF_a
    
    ! Isotope content
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_a_prev
    REAL(dp), DIMENSION(:,:  ), POINTER     :: IsoRef
    REAL(dp), DIMENSION(:,:  ), POINTER     :: IsoSurf
    REAL(dp), DIMENSION(:,:  ), POINTER     :: IsoIce
    REAL(dp), DIMENSION(:,:  ), POINTER     :: MB_iso
    INTEGER :: wHi_a_prev, wIsoRef, wIsoSurf, wIsoIce, wMB_iso
    
    ! ELRA GIA model
    INTEGER,                    POINTER     :: flex_prof_rad
    REAL(dp), DIMENSION(:,:  ), POINTER     :: flex_prof
    REAL(dp), DIMENSION(:,:  ), POINTER     :: surface_load_topo
    REAL(dp), DIMENSION(:,:  ), POINTER     :: surface_load
    REAL(dp), DIMENSION(:,:  ), POINTER     :: surface_load_rel
    REAL(dp), DIMENSION(:,:  ), POINTER     :: surface_load_rel_ext
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_eq
    INTEGER :: wflex_prof_rad, wflex_prof, wsurface_load_topo, wsurface_load, wsurface_load_rel, wsurface_load_rel_ext, wdHb_eq
    
  END TYPE type_ice_model
  
  TYPE type_zeta_coefficients
    ! Coefficients used in the scaled vertical coordinate transformation
    ! NOTE: local instead of shared memory, since these are not indexed from 1
    !       and that's not possible (or at least not convenient) with shared memory.

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: a_k
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: b_k
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: c_k
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: d_k
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: a_zeta
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: b_zeta
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: c_zeta
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: a_zetazeta
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: b_zetazeta
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: c_zetazeta

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: z_zeta_minus
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: a_zeta_minus
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: b_zeta_minus

    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: b_zeta_plus
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: c_zeta_plus
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: d_zeta_plus
    
  END TYPE type_zeta_coefficients
  
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
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_WE                       ! Monthly mean west-east wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_SN                       ! Monthly mean south_north wind speed (m/s)
    INTEGER :: wT2m, wPrecip, wHs_ref, wWind_WE, wWind_SN ! MPI windows to all these memory spaces
    
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
    
    ! The GCM snapshots.
    TYPE(type_subclimate_global)            :: GCM_PI
    TYPE(type_subclimate_global)            :: GCM_warm
    TYPE(type_subclimate_global)            :: GCM_cold
    
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
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb                            ! Bed topography that was used to force the GCM (m w.r.t. PD sea-level)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs                            ! Orography      that was used to force the GCM (m w.r.t. PD sea-level) (fine   ISM resolution)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_WE                       ! Monthly mean west-east   wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_SN                       ! Monthly mean south-north wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_LR                       ! Monthly mean wind speed in the x-direction (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_DU                       ! Monthly mean wind speed in the y-direction (m/s)
    INTEGER :: wT2m, wPrecip, wHs_ref, wHi, wHb, wHs, wWind_WE, wWind_SN, wWind_LR, wWind_DU ! MPI windows to all these memory spaces
    
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
    TYPE(type_subclimate_region)            :: GCM_PI         ! Pre-industrial  (e.g. HadCM3, Singarayer & Valdes, 2010), for bias correction
    TYPE(type_subclimate_region)            :: GCM_warm       ! Warm            (e.g. HadCM3, Singarayer & Valdes, 2010)
    TYPE(type_subclimate_region)            :: GCM_cold       ! Cold            (e.g. HadCM3, Singarayer & Valdes, 2010)
    TYPE(type_subclimate_region)            :: applied        ! Final applied climate
    
    ! GCM bias
    REAL(dp), DIMENSION(:,:,:), POINTER     :: GCM_bias_T2m
    REAL(dp), DIMENSION(:,:,:), POINTER     :: GCM_bias_Precip
    INTEGER :: wGCM_bias_T2m, wGCM_bias_Precip
          
  END TYPE type_climate_model
  
  TYPE type_SMB_model
    ! The different SMB components, calculated from the prescribed climate
    
    ! Tuning parameters for the IMAU-ITM SMB model (different for each region, set from config)
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
    INTEGER :: wHi, wHb
              
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
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb_raw
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SL_raw
    REAL(dp), DIMENSION(:,:,:), POINTER     :: FirnDepth_raw
    REAL(dp), DIMENSION(:,:  ), POINTER     :: MeltPreviousYear_raw
    INTEGER :: wnx, wny, wnz, wnt, wx, wy, wzeta, wtime, wHi_raw, wHb_raw, wHs_raw, wTi_raw, wdHb_raw, wSL_raw, wFirnDepth_raw, wMeltPreviousYear_raw
    
    ! The data mapped to the model grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Ti
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHb
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SL
    REAL(dp), DIMENSION(:,:,:), POINTER     :: FirnDepth
    REAL(dp), DIMENSION(:,:  ), POINTER     :: MeltPreviousYear
    INTEGER :: wHi, wHb, wTi, wdHb, wSL, wFirnDepth, wMeltPreviousYear
              
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
    INTEGER,                    POINTER     :: ghf_nlon, ghf_nlat
    REAL(dp), DIMENSION(:    ), POINTER     :: ghf_lon, ghf_lat
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ghf_ghf
    INTEGER :: wghf_nlon, wghf_nlat, wghf_lon, wghf_lat, wghf_ghf

    ! External forcing: climate/SMB forcing data
    TYPE(type_netcdf_climate_forcing)       :: netcdf_clim
    INTEGER,                    POINTER     :: clim_nlon, clim_nlat
    INTEGER,                    POINTER     :: clim_nx, clim_ny
    INTEGER,                    POINTER     :: clim_nyears
    REAL(dp), DIMENSION(:    ), POINTER     :: clim_lon, clim_lat
    REAL(dp), DIMENSION(:    ), POINTER     :: clim_x, clim_y
    REAL(dp), DIMENSION(:    ), POINTER     :: clim_time
    REAL(dp),                   POINTER     :: clim_t0, clim_t1
    REAL(dp), DIMENSION(:,:,:), POINTER     :: clim_T2m0, clim_T2m1, clim_T2m2, clim_Precip0,   clim_Precip1,   clim_Precip2     ! Monthly
    REAL(dp), DIMENSION(:,:  ), POINTER     :: clim_SMB0, clim_SMB1, clim_SMB2, clim_T2m_year0, clim_T2m_year1, clim_T2m_year2   ! Yearly
    INTEGER :: wclim_nyears, wclim_nlat, wclim_nlon, wclim_time, wclim_lat, wclim_lon, wclim_t0, wclim_t1, wclim_T2m0, wclim_T2m1, wclim_T2m2
    INTEGER :: wclim_Precip0, wclim_Precip1, wclim_Precip2, wclim_SMB0, wclim_SMB1, wclim_SMB2, wclim_nx, wclim_ny, wclim_x, wclim_y
    INTEGER :: wclim_T2m_year0, wclim_T2m_year1, wclim_T2m_year2
    
  END TYPE type_forcing_data
  
  TYPE type_SELEN_mesh
    ! The global SELEN mesh, in UFEMISM mesh format
    
    INTEGER,                    POINTER     :: nV                            ! Number of vertices
    INTEGER,                    POINTER     :: nTri                          ! Number of triangles
    INTEGER,                    POINTER     :: nC_mem                        ! Maximum allowed number of connections per vertex
    INTEGER :: wnV, wnTri, wnC_mem

    ! Actual mesh data
    REAL(dp), DIMENSION(:,:  ), POINTER     :: V                             ! The X and Y coordinates of all the vertices
    INTEGER,  DIMENSION(:,:  ), POINTER     :: Tri                           ! The triangle array: Tri(ti) = [vi1, vi2, vi3] (vertices ordered counter-clockwise)
    INTEGER,  DIMENSION(:    ), POINTER     :: nC                            ! The number of other vertices this vertex is connected to
    INTEGER,  DIMENSION(:,:  ), POINTER     :: C                             ! The list   of other vertices this vertex is connected to (ordered counter-clockwise, from edge to edge for edge vertices)
    INTEGER,  DIMENSION(:    ), POINTER     :: niTri                         ! The number of triangles this vertex is a part of
    INTEGER,  DIMENSION(:,:  ), POINTER     :: iTri                          ! The list   of triangles this vertex is a part of (ordered counter-clockwise)
    INTEGER :: wV, wTri, wnC, wC, wniTri, wiTri
    
    ! Lat/lon
    REAL(dp), DIMENSION(:    ), POINTER     :: lat                           ! Latitude  (degrees north)
    REAL(dp), DIMENSION(:    ), POINTER     :: lon                           ! Longitude (degrees east)
    INTEGER :: wlat, wlon
    
    ! "anchor points" (used in some places to reduce spherical harmonics calculations)
    INTEGER,                    POINTER     :: nanc                          ! Number of anchor points
    REAL(dp), DIMENSION(:    ), POINTER     :: ancplist_lat                  ! Latitude of each anchor point
    INTEGER,  DIMENSION(:,:  ), POINTER     :: ancplist_isl                  ! Range of global grid pixels for each anchor point
    INTEGER,  DIMENSION(:    ), POINTER     :: ianc                          ! Anchor point index of each global grid pixel
    INTEGER :: wnanc, wancplist_lat, wancplist_isl, wianc
    
  END TYPE type_SELEN_mesh
  
  TYPE type_SELEN_global
    ! Global SELEN data fields, spherical harmonics, etc.
    
    ! Timers to trigger solving the SLE in the main program
    REAL(dp)                                :: t0_SLE, t1_SLE
    
    ! NetCDF global topography input
    TYPE(type_netcdf_SELEN_global_topo)     :: netcdf_topo
    
    ! NetCDF output
    TYPE(type_netcdf_SELEN_output)          :: output
    
    ! The irregular global mesh
    TYPE( type_SELEN_mesh)                  :: mesh
    
    ! Global reference fields
    REAL(dp), DIMENSION(:    ), POINTER     :: topo_ref               ! Reference topography
    REAL(dp), DIMENSION(:    ), POINTER     :: load_ref               ! Reference load
    INTEGER,  DIMENSION(:    ), POINTER     :: of_ref                 ! Reference ocean function
    INTEGER,  DIMENSION(:,:  ), POINTER     :: icemodel_region        ! To which ice model region a SELEN pixel belongs (if any), and their ice model subgrid index
    INTEGER,                    POINTER     :: nel_icemodel           ! Number of pixels that belong to ice model regions
    INTEGER,  DIMENSION(:,:  ), POINTER     :: isl_icemodel           ! List   of pixels that belong to ice model regions
    INTEGER :: wtopo_ref, wload_ref, wof_ref, wicemodel_region, wnel_icemodel, wisl_icemodel
    
    ! Global end results (ice loading, geoid, bed deformation, geoid rate, bed deformation rate)
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi_glob
    REAL(dp), DIMENSION(:    ), POINTER     :: Hi_rel_glob
    REAL(dp), DIMENSION(:    ), POINTER     :: N_glob
    REAL(dp), DIMENSION(:    ), POINTER     :: U_glob
    INTEGER,  DIMENSION(:    ), POINTER     :: of_glob
    INTEGER :: wHi_glob, wHi_rel_glob, wN_glob, wU_glob, wof_glob
    
    ! MJ and LJ values
    INTEGER,  DIMENSION(:    ), POINTER     ::  MJ_VAL,  LJ_VAL
    INTEGER                                 :: wMJ_VAL, wLJ_VAL

    ! Reference data for SELEN in SH
    REAL(dp),   DIMENSION(:,:), POINTER     :: ALF
    COMPLEX*16, DIMENSION(:  ), POINTER     :: dLONG_TABLE
    COMPLEX*16, DIMENSION(:,:), POINTER     :: LONG_TABLE
    INTEGER :: wALF, wdLONG_TABLE, wLONG_TABLE
    
    ! Ice loading history on the global grid on the irregular moving time window
    REAL(dp),   DIMENSION(:    ), POINTER     :: dice_loading_history_irreg_glob
    REAL(dp),   DIMENSION(:,:  ), POINTER     :: ice_loading_history_irreg_glob
    INTEGER :: wdice_loading_history_irreg_glob, wice_loading_history_irreg_glob
    
    ! Memory of SELEN
    COMPLEX*16, DIMENSION(:    ), POINTER     :: MEM_S
    COMPLEX*16, DIMENSION(:    ), POINTER     :: MEM_U
    INTEGER                                   :: wMEM_S, wMEM_U
    
    ! SELEN internal data
    INTEGER,    DIMENSION(:    ), POINTER     :: DM                  ! INTEGER,    DIMENSION( C_SLE%JMAX                )     :: DM
    REAL(dp),   DIMENSION(:    ), POINTER     :: X                   ! REAL(dp),   DIMENSION( C_SLE%NP,  0:C_SLE%TCONV  )     :: X
    INTEGER :: wDM, wX
   
    COMPLEX*16, DIMENSION(:    ), POINTER     :: ZE                  ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: ZE                  ! Eustatic Z array for water loading
    COMPLEX*16, DIMENSION(:    ), POINTER     :: SE                  ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: SE                  ! Eustatic S array for ice loading
    COMPLEX*16, DIMENSION(:    ), POINTER     :: AAAA                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: AAAA
    COMPLEX*16, DIMENSION(:    ), POINTER     :: AAAA_MOD            ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: AAAA_MOD
    COMPLEX*16, DIMENSION(:    ), POINTER     :: BBBB                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: BBBB
    COMPLEX*16, DIMENSION(:    ), POINTER     :: BBBB_MOD            ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: BBBB_MOD
    COMPLEX*16, DIMENSION(:    ), POINTER     :: HHHH                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: HHHH
    COMPLEX*16, DIMENSION(:    ), POINTER     :: KKKK                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: KKKK
    COMPLEX*16, DIMENSION(:    ), POINTER     :: TTTT                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: IIII                ! Ice loading in SH
    COMPLEX*16, DIMENSION(:    ), POINTER     :: IIII                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: TTTT                ! Relative ice loading in SH
    INTEGER :: wZE, wSE, wAAAA, wAAAA_MOD, wBBBB, wBBBB_MOD, wHHHH, wKKKK, wTTTT, wIIII
    
    ! New arrays for the rotational feedback
    COMPLEX*16, DIMENSION(:    ), POINTER     :: load_ice_and_water  ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: load_ice_and_water  ! Ice and water loading in SH for rotational feedback
    COMPLEX*16, DIMENSION(:    ), POINTER     :: ZROT                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: ZROT
    COMPLEX*16, DIMENSION(:    ), POINTER     :: ZROT_MOD            ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: ZROT_MOD
    COMPLEX*16, DIMENSION(:    ), POINTER     :: SROT                ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: SROT
    COMPLEX*16, DIMENSION(:    ), POINTER     :: SROT_MOD            ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: SROT_MOD
    COMPLEX*16, DIMENSION(:    ), POINTER     :: ZROTVV              ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: ZROTVV
    COMPLEX*16, DIMENSION(:    ), POINTER     :: SROTVV              ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: SROTVV
    INTEGER :: wload_ice_and_water, wZROT, wZROT_MOD, wSROT, wSROT_MOD, wZROTVV, wSROTVV
 
    ! Results
    COMPLEX*16, DIMENSION(:    ), POINTER     :: S                   ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: S                   ! relative sea level
    COMPLEX*16, DIMENSION(:    ), POINTER     :: U                   ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: U                   ! Solid Earth deformation (N = S+U) 
    COMPLEX*16, DIMENSION(:    ), POINTER     :: Z                   ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV,0:C%selen_max_iteration)     :: Z
    INTEGER :: wS, wU, wZ
 
    ! New arrays for tdof
    REAL(dp),   DIMENSION(:    ), POINTER     :: newalf              ! REAL(dp),   DIMENSION( C_SLE%JMAX,  C_SLE%NANCH  )     :: newalf              ! New values of alf
    REAL(dp),   DIMENSION(:    ), POINTER     :: slc                 ! REAL(dp),   DIMENSION( C_SLE%NP,  0:C_SLE%TCONV  )     :: slc                 ! sea level forward in time for all elements
    INTEGER,    DIMENSION(:    ), POINTER     :: WET                 ! INTEGER,    DIMENSION( C_SLE%NP,  0:C_SLE%TCONV  )     :: WET                 ! new ocean function (0 or 1) forward in time
    REAL(dp),   DIMENSION(:    ), POINTER     :: newtopo             ! REAL(dp),   DIMENSION( C_SLE%NP,  0:C_SLE%TCONV  )     :: newtopo             ! new topograhpy forward in time
    INTEGER :: wnewalf, wslc, wWET, wnewtopo
 
    INTEGER,    DIMENSION(:    ), POINTER     :: fj                  ! INTEGER,    DIMENSION( C_SLE%NP,  0:C_SLE%TCONV  )     :: fj                  ! floating function for ice thickness (float or not)
    INTEGER,    DIMENSION(:    ), POINTER     :: xcf                 ! INTEGER,    DIMENSION( C_SLE%NP,  0:C_SLE%TCONV+1)     :: xcf
    INTEGER :: wfj, wxcf
 
    INTEGER,    DIMENSION(:    ), POINTER     :: newet2              ! INTEGER,    DIMENSION( C_SLE%NP                  )     :: newet2              ! new ocean function after iteration
    REAL(dp),   DIMENSION(:    ), POINTER     :: newtopo2            ! REAL(dp),   DIMENSION( C_SLE%NP                  )     :: newtopo2            ! new topograhpy after iteration
    INTEGER :: wnewet2, wnewtopo2
 
    REAL(dp),   DIMENSION(:    ), POINTER     :: S_global, U_global  ! REAL(dp),   DIMENSION( C_SLE%NP                  )     :: S_global, U_global  ! Global fields of S and U for output
    INTEGER :: wS_global, wU_global
 
    COMPLEX*16, DIMENSION(:    ), POINTER     :: MSint               ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: MSint               ! Interpolated memory of relative sea level
    COMPLEX*16, DIMENSION(:    ), POINTER     :: MUint               ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: MUint               ! Interpolationd memory of solid earth deformation (N is geoid)
    INTEGER :: wMSint, wMUint
 
    COMPLEX*16, DIMENSION(:    ), POINTER     :: MSread              ! COMPLEX*16, DIMENSION( C_SLE%JMAX                )     :: MSread              ! Initial memory of relative sea level
    INTEGER :: wMSread
 
    ! New variables for spline interpolation
    REAL(dp),   DIMENSION(:    ), POINTER     :: int_time            ! REAL(dp),   DIMENSION(            0:C_SLE%TCONV  )     :: int_time            ! time points (depending on DELTA) for interpolation
    COMPLEX*16, DIMENSION(:    ), POINTER     :: sp1, spn, preSint   ! COMPLEX*16, DIMENSION( C_SLE%JMAX)                     :: sp1, spn, preSint   ! spline in- and output for S
    COMPLEX*16, DIMENSION(:    ), POINTER     :: up1, upn, preUint   ! COMPLEX*16, DIMENSION( C_SLE%JMAX)                     :: up1, upn, preUint   ! spline in- and output for U
    COMPLEX*16, DIMENSION(:    ), POINTER     :: s2,u2               ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: s2,u2               ! spline output
    INTEGER :: wint_time, wsp1, wspn, wpreSint, wup1, wupn, wpreUint, ws2, wu2
 
    COMPLEX*16, DIMENSION(:    ), POINTER     :: varpreoc            ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: varpreoc            ! pre determined SH coefficients ocean function
    COMPLEX*16, DIMENSION(:    ), POINTER     :: varoc               ! COMPLEX*16, DIMENSION( C_SLE%JMAX,0:C_SLE%TCONV  )     :: varoc               ! variable SH coefficient of ocean function
    COMPLEX*16, DIMENSION(:    ), POINTER     :: varoc_inv           ! COMPLEX*16, DIMENSION(            0:C_SLE%TCONV  )     :: varoc_inv           ! 1 / varoc(1,:) - the first degree
    INTEGER :: wvarpreoc, wvaroc, wvaroc_inv
 
  END TYPE type_SELEN_global
  
  TYPE type_SELEN_regional
    ! SELEN data for this model region
    
    ! Mapping between the icemodel square grid and the SELEN global grid
    REAL(dp),                   POINTER :: rad                    ! Radius of the circular disc describing a square grid element (used for mapping)
    INTEGER,  DIMENSION(:,:  ), POINTER :: map_nisl               ! For each ice model square grid cell, the number of SELEN global grid cells to which it contributes
    INTEGER,  DIMENSION(:,:,:), POINTER :: map_isl                ! For each ice model square grid cell, the list   of SELEN global grid cells to which it contributes
    REAL(dp), DIMENSION(:,:,:), POINTER :: map_w                  ! The weights of those contributions
    INTEGER,                    POINTER :: nr                     ! The number of SELEN global grid pixels that lie inside this ice model region (accounting for region overlap)
    INTEGER                             :: ir1, ir2               ! Parallelisation process domains
    INTEGER,  DIMENSION(:    ), POINTER :: map_isl_region2glob    ! For each SELEN global grid pixel inside this ice model region, the corresponding global grid pixel index
    INTEGER :: wrad, wmap_nisl, wmap_isl, wmap_w, wnslr, wmap_isl_region2glob
    
    ! Reference fields
    REAL(dp), DIMENSION(:,:  ), POINTER :: load_ref               ! Reference loading field    (on the regional square grid)
    REAL(dp), DIMENSION(:,:  ), POINTER :: topo_ref               ! Reference topography field (on the regional square grid)
    INTEGER :: wload_ref, wtopo_ref
    
    ! Path to directory containing spherical harmonics binary files for this region
    CHARACTER(LEN=256)                  :: sh_foldername
    
    ! Ice loading history (on the regional square grid, both on the regular and irregular moving time window)
    REAL(dp), DIMENSION(:,:,:), POINTER :: ice_loading_history_reg_sq
    REAL(dp), DIMENSION(:,:,:), POINTER :: ice_loading_history_irreg_sq
    INTEGER :: wice_loading_history_reg_sq, wice_loading_history_irreg_sq
    
    ! SELEN results
    REAL(dp), DIMENSION(:,:  ), POINTER :: dHb_t_grid_GIA
    REAL(dp), DIMENSION(:,:  ), POINTER :: dHb_tplusdt_grid_GIA
    REAL(dp), DIMENSION(:,:  ), POINTER :: SL_t_grid_GIA
    REAL(dp), DIMENSION(:,:  ), POINTER :: SL_tplusdt_grid_GIA
    REAL(dp), DIMENSION(:,:  ), POINTER :: dHb_t
    REAL(dp), DIMENSION(:,:  ), POINTER :: dHb_tplusdt
    REAL(dp), DIMENSION(:,:  ), POINTER :: SL_t
    REAL(dp), DIMENSION(:,:  ), POINTER :: SL_tplusdt
    INTEGER :: wdHb_t_grid_GIA, wdHb_tplusdt_grid_GIA, wSL_t_grid_GIA, wSL_tplusdt_grid_GIA
    INTEGER :: wdHb_t,          wdHb_tplusdt,          wSL_t,          wSL_tplusdt
    
  END TYPE type_SELEN_regional
  
  TYPE type_model_region
    ! Contains all the different data structures, organised by sub-model (ice, climate)
    
    ! Metadata
    CHARACTER(LEN=3)                        :: name                  ! NAM, EAS, GRL, ANT
    CHARACTER(LEN=256)                      :: long_name             ! North America, Eurasia, Greenland, Antarctica
    
    ! The current time (step) of this particular region.
    REAL(dp), POINTER                       :: time
    REAL(dp), POINTER                       :: dt
    REAL(dp), POINTER                       :: dt_prev
    INTEGER :: wtime, wdt, wdt_prev
    
    ! Timers and switches for determining which modules need to be called at what points in time during the simulation
    REAL(dp), POINTER                       :: dt_crit_SIA
    REAL(dp), POINTER                       :: dt_crit_SSA
    REAL(dp), POINTER                       :: dt_crit_ice, dt_crit_ice_prev
    REAL(dp), POINTER                       :: t_last_SIA,     t_next_SIA
    REAL(dp), POINTER                       :: t_last_SSA,     t_next_SSA
    REAL(dp), POINTER                       :: t_last_DIVA,    t_next_DIVA
    REAL(dp), POINTER                       :: t_last_thermo,  t_next_thermo
    REAL(dp), POINTER                       :: t_last_output,  t_next_output
    REAL(dp), POINTER                       :: t_last_climate, t_next_climate
    REAL(dp), POINTER                       :: t_last_SMB,     t_next_SMB
    REAL(dp), POINTER                       :: t_last_BMB,     t_next_BMB
    REAL(dp), POINTER                       :: t_last_ELRA,    t_next_ELRA
    LOGICAL,  POINTER                       :: do_SIA
    LOGICAL,  POINTER                       :: do_SSA
    LOGICAL,  POINTER                       :: do_DIVA
    LOGICAL,  POINTER                       :: do_thermo
    LOGICAL,  POINTER                       :: do_climate
    LOGICAL,  POINTER                       :: do_SMB
    LOGICAL,  POINTER                       :: do_BMB
    LOGICAL,  POINTER                       :: do_output
    LOGICAL,  POINTER                       :: do_ELRA
    INTEGER :: wdt_crit_SIA, wdt_crit_SSA, wdt_crit_ice, wdt_crit_ice_prev
    INTEGER :: wt_last_SIA, wt_last_SSA, wt_last_DIVA, wt_last_thermo, wt_last_output, wt_last_climate, wt_last_SMB, wt_last_BMB, wt_last_ELRA
    INTEGER :: wt_next_SIA, wt_next_SSA, wt_next_DIVA, wt_next_thermo, wt_next_output, wt_next_climate, wt_next_SMB, wt_next_BMB, wt_next_ELRA
    INTEGER ::     wdo_SIA,     wdo_SSA,     wdo_DIVA,     wdo_thermo,     wdo_output,     wdo_climate,     wdo_SMB,     wdo_BMB,     wdo_ELRA
    
    ! The region's ice sheet's volume and volume above flotation (in mSLE, so the second one is the ice sheets GMSL contribution)
    REAL(dp), POINTER                       :: ice_area
    REAL(dp), POINTER                       :: ice_volume
    REAL(dp), POINTER                       :: ice_volume_PD
    REAL(dp), POINTER                       :: ice_volume_above_flotation
    REAL(dp), POINTER                       :: ice_volume_above_flotation_PD
    REAL(dp), POINTER                       :: GMSL_contribution
    INTEGER :: wice_area, wice_volume, wice_volume_PD, wice_volume_above_flotation, wice_volume_above_flotation_PD, wGMSL_contribution
    
    ! Regionally integrated mass balance components
    REAL(dp), POINTER                       :: int_T2m
    REAL(dp), POINTER                       :: int_snowfall
    REAL(dp), POINTER                       :: int_rainfall
    REAL(dp), POINTER                       :: int_melt
    REAL(dp), POINTER                       :: int_refreezing
    REAL(dp), POINTER                       :: int_runoff
    REAL(dp), POINTER                       :: int_SMB
    REAL(dp), POINTER                       :: int_BMB
    REAL(dp), POINTER                       :: int_MB
    INTEGER :: wint_T2m, wint_snowfall, wint_rainfall, wint_melt, wint_refreezing, wint_runoff, wint_SMB, wint_BMB, wint_MB
    
    ! Variables related to the englacial isotope content
    REAL(dp), POINTER                       :: mean_isotope_content
    REAL(dp), POINTER                       :: mean_isotope_content_PD
    REAL(dp), POINTER                       :: d18O_contribution
    REAL(dp), POINTER                       :: d18O_contribution_PD
    INTEGER :: wmean_isotope_content, wmean_isotope_content_PD, wd18O_contribution, wd18O_contribution_PD
        
    ! Reference data fields
    TYPE(type_PD_data_fields)               :: PD                                        ! The present-day  data fields for this model region, on a high-res Cartesian grid
    TYPE(type_PD_data_fields)               :: topo                                      ! The (paleo-)topo data fields for this model region, on a high-res Cartesian grid
    TYPE(type_init_data_fields)             :: init                                      ! The initial      data fields for this model region, on a high-res Cartesian grid
    
    ! Mask where ice is not allowed to form (so Greenland is not included in NAM and EAS, and Ellesmere is not included in GRL)
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_noice
    INTEGER                                 :: wmask_noice
        
    ! Sub-models
    TYPE(type_grid)                         :: grid                                      ! This region's x,y grid
    TYPE(type_grid)                         :: grid_GIA                                  ! Square grid used for GIA (ELRA or SELEN) so that it can use a different resolution from the ice model one
    TYPE(type_ice_model)                    :: ice                                       ! All the ice model data for this model region
    TYPE(type_climate_model)                :: climate                                   ! All the climate data for this model region
    TYPE(type_SMB_model)                    :: SMB                                       ! The different SMB components for this model region
    TYPE(type_BMB_model)                    :: BMB                                       ! The different BMB components for this model region
    TYPE(type_SELEN_regional)               :: SELEN                                     ! SELEN input and output data for this model region
    
    ! Output netcdf files
    TYPE(type_netcdf_restart)               :: restart
    TYPE(type_netcdf_help_fields)           :: help_fields
    TyPE(type_netcdf_scalars_regional)      :: scalars
    
    ! Computation times for this region
    REAL(dp), POINTER                       :: tcomp_total
    REAL(dp), POINTER                       :: tcomp_ice
    REAL(dp), POINTER                       :: tcomp_thermo
    REAL(dp), POINTER                       :: tcomp_climate
    REAL(dp), POINTER                       :: tcomp_GIA
    INTEGER :: wtcomp_total, wtcomp_ice, wtcomp_thermo, wtcomp_climate, wtcomp_GIA
    
  END TYPE type_model_region
  
  TYPE type_global_scalar_data
    ! Structure containing some global scalar values: sea level, CO2, d18O components, computation times, etc.
    
    ! Netcdf file
    TYPE(type_netcdf_scalars_global)        :: netcdf
    
    ! Sea level
    REAL(dp), POINTER                       :: GMSL                                      ! Global mean sea level change
    REAL(dp), POINTER                       :: GMSL_NAM                                  ! Global mean sea level change (contribution from ice in North America)
    REAL(dp), POINTER                       :: GMSL_EAS                                  ! Global mean sea level change (contribution from ice in Eurasia)
    REAL(dp), POINTER                       :: GMSL_GRL                                  ! Global mean sea level change (contribution from ice in Greenland)
    REAL(dp), POINTER                       :: GMSL_ANT                                  ! Global mean sea level change (contribution from ice in Antarctica)
    INTEGER :: wGMSL, wGMSL_NAM, wGMSL_EAS, wGMSL_GRL, wGMSL_ANT
    
    ! CO2
    REAL(dp), POINTER                       :: CO2_obs                                   ! Observed atmospheric CO2
    REAL(dp), POINTER                       :: CO2_mod                                   ! Modelled atmospheric CO2
    INTEGER :: wCO2_obs, wCO2_mod
    
    ! d18O
    REAL(dp), POINTER                       :: d18O_obs                                  ! Observed benthic d18O
    REAL(dp), POINTER                       :: d18O_mod                                  ! Modelled benthic d18O
    REAL(dp), POINTER                       :: d18O_ice                                  ! Contribution to benthic d18O from global ice volume
    REAL(dp), POINTER                       :: d18O_Tdw                                  ! Contribution to benthic d18O from deep-water temperature
    REAL(dp), POINTER                       :: d18O_NAM                                  ! Contribution to benthic d18O from ice in North America
    REAL(dp), POINTER                       :: d18O_EAS                                  ! Contribution to benthic d18O from ice in Eurasia
    REAL(dp), POINTER                       :: d18O_GRL                                  ! Contribution to benthic d18O from ice in Greenland
    REAL(dp), POINTER                       :: d18O_ANT                                  ! Contribution to benthic d18O from ice in Antarctica
    INTEGER :: wd18O_obs, wd18O_mod, wd18O_ice, wd18O_Tdw, wd18O_NAM, wd18O_EAS, wd18O_GRL, wd18O_ANT
    
    ! Temperature
    REAL(dp), POINTER                       :: dT_glob                                   ! Global mean annual surface temperature change
    REAL(dp), POINTER                       :: dT_dw                                     ! Deep-water temperature change
    INTEGER :: wdT_glob, wdT_dw
    
    ! Computation times for all regions combined
    REAL(dp), POINTER                       :: tcomp_total
    REAL(dp), POINTER                       :: tcomp_ice
    REAL(dp), POINTER                       :: tcomp_thermo
    REAL(dp), POINTER                       :: tcomp_climate
    REAL(dp), POINTER                       :: tcomp_GIA
    INTEGER :: wtcomp_total, wtcomp_ice, wtcomp_thermo, wtcomp_climate, wtcomp_GIA
  
  END TYPE type_global_scalar_data
  
CONTAINS

END MODULE data_types_module
