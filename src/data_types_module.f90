MODULE data_types_module

  ! Contains all the different TYPEs for storing data. Put all together in a separate module so that
  ! all subroutines can use all types without interdependency conflicts, and also to make the
  ! modules with the actual physics code more readable.
  ! If only Types could be collapsed in BBEdit...

  USE configuration_module,        ONLY: dp, C
  USE data_types_netcdf_module

  IMPLICIT NONE

  TYPE type_grid
    ! A regular square grid

    INTEGER,                    POINTER     :: nx, ny
    REAL(dp),                   POINTER     :: dx
    REAL(dp), DIMENSION(:    ), POINTER     :: x, y
    REAL(dp),                   POINTER     :: xmin, xmax, ymin, ymax
    INTEGER :: wnx, wny, wdx, wx, wy, wxmin, wxmax, wymin, wymax
    INTEGER                                 :: i1, i2, j1, j2 ! Parallelisation by domain decomposition

    REAL(dp),                   POINTER     :: lambda_M
    REAL(dp),                   POINTER     :: phi_M
    REAL(dp),                   POINTER     :: beta_stereo
    REAL(dp), DIMENSION(:,:  ), POINTER     :: lat, lon
    INTEGER :: wlambda_M, wphi_M, wbeta_stereo, wlat, wlon

  END TYPE type_grid

  TYPE type_grid_lonlat
    ! A regular global lon/lat-grid

    INTEGER,                    POINTER     :: nlon, nlat
    REAL(dp), DIMENSION(:    ), POINTER     :: lon, lat
    INTEGER :: wnlon, wnlat, wlon, wlat
    INTEGER                                 :: i1, i2, j1, j2 ! Parallelisation by domain decomposition

  END TYPE type_grid_lonlat

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

    ! Sub-grid bedrock
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Hb_CDF_a
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Hb_CDF_b
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Hb_CDF_cx
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Hb_CDF_cy
    INTEGER :: wHb_CDF_a, wHb_CDF_b, wHb_CDF_cx, wHb_CDF_cy

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
    REAL(dp), DIMENSION(:,:  ), POINTER     :: w_surf_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_surf_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_surf_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: uabs_surf_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_base_a              ! Ice velocity at the base [m yr^-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_base_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: w_base_a
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_base_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_base_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: uabs_base_a
    REAL(dp), DIMENSION(:,:,:), POINTER     :: u_3D_SIA_cx           ! Separate fields for the SIA/SSA components, required for the old hybrid method
    REAL(dp), DIMENSION(:,:,:), POINTER     :: v_3D_SIA_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: u_SSA_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: v_SSA_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: R_shear               ! Shearing ratio; 1 = full shearing, 0 = full sliding
    INTEGER :: wu_3D_a, wv_3D_a, wu_3D_cx, wv_3D_cy, ww_3D_a
    INTEGER :: wu_vav_a,  wv_vav_a,  wu_vav_cx,  wv_vav_cy,  wuabs_vav_a
    INTEGER :: wu_surf_a, wv_surf_a, wu_surf_cx, wv_surf_cy, wuabs_surf_a, ww_surf_a
    INTEGER :: wu_base_a, wv_base_a, wu_base_cx, wv_base_cy, wuabs_base_a, ww_base_a
    INTEGER :: wu_3D_SIA_cx, wv_3D_SIA_cy, wu_SSA_cx, wv_SSA_cy
    INTEGER :: wR_shear

    ! Different masks
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_land_a           ! Land touching air or ice
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_ocean_a          ! Land covered by ocean (possibly covered by floating ice)
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_lake_a           ! Land covered by lake  (possibly covered by floating ice)
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_ice_a            ! Any      ice
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_sheet_a          ! Grounded ice
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_shelf_a          ! Floating ice
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_coast_a          ! On A-grid: land bordering ocean
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_margin_a         ! On A-grid: ice  bordering non-ice
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_gl_a             ! On A-grid: sheet bordering shelf
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_cf_a             ! On A-grid: shelf bordering ocean
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_a                ! Multi-info mask, only used for writing to output
    REAL(dp), DIMENSION(:,:  ), POINTER     :: f_grnd_a              ! Grounded fraction (used to determine basal friction in DIVA)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: f_grnd_cx
    REAL(dp), DIMENSION(:,:  ), POINTER     :: f_grnd_cy
    REAL(dp), DIMENSION(:,:  ), POINTER     :: f_grnd_b
    INTEGER,  DIMENSION(:,:  ), POINTER     :: basin_ID              ! The drainage basin to which each grid cell belongs
    INTEGER,                    POINTER     :: nbasins               ! Total number of basins defined for this region
    INTEGER :: wmask_land_a, wmask_ocean_a, wmask_lake_a, wmask_ice_a, wmask_sheet_a, wmask_shelf_a
    INTEGER :: wmask_coast_a, wmask_margin_a, wmask_gl_a, wmask_cf_a, wmask_a
    INTEGER :: wf_grnd_a, wf_grnd_cx, wf_grnd_cy, wf_grnd_b
    INTEGER :: wbasin_ID, wnbasins

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

    ! Ice dynamics - driving stress
    REAL(dp), DIMENSION(:,:  ), POINTER     :: taudx_cx              ! Driving stress taud in the x-direction (on the cx-grid)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: taudy_cy              !      "    "      "     "   y-direction ( "  "  cy-grid)
    INTEGER :: wtaudx_cx, wtaudy_cy

    ! Ice dynamics - basal hydrology
    REAL(dp), DIMENSION(:,:  ), POINTER     :: overburden_pressure_a ! Overburden pressure ( = H * rho_i * g) [Pa]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: pore_water_pressure_a ! Pore water pressure (determined by basal hydrology model) [Pa]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Neff_a                ! Effective pressure ( = overburden pressure - pore water pressure) [Pa]
    INTEGER :: woverburden_pressure_a, wpore_water_pressure_a, wNeff_a

    ! Ice dynamics - basal roughness / friction
    REAL(dp), DIMENSION(:,:  ), POINTER     :: phi_fric_a            ! Till friction angle (degrees)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: tauc_a                ! Till yield stress tauc   (used when choice_sliding_law = "Coloumb" or "Coulomb_regularised")
    REAL(dp), DIMENSION(:,:  ), POINTER     :: alpha_sq_a            ! Coulomb-law friction coefficient [unitless]         (used when choice_sliding_law =             "Tsai2015", or "Schoof2005")
    REAL(dp), DIMENSION(:,:  ), POINTER     :: beta_sq_a             ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3] (used when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005")
    INTEGER :: wphi_fric_a, wtauc_a, walpha_sq_a, wbeta_sq_a

    ! Ice dynamics - basal inversion
    REAL(dp), DIMENSION(:,:  ), POINTER     :: BIV_uabs_surf_target  ! Target surface velocity for the basal inversion [m/yr]
    INTEGER :: wBIV_uabs_surf_target

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
    INTEGER                                 :: DIVA_SOR_nit
    REAL(dp)                                :: DIVA_SOR_tol
    REAL(dp)                                :: DIVA_SOR_omega
    REAL(dp)                                :: DIVA_PETSc_rtol
    REAL(dp)                                :: DIVA_PETSc_abstol
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
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_eff_cf_a           ! Effective ice thickness at calving front pixels (= Hi of thinnest non-calving-front neighbour)
    INTEGER :: wfloat_margin_frac_a, wHi_eff_cf_a

    ! Ice dynamics - prescribed retreat mask
    REAL(dp),                   POINTER     :: ice_fraction_retreat_mask_t0  ! Time of first     prescribed retreat mask timeframe ( <= region%time)
    REAL(dp),                   POINTER     :: ice_fraction_retreat_mask_t1  ! Time of second    prescribed retreat mask timeframe ( >= region%time)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ice_fraction_retreat_mask0    !         First     prescribed retreat mask timeframe
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ice_fraction_retreat_mask1    !         Second    prescribed retreat mask timeframe
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ice_fraction_retreat_mask     ! Time-interpolated prescribed retreat mask
    REAL(dp), DIMENSION(:,:  ), POINTER     :: retreat_mask_Hi_ref
    INTEGER :: wice_fraction_retreat_mask_t0, wice_fraction_retreat_mask_t1, wice_fraction_retreat_mask0
    INTEGER :: wice_fraction_retreat_mask1, wice_fraction_retreat_mask, wretreat_mask_Hi_ref

    ! Ice dynamics - predictor/corrector ice thickness update
    REAL(dp),                   POINTER     :: pc_zeta
    REAL(dp), DIMENSION(:,:  ), POINTER     :: pc_tau
    REAL(dp),                   POINTER     :: pc_eta
    REAL(dp),                   POINTER     :: pc_eta_prev
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHidt_Hnm1_unm1
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHidt_Hn_un
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dHidt_Hstarnp1_unp1
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_old
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_pred
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi_corr
    INTEGER :: wpc_zeta, wpc_tau, wpc_eta, wpc_eta_prev
    INTEGER :: wdHidt_Hnm1_unm1, wdHidt_Hn_un, wdHidt_Hstarnp1_unp1, wHi_old, wHi_pred, wHi_corr

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

  TYPE type_climate_snapshot
    ! A single climate snapshot

    CHARACTER(LEN=256)                      :: name                          ! 'ERA40', 'HadCM3_PI', etc.

    ! Metadata
    REAL(dp),                   POINTER     :: CO2
    REAL(dp),                   POINTER     :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp),                   POINTER     :: orbit_ecc                     ! Orbital parameters
    REAL(dp),                   POINTER     :: orbit_obl
    REAL(dp),                   POINTER     :: orbit_pre
    REAL(dp),                   POINTER     :: sealevel
    INTEGER :: wCO2, worbit_time, worbit_ecc, worbit_obl, worbit_pre, wsealevel

    ! Climate data
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs                            ! Orography (m w.r.t. PD sea level)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T2m                           ! Monthly mean 2m air temperature (K)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Precip                        ! Monthly mean precipitation (m)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_WE                       ! Monthly mean west-east   wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_SN                       ! Monthly mean south-north wind speed (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_LR                       ! Monthly mean wind speed in the x-direction (m/s)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_DU                       ! Monthly mean wind speed in the y-direction (m/s)
    INTEGER :: wHs, wT2m, wPrecip, wHs_ref, wWind_WE, wWind_SN, wWind_LR, wWind_DU

    ! Spatially variable lapse rate for GCM snapshots (see Berends et al., 2018)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: lambda
    INTEGER :: wlambda

    ! Reference absorbed insolation (for GCM snapshots), or insolation at model time for the applied climate
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Q_TOA                         ! Monthly mean insolation at the top of the atmosphere (W/m2) (taken from the prescribed insolation solution at orbit_time)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Albedo                        ! Monthly mean surface albedo (calculated using our own SMB scheme for consistency)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: I_abs                         ! Total yearly absorbed insolation, used in the climate matrix for interpolation
    INTEGER :: wQ_TOA, wAlbedo, wI_abs

  END TYPE type_climate_snapshot

  TYPE type_climate_model_PD_obs
    ! The "PD_obs" climate model option: fixed present-day climate

    ! The climate snapshot
    TYPE(type_climate_snapshot)             :: snapshot

  END TYPE type_climate_model_PD_obs

  TYPE type_climate_model_direct
    ! The "direct_climate" climate model option: directly prescribed climate

    ! Timestamps of the two timeframes
    REAL(dp)                                :: t0, t1

    ! The two timeframes enveloping the model time
    TYPE( type_climate_snapshot)            :: timeframe0
    TYPE( type_climate_snapshot)            :: timeframe1

  END TYPE type_climate_model_direct

  TYPE type_climate_model_direct_SMB_snapshot
    ! A single timeframe of monthly SMB and monthly 2-m air temperature on the model grid

    REAL(dp), DIMENSION(:,:  ), POINTER     :: SMB
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T2m
    INTEGER :: wSMB, wT2m

  END TYPE type_climate_model_direct_SMB_snapshot

  TYPE type_climate_model_direct_SMB
    ! The "direct_SMB" climate model option: directly prescribed SMB + temperature

    ! Timestamps of the two timeframes
    REAL(dp)                                :: t0, t1

    ! The two timeframes enveloping the model time
    TYPE( type_climate_model_direct_SMB_snapshot) :: timeframe0
    TYPE( type_climate_model_direct_SMB_snapshot) :: timeframe1

  END TYPE type_climate_model_direct_SMB

  TYPE type_climate_model_ISMIP_style_snapshot
    ! A single timeframe of data fields for the ISMIP-style (SMB + aSMB + dSMBdz + ST + aST + dSTdz) forcing

    ! ISMIP-style forcing data
    REAL(dp), DIMENSION(:,:  ), POINTER     :: aSMB
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dSMBdz
    REAL(dp), DIMENSION(:,:  ), POINTER     :: aST
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dSTdz
    INTEGER :: waSMB, wdSMBdz, waST, wdSTdz

  END TYPE type_climate_model_ISMIP_style_snapshot

  TYPE type_climate_model_ISMIP_style
    ! The "ISMIP-style" climate model option: aSMB + dSMBdz + aST + dSTdz

    ! The baseline elevastion, SMB, and surface temperature
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs_ref
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SMB_ref
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ST_ref
    INTEGER :: wHs_ref, wSMB_ref, wST_ref

    ! Timestamps of the two timeframes
    REAL(dp)                                :: t0, t1

    ! The two timeframes enveloping the model time
    TYPE( type_climate_model_ISMIP_style_snapshot) :: timeframe0
    TYPE( type_climate_model_ISMIP_style_snapshot) :: timeframe1

    ! Time-interpolated values
    REAL(dp), DIMENSION(:,:  ), POINTER     :: aSMB
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dSMBdz
    REAL(dp), DIMENSION(:,:  ), POINTER     :: aST
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dSTdz
    INTEGER :: waSMB, wdSMBdz, waST, wdSTdz

   ! Elevation-corrected temperature and SMB
    REAL(dp), DIMENSION(:,:  ), POINTER     :: SMB
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ST
    INTEGER :: wSMB, wST

  END TYPE type_climate_model_ISMIP_style

  TYPE type_climate_model_matrix
    ! The "matrix" climate model option: three GCM snapshots (warm, cold, and PI), and a PD reanalysis snapshot to use for bias correction

    ! The three GCM snapshots
    TYPE(type_climate_snapshot)             :: GCM_PI
    TYPE(type_climate_snapshot)             :: GCM_warm
    TYPE(type_climate_snapshot)             :: GCM_cold

    ! The present-day climate
    TYPE(type_climate_snapshot)             :: PD_obs

    ! GCM bias
    REAL(dp), DIMENSION(:,:,:), POINTER     :: GCM_bias_T2m
    REAL(dp), DIMENSION(:,:,:), POINTER     :: GCM_bias_Precip
    INTEGER :: wGCM_bias_T2m, wGCM_bias_Precip

    ! Total yearly absorbed insolation, used in the climate matrix for interpolation
    REAL(dp), DIMENSION(:,:  ), POINTER     :: I_abs
    INTEGER :: wI_abs

  END TYPE type_climate_model_matrix

  TYPE type_climate_model
    ! All the different climate model options, and the applied climate

    ! All the different climate model options
    TYPE(type_climate_model_PD_obs)         :: PD_obs             ! The "PD_obs"          climate model option: fixed present-day climate
    TYPE(type_climate_model_direct)         :: direct             ! The "direct_climate"  climate model option: directly prescribed climate
    TYPE(type_climate_model_direct_SMB)     :: direct_SMB         ! The "direct_SMB"      climate model option: directly prescribed SMB + temperature
    TYPE(type_climate_model_ISMIP_style)    :: ISMIP_style        ! The "ISMIP-style"     climate model option: aSMB + dSMBdz + aST + dSTdz
    TYPE(type_climate_model_matrix)         :: matrix             ! The "matrix"          climate model option: three GCM snapshots (warm, cold, and PI), and a PD reanalysis snapshot to use for bias correction

    ! The applied climate
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Q_TOA              ! Monthly mean insolation at the top of the atmosphere [W/m2]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T2m                ! Monthly mean 2-m air temperature                     [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Precip             ! Monthly total precipitation                          [m]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_LR            ! Monthly mean wind in the x-direction                 [m/s]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Wind_DU            ! Monthly mean wind in the y-direction                 [m/s]
    INTEGER :: wQ_TOA, wT2m, wPrecip, wWind_LR, wWind_DU

  END TYPE type_climate_model

  TYPE type_ocean_snapshot_global
    ! Global ocean snapshot, either from present-day observations (e.g. WOA18) or from a GCM ocean snapshot.

    CHARACTER(LEN=256)                      :: name                          ! 'WOA', 'COSMOS_LGM', etc.

    ! NetCDF file containing the data
    TYPE(type_netcdf_global_ocean_data)     :: netcdf

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
    INTEGER :: wCO2, worbit_time, worbit_ecc, worbit_obl, worbit_pre

    ! Global 3-D ocean data on the original vertical grid
    REAL(dp),                   POINTER     :: T_ocean_mean                  ! Regional mean ocean temperature (used for basal melt when no ocean temperature data is provided)
    REAL(dp), DIMENSION(:    ), POINTER     :: z_ocean_raw                   ! Vertical coordinate of the 3-D ocean fields [m below sea surface]
    INTEGER,                    POINTER     :: nz_ocean_raw                  ! Number of vertical layers in the 3-D ocean fields
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T_ocean_raw                   ! 3-D annual mean ocean temperature [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: S_ocean_raw                   ! 3-D annual mean ocean salinity    [PSU]
    INTEGER :: wT_ocean_mean, wz_ocean_raw, wnz_ocean_raw, wT_ocean_raw, wS_ocean_raw

    ! Global 3-D ocean data on the ice-model vertical grid
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T_ocean                       ! 3-D annual mean ocean temperature [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: S_ocean                       ! 3-D annual mean ocean salinity    [PSU]
    INTEGER :: wT_ocean, wS_ocean

    ! Paralelisation
    INTEGER                                 :: i1, i2                        ! Grid domain (:,i1:i2) of each process

  END TYPE type_ocean_snapshot_global

  TYPE type_ocean_matrix_global
    ! The ocean matrix data structure. Contains the PD observations and all the different global GCM snapshots.

    ! The present-day observed ocean (e.g. WOA18)
    TYPE(type_ocean_snapshot_global)        :: PD_obs

    ! The GCM ocean snapshots for climate and ocean.
    TYPE(type_ocean_snapshot_global)        :: GCM_PI
    TYPE(type_ocean_snapshot_global)        :: GCM_warm
    TYPE(type_ocean_snapshot_global)        :: GCM_cold

  END TYPE type_ocean_matrix_global

  TYPE type_ocean_snapshot_regional
    ! Global ocean snapshot, either from present-day observations (e.g. WOA18) or from a GCM snapshot, projected onto the regional model grid.

    CHARACTER(LEN=256)                      :: name                          ! 'ERA40', 'HadCM3_PI', etc.

    ! General forcing info (not relevant for PD observations)
    REAL(dp),                   POINTER     :: CO2                           ! CO2 concentration in ppm that was used to force the GCM
    REAL(dp),                   POINTER     :: orbit_time                    ! The time (in ky ago) for the orbital forcing (Q_TOA can then be read from Laskar data)
    REAL(dp),                   POINTER     :: orbit_ecc                     ! Orbital parameters that were used to force the GCM
    REAL(dp),                   POINTER     :: orbit_obl
    REAL(dp),                   POINTER     :: orbit_pre
    REAL(dp),                   POINTER     :: sealevel
    INTEGER :: wCO2, worbit_time, worbit_ecc, worbit_obl, worbit_pre, wsealevel

    ! Ocean data
    REAL(dp),                   POINTER     :: T_ocean_mean                  ! Regional mean ocean temperature (used for basal melt when no ocean temperature data is provided)
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T_ocean                       ! 3-D annual mean ocean temperature [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: S_ocean                       ! 3-D annual mean ocean salinity    [PSU]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T_ocean_ext                   ! 3-D annual mean ocean temperature, extrapolated beneath ice shelves [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: S_ocean_ext                   ! 3-D annual mean ocean salinity   , extrapolated beneath ice shelves [PSU]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T_ocean_corr_ext              ! Bias-corrected 3-D annual mean ocean temperature, extrapolated beneath ice shelves [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: S_ocean_corr_ext              ! Bias-corrected 3-D annual mean ocean salinity,    extrapolated beneath ice shelves [PSU]
    INTEGER :: wT_ocean_mean, wT_ocean, wS_ocean, wT_ocean_ext, wS_ocean_ext, wT_ocean_corr_ext, wS_ocean_corr_ext

    ! History of the weighing fields
    REAL(dp), DIMENSION(:,:,:), POINTER     :: w_tot_history
    INTEGER,                    POINTER     :: nw_tot_history
    INTEGER :: ww_tot_history, wnw_tot_history

  END TYPE type_ocean_snapshot_regional

  TYPE type_ocean_matrix_regional
    ! All the relevant ocean data fields (PD observations, GCM snapshots, and final, applied ocean) on the model region grid

    TYPE(type_ocean_snapshot_regional)      :: PD_obs                        ! PD observations (e.g. WOA18)
    TYPE(type_ocean_snapshot_regional)      :: GCM_PI                        ! Pre-industrial reference snapshot, for bias correction
    TYPE(type_ocean_snapshot_regional)      :: GCM_warm                      ! Warm snapshot
    TYPE(type_ocean_snapshot_regional)      :: GCM_cold                      ! Cold snapshot
    TYPE(type_ocean_snapshot_regional)      :: applied                       ! Final applied ocean

  END TYPE type_ocean_matrix_regional

  TYPE type_highres_ocean_data
    ! High-resolution (extrapolated) regional ocean data

    ! NetCDF files
    TYPE( type_netcdf_reference_geometry)      :: netcdf_geo
    TYPE( type_netcdf_extrapolated_ocean_data) :: netcdf

    ! Grid
    TYPE( type_grid)                        :: grid

    ! Geometry data
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi                            ! Ice thickness [m]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb                            ! Bedrock elevation [m]
    INTEGER :: wHi, wHb

    ! Ice basins
    INTEGER,  DIMENSION(:,:  ), POINTER     :: basin_ID                      ! The drainage basin to which each grid cell belongs
    INTEGER,                    POINTER     :: nbasins                       ! Total number of basins defined for this region
    INTEGER :: wbasin_ID, wnbasins

    ! Raw and extrapolated ocean data
    REAL(dp), DIMENSION(:,:,:), POINTER     :: T_ocean                       ! 3-D annual mean ocean temperature [K]
    REAL(dp), DIMENSION(:,:,:), POINTER     :: S_ocean                       ! 3-D annual mean ocean salinity    [PSU]
    INTEGER ::wT_ocean, wS_ocean

  END TYPE type_highres_ocean_data

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

    ! General data fields
    REAL(dp), DIMENSION(:,:  ), POINTER     :: BMB                           ! The basal mass balance (same as SMB: negative means ice loss, positive means ice gain!) [m/yr]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: BMB_sheet                     ! The basal mass balance underneath the land-based ice sheet [m/yr]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: BMB_shelf                     ! The basal mass balance underneath the floating   ice shelf [m/yr]

    ! The linear/quadratic models from Favier et al. (2019)
    ! =====================================================

    REAL(dp), DIMENSION(:,:  ), POINTER     :: T_ocean_base                  ! Ocean temperature    at the ice shelf base
    REAL(dp), DIMENSION(:,:  ), POINTER     :: T_ocean_freeze_base           ! Ocean freezing point at the ice shelf base (depends on pressure and salinity)
    INTEGER :: wT_ocean_base, wT_ocean_freeze_base

    ! The Lazeroms (2018) plume model
    ! ===============================

    ! NOTE: also uses T_ocean_base!

    INTEGER,  DIMENSION(:,:  ), POINTER     :: search_directions             ! The 16 search directions
    REAL(dp), DIMENSION(:,:  ), POINTER     :: eff_plume_source_depth        ! Effective plume source depth (average source depth over all valid plume paths)
    REAL(dp), DIMENSION(:,:  ), POINTER     :: eff_basal_slope               ! Effective basal slope        (average slope        over all valid plume paths)
    INTEGER :: wsearch_directions, weff_plume_source_depth, weff_basal_slope

    ! The PICO model
    ! ==============

    ! NOTE: also uses search_directions!

    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_d_GL                     ! Distance to grounding line [m]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_d_IF                     ! Distance to ice front      [m]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_r                        ! Relative distance to grounding line [0-1]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_A                        ! Area covered by each ocean box in each basin [n_basins x n_boxes]
    INTEGER,  DIMENSION(:    ), POINTER     :: PICO_n_D                      ! Number of ocean boxes for each ice basin
    INTEGER,  DIMENSION(:,:  ), POINTER     :: PICO_k                        ! PICO ocean box number to which the shelf grid cells belong
    INTEGER :: wPICO_d_GL, wPICO_d_IF, wPICO_r, wPICO_A, wPICO_n_D, wPICO_k

    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_T                        ! 2-D     ambient temperature [K]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_Tk                       ! Average ambient temperature within each basin-box
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_S                        ! 2-D     ambient salinity    [PSU]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_Sk                       ! Average ambient salinity    within each basin-box
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_p                        ! 2-D     basal pressure      [Pa]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_pk                       ! Average basal pressure      within each basin-box
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_m                        ! 2-D     melt rate           [m/yr]
    REAL(dp), DIMENSION(:,:  ), POINTER     :: PICO_mk                       ! Average melt rate           within each basin-box
    INTEGER :: wPICO_T, wPICO_Tk, wPICO_S, wPICO_Sk, wPICO_p, wPICO_pk, wPICO_m, wPICO_mk

    ! The ANICE_legacy BMB model
    ! ==========================

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

    ! Additional data fields
    REAL(dp), DIMENSION(:,:  ), POINTER     :: sub_angle                     ! "subtended angle"      for the sub-shelf melt parameterisation
    REAL(dp), DIMENSION(:,:  ), POINTER     :: dist_open                     ! distance to open ocean for the sub-shelf melt parameterisation
    INTEGER :: wBMB, wBMB_sheet, wBMB_shelf, wsub_angle, wdist_open

  END TYPE type_BMB_model

  TYPE type_reference_geometry
    ! Data structure containing a reference ice-sheet geometry

    ! NetCDF file containing the data
    TYPE(type_netcdf_reference_geometry)    :: netcdf

    ! Data on the model grid
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hi
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hb
    REAL(dp), DIMENSION(:,:  ), POINTER     :: Hs
    INTEGER :: wHi, wHb, wHs

    ! Sub-grid bedrock
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Hb_CDF_a
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Hb_CDF_b
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Hb_CDF_cx
    REAL(dp), DIMENSION(:,:,:), POINTER     :: Hb_CDF_cy
    INTEGER :: wHb_CDF_a, wHb_CDF_b, wHb_CDF_cx, wHb_CDF_cy

  END TYPE type_reference_geometry

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

    ! CO2 record
    REAL(dp), DIMENSION(:    ), POINTER     :: CO2_time
    REAL(dp), DIMENSION(:    ), POINTER     :: CO2_record
    REAL(dp),                   POINTER     :: CO2_obs
    INTEGER :: wCO2_time, wCO2_record, wCO2_obs

    ! d18O record
    REAL(dp), DIMENSION(:    ), POINTER     :: d18O_time
    REAL(dp), DIMENSION(:    ), POINTER     :: d18O_record
    REAL(dp),                   POINTER     :: d18O_obs
    REAL(dp),                   POINTER     :: d18O_obs_PD
    INTEGER :: wd18O_time, wd18O_record, wd18O_obs, wd18O_obs_PD

    ! Insolation reconstruction
    TYPE(type_netcdf_insolation)            :: netcdf_ins
    INTEGER,                    POINTER     :: ins_nyears
    INTEGER,                    POINTER     :: ins_nlat
    REAL(dp), DIMENSION(:    ), POINTER     :: ins_time
    REAL(dp), DIMENSION(:    ), POINTER     :: ins_lat
    REAL(dp),                   POINTER     :: ins_t0, ins_t1
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ins_Q_TOA0, ins_Q_TOA1
    INTEGER :: wins_nyears, wins_nlat, wins_time, wins_lat, wins_t0, wins_t1, wins_Q_TOA0, wins_Q_TOA1

    ! Geothermal heat flux
    TYPE(type_netcdf_geothermal_heat_flux)  :: netcdf_ghf
    INTEGER,                    POINTER     :: ghf_nlon, ghf_nlat
    REAL(dp), DIMENSION(:    ), POINTER     :: ghf_lon, ghf_lat
    REAL(dp), DIMENSION(:,:  ), POINTER     :: ghf_ghf
    INTEGER :: wghf_nlon, wghf_nlat, wghf_lon, wghf_lat, wghf_ghf

    ! External forcing: sea level record
    REAL(dp), DIMENSION(:    ), POINTER     :: sealevel_time
    REAL(dp), DIMENSION(:    ), POINTER     :: sealevel_record
    REAL(dp),                   POINTER     :: sealevel_obs
    INTEGER :: wsealevel_time, wsealevel_record, wsealevel_obs

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
    REAL(dp), POINTER                       :: t_last_ocean,   t_next_ocean
    REAL(dp), POINTER                       :: t_last_SMB,     t_next_SMB
    REAL(dp), POINTER                       :: t_last_BMB,     t_next_BMB
    REAL(dp), POINTER                       :: t_last_ELRA,    t_next_ELRA
    REAL(dp), POINTER                       :: t_last_BIV,     t_next_BIV
    LOGICAL,  POINTER                       :: do_SIA
    LOGICAL,  POINTER                       :: do_SSA
    LOGICAL,  POINTER                       :: do_DIVA
    LOGICAL,  POINTER                       :: do_thermo
    LOGICAL,  POINTER                       :: do_climate
    LOGICAL,  POINTER                       :: do_ocean
    LOGICAL,  POINTER                       :: do_SMB
    LOGICAL,  POINTER                       :: do_BMB
    LOGICAL,  POINTER                       :: do_output
    LOGICAL,  POINTER                       :: do_ELRA
    LOGICAL,  POINTER                       :: do_BIV
    INTEGER :: wdt_crit_SIA, wdt_crit_SSA, wdt_crit_ice, wdt_crit_ice_prev
    INTEGER :: wt_last_SIA, wt_last_SSA, wt_last_DIVA, wt_last_thermo, wt_last_output, wt_last_climate, wt_last_ocean, wt_last_SMB, wt_last_BMB, wt_last_ELRA, wt_last_BIV
    INTEGER :: wt_next_SIA, wt_next_SSA, wt_next_DIVA, wt_next_thermo, wt_next_output, wt_next_climate, wt_next_ocean, wt_next_SMB, wt_next_BMB, wt_next_ELRA, wt_next_BIV
    INTEGER ::     wdo_SIA,     wdo_SSA,     wdo_DIVA,     wdo_thermo,     wdo_output,     wdo_climate,     wdo_ocean,     wdo_SMB,     wdo_BMB,     wdo_ELRA,     wdo_BIV

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

    ! Reference geometries
    TYPE(type_reference_geometry)           :: refgeo_init                               ! Initial         ice-sheet geometry
    TYPE(type_reference_geometry)           :: refgeo_PD                                 ! Present-day     ice-sheet geometry
    TYPE(type_reference_geometry)           :: refgeo_GIAeq                              ! GIA equilibrium ice-sheet geometry

    ! Mask where ice is not allowed to form (so Greenland is not included in NAM and EAS, and Ellesmere is not included in GRL)
    INTEGER,  DIMENSION(:,:  ), POINTER     :: mask_noice
    INTEGER                                 :: wmask_noice

    ! Sub-models
    TYPE(type_grid)                         :: grid                                      ! This region's x,y grid
    TYPE(type_grid)                         :: grid_GIA                                  ! Square grid used for GIA (ELRA or SELEN) so that it can use a different resolution from the ice model one
    TYPE(type_ice_model)                    :: ice                                       ! All the ice model data for this model region
    TYPE(type_climate_model)                :: climate                                   ! All the climate data for this model region
    TYPE(type_ocean_matrix_regional)        :: ocean_matrix                              ! All the ocean data for this model region
    TYPE(type_SMB_model)                    :: SMB                                       ! The different SMB components for this model region
    TYPE(type_BMB_model)                    :: BMB                                       ! The different BMB components for this model region
    TYPE(type_SELEN_regional)               :: SELEN                                     ! SELEN input and output data for this model region

    ! Output netcdf files
    TYPE(type_netcdf_restart)               :: restart
    TYPE(type_netcdf_help_fields)           :: help_fields
    TYPE(type_netcdf_scalars_regional)      :: scalars

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
