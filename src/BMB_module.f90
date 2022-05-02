MODULE BMB_module

  ! Contains all the routines for calculating the basal mass balance.

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_ocean_snapshot_regional, type_BMB_model, type_reference_geometry
  USE netcdf_module,                   ONLY: debug, write_to_debug_file 
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             interpolate_ocean_depth, interp_bilin_2D
  USE forcing_module,                  ONLY: forcing, get_insolation_at_time_month_and_lat
  IMPLICIT NONE
    
CONTAINS

! == The main routines that should be called from the main ice model/program
! ==========================================================================

  SUBROUTINE run_BMB_model( grid, ice, ocean, BMB, region_name, time, refgeo_init)
    ! Run the selected BMB model
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                      INTENT(IN)    :: grid 
    TYPE(type_ice_model),                 INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),   INTENT(INOUT) :: ocean
    TYPE(type_BMB_model),                 INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                     INTENT(IN)    :: region_name
    REAL(dp),                             INTENT(IN)    :: time
    TYPE(type_reference_geometry),        INTENT(IN)    :: refgeo_init
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: amplification_factor
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Run the selected shelf BMB model
    IF     (C%choice_BMB_shelf_model == 'uniform') THEN
      BMB%BMB_shelf( :,grid%i1:grid%i2) = C%BMB_shelf_uniform
      CALL sync
    ELSEIF (C%choice_BMB_shelf_model == 'idealised') THEN
      CALL run_BMB_model_idealised(            grid, ice,        BMB, time)
    ELSEIF (C%choice_BMB_shelf_model == 'ANICE_legacy') THEN
      CALL run_BMB_model_ANICE_legacy(         grid, ice,        BMB, region_name, time)
    ELSEIF (C%choice_BMB_shelf_model == 'Favier2019_lin') THEN
      CALL run_BMB_model_Favier2019_linear(    grid, ice, ocean, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'Favier2019_quad') THEN
      CALL run_BMB_model_Favier2019_quadratic( grid, ice, ocean, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'Favier2019_Mplus') THEN
      CALL run_BMB_model_Favier2019_Mplus(     grid, ice, ocean, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'Lazeroms2018_plume') THEN
      CALL run_BMB_model_Lazeroms2018_plume(   grid, ice, ocean, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'PICO') THEN
      CALL run_BMB_model_PICO(                 grid, ice, ocean, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'PICOP') THEN
      CALL run_BMB_model_PICOP(                grid, ice, ocean, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'inverse_shelf_geometry') THEN
      CALL inverse_BMB_shelf_geometry(         grid, ice,        BMB, refgeo_init)
    ELSE
      CALL crash('unknown choice_BMB_shelf_model "' // TRIM(C%choice_BMB_shelf_model) // '"!')
    END IF
    
    ! Run the selected sheet BMB model
    IF     (C%choice_BMB_sheet_model == 'uniform') THEN
      BMB%BMB_sheet( :,grid%i1:grid%i2) = C%BMB_sheet_uniform
    ELSE
      CALL crash('unknown choice_BMB_sheet_model "' // TRIM(C%choice_BMB_sheet_model) // '"!')
    END IF
    
    ! Extrapolate melt field from the regular (FCMP) mask to the (extended) PMP mask
    IF (C%choice_BMB_subgrid == 'PMP') CALL extrapolate_melt_from_FCMP_to_PMP( grid, ice, BMB)
    
    
    ! If desired, the BMB can be tuned locally by changing the amplification factor
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      amplification_factor = 1._dp
    
      IF (C%choice_BMB_shelf_amplification == 'uniform') THEN
        amplification_factor = 1._dp
      ELSE IF (C%choice_BMB_shelf_amplification == 'basin') THEN
        IF (region_name == 'ANT') THEN
          amplification_factor = C%basin_BMB_amplification_factor_ANT( ice%basin_ID( j,i))
        ELSE IF (region_name == 'GRL') THEN
          amplification_factor = C%basin_BMB_amplification_factor_GRL( ice%basin_ID( j,i))
        ELSE
          amplification_factor = 1._dp
        END IF 
      ELSE
        CALL crash('unknown choice_BMB_shelf_amplification "' // TRIM(C%choice_BMB_shelf_amplification) // '"!')
      END IF       
    
      BMB%BMB_shelf( j,i) = amplification_factor * BMB%BMB_shelf( j,i)
    
    END DO
    END DO
    CALL sync    
    
    ! Add sheet and shelf melt rates together, applying the selected scheme for sub-grid shelf melt
    ! (see Leguy et al. 2021 for explanations of the three schemes)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! No sub-grid scaling for sub-sheet melt yet
      BMB%BMB( j,i) = 0._dp
      IF (ice%mask_sheet_a( j,i) == 1._dp) BMB%BMB( j,i) = BMB%BMB_sheet( j,i)
      
      ! Different sub-grid schemes for sub-shelf melt
      IF     (C%choice_BMB_subgrid == 'FCMP') THEN
        IF (ice%mask_shelf_a( j,i) == 1) BMB%BMB( j,i) = BMB%BMB( j,i) + BMB%BMB_shelf( j,i)
      ELSEIF (C%choice_BMB_subgrid == 'PMP') THEN
        BMB%BMB( j,i) = BMB%BMB( j,i) + (1._dp - ice%f_grnd_a( j,i)) * BMB%BMB_shelf( j,i)
      ELSEIF (C%choice_BMB_subgrid == 'NMP') THEN
        IF (ice%f_grnd_a( j,i) == 0._dp) BMB%BMB( j,i) = BMB%BMB( j,i) + BMB%BMB_shelf( j,i)
      ELSE
        CALL crash('unknown choice_BMB_subgrid "' // TRIM(C%choice_BMB_subgrid) // '"!')
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Limit basal melt
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      BMB%BMB( j,i) = MIN( C%BMB_max, MAX( -C%BMB_max, BMB%BMB( j,i) ))
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( BMB%BMB_sheet, 'BMB%BMB_sheet')
    CALL check_for_NaN_dp_2D( BMB%BMB_shelf, 'BMB%BMB_shelf')
    CALL check_for_NaN_dp_2D( BMB%BMB,       'BMB%BMB'      )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE run_BMB_model
  SUBROUTINE initialise_BMB_model( grid, ice, BMB, region_name)
    ! Allocate memory for the data fields of the SMB model.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (par%master) WRITE (0,*) '  Initialising BMB model: sheet = "', TRIM(C%choice_BMB_sheet_model), '", shelf = "', TRIM(C%choice_BMB_shelf_model), '"...'
    
    ! General
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%BMB      , BMB%wBMB      )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%BMB_shelf, BMB%wBMB_shelf)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%BMB_sheet, BMB%wBMB_sheet)
    
    ! Shelf
    IF     (C%choice_BMB_shelf_model == 'uniform' .OR. &
            C%choice_BMB_shelf_model == 'idealised') THEN
      ! Nothing else needs to be done
    ELSEIF (C%choice_BMB_shelf_model == 'ANICE_legacy') THEN
      CALL initialise_BMB_model_ANICE_legacy( grid, BMB, region_name)
    ELSEIF (C%choice_BMB_shelf_model == 'Favier2019_lin' .OR. &
            C%choice_BMB_shelf_model == 'Favier2019_quad' .OR. &
            C%choice_BMB_shelf_model == 'Favier2019_Mplus') THEN
      CALL initialise_BMB_model_Favier2019( grid, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'Lazeroms2018_plume') THEN
      CALL initialise_BMB_model_Lazeroms2018_plume( grid, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'PICO') THEN
      CALL initialise_BMB_model_PICO(  grid, ice, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'PICOP') THEN
      CALL initialise_BMB_model_PICOP( grid, ice, BMB)
    ELSEIF (C%choice_BMB_shelf_model == 'inverse_shelf_geometry') THEN
      BMB%BMB_shelf( :,grid%i1:grid%i2) = C%BMB_shelf_uniform
      CALL sync
    ELSE ! IF     (C%choice_BMB_shelf_model == 'uniform') THEN
      CALL crash('unknown choice_BMB_shelf_model "' // TRIM(C%choice_BMB_shelf_model) // '"!')
    END IF
    
    ! Sheet
    IF     (C%choice_BMB_sheet_model == 'uniform') THEN
      ! Nothing else needs to be done
    ELSE ! IF     (C%choice_BMB_sheet_model == 'uniform') THEN
      CALL crash('unknown choice_BMB_sheet_model "' // TRIM(C%choice_BMB_sheet_model) // '"!')
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
      
  END SUBROUTINE initialise_BMB_model
  
! == Idealised BMB schemes
! ========================

  SUBROUTINE run_BMB_model_idealised( grid, ice, BMB, time)
    ! Idealised BMB schemes
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_idealised'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF     (C%choice_idealised_BMB_shelf == 'MISMIP+') THEN
      ! The schematic basal melt used in the MISMIPplus experiments
      
      CALL run_BMB_model_MISMIPplus( grid, ice, BMB, time)
      
    ELSE
      CALL crash('unknown choice_idealised_BMB_shelf "' // TRIM(C%choice_idealised_BMB_shelf) // '"!')
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE run_BMB_model_idealised
  SUBROUTINE run_BMB_model_MISMIPplus( grid, ice, BMB, time)
    ! The schematic basal melt used in the MISMIPplus experiments
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_MISMIPplus'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: zd, cavity_thickness
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF     (C%MISMIPplus_scenario == 'ice0') THEN
      ! The reference scenario; no basal melt ever
      
      BMB%BMB(       :,grid%i1:grid%i2) = 0._dp
      BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
      BMB%BMB_shelf( :,grid%i1:grid%i2) = 0._dp
      CALL sync
      
    ELSEIF (C%MISMIPplus_scenario == 'ice1r' .OR. C%MISMIPplus_scenario == 'ice1ra') THEN
      ! Increased melt for 100 yr, followed by zero melt for 100 yr
      
      IF (time < 0._dp) THEN
        ! t = -10,000 to t = 0: spin-up, no melt
        
        BMB%BMB(       :,grid%i1:grid%i2) = 0._dp
        BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
        BMB%BMB_shelf( :,grid%i1:grid%i2) = 0._dp
        CALL sync
        
      ELSEIF (time < 100._dp) THEN
        ! t = 0 to t = 100: melt
        
        BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          
          zd = ice%Hs_a( j,i) - ice%Hi_a( j,i)
          cavity_thickness = MAX( 0._dp, zd - ice%Hb_a( j,i))
          
          ! Cornford et al. (2020), Eq. 7
          BMB%BMB_shelf( j,i) = -0.2_dp * TANH( cavity_thickness / 75._dp) * MAX( -100._dp - zd, 0._dp)
          
        END DO
        END DO
        CALL sync
        
      ELSE ! IF (time < 0._dp) THEN
        ! After t = 100: no melt
        
        BMB%BMB(       :,grid%i1:grid%i2) = 0._dp
        BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
        BMB%BMB_shelf( :,grid%i1:grid%i2) = 0._dp
        CALL sync
        
      END IF ! IF (time < 0._dp) THEN
      
    ELSEIF (C%MISMIPplus_scenario == 'ice1r' .OR. C%MISMIPplus_scenario == 'ice1rr') THEN
      ! Increased melt forever
      
      IF (time < 0._dp) THEN
        ! t = -10,000 to t = 0: spin-up, no melt
        
        BMB%BMB(       :,grid%i1:grid%i2) = 0._dp
        BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
        BMB%BMB_shelf( :,grid%i1:grid%i2) = 0._dp
        CALL sync
        
      ELSE ! IF (time < 0._dp) THEN
        ! t > 0: melt
        
        BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          
          zd = ice%Hs_a( j,i) - ice%Hi_a( j,i)
          cavity_thickness = MAX( 0._dp, zd - ice%Hb_a( j,i))
          
          ! Cornford et al. (2020), Eq. 7
          BMB%BMB_shelf( j,i) = -0.2_dp * TANH( cavity_thickness / 75._dp) * MAX( -100._dp - zd, 0._dp)
          
        END DO
        END DO
        CALL sync
        
      END IF ! IF (time < 0._dp) THEN
      
    ELSEIF (C%MISMIPplus_scenario == 'ice2r' .OR. C%MISMIPplus_scenario == 'ice2ra') THEN
      ! Increased "calving" for 100 yr, followed by zero "calving" for 100 yr
      
      IF (time < 0._dp) THEN
        ! t = -10,000 to t = 0: spin-up, no "calving"
        
        BMB%BMB(       :,grid%i1:grid%i2) = 0._dp
        BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
        BMB%BMB_shelf( :,grid%i1:grid%i2) = 0._dp
        CALL sync
        
      ELSEIF (time < 100._dp) THEN
        ! t = 0 to t = 100: "calving"
        
        BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          IF (grid%x( i) > 80000._dp) THEN ! Actually the border is at x = 480 km, but our coordinate system is shifted...
            BMB%BMB_shelf( j,i) = -100._dp
          ELSE
            BMB%BMB_shelf( j,i) = 0._dp
          END IF
        END DO
        END DO
        CALL sync
        
      ELSE ! IF (time < 0._dp) THEN
        ! After t = 100: no "calving"
        
        BMB%BMB(       :,grid%i1:grid%i2) = 0._dp
        BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
        BMB%BMB_shelf( :,grid%i1:grid%i2) = 0._dp
        CALL sync
        
      END IF ! IF (time < 0._dp) THEN
      
    ELSEIF (C%MISMIPplus_scenario == 'ice2r' .OR. C%MISMIPplus_scenario == 'ice2rr') THEN
      ! Increased "calving" forever
      
      IF (time < 0._dp) THEN
        ! t = -10,000 to t = 0: spin-up, no "calving"
        
        BMB%BMB(       :,grid%i1:grid%i2) = 0._dp
        BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
        BMB%BMB_shelf( :,grid%i1:grid%i2) = 0._dp
        CALL sync
        
      ELSE
        ! t > 0: "calving"
        
        BMB%BMB_sheet( :,grid%i1:grid%i2) = 0._dp
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          IF (grid%x( i) > 80000._dp) THEN ! Actually the border is at x = 480 km, but our coordinate system is shifted...
            BMB%BMB_shelf( j,i) = -100._dp
          ELSE
            BMB%BMB_shelf( j,i) = 0._dp
          END IF
        END DO
        END DO
        CALL sync
        
      END IF ! IF (time < 0._dp) THEN
      
    ELSE
      CALL crash('unknown MISMIPplus_scenario "' // TRIM(C%MISMIPplus_scenario) // '"!')
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE run_BMB_model_MISMIPplus

! == The ANICE_legacy sub-shelf melt model
! ========================================

  SUBROUTINE run_BMB_model_ANICE_legacy( grid, ice, BMB, region_name, time)
    ! Calculate sub-shelf melt with the ANICE_legacy model, which is based on the glacial-interglacial
    ! parameterisation by Pollard & DeConto (2012), the distance-to-open-ocean + subtended-angle parameterisation
    ! by Pollard & DeConto (2012), and the linear temperature-melt relation by Martin et al. (2011).
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_ANICE_legacy'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: T_ocean_mean, Q_TOA_jun_65N, Q_TOA_jan_80S
    REAL(dp)                                           :: BMB_shelf                             ! Sub-shelf melt rate for non-exposed shelf  [m/year]
    REAL(dp)                                           :: BMB_shelf_exposed                     ! Sub-shelf melt rate for exposed shelf      [m/year]
    REAL(dp)                                           :: BMB_deepocean                         ! Sub-shelf melt rate for deep-ocean areas   [m/year]
    REAL(dp)                                           :: w_ins, w_PD, w_warm, w_cold, w_deep, w_expo
    REAL(dp)                                           :: T_freeze                              ! Freezing temperature at the base of the shelf (Celcius)
    REAL(dp)                                           :: water_depth
    REAL(dp), PARAMETER                                :: cp0        = 3974._dp                 ! specific heat capacity of the ocean mixed layer (J kg-1 K-1) 
    REAL(dp), PARAMETER                                :: gamma_T    = 1.0E-04_dp               ! Thermal exchange velocity (m s-1)
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    ! Initialise
    BMB%sub_angle( :,grid%i1:grid%i2) = 360._dp
    BMB%dist_open( :,grid%i1:grid%i2) = 0._dp
    w_ins                             = 0._dp
    w_PD                              = 0._dp
    w_warm                            = 0._dp
    w_cold                            = 0._dp
    w_deep                            = 0._dp
    w_expo                            = 0._dp
    BMB_shelf                         = 0._dp
    BMB_shelf_exposed                 = 0._dp
    BMB_deepocean                     = 0._dp
    CALL sync
    
    ! Find the "subtended angle" and distance-to-open-ocean of all shelf pixels
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%mask_shelf_a( j,i) == 1) THEN
        BMB%sub_angle( j,i) = subtended_angle(     grid, i, j, ice%mask_land_a, ice%mask_ocean_a, ice%mask_ice_a, ice%mask_sheet_a, ice%mask_shelf_a)
        BMB%dist_open( j,i) = distance_open_ocean( grid, i, j, ice%mask_land_a, ice%mask_ocean_a, ice%mask_ice_a,                   ice%mask_shelf_a)
      END IF
    END DO
    END DO
    CALL sync
    
    ! Find the weight from insolation
    IF (region_name == 'NAM' .OR. region_name == 'EAS' .OR. region_name == 'GRL') THEN
      CALL get_insolation_at_time_month_and_lat( time, 6, 65._dp, Q_TOA_jun_65N)
      w_ins = MAX(0._dp, (Q_TOA_jun_65N - 462.29_dp) / 40._dp)
    ELSEIF (region_name == 'ANT') THEN
      CALL get_insolation_at_time_month_and_lat( time, 1, -80._dp, Q_TOA_jan_80S)
      w_ins = MAX(0._dp, (Q_TOA_jan_80S - 532.19_dp) / 40._dp)
    END IF
    
    ! Initialise to prevent compiler warnings
    T_ocean_mean = 0._dp
    
    ! Determine mean ocean temperature and basal melt rates for deep ocean and exposed shelves
    IF (C%choice_forcing_method == 'CO2_direct') THEN
    
      ! Use the prescribed CO2 record as a glacial index
      IF (forcing%CO2_obs > 280._dp) THEN
        ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 400 ppmv
        w_PD   = MIN(1.25_dp, MAX(-0.25_dp, (400._dp - forcing%CO2_obs) / (400._dp - 280._dp) )) - w_ins
        w_warm = 1._dp - w_PD
        w_cold = 0._dp
      ELSE
        ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
        w_PD   = MIN(1.25_dp, MAX(-0.25_dp, (forcing%CO2_obs - 190._dp) / (280._dp - 190._dp) )) + w_ins
        w_cold = 1._dp - w_PD
        w_warm = 0._dp
      END IF
      
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
    
      ! Use modelled CO2 as a glacial index
      IF (forcing%CO2_mod > 280._dp) THEN
        ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 400 ppmv
        w_PD   = MIN(1.25_dp, MAX(-0.25_dp, (400._dp - forcing%CO2_mod) / (400._dp - 280._dp) )) - w_ins
        w_warm = 1._dp - w_PD
        w_cold = 0._dp
      ELSE
        ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
        w_PD   = MIN(1.25_dp, MAX(-0.25_dp, (forcing%CO2_mod - 190._dp) / (280._dp - 190._dp) )) + w_ins
        w_cold = 1._dp - w_PD
        w_warm = 0._dp
      END IF
      
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
    
      ! Use modelled global mean annual temperature change as a glacial index
      IF (forcing%dT_glob_inverse > 0._dp) THEN
        ! Warmer than present, interpolate between "PD" and "warm", assuming "warm" means 405 ppmv
        w_warm = MIN(1.25_dp, MAX(-0.25_dp, forcing%dT_glob_inverse / 12._dp )) - w_ins
        w_PD   = 1._dp - w_warm
        w_cold = 0._dp
      ELSE
        ! Colder than present, interpolate between "PD" and "cold", assuming "cold" means 190 ppmv
        w_cold = MIN(1.25_dp, MAX(-0.25_dp, -forcing%dT_glob_inverse / 12._dp )) + w_ins
        w_PD   = 1._dp - w_cold
        w_warm = 0._dp
      END IF
      
    ELSEIF (C%choice_forcing_method == 'climate_direct' .OR. C%choice_forcing_method == 'SMB_direct') THEN
      ! In this case, no CO2/d18O forcing is used; just assume PD weights
      
      w_warm = 0._dp
      w_cold = 0._dp
      w_PD   = 1._dp
      
    ELSE ! IF (C%choice_forcing_method == 'CO2_direct') THEN
      CALL crash('unknown choice_forcing_method "' // TRIM(C%choice_forcing_method) // '"!')
    END IF ! IF (C%choice_forcing_method == 'CO2_direct') THEN
    
    T_ocean_mean      = w_PD * BMB%T_ocean_mean_PD      + w_warm * BMB%T_ocean_mean_warm      + w_cold * BMB%T_ocean_mean_cold
    BMB_deepocean     = w_PD * BMB%BMB_deepocean_PD     + w_warm * BMB%BMB_deepocean_warm     + w_cold * BMB%BMB_deepocean_cold
    BMB_shelf_exposed = w_PD * BMB%BMB_shelf_exposed_PD + w_warm * BMB%BMB_shelf_exposed_warm + w_cold * BMB%BMB_shelf_exposed_cold
    CALL sync
    
    ! Use the (interpolated, spatially uniform) ocean temperature and the subtended angle + distance-to-open-ocean
    ! to calculate sub-shelf melt rates using the parametrisation from Martin et al., 2011
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      IF (ice%mask_shelf_a( j,i) == 1) THEN
        ! Sub-shelf melt

        ! Freezing temperature at the bottom of the ice shelves, scaling with depth below water level
        T_freeze = 0.0939_dp - 0.057_dp * 35._dp - 7.64E-04_dp * ice%Hi_a( j,i) * ice_density / seawater_density

        ! Sub-shelf melt rate for non-exposed shelves (Martin, TC, 2011) - melt values, when T_ocean > T_freeze.
        BMB_shelf   = seawater_density * cp0 * sec_per_year * gamma_T * BMB%subshelf_melt_factor * &
                   (T_ocean_mean - T_freeze) / (L_fusion * ice_density)

      ELSE
        BMB_shelf = 0._dp
      END IF

      IF (ice%mask_shelf_a( j,i) == 1 .OR. ice%mask_ocean_a( j,i) == 1) THEN
      
        water_depth = ice%SL_a( j,i) - ice%Hb_a( j,i)
        w_deep = MAX(0._dp, MIN(1._dp, (water_depth - BMB%deep_ocean_threshold_depth) / 200._dp))
        w_expo = MAX(0._dp, MIN(1._dp, (BMB%sub_angle( j,i) - 80._dp)/30._dp)) * EXP(-BMB%dist_open( j,i) / 100000._dp)
        
        BMB%BMB_shelf( j,i) = (w_deep * BMB_deepocean) + (1._dp - w_deep) * (w_expo * BMB_shelf_exposed + (1._dp - w_expo) * BMB_shelf)
        
      ELSE  
        BMB%BMB_shelf( j,i) = 0._dp
      END IF

    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( BMB%sub_angle, 'BMB%sub_angle')
    CALL check_for_NaN_dp_2D( BMB%dist_open, 'BMB%dist_open')
    CALL check_for_NaN_dp_2D( BMB%BMB_shelf, 'BMB%BMB_shelf')
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
          
  END SUBROUTINE run_BMB_model_ANICE_legacy
  FUNCTION distance_open_ocean( grid, i_shelf, j_shelf, mask_land, mask_ocean, mask_ice, mask_shelf) RESULT( open_distance)
    ! Determine the distance to the open ocean (expressed in number of grid cells). A  solution is found for all 16 directions
    ! also when encountering land points, keep searching to find a grid point. 
    ! The minimum value is the solution

    ! The 16 directions
    !
    !              12 13 14
    !             11 ---- 15
    !            10 --  -- 16
    !           09 -- ij -- 01
    !            08 --  -- 02 
    !             07 ---- 03
    !              06 05 04 
    !

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    INTEGER,                             INTENT(IN)    :: j_shelf, i_shelf
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_land
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_ocean
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_ice
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_shelf
    REAL(dp)                                           :: open_distance
    
    ! Local variables:
    REAL(dp), DIMENSION(16)                            :: distance             ! list of all distances for 16 directions (m)
    INTEGER                                            :: iloop,imax           ! loop integer and maximum value (max of NX/NY)
    INTEGER                                            :: ii,jj
    
    imax          = MAX( grid%nx, grid%ny)
    distance      = 0
    
    IF (mask_shelf( j_shelf,i_shelf) == 0) THEN
      IF     (mask_land(  j_shelf,i_shelf) == 1) THEN
        open_distance = REAL(grid%nx+grid%ny,dp) * grid%dx
      ELSE
        open_distance = 0._dp
      END IF
      RETURN
    END IF
   
    ! direction 01: (should work fine with imax, no possibility that jj or ii will exceed nx or ny)
    ! but what about when jj or ii are lower than 1..
    DO iloop = 1,imax
      ii = i_shelf
      jj = j_shelf + iloop
       
        ! first check if the boundaries of the grid are reached
        IF (jj > grid%ny) THEN
          distance(1) = 100._dp
          EXIT
        END IF  
        
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(1) = real(jj) - real(j_shelf)
          EXIT
        END IF    
      END DO

      ! direction 02:
      DO iloop = 1,imax   
        ii = i_shelf - iloop
        jj = j_shelf + 2*iloop
        IF (jj > grid%ny .OR. ii < 1) THEN
          distance(2) = 100._dp
          EXIT
        END IF  
        
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(2) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 03:
      DO iloop = 1,imax   
        ii = i_shelf - iloop
        jj = j_shelf + iloop
        IF (jj > grid%ny .OR. ii < 1) THEN
          distance(3) = 100._dp
          EXIT
        END IF  

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(3) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 04:
      DO iloop = 1,imax   
        ii = i_shelf - 2*iloop
        jj = j_shelf + iloop
        IF (jj > grid%ny .OR. ii < 1) THEN
          distance(4) = 100._dp
          EXIT
        END IF  
        
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(4) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 05:
      DO iloop = 1,imax   
        ii = i_shelf - iloop
        jj = j_shelf
        IF ( ii < 1) THEN
          distance(5) = 100._dp
          EXIT
        END IF  

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(5) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 06:
      DO iloop = 1,imax   
        ii = i_shelf - 2*iloop
        jj = j_shelf - iloop
        IF (jj < 1 .OR. ii < 1) THEN
          distance(6) = 100._dp
          EXIT
        END IF
        
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(6) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 07:
      DO iloop = 1,imax   
        ii = i_shelf - iloop
        jj = j_shelf - iloop
        IF (jj < 1 .OR. ii < 1) THEN
          distance(7) = 100._dp
          EXIT
        END IF 
                 
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(7) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 08:
      DO iloop = 1,imax   
        ii = i_shelf - iloop
        jj = j_shelf - 2*iloop
        IF (jj < 1 .OR. ii < 1) THEN
          distance(8) = 100._dp
          EXIT
        END IF
                  
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(8) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 09:
      DO iloop = 1,imax   
        ii = i_shelf 
        jj = j_shelf - iloop
        IF (jj < 1) THEN
          distance(9) = 100._dp
          EXIT
        END IF
                  
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(9) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 10:
      DO iloop = 1,imax   
        ii = i_shelf + iloop 
        jj = j_shelf - 2*iloop
        IF (jj < 1 .OR. ii > grid%nx) THEN
          distance(10) = 100._dp
          EXIT
        END IF
                  
        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(10) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 11:
      DO iloop = 1,imax   
        ii = i_shelf + iloop 
        jj = j_shelf - iloop
        IF (jj < 1 .OR. ii > grid%nx) THEN
          distance(11) = 100._dp
          EXIT
        END IF

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(11) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 12:
      DO iloop = 1,imax   
        ii = i_shelf + 2*iloop 
        jj = j_shelf - iloop
        IF (jj < 1 .OR. ii > grid%nx) THEN
          distance(12) = 100._dp
          EXIT
        END IF

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(12) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 13:
      DO iloop = 1,imax   
        ii = i_shelf + iloop 
        jj = j_shelf
        IF (ii > grid%nx) THEN
          distance(13) = 100._dp
          EXIT
        END IF

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(13) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 14:
      DO iloop = 1,imax   
        ii = i_shelf + 2*iloop 
        jj = j_shelf + iloop
        IF (jj > grid%ny .OR. ii > grid%nx) THEN
          distance(14) = 100._dp
          EXIT
        END IF

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(14) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 15:
      DO iloop = 1,imax   
        ii = i_shelf + iloop 
        jj = j_shelf + iloop
        IF (jj > grid%ny .OR. ii > grid%nx) THEN
          distance(15) = 100._dp
          EXIT
        END IF

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(15) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! direction 16:
      DO iloop = 1,imax   
        ii = i_shelf + iloop 
        jj = j_shelf + 2*iloop
        IF (jj > grid%ny .OR. ii > grid%nx) THEN
          distance(16) = 100._dp
          EXIT
        END IF

        IF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
          distance(16) = SQRT( (real(jj) - real(j_shelf))**2 + (real(ii) - real(i_shelf))**2 )
          EXIT
        END IF    
      END DO

      ! Calculate minimum distance
      open_distance = MINVAL(distance) * grid%dx

    END FUNCTION distance_open_ocean
  FUNCTION subtended_angle( grid, i_shelf, j_shelf, mask_land, mask_ocean, mask_ice, mask_sheet, mask_shelf) RESULT(angle_sub)
    ! Determine the subtended angle to the open ocean. When encountering land points
    ! the angle is zet to zero (i.e. a loss of 1/16th of the fully possible 360). 
    ! The minimum value is the solution
    
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    INTEGER,                             INTENT(IN)    :: j_shelf, i_shelf
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_land
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_ocean
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_ice
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_sheet
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_shelf
    REAL(dp)                                           :: angle_sub
    
    ! Local variables:
    REAL(dp), DIMENSION(16)                            :: angle                !  list of open-ocean angles for 16 directions (degrees)
    INTEGER                                            :: iloop,imax           ! loop integer and maximum value (max of NX/NY)
    INTEGER                                            :: ii,jj

    imax      = MAX(grid%nx,grid%ny)
    angle     = 1._dp
    angle_sub = 360._dp
    
    IF (mask_shelf( j_shelf,i_shelf) == 0) THEN
      IF     (mask_land(  j_shelf,i_shelf) == 1) THEN
        angle_sub = 0._dp
      ELSE
        angle_sub = 360._dp
      END IF
      RETURN
    END IF
    
    ! direction 01:
    DO iloop = 1,imax   
      ii = i_shelf
      jj = j_shelf + iloop
      
      ! first check if the boundaries of the grid are reached
      IF (jj > grid%ny) THEN
        angle(1) = 0._dp
        EXIT
      END IF  
      
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(1) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(1) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 02:
    DO iloop = 1,imax   
      ii = i_shelf - iloop
      jj = j_shelf + 2*iloop
      IF (jj > grid%ny .OR. ii < 1) THEN
        angle(2) = 0._dp
        EXIT
      END IF  
      
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(2) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(2) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 03:
    DO iloop = 1,imax   
      ii = i_shelf - iloop
      jj = j_shelf + iloop
      IF (jj > grid%ny .OR. ii < 1) THEN
        angle(3) = 0._dp
        EXIT
      END IF  

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(3) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(3) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 04:
    DO iloop = 1,imax   
      ii = i_shelf - 2*iloop
      jj = j_shelf + iloop
      IF (jj > grid%ny .OR. ii < 1) THEN
        angle(4) = 0._dp
        EXIT
      END IF  
      
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(4) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(4) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 05:
    DO iloop = 1,imax   
      ii = i_shelf - iloop
      jj = j_shelf
      IF ( ii < 1) THEN
        angle(5) = 0._dp
        EXIT
      END IF  

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(5) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(5) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 06:
    DO iloop = 1,imax   
      ii = i_shelf - 2*iloop
      jj = j_shelf - iloop
      IF (jj < 1 .OR. ii < 1) THEN
        angle(6) = 0._dp
        EXIT
      END IF
      
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(6) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(6) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 07:
    DO iloop = 1,imax   
      ii = i_shelf - iloop
      jj = j_shelf - iloop
      IF (jj < 1 .OR. ii < 1) THEN
        angle(7) = 0._dp
        EXIT
      END IF 
               
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(7) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(7) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 08:
    DO iloop = 1,imax   
      ii = i_shelf - iloop
      jj = j_shelf - 2*iloop
      IF (jj < 1 .OR. ii < 1) THEN
        angle(8) = 0._dp
        EXIT
      END IF
                
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(8) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(8) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 09:
    DO iloop = 1,imax   
      ii = i_shelf 
      jj = j_shelf - iloop
      IF (jj < 1) THEN
        angle(9) = 0._dp
        EXIT
      END IF
                
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(9) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(9) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 10:
    DO iloop = 1,imax   
      ii = i_shelf + iloop 
      jj = j_shelf - 2*iloop
      IF (jj < 1 .OR. ii > grid%nx) THEN
        angle(10) = 0._dp
        EXIT
      END IF
                
      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(10) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(10) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 11:
    DO iloop = 1,imax   
      ii = i_shelf + iloop 
      jj = j_shelf - iloop
      IF (jj < 1 .OR. ii > grid%nx) THEN
        angle(11) = 0._dp
        EXIT
      END IF

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(11) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(11) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 12:
    DO iloop = 1,imax   
      ii = i_shelf + 2*iloop 
      jj = j_shelf - iloop
      IF (jj < 1 .OR. ii > grid%nx) THEN
        angle(12) = 0._dp
        EXIT
      END IF

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(12) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(12) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 13:
    DO iloop = 1,imax   
      ii = i_shelf + iloop 
      jj = j_shelf
      IF (ii > grid%nx) THEN
        angle(13) = 0._dp
        EXIT
      END IF

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(13) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(13) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 14:
    DO iloop = 1,imax   
      ii = i_shelf + 2*iloop 
      jj = j_shelf + iloop
      IF (jj > grid%ny .OR. ii > grid%nx) THEN
        angle(14) = 0._dp
        EXIT
      END IF

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(14) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(14) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 15:
    DO iloop = 1,imax   
      ii = i_shelf + iloop 
      jj = j_shelf + iloop
      IF (jj > grid%ny .OR. ii > grid%nx) THEN
        angle(15) = 0._dp
        EXIT
      END IF

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(15) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(15) = 1._dp
        EXIT
      END IF    
    END DO

    ! direction 16:
    DO iloop = 1,imax   
      ii = i_shelf + iloop 
      jj = j_shelf + 2*iloop
      IF (jj > grid%ny .OR. ii > grid%nx) THEN
        angle(16) = 0._dp
        EXIT
      END IF

      IF (mask_sheet( jj,ii) == 1 .OR. mask_land( jj,ii) == 1) THEN
        angle(16) = 0._dp
        EXIT
      ELSEIF (mask_ocean( jj,ii) == 1 .AND. mask_ice( jj,ii) == 0) THEN
        angle(16) = 1._dp
        EXIT
      END IF    
    END DO

    ! Calculate total subtanded angle
    angle_sub = SUM(angle(1:16)) * 360._dp / 16._dp

  END FUNCTION subtended_angle
  SUBROUTINE initialise_BMB_model_ANICE_legacy( grid, BMB, region_name)
    ! Allocate memory for the data fields of the ANICE_legacy shelf BMB model.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model_ANICE_legacy'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Variables
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%BMB_shelf, BMB%wBMB_shelf)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%sub_angle, BMB%wsub_angle)
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%dist_open, BMB%wdist_open)
    
    ! Tuning parameters
    CALL allocate_shared_dp_0D( BMB%T_ocean_mean_PD,            BMB%wT_ocean_mean_PD           )
    CALL allocate_shared_dp_0D( BMB%T_ocean_mean_cold,          BMB%wT_ocean_mean_cold         )
    CALL allocate_shared_dp_0D( BMB%T_ocean_mean_warm,          BMB%wT_ocean_mean_warm         )
    CALL allocate_shared_dp_0D( BMB%BMB_deepocean_PD,           BMB%wBMB_deepocean_PD          )
    CALL allocate_shared_dp_0D( BMB%BMB_deepocean_cold,         BMB%wBMB_deepocean_cold        )
    CALL allocate_shared_dp_0D( BMB%BMB_deepocean_warm,         BMB%wBMB_deepocean_warm        )
    CALL allocate_shared_dp_0D( BMB%BMB_shelf_exposed_PD,       BMB%wBMB_shelf_exposed_PD      )
    CALL allocate_shared_dp_0D( BMB%BMB_shelf_exposed_cold,     BMB%wBMB_shelf_exposed_cold    )
    CALL allocate_shared_dp_0D( BMB%BMB_shelf_exposed_warm,     BMB%wBMB_shelf_exposed_warm    )
    CALL allocate_shared_dp_0D( BMB%subshelf_melt_factor,       BMB%wsubshelf_melt_factor      )
    CALL allocate_shared_dp_0D( BMB%deep_ocean_threshold_depth, BMB%wdeep_ocean_threshold_depth)
    
    IF (region_name == 'NAM') THEN
      BMB%T_ocean_mean_PD            = C%T_ocean_mean_PD_NAM
      BMB%T_ocean_mean_cold          = C%T_ocean_mean_cold_NAM
      BMB%T_ocean_mean_warm          = C%T_ocean_mean_warm_NAM
      BMB%BMB_deepocean_PD           = C%BMB_deepocean_PD_NAM
      BMB%BMB_deepocean_cold         = C%BMB_deepocean_cold_NAM
      BMB%BMB_deepocean_warm         = C%BMB_deepocean_warm_NAM
      BMB%BMB_shelf_exposed_PD       = C%BMB_shelf_exposed_PD_NAM
      BMB%BMB_shelf_exposed_cold     = C%BMB_shelf_exposed_cold_NAM
      BMB%BMB_shelf_exposed_warm     = C%BMB_shelf_exposed_warm_NAM
      BMB%subshelf_melt_factor       = C%subshelf_melt_factor_NAM
      BMB%deep_ocean_threshold_depth = C%deep_ocean_threshold_depth_NAM
    ELSEIF (region_name == 'EAS') THEN
      BMB%T_ocean_mean_PD            = C%T_ocean_mean_PD_EAS
      BMB%T_ocean_mean_cold          = C%T_ocean_mean_cold_EAS
      BMB%T_ocean_mean_warm          = C%T_ocean_mean_warm_EAS
      BMB%BMB_deepocean_PD           = C%BMB_deepocean_PD_EAS
      BMB%BMB_deepocean_cold         = C%BMB_deepocean_cold_EAS
      BMB%BMB_deepocean_warm         = C%BMB_deepocean_warm_EAS
      BMB%BMB_shelf_exposed_PD       = C%BMB_shelf_exposed_PD_EAS
      BMB%BMB_shelf_exposed_cold     = C%BMB_shelf_exposed_cold_EAS
      BMB%BMB_shelf_exposed_warm     = C%BMB_shelf_exposed_warm_EAS
      BMB%subshelf_melt_factor       = C%subshelf_melt_factor_EAS
      BMB%deep_ocean_threshold_depth = C%deep_ocean_threshold_depth_EAS
    ELSEIF (region_name == 'GRL') THEN
      BMB%T_ocean_mean_PD            = C%T_ocean_mean_PD_GRL
      BMB%T_ocean_mean_cold          = C%T_ocean_mean_cold_GRL
      BMB%T_ocean_mean_warm          = C%T_ocean_mean_warm_GRL
      BMB%BMB_deepocean_PD           = C%BMB_deepocean_PD_GRL
      BMB%BMB_deepocean_cold         = C%BMB_deepocean_cold_GRL
      BMB%BMB_deepocean_warm         = C%BMB_deepocean_warm_GRL
      BMB%BMB_shelf_exposed_PD       = C%BMB_shelf_exposed_PD_GRL
      BMB%BMB_shelf_exposed_cold     = C%BMB_shelf_exposed_cold_GRL
      BMB%BMB_shelf_exposed_warm     = C%BMB_shelf_exposed_warm_GRL
      BMB%subshelf_melt_factor       = C%subshelf_melt_factor_GRL
      BMB%deep_ocean_threshold_depth = C%deep_ocean_threshold_depth_GRL
    ELSEIF (region_name == 'ANT') THEN
      BMB%T_ocean_mean_PD            = C%T_ocean_mean_PD_ANT
      BMB%T_ocean_mean_cold          = C%T_ocean_mean_cold_ANT
      BMB%T_ocean_mean_warm          = C%T_ocean_mean_warm_ANT
      BMB%BMB_deepocean_PD           = C%BMB_deepocean_PD_ANT
      BMB%BMB_deepocean_cold         = C%BMB_deepocean_cold_ANT
      BMB%BMB_deepocean_warm         = C%BMB_deepocean_warm_ANT
      BMB%BMB_shelf_exposed_PD       = C%BMB_shelf_exposed_PD_ANT
      BMB%BMB_shelf_exposed_cold     = C%BMB_shelf_exposed_cold_ANT
      BMB%BMB_shelf_exposed_warm     = C%BMB_shelf_exposed_warm_ANT
      BMB%subshelf_melt_factor       = C%subshelf_melt_factor_ANT
      BMB%deep_ocean_threshold_depth = C%deep_ocean_threshold_depth_ANT
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
      
  END SUBROUTINE initialise_BMB_model_ANICE_legacy
  
! == The Favier et al. (2019) sub-shelf melt parameterisations
! ============================================================

  SUBROUTINE run_BMB_model_Favier2019_linear(         grid, ice, ocean, BMB)
    ! Calculate sub-shelf melt with Favier et al. (2019) linear parameterisation
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_Favier2019_linear'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: dT
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Calculate ocean temperature and freezing point at the base of the shelf
    CALL calc_ocean_temperature_at_shelf_base(    grid, ice, ocean, BMB)
    CALL calc_ocean_freezing_point_at_shelf_base( grid, ice, ocean, BMB)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! Initialise
      BMB%BMB_shelf( j,i) = 0._dp
      
      IF (ice%mask_shelf_a( j,i) == 1) THEN
        
        ! Temperature forcing
        dT = BMB%T_ocean_base( j,i) - BMB%T_ocean_freeze_base( j,i)
        
        ! Favier et al. (2019), Eq. 2
        BMB%BMB_shelf( j,i) = -sec_per_year * C%BMB_Favier2019_lin_GammaT * seawater_density * cp_ocean * dT / (ice_density * L_fusion)
        
      END IF ! IF (ice%mask_shelf_a( j,i) == 1) THEN
      
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
          
  END SUBROUTINE run_BMB_model_Favier2019_linear
  SUBROUTINE run_BMB_model_Favier2019_quadratic(      grid, ice, ocean, BMB)
    ! Calculate sub-shelf melt with Favier et al. (2019) quadratic parameterisation
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_Favier2019_quadratic'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: dT
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Calculate ocean temperature and freezing point at the base of the shelf
    CALL calc_ocean_temperature_at_shelf_base(    grid, ice, ocean, BMB)
    CALL calc_ocean_freezing_point_at_shelf_base( grid, ice, ocean, BMB)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! Initialise
      BMB%BMB_shelf( j,i) = 0._dp
      
      IF (ice%mask_shelf_a( j,i) == 1) THEN
        
        ! Temperature forcing
        dT = BMB%T_ocean_base( j,i) - BMB%T_ocean_freeze_base( j,i)
        
        ! Favier et al. (2019), Eq. 4
        ! BMB%BMB_shelf( j,i) = -sec_per_year * C%BMB_Favier2019_quad_GammaT * (seawater_density * cp_ocean * dT / (ice_density * L_fusion))**2._dp
        ! Altered to allow for negative basal melt (i.e. refreezing) when dT < 0 
        BMB%BMB_shelf( j,i) = -sec_per_year * C%BMB_Favier2019_quad_GammaT * SIGN(dT,1._dp) * (seawater_density * cp_ocean * dT / (ice_density * L_fusion))**2._dp
        
      END IF ! IF (ice%mask_shelf_a( j,i) == 1) THEN
      
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
          
  END SUBROUTINE run_BMB_model_Favier2019_quadratic
  SUBROUTINE run_BMB_model_Favier2019_Mplus(          grid, ice, ocean, BMB)
    ! Calculate sub-shelf melt with Favier et al. (2019) M+ parameterisation
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_Favier2019_Mplus'
    INTEGER                                            :: i,j,n_shelf,basin_i
    REAL(dp)                                           :: dT,dT_av
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Calculate ocean temperature and freezing point at the base of the shelf
    CALL calc_ocean_temperature_at_shelf_base(    grid, ice, ocean, BMB)
    CALL calc_ocean_freezing_point_at_shelf_base( grid, ice, ocean, BMB)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! Initialise
      BMB%BMB_shelf( j,i) = 0._dp
      
      IF (ice%mask_shelf_a( j,i) == 1) THEN
        
        ! Temperature forcing
        dT = BMB%T_ocean_base( j,i) - BMB%T_ocean_freeze_base( j,i)
        
        ! Favier et al. (2019), Eq. 5 (part 1, without the basin-averaged term, that comes later)
        BMB%BMB_shelf( j,i) = -sec_per_year * C%BMB_Favier2019_Mplus_GammaT * (seawater_density * cp_ocean / (ice_density * L_fusion))**2._dp * dT
        
      END IF ! IF (ice%mask_shelf_a( j,i) == 1) THEN
      
    END DO
    END DO
    CALL sync
    
    ! Calculate and apply basin-averaged temperature forcing
    DO basin_i = 1, ice%nbasins
    
      dT_av   = 0._dp
      n_shelf = 0
    
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (ice%basin_ID( j,i) == basin_i .AND. ice%mask_shelf_a( j,i) == 1) THEN
          dT = BMB%T_ocean_base( j,i) - BMB%T_ocean_freeze_base( j,i)
          dT_av   = dT_av   + dT
          n_shelf = n_shelf + 1
        END IF
      END DO
      END DO
        
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, dT_av,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_shelf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        
      ! Safety
      IF (n_shelf == 0) THEN
        dT_av = 0._dp
      ELSE
        dT_av = dT_av / n_shelf
      END IF
      
      ! Add last term to the equation
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (ice%basin_ID( j,i) == basin_i .AND. ice%mask_shelf_a( j,i) == 1) THEN
          BMB%BMB_shelf( j,i) = BMB%BMB_shelf( j,i) * dT_av
        END IF
      END DO
      END DO
      CALL sync
    
    END DO ! DO basin_i = 1, ice%nbasins
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
          
  END SUBROUTINE run_BMB_model_Favier2019_Mplus
  SUBROUTINE calc_ocean_temperature_at_shelf_base(    grid, ice, ocean, BMB)
    ! Calculate ocean temperature at the base of the shelf by interpolating
    ! the 3-D ocean temperature field in the vertical column
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_ocean_temperature_at_shelf_base'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: depth
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Initialise at zero
      BMB%T_ocean_base( j,i) = 0._dp
      
      IF (ice%mask_shelf_a( j,i) == 1) THEN
      
        ! Calculate depth
        depth = MAX( 0.1_dp, ice%Hi_a( j,i) * ice_density / seawater_density)
        
        ! Find ocean temperature at this depth
        CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( :,j,i), depth, BMB%T_ocean_base( j,i))
        
      END IF ! IF (ice%mask_shelf_a( j,i) == 1) THEN
      
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
          
  END SUBROUTINE calc_ocean_temperature_at_shelf_base
  SUBROUTINE calc_ocean_freezing_point_at_shelf_base( grid, ice, ocean, BMB)
    ! Calculate the ocean freezing point at the base of the shelf, needed to calculate
    ! basal melt in the different parameterisations from Favier et al. (2019)
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_ocean_freezing_point_at_shelf_base'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: depth
    REAL(dp)                                           :: S0                   ! Practical salinity [PSU]
    REAL(dp), PARAMETER                                :: lambda1 = -0.0575_dp ! Liquidus slope                [degC PSU^-1] (Favier et al. (2019), Table 2)
    REAL(dp), PARAMETER                                :: lambda2 = 0.0832_dp  ! Liquidus intercept            [degC]        (Favier et al. (2019), Table 2)
    REAL(dp), PARAMETER                                :: lambda3 = 7.59E-4_dp ! Liquidus pressure coefficient [degC m^-1]   (Favier et al. (2019), Table 2)
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Initialise at zero
      BMB%T_ocean_freeze_base( j,i) = 0._dp
      
      IF (ice%mask_shelf_a( j,i) == 1) THEN
      
        ! Calculate depth
        depth = MAX( 0.1_dp, ice%Hi_a( j,i) * ice_density / seawater_density)
        
        ! Find salinity at this depth
        CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S_ocean_corr_ext( :,j,i), depth, S0)
        
        ! Calculate ocean freezing temperature (Favier et al. (2019), Eq. 3) in degrees Celsius
        BMB%T_ocean_freeze_base( j,i) = lambda1 * S0 + lambda2 - lambda3 * depth
        
      END IF ! IF (ice%mask_shelf_a( j,i) == 1) THEN
      
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
          
  END SUBROUTINE calc_ocean_freezing_point_at_shelf_base
  SUBROUTINE initialise_BMB_model_Favier2019( grid, BMB)
    ! Allocate memory for the data fields of the Favier et al. (2019) shelf BMB parameterisations.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model_Favier2019'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Variables
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%T_ocean_base,        BMB%wT_ocean_base       )
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, BMB%T_ocean_freeze_base, BMB%wT_ocean_freeze_base)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
      
  END SUBROUTINE initialise_BMB_model_Favier2019
  
! == The Lazeroms et al. (2018) quasi-2-D plume parameterisation
! ==============================================================

  SUBROUTINE run_BMB_model_Lazeroms2018_plume( grid, ice, ocean, BMB)
    ! Calculate basal melt using the quasi-2-D plume parameterisation by Lazeroms et al. (2018), following the equations in Appendix A
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_Lazeroms2018_plume'
    INTEGER                                            :: i,j,basin_i
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  list_GL
    INTEGER                                            :: wlist_GL
    INTEGER,                    POINTER                ::  n_GL
    INTEGER                                            :: wn_GL
    REAL(dp)                                           :: depth                               ! Ice-shelf base depth ("draft") [m]
    REAL(dp)                                           :: Ta                                  ! Ambient temperature at the ice-shelf base [degC]
    REAL(dp)                                           :: Sa                                  ! Ambient salinity    at the ice-shelf base [PSU]
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    ! Initialise
    BMB%eff_plume_source_depth( :,grid%i1:grid%i2) = 0._dp
    BMB%eff_basal_slope(        :,grid%i1:grid%i2) = 0._dp
    BMB%BMB_shelf(              :,grid%i1:grid%i2) = 0._dp
    CALL sync
    
    ! Calculate ocean temperature at the base of the shelf
    CALL calc_ocean_temperature_at_shelf_base( grid, ice, ocean, BMB)
    
    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%nx*grid%ny, 2, list_GL, wlist_GL)
    CALL allocate_shared_int_0D(                     n_GL,     wn_GL    )
    
    ! The new look-for-plume-source approach needs to be run separately for each basin!
    DO basin_i = 1, ice%nbasins
      
      ! List all grounding-line and ice-front grid cells
      CALL list_GL_grid_cells_basin( grid, ice, ice%mask_sheet_a, ice%mask_shelf_a, basin_i, list_GL, n_GL)
      
      ! Determine the effective plume path to find the effective plume source depth and basal slope fo all shelf grid cells in this basin
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (ice%mask_shelf_a( j,i) == 1 .AND. ice%basin_ID( j,i) == basin_i) THEN
          CALL find_effective_plume_path( grid, ice, BMB, list_GL, n_GL, i,j, BMB%eff_plume_source_depth( j,i), BMB%eff_basal_slope( j,i))
        END IF
      END DO
      END DO
      CALL sync
      
    END DO ! DO basin_i = 1, ice%nbasins
      
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      IF (ice%mask_shelf_a( j,i) == 1) THEN
        
        ! Calculate the depth of the shelf base
        depth = MAX( 0.1_dp, ice%Hi_a( j,i) - ice%Hs_a( j,i))   ! Depth is positive when below the sea surface!
          
        ! Find ambient temperature and salinity at the ice-shelf base
        IF (C%choice_BMB_shelf_model == 'Lazeroms2018_plume') THEN
          ! Use the extrapolated ocean temperature+salinity fields
        
          CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( :,j,i), depth, Ta)
          CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S_ocean_corr_ext( :,j,i), depth, Sa)
          
        ELSEIF (C%choice_BMB_shelf_model == 'PICOP') THEN
          ! Use the results from the PICO ocean box model
          
          Ta = BMB%PICO_T( j,i)
          Sa = BMB%PICO_S( j,i)
          
        ELSE
          CALL crash('finding ambient temperature+salinity only defined for choice_BMB_shelf_model = "Lazeroms2018_plume" or "PICOP"!')
        END IF
        
        ! Apply the plume model
        CALL plume_model_Lazeroms2018( BMB%eff_plume_source_depth( j,i), depth, Ta, Sa, BMB%eff_basal_slope( j,i), BMB%BMB_shelf( j,i))
        
      END IF ! IF (ice%mask_shelf_a( j,i) == 1) THEN
      
    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wlist_GL)
    CALL deallocate_shared( wn_GL    )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE run_BMB_model_Lazeroms2018_plume
  SUBROUTINE plume_model_Lazeroms2018( plume_source_depth, depth, Ta, Sa, basal_slope, melt_rate)
    ! Calculate basal melt using the quasi-2-D plume parameterisation by Lazeroms et al. (2018), following the equations in Appendix A
    
    IMPLICIT NONE
    
    ! In/output variables
    REAL(dp),                            INTENT(IN)    :: plume_source_depth    ! Water depth at the grounding-line plume source   [m]
    REAL(dp),                            INTENT(IN)    :: depth                 ! Depth of the ice shelf base                      [m]
    REAL(dp),                            INTENT(IN)    :: Ta                    ! Ambient temperature at the ice-shelf base        [degC]
    REAL(dp),                            INTENT(IN)    :: Sa                    ! Ambient salinity    at the ice-shelf base        [PSU]
    REAL(dp),                            INTENT(IN)    :: basal_slope           ! Basal slope                                      [m/m]
    REAL(dp),                            INTENT(OUT)   :: melt_rate             ! Basal melt rate                                  [m/yr]
    
    ! Local variables:
    REAL(dp)                                           :: Tf_GL                               ! Freezing temperature at effective grounding-line plume source
    REAL(dp)                                           :: alpha                               ! Effective local angle ( = atan(slope))
    REAL(dp)                                           :: sinalpha                            ! sin( alpha) (appears so often that it's better to give it its own variable)
    REAL(dp)                                           :: GammaTS                             ! Effective heat exchange coefficient
    REAL(dp)                                           :: g_alpha                             ! Geometry term
    REAL(dp)                                           :: g_alpha_term1, g_alpha_term2, g_alpha_term3, sqrtCd_GammaTS
    REAL(dp)                                           :: M                                   ! Empirically derived melt-rate scale
    REAL(dp)                                           :: l_geo                               ! Geometry- and temperature-dependent length scale
    REAL(dp)                                           :: Xhat                                ! Dimensionless scaled distance along plume path
    REAL(dp)                                           :: Mhat                                ! Dimensionless melt curve
    
    ! Constant parameters (Lazeroms et al. (2018), Table 1)
    REAL(dp), PARAMETER                                :: E0              =  3.6E-2_dp        ! Entrainment coefficient             [unitless]
    REAL(dp), PARAMETER                                :: Cd              =  2.5E-3_dp        ! Drag coefficient                    [unitless]
    REAL(dp), PARAMETER                                :: lambda1         = -0.0573_dp        ! Freezing point salinity coefficient [K]
    REAL(dp), PARAMETER                                :: lambda2         =  0.0832_dp        ! Freezing point offset               [K]
    REAL(dp), PARAMETER                                :: lambda3         =  7.61E-4_dp       ! Freezing point depth coefficient    [K m^-1]
    REAL(dp), PARAMETER                                :: M0              =  10._dp           ! Melt-rate parameter                 [m yr^1 degC^-2]
    REAL(dp), PARAMETER                                :: sqrtCd_GammaTS0 =  6.0E-4_dp        ! Heat exchange parameter             [unitless]
    REAL(dp), PARAMETER                                :: gamma1          =  0.545_dp         ! Heat exchange parameter             [unitless]
    REAL(dp), PARAMETER                                :: gamma2          =  3.5E-5_dp        ! Heat exchange parameter             [m^-1]
    REAL(dp), PARAMETER                                :: x0              =  0.56_dp          ! Empirically derived dimensionless scaling factor
        
    ! Calculate alpha based on the effective basal slope
    alpha    = ATAN( basal_slope)
    sinalpha = SIN( alpha)
    
    ! Calculate freezing temperature at effective grounding-line plume source (Lazeroms et al., 2018, Eq. A7)
    Tf_GL = lambda1 * Sa + lambda2 + lambda3 * plume_source_depth
    
    ! Calculate the effective heat exchange coefficient (Lazeroms et al., 2018, Eq. A8)
    GammaTS = C%BMB_Lazeroms2018_GammaT * (gamma1 + gamma2 * ((Ta - Tf_GL) / lambda3) * ((E0 * sinalpha) / (sqrtCD_GammaTS0 + E0 * sinalpha)))
    sqrtCd_GammaTS = SQRT( Cd) * GammaTS
    
    ! Calculate the geometrical factor in the melt-rate expression (Lazeroms et al., 2018, Eq. A3)
    g_alpha_term1 = SQRT(      sinalpha  / (Cd             + E0 * sinalpha))
    g_alpha_term2 = SQRT( sqrtCd_GammaTS / (sqrtCd_GammaTS + E0 * sinalpha))
    g_alpha_term3 =     ( E0 * sinalpha  / (sqrtCd_GammaTS + E0 * sinalpha))
    g_alpha = g_alpha_term1 * g_alpha_term2 * g_alpha_term3
    
    ! Calculate the empirically derived melt-rate scale (Lazeroms et al., 2018, Eq. A9)
    M = M0 * g_alpha * (Ta - Tf_GL)**2
    
    ! Calculate the geometry- and temperature-dependent length scale (Lazeroms et al., 2018, Eq. A10)
    l_geo = ((Ta - Tf_GL) / lambda3) * (x0 * sqrtCd_GammaTS + E0 * sinalpha) / (x0 * (sqrtCd_GammaTS + E0 * sinalpha))
    
    ! Calculate the dimensionless scaled distance along plume path (Lazeroms et al., 2018, Eq. A11),
    !   and evaluate the dimensionless melt curve at that point
    Xhat = MIN(1._dp, MAX(0._dp, (depth - plume_source_depth) / l_geo ))
    Mhat = Lazeroms2018_dimensionless_melt_curve( Xhat)
    
    ! Finally, calculate the basal melt rate (Lazeroms et al., 2018, Eq. A12)
    melt_rate = -M * Mhat
    
  END SUBROUTINE plume_model_Lazeroms2018
  SUBROUTINE find_effective_plume_path( grid, ice, BMB, stack_GL, n_GL, i,j, eff_plume_source_depth, eff_basal_slope)
    ! Find the effective plume source depth and basal slope for shelf grid cell [i,j]
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: stack_GL
    INTEGER,                             INTENT(IN)    :: n_GL
    INTEGER,                             INTENT(IN)    :: i,j
    REAL(dp),                            INTENT(OUT)   :: eff_plume_source_depth
    REAL(dp),                            INTENT(OUT)   :: eff_basal_slope
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_effective_plume_path'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF     (C%BMB_Lazeroms2018_find_GL_scheme == 'GL_average') THEN
      CALL find_effective_plume_path_GL_average(      grid, ice, BMB,            i,j, eff_plume_source_depth, eff_basal_slope)
    ELSEIF (C%BMB_Lazeroms2018_find_GL_scheme == 'along_ice_flow') THEN
      CALL find_effective_plume_path_along_ice_flow(  grid, ice,                 i,j, eff_plume_source_depth, eff_basal_slope)
    ELSEIF (C%BMB_Lazeroms2018_find_GL_scheme == 'GL_average_Tijn') THEN
      CALL find_effective_plume_path_GL_average_Tijn( grid, ice, stack_GL, n_GL, i,j, eff_plume_source_depth, eff_basal_slope)
    ELSE
      CALL crash('unknown BMB_Lazeroms2018_find_GL_scheme "' // TRIM(C%BMB_Lazeroms2018_find_GL_scheme) // '"!')
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE find_effective_plume_path
  SUBROUTINE find_effective_plume_path_GL_average( grid, ice, BMB, i,j, eff_plume_source_depth, eff_basal_slope)
    ! Find the effective plume source depth and basal slope for shelf grid cell [i,j], following
    ! the approach outlined in Lazeroms et al. (2018), Sect. 2.3
    !
    ! Based loosely on code from IMAU-ICE v1.0 by Heiko Goelzer (2019?)
    !
    ! The straight line between the shelf point and a grounding line point indicates the 1-D path that a meltwater plume from the grounding line may 
    ! have taken to reach the upper point. From this path, we determine the local slope beneath the ice shelf needed for the UPP melt parametrization. 
    ! The average grounding line depth and average slope are a measure of the combined effect of plumes from multiple directions.
    !
    ! ***NOTE*** Only grounding-line points lying DEEPER than the shelf point 
    ! can be a valid physical source of the meltwater plume and, consequently,
    ! valid input for the UPP parametrization. Points lying above depth(i,j)
    ! are ignored. Furthermore, we only consider directions in which the LOCAL
    ! slope at (i,j) is positive.

    ! Based on distance_open_ocean by Bas de Boer:
    ! searches for grounding line points in 16 directions:
    !
    !             12 13 14
    !            11 ---- 15
    !           10 --  -- 16
    !          09 -- ij -- 01
    !           08 --  -- 02
    !            07 ---- 03
    !             06 05 04
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    INTEGER,                             INTENT(IN)    :: i,j
    REAL(dp),                            INTENT(OUT)   :: eff_plume_source_depth
    REAL(dp),                            INTENT(OUT)   :: eff_basal_slope
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_effective_plume_path_GL_average'
    INTEGER                                            :: n
    INTEGER                                            :: dpi,dpj,ip1,jp1,ip2,jp2
    REAL(dp)                                           :: basal_slope
    REAL(dp)                                           :: zb_shelf, dist, zb_dp
    LOGICAL                                            :: reached_end, found_source
    REAL(dp)                                           :: TAF1,TAF2,lambda_GL
    REAL(dp)                                           :: Hb1,Hb2,plume_source_depth
    INTEGER                                            :: n_valid_plumes
    REAL(dp)                                           :: sum_plume_source_depths, sum_basal_slopes
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Initialise
    n_valid_plumes          = 0
    sum_plume_source_depths = 0._dp
    sum_basal_slopes        = 0._dp
    
    ! Calculate the shelf base depth (draft)
    zb_shelf = ice%Hs_a( j,i) - ice%Hi_a( j,i)
    
    ! Investigate all 16 search directions
    DO n = 1, 16
      
      ! The search direction vector
      dpi = BMB%search_directions( n,1)
      dpj = BMB%search_directions( n,2)
      
      ! Initialise the search pointer at the shelf grid cell
      ip1 = i
      jp1 = j
      ip2 = i + dpi
      jp2 = j + dpj
        
      ! If the search direction already points out of the model domain, don't bother
      IF (ip2 < 1 .OR. ip2 > grid%nx .OR. jp2 < 1 .OR. jp2 > grid%ny) CYCLE
      
      ! Calculate the basal slope in this direction (Lazeroms et al. (2018), Eq. 12)
      zb_dp = ice%Hs_a( jp2,ip2) - ice%Hi_a( jp2,ip2)
      dist  = SQRT( REAL(dpi,dp)**2 + REAL(dpj,dp)**2) * grid%dx
      basal_slope = (zb_shelf - zb_dp) / dist
      
      ! If this slope is negative, this plume is not valid
      IF (basal_slope < 0._dp) CYCLE
      
      ! Search in this direction
      reached_end  = .FALSE.
      found_source = .FALSE.
      DO WHILE (.NOT. reached_end)
        
        ! If the pointer exits the model domain, stop the search
        IF (ip2 < 1 .OR. ip2 > grid%nx .OR. jp2 < 1 .OR. jp2 > grid%ny) THEN
          reached_end  = .TRUE.
          found_source = .FALSE.
          EXIT
        END IF
        
        ! If the pointer encounters open ocean, stop the search
        IF (ice%mask_ocean_a( jp2,ip2) == 1 .AND. ice%mask_shelf_a( jp2,ip2) == 0) THEN
          reached_end  = .TRUE.
          found_source = .FALSE.
          EXIT
        END IF
        
        ! If the pointer encounters grounded ice, there must be a grounding line here
        IF (ice%mask_sheet_a( jp2,ip2) == 1) THEN
          
          reached_end  = .TRUE.
          found_source = .TRUE.
          
          ! Interpolate the thickness above flotation to find the sub-grid grounding-line depth ( = plume source depth)
          TAF1 = ice%TAF_a( jp1,ip1)
          TAF2 = ice%TAF_a( jp2,ip2)
          lambda_GL = TAF1 / (TAF1 - TAF2)
          
          Hb1  = ice%Hb_a( jp1,ip1)
          Hb2  = ice%Hb_a( jp2,ip2)
          plume_source_depth = (1._dp - lambda_GL) * Hb1 + lambda_GL * Hb2
          
          ! If this grounding line is less deep than the shelf, don't use it
          IF (plume_source_depth > zb_shelf) THEN
            found_source = .FALSE.
          END IF
          
        END IF ! IF (ice%mask_sheet_a( jp2,ip2) == 1) THEN
        
        ! If none of these exceptions were triggered, advance the pointer in the search direction
        ip1 = ip2
        jp1 = jp2
        ip2 = ip1 + dpi
        jp2 = jp1 + dpj
        
      END DO ! DO WHILE (.NOT. reached_end)
      
      ! If a valid plume source was found, add it to the sum for averaging
      IF (found_source) THEN
        n_valid_plumes          = n_valid_plumes          + 1
        sum_plume_source_depths = sum_plume_source_depths + plume_source_depth
        sum_basal_slopes        = sum_basal_slopes        + basal_slope
      END IF
      
    END DO ! DO n = 1, 16
    
    ! Define the effective plume source depth and basal slope as the average
    ! of those values for all valid plume paths
    
    IF (n_valid_plumes > 0) THEN
      eff_plume_source_depth = sum_plume_source_depths / REAL(n_valid_plumes,dp)
      eff_basal_slope        = sum_basal_slopes        / REAL(n_valid_plumes,dp)
    ELSE
      ! Exception for when no valid plume sources were found
      eff_plume_source_depth = zb_shelf
      eff_basal_slope        = 1E-10_dp   ! Because the melt parameterisation yields NaN for a zero slope
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE find_effective_plume_path_GL_average
  SUBROUTINE find_effective_plume_path_along_ice_flow( grid, ice, i,j, eff_plume_source_depth, eff_basal_slope)
    ! Find the effective plume source depth and basal slope for shelf grid cell [i,j] by
    ! assuming plumes follow the same horizontal flow field as the ice shelf
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    INTEGER,                             INTENT(IN)    :: i,j
    REAL(dp),                            INTENT(OUT)   :: eff_plume_source_depth
    REAL(dp),                            INTENT(OUT)   :: eff_basal_slope
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_effective_plume_path_along_ice_flow'
    REAL(dp), DIMENSION(2)                             :: t1 ,t2 ,uv
    INTEGER                                            :: nit, nit_max
    REAL(dp)                                           :: u1, v1, TAF1, TAF2, depth1, Hi2, Hs2, depth2, lambda, GLx, GLy
    LOGICAL                                            :: found_GL, got_stuck
    REAL(dp)                                           :: zb_shelf, dist
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Calculate the shelf base depth (draft)
    zb_shelf = ice%Hs_a( j,i) - ice%Hi_a( j,i)
    
    ! Track a tracer upstream along the ice flow field until grounded ice is found
    t1     = [grid%x( i), grid%y( j)]
    t2     = t1
    TAF1   = ice%TAF_a( j,i)
    depth1 = zb_shelf
    
    found_GL  = .FALSE.
    got_stuck = .FALSE.
    nit       = 0
    nit_max   = grid%nx + grid%ny
    DO WHILE (.NOT. found_GL)
      
      ! Safety
      nit = nit + 1
      IF (nit > nit_max) THEN
        got_stuck = .TRUE.
        EXIT
      END IF
      
      ! Interpolate velocity field to exact tracer location
      u1 = interp_bilin_2D( ice%u_vav_a, grid%x, grid%y, t1( 1), t1( 2))
      v1 = interp_bilin_2D( ice%v_vav_a, grid%x, grid%y, t1( 1), t1( 2))
      
      ! Move the tracer upstream
      uv = [u1, v1]
      uv = uv / NORM2( uv)
      t2 = t1 - uv * grid%dx
      
      ! Safety
      IF (t2( 1) < grid%xmin .OR. t2( 1) > grid%xmax .OR. t2( 2) < grid%ymin .OR. t2( 2) > grid%ymax) THEN
        got_stuck = .TRUE.
        EXIT
      END IF
      
      ! Interpolate data to new tracer location
      TAF2   = interp_bilin_2D( ice%TAF_a, grid%x, grid%y, t2( 1), t2( 2))
      Hi2    = interp_bilin_2D( ice%Hi_a,  grid%x, grid%y, t2( 1), t2( 2))
      Hs2    = interp_bilin_2D( ice%Hs_a,  grid%x, grid%y, t2( 1), t2( 2))
      depth2 = Hs2 - Hi2
      
      ! Check if we've found grounded ice
      IF (TAF2 > 0._dp) THEN
      
        found_GL = .TRUE.
        
        ! Find exact GL depth
        lambda = TAF1 / (TAF1 - TAF2)
        eff_plume_source_depth = lambda * depth2 + (1._dp - lambda) * depth1
        GLx                    = lambda * t2( 1) + (1._dp - lambda) * t1( 1)
        GLy                    = lambda * t2( 2) + (1._dp - lambda) * t1( 2)
      
        ! Calculate the plume source depth and basal slope
        dist  = SQRT( (GLx - grid%x( i))**2 + (GLy - grid%y( j))**2 )
        eff_basal_slope = (zb_shelf - eff_plume_source_depth) / dist
        
      ELSE
        ! Cycle the tracer
        
        t1     = t2
        TAF1   = TAF2
        depth1 = depth2
        
      END IF ! IF (TAF2 > 0._dp) THEN
      
    END DO ! DO WHILE (.NOT. found_GL)
    
    ! Exception for when the GL source is less deep than the shelf base
    IF (eff_plume_source_depth >= zb_shelf .OR. got_stuck) THEN
      eff_plume_source_depth = zb_shelf
      eff_basal_slope        = 1E-10_dp   ! Because the melt parameterisation yields NaN for a zero slope
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE find_effective_plume_path_along_ice_flow
  SUBROUTINE find_effective_plume_path_GL_average_Tijn( grid, ice, stack_GL, n_GL, i,j, eff_plume_source_depth, eff_basal_slope)
    ! Find the effective plume source depth and basal slope for shelf grid cell [i,j]
    ! 
    ! New approach by Tijn Berends (March 2022): similar to the original approach by Lazeroms et al.,
    ! but instead of averaging over GL points in only the 16 directions, which often leads to sharp discontinuities
    ! in the resulting melt fields, just average over ALL (valid) grounding-line pixels.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: stack_GL
    INTEGER,                             INTENT(IN)    :: n_GL
    INTEGER,                             INTENT(IN)    :: i,j
    REAL(dp),                            INTENT(OUT)   :: eff_plume_source_depth
    REAL(dp),                            INTENT(OUT)   :: eff_basal_slope
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'find_effective_plume_path_GL_average_Tijn'
    INTEGER                                            :: n,ii,jj
    REAL(dp)                                           :: zb_shelf, plume_source_depth, dist, basal_slope, w
    REAL(dp)                                           :: sum_plume_source_depths, sum_basal_slopes, sum_w
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Initialise
    sum_plume_source_depths = 0._dp
    sum_basal_slopes        = 0._dp
    sum_w                   = 0._dp
    
    ! Calculate the shelf base depth (draft)
    zb_shelf = ice%Hs_a( j,i) - ice%Hi_a( j,i)
    
    ! Average over all grounding-line cells
    DO n = 1, n_GL
      
      ! The grounding-line grid cell
      ii = stack_GL( n,1)
      jj = stack_GL( n,2)
      
      ! Calculate the plume source depth and cavity slope
      CALL calc_GL_depth( grid, ice, ii, jj, plume_source_depth)
      !plume_source_depth = ice%Hb_a( jj,ii)
      dist = MAX( grid%dx, SQRT( REAL( ii - i, dp)**2 + REAL( jj - j, dp)**2) * grid%dx)
      basal_slope = (zb_shelf - plume_source_depth) / dist
      
      ! If this slope is negative, this plume is not valid
      IF (basal_slope < 0._dp) CYCLE
      
      ! Weight factor based on inverse of square of distance and plume source depth
      w = -plume_source_depth / dist**2
      
      ! Add to the sum
      sum_w                   = sum_w                   + w
      sum_plume_source_depths = sum_plume_source_depths + w * plume_source_depth
      sum_basal_slopes        = sum_basal_slopes        + w * basal_slope
      
    END DO ! DO n = 1, 16
    
    ! Define the effective plume source depth and basal slope as the average
    ! of those values for all valid plume paths
    
    IF (sum_w > 0._dp) THEN
      eff_plume_source_depth = sum_plume_source_depths / sum_w
      eff_basal_slope        = sum_basal_slopes        / sum_w
    ELSE
      ! Exception for when no valid plume sources were found
      eff_plume_source_depth = zb_shelf
      eff_basal_slope        = 1E-10_dp   ! Because the melt parameterisation yields NaN for a zero slope
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE find_effective_plume_path_GL_average_Tijn
  FUNCTION Lazeroms2018_dimensionless_melt_curve( xhat) RESULT( Mhat)
    ! The dimensionless melt curve from Lazeroms et al. (2018), Appendix A, Eq. A13
    
    IMPLICIT NONE
  
    REAL(dp), INTENT(IN)                                 :: xhat                        ! Scaled distance  from grounding line
    REAL(dp)                                             :: Mhat                        ! Scaled melt rate from polynomial fit
    
    ! L variables: polynomial coefficients
    REAL(dp), PARAMETER                                  :: p11 =  6.387953795485420E4_dp
    REAL(dp), PARAMETER                                  :: p10 = -3.520598035764990E5_dp
    REAL(dp), PARAMETER                                  :: p9  =  8.466870335320488E5_dp
    REAL(dp), PARAMETER                                  :: p8  = -1.166290429178556E6_dp
    REAL(dp), PARAMETER                                  :: p7  =  1.015475347943186E6_dp
    REAL(dp), PARAMETER                                  :: p6  = -5.820015295669482E5_dp
    REAL(dp), PARAMETER                                  :: p5  =  2.218596970948727E5_dp
    REAL(dp), PARAMETER                                  :: p4  = -5.563863123811898E4_dp
    REAL(dp), PARAMETER                                  :: p3  =  8.927093637594877E3_dp
    REAL(dp), PARAMETER                                  :: p2  = -8.951812433987858E2_dp
    REAL(dp), PARAMETER                                  :: p1  =  5.527656234709359E1_dp
    REAL(dp), PARAMETER                                  :: p0  =  0.1371330075095435_dp
    
    Mhat = p11 * xhat**11 + &
           p10 * xhat**10 + &
           p9  * xhat**9  + &
           p8  * xhat**8  + &
           p7  * xhat**7  + &
           p6  * xhat**6  + &
           p5  * xhat**5  + &
           p4  * xhat**4  + &
           p3  * xhat**3  + &
           p2  * xhat**2  + &
           p1  * xhat + p0
    
  END FUNCTION Lazeroms2018_dimensionless_melt_curve
  SUBROUTINE initialise_BMB_model_Lazeroms2018_plume( grid, BMB)
    ! Allocate memory for the data fields of the Favier et al. (2019) shelf BMB parameterisations.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model_Lazeroms2018_plume'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Variables
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx, BMB%T_ocean_base,           BMB%wT_ocean_base          )
    CALL allocate_shared_int_2D( 16,      2,       BMB%search_directions,      BMB%wsearch_directions     )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx, BMB%eff_plume_source_depth, BMB%weff_plume_source_depth)
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx, BMB%eff_basal_slope,        BMB%weff_basal_slope       )
    
    ! Define the 16 search directions
    IF (par%master) THEN
      BMB%search_directions(  1,:) = (/  1,  0 /)
      BMB%search_directions(  2,:) = (/  2, -1 /)
      BMB%search_directions(  3,:) = (/  1, -1 /)
      BMB%search_directions(  4,:) = (/  1, -2 /)
      BMB%search_directions(  5,:) = (/  0, -1 /)
      BMB%search_directions(  6,:) = (/ -1, -2 /)
      BMB%search_directions(  7,:) = (/ -1, -1 /)
      BMB%search_directions(  8,:) = (/ -2, -1 /)
      BMB%search_directions(  9,:) = (/ -1,  0 /)
      BMB%search_directions( 10,:) = (/ -2,  1 /)
      BMB%search_directions( 11,:) = (/ -1,  1 /)
      BMB%search_directions( 12,:) = (/ -1,  2 /)
      BMB%search_directions( 13,:) = (/  0,  1 /)
      BMB%search_directions( 14,:) = (/  1,  2 /)
      BMB%search_directions( 15,:) = (/  1,  1 /)
      BMB%search_directions( 16,:) = (/  2,  1 /)
    END IF
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
      
  END SUBROUTINE initialise_BMB_model_Lazeroms2018_plume
  
! == The PICO ocean box model
! ===========================

  SUBROUTINE run_BMB_model_PICO( grid, ice, ocean, BMB)
    ! Calculate basal melt using the PICO ocean box model
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_PICO'
    INTEGER                                            :: basin_i
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Assign PICO ocean boxes to all shelf grid cells in all basins
    CALL PICO_assign_ocean_boxes( grid, ice, BMB)
    
    ! Run PICO for all basins
    DO basin_i = 1, ice%nbasins
      CALL run_BMB_model_PICO_basin( grid, ice, ocean, BMB, basin_i)
    END DO
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE run_BMB_model_PICO
  SUBROUTINE PICO_assign_ocean_boxes( grid, ice, BMB)
    ! Assign PICO ocean boxes to shelf grid cells using the distance-to-grounding-line / distance-to-calving-front
    ! approach outlined in Reese et al. (2018)
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'PICO_assign_ocean_boxes'
    INTEGER                                            :: i,j,k,basin_i
    REAL(dp), DIMENSION(:    ), POINTER                ::  d_GL_D
    INTEGER                                            :: wd_GL_D
    REAL(dp)                                           :: d_max
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: n_cells_per_box
    LOGICAL                                            :: do_reduce_n_D
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! For each shelf grid cell, calculate the shortest linear distance to the grounding line d_GL,
    ! and to the ice front d_IF, and the scaled flowline distance r (Reese et al. (2018), Eq. 10)
    CALL PICO_calc_r( grid, ice, BMB%PICO_d_GL, BMB%PICO_d_IF, BMB%PICO_r)
    
    ! Calculate maximum distance to grounding line within each basin
    CALL allocate_shared_dp_1D(  ice%nbasins, d_GL_D, wd_GL_D)
    DO basin_i = 1, ice%nbasins
      d_max = 0._dp
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        IF (ice%basin_ID( j,i) == basin_i) THEN
          d_max = MAX( d_max, BMB%PICO_d_GL( j,i))
        END IF
      END DO
      END DO
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, d_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      IF (par%master) d_GL_D( basin_i) = d_max
      CALL sync
    END DO ! DO basin_i = 1, ice%nbasins
    
    ! Calculate maximum distance to grounding line within the entire model doman
    IF (par%master) d_max = MAXVAL( d_GL_D)
    CALL MPI_BCAST( d_max, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    
    ! Determine number of PICO boxes for each basin
    IF (par%master) THEN
      DO basin_i = 1, ice%nbasins
        ! Reese et al. (2018), Eq. 9
        BMB%PICO_n_D( basin_i) = 1 + INT( SQRT( d_GL_D( basin_i) / d_max) * REAL( C%BMB_PICO_nboxes - 1,dp))
      END DO
    END IF
    CALL sync
    
  ! Assign PICO boxes to all shelf grid cells in all basins
  ! =======================================================
    
    ALLOCATE( n_cells_per_box( C%BMB_PICO_nboxes))
    
    DO basin_i = 1, ice%nbasins
    
      ! If necessary, reduce the maximum number of boxes for a basin
      ! until every box is assigned at least one grid cell
      
      n_cells_per_box = 0
      
      DO WHILE (.TRUE.)
        
        ! Assign ocean boxes according to Reese et al. (2018), Eq. 11
        DO i = grid%i1, grid%i2
        DO j = 1, grid%ny
          IF (ice%basin_ID( j,i) == basin_i) THEN
            BMB%PICO_k( j,i) = 0
            IF (ice%mask_shelf_a( j,i) == 1) THEN
              DO k = 1, BMB%PICO_n_D( basin_i)
                  IF (1._dp - SQRT( REAL(BMB%PICO_n_D( basin_i) - k + 1, dp) / REAL( BMB%PICO_n_D( basin_i), dp)) <= BMB%PICO_r( j,i) .AND. &
                      1._dp - SQRT( REAL(BMB%PICO_n_D( basin_i) - k    , dp) / REAL( BMB%PICO_n_D( basin_i), dp)) >= BMB%PICO_r( j,i)) THEN
                  BMB%PICO_k( j,i) = k
                  n_cells_per_box( k) = n_cells_per_box( k) + 1
                END IF
              END DO ! DO k = 1, BMB%PICO_n_D( basin_i)
            END IF ! IF (ice%mask_shelf_a( j,i) == 1) THEN
          END IF ! IF (ice%basin_ID( j,i) == basin_i) THEN
        END DO
        END DO
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, n_cells_per_box, C%BMB_PICO_nboxes, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        IF (par%master) BMB%PICO_A( basin_i,:) = n_cells_per_box * grid%dx**2
        CALL sync
        
        ! If any of the boxes have zero grid cells assigned, reduce the maximum number of boxes and try again
        do_reduce_n_D = .FALSE.
        DO k = 1, BMB%PICO_n_D( basin_i)
          IF (n_cells_per_box( k) == 0) THEN
            do_reduce_n_D = .TRUE.
          END IF
        END DO
        IF (do_reduce_n_D) THEN
          IF (par%master) BMB%PICO_n_D( basin_i) = BMB%PICO_n_D( basin_i) - 1
          CALL sync
        ELSE
          EXIT
        END IF
        
      END DO ! DO WHILE (.TRUE.)
      
    END DO ! DO basin_i = 1, ice%nbasins
    
    ! Clean up after yourself
    DEALLOCATE( n_cells_per_box)
    CALL deallocate_shared( wd_GL_D)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE PICO_assign_ocean_boxes
  SUBROUTINE run_BMB_model_PICO_basin( grid, ice, ocean, BMB, basin_i)
    ! Calculate basal melt for ice basin i using the PICO ocean box model
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    INTEGER,                             INTENT(IN)    :: basin_i
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_PICO_basin'
    INTEGER                                            :: i,j,k
    REAL(dp)                                           :: Tk0,Sk0
    REAL(dp)                                           :: nu,lambda
    REAL(dp)                                           :: g1,g2,s,Crbsa,Tstar,x,y
    REAL(dp)                                           :: q
    
    ! Constants (Reese et al. (2018), Table 1)
    REAL(dp), PARAMETER                                :: aa          = -0.0572_dp     ! Salinity coefficient of freezing equation         [degC PSU^-1]
    REAL(dp), PARAMETER                                :: bb          =  0.0788_dp     ! Constant coefficient of freezing equation         [degC]
    REAL(dp), PARAMETER                                :: cc          = 7.77E-8_dp     ! Pressure coefficient of freezing equation         [degC Pa^-1]
    REAL(dp), PARAMETER                                :: alpha       =  7.5E-5_dp     ! Thermal expansion coefficient in EOS              [degC^-1]
    REAL(dp), PARAMETER                                :: beta        =  7.7E-4_dp     ! Salt contraction coefficient in EOS               [PSU^-1]
    REAL(dp), PARAMETER                                :: rhostar     = 1033_dp        ! Reference density in EOS                          [kg m^-3]
    REAL(dp), PARAMETER                                :: C_overturn  = 1.0E6_dp       ! Overturning strength                              [m^6 s^-1 kg^-1]
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Initialise
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%basin_ID( j,i) == basin_i) THEN
        BMB%PICO_T( j,i) = 0._dp
        BMB%PICO_S( j,i) = 0._dp
        BMB%PICO_m( j,i) = 0._dp
      END IF
    END DO
    END DO
    CALL sync
    
    ! Some intermediary constants (Reese et al. (2018), just after Eq. A2)
    nu     = ice_density / seawater_density
    lambda = L_fusion / cp_ocean
    
    ! Calculate temperature and salinity in box B0 for this basin
    CALL PICO_calc_T0_S0( grid, ice, ocean, basin_i, Tk0, Sk0)
    
    ! Calculate 2-D + box-averaged basal pressures
    BMB%PICO_p( :,grid%i1:grid%i2) = ice_density * grav * ice%Hi_a( :,grid%i1:grid%i2)
    DO k = 1, BMB%PICO_n_D( basin_i)
      CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_p, basin_i, k, BMB%PICO_pk( basin_i, k))
    END DO
    
  ! Calculate solution for box 1
  ! ============================
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      IF (ice%basin_ID( j,i) == basin_i .AND. BMB%PICO_k( j,i) == 1) THEN
        
        ! Reese et al. (2018), just before Eq. A6
        g1 = BMB%PICO_A( basin_i, 1) * C%BMB_PICO_GammaTstar
        g2 = g1 / (nu * lambda)
        Tstar = aa * Sk0 + bb - cc * BMB%PICO_p( j,i) - Tk0
        
        ! Reese et al. (2018), just after Eq. A11
        s = Sk0 / (nu * lambda)
        
        ! Intermediary constants
        Crbsa = C_overturn * rhostar * (beta * s - alpha)
        
        ! Reese et al. (2018), Eq. A12
        x = -g1 / (2._dp * Crbsa) + SQRT( (g1 / (2._dp * Crbsa))**2 - (g1 * Tstar / Crbsa))
        
        ! Reese et al. (2018), Eq. A8
        y = Sk0 * x / (nu * lambda)
        
        BMB%PICO_T( j,i) = Tk0 - x
        BMB%PICO_S( j,i) = Sk0 - y
        
        ! Reese et al. (2019), Eq. 13
        BMB%PICO_m( j,i) = sec_per_year * C%BMB_PICO_GammaTstar / (nu*lambda) * (aa * BMB%PICO_S( j,i) + bb - cc * BMB%PICO_p( j,i) - BMB%PICO_T( j,i))
        
      END IF ! IF (BMB%PICO_k( j,i) == 1) THEN
      
    END DO
    END DO
    CALL sync
    
    ! Calculate box-averaged values
    CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_T, basin_i, 1, BMB%PICO_Tk( basin_i, 1))
    CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_S, basin_i, 1, BMB%PICO_Sk( basin_i, 1))
    CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_m, basin_i, 1, BMB%PICO_mk( basin_i, 1))
    
    ! Calculate overturning strength (Reese et al. (2018), Eq. A9)
    q = C_overturn * rhostar * (beta * (Sk0 - BMB%PICO_Sk( basin_i, 1)) - alpha * (Tk0 - BMB%PICO_Tk( basin_i, 1))) 
    
  ! Calculate solutions for subsequent boxes
  ! ========================================
    
    DO k = 2, BMB%PICO_n_D( basin_i)
      
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        
        IF (ice%basin_ID( j,i) == basin_i .AND. BMB%PICO_k( j,i) == k) THEN
        
          ! Reese et al. (2018), just before Eq. A6
          g1 = BMB%PICO_A( basin_i, k) * C%BMB_PICO_GammaTstar
          g2 = g1 / (nu * lambda)
          Tstar = aa * BMB%PICO_Sk( basin_i, k-1) + bb - cc * BMB%PICO_p( j,i) - BMB%PICO_Tk( basin_i, k-1)
          
          ! Reese et al. (2018), Eq. A13
          x = -g1 * Tstar / (q + g1 - g2 * aa * BMB%PICO_Sk( basin_i, k-1))
        
          ! Reese et al. (2018), Eq. A8
          y = BMB%PICO_Sk( basin_i, k-1) * x / (nu * lambda)
        
          BMB%PICO_T( j,i) = BMB%PICO_Tk( basin_i, k-1) - x
          BMB%PICO_S( j,i) = BMB%PICO_Sk( basin_i, k-1) - y
        
          ! Reese et al. (2019), Eq. 13
          BMB%PICO_m( j,i) = sec_per_year * C%BMB_PICO_GammaTstar / (nu*lambda) * (aa * BMB%PICO_S( j,i) + bb - cc * BMB%PICO_p( j,i) - BMB%PICO_T( j,i))
          
        END IF ! IF (BMB%PICO_k( j,i) == k) THEN
        
      END DO
      END DO
      CALL sync
    
      ! Calculate box-averaged values
      CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_T, basin_i, k, BMB%PICO_Tk( basin_i, k))
      CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_S, basin_i, k, BMB%PICO_Sk( basin_i, k))
      CALL PICO_calc_box_average( grid, ice, BMB, BMB%PICO_m, basin_i, k, BMB%PICO_mk( basin_i, k))
      
    END DO ! DO k = 2, BMB%PICO_n_D( basin_i)
    
    ! Copy melt rates to final data field      
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%basin_ID( j,i) == basin_i) BMB%BMB_shelf( j,i) = BMB%PICO_m( j,i)
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE run_BMB_model_PICO_basin
  SUBROUTINE PICO_calc_r( grid, ice, d_GL, d_IF, r)
    ! For each shelf grid cell, calculate the shortest linear distance to the grounding line d_GL,
    ! and to the ice front d_IF, and the scaled flowline distance r (Reese et al. (2018), Eq. 10)
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: d_GL, d_IF, r
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'PICO_calc_r'
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  mask_shelf
    INTEGER                                            :: wmask_shelf
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  list_GL,  list_IF
    INTEGER                                            :: wlist_GL, wlist_IF
    INTEGER,                    POINTER                ::  n_GL,  n_IF
    INTEGER                                            :: wn_GL, wn_IF
    INTEGER                                            :: basin_i,i,j
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx,   mask_shelf, wmask_shelf)
    CALL allocate_shared_int_2D( grid%nx*grid%ny, 2, list_GL   , wlist_GL   )
    CALL allocate_shared_int_2D( grid%nx*grid%ny, 2, list_IF   , wlist_IF   )
    CALL allocate_shared_int_0D(                     n_GL      , wn_GL      )
    CALL allocate_shared_int_0D(                     n_IF      , wn_IF      )
    
    ! Make temporary copy of the shelf mask
    mask_shelf( :,grid%i1:grid%i2) = ice%mask_shelf_a( :,grid%i1:grid%i2)
    
    ! Fill "holes" in the shelf mask (small open ocean pockets surrounded by shelf)
    CALL fill_shelf_holes( grid, ice%mask_ocean_a, mask_shelf)
    
    ! Treat every basin separately
    DO basin_i = 1, ice%nbasins
      
      ! List all grounding-line and ice-front grid cells
      CALL list_GL_grid_cells_basin( grid, ice, ice%mask_sheet_a,     mask_shelf  , basin_i, list_GL, n_GL)
      CALL list_IF_grid_cells_basin( grid, ice,     mask_shelf  , ice%mask_ocean_a, basin_i, list_IF, n_IF)
      
      ! Calculate shortest linear distance to the grounding line and to
      ! the ice front for all shelf grid cells within the specified basin
      
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        
        IF (ice%mask_shelf_a( j,i) == 1 .AND. ice%basin_ID( j,i) == basin_i) THEN
          
          ! Calculate shortest linear distance to the grounding line and the ice front
          CALL calc_shortest_distance_to_GL( grid, list_GL, n_GL, i,j, d_GL( j,i))
          CALL calc_shortest_distance_to_IF( grid, list_IF, n_IF, i,j, d_IF( j,i))
          
          ! Regularisation (to prevent divide-by-zero errors)
          d_GL( j,i) = MAX( d_GL( j,i), grid%dx / 2._dp)
          d_IF( j,i) = MAX( d_IF( j,i), grid%dx / 2._dp)
    
          ! Calculate the scaled flowline distance r (Reese et al. (2018), Eq. 10)
          r( j,i) = d_GL( j,i) / (d_GL( j,i) + d_IF( j,i))
          
        ELSEIF (ice%mask_shelf_a( j,i) == 0 .AND. ice%mask_ocean_a( j,i) == 1 .AND. ice%basin_ID( j,i) == basin_i) THEN
          
          r( j,i) = 1._dp
          
        END IF ! IF (ice%mask_shelf_a( j,i) == 1 .AND. ice%basin_ID( j,i) == basin_i) THEN
        
      END DO
      END DO
      CALL sync
      
    END DO ! DO basin_i = 1, ice%nbasins
    
    ! Clean up after yourself
    CALL deallocate_shared( wmask_shelf)
    CALL deallocate_shared( wlist_GL   )
    CALL deallocate_shared( wlist_IF   )
    CALL deallocate_shared( wn_GL      )
    CALL deallocate_shared( wn_IF      )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE PICO_calc_r
  SUBROUTINE PICO_calc_T0_S0( grid, ice, ocean, basin_i, Tk0, Sk0)
    ! Find temperature and salinity in box B0 (defined as mean ocean-floor value at the calving front)
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    INTEGER,                             INTENT(IN)    :: basin_i
    REAL(dp),                            INTENT(OUT)   :: Tk0,Sk0
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'PICO_calc_T0_S0'
    INTEGER                                            :: i,j,n
    REAL(dp)                                           :: depth, T_floor, S_floor, depth_max
    INTEGER                                            :: ii,jj
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Average ocean-floor temperature and salinity over this basin's ocean-next-to-floating-ice pixels
    n   = 0
    Tk0 = 0._dp
    Sk0 = 0._dp
    
    DO i = MAX(2,grid%i1), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
      
      IF (ice%basin_ID( j,i) == basin_i .AND. ice%mask_ocean_a( j,i) == 1 .AND. ice%mask_ice_a( j,i) == 0) THEN
        IF (ice%mask_shelf_a( j-1,i-1) == 1 .OR. &
            ice%mask_shelf_a( j-1,i  ) == 1 .OR. &
            ice%mask_shelf_a( j-1,i+1) == 1 .OR. &
            ice%mask_shelf_a( j  ,i-1) == 1 .OR. &
            ice%mask_shelf_a( j  ,i+1) == 1 .OR. &
            ice%mask_shelf_a( j+1,i-1) == 1 .OR. &
            ice%mask_shelf_a( j+1,i  ) == 1 .OR. &
            ice%mask_shelf_a( j+1,i+1) == 1) THEN
          ! This pixel is open ocean next to floating ice
          
          ! Find ocean-floor temperature and salinity
          depth = MAX( 0.1_dp, -ice%Hb_a( j,i))
          CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( :,j,i), depth, T_floor)
          CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S_ocean_corr_ext( :,j,i), depth, S_floor)
          
          ! Add to sum
          n   = n   + 1
          Tk0 = Tk0 + T_floor
          Sk0 = Sk0 + S_floor
          
        END IF
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Combine results from processes
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n,   1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, Tk0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, Sk0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    Tk0 = Tk0 / REAL(n,dp)
    Sk0 = Sk0 / REAL(n,dp)
    
    ! Safety
    IF (n == 0) THEN
      ! No ocean-next-to-shelf grid cells could be found within this basin;
      ! instead, just take the value from the deepest ocean floor grid cell (which must by definition be ice-covered...)
      
      IF (par%master) THEN
      
        depth_max = 0._dp
        ii = 0
        jj = 0
        
        DO i = 1, grid%nx
        DO j = 1, grid%ny
          IF (ice%mask_ocean_a( j,i) == 1) THEN
            depth = -ice%Hb_a( j,i)
            IF (depth > depth_max) THEN
              depth_max = depth
              ii = i
              jj = j
            END IF
          END IF
        END DO
        END DO
        
        IF (ii == 0 .OR. jj == 0) THEN
          CALL crash('couldnt find deepest ocean floor grid cell!')
        END IF
        
        ! Find ocean-floor temperature and salinity
        CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( :,jj,ii), depth_max, Tk0)
        CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S_ocean_corr_ext( :,jj,ii), depth_max, Sk0)
        
      END IF ! IF (par%master) THEN
      
      CALL MPI_BCAST( Tk0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST( Sk0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      
    END IF ! IF (n == 0) THEN
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE PICO_calc_T0_S0
  SUBROUTINE PICO_calc_box_average( grid, ice, BMB, d, basin_i, k, d_av)
    ! Calculate the average d_av of field d over ocean box k in basin i
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d
    INTEGER,                             INTENT(IN)    :: basin_i
    INTEGER,                             INTENT(IN)    :: k
    REAL(dp),                            INTENT(OUT)   :: d_av
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'PICO_calc_box_average'
    INTEGER                                            :: i,j,n
    REAL(dp)                                           :: d_sum
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    n     = 0
    d_sum = 0._dp
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%basin_ID( j,i) == basin_i .AND. BMB%PICO_k( j,i) == k) THEN
        n     = n     + 1
        d_sum = d_sum + d( j,i)
      END IF
    END DO
    END DO
    
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, n,     1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, d_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    d_av = d_sum / REAL(n,dp)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE PICO_calc_box_average
  SUBROUTINE initialise_BMB_model_PICO( grid, ice, BMB)
    ! Allocate memory for the data fields of the PICO ocean box model
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model_PICO'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Variables
    CALL allocate_shared_int_2D( 16,      2,        BMB%search_directions,   BMB%wsearch_directions  )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_d_GL,           BMB%wPICO_d_GL          )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_d_IF,           BMB%wPICO_d_IF          )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_r,              BMB%wPICO_r             )
    CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_A, BMB%wPICO_A             )
    CALL allocate_shared_int_1D( ice%nbasins,       BMB%PICO_n_D,            BMB%wPICO_n_D           )
    CALL allocate_shared_int_2D( grid%ny, grid%nx,  BMB%PICO_k,              BMB%wPICO_k             )
    
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_T,              BMB%wPICO_T             )
    CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_Tk, BMB%wPICO_Tk           )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_S,              BMB%wPICO_S             )
    CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_Sk, BMB%wPICO_Sk           )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_p,              BMB%wPICO_p             )
    CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_pk, BMB%wPICO_pk           )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_m,              BMB%wPICO_m             )
    CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_mk, BMB%wPICO_mk           )
    
    ! Define the 16 search directions
    IF (par%master) THEN
      BMB%search_directions(  1,:) = (/  1,  0 /)
      BMB%search_directions(  2,:) = (/  2, -1 /)
      BMB%search_directions(  3,:) = (/  1, -1 /)
      BMB%search_directions(  4,:) = (/  1, -2 /)
      BMB%search_directions(  5,:) = (/  0, -1 /)
      BMB%search_directions(  6,:) = (/ -1, -2 /)
      BMB%search_directions(  7,:) = (/ -1, -1 /)
      BMB%search_directions(  8,:) = (/ -2, -1 /)
      BMB%search_directions(  9,:) = (/ -1,  0 /)
      BMB%search_directions( 10,:) = (/ -2,  1 /)
      BMB%search_directions( 11,:) = (/ -1,  1 /)
      BMB%search_directions( 12,:) = (/ -1,  2 /)
      BMB%search_directions( 13,:) = (/  0,  1 /)
      BMB%search_directions( 14,:) = (/  1,  2 /)
      BMB%search_directions( 15,:) = (/  1,  1 /)
      BMB%search_directions( 16,:) = (/  2,  1 /)
    END IF
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
      
  END SUBROUTINE initialise_BMB_model_PICO
  
! == The PICOP ocean box + plume model
! ====================================

  SUBROUTINE run_BMB_model_PICOP( grid, ice, ocean, BMB)
    ! Calculate basal melt using the PICOP ocean box + plume model
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_ocean_snapshot_regional),  INTENT(IN)    :: ocean
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_BMB_model_PICOP'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! First run the PICO ocean box model to determine the temperature and salinity in the cavity
    CALL run_BMB_model_PICO( grid, ice, ocean, BMB)
    
    ! Then run the Lazeroms (2018) plume parameterisation to calculate melt rates
    CALL run_BMB_model_Lazeroms2018_plume( grid, ice, ocean, BMB)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE run_BMB_model_PICOP
  SUBROUTINE initialise_BMB_model_PICOP( grid, ice, BMB)
    ! Allocate memory for the data fields of the PICOP ocean box + plume model
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_BMB_model_PICOP'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Variables
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%T_ocean_base,           BMB%wT_ocean_base          )
    CALL allocate_shared_int_2D( 16,      2,        BMB%search_directions,      BMB%wsearch_directions     )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%eff_plume_source_depth, BMB%weff_plume_source_depth)
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%eff_basal_slope,        BMB%weff_basal_slope       )
    
    CALL allocate_shared_int_2D( 16,      2,        BMB%search_directions,   BMB%wsearch_directions  )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_d_GL,           BMB%wPICO_d_GL          )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_d_IF,           BMB%wPICO_d_IF          )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_r,              BMB%wPICO_r             )
    CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_A, BMB%wPICO_A             )
    CALL allocate_shared_int_1D( ice%nbasins,       BMB%PICO_n_D,            BMB%wPICO_n_D           )
    CALL allocate_shared_int_2D( grid%ny, grid%nx,  BMB%PICO_k,              BMB%wPICO_k             )
    
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_T,              BMB%wPICO_T             )
    CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_Tk, BMB%wPICO_Tk           )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_S,              BMB%wPICO_S             )
    CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_Sk, BMB%wPICO_Sk           )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_p,              BMB%wPICO_p             )
    CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_pk, BMB%wPICO_pk           )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx,  BMB%PICO_m,              BMB%wPICO_m             )
    CALL allocate_shared_dp_2D(  ice%nbasins, C%BMB_PICO_nboxes, BMB%PICO_mk, BMB%wPICO_mk           )
    
    ! Define the 16 search directions
    IF (par%master) THEN
      BMB%search_directions(  1,:) = (/  1,  0 /)
      BMB%search_directions(  2,:) = (/  2, -1 /)
      BMB%search_directions(  3,:) = (/  1, -1 /)
      BMB%search_directions(  4,:) = (/  1, -2 /)
      BMB%search_directions(  5,:) = (/  0, -1 /)
      BMB%search_directions(  6,:) = (/ -1, -2 /)
      BMB%search_directions(  7,:) = (/ -1, -1 /)
      BMB%search_directions(  8,:) = (/ -2, -1 /)
      BMB%search_directions(  9,:) = (/ -1,  0 /)
      BMB%search_directions( 10,:) = (/ -2,  1 /)
      BMB%search_directions( 11,:) = (/ -1,  1 /)
      BMB%search_directions( 12,:) = (/ -1,  2 /)
      BMB%search_directions( 13,:) = (/  0,  1 /)
      BMB%search_directions( 14,:) = (/  1,  2 /)
      BMB%search_directions( 15,:) = (/  1,  1 /)
      BMB%search_directions( 16,:) = (/  2,  1 /)
    END IF
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
      
  END SUBROUTINE initialise_BMB_model_PICOP
  
! == Some generally useful tools
! ==============================
  
  ! Routine for extrapolating melt field from the regular (FCMP) mask
  ! to the (extended) PMP mask
  SUBROUTINE extrapolate_melt_from_FCMP_to_PMP( grid, ice, BMB)
    ! All the BMB parameterisations are implicitly run using the FCMP sub-grid scheme
    ! (i.e. they are only applied to grid cells whose centre is floating).
    ! Calculating melt rates for partially-floating-but-grounded-at-the-centre cells (which
    ! is needed in the PMP scheme) is not straightforward; instead, just extrapolate values into
    ! the grid cells where this is the case.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extrapolate_melt_from_FCMP_to_PMP'
    INTEGER                                            :: i,j,i1,i2,j1,j2,ii,jj,n,n_ext
    INTEGER,  DIMENSION(:,:  ), POINTER                :: mask_FCMP, mask_PMP
    REAL(dp), DIMENSION(:,:  ), POINTER                :: BMB_shelf_extra
    REAL(dp), DIMENSION(:,:  ), POINTER                :: dBMBdx, dBMBdy
    INTEGER                                            :: wmask_FCMP, wmask_PMP, wBMB_shelf_extra, wdBMBdx, wdBMBdy
    REAL(dp)                                           :: BMB_av
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx, mask_FCMP,       wmask_FCMP      )
    CALL allocate_shared_int_2D( grid%ny, grid%nx, mask_PMP,        wmask_PMP       )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx, BMB_shelf_extra, wBMB_shelf_extra)
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx, dBMBdx,          wdBMBdx         )
    CALL allocate_shared_dp_2D(  grid%ny, grid%nx, dBMBdy,          wdBMBdy         )
    
    ! Define the two masks
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      mask_FCMP( j,i) = 0
      mask_PMP(  j,i) = 0
      IF (ice%mask_shelf_a( j,i) == 1) mask_FCMP( j,i) = 1
      IF (ice%f_grnd_a( j,i) < 1._dp .AND. ice%mask_ice_a( j,i) == 1) mask_PMP(  j,i) = 1
    END DO
    END DO
    CALL sync
    
    ! Calculate spatial derivatives of the melt field, accounting for masked grid cells
    DO i = MAX(2,grid%i2), MIN(grid%nx-1,grid%i2)
    DO j = 2, grid%ny-1
    
      ! d/dx
      IF (mask_FCMP( j,i) == 0) THEN
        dBMBdx( j,i) = 0._dp
      ELSE
        IF     (mask_FCMP( j,i-1) == 1 .AND. mask_FCMP( j,i+1) == 1) THEN
          ! Central differencing
          dBMBdx( j,i) = (BMB%BMB_shelf( j,i+1) - BMB%BMB_shelf( j,i-1)) / (2._dp * grid%dx)
        ELSEIF (mask_FCMP( j,i-1) == 0 .AND. mask_FCMP( j,i+1) == 1) THEN
          ! One-sided differencing
          dBMBdx( j,i) = (BMB%BMB_shelf( j,i+1) - BMB%BMB_shelf( j,i  )) / grid%dx
        ELSEIF (mask_FCMP( j,i-1) == 1 .AND. mask_FCMP( j,i+1) == 0) THEN
          ! One-sided differencing
          dBMBdx( j,i) = (BMB%BMB_shelf( j,i  ) - BMB%BMB_shelf( j,i-1)) / grid%dx
        ELSE
          dBMBdx( j,i) = 0._dp
        END IF
      END IF
      
      ! d/dy
      IF (mask_FCMP( j,i) == 0) THEN
        dBMBdy( j,i) = 0._dp
      ELSE
        IF     (mask_FCMP( j-1,i) == 1 .AND. mask_FCMP( j+1,i) == 1) THEN
          ! Central differencing
          dBMBdy( j,i) = (BMB%BMB_shelf( j+1,i) - BMB%BMB_shelf( j-1,i)) / (2._dp * grid%dx)
        ELSEIF (mask_FCMP( j-1,i) == 0 .AND. mask_FCMP( j+1,i) == 1) THEN
          ! One-sided differencing
          dBMBdy( j,i) = (BMB%BMB_shelf( j+1,i) - BMB%BMB_shelf( j  ,i)) / grid%dx
        ELSEIF (mask_FCMP( j-1,i) == 1 .AND. mask_FCMP( j+1,i) == 0) THEN
          ! One-sided differencing
          dBMBdy( j,i) = (BMB%BMB_shelf( j  ,i) - BMB%BMB_shelf( j-1,i)) / grid%dx
        ELSE
          dBMBdy( j,i) = 0._dp
        END IF
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Extrapolate melt from the FCMP to the PMP mask by taking the average over all FCMP neighbours
    n_ext = 1
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Initialise
      BMB_shelf_extra( j,i) = 0._dp
      
      IF (mask_PMP( j,i) == 1) THEN
        IF (mask_FCMP( j,i) == 1) THEN
          ! Simply copy the FCMP value
          
          BMB_shelf_extra( j,i) = BMB%BMB_shelf( j,i)
          
        ELSE
          ! Calculate average melt rate over all FCMP neighbours
          
          n      = 0
          BMB_av = 0._dp
          n_ext  = 0
          
          DO WHILE (n == 0)
            
            n_ext = n_ext + 1
            
            i1 = MAX( 1      , i - n_ext)
            i2 = MIN( grid%nx, i + n_ext)
            j1 = MAX( 1      , j - n_ext)
            j2 = MIN( grid%ny, j + n_ext)
            
            DO ii = i1, i2
            DO jj = j1, j2
              IF (mask_FCMP( jj,ii) == 1) THEN
                n = n + 1
                BMB_av = BMB_av + BMB%BMB_shelf( jj,ii) + REAL(ii-i,dp) * grid%dx * dBMBdx( j,i) &
                                                        + REAL(jj-j,dp) * grid%dx * dBMBdy( j,i)
              END IF
            END DO
            END DO
            
          END DO ! DO WHILE (n == 0)
          
          BMB_av = BMB_av / REAL(n,dp)
          BMB_shelf_extra( j,i) = BMB_av
          
        END IF ! IF (mask_FCMP( j,i) == 1) THEN
      END IF ! IF (mask_PMP( j,i) == 1) THEN
      
    END DO
    END DO
    CALL sync
    
    ! Copy results back to original array
    BMB%BMB_shelf( :,grid%i1:grid%i2) = BMB_shelf_extra( :,grid%i1:grid%i2)
    
    ! Clean up after yourself
    CALL deallocate_shared( wmask_FCMP      )
    CALL deallocate_shared( wmask_PMP       )
    CALL deallocate_shared( wBMB_shelf_extra)
    CALL deallocate_shared( wdBMBdx         )
    CALL deallocate_shared( wdBMBdy         )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE extrapolate_melt_from_FCMP_to_PMP
  SUBROUTINE fill_shelf_holes( grid, mask_ocean, mask_shelf)
    ! Fill "holes" in masks (small open ocean pockets surrounded by shelf)
    !
    ! Only ocean that is connected via other open ocean pixels to the domain border
    ! count as "true" open ocean
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_ocean
    INTEGER,  DIMENSION(:,:  ),          INTENT(INOUT) :: mask_shelf
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'fill_shelf_holes'
    INTEGER                                            :: i,j,ii,jj
    INTEGER,  DIMENSION(:,:  ), POINTER                ::  map
    INTEGER                                            :: wmap
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE            :: stack
    INTEGER                                            :: stackN
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx, map, wmap)
    
    ! Let the master do the work
    IF (par%master) THEN
      
      ! Allocate memory for the stack
      ALLOCATE( stack( grid%nx*grid%ny,2))
      
      ! Initialise the map and stack
      map    = 0
      stack  = 0
      stackN = 0
      
      ! West
      i = 1
      DO j = 1, grid%ny
        stackN = stackN + 1
        stack( stackN,:) = [i,j]
        map( j,i) = 1
      END DO
      
      ! East
      i = grid%nx
      DO j = 1, grid%ny
        stackN = stackN + 1
        stack( stackN,:) = [i,j]
        map( j,i) = 1
      END DO
      
      ! South
      j = 1
      DO i = 1, grid%nx
        stackN = stackN + 1
        stack( stackN,:) = [i,j]
        map( j,i) = 1
      END DO
      
      ! North
      j = grid%ny
      DO i = 1, grid%nx
        stackN = stackN + 1
        stack( stackN,:) = [i,j]
        map( j,i) = 1
      END DO
      
      ! Perform a flood fill
      DO WHILE (stackN > 0)
        
        ! Take the last grid cell in the stack
        i = stack( stackN,1)
        j = stack( stackN,2)
        map( j,i) = -1
        stackN = stackN - 1
        
        IF (mask_ocean( j,i) == 1 .AND. mask_shelf( j,i) == 0) THEN
          ! This grid cell is open ocean; mark it, and add its non-stacked neighbours to the stack
          
          map( j,i) = 2
          DO ii = MAX( 1,i-1), MIN( grid%nx,i+1)
          DO jj = MAX( 1,j-1), MIN( grid%ny,j+1)
            IF (map( jj,ii) == 0) THEN
              map( jj,ii) = 1
              stackN = stackN + 1
              stack( stackN,:) = [ii,jj]
            END IF
          END DO
          END DO
          
        END IF
        
      END DO ! DO WHILE (stackN > 0)
      
      ! Clean up after yourself
      DEALLOCATE( stack)
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Remove ocean pockets
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      IF (mask_ocean( j,i) == 1 .AND. mask_shelf( j,i) == 0 .AND. map( j,i) == 0) THEN
        ! This grid cell is technically ice-free, but since it has no connection to the
        ! open ocean, let PICO pretend it's really shelf
        mask_shelf( j,i) = 1
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wmap)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
    
  END SUBROUTINE fill_shelf_holes
  SUBROUTINE list_GL_grid_cells_basin( grid, ice, mask_sheet, mask_shelf, basin_i, list_GL, n_GL)
    ! List the indices of all grounding-line grid cells within the specified basin
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_sheet
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_shelf
    INTEGER,                             INTENT(IN)    :: basin_i
    INTEGER,  DIMENSION(:,:  ),          INTENT(INOUT) :: list_GL
    INTEGER,                             INTENT(INOUT) :: n_GL
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'list_GL_grid_cells_basin'
    INTEGER                                            :: i,j,ii,jj
    LOGICAL                                            :: is_GL
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (par%master) THEN
    
      ! Initialise
      list_GL = 0
      n_GL    = 0
      
      ! Find and list all grounding-line grid cells within the specified basin
      DO i = 1, grid%nx
      DO j = 1, grid%ny
        
        IF (mask_shelf( j,i) == 1 .AND. ice%basin_ID( j,i) == basin_i) THEN
          
          is_GL = .FALSE.
          
          DO ii = MAX( 1,i-1), MIN( grid%nx,i+1)
          DO jj = MAX( 1,j-1), MIN( grid%ny,j+1)
            IF (mask_sheet( jj,ii) == 1) THEN
              ! This is a grounding-line grid cell
              is_GL = .TRUE.
              EXIT
            END IF
          END DO
          IF (is_GL) EXIT
          END DO
          
          IF (is_GL) THEN
            n_GL = n_GL + 1
            list_GL( n_GL,:) = [i,j]
          END IF
          
        END IF ! IF (mask_shelf( j,i) == 1 .AND. ice%basin_ID( j,i) == basin_i) THEN
        
      END DO
      END DO
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE list_GL_grid_cells_basin
  SUBROUTINE list_IF_grid_cells_basin( grid, ice, mask_shelf, mask_ocean, basin_i, list_IF, n_IF)
    ! List the indices of all (floating) ice-front grid cells within the specified basin
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_shelf
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: mask_ocean
    INTEGER,                             INTENT(IN)    :: basin_i
    INTEGER,  DIMENSION(:,:  ),          INTENT(INOUT) :: list_IF
    INTEGER,                             INTENT(INOUT) :: n_IF
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'list_IF_grid_cells_basin'
    INTEGER                                            :: i,j,ii,jj
    LOGICAL                                            :: is_IF
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (par%master) THEN
    
      ! Initialise
      list_IF = 0
      n_IF    = 0
      
      ! Find and list all (floating) ice-front grid cells within the specified basin
      DO i = 1, grid%nx
      DO j = 1, grid%ny
        
        IF (mask_shelf( j,i) == 1 .AND. ice%basin_ID( j,i) == basin_i) THEN
          
          is_IF = .FALSE.
          
          DO ii = MAX( 1,i-1), MIN( grid%nx,i+1)
          DO jj = MAX( 1,j-1), MIN( grid%ny,j+1)
            IF (mask_shelf( jj,ii) == 0 .AND. mask_ocean( jj,ii) == 1) THEN
              ! This is an ice-front grid cell
              is_IF = .TRUE.
              EXIT
            END IF
          END DO
          IF (is_IF) EXIT
          END DO
          
          IF (is_IF) THEN
            n_IF = n_IF + 1
            list_IF( n_IF,:) = [i,j]
          END IF
          
        END IF ! IF (mask_shelf( j,i) == 1 .AND. ice%basin_ID( j,i) == basin_i) THEN
        
      END DO
      END DO
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE list_IF_grid_cells_basin
  SUBROUTINE calc_shortest_distance_to_GL( grid, list_GL, n_GL, i,j, d_GL)
    ! Calculate the shortest linear distance from grid cell i,j to the grounding line
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: list_GL
    INTEGER,                             INTENT(IN)    :: n_GL
    INTEGER,                             INTENT(IN)    :: i,j
    REAL(dp),                            INTENT(OUT)   :: d_GL
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_shortest_distance_to_GL'
    INTEGER                                            :: n,ii,jj
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Initialise
    d_GL = REAL( MAX( grid%ny, grid%nx), dp) * grid%dx
    
    ! Find the shortest linear distance from grid cell i,j to the grounding line
    DO n = 1, n_GL
      ii = list_GL( n,1)
      jj = list_GL( n,2)
      d_GL = MIN( d_GL, SQRT( REAL( ii-i,dp)**2 + REAL( jj-j,dp)**2) * grid%dx)
    END DO
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE calc_shortest_distance_to_GL
  SUBROUTINE calc_shortest_distance_to_IF( grid, list_IF, n_IF, i,j, d_IF)
    ! Calculate the shortest linear distance from grid cell i,j to the ice front
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)    :: list_IF
    INTEGER,                             INTENT(IN)    :: n_IF
    INTEGER,                             INTENT(IN)    :: i,j
    REAL(dp),                            INTENT(OUT)   :: d_IF
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_shortest_distance_to_IF'
    INTEGER                                            :: n,ii,jj
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Initialise
    d_IF = REAL( MAX( grid%ny, grid%nx), dp) * grid%dx
    
    ! Find the shortest linear distance from grid cell i,j to the ice front
    DO n = 1, n_IF
      ii = list_IF( n,1)
      jj = list_IF( n,2)
      d_IF = MIN( d_IF, SQRT( REAL( ii-i,dp)**2 + REAL( jj-j,dp)**2) * grid%dx)
    END DO
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE calc_shortest_distance_to_IF
  SUBROUTINE calc_GL_depth( grid, ice, i, j, Hb_GL)
    ! Interpolate the water depth at the nearest exact GL position
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    INTEGER,                             INTENT(IN)    :: i,j
    REAL(dp),                            INTENT(OUT)   :: Hb_GL
    
    ! Local variables:
    REAL(dp)                                           :: Hb_GL_sum
    REAL(dp)                                           :: w_sum
    REAL(dp)                                           :: TAF1, TAF2, Hb1, Hb2, lambda_GL
    
    ! Safety
    IF (ice%TAF_a( j,i) > 0._dp) THEN
      CALL crash('only defined for shelf cells!')
    END IF
    
    ! Initialise
    Hb_GL_sum = 0._dp
    w_sum     = 0._dp
    
    TAF1 = ice%TAF_a( j,i)
    Hb1  = ice%Hb_a(  j,i)
    
    ! West
    IF (i > 1) THEN
      IF (ice%TAF_a( j,i-1) > 0._dp) THEN
        
        TAF2 = ice%TAF_a( j,i-1)
        Hb2  = ice%Hb_a(  j,i-1)
        
        lambda_GL = TAF1 / (TAF1 - TAF2)
        Hb_GL_sum = Hb_GL_sum + (1._dp - lambda_GL) * Hb1 + lambda_GL * Hb2
        w_sum     = w_sum + 1._dp
        
      END IF
    END IF
    
    ! East
    IF (i < grid%nx) THEN
      IF (ice%TAF_a( j,i+1) > 0._dp) THEN
        
        TAF2 = ice%TAF_a( j,i+1)
        Hb2  = ice%Hb_a(  j,i+1)
        
        lambda_GL = TAF1 / (TAF1 - TAF2)
        Hb_GL_sum = Hb_GL_sum + (1._dp - lambda_GL) * Hb1 + lambda_GL * Hb2
        w_sum     = w_sum + 1._dp
        
      END IF
    END IF
    
    ! South
    IF (j > 1) THEN
      IF (ice%TAF_a( j-1,i) > 0._dp) THEN
        
        TAF2 = ice%TAF_a( j-1,i)
        Hb2  = ice%Hb_a(  j-1,i)
        
        lambda_GL = TAF1 / (TAF1 - TAF2)
        Hb_GL_sum = Hb_GL_sum + (1._dp - lambda_GL) * Hb1 + lambda_GL * Hb2
        w_sum     = w_sum + 1._dp
        
      END IF
    END IF
    
    ! North
    IF (j < grid%ny) THEN
      IF (ice%TAF_a( j+1,i) > 0._dp) THEN
        
        TAF2 = ice%TAF_a( j+1,i)
        Hb2  = ice%Hb_a(  j+1,i)
        
        lambda_GL = TAF1 / (TAF1 - TAF2)
        Hb_GL_sum = Hb_GL_sum + (1._dp - lambda_GL) * Hb1 + lambda_GL * Hb2
        w_sum     = w_sum + 1._dp
        
      END IF
    END IF
    
    ! Southwest
    IF (i > 1 .AND. j > 1) THEN
      IF (ice%TAF_a( j-1,i-1) > 0._dp) THEN
        
        TAF2 = ice%TAF_a( j-1,i-1)
        Hb2  = ice%Hb_a(  j-1,i-1)
        
        lambda_GL = TAF1 / (TAF1 - TAF2)
        Hb_GL_sum = Hb_GL_sum + (1._dp - lambda_GL) * Hb1 + lambda_GL * Hb2
        w_sum     = w_sum + 1._dp
        
      END IF
    END IF
    
    ! Southeast
    IF (i < grid%nx .AND. j > 1) THEN
      IF (ice%TAF_a( j-1,i+1) > 0._dp) THEN
        
        TAF2 = ice%TAF_a( j-1,i+1)
        Hb2  = ice%Hb_a(  j-1,i+1)
        
        lambda_GL = TAF1 / (TAF1 - TAF2)
        Hb_GL_sum = Hb_GL_sum + (1._dp - lambda_GL) * Hb1 + lambda_GL * Hb2
        w_sum     = w_sum + 1._dp
        
      END IF
    END IF
    
    ! Northwest
    IF (i > 1 .AND. j < grid%ny) THEN
      IF (ice%TAF_a( j+1,i-1) > 0._dp) THEN
        
        TAF2 = ice%TAF_a( j+1,i-1)
        Hb2  = ice%Hb_a(  j+1,i-1)
        
        lambda_GL = TAF1 / (TAF1 - TAF2)
        Hb_GL_sum = Hb_GL_sum + (1._dp - lambda_GL) * Hb1 + lambda_GL * Hb2
        w_sum     = w_sum + 1._dp
        
      END IF
    END IF
    
    ! Northeast
    IF (i < grid%nx .AND. j < grid%ny) THEN
      IF (ice%TAF_a( j+1,i+1) > 0._dp) THEN
        
        TAF2 = ice%TAF_a( j+1,i+1)
        Hb2  = ice%Hb_a(  j+1,i+1)
        
        lambda_GL = TAF1 / (TAF1 - TAF2)
        Hb_GL_sum = Hb_GL_sum + (1._dp - lambda_GL) * Hb1 + lambda_GL * Hb2
        w_sum     = w_sum + 1._dp
        
      END IF
    END IF
    
    ! Average
    IF (w_sum > 0._dp) THEN
      Hb_GL = Hb_GL_sum / w_sum
    ELSE
      CALL crash('whaa!')
    END IF
    
  END SUBROUTINE calc_GL_depth
  
! == Invert for BMB to achieve correct shelf geometry
! ===================================================
  
  SUBROUTINE inverse_BMB_shelf_geometry( grid, ice, BMB, refgeo_init)
    ! Invert for the BMB that is needed to achieve a steady-state shelf geometry
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid 
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_BMB_model),                INTENT(INOUT) :: BMB
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'inverse_BMB_shelf_geometry'
    INTEGER                                            :: i,j
    REAL(dp), PARAMETER                                :: H0    = 10._dp
    REAL(dp), PARAMETER                                :: dHdt0 = 10._dp
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      IF (ice%mask_shelf_a( j,i) == 1) THEN
        
        BMB%BMB_shelf( j,i) = BMB%BMB_shelf( j,i) - ((ice%Hi_a( j,i) - refgeo_init%Hi( j,i)) / H0) - (ice%dHi_dt_a( j,i) / dHdt0)
        
      END IF
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE inverse_BMB_shelf_geometry
  
END MODULE BMB_module
