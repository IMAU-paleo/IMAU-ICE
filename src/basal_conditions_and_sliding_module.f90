MODULE basal_conditions_and_sliding_module

  ! Contains all the routines for calculating the basal conditions underneath the ice.

  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_ice_model, type_reference_geometry, type_BIV_target_velocity, &
                                             type_BIV_bed_roughness
  USE netcdf_module,                   ONLY: debug, write_to_debug_file, inquire_BIV_target_velocity, &
                                             read_BIV_target_velocity, create_BIV_bed_roughness_file, inquire_BIV_bed_roughness_file, &
                                             read_BIV_bed_roughness_file
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             SSA_Schoof2006_analytical_solution, smooth_Gaussian_2D, &
                                             extrapolate_Gaussian_floodfill, transpose_dp_2D, &
                                             map_square_to_square_cons_2nd_order_2D
  USE general_ice_model_data_module,   ONLY: surface_elevation

  IMPLICIT NONE
  
CONTAINS

  ! The main routine, to be called from the ice_velocity_module
  SUBROUTINE calc_basal_conditions( grid, ice)
    ! Determine the basal conditions underneath the ice

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Basal hydrology
    CALL calc_basal_hydrology( grid, ice)
    
    ! Bed roughness
    CALL calc_bed_roughness( grid, ice)
    
  END SUBROUTINE calc_basal_conditions
  SUBROUTINE initialise_basal_conditions( grid, ice)
    ! Allocation and initialisation

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Basal hydrology
    CALL initialise_basal_hydrology( grid, ice)
    
    ! Bed roughness
    CALL initialise_bed_roughness( grid, ice)
    
    ! Basal inversion
    IF (C%do_BIVgeo) THEN
      CALL initialise_basal_inversion( grid, ice)
    END IF
    
  END SUBROUTINE initialise_basal_conditions

! == Basal hydrology
! ==================

  SUBROUTINE calc_basal_hydrology( grid, ice)
    ! Calculate the pore water pressure and effective basal pressure

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    ! Calculate pore water pressure using the chosen basal hydrology model
    ! ====================================================================
    
    IF     (C%choice_basal_hydrology == 'saturated') THEN
      ! Assume all marine till is saturated (i.e. pore water pressure is equal to water pressure at depth everywhere)
      CALL calc_pore_water_pressure_saturated( grid, ice)
    ELSEIF (C%choice_basal_hydrology == 'Martin2011') THEN
      ! The Martin et al. (2011) parameterisation of pore water pressure
      CALL calc_pore_water_pressure_Martin2011( grid, ice)
    ELSE
      IF (par%master) WRITE(0,*) 'calc_basal_hydrology - ERROR: unknown choice_basal_hydrology "', TRIM(C%choice_basal_hydrology), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Calculate overburden and effective pressure
    ! ===========================================
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%overburden_pressure_a( j,i) = ice_density * grav * ice%Hi_a( j,i)
      ice%Neff_a(                j,i) = MAX(0._dp, ice%overburden_pressure_a( j,i) - ice%pore_water_pressure_a( j,i))
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calc_basal_hydrology
  SUBROUTINE initialise_basal_hydrology( grid, ice)
    ! Allocation and initialisation

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Allocate shared memory
    IF     (C%choice_basal_hydrology == 'saturated') THEN
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%pore_water_pressure_a, ice%wpore_water_pressure_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%overburden_pressure_a, ice%woverburden_pressure_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%Neff_a               , ice%wNeff_a               )
    ELSEIF (C%choice_basal_hydrology == 'Martin2011') THEN
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%pore_water_pressure_a, ice%wpore_water_pressure_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%overburden_pressure_a, ice%woverburden_pressure_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%Neff_a               , ice%wNeff_a               )
    ELSE
      IF (par%master) WRITE(0,*) 'initialise_basal_hydrology - ERROR: unknown choice_basal_hydrology "', TRIM(C%choice_basal_hydrology), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE initialise_basal_hydrology
  
  SUBROUTINE calc_pore_water_pressure_saturated( grid, ice)
    ! Calculate the pore water pressure
    !
    ! Assume all till is saturated, i.e. pore water pressure = -rho_w * g * Hb

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%pore_water_pressure_a( j,i) = -seawater_density * grav * (ice%SL_a( j,i) - ice%Hb_a( j,i))
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calc_pore_water_pressure_saturated
  SUBROUTINE calc_pore_water_pressure_Martin2011( grid, ice)
    ! Calculate the pore water pressure
    !
    ! Use the parameterisation from Martin et al. (2011)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: lambda_p
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      ! Pore water pressure scaling factor (Martin et al., 2011, Eq. 12)
      lambda_p = MIN( 1._dp, MAX( 0._dp, 1._dp - (ice%Hb_a( j,i) - ice%SL_a( j,i) - C%Martin2011_hydro_Hb_min) / (C%Martin2011_hydro_Hb_max - C%Martin2011_hydro_Hb_min) ))
  
      ! Pore water pressure (Martin et al., 2011, Eq. 11)
      ice%pore_water_pressure_a( j,i) = 0.96_dp * ice_density * grav * ice%Hi_a( j,i) * lambda_p 
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calc_pore_water_pressure_Martin2011
  
! == Bed roughness
! ================

  SUBROUTINE calc_bed_roughness( grid, ice)
    ! Calculate the bed roughness

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! In case of no sliding or "idealised" sliding (e.g. ISMIP-HOM experiments), no bed roughness is required
    IF (C%choice_sliding_law == 'no_sliding' .OR. &
        C%choice_sliding_law == 'idealised') RETURN
    
    IF (C%choice_basal_roughness == 'uniform') THEN
      ! Apply a uniform bed roughness
      
      IF     (C%choice_sliding_law == 'no_sliding') THEN
        IF (par%master) WRITE(0,*) 'initialise_basal_inversion - ERROR: cannot do basal inversion when no sliding is allowed!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (C%choice_sliding_law == 'Weertman') THEN
        ice%beta_sq_a(   :,grid%i1:grid%i2) = C%uniform_Weertman_beta_sq
      ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
              C%choice_sliding_law == 'Coulomb_regularised') THEN
        ice%phi_fric_a(  :,grid%i1:grid%i2) = C%uniform_Coulomb_phi_fric
      ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
        ice%alpha_sq_a(  :,grid%i1:grid%i2) = C%uniform_Tsai2015_alpha_sq
        ice%beta_sq_a(   :,grid%i1:grid%i2) = C%uniform_Tsai2015_beta_sq
      ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
        ice%alpha_sq_a(  :,grid%i1:grid%i2) = C%uniform_Schoof2005_alpha_sq
        ice%beta_sq_a(   :,grid%i1:grid%i2) = C%uniform_Schoof2005_beta_sq
      ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
        ice%phi_fric_a(  :,grid%i1:grid%i2) = C%uniform_Coulomb_phi_fric
      ELSE
        IF (par%master) WRITE(0,*) 'initialise_basal_inversion - ERROR: unknown choice_sliding_law "', TRIM(C%choice_sliding_law), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    
    ELSEIF (C%choice_basal_roughness == 'parameterised') THEN
      ! Apply the chosen parameterisation of bed roughness
      
      IF     (C%choice_param_basal_roughness == 'Martin2011') THEN
        ! The Martin et al. (2011) parameterisation of basal roughness (specifically the till friction angle and till yield stress)
        CALL calc_bed_roughness_Martin2011( grid, ice)
      ELSEIF (C%choice_param_basal_roughness == 'SSA_icestream') THEN
        ! The basal roughness parameterisation in the SSA_icestream idealised-geometry experiment
        CALL calc_bed_roughness_SSA_icestream( grid, ice)
      ELSEIF (C%choice_param_basal_roughness == 'MISMIPplus') THEN
        ! The basal roughness parameterisation in the MISMIP+ idealised-geometry experiment
        CALL calc_bed_roughness_MISMIPplus( grid, ice)
      ELSEIF (C%choice_param_basal_roughness == 'BIVMIP_A') THEN
        ! The basal roughness parameterisation in the BIVMIP_A idealised-geometry experiment
        CALL calc_bed_roughness_BIVMIP_A( grid, ice)
      ELSEIF (C%choice_param_basal_roughness == 'BIVMIP_B') THEN
        ! The basal roughness parameterisation in the BIVMIP_B idealised-geometry experiment
        CALL calc_bed_roughness_BIVMIP_B( grid, ice)
      ELSEIF (C%choice_param_basal_roughness == 'BIVMIP_C') THEN
        ! The basal roughness parameterisation in the BIVMIP_C idealised-geometry experiment
        CALL calc_bed_roughness_BIVMIP_C( grid, ice)
      ELSE
        IF (par%master) WRITE(0,*) 'calc_bed_roughness - ERROR: unknown choice_param_basal_roughness "', TRIM(C%choice_param_basal_roughness), '"!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    ELSEIF (C%choice_basal_roughness == 'prescribed') THEN
      ! Basal roughness has been initialised from an external file; no need to do anything
      
    ELSEIF (C%choice_basal_roughness == 'inversion') THEN
      ! Basal roughness is updated by the inversion routines; no need to do anything
      
    ELSE
      IF (par%master) WRITE(0,*) 'calc_bed_roughness - ERROR: unknown choice_basal_roughness "', TRIM(C%choice_basal_roughness), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_bed_roughness
  SUBROUTINE initialise_bed_roughness( grid, ice)
    ! Allocation and initialisation

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Allocate shared memory
    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No sliding allowed
    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! Sliding laws for some idealised experiments
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%beta_sq_a , ice%wbeta_sq_a )
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Regularised Coulomb-type sliding law
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%phi_fric_a, ice%wphi_fric_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%tauc_a    , ice%wtauc_a    )
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%alpha_sq_a, ice%walpha_sq_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%beta_sq_a , ice%wbeta_sq_a )
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%alpha_sq_a, ice%walpha_sq_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%beta_sq_a , ice%wbeta_sq_a )
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%phi_fric_a, ice%wphi_fric_a)
      CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%tauc_a    , ice%wtauc_a    )
    ELSE
      IF (par%master) WRITE(0,*) 'initialise_bed_roughness - ERROR: unknown choice_sliding_law "', TRIM(C%choice_sliding_law), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! If bed roughness is prescribed, read it from the provided NetCDF file
    IF (C%choice_basal_roughness == 'prescribed') THEN
      CALL initialise_bed_roughness_from_file( grid, ice)
    END IF
    
  END SUBROUTINE initialise_bed_roughness
  SUBROUTINE initialise_bed_roughness_from_file( grid, ice)
    ! Initialise bed roughness with data from an external NetCDF file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    IF     (C%choice_sliding_law == 'no_sliding' .OR. &
            C%choice_sliding_law == 'idealised') THEN
      ! No sliding allowed / sliding laws for some idealised experiments
      IF (par%master) WRITE(0,*) 'initialise_bed_roughness_from_file - ERROR: not defined for choice_sliding_law "', TRIM(C%choice_sliding_law), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      CALL initialise_bed_roughness_from_file_Weertman( grid, ice)
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Coulomb-type sliding law
      CALL initialise_bed_roughness_from_file_Coulomb( grid, ice)
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      CALL initialise_bed_roughness_from_file_Tsai2015( grid, ice)
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      CALL initialise_bed_roughness_from_file_Schoof2005( grid, ice)
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      CALL initialise_bed_roughness_from_file_ZoetIverson( grid, ice)
    ELSE
      IF (par%master) WRITE(0,*) 'initialise_bed_roughness_from_file - ERROR: unknown choice_sliding_law "', TRIM(C%choice_sliding_law), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
  
  END SUBROUTINE initialise_bed_roughness_from_file
  SUBROUTINE initialise_bed_roughness_from_file_Weertman( grid, ice)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Weertman-type sliding law: bed roughness described by beta_sq

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    TYPE(type_BIV_bed_roughness)                       :: BIV
    
    ! Determine filename
    BIV%netcdf%filename = C%basal_roughness_filename
    
    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'
      
    ! Inquire grid data from the NetCDF file
    CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
    CALL allocate_shared_int_0D( BIV%ny, BIV%wny)
    
    IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
    CALL sync
    
    ! Allocate memory - grid
    CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
    CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%beta_sq , BIV%wbeta_sq )
    
    ! Read grid & bed roughness data from file
    IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
    CALL sync
  
    ! Safety
    CALL check_for_NaN_dp_2D( BIV%beta_sq,  'BIV%beta_sq',  'initialise_bed_roughness_from_file_Weertman')
    
    ! Since we want data represented as [j,i] internally, transpose the data we just read.
    CALL transpose_dp_2D( BIV%beta_sq,  BIV%wbeta_sq )
    
    ! Map (transposed) raw data to the model grid
    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, grid%nx, grid%ny, grid%x, grid%y, BIV%beta_sq , ice%beta_sq_a )
    
    ! Deallocate raw data
    CALL deallocate_shared( BIV%wnx      )
    CALL deallocate_shared( BIV%wny      )
    CALL deallocate_shared( BIV%wx       )
    CALL deallocate_shared( BIV%wy       )
    CALL deallocate_shared( BIV%wbeta_sq )
  
  END SUBROUTINE initialise_bed_roughness_from_file_Weertman
  SUBROUTINE initialise_bed_roughness_from_file_Coulomb( grid, ice)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Coulomb-type sliding law: bed roughness described by phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    TYPE(type_BIV_bed_roughness)                       :: BIV
    
    ! Determine filename
    BIV%netcdf%filename = C%basal_roughness_filename
    
    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'
      
    ! Inquire grid data from the NetCDF file
    CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
    CALL allocate_shared_int_0D( BIV%ny, BIV%wny)
    
    IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
    CALL sync
    
    ! Allocate memory - grid
    CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
    CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%phi_fric, BIV%wphi_fric)
    
    ! Read grid & bed roughness data from file
    IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
    CALL sync
  
    ! Safety
    CALL check_for_NaN_dp_2D( BIV%phi_fric, 'BIV%phi_fric', 'initialise_bed_roughness_from_file_Coulomb')
    
    ! Since we want data represented as [j,i] internally, transpose the data we just read.
    CALL transpose_dp_2D( BIV%phi_fric, BIV%wphi_fric)
    
    ! Map (transposed) raw data to the model grid
    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, grid%nx, grid%ny, grid%x, grid%y, BIV%phi_fric, ice%phi_fric_a)
    
    ! Deallocate raw data
    CALL deallocate_shared( BIV%wnx      )
    CALL deallocate_shared( BIV%wny      )
    CALL deallocate_shared( BIV%wx       )
    CALL deallocate_shared( BIV%wy       )
    CALL deallocate_shared( BIV%wphi_fric)
  
  END SUBROUTINE initialise_bed_roughness_from_file_Coulomb
  SUBROUTINE initialise_bed_roughness_from_file_Tsai2015( grid, ice)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Tsai 2015 sliding law: bed roughness described by alpha_sq & beta_sq

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    TYPE(type_BIV_bed_roughness)                       :: BIV
    
    ! Determine filename
    BIV%netcdf%filename = C%basal_roughness_filename
    
    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'
      
    ! Inquire grid data from the NetCDF file
    CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
    CALL allocate_shared_int_0D( BIV%ny, BIV%wny)
    
    IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
    CALL sync
    
    ! Allocate memory - grid
    CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
    CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%alpha_sq, BIV%walpha_sq)
    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%beta_sq , BIV%wbeta_sq )
    
    ! Read grid & bed roughness data from file
    IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
    CALL sync
  
    ! Safety
    CALL check_for_NaN_dp_2D( BIV%alpha_sq, 'BIV%alpha_sq', 'initialise_bed_roughness_from_file_Tsai2015')
    CALL check_for_NaN_dp_2D( BIV%beta_sq,  'BIV%beta_sq',  'initialise_bed_roughness_from_file_Tsai2015')
    
    ! Since we want data represented as [j,i] internally, transpose the data we just read.
    CALL transpose_dp_2D( BIV%alpha_sq, BIV%walpha_sq)
    CALL transpose_dp_2D( BIV%beta_sq,  BIV%wbeta_sq )
    
    ! Map (transposed) raw data to the model grid
    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, grid%nx, grid%ny, grid%x, grid%y, BIV%alpha_sq, ice%alpha_sq_a)
    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, grid%nx, grid%ny, grid%x, grid%y, BIV%beta_sq , ice%beta_sq_a )
    
    ! Deallocate raw data
    CALL deallocate_shared( BIV%wnx      )
    CALL deallocate_shared( BIV%wny      )
    CALL deallocate_shared( BIV%wx       )
    CALL deallocate_shared( BIV%wy       )
    CALL deallocate_shared( BIV%walpha_sq)
    CALL deallocate_shared( BIV%wbeta_sq )
  
  END SUBROUTINE initialise_bed_roughness_from_file_Tsai2015
  SUBROUTINE initialise_bed_roughness_from_file_Schoof2005( grid, ice)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Schoof 2005 sliding law: bed roughness described by alpha_sq & beta_sq

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    TYPE(type_BIV_bed_roughness)                       :: BIV
    
    ! Determine filename
    BIV%netcdf%filename = C%basal_roughness_filename
    
    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'
      
    ! Inquire grid data from the NetCDF file
    CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
    CALL allocate_shared_int_0D( BIV%ny, BIV%wny)
    
    IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
    CALL sync
    
    ! Allocate memory - grid
    CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
    CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%alpha_sq, BIV%walpha_sq)
    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%beta_sq , BIV%wbeta_sq )
    
    ! Read grid & bed roughness data from file
    IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
    CALL sync
  
    ! Safety
    CALL check_for_NaN_dp_2D( BIV%alpha_sq, 'BIV%alpha_sq', 'initialise_bed_roughness_from_file_Schoof2005')
    CALL check_for_NaN_dp_2D( BIV%beta_sq,  'BIV%beta_sq',  'initialise_bed_roughness_from_file_Schoof2005')
    
    ! Since we want data represented as [j,i] internally, transpose the data we just read.
    CALL transpose_dp_2D( BIV%alpha_sq, BIV%walpha_sq)
    CALL transpose_dp_2D( BIV%beta_sq,  BIV%wbeta_sq )
    
    ! Map (transposed) raw data to the model grid
    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, grid%nx, grid%ny, grid%x, grid%y, BIV%alpha_sq, ice%alpha_sq_a)
    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, grid%nx, grid%ny, grid%x, grid%y, BIV%beta_sq , ice%beta_sq_a )
    
    ! Deallocate raw data
    CALL deallocate_shared( BIV%wnx      )
    CALL deallocate_shared( BIV%wny      )
    CALL deallocate_shared( BIV%wx       )
    CALL deallocate_shared( BIV%wy       )
    CALL deallocate_shared( BIV%walpha_sq)
    CALL deallocate_shared( BIV%wbeta_sq )
  
  END SUBROUTINE initialise_bed_roughness_from_file_Schoof2005
  SUBROUTINE initialise_bed_roughness_from_file_ZoetIverson( grid, ice)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Zoet-Iverson sliding law: bed roughness described by phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    TYPE(type_BIV_bed_roughness)                       :: BIV
    
    ! Determine filename
    BIV%netcdf%filename = C%basal_roughness_filename
    
    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'
      
    ! Inquire grid data from the NetCDF file
    CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
    CALL allocate_shared_int_0D( BIV%ny, BIV%wny)
    
    IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
    CALL sync
    
    ! Allocate memory - grid
    CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
    CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%phi_fric, BIV%wphi_fric)
    
    ! Read grid & bed roughness data from file
    IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
    CALL sync
  
    ! Safety
    CALL check_for_NaN_dp_2D( BIV%phi_fric, 'BIV%phi_fric', 'initialise_bed_roughness_from_file_ZoetIverson')
    
    ! Since we want data represented as [j,i] internally, transpose the data we just read.
    CALL transpose_dp_2D( BIV%phi_fric, BIV%wphi_fric)
    
    ! Map (transposed) raw data to the model grid
    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, grid%nx, grid%ny, grid%x, grid%y, BIV%phi_fric, ice%phi_fric_a)
    
    ! Deallocate raw data
    CALL deallocate_shared( BIV%wnx      )
    CALL deallocate_shared( BIV%wny      )
    CALL deallocate_shared( BIV%wx       )
    CALL deallocate_shared( BIV%wy       )
    CALL deallocate_shared( BIV%wphi_fric)
  
  END SUBROUTINE initialise_bed_roughness_from_file_ZoetIverson
  
  ! The Martin et al. (2011) till parameterisation
  SUBROUTINE calc_bed_roughness_Martin2011( grid, ice)
    ! Calculate the till friction angle phi_fric and till yield stress tauc,
    ! using the till model by Martin et al. (2011).
    ! 
    ! Only applicable when choice_sliding_law = "Coulomb" or "Coulomb_regularised"

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: w_Hb
    
    ! Safety
    IF (.NOT. (C%choice_sliding_law == 'Coulomb' .OR. C%choice_sliding_law == 'Coulomb_regularised' .OR. C%choice_sliding_law == 'Zoet-Iverson')) THEN
      IF (par%master) WRITE(0,*) 'calc_bed_roughness_Martin2011 - ERROR: only applicable when choice_sliding_law = "Coulomb", "Coulomb_regularised", or "Zoet-Iverson"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
  
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
  
      ! Martin et al. (2011) Eq. 10
      w_Hb = MIN( 1._dp, MAX( 0._dp, (ice%Hb_a( j,i) - C%Martin2011till_phi_Hb_min) / (C%Martin2011till_phi_Hb_max - C%Martin2011till_phi_Hb_min) ))
      ice%phi_fric_a( j,i) = (1._dp - w_Hb) * C%Martin2011till_phi_min + w_Hb * C%Martin2011till_phi_max
    
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( ice%phi_fric_a, 'ice%phi_fric_a', 'Martin_2011_till_model')
    
  END SUBROUTINE calc_bed_roughness_Martin2011
  
  ! Idealised cases
  SUBROUTINE calc_bed_roughness_SSA_icestream( grid, ice)
    ! Determine the basal conditions underneath the ice
    ! 
    ! Idealised case: SSA_icestream (i.e. the Schoof 2006 analytical solution)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: dummy1
  
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      CALL SSA_Schoof2006_analytical_solution( ABS(ice%dHs_dx_a( j,i)), ice%Hi_a( j,i), ice%A_flow_vav_a( j,i), grid%y(j), dummy1, ice%tauc_a( j,i))
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE calc_bed_roughness_SSA_icestream
  SUBROUTINE calc_bed_roughness_MISMIPplus( grid, ice)
    ! Determine the basal conditions underneath the ice
    ! 
    ! Idealised case: MISMIP+ (see Asay-Davis et al., 2016)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), PARAMETER                                :: MISMIPplus_alpha_sq = 0.5_dp   ! Coulomb-law friction coefficient [unitless];         see Asay-Davis et al., 2016
    REAL(dp), PARAMETER                                :: MISMIPplus_beta_sq  = 1.0E4_dp ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3]; idem dito
  
    IF     (C%choice_sliding_law == 'Weertman') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the first (Weertman) sliding law option
       
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%beta_sq_a( j,i) = MISMIPplus_beta_sq
      END DO
      END DO
      CALL sync
          
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the second (Tsai et al., 2015) sliding law option
       
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%alpha_sq_a( j,i) = MISMIPplus_alpha_sq
        ice%beta_sq_a(  j,i) = MISMIPplus_beta_sq
      END DO
      END DO
      CALL sync
      
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the third (Schoof, 2005) sliding law option
       
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%alpha_sq_a( j,i) = MISMIPplus_alpha_sq
        ice%beta_sq_a(  j,i) = MISMIPplus_beta_sq
      END DO
      END DO
      CALL sync
      
    ELSE
      IF (par%master) WRITE(0,*) 'calc_bed_roughness_MISMIPplus - ERROR: only defined when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_bed_roughness_MISMIPplus
  SUBROUTINE calc_bed_roughness_BIVMIP_A( grid, ice)
    ! Determine the basal conditions underneath the ice
    ! 
    ! Idealised case: BIVMIP experiment A (uniform bed plus single southward ice stream)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    REAL(dp), PARAMETER                                :: xc      =    0.0E3_dp  ! Ice-stream midpoint (x-coordinate)
    REAL(dp), PARAMETER                                :: yc      = -700.0E3_dp  ! Ice-stream midpoint (y-coordinate)
    REAL(dp), PARAMETER                                :: sigma_x =   50.0E3_dp  ! Ice-stream half-width in x-direction
    REAL(dp), PARAMETER                                :: sigma_y =  400.0E3_dp  ! Ice-stream half-width in y-direction
    REAL(dp), PARAMETER                                :: phi_max =  5.0_dp      ! Uniform till friction angle outside of ice stream
    REAL(dp), PARAMETER                                :: phi_min =  0.0_dp      ! Lowest  till friction angle at ice-stream midpoint
    
    IF     (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised' .OR. &
            C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! (Regularised) Coulomb-type / Zoet-Iverson sliding law: parameterise phi_fric
  
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        
        ! Calculate till friction angle
        ice%phi_fric_a( j,i) = phi_max - (phi_max - phi_min) * EXP( -0.5_dp * &
          (((grid%x( i) - xc) / sigma_x)**2 + &
           ((grid%y( j) - yc) / sigma_y)**2))
        
      END DO
      END DO
      CALL sync
    
    ELSE
      IF (par%master) WRITE(0,*) 'calc_bed_roughness_BIVMIP_A - ERROR: choice_sliding_law "', TRIM(C%choice_sliding_law), ' not implemented"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_bed_roughness_BIVMIP_A
  SUBROUTINE calc_bed_roughness_BIVMIP_B( grid, ice)
    ! Determine the basal conditions underneath the ice
    ! 
    ! Idealised case: BIVMIP experiment B (uniform bed plus single southward ice stream)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    REAL(dp), PARAMETER                                :: xc      =    0.0E3_dp  ! Ice-stream midpoint (x-coordinate)
    REAL(dp), PARAMETER                                :: yc      = -550.0E3_dp  ! Ice-stream midpoint (y-coordinate)
    REAL(dp), PARAMETER                                :: sigma_x =   75.0E3_dp  ! Ice-stream half-width in x-direction
    REAL(dp), PARAMETER                                :: sigma_y =  300.0E3_dp  ! Ice-stream half-width in y-direction
    REAL(dp), PARAMETER                                :: phi_max =  5.0_dp      ! Uniform till friction angle outside of ice stream
    REAL(dp), PARAMETER                                :: phi_min =  0.0_dp      ! Lowest  till friction angle at ice-stream midpoint
    
    IF     (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised' .OR. &
            C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! (Regularised) Coulomb-type / Zoet-Iverson sliding law: parameterise phi_fric
  
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        
        ! Calculate till friction angle
        ice%phi_fric_a( j,i) = phi_max - (phi_max - phi_min) * EXP( -0.5_dp * &
          (((grid%x( i) - xc) / sigma_x)**2 + &
           ((grid%y( j) - yc) / sigma_y)**2))
        
      END DO
      END DO
      CALL sync
    
    ELSE
      IF (par%master) WRITE(0,*) 'calc_bed_roughness_BIVMIP_B - ERROR: choice_sliding_law "', TRIM(C%choice_sliding_law), ' not implemented"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_bed_roughness_BIVMIP_B
  SUBROUTINE calc_bed_roughness_BIVMIP_C( grid, ice)
    ! Determine the basal conditions underneath the ice
    ! 
    ! Idealised case: BIVMIP experiment C (parameterised bed following Martin et al (2011), plus a narrow southward ice-stream region)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    REAL(dp), PARAMETER                                :: xc      =     0.0E3_dp  ! Ice-stream midpoint (x-coordinate)
    REAL(dp), PARAMETER                                :: yc      = -1000.0E3_dp  ! Ice-stream midpoint (y-coordinate)
    REAL(dp), PARAMETER                                :: sigma_x =    75.0E3_dp  ! Ice-stream half-width in x-direction
    REAL(dp), PARAMETER                                :: sigma_y =   500.0E3_dp  ! Ice-stream half-width in y-direction
    
    IF     (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised' .OR. &
            C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! (Regularised) Coulomb-type / Zoet-Iverson sliding law: parameterise phi_fric
    
      ! Start with the Martin et al. (2011) till parameterisation
      CALL calc_bed_roughness_Martin2011( grid, ice)
     
      ! Then add the ice-stream region
      DO i = grid%i1, grid%i2
      DO j = 1, grid%ny
        ice%phi_fric_a( j,i) = ice%phi_fric_a( j,i) * (1._dp - EXP( -0.5_dp * &
          (((grid%x( i) - xc) / sigma_x)**2 + &
           ((grid%y( j) - yc) / sigma_y)**2)))
      END DO
      END DO
      CALL sync
    
    ELSE
      IF (par%master) WRITE(0,*) 'calc_bed_roughness_BIVMIP_C - ERROR: choice_sliding_law "', TRIM(C%choice_sliding_law), ' not implemented"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_bed_roughness_BIVMIP_C
  
! == Sliding laws
! ===============

  SUBROUTINE calc_sliding_law( grid, ice, u_a, v_a, beta_a)
    ! Calculate the sliding term beta in the SSA/DIVA using the specified sliding law
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
      
    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No sliding allowed (choice of beta is trivial)
      beta_a( :,grid%i1:grid%i2) = 0._dp
      CALL sync
    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! Sliding laws for some idealised experiments
      CALL calc_sliding_law_idealised(           grid, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      CALL calc_sliding_law_Weertman(            grid, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Coulomb') THEN
      ! Coulomb-type sliding law
      CALL calc_sliding_law_Coulomb(             grid, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Regularised Coulomb-type sliding law
      CALL calc_sliding_law_Coulomb_regularised( grid, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      CALL calc_sliding_law_Tsai2015(            grid, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      CALL calc_sliding_law_Schoof2005(          grid, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      CALL calc_sliding_law_ZoetIverson(         grid, ice, u_a, v_a, beta_a)
    ELSE
      IF (par%master) WRITE(0,*) 'calc_sliding_law - ERROR: unknown choice_sliding_law "', TRIM(C%choice_sliding_law), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_sliding_law

  SUBROUTINE calc_sliding_law_Weertman( grid, ice, u_a, v_a, beta_a)
    ! Weertman-type ("power law") sliding law
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)
      
      ! Asay-Davis et al. (2016), Eq. 6
      beta_a( j,i) = ice%beta_sq_a( j,i) * uabs_a ** (1._dp / C%slid_Weertman_m - 1._dp)
      
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_Weertman')
    
  END SUBROUTINE calc_sliding_law_Weertman
  SUBROUTINE calc_sliding_law_Coulomb( grid, ice, u_a, v_a, beta_a)
    ! Coulomb-type sliding law
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Calculate the till yield stress from the till friction angle and the effective pressure
      ice%tauc_a( j,i) = TAN((pi / 180._dp) * ice%phi_fric_a( j,i)) * ice%Neff_a( j,i)
      
      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)
      
      beta_a( j,i) = ice%tauc_a( j,i) / uabs_a
      
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_Coulomb')
    
  END SUBROUTINE calc_sliding_law_Coulomb
  SUBROUTINE calc_sliding_law_Coulomb_regularised( grid, ice, u_a, v_a, beta_a)
    ! Regularised Coulomb-type sliding law
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Calculate the till yield stress from the till friction angle and the effective pressure
      ice%tauc_a( j,i) = TAN((pi / 180._dp) * ice%phi_fric_a( j,i)) * ice%Neff_a( j,i)
      
      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)
      
      beta_a( j,i) = ice%tauc_a( j,i) * uabs_a ** (C%slid_Coulomb_reg_q_plastic - 1._dp) / (C%slid_Coulomb_reg_u_threshold ** C%slid_Coulomb_reg_q_plastic)
      
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_Coulomb_regularised')
    
  END SUBROUTINE calc_sliding_law_Coulomb_regularised
  SUBROUTINE calc_sliding_law_Tsai2015( grid, ice, u_a, v_a, beta_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Davis et al., 2016)
    ! 
    ! Asay-Davis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    ! 
    ! Tsai et al.: Marine ice-sheet profiles and stability under Coulomb basal conditions,
    ! Journal of Glaciology 61, 205–215, doi:10.3189/2015JoG14J221, 2015.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)
      
      ! Asay-Davis et al. (2016), Eq. 7
      beta_a( j,i) = MIN( ice%alpha_sq_a( j,i) * ice%Neff_a( j,i), ice%beta_sq_a( j,i) * uabs_a ** (1._dp / C%slid_Weertman_m)) * uabs_a**(-1._dp)
      
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_Tsai2015')
    
  END SUBROUTINE calc_sliding_law_Tsai2015
  SUBROUTINE calc_sliding_law_Schoof2005( grid, ice, u_a, v_a, beta_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Davis et al., 2016)
    ! 
    ! Asay-Davis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    ! 
    ! Schoof: The effect of cavitation on glacier sliding, P. Roy. Soc. A-Math. Phy., 461, 609–627, doi:10.1098/rspa.2004.1350, 2005
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a, alpha_sq, beta_sq
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)
      
      ! Abbreviations for shorter code
      alpha_sq = ice%alpha_sq_a( j,i)
      beta_sq  = ice%beta_sq_a(  j,i)
      
      ! Asay-Davis et al. (2016), Eq. 11
      beta_a( j,i) = ((beta_sq * uabs_a**(1._dp / C%slid_Weertman_m) * alpha_sq * ice%Neff_a( j,i)) / &
        ((beta_sq**C%slid_Weertman_m * uabs_a + (alpha_sq * ice%Neff_a( j,i))**C%slid_Weertman_m)**(1._dp / C%slid_Weertman_m))) * uabs_a**(-1._dp)
        
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_Schoof2005')
    
  END SUBROUTINE calc_sliding_law_Schoof2005
  SUBROUTINE calc_sliding_law_ZoetIverson( grid, ice, u_a, v_a, beta_a)
    ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: uabs_a
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs_a = SQRT( C%slid_delta_v**2 + u_a( j,i)**2 + v_a( j,i)**2)
      
      ! Zoet & Iverson (2020), Eq. (3) (divided by u to give beta = tau_b / u)
      beta_a( j,i) = ice%Neff_a( j,i) * TAN((pi / 180._dp) * ice%phi_fric_a( j,i)) * (uabs_a**(1._dp / C%slid_ZI_p - 1._dp)) * ((uabs_a + C%slid_ZI_ut)**(-1._dp / C%slid_ZI_p))
      
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_ZoetIverson')
    
  END SUBROUTINE calc_sliding_law_ZoetIverson
  
  SUBROUTINE calc_sliding_law_idealised( grid, ice, u_a, v_a, beta_a)
    ! Sliding laws for some idealised experiments
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    REAL(dp) :: dummy_dp
    
    ! To prevent compiler warnings...
    dummy_dp = u_a( 1,1)
    dummy_dp = v_a( 1,1)
      
    IF     (C%choice_idealised_sliding_law == 'ISMIP_HOM_C') THEN
      ! ISMIP-HOM experiment C
      
      CALL calc_sliding_law_idealised_ISMIP_HOM_C( grid, beta_a)
      
    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_D') THEN
      ! ISMIP-HOM experiment D
      
      CALL calc_sliding_law_idealised_ISMIP_HOM_D( grid, beta_a)
      
    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_E') THEN
      ! ISMIP-HOM experiment E
      
      CALL calc_sliding_law_idealised_ISMIP_HOM_E( grid, beta_a)
      
    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_F') THEN
      ! ISMIP-HOM experiment F
      
      CALL calc_sliding_law_idealised_ISMIP_HOM_F( grid, ice, beta_a)
      
    ELSE
      IF (par%master) WRITE(0,*) 'calc_sliding_law_idealised - ERROR: unknown choice_idealised_sliding_law "', TRIM(C%choice_idealised_sliding_law), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_sliding_law_idealised
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C( grid, beta_a)
    ! Sliding laws for some idealised experiments
    ! 
    ! ISMIP-HOM experiment C
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      beta_a( j,i) = 1000._dp + 1000._dp * SIN( 2._dp * pi * grid%x( i) / C%ISMIP_HOM_L) * SIN( 2._dp * pi * grid%y( j) / C%ISMIP_HOM_L)
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_idealised_ISMIP_HOM_C')
    
  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D( grid, beta_a)
    ! Sliding laws for some idealised experiments
    ! 
    ! ISMIP-HOM experiment D
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      beta_a( j,i) = 1000._dp + 1000._dp * SIN( 2._dp * pi * grid%x( i) / C%ISMIP_HOM_L)
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_idealised_ISMIP_HOM_D')
    
  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_E( grid, beta_a)
    ! Sliding laws for some idealised experiments
    ! 
    ! ISMIP-HOM experiment E: use the externally prescribed slip zone in Haut Glacier d'Arolla
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp)                                           :: dummy1, dummy2, dummy3
    INTEGER                                            :: ios,slides
    
    IF (par%master) THEN
      OPEN(   UNIT = 1337, FILE=C%ISMIP_HOM_E_Arolla_filename, ACTION='READ')
      DO i = 1, 51
        READ( UNIT = 1337, FMT=*, IOSTAT=ios) dummy1, dummy2, dummy3, slides
        IF (slides == 1) THEN
          DO j = 1, grid%ny
            beta_a( j,i) = 0._dp
          END DO
        ELSE
          DO j = 1, grid%ny
            beta_a( j,i) = 1.0E30_dp
          END DO
        END IF
        IF (ios /= 0) THEN
          IF (par%master) WRITE(0,*) ' calc_sliding_law_idealised_ISMIP_HOM_E - ERROR: length of text file "', TRIM(C%ISMIP_HOM_E_Arolla_filename), '" should be 51 lines!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END DO
      CLOSE( UNIT  = 1337)
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_idealised_ISMIP_HOM_E')
    
  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_E
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F( grid, ice, beta_a)
    ! Sliding laws for some idealised experiments
    ! 
    ! ISMIP-HOM experiment F
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: beta_a
    
    ! Local variables:
    INTEGER                                            :: i,j
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      beta_a( j,i) = (ice%A_flow_vav_a( j,i) * 1000._dp)**(-1._dp)
    END DO
    END DO
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( beta_a, 'beta_a', 'calc_sliding_law_idealised_ISMIP_HOM_F')
    
  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F
  
! == Basal inversion routines
! ===========================

  SUBROUTINE basal_inversion_geo( grid, ice, refgeo_init, dt)
    ! Invert for bed roughness based on modelled vs. target geometry

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    REAL(dp),                            INTENT(IN)    :: dt
    
    IF     (C%choice_BIVgeo_method == 'PDC2012') THEN
      ! Update the bed roughness according to Pollard & DeConto (2012)
      
      CALL update_bed_roughness_PDC2012( grid, ice, refgeo_init)
      
    ELSEIF (C%choice_BIVgeo_method == 'Lipscomb2021') THEN
      ! Update the bed roughness according to Lipscomb et al. (2021)
      
      CALL update_bed_roughness_Lipscomb2021( grid, ice, refgeo_init, dt)
      
    ELSEIF (C%choice_BIVgeo_method == 'CISM+') THEN
      ! Update the bed roughness according to the geometry+velocity-based CISM+ approach
      
      CALL update_bed_roughness_CISMplus( grid, ice, refgeo_init, dt)
      
    ELSE
      IF (par%master) WRITE(0,*) 'basal_inversion_geo - ERROR: unknown choice_BIVgeo_method "', TRIM(C%choice_BIVgeo_method), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE basal_inversion_geo
  
  SUBROUTINE update_bed_roughness_PDC2012( grid, ice, refgeo_init)
    ! Update the bed roughness according to Pollard & DeConto (2012)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    
    IF (C%choice_sliding_law == 'Coulomb' .OR. &
        C%choice_sliding_law == 'Coulomb_regularised') THEN
      
      CALL update_bed_roughness_PDC2012_Coulomb( grid, ice, refgeo_init)
      
    ELSE
      IF (par%master) WRITE(0,*) 'update_bed_roughness_PDC2012 - ERROR: not implemented for sliding law "', TRIM(C%choice_sliding_law), '"'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE update_bed_roughness_PDC2012
  SUBROUTINE update_bed_roughness_PDC2012_Coulomb( grid, ice, refgeo_init)
    ! Update the bed roughness according to Pollard & DeConto (2012)
    ! 
    ! For the case of a (regularised) Coulomb sliding law: update phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dz
    INTEGER                                            :: wdz
    REAL(dp), PARAMETER                                :: sigma = 500._dp
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dz, wdz)
    
    ! Calculate dz (Pollard % DeConto (2012), just after Eq. 5)
    DO i = grid%i1, grid%I2
    DO j = 1, grid%ny
    
      IF (ice%mask_sheet_a( j,i) == 1) THEN
      
        dz( j,i) = MAX( -1.5_dp, MIN( 1.5_dp, (ice%Hs_a( j,i) - refgeo_init%Hs( j,i)) / C%BIVgeo_PDC2012_hinv ))
        
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Extrapolate new values outside the ice sheet
    CALL extrapolate_updated_bed_roughness( grid, ice, dz)
    
    ! Apply regularisation
    CALL smooth_Gaussian_2D( grid, dz, sigma)
    
    ! Update bed roughness
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ! Pollard % DeConto (2012), Eq. 5
      ice%phi_fric_a( j,i) = MAX( 0.001_dp, MIN( 30._dp, ice%phi_fric_a( j,i) * (10._dp**(-dz( j,i))) ))
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE update_bed_roughness_PDC2012_Coulomb
  
  SUBROUTINE update_bed_roughness_Lipscomb2021( grid, ice, refgeo_init, dt)
    ! Update the bed roughness according to Lipscomb et al. (2021)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    REAL(dp),                            INTENT(IN)    :: dt
    
    IF (C%choice_sliding_law == 'Coulomb' .OR. &
        C%choice_sliding_law == 'Coulomb_regularised') THEN
      
      CALL update_bed_roughness_Lipscomb2021_Coulomb( grid, ice, refgeo_init, dt)
      
    ELSE
      IF (par%master) WRITE(0,*) 'update_bed_roughness_Lipscomb2021 - ERROR: not implemented for sliding law "', TRIM(C%choice_sliding_law), '"'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE update_bed_roughness_Lipscomb2021
  SUBROUTINE update_bed_roughness_Lipscomb2021_Coulomb( grid, ice, refgeo_init, dt)
    ! Update the bed roughness according to Lipscomb et al. (2021)
    ! 
    ! For the case of a (regularised) Coulomb sliding law: update phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    REAL(dp),                            INTENT(IN)    :: dt
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dCdt
    INTEGER                                            :: wdCdt
    REAL(dp), PARAMETER                                :: sigma = 500._dp
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dCdt, wdCdt)
    
    ! Calculate dCdt
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      IF (ice%mask_sheet_a( j,i) == 1) THEN
      
        ! Lipscomb et al. (2021), Eq. 8 (slightly altered for Coulomb friction)
        dCdt( j,i) = -ice%phi_fric_a( j,i) / C%BIVgeo_Lipscomb2021_H0 * ( (ice%Hi_a( j,i) - refgeo_init%Hi( j,i)) / C%BIVgeo_Lipscomb2021_tauc + 2._dp * ice%dHi_dt_a( j,i))
        
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Extrapolate new values outside the ice sheet
    CALL extrapolate_updated_bed_roughness( grid, ice, dCdt)
    
    ! Apply regularisation
    CALL smooth_Gaussian_2D( grid, dCdt, sigma)
    
    ! Update bed roughness
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%phi_fric_a( j,i) = MAX( 1E-6_dp, MIN( 30._dp, ice%phi_fric_a( j,i) + dCdt( j,i) * dt ))
    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wdCdt)
    
  END SUBROUTINE update_bed_roughness_Lipscomb2021_Coulomb
  
  SUBROUTINE update_bed_roughness_CISMplus( grid, ice, refgeo_init, dt)
    ! Update the bed roughness according to the geometry+velocity-based CISM+ approach

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    REAL(dp),                            INTENT(IN)    :: dt
    
    IF (C%choice_sliding_law == 'Coulomb' .OR. &
        C%choice_sliding_law == 'Coulomb_regularised' .OR. &
        C%choice_sliding_law == 'Zoet-Iverson') THEN
      
      CALL update_bed_roughness_CISMplus_Coulomb( grid, ice, refgeo_init, dt)
      
    ELSE
      IF (par%master) WRITE(0,*) 'update_bed_roughness_Lipscomb2021 - ERROR: not implemented for sliding law "', TRIM(C%choice_sliding_law), '"'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE update_bed_roughness_CISMplus
  SUBROUTINE update_bed_roughness_CISMplus_Coulomb( grid, ice, refgeo_init, dt)
    ! Update the bed roughness according to the geometry+velocity-based CISM+ approach
    ! 
    ! For the case of a (regularised) Coulomb / Zoet-Iverson sliding law: update phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    REAL(dp),                            INTENT(IN)    :: dt
    
    ! Local variables:
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  dCdt
    INTEGER                                            :: wdCdt
    REAL(dp), PARAMETER                                :: sigma = 500._dp
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, dCdt, wdCdt)
    
    ! Calculate dCdt
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
    
      IF (ice%mask_sheet_a( j,i) == 1) THEN
      
        ! CISM+ inversion equation
        dCdt( j,i) = -ice%phi_fric_a( j,i) / C%BIVgeo_CISMplus_tauc * &
          ( C%BIVgeo_CISMplus_wH * (ice%Hi_a(        j,i) - refgeo_init%Hi(           j,i)) / ( C%BIVgeo_CISMplus_H0) + &
            C%BIVgeo_CISMplus_wu * (ice%uabs_surf_a( j,i) - ice%BIV_uabs_surf_target( j,i)) / (-C%BIVgeo_CISMplus_u0) )
        
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Extrapolate new values outside the ice sheet
    CALL extrapolate_updated_bed_roughness( grid, ice, dCdt)
    
    ! Apply regularisation
    CALL smooth_Gaussian_2D( grid, dCdt, sigma)
    
    ! Update bed roughness
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%phi_fric_a( j,i) = MAX( 1E-6_dp, MIN( 30._dp, ice%phi_fric_a( j,i) + dCdt( j,i) * dt ))
    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wdCdt)
    
  END SUBROUTINE update_bed_roughness_CISMplus_Coulomb
  
  SUBROUTINE extrapolate_updated_bed_roughness( grid, ice, d)
    ! The different geometry-based basal inversion routine only yield values
    ! beneath grounded ice; extrapolate new values outward to cover the entire domain.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d
    
    ! Local variables:
    INTEGER                                            :: i,j
    INTEGER,  DIMENSION(:,:  ), POINTER                :: mask,  mask_filled
    INTEGER                                            :: wmask, wmask_filled
    REAL(dp)                                           :: sigma
    
    ! Allocate shared memory
    CALL allocate_shared_int_2D( grid%ny, grid%nx, mask       , wmask       )
    CALL allocate_shared_int_2D( grid%ny, grid%nx, mask_filled, wmask_filled)
    
    ! Fill the mask used for the extrapolation
    ! 
    ! NOTE: not parallelised!
    ! 
    ! Note about the mask:
    !    2 = data provided
    !    1 = no data provided, fill allowed
    !    0 = no fill allowed
    ! (so basically this routine extrapolates data from the area
    !  where mask == 2 into the area where mask == 1)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      IF (ice%mask_sheet_a( j,i) == 1) THEN
        mask( j,i) = 2
      ELSE
        mask( j,i) = 1
      END IF
    END DO
    END DO
    CALL sync
    
    ! Perform the extrapolation
    sigma = grid%dx / 2._dp
    IF (par%master) CALL extrapolate_Gaussian_floodfill( grid, mask, d, sigma, mask_filled)
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wmask)
    CALL deallocate_shared( wmask_filled)
    
  END SUBROUTINE extrapolate_updated_bed_roughness
  SUBROUTINE initialise_basal_inversion( grid, ice)
    ! Fill in the initial guess for the bed roughness

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Initial guess for bed roughness: just a uniform value
    ! =====================================================
    
    IF     (C%choice_sliding_law == 'no_sliding') THEN
      IF (par%master) WRITE(0,*) 'initialise_basal_inversion - ERROR: cannot do basal inversion when no sliding is allowed!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ice%beta_sq_a(   :,grid%i1:grid%i2) = C%uniform_Weertman_beta_sq
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      ice%phi_fric_a(  :,grid%i1:grid%i2) = C%uniform_Coulomb_phi_fric
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ice%alpha_sq_a(  :,grid%i1:grid%i2) = C%uniform_Tsai2015_alpha_sq
      ice%beta_sq_a(   :,grid%i1:grid%i2) = C%uniform_Tsai2015_beta_sq
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ice%alpha_sq_a(  :,grid%i1:grid%i2) = C%uniform_Schoof2005_alpha_sq
      ice%beta_sq_a(   :,grid%i1:grid%i2) = C%uniform_Schoof2005_beta_sq
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ice%phi_fric_a(  :,grid%i1:grid%i2) = C%uniform_Coulomb_phi_fric
    ELSE
      IF (par%master) WRITE(0,*) 'initialise_basal_inversion - ERROR: unknown choice_sliding_law "', TRIM(C%choice_sliding_law), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! If needed, initialise target velocity fields
    ! ============================================
    
    IF     (C%choice_BIVgeo_method == 'PDC2012' .OR. &
            C%choice_BIVgeo_method == 'Lipscomb2021') THEN
      ! Not needed in these methods
    ELSEIF (C%choice_BIVgeo_method == 'CISM+') THEN
      ! Needed in these methods
      
      CALL initialise_basal_inversion_target_velocity( grid, ice)
      
    ELSE
      IF (par%master) WRITE(0,*) 'initialise_basal_inversion - ERROR: unknown choice_BIVgeo_method "', TRIM(C%choice_BIVgeo_method), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE initialise_basal_inversion
  SUBROUTINE initialise_basal_inversion_target_velocity( grid, ice)
    ! Initialise the target velocity fields used in a velocity-based basal inversion routine

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    TYPE(type_BIV_target_velocity)                     :: BIV_target
    INTEGER                                            :: i,j,i1,i2
    
    ! Determine filename
    IF (C%choice_BIVgeo_method == 'CISM+') THEN
      BIV_target%netcdf%filename = C%BIVgeo_CISMplus_target_filename
    ELSE
      IF (par%master) WRITE(0,*) 'initialise_basal_inversion_target_velocity - ERROR: unknown choice_BIVgeo_method "', TRIM(C%choice_BIVgeo_method), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    IF (par%master) WRITE(0,*) '  Initialising basal inversion target velocity from file ', TRIM( BIV_target%netcdf%filename), '...'
      
    ! Inquire grid data from the NetCDF file
    CALL allocate_shared_int_0D( BIV_target%nx, BIV_target%wnx)
    CALL allocate_shared_int_0D( BIV_target%ny, BIV_target%wny)
    
    IF (par%master) CALL inquire_BIV_target_velocity( BIV_target)
    CALL sync
    
    ! Allocate memory
    CALL allocate_shared_dp_1D( BIV_target%nx,                BIV_target%x        , BIV_target%wx        )
    CALL allocate_shared_dp_1D(                BIV_target%ny, BIV_target%y        , BIV_target%wy        )
    CALL allocate_shared_dp_2D( BIV_target%nx, BIV_target%ny, BIV_target%u_surf   , BIV_target%wu_surf   )
    CALL allocate_shared_dp_2D( BIV_target%nx, BIV_target%ny, BIV_target%v_surf   , BIV_target%wv_surf   )
    CALL allocate_shared_dp_2D( BIV_target%nx, BIV_target%ny, BIV_target%uabs_surf, BIV_target%wuabs_surf)
    
    IF (par%master) CALL read_BIV_target_velocity( BIV_target)
    CALL sync
  
    ! Safety
    CALL check_for_NaN_dp_2D( BIV_target%u_surf, 'BIV_target%u_surf', 'initialise_basal_inversion_target_velocity')
    CALL check_for_NaN_dp_2D( BIV_target%v_surf, 'BIV_target%v_surf', 'initialise_basal_inversion_target_velocity')
    
    ! Get absolute velocity
    CALL partition_list( BIV_target%nx, par%i, par%n, i1, i2)
    DO i = i1, i2
    DO j = 1, BIV_target%ny
      BIV_target%uabs_surf( i,j) = SQRT( BIV_target%u_surf( i,j)**2 + BIV_target%v_surf( i,j)**2)
    END DO
    END DO
    CALL sync
    
    ! Since we want data represented as [j,i] internally, transpose the data we just read.
    CALL transpose_dp_2D( BIV_target%uabs_surf, BIV_target%wuabs_surf)
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, ice%BIV_uabs_surf_target, ice%wBIV_uabs_surf_target)
    
    ! Map (transposed) raw data to the model grid
    CALL map_square_to_square_cons_2nd_order_2D( BIV_target%nx, BIV_target%ny, BIV_target%x, BIV_target%y, grid%nx, grid%ny, grid%x, grid%y, BIV_target%uabs_surf, ice%BIV_uabs_surf_target)
    
    ! Deallocate raw data
    CALL deallocate_shared( BIV_target%wnx       )
    CALL deallocate_shared( BIV_target%wny       )
    CALL deallocate_shared( BIV_target%wx        )
    CALL deallocate_shared( BIV_target%wy        )
    CALL deallocate_shared( BIV_target%wu_surf   )
    CALL deallocate_shared( BIV_target%wv_surf   )
    CALL deallocate_shared( BIV_target%wuabs_surf)
    
  END SUBROUTINE initialise_basal_inversion_target_velocity
  SUBROUTINE write_inverted_bed_roughness_to_file( grid, ice)
    ! Create a new NetCDF file and write the inverted bed roughness to it

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    CALL create_BIV_bed_roughness_file( grid, ice)
    
  END SUBROUTINE write_inverted_bed_roughness_to_file
  
END MODULE basal_conditions_and_sliding_module
