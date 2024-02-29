MODULE isotopes_module

  ! Contains all the routines for calculating the isotope content of the ice sheet.

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_grid, type_model_region
  USE parameters_module
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             surface_elevation
  USE forcing_module,                  ONLY: forcing
  USE netcdf_input_module,             ONLY: read_field_from_file_2D

  USE netcdf_debug_module,             ONLY: save_variable_as_netcdf_int_1D, save_variable_as_netcdf_int_2D, save_variable_as_netcdf_int_3D, &
                                             save_variable_as_netcdf_dp_1D,  save_variable_as_netcdf_dp_2D,  save_variable_as_netcdf_dp_3D



  IMPLICIT NONE

CONTAINS

! == Run the isotopes model (the main routine called from run_model_region)
  SUBROUTINE run_isotopes_model( region)
    ! Run the isotopes model (the main routine called from run_model_region)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_isotopes_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_isotopes_model == 'none') THEN
      ! Do nothing
    ELSEIF (C%choice_ice_isotopes_model == 'uniform') THEN
      ! Assign a uniform d18O value to all glacial ice
      region%ice%IsoIce( :,region%grid%i1:region%grid%i2) = C%uniform_ice_d18O
      CALL sync
    ELSEIF (C%choice_ice_isotopes_model == 'ANICE_legacy') THEN
      ! The ANICE_legacy englacial isotopes model
      CALL run_isotopes_model_ANICE_legacy( region)
    ELSE
      CALL crash('unknown choice_ice_isotopes_model "' // TRIM(C%choice_ice_isotopes_model) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_isotopes_model
  SUBROUTINE run_isotopes_model_ANICE_legacy( region)

    ! Based on the ANICE routines by Bas de Boer (November 2010).
    !
    ! using the implicit ice fluxes (vertical averaged) and applied mass balance fields from the
    ! ice thickness subroutine to calculate the advected isotope flux from all 4 directions.
    !
    ! Because of the dynamic shelf, we can also calculate IsoIce over the shelf area and threat
    ! all ice covered gridpoints the same since we are using vertically averaged velocities
    !
    ! IsoIce_adv should be zero for no ice (or at least have a small value, just like Hi)

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_isotopes_model_ANICE_legacy'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: Ts, Ts_ref, Hs, Hs_ref
    REAL(dp)                                           :: IsoMin,IsoMax    ! minimum and maximum value of IsoIce

    REAL(dp)                                           :: dIso_xl, dIso_xr, dIso_yd, dIso_yu, VIso
    REAL(dp), DIMENSION(:,:), POINTER                  ::  dIso_dt,  IsoIce_new
    REAL(dp)                                           :: left, right, up, down
    INTEGER                                            :: wdIso_dt, wIsoIce_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise
    region%ice%IsoSurf( :,region%grid%i1:region%grid%i2) = 0._dp

    ! Calculate the isotope content of annual mean precipitation
    ! ==========================================================

    IsoMax = -1E9_dp
    IsoMin =  1E9_dp

    region%ice%IsoSurf( :,region%grid%i1:region%grid%i2) = 0._dp

    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny

      IF (region%ice%mask_ice_a( j,i) == 1) THEN

        Ts     = SUM( region%climate%T2m( :,j,i)) / 12._dp
        Ts_ref = SUM( region%climate%matrix%PD_obs%T2m(  :,j,i)) / 12._dp
        Hs     = region%ice%Hs_a( j,i)
        Hs_ref = region%climate%matrix%PD_obs%Hs( j,i)

        region%ice%IsoSurf( j,i) = region%ice%IsoRef( j,i)                &
                                 + 0.35_dp              * (Ts - Ts_ref    &
                                 - C%constant_lapserate * (Hs - Hs_ref))  &
                                 - 0.0062_dp            * (Hs - Hs_ref)   ! from Clarke et al., 2005

        IsoMax = MAX( IsoMax, region%ice%IsoSurf( j,i))
        IsoMin = MIN( IsoMin, region%ice%IsoSurf( j,i))
      END IF

    END DO
    END DO
    CALL sync

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, IsoMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, IsoMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! Some regions may have undergone calving, but these regions contain no ice, but do contain Iso Ice.
    ! That should not cause problems, but still it is better to make sure IsoIce is 0. 
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
      IF (region%ice%mask_ice_a( j,i) == 0) THEN
        region%ice%IsoIce( j,i) = 0._dp ! Remove Iso Ice were there is no ice anymore
      END IF
    END DO
    END DO

    ! Calculate the mass gain/loss of d18O
    region%ice%MB_iso( :,region%grid%i1:region%grid%i2) = 0._dp

    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny

      IF ( region%SMB%SMB_year( j,i) > 0._dp ) THEN
        ! Surface accumulation has the isotope content of precipitation
        region%ice%MB_iso( j,i) = region%ice%MB_iso( j,i) + region%SMB%SMB_year( j,i) * region%ice%float_margin_frac_a( j,i) * region%ice%IsoSurf( j,i) * region%grid%dx * region%grid%dx    ! (applied MB) mass gain, so d18O from precipitation
      ELSE
        ! Surface melt has the isotope content of the ice itself
        region%ice%MB_iso( j,i) = region%ice%MB_iso( j,i) + region%SMB%SMB_year( j,i) * region%ice%float_margin_frac_a( j,i) * region%ice%IsoIce(  j,i) * region%grid%dx * region%grid%dx    ! (applied MB) mass loss, so d18O from ice
      END IF

      ! Both basal melt and basal freezing have the isotope content of the ice itself (the latter
      ! is not really true, but it's the best we can do for now)
      region%ice%MB_iso(   j,i) = region%ice%MB_iso( j,i) + region%BMB%BMB( j,i) * region%ice%float_margin_frac_a( j,i) * region%ice%IsoIce(  j,i) * region%grid%dx * region%grid%dx

      ! In case more ice is removed by melt than physically possible, set IsoIce to 0 and Mb to 0 just to be sure.
      ! Therefore, all Iso is removed. All ice has melted, but some may have been transported from neighbouring cells.
      ! To reflect this, set MB_iso to the remaining chunk of Iso
      IF (((region%BMB%BMB(     j,i) + region%SMB%SMB_year(            j,i))  * region%ice%float_margin_frac_a( j,i) * region%dt)  < -region%ice%Hi_a( j,i)) THEN
           region%ice%IsoIce( j,i) = 0._dp ! Reset Iso Ice to be sure
           region%ice%MB_iso( j,i) = 0._dp ! After resetting Iso Ice there does not need to be a MB flux for this cell
      END IF
      
    END DO
    END DO
    CALL sync

    ! Calculate the new d18O_ice from the ice fluxes and applied mass balance
    ! =======================================================================

    CALL allocate_shared_dp_2D( region%grid%ny, region%grid%nx, IsoIce_new, wIsoIce_new)
    CALL allocate_shared_dp_2D( region%grid%ny, region%grid%nx, dIso_dt,    wdIso_dt   )

    IsoIce_new( :,region%grid%i1:region%grid%i2) = 0._dp
    dIso_dt(    :,region%grid%i1:region%grid%i2) = 0._dp

    DO i = MAX(2,region%grid%i1), MIN(region%grid%nx-1,region%grid%i2)
    DO j = 2, region%grid%ny-1
      IF (region%ice%mask_ice_a( j,i) == 1) THEN

        ! Calculate advective isotope fluxes
        IF (region%ice%U_vav_cx( j  ,i-1) > 0._dp) THEN
          dIso_xl = region%ice%U_vav_cx( j  ,i-1) * region%ice%Hi_a( j  ,i-1) * region%ice%IsoIce( j  ,i-1) * region%grid%dx
        ELSE
          dIso_xl = region%ice%U_vav_cx( j  ,i-1) * region%ice%Hi_a( j  ,i  ) * region%ice%IsoIce( j  ,i  ) * region%grid%dx
        END IF
      
        IF (region%ice%U_vav_cx( j  ,i  ) > 0._dp) THEN
          dIso_xr = region%ice%U_vav_cx( j  ,i  ) * region%ice%Hi_a( j  ,i  ) * region%ice%IsoIce( j  ,i  ) * region%grid%dx
        ELSE
          dIso_xr = region%ice%U_vav_cx( j  ,i  ) * region%ice%Hi_a( j  ,i+1) * region%ice%IsoIce( j  ,i+1) * region%grid%dx
        END IF

        IF (region%ice%V_vav_cy( j-1,i  ) > 0._dp) THEN
          dIso_yd = region%ice%V_vav_cy( j-1,i  ) * region%ice%Hi_a( j-1,i  ) * region%ice%IsoIce( j-1,i  ) * region%grid%dx
        ELSE
          dIso_yd = region%ice%V_vav_cy( j-1,i  ) * region%ice%Hi_a( j  ,i  ) * region%ice%IsoIce( j  ,i  ) * region%grid%dx
        END IF
        
        IF (region%ice%V_vav_cy( j  ,i  ) > 0._dp) THEN
          dIso_yu = region%ice%V_vav_cy( j  ,i  ) * region%ice%Hi_a( j  ,i  ) * region%ice%IsoIce( j  ,i  ) * region%grid%dx
        ELSE
          dIso_yu = region%ice%V_vav_cy( j  ,i  ) * region%ice%Hi_a( j+1,i  ) * region%ice%IsoIce( j+1,i  ) * region%grid%dx
        END IF

        ! Calculate total isotope balance
        dIso_dt( j,i) = region%ice%MB_iso( j,i) + dIso_xl - dIso_xr + dIso_yd - dIso_yu

        ! Update vertically averaged ice isotope content
        VIso = region%ice%IsoIce( j,i) * region%ice%Hi_a( j,i) * region%grid%dx * region%grid%dx
        
        ! First, scale dIso_dt so more cannot leave than 
        VIso = VIso + dIso_dt( j,i) * region%dt

        ! Calcualte the new IsoIce and limit it between IsoMin and IsoMax
        IsoIce_new( j,i) = MIN( IsoMax, MAX( IsoMin, VIso / (region%grid%dx * region%grid%dx * region%ice%Hi_tplusdt_a( j,i)) ))
       
        ! FAIL SAVE: For very thin ice this calculation does not work. So assume IsoIce is 0.
        IF (region%ice%Hi_tplusdt_a( j,i) < 0.1_dp) THEN
          IsoIce_new( j,i) = 0._dp
        END IF

        ! FAIL SAVE: Calving front take the IsoIce from the source of advected ice
        IF (region%ice%float_margin_frac_a( j,i) < 0.99_dp) THEN
           
           ! Check the source areas:
           left  = 0._dp
           right = 0._dp
           up    = 0._dp
           down  = 0._dp

           ! Compare the ice that is being added to calculate a weighted average
           IF (region%ice%mask_ice_a( j,i-1) == 1) left   = ABS(MAX(0._dp,region%ice%U_vav_cx( j  ,i-1) * region%ice%Hi_a( j,  i-1)))
           IF (region%ice%mask_ice_a( j,i+1) == 1) right  = ABS(MIN(0._dp,region%ice%U_vav_cx( j  ,i  ) * region%ice%Hi_a( j,  i+1)))
           IF (region%ice%mask_ice_a( j+1,i) == 1) up     = ABS(MIN(0._dp,region%ice%V_vav_cy( j  ,i  ) * region%ice%Hi_a( j+1,i  )))
           IF (region%ice%mask_ice_a( j-1,i) == 1) down   = ABS(MAX(0._dp,region%ice%V_vav_cy( j-1,i  ) * region%ice%Hi_a( j-1,i  )))

           ! Calculate the weighted average of the ice sources
           IsoIce_new( j,i) = (region%ice%IsoIce( j,i-1) * left + region%ice%IsoIce( j,i+1) * right + &
                               region%ice%IsoIce( j+1,i) * up   + region%ice%IsoIce( j-1,i) * down) / &
                               (left + right + up + down)

           ! Only a very small amount of ice is being transported to the calving front. Use IsoSurf instead 
           IF ((left + right + up + down) < 0.001) THEN
             IsoIce_new( j,i) = region%ice%IsoSurf( j,i)
           END IF 

        END IF !  (region%ice%float_margin_frac_a( j,i) < 0.99_dp)
      END IF ! IF (ice%mask_ice_a( j,i) == 1) THEN
    END DO
    END DO
    CALL sync

    ! Update field
    region%ice%IsoIce( :,region%grid%i1:region%grid%i2) = IsoIce_new( :,region%grid%i1:region%grid%i2)

    ! Calculate mean isotope content of the whole ice sheet
    CALL calculate_isotope_content( region%grid, region%ice%Hi_tplusdt_a, region%ice%IsoIce, region%mean_isotope_content, region%d18O_contribution)

    ! Clean up after yourself
    CALL deallocate_shared( wdIso_dt   )
    CALL deallocate_shared( wIsoIce_new)

    ! Safety
    CALL check_for_NaN_dp_2D( region%ice%IsoSurf, 'region%ice%IsoSurf')
    CALL check_for_NaN_dp_2D( region%ice%MB_iso , 'region%ice%MB_iso' )
    CALL check_for_NaN_dp_2D( region%ice%IsoIce , 'region%ice%IsoIce' )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_isotopes_model_ANICE_legacy

! == Calculate the mean isotope content and benthic d18O contribution of this region's glacial ice
  SUBROUTINE calculate_isotope_content( grid, Hi, IsoIce, mean_isotope_content, d18O_contribution)
    ! Calculate mean isotope content of the whole ice sheet

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:),            INTENT(IN)    :: Hi
    REAL(dp), DIMENSION(:,:),            INTENT(IN)    :: IsoIce
    REAL(dp),                            INTENT(OUT)   :: mean_isotope_content
    REAL(dp),                            INTENT(OUT)   :: d18O_contribution

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calculate_isotope_content'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: Hi_msle
    REAL(dp)                                           :: total_isotope_content
    REAL(dp)                                           :: total_ice_volume_msle

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate total isotope content
    ! ===============================

    total_isotope_content = 0._dp
    total_ice_volume_msle = 0._dp

    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny

      IF (Hi( j,i) > 0._dp) THEN
        Hi_msle = Hi( j,i) * grid%dx * grid%dx * ice_density / (seawater_density * ocean_area)
        total_isotope_content = total_isotope_content + Hi_msle * IsoIce( j,i)
        total_ice_volume_msle = total_ice_volume_msle + Hi_msle
      END IF

    END DO
    END DO
    CALL sync

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_isotope_content, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_ice_volume_msle, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Weighted average of isotope content with Hi
    IF (par%master) THEN
      IF (total_ice_volume_msle > 0._dp) THEN
        mean_isotope_content  = total_isotope_content / total_ice_volume_msle
      ELSE
        mean_isotope_content  = 0._dp
      END IF
    END IF
    CALL sync

    ! Contribution to benthic d18O
    d18O_contribution = -1._dp * mean_isotope_content * total_ice_volume_msle / mean_ocean_depth

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_isotope_content

! == Initialise the isotopes model (allocating shared memory)
  SUBROUTINE initialise_isotopes_model( region)
    ! Allocate memory for the data fields of the isotopes model.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_isotopes_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_ice_isotopes_model == 'none') THEN
      ! Do nothing
    ELSEIF (C%choice_ice_isotopes_model == 'uniform') THEN
      ! Assign a uniform d18O value to all glacial ice

      ! Allocate shared memory
      CALL allocate_shared_dp_2D( region%grid%ny, region%grid%nx, region%ice%IsoIce, region%ice%wIsoIce)

      ! Assign value
      region%ice%IsoIce( :,region%grid%i1:region%grid%i2) = C%uniform_ice_d18O
      CALL sync

    ELSEIF (C%choice_ice_isotopes_model == 'ANICE_legacy') THEN
      ! The ANICE_legacy englacial isotopes model

      CALL initialise_isotopes_model_ANICE_legacy( region)

    ELSE
      CALL crash('unknown choice_ice_isotopes_model "' // TRIM(C%choice_ice_isotopes_model) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_isotopes_model
  SUBROUTINE initialise_isotopes_model_ANICE_legacy( region)
    ! Allocate memory for the data fields of the isotopes model.

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_isotopes_model_ANICE_legacy'
    INTEGER                                            :: i,j
    REAL(dp)                                           :: Ts, Ts_ref, Hs, Hs_ref
    CHARACTER(LEN=256)                                 :: filename
    CHARACTER(LEN=256)                                 :: choice_d18O_inverse_init
    REAL(dp)                                           :: time_to_restart_from

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE (0,*) '  Initialising ANICE_legacy isotopes model ...'

    IF (.NOT.(C%choice_climate_model == 'matrix')) THEN
      CALL crash('Can only be called if choice_climate_model_config == matrix for now') 
    END IF

    ! Determine which type of d18O initialization, and which filename
    IF  (C%do_NAM) THEN
      choice_d18O_inverse_init = C%choice_d18O_inverse_init_NAM
      filename                 = C%filename_d18O_inverse_init_NAM
      time_to_restart_from     = C%time_to_restart_from_NAM
    ELSEIF (C%do_EAS) THEN
      choice_d18O_inverse_init = C%choice_d18O_inverse_init_EAS
      filename                 = C%filename_d18O_inverse_init_EAS
      time_to_restart_from     = C%time_to_restart_from_EAS
    ELSEIF (C%do_GRL) THEN
      choice_d18O_inverse_init = C%choice_d18O_inverse_init_GRL
      filename                 = C%filename_d18O_inverse_init_GRL
      time_to_restart_from     = C%time_to_restart_from_GRL
    ELSEIF (C%do_ANT) THEN
      choice_d18O_inverse_init = C%choice_d18O_inverse_init_ANT
      filename                 = C%filename_d18O_inverse_init_ANT
      time_to_restart_from     = C%time_to_restart_from_ANT
    END IF

    ! Allocate memory
    CALL allocate_shared_dp_2D( region%grid%ny, region%grid%nx, region%ice%IsoRef    , region%ice%wIsoRef    )
    CALL allocate_shared_dp_2D( region%grid%ny, region%grid%nx, region%ice%IsoSurf   , region%ice%wIsoSurf   )
    CALL allocate_shared_dp_2D( region%grid%ny, region%grid%nx, region%ice%IsoIce    , region%ice%wIsoIce    )
    CALL allocate_shared_dp_2D( region%grid%ny, region%grid%nx, region%ice%MB_iso    , region%ice%wMB_iso    )

    ! Calculate reference field of d18O of precipitation
    ! (Zwally, H. J. and Giovinetto, M. B.: Areal distribution of the oxygen-isotope ratio in Greenland, Annals of Glaciology 25, 208-213, 1997)
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
      region%ice%IsoRef( j,i) = 0.691_dp * SUM( region%climate%matrix%PD_obs%T2m( :,j,i) / 12._dp) - 202.172_dp
    END DO
    END DO
    CALL sync

    ! ===== Present-day =====
    ! =======================

    ! Initialise ice sheet isotope content with the isotope content of annual mean precipitation
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny

      IF (region%refgeo_PD%Hi( j,i) > 0._dp) THEN

        Ts     = SUM( region%climate%matrix%PD_obs%T2m( :,j,i)) / 12._dp
        Ts_ref = SUM( region%climate%matrix%PD_obs%T2m( :,j,i)) / 12._dp
        Hs     = surface_elevation( region%refgeo_PD%Hi( j,i), region%refgeo_PD%Hb( j,i), 0._dp)
        Hs_ref = region%climate%matrix%PD_obs%Hs( j,i)

        region%ice%IsoIce( j,i) = region%ice%IsoRef( j,i)                &
                                + 0.35_dp              * (Ts - Ts_ref    &
                                - C%constant_lapserate * (Hs - Hs_ref))  &
                                - 0.0062_dp            * (Hs - Hs_ref)   ! from Clarke et al., 2005

      ELSE
        region%ice%IsoIce( j,i) = 0._dp ! = No ice
      END IF
    
    END DO
    END DO
    CALL sync

    ! Calculate mean isotope content of the whole ice sheet at present-day
    CALL calculate_isotope_content( region%grid, region%refgeo_PD%Hi, region%ice%IsoIce, region%mean_isotope_content_PD, region%d18O_contribution_PD)

    ! ===== Initial =====
    ! ===================
    IF (choice_d18O_inverse_init == 'init') THEN

      ! Initialise ice sheet isotope content with the isotope content of annual mean precipitation
      DO i = region%grid%i1, region%grid%i2
      DO j = 1, region%grid%ny

        IF (region%ice%mask_ice_a( j,i) == 1) THEN

          Ts     = SUM( region%climate%T2m( :,j,i)) / 12._dp
          Ts_ref = SUM( region%climate%matrix%PD_obs%T2m(  :,j,i)) / 12._dp
          Hs     = region%ice%Hs_a( j,i)
          Hs_ref = region%climate%matrix%PD_obs%Hs( j,i)

          region%ice%IsoIce( j,i) = region%ice%IsoRef( j,i)                &
                                  + 0.35_dp              * (Ts - Ts_ref    &
                                  - C%constant_lapserate * (Hs - Hs_ref))  &
                                  - 0.0062_dp            * (Hs - Hs_ref)   ! from Clarke et al., 2005
                                  

        ELSE
          region%ice%IsoIce( j,i) = 0._dp ! = No ice
        END IF
 
      END DO
      END DO
      CALL sync
     
    ELSEIF (choice_d18O_inverse_init == 'restart') THEN
      ! Read data from input file
      CALL read_field_from_file_2D(         filename, 'IsoIce', region%grid, region%ice%IsoIce, 'N/A', time_to_restart_from)
    
    ELSE
       CALL crash('unknown choice_forcing_method "' // TRIM(choice_d18O_inverse_init) // '"!')
    END IF

    ! Calculate mean isotope content of the whole ice sheet at the start of the simulation
    CALL calculate_isotope_content( region%grid, region%ice%Hi_a, region%ice%IsoIce, region%mean_isotope_content, region%d18O_contribution)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_isotopes_model_ANICE_legacy

END MODULE isotopes_module
