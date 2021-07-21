PROGRAM IMAU_ICE_program
  ! IMAU-ICE: the ice model developed by the Institute for Marine and Atmospheric research Utrecht (IMAU)
  !
  ! e-mail: c.j.berends@uu.nl
  !
  ! After some initialisation work (starting the program on multiple cores using MPI_INIT,
  ! reading the config file, creating an output folder, etc.), this program runs the four
  ! copies of the ice-sheet model (North America, Eurasia, Greenland and Antarctica).
  ! These are all run individually for 100 years (the "coupling interval"), after which
  ! control is passed back to this program. At this point, the sea-level model SELEN can be
  ! called, some global output data is calculated and written to the output file, and the
  ! coupling loop is run again.
  !
  ! The four ice-sheet models are four instances of the "model_region" data type (declared in
  ! the data_types module), which is accepted as an argument by the "run_model" subroutine.
  ! Model data is arranged into several large structures, all of which are declared in the 
  ! data_types module, and used by the different model subroutines. This prevents dependency
  ! problems during compiling, and makes the modules very clean and easy to read.

  USE mpi
  USE parallel_module,                 ONLY: par, sync, cerr, ierr
  USE configuration_module,            ONLY: dp, C, initialise_model_configuration
  USE data_types_module,               ONLY: type_model_region, type_climate_matrix, type_SELEN_global
  USE petsc_module,                    ONLY: initialise_petsc, finalise_petsc
  USE forcing_module,                  ONLY: forcing, initialise_insolation_data, update_insolation_data, initialise_CO2_record, update_CO2_at_model_time, &
                                             initialise_d18O_record, update_d18O_at_model_time, initialise_d18O_data, update_global_mean_temperature_change_history, &
                                             calculate_modelled_d18O, initialise_inverse_routine_data, inverse_routine_global_temperature_offset, inverse_routine_CO2, &
                                             initialise_geothermal_heat_flux
  USE climate_module,                  ONLY: initialise_climate_matrix
  USE global_text_output_module,       ONLY: create_text_output_file, write_text_output
  USE derivatives_and_grids_module,    ONLY: initialise_zeta_discretisation
  USE IMAU_ICE_main_model,             ONLY: initialise_model, run_model
  USE SELEN_main_module,               ONLY: initialise_SELEN, run_SELEN

  IMPLICIT NONE
  
  CHARACTER(LEN=256), PARAMETER          :: version_number = '2.0_dev'
  
  INTEGER                                :: process_rank, number_of_processes
  
  ! The four model regions
  TYPE(type_model_region)                :: NAM, EAS, GRL, ANT
  
  ! The global climate matrix
  TYPE(type_climate_matrix)              :: matrix
  
  ! SELEN
  TYPE(type_SELEN_global)                :: SELEN
  REAL(dp)                               :: ocean_area
  REAL(dp)                               :: ocean_depth
  
  ! Coupling timer and global mean sea level calculation
  REAL(dp)                               :: t_coupling, t_end_models
  REAL(dp)                               :: GMSL_NAM, GMSL_EAS, GMSL_GRL, GMSL_ANT, GMSL_glob
  
  ! Computation time tracking
  REAL(dp)                               :: tstart, tstop, dt
  INTEGER                                :: nr, ns, nm, nh, nd
  
  ! ======================================================================================
  
  ! MPI Initialisation
  ! ==================
  
  ! Use MPI to create copies of the program on all the processors, so the model can run in parallel.
  CALL MPI_INIT(ierr)
  
  ! Get rank of current process and total number of processes
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, process_rank, ierr)
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, number_of_processes, ierr)

  par%i      = process_rank
  par%n      = number_of_processes  
  par%master = (par%i == 0)
  
  IF (par%master) WRITE(0,*) ''
  IF (par%master) WRITE(0,*) '===================================================='
  IF (par%master) WRITE(0,'(A,A,A,I3,A)') ' ===== Running IMAU-ICE v', TRIM(version_number), ' on ', number_of_processes, ' cores ====='
  IF (par%master) WRITE(0,*) '===================================================='
  
  tstart = MPI_WTIME()
  
  ! PETSc Initialisation
  ! ====================
  
  ! Basically just a call to PetscInitialize
  CALL initialise_petsc
    
  ! Set up the model configuration from the provided config file(s) and create an output directory
  ! ==============================================================================================
  
  CALL initialise_model_configuration( version_number)
    
  ! ===== Initialise parameters for the vertical scaled coordinate transformation =====
  ! (the same for all model regions, so stored in the "C" structure)
  ! ===================================================================================
  
  CALL initialise_zeta_discretisation
  
  ! ===== Initialise forcing data =====
  ! ===================================
  
  CALL initialise_d18O_data
  CALL initialise_insolation_data
  CALL initialise_CO2_record
  CALL initialise_d18O_record
  CALL initialise_inverse_routine_data
  CALL initialise_geothermal_heat_flux
  
  ! ===== Initialise the climate matrix =====
  ! =========================================
  
  CALL initialise_climate_matrix( matrix)

  ! ===== Initialise the model regions ======
  ! =========================================
  
  IF (C%do_NAM) CALL initialise_model( NAM, 'NAM', matrix)
  IF (C%do_EAS) CALL initialise_model( EAS, 'EAS', matrix)
  IF (C%do_GRL) CALL initialise_model( GRL, 'GRL', matrix)
  IF (C%do_ANT) CALL initialise_model( ANT, 'ANT', matrix)
    
  ! Determine GMSL contributions of all simulated ice sheets
  GMSL_NAM = 0._dp
  GMSL_EAS = 0._dp
  GMSL_GRL = 0._dp
  GMSL_ANT = 0._dp
  IF (C%do_NAM) GMSL_NAM = NAM%GMSL_contribution
  IF (C%do_EAS) GMSL_EAS = EAS%GMSL_contribution
  IF (C%do_GRL) GMSL_GRL = GRL%GMSL_contribution
  IF (C%do_ANT) GMSL_ANT = ANT%GMSL_contribution
  GMSL_glob = GMSL_NAM + GMSL_EAS + GMSL_GRL + GMSL_ANT 
  
  ! Determine d18O contributions of all simulated ice sheets
  CALL update_global_mean_temperature_change_history( NAM, EAS, GRL, ANT)
  CALL calculate_modelled_d18O( NAM, EAS, GRL, ANT)
  
  ! ===== Initialise SELEN =====
  ! ============================
  
  IF (C%choice_GIA_model == 'SELEN' .OR. C%choice_sealevel_model == 'SELEN') THEN
    CALL initialise_SELEN( SELEN, NAM, EAS, GRL, ANT, version_number)
  END IF
    
  ! Timers and switch
  IF (C%SELEN_run_at_t_start) THEN
    SELEN%t0_SLE = C%start_time_of_run - C%dt_SELEN
    SELEN%t1_SLE = C%start_time_of_run
  ELSE
    SELEN%t0_SLE = C%start_time_of_run
    SELEN%t1_SLE = C%start_time_of_run + C%dt_SELEN
  END IF
    
  ! Create the global text output file
  IF (par%master) CALL create_text_output_file
  
  ! Write global data at t=0 to output file
  IF (par%master) CALL write_text_output( &
                          C%start_time_of_run,               &  ! time
                          GMSL_glob,                         &  ! global mean sea level
                          forcing%CO2_obs,                   &  ! observed CO2  from prescribed record (if any)
                          forcing%CO2_mod,                   &  ! modelled CO2                         (if any)
                          forcing%d18O_obs,                  &  ! observed d18O from prescribed record (if any)
                          forcing%d18O_mod,                  &  ! modelled d18O                        (always)
                          forcing%d18O_from_ice_volume_mod,  &  ! contribution to modelled d18O from ice volume
                          forcing%d18O_from_temperature_mod, &  !     ""            ""          ""   deep-sea temperature change
                          GMSL_NAM,                          &  ! contribution to GMSL from North America
                          GMSL_EAS,                          &  ! contribution to GMSL from Eurasia
                          GMSL_GRL,                          &  ! contribution to GMSL from Greenland
                          GMSL_ANT,                          &  ! contribution to GMSL from Antarctica
                          forcing%d18O_NAM,                  &  ! mean isotope content of North America
                          forcing%d18O_EAS,                  &  ! mean isotope content of Eurasia
                          forcing%d18O_GRL,                  &  ! mean isotope content of Greenland
                          forcing%d18O_ANT,                  &  ! mean isotope content of Antarctica
                          forcing%dT_glob,                   &  ! global mean surface temperature anomaly
                          forcing%dT_deepwater               )  ! deep-water temperature anomaly
  
! =============================
! ===== The big time loop =====
! =============================

  t_coupling = C%start_time_of_run
  
  DO WHILE (t_coupling < C%end_time_of_run)
  
    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,'(A,F9.3,A)') ' Coupling model: t = ', t_coupling/1000._dp, ' kyr'
    
    ! Solve the SLE
    IF (t_coupling >= SELEN%t1_SLE .AND. (C%choice_GIA_model == 'SELEN' .OR. C%choice_sealevel_model == 'SELEN')) THEN
      CALL run_SELEN( SELEN, NAM, EAS, GRL, ANT, t_coupling, ocean_area, ocean_depth)
      SELEN%t0_SLE = t_coupling
      SELEN%t1_SLE = t_coupling + C%dt_SELEN
    END IF
  
    ! Update global insolation forcing, CO2, and d18O at the current model time
    CALL update_insolation_data(    t_coupling)
    CALL update_CO2_at_model_time(  t_coupling)
    CALL update_d18O_at_model_time( t_coupling)
    
    ! Update regional sea level (needs to be moved to separate subroutine at some point!)
    IF (C%choice_sealevel_model == 'fixed') THEN
      IF (C%do_NAM) NAM%ice%SL_a( :,NAM%grid%i1:NAM%grid%i2) = C%fixed_sealevel
      IF (C%do_EAS) EAS%ice%SL_a( :,EAS%grid%i1:EAS%grid%i2) = C%fixed_sealevel
      IF (C%do_GRL) GRL%ice%SL_a( :,GRL%grid%i1:GRL%grid%i2) = C%fixed_sealevel
      IF (C%do_ANT) ANT%ice%SL_a( :,ANT%grid%i1:ANT%grid%i2) = C%fixed_sealevel
    ELSEIF (C%choice_sealevel_model == 'eustatic') THEN
      IF (C%do_NAM) NAM%ice%SL_a( :,NAM%grid%i1:NAM%grid%i2) = GMSL_glob
      IF (C%do_EAS) EAS%ice%SL_a( :,EAS%grid%i1:EAS%grid%i2) = GMSL_glob
      IF (C%do_GRL) GRL%ice%SL_a( :,GRL%grid%i1:GRL%grid%i2) = GMSL_glob
      IF (C%do_ANT) ANT%ice%SL_a( :,ANT%grid%i1:ANT%grid%i2) = GMSL_glob
    ELSEIF (C%choice_sealevel_model == 'SELEN') THEN
      ! Sea level fields are filled in the SELEN routines
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_sealevel_model "', TRIM(C%choice_sealevel_model), '" not implemented in IMAU_ICE_program!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Run all four model regions for 100 years
    t_end_models = MIN(C%end_time_of_run, t_coupling + C%dt_coupling)
    
    IF (C%do_NAM) CALL run_model( NAM, t_end_models)
    IF (C%do_EAS) CALL run_model( EAS, t_end_models)
    IF (C%do_GRL) CALL run_model( GRL, t_end_models)
    IF (C%do_ANT) CALL run_model( ANT, t_end_models)
    
    ! Advance coupling time
    t_coupling = t_end_models
    
    ! Determine GMSL contributions of all simulated ice sheets
    GMSL_NAM = 0._dp
    GMSL_EAS = 0._dp
    GMSL_GRL = 0._dp
    GMSL_ANT = 0._dp
    IF (C%do_NAM) GMSL_NAM = NAM%GMSL_contribution
    IF (C%do_EAS) GMSL_EAS = EAS%GMSL_contribution
    IF (C%do_GRL) GMSL_GRL = GRL%GMSL_contribution
    IF (C%do_ANT) GMSL_ANT = ANT%GMSL_contribution
    GMSL_glob = GMSL_NAM + GMSL_EAS + GMSL_GRL + GMSL_ANT
  
    ! Calculate contributions to global mean sea level and benthic d18O from the different ice sheets
    CALL update_global_mean_temperature_change_history( NAM, EAS, GRL, ANT)
    CALL calculate_modelled_d18O( NAM, EAS, GRL, ANT)
    
    ! If applicable, call the inverse routine to update the climate forcing parameter
    IF     (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      CALL inverse_routine_global_temperature_offset
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CALL inverse_routine_CO2
    ELSEIF (C%choice_forcing_method == 'CO2_direct') THEN
      ! No inverse routine is used in these forcing methods
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in IMAU_ICE_program!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Write global data to output file
    IF (par%master) CALL write_text_output( &
                            t_coupling,                        &  ! time
                            GMSL_glob,                         &  ! global mean sea level
                            forcing%CO2_obs,                   &  ! observed CO2  from prescribed record (if any)
                            forcing%CO2_mod,                   &  ! modelled CO2                         (if any)
                            forcing%d18O_obs,                  &  ! observed d18O from prescribed record (if any)
                            forcing%d18O_mod,                  &  ! modelled d18O                        (always)
                            forcing%d18O_from_ice_volume_mod,  &  ! contribution to modelled d18O from ice volume
                            forcing%d18O_from_temperature_mod, &  !     ""            ""          ""   deep-sea temperature change
                            GMSL_NAM,                          &  ! contribution to GMSL from North America
                            GMSL_EAS,                          &  ! contribution to GMSL from Eurasia
                            GMSL_GRL,                          &  ! contribution to GMSL from Greenland
                            GMSL_ANT,                          &  ! contribution to GMSL from Antarctica
                            forcing%d18O_NAM,                  &  ! mean isotope content of North America
                            forcing%d18O_EAS,                  &  ! mean isotope content of Eurasia
                            forcing%d18O_GRL,                  &  ! mean isotope content of Greenland
                            forcing%d18O_ANT,                  &  ! mean isotope content of Antarctica
                            forcing%dT_glob,                   &  ! global mean surface temperature anomaly
                            forcing%dT_deepwater               )  ! deep-water temperature anomaly
      
  END DO ! DO WHILE (t_coupling < C%end_time_of_run)
  
! ====================================
! ===== End of the big time loop =====
! ====================================
  
  ! Write total elapsed time to screen
  ! ==================================
  
  tstop = MPI_WTIME()
  
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
  
  IF (par%master) WRITE(0,'(A)') ''
  IF (par%master) WRITE(0,'(A)') ' ================================================================================'
  IF (par%master) WRITE(0,'(A,I2,A,I2,A,I2,A,I2,A)') ' ===== Simulation finished in ', nd, ' days, ', nh, ' hours, ', nm, ' minutes and ', ns, ' seconds! ====='
  IF (par%master) WRITE(0,'(A)') ' ================================================================================'
  IF (par%master) WRITE(0,'(A)') ''
  
  ! Finalise MPI and PETSc
  CALL finalise_petsc
  CALL MPI_FINALIZE( ierr) 
    
END PROGRAM IMAU_ICE_program
