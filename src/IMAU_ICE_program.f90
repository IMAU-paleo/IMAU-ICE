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
  USE configuration_module,            ONLY: dp, C, initialise_model_configuration, write_total_model_time_to_screen
  USE data_types_module,               ONLY: type_model_region, type_climate_matrix, type_SELEN_global, type_global_scalar_data
  USE petsc_module,                    ONLY: initialise_petsc, finalise_petsc
  USE forcing_module,                  ONLY: forcing, initialise_insolation_data, update_insolation_data, initialise_CO2_record, update_CO2_at_model_time, &
                                             initialise_d18O_record, update_d18O_at_model_time, initialise_d18O_data, update_global_mean_temperature_change_history, &
                                             calculate_modelled_d18O, initialise_inverse_routine_data, inverse_routine_global_temperature_offset, inverse_routine_CO2, &
                                             initialise_geothermal_heat_flux, initialise_climate_SMB_forcing_data, update_climate_SMB_forcing_data
  USE climate_module,                  ONLY: initialise_climate_matrix
  USE derivatives_and_grids_module,    ONLY: initialise_zeta_discretisation
  USE IMAU_ICE_main_model,             ONLY: initialise_model, run_model
  USE SELEN_main_module,               ONLY: initialise_SELEN, run_SELEN
  USE scalar_data_output_module,       ONLY: initialise_global_scalar_data, write_global_scalar_data
  USE general_ice_model_data_module,   ONLY: MISMIPplus_adapt_flow_factor

  IMPLICIT NONE
  
  CHARACTER(LEN=256), PARAMETER          :: version_number = '2.0'
  
  INTEGER                                :: process_rank, number_of_processes
  
  ! The four model regions
  TYPE(type_model_region)                :: NAM, EAS, GRL, ANT
  
  ! The global climate matrix
  TYPE(type_climate_matrix)              :: matrix
  
  ! SELEN
  TYPE(type_SELEN_global)                :: SELEN
  REAL(dp)                               :: ocean_area
  REAL(dp)                               :: ocean_depth
  
  ! Global scalar data (sea level, CO2, d18O, etc.)
  TYPE(type_global_scalar_data)          :: global_data
  
  ! Coupling timer
  REAL(dp)                               :: t_coupling, t_end_models
  
  ! Computation time tracking
  REAL(dp)                               :: tstart, tstop
  
  ! MISMIPplus flow factor tuning
  REAL(dp)                               :: Hprev, Hcur
  
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
    
  ! ===== Create the global scalar output file =====
  ! ================================================
  
  CALL initialise_global_scalar_data( global_data)
  
  ! ===== Initialise the climate matrix =====
  ! =========================================
  
  IF ((C%choice_forcing_method == 'SMB_direct') .OR. (C%choice_forcing_method == 'climate_direct')) THEN
    CALL initialise_climate_SMB_forcing_data
  ELSE
    CALL initialise_climate_matrix( matrix)
  END IF
    
  ! ===== Initialise the model regions ======
  ! =========================================
  
  IF (C%do_NAM) CALL initialise_model( NAM, 'NAM', matrix)
  IF (C%do_EAS) CALL initialise_model( EAS, 'EAS', matrix)
  IF (C%do_GRL) CALL initialise_model( GRL, 'GRL', matrix)
  IF (C%do_ANT) CALL initialise_model( ANT, 'ANT', matrix)
    
  ! Set GMSL contributions of all simulated ice sheets
  IF (par%master) THEN
    global_data%GMSL_NAM = 0._dp
    global_data%GMSL_EAS = 0._dp
    global_data%GMSL_GRL = 0._dp
    global_data%GMSL_ANT = 0._dp
    IF (C%do_NAM) global_data%GMSL_NAM = NAM%GMSL_contribution
    IF (C%do_EAS) global_data%GMSL_EAS = EAS%GMSL_contribution
    IF (C%do_GRL) global_data%GMSL_GRL = GRL%GMSL_contribution
    IF (C%do_ANT) global_data%GMSL_ANT = ANT%GMSL_contribution
  END IF ! IF (par%master) THEN
  CALL sync
  
  ! Determine global mean sea level
  IF     (C%choice_sealevel_model == 'fixed') THEN
    IF (par%master) global_data%GMSL = C%fixed_sealevel
  ELSEIF (C%choice_sealevel_model == 'eustatic' .OR. C%choice_sealevel_model == 'SELEN') THEN
    IF (par%master) global_data%GMSL = global_data%GMSL_NAM + global_data%GMSL_EAS + global_data%GMSL_GRL + global_data%GMSL_ANT 
  ELSE
    IF (par%master) WRITE(0,*) '  ERROR: choice_sealevel_model "', TRIM(C%choice_sealevel_model), '" not implemented in IMAU_ICE_program!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
  END IF
  CALL sync
  
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
    
  ! Write global scalar data at the start of the simulation to output file
  CALL write_global_scalar_data( global_data, NAM, EAS, GRL, ANT, forcing, C%start_time_of_run)
  
! =============================
! ===== The big time loop =====
! =============================

  t_coupling = C%start_time_of_run
  Hcur       = 0._dp
  
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
    IF ((C%choice_forcing_method == 'SMB_direct') .OR. (C%choice_forcing_method == 'climate_direct')) CALL update_climate_SMB_forcing_data (t_coupling)
    
    ! Update regional sea level
    IF (C%choice_sealevel_model == 'fixed' .OR. C%choice_sealevel_model == 'eustatic') THEN
      ! Local sea level is equal to the eustatic signal
      IF (C%do_NAM) NAM%ice%SL_a( :,NAM%grid%i1:NAM%grid%i2) = global_data%GMSL
      IF (C%do_EAS) EAS%ice%SL_a( :,EAS%grid%i1:EAS%grid%i2) = global_data%GMSL
      IF (C%do_GRL) GRL%ice%SL_a( :,GRL%grid%i1:GRL%grid%i2) = global_data%GMSL
      IF (C%do_ANT) ANT%ice%SL_a( :,ANT%grid%i1:ANT%grid%i2) = global_data%GMSL
      CALL sync
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
    
    ! Set GMSL contributions of all simulated ice sheets
    global_data%GMSL_NAM = 0._dp
    global_data%GMSL_EAS = 0._dp
    global_data%GMSL_GRL = 0._dp
    global_data%GMSL_ANT = 0._dp
    IF (C%do_NAM) global_data%GMSL_NAM = NAM%GMSL_contribution
    IF (C%do_EAS) global_data%GMSL_EAS = EAS%GMSL_contribution
    IF (C%do_GRL) global_data%GMSL_GRL = GRL%GMSL_contribution
    IF (C%do_ANT) global_data%GMSL_ANT = ANT%GMSL_contribution
    
    ! Determine global mean sea level
    IF     (C%choice_sealevel_model == 'fixed') THEN
      global_data%GMSL = C%fixed_sealevel
    ELSEIF (C%choice_sealevel_model == 'eustatic' .OR. C%choice_sealevel_model == 'SELEN') THEN
      global_data%GMSL = global_data%GMSL_NAM + global_data%GMSL_EAS + global_data%GMSL_GRL + global_data%GMSL_ANT 
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_sealevel_model "', TRIM(C%choice_sealevel_model), '" not implemented in IMAU_ICE_program!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
  
    ! Calculate contributions to global mean sea level and benthic d18O from the different ice sheets
    CALL update_global_mean_temperature_change_history( NAM, EAS, GRL, ANT)
    CALL calculate_modelled_d18O( NAM, EAS, GRL, ANT)
    
    ! If applicable, call the inverse routine to update the climate forcing parameter
    IF     (C%choice_forcing_method == 'd18O_inverse_dT_glob') THEN
      CALL inverse_routine_global_temperature_offset
    ELSEIF (C%choice_forcing_method == 'd18O_inverse_CO2') THEN
      CALL inverse_routine_CO2
    ELSEIF ((C%choice_forcing_method == 'CO2_direct') .OR. (C%choice_forcing_method == 'SMB_direct') .OR. (C%choice_forcing_method == 'climate_direct')) THEN
      ! No inverse routine is used in these forcing methods
    ELSE
      IF (par%master) WRITE(0,*) '  ERROR: choice_forcing_method "', TRIM(C%choice_forcing_method), '" not implemented in IMAU_ICE_program!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Write global scalar data to output file
    CALL write_global_scalar_data( global_data, NAM, EAS, GRL, ANT, forcing, t_coupling)
    
    ! MISMIP+ flow factor tuning for GL position
    IF (C%MISMIPplus_do_tune_A_for_GL) THEN
      Hprev = Hcur
      Hcur  = ANT%ice%Hs_a( CEILING( REAL(ANT%grid%ny,dp)/2._dp), 1)
      IF (par%master) WRITE(0,*) 'Hprev = ', Hprev, ', Hcur = ', Hcur
      IF (ABS(1._dp - Hcur/Hprev) < 5.0E-3_dp) THEN
        ! The model has converged to a steady state; adapt the flow factor
        CALL MISMIPplus_adapt_flow_factor( ANT%grid, ANT%ice)
      END IF
    END IF
      
  END DO ! DO WHILE (t_coupling < C%end_time_of_run)
  
! ====================================
! ===== End of the big time loop =====
! ==================================== 
  
  ! Write total elapsed time to screen
  tstop = MPI_WTIME()
  IF (par%master) CALL write_total_model_time_to_screen( tstart, tstop)
  CALL sync
  
  ! Finalise MPI and PETSc
  CALL finalise_petsc
  CALL MPI_FINALIZE( ierr)
    
END PROGRAM IMAU_ICE_program
