MODULE auto_tuner_module

  ! Contains all routines for the auto-tuner
  USE mpi
  USE configuration_module,            ONLY: dp, C           
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared, partition_list
  USE data_types_module,               ONLY: type_model_region, type_climate_matrix, type_auto_tuner, type_global_scalar_data, type_forcing_data
  USE netcdf_module,                   ONLY: create_restart_file, create_help_fields_file, write_to_restart_file, write_to_help_fields_file, &
                                             create_regional_scalar_output_file, create_global_scalar_output_file

  IMPLICIT NONE
  
CONTAINS


  ! The initialisation routine
  SUBROUTINE initialise_auto_tuner( auto_tuner, GMSL)
    ! Initialize the auto_tuner
    
    IMPLICIT NONE
    
    TYPE(type_auto_tuner), INTENT(INOUT)   :: auto_tuner
    REAL(dp),  INTENT(IN)                  :: GMSL

    ! Initialization
    CALL allocate_shared_dp_0D( auto_tuner%NAM_peak_SLC          , auto_tuner%wNAM_peak_SLC)
    CALL allocate_shared_dp_0D( auto_tuner%EAS_peak_SLC          , auto_tuner%wEAS_peak_SLC)
    CALL allocate_shared_dp_0D( auto_tuner%GRL_peak_SLC          , auto_tuner%wGRL_peak_SLC)
    CALL allocate_shared_dp_0D( auto_tuner%ANT_peak_SLC          , auto_tuner%wANT_peak_SLC)

    CALL allocate_shared_dp_0D( auto_tuner%Glob_peak_SLC         , auto_tuner%wGlob_peak_SLC)
    CALL allocate_shared_dp_0D( auto_tuner%GMSL_initial          , auto_tuner%wGMSL_initial)
    CALL allocate_shared_int_0D( auto_tuner%n_iteration           , auto_tuner%wn_iteration)

    ! Initialise values for the auto-tuner
    auto_tuner%NAM_peak_SLC  = 0._dp
    auto_tuner%EAS_peak_SLC  = 0._dp
    auto_tuner%GRL_peak_SLC  = 0._dp
    auto_tuner%ANT_peak_SLC  = 0._dp
    auto_tuner%Glob_peak_SLC = 0._dp

    ! Initialise the counter that keeps track of the current auto-tuner iteration
    auto_tuner%n_iteration   = 1._dp

    ! Storing the initial GMSL into the auto-tuner
    auto_tuner%GMSL_initial = GMSL

  END SUBROUTINE initialise_auto_tuner
  
  
  SUBROUTINE auto_tuner_main(NAM, EAS, GRL, ANT, auto_tuner, t_coupling, GMSL_NAM, GMSL_EAS, GMSL_GRL, GMSL_ANT, GMSL_Glob,forcing, global_data)
        ! This routine of the auto-tuner serves two purposes. First it keeps track of the higherst SLC (either positive or negative).
        ! Secondly, it checks if the stopping condition is met. This stopping condition is met either when the peak SLC has decreased by
        ! a certain fraction, or if the t_coupling is below C%at_latest_stopping_time.
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_model_region),                     INTENT(INOUT)       :: NAM
    TYPE(type_model_region),                     INTENT(INOUT)       :: EAS
    TYPE(type_model_region),                     INTENT(INOUT)       :: GRL
    TYPE(type_model_region),                     INTENT(INOUT)       :: ANT
    TYPE(type_auto_tuner),                       INTENT(INOUT)       :: auto_tuner
    TYPE(type_forcing_data),                     INTENT(INOUT)       :: forcing
    TYPE(type_global_scalar_data),               INTENT(INOUT)       :: global_data

    REAL(dp),                                    INTENT(INOUT)       :: GMSL_NAM
    REAL(dp),                                    INTENT(INOUT)       :: GMSL_EAS
    REAL(dp),                                    INTENT(INOUT)       :: GMSL_GRL
    REAL(dp),                                    INTENT(INOUT)       :: GMSL_ANT    
    REAL(dp),                                    INTENT(INOUT)       :: GMSL_Glob

    REAL(dp),                                    INTENT(INOUT)       :: t_coupling
 
    ! Local variables:
    LOGICAL                                            :: at_SLC_final_stopping_condition     ! Keeps track of the stopping condition for the full run
    INTEGER                                            :: i,j
    CHARACTER(LEN=1)                                   :: number_filename

    ! ===== Check for auto-tuner new iteration conditions ===== !
        
    ! Keep track of the peak volume
    IF (C%do_NAM .AND. (ABS(auto_tuner%NAM_peak_SLC) < ABS(GMSL_NAM))) auto_tuner%NAM_peak_SLC = GMSL_NAM
    IF (C%do_EAS .AND. (ABS(auto_tuner%EAS_peak_SLC) < ABS(GMSL_EAS))) auto_tuner%EAS_peak_SLC = GMSL_EAS
    IF (C%do_GRL .AND. (ABS(auto_tuner%GRL_peak_SLC) < ABS(GMSL_GRL))) auto_tuner%GRL_peak_SLC = GMSL_GRL
    IF (C%do_ANT .AND. (ABS(auto_tuner%ANT_peak_SLC) < ABS(GMSL_ANT))) auto_tuner%ANT_peak_SLC = GMSL_ANT
    IF (ABS(auto_tuner%Glob_peak_SLC) < ABS(GMSL_Glob)) auto_tuner%Glob_peak_SLC = GMSL_Glob
            
    ! Check for the stopping condition of the current iteration
    IF ((t_coupling >= C%at_earliest_stopping_time .AND.  (ABS(GMSL_glob) < ABS(auto_tuner%Glob_peak_SLC*C%SLC_for_stopping_condition))) .OR. (t_coupling > C%at_latest_stopping_time))     THEN
            ! New iteration condition is met if either: 1. The current model time is after the earliest stopping time AND the peak SLC is below the maximum with a certain percentage
            !  2. OR the current model time exceeds the C%at_SLC_iteration_stopping_condition and a stopping condition is forced
 
            ! Tell the user a new iteration has started
            IF (par%master) WRITE(0,*) '---------------------------------------------------------'
            IF (par%master) WRITE(0,'(A,I3,A,I3,A)') 'Auto-tuner iteration ', (auto_tuner%n_iteration),' finished, ', (C%n_at_iterations-1),' remaining'
            IF (par%master) WRITE(0,*) ' '
  
            ! ===== Check for IMAU-ICE program stopping conditions ===== !
            ! Condition for a new iteration has been met
            IF (par%master) C%n_at_iterations      = C%n_at_iterations - 1      ! Once it reaches 0, the program is forced to stop
            IF (par%master) auto_tuner%n_iteration = auto_tuner%n_iteration + 1 ! Keeping track of the current number of iteration

            ! Check if IMAU_ICE needs to stop
            at_SLC_final_stopping_condition = .TRUE.

            IF (C%n_at_iterations>0) THEN
                    IF (C%do_NAM .AND. (ABS(auto_tuner%NAM_peak_SLC - C%SLC_ref_NAM) > ABS(C%SLC_ref_NAM*(1._dp - C%SLC_for_final_stopping_condition)))) at_SLC_final_stopping_condition = .FALSE. ! Checks if NAM SLC decreased enough compared to the peak SLC, otherwise no stopping condition
                    IF (C%do_EAS .AND. (ABS(auto_tuner%EAS_peak_SLC - C%SLC_ref_EAS) > ABS(C%SLC_ref_EAS*(1._dp - C%SLC_for_final_stopping_condition)))) at_SLC_final_stopping_condition = .FALSE. ! Checks if EAS SLC decreased enough compared to the peak SLC, otherwise no stopping condition
                    IF (C%do_GRL .AND. (ABS(auto_tuner%GRL_peak_SLC - C%SLC_ref_GRL) > ABS(C%SLC_ref_GRL*(1._dp - C%SLC_for_final_stopping_condition)))) at_SLC_final_stopping_condition = .FALSE. ! Checks if GRL SLC decreased enough compared to the peak SLC, otherwise no stopping condition
                    IF (C%do_ANT .AND. (ABS(auto_tuner%ANT_peak_SLC - C%SLC_ref_ANT) > ABS(C%SLC_ref_ANT*(1._dp - C%SLC_for_final_stopping_condition)))) at_SLC_final_stopping_condition = .FALSE. ! Checks if ANT SLC decreased enough compared to the peak SLC, otherwise no stopping condition
            END IF

            ! Check if the stopping condition for IMAU-ICE has been met
            IF (at_SLC_final_stopping_condition) THEN ! Stop the Program, the run is finished or was forced to finish

                    IF (C%n_at_iterations<=0 .AND. par%master) THEN
                            WRITE(0,*) 'Auto-tuner finished all iterations'
                    ELSEIF (par%master) THEN
                            WRITE(0,*) 'Auto-tuner has finished early, found SLC close enough to the reference'
                    END IF
                
                    IF (par%master) WRITE(0,*) '---------------------------------------------------------'

                    ! The program stops either when the final stopping condition is met (e.g. the tuned volume is approached)
                    ! OR if the number of iterations given by the user is equal to the number of finalized iterations
                    t_coupling = C%end_time_of_run ! Stops the program


            ELSE ! Prepare for a new iteration

                    IF (par%master) WRITE(0,*) 'Preparing for new iteration...'
               
                    ! Get the new SMB constant
                    IF (par%master) THEN
                            IF (C%do_NAM) NAM%SMB%C_abl_constant = NAM%SMB%C_abl_constant + ((auto_tuner%NAM_peak_SLC - C%SLC_ref_NAM) * C%multiplication_factor_NAM)
                            IF (C%do_EAS) EAS%SMB%C_abl_constant = EAS%SMB%C_abl_constant + ((auto_tuner%EAS_peak_SLC - C%SLC_ref_EAS) * C%multiplication_factor_EAS)
                            IF (C%do_GRL) GRL%SMB%C_abl_constant = GRL%SMB%C_abl_constant + ((auto_tuner%GRL_peak_SLC - C%SLC_ref_GRL) * C%multiplication_factor_GRL)
                            IF (C%do_ANT) ANT%SMB%C_abl_constant = ANT%SMB%C_abl_constant + ((auto_tuner%ANT_peak_SLC - C%SLC_ref_ANT) * C%multiplication_factor_ANT)
                    END IF

                    ! Synchronise, this is to prevent the peak_SLC to be overwritten
                    ! before it can change the SMB
                    CALL sync

                    ! Reset the peak sea level contribution
                    auto_tuner%NAM_peak_SLC  = 0._dp
                    auto_tuner%EAS_peak_SLC  = 0._dp
                    auto_tuner%GRL_peak_SLC  = 0._dp
                    auto_tuner%ANT_peak_SLC  = 0._dp
                    auto_tuner%Glob_peak_SLC = 0._dp
            
                    IF (par%master) WRITE(0,*) 'Resetting all values to the initial conditions...'

                    ! Reset the 2D and 3D fields that need to be converted back to initial values, create new files and filenames
                    IF (C%do_NAM) CALL auto_tuner_reset_to_initial_values(NAM,t_coupling,auto_tuner,forcing, number_filename)
                    IF (C%do_EAS) CALL auto_tuner_reset_to_initial_values(EAS,t_coupling,auto_tuner,forcing, number_filename)
                    IF (C%do_GRL) CALL auto_tuner_reset_to_initial_values(GRL,t_coupling,auto_tuner,forcing, number_filename)
                    IF (C%do_ANT) CALL auto_tuner_reset_to_initial_values(ANT,t_coupling,auto_tuner,forcing, number_filename)

                    IF (par%master) CALL create_global_scalar_output_file( global_data%netcdf, number_filename)

                    IF (par%master) WRITE(0,*) '---------------------------------------------------------'
                    IF (par%master) WRITE(0,*) ' '
  
            END IF
    END IF

  END SUBROUTINE auto_tuner_main
  
  
  SUBROUTINE auto_tuner_reset_to_initial_values(region,t_coupling,auto_tuner,forcing,number_filename)
    ! auto_tuner_reset_to_initial_values: Initialize the auto_tuner variables
    ! This subroutine is called when the auto-tuner is tuned on, and the simulation has completed a 
    ! cycle. It resets all relavent variables to the initial conditions so a new cycle can be started.
    
    IMPLICIT NONE
   
    TYPE(type_model_region),                     INTENT(INOUT)       :: region
    TYPE(type_auto_tuner),                       INTENT(INOUT)       :: auto_tuner
    REAL(dp),                                    INTENT(INOUT)       :: t_coupling
    TYPE(type_forcing_data),                     INTENT(INOUT)       :: forcing
    CHARACTER(LEN=1),                            INTENT(INOUT)       :: number_filename
    
    ! Local variables:
    CHARACTER(LEN=20)                                  :: short_filename

    INTEGER                                            :: i,j,n
 
    ! Resetting topography to initial conditions
    DO i = region%grid%i1, region%grid%i2
    DO j = 1, region%grid%ny
     
        region%ice%Hi_a(        j,i) = region%init%Hi( j,i)
        region%ice%Hb_a(        j,i) = region%init%Hb( j,i)
        region%ice%Hs_a(        j,i) = MAX( region%ice%SL_a( j,i), (region%ice%Hb_a( j,i) + region%ice%Hi_a( j,i)))
    
    END DO
    END DO   
    
    CALL sync
      
    ! Time travel to the initial run-time
    t_coupling        = C%start_time_of_run ! Model Time
    region%time       = C%start_time_of_run ! Regional Time (which keeps track of the previous time-step)
   
    ! Time for each feature in IMAU-ICE that keeps track of a time-step
    region%t_next_SIA  = C%start_time_of_run
    region%t_next_SSA  = C%start_time_of_run
    region%t_next_DIVA = C%start_time_of_run

    region%t_next_thermo  = C%start_time_of_run
    region%t_next_climate = C%start_time_of_run
    region%t_next_SMB     = C%start_time_of_run
    region%t_next_BMB     = C%start_time_of_run
    region%t_next_ELRA    = C%start_time_of_run
    region%t_next_output  = C%start_time_of_run 

    ! Reset the sea level contribution
    region%ice%SL_a( :,region%grid%i1:region%grid%i2) = auto_tuner%GMSL_initial

    ! Get the number to fill into the file name
    IF (auto_tuner%n_iteration == 1) THEN
        number_filename = '1'
    ELSEIF  (auto_tuner%n_iteration == 2) THEN
        number_filename = '2'
    ELSEIF  (auto_tuner%n_iteration == 3) THEN
        number_filename = '3'
    ELSEIF  (auto_tuner%n_iteration == 4) THEN
        number_filename = '4'
    ELSEIF  (auto_tuner%n_iteration == 5) THEN
        number_filename = '5'
    ELSEIF  (auto_tuner%n_iteration == 6) THEN
        number_filename = '6'
    ELSEIF  (auto_tuner%n_iteration == 7) THEN
        number_filename = '7'
    ELSEIF  (auto_tuner%n_iteration == 8) THEN
        number_filename = '8'
    ELSEIF  (auto_tuner%n_iteration == 9) THEN
        number_filename = '9'
    END IF

    ! Set output filenames for this region
    short_filename = 'restart_NAM__.nc'
    short_filename(9:11)  = region%name
    short_filename(13:13) = number_filename
    
    DO n = 1, 256
      region%restart%filename(n:n) = ' '
    END DO
    region%restart%filename = TRIM(C%output_dir)//TRIM(short_filename)

    short_filename = 'help_fields_NAM__.nc'
    short_filename(13:15) = region%name

    ! At the moment this is a quick fix:
    short_filename(17:17) = number_filename
    
    DO n = 1, 256
      region%help_fields%filename(n:n) = ' '
    END DO
    
    region%help_fields%filename = TRIM(C%output_dir)//TRIM(short_filename)

    ! Create the new (empty) netcdf files
    IF (par%master)                             CALL create_restart_file( region, forcing)
    IF (par%master)                             CALL create_help_fields_file( region)
    ! IF (par%master .AND. C%do_write_debug_data) CALL create_debug_file( region)
    
    ! Create new output names for the scalar variables
    ! IF (par%master)                             CALL create_regional_scalar_output_file( region, number_filename)

  END SUBROUTINE auto_tuner_reset_to_initial_values
  
END MODULE auto_tuner_module
