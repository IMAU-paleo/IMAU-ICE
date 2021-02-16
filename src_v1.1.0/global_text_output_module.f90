MODULE global_text_output_module

  USE mpi
  USE parallel_module,                 ONLY: par, sync
  USE configuration_module,            ONLY: dp, C

  IMPLICIT NONE

CONTAINS
  
  SUBROUTINE create_text_output_file
    
    IMPLICIT NONE
    
    CHARACTER(LEN=256)                            :: filename
    
    filename = TRIM(C%output_dir) // 'aa_general_output.txt'
    OPEN(UNIT  = 1337, FILE = filename, STATUS = 'NEW')
    
    WRITE(UNIT = 1337, FMT = '(A)') '% IMAU-ICE global output data'
    WRITE(UNIT = 1337, FMT = '(A)') '%'
    WRITE(UNIT = 1337, FMT = '(A)') '% Time     : in yr, so LGM occurs at -21000'
    WRITE(UNIT = 1337, FMT = '(A)') '% sealevel : global mean sea level in m w.r.t. PD, so a sea-level drop shows up as a negative number'
    WRITE(UNIT = 1337, FMT = '(A)') '% CO2_obs  : observed CO2          from a prescribed record (set to zero if no record is prescribed)'
    WRITE(UNIT = 1337, FMT = '(A)') '% CO2_mod  : modelled CO2          from the inverse routine (set to zero if the inverse routine is not used)'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_obs : observed benthic d18O from a prescribed record (set to zero if no record is prescribed)'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_mod : modelled benthic d18O from the inverse routine (set to zero if the inverse routine is not used)'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_ice : contribution to benthic d18O from global ice volume'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_Tdw : contribution to benthic d18O from deep-water temperature change'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_NAM   : global mean sea level contribution from the North American ice sheet (= ice volume above flotation / ocean area)'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_EAS   : global mean sea level contribution from the Eurasian       ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_GRL   : global mean sea level contribution from the Greenland      ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% SL_ANT   : global mean sea level contribution from the Antarctic      ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_NAM : contribution to benthic d18O       from the North American ice sheet (= mean isotope content * sea-level equivalent ice volume / mean ocean depth)'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_EAS : contribution to benthic d18O       from the Eurasian       ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_GRL : contribution to benthic d18O       from the Greenland      ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% d18O_ANT : contribution to benthic d18O       from the Antarctic      ice sheet'
    WRITE(UNIT = 1337, FMT = '(A)') '% dT_glob  : global mean annual surface temperature change (scaled to sea-level)'
    WRITE(UNIT = 1337, FMT = '(A)') '% dT_dw    : deep-water temperature anomaly'
    WRITE(UNIT = 1337, FMT = '(A)') ''
    WRITE(UNIT = 1337, FMT = '(A,A)') '     Time     sealevel    CO2_obs    CO2_mod   d18O_obs   d18O_mod   d18O_ice   d18O_Tdw     ', &
                                      'SL_NAM     SL_EAS     SL_GRL     SL_ANT   d18O_NAM   d18O_EAS   d18O_GRL   d18O_ANT    dT_glob      dT_dw'
    
    CLOSE(UNIT = 1337)
    
  END SUBROUTINE create_text_output_file
  SUBROUTINE write_text_output( time, SL_glob, CO2_obs, CO2_mod, d18O_obs, d18O_mod, d18O_ice, d18O_Tdw, &
    SL_NAM, SL_EAS, SL_GRL, SL_ANT, d18O_NAM, d18O_EAS, d18O_GRL, d18O_ANT, dT_glob, dT_deepwater)
    ! Write data to global output file
  
    IMPLICIT NONE  
    
    ! In/output variables:
    REAL(dp),                   INTENT(IN)        :: time, SL_glob, CO2_obs, CO2_mod, d18O_obs, d18O_mod, d18O_ice, d18O_Tdw
    REAL(dp),                   INTENT(IN)        :: SL_NAM, SL_EAS, SL_GRL, SL_ANT
    REAL(dp),                   INTENT(IN)        :: d18O_NAM, d18O_EAS, d18O_GRL, d18O_ANT
    REAL(dp),                   INTENT(IN)        :: dT_glob, dT_deepwater
    
    ! Local variables:
    CHARACTER(LEN=256)                            :: filename
    
!    WRITE(0,*) ' time         = ', time
!    WRITE(0,*) ' SL_glob      = ', SL_glob
!    WRITE(0,*) ' CO2_obs      = ', CO2_obs
!    WRITE(0,*) ' CO2_mod      = ', CO2_mod
!    WRITE(0,*) ' d18O_obs     = ', d18O_obs
!    WRITE(0,*) ' d18O_mod     = ', d18O_mod
!    WRITE(0,*) ' d18O_ice     = ', d18O_ice
!    WRITE(0,*) ' d18O_Tdw     = ', d18O_Tdw
!    WRITE(0,*) ' SL_NAM       = ', SL_NAM
!    WRITE(0,*) ' SL_EAS       = ', SL_EAS
!    WRITE(0,*) ' SL_GRL       = ', SL_GRL
!    WRITE(0,*) ' SL_ANT       = ', SL_ANT
!    WRITE(0,*) ' d18O_NAM     = ', d18O_NAM
!    WRITE(0,*) ' d18O_EAS     = ', d18O_EAS
!    WRITE(0,*) ' d18O_GRL     = ', d18O_GRL
!    WRITE(0,*) ' d18O_ANT     = ', d18O_ANT
!    WRITE(0,*) ' dT_glob      = ', dT_glob
!    WRITE(0,*) ' dT_deepwater = ', dT_deepwater
            
    filename = TRIM(C%output_dir) // 'aa_general_output.txt'
    OPEN(UNIT  = 1337, FILE = filename, ACCESS = 'APPEND')
    
    WRITE(UNIT = 1337, FMT = '(18F11.2)') &
      time,         &  ! 1
      SL_glob,      &  ! 2
      CO2_obs,      &  ! 3
      CO2_mod,      &  ! 4
      d18O_obs,     &  ! 5
      d18O_mod,     &  ! 6
      d18O_ice,     &  ! 7
      d18O_Tdw,     &  ! 8
      SL_NAM,       &  ! 9
      SL_EAS,       &  ! 10
      SL_GRL,       &  ! 11
      SL_ANT,       &  ! 12
      d18O_NAM,     &  ! 13
      d18O_EAS,     &  ! 14
      d18O_GRL,     &  ! 15
      d18O_ANT,     &  ! 16
      dT_glob,      &  ! 17
      dT_deepwater     ! 18
    
    CLOSE(UNIT = 1337)
    
  END SUBROUTINE write_text_output

END MODULE global_text_output_module
