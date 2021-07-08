MODULE SELEN_mapping_module
! Routines for mapping data between the UFEMISM square grid and the SELEN grid
  
  USE mpi
  USE configuration_module,          ONLY: dp, C
  USE parameters_module
  USE parallel_module,               ONLY: par, sync, cerr, ierr, &
                                           allocate_shared_int_0D, allocate_shared_dp_0D, &
                                           allocate_shared_int_1D, allocate_shared_dp_1D, &
                                           allocate_shared_int_2D, allocate_shared_dp_2D, &
                                           allocate_shared_int_3D, allocate_shared_dp_3D, &
                                           deallocate_shared
  USE data_types_module,             ONLY: type_SELEN_global, type_model_region

  IMPLICIT NONE

CONTAINS
  SUBROUTINE create_GIA_grid_to_SELEN_maps( SELEN, region, region_label)
    ! Create mapping arrays between the (square) ice model GIA grid,
    ! and the SELEN global grid. SELEN global grid pixels that are
    ! already covered by another ice model region are skipped, to make
    ! sure that no areas are double-counted.
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(INOUT) :: SELEN
    TYPE(type_model_region),             INTENT(INOUT) :: region
    INTEGER,                             INTENT(IN)    :: region_label ! 1 = NAM, 2 = EAS, 3 = GRL, 4 = ANT
        
    INTEGER                                            :: nmax
    INTEGER                                            :: i,j,n
    REAL(dp)                                           :: lon_gr, lat_gr
    INTEGER                                            :: iapl,iapu,isll,islu,isl
    REAL(dp)                                           :: lon_sl, lat_sl
    REAL(dp)                                           :: deg, intarea, pixarea, w
    
    IF (par%master) WRITE(0,*) '   Creating icemodel-to-SELEN mapping arrays for region ', region%name, '...'
  
    ! Determine the radius of the circular discs representing this region's square grid pixels
    CALL allocate_shared_dp_0D( region%SELEN%rad, region%SELEN%wrad)
    IF (par%master) region%SELEN%rad = C%dx_GIA / (SQRT(pi) * earth_radius)
    CALL sync
    
    ! The area of a SELEN global grid pixel
    pixarea = (2._dp * pi * earth_radius**2 * (1._dp - COS( C%SELEN_alfa)))
    
    ! Estimate the maximum number of SELEN global grid pixels to which each square grid pixel can contribute,
    ! so we know how much memory to allocate
    nmax = 8 * MAX( CEILING( region%SELEN%rad / C%SELEN_alfa), CEILING( C%SELEN_alfa / region%SELEN%rad))
    
    ! Allocate and initialise shared memory for mapping arrays
    CALL allocate_shared_int_2D(       region%grid_GIA%ny, region%grid_GIA%nx, region%SELEN%map_nisl, region%SELEN%wmap_nisl)
    CALL allocate_shared_int_3D( nmax, region%grid_GIA%ny, region%grid_GIA%nx, region%SELEN%map_isl,  region%SELEN%wmap_isl )
    CALL allocate_shared_dp_3D(  nmax, region%grid_GIA%ny, region%grid_GIA%nx, region%SELEN%map_w,    region%SELEN%wmap_w   )
        
    region%SELEN%map_nisl(   :,region%grid_GIA%i1:region%grid_GIA%i2) = 0
    region%SELEN%map_isl(  :,:,region%grid_GIA%i1:region%grid_GIA%i2) = 0
    region%SELEN%map_w(    :,:,region%grid_GIA%i1:region%grid_GIA%i2) = 0._dp
    
    ! Go over all ice model square grid points
    DO i = region%grid_GIA%i1, region%grid_GIA%i2
    DO j = 1, region%grid_GIA%ny
    
      ! The coordinates of this GIA grid point
      lat_gr = region%grid_GIA%lat( j,i)
      lon_gr = region%grid_GIA%lon( j,i)
      
      ! Find the range of SELEN grid points that lie within the same latitude band
      iapl = 1
      DO WHILE (SELEN%mesh%ancplist_lat( iapl) < lat_gr - ((180._dp / pi) * (region%SELEN%rad + C%SELEN_alfa)))
        iapl = iapl + 1
      END DO
      iapl = MAX(1, iapl - 1)
      isll = SELEN%mesh%ancplist_isl( iapl,1)
      
      iapu = SELEN%mesh%nanc
      DO WHILE (SELEN%mesh%ancplist_lat( iapu) > lat_gr + ((180._dp / pi) * (region%SELEN%rad + C%SELEN_alfa)))
        iapu = iapu - 1
      END DO
      iapu = MIN(SELEN%mesh%nanc, iapu + 1)
      islu = SELEN%mesh%ancplist_isl( iapu,2)
      
      ! Inspect all SELEN global grid pixels that lie within the latitude band      
      DO isl = isll, islu
      
        ! Don't inspect SELEN global grid pixels that already belong to another ice model region
        IF (SELEN%icemodel_region( isl,1) > 0 .AND. SELEN%icemodel_region( isl,1) /= region_label) CYCLE
        
        ! The coordinates of this SELEN global grid pixels
        lon_sl = SELEN%mesh%lon( isl)
        lat_sl = SELEN%mesh%lat( isl)
        
        ! The angular distance between this ice model square grid point and this SELEN global grid point
        deg = angdist( lon_gr, lat_gr, lon_sl, lat_sl)
        
        w = 0._dp
        IF ( (180._dp / pi) * MIN( C%SELEN_alfa, region%SELEN%rad) <= ((180._dp / pi) * MAX( C%SELEN_alfa, region%SELEN%rad) - deg) ) THEN
          ! This ice model square grid point's circular disc is fully enclosed by this SELEN global grid pixel
          intarea  = 2._dp * pi * earth_radius**2 * (1._dp - COS( MIN( C%SELEN_alfa, region%SELEN%rad)))
          w = intarea / pixarea
        ELSEIF ( (180._dp / pi) * (C%SELEN_alfa + region%SELEN%rad) > deg ) THEN
          ! This ice model square grid point's circular disc is partly enclosed by this SELEN global grid pixel
          intarea  = optintersect( C%SELEN_alfa, region%SELEN%rad, deg)
          w = intarea / pixarea
        END IF
        
        IF (w > 0._dp) THEN
          ! List this contribution
          n = region%SELEN%map_nisl( j,i) + 1
          region%SELEN%map_nisl(     j,i) = n
          region%SELEN%map_isl(    n,j,i) = isl
          region%SELEN%map_w(      n,j,i) = w
          ! Mark this SELEN global grid pixel as part of this ice model region
          SELEN%icemodel_region( isl,1) = region_label
        END IF
        
      END DO ! DO isl = isll, islu
      
    END DO ! DO j = 1, region%grid_GIA%ny
    END DO ! DO i = region%grid_GIA%i1, region%grid_GIA%i2
    CALL sync
       
  END SUBROUTINE create_GIA_grid_to_SELEN_maps
  SUBROUTINE map_GIA_grid_to_SELEN( SELEN, region, d_gr, d_sl, do_scale)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(IN)    :: SELEN
    TYPE(type_model_region),             INTENT(IN)    :: region
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: d_gr
    REAL(dp), DIMENSION(:    ),          INTENT(INOUT) :: d_sl
    LOGICAL,                             INTENT(IN)    :: do_scale
    
    ! Local variables:
    INTEGER                                            :: ierr
    INTEGER                                            :: i,j,n,isl
    REAL(dp)                                           :: w
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: d_sl_temp
    REAL(dp)                                           :: int_gr, int_sl, R
    
    ALLOCATE( d_sl_temp( SELEN%mesh%nV))
    d_sl_temp = 0._dp
    
    DO i = region%grid_GIA%i1, region%grid_GIA%i2
    DO j = 1, region%grid_GIA%ny
      
      DO n = 1, region%SELEN%map_nisl(   j,i)
        isl = region%SELEN%map_isl(    n,j,i)
        w   = region%SELEN%map_w(      n,j,i)
        d_sl_temp( isl) = d_sl_temp( isl) + (w * d_gr( j,i))
      END DO ! DO n = 1, region%SELEN%map_nisl( i,j)
      
    END DO ! DO j = 1, region%ny_square
    END DO ! DO i = region%i1_square, region%i2_square
    
    CALL MPI_REDUCE( d_sl_temp, d_sl, SELEN%mesh%nV, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
    IF (do_scale) THEN
      ! Guarantee conservation of mass when mapping ice thickness
      
      int_gr = SUM( d_gr( :,region%grid_GIA%i1:region%grid_GIA%i2)) * C%dx_GIA**2
      int_sl = SUM( d_sl( C%SELEN_i1:C%SELEN_i2)) * 2._dp * pi * earth_radius**2 * (1._dp - COS(C%SELEN_alfa))
      
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, int_gr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, int_sl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      R = 1._dp
      IF (int_gr > 0._dp .AND. int_sl > 0._dp) THEN
        R = int_sl / int_gr
      END IF
      
      d_sl( C%SELEN_i1:C%SELEN_i2) = d_sl( C%SELEN_i1:C%SELEN_i2) / R
      CALL sync
      
    END IF ! IF (do_scale)
    
  END SUBROUTINE map_GIA_grid_to_SELEN
  
  FUNCTION angdist( lonp, latp, lonq, latq) RESULT( deg)
    ! Angular distance between points p and q
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp) :: lonp, latp, lonq, latq
    REAL(dp) :: deg
    
    deg =  (180._dp / pi) * ACOS( (COS( (pi / 180._dp) * (90._dp - latp)) * COS( (pi / 180._dp) * (90._dp - latq)) ) + &
                                  (SIN( (pi / 180._dp) * (90._dp - latp)) * SIN( (pi / 180._dp) * (90._dp - latq))   * & 
                                   COS( (pi / 180._dp) * (lonp - lonq)) ))
  END FUNCTION angdist  
  FUNCTION optintersect( radp, radi, deg) RESULT( intarea)
    
    IMPLICIT NONE
  
    ! In/output variables:
    REAL(dp)           :: radp, radi, deg
    REAL(dp)           :: intarea
    
    ! Local variables:
    INTEGER, PARAMETER :: a = 0
    INTEGER, PARAMETER :: b = 1
    REAL(dp)           :: x         ! = (1-(deg-abs(radp-radi))/(radp+radi-abs(radp-radi)))
    REAL(dp)           :: smooth    ! = ( (-2*((x-a/b-a))**3) + (3*((x-a)/(b-a))**2) ) 
    
    x       = 1._dp - (deg- (180._dp / pi) * abs(radp-radi)) / ((180._dp / pi) * (radp+radi-abs(radp-radi)))
    smooth  = (-2._dp * ( (x-dble(a)/dble(b-a)))**3) + (3._dp*((x-dble(a))/(dble(b-a)))**2)
    intarea = 2._dp * pi * earth_radius**2 * (1._dp - COS(MIN(radp,radi))) * smooth
    
  END FUNCTION optintersect

END MODULE SELEN_mapping_module

