MODULE parallel_module

  ! A collection of different routines that make parallel programming a lot easier.

  USE mpi
  USE configuration_module,        ONLY: dp, n_MPI_windows
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER
  
  IMPLICIT NONE
    
  TYPE parallel_info
    
    INTEGER :: i        ! ID of this process (0 = master, >0 = slave)
    INTEGER :: n        ! Total number of processes (1 = single-core, >1 = master+slaves)
    LOGICAL :: master   ! Whether or not the current process is the master process
    
  END TYPE parallel_info
    
  TYPE(parallel_info), SAVE :: par
  INTEGER                   :: cerr, ierr

CONTAINS

  ! Synchronise the different processes
  SUBROUTINE sync
    ! Use MPI_BARRIER to synchronise all the processes         
    
    IMPLICIT NONE
    
    INTEGER                                           :: ierr
    
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
  
  END SUBROUTINE sync
    
  ! Allocate some shared memory space on all the processes (cannot be done on only one process).
  ! Only the master process allocates actual space, the slaves allocate zero bytes.
  ! All processes then associate their own pointers with the master's memory space,
  ! so that all processes can access the same physical memory and work in parallel.
  ! All processes must save the "window" to their own memory space, so that it can be deallocated if needed.
  SUBROUTINE allocate_shared_int_0D(                 p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                    POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = 4_MPI_ADDRESS_KIND
      disp_unit   = 4
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN  
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p)
    
    ! Initialise memory with zeros
    IF (par%master) p = 0
    CALL sync
    
    ! Update the max_window memory leak tracker
    n_MPI_windows = n_MPI_windows + 1
  
  END SUBROUTINE allocate_shared_int_0D  
  SUBROUTINE allocate_shared_int_1D(     n1,         p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1         ! Dimension(s) of memory to be allocated
    INTEGER,  DIMENSION(:    ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*4_MPI_ADDRESS_KIND
      disp_unit   = 4
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN  
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1])
    
    ! Initialise memory with zeros
    IF (par%master) p = 0
    CALL sync
    
    ! Update the max_window memory leak tracker
    n_MPI_windows = n_MPI_windows + 1
  
  END SUBROUTINE allocate_shared_int_1D  
  SUBROUTINE allocate_shared_int_2D(     n1, n2,     p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2     ! Dimension(s) of memory to be allocated
    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*n2*4_MPI_ADDRESS_KIND
      disp_unit   = n1*4
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN 
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2])
    
    ! Initialise memory with zeros
    IF (par%master) p = 0
    CALL sync
    
    ! Update the max_window memory leak tracker
    n_MPI_windows = n_MPI_windows + 1
  
  END SUBROUTINE allocate_shared_int_2D  
  SUBROUTINE allocate_shared_int_3D(     n1, n2, n3, p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2, n3 ! Dimension(s) of memory to be allocated
    INTEGER,  DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*n2*n3*4_MPI_ADDRESS_KIND
      disp_unit   = n1*n2*4
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN   
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2, n3])
    
    ! Initialise memory with zeros
    IF (par%master) p = 0
    CALL sync
    
    ! Update the max_window memory leak tracker
    n_MPI_windows = n_MPI_windows + 1
  
  END SUBROUTINE allocate_shared_int_3D  
  SUBROUTINE allocate_shared_dp_0D(                  p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    REAL(dp),                   POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = 8_MPI_ADDRESS_KIND
      disp_unit   = 8
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p)
    
    ! Initialise memory with zeros
    IF (par%master) p = 0
    CALL sync
    
    ! Update the max_window memory leak tracker
    n_MPI_windows = n_MPI_windows + 1
  
  END SUBROUTINE allocate_shared_dp_0D  
  SUBROUTINE allocate_shared_dp_1D(      n1,         p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1         ! Dimension(s) of memory to be allocated
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*8_MPI_ADDRESS_KIND
      disp_unit   = 8
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN  
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    CALL sync
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1])
    
    ! Initialise memory with zeros
    IF (par%master) p = 0._dp
    CALL sync
    
    ! Update the max_window memory leak tracker
    n_MPI_windows = n_MPI_windows + 1
  
  END SUBROUTINE allocate_shared_dp_1D  
  SUBROUTINE allocate_shared_dp_2D(      n1, n2,     p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2     ! Dimension(s) of memory to be allocated
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*n2*8_MPI_ADDRESS_KIND
      disp_unit   = n1*8
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN    
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    CALL sync
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2])
    
    ! Initialise memory with zeros
    IF (par%master) p = 0._dp
    CALL sync
    
    ! Update the max_window memory leak tracker
    n_MPI_windows = n_MPI_windows + 1
  
  END SUBROUTINE allocate_shared_dp_2D  
  SUBROUTINE allocate_shared_dp_3D(      n1, n2, n3, p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2, n3 ! Dimension(s) of memory to be allocated
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*n2*n3*8_MPI_ADDRESS_KIND
      disp_unit   = n1*n2*8
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN   
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    CALL sync
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2, n3])
    
    ! Initialise memory with zeros
    IF (par%master) p = 0._dp
    CALL sync
    
    ! Update the max_window memory leak tracker
    n_MPI_windows = n_MPI_windows + 1
  
  END SUBROUTINE allocate_shared_dp_3D  
  SUBROUTINE allocate_shared_bool_0D(                p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    LOGICAL,                    POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = 4_MPI_ADDRESS_KIND
      disp_unit   = 4
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN   
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p)
    
    ! Update the max_window memory leak tracker
    n_MPI_windows = n_MPI_windows + 1
  
  END SUBROUTINE allocate_shared_bool_0D
  SUBROUTINE allocate_shared_complex_1D( n1,         p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1         ! Dimension(s) of memory to be allocated
    COMPLEX*16, DIMENSION(:    ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*16_MPI_ADDRESS_KIND
      disp_unit   = 16
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN  
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1])
    
    ! Initialise memory with zeros
    IF (par%master) p = (0.,0.)
    CALL sync
    
    ! Update the max_window memory leak tracker
    n_MPI_windows = n_MPI_windows + 1
  
  END SUBROUTINE allocate_shared_complex_1D  
  SUBROUTINE allocate_shared_complex_2D( n1, n2,     p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2     ! Dimension(s) of memory to be allocated
    COMPLEX*16, DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*n2*16_MPI_ADDRESS_KIND
      disp_unit   = n1*16
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN    
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2])
    
    ! Initialise memory with zeros
    IF (par%master) p = (0.,0.)
    CALL sync
    
    ! Update the max_window memory leak tracker
    n_MPI_windows = n_MPI_windows + 1
  
  END SUBROUTINE allocate_shared_complex_2D 
  SUBROUTINE allocate_shared_complex_3D( n1, n2, n3, p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2, n3 ! Dimension(s) of memory to be allocated
    COMPLEX*16, DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*n2*n3*16_MPI_ADDRESS_KIND
      disp_unit   = n1*n2*16
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN   
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2, n3])
    
    ! Initialise memory with zeros
    IF (par%master) p = (0.,0.)
    CALL sync
    
    ! Update the max_window memory leak tracker
    n_MPI_windows = n_MPI_windows + 1
  
  END SUBROUTINE allocate_shared_complex_3D   
    
  ! Use the "window" to the allocated memory space (which consists of zero bytes for the slave)
  ! to deallocate that memory.
  SUBROUTINE deallocate_shared( win)
    ! Use MPI_WIN_FREE to deallocate shared memory space for an array.
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(INOUT) :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    
    CALL MPI_WIN_FREE( win, ierr)
    
    n_MPI_windows = n_MPI_windows - 1
  
  END SUBROUTINE deallocate_shared
    
  SUBROUTINE partition_list( ntot, i, n, i1, i2)
    ! Partition a list into parallel ranges (e.g. vertex domains)
  
    ! In/output variables:
    INTEGER,                    INTENT(IN)        :: ntot, i, n
    INTEGER,                    INTENT(OUT)       :: i1, i2
    
    ! Local variables:
    INTEGER                                       :: vi
    INTEGER, DIMENSION(:    ), ALLOCATABLE        :: to_which_process_do_I_belong
    
    ALLOCATE( to_which_process_do_I_belong( ntot))
    
    DO vi = 1, ntot
      to_which_process_do_i_belong( vi) = NINT( REAL(vi-1,dp) * REAL(n-1,dp) / REAL(ntot-1,dp) )
    END DO
    
    i1 = 1
    DO WHILE (to_which_process_do_i_belong( i1) /= i)
      i1 = i1+1
      IF (i1 > ntot) THEN
        i1 =  0
        i2 = -1
        RETURN
      END IF
    END DO
    
    i2 = ntot
    DO WHILE (to_which_process_do_i_belong( i2) /= i)
      i2 = i2-1
      IF (i2 < 1) THEN
        i1 =  0
        i2 = -1
        RETURN
      END IF
    END DO
    
  END SUBROUTINE partition_list
  
END MODULE parallel_module
