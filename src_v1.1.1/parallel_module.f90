MODULE parallel_module
  ! A collection of different routines that make parallel programming in UFEMISM a lot easier.

  USE mpi
  USE configuration_module,        ONLY: dp
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER
  
  IMPLICIT NONE
    
  TYPE parallel_info
    
    INTEGER :: i        ! ID of this process (0 = master, >0 = slave)
    INTEGER :: n        ! Total number of processes (1 = single-core, >1 = master+slaves)
    LOGICAL :: master   ! Whether or not the current process is the master process
    
  END TYPE parallel_info
    
  TYPE(parallel_info), SAVE :: par
  
  LOGICAL :: debug_check_for_memory_leaks = .FALSE.

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
  SUBROUTINE allocate_shared_int_0D(              p, win)
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
    ELSE
      ! Initialise memory with zero
      p = 0
    END IF
    CALL sync
    
    IF (par%master .AND. debug_check_for_memory_leaks) WRITE(0,*) '     allocate_shared_int_0D: win = ', win
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p)
  
  END SUBROUTINE allocate_shared_int_0D  
  SUBROUTINE allocate_shared_int_1D(  n1,         p, win)
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
    ELSE
      ! Initialise memory with zero
      p = 0
    END IF
    CALL sync
    
    IF (par%master .AND. debug_check_for_memory_leaks) WRITE(0,*) '     allocate_shared_int_1D: win = ', win
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1])
  
  END SUBROUTINE allocate_shared_int_1D  
  SUBROUTINE allocate_shared_int_2D(  n1, n2,     p, win)
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
    ELSE
      ! Initialise memory with zero
      p = 0
    END IF
    CALL sync
    
    IF (par%master .AND. debug_check_for_memory_leaks) WRITE(0,*) '     allocate_shared_int_2D: win = ', win
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2])
  
  END SUBROUTINE allocate_shared_int_2D  
  SUBROUTINE allocate_shared_int_3D(  n1, n2, n3, p, win)
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
    ELSE
      ! Initialise memory with zero
      p = 0
    END IF
    CALL sync
    
    IF (par%master .AND. debug_check_for_memory_leaks) WRITE(0,*) '     allocate_shared_int_3D: win = ', win
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2, n3])
  
  END SUBROUTINE allocate_shared_int_3D  
  SUBROUTINE allocate_shared_dp_0D(               p, win)
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
    ELSE
      ! Initialise memory with zero
      p = 0._dp
    END IF
    CALL sync
    
    IF (par%master .AND. debug_check_for_memory_leaks) WRITE(0,*) '     allocate_shared_dp_0D:  win = ', win
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p)
  
  END SUBROUTINE allocate_shared_dp_0D  
  SUBROUTINE allocate_shared_dp_1D(   n1,         p, win)
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
    ELSE
      ! Initialise memory with zero
      p = 0._dp
    END IF
    CALL sync
    
    IF (par%master .AND. debug_check_for_memory_leaks) WRITE(0,*) '     allocate_shared_dp_1D:  win = ', win
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1])
  
  END SUBROUTINE allocate_shared_dp_1D  
  SUBROUTINE allocate_shared_dp_2D(   n1, n2,     p, win)
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
    ELSE
      ! Initialise memory with zero
      p = 0._dp
    END IF
    CALL sync
    
    IF (par%master .AND. debug_check_for_memory_leaks) WRITE(0,*) '     allocate_shared_dp_2D:  win = ', win
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2])
  
  END SUBROUTINE allocate_shared_dp_2D  
  SUBROUTINE allocate_shared_dp_3D(   n1, n2, n3, p, win)
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
    ELSE
      ! Initialise memory with zero
      p = 0._dp
    END IF
    CALL sync
    
    IF (par%master .AND. debug_check_for_memory_leaks) WRITE(0,*) '     allocate_shared_dp_3D:  win = ', win
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2, n3])
  
  END SUBROUTINE allocate_shared_dp_3D  
  SUBROUTINE allocate_shared_bool_1D( n1,         p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1         ! Dimension(s) of memory to be allocated
    LOGICAL,  DIMENSION(:    ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
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
    ELSE
      ! Initialise memory with zero
      p = .FALSE.
    END IF
    CALL sync
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1])
  
  END SUBROUTINE allocate_shared_bool_1D
    
  ! Use the "window" to the allocated memory space (which consists of zero bytes for the slave)
  ! to deallocate that memory.
  SUBROUTINE deallocate_shared( win)
    ! Use MPI_WIN_FREE to deallocate shared memory space for an array.
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(INOUT) :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    
    IF (par%master .AND. debug_check_for_memory_leaks) WRITE(0,*) '     deallocate_shared:      win = ', win
    
    CALL MPI_WIN_FREE( win, ierr)
  
  END SUBROUTINE deallocate_shared
    
  SUBROUTINE partition_list( ntot, i, n, i1, i2)
    ! Partition a list into parallel ranges (e.g. vertex domains)
  
    ! In/output variables:
    INTEGER,                    INTENT(IN)        :: ntot, i, n
    INTEGER,                    INTENT(OUT)       :: i1, i2
    
    IF (ntot > n*2) THEN
      i1 = MAX(1,    FLOOR(REAL(ntot *  i      / n)) + 1)
      i2 = MIN(ntot, FLOOR(REAL(ntot * (i + 1) / n)))
    ELSE
      IF (i==0) THEN
        i1 = 1
        i2 = ntot
      ELSE
        i1 = 1
        i2 = 0
      END IF
    END IF
    
  END SUBROUTINE partition_list
  
END MODULE parallel_module
