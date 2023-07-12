MODULE petsc_module

  ! Contains routines that use the PETSc matrix solvers

#include <petsc/finclude/petscksp.h>

  USE mpi
  USE petscksp
  USE parallel_module,                 ONLY: par, sync, cerr, ierr, partition_list, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE data_types_module,               ONLY: type_sparse_matrix_CSR

  IMPLICIT NONE

CONTAINS

! == Solve a square CSR matrix equation with PETSc
  SUBROUTINE solve_matrix_equation_CSR_PETSc( CSR, PETSc_rtol, PETSc_abstol)
    ! Solve the matrix equation using a Krylov solver from PETSc

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: CSR
    REAL(dp),                            INTENT(IN)    :: PETSc_rtol
    REAL(dp),                            INTENT(IN)    :: PETSc_abstol

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_matrix_equation_CSR_PETSc'
    TYPE(PetscErrorCode)                               :: perr
    TYPE(PetscInt)                                     :: its
    TYPE(tMat)                                         :: A
    TYPE(tVec)                                         :: b
    TYPE(tVec)                                         :: x
    TYPE(tKSP)                                         :: KSP_solver
    TYPE(PetscReal)                                    :: rtol, abstol

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Convert the CSR matrix from the native IMAU-ICE storage structure to a PETSc matrix
    CALL convert_CSR_to_petsc( CSR, A, b, x)

    ! Set up the KSP solver
    CALL KSPcreate( PETSC_COMM_WORLD, KSP_solver, perr)

    ! Set operators. Here the matrix that defines the linear system
    ! also serves as the preconditioning matrix.
    CALL KSPSetOperators( KSP_solver, A, A, perr)

    ! Iterative solver tolerances
    rtol   = PETSc_rtol     ! The relative convergence tolerance; relative decrease in the (possibly preconditioned) residual norm
    abstol = PETSc_abstol   ! The absolute convergence tolerance; absolute size of the (possibly preconditioned) residual norm
    CALL KSPSetTolerances( KSP_solver, rtol, abstol, PETSC_DEFAULT_REAL, PETSC_DEFAULT_INTEGER, perr)

    ! Set runtime options, e.g.,
    !     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    ! These options will override those specified above as long as
    ! KSPSetFromOptions() is called _after_ any other customization routines.
    CALL KSPSetFromOptions( KSP_solver, perr)

    ! Solve the linear system
    CALL KSPSolve( KSP_solver, b, x, perr)

    ! Find out how many iterations it took
    CALL KSPGetIterationNumber( KSP_solver, its, perr)
    !IF (par%master) WRITE(0,*) '   PETSc solved Ax=b in ', its, ' iterations'

    ! Get the solution back to the native IMAU-ICE storage structure
    CALL convert_petsc_solution_to_CSR( CSR, x)
    CALL sync

    ! Clean up after yourself
    CALL KSPDestroy( KSP_solver, perr)
    CALL VecDestroy( x, perr)
    CALL VecDestroy( b, perr)
    CALL MatDestroy( A, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_matrix_equation_CSR_PETSc

! == Convert native IMAU-ICE CSR matrix to PETSc matrix
  SUBROUTINE convert_CSR_to_petsc( CSR, A, b, x)
    ! Convert a native IMAU-ICE Compressed Sparse Row (CSR) type matrix
    ! to a PETSc matrix.
    !
    ! Mostly copied from: https://petsc.org/release/documentation/manual/getting_started/#parallel-programming

    IMPLICIT NONE

    ! In- and output variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'convert_CSR_to_petsc'
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: CSR
    TYPE(tMat),                          INTENT(INOUT) :: A
    TYPE(tVec),                          INTENT(INOUT) :: b
    TYPE(tVec),                          INTENT(INOUT) :: x

    ! Local variables:
    TYPE(PetscErrorCode)                               :: perr
    INTEGER                                            :: k,i,j,istart,iend
    TYPE(PetscReal)                                    :: v

    ! Add routine to path
    CALL init_routine( routine_name)

  ! == Matrix A ==
  ! ==============

    ! Initialise the matrix object
    CALL MatCreate( PETSC_COMM_WORLD, A, perr)

    ! Set the matrix type to parallel (MPI) Aij
    CALL MatSetType( A, 'mpiaij', perr)

    ! Set the size, let PETSc automatically determine parallelisation domains
    CALL MatSetSizes( A, PETSC_DECIDE, PETSC_DECIDE, CSR%m, CSR%n, perr)

    ! Not entirely sure what this one does, but apparently it's really important
    CALL MatSetFromOptions( A, perr)

    ! Tell PETSc how much memory needs to be allocated
    CALL MatMPIAIJSetPreallocation( A, CSR%nnz_per_row_max+1, PETSC_NULL_INTEGER, CSR%nnz_per_row_max+1, PETSC_NULL_INTEGER, perr)

    ! Get parallelisation domains ("ownership ranges")
    CALL MatGetOwnershipRange( A, istart, iend, perr)

    ! Fill in matrix values
    DO i = istart+1,iend ! +1 because PETSc indexes from 0
      DO k = CSR%A_ptr( i), CSR%A_ptr( i+1)-1

        j = CSR%A_index( k)
        v = CSR%A_val(   k)

        CALL MatSetValues( A, 1, i-1, 1, j-1, v, INSERT_VALUES, perr)

      END DO
    END DO
    CALL sync

  ! == Vectors b,x ==
  ! =================

    ! Create parallel vectors.
    CALL VecCreate( PETSC_COMM_WORLD, x, perr)
    CALL VecSetSizes( x, PETSC_DECIDE, CSR%n, perr)
    CALL VecSetFromOptions( x, perr)

    CALL VecCreate( PETSC_COMM_WORLD, b, perr)
    CALL VecSetSizes( b, PETSC_DECIDE, CSR%m, perr)
    CALL VecSetFromOptions( b, perr)

    ! Fill in vector values
    DO i = istart+1,iend ! +1 because PETSc indexes from 0

      v = CSR%b( i)
      CALL VecSetValues( b, 1, i-1, v, INSERT_VALUES, perr)

      v = CSR%x( i)
      CALL VecSetValues( x, 1, i-1, v, INSERT_VALUES, perr)

    END DO
    CALL sync

    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.

    CALL MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY, perr)
    CALL VecAssemblyBegin( b, perr)
    CALL VecAssemblyBegin( x, perr)

    CALL MatAssemblyEnd(   A, MAT_FINAL_ASSEMBLY, perr)
    CALL VecAssemblyEnd(   b, perr)
    CALL VecAssemblyEnd(   x, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE convert_CSR_to_petsc
  SUBROUTINE convert_petsc_solution_to_CSR( CSR, x)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: CSR
    TYPE(tVec),                          INTENT(IN)    :: x

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'convert_petsc_solution_to_CSR'
    TYPE(PetscErrorCode)                               :: perr
    TYPE(PetscInt)                                     :: istart,iend,i
    TYPE(PetscInt),    DIMENSION(1)                    :: ix
    TYPE(PetscScalar), DIMENSION(1)                    :: v

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get parallelisation domains ("ownership ranges")
    CALL VecGetOwnershipRange( x, istart, iend, perr)

    ! Get values
    DO i = istart+1,iend
      ix(1) = i-1
      CALL VecGetValues( x, 1, ix, v, perr)
      CSR%x( i) = v(1)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE convert_petsc_solution_to_CSR

! == General PETSc initialisation and finalisation
  SUBROUTINE initialise_petsc
    ! Initialise PETSc

    IMPLICIT NONE

    ! Local variables:
    TYPE(PetscErrorCode)                               :: perr

    ! Initialise PETSc MPI stuff
    CALL PetscInitialize( PETSC_NULL_CHARACTER, perr)

  END SUBROUTINE initialise_petsc
  SUBROUTINE finalise_petsc
    ! Finalise PETSc

    IMPLICIT NONE

    TYPE(PetscErrorCode)                               :: perr

    ! Finalise PETSc MPI stuff
    CALL PetscFinalize( perr)

  END SUBROUTINE finalise_petsc

END MODULE petsc_module
