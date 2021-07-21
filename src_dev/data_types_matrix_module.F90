MODULE data_types_matrix_module
  ! Contains the TYPE representing a matrix equation (in regular Fortran data, as well as in PETSc types)
  
#include <petsc/finclude/petscksp.h>

  USE configuration_module,        ONLY: dp, C
  USE petscksp

  IMPLICIT NONE
  
  TYPE type_sparse_matrix_CSR
    ! A matrix equation Ax=b, represented in the Compressed Sparse Row (CSR) format
    
    INTEGER,                    POINTER     :: m,n             ! A = [m-by-n]
    INTEGER,                    POINTER     :: nnz_per_row_max ! Maximum number of non-zero entries per row in A (determines how much memory is allocated)
    INTEGER,                    POINTER     :: nnz_max         ! Maximum number of non-zero entries         in A (determines how much memory is allocated)
    INTEGER,  DIMENSION(:    ), POINTER     :: A_ptr
    INTEGER,  DIMENSION(:    ), POINTER     :: A_index
    REAL(dp), DIMENSION(:    ), POINTER     :: A_val
    REAL(dp), DIMENSION(:    ), POINTER     :: x
    REAL(dp), DIMENSION(:    ), POINTER     :: b
    INTEGER :: wm, wn, wnnz_per_row_max, wnnz_max, wA_ptr, wA_index, wA_val, wx, wb
    
    ! PETSc version of the same, including the KSP solver
    TYPE(tMat)                              :: pA
    TYPE(tVec)                              :: pb
    TYPE(tVec)                              :: px
    TYPE(tKSP)                              :: KSP_solver
    
  END TYPE type_sparse_matrix_CSR
  
CONTAINS

END MODULE data_types_matrix_module
