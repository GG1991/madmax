#include "macro.h"

int alloc_memory(void){

  int ierr;

  if(params.calc_mode == CALC_MODE_EIGEN || params.calc_mode == CALC_MODE_NORMAL){

    int nnz = dim * MAX_ADJ_NODES;

    ierr = MatCreate(MACRO_COMM, &A);
    ierr = MatSetSizes(A, dim*nmynods, dim*nmynods, dim*ntotnod, dim*ntotnod);
    ierr = MatSetType(A, MATAIJ);
    ierr = MatSeqAIJSetPreallocation(A, nnz, NULL);
    ierr = MatMPIAIJSetPreallocation(A, nnz, NULL, nnz, NULL );
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    ierr = MatSetUp(A);
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    int *ghost_index = malloc(nghost*dim*sizeof(int));

    for(int i = 0 ; i < nghost ; i++){
      for(int d = 0 ; d < dim ; d++)
	ghost_index[i*dim + d] = loc2petsc[nmynods + i]*dim + d;
    }

    ierr = VecCreateGhost(MACRO_COMM, dim*nmynods, dim*ntotnod, nghost*dim, ghost_index, &x);

    if(params.calc_mode == CALC_MODE_EIGEN){

      ierr = MatCreate(MACRO_COMM, &M);
      ierr = MatSetSizes(M, dim*nmynods, dim*nmynods, dim*ntotnod, dim*ntotnod);
      ierr = MatSetType(M, MATAIJ);
      ierr = MatSeqAIJSetPreallocation(M, nnz, NULL);
      ierr = MatMPIAIJSetPreallocation(M, nnz, NULL, nnz, NULL);
      ierr = MatSetOption(M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      ierr = MatSetUp(M);
      ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
      ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);

    }
    else if(params.calc_mode == CALC_MODE_NORMAL){

      ierr = VecDuplicate(x, &dx);
      ierr = VecDuplicate(x, &b);

    }
  }

  if(ierr == 0) flags.allocated = true;

  return ierr;
}
