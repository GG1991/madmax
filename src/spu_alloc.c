/*
 *  General allocation routines 
 *
 */

#include "sputnik.h"

int AllocMatrixVector(MPI_Comm *comm)
{

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create parallel matrix, specifying only its global dimensions.
     When using MatCreate(), the matrix format can be specified at
     runtime. Also, the parallel partitioning of the matrix is
     determined by PETSc at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.
  */

  int Istart, Iend;
  int ierr;

  MPI_Comm_size(*comm, &nproc);
  MPI_Comm_rank(*comm, &rank);

  ierr = MatCreate(*comm,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,NMyNod*3,NMyNod*3,NTotalNod*3,NTotalNod*3);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,81,NULL,81,NULL);CHKERRQ(ierr);

  /*
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned.
  */
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  if( Istart != StartIndexRank[rank]*3 ){
    printf("AllocMatrixVector: error on indeces set for matrix and vector.\n");
  }
  if(rank<nproc-1){
    if( Iend != StartIndexRank[rank+1]*3 ){
      printf("AllocMatrixVector: error on indeces set for matrix and vector.\n");
    }
  }
  else{
    if( Iend != NTotalNod*3 ){
      printf("AllocMatrixVector: error on indeces set for matrix and vector.\n");
    }
  }


  /*
     Create parallel vectors.
      - We form 1 vector from scratch and then duplicate as needed.
      - When using VecCreate(), VecSetSizes and VecSetFromOptions()
        in this example, we specify only the
        vector's global dimension; the parallel partitioning is determined
        at runtime.
      - When solving a linear system, the vectors and matrices MUST
        be partitioned accordingly.  PETSc automatically generates
        appropriately partitioned matrices and vectors when MatCreate()
        and VecCreate() are used with the same communicator.
      - The user can alternatively specify the local vector and matrix
        dimensions when more sophisticated partitioning is needed
        (replacing the PETSC_DECIDE argument in the VecSetSizes() statement
        below).
  */
  ierr = VecCreate(*comm,&u);CHKERRQ(ierr);
  ierr = VecSetSizes(u,NMyNod*3,m*n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&b);CHKERRQ(ierr);
  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);


  return 0;
}
