/*
 *  General allocation routines 
 *
 */

#include "sputnik.h"

int AllocMatrixVector(MPI_Comm comm, int nlocal, int ntotal, Mat *A, Vec *x, Vec *b)
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

     Input:
     a) comm : MPI communicator
     b) nlocal : # of local componenets
     c) ntotal : # of total componenets
     
     Output:
     a) A : the PETSc matrix
     b) x,b : the PETSc vectors
  */

  int rank, nproc, ierr;

  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  ierr = MatCreate(comm,A);CHKERRQ(ierr);
  ierr = MatSetSizes(*A,nlocal,nlocal,ntotal,ntotal);CHKERRQ(ierr);
  ierr = MatSetFromOptions(*A);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(*A,81,NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(*A,81,NULL,81,NULL);CHKERRQ(ierr);

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
  ierr = VecCreate(comm,x);CHKERRQ(ierr);
  ierr = VecSetSizes(*x,nlocal,ntotal);CHKERRQ(ierr);
  ierr = VecSetFromOptions(*x);CHKERRQ(ierr);
  ierr = VecDuplicate(*x,b);CHKERRQ(ierr);

  return 0;
}
