/*
 *  General allocation routines 
 *
 */

#include "sputnik.h"
#include "macro.h"

int MacroAllocMatrixVector(MPI_Comm MACRO_COMM, int nlocal, int ntotal)
{
  /*
     Allocates memory for Mat A, Vec x, Vec dx, Vec b

     Input>
     a) MACRO_COMM > MPI communicator
     b) nlocal > # of local componenets
     c) ntotal > # of total componenets

     Output>
     a) A      > the PETSc matrix
     b) x,dx,b > the PETSc vectors

   */

  int rank, nproc, ierr;

  MPI_Comm_size(MACRO_COMM, &nproc);
  MPI_Comm_rank(MACRO_COMM, &rank);

  ierr = MatCreate(MACRO_COMM,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,ntotal,ntotal);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(A,81,NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,81,NULL,81,NULL);CHKERRQ(ierr);

  /*
     Create parallel vectors.
     - We form 1 vector from scratch and then duplicate as needed.
     This vector has ghost pathing in order to Get & and Set 
     ghost values in an easy way. This is a great PETSc tool.
     - When solving a linear system, the vectors and matrices MUST
     be partitioned accordingly.  PETSc automatically generates
     appropriately partitioned matrices and vectors when MatCreate()
     and VecCreate() are used with the same communicator.
   */
  int i, d, *ghostsIndex;
  ghostsIndex = malloc(NMyGhost*3* sizeof(int));

  for(i=0;i<NMyGhost;i++){
    for(d=0;d<3;d++){
      ghostsIndex[i*3+d] = loc2petsc[NMyNod + i]*3+d;
    }
  }

  //  ierr = VecCreate(comm,x);CHKERRQ(ierr);
  //  ierr = VecSetSizes(x,NMyNod,NTotalNod);CHKERRQ(ierr);
  ierr = VecCreateGhost(MACRO_COMM, NMyNod*3, NTotalNod*3, NMyGhost*3, ghostsIndex, &x); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&dx);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

  free(ghostsIndex);

  /*
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned.
   */
  int Istart, Iend;
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  if( Istart != StartIndexRank[rank_mac]*3 )
    SETERRQ(MACRO_COMM,1,"error on indeces set for matrix and vector.");
      
  if(rank_mac<nproc_mac-1){
    if( Iend != StartIndexRank[rank_mac+1]*3 ){
      SETERRQ(MACRO_COMM,1,"error on indeces set for matrix and vector.");
    }
  }
  else{
    if( Iend != NTotalNod*3 ){
      SETERRQ(MACRO_COMM,1,"error on indeces set for matrix and vector.");
    }
  }

  return 0;
}
