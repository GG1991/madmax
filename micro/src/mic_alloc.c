/*
 *  General allocation routines 
 *
 */

#include "sputnik.h"
#include "micro.h"

int MicroAllocMatrixVector(MPI_Comm comm, int nlocal, int ntotal)
{
  /*
     Allocates memory for Mat A, Vec x, Vec dx, Vec b

     Input>
     a) comm > MPI communicator
     b) nlocal > # of local componenets
     c) ntotal > # of total componenets

     Output>
     a) A      > the PETSc matrix
     b) x,dx,b > the PETSc vectors

   */

  int rank, nproc, ierr;

  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  ierr = MatCreate(comm,&A);CHKERRQ(ierr);
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
  ierr = VecCreateGhost(comm, NMyNod*3, NTotalNod*3, NMyGhost*3, ghostsIndex, &x); CHKERRQ(ierr);
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
  if( Istart != StartIndexRank[rank_mic]*3 ){
    printf("mic_main: error on indeces set for matrix and vector.\n");
    return 1;
  }
  if(rank_mic<nproc_mic-1){
    if( Iend != StartIndexRank[rank_mic+1]*3 ){
      printf("mic_main: error on indeces set for matrix and vector.\n");
      return 1;
    }
  }
  else{
    if( Iend != NTotalNod*3 ){
      printf("mic_main: error on indeces set for matrix and vector.\n");
      return 1;
    }
  }

  return 0;
}