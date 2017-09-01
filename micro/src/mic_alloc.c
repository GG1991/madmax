/*
 *  General allocation routines 
 *
 */

#include "sputnik.h"
#include "micro.h"

int mic_alloc(MPI_Comm MICRO_COMM)
{
  /*
     Allocates memory for Mat A, Vec x, Vec dx, Vec b

     Input>
     a) MICRO_COMM > MPI communicator
     b) nlocal > # of local componenets
     c) ntotal > # of total componenets

     Output>
     a) A      > the PETSc matrix
     b) x,dx,b > the PETSc vectors

   */

  int rank, nproc, ierr, nlocal, ntotal;

  if(flag_reactions == PETSC_TRUE){
    nlocal = 3*nmynods + 3*nmybcnods;
    ntotal = 3*NTotalNod+ 3*nallbcnods;
  }else{
    nlocal = 3*nmynods;
    ntotal = 3*NTotalNod;
  }


  MPI_Comm_size(MICRO_COMM, &nproc);
  MPI_Comm_rank(MICRO_COMM, &rank);

  ierr = MatCreate(MICRO_COMM,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,ntotal,ntotal);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSeqAIJSetPreallocation(A,117,NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,117,NULL,117,NULL);CHKERRQ(ierr);
  ierr = MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);

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
  ghostsIndex = malloc(nghost*3* sizeof(int));

  for(i=0;i<nghost;i++){
    for(d=0;d<3;d++){
      ghostsIndex[i*3+d] = loc2petsc[nmynods + i]*3+d;
    }
  }

  //  ierr = VecCreate(comm,x);CHKERRQ(ierr);
  //  ierr = VecSetSizes(x,nmynods,NTotalNod);CHKERRQ(ierr);
  ierr = VecCreateGhost(MICRO_COMM, nmynods*3, NTotalNod*3, nghost*3, ghostsIndex, &x); CHKERRQ(ierr);
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
  if( Istart != StartIndexRank[rank_mic]*3 )
    SETERRQ(MICRO_COMM,1,"error on indeces set for matrix and vector.");
      
  if(rank_mic<nproc_mic-1){
    if( Iend != StartIndexRank[rank_mic+1]*3 ){
      SETERRQ(MICRO_COMM,1,"error on indeces set for matrix and vector.");
    }
  }
  else{
    if( Iend != NTotalNod*3 ){
      SETERRQ(MICRO_COMM,1,"error on indeces set for matrix and vector.");
    }
  }

  return 0;
}
