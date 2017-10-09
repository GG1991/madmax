/*
    Allocation routines 
    The sizes of the matrices depends on the homogenization technique 
    that it is going to be used
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

  if(homo_type==UNIF_STRAINS){

    int rank, nproc, ierr, nlocal, ntotal;

    nlocal = dim*nmynods;
    ntotal = dim*ntotnod;

    MPI_Comm_size(MICRO_COMM, &nproc);
    MPI_Comm_rank(MICRO_COMM, &rank);

    ierr = MatCreate(MICRO_COMM,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,nlocal,nlocal,ntotal,ntotal);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A,117,NULL);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(A,117,NULL,117,NULL);CHKERRQ(ierr);
    ierr = MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);

    int i, d, *ghostsIndex;
    ghostsIndex = malloc(nghost*dim* sizeof(int));

    for(i=0;i<nghost;i++){
      for(d=0;d<dim;d++){
	ghostsIndex[i*dim+d] = loc2petsc[nmynods + i]*dim+d;
      }
    }

    ierr = VecCreateGhost(MICRO_COMM, nmynods*dim, ntotnod*dim, nghost*dim, ghostsIndex, &x); CHKERRQ(ierr);
    ierr = VecDuplicate(x,&dx);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&b1);CHKERRQ(ierr);

    free(ghostsIndex);

  }

  return 0;
}
