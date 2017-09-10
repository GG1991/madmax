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

  if(homo.type==LD){

    int rank, nproc, ierr, nlocal, ntotal;

    nlocal = dim*nmynods;
    ntotal = dim*NTotalNod;

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

    ierr = VecCreateGhost(MICRO_COMM, nmynods*dim, NTotalNod*dim, nghost*dim, ghostsIndex, &x); CHKERRQ(ierr);
    ierr = VecDuplicate(x,&dx);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&b1);CHKERRQ(ierr);

    free(ghostsIndex);

  }else if(homo.type==LD_LAGRAN_SEQ){

    /*
       Linear displacements with Lagrangian Multiplier approach
     */

    int n,ierr;

    n = (nmynods + ((homog_ld_lagran_t*)homo.st)->nnods_bc) * dim;

    ierr = MatCreateSeqAIJ(MICRO_COMM,n,n,177,PETSC_NULL,&A);CHKERRQ(ierr);
    ierr = MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);

    ierr = VecCreateSeq(MICRO_COMM, n, &x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&dx);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&b1);CHKERRQ(ierr);

  }else if(homo.type==LD_LAGRAN_PAR){

//    nlocal_ext = nlocal + ((homog_ld_lagran_t*)homo.st)->nnods_bc * 3;
//
//    ierr = MPI_Allreduce(&nlocal_ext, &ntotal_ext, 1, MPI_INT, MPI_SUM, MICRO_COMM);
//    
//    ierr = MatCreate(MICRO_COMM,&J);CHKERRQ(ierr); // this communicator should have only 1 process here
//    ierr = MatSetSizes(J,nlocal_ext,nlocal_ext,ntotal_ext,ntotal_ext);CHKERRQ(ierr);
//    ierr = MatSetFromOptions(J);CHKERRQ(ierr);
//    ierr = MatSeqAIJSetPreallocation(J,117,NULL);CHKERRQ(ierr);
//    ierr = MatMPIAIJSetPreallocation(J,117,NULL,117,NULL);CHKERRQ(ierr);
//    ierr = MatSetOption(J,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);
//
//    ierr = VecCreate(MICRO_COMM, &xe);CHKERRQ(ierr);
//    ierr = VecSetSizes(xe, nlocal_ext, ntotal_ext);CHKERRQ(ierr);

  }

  return 0;
}
