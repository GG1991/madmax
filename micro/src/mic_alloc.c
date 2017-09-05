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

  nlocal = 3*nmynods;
  ntotal = 3*NTotalNod;

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

  ierr = VecCreateGhost(MICRO_COMM, nmynods*3, NTotalNod*3, nghost*3, ghostsIndex, &x); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&dx);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

  free(ghostsIndex);

  if(homo.type==LD_LAGRAN_SEQ){

    /*
       Linear displacements with Lagrangian Multiplier approach
    */

    int nlocal_ext, ntotal_ext;

    nlocal_ext = nlocal + ((homog_ld_lagran_t*)homo.st)->nnods_bc * 3;

    ierr = MPI_Allreduce(&nlocal_ext, &ntotal_ext, 1, MPI_INT, MPI_SUM, MICRO_COMM);

    ierr = MatCreate(MICRO_COMM,&J);CHKERRQ(ierr);
    ierr = MatSetSizes(J,nlocal_ext,nlocal_ext,ntotal_ext,ntotal_ext);CHKERRQ(ierr);
    ierr = MatSetFromOptions(J);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(J,117,NULL);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(J,117,NULL,117,NULL);CHKERRQ(ierr);
    ierr = MatSetOption(J,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ(ierr);

    ierr = VecCreate(MICRO_COMM, &xe);CHKERRQ(ierr);
    ierr = VecSetSizes(xe, nlocal_ext, ntotal_ext);CHKERRQ(ierr);
  }

  return 0;
}
