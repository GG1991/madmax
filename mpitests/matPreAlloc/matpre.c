#include "petscksp.h"
#include "stdio.h"

int main(int argc, char **argv)
{

  int    Istart, Iend, ierr;
  Mat    A;
  ierr = PetscInitialize(&argc,&argv,(char*)0,NULL);
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,1000000,1000000);CHKERRQ(ierr);
//  ierr = MatSetType(A,MATMPIAIJ);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,81,NULL,81,NULL);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  printf("Istart %d, Iend %d\n",Istart,Iend);

  ierr = PetscFinalize();
  return 0;
}
