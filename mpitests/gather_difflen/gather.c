#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define nstrings 5

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv); 

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(size!=3){
    MPI_Finalize();
    printf("execution with np = 3 only\n");
    return 1;
  }

  int a[5],b[5],c[5],len_a,len_b,len_c;

  if(rank == 0){ 
    a[0] = 0 ; 
    b[0] = 30; 
    c[0] = 60; 
    len_a = len_b = len_c = 1;
    //a,b,c can have different lengths 
  }
  else if(rank == 1){
    a[0] = 1 ; a[1] = 11; 
    b[0] = 31; b[1] = 41; 
    c[0] = 61; c[1] = 71; 
    len_a = len_b = len_c = 2;
    //a,b,c can have different lengths 
  }
  else if(rank == 2){
    a[0] = 2 ; a[1] = 12; a[2] = 22;
    b[0] = 32; b[1] = 42; b[2] = 52;
    c[0] = 62; c[1] = 72; c[2] = 82;
    len_a = len_b = len_c = 3;
    //a,b,c can have different lengths 
  }

  /*

     At the end I should have

     rank 0 : a = { 0  | 1  11  | 2  12 22 }
     rank 1 : b = { 30 | 31 41  | 32 42 52 }
     rank 2 : c = { 60 | 61 71  | 62 72 82 }

   */

  int lengths[3];
  int tot_length;

  /*
   * we Gather the lengths for every process as root 
   */

  int mylen,i,root = 0;

  for(i=0;i<3;i++){
    root = i;

    mylen = rank + 1;

    MPI_Gather(&mylen, 
	1, 
	MPI_INT,
	lengths, 
	1,
	MPI_INT,
	root,
	MPI_COMM_WORLD);

  }

  /*
     calculate tot_length of my elem array 
     and determine the disp vector
   */

  int displs[3];
  displs[0] = 0;
  tot_length = lengths[0];

  for(i=1;i<3;i++){
    displs[i] = displs[i-1] + lengths[i-1];
    tot_length += lengths[i];
  }

  printf("rank %d lengths = %d %d %d displs = %d %d %d \n", rank, \
      lengths[0],lengths[1],lengths[2],                           \
      displs[0],displs[1],displs[2]);

  int *elem = malloc(tot_length * sizeof(int));
  int  myvec[5];
  int  j;

  for(i=0;i<3;i++){

    root = i;
    if(root == 0){
      for(j=0;j<5;j++)
	myvec[j]=a[j];
    }
    else if(root == 1){
      for(j=0;j<5;j++)
	myvec[j]=b[j];
    }
    else if(root == 2){
      for(j=0;j<5;j++)
	myvec[j]=c[j];
    }

    MPI_Gatherv(myvec,
	mylen,
	MPI_INT,
	elem,
	lengths,
	displs,
	MPI_INT,
	root,
	MPI_COMM_WORLD);

  }

  printf("rank %d elem : ",rank);
  for(i=0;i<tot_length;i++){
    printf("%d ", elem[i]);
  }
  printf("\n");


  MPI_Finalize();
  return 0;
}
