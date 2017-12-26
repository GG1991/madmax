#include "mesh.h"
#include <stdbool.h>

#define NRM  "\x1B[0m"
#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"

#define PRINT_ARRAY(array, length) {\
  for(int i = 0 ; i < length ; i++){ \
    printf("%d ", array[i]); \
  } \
  printf("\n"); \
}

#define COMP_ARRAY(arr_1, len_1, arr_2, len_2, flag) {\
  flag = true;\
  if(len_1 != len_2){\
    flag = false;\
  }else{\
    int i = 0;\
    while (i < len_1){\
      if(arr_1[i] != arr_2[i]){\
	flag = false;\
	break;\
      }\
      i++;\
    }\
  }\
}

int array_in_1[] = {1, 3, 4, 4, 3, 1, 1, 6, 5};
int array_out_1[] = {1, 3, 4, 5, 6};
mesh_t mesh;

int main(void){

  int number_of_test = 1;
  bool flag_pass = true;

  MPI_Init(NULL, NULL);

  int rank, nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(nproc != 2){
    if(rank == 0){
      for(int i = 0 ; i < number_of_test ; i++)
	printf(NRM "test %d %s" NRM "\n", i, RED"FAIL");
    }
    goto end;
  } 
  
  mesh.dim = 2;
  mesh.nelm_local = 2;
  mesh.nelm_total = 4;
  mesh.eptr = malloc(3*sizeof(int));
  mesh.eptr[0] = 0;
  mesh.eptr[1] = 4;
  mesh.eptr[2] = 8;

  mesh.eind = malloc(4*2*sizeof(int));
  if(rank == 0){
    mesh.eind[0] = 0; mesh.eind[1] = 1; mesh.eind[2] = 4; mesh.eind[3] = 3;
    mesh.eind[4] = 1; mesh.eind[5] = 2; mesh.eind[6] = 4; mesh.eind[7] = 5;
  }else if(rank == 1){
    mesh.eind[0] = 3; mesh.eind[1] = 4; mesh.eind[2] = 6; mesh.eind[3] = 7;
    mesh.eind[4] = 4; mesh.eind[5] = 5; mesh.eind[6] = 7; mesh.eind[7] = 8;
  }

  mesh.coord = malloc(2*9*sizeof(double));
  for(int i = 0 ; i < 3 ; i++){
    for(int j = 0 ; j < 3 ; j++){
      mesh.coord[(i*3 + j)*2 + 0] = i;
      mesh.coord[(i*3 + j)*2 + 1] = j;
    }
  }

  mesh.elm_centroid = malloc(2*2*sizeof(double));
  for(int i = 0 ; i < 2 ; i++){
    for(int d = 0 ; d < 2 ; d++){
      mesh.elm_centroid[i*2 + d] = 0.0;
      for(int n = 0 ; n < 4 ; n++)
	mesh.elm_centroid[i*2 + d] += mesh.coord[(mesh.eind[mesh.eptr[i] + n])*2 + d];
      mesh.elm_centroid[i*2 + d] /= 4;
    }
  }


//  printf(NRM "test 1 %s" NRM "\n", (flag_pass == true) ? GRN"OK" : RED"FAIL");

end:

  MPI_Finalize();
  return 0;
}
