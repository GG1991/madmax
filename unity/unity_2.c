#include "gmsh.h"

#define NRM  "\x1B[0m"
#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"

#define PRINT_ARRAY(array, length) {\
if( rank == 0 ){ \
for(int i = 0 ; i < length ; i++){ \
    printf("%d ", array[i]); \
  } \
  printf("\n"); \
}\
}


gmsh_mesh_t gmsh_mesh;

int main(void){

  const char file_name[] = "data/cube_2d.msh";
  int number_of_test = 2;

  MPI_Init(NULL, NULL);

  int rank, nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(nproc > 3){
    if(rank == 0){
      for(int i = 0 ; i < number_of_test ; i++)
	printf(NRM "test %d %s" NRM "\n", i, RED"FAIL");
    }
    goto end;
  } 

  gmsh_mesh.dim = 2;
  gmsh_read_vol_elms_csr_format_parall(MPI_COMM_WORLD, file_name, &gmsh_mesh);
  if(rank == 0){
    printf(NRM "test 1 %s" NRM "\n", (gmsh_mesh.num_vol_elems == 81) ? GRN"OK" : RED"FAIL");
    switch(nproc){
      case 1:
	printf(NRM "test 2 %s" NRM "\n", (gmsh_mesh.num_vol_elems_local == 81) ? GRN"OK" : RED"FAIL");
	break;
      case 2:
	printf(NRM "test 2 %s" NRM "\n", (gmsh_mesh.num_vol_elems_local == 41) ? GRN"OK" : RED"FAIL");
	break;
      case 3:
	printf(NRM "test 2 %s" NRM "\n", (gmsh_mesh.num_vol_elems_local == 27) ? GRN"OK" : RED"FAIL");
	break;
    }
  }

  PRINT_ARRAY(gmsh_mesh.elem_per_proc, nproc);

//  TEST("$Elements", 1, 4202)
//  TEST("$Nodes", 2, 152)
//  TEST("$PhysicalNames", 3, 35)

end:

  MPI_Finalize();
  return 0;
}
