#include "gmsh.h"
#include <stdbool.h>

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

int array_test_5_1proc[] = {100, 28, 1, 29};
int array_test_5_2proc[] = {64, 72, 73, 65};
int array_test_5_3proc[] = {52, 60, 34, 35};

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

  bool flag_pass = true;
  if(rank == 0){

    switch(nproc){

      case 1:
	printf(NRM "test 1 %s" NRM "\n", (gmsh_mesh.num_vol_elems == 81) ? GRN"OK" : RED"FAIL");
	printf(NRM "test 2 %s" NRM "\n", (gmsh_mesh.num_vol_elems_local == 81) ? GRN"OK" : RED"FAIL");
	printf(NRM "test 3 %s" NRM "\n", (gmsh_mesh.num_surf_elems == 39) ? GRN"OK" : RED"FAIL");
	printf(NRM "test 4 %s" NRM "\n", (gmsh_mesh.eptr[1] - gmsh_mesh.eptr[0] == 4) ? GRN"OK" : RED"FAIL");
	for(int i = 0 ; i < 4 ; i++)
	  if(gmsh_mesh.eind[gmsh_mesh.eptr[gmsh_mesh.num_vol_elems_local-1] + i] != array_test_5_1proc[i])
	    flag_pass = false;
	printf(NRM "test 5 %s" NRM "\n", (flag_pass == true) ? GRN"OK" : RED"FAIL");
	break;

      case 2:
	printf(NRM "test 1 %s" NRM "\n", (gmsh_mesh.num_vol_elems == 81) ? GRN"OK" : RED"FAIL");
	printf(NRM "test 2 %s" NRM "\n", (gmsh_mesh.num_vol_elems_local == 41) ? GRN"OK" : RED"FAIL");
	printf(NRM "test 3 %s" NRM "\n", (gmsh_mesh.num_surf_elems == 39) ? GRN"OK" : RED"FAIL");
	printf(NRM "test 4 %s" NRM "\n", (gmsh_mesh.eptr[1] - gmsh_mesh.eptr[0] == 4) ? GRN"OK" : RED"FAIL");
	for(int i = 0 ; i < 4 ; i++)
	  if(gmsh_mesh.eind[gmsh_mesh.eptr[gmsh_mesh.num_vol_elems_local-1] + i] != array_test_5_2proc[i])
	    flag_pass = false;
	printf(NRM "test 5 %s" NRM "\n", (flag_pass == true) ? GRN"OK" : RED"FAIL");
	break;

      case 3:
	printf(NRM "test 1 %s" NRM "\n", (gmsh_mesh.num_vol_elems == 81) ? GRN"OK" : RED"FAIL");
	printf(NRM "test 2 %s" NRM "\n", (gmsh_mesh.num_vol_elems_local == 27) ? GRN"OK" : RED"FAIL");
	printf(NRM "test 3 %s" NRM "\n", (gmsh_mesh.num_surf_elems == 39) ? GRN"OK" : RED"FAIL");
	printf(NRM "test 4 %s" NRM "\n", (gmsh_mesh.eptr[1] - gmsh_mesh.eptr[0] == 4) ? GRN"OK" : RED"FAIL");
	for(int i = 0 ; i < 4 ; i++)
	  if(gmsh_mesh.eind[gmsh_mesh.eptr[gmsh_mesh.num_vol_elems_local-1] + i] != array_test_5_3proc[i])
	    flag_pass = false;
	printf(NRM "test 5 %s" NRM "\n", (flag_pass == true) ? GRN"OK" : RED"FAIL");
	break;
    }
  }

end:

  MPI_Finalize();
  return 0;
}
