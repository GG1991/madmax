/*

   MICRO main function

   Program for solving the microscopic problem 
   for multi-scale approach.

   Author: Guido Giuntoli

 */


static char help[] = "Solves an RVE problem.\n\n";



#include "micro.h"


int main(int argc, char **argv)
{

    int        ierr;

    world_comm = MPI_COMM_WORLD;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(world_comm, &nproc_wor);
    ierr = MPI_Comm_rank(world_comm, &rank_wor);
    
    if(argc>1){
      strcpy(input_n,argv[1]);
    }
    else{
       printf("mac_main.c:no input file has been given\n");
    }

    spu_parse_scheme(input_n);

    /* 
       Stablish a new local communicator and a set of 
       intercommunicators with micro programs 
     */
    mic_comm_init();
    	
    spu_parse_mesh(input_n);
   
    //************************************************************ 
    // Set PETSc communicator to macro_comm
    PETSC_COMM_WORLD = micro_comm;
    ierr = PetscInitialize(&argc,&argv,(char*)0,help);

    //
    // read mesh
    //    
    strcpy(mesh_f,"gmsh");
    read_mesh(&micro_comm, mesh_n, mesh_f, &elmdist, &eptr, &eind);
    nelm = elmdist[rank_mic+1] - elmdist[rank_mic];
    part = (int*)malloc(nelm * sizeof(int));
    part_mesh_PARMETIS(&micro_comm, elmdist, eptr, eind, part, NULL, PARMETIS_MESHKWAY );

    ierr = PetscFinalize();
    ierr = MPI_Finalize();

    return 0;
}
