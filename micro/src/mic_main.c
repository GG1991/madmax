/*

   MICRO main function

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
   
    //************************************************************ 
    // Set PETSc communicator to macro_comm
    PETSC_COMM_WORLD = macro_comm;
    ierr = PetscInitialize(&argc,&argv,(char*)0,help);

    //
    // read mesh
    //    
//    read_mesh(macro_comm, mesh_n, mesh_f, &elmdist, &eptr, &eind);

    ierr = PetscFinalize();
    ierr = MPI_Finalize();

    return 0;
}
