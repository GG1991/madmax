/*

   MACRO main function

   Program for solving the displacement field inside a solid 
   structure representing the macrostructure

   Author: Guido Giuntoli

 */


static char help[] = "Solves the displacement field inside a solid structure. \
	    It has the capability of being couple with MICRO, a code for solving and RVE problem.n\n";


#include "macro.h"


int main(int argc, char **argv)
{

    int        ierr;
    bool       set;


    world = MPI_COMM_WORLD;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(world, &nproc_wor);
    ierr = MPI_Comm_rank(world, &rank_wor);
    
    spu_parse_scheme(input_n);

    /* 
       Stablish a new local communicator and a set of 
       intercommunicators with micro programs 
     */
    mac_comm_init();
    ierr = MPI_Comm_size(world, &nproc_mac);
   
    //************************************************************ 
    // Set PETSc communicator to macro_comm
    PETSC_COMM_WORLD = macro_comm;
    ierr = PetscInitialize(&argc,&argv,(char*)0,help);

    //
    // read mesh
    //    
    read_mesh(macro_comm, mesh_n, mesh_f, &elmdist, &eptr, &eind);

    ierr = PetscFinalize();
    ierr = MPI_Finalize();

    return 0;
}
