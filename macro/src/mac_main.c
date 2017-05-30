
/*
   Peuge main function

   Program for solving the displacement field inside a solid structure


 */


static char help[] = "Solves the displacement field inside a solid structure. \
	    It has the capability of being couple with Giant, a code for solving and RVE problem.n\n";


#include "peuge.h"


int main(int argc, char **argv)
{

    char       mesh_n[64];    // Mesh file name
    char       mesh_f[4];     // Mesh format name
    int        color;         // color of this process
    int        nproc_mac;     // number of macroscopic process (used to determine color)
    int        nsubs_mac;     // number of macroscopic subdomains to performe the partition
    int        nkind_mic;     // number of microscopic kinds 
    int      * mpinfo_mic;    // array specifying how many process for each micro-kind are going to be used
    int        ierr;

    bool       set;
    MPI_Comm   world = MPI_COMM_WORLD;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(world, &nproc);
    ierr = MPI_Comm_rank(world, &rank);
    
    parse_mpi("mpi.dat", &nproc_mac, &nsubs_mac, &nkind_mic, &mpinfo_mic);

    /* We change our PETSc communicator, macroprocess will solve a distributed problem
       will all microkinds per macroprocess will solve another one 
     */
    
//    PETSC_COMM_WORLD = macro;
//
//    ierr = PetscInitialize(&argc,&args,(char*)0,help);
//
//    ierr = PetscPrintf(PETSC_COMM_WORLD,"processors  : %d\n",nproc);CHKERRQ(ierr);

    // read options from command line 
    //    
    // couple flag ( 1:coupled 0:not coupled) -> -c 1
    //
//    ierr = PetscOptionsGetBool(NULL,NULL,"-c",&couple_fl,&set);CHKERRQ(ierr);
//    if(set == PETSC_FALSE){
//	couple_fl = PETSC_FALSE;
//    }
//    ierr = PetscPrintf(PETSC_COMM_WORLD,"coupling    : %d\n",couple_fl);CHKERRQ(ierr);

    //
    // mesh file format  -> -mfor gmsh
    //
//    ierr = PetscOptionsGetString(NULL,NULL,"-mfor",mesh_f,16,&set);CHKERRQ(ierr);
//    if(set == PETSC_FALSE){
//	strcpy(mesh_f,"gmsh");
//    }
//    ierr = PetscPrintf(PETSC_COMM_WORLD,"mesh format : %s\n",mesh_f);CHKERRQ(ierr);

    //
    // mesh file name    -> -mesh rve/cube_unif/cube.msh
    //
//    ierr = PetscOptionsGetString(NULL,NULL,"-mesh",mesh_n,16,&set);CHKERRQ(ierr);
//    if(set == PETSC_FALSE){
//	strcpy(mesh_n,"rve/cube_unif/cube.msh");
//    }
//    ierr = PetscPrintf(PETSC_COMM_WORLD,"mesh file   : %s\n",mesh_n);CHKERRQ(ierr);
//
//
    // read mesh
    //    
    //    ierr = peu_rmsh(mesh_n, mesh_f);
//
//    ierr = PetscFinalize();
    ierr = MPI_Finalize();

    return 0;
}
