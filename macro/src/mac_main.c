/*

   MACRO main function

   Program for solving the displacement field inside a solid 
   structure representing the macrostructure

 */


static char help[] = "Solves the displacement field inside a solid structure. \
	    It has the capability of being couple with MICRO, a code for solving and RVE problem.n\n";


#include "macro.h"


int main(int argc, char **argv)
{

    char       mesh_n[MESH_N_LENGTH]; // Mesh file name
    char       mesh_f[4];             // Mesh format name
    int        nproc_tot;             // number of total process 
    int        nproc_mac;              
    int        ierr;
    bool       couple_fl;
    bool       set;

    MPI_Comm     macro_comm;
    MPI_Comm   * macmic_comm;
    MPI_Comm     world;
    spu_comm_t   spu_comm;

    world = MPI_COMM_WORLD;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(world, &nproc_tot);
    ierr = MPI_Comm_rank(world, &rank);
//    ierr = PetscPrintf(PETSC_COMM_WORLD,"sizeof(world)  : %d\n",sizeof(MPI_Comm));CHKERRQ(ierr);
    
    parse_mpi("mpi.dat", &spu_comm);

    /* stablish a new local communicator and a set of intercommunicators 
       with micro programs 
     */
    mac_color(world, spu_comm, &macro_comm, &macmic_comm);
    ierr = MPI_Comm_size(world, &nproc_mac);
   
    // Set PETSc communicator to macro_comm
    PETSC_COMM_WORLD = macro_comm;
    ierr = PetscInitialize(&argc,&argv,(char*)0,help);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"processors  : %d\n",nproc_mac);CHKERRQ(ierr);
    //
    // read options from command line 
    //    
    // couple flag ( 1:coupled 0:not coupled) -> -c 1
    //
    ierr = PetscOptionsGetBool(NULL,NULL,"-c",(PetscBool*)&couple_fl,(PetscBool*)&set);CHKERRQ(ierr);
    if(set == PETSC_FALSE){
	couple_fl = PETSC_FALSE;
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"coupling    : %d\n",couple_fl);CHKERRQ(ierr);

    //
    // mesh file format  -> -mfor gmsh
    //
    ierr = PetscOptionsGetString(NULL,NULL,"-mfor",mesh_f,16,(PetscBool*)&set);CHKERRQ(ierr);
    if(set == PETSC_FALSE){
	strcpy(mesh_f,"gmsh");
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"mesh format : %s\n",mesh_f);CHKERRQ(ierr);

    //
    // mesh file name    -> -mesh rve/cube_unif/cube.msh
    //
    ierr = PetscOptionsGetString(NULL,NULL,"-mesh",mesh_n,MESH_N_LENGTH,(PetscBool*)&set);CHKERRQ(ierr);
    if(set == PETSC_FALSE){
	strcpy(mesh_n,"../meshes/cube_unif/cube.msh");
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"mesh file   : %s\n",mesh_n);CHKERRQ(ierr);

    //
    // read mesh
    //    
    read_mesh(macro_comm, mesh_n, mesh_f, &elmdist, &eptr, &eind);

    ierr = PetscFinalize();
    ierr = MPI_Finalize();

    return 0;
}
