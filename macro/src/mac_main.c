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

    int        i, ierr;
    char       *myname = strdup("macro");

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

    print_flag = false;
    for(i=2;i<argc;i++){
      if(strcmp(argv[i],"-p")==0){
	  print_flag = true;
      }
    }

    spu_parse_scheme(input_n);
    
    /* 
       Stablish a new local communicator and a set of 
       intercommunicators with micro programs 
     */
    mac_comm_init();
    
    // file for measuring time in "main" rutine
    time_fl = fopen("time_mac.dat","w");

    if(rank_mac == 0){
      fprintf(time_fl, "%-20s", "ranks");
      for(i=0;i<nproc_mac;i++){
	fprintf(time_fl," %-12d",i);
      }
      fprintf(time_fl,"\n");
    }

    // here we collect time from all processes
    if(rank_mac == 0){
      time_vec = calloc(nproc_mac, sizeof(double));
    }

    t0 = MPI_Wtime();      /* ON time lapse */

    spu_parse_mesh(input_n);

    t1 = MPI_Wtime() - t0;
    ierr = MPI_Gather(&t1, 1, MPI_DOUBLE, time_vec, 1, MPI_DOUBLE, 0, macro_comm);
    if(ierr){
      return 1;
    }
    if(rank_mac == 0){
      fprintf(time_fl, "%-20s", "spu_parse_mesh");
      for(i=0;i<nproc_mac;i++){
	fprintf(time_fl," %e",time_vec[i]);
      }
      fprintf(time_fl,"\n");
    }                       /* OFF time lapse */
    

    //************************************************************ 
    // Set PETSc communicator to macro_comm
    PETSC_COMM_WORLD = macro_comm;
    ierr = PetscInitialize(&argc,&argv,(char*)0,help);

    //
    // read mesh
    //    
    t0 = MPI_Wtime();      /* ON time lapse */

    strcpy(mesh_f,"gmsh");
    read_mesh_elmv(&macro_comm, myname, mesh_n, mesh_f);

    t1 = MPI_Wtime() - t0;
    ierr = MPI_Gather(&t1, 1, MPI_DOUBLE, time_vec, 1, MPI_DOUBLE, 0, macro_comm);
    if(ierr){
      return 1;
    }
    if(rank_mac == 0){
      fprintf(time_fl, "%-20s", "read_mesh");
      for(i=0;i<nproc_mac;i++){
	fprintf(time_fl," %e",time_vec[i]);
      }
      fprintf(time_fl,"\n");
    }                       /* OFF time lapse */

    nelm = elmdist[rank_mac+1] - elmdist[rank_mac];
    part = (int*)malloc(nelm * sizeof(int));

    t0 = MPI_Wtime();      /* ON time lapse */

    // partition the mesh
    part_mesh_PARMETIS(&macro_comm, time_fl, myname, NULL, PARMETIS_MESHKWAY );

    // We delete repeated nodes and save the <nnod_glo> values on <nod_glo> in order
    clean_vector_qsort(&macro_comm, myname, eptr[nelm], eind, &nod_glo, &nnod_glo);

    // calculate <*ghosts> and <nghosts> 
    calculate_ghosts(&macro_comm, myname);

    t1 = MPI_Wtime() - t0;
    ierr = MPI_Gather(&t1, 1, MPI_DOUBLE, time_vec, 1, MPI_DOUBLE, 0, macro_comm);
    if(ierr){
      return 1;
    }
    if(rank_mac == 0){
      fprintf(time_fl, "%-20s", "part_mesh_PARMETIS");
      for(i=0;i<nproc_mac;i++){
	fprintf(time_fl," %e",time_vec[i]);
      }
      fprintf(time_fl,"\n");
    }                       /* OFF time lapse */
    
    spu_vtk_partition( myname, mesh_n, &macro_comm );

    fclose(time_fl);

    ierr = PetscFinalize();
    ierr = MPI_Finalize();

    return 0;
}
