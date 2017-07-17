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

    spu_parse_scheme( input_n );
    SpuParseMaterials( input_n );
    
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

    /******************/
    /* ON time lapse */
    t0 = MPI_Wtime();
    spu_parse_mesh(input_n);
    t1 = MPI_Wtime() - t0;
    save_time(&MACRO_COMM, "spu_parse_mesh", time_fl, t1);
    /* OFF time lapse */
    /******************/
    

    //************************************************************ 
    // Set PETSc communicator to MACRO_COMM
    PETSC_COMM_WORLD = MACRO_COMM;
    ierr = PetscInitialize(&argc,&argv,(char*)0,help);

    //
    // read mesh
    //    
    /******************/
    /* ON time lapse */
    t0 = MPI_Wtime();
    strcpy(mesh_f,"gmsh");
    read_mesh_elmv(&MACRO_COMM, myname, mesh_n, mesh_f);
    t1 = MPI_Wtime() - t0;
    save_time(&MACRO_COMM, "read_mesh", time_fl, t1);
    /* OFF time lapse */
    /******************/

    nelm = elmdist[rank_mac+1] - elmdist[rank_mac];
    part = (int*)malloc(nelm * sizeof(int));

    // partition the mesh

    /******************/
    /* ON time lapse */
    t0 = MPI_Wtime();
    part_mesh_PARMETIS(&MACRO_COMM, time_fl, myname, NULL, PARMETIS_MESHKWAY );
    t1 = MPI_Wtime() - t0;
    save_time(&MACRO_COMM, "part_mesh_PARMETIS", time_fl, t1);
    /* OFF time lapse */
    /******************/

    // We delete repeated nodes and save the <NAllMyNod> values on <AllMyNodOrig> in order
    /******************/
    /* ON time lapse */
    clean_vector_qsort(&MACRO_COMM, myname, eptr[nelm], eind, &AllMyNodOrig, &NAllMyNod);
    t1 = MPI_Wtime() - t0;
    save_time(&MACRO_COMM, "AllMyNodOrig calc", time_fl, t1);
    /* OFF time lapse */
    /******************/

    // calculate <*ghosts> and <nghosts> 

    /******************/
    /* ON time lapse */
    t0 = MPI_Wtime();
    calculate_ghosts(&MACRO_COMM, myname);
    t1 = MPI_Wtime() - t0;
    save_time(&MACRO_COMM, "ghosts", time_fl, t1);
    /* OFF time lapse */
    /******************/

    // reenumerate nodes

    /******************/
    /* ON time lapse */
    t0 = MPI_Wtime();
    reenumerate_PETSc(&MACRO_COMM);
    t1 = MPI_Wtime() - t0;
    save_time(&MACRO_COMM, "reenumerate", time_fl, t1);
    /* OFF time lapse */
    /******************/

    /******************/
    /* ON time lapse */
    t0 = MPI_Wtime();
    read_mesh_coord(&MACRO_COMM, myname, mesh_n, mesh_f);
    t1 = MPI_Wtime() - t0;
    save_time(&MACRO_COMM, "read coord", time_fl, t1);
    /* OFF time lapse */
    /******************/

    /******************/
    /* ON time lapse */
    t0 = MPI_Wtime();
    char  vtkfile_n[NBUF];

    sprintf(vtkfile_n,"%s_part_%d.vtk",myname,rank_mac);
    spu_vtk_partition( vtkfile_n, &MACRO_COMM );
    t1 = MPI_Wtime() - t0;
    save_time(&MACRO_COMM, "vtk_partition", time_fl, t1);
    /* OFF time lapse */
    /******************/

    AllocMatrixVector( MACRO_COMM, NMyNod*3, NTotalNod*3, &A, &x, &b);
    /*
       Currently, all PETSc parallel matrix formats are partitioned by
       contiguous chunks of rows across the processors.  Determine which
       rows of the matrix are locally owned.
     */
    int Istart, Iend;
    ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
    if( Istart != StartIndexRank[rank_mac]*3 ){
      printf("AllocMatrixVector: error on indeces set for matrix and vector.\n");
      return 1;
    }
    if(rank_mac<nproc_mac-1){
      if( Iend != StartIndexRank[rank_mac+1]*3 ){
	printf("AllocMatrixVector: error on indeces set for matrix and vector.\n");
	return 1;
      }
    }
    else{
      if( Iend != NTotalNod*3 ){
	printf("AllocMatrixVector: error on indeces set for matrix and vector.\n");
	return 1;
      }
    }


    fclose(time_fl);

    ierr = PetscFinalize();
    ierr = MPI_Finalize();

    return 0;
}
