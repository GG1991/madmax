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

    int        i, ierr;
    char       *myname = strdup("micro");

    world_comm = MPI_COMM_WORLD;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(world_comm, &nproc_wor);
    ierr = MPI_Comm_rank(world_comm, &rank_wor);
    
    if(argc>1){
      strcpy(input_n,argv[1]);
    }
    else{
       printf("mic_main.c:no input file has been given\n");
    }

    spu_parse_scheme(input_n);

    /* 
       Stablish a new local communicator and a set of 
       intercommunicators with micro programs 
     */
    mic_comm_init();

    // file for measuring time in "main" rutine
    time_fl = fopen("time_mic.dat","w");

    if(rank_mic == 0){
      fprintf(time_fl, "%-20s", "ranks");
      for(i=0;i<nproc_mic;i++){
	fprintf(time_fl," %-12d",i);
      }
      fprintf(time_fl,"\n");
    }

    // here we collect time from all processes
    if(rank_mic == 0){
      time_vec = calloc(nproc_mic, sizeof(double));
    }

    t0 = MPI_Wtime();
    spu_parse_mesh(input_n);
    t1 = MPI_Wtime() - t0;
    ierr = MPI_Gather(&t1, 1, MPI_DOUBLE, time_vec, 1, MPI_DOUBLE, 0, micro_comm);
    if(ierr){
      return 1;
    }
    if(rank_mic == 0){
      fprintf(time_fl, "%-20s", "spu_parse_mesh");
      for(i=0;i<nproc_mic;i++){
	fprintf(time_fl," %e",time_vec[i]);
      }
      fprintf(time_fl,"\n");
    }
    

    //************************************************************ 
    // Set PETSc communicator to micro_comm
    PETSC_COMM_WORLD = micro_comm;
    ierr = PetscInitialize(&argc,&argv,(char*)0,help);

    //
    // read mesh
    //    
    t0 = MPI_Wtime();

    strcpy(mesh_f,"gmsh");
    read_mesh_elmv(&micro_comm, myname, mesh_n, mesh_f);

    t1 = MPI_Wtime() - t0;
    ierr = MPI_Gather(&t1, 1, MPI_DOUBLE, time_vec, 1, MPI_DOUBLE, 0, micro_comm);
    if(ierr){
      return 1;
    }
    if(rank_mic == 0){
      fprintf(time_fl, "%-20s", "read_mesh");
      for(i=0;i<nproc_mic;i++){
	fprintf(time_fl," %e",time_vec[i]);
      }
      fprintf(time_fl,"\n");
    }

    nelm = elmdist[rank_mic+1] - elmdist[rank_mic];
    part = (int*)malloc(nelm * sizeof(int));

    t0 = MPI_Wtime();
    part_mesh_PARMETIS(&micro_comm, time_fl, myname, NULL, PARMETIS_MESHKWAY );
    t1 = MPI_Wtime() - t0;
    ierr = MPI_Gather(&t1, 1, MPI_DOUBLE, time_vec, 1, MPI_DOUBLE, 0, micro_comm);
    if(ierr){
      return 1;
    }
    if(rank_mic == 0){
      fprintf(time_fl, "%-20s", "part_mesh_PARMETIS");
      for(i=0;i<nproc_mic;i++){
	fprintf(time_fl," %e",time_vec[i]);
      }
      fprintf(time_fl,"\n");
    }

    fclose(time_fl);

    ierr = PetscFinalize();
    ierr = MPI_Finalize();

    return 0;
}
