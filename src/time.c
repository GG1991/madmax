#include "sputnik.h"

int save_time(MPI_Comm *comm, const char *string, FILE *file, double dt)
{

  int    rank, nproc, i, ierr;

  double *time_vec = NULL;

  MPI_Comm_rank(*comm, &rank);
  MPI_Comm_size(*comm, &nproc);

  if(rank==0){
    time_vec = malloc(nproc*sizeof(double));
  }

  ierr = MPI_Gather(&dt, 1, MPI_DOUBLE, time_vec, 1, MPI_DOUBLE, 0, *comm);
  if(ierr){
    return 1;
  }
  if(rank == 0){
    fprintf(file,"%-20s", string);
    for(i=0;i<nproc;i++){
      fprintf(file," %e",time_vec[i]);
    }
    fprintf(file,"\n");
  }

  return 0;
}
