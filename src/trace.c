#include "trace.h"


int init_trace(MPI_Comm COMM, const char *file_name) {

  int rank, size;
  MPI_Comm_rank(COMM, &rank);
  MPI_Comm_size(COMM, &size);

  trace.t0 = MPI_Wtime();

  if (rank == 0) {
    trace.file = fopen(file_name , "w");
    fprintf(trace.file, "%-10s", "rank");
    for (int i = 0 ; i < size ; i++)
      fprintf(trace.file, "%-8d ", i);
    fprintf(trace.file, "\n");
  }

  return 0;

}


int save_event(MPI_Comm COMM, const char *event) {

  int ierr, rank, size;
  MPI_Comm_rank( COMM, &rank );
  MPI_Comm_size( COMM, &size );

  double time = MPI_Wtime();

  if (rank == 0) {

    double *times = malloc( size * sizeof(double));
    ierr = MPI_Gather( &time, 1, MPI_DOUBLE, times, 1, MPI_DOUBLE, 0, COMM); if (ierr) return 1;
    fprintf( trace.file, "%-10s", event);
    for (int i = 0 ; i < size ; i++)
      fprintf( trace.file, "%1.4lf   ", (times[i] - trace.t0));
    fprintf( trace.file, "\n");
  }
  else
    ierr = MPI_Gather(&time, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0, COMM); if (ierr) return 1;

  return 0;
}


int end_trace(MPI_Comm COMM) {

  int rank;
  MPI_Comm_rank(COMM, &rank);

  if (rank == 0)
    fclose(trace.file);

  return 0;
}
