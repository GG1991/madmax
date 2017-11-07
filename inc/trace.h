#include <mpi.h>

#ifndef _TRACE_H_
#define _TRACE_H_

typedef struct _trace_t{

  FILE   *file;
  double t0;

}trace_t;


trace_t trace;

int init_trace( MPI_Comm COMM, const char *file_name );
int save_event( MPI_Comm COMM, const char *event );
int end_trace ( MPI_Comm COMM );

#endif
