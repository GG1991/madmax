#ifndef _COMM_H_
#define _COMM_H_


#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>

#define MAX_VOIGT 6

#define ACTION_NULL 0
#define ACTION_MICRO_CALC_STRESS 1
#define ACTION_MICRO_CALC_C_TANGENT 2
#define ACTION_MICRO_CALC_RHO 3
#define ACTION_MICRO_END 4

typedef struct {

  int action;
  int num_voigt;
  double strain_mac[MAX_VOIGT];
  double stress_ave[MAX_VOIGT];
  double c_tangent_ave[MAX_VOIGT*MAX_VOIGT];
  double rho;

}message_t;

message_t message;
MPI_Datatype mpi_message_t;

MPI_Comm   WORLD_COMM;
MPI_Comm   MICRO_COMM;
MPI_Comm   MACRO_COMM;

#define COLOR_MACRO 1
#define COLOR_MICRO 2

typedef struct{

  int color;
  int macro_leader;
  int micro_slave;

}comm_t;

extern comm_t comm;

int comm_init_message(message_t *message);
int comm_macro_send(message_t *message, comm_t *comm);
int comm_macro_recv(message_t *message, comm_t *comm);
int comm_micro_send(message_t *message, comm_t *comm);
int comm_micro_recv(message_t *message, comm_t *comm);
int comm_finalize_message(void);
int comm_coloring(MPI_Comm WORLD_COMM, comm_t *comm, MPI_Comm *LOCAL_COMM);


#endif
