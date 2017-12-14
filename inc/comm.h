#ifndef _COMM_H_
#define _COMM_H_

#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>

#define  COLOR_MACRO 1
#define  COLOR_MICRO 2
#define  MAX_VOIGT 6

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

int comm_init_message(message_t *message);
int comm_macro_send(message_t *message);
int comm_macro_recv(message_t *message);
int comm_micro_send(message_t *message);
int comm_micro_recv(message_t *message);
int comm_finalize_message(void);


#define  COUP_NULL          0
#define  COUP_1             1

int       coup_type;

MPI_Comm   WORLD_COMM;
MPI_Comm   MICRO_COMM;
MPI_Comm   MACRO_COMM;

int        color;
int        rank_wor;
int        nproc_wor;

typedef struct mac_coup_1_t_{

  int     mic_rank;     // rank of micro worker
  double  homo_cij[36]; // tangent constitutive tensor of the micro-structure

}mac_coup_1_t;

typedef struct mic_coup_1_t_{

  int   mac_rank;      // rank of macro leader
  int   im_leader;     // 1 if im the leader 0 if not

}mic_coup_1_t;

typedef struct coupling_t_{

  int   type;
  void  *coup;

}coupling_t;

coupling_t macmic;

int macmic_coloring(MPI_Comm WORLD_COMM, int *color, coupling_t *macmic, MPI_Comm *LOCAL_COMM, bool flag_coupling);

#endif
