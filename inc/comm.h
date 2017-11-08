#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>

#ifndef _COMM_H_
#define _COMM_H_

#define  MACRO              1     // MACRO IDs and colors
#define  MICRO              2     // MICRO IDs and colors
#define  COUP_NULL          0
#define  COUP_1             1

#define  MIC_END            1
#define  MAC2MIC_STRAIN     2
#define  C_HOMO             3
#define  RHO                4

int       coup_type;

/*
   COUP_1 > linear analisis only one kind of micro structure with homogeneum material
   ____________________                      __________
   |microstructure     | -> macro-structure  |         | -> micro-structure
   |      _____________|                     | ||||||| |
   |      |homogeneous |                     |         |
   |______|____________|                     |_________|

   We store the constitutive tensor on the <macro> processes, they 
   are calculated only once by their respectively <micro> workers.
   All the calculations should give the same values of <homo_cij>.

*/

MPI_Comm   WORLD_COMM;
MPI_Comm   MICRO_COMM;
MPI_Comm   MACRO_COMM;

int        color;
int        rank_wor;
int        nproc_wor;
bool       flag_coupling;

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

int macmic_coloring(MPI_Comm WORLD_COMM, int *color, coupling_t *macmic, MPI_Comm *LOCAL_COMM);

int mic_recv_signal(MPI_Comm WORLD_COMM, int *signal);
int mic_recv_strain(MPI_Comm WORLD_COMM, double strain[6]);
int mic_recv_macro_gp(MPI_Comm WORLD_COMM, int *macro_gp);
int mic_send_strain(MPI_Comm WORLD_COMM, double strain[6]);
int mic_send_stress(MPI_Comm WORLD_COMM, double stress[6]);
int mic_send_c_homo(MPI_Comm WORLD_COMM, int nvoi, double c_homo[36]);
int mic_send_rho(MPI_Comm WORLD_COMM, double *rho);

int mac_send_signal(MPI_Comm WORLD_COMM, int signal);
int mac_send_strain(MPI_Comm WORLD_COMM, double strain[6]);
int mac_recv_stress(MPI_Comm WORLD_COMM, double stress[6]);
int mac_recv_rho(MPI_Comm WORLD_COMM, double *rho);
int mac_recv_c_homo(MPI_Comm WORLD_COMM, int nvoi, double c_homo[36]);
int mac_send_macro_gp(MPI_Comm WORLD_COMM, int *macro_gp);

#endif