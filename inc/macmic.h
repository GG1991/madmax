/*

   Header to define data structures for <macro> & <micro> 
   programs

   Author: Guido Giuntoli
   Date: 28-07-2017

*/

#include <stdlib.h>
#include <stdbool.h>
#include "petscksp.h"

#ifndef MACMIC_H
#define MACMIC_H

#define MACRO         1     // MACRO IDs and colors
#define MICRO         2     // MICRO IDs and colors

#define MIC_END            1
#define MAC2MIC_STRAIN     2

/*
   This structure represents a Gauss points
   in this case it has an element called 
   <param_d> to store those variable that
   should be represented by double precision
 */

typedef struct gauss_t_{

  double *param_d;

}gauss_t;

gauss_t * gauss;

#define  COUP_NULL  0
#define  COUP_1     1

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

/*
   Global Variables
*/

int          *remote_ranks;     //  remote ranks if micro processes
int          color;

int          nstruc_mic;        // number of micro structures
int          *nproc_per_mic;    // number of processes per micro structure ( size = nstruc_mic )
int          nproc_mic_group;   // number of micro process in a group = sum_i nproc_per_mic[i]
int          nmic_worlds;       // number of micro worlds nproc_mic / nproc_mic_group
int          scheme;            // communication approach

PetscBool    flag_coupling;

// Matrices and vectors

Mat           A;                    // Jacobian Matrix          
Vec           x, dx, b;             // Vectors unknowns and RHS 
KSP           ksp;                  // linear solver context    
int           reason;

double  *stress, *strain, *energy;  // Averange strain, stress and energy on each element

/*
   Common functions
*/

int macmic_coloring(MPI_Comm WORLD_COMM, int *color, coupling_t *macmic, MPI_Comm *LOCAL_COMM);
int mic_recv_signal(MPI_Comm WORLD_COMM, int *signal);
int mic_recv_strain(MPI_Comm WORLD_COMM, double strain[6]);
int mic_send_strain(MPI_Comm WORLD_COMM, double strain[6]);
int mic_send_stress(MPI_Comm WORLD_COMM, double stress[6]);
int mic_sent_ttensor(MPI_Comm WORLD_COMM, double ttensor[36]);

int mac_send_signal(MPI_Comm WORLD_COMM, int signal);
int mac_send_strain(MPI_Comm WORLD_COMM, double strain[6]);
int mac_recv_stress(MPI_Comm WORLD_COMM, double stress[6]);

int mac_calc_homo_cij(double homo_cij[36]);

#endif
