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

#define FLAG_VTK_NONE 0
#define FLAG_VTK_PART 1
#define FLAG_VTK_DISP 2

#define MACMIC_START  1
#define MACMIC_END    2
#define SIGNAL_NULL         -1
#define SIGNAL_MICRO_CALC    3
#define SIGNAL_MICRO_END     4

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

typedef struct coupMac_1_t_{

  int   mic_rank; /* rank of micro worker */

}coupMac_1_t;

typedef struct coupMic_1_t_{

  int   mac_rank; /* rank of macro leader */
  int   im_leader;/* 1 if im the leader 0 if not */

}coupMic_1_t;

typedef struct coupling_t_{

  int   type;
  void  *coup;

}coupling_t;

coupling_t macmic;

/*
   Global Variables
*/
MPI_Comm     MICRO_COMM;
MPI_Comm     MACRO_COMM;

int          *remote_ranks;     //  remote ranks if micro processes

int          nev;

MPI_Comm     WORLD_COMM;
int          color;
int          rank_wor;                 //  rank on world comm
int          nproc_wor;                //  # of processes (WORLD_COMM)

int          nstruc_mic;               // number of micro structures
int          *nproc_per_mic;           // number of processes per micro structure ( size = nstruc_mic )
int          nproc_mic_group;          // number of micro process in a group = sum_i nproc_per_mic[i]
int          nmic_worlds;              // number of micro worlds nproc_mic / nproc_mic_group
int          scheme;                   // communication approach

PetscBool    flag_coupling;
int          flag_print_vtk;
// Matrices and vectors

Mat           A;                    /* Jacobian Matrix          */
Vec           x, dx, b;             /* Vectors unknowns and RHS */
KSP           ksp;                  /* linear solver context    */
KSPConvergedReason  reason;

double        *stress, *strain;     // Averange strain and stress on each element

/*
   Common functions
*/

int MacMicInitGaussStructure(int *eptr, int nelm);
int MacMicParseScheme( char *input );
int MacMicColoring(MPI_Comm WORLD_COMM, int *color, coupling_t *macmic, MPI_Comm *LOCAL_COMM);
int MicCommWaitSignal( MPI_Comm WORLD_COMM, int *signal );
int MicCommWaitStartSignal( MPI_Comm WORLD_COMM );
int MicCommRecvStrain( MPI_Comm WORLD_COMM );
int MicCommRecvGPnum( MPI_Comm WORLD_COMM );
int MicCommSendAveStressAndTanTensor( MPI_Comm WORLD_COMM );
int MicCommSendAveStress( MPI_Comm WORLD_COMM );
int MicCommSendAveTTensor( MPI_Comm WORLD_COMM );

int MacCommSendSignal( MPI_Comm WORLD_COMM, int signal );

#endif
