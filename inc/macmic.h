/*

   Header to define data structures for <macro> & <micro> 
   programs

   Author: Guido Giuntoli
   Date: 28-07-2017

*/

#include <stdlib.h>
#include "petscksp.h"

#ifndef MACMIC_H
#define MACMIC_H

#define MACRO         1     // MACRO IDs and colors
#define MICRO         2     // MICRO IDs and colors

#define FLAG_VTK_NONE 0
#define FLAG_VTK_PART 1
#define FLAG_VTK_DISP 2

#define MPI_MICRO_START 1

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

  int   myMicWorker;

}coupMac_1_t;

typedef struct coupMic_1_t_{

  int   myMacLeader;
  int   imMicLeader;

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

int          nproc_mac;         //  # of macro processes (WORLD_COMM)  
int          nproc_mic;         //  # of micro processes (MICRO_COMM)
int          rank_mic;          //  rank on macro comm
int          rank_mac;          //  rank on macro comm
int          *remote_ranks;     //  remote ranks if micro processes

int          nev;

MPI_Comm     WORLD_COMM;
int          rank_wor;                 //  rank on world comm
int          nproc_wor;                //  # of processes (WORLD_COMM)
int          flag_print_vtk;
int          *id_vec;                  // ID vector size = #proc (info of which ID has each rank
int          nproc_mac_tot;            // number of macro processes total (inside WORLD_COMM)
int          nproc_mic_tot;            // number of micro processes total (inside WORLD_COMM)
int          nstruc_mic;               // number of micro structures
int          *nproc_per_mic;           // number of processes per micro structure ( size = nstruc_mic )
int          nproc_mic_group;          // number of micro process in a group = sum_i nproc_per_mic[i]
int          nmic_worlds;              // number of micro worlds nproc_mic / nproc_mic_group
int          scheme;                   // communication approach
PetscBool    flag_coupling;

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
int MacMicColoring(int id);

#endif
