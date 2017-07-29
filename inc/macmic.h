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

#define  FLAG_VTK_NONE 0
#define  FLAG_VTK_PART 1
#define  FLAG_VTK_DISP 2

#define  MPI_MICRO_START 1

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

#define  COMM_NULL        0
#define  COMM_MACRO_MICRO 1

typedef struct CoupMac_1_t_{

  int   my_micro_worker;

}CoupMac_1_t;

typedef struct CoupMic_1_t_{

  int   my_macro_leader;
  int   im_the_micro_leader;

}coupMic_1_t;

typedef struct coupling_t_{

  int   type;
  int   id;
  void  *coup;

}coupling_t;

coupling_t macmic;

/*
   Global Variables
*/

int           flag_print_vtk;
PetscBool     flag_coupling;

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

#endif