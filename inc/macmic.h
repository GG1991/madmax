/*

   Header to define data structures for <macro> & <micro> 
   programs

   Author: Guido Giuntoli
   Date: 28-07-2017

*/

#include <stdlib.h>

#ifndef MACMIC_H
#define MACMIC_H

#define  FLAG_VTK_NONE 0
#define  FLAG_VTK_PART 1
#define  FLAG_VTK_DISP 2

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

/*
   Global Variables
*/
int           flag_print_vtk;

/*
   Common functions
*/

int MacMicInitGaussStructure(int *eptr, int nelm);

#endif
