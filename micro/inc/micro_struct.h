#ifndef MICRO_STRUCT_H
#define MICRO_STRUCT_H

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

#define FIBER_CILIN  1
#define FIBER_PLANAR 2

#define MATRIX_ID 0
#define FIBER_ID  1

/*
   represents a square with a circular fibers in the middle

   from command line this should be defined for 2d

   -micro_struct "fiber_cilin <size[0]> <size[1]> <nx> <ny> <radio> <desv[0]> <desv[1]>"

*/

typedef struct _fiber_cilin_t
{

  double  *size;        // size of the rve
  double  *desv;        // desviation of the fiber from the center
  double   radio;       // radius of the fiber
  int      nx;          // number of fibers in x
  int      ny;          // number of fibers in y

}fiber_cilin_t;


/*
   represents a square with a planar fibers in the middle

   from command line this should be defined for 2d

   -micro_struct "fiber_planar <size[0]> <size[1]> <ntype> <angle[0]> <seps[0]> <width[0]> <desv[0]> <numb[0]>
                                                           <angle[1]> <seps[1]> <width[1]> <desv[1]> <numb[1]>
      						     ..."

 */

typedef struct _fiber_planar_t
{

  double  *size;         // size of the rve
  int      ntype;        // number of types 
  double  *angles;       // angle for each type ( 1 value in 2d (theta), 2 values in 3d (theta and phy))
  double  *seps;         // separation for each type
  double  *width;        // width for each type
  double  *desv;         // desviation from the center for each type ( 1 value in 2d and 3d ("y" displacement) )
  int     *numb;         // total number of fibers for each type

}fiber_planar_t;


/* micro structure */

typedef struct _micro_struct_t
{

  int      type;
  void    *data;

}micro_struct_t;

int init_micro_structure( int dim, micro_struct_t * micro_struct_t, char * format );
int get_elem_id( int dim, micro_struct_t *micro_struct_t, double *elem_centroid, int *elem_id );

#endif
