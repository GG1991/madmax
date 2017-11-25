#ifndef MICRO_STRUCT_H
#define MICRO_STRUCT_H

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

#define FIBER_CILIN  1
#define FIBER_PLANAR 2

#define ID_MATRIX 0
#define ID_FIBER  1

/*
   represents a square with a circular fibers in the middle

   from command line this should be defined for 2d

   -micro_struct "fiber_cilin <size[0]> <size[1]> <nx> <ny> <radio> <desv[0]> <desv[1]>"

*/

typedef struct _fiber_cilin_t
{

  double  *desv;        // desviation of the fiber from the center
  double   radio;       // radius of the fiber
  int      nx_fib;      // number of fibers in x
  int      ny_fib;      // number of fibers in y

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

  int      type;         // micro structure type
  double  *size;         // size of the rve
  void    *data;         // data of the micro structure

}micro_struct_t;

int micro_struct_init( int dim, micro_struct_t *micro_struct, char *format );
int micro_struct_get_elem_id( int dim, micro_struct_t *micro_struct, double *elem_centroid, int *elem_id );

int micro_struct_init_elem_type(
    micro_struct_t *micro_struct,
    int dim,
    int nelm,
    int (*get_centroid)( int e, int dim, double *elem_centroid ),
    int *elem_type );

#endif
