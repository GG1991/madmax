#ifndef MICRO_STRUCT_H
#define MICRO_STRUCT_H

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "geom.h"

#define FIBER_CILIN  1
#define FIBER_LINE   2

#define ID_MATRIX 0
#define ID_FIBER  1

/* Represents a square with a circular fibers in the middle
 *
 * note : for 2d and 3d
 */

typedef struct _fiber_cilin_t
{

  double  *desv;         // desviation of the fiber from the center
  double   radio;        // radius of the fiber
  int      nx_fib;       // number of fibers in x
  int      ny_fib;       // number of fibers in y

}fiber_cilin_t;


/* Represents a square with a linear fibers
 *
 * note : only for 2d
 */

typedef struct _fiber_line_t
{

  int      ntype;        // number of types 
  int     *nfib;         // fiber number               (for each type)
  double  *theta;        // fiber angle                (for each type)
  double  *seps;         // fiber separation           (for each type)
  double  *width;        // fiber width                (for each type)
  double  *desv;         // fiber central desviation   (for each type)

}fiber_line_t;


/* Main micro-structure for representing the arrange of material on the domain */

typedef struct _micro_struct_t
{

  int      type;         // micro structure type
  double  *size;         // size of the rve
  void    *data;         // data of the micro structure

}micro_struct_t;


int micro_struct_init( int dim, const char *string,
    micro_struct_t *micro_struct );

int micro_struct_get_elem_id( int dim, micro_struct_t *micro_struct,
    double *elem_centroid, int *elem_id );

int micro_struct_init_elem_type( micro_struct_t *micro_struct,
    int dim, int nelm,
    int (*get_centroid)( int e, int dim, double *elem_centroid ),
    int *elem_type );

#endif
