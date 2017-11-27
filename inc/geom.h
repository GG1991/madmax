#ifndef GEOM_H
#define GEOM_H

#define GEOM_TOL 1.0e-10

#include <math.h>

int geom_2d_line_side( const double n_line[2], const double p_line[2], const double point[2] );

#endif
