/*
 * Functions to calculate geometrical properties
 *
 *
 * author : Guido Giuntoli
 * date   : 26 - 11 - 2017
 *
 */

#include "geometry.h"

int geom_2d_line_side( const double n_line[2], const double p_line[2], const double point[2] )
{
  /* Determines the side in where a point is respect to a line
   *
   * The side is defined respect to the sign of the dot product
   * between the normal to line "n_line" and the vector that join
   * the point "point" with a point that belongs to the line
   * "p_line". side = sign ( < n_line , point - p_line > )
   *
   *
   * returns -2 is error
   * returns -1 for side 1
   * returns  1 for side 2
   * returns  0 if the point is in the line
   */

  if( !n_line || !p_line || !point ) return -2;

  int    i;
  double side = 0;

  for( i=0 ; i<2 ; i++ )
    side += n_line[i] * (p_line[i] - point[i]);

  if( fabs(side) < GEOM_TOL ) return  0; /* is in the line ? */
  else if( side > 0 )         return  1;
  else                        return -1;

  return -2; 
}
