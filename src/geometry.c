#include "geometry.h"

int geometry_2d_line_side(double n_line[2], double p_line[2], double point[2])
{
  if (n_line == NULL || p_line == NULL || point == NULL) return -2;

  double side = 0.0;
  for (int i = 0; i < 2; i++)
    side += n_line[i] * (p_line[i] - point[i]);

  if (fabs(side) < GEOM_TOL)
    return  0;

  else if (side > 0)
    return  1;

  else
    return -1;

  return -2;
}
