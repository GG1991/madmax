/* Materials for solid structures calculation
 * 
 */

#ifndef MATERIAL_H
#define MATERIAL_H

#include "list.h"
#include "myio.h"

#define MAT_ELASTIC    0
#define MAT_MICRO      1

#define MAX_NUM_OF_MATERIALS 4

typedef struct material_t_{

  char  *name;
  int   type_id;
  int   id;
  void  *type;

}material_t;

/* Linear Elastic Material */
typedef struct _type_0{

  double young;
  double poisson;
  double lambda;
  double mu;
  double rho;

}type_0;

list_t material_list;

int material_get_stress(material_t *mat, int dim, double *strain, double *stress);
int material_get_c_tang(material_t *mat, int dim, double *strain, double *c_tang);
int material_get_rho(material_t *mat, int dim, double *rho);
int material_are_linear(list_t * material_list);
int material_fill_list_from_command_line(int argc, char **argv, list_t *material_list);

#endif
