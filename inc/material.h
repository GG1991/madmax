#ifndef MATERIAL_H
#define MATERIAL_H

#include <stdbool.h>
#include "list.h"
#include "myio.h"

#define MAT_ELASTIC    0
#define MAT_MICRO      1

#define MAX_NUM_OF_MATERIALS 4

typedef struct{

  char *name;
  int type_id;
  int id;
  void *type;
  bool is_linear;

}material_t;

typedef struct _type_0{ /* Linear Elastic Material */

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
int material_fill_list_from_command_line(command_line_t *command_line, list_t *material_list);
bool material_are_all_linear(list_t * material_list);
bool material_check_in_list(list_t *material_list, char *name_to_check);

#endif
