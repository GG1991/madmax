/*
 * Materials definition types
 *
 * material_t material;
 * lambda = ((type_00*)material.type)->lambda
 *
 */

#include "list.h"

#ifndef _MATERIALH_
#define _MATERIALH_

typedef struct material_t_{

  int   typeID;
  void  *type;

}material_t;

/*
 * Linear Elastic Material
 *
 */
typedef struct _type_00{

  double young;
  double poisson;
  double lambda;
  double mu;

}type_00;

list_t material_list;

#endif
