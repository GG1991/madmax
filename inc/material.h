/*
   Materials definition types
  
   material_t material;
   lambda = ((type_0*)material.type)->lambda
  
 */

#include "list.h"

#define TYPE_0    0
#define TYPE01    1

#define MICRO  2

#ifndef _MATERIALH_
#define _MATERIALH_

typedef struct material_t_{

  char  *name;
  int   type_id;
  int   GmshID;
  void  *type;

}material_t;

/*
   Linear Elastic Material
 */
typedef struct _type_0{

  double young;
  double poisson;
  double lambda;
  double mu;
  double rho;

}type_0;

list_t material_list;

#endif
