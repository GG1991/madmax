/*
   Materials definition types
  
   material_t material;
   lambda = ((type_00*)material.type)->lambda
  
 */

#include "list.h"

#define TYPE00    0
#define TYPE01    1

#define MICRO  2

#ifndef _MATERIALH_
#define _MATERIALH_

typedef struct material_t_{

  char  *name;
  int   typeID;
  int   GmshID;
  void  *type;

}material_t;

/*
   Linear Elastic Material
 */
typedef struct _type_00{

  double young;
  double poisson;
  double lambda;
  double mu;
  double rho;

}type_00;

list_t material_list;

#endif
