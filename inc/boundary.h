/*
 * Boundary Conditions Structures
 *
 */

#include "list.h"
#include "fun.h"

#ifndef _BOUNDARYH_
#define _BOUNDARYH_

typedef struct boundary_t_{

  char   *name;
  int    typeID;
  int    order;
  int    nfx;
  int    nfy;
  int    nfz;
  f1d_t  *fx; 
  f1d_t  *fy;
  f1d_t  *fz;
  int    GmshID;
  int    NNod;
  int    *MyNod;

}boundary_t;

list_t boundary_list;

#endif
