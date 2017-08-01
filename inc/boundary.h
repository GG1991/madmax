/*
 * Boundary Conditions Structures
 *
 */

#include "list.h"
#include "fun.h"

#ifndef _BOUNDARYH_
#define _BOUNDARYH_

typedef struct mac_boundary_t_{

  int      kind;
  int      order;
  int      nfx;
  int      nfy;
  int      nfz;
  f1d_t    *fx; 
  f1d_t    *fy;
  f1d_t    *fz;
  int      GmshID;
  int      NNods;
  int      *Nods;

  /* Staff to set boundary conditions */
  int      NDirPerNode;
  int      NNeuPerNode;
  int      *indeces;
  int      NDirIndeces;
  int      *DirichletIndeces;
  int      NNeuIndeces;
  int      *NeumannIndeces;
  double   *values;
  double   *DirichletValues; 
  double   *NeumannValues;

}mac_boundary_t;

typedef struct mic_boundary_t_{

  void * voidv;

}mic_boundary_t;

/* 
   All we want to know about a generic boundary is its 
   GmshId, name and the nodes that it has
*/
typedef struct boundary_t_{

  char    *name;
  int     GmshID;
  list_t  Nods;
  void    *bvoid;  //this can be <mac_boundary_t> or <mic_boundary_t>

}boundary_t;

list_t boundary_list;      // it is used to store things particular for each problem

#endif
