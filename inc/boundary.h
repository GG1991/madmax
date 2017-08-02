/*
   Boundary Conditions Structures
  
   Author > Guido Giuntoli
   Date   > 02-08-2017
 */

#include "list.h"
#include "fun.h"

#ifndef _BOUNDARYH_
#define _BOUNDARYH_

typedef struct mac_boundary_t_{

  int      kind;
  int      order;
  int      nfx, nfy, nfz;
  f1d_t    *fx, *fy, *fz;
  int      nnod, *nods;
  int      ndir_pn, ndir, *dir_idx;
  double   *dir_val;
  int      nneu_pn, nneu, *neu_idx;
  double   *neu_val;

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
