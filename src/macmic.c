/*
   Common functions for <macro> & <micro> programs
*/

#include "macmic.h"

int MacMicInitGaussStructure(int *eptr, int nelm)
{
  /*
     Allocates memory for Gauss point global structure <gauss>
     of type <gauss_t>
  */
  int i, ngp, ngauss = 0;
  for(i=1;i<(nelm+1);i++){
    // we assume that each element has the same number of gauss points as vertices 
    ngp = eptr[i] - eptr[i-1]; 
    ngauss += ngp; 
  }
  gauss = malloc(ngauss * sizeof(gauss_t)); if(!gauss) return 1;

  for(i=0;i<ngauss;i++){
    gauss[i].param_d = NULL;
  }

  return 0;
}
