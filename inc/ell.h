#include <stdlib.h> 

#ifndef ELL_H_
#define ELL_H_

typedef struct ell_matrix_ {
  int nrow;
  int nnz;
  int *cols;
  double *vals;
} ell_matrix;

int ell_init (ell_matrix * m, int nrow, int nnz);

#endif
