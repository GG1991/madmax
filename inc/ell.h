#include <stdlib.h> 
#include <stdio.h>

#ifndef ELL_H_
#define ELL_H_

typedef struct ell_matrix_ {
  int nrow;
  int ncol;
  int nnz;
  int *cols;
  double *vals;
} ell_matrix;

int ell_init (ell_matrix * m, int nrow, int ncol, int nnz);
int ell_print_full (ell_matrix * m);
int ell_search_val (ell_matrix * m, int row, int col, double *val);

#endif
