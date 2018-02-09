#ifndef ELL_H_
#define ELL_H_

#include <stdlib.h> 
#include <stdio.h>

#define NRM  "\x1B[0m"
#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"

typedef struct ell_matrix_ {
  int nrow;
  int ncol;
  int nnz;
  int *cols;
  double *vals;
} ell_matrix;

int ell_init (ell_matrix * m, int nrow, int ncol, int nnz);
int ell_set_val (ell_matrix * m, int row, int col, double val);
int ell_add_val (ell_matrix * m, int row, int col, double val);
int ell_mvp (ell_matrix * m, double *x, double *y);
int ell_print_full (ell_matrix * m);
int ell_search_val (ell_matrix * m, int row, int col, double *val);

#endif
