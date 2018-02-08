#include "ell.h"

int ell_init (ell_matrix * m, int nrow, int nnz)
{
  if (m == NULL) return 1;
  m->nnz = nnz;
  m->nrow = nrow;
  m->cols = malloc((nrow*nnz) * sizeof(int));
  m->vals = malloc((nrow*nnz) * sizeof(double));
  return 0;
}
