#include "ell.h"

int ell_init (ell_matrix * m, int nrow, int ncol, int nnz)
{
  if (m == NULL) return 1;
  m->nnz = nnz;
  m->nrow = nrow;
  m->ncol = ncol;
  m->cols = malloc((nrow*nnz) * sizeof(int));
  m->vals = malloc((nrow*nnz) * sizeof(double));
  if (m->vals == NULL || m->cols == NULL) return 2;
  for (int i = 0 ; i < (nrow*nnz) ; i++) m->cols[i] = -1;
  for (int i = 0 ; i < (nrow*nnz) ; i++) m->vals[i] = +0;
  return 0;
}

int ell_set_val (ell_matrix * m, int row, int col, double val)
{
  if (row >= m->nrow || col >= m->ncol) {
    printf(RED "ell error: row %d or col %d greater than the dimension of the matrix\n" NRM, row, col);
    return 1;
  }
  if (row < 0 || col < 0) {
    printf(RED "ell error: negative values in row %d or col %d\n" NRM, row, col);
    return 2;
  }
  int j = 0;
  while (j < m->nnz) {
    if (m->cols[(row*m->nnz) + j] == -1) {
      m->cols[(row*m->nnz) + j] = col;
      m->vals[(row*m->nnz) + j] = val;
      return 0;
    } else if (m->cols[(row*m->nnz) + j] == col) {
      m->vals[(row*m->nnz) + j] = val;
      return 0;
    }
    j++;
  }
  if (j == m->nnz) {
    printf(RED "ell error: not enought space to store value in row %d and col %d\n" NRM, row, col);
    return 3;
  }
  return 4;
}

int ell_add_val (ell_matrix * m, int row, int col, double val)
{
  if (row >= m->nrow || col >= m->ncol) {
    printf(RED "ell error: row %d or col %d greater than the dimension of the matrix\n" NRM, row, col);
    return 1;
  }
  if (row < 0 || col < 0) {
    printf(RED "ell error: negative values in row %d or col %d\n" NRM, row, col);
    return 2;
  }
  int j = 0;
  while (j < m->nnz) {
    if (m->cols[(row*m->nnz) + j] == -1) {
      m->cols[(row*m->nnz) + j] = col;
      m->vals[(row*m->nnz) + j] = val;
      return 0;
    } else if (m->cols[(row*m->nnz) + j] == col) {
      m->vals[(row*m->nnz) + j] += val;
      return 0;
    }
    j++;
  }
  if (j == m->nnz) {
    printf(RED "ell error: not enought space to add value in row %d and col %d\n" NRM, row, col);
    return 3;
  }
  return 4;
}

int ell_mvp (ell_matrix * m, double *x, double *y)
{
  if (x == NULL || y == NULL) return 1;
  for (int i = 0 ; i < m->nrow ; i++) {
    y[i] = 0;
    int j = 0;
    while (j < m->nnz) {
      if (m->cols[(i*m->nnz) + j] == -1) break;
      y[i] += m->vals[(i*m->nnz) + j] * x[m->cols[(i*m->nnz) + j]];
      j++;
    }
  }
  return 0;
}

int ell_print_full (ell_matrix * m)
{
  if (m == NULL) return 1;
  if (m->vals == NULL || m->cols == NULL) return 2;
  double val;
  for (int i = 0 ; i < m->nrow ; i++) {
    for (int j = 0 ; j < m->ncol ; j++) {
      printf("%lf%s",(ell_search_val(m, i, j, &val) == 0)?val:0.0, (j == m->ncol - 1) ? "\n" : " ");
    }
  }
  return 0;
}

int ell_search_val (ell_matrix * m, int row, int col, double *val)
{
  if (row >= m->nrow || col >= m->ncol) return 1;
  int j = 0;
  while (j < m->nnz) {
    if (m->cols[(row * m->nnz) + j] == col) break;
    j++;
  }
  if (j == m->nnz) return -1;
  *val = m->vals[(row * m->nnz) + j];
  return 0;
}
