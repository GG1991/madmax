#include "ell.h"
#include <stdio.h>

#define NRM  "\x1B[0m"
#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"

int main(void)
{
  printf("ell.c / ell.h test\n");

  double x[4] = {0.029412, 0.029412, 0.029412, 0.029412};
  double y[4] = {1.0, 2.0, 3.0, 4.0};
  double x1[4] = {0.0, 0.0, 0.0, 0.0};
  double b[4] = {1.0, 1.0, 1.0, 1.0};
  ell_matrix m;
  ell_solver solver;
  ell_init(&m, 4, 4, 4);
  ell_set_val(&m, 0, 0, 1); ell_set_val(&m, 0, 1, 0); ell_set_val(&m, 0, 2, 0); ell_set_val(&m, 0, 3, 0);
  ell_set_val(&m, 1, 0,-1); ell_set_val(&m, 1, 1, 2); ell_set_val(&m, 1, 2,-1); ell_set_val(&m, 1, 3, 0);
  ell_set_val(&m, 2, 0, 0); ell_set_val(&m, 2, 1,-1); ell_set_val(&m, 2, 2, 2); ell_set_val(&m, 2, 3,-1);
  ell_set_val(&m, 3, 0, 0); ell_set_val(&m, 3, 1, 0); ell_set_val(&m, 3, 2, 0); ell_set_val(&m, 3, 3, 1);
  printf("\nm=\n");
  ell_print_full(&m);
  ell_mvp(&m, x, y);
  ell_solve_jacobi(&solver, &m, b, x1);

  printf("\ny=\n");
  for (int i = 0 ; i < 4 ; i++) {
    printf("%lf\n", y[i]);
  }

  printf("\nx1=\n");
  for (int i = 0 ; i < 4 ; i++) {
    printf("%lf\n", x1[i]);
  }
  printf("\nerror = %lf\n", solver.err);
  printf("its = %d\n", solver.its);

  return 0;
}
