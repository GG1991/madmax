#include "ell.h"
#include <stdio.h>

#define NRM  "\x1B[0m"
#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"

int main(void)
{
  printf("ell.c / ell.h test\n");

  double x[4] = {1.0, 2.0, 3.0, 4.0};
  double y[4] = {1.0, 2.0, 3.0, 4.0};
  ell_matrix m;
  ell_init(&m, 4, 4, 4);
  ell_set_val(&m, 1, 1, 2.71);
  ell_set_val(&m, 1, 2, 2.69);
  ell_set_val(&m, 1, 0, 2.69);
  ell_set_val(&m, 1, 3, 2.69);
  ell_add_val(&m, 1, 1,-2.00);
  ell_add_val(&m, 3, 1,-2.00);
  ell_print_full(&m);
  ell_mvp(&m, x, y);

  printf("\n");
  for (int i = 0 ; i < 4 ; i++) {
    printf("%lf\n", y[i]);
  }

  return 0;
}
