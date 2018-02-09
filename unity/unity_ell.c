#include "ell.h"
#include <stdio.h>

#define NRM  "\x1B[0m"
#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"

int main(void)
{
  printf("ell.c / ell.h test\n");

  ell_matrix m;
  ell_init(&m, 4, 4, 4);
  ell_set_val(&m, 1, 1, 2.71);
  ell_set_val(&m, 1, 2, 2.69);
  ell_set_val(&m, 1, 0, 2.69);
  ell_set_val(&m, 1, 3, 2.69);
  ell_set_val(&m, 1, 0, 1.00);
  ell_set_val(&m, 3, 3, 3.14);
  ell_add_val(&m, 1, 1,-2.00);
  ell_print_full(&m);

  return 0;
}
