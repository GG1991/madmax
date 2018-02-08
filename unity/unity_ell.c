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

  return 0;
}
