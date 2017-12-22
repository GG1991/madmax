#include "util.h"
#include <stdbool.h>

#define NRM  "\x1B[0m"
#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"

#define PRINT_ARRAY(array, length) {\
  for(int i = 0 ; i < length ; i++){ \
    printf("%d ", array[i]); \
  } \
  printf("\n"); \
}

#define COMP_ARRAY(arr_1, len_1, arr_2, len_2, flag) {\
  flag = true;\
  if(len_1 != len_2){\
    flag = false;\
  }else{\
    int i = 0;\
    while (i < len_1){\
      if(arr_1[i] != arr_2[i]){\
	flag = false;\
	break;\
      }\
      i++;\
    }\
  }\
}

int array_in_1[] = {1, 3, 4, 4, 3, 1, 1, 6, 5};
int array_out_1[] = {1, 3, 4, 5, 6};

int main(void){

  int number_of_test = 2;
  bool flag_pass = true;

  int *array_out, n_out;
  util_clean_and_sort_vector(array_in_1, 9, &array_out, &n_out);
  printf(NRM "test 1 %s" NRM "\n", (n_out == 5) ? GRN"OK" : RED"FAIL");
  COMP_ARRAY(array_out_1, 5, array_out, n_out, flag_pass);
  printf(NRM "test 1 %s" NRM "\n", (flag_pass == true) ? GRN"OK" : RED"FAIL");

  return 0;
}
