#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define ARRAY_SET_TO_ZERO(array, length) {for(int i = 0 ; i < length ; i++) array[i] = 0.0;}
#define ARRAY_COPY(array_1, array_2, length) {for(int i = 0 ; i < length ; i++) array_1[i] = array_2[i];}

int strbin2dec(char *str);
int util_is_in_vector(int val, int *vector, int size);
int util_clean_and_sort_vector(int *in_vector, int n_in, int **out_vector, int *n_out);
int util_cmpfunc(const void * a, const void * b);
int util_sort_vector_intersec(int *array_1, int n1, int *array_2, int n2, int **inter_arr, int *n_inter);

#endif
