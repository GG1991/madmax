#ifndef UTIL_H
#define UTIL_H

#define ARRAY_SET_TO_ZERO(array, length) {for(int i = 0 ; i < length ; i++) array[i] = 0.0;}
#define ARRAY_COPY(array_1, array_2, length) {for(int i = 0 ; i < length ; i++) array_1[i] = array_2[i];}

#ifdef PETSC
#include "petscksp.h"
int print_petsc_ksp_info( MPI_Comm COMM, KSP ksp);
#endif

int strbin2dec(char *str);
int util_is_in_vector(int val, int *vector, int size);
int util_clean_and_sort_vector(int *in_vector, int n_in, int **out_vector, int *n_out);
int util_cmpfunc(const void * a, const void * b);

#endif
