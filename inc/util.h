#ifndef UTIL_H
#define UTIL_H

#define ARRAY_SET_TO_ZERO(array, length) {for(int i = 0 ; i < length ; i++) array[i] = 0.0;}
#define ARRAY_COPY(array_1, array_2, length) {for(int i = 0 ; i < length ; i++) array_1[i] = array_2[i];}

#ifdef PETSC
#include "petscksp.h"
int print_petsc_ksp_info( MPI_Comm COMM, KSP ksp);
#endif

int strbin2dec( char *str );

#endif
