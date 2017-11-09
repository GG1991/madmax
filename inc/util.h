#ifndef UTIL_H
#define UTIL_H

#ifdef PETSC
#include "petscksp.h"
int print_ksp_info( MPI_Comm COMM, KSP ksp);
#endif

int strbin2dec( char *str );

#endif
