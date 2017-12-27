#ifndef SOLVERS_H
#define SOLVERS_H

#define SOLVER_NULL 0
#define SOLVER_PETSC 1

typedef struct{

  int type;

}solver_t;

extern solver_t solver;

#ifdef PETSC
#include "petscksp.h"
int solvers_print_petsc_ksp_info(MPI_Comm COMM, KSP ksp);
#endif

#ifdef SLEPC
#include "slepceps.h"
#endif

#endif
