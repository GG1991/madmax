#ifndef SOLVERS_H
#define SOLVERS_H


#define SOLVER_NULL 0
#define SOLVER_PETSC 1

typedef struct{

  int type;

}solver_t;

extern solver_t solver;


#endif
