#include "micro.h"

void init_variables(void){

  params.homog_method = HOMOG_METHOD_NULL;
  params.non_linear_max_its = 2;
  params.non_linear_min_norm_tol = 1.0e-4;

  solver.type = SOLVER_PETSC;

  flags.coupled = false;
  flags.linear_materials = false;
  flags.allocated = false;
  flags.c_linear_calculated = false;
  flags.print_pvtu = false;
  flags.print_vectors = false;
  flags.print_matrices = false;
  
  MPI_Init(&command_line.argc, &command_line.argv);

  WORLD_COMM = MPI_COMM_WORLD;
  MICRO_COMM = MPI_COMM_WORLD;

  comm_init_message(&message);

  A   = NULL;
  b   = NULL;
  x   = NULL;
  dx  = NULL;

  return;
}
