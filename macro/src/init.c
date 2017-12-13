#include "macro.h"

void init_variables(void){

  params.calc_mode = CALC_MODE_NULL;
  params.non_linear_max_its = 1;
  params.non_linear_min_norm_tol = 1.0;
  params.num_eigen_vals = 1;
  params.tf = 0.0;
  params.dt = 0.0;
  params.t = 0.0;
  params.ts = 0;
  params.energy_stored = 1.0;

  solver.type = SOLVER_PETSC;

  flags.coupled = false;
  flags.allocated = false;
  flags.print_vectors = false;
  flags.print_matrices = false;

  comm_init_message(&message);

  A = NULL;
  b = NULL;
  x = NULL;
  dx = NULL;

  return;
}
