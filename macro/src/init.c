#include "macro.h"

void init_variables(void) {

  MPI_Init(&command_line.argc, &command_line.argv);

  WORLD_COMM = MPI_COMM_WORLD;
  MACRO_COMM = MPI_COMM_WORLD;

  MPI_Comm_size(WORLD_COMM, &nproc_wor);
  MPI_Comm_rank(WORLD_COMM, &rank_wor);

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
  flags.print_pvtu = false;
  flags.print_vectors = false;
  flags.print_matrices = false;

  mesh.partition = PARMETIS_GEOM;

  comm_init_message(&message);

  A = NULL;
  b = NULL;
  x = NULL;
  dx = NULL;

  return;
}
