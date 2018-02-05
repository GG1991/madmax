#include "micro.h"

int init_variables_1(void)
{
  params.multis_method = MULTIS_NULL;
  params.fe2_bc = BC_NULL;
  params.non_linear_max_its = 2;
  params.non_linear_min_norm_tol = 1.0e-4;
  params.rho = 1.0e7;

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

  A = NULL;
  b = NULL;
  x = NULL;
  dx = NULL;

  return 0;
}

int init_variables_2(void)
{

  for (int gp = 0; gp < ngp ; gp++) {
    for (int is = 0; is < mesh_struct.npe ; is++) {
      if (dim == 2) {
	struct_bmat[0][is*dim + 0][gp] = struct_dsh[is][0][gp];
	struct_bmat[0][is*dim + 1][gp] = 0;
	struct_bmat[1][is*dim + 0][gp] = 0;
	struct_bmat[1][is*dim + 1][gp] = struct_dsh[is][1][gp];
	struct_bmat[2][is*dim + 0][gp] = struct_dsh[is][1][gp];
	struct_bmat[2][is*dim + 1][gp] = struct_dsh[is][0][gp];
      }
    }
  }

  int ierr = 0;
  if (material_are_all_linear(&material_list) == true){
    MIC_PRINTF_0("calc ctan around for linear micro-structure...\n")
    ierr = homog_calculate_c_tangent_around_zero(params.c_tangent_linear);
    flags.linear_materials = true;
    flags.c_linear_calculated = true;
  }

  return ierr;
}
