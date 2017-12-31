#include "micro.h"

int finalize(void)
{
  int ierr;

  if (flags.allocated == true) {

    free(elem_index);
    free(elem_nods);
    free(elem_disp);
    free(stress_gp);
    free(strain_gp);
    free(elem_strain);
    free(elem_stress);
    free(elem_energy);
    free(elem_type);

    for (int i = 0; i < nvoi; i++) {
      for (int j = 0; j < mesh_struct.npe*dim; j++)
	free(struct_bmat[i][j]);
      free(struct_bmat[i]);
    }
    free(struct_bmat);

    if (params.homog_method == HOMOG_METHOD_UNIF_STRAINS) {
      if (solver.type == SOLVER_PETSC) {

	ierr = PetscFinalize();
      }
    }
  }

  comm_finalize_message();

  ierr = MPI_Finalize();

  return ierr;
}
