#include "micro.h"

int alloc_memory(void) {

  int ierr = 0;

  if (params.homog_method == HOMOG_METHOD_UNIF_STRAINS) {

    if (solver.type == SOLVER_PETSC) {

#ifdef PETSC
      PETSC_COMM_WORLD = MICRO_COMM;
      PetscInitialize(&command_line.argc, &command_line.argv, (char*)0, NULL);

      int nnz = (mesh_struct.dim == 2)? 18:81;

      MatCreate(MICRO_COMM,&A);
      MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, mesh_struct.nn * mesh_struct.dim, mesh_struct.nn * mesh_struct.dim);
      MatSetFromOptions(A);
      MatSeqAIJSetPreallocation(A, nnz, NULL);
      MatMPIAIJSetPreallocation(A, nnz, NULL, nnz, NULL);
      MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

      VecCreate(MICRO_COMM, &x);
      VecSetSizes(x, PETSC_DECIDE, mesh_struct.nn * mesh_struct.dim);
      VecSetFromOptions(x);
      VecDuplicate(x, &dx);
      VecDuplicate(x, &b);

      KSPCreate(PETSC_COMM_WORLD, &ksp);
      KSPSetFromOptions(ksp);

#else
      return 1;
#endif

    }
  }

  elem_index = malloc(dim*mesh_struct.npe*sizeof(int));
  elem_nods = malloc(mesh_struct.npe*sizeof(int));
  elem_disp = malloc(dim*mesh_struct.npe*sizeof(double));
  stress_gp = malloc(nvoi*sizeof(double));
  strain_gp = malloc(nvoi*sizeof(double));
  elem_strain = malloc(mesh_struct.nelm*nvoi*sizeof(double));
  elem_stress = malloc(mesh_struct.nelm*nvoi*sizeof(double));
  elem_energy = malloc(mesh_struct.nelm*sizeof(double));
  elem_type = malloc(mesh_struct.nelm*sizeof(int));

  struct_bmat = malloc(nvoi*sizeof(double**));
  for (int i = 0 ; i < nvoi ; i++) {
    struct_bmat[i] = malloc(mesh_struct.npe*dim*sizeof(double*));
    for (int j = 0 ; j < mesh_struct.npe*dim ; j++)
      struct_bmat[i][j] = malloc(ngp*sizeof(double));
  }

  flags.allocated = true;

  return ierr;
}
