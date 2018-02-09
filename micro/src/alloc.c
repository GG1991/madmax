#include "micro.h"

int alloc_memory(void) {

  int ierr = 0;
  int nnz = (mesh_struct.dim == 2)? 18:81;

  if (params.multis_method == MULTIS_FE2) {

    int nrows = 0;
    switch (params.fe2_bc)  {
      case BC_USTRAIN:
	nrows = mesh_struct.nn * mesh_struct.dim;
	break;
      case BC_PERIODIC:
	nrows = (mesh_struct.nn + (mesh_struct.nx - 2) + (mesh_struct.ny - 2)) * mesh_struct.dim;
	break;
      case BC_USTRESS:
	nrows = mesh_struct.nn*mesh_struct.dim + nvoi;
	break;
    }

    if (params.solver == SOL_PETSC) {

#ifdef PETSC
      PETSC_COMM_WORLD = MICRO_COMM;
      PetscInitialize(&command_line.argc, &command_line.argv, (char*)0, NULL);

      MatCreate(MICRO_COMM,&A);
      MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nrows, nrows);
      MatSetFromOptions(A);
      MatSeqAIJSetPreallocation(A, nnz, NULL);
      MatMPIAIJSetPreallocation(A, nnz, NULL, nnz, NULL);
      MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

      VecCreate(MICRO_COMM, &x);
      VecSetSizes(x, PETSC_DECIDE, nrows);
      VecSetFromOptions(x);
      VecDuplicate(x, &dx);
      VecDuplicate(x, &b);

      KSPCreate(PETSC_COMM_WORLD, &ksp);
      KSPSetFromOptions(ksp);

#else
      return 1;
#endif

    } else if (params.solver == SOL_ELL) {

      x_ell = malloc(nrows*sizeof(double));
      dx_ell = malloc(nrows*sizeof(double));
      b_ell = malloc(nrows*sizeof(double));
      ell_init(&jac_ell, nrows, nrows, nnz);
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
