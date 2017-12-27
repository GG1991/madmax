#include "macro.h"

int alloc_memory(void) {

  int ierr;

  if (solver.type == SOLVER_PETSC) {

#ifdef SLEPC
    PETSC_COMM_WORLD = MACRO_COMM;
    SlepcInitialize(&command_line.argc, &command_line.argv, (char*)0, NULL);
#elif  PETSC
    PETSC_COMM_WORLD = MACRO_COMM;
    PetscInitialize(&command_line.argc, &command_line.argv, (char*)0, NULL);
#else
    return 1;
#endif

    if (params.calc_mode == CALC_MODE_EIGEN || params.calc_mode == CALC_MODE_NORMAL) {

      int nnz = dim * MAX_ADJ_NODES;

      ierr = MatCreate(MACRO_COMM, &A);
      ierr = MatSetSizes(A, dim*mesh.nnods_local, dim*mesh.nnods_local, dim*mesh.nnods_total, dim*mesh.nnods_total);
      ierr = MatSetType(A, MATAIJ);
      ierr = MatSeqAIJSetPreallocation(A, nnz, NULL);
      ierr = MatMPIAIJSetPreallocation(A, nnz, NULL, nnz, NULL);
      ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      ierr = MatSetUp(A);
      ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
      ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

      int *ghost_index = malloc(mesh.nnods_ghost * dim *sizeof(int));

      for (int i = 0 ; i < mesh.nnods_ghost ; i++) {
	for (int d = 0 ; d < dim ; d++)
	  ghost_index[i*dim + d] = mesh.local_to_global[mesh.nnods_local + i]*dim + d;
      }

      ierr = VecCreateGhost(MACRO_COMM, dim*mesh.nnods_local, dim*mesh.nnods_total, dim*mesh.nnods_ghost, ghost_index, &x);

      if (params.calc_mode == CALC_MODE_EIGEN) {

	ierr = MatCreate(MACRO_COMM, &M);
	ierr = MatSetSizes(M, dim*mesh.nnods_local, dim*dim*mesh.nnods_local, dim*mesh.nnods_total, dim*mesh.nnods_total);
	ierr = MatSetType(M, MATAIJ);
	ierr = MatSeqAIJSetPreallocation(M, nnz, NULL);
	ierr = MatMPIAIJSetPreallocation(M, nnz, NULL, nnz, NULL);
	ierr = MatSetOption(M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
	ierr = MatSetUp(M);
	ierr = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
	ierr = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);

      }
      else if (params.calc_mode == CALC_MODE_NORMAL) {

	ierr = VecDuplicate(x, &dx);
	ierr = VecDuplicate(x, &b);

      }
    }
  }

  loc_elem_index = malloc(npe_max*dim*sizeof(int));
  glo_elem_index = malloc(npe_max*dim*sizeof(int));
  elem_disp = malloc(npe_max*dim*sizeof(double));
  elem_coor = malloc(npe_max*dim*sizeof(double));
  k_elem = malloc(npe_max*dim*npe_max*dim*sizeof(double));
  m_elem = malloc(npe_max*dim*npe_max*dim*sizeof(double));
  res_elem = malloc(npe_max*dim*sizeof(double));
  stress_gp = malloc(nvoi*sizeof(double));
  strain_gp = malloc(nvoi*sizeof(double));
  c = malloc(nvoi*nvoi*sizeof(double));
  elem_strain = malloc(mesh.nelm_local * nvoi*sizeof(double));
  elem_stress = malloc(mesh.nelm_local * nvoi*sizeof(double));
  elem_energy = malloc(mesh.nelm_local * sizeof(double));
  elem_type = malloc(mesh.nelm_local * sizeof(int));
  flag_neg_detj  = 0;

  bmat = malloc(nvoi*sizeof(double**));
  for (int i = 0 ; i < nvoi  ; i++) {
    bmat[i] = malloc(npe_max*dim*sizeof(double*));
    for (int j = 0 ; j < npe_max*dim ; j++)
      bmat[i][j] = malloc(ngp_max*sizeof(double));
  }

  dsh = malloc(npe_max*sizeof(double**));
  for (int i = 0 ; i < npe_max ; i++) {
    dsh[i] = malloc( dim * sizeof(double*));
    for (int j = 0 ; j < dim ; j++)
      dsh[i][j] = malloc( ngp_max * sizeof(double));
  }

  jac = malloc(dim*sizeof(double*));
  for (int k = 0 ; k < dim ; k++)
    jac[k] = malloc( dim * sizeof(double));

  jac_inv = malloc(dim*sizeof(double*));
  for (int i = 0 ; i < dim ; i++)
    jac_inv[i] = malloc( dim * sizeof(double));

  detj = malloc(ngp_max*sizeof(double));

  if (ierr == 0) flags.allocated = true;

  return ierr;
}
