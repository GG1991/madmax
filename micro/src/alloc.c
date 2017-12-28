#include "micro.h"

int alloc_memory(void) {

  int ierr = 0;

  if (solver.type == SOLVER_PETSC) {

#ifdef PETSC
    PETSC_COMM_WORLD = MICRO_COMM;
    PetscInitialize(&command_line.argc, &command_line.argv, (char*)0, NULL);
#else
    return 1;
#endif

    int nnz = (dim==2)? 18:81;

    MatCreate(MICRO_COMM,&A);
    MatSetSizes(A, nl*dim, nl*dim, nn*dim, nn*dim);
    MatSetFromOptions(A);
    MatSeqAIJSetPreallocation(A, nnz, NULL);
    MatMPIAIJSetPreallocation(A, nnz, NULL, nnz, NULL);
    MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  }

  loc_elem_index = malloc(dim*npe*sizeof(int));
  glo_elem_index = malloc(dim*npe*sizeof(int));
  elem_disp = malloc(dim*npe*sizeof(double));
  stress_gp = malloc(nvoi*sizeof(double));
  strain_gp = malloc(nvoi*sizeof(double));
  elem_strain = malloc(nelm*nvoi*sizeof(double));
  elem_stress = malloc(nelm*nvoi*sizeof(double));
  elem_energy = malloc(nelm*sizeof(double));
  elem_type = malloc(nelm*sizeof(int));

  struct_bmat = malloc(nvoi*sizeof(double**));
  for (int i = 0 ; i < nvoi ; i++) {
    struct_bmat[i] = malloc(npe*dim*sizeof(double*));
    for (int j = 0 ; j < npe*dim ; j++)
      struct_bmat[i][j] = malloc(ngp*sizeof(double));
  }

  int *ghost_index;
  if (nproc_mic == 1)
    ngho = 0;
  else{
    if (rank_mic == 0)
      ngho = 0;
    else
      ngho = ( dim == 2 ) ? nx : nx*nz;
  }

  ghost_index = malloc( ngho*dim *sizeof(int) );

  int i;
  for ( i = 0 ; i < ngho*dim  ; i++ )
    ghost_index[i] = istart - (( dim == 2 )? nx : nx*nz)*dim + i;

  VecCreateGhost(MICRO_COMM, nl*dim, nn*dim, ngho*dim, ghost_index, &x);
  VecDuplicate(x,&dx);
  VecDuplicate(x,&b);

  free(ghost_index);

  /* alloc arrays for boundary condition setting */
  if (nproc_mic == 1) {
    ndir_ix = 2*nx + 2*(nyl-2);
  }
  else{
    ndir_ix = ( rank_mic == 0 || rank_mic == (nproc_mic-1) ) ? (nx + (nyl-1)*2) : (nyl * 2);
  }
  ndir_ix    *= dim;
  dir_ix_loc  = malloc ( ndir_ix * sizeof(int));
  dir_ix_glo  = malloc ( ndir_ix * sizeof(int));
  coor_dir    = malloc ( ndir_ix * sizeof(double));

  /* fill the array of dirichlet indeces "dir_ix_loc" and the coordinates "coor_dir" */
  int n, d, c = 0;
  if (rank_mic == 0) {

    /* y = 0 */
    for ( n = 0 ; n < nx ; n++ ) {
      for ( d = 0 ; d < dim ; d++ )
	dir_ix_loc[c*dim + d] = n*dim + d;
      coor_dir[c*dim + 0] = n*hx;
      coor_dir[c*dim + 1] = 0.0;
      c++;
    }
  }
  if ( rank_mic == (nproc_mic - 1) ) {

    /* y = ly */
    for ( n = 0 ; n < nx ; n++ ) {
      for ( d = 0 ; d < dim ; d++ )
	dir_ix_loc[c*dim + d] = ((nyl-1)*nx + n)*dim + d;
      coor_dir[c*dim + 0] = n*hx;
      coor_dir[c*dim + 1] = ly;
      c++;
    }
  }
  if ( nproc_mic > 1)
  {
    if ( rank_mic == 0 )
    {
      /* x = 0 */
      for ( n = 0 ; n < (nyl - 1) ; n++ ) {
	for ( d = 0 ; d < dim ; d++ )
	  dir_ix_loc[c*dim + d] = (n+1)*nx*dim + d;
	coor_dir[c*dim + 0] = 0;
	coor_dir[c*dim + 1] = (ny_inf + n+1)*hy;
	c++;
      }

      /* x = lx */
      for ( n = 0 ; n < (nyl - 1) ; n++ ) {
	for ( d = 0 ; d < dim ; d++ )
	  dir_ix_loc[c*dim + d] = (2*nx-1)*dim + n*nx*dim + d;
	coor_dir[c*dim + 0] = lx;
	coor_dir[c*dim + 1] = (ny_inf + n + 1)*hy;
	c++;
      }
    }
    else if ( rank_mic == (nproc_mic - 1) )
    {
      /* x = 0 */
      for ( n = 0 ; n < (nyl - 1) ; n++ ) {
	for ( d = 0 ; d < dim ; d++ )
	  dir_ix_loc[c*dim + d] = n*nx*dim + d;
	coor_dir[c*dim + 0] = 0;
	coor_dir[c*dim + 1] = ( ny_inf + n )*hy;
	c++;
      }

      /* x = lx */
      for ( n = 0 ; n < (nyl - 1) ; n++ ) {
	for ( d = 0 ; d < dim ; d++ )
	  dir_ix_loc[c*dim + d] = (nx-1)*dim + n*nx*dim + d;
	coor_dir[c*dim + 0] = lx;
	coor_dir[c*dim + 1] = ( ny_inf + n )*hy;
	c++;
      }
    }
    else
    {
      /* 0 < rank_mic < (nproc - 1) */
      /* x = 0 */
      for ( n = 0 ; n < nyl ; n++ ) {
	for ( d = 0 ; d < dim ; d++ )
	  dir_ix_loc[c*dim + d] = n*nx*dim + d;
	coor_dir[c*dim + 0] = 0;
	coor_dir[c*dim + 1] = ( ny_inf + n )*hy;
	c++;
      }

      /* x = lx */
      for ( n = 0 ; n < nyl ; n++ ) {
	for ( d = 0 ; d < dim ; d++ )
	  dir_ix_loc[c*dim + d] = (nx-1)*dim + n*nx*dim + d;
	coor_dir[c*dim + 0] = lx;
	coor_dir[c*dim + 1] = ( ny_inf + n )*hy;
	c++;
      }
    }
  }
  else{

    /* x = 0 */
    for ( n = 0 ; n < nyl-2; n++ ) {
      for ( d = 0 ; d < dim ; d++ )
	dir_ix_loc[c*dim + d] = (n+1)*nx*dim + d;
      coor_dir[c*dim + 0] = 0;
      coor_dir[c*dim + 1] = (ny_inf + n+1)*hy;
      c++;
    }

    /* x = lx */
    for ( n = 0 ; n < nyl-2; n++ ) {
      for ( d = 0 ; d < dim ; d++ )
	dir_ix_loc[c*dim + d] = ((n+2)*nx-1)*dim + d;
      coor_dir[c*dim + 0] = lx;
      coor_dir[c*dim + 1] = (ny_inf + n + 1)*hy;
      c++;
    }

  }

  for ( i = 0; i < ndir_ix ; i++ )
    dir_ix_glo[i] = local_to_global_index( dir_ix_loc[i] );

  flags.allocated = true;

  return ierr;
}
