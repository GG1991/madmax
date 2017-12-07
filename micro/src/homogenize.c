#include "micro.h"


int mic_homogenize_taylor(double *strain_mac, double *strain_ave, double *stress_ave){
    
  int ierr;
  int ne_i = 0, ne_m = 0;
  double *c_i = malloc( nvoi*nvoi * sizeof(double));
  double *c_m = malloc( nvoi*nvoi * sizeof(double));
  double *c   = malloc( nvoi*nvoi * sizeof(double));
  double vol_i = 0.0, vol_ia = 0.0;
  double vol_m = 0.0, vol_ma = 0.0;

  if(params.flag_first_homogenization){

    for(int e = 0 ; e < nelm ; e++){

      if(elem_type[e] == 1){
	vol_ia += vol_elem;
	ne_i++;
      }
      else{
	vol_ma += vol_elem;
	ne_m++;
      }

    }
    ierr = MPI_Allreduce( &vol_ia, &vol_i, 1, MPI_DOUBLE, MPI_SUM, MICRO_COMM); if(ierr) return 1;
    ierr = MPI_Allreduce( &vol_ma, &vol_m, 1, MPI_DOUBLE, MPI_SUM, MICRO_COMM); if(ierr) return 1;
    vi = vol_i / vol_tot;
    vm = vol_m / vol_tot;
    myio_printf(&MICRO_COMM, "vi = %lf \n", vi );
    myio_printf(&MICRO_COMM, "vm = %lf \n", vm );
  }
  get_c_tan("FIBER" , -1, -1, NULL, c_i);
  get_c_tan("MATRIX", -1, -1, NULL, c_m);
 
  if(params.homog_method == HOMOG_METHOD_TAYLOR_PARALLEL){

    for(int i = 0 ; i < nvoi ; i++){
      for(int j = 0 ; j < nvoi ; j++)
	c[i*nvoi + j] = vi * c_i[i*nvoi + j] + vm * c_m[i*nvoi + j];
    }

  }
  else if(params.homog_method == HOMOG_METHOD_TAYLOR_SERIAL){

    int s;
    double *c_mi = malloc(nvoi*nvoi*sizeof(double));
    double *c_ii = malloc(nvoi*nvoi*sizeof(double));
    gsl_matrix_view  gsl_c_m , gsl_c_i ;
    gsl_matrix_view  gsl_c_mi, gsl_c_ii;

    gsl_permutation  *p;

    p  = gsl_permutation_alloc(nvoi);

    gsl_c_i  = gsl_matrix_view_array( c_i , nvoi, nvoi );
    gsl_c_m  = gsl_matrix_view_array( c_m , nvoi, nvoi );
    gsl_c_ii = gsl_matrix_view_array( c_ii, nvoi, nvoi );
    gsl_c_mi = gsl_matrix_view_array( c_mi, nvoi, nvoi );

    gsl_linalg_LU_decomp (&gsl_c_m.matrix, p, &s);
    gsl_linalg_LU_invert (&gsl_c_m.matrix, p, &gsl_c_mi.matrix);
    gsl_linalg_LU_decomp (&gsl_c_i.matrix, p, &s);
    gsl_linalg_LU_invert (&gsl_c_i.matrix, p, &gsl_c_ii.matrix);

    for(int i = 0; i < nvoi*nvoi ; i++)
      c_i[i] = vi * c_ii[i] + vm * c_mi[i];

    gsl_linalg_LU_decomp (&gsl_c_i.matrix, p, &s);
    gsl_linalg_LU_invert (&gsl_c_i.matrix, p, &gsl_c_mi.matrix);

    for(int i = 0; i < nvoi; i++){
      for(int j = 0; j < nvoi; j++)
	c[i*nvoi + j] = gsl_matrix_get( &gsl_c_mi.matrix, i, j );
    }

    gsl_permutation_free (p);
    free(c_mi);
    free(c_ii);
  }

  for(int i = 0; i < nvoi; i++){
    strain_ave[i] = strain_mac[i];
    stress_ave[i] = 0.0;
    for(int j = 0; j < nvoi; j++)
      stress_ave[i] += c[i*nvoi + j] * strain_mac[j];
  }

  return 0;
}


int mic_homog_us(double *strain_mac, double *strain_ave, double *stress_ave){

  if(params.flag_have_allocated == false){

    params.flag_have_allocated = true;

    int nnz = (dim==2)? 18:81;

    MatCreate(MICRO_COMM,&A);
    MatSetSizes(A, nl*dim, nl*dim, nn*dim, nn*dim);
    MatSetFromOptions(A);
    MatSeqAIJSetPreallocation(A, nnz, NULL);
    MatMPIAIJSetPreallocation(A, nnz, NULL, nnz, NULL);
    MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    MatGetOwnershipRange(A,&istart,&iend);
    nstart = istart / dim;
    nend   = iend   / dim;

    int *ghost_index;
    if(nproc_mic == 1)
      ngho = 0;
    else{
      if(rank_mic == 0)
	ngho = 0;
      else
	ngho = ( dim == 2 ) ? nx : nx*nz;
    }

    ghost_index = malloc( ngho*dim *sizeof(int) );

    int i;
    for( i = 0 ; i < ngho*dim  ; i++ )
      ghost_index[i] = istart - (( dim == 2 )? nx : nx*nz)*dim + i;

    VecCreateGhost(MICRO_COMM, nl*dim, nn*dim, ngho*dim, ghost_index, &x);
    VecDuplicate(x,&dx);
    VecDuplicate(x,&b);

    free(ghost_index);

    /* alloc arrays for boundary condition setting */
    if(nproc_mic == 1){
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
    if(rank_mic == 0){

      /* y = 0 */
      for( n = 0 ; n < nx ; n++ ){
	for( d = 0 ; d < dim ; d++ )
	  dir_ix_loc[c*dim + d] = n*dim + d;
	coor_dir[c*dim + 0] = n*hx;
	coor_dir[c*dim + 1] = 0.0;
	c++;
      }
    }
    if( rank_mic == (nproc_mic - 1) ){

      /* y = ly */
      for( n = 0 ; n < nx ; n++ ){
	for( d = 0 ; d < dim ; d++ )
	  dir_ix_loc[c*dim + d] = ((nyl-1)*nx + n)*dim + d;
	coor_dir[c*dim + 0] = n*hx;
	coor_dir[c*dim + 1] = ly;
	c++;
      }
    }
    if( nproc_mic > 1)
    {
      if( rank_mic == 0 )
      {
	/* x = 0 */
	for( n = 0 ; n < (nyl - 1) ; n++ ){
	  for( d = 0 ; d < dim ; d++ )
	    dir_ix_loc[c*dim + d] = (n+1)*nx*dim + d;
	  coor_dir[c*dim + 0] = 0;
	  coor_dir[c*dim + 1] = (ny_inf + n+1)*hy;
	  c++;
	}

	/* x = lx */
	for( n = 0 ; n < (nyl - 1) ; n++ ){
	  for( d = 0 ; d < dim ; d++ )
	    dir_ix_loc[c*dim + d] = (2*nx-1)*dim + n*nx*dim + d;
	  coor_dir[c*dim + 0] = lx;
	  coor_dir[c*dim + 1] = (ny_inf + n + 1)*hy;
	  c++;
	}
      }
      else if( rank_mic == (nproc_mic - 1) )
      {
	/* x = 0 */
	for( n = 0 ; n < (nyl - 1) ; n++ ){
	  for( d = 0 ; d < dim ; d++ )
	    dir_ix_loc[c*dim + d] = n*nx*dim + d;
	  coor_dir[c*dim + 0] = 0;
	  coor_dir[c*dim + 1] = ( ny_inf + n )*hy;
	  c++;
	}

	/* x = lx */
	for( n = 0 ; n < (nyl - 1) ; n++ ){
	  for( d = 0 ; d < dim ; d++ )
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
	for( n = 0 ; n < nyl ; n++ ){
	  for( d = 0 ; d < dim ; d++ )
	    dir_ix_loc[c*dim + d] = n*nx*dim + d;
	  coor_dir[c*dim + 0] = 0;
	  coor_dir[c*dim + 1] = ( ny_inf + n )*hy;
	  c++;
	}

	/* x = lx */
	for( n = 0 ; n < nyl ; n++ ){
	  for( d = 0 ; d < dim ; d++ )
	    dir_ix_loc[c*dim + d] = (nx-1)*dim + n*nx*dim + d;
	  coor_dir[c*dim + 0] = lx;
	  coor_dir[c*dim + 1] = ( ny_inf + n )*hy;
	  c++;
	}
      }
    }
    else{

      /* x = 0 */
      for( n = 0 ; n < nyl-2; n++ ){
	for( d = 0 ; d < dim ; d++ )
	  dir_ix_loc[c*dim + d] = (n+1)*nx*dim + d;
	coor_dir[c*dim + 0] = 0;
	coor_dir[c*dim + 1] = (ny_inf + n+1)*hy;
	c++;
      }

      /* x = lx */
      for( n = 0 ; n < nyl-2; n++ ){
	for( d = 0 ; d < dim ; d++ )
	  dir_ix_loc[c*dim + d] = ((n+2)*nx-1)*dim + d;
	coor_dir[c*dim + 0] = lx;
	coor_dir[c*dim + 1] = (ny_inf + n + 1)*hy;
	c++;
      }

    }

    for( i = 0; i < ndir_ix ; i++ )
      dir_ix_glo[i] = local_to_global_index( dir_ix_loc[i] );

  } // first time for allocation

  /* Set the displacement boundary conditions "u = E . X" */

  VecZeroEntries(x);
  VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

  Vec     x_loc;
  double  *x_arr;
  VecGhostGetLocalForm(x, &x_loc);
  VecGetArray(x_loc, &x_arr);

  if(dim == 2){

    double  displ[2]; // (ux,uy) displacement

    for(int n = 0 ; n < ndir_ix/dim ; n++ ){
      strain_x_coord(strain_mac, &coor_dir[n*dim], displ);
      for(int d = 0 ; d < dim ; d++)
	x_arr[dir_ix_loc[n*dim + d]] = displ[d];
    }
  }

  VecRestoreArray(x_loc, &x_arr);
  VecGhostRestoreLocalForm(x, &x_loc);
  VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

  int     i;
  int     nr_its = 0; 
  double  *b_arr;
  double  norm = params.non_linear_min_norm_tol*10;

  VecNorm( x , NORM_2 , &norm );
  PRINTF2("|x| = %lf\n", norm);

  while(nr_its < params.non_linear_max_its && norm > params.non_linear_min_norm_tol){

    save_event(MICRO_COMM, "ass_0");

    assembly_b();
    
    VecGetArray(b, &b_arr);
    for( i = 0; i < ndir_ix ; i++ )
      b_arr[dir_ix_loc[i]] = 0.0;
    VecRestoreArray(b, &b_arr);
    VecGhostUpdateBegin(b, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(b, INSERT_VALUES, SCATTER_FORWARD);

    VecNorm(b, NORM_2, &norm);
    PRINTF2("|b| = %lf \n", norm);

    if(norm < params.non_linear_min_norm_tol) break;

    VecScale( b, -1.0 );

    assembly_A();

    MatZeroRowsColumns(A, ndir_ix, dir_ix_glo, 1.0, NULL, NULL);
    save_event(MICRO_COMM, "ass_1");

    save_event(MICRO_COMM, "sol_0");
    KSPSetOperators(ksp, A, A);
    KSPSolve(ksp, b, dx);
    save_event(MICRO_COMM, "sol_1");

    VecAXPY(x, 1.0, dx);
    VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

    nr_its ++;
  }
  save_event(MICRO_COMM, "ass_1");

  get_averages(strain_ave, stress_ave);

  if( flag_print & (1<<PRINT_PETSC) ){
    PetscViewer  viewer;
    PetscViewerASCIIOpen(MICRO_COMM,"A.dat" ,&viewer); MatView(A ,viewer);
    PetscViewerASCIIOpen(MICRO_COMM,"b.dat" ,&viewer); VecView(b ,viewer);
    PetscViewerASCIIOpen(MICRO_COMM,"x.dat" ,&viewer); VecView(x ,viewer);
    PetscViewerASCIIOpen(MICRO_COMM,"dx.dat",&viewer); VecView(dx,viewer);
  }

  return 0;
}


int homogenize_init(void){

  params.have_linear_materials = (material_are_all_linear(&material_list) == true) ? true : false;

  int ierr = 0;

  if(params.have_linear_materials == true){

    ierr = homogenize_calculate_c_tangent_around_zero(params.c_tangent_linear);

    params.c_tangent_linear_calculated = true;

  }

  return ierr;
}


int homogenize_get_strain_stress(double *strain_mac, double *strain_ave, double *stress_ave){

  if(params.have_linear_materials == true && params.c_tangent_linear_calculated == true){

    for(int i = 0 ; i < nvoi ; i++){
      strain_ave[i] = strain_mac[i];
      stress_ave[i] = 0.0;
      for(int j = 0 ; j < nvoi ; j++)
	stress_ave[i] += params.c_tangent_linear[i*nvoi + j] * strain_mac[j];
    }

  }
  else{

    int ierr = homogenize_get_strain_stress_non_linear(strain_mac, strain_ave, stress_ave);
    if(ierr) return 1;

  }
  return 0;
}


int homogenize_get_strain_stress_non_linear(double *strain_mac, double *strain_ave, double *stress_ave){

  int ierr = 0;

  if(params.homog_method == HOMOG_METHOD_TAYLOR_PARALLEL || params.homog_method == HOMOG_METHOD_TAYLOR_SERIAL){

    ierr = mic_homogenize_taylor(strain_mac, strain_ave, stress_ave);
  }
  else if(params.homog_method == HOMOG_METHOD_UNIF_STRAINS){

    ierr = mic_homog_us(strain_mac, strain_ave, stress_ave);
  }

  return ierr;
}


int homogenize_get_c_tangent(double *strain_mac, double **c_tangent){

  int ierr = 0;

  if(params.have_linear_materials == true && params.have_linear_materials == true){

    (*c_tangent) = params.c_tangent_linear;

  }
  else if(params.have_linear_materials == false){

    ierr = homogenize_calculate_c_tangent(strain_mac, params.c_tangent);
    (*c_tangent) = params.c_tangent;

  }

  return ierr;
}


int homogenize_calculate_c_tangent_around_zero(double *c_tangent){

  double strain[MAX_NVOIGT];
  ARRAY_SET_TO_ZERO(strain, nvoi);

  int ierr = homogenize_calculate_c_tangent(strain, c_tangent);

  return ierr;
}


int homogenize_calculate_c_tangent(double *strain_mac, double *c_tangent){

  int ierr = 0;
  double strain_1[MAX_NVOIGT], strain_2[MAX_NVOIGT];
  double stress_1[MAX_NVOIGT], stress_2[MAX_NVOIGT];
  double strain_aux[MAX_NVOIGT];

  ARRAY_COPY(strain_1, strain_mac, nvoi);

  ierr = homogenize_get_strain_stress(strain_1, strain_aux, stress_1);

  for(int i = 0 ; i < nvoi ; i++){

    ARRAY_COPY(strain_2, strain_mac, nvoi);

    strain_2[i] = strain_2[i] + HOMOGENIZE_DELTA_STRAIN;

    ierr = homogenize_get_strain_stress(strain_2, strain_aux, stress_2);

    for(int j = 0 ; j < nvoi ; j++)
      c_tangent[j*nvoi + i] = (stress_2[j] - stress_1[j]) / (strain_2[i] - strain_1[i]);

  }

  return ierr;
}


int strain_x_coord( double *strain , double *coord , double *u ){

  if(dim == 2){
    u[0] = strain[0]   * coord[0] + strain[2]/2 * coord[1] ;
    u[1] = strain[2]/2 * coord[0] + strain[1]   * coord[1] ;
  }

  return 0;
}


int assembly_b(void){

  VecZeroEntries(b);
  VecGhostUpdateBegin(b ,INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(b ,INSERT_VALUES, SCATTER_FORWARD);

  Vec b_loc;
  double *b_arr;

  VecGhostGetLocalForm(b, &b_loc);
  VecGetArray(b_loc, &b_arr);

  double *res_elem  = malloc(dim*npe*sizeof(double));

  for(int e = 0 ; e < nelm ; e++){

    ARRAY_SET_TO_ZERO(res_elem, npe*dim)

    get_local_elem_index(e, loc_elem_index);

    for(int gp = 0; gp < ngp ; gp++){

      get_strain(e, gp, strain_gp);
      get_stress(e, gp, strain_gp, stress_gp);

      for(int i = 0 ; i < npe*dim ; i++){
	for(int j = 0; j < nvoi ; j++)
	  res_elem[i] += struct_bmat[j][i][gp] * stress_gp[j] * struct_wp[gp];
      }

    }

    for(int i = 0 ; i < npe*dim ; i++ )
      b_arr[loc_elem_index[i]] += res_elem[i];

  }

  VecRestoreArray(b_loc, &b_arr);
  VecGhostRestoreLocalForm(b, &b_loc);

  VecGhostUpdateBegin(b, ADD_VALUES, SCATTER_REVERSE);
  VecGhostUpdateEnd(b, ADD_VALUES, SCATTER_REVERSE);
  VecGhostUpdateBegin(b, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(b, INSERT_VALUES, SCATTER_FORWARD);

  return 0;
}


int assembly_A(void){

  MatZeroEntries(A);

  int    e, gp;
  double *k_elem    = malloc( dim*npe*dim*npe * sizeof(double));
  double *c         = malloc( nvoi*nvoi       * sizeof(double));
  int    i, j, k, h;

  for( e = 0 ; e < nelm ; e++ ){

    /* set to 0 res_elem */
    for( i = 0 ; i < npe*dim*npe*dim ; i++ )
      k_elem[i] = 0.0;

    get_global_elem_index(e, glo_elem_index);

    for( gp = 0; gp < ngp ; gp++ ){

      get_strain( e , gp, strain_gp );
      get_c_tan( NULL , e , gp , strain_gp , c );

      for( i = 0 ; i < npe*dim ; i++ ){
	for( j = 0 ; j < npe*dim ; j++ ){
	  for( k = 0; k < nvoi ; k++ ){
	    for( h = 0; h < nvoi ; h++ )
	      k_elem[ i*npe*dim + j] += \
	      struct_bmat[h][i][gp] * c[ h*nvoi + k ] * struct_bmat[k][j][gp] * struct_wp[gp];
	  }
	}
      }

    }
    MatSetValues( A, npe*dim, glo_elem_index, npe*dim, glo_elem_index, k_elem, ADD_VALUES );

  }

  MatAssemblyBegin( A , MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd(   A , MAT_FINAL_ASSEMBLY );

  return 0;
}


int get_averages(double *strain_ave, double *stress_ave){

  int    i, e, gp;
  int    ierr;

  double * strain_part = malloc( nvoi * sizeof(double) );
  double * stress_part = malloc( nvoi * sizeof(double) );

  for( i = 0 ; i < nvoi ; i++ )
    strain_part[i] = stress_part[i] = 0.0;

  for( e = 0 ; e < nelm ; e++ )
  {
    for( gp = 0; gp < ngp ; gp++ )
    {

      /* calc strain gp */
      get_strain( e , gp, strain_gp );

      /* we get stress = f(strain) */
      get_stress( e , gp , strain_gp , stress_gp );

      for( i = 0; i < nvoi ; i++ )
      {
	stress_part[i] += stress_gp[i] * struct_wp[gp];
	strain_part[i] += strain_gp[i] * struct_wp[gp];
      }
    }
  }

  ierr = MPI_Allreduce( stress_part, stress_ave, nvoi, MPI_DOUBLE, MPI_SUM, MICRO_COMM );
  ierr = MPI_Allreduce( strain_part, strain_ave, nvoi, MPI_DOUBLE, MPI_SUM, MICRO_COMM );

  for( i = 0; i < nvoi ; i++ )
  {
    stress_ave[i] /= vol_tot;
    strain_ave[i] /= vol_tot;
  }
  return ierr;
}


int get_elem_properties(void){

  for(int e = 0 ; e < nelm ; e++){

    double strain_aux[MAX_NVOIGT];
    double stress_aux[MAX_NVOIGT];
    for(int v = 0 ; v < nvoi ; v++)
      strain_aux[v] = stress_aux[v] = 0.0;

    for (int gp = 0 ; gp < ngp ; gp++){

      get_strain( e , gp, strain_gp );
      get_stress( e , gp, strain_gp, stress_gp );
      for(int v = 0 ; v < nvoi ; v++ ){
	strain_aux[v] += strain_gp[v] * struct_wp[gp];
	stress_aux[v] += stress_gp[v] * struct_wp[gp];
      }

    }
    for (int v = 0 ; v < nvoi ; v++){
      elem_strain[ e*nvoi + v ] = strain_aux[v] / vol_elem;
      elem_stress[ e*nvoi + v ] = stress_aux[v] / vol_elem;
    }

  }

  return 0;
}


int get_strain(int e, int gp, double *strain_gp){

  Vec     x_loc;
  double  *x_arr;

  VecGhostGetLocalForm( x    , &x_loc );
  VecGetArray         ( x_loc, &x_arr );

  int    i , v;

  get_local_elem_index(e, loc_elem_index);

  for( i = 0 ; i < npe*dim ; i++ )
    elem_disp[i] = x_arr[loc_elem_index[i]];

  for( v = 0; v < nvoi ; v++ ){
    strain_gp[v] = 0.0;
    for( i = 0 ; i < npe*dim ; i++ )
      strain_gp[v] += struct_bmat[v][i][gp] * elem_disp[i];
  }

  return 0;
}


int get_stress(int e, int gp, double *strain_gp, double *stress_gp){

  char *word_to_search;

  switch(elem_type[e]){

    case ID_FIBER:
      word_to_search = strdup("FIBER");
      break;

    case ID_MATRIX:
      word_to_search = strdup("MATRIX");
      break;

    default:
      return 1;
  }

  material_t  *mat_p;
  node_list_t *pm = material_list.head;

  while(pm != NULL){
    mat_p = (material_t *)pm->data;
    if(strcmp(mat_p->name, word_to_search) == 0) break;
    pm = pm->next;
  }

  if(pm == NULL) return 1;

  int ierr = material_get_stress(mat_p, dim, strain_gp, stress_gp);

  return ierr;
}


int get_c_tan( const char *name, int e , int gp, double *strain_gp , double *c_tan )
{

  char *word_to_search;

  if(name != NULL){
    word_to_search = strdup(name);
  }
  else if( elem_type[e] == ID_FIBER )
  {
    word_to_search = strdup("FIBER");
  }
  else if( elem_type[e] == ID_MATRIX )
  {
    word_to_search = strdup("MATRIX");
  }

  material_t  *mat_p;
  node_list_t *pm = material_list.head;

  while( pm )
  {
    mat_p = (material_t *)pm->data;
    if( strcmp ( mat_p->name , word_to_search ) == 0 ) break;
    pm = pm->next;
  }
  if(!pm) return 1;

  material_get_c_tang( mat_p, dim, strain_gp, c_tan );

  return 0;
}


int get_elem_centroid(int e, int dim, double *centroid){

  if(dim == 2){
    if( rank_mic == 0 ){
      centroid[0] = ( e%nex + 0.5 )*hx;
      centroid[1] = ( e/nex + 0.5 )*hy;
    }
    else{
      centroid[0] = ( e%nex + 0.5              )*hx;
      centroid[1] = ( e/nex + 0.5 + ny_inf - 1 )*hy;
    }
  }

  return 0;
}


int get_local_elem_index(int e, int *loc_elem_index){

  if(dim == 2){
    if(  rank_mic == 0 ){
      for(int d = 0 ; d < dim ; d++){
	int n0 = (e%nex) + (e/nex)*nx;
	loc_elem_index[ 0*dim + d ] = ( n0          ) * dim + d ;
	loc_elem_index[ 1*dim + d ] = ( n0 + 1      ) * dim + d ;
	loc_elem_index[ 2*dim + d ] = ( n0 + nx + 1 ) * dim + d ;
	loc_elem_index[ 3*dim + d ] = ( n0 + nx + 0 ) * dim + d ;
      }
    }
    else if(e >= nex ){
      for(int d = 0 ; d < dim ; d++){
	int n0 = (e%nex) + (e/nex-1)*nx;
	loc_elem_index[ 0*dim + d ] = ( n0          ) * dim + d ;
	loc_elem_index[ 1*dim + d ] = ( n0 + 1      ) * dim + d ;
	loc_elem_index[ 2*dim + d ] = ( n0 + nx + 1 ) * dim + d ;
	loc_elem_index[ 3*dim + d ] = ( n0 + nx     ) * dim + d ;
      }

    }
    else{
      for(int d = 0 ; d < dim ; d++){
	int n0 = (e%nex) + (e/nex)*nx;
	loc_elem_index[ 0*dim + d ] = ( n0 + nl     ) * dim + d ; // is a ghost
	loc_elem_index[ 1*dim + d ] = ( n0 + nl + 1 ) * dim + d ; // is a ghost
	loc_elem_index[ 2*dim + d ] = ( n0 + 1      ) * dim + d ;
	loc_elem_index[ 3*dim + d ] = ( n0          ) * dim + d ;
      }
    }
  }
  return 0;
}


int get_local_elem_node(int e , int *n_loc){

  if(dim == 2){
    if(rank_mic == 0){
      int n0 = (e%nex) + (e/nex)*nx;
      n_loc[0] = n0          ;
      n_loc[1] = n0 + 1      ;
      n_loc[2] = n0 + nx + 1 ;
      n_loc[3] = n0 + nx     ;
    }
    else if(e >= nex ){
      int n0 = (e%nex) + (e/nex-1)*nx;
      n_loc[0] = n0          ;
      n_loc[1] = n0 + 1      ;
      n_loc[2] = n0 + nx + 1 ;
      n_loc[3] = n0 + nx     ;
    }
    else{
      int n0 = (e%nex) + (e/nex)*nx;
      n_loc[0] = n0 + nl     ; // is a ghost
      n_loc[1] = n0 + nl + 1 ; // is a ghost
      n_loc[2] = n0 + 1      ;
      n_loc[3] = n0          ;
    }
  }
  return 0;
}


int get_global_elem_index(int e, int *glo_elem_index){
  
  if(dim == 2){
    if(rank_mic == 0){
      for(int d = 0 ; d < dim ; d++ ){
	int n0 = (e%nex) + (e/nex)*nx;
	glo_elem_index[ 0*dim + d ] = ( n0          ) * dim + d ;
	glo_elem_index[ 1*dim + d ] = ( n0 + 1      ) * dim + d ;
	glo_elem_index[ 2*dim + d ] = ( n0 + nx + 1 ) * dim + d ;
	glo_elem_index[ 3*dim + d ] = ( n0 + nx     ) * dim + d ;
      }
    }
    else{
      for(int d = 0 ; d < dim ; d++){
	int n0 = (e%nex) + (e/nex)*nx + (ny_inf-1)*nx;
	glo_elem_index[ 0*dim + d ] = ( n0          ) * dim + d ;
	glo_elem_index[ 1*dim + d ] = ( n0 + 1      ) * dim + d ;
	glo_elem_index[ 2*dim + d ] = ( n0 + nx + 1 ) * dim + d ;
	glo_elem_index[ 3*dim + d ] = ( n0 + nx     ) * dim + d ;
      }
    }
  }

  return 0;
}


int local_to_global_index( int local )
{

  if ( rank_mic == 0 ) return local;

  if( dim == 2 ){

    return ny_inf * nx * dim + local;

  }

  return 0;
}


int get_node_local_coor( int n , double * coord )
{
  if(dim == 2){
    if( rank_mic == 0 ){
      coord[0] = (n % nx)*hx;
      coord[1] = (n / nx)*hy;
    }
    else{
      coord[0] = (n % nx)*hx;
      coord[1] = (n / nx + ny_inf)*hy;
    }
  }
  return 0;
}


int get_node_ghost_coor( int n , double * coord )
{
  if( dim == 2 ){
    coord[0] = (n % nx)*hx;
    coord[1] = (ny_inf + n / nx - 1)*hy;
  }
  return 0;
}


int init_shapes( double ***sh, double ****dsh, double **wp )
{
  int    nsh = ( dim == 2 ) ? 4 : 8;
  int    gp;
  double *xp = malloc( ngp*dim * sizeof(double));

  if( dim == 2 )
  {
    xp[0] = -0.577350269189626;   xp[1]= -0.577350269189626;
    xp[2] = +0.577350269189626;   xp[3]= -0.577350269189626;
    xp[4] = +0.577350269189626;   xp[5]= +0.577350269189626;
    xp[6] = -0.577350269189626;   xp[7]= +0.577350269189626;
  }

  int is, d;

  *sh = malloc( nsh * sizeof(double*));
  for( is = 0 ; is < nsh ; is++ ){
    (*sh)[is] = malloc( ngp * sizeof(double));
  }
  
  *dsh  = malloc( nsh * sizeof(double**));
  for( is = 0 ; is < nsh ; is++ ){
    (*dsh)[is] = malloc( dim * sizeof(double*));
    for( d = 0 ; d < dim ; d++ ){
      (*dsh)[is][d] = malloc( ngp * sizeof(double));
    }
  }

  *wp   = malloc( ngp * sizeof(double));

  if( dim == 2 )
  {
    
    for( gp = 0 ; gp < ngp ; gp++ ){
      (*sh)[0][gp] = (1 - xp[2*gp]) * (1 - xp[2*gp+1])/4;
      (*sh)[1][gp] = (1 + xp[2*gp]) * (1 - xp[2*gp+1])/4;
      (*sh)[2][gp] = (1 + xp[2*gp]) * (1 + xp[2*gp+1])/4;
      (*sh)[3][gp] = (1 - xp[2*gp]) * (1 + xp[2*gp+1])/4;
    }

    for( gp = 0 ; gp < ngp ; gp++ ){
      (*dsh)[0][0][gp] = -1 * (1 - xp[2*gp+1]) /4 * 2/hx; // d phi / d x
      (*dsh)[1][0][gp] = +1 * (1 - xp[2*gp+1]) /4 * 2/hx;
      (*dsh)[2][0][gp] = +1 * (1 + xp[2*gp+1]) /4 * 2/hx;
      (*dsh)[3][0][gp] = -1 * (1 + xp[2*gp+1]) /4 * 2/hx;
      (*dsh)[0][1][gp] = -1 * (1 - xp[2*gp+0]) /4 * 2/hy; // d phi / d y
      (*dsh)[1][1][gp] = -1 * (1 + xp[2*gp+0]) /4 * 2/hy;
      (*dsh)[2][1][gp] = +1 * (1 + xp[2*gp+0]) /4 * 2/hy;
      (*dsh)[3][1][gp] = +1 * (1 - xp[2*gp+0]) /4 * 2/hy;
    }

    for( gp = 0 ; gp < ngp ; gp++ )
      (*wp)[gp] = vol_elem / ngp;

  }

  free(xp);

  return 0;
}
