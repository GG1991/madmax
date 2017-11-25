/*
   Routines for performing homogenization on RVE

   Author > Guido Giuntoli
   Date   > 18-08-2017
 */

#include "micro.h"

int mic_homogenize_taylor( MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6] )
{
  /* mixture theory */
    
  int          i, j, e, ierr;
  int          ne_i = 0, ne_m = 0;
  double       *c_i = malloc( nvoi*nvoi * sizeof(double));
  double       *c_m = malloc( nvoi*nvoi * sizeof(double));
  double       *c   = malloc( nvoi*nvoi * sizeof(double));
  double       vol_i = 0.0, vol_ia = 0.0;  // inclusion volume
  double       vol_m = 0.0, vol_ma = 0.0;  // matrix volume

  if(first_time_homo){

    for( e = 0 ; e < nelm ; e++ ){

      if( elem_type[e] == 1 ){
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
    vi   = vol_i / vol_tot;              // inclusion fraction
    vm   = vol_m / vol_tot;              // matrix fraction
    printf_p( &MICRO_COMM, "vi = %lf \n", vi );
    printf_p( &MICRO_COMM, "vm = %lf \n", vm );
  }
  get_c_tan("FIBER" , 0, 0, NULL, c_i);  //returns c_i of FIBER
  get_c_tan("MATRIX", 0, 0, NULL, c_m);  //returns c_m of MATRIX
 
  if( homo_type == TAYLOR_P ){

    /* PARALLEL THEORY */
    for( i = 0 ; i < nvoi ; i++ ){
      for( j = 0 ; j < nvoi ; j++ )
	c[i*nvoi + j] = vi * c_i[i*nvoi + j] + vm * c_m[i*nvoi + j];
    }

  }
  else if( homo_type == TAYLOR_S ){

    /* SERIAL THEORY */
    int              s;
    double *c_mi = malloc( nvoi*nvoi * sizeof(double));
    double *c_ii = malloc( nvoi*nvoi * sizeof(double));
    gsl_matrix_view  gsl_c_m , gsl_c_i ;  // gsl arrays of c_i and c_m
    gsl_matrix_view  gsl_c_mi, gsl_c_ii;  // gsl arrays of c_i and c_m (inverted)

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

    /* we reuse c_i */
    for( i = 0; i < nvoi*nvoi ; i++)
      c_i[i] = vi * c_ii[i] + vm * c_mi[i];

    gsl_linalg_LU_decomp (&gsl_c_i.matrix, p, &s);
    gsl_linalg_LU_invert (&gsl_c_i.matrix, p, &gsl_c_mi.matrix);

    for( i = 0; i < nvoi; i++ ){
      for( j = 0; j < nvoi; j++ )
	c[i*nvoi + j] = gsl_matrix_get( &gsl_c_mi.matrix, i, j );
    }

    gsl_permutation_free (p);
    free(c_mi);
    free(c_ii);
  }

  for( i = 0; i < nvoi; i++ ){
    strain_ave[i] = strain_mac[i];
    stress_ave[i] = 0.0;
    for( j = 0; j < nvoi; j++ )
      stress_ave[i] += c[i*nvoi + j] * strain_mac[j];
  }

  return 0;
}

/****************************************************************************************************/

int mic_homog_us(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  /* 
     homogenization routine with uniform strain BC
     for structured mesh 

     a) alloc matrix/vector if it is first time
     b) assembly residue -> evaluate norm
     c) assembly jacobian
     d) solve system 
     e) iterate to b) with NR method
     f) get average properties

     partimos el dominio en tandas por el eje y

   */

  if( flag_first_alloc == true ){

    flag_first_alloc = false;

    int nnz = (dim==2)? 18:81;                        // nonzeros per row

    MatCreate(MICRO_COMM,&A);
    MatSetSizes(A, nl*dim, nl*dim, nn*dim, nn*dim);
    MatSetFromOptions(A);
    MatSeqAIJSetPreallocation(A, nnz, NULL);
    MatMPIAIJSetPreallocation(A, nnz, NULL, nnz, NULL);
    MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    MatGetOwnershipRange(A,&istart,&iend);
    nstart = istart / dim;
    nend   = iend   / dim;
    //    ny_inf = ( dim == 2 ) ? nstart / nx : nstart / (nx*nz);

    int *ghost_index;
    if( nproc_mic == 1 )
      ngho = 0;
    else{
      if( rank_mic == 0 )
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

  VecZeroEntries( x );
  VecGhostUpdateBegin( x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd  ( x, INSERT_VALUES, SCATTER_FORWARD);

  Vec     x_loc;
  double  *x_arr;
  VecGhostGetLocalForm( x    , &x_loc );
  VecGetArray         ( x_loc, &x_arr );

  int  n , d;
  if( dim == 2 ){

    double  displ[2]; // (ux,uy) displacement

    for( n = 0 ; n < ndir_ix/dim ; n++ )
    {
      /* calc displ on the node */
      strain_x_coord( strain_mac , &coor_dir[n*dim] , displ );
      for( d = 0 ; d < dim ; d++ )
	x_arr[dir_ix_loc[n*dim + d]] = displ[d];
    }
  }

  VecRestoreArray         ( x_loc , &x_arr );
  VecGhostRestoreLocalForm( x     , &x_loc );
  VecGhostUpdateBegin( x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd  ( x, INSERT_VALUES, SCATTER_FORWARD);

  int     i;
  int     nr_its = 0; 
  double  *b_arr;
  double  norm = nr_norm_tol*10;

  VecNorm( x , NORM_2 , &norm );
  if(!flag_coupling)
    printf_p( &MICRO_COMM,"|x| = %lf \n",norm);

  while( nr_its < nr_max_its && norm > nr_norm_tol )
  {
    save_event( MICRO_COMM, "ass_0" );

    /* assembly "b" (residue) using "x" (displacement) */
    assembly_b();
    
    VecGetArray( b, &b_arr );
    for( i = 0; i < ndir_ix ; i++ )
      b_arr[dir_ix_loc[i]] = 0.0;
    VecRestoreArray( b, &b_arr );
    VecGhostUpdateBegin( b, INSERT_VALUES, SCATTER_FORWARD );
    VecGhostUpdateEnd  ( b, INSERT_VALUES, SCATTER_FORWARD );

    VecNorm( b , NORM_2 , &norm );

    if(!flag_coupling)
      printf_p( &MICRO_COMM,"|b| = %lf \n",norm);

    if( !(norm > nr_norm_tol) ) break;

    VecScale( b, -1.0 );

    /* assembly "A" (jacobian) using "x" (displacement) */
    assembly_A();

    MatZeroRowsColumns( A, ndir_ix, dir_ix_glo, 1.0, NULL, NULL );
    save_event( MICRO_COMM, "ass_1" );

    save_event( MICRO_COMM, "sol_0" );
    KSPSetOperators( ksp, A, A );
    KSPSolve( ksp, b, dx );
    save_event( MICRO_COMM, "sol_1" );

    VecAXPY( x, 1.0, dx );
    VecGhostUpdateBegin( x, INSERT_VALUES, SCATTER_FORWARD );
    VecGhostUpdateEnd  ( x, INSERT_VALUES, SCATTER_FORWARD );

    nr_its ++;
  }
  save_event( MICRO_COMM, "ass_1" );

  /* get the integrals */
  get_averages( strain_ave, stress_ave );

  if( flag_print & (1<<PRINT_PETSC) ){
    PetscViewer  viewer;
    PetscViewerASCIIOpen(MICRO_COMM,"A.dat" ,&viewer); MatView(A ,viewer);
    PetscViewerASCIIOpen(MICRO_COMM,"b.dat" ,&viewer); VecView(b ,viewer);
    PetscViewerASCIIOpen(MICRO_COMM,"x.dat" ,&viewer); VecView(x ,viewer);
    PetscViewerASCIIOpen(MICRO_COMM,"dx.dat",&viewer); VecView(dx,viewer);
  }

  return 0;
}

/****************************************************************************************************/

int strain_x_coord( double * strain , double * coord , double * u )
{
  /* b = mat . a */

  if( dim == 2 ){
    u[0] = strain[0]   * coord[0] + strain[2]/2 * coord[1] ;
    u[1] = strain[2]/2 * coord[0] + strain[1]   * coord[1] ;
  }

  return 0;
}

/****************************************************************************************************/

int assembly_b(void)
{

  VecZeroEntries( b);
  VecGhostUpdateBegin( b , INSERT_VALUES , SCATTER_FORWARD );
  VecGhostUpdateEnd  ( b , INSERT_VALUES , SCATTER_FORWARD );

  Vec      b_loc;
  double  *b_arr;

  VecGhostGetLocalForm( b , &b_loc );
  VecGetArray         ( b_loc, &b_arr );

  int    e, gp;
  int    i, j;
  double *res_elem  = malloc( dim*npe * sizeof(double));

  for( e = 0 ; e < nelm ; e++ ){

    /* set to 0 res_elem */
    for( i = 0 ; i < npe*dim ; i++ )
      res_elem[i] = 0.0;

    /* get the local indeces of the element vertex nodes */
    get_local_elem_index(e, loc_elem_index);

    for( gp = 0; gp < ngp ; gp++ ){

      /* calc strain gp */
      get_strain( e , gp, strain_gp );

      /* we get stress = f(strain) */
      get_stress( e , gp , strain_gp , stress_gp );

      for( i = 0 ; i < npe*dim ; i++ ){
	for( j = 0; j < nvoi ; j++ )
	  res_elem[i] += struct_bmat[j][i][gp] * stress_gp[j] * struct_wp[gp];
      }

    }

    for( i = 0 ; i < npe*dim ; i++ )
      b_arr[loc_elem_index[i]] += res_elem[i];

  }

  VecRestoreArray         ( b_loc, &b_arr );
  VecGhostRestoreLocalForm( b    , &b_loc );

  /* from the local and ghost part with add to all processes */
  VecGhostUpdateBegin( b, ADD_VALUES   , SCATTER_REVERSE );
  VecGhostUpdateEnd  ( b, ADD_VALUES   , SCATTER_REVERSE );
  VecGhostUpdateBegin( b, INSERT_VALUES, SCATTER_FORWARD );
  VecGhostUpdateEnd  ( b, INSERT_VALUES, SCATTER_FORWARD );

  return 0;
}

/****************************************************************************************************/

int assembly_A( void )
{

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

      /* calc strain gp */
      get_strain( e , gp, strain_gp );

      /* we get stress = f(strain) */
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

  /* communication between processes */
  MatAssemblyBegin( A , MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd(   A , MAT_FINAL_ASSEMBLY );

  return 0;
}

/****************************************************************************************************/

int get_averages( double * strain_ave, double * stress_ave )
{
  /* Calculate averange strain and stress tensors on the hole domain
     the operation is an All_Reduce on all processes */

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

/****************************************************************************************************/

int get_elem_properties( void )
{

  /* fills *elem_strain, *elem_stress, *elem_type, *elem_energy */

  int      e, v, gp;
  double  *strain_aux = malloc( nvoi * sizeof(double) );
  double  *stress_aux = malloc( nvoi * sizeof(double) );

  for ( e = 0 ; e < nelm ; e++ ){

    for ( v = 0 ; v < nvoi ; v++ )
      strain_aux[v] = stress_aux[v] = 0.0;

    for ( gp = 0 ; gp < ngp ; gp++ ){

      get_strain( e , gp, strain_gp );
      get_stress( e , gp, strain_gp, stress_gp );
      for ( v = 0 ; v < nvoi ; v++ ){
	strain_aux[v] += strain_gp[v] * struct_wp[gp];
	stress_aux[v] += stress_gp[v] * struct_wp[gp];
      }

    }
    for ( v = 0 ; v < nvoi ; v++ ){
      elem_strain[ e*nvoi + v ] = strain_aux[v] / vol_elem;
      elem_stress[ e*nvoi + v ] = stress_aux[v] / vol_elem;
    }

  }

  return 0;
}

/****************************************************************************************************/

int get_strain( int e , int gp, double *strain_gp )
{

  Vec     x_loc;
  double  *x_arr;

  VecGhostGetLocalForm( x    , &x_loc );
  VecGetArray         ( x_loc, &x_arr );

  int    i , v;

  /* get the local indeces of the element vertex nodes */
  get_local_elem_index(e, loc_elem_index);

  /* get the elemental displacements */
  for( i = 0 ; i < npe*dim ; i++ )
    elem_disp[i] = x_arr[loc_elem_index[i]];

  /* calc strain gp */
  for( v = 0; v < nvoi ; v++ ){
    strain_gp[v] = 0.0;
    for( i = 0 ; i < npe*dim ; i++ )
      strain_gp[v] += struct_bmat[v][i][gp] * elem_disp[i];
  }

  return 0;
}

/****************************************************************************************************/

int get_stress( int e , int gp, double *strain_gp , double *stress_gp )
{

  /* returns the stress according to the elemet type */

  char *word_to_search;

  if( elem_type[e] == ID_FIBER )
  {
    word_to_search = strdup( "FIBER" );
  }
  else if( elem_type[e] == ID_MATRIX )
  {
    word_to_search = strdup( "MATRIX" );
  }

  material_t  *mat_p;
  node_list_t *pm = material_list.head;

  while( pm )
  {
    mat_p = (material_t *)pm->data;
    if( strcmp ( mat_p->name , word_to_search ) == 0 ) break;
    pm = pm->next;
  }
  if( !pm ) return 1;

  /*
     now that we now the material (mat_p) we calculate
     stress = f(strain)
   */
  mat_get_stress( mat_p, dim, strain_gp, stress_gp );

  return 0;
}

/****************************************************************************************************/

int get_c_tan( const char *name, int e , int gp, double *strain_gp , double *c_tan )
{

  /*
     returns the c_tan (tangent constitutive tensor)
     according to the elemet type
   */

  char *word_to_search;

  if( elem_type[e] == ID_FIBER )
  {
    word_to_search = strdup( "FIBER" );
  }
  else if( elem_type[e] == ID_MATRIX )
  {
    word_to_search = strdup( "MATRIX" );
  }

  material_t  *mat_p;
  node_list_t *pm = material_list.head;

  while( pm )
  {
    mat_p = (material_t *)pm->data;
    if( strcmp ( mat_p->name , word_to_search ) == 0 ) break;
    pm = pm->next;
  }
  if( !pm ) return 1;

  /*
     now that we now the material (mat_p) we calculate
     stress = f(strain)
   */
  mat_get_c_tang( mat_p, dim, strain_gp, c_tan );

  return 0;
}

/****************************************************************************************************/

int get_elem_centroid( int e, int dim, double *centroid )
{

  /* formula only valid for sequencial now */

  if( dim == 2 ){
    if( rank_mic == 0 ){
      centroid[0] = ( e % nex + 0.5 ) * hx;
      centroid[1] = ( e / nex + 0.5 ) * hy;
    }
    else{
      centroid[0] = ( e % nex + 0.5              ) * hx;
      centroid[1] = ( e / nex + 0.5 + ny_inf - 1 ) * hy;
    }
  }

  return 0;
}

/****************************************************************************************************/

int get_local_elem_index( int e, int *loc_elem_index )
{

  /* 
     Returns the local position in the distributed vector of the 
     indeces corresponding to an element vertex 
   */
  
  int d ;
  if( dim == 2 )
  {
    int n0;
    if(  rank_mic == 0 ){
      for( d = 0 ; d < dim ; d++ ){
	n0 = (e%nex) + (e/nex)*nx;
	loc_elem_index[ 0*dim + d ] = ( n0          ) * dim + d ;
	loc_elem_index[ 1*dim + d ] = ( n0 + 1      ) * dim + d ;
	loc_elem_index[ 2*dim + d ] = ( n0 + nx + 1 ) * dim + d ;
	loc_elem_index[ 3*dim + d ] = ( n0 + nx + 0 ) * dim + d ;
      }
    }
    else if(e >= nex ){
      for( d = 0 ; d < dim ; d++ ){
	n0 = (e%nex) + (e/nex-1)*nx;
	loc_elem_index[ 0*dim + d ] = ( n0          ) * dim + d ;
	loc_elem_index[ 1*dim + d ] = ( n0 + 1      ) * dim + d ;
	loc_elem_index[ 2*dim + d ] = ( n0 + nx + 1 ) * dim + d ;
	loc_elem_index[ 3*dim + d ] = ( n0 + nx     ) * dim + d ;
      }
    }
    else{
      for( d = 0 ; d < dim ; d++ ){
	n0 = (e%nex) + (e/nex)*nx;
	loc_elem_index[ 0*dim + d ] = ( n0 + nl     ) * dim + d ; // is a ghost
	loc_elem_index[ 1*dim + d ] = ( n0 + nl + 1 ) * dim + d ; // is a ghost
	loc_elem_index[ 2*dim + d ] = ( n0 + 1      ) * dim + d ;
	loc_elem_index[ 3*dim + d ] = ( n0          ) * dim + d ;
      }
    }
  }
  return 0;
}

/****************************************************************************************************/

int get_local_elem_node( int e , int *n_loc )
{
  if( dim == 2 )
  {
    int n0;
    if(  rank_mic == 0 ){
      n0 = (e%nex) + (e/nex)*nx;
      n_loc[0] = n0          ;
      n_loc[1] = n0 + 1      ;
      n_loc[2] = n0 + nx + 1 ;
      n_loc[3] = n0 + nx     ;
    }
    else if(e >= nex ){
      n0 = (e%nex) + (e/nex-1)*nx;
      n_loc[0] = n0          ;
      n_loc[1] = n0 + 1      ;
      n_loc[2] = n0 + nx + 1 ;
      n_loc[3] = n0 + nx     ;
    }
    else{
      n0 = (e%nex) + (e/nex)*nx;
      n_loc[0] = n0 + nl     ; // is a ghost
      n_loc[1] = n0 + nl + 1 ; // is a ghost
      n_loc[2] = n0 + 1      ;
      n_loc[3] = n0          ;
    }
  }
  return 0;
}

/****************************************************************************************************/

int get_global_elem_index( int e, int *glo_elem_index )
{

  /* 
     returns the local position in the distributed vector of the 
     indeces corresponding to an element vertex 
   */
  
  int d;
  if( dim == 2 )
  {
    int n0;
    if(  rank_mic == 0 ){
      for( d = 0 ; d < dim ; d++ ){
	n0 = (e%nex) + (e/nex)*nx;
	glo_elem_index[ 0*dim + d ] = ( n0          ) * dim + d ;
	glo_elem_index[ 1*dim + d ] = ( n0 + 1      ) * dim + d ;
	glo_elem_index[ 2*dim + d ] = ( n0 + nx + 1 ) * dim + d ;
	glo_elem_index[ 3*dim + d ] = ( n0 + nx     ) * dim + d ;
      }
    }
    else{
      for( d = 0 ; d < dim ; d++ ){
	n0 = (e%nex) + (e/nex)*nx + (ny_inf-1)*nx;
	glo_elem_index[ 0*dim + d ] = ( n0          ) * dim + d ;
	glo_elem_index[ 1*dim + d ] = ( n0 + 1      ) * dim + d ;
	glo_elem_index[ 2*dim + d ] = ( n0 + nx + 1 ) * dim + d ;
	glo_elem_index[ 3*dim + d ] = ( n0 + nx     ) * dim + d ;
      }
    }
  }
  return 0;
}

/****************************************************************************************************/

int local_to_global_index( int local )
{

  if ( rank_mic == 0 ) return local;

  if( dim == 2 ){

    return ny_inf * nx * dim + local;

  }

  return 0;
}

/****************************************************************************************************/

int mic_homogenize(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  /*
     Performs linear homogenization of the RVE 

     UNIF_STRAINS            > u = E . x

     HOMO_TAYLOR_S
   */
  int ierr;

  if(homo_type==TAYLOR_P || homo_type==TAYLOR_S){

    ierr = mic_homogenize_taylor(MICRO_COMM, strain_mac, strain_ave, stress_ave);
    if(ierr){
      return 1;
    }
  }
  else if(homo_type==UNIF_STRAINS){

    ierr = mic_homog_us(MICRO_COMM, strain_mac, strain_ave, stress_ave);
    if(ierr) return 1;

  }
  if(first_time_homo) first_time_homo = 0;

  return 0;
}

/****************************************************************************************************/

int mic_calc_c_homo(MPI_Comm MICRO_COMM, double strain_mac[6], double c_homo[36])
{

  /* 
     Si la micro estructura está integramente conformada por materiales
     lineales entonces este tensor será siempre el mismo para cada punto 
     de gauss en la macro escala entonces es eficiente almacenar c_homo_linear
  */

  int i, ierr;

  if(flag_linear_micro){

    if(first_time_homo){
      ierr = mic_calc_c_homo_lineal(MICRO_COMM, c_homo_lineal);
      if(ierr){
	return 1;
      }
    }
    for(i=0;i<nvoi*nvoi;i++){
      c_homo[i] = c_homo_lineal[i];
    }

  }
  else{
    return 1;
  }
  return 0;
}

/****************************************************************************************************/

int mic_calc_c_homo_lineal(MPI_Comm MICRO_COMM, double c_homo_lineal[36])
{

  int       i, j;
  int       ierr;
  double    strain[6];
  double    strain_ave[6];
  double    stress_ave[6];

  for(i=0;i<nvoi*nvoi;i++){
    c_homo_lineal[i]=0.0;
  }

  for(i=0;i<nvoi;i++){

    for(j=0;j<nvoi;j++){
      strain[j]=0.0;
    }
    strain[i]=0.005;

    ierr = mic_homogenize(MICRO_COMM, strain, strain_ave, stress_ave);
    if(ierr){
      return 1;
    }

    for(j=0;j<nvoi;j++){
      c_homo_lineal[j*nvoi+i] = stress_ave[j] / strain_ave[i];
    }

  }

  return 0;
}

/****************************************************************************************************/

int mic_calc_stress_ave(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  /* 
     Si la micro estructura está integramente conformada por materiales
     lineales entonces este tensor será siempre el mismo para cada punto 
     de gauss en la macro escala entonces es eficiente almacenar c_homo_linear
  */

  int i, j, ierr;

  if(flag_linear_micro){

    for(i=0;i<nvoi;i++){
      strain_ave[i] = strain_mac[i];
      stress_ave[i] = 0.0;
      for( j = 0 ; j < nvoi ; j++ )
	stress_ave[i] += c_homo_lineal[i*nvoi+j] * strain_mac[j];
    }
    return 0;

  }
  else{

    ierr = mic_homogenize(MICRO_COMM, strain_mac, strain_ave, stress_ave);
    if(ierr){
      return 1;
    }
    return 0;

  }
  return 0;
}

/****************************************************************************************************/

int mic_check_linear_material(void)
{

  flag_linear_micro = 1;

  return 0;
}

/****************************************************************************************************/

int get_node_local_coor( int n , double * coord )
{
  if( dim == 2 ){
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

/****************************************************************************************************/

int get_node_ghost_coor( int n , double * coord )
{
  if( dim == 2 ){
    coord[0] = (n % nx)*hx;
    coord[1] = (ny_inf + n / nx - 1)*hy;
  }
  return 0;
}

/****************************************************************************************************/

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

/****************************************************************************************************/
