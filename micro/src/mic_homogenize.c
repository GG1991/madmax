/*
   Routines for performing homogenization on RVE

   Author > Guido Giuntoli
   Date   > 18-08-2017

 */

#include "micro.h"

int get_local_elem_index( int e, int *loc_index );
int get_global_elem_index( int e, int *glo_elem_index );
int assembly_residual_struct( void );
int assembly_jacobian_struct( void );
int get_stress( int e , int gp, double *strain_gp , double *stress_gp );
int get_strain( int e , int gp, double *strain_gp );
int get_c_tan( int e , int gp, double *strain_gp , double *c_tan );
int get_centroid_struct( int e, double *centroid );
int is_in_fiber( int e );
int strain_x_coord( double * strain , double * coord , double * u );
int get_node_local_coor( int n , double * coord );
int get_node_ghost_coor( int n , double * coord );
int get_local_elem_node( int e , int * n_loc );
int local_to_global_index( int local );
int init_shapes(double ***sh, double ****dsh, double **wp, double *h, int dim);

int mic_homogenize_taylor(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{
  /* 
     This work only if inclusion names starts with FIBER and matrix name start with MATRIX
     we assume only on kind of them 
    
     This kind of homogenization do not need to perform 3 (dim=2) or 6 (dim=3) experiments
  */
  int          i, j, e, ierr;
  double       c_i[6][6];
  double       c_m[6][6];
  double       c[6][6];
  double       vol_i = 0.0, vol_ia = 0.0;  // inclusion volume
  double       vol_m = 0.0, vol_ma = 0.0;  // matrix volume
  double       vol_t = 0.0;  // total volume
  double       vol_e = 0.0;  // element volume
  material_t   *mat = NULL; 

  if(first_time_homo){
    for(e=0;e<nelm;e++){

      ierr = get_mat_from_elem(e, &mat);
      ierr = get_elem_vol(e, &vol_e);
      if(!strncmp(mat->name,"FIBER",5)){
	vol_ia += vol_e;
      }
      else if(!strncmp(mat->name,"MATRIX",6)){
	vol_ma += vol_e;
      }
      else{
	PetscPrintf(MICRO_COMM,	"not possible to compute a Taylor homogenization with a material of type\n",mat->name);
	return 1;
      }

    }
    ierr = MPI_Allreduce(&vol_ia, &vol_i, 1, MPI_DOUBLE, MPI_SUM, MICRO_COMM);CHKERRQ(ierr);
    ierr = MPI_Allreduce(&vol_ma, &vol_m, 1, MPI_DOUBLE, MPI_SUM, MICRO_COMM);CHKERRQ(ierr);
    vol_t = vol_i + vol_m;
    vi = vol_i / vol_t; // inclusion fraction
    vm = vol_m / vol_t; // matrix fraction
  }
  get_c("FIBER" , 0, 0, strain_mac, c_i);  //returns c_i of FIBER
  get_c("MATRIX", 0, 0, strain_mac, c_m);  //returns c_m of MATRIX
 
  if(homo_type==TAYLOR_P){
    for(i=0;i<nvoi;i++){
      for(j=0;j<nvoi;j++){
	c[i][j] = vi * c_i[i][j] + vm * c_m[i][j];
      }
    }
  }
  else if(homo_type==TAYLOR_S){

    int              s;
    double           c_ia[36], c_ma[36]; // matrices in array
    double           c_ii[36], c_mi[36]; // inverted matrices
    double           c_a[36] , c_ai[36];
    gsl_matrix_view  m, mi;
    gsl_permutation  *p;

    p  = gsl_permutation_alloc(nvoi);

    for(i = 0; i < nvoi; ++i){
      for(j = 0; j < nvoi; ++j){
	c_ia[i*nvoi + j] = c_i[i][j];
	c_ma[i*nvoi + j] = c_m[i][j];
      }
    }
    
    m  = gsl_matrix_view_array(c_ia,nvoi,nvoi);
    mi = gsl_matrix_view_array(c_ii,nvoi,nvoi);
    
    gsl_linalg_LU_decomp (&m.matrix, p, &s);    
    gsl_linalg_LU_invert (&m.matrix, p, &mi.matrix);
    
    for(i = 0; i < nvoi; ++i){
      for(j = 0; j < nvoi; ++j){
	c_ii[i*nvoi+j] = gsl_matrix_get(&mi.matrix,i,j);
      }
    }

    m  = gsl_matrix_view_array(c_ma,nvoi,nvoi);
    mi = gsl_matrix_view_array(c_mi,nvoi,nvoi);
    
    gsl_linalg_LU_decomp (&m.matrix, p, &s);    
    gsl_linalg_LU_invert (&m.matrix, p, &mi.matrix);
    
    for(i = 0; i < nvoi; ++i){
      for(j = 0; j < nvoi; ++j){
	c_mi[i*nvoi+j] = gsl_matrix_get(&mi.matrix,i,j);
      }
    }
    for(i = 0; i < nvoi*nvoi; ++i){
      c_a[i] = vi * c_ii[i] + vm * c_mi[i];
    }

    m  = gsl_matrix_view_array(c_a ,nvoi,nvoi);
    mi = gsl_matrix_view_array(c_ai,nvoi,nvoi);
    
    gsl_linalg_LU_decomp (&m.matrix, p, &s);    
    gsl_linalg_LU_invert (&m.matrix, p, &mi.matrix);
    
    for(i = 0; i < nvoi; ++i){
      for(j = 0; j < nvoi; ++j){
	c[i][j] = gsl_matrix_get(&mi.matrix,i,j);
      }
    }
     
    gsl_permutation_free (p);
  }
  for(i=0;i<nvoi;i++){
    strain_ave[i] = strain_mac[i];
    stress_ave[i] = 0.0;
    for(j=0;j<nvoi;j++){
      stress_ave[i] += c[i][j] * strain_mac[j];
    }
  }

  return 0;
}
/****************************************************************************************************/
int mic_homogenize_unif_strains(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{
  /*
     ub = E . x
  */

  int    ierr, nr_its=-1;
  int    nnods_bc, *nods_bc;
  int    i, d, k;
  int    *ix_loc, *ix_glo;
  double *ub;
  double norm=2*nr_norm_tol;
  double strain_matrix[3][3];
  double *x_arr, *b_arr;

  if(dim==2){
    strain_matrix[0][0]=strain_mac[0]    ; strain_matrix[0][1]=strain_mac[2]*0.5;
    strain_matrix[1][0]=strain_mac[2]*0.5; strain_matrix[1][1]=strain_mac[1]    ;
  }
  else if(dim==3){
    strain_matrix[0][0]=strain_mac[0]; strain_matrix[0][1]=strain_mac[3]; strain_matrix[0][2]=strain_mac[5];
    strain_matrix[1][0]=strain_mac[3]; strain_matrix[1][1]=strain_mac[1]; strain_matrix[1][2]=strain_mac[4];
    strain_matrix[2][0]=strain_mac[5]; strain_matrix[2][1]=strain_mac[4]; strain_matrix[2][2]=strain_mac[2];
  }

  /* get array of nods in al boundary without repetition */
  ierr = get_nods_bc( &nods_bc, &nnods_bc);
  if(ierr){
    return 1;
  }
  ix_loc = malloc(nnods_bc*dim*sizeof(int));    /* indeces to search for local coordinates */
  ix_glo = malloc(nnods_bc*dim*sizeof(int));    /* indeces to be set at boundary */
  ub     = malloc(nnods_bc*dim*sizeof(double)); /* values  to be set at boundary */
  ierr   = get_nods_index( nods_bc, nnods_bc, ix_loc, ix_glo);

  // first we fill the value ub_val in <homog_ld_lagran_t> structure
  for(i=0; i<nnods_bc; i++){
    for(d=0;d<dim;d++){
      ub[i*dim+d] = 0.0;
      for(k=0;k<dim;k++){
	ub[ i * dim + d] = ub[ i * dim + d] +  strain_matrix[d][k] * coord[ix_loc[i* dim + k]];
      }
    }
  }

  // first we fill the value ub_val in <homog_ld_lagran_t> structure
  ierr = VecGetArray(x, &x_arr);CHKERRQ(ierr);
  for(i=0; i<nnods_bc; i++){
    for(d=0;d<dim;d++){
      x_arr[ix_loc[i*dim+d]] = ub[i*dim+d]; 
    }
  }
  ierr = VecRestoreArray(x, &x_arr);CHKERRQ(ierr);

  while( nr_its < nr_max_its && norm > nr_norm_tol )
  {
    ierr = assembly_residual_sd(&x, &b);CHKERRQ(ierr);
    ierr = VecGetArray(b, &b_arr);CHKERRQ(ierr);
    for(i=0; i<nnods_bc; i++){
      for(d=0;d<dim;d++){
	b_arr[ix_loc[i*dim+d]] = 0.0;
      }
    }
    ierr = VecRestoreArray(b, &b_arr);CHKERRQ(ierr);
    ierr = VecScale(b,-1.0); CHKERRQ(ierr);

    ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
    if(!flag_coupling)
      PetscPrintf(MICRO_COMM,"|b| = %lf \n",norm);

    if( !(norm > nr_norm_tol) )break;

    if(first_time_homo || !flag_linear_micro){
      ierr = assembly_jacobian_sd(&A);
    }
    ierr = MatZeroRowsColumns(A, nnods_bc*dim, ix_glo, 1.0, NULL, NULL); CHKERRQ(ierr);
    ierr = KSPSolve(ksp,b,dx);CHKERRQ(ierr);

    ierr = VecAXPY( x, 1.0, dx); CHKERRQ(ierr);
    ierr = MatMult(A, dx, b1);CHKERRQ(ierr);
    
    if( flag_print & (1<<PRINT_PETSC) ){
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"A.dat",&viewer); CHKERRQ(ierr);
      ierr = MatView(A,viewer); CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"b.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(b,viewer); CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"b1.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(b1,viewer); CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(x,viewer); CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"dx.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(x,viewer); CHKERRQ(ierr);
    }

    nr_its ++;
  }

  ierr = calc_ave_strain_stress(MICRO_COMM, &x, strain_ave, stress_ave);CHKERRQ(ierr);

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
    ny_inf = ( dim == 2 ) ? nstart / nx : nstart / (nx*nz);

    int *ghost_index;;
    if( nproc_mic == 1 )
      ngho = 0;
    else{
      if( rank_mic == 0 )
	ngho = 0;
      else
	ngho = ( dim == 2 ) ? nx : nx*nz;
    }

    ghost_index = malloc( ngho*dim *sizeof(int) );

    int i, j;
    for( i = 0 ; i < ngho*dim  ; i++ )
      ghost_index[i] = istart - (( dim == 2 )? nx : nx*nz)*dim + i;

    VecCreateGhost(MICRO_COMM, nl*dim, nn*dim, ngho*dim, ghost_index, &x);
    VecDuplicate(x,&dx);
    VecDuplicate(x,&b);

    free(ghost_index);

    /* Initilize shape functions, derivatives, jacobian, b_matrix */

    double h[3]; h[0] = hx; h[1] = hx; h[2] = hz;

    init_shapes( &struct_sh, &struct_dsh, &struct_wp, h, dim);
    
    /* alloc the B matrix */
    struct_bmat = malloc( nvoi * sizeof(double**));
    for( i = 0 ; i < nvoi  ; i++ ){
      struct_bmat[i] = malloc( npe*dim * sizeof(double*));
      for( j = 0 ; j < npe*dim ; j++ )
	struct_bmat[i][j] = malloc( ngp * sizeof(double));
    }

    /* calc B matrix */
    int is, gp;
    for( gp = 0; gp < ngp ; gp++ ){
      for( is = 0; is < npe ; is++ ){
	if( dim == 2 ){
	  struct_bmat[0][is*dim + 0][gp] = struct_dsh[is][0][gp];
	  struct_bmat[0][is*dim + 1][gp] = 0;
	  struct_bmat[1][is*dim + 0][gp] = 0;
	  struct_bmat[1][is*dim + 1][gp] = struct_dsh[is][1][gp];
	  struct_bmat[2][is*dim + 0][gp] = struct_dsh[is][1][gp];
	  struct_bmat[2][is*dim + 1][gp] = struct_dsh[is][0][gp];
	}
      }
    }

    /* alloc the local and global index vector for assembly */
    loc_elem_index = malloc( dim*npe * sizeof(int));
    glo_elem_index = malloc( dim*npe * sizeof(int));
    elem_disp = malloc( dim*npe * sizeof(double));
    stress_gp = malloc( nvoi * sizeof(double) );
    strain_gp = malloc( nvoi * sizeof(double) );

    /* alloc arrays for boundary condition setting */
    if(nproc_mic == 1){
      ndir_ix = 2*nx + 2*(nyl-2);
    }
    else{
      if( rank_mic == 0 || rank_mic == (nproc_mic - 1) ){
	ndir_ix = nx + (nyl - 1)*2;
      }
      else{ 
	ndir_ix = nyl * 2;
      }
    }
    ndir_ix    *= dim;
    dir_ix_loc  = malloc ( ndir_ix * sizeof(int));
    dir_ix_glo  = malloc ( ndir_ix * sizeof(int));
    coor_dir    = malloc ( ndir_ix * sizeof(double));

    int n, d, c;

    c = 0;
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
      if( rank_mic == 0 || rank_mic == (nproc_mic - 1) )
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
      else
      {
	/* 0 < rank_mic < (nproc - 1) */
	/* x = 0 */
	for( n = 0 ; n < nyl ; n++ ){
	  for( d = 0 ; d < dim ; d++ )
	    dir_ix_loc[c*dim + d] = n*nx*dim + d;
	  coor_dir[c*dim + 0] = 0;
	  coor_dir[c*dim + 1] = (ny_inf + n + 1)*hy;
	  c++;
	}

	/* x = lx */
	for( n = 0 ; n < nyl ; n++ ){
	  for( d = 0 ; d < dim ; d++ )
	    dir_ix_loc[c*dim + d] = (nx-1)*dim + n*nx*dim + d;
	  coor_dir[c*dim + 0] = lx;
	  coor_dir[c*dim + 1] = (ny_inf + n + 1)*hy;
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

  Vec     x_loc;
  double  *x_arr;

  VecZeroEntries( x );
  VecGhostGetLocalForm( x , &x_loc );
  VecGetArray( x_loc, &x_arr );

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

  VecRestoreArray( x_loc , &x_arr );
  VecGhostRestoreLocalForm( x , &x_loc );
  VecGhostUpdateBegin ( x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd   ( x, INSERT_VALUES, SCATTER_FORWARD);

  int     i;
  int     nr_its = 0; 
  double  *b_arr;
  double  norm = nr_norm_tol*10;

  VecNorm( x , NORM_2 , &norm );
  if(!flag_coupling)
    PetscPrintf(MICRO_COMM,"|x| = %lf \n",norm);

  while( nr_its < nr_max_its && norm > nr_norm_tol )
  {
    assembly_residual_struct(); // assembly "b" (residue) using "x" (displacement)

    VecGetArray( b, &b_arr );
    for( i = 0; i < ndir_ix ; i++ )
      b_arr[dir_ix_loc[i]] = 0.0;
    VecRestoreArray( b, &b_arr );
    VecGhostUpdateBegin ( b, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd   ( b, INSERT_VALUES, SCATTER_FORWARD);
    VecNorm( b , NORM_2 , &norm );
    if(!flag_coupling)
      PetscPrintf(MICRO_COMM,"|b| = %lf \n",norm);

    if( !(norm > nr_norm_tol) ) break;

    VecScale( b, -1.0 );
    assembly_jacobian_struct(); // assembly "A" (jacobian) using "x" (displacement)
    MatZeroRowsColumns(A, ndir_ix, dir_ix_glo, 1.0, NULL, NULL);

    KSPSetOperators(ksp,A,A);
    KSPSolve( ksp, b, dx );
    VecAXPY( x, 1.0, dx);
    VecGhostUpdateBegin ( x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd   ( x, INSERT_VALUES, SCATTER_FORWARD);

    nr_its ++;
  }

  if( flag_print & (1<<PRINT_PETSC) ){
    PetscViewerASCIIOpen(MICRO_COMM,"A.dat",&viewer);
    MatView(A,viewer);
    PetscViewerASCIIOpen(MICRO_COMM,"b.dat",&viewer);
    VecView(b,viewer);
    PetscViewerASCIIOpen(MICRO_COMM,"x.dat",&viewer);
    VecView(x,viewer);
    PetscViewerASCIIOpen(MICRO_COMM,"dx.dat",&viewer);
    VecView(x,viewer);
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
int assembly_residual_struct(void)
{

  VecZeroEntries(b);

  /* from the owning processes to the ghosts in the others processes */
  VecGhostUpdateBegin( b , INSERT_VALUES , SCATTER_FORWARD );
  VecGhostUpdateEnd(   b , INSERT_VALUES , SCATTER_FORWARD );

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

  VecRestoreArray( b_loc , &b_arr);
  VecGhostRestoreLocalForm( b , &b_loc );

  /* from the local and ghost part with add to all processes */
  VecGhostUpdateBegin( b, ADD_VALUES, SCATTER_REVERSE );
  VecGhostUpdateEnd  ( b, ADD_VALUES, SCATTER_REVERSE );

  return 0;
}
/****************************************************************************************************/
int assembly_jacobian_struct( void )
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
      get_c_tan( e , gp , strain_gp , c );

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
int calc_stress_strain_energy( double * stress, double * strain, double * energy )
{

  int e, v, gp;

  for ( e = 0 ; e < nelm ; e++ ){

    for ( v = 0 ; v < nvoi ; v++ )
      strain[ e*nvoi + v ] = stress[ e*nvoi + v ] = 0.0;

    for ( gp = 0 ; gp < ngp ; gp++ ){

      get_strain( e , gp, strain_gp );
      get_stress( e , gp, strain_gp, stress_gp );
      for ( v = 0 ; v < nvoi ; v++ ){
	strain[ e*nvoi + v ] += strain_gp[v] * struct_wp[gp];
	stress[ e*nvoi + v ] += stress_gp[v] * struct_wp[gp];
      }

    }
  }

  return 0;
}
/****************************************************************************************************/
int get_strain( int e , int gp, double *strain_gp )
{

  Vec     x_loc;
  double  *x_arr;

  VecGhostGetLocalForm( x , &x_loc );
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

  material_t  *mat_p;
  node_list_t *pm;

  switch( micro_type )
  {
    case CIRCULAR_FIBER:
      if(is_in_fiber( e ))
      { /* is in the fiber */
	pm = material_list.head;
	while( pm ){
          /* search FIBER */
	  mat_p = (material_t *)pm->data;
	  if( strcmp ( mat_p->name , "FIBER" ) == 0 ) break;
	  pm = pm->next;
	}
	if( !pm ) return 1;
      }
      else
      { /* is in the matrix */
	pm = material_list.head;
	while( pm ){
          /* search MATRIX */
	  mat_p = (material_t *)pm->data;
	  if( strcmp ( mat_p->name , "MATRIX" ) == 0 ) break;
	  pm = pm->next;
	}
	if( !pm ) return 1;
      }

      break;

    default:
      return 1;

  }

  /* now that we now the material (mat_p) we calculate stress = f(strain) */

  if( mat_p->type_id == TYPE_0 ){

    /* is a linear material stress = C * strain */

    double  young   = ((type_0*)mat_p->type)->young;
    double  poisson = ((type_0*)mat_p->type)->poisson;
    double  c[3][3];
    int     i , j;

    if(dim==2){

      /* Plane strain ( long fibers case ) */
      c[0][0]=1.0-poisson; c[0][1]=poisson    ; c[0][2]=0.0            ;
      c[1][0]=poisson    ; c[1][1]=1.0-poisson; c[1][2]=0.0            ;
      c[2][0]=0.0        ; c[2][1]=0.0        ; c[2][2]=(1-2*poisson)/2;

      for( i = 0; i < nvoi ; i++ ){
	for( j = 0 ; j < nvoi ; j++ )
	  c[i][j] *= young/((1+poisson)*(1-2*poisson));
      }
      for( i = 0; i < nvoi ; i++ ){
	stress_gp[i] = 0.0;
	for( j = 0 ; j < nvoi ; j++ )
	  stress_gp[i] += c[i][j] * strain_gp[j];
      }

    }
  }

  return 0;
}
/****************************************************************************************************/
int get_c_tan( int e , int gp, double *strain_gp , double *c_tan )
{

  /* returns the stress according to the elemet type */

  material_t  *mat_p;
  node_list_t *pm;

  switch( micro_type ){

    case CIRCULAR_FIBER:

      if(is_in_fiber( e ))
      {	/* is in the fiber */
	pm = material_list.head;
	while( pm ){
          /* search FIBER */
	  mat_p = (material_t *)pm->data;
	  if( strcmp ( mat_p->name , "FIBER" ) == 0 ) break;
	  pm = pm->next;
	}
	if( !pm ) return 1;
      }
      else
      { /* is in the matrix */
	pm = material_list.head;
	while( pm ){
          /* search MATRIX */
	  mat_p = (material_t *)pm->data;
	  if( strcmp ( mat_p->name , "MATRIX" ) == 0 ) break;
	  pm = pm->next;
	}
	if( !pm ) return 1;
      }
      break;

    default:
      return 1;

  }

  /* now that we now the material (mat_p) we calculate stress = f(strain) */

  if( mat_p->type_id == TYPE_0 ){

    /* is a linear material stress = C * strain */
    double  young   = ((type_0*)mat_p->type)->young;
    double  poisson = ((type_0*)mat_p->type)->poisson;
    int     i , j;

    if(dim==2){

      /* Plane strain ( long fibers case ) */
      c_tan[0*nvoi+0]=1.0-poisson; c_tan[0*nvoi+1]=poisson    ; c_tan[0*nvoi+2]=0.0            ;
      c_tan[1*nvoi+0]=poisson    ; c_tan[1*nvoi+1]=1.0-poisson; c_tan[1*nvoi+2]=0.0            ;
      c_tan[2*nvoi+0]=0.0        ; c_tan[2*nvoi+1]=0.0        ; c_tan[2*nvoi+2]=(1-2*poisson)/2;

      for( i = 0; i < nvoi ; i++ ){
	for( j = 0 ; j < nvoi ; j++ )
	  c_tan[i*nvoi + j] *= young/((1+poisson)*(1-2*poisson));
      }

    }
  }

  return 0;
}
/****************************************************************************************************/
int is_in_fiber( int e )
{

  int    i, j, d;
  double l = -1, centroid[3];
  double deviation[2];

  get_centroid_struct( e, centroid );

  for( i = 0 ; i < nx_fibers ; i++ ){
    for( j = 0 ; j < ny_fibers ; j++ ){
      deviation[0] = fiber_cilin_center_devi[0] - LX/2 + (lx/nx_fibers)/2 + i*(lx/nx_fibers);
      deviation[1] = fiber_cilin_center_devi[1] - LY/2 + (ly/ny_fibers)/2 + j*(ly/ny_fibers);
      l = 0.0;
      for(d=0;d<2;d++){
	l = l + pow(centroid[d]-(center_domain[d]+deviation[d]),2);
      }
      l = sqrt(l);
      if (l<=fiber_cilin_r) return 1;
    }
  }
  return 0;

}
/****************************************************************************************************/
int get_centroid_struct( int e, double *centroid )
{

  /* formula only valid for sequencial now */

  if( dim == 2 ){
    centroid[0] = ( e % nex + 0.5         ) * hx;
    centroid[1] = ( e / nex + 0.5 + ny_inf) * hy;
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
    if(  rank_mic == 0 ){
      for( d = 0 ; d < dim ; d++ ){
	loc_elem_index[ 0*dim + d ] = ( (e%nex) + (e/nex)*nx + 0 )      * dim + d ;
	loc_elem_index[ 1*dim + d ] = ( (e%nex) + (e/nex)*nx + 1 )      * dim + d ;
	loc_elem_index[ 2*dim + d ] = ( (e%nex) + (e/nex)*nx + nx + 1 ) * dim + d ;
	loc_elem_index[ 3*dim + d ] = ( (e%nex) + (e/nex)*nx + nx + 0 ) * dim + d ;
      }
    }
    else if(e >= nex ){
      for( d = 0 ; d < dim ; d++ ){
	loc_elem_index[ 0*dim + d ] = ( (e%nex) + (e/nex-1)*nx + 0 )      * dim + d ;
	loc_elem_index[ 1*dim + d ] = ( (e%nex) + (e/nex-1)*nx + 1 )      * dim + d ;
	loc_elem_index[ 2*dim + d ] = ( (e%nex) + (e/nex-1)*nx + nx + 1 ) * dim + d ;
	loc_elem_index[ 3*dim + d ] = ( (e%nex) + (e/nex-1)*nx + nx + 0 ) * dim + d ;
      }
    }
    else{
      for( d = 0 ; d < dim ; d++ ){
	loc_elem_index[ 0*dim + d ] = ( (e%nex) + (e/nex)*nx + 0 + nl ) * dim + d ; // is a ghost
	loc_elem_index[ 1*dim + d ] = ( (e%nex) + (e/nex)*nx + 1 + nl ) * dim + d ; // is a ghost
	loc_elem_index[ 2*dim + d ] = ( (e%nex) + (e/nex)*nx + 1 )      * dim + d ;
	loc_elem_index[ 3*dim + d ] = ( (e%nex) + (e/nex)*nx + 0 )      * dim + d ;
      }
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
    if(  rank_mic == 0 ){
      for( d = 0 ; d < dim ; d++ ){
	glo_elem_index[ 0*dim + d ] = ( (e%nex) + (e/nex)*nx + 0      ) * dim + d ;
	glo_elem_index[ 1*dim + d ] = ( (e%nex) + (e/nex)*nx + 1      ) * dim + d ;
	glo_elem_index[ 2*dim + d ] = ( (e%nex) + (e/nex)*nx + nx + 1 ) * dim + d ;
	glo_elem_index[ 3*dim + d ] = ( (e%nex) + (e/nex)*nx + nx + 0 ) * dim + d ;
      }
    }
    else{
      for( d = 0 ; d < dim ; d++ ){
	glo_elem_index[ 0*dim + d ] = ( (e%nex) + (e/nex)*nx + (ny_inf-1)*nx + 0      ) * dim + d ;
	glo_elem_index[ 1*dim + d ] = ( (e%nex) + (e/nex)*nx + (ny_inf-1)*nx + 1      ) * dim + d ;
	glo_elem_index[ 2*dim + d ] = ( (e%nex) + (e/nex)*nx + (ny_inf-1)*nx + nx + 1 ) * dim + d ;
	glo_elem_index[ 3*dim + d ] = ( (e%nex) + (e/nex)*nx + (ny_inf-1)*nx + nx + 0 ) * dim + d ;
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

    return ny_inf*nx*dim + local;

  }

  return 0;
}
/****************************************************************************************************/
int get_local_elem_node( int e , int *n_loc )
{
  if( dim == 2 )
  {
    if( e >= nex || rank_mic == 0 ){
      /* does not have ghost values */
      n_loc[0] = ( (e%nex) + (e/nex)*nx + 0)      ;
      n_loc[1] = ( (e%nex) + (e/nex)*nx + 1)      ;
      n_loc[2] = ( (e%nex) + (e/nex)*nx + nx + 1) ;
      n_loc[3] = ( (e%nex) + (e/nex)*nx + nx + 0) ;
    }
    else{
      /* has ghost values in the bottom */
      n_loc[0] = ( (e%nex) + (e/nex)*nx + 0) + nl ;
      n_loc[1] = ( (e%nex) + (e/nex)*nx + 1) + nl ;
      n_loc[2] = ( (e%nex) + (e/nex)*nx + 1) ;
      n_loc[3] = ( (e%nex) + (e/nex)*nx + 0) ;
    }
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
int micro_pvtu(MPI_Comm PROBLEM_COMM, char *name, double *strain, double *stress, double *energy)
{

  FILE    *fm;
  char    file_name[NBUF];
  double  *xvalues;
  Vec     xlocal;

  int rank, nproc;
  MPI_Comm_size(PROBLEM_COMM, &nproc);
  MPI_Comm_rank(PROBLEM_COMM, &rank);

  if( rank == 0 ){

    /* rank 0 writes the .pvtu file first */
    strcpy(file_name,name);
    strcat(file_name,".pvtu");
    fm = fopen(file_name,"w");

    fprintf(fm, "<?xml version=\"1.0\"?>\n"
	"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
	"<PUnstructuredGrid GhostLevel=\"0\">\n"
	"<PPoints>\n"
	"<PDataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\"/>\n"
	"</PPoints>\n"
	"<PCells>\n"
	"<PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>\n"
	"</PCells>\n" 

	"<PPointData Vectors=\"displ\">\n" 
	"<PDataArray type=\"Float64\" Name=\"displ\"    NumberOfComponents=\"3\" />\n"
	"<PDataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" />\n"
	"</PPointData>\n"

	"<PCellData>\n"
	"<PDataArray type=\"Int32\"   Name=\"part\" NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\"/>\n"
	"</PCellData>\n" , nvoi , nvoi);

    int i;
    for( i = 0 ; i < nproc ; i++ ){
      sprintf(file_name,"%s_%d",name,i);
      fprintf(fm,	"<Piece Source=\"%s.vtu\"/>\n",file_name);
    }
    fprintf(fm,	"</PUnstructuredGrid>\n" 
      "</VTKFile>\n" );

    fclose(fm);

  } // rank = 0

  sprintf(file_name,"%s_%d.vtu",name,rank);
  fm = fopen(file_name,"w"); 
  if(!fm){
    PetscPrintf(PETSC_COMM_WORLD,"Problem trying to opening file %s for writing\n", file_name);
    return 1;
  }

  fprintf(fm, 
      "<?xml version=\"1.0\"?>\n"
      "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      "<UnstructuredGrid>\n");
  fprintf(fm,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nl + ngho, nelm);
  fprintf(fm,"<Points>\n");
  fprintf(fm,"<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");

  double *coord = malloc( dim * sizeof(double));
  int    n , d;

  for( n = 0 ; n < nl ; n++ ){
    get_node_local_coor( n , coord );
    for( d = 0 ; d < dim ; d++ )
      fprintf(fm,"%e ",  coord[d] );
    for( d = dim ; d < 3 ; d++ )
      fprintf(fm,"%e ",0.0);
    fprintf(fm,"\n");
  }
  for( n = 0 ; n < ngho ; n++ ){
    get_node_ghost_coor( n , coord );
    for( d = 0 ; d < dim ; d++ )
      fprintf(fm,"%e ",  coord[d] );
    for( d = dim ; d < 3 ; d++ )
      fprintf(fm,"%e ",0.0);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");
  fprintf(fm,"</Points>\n");
  fprintf(fm,"<Cells>\n");

  int e;
  int * n_loc = malloc( npe * sizeof(int));

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for ( e = 0 ; e < nelm ; e++ ){
    get_local_elem_node( e , n_loc );
    for ( n = 0 ; n < npe ; n++ )
      fprintf(fm,"%d ", n_loc[n]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  int ce = npe;
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for ( e = 0 ; e < nelm ; e++ ){
    fprintf(fm,"%d ", ce);
    ce += npe;
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for ( e = 0 ; e < nelm ; e++ )
    fprintf(fm, "%d ",vtkcode( dim , npe ) );  
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"</Cells>\n");
  
  fprintf(fm,"<PointData Vectors=\"displ\">\n"); // Vectors inside is a filter we should not use this here

  /* <displ> */
  VecGhostUpdateBegin( x , INSERT_VALUES , SCATTER_FORWARD);
  VecGhostUpdateEnd(   x , INSERT_VALUES , SCATTER_FORWARD);
  VecGhostGetLocalForm(x , &xlocal );

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
  VecGetArray( xlocal , &xvalues );
  for( n = 0 ; n < (nl + ngho) ; n++ ){
    for( d = 0 ; d < dim ; d++ )
      fprintf(fm, "%lf ", xvalues[ n * dim + d ]);
    for( d = dim ; d < 3 ; d++ )
      fprintf(fm,"%lf ",0.0);
    fprintf(fm,"\n");
  }
  VecRestoreArray( xlocal , &xvalues );
  fprintf(fm,"</DataArray>\n");

  /* <residual> */
  VecGhostUpdateBegin( b , INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(   b , INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(b , &xlocal);

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
  VecGetArray(xlocal, &xvalues);
  for( n = 0 ; n < (nl + ngho) ; n++ ){
    for( d = 0 ; d < dim ; d++ )
      fprintf(fm, "%lf ", xvalues[ n * dim + d ]);
    for( d = dim ; d < 3 ; d++ )
      fprintf(fm, "%lf ", 0.0);
    fprintf(fm,"\n");
  }
  VecRestoreArray( xlocal , &xvalues );
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"</PointData>\n");
  fprintf(fm,"<CellData>\n");

  /* <part> */
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"part\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for( e = 0; e < nelm ; e++ )
    fprintf(fm,"%d ",rank);  
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  int v;

  /* <strain> */
  fprintf(fm,"<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for( e = 0; e < nelm ; e++ ){
    for( v = 0 ; v < nvoi ; v++ )
      fprintf(fm, "%lf ", strain[ e*nvoi + v ]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  /* <stress> */
  fprintf(fm,"<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for( e = 0; e < nelm ; e++ ){
    for( v = 0 ; v < nvoi ; v++ )
      fprintf(fm, "%lf ", stress[ e*nvoi + v ]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  /* <elm_id> */
//  fprintf(fm,
//      "<DataArray type=\"Int32\" Name=\"elm_id\" NumberOfComponents=\"1\" >\n");
//  for (i=0;i<nelm;i++){
//    fprintf(fm, "%d ", elm_id[i]);
//  }
//  fprintf(fm,"\n");
//  fprintf(fm,"</DataArray>\n");

  /* <energy> */
//  fprintf(fm,"<DataArray type=\"Float64\" Name=\"energy\" NumberOfComponents=\"1\" format=\"ascii\">\n");
//  for (i=0;i<nelm;i++){
//    fprintf(fm, "%lf ", energy[i]);
//  }
//  fprintf(fm,"\n");
//  fprintf(fm,"</DataArray>\n");

  fprintf(fm,
      "</CellData>\n"
      "</Piece>\n"
      "</UnstructuredGrid>\n"
      "</VTKFile>\n");

  fclose(fm);
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
void micro_print_info( void ){

  FILE *fm = fopen("mic_info.dat","w");

  int    *i_data;
  int    i , ierr;
  
  if( rank_mic == 0 ){
    fprintf(fm,"-----------\n");
    fprintf(fm,"nproc %d\n", nproc_mic);
    fprintf(fm,"-----------\n");
    fprintf(fm,"nx    ny    nz\n");
    fprintf(fm,"%2d    %2d    %2d\n", nx, ny, nz);
    fprintf(fm,"lx    ly    lz\n");
    fprintf(fm,"%lf    %lf    %lf\n", lx, ly, lz);
    fprintf(fm,"hx    hy    hz\n");
    fprintf(fm,"%lf    %lf    %lf\n", hx, hy, hz);
    fprintf(fm,"nex   ney   nez\n");
    fprintf(fm,"%2d    %2d    %2d\n", nex, ny-1, nez);
    fprintf(fm,"-----------\n");
  }

  if( rank_mic == 0 ){
    i_data = malloc(nproc_mic*sizeof(double));
    fprintf(fm,"%-20s","rank ");
    for( i = 0 ; i < nproc_mic ; i++ )
      fprintf(fm,"%d ", i);
    fprintf(fm,"\n");
    fprintf(fm,"%-20s","");
    for( i = 0 ; i < nproc_mic ; i++ )
      fprintf(fm,"--");
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&ney, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return;

  if( rank_mic == 0 ){
    fprintf(fm,"%-20s","eyl");
    for( i = 0 ; i < nproc_mic ; i++ )
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&nyl, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return;

  if( rank_mic == 0 ){
    fprintf(fm,"%-20s","nyl");
    for( i = 0 ; i < nproc_mic ; i++ )
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&ny_inf, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return;

  if( rank_mic == 0 ){
    fprintf(fm,"%-20s","ny_inf");
    for( i = 0 ; i < nproc_mic ; i++ )
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&ngho, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return;

  if( rank_mic == 0 ){
    fprintf(fm,"%-20s","ngho");
    for( i = 0 ; i < nproc_mic ; i++ )
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  double norm;

  VecNorm( b, NORM_2, &norm );
  if( rank_mic == 0 )
    fprintf(fm,"|b| = %lf \n", norm);

  VecNorm( x, NORM_2, &norm );
  if( rank_mic == 0 )
    fprintf(fm,"|x| = %lf \n", norm);

  MatNorm( A , NORM_FROBENIUS , &norm );
  if( rank_mic == 0 )
    fprintf(fm,"|A| = %lf \n",norm);

  fclose(fm);
  return;
}
/****************************************************************************************************/
int init_shapes(double ***sh, double ****dsh, double **wp, double *h, int dim)
{
  int nsh = ( dim == 2 ) ? 4 : 8;
  int ngp = ( dim == 2 ) ? 4 : 8;
  int gp;
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

    double hx = h[0], hy = h[1];
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

    double vol = hx*hy;
    for( gp = 0 ; gp < ngp ; gp++ )
      (*wp)[gp] = vol / ngp;

  }

  free(xp);

  return 0;
}
