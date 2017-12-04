#include "macro.h"

static char help[] = 
"MACRO MULTISCALE CODE                                                                             \n"
"Solves the displacement field inside a solid structure.                                           \n"
"-coupl       : coupled with \"micro\" code for solving multiscale problem                         \n"
"-normal      : normal execution, solves a time dependent boundary condition problem               \n"
"-testcomm    : communication testing with the \"micro\" code                                      \n"
"-eigen       : calculates the eigensystem Mx = -(1/omega)Kx                                       \n"
"-print_petsc : prints petsc structures on files such as Mat and Vec objects                       \n"
"-print_vtu   : prints solutions on .vtu and .pvtu files                                           \n";

params_t params;

int main(int argc, char **argv)
{

#define CHECK_AND_GOTO(error){if(error){myio_printf(&MACRO_COMM, "error line %d at %s\n", __LINE__, __FILE__);goto end;}}
#define CHECK_INST_ELSE_GOTO(cond,instr){if(cond){instr}else{myio_printf(&MACRO_COMM, "error line %d at %s\n", __LINE__, __FILE__);goto end;}}

  int        ierr;
  int        i, j;

  myname = strdup("macro");

  myio_comm_line_init(argc, argv, &command_line);

  WORLD_COMM = MPI_COMM_WORLD;
  MPI_Init( &argc, &argv );
  MPI_Comm_size( WORLD_COMM, &nproc_wor );
  MPI_Comm_rank( WORLD_COMM, &rank_wor );

  flag_coupling = PETSC_FALSE;

  myio_comm_line_search_option(&command_line, "-coupl");
  if(command_line.found){
    flag_coupling = PETSC_TRUE;
  }

  macmic.type = COUP_1;
  color = MACRO;
  ierr = macmic_coloring( WORLD_COMM, &color, &macmic, &MACRO_COMM );
  if( ierr ){
    myio_printf(&PETSC_COMM_WORLD, "problem in coloring\n" );
    goto end;
  }

  MPI_Comm_size(MACRO_COMM, &nproc_mac);
  MPI_Comm_rank(MACRO_COMM, &rank_mac);

  PETSC_COMM_WORLD = MACRO_COMM;

#ifdef SLEPC
  SlepcInitialize(&argc,&argv,(char*)0,help);
#elif  PETSC
  PetscInitialize(&argc,&argv,(char*)0,help);
#endif

  myio_printf(&MACRO_COMM,
      "--------------------------------------------------\n"
      "  MACRO: COMPOSITE MATERIAL MULTISCALE CODE\n"
      "--------------------------------------------------\n");

  init_variables(&params);
  myio_comm_line_search_option(&command_line, "-normal");

  if(command_line.found){

    params.calc_mode = CALC_MODE_NORMAL;

    myio_comm_line_get_double(&command_line, "-tf");
    CHECK_INST_ELSE_GOTO(command_line.found, params.final_time = command_line.double_val;)

    myio_comm_line_get_double(&command_line, "-dt");
    CHECK_INST_ELSE_GOTO(command_line.found, params.delta_time = command_line.double_val;)
  }

  myio_comm_line_search_option(&command_line, "-testcomm");
  if(command_line.found){
    params.calc_mode = CALC_MODE_TEST;
  }

  myio_comm_line_search_option(&command_line, "-eigen");
  if(command_line.found){
    params.calc_mode = CALC_MODE_EIGEN;
    myio_comm_line_get_double(&command_line, "-energy_stored");
    if(command_line.found) params.energy_stored = command_line.double_val;
  }

  mesh_f = FORMAT_GMSH;

  myio_comm_line_get_string(&command_line, "-mesh");
  if(command_line.found){
    mesh_n = strdup(command_line.str);
  }
  else{
    myio_printf(&MACRO_COMM,"mesh file not given on command line.\n");
    goto end;
  }

  FILE *fm = fopen( mesh_n, "r");
  if( fm == NULL ){
    myio_printf(&MACRO_COMM,"mesh file not found.\n");
    goto end;
  }
  myio_comm_line_get_int(&command_line, "-dim");
  CHECK_INST_ELSE_GOTO(command_line.found, dim = command_line.int_val;)

  nvoi    = (dim == 2) ? 3 : 6;
  npe_max = (dim == 2) ? 4 : 8;
  ngp_max = npe_max;

  flag_print = 0;
  myio_comm_line_search_option(&command_line, "-print_petsc");
  if(command_line.found) flag_print = flag_print | (1<<PRINT_PETSC);

  myio_comm_line_search_option(&command_line, "-print_vtu");
  if(command_line.found) flag_print = flag_print | (1<<PRINT_VTU);

  myio_comm_line_get_int(&command_line, "-nl_max_its");
  if(command_line.found) params.non_linear_max_its = command_line.int_val;

  myio_comm_line_get_int(&command_line, "-nl_min_norm_tol");
  if(command_line.found) params.non_linear_min_norm_tol = command_line.double_val;

  ierr = function_fill_list_from_command_line(&command_line, &function_list);
  CHECK_AND_GOTO(ierr);

  ierr = mesh_fill_boundary_list_from_command_line(&command_line, &boundary_list);
  CHECK_AND_GOTO(ierr);

  ierr = material_fill_list_from_command_line(&command_line, &material_list);
  CHECK_AND_GOTO(ierr);

  partition_algorithm = PARMETIS_GEOM;
  myio_comm_line_search_option(&command_line, "-part_kway");
  if(command_line.found) partition_algorithm = PARMETIS_MESHKWAY;

  myio_comm_line_search_option(&command_line, "-part_geom");
  if(command_line.found) partition_algorithm = PARMETIS_GEOM;

  myio_printf(&MACRO_COMM, "reading mesh elements\n" );
  ierr = read_mesh_elmv(MACRO_COMM, myname, mesh_n, mesh_f);
  CHECK_AND_GOTO(ierr);

  myio_printf(&MACRO_COMM, "partitioning and distributing mesh\n");
  ierr = part_mesh(MACRO_COMM, myname, NULL);
  CHECK_AND_GOTO(ierr);

  myio_printf(&MACRO_COMM, "calculating ghost nodes\n");
  ierr = calc_local_and_ghost(MACRO_COMM, nallnods, allnods, &ntotnod, &nmynods, &mynods, &nghost, &ghost );
  CHECK_AND_GOTO(ierr);

  myio_printf(&MACRO_COMM, "reenumering nodes\n");
  ierr = reenumerate_PETSc( MACRO_COMM );
  CHECK_AND_GOTO(ierr);

  myio_printf(&MACRO_COMM, "reading Coordinates\n");
  ierr = read_coord(mesh_n, nmynods, mynods, nghost , ghost, &coord );
  CHECK_AND_GOTO(ierr);

  ierr = read_bc();
  CHECK_AND_GOTO(ierr);

  list_init(&physical_list, sizeof(physical_t), NULL );
  gmsh_get_physical_list(mesh_n, &physical_list);

  myio_printf(&MACRO_COMM, "allocating ");

  A   = NULL;
  b   = NULL;
  x   = NULL;
  dx  = NULL;

  int ixpe = npe_max * dim;
  loc_elem_index = malloc( ixpe * sizeof(int));
  glo_elem_index = malloc( ixpe * sizeof(int));
  elem_disp      = malloc( ixpe * sizeof(double));
  elem_coor      = malloc( ixpe * sizeof(double));
  k_elem         = malloc( ixpe*ixpe * sizeof(double));
  m_elem         = malloc( ixpe*ixpe * sizeof(double));
  res_elem       = malloc( ixpe * sizeof(double));
  stress_gp      = malloc( nvoi * sizeof(double));
  strain_gp      = malloc( nvoi * sizeof(double));
  c              = malloc( nvoi*nvoi * sizeof(double));
  if( flag_print & ( 1 << PRINT_VTU ) ){
    elem_strain  = malloc( nelm*nvoi * sizeof(double));
    elem_stress  = malloc( nelm*nvoi * sizeof(double));
    elem_energy  = malloc( nelm      * sizeof(double));
    elem_type    = malloc( nelm      * sizeof(int));
  }
  flag_neg_detj  = 0;

  bmat = malloc( nvoi * sizeof(double**));
  for( i = 0 ; i < nvoi  ; i++ ){
    bmat[i] = malloc( ixpe * sizeof(double*));
    for( j = 0 ; j < ixpe ; j++ )
      bmat[i][j] = malloc( ngp_max * sizeof(double));
  }

  dsh  = malloc( npe_max * sizeof(double**));
  for( i = 0 ; i < npe_max ; i++ ){
    dsh[i] = malloc( dim * sizeof(double*));
    for( j = 0 ; j < dim ; j++ )
      dsh[i][j] = malloc( ngp_max * sizeof(double));
  }

  jac = malloc( dim * sizeof(double*));
  for( int k = 0 ; k < dim ; k++ )
    jac[k] = malloc( dim * sizeof(double));

  jac_inv = malloc( dim * sizeof(double*));
  for( i = 0 ; i < dim ; i++ )
    jac_inv[i] = malloc( dim * sizeof(double));

  detj = malloc( ngp_max * sizeof(double));

  myio_printf(&MACRO_COMM, "ok\n");

  ierr = fem_inigau();

  int      nr_its = -1;
  double   norm = -1.0;
  double   limit[6];

  ierr = get_bbox_local_limits(coord, nallnods, &limit[0], &limit[2], &limit[4]);
  myio_printf(&MACRO_COMM,"Limit = ");
  for( i = 0 ; i < nvoi ; i++ )
    myio_printf(&MACRO_COMM,"%lf ", limit[i]);
  myio_printf(&MACRO_COMM,"\n");

  if(params.calc_mode == CALC_MODE_EIGEN){

    EPS     eps;
    int     nnz = ( dim == 2 ) ? dim*9 : dim*27;

    MatCreate( MACRO_COMM, &A );
    MatSetSizes( A, dim*nmynods, dim*nmynods, dim*ntotnod, dim*ntotnod );
    MatSetType( A, MATAIJ );
    MatSeqAIJSetPreallocation( A, nnz, NULL );
    MatMPIAIJSetPreallocation( A, nnz, NULL, nnz, NULL );
    MatSetOption( A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );
    MatSetUp( A );
    MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd  ( A, MAT_FINAL_ASSEMBLY );

    MatCreate( MACRO_COMM, &M );
    MatSetSizes( M, dim*nmynods, dim*nmynods, dim*ntotnod, dim*ntotnod );
    MatSetType( M, MATAIJ );
    MatSeqAIJSetPreallocation( M, nnz, NULL );
    MatMPIAIJSetPreallocation( M, nnz, NULL, nnz, NULL );
    MatSetOption( M, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );
    MatSetUp( M );
    MatAssemblyBegin( M, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd  ( M, MAT_FINAL_ASSEMBLY );

    int *ghost_index = malloc( nghost*dim *sizeof(int) );

    int d;
    for( i = 0 ; i < nghost ; i++ ){
      for( d = 0 ; d < dim ; d++ )
	ghost_index[ i * dim + d ] = loc2petsc[ nmynods + i ] * dim + d;
    }

    VecCreateGhost( MACRO_COMM, dim*nmynods, dim*ntotnod, nghost*dim, ghost_index, &x );
    VecZeroEntries( x );

    VecGhostUpdateBegin( x , INSERT_VALUES , SCATTER_FORWARD );
    VecGhostUpdateEnd  ( x , INSERT_VALUES , SCATTER_FORWARD );

    ierr = assembly_AM();
    if( ierr ){
      myio_printf(&MACRO_COMM,"problem during matrix assembly\n");
      goto end;
    }

    node_list_t *pn = boundary_list.head;
    while( pn )
    {
      mesh_boundary_t * bou = (mesh_boundary_t * )pn->data;
      MatZeroRowsColumns( M, bou->ndirix, bou->dir_glo_ixs, 1.0, NULL, NULL );
      MatZeroRowsColumns( A, bou->ndirix, bou->dir_glo_ixs, 1.0, NULL, NULL );
      pn = pn->next;
    }
    MatAssemblyBegin( M, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd  ( M, MAT_FINAL_ASSEMBLY );
    MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd  ( A, MAT_FINAL_ASSEMBLY );

    if(flag_print & ( 1<<PRINT_PETSC )){
      PetscViewer  viewer;
      PetscViewerASCIIOpen(MACRO_COMM,"M.dat" ,&viewer); MatView(M ,viewer);
      PetscViewerASCIIOpen(MACRO_COMM,"A.dat" ,&viewer); MatView(A ,viewer);
    }

    int    nconv;
    double error;

    EPSCreate(MACRO_COMM, &eps);
    EPSSetOperators(eps, M, A);
    EPSSetProblemType(eps, EPS_GHEP);
    EPSSetFromOptions(eps);
    EPSGetDimensions(eps, &params.num_eigen_vals, NULL, NULL);
    params.eigen_vals = malloc( params.num_eigen_vals*sizeof(double));
    myio_printf(&MACRO_COMM,"Number of requested eigenvalues: %d\n", params.num_eigen_vals);

    EPSSolve(eps);
    EPSGetConverged(eps, &nconv);
    myio_printf(&MACRO_COMM,"Number of converged eigenpairs: %d\n",nconv);

    for( i=0 ; i<params.num_eigen_vals ; i++ ){

      EPSGetEigenpair( eps, i, &params.eigen_vals[i], NULL, x, NULL );
      EPSComputeError( eps, i, EPS_ERROR_RELATIVE, &error );
      myio_printf(&MACRO_COMM, "omega %d = %e   error = %e\n", i, params.eigen_vals[i], error);

      if(flag_print & (1<<PRINT_VTU))
      { 
	get_elem_properties();
	sprintf( filename, "macro_eigen_%d", i);
	macro_pvtu( filename );
      }

    }

    EPSDestroy(&eps);

  }
  else if(params.calc_mode == CALC_MODE_NORMAL){

    int     nnz = ( dim == 2 ) ? dim*9 : dim*27;
    KSP     ksp;

    MatCreate( MACRO_COMM, &A );
    MatSetSizes( A, dim*nmynods, dim*nmynods, dim*ntotnod, dim*ntotnod );
    MatSetType( A, MATAIJ );
    MatSeqAIJSetPreallocation( A, nnz, NULL );
    MatMPIAIJSetPreallocation( A, nnz, NULL, nnz, NULL );
    MatSetOption( A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE );
    MatSetUp( A );
    MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd  ( A, MAT_FINAL_ASSEMBLY );

    int * ghost_index = malloc( nghost*dim * sizeof(int) );

    int d;
    for( i = 0 ; i < nghost ; i++ ){
      for( d = 0 ; d < dim ; d++ )
	ghost_index[ i * dim + d ] = loc2petsc[ nmynods + i ] * dim + d;
    }

    VecCreateGhost( MACRO_COMM, dim*nmynods, dim*ntotnod, nghost*dim, ghost_index, &x );
    VecDuplicate( x, &dx );
    VecDuplicate( x, &b );

    KSPCreate( MACRO_COMM, &ksp );
    KSPSetFromOptions( ksp );

    params.time = 0.0;
    params.time_step = 0;

    VecZeroEntries( x );
    VecGhostUpdateBegin( x, INSERT_VALUES, SCATTER_FORWARD );
    VecGhostUpdateEnd  ( x, INSERT_VALUES, SCATTER_FORWARD );

    while( params.time < (params.final_time + 1.0e-10) ){

      myio_printf(&MACRO_COMM,"\ntime step %-3d %-e seg\n", params.time_step, params.time);

      update_boundary( params.time , &function_list, &boundary_list );

      Vec      x_loc,  b_loc;
      double  *x_arr, *b_arr;

      VecGhostGetLocalForm( x    , &x_loc );
      VecGetArray(          x_loc, &x_arr );

      node_list_t * pn = boundary_list.head;
      while( pn )
      {
	mesh_boundary_t *bou = ( mesh_boundary_t * )pn->data;
	for( i = 0 ; i < bou->ndirix ; i++ )
	  x_arr[bou->dir_loc_ixs[i]] = bou->dir_val[i];
	pn = pn->next;
      }

      VecRestoreArray         ( x_loc, &x_arr );
      VecGhostRestoreLocalForm( x    , &x_loc );

      VecGhostUpdateBegin( x, INSERT_VALUES, SCATTER_FORWARD );
      VecGhostUpdateEnd  ( x, INSERT_VALUES, SCATTER_FORWARD );

      nr_its = 0; norm = 2 * params.non_linear_min_norm_tol;
      while( nr_its < params.non_linear_max_its && norm > params.non_linear_min_norm_tol )
      {

	myio_printf(&MACRO_COMM, "MACRO: assembling residual\n" );
	assembly_b();

	if( flag_neg_detj == 1)
	  myio_printf(&MACRO_COMM, "MACRO: warning negative jacobian detected\n");

	VecGhostGetLocalForm( b    , &b_loc );
	VecGetArray         ( b_loc, &b_arr );
	pn = boundary_list.head;
	while( pn )
	{
	  mesh_boundary_t *bou = ( mesh_boundary_t * )pn->data;
	  for( i = 0 ; i < bou->ndirix ; i++ )
	    b_arr[bou->dir_loc_ixs[i]] = 0.0;
	  pn = pn->next;
	}
	VecRestoreArray         ( b_loc, &b_arr );
	VecGhostRestoreLocalForm( b    , &b_loc );
	VecGhostUpdateBegin( b, INSERT_VALUES, SCATTER_FORWARD );
	VecGhostUpdateEnd  ( b, INSERT_VALUES, SCATTER_FORWARD );

	VecNorm( b, NORM_2, &norm );
	VecScale( b, -1.0 );

	myio_printf(&MACRO_COMM,"MACRO: |b| = %e\n", norm );

	if( norm < params.non_linear_min_norm_tol ) break;

	myio_printf(&MACRO_COMM, "MACRO: assembling jacobian\n");
	assembly_A();

	node_list_t *pn = boundary_list.head;
	while( pn )
	{
	  mesh_boundary_t * bou = (mesh_boundary_t * )pn->data;
	  MatZeroRowsColumns( A, bou->ndirix, bou->dir_glo_ixs, 1.0, NULL, NULL );
	  pn = pn->next;
	}
	MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd  ( A, MAT_FINAL_ASSEMBLY );

	myio_printf(&MACRO_COMM, "MACRO: solving system\n" );
	KSPSetOperators( ksp, A, A );
	KSPSolve( ksp, b, dx );
	print_ksp_info( MACRO_COMM, ksp);
	myio_printf(&MACRO_COMM, "\n");

	VecAXPY( x, 1.0, dx );
	VecGhostUpdateBegin( x, INSERT_VALUES, SCATTER_FORWARD );
	VecGhostUpdateEnd  ( x, INSERT_VALUES, SCATTER_FORWARD );

	nr_its ++;
      }

      if(flag_print & (1<<PRINT_PETSC)){
	PetscViewer  viewer;
	PetscViewerASCIIOpen(MACRO_COMM,"A.dat" ,&viewer); MatView(A ,viewer);
	PetscViewerASCIIOpen(MACRO_COMM,"b.dat" ,&viewer); VecView(b ,viewer);
	PetscViewerASCIIOpen(MACRO_COMM,"dx.dat",&viewer); VecView(dx,viewer);
	PetscViewerASCIIOpen(MACRO_COMM,"x.dat" ,&viewer); VecView(x ,viewer);
      }

      if( flag_print & (1<<PRINT_VTU) )
      { 
	get_elem_properties();
	sprintf( filename, "macro_t_%d", params.time_step );
	macro_pvtu( filename );
      }

      params.time += params.delta_time;
      params.time_step ++;
    }
    KSPDestroy(&ksp);
  }
  else if(params.calc_mode == CALC_MODE_TEST){

    double   strain_mac[6] = {0.1, 0.1, 0.2, 0.0, 0.0, 0.0};
    double   stress_mac[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for( i = 0 ; i < nvoi ; i++ ){
      for( j = 0 ; j < nvoi ; j++ )
	strain_mac[j] = 0.0;
      strain_mac[i] = 0.005;
      ierr = mac_send_signal(WORLD_COMM, MAC2MIC_STRAIN);
      ierr = mac_send_strain(WORLD_COMM, strain_mac    );
      ierr = mac_recv_stress(WORLD_COMM, stress_mac    );
      myio_printf(&MACRO_COMM,"\nstress_ave = ");
      for( j = 0 ; j < nvoi ; j++ )
	myio_printf(&MACRO_COMM,"%e ",stress_mac[j]);
      myio_printf(&MACRO_COMM,"\n");
    }

  }

  free(loc_elem_index); 
  free(glo_elem_index); 
  free(elem_disp     ); 
  free(stress_gp     ); 
  free(strain_gp     ); 
  if( flag_print & ( 1 << PRINT_VTU ) ){
    free(elem_strain); 
    free(elem_stress); 
    free(elem_energy); 
    free(elem_type); 
  }

  for( i = 0 ; i < nvoi  ; i++ ){
    for( j = 0 ; j < ixpe ; j++ )
      free(bmat[i][j]);
    free(bmat[i]);
  }
  free(bmat);

  for( i = 0 ; i < npe_max ; i++ ){
    for( j = 0 ; j < dim ; j++ )
      free(dsh[i][j]);
    free(dsh[i]);
  }
  free(dsh);

end:

  if(flag_coupling){
    ierr = mac_send_signal(WORLD_COMM, MIC_END);
    if(ierr){
      myio_printf(&PETSC_COMM_WORLD, "macro: problem sending MIC_END to micro\n");
      return 1;
    }
  }

  list_clear(&material_list);
  list_clear(&physical_list);
  list_clear(&function_list);

  MatDestroy(&A);
  VecDestroy(&x);
  VecDestroy(&b);

  myio_printf(&MACRO_COMM,
      "--------------------------------------------------\n"
      "  MACRO: FINISH COMPLETE\n"
      "--------------------------------------------------\n");

#ifdef SLEPC
  SlepcFinalize();
#elif  PETSC
  PetscFinalize();
#endif

  ierr = MPI_Finalize();

  return 0;
}


int read_bc()
{

  int        *ix, i, d, da, n;
  int         ierr;
  mesh_boundary_t    *bou;

  node_list_t *pn = boundary_list.head;
  while( pn )
  {
    bou = ( mesh_boundary_t * )pn->data;
    ierr = gmsh_get_node_index( mesh_n, bou->name, nmynods, mynods, dim, &n, &ix );
    if( ierr ){
      myio_printf(&MACRO_COMM, "problem finding nodes of boundary %s on msh file\n", bou->name );
      return 1;
    }
    bou->ndir        = n;
    bou->ndirix      = bou->ndir * bou->ndirpn;
    bou->dir_val     = malloc( bou->ndirix * sizeof(double));
    bou->dir_loc_ixs = malloc( bou->ndirix * sizeof(int));
    bou->dir_glo_ixs = malloc( bou->ndirix * sizeof(int));
    for( i = 0 ; i < n ; i++ ){
      da = 0;
      int * p = bsearch( &ix[i], mynods, nmynods, sizeof(int), cmpfunc );
      for( d = 0 ; d < dim ; d++ )
	if( bou->kind & (1<<d) ) {
	  bou->dir_loc_ixs[i* (bou->ndirpn) + da] = (p - mynods) * dim + d;
	  bou->dir_glo_ixs[i* (bou->ndirpn) + da] = loc2petsc[(p - mynods)] * dim + d;
	  da++;
	}
    }
    free(ix);
    pn = pn->next;
  }

  return 0;
}


int read_coord( char *mesh_n, int nmynods, int *mynods, int nghost , int *ghost, double **coord )
{

  (*coord) = malloc( ( nmynods + nghost )*dim * sizeof(double));

  int ierr = gmsh_read_coord_parall( mesh_n, dim, nmynods, mynods, nghost , ghost, *coord );

  return ierr;
}


int assembly_b( void )
{
  int      npe, ngp;
  int      e, gp, i, j;
  double  *wp;

  double  *b_arr;
  Vec      b_loc;

  VecGhostGetLocalForm( b    , &b_loc );
  VecGetArray         ( b_loc, &b_arr );

  for( i = 0 ; i < nallnods*dim ; i++ )
    b_arr[i] = 0.0;

  for( e = 0 ; e < nelm ; e++ )
  {
    get_local_elem_index( e, loc_elem_index );
    ngp = npe = eptr[e+1] - eptr[e];

    for( i = 0 ; i < npe*dim ; i++ )
      res_elem[i] = 0.0;

    get_dsh( e, loc_elem_index, dsh, detj);
    get_bmat( e, dsh, bmat );
    get_wp( dim, npe, &wp );

    for( gp = 0; gp < ngp ; gp++ ){

      if( detj[gp] < 0.0 ) flag_neg_detj = 1;
      detj[gp] = fabs( detj[gp] );

      /* calc strain at gp */
      get_strain( e , gp, loc_elem_index, dsh, bmat, strain_gp );

      /* calc stress at gp */
      get_stress( e , gp, strain_gp, stress_gp );

      for( i = 0 ; i < npe*dim ; i++ ){
	for( j = 0; j < nvoi ; j++ )
	  res_elem[i] += bmat[j][i][gp] * stress_gp[j] * wp[gp] * detj[gp];
      }
    }

#ifdef ZERO
    for( i = 0 ; i < (npe * dim) ; i++ ) 
      res_elem[i] = ( fabs(res_elem[i]) < 1.0e-6 ) ? 0.0 : res_elem[i];
#endif

    for( i = 0 ; i < ( npe * dim ) ; i++ )
      b_arr[ loc_elem_index[i] ] += res_elem[i];

  }

  VecRestoreArray         ( b_loc , &b_arr );
  VecGhostRestoreLocalForm( b     , &b_loc );

  VecGhostUpdateBegin( b, ADD_VALUES   , SCATTER_REVERSE );
  VecGhostUpdateEnd  ( b, ADD_VALUES   , SCATTER_REVERSE );
  VecGhostUpdateBegin( b, INSERT_VALUES, SCATTER_FORWARD );
  VecGhostUpdateEnd  ( b, INSERT_VALUES, SCATTER_FORWARD );

  return 0;
}


int assembly_AM( void )
{

  MatZeroEntries(A);
  MatZeroEntries(M);

  int       e, gp;
  int       i, j, d, k, h;
  int       npe, ngp;
  int       ierr;
  double    rho_gp;
  double  **sh;
  double   *wp;

  for( e = 0 ; e < nelm ; e++ )
  {
    ngp = npe = eptr[e+1] - eptr[e];

    for( i = 0 ; i < npe*dim*npe*dim ; i++ ){
      m_elem[i] = 0.0;
      k_elem[i] = 0.0;
    }

    get_local_elem_index ( e, loc_elem_index );
    get_global_elem_index( e, glo_elem_index );

    get_sh( dim, npe, &sh );
    get_dsh( e, loc_elem_index, dsh, detj);
    get_bmat( e, dsh, bmat );
    get_wp( dim, npe, &wp );

    for( gp = 0; gp < ngp ; gp++ ){

      detj[gp] = fabs( detj[gp] );

      get_strain( e , gp, loc_elem_index, dsh, bmat, strain_gp );

      ierr = get_c_tan( NULL , e , gp , strain_gp , c ); if( ierr ) return 1;
      ierr = get_rho  ( NULL , e , &rho_gp );            if( ierr ) return 1;

      for( i = 0 ; i < npe*dim ; i++ ){
	for( j = 0 ; j < npe*dim ; j++ ){
	  for( k = 0; k < nvoi ; k++ ){
	    for( h = 0; h < nvoi ; h++ )
	      k_elem[ i*npe*dim + j] += \
	      bmat[h][i][gp] * c[ h*nvoi + k ] * bmat[k][j][gp] * wp[gp] * detj[gp] ;
	  }
	}
      }

      for( d = 0 ; d < dim ; d++ ){
	for( i = 0 ; i < npe; i++ ){
	  for( j = 0 ; j < npe; j++ )
	    m_elem[ (i*dim)*(npe*dim) + j*dim + (d*dim*npe + d)] += \
	    rho_gp * sh[i][gp] * sh[j][gp] * wp[gp] * detj[gp] ;
	}
      }

    }
    MatSetValues( A, npe*dim, glo_elem_index, npe*dim, glo_elem_index, k_elem, ADD_VALUES );
    MatSetValues( M, npe*dim, glo_elem_index, npe*dim, glo_elem_index, m_elem, ADD_VALUES );

  }

  MatAssemblyBegin( A , MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd  ( A , MAT_FINAL_ASSEMBLY );
  MatAssemblyBegin( M , MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd  ( M , MAT_FINAL_ASSEMBLY );

  return 0;
}


int assembly_A( void )
{

  MatZeroEntries(A);

  int      e, gp;
  int      i, j, k, h;
  int      ierr;
  double  *wp;

  for( e = 0 ; e < nelm ; e++ )
  {
    int npe = eptr[e+1] - eptr[e];
    int ngp = npe;

    for( i = 0 ; i < npe*dim*npe*dim ; i++ )
      k_elem[i] = 0.0;
    
    get_local_elem_index (e, loc_elem_index);
    get_global_elem_index(e, glo_elem_index);

    get_dsh( e, loc_elem_index, dsh, detj);
    get_bmat( e, dsh, bmat );
    get_wp( dim, npe, &wp );

    for( gp = 0; gp < ngp ; gp++ ){

      detj[gp] = fabs( detj[gp] );

      get_strain( e , gp, loc_elem_index, dsh, bmat, strain_gp );

      ierr = get_c_tan( NULL , e , gp , strain_gp , c ); if( ierr ) return 1;

      for( i = 0 ; i < npe*dim ; i++ ){
	for( j = 0 ; j < npe*dim ; j++ ){
	  for( k = 0; k < nvoi ; k++ ){
	    for( h = 0; h < nvoi ; h++ )
	      k_elem[ i*npe*dim + j] += \
	      bmat[h][i][gp] * c[ h*nvoi + k ] * bmat[k][j][gp] * wp[gp] * detj[gp] ;
	  }
	}
      }

    }
    MatSetValues( A, npe*dim, glo_elem_index, npe*dim, glo_elem_index, k_elem, ADD_VALUES );

  }

  MatAssemblyBegin( A , MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd  ( A , MAT_FINAL_ASSEMBLY );

  return 0;
}


int get_strain( int e , int gp, int *loc_elem_index, double ***dsh_gp,  double ***bmat, double *strain_gp )
{

  double  *x_arr; 
  Vec      x_loc; 
  VecGhostGetLocalForm( x    , &x_loc );
  VecGetArray         ( x_loc, &x_arr );

  int  i , v;
  int  npe = eptr[e+1] - eptr[e];
  for( i = 0 ; i < ( npe * dim ) ; i++ )
    elem_disp[i] = x_arr[ loc_elem_index[i] ];

  VecRestoreArray         ( x_loc , &x_arr);
  VecGhostRestoreLocalForm( x     , &x_loc);

  for( v = 0; v < nvoi ; v++ ){
    strain_gp[v] = 0.0;
    for( i = 0 ; i < npe*dim ; i++ )
      strain_gp[v] += bmat[v][i][gp] * elem_disp[i];
    strain_gp[v] = ( fabs(strain_gp[v]) < 1.0e-6 ) ? 0.0 : strain_gp[v];
  }

  return 0;
}


int get_stress( int e , int gp, double *strain_gp , double *stress_gp )
{

  char        name_s[64];
  int         ierr;
  int         macro_gp = 0;
  material_t  *mat_p;
  get_mat_name( elm_id[e], name_s );

  node_list_t *pn = material_list.head;
  while( pn ){
    mat_p = ( material_t * )pn->data;
    if( strcmp( name_s, mat_p->name ) == 0 ) break;
    pn = pn->next;
  }
  if( pn == NULL ){
    myio_printf(&MACRO_COMM, "Material %s corresponding to element %d not found on material list\n", name_s, e );
    return 1;
  }

  if( mat_p->type_id == MAT_MICRO )
  {
    ierr = mac_send_signal(WORLD_COMM, MAC2MIC_STRAIN); if(ierr) return 1;
    ierr = mac_send_strain(WORLD_COMM, strain_gp);
    ierr = mac_send_macro_gp(WORLD_COMM, &macro_gp);
    ierr = mac_recv_stress(WORLD_COMM, stress_gp);
  }
  else
    material_get_stress( mat_p, dim, strain_gp, stress_gp );

#ifdef ZERO
  int v;
  for( v = 0 ; v < nvoi ; v++ )
    stress_gp[v] = ( fabs(stress_gp[v]) < 1.0e-6 ) ? 0.0 : stress_gp[v];
#endif

  return 0;
}


int get_c_tan( const char * name, int e , int gp , double * strain_gp , double * c_tan )
{

  char        name_s[64];
  int         ierr;
  int         macro_gp = 0;
  material_t  *mat_p;
  get_mat_name( elm_id[e], name_s );

  node_list_t *pn = material_list.head;
  while( pn ){
    mat_p = ( material_t * )pn->data;
    if( strcmp( name_s, mat_p->name ) == 0 ) break;
    pn = pn->next;
  }
  if( pn == NULL ){
    myio_printf(&MACRO_COMM, "Material %s corresponding to element %d not found on material list\n", name_s, e );
    return 1;
  }

  if( mat_p->type_id == MAT_MICRO )
  {
    ierr = mac_send_signal(WORLD_COMM, C_HOMO); if(ierr) return 1;
    ierr = mac_send_strain(WORLD_COMM, strain_gp);
    ierr = mac_send_macro_gp(WORLD_COMM, &macro_gp);
    ierr = mac_recv_c_homo(WORLD_COMM, nvoi, c_tan);
  }
  else
    material_get_c_tang( mat_p, dim, strain_gp, c_tan );

  return 0;
}


int get_rho( const char * name, int e , double * rho )
{

  material_t  *mat_p;
  char         name_s[64];
  int          ierr;

  get_mat_name( elm_id[e], name_s );

  node_list_t *pn = material_list.head;
  while( pn ){
    mat_p = ( material_t * )pn->data;
    if( strcmp( name_s, mat_p->name ) == 0 ) break;
    pn = pn->next;
  }
  if( pn == NULL ){
    myio_printf(&MACRO_COMM, "Material %s corresponding to element %d not found on material list\n", name_s, e );
    return 1;
  }

  if( mat_p->type_id == MAT_MICRO )
  {
    ierr = mac_send_signal(WORLD_COMM, RHO); if(ierr) return 1;
    ierr = mac_recv_rho(WORLD_COMM, rho);
  }
  else
    material_get_rho( mat_p, dim, rho );

  return 0;
}


int get_mat_name( int id , char * name_s )
{

  node_list_t *pn;
  physical_t  *phy_p;
  
  pn = physical_list.head;
  while( pn ){
    phy_p = ( physical_t * )pn->data;
    if( id == phy_p->id ) break;
    pn = pn->next;
  }
  if( pn == NULL ) return 1;

  strcpy( name_s, phy_p->name );

  return 0;
}


int get_global_elem_index( int e, int * glo_elem_index )
{

  int  n, d;
  int  npe = eptr[e+1] - eptr[e];

  for( n = 0 ; n < npe ; n++ ){
    for( d = 0 ; d < dim ; d++ )
      glo_elem_index[ n * dim + d ] = loc2petsc[ eind[ eptr[e] + n ] ] * dim + d;
  }
  return 0;
}


int get_local_elem_index( int e, int * loc_elem_index )
{

  int  n, d;
  int  npe = eptr[e+1] - eptr[e];

  for( n = 0 ; n < npe ; n++ ){
    for( d = 0 ; d < dim ; d++ )
      loc_elem_index[ n * dim + d ] = eind[ eptr[e] + n ] * dim + d;
  }
  return 0;
}


int get_dsh( int e, int *loc_elem_index, double ***dsh, double *detj )
{

  double ***dsh_master;
  int       i, gp;
  int       npe = eptr[e+1] - eptr[e];
  int       ngp = npe;

  for( i = 0 ; i < npe*dim ; i++ )
    elem_coor[i] = coord[loc_elem_index[i]];

  for( gp = 0; gp < ngp ; gp++ ){

    fem_get_dsh_master( npe, dim, &dsh_master );

    fem_calc_jac( dim, npe, gp, elem_coor, dsh_master, jac );
    fem_invjac( dim, jac, jac_inv, &detj[gp] );
    fem_trans_dsh( dim, npe, gp, jac_inv, dsh_master, dsh );

  }

  return 0;
}


int get_bmat( int e, double ***dsh, double ***bmat )
{


  int       i, gp;
  int       npe = eptr[e+1] - eptr[e];
  int       ngp = npe;

  if( dim == 2 ){
    for( i = 0 ; i < npe ; i++ ){
      for( gp = 0; gp < ngp ; gp++ ){
	bmat[0][i*dim + 0][gp] = dsh[i][0][gp];
	bmat[0][i*dim + 1][gp] = 0             ;
	bmat[1][i*dim + 0][gp] = 0             ;
	bmat[1][i*dim + 1][gp] = dsh[i][1][gp];
	bmat[2][i*dim + 0][gp] = dsh[i][1][gp];
	bmat[2][i*dim + 1][gp] = dsh[i][0][gp];
      }
    }
  }

  return 0;
}


int get_sh( int dim, int npe, double ***sh )
{

  fem_get_sh( npe, dim, sh );

  return 0;
}


int get_wp( int dim, int npe, double **wp )
{

  fem_get_wp( npe, dim, wp );

  return 0;
}


int get_elem_properties( void )
{

  int      e, v, gp;
  double  *strain_aux = malloc( nvoi * sizeof(double) );
  double  *stress_aux = malloc( nvoi * sizeof(double) );
  double  *wp;

  for ( e = 0 ; e < nelm ; e++ ){

    int     npe = eptr[e+1] - eptr[e];
    int     ngp = npe;
    double  vol_elem = 0.0;

    for ( v = 0 ; v < nvoi ; v++ ) 
      strain_aux[v] = stress_aux[v] = 0.0;

    get_local_elem_index (e, loc_elem_index);

    get_dsh( e, loc_elem_index, dsh, detj);
    get_bmat( e, dsh, bmat );
    get_wp( dim, npe, &wp );

    for ( gp = 0 ; gp < ngp ; gp++ ){

      detj[gp] = fabs( detj[gp] );

      get_strain( e , gp, loc_elem_index, dsh, bmat, strain_gp );
      get_stress( e , gp, strain_gp, stress_gp );
      for ( v = 0 ; v < nvoi ; v++ ){
	strain_aux[v] += strain_gp[v] * detj[gp] * wp[gp];
	stress_aux[v] += stress_gp[v] * detj[gp] * wp[gp];
      }
      vol_elem += detj[gp] * wp[gp];
    }
    for ( v = 0 ; v < nvoi ; v++ ){
      elem_strain[ e*nvoi + v ] = strain_aux[v] / vol_elem;
      elem_stress[ e*nvoi + v ] = stress_aux[v] / vol_elem;
    }

    physical_t * phy;
    node_list_t * pn = physical_list.head;
    while ( pn )
    {
      phy = pn->data;
      if( phy->id == elm_id[e] ) break;
      pn = pn->next;
    }
    if( !pn ) return 1;

    int type = 0;
    pn = material_list.head;
    while ( pn )
    {
      material_t *mat = pn->data;
      if( strcmp( phy->name , mat->name) == 0 ) break;
      pn = pn->next;
      type ++;
    }
    if( !pn ) return 1;

    elem_type[e] = type;
  }

  return 0;
}


int update_boundary( double t , list_t * function_list, list_t * boundary_list )
{

  node_list_t * pn = boundary_list->head;
  while( pn )
  {
    mesh_boundary_t * bou = ( mesh_boundary_t * ) pn->data;
    function_t   * function = NULL;
    int i, d;
    for( d = 0 ; d < dim ; d++ ){
      function_get_from_list( bou->fnum[d] , function_list , &function );
      double val;
      function_eval( t , function , &val );
      for( i = 0 ; i < bou->ndir ; i++ )
	bou->dir_val[ i* (bou->ndirpn) + d ] = val;
    }
    pn = pn->next;
  }

  return 0;
}


int macro_pvtu( char *name )
{

  FILE    *fm;
  char    file_name[NBUF];
  double  *xvalues;
  Vec     xlocal;

  if( rank_mac == 0 ){

    strcpy(file_name,name);
    strcat(file_name,".pvtu");
    fm = fopen(file_name,"w");

    fprintf(fm, "<?xml version=\"1.0\"?>\n"
	"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
	"<PUnstructuredGrid GhostLevel=\"0\">\n"
	"<PPoints>\n"
	"<PDataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\"/>\n"
	"</PPoints>\n"
	"<PCells>\n"
	"<PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>\n"
	"</PCells>\n");

    fprintf(fm, "<PPointData Vectors=\"displ\">\n");
    if( x != NULL )
    fprintf(fm, "<PDataArray type=\"Float64\" Name=\"displ\"    NumberOfComponents=\"3\" />\n");
    if( b != NULL )
      fprintf(fm, "<PDataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" />\n");

    fprintf(fm, "</PPointData>\n"
	"<PCellData>\n"
	"<PDataArray type=\"Int32\"   Name=\"part\"   NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\"/>\n"
	"<PDataArray type=\"Int32\"   Name=\"elem_type\" NumberOfComponents=\"1\"/>\n"
	"</PCellData>\n" , nvoi , nvoi);

    int i;
    for( i = 0 ; i < nproc_mac ; i++ ){
      sprintf(file_name,"%s_%d", name, i );
      fprintf(fm, "<Piece Source=\"%s.vtu\"/>\n", file_name );
    }
    fprintf(fm,	"</PUnstructuredGrid>\n</VTKFile>\n" );

    fclose(fm);

  } // rank = 0

  sprintf( file_name, "%s_%d.vtu", name, rank_mac);
  fm = fopen(file_name,"w"); 
  if(!fm){
    myio_printf(&PETSC_COMM_WORLD,"Problem trying to opening file %s for writing\n", file_name);
    return 1;
  }

  fprintf(fm, 
      "<?xml version=\"1.0\"?>\n"
      "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      "<UnstructuredGrid>\n");
  fprintf(fm,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nallnods, nelm);
  fprintf(fm,"<Points>\n");
  fprintf(fm,"<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");

  int    n , d;

  for( n = 0 ; n < nallnods ; n++ ){
    for( d = 0 ; d < dim ; d++ )
      fprintf( fm,"% 01.6e ",  coord[n*dim + d] );
    for( d = dim ; d < 3 ; d++ )
      fprintf( fm, "% 01.6e ", 0.0 );
    fprintf( fm, "\n" );
  }
  fprintf(fm,"</DataArray>\n");
  fprintf(fm,"</Points>\n");
  fprintf(fm,"<Cells>\n");

  int npe, e;

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for ( e = 0 ; e < nelm ; e++ ){
    npe = eptr[e+1] - eptr[e];
    for ( n = 0 ; n < npe ; n++ )
      fprintf(fm,"%-6d ", eind[eptr[e]+n]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  int ce = 0.0;
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for ( e = 0 ; e < nelm ; e++ ){
    npe = eptr[e+1] - eptr[e];
    ce += npe;
    fprintf(fm,"%d ", ce);
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for ( e = 0 ; e < nelm ; e++ )
    fprintf(fm, "%-3d ", vtkcode( dim , npe ) );
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"</Cells>\n");

  fprintf(fm,"<PointData Vectors=\"displ\">\n"); // Vectors inside is a filter we should not use this here

  if( x != NULL ){
    VecGhostUpdateBegin( x , INSERT_VALUES , SCATTER_FORWARD);
    VecGhostUpdateEnd(   x , INSERT_VALUES , SCATTER_FORWARD);
    VecGhostGetLocalForm(x , &xlocal );

    fprintf(fm,"<DataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
    VecGetArray( xlocal , &xvalues );
    for( n = 0 ; n < nallnods ; n++ ){
      for( d = 0 ; d < dim ; d++ )
	fprintf(fm, "% 01.6e ", xvalues[ n * dim + d ]);
      for( d = dim ; d < 3 ; d++ )
	fprintf(fm,"% 01.6e ",0.0);
      fprintf(fm,"\n");
    }
    VecRestoreArray( xlocal , &xvalues );
    fprintf(fm,"</DataArray>\n");
  }

  if( b != NULL ){
    VecGhostUpdateBegin( b , INSERT_VALUES,SCATTER_FORWARD);
    VecGhostUpdateEnd  ( b , INSERT_VALUES,SCATTER_FORWARD);
    VecGhostGetLocalForm(b , &xlocal);

    fprintf(fm,"<DataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
    VecGetArray(xlocal, &xvalues);
    for( n = 0 ; n < nallnods ; n++ ){
      for( d = 0 ; d < dim ; d++ )
	fprintf(fm, "% 01.6e ", xvalues[ n * dim + d ]);
      for( d = dim ; d < 3 ; d++ )
	fprintf(fm, "% 01.6e ", 0.0);
      fprintf(fm,"\n");
    }
    VecRestoreArray( xlocal , &xvalues );
    fprintf(fm,"</DataArray>\n");

  }
  fprintf(fm,"</PointData>\n");
  fprintf(fm,"<CellData>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"part\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for( e = 0; e < nelm ; e++ )
    fprintf( fm, "%d ", rank_mac );  
  fprintf( fm, "\n");
  fprintf( fm, "</DataArray>\n");

  int v;

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for( e = 0 ; e < nelm ; e++ ){
    for( v = 0 ; v < nvoi ; v++ )
      fprintf( fm, "% 01.6e ", elem_strain[ e*nvoi + v ]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for( e = 0; e < nelm ; e++ ){
    for( v = 0 ; v < nvoi ; v++ )
      fprintf(fm, "% 01.6e ", elem_stress[ e*nvoi + v ]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"elem_type\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for( e = 0; e < nelm ; e++ )
    fprintf( fm, "%d ", elem_type[e] );
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,
      "</CellData>\n"
      "</Piece>\n"
      "</UnstructuredGrid>\n"
      "</VTKFile>\n");

  fclose(fm);
  return 0;
}
