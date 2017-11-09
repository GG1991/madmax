/*
   MACRO main function

   Program for solving the displacement field inside a solid 
   structure representing the macrostructure

   Author> Guido Giuntoli
   Date> 28-07-2017
 */

#include "macro.h"

static char help[] = 
"MACRO MULTISCALE CODE                                                                             \n"
"Solves the displacement field inside a solid structure.                                           \n"
"-coupl       : coupled with \"micro\" code for solving multiscale problem                         \n"
"-normal      : normal execution, solves a time dependent boundary condition problem               \n"
"-testcomm    : communication testing with the \"micro\" code                                      \n"
"-eigensys    : calculates the eigensystem Mx = -(1/omega)Kx                                       \n"
"-print_petsc : prints petsc structures on files such as Mat and Vec objects                       \n"
"-print_vtu   : prints solutions on .vtu and .pvtu files                                           \n";

int read_bc     ( void     );
int assembly_A  ( void     );
int assembly_M  ( void     );
int assembly_b  ( void     );
int update_bound( double t );
int get_global_elem_index(int e, int * glo_elem_index);
int get_local_elem_index (int e, int * loc_elem_index);

int main(int argc, char **argv)
{

  int        ierr, ierr_1 = 0;
  int        i, j;
  double     tf, dt;
  PetscBool  set;

  myname = strdup("macro");

  WORLD_COMM = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(WORLD_COMM, &nproc_wor);
  MPI_Comm_rank(WORLD_COMM, &rank_wor);

  /* 
     We start PETSc before coloring here for using command line reading tools only
     Then we finalize it
  */
  PETSC_COMM_WORLD = WORLD_COMM;
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);

  /* Coupling Options */
  flag_coupling = PETSC_FALSE;
  PetscOptionsHasName(NULL,NULL,"-coupl",&set);
  macmic.type = 0;
  if( set == PETSC_TRUE ){
    flag_coupling = PETSC_TRUE;
    macmic.type = COUP_1;
  }

  /* Stablish a new local communicator */
  color = MACRO;
  ierr = macmic_coloring(WORLD_COMM, &color, &macmic, &MACRO_COMM);
  if( ierr ) ierr_1 = 1; goto end_mac_0;

  MPI_Comm_size(MACRO_COMM, &nproc_mac);
  MPI_Comm_rank(MACRO_COMM, &rank_mac);
  
end_mac_0:
  PetscFinalize();
  if(ierr_1) goto end_mac_2;

  /*
     Set PETSc communicator to MACRO_COMM
     and start again
  */
  PETSC_COMM_WORLD = MACRO_COMM;

#ifdef SLEPC
  SlepcInitialize(&argc,&argv,(char*)0,help);
#elif  PETSC
  PetscInitialize(&argc,&argv,(char*)0,help);
#endif

  PetscPrintf(MACRO_COMM,
      "--------------------------------------------------\n"
      "  MACRO: COMPOSITE MATERIAL MULTISCALE CODE\n"
      "--------------------------------------------------\n");

  /* execution mode */
  macro_mode  = NORMAL;
  PetscOptionsHasName( NULL, NULL, "-normal", &set );
  if( set == PETSC_TRUE ){
    macro_mode = NORMAL;
    PetscPrintf( MACRO_COMM, "MACRO MODE : NORMAL\n" );
  }
  PetscOptionsHasName( NULL, NULL, "-testcomm", &set );
  if( set == PETSC_TRUE ){
    macro_mode = TEST_COMM;
    PetscPrintf( MACRO_COMM, "MACRO MODE : TEST_COMM\n" );
  }
  PetscOptionsHasName( NULL,NULL,"-eigensys",&set);
  if( set == PETSC_TRUE ){
    macro_mode = EIGENSYSTEM;
    PetscPrintf(MACRO_COMM,"MACRO MODE : EIGENSYSTEM\n");
  }

  /* Mesh and Input Options */
  mesh_f = FORMAT_NULL;
  PetscOptionsHasName( NULL, NULL, "-mesh_gmsh",&set);
  if( set == PETSC_TRUE ) mesh_f = FORMAT_GMSH;
  PetscOptionsHasName(NULL,NULL,"-mesh_alya",&set);
  if( set == PETSC_TRUE ) mesh_f = FORMAT_ALYA;
  if( mesh_f == FORMAT_NULL ){
    PetscPrintf(MPI_COMM_SELF,"mesh format not given on command line.\n");
    goto end_mac_1;
  }
  PetscOptionsGetString(NULL, NULL, "-mesh", mesh_n, 128, &set);
  if( set == PETSC_FALSE ){
    PetscPrintf(MPI_COMM_SELF,"mesh file not given on command line.\n");
    goto end_mac_1;
  }
  PetscOptionsGetInt(NULL, NULL, "-dim", &dim, &set);
  if( set == PETSC_FALSE ){
    PetscPrintf(MPI_COMM_SELF,"dimension (-dim <dim>) not given\n");
    goto end_mac_1;
  }
  nvoi    = (dim == 2) ? 3 : 6;
  npe_max = (dim == 2) ? 4 : 8;

  /* Printing Options */
  flag_print = 0;
  PetscOptionsHasName(NULL,NULL,"-print_petsc",&set);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_PETSC);
  PetscOptionsHasName(NULL,NULL,"-print_vtu",&set);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_VTU);
  
  /* Newton-Raphson solver options */
  PetscOptionsGetInt(NULL, NULL,  "-nr_max_its", &nr_max_its, &set);
  if(set==PETSC_FALSE) nr_max_its=5;
  PetscOptionsGetReal(NULL, NULL, "-nr_norm_tol", &nr_norm_tol, &set);
  if(set==PETSC_FALSE) nr_norm_tol=1.0e-7;

  /* structured grid interp */
  PetscOptionsGetInt(NULL, NULL, "-nx_interp", &nx_interp, &set);
  if(set==PETSC_FALSE) nx_interp = 2;
  PetscOptionsGetInt(NULL, NULL, "-ny_interp", &ny_interp, &set);
  if(set==PETSC_FALSE) ny_interp = 2;
  PetscOptionsGetInt(NULL, NULL, "-nz_interp", &nz_interp, &set);
  if(set==PETSC_FALSE) nz_interp = 2;

  /* read function */
  { 
    int     nval = 16;
    double  values[nval];
    PetscOptionsGetRealArray( NULL, NULL, "-function", values, &nval, &set);
    if( set == PETSC_TRUE )
    {
      if( nval % 2 != 0 ){
	PetscPrintf(MPI_COMM_SELF,"odd number of argument for -function\n");
	ierr_1 = 1;
	goto end_mac_0;
      }
      func_bc.fnum = 1;
      func_bc.n    = nval/2;
      func_bc.x    = malloc( func_bc.n * sizeof(double));
      func_bc.y    = malloc( func_bc.n * sizeof(double));
      for( i = 0 ; i < func_bc.n ; i++ ){
	func_bc.x[i] = values[2*i+0];
	func_bc.y[i] = values[2*i+1];
      }
    }
    else if( macro_mode == NORMAL ){
      PetscPrintf(MPI_COMM_SELF,"-function is request to impose non trivial BC.\n");
      ierr_1 = 1;
      goto end_mac_0;
    }
  }

  /* read boundary elements */
  { 

    int      nval = 4;
    char     *string[nval];
    char     *data;
    bound_t  bou;

    list_init( &boundary_list, sizeof(bound_t), NULL );
    PetscOptionsGetStringArray( NULL, NULL, "-boundary", string, &nval, &set );
    if( set == PETSC_TRUE )
    {
      for( i = 0 ; i < nval ; i++ ){
	data = strtok(string[i], " \n");
	bou.name     = strdup(data); 
	data = strtok(NULL, " \n");
	bou.kind     = strbin2dec(data);
	if(dim == 2){
	  if( bou.kind == 1 || bou.kind == 2 ) bou.ndirpn =  1;
	  if( bou.kind == 3 )                  bou.ndirpn =  2;
	  if( bou.kind == 0 )                  bou.ndirpn =  0;
	}
	bou.nneupn   = dim - bou.ndirpn;
	bou.fnum     = malloc(dim * sizeof(int));
	for( j = 0 ; j < dim ; j++ ){
	  data = strtok(NULL, " \n");
	  bou.fnum[j] = atoi(data);
	}
	bou.dir_ixs = NULL;
	bou.dir_val = NULL;
	list_insertlast( &boundary_list, &bou );
      }
    }
    else if( macro_mode == NORMAL ){
      PetscPrintf(MACRO_COMM,"-boundary should be set.\n");
      ierr_1 = 1;
      goto end_mac_0;
    }
  }

  /* flow execution variables for NORMAL mode */
  if( macro_mode == NORMAL ){
    PetscOptionsGetReal(NULL,NULL,"-tf",&tf,&set);
    if(set == PETSC_FALSE){
      PetscPrintf(MPI_COMM_SELF,"-tf not given.\n");
      goto end_mac_1;
    }
    PetscOptionsGetReal(NULL,NULL,"-dt",&dt,&set);
    if(set == PETSC_FALSE){
      PetscPrintf(MPI_COMM_SELF,"-dt not given.\n");
      goto end_mac_1;
    }
  }

  /* Mesh partition algorithms */
  partition_algorithm = PARMETIS_GEOM;
  PetscOptionsHasName(NULL,NULL,"-part_meshkway",&set);
  if( set == PETSC_TRUE ) partition_algorithm = PARMETIS_MESHKWAY;
  PetscOptionsHasName(NULL,NULL,"-part_geom",&set);
  if( set == PETSC_TRUE ) partition_algorithm = PARMETIS_GEOM;

  if(flag_coupling)
    PetscPrintf(MACRO_COMM,"MACRO: COUPLING\n");
  else
    PetscPrintf(MACRO_COMM,"MACRO: STANDALONE\n");

  /* read mesh */    
  PetscPrintf( MACRO_COMM, "Reading mesh elements\n" );
  read_mesh_elmv(MACRO_COMM, myname, mesh_n, mesh_f);

  /* partition the mesh */
  PetscPrintf(MACRO_COMM,"Partitioning and distributing mesh\n");
  part_mesh_PARMETIS(&MACRO_COMM, myname, NULL);

  /* calculate <*ghosts> and <nghosts> */
  PetscPrintf(MACRO_COMM,"Calculating Ghost Nodes\n");
  calculate_ghosts(&MACRO_COMM, myname);

  /* reenumerate Nodes */
  PetscPrintf(MACRO_COMM,"Reenumering nodes\n");
  reenumerate_PETSc(MACRO_COMM);

  /* coordinate Reading */
  PetscPrintf(MACRO_COMM,"Reading Coordinates\n");
  read_mesh_coord(MACRO_COMM, mesh_n, mesh_f);

  /* read boundaries */
  read_bc();
  

  /**************************************************/
  /* alloc and init variables */
  PetscPrintf(MACRO_COMM, "allocating matrices & vectors\n");
  mac_alloc(MACRO_COMM);

  A   = NULL;
  b   = NULL;
  x   = NULL;
  dx  = NULL;

  /* alloc variables*/
  int ixpe = npe_max * dim;  // number of indeces per element
  loc_elem_index = malloc( ixpe * sizeof(int));
  glo_elem_index = malloc( ixpe * sizeof(int));
  elem_disp      = malloc( ixpe * sizeof(double));
  k_elem         = malloc( ixpe*ixpe * sizeof(double));
  stress_gp      = malloc( nvoi * sizeof(double));
  strain_gp      = malloc( nvoi * sizeof(double));
  c              = malloc( nvoi*nvoi * sizeof(double));
  if( flag_print & ( 1 << PRINT_VTU ) ){
    elem_strain  = malloc( nelm*nvoi * sizeof(double));
    elem_stress  = malloc( nelm*nvoi * sizeof(double));
    elem_energy  = malloc( nelm      * sizeof(double));
    elem_type    = malloc( nelm      * sizeof(int));
  }

  /* alloc the B matrix */
  bmat = malloc( nvoi * sizeof(double*));
  for( i = 0 ; i < nvoi  ; i++ )
    bmat[i] = malloc( ixpe * sizeof(double));


  /**************************************************/

  /* Setting solver options */
  KSPCreate(MACRO_COMM,&ksp);
  KSPSetFromOptions(ksp);

  /* Init Gauss point shapes functions and derivatives */
  ierr = fem_inigau();

  int      nr_its = -1;
  double   norm = -1.0;
  double   limit[6];

  ierr = get_bbox_local_limits(coord, nallnods, &limit[0], &limit[2], &limit[4]);
  PetscPrintf(MACRO_COMM,"Limit = ");
  for( i = 0 ; i < nvoi ; i++ )
    PetscPrintf(MACRO_COMM,"%lf ",limit[i]);
  PetscPrintf(MACRO_COMM,"\n");

  if(macro_mode == TEST_COMM){

    double   strain_mac[6] = {0.1, 0.1, 0.2, 0.0, 0.0, 0.0};
    double   stress_mac[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    /* TEST> Sends a calculating strain to micro and obtain the stress */
    for( i = 0 ; i < nvoi ; i++ ){
      for( j = 0 ; j < nvoi ; j++ )
	strain_mac[j] = 0.0;
      strain_mac[i] = 0.005;
      ierr = mac_send_signal(WORLD_COMM, MAC2MIC_STRAIN);
      ierr = mac_send_strain(WORLD_COMM, strain_mac    );
      ierr = mac_recv_stress(WORLD_COMM, stress_mac    );
      PetscPrintf(MACRO_COMM,"\nstress_ave = ");
      for( j = 0 ; j < nvoi ; j++ )
	PetscPrintf(MACRO_COMM,"%e ",stress_mac[j]);
      PetscPrintf(MACRO_COMM,"\n");
    }

  }
  else if( macro_mode == EIGENSYSTEM ){

    int     nlocal, ntotal;
    double  omega;
    EPS     eps;

    nlocal = dim * nmynods;
    ntotal = dim * ntotnod;

    MatCreate(MACRO_COMM,&M);
    MatSetSizes(M,nlocal,nlocal,ntotal,ntotal);
    MatSetFromOptions(M);
    MatSeqAIJSetPreallocation(M,117,NULL);
    MatMPIAIJSetPreallocation(M,117,NULL,117,NULL);
    MatSetOption(M,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    assembly_M();
    assembly_A();

    /* set dirichlet bc */
    node_list_t *pn;
    bound_t     *bou;
    pn = boundary_list.head;
    while( pn )
    {
      bou  = (bound_t * )pn->data;
      MatZeroRowsColumns(M, bou->ndirix, bou->dir_ixs, 1.0, NULL, NULL);
      MatZeroRowsColumns(A, bou->ndirix, bou->dir_ixs, 1.0, NULL, NULL);
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
      PetscViewerASCIIOpen(MACRO_COMM,"b.dat" ,&viewer); VecView(b ,viewer);
      PetscViewerASCIIOpen(MACRO_COMM,"dx.dat",&viewer); VecView(dx,viewer);
      PetscViewerASCIIOpen(MACRO_COMM,"x.dat" ,&viewer); VecView(x ,viewer);
    }

    int    nev;   // number of request eigenpairs
    int    nconv; // number of converged eigenpairs
    double error;

    EPSCreate(MACRO_COMM,&eps);
    EPSSetOperators(eps,M,A);
    EPSSetProblemType(eps,EPS_GHEP);
    EPSSetFromOptions(eps);
    EPSGetDimensions(eps,&nev,NULL,NULL);
    PetscPrintf(MACRO_COMM,"Number of requested eigenvalues: %D\n",nev);

    EPSSolve(eps);
    EPSGetConverged(eps,&nconv);
    PetscPrintf(MACRO_COMM,"Number of converged eigenpairs: %D\n",nconv);

    for( i = 0 ; i < nev ; i++ ){

      EPSGetEigenpair(eps,i,&omega,NULL,x,NULL);
      EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);
      PetscPrintf(MACRO_COMM, "omega %d = %e   error = %e\n", i, omega, error);

      if(flag_print & (1<<PRINT_VTU))
      { 
	ierr = calc_strain_stress_energy(&x, strain, stress, energy);
	ierr = interpolate_structured_2d(limit, nx_interp, ny_interp, energy, energy_interp);
	sprintf( filename, "macro_eigen_%d", i);
	ierr = write_vtu(MACRO_COMM, filename, &x, &b, strain, stress, energy);
      }

    }

    EPSDestroy(&eps);

  }
  else if( macro_mode == NORMAL ){

    /* Begin time dependent loop */
    
    double   t = 0.0;
    int      time_step = 0;

    while( t < (tf + 1.0e-10) ){

      PetscPrintf(MACRO_COMM,"\ntime step %3d %e seg\n", time_step, t);

      /* Setting Displacement on Dirichlet Indeces on <x> */

      nr_its = 0; norm = 2*nr_norm_tol;
      while( nr_its < nr_max_its && norm > nr_norm_tol )
      {

	/* Assemblying residual */
	PetscPrintf(MACRO_COMM, "Assembling Residual\n");
	assembly_b();
	VecNorm(b,NORM_2,&norm);
	PetscPrintf(MACRO_COMM,"|b| = %e\n",norm);
	VecScale( b, -1.0 );

	if( norm < nr_norm_tol )break;

	/* Assemblying Jacobian */
	PetscPrintf(MACRO_COMM, "Assembling Jacobian\n");
	assembly_A();

	/* Solving Problem */
	PetscPrintf(MACRO_COMM, "Solving Linear System\n ");
	KSPSolve(ksp,b,dx);
	print_ksp_info( MACRO_COMM, ksp );
	VecAXPY(x, 1.0, dx);


	if(flag_print & (1<<PRINT_PETSC)){
	  PetscViewer  viewer;
	  PetscViewerASCIIOpen(MACRO_COMM,"M.dat" ,&viewer); MatView(M ,viewer);
	  PetscViewerASCIIOpen(MACRO_COMM,"A.dat" ,&viewer); MatView(A ,viewer);
	  PetscViewerASCIIOpen(MACRO_COMM,"b.dat" ,&viewer); VecView(b ,viewer);
	  PetscViewerASCIIOpen(MACRO_COMM,"dx.dat",&viewer); VecView(dx,viewer);
	  PetscViewerASCIIOpen(MACRO_COMM,"x.dat" ,&viewer); VecView(x ,viewer);
	}

	nr_its ++;
      }

      if( flag_print & (1<<PRINT_VTU) )
      { 
	ierr = calc_strain_stress_energy( &x, strain, stress, energy);
	sprintf( filename, "macro_t_%d", time_step );
	write_vtu( MACRO_COMM, filename, &x, &b, strain, stress, energy );
      }

      t += dt;
      time_step ++;
    }
  }

  /**************************************************/
  /* free variables*/
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

  /* free the B matrix */
  for( i = 0 ; i < nvoi  ; i++ )
    free(bmat[i]);
  free(bmat);
  /**************************************************/

end_mac_1:

  /* Stop signal to micro if it is coupled */
  if(flag_coupling){
    ierr = mac_send_signal(WORLD_COMM, MIC_END);
    if(ierr){
      PetscPrintf(PETSC_COMM_WORLD, "macro: problem sending MIC_END to micro\n");
      return 1;
    }
  }


  list_clear(&material_list);
  list_clear(&physical_list);

  MatDestroy(&A);
  VecDestroy(&x);
  VecDestroy(&b);
  KSPDestroy(&ksp);

  PetscPrintf(MACRO_COMM,
      "--------------------------------------------------\n"
      "  MACRO: FINISH COMPLETE\n"
      "--------------------------------------------------\n");

#ifdef SLEPC
  SlepcFinalize();
#elif  PETSC
  PetscFinalize();
#endif

end_mac_2:
  ierr = MPI_Finalize();

  return 0;
}

/****************************************************************************************************/

int read_bc()
{

  /* 
     completes the field of each "bound_t" in boundary list 
     int     nix;
     int     *disp_ixs;
     double  *disp_val;
   */

  node_list_t *pn;
  bound_t     *bou;
  int         i, d, da;

  pn = boundary_list.head;
  while( pn ){
    bou = ( bound_t * )pn->data;
    int *ix, n;
    gmsh_get_node_index( mesh_n, bou->name, nmynods, mynods, dim, &n, &ix );
    bou->ndirix  = n * bou->ndirpn;
    bou->dir_ixs = malloc( bou->ndirix * sizeof(int)); 
    bou->dir_val = malloc( bou->ndirix * sizeof(double)); 
    for( i = 0 ; i < n ; i++ ){
      da = 0;
      for( d = 0 ; d < dim ; d++ )
	if( bou->kind & (1<<d) ) 
	  bou->dir_ixs[i* (bou->ndirpn) + da++] = ix[i] * dim + d;
    }
    free(ix);
    pn = pn->next;
  }

  return 0;
}

/****************************************************************************************************/

int assembly_A( void )
{

  MatZeroEntries(A);

  int     e, gp;
  int     i, j, k, h;
  int     npe, ngp;

  for( e = 0 ; e < nelm ; e++ ){

    ngp = npe = eptr[e+1] - eptr[e];

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
	      bmat[h][i] * c[ h*nvoi + k ] * struct_bmat[k][j] * wp[gp];
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

int get_global_elem_index(int e, int * glo_elem_index)
{
  return 0;
}

/****************************************************************************************************/

int get_local_elem_index(int e, int * loc_elem_index)
{
  return 0;
}

/****************************************************************************************************/
