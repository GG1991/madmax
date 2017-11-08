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

int main(int argc, char **argv)
{

  int        ierr, ierr_1 = 0;
  double     tf, dt;
  char       mesh_n[NBUF];
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
  PetscOptionsHasName(NULL,NULL,"-testcomm",&set);
  if(set == PETSC_TRUE){
    macro_mode = TEST_COMM;
    PetscPrintf(MACRO_COMM,"MACRO MODE : TEST_COMM\n");
  }
  PetscOptionsHasName(NULL,NULL,"-eigensys",&set);
  if(set == PETSC_TRUE){
    macro_mode = EIGENSYSTEM;
    PetscPrintf(MACRO_COMM,"MACRO MODE : EIGENSYSTEM\n");
  }

  /* Mesh and Input Options */
  mesh_f = FORMAT_NULL;
  PetscOptionsHasName(NULL,NULL,"-mesh_gmsh",&set);
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
  nvoi = (dim == 2) ? 3 : 6;


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
  if(set==PETSC_FALSE) nx_interp=2;
  PetscOptionsGetInt(NULL, NULL, "-ny_interp", &ny_interp, &set);
  if(set==PETSC_FALSE) ny_interp=2;
  PetscOptionsGetInt(NULL, NULL, "-nz_interp", &nz_interp, &set);
  if(set==PETSC_FALSE) nz_interp=2;

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
  partition_algorithm = PARMETIS_MESHKWAY;
  PetscOptionsHasName(NULL,NULL,"-part_meshkway",&set);
  if( set == PETSC_TRUE ) partition_algorithm = PARMETIS_MESHKWAY;
  PetscOptionsHasName(NULL,NULL,"-part_geom",&set);
  if( set == PETSC_TRUE ) partition_algorithm = PARMETIS_GEOM;

  file_out = NULL;
  if(rank_mac==0) file_out = fopen("macro_structures.dat","w");

  if(flag_coupling)
    PetscPrintf(MACRO_COMM,"MACRO: COUPLING\n");
  else
    PetscPrintf(MACRO_COMM,"MACRO: STANDALONE\n");

  /* read mesh */    
  PetscPrintf(MACRO_COMM,"Reading mesh elements\n");
  ierr = read_mesh_elmv(MACRO_COMM, myname, mesh_n, mesh_f);
  if(ierr){
    goto end_mac_1;
  }

  /* partition the mesh */
  PetscPrintf(MACRO_COMM,"Partitioning and distributing mesh\n");
  ierr = part_mesh_PARMETIS(&MACRO_COMM, time_fl, myname, NULL);
  if(ierr){
    goto end_mac_1;
  }

  /* Calculate <*ghosts> and <nghosts> */
  PetscPrintf(MACRO_COMM,"Calculating Ghost Nodes\n");
  ierr = calculate_ghosts(&MACRO_COMM, myname);
  if(ierr){
    goto end_mac_1;
  }

  /* Reenumerate Nodes */
  PetscPrintf(MACRO_COMM,"Reenumering nodes\n");
  ierr = reenumerate_PETSc(MACRO_COMM);
  if(ierr){
    goto end_mac_1;
  }

  /* Coordinate Reading */
  PetscPrintf(MACRO_COMM,"Reading Coordinates\n");
  ierr = read_mesh_coord(MACRO_COMM, mesh_n, mesh_f);
  if(ierr){
    goto end_mac_1;
  }

  list_init(&physical_list, sizeof(physical_t), NULL);
  list_init(&function_list, sizeof(physical_t), NULL);

  /* Read Physical entities */
  ierr = read_physical_entities(MACRO_COMM, mesh_n, mesh_f);
  if(ierr){
    PetscPrintf(MACRO_COMM,"Problem parsing physical entities from mesh file\n");
    goto end_mac_1;
  }

  /* Read boundaries */
  ierr = set_id_on_material_and_boundary(MACRO_COMM);
  if(ierr){
    PetscPrintf(MACRO_COMM,"Problem determing ids on materials and boundaries\n");
    goto end_mac_1;
  }
  ierr = check_elm_id();
  if(ierr){
    PetscPrintf(MACRO_COMM,"Problem checking physical ids\n");
    goto end_mac_1;
  }
  
  read_boundary(MACRO_COMM, mesh_n, mesh_f);
  mac_init_boundary(MACRO_COMM, &boundary_list);

  /* Allocate matrices & vectors */ 
  PetscPrintf(MACRO_COMM, "allocating matrices & vectors\n");
  mac_alloc(MACRO_COMM);

  /* Setting solver options */
  KSPCreate(MACRO_COMM,&ksp);
  KSPSetType(ksp,KSPCG);
  KSPSetFromOptions(ksp);
  KSPSetOperators(ksp,A,A);

  /* Init Gauss point shapes functions and derivatives */
  ierr = fem_inigau();

  int      i, j, nr_its = -1;
  double   norm = -1.0;
  double   limit[6];

  ierr = get_bbox_local_limits(coord, nallnods, &limit[0], &limit[2], &limit[4]);
  PetscPrintf(MACRO_COMM,"Limit = ");
  for( i = 0 ; i < nvoi ; i++ )
    PetscPrintf(MACRO_COMM,"%lf ",limit[i]);
  PetscPrintf(MACRO_COMM,"\n");

  // Initial condition <x> = 0
  VecZeroEntries(x);

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
      for( j = 0 ; j < nvoi ; j++ ){
	PetscPrintf(MACRO_COMM,"%e ",stress_mac[j]);
      }
      PetscPrintf(MACRO_COMM,"\n");
    }

  }
  else if( macro_mode == EIGENSYSTEM ){

    int    nlocal, ntotal;
    double omega;
    EPS   eps;

    nlocal = dim * nmynods;
    ntotal = dim * ntotnod;

    MatCreate(MACRO_COMM,&M);
    MatSetSizes(M,nlocal,nlocal,ntotal,ntotal);
    MatSetFromOptions(M);
    MatSeqAIJSetPreallocation(M,117,NULL);
    MatMPIAIJSetPreallocation(M,117,NULL,117,NULL);
    MatSetOption(M,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    ierr = assembly_mass(&M);
    if(ierr){
      PetscPrintf(MACRO_COMM, "problem assembling mass matrix\n");
      goto end_mac_1;
    }
    ierr = assembly_jacobian_sd(&A);

    int  *dir_idx;
    int  ndir;
    node_list_t    *pBound;
    mac_boundary_t *mac_boundary;

    pBound = boundary_list.head;
    while(pBound){
      mac_boundary  = ((boundary_t*)pBound->data)->bvoid;
      dir_idx = mac_boundary->dir_idx;
      ndir = mac_boundary->ndir;
      MatZeroRowsColumns(M, ndir, dir_idx, 1.0, NULL, NULL);
      MatZeroRowsColumns(A, ndir, dir_idx, 1.0, NULL, NULL);
      pBound = pBound->next;
    }
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd  (M, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);

    if(flag_print & (1<<PRINT_PETSC)){
      PetscViewer  viewer;
      PetscViewerASCIIOpen(MACRO_COMM,"M.dat" ,&viewer); MatView(M ,viewer);
      PetscViewerASCIIOpen(MACRO_COMM,"A.dat" ,&viewer); MatView(A ,viewer);
      PetscViewerASCIIOpen(MACRO_COMM,"b.dat" ,&viewer); VecView(b ,viewer);
      PetscViewerASCIIOpen(MACRO_COMM,"dx.dat",&viewer); VecView(dx,viewer);
      PetscViewerASCIIOpen(MACRO_COMM,"x.dat" ,&viewer); VecView(x ,viewer);
    }

    int nev;   // number of request eigenpairs
    int nconv; // number of converged eigenpairs
    double error;

    EPSCreate(MACRO_COMM,&eps);
    EPSSetOperators(eps,M,A);
    EPSSetProblemType(eps,EPS_GHEP);
    EPSSetFromOptions(eps);
    EPSGetDimensions(eps,&nev,NULL,NULL);
    PetscPrintf(MACRO_COMM,"Number of requested eigenvalues: %D\n",nev);

    ierr = EPSSolve(eps);
    ierr = EPSGetConverged(eps,&nconv);
    PetscPrintf(MACRO_COMM,"Number of converged eigenpairs: %D\n",nconv);

    for( i = 0 ; i < nev ; i++ ){

      EPSGetEigenpair(eps,i,&omega,NULL,x,NULL);
      EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);
      PetscPrintf(MACRO_COMM, "omega %d = %e   error = %e\n", i, omega, error);

      if(flag_print & (1<<PRINT_VTU)){ 

	strain = malloc(nelm*nvoi*sizeof(double));
	stress = malloc(nelm*nvoi*sizeof(double));
	energy = malloc(nelm*sizeof(double));
	energy_interp = malloc(nelm*sizeof(double));
	ierr = assembly_residual_sd(&x, &b);
	ierr = calc_strain_stress_energy(&x, strain, stress, energy);
	ierr = interpolate_structured_2d(limit, nx_interp, ny_interp, energy, energy_interp);

	if(flag_print & (1<<PRINT_VTU)){ 
	  sprintf(filename,"%s_eigen_%d", myname, i);
	  ierr = write_vtu(MACRO_COMM, filename, &x, &b, strain, stress, energy);
	}
	free(stress); free(strain); free(energy);

      }

    }

    EPSDestroy(&eps);

  }
  else if( macro_mode == NORMAL ){

    /* Begin time dependent loop */
    
    double   t = t0;
    int      time_step = 0;

    while( t < (tf + 1.0e-10) ){

      PetscPrintf(MACRO_COMM,"\ntime step %3d %e seg\n", time_step, t);

      /* Setting Displacement on Dirichlet Indeces on <x> */
      MacroSetDisplacementOnBoundary(t, &x);

      nr_its = 0; norm = 2*nr_norm_tol;
      while( nr_its < nr_max_its && norm > nr_norm_tol )
      {

	/* Assemblying residual */
	PetscPrintf(MACRO_COMM, "Assembling Residual\n");
	ierr = assembly_residual_sd(&x, &b);
	if(ierr){
	  PetscPrintf(MACRO_COMM, "problem assembling residual\n");
	  goto end_mac_1;
	}
	MacroSetBoundaryOnResidual(&b);
	VecNorm(b,NORM_2,&norm);
	PetscPrintf(MACRO_COMM,"|b| = %e\n",norm);
	VecScale( b, -1.0 );

	if( norm < nr_norm_tol )break;

	/* Assemblying Jacobian */
	PetscPrintf(MACRO_COMM, "Assembling Jacobian\n");
	ierr = assembly_jacobian_sd(&A);
	ierr = MacroSetBoundaryOnJacobian(&A); CHKERRQ(ierr);

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
	strain = malloc( nelm*nvoi * sizeof(double));
	stress = malloc( nelm*nvoi * sizeof(double));
	energy = malloc( nelm      * sizeof(double));
	energy_interp = malloc(nelm* sizeof(double));
	ierr = assembly_residual_sd(&x, &b);
	ierr = calc_strain_stress_energy(&x, strain, stress, energy);
	ierr = interpolate_structured_2d(limit, nx_interp, ny_interp, energy, energy_interp);
	sprintf( filename, "macro_t_%d", time_step );
	write_vtu( MACRO_COMM, filename, &x, &b, strain, stress, energy );
	free(stress); free(strain); free(energy);
      }

      t += dt;
      time_step ++;
    }
  }

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
