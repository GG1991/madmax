/*
   MACRO main function

   Program for solving the displacement field inside a solid 
   structure representing the macrostructure

   Author> Guido Giuntoli
   Date> 28-07-2017
 */

static char help[] = 
"MACRO MULTISCALE CODE\n"
"Solves the displacement field inside a solid structure. \n"
"-coupl    [0 (no coupling ) | 1 (coupling with micro)]\n"
"-testcomm [0 (no test) | 1 (sends a strain value and receive a stress calculated from micro)]\n"
"-print_petsc prints petsc structures on files such as Mat and Vec objects\n"
"-print_vtu prints solutions on .vtu and .pvtu files\n";

#include "macro.h"

int main(int argc, char **argv)
{

  int        ierr, ierr_1=0;
  char       *myname = strdup("macro");
  char       vtkfile_n[NBUF];
  double     t0=0.0, tf, dt;
  PetscBool  set;

  WORLD_COMM = MPI_COMM_WORLD;
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(WORLD_COMM, &nproc_wor);
  ierr = MPI_Comm_rank(WORLD_COMM, &rank_wor);

  /* 
     We start PETSc before coloring here for using command line reading tools only
     Then we finalize it
  */
  PETSC_COMM_WORLD = WORLD_COMM;
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);

  /*
     Printing Options
  */
  flag_print = 0;
  ierr = PetscOptionsHasName(NULL,NULL,"-print_petsc",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_PETSC);
  ierr = PetscOptionsHasName(NULL,NULL,"-print_vtk",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_VTK);
  ierr = PetscOptionsHasName(NULL,NULL,"-print_part",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_VTKPART);
  ierr = PetscOptionsHasName(NULL,NULL,"-print_vtu",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_VTU);
  ierr = PetscOptionsHasName(NULL,NULL,"-print_all",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_ALL);

  /*
     Coupling Options
  */
  flag_coupling = PETSC_FALSE;
  ierr = PetscOptionsHasName(NULL,NULL,"-coupl",&set);CHKERRQ(ierr);
  macmic.type = 0;
  if(set == PETSC_TRUE){
    flag_coupling = PETSC_TRUE;
    macmic.type = COUP_1;
  }
  flag_testcomm  = TESTCOMM_NULL;
  ierr = PetscOptionsHasName(NULL,NULL,"-testcomm",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) flag_testcomm  = TESTCOMM_STRAIN;
  /*
     Mesh and Input Options
  */
  mesh_f = FORMAT_NULL;
  ierr = PetscOptionsHasName(NULL,NULL,"-mesh_gmsh",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) mesh_f = FORMAT_GMSH;
  ierr = PetscOptionsHasName(NULL,NULL,"-mesh_alya",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) mesh_f = FORMAT_ALYA;
  if(mesh_f == FORMAT_NULL)SETERRQ(MACRO_COMM,1,"mesh format not given on command line.");
  ierr = PetscOptionsGetString(NULL, NULL, "-mesh", mesh_n, 128, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE) SETERRQ(MACRO_COMM,1,"mesh file not given on command line.");
  ierr = PetscOptionsGetString(NULL, NULL, "-input", input_n, 128, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE) SETERRQ(MACRO_COMM,1,"input file not given.");
  ierr = PetscOptionsGetInt(NULL, NULL, "-dim", &dim, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE){
    PetscPrintf(MPI_COMM_SELF,"dimension (-dim <dim>) not given\n");
    ierr_1 = 1;
    goto end_mac_0;
  }
  switch(dim){
    case 2:
      nvoi=3;
      break;
    case 3:
      nvoi=6;
      break;
    default:
      PetscPrintf(MPI_COMM_SELF,"dimension number %d not allowded\n", dim);
      ierr_1 = 1;
      goto end_mac_1;
  }
  /* Solver Options */
  ierr = PetscOptionsGetInt(NULL, NULL, "-nr_max_its", &nr_max_its, &set);
  if(set==PETSC_FALSE) nr_max_its=5;
  ierr = PetscOptionsGetReal(NULL, NULL, "-nr_norm_tol", &nr_norm_tol, &set);
  if(set==PETSC_FALSE) nr_norm_tol=1.0e-7;

  /* flow execution variables */
  ierr = PetscOptionsGetReal(NULL,NULL,"-tf",&tf,&set);CHKERRQ(ierr);
  if(set == PETSC_FALSE) SETERRQ(MACRO_COMM,1,"-tf not given.");
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,&set);CHKERRQ(ierr);
  if(set == PETSC_FALSE) SETERRQ(MACRO_COMM,1,"-dt not given.");

  ierr = PetscOptionsHasName(NULL,NULL,"-reactions",&set);CHKERRQ(ierr);
  flag_reactions = (set==PETSC_TRUE) ? PETSC_TRUE : PETSC_FALSE;

  /* Mesh partition algorithms */
  partition_algorithm = PARMETIS_MESHKWAY;
  ierr = PetscOptionsHasName(NULL,NULL,"-part_meshkway",&set);CHKERRQ(ierr);
  if(set==PETSC_TRUE) partition_algorithm = PARMETIS_MESHKWAY;
  ierr = PetscOptionsHasName(NULL,NULL,"-part_geom",&set);CHKERRQ(ierr);
  if(set==PETSC_TRUE) partition_algorithm = PARMETIS_GEOM;



  /* 
     Stablish a new local communicator
  */
  color = MACRO;
  ierr = macmic_coloring(WORLD_COMM, &color, &macmic, &MACRO_COMM);
  ierr = MPI_Comm_size(MACRO_COMM, &nproc_mac);
  ierr = MPI_Comm_rank(MACRO_COMM, &rank_mac);

  
end_mac_0:
  ierr = PetscFinalize();CHKERRQ(ierr);
  if(ierr_1) goto end_mac_2;

  /*
     Set PETSc communicator to MACRO_COMM
     and start again
  */
  PETSC_COMM_WORLD = MACRO_COMM;
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);

  PetscPrintf(MACRO_COMM,
      "--------------------------------------------------\n"
      "  MACRO: COMPOSITE MATERIAL MULTISCALE CODE\n"
      "--------------------------------------------------\n");

  file_out = NULL;
  if(rank_mac==0) file_out = fopen("macro_structures.dat","w");

#if defined(PETSC_USE_LOG)
  ierr = PetscLogEventRegister("Read_Elems_of_Mesh"    ,PETSC_VIEWER_CLASSID,&EVENT_READ_MESH_ELEM);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Partition_Mesh"        ,PETSC_VIEWER_CLASSID,&EVENT_PART_MESH);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Calculate_Ghosts_Nodes",PETSC_VIEWER_CLASSID,&EVENT_CALC_GHOSTS);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Reenumerates_Nodes"    ,PETSC_VIEWER_CLASSID,&EVENT_REENUMERATE);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Read_Coordinates"      ,PETSC_VIEWER_CLASSID,&EVENT_READ_COORD);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Init_Gauss_Points"     ,PETSC_VIEWER_CLASSID,&EVENT_INIT_GAUSS);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Allocate_Mat_and_Vec"  ,PETSC_VIEWER_CLASSID,&EVENT_ALLOC_MATVEC);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Set_Displ_on_Bou "     ,PETSC_VIEWER_CLASSID,&EVENT_SET_DISP_BOU);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Assembly_Jacobian"     ,PETSC_VIEWER_CLASSID,&EVENT_ASSEMBLY_JAC);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Assembly_Residual"     ,PETSC_VIEWER_CLASSID,&EVENT_ASSEMBLY_RES);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Solve_Linear_System"   ,PETSC_VIEWER_CLASSID,&EVENT_SOLVE_SYSTEM);CHKERRQ(ierr);
#endif

  if(flag_coupling){
    PetscPrintf(MACRO_COMM,
	"--------------------------------------------------\n"
	"  MACRO: COUPLING \n"
	"--------------------------------------------------\n");
  }
  else{
    PetscPrintf(MACRO_COMM,
	"--------------------------------------------------\n"
	"  MACRO: STANDALONE \n"
	"--------------------------------------------------\n");
  }

  /*
     Register various stages for profiling
  */
  ierr = PetscLogStageRegister("Read Mesh Elements",&stages[0]);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Linear System 1",&stages[1]);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Linear System 2",&stages[2]);CHKERRQ(ierr);

  /*
     read mesh
  */    
  ierr = PetscLogStagePush(stages[0]);CHKERRQ(ierr);

  ierr = PetscLogEventBegin(EVENT_READ_MESH_ELEM,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscPrintf(MACRO_COMM,"Reading mesh elements\n");CHKERRQ(ierr);
  ierr = read_mesh_elmv(MACRO_COMM, myname, mesh_n, mesh_f);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_READ_MESH_ELEM,0,0,0,0);CHKERRQ(ierr);

  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /*
     partition the mesh
  */
  ierr = PetscPrintf(MACRO_COMM,"Partitioning and distributing mesh\n");CHKERRQ(ierr);
  ierr = PetscLogEventBegin(EVENT_PART_MESH,0,0,0,0);CHKERRQ(ierr);
  ierr = part_mesh_PARMETIS(&MACRO_COMM, time_fl, myname, NULL);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_PART_MESH,0,0,0,0);CHKERRQ(ierr);

  /*
     Calculate <*ghosts> and <nghosts> 
  */
  ierr = PetscPrintf(MACRO_COMM,"Calculating Ghost Nodes\n");CHKERRQ(ierr);
  ierr = PetscLogEventBegin(EVENT_CALC_GHOSTS,0,0,0,0);CHKERRQ(ierr);
  ierr = calculate_ghosts(&MACRO_COMM, myname);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_CALC_GHOSTS,0,0,0,0);CHKERRQ(ierr);

  /*
     Reenumerate Nodes
  */
  ierr = PetscPrintf(MACRO_COMM,"Reenumering nodes\n");CHKERRQ(ierr);
  ierr = PetscLogEventBegin(EVENT_REENUMERATE,0,0,0,0);CHKERRQ(ierr);
  ierr = reenumerate_PETSc(MACRO_COMM);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_REENUMERATE,0,0,0,0);CHKERRQ(ierr);

  /*
     Coordinate Reading
  */
  ierr = PetscPrintf(MACRO_COMM,"Reading Coordinates\n");CHKERRQ(ierr);
  ierr = PetscLogEventBegin(EVENT_READ_COORD,0,0,0,0);CHKERRQ(ierr);
  ierr = read_mesh_coord(MACRO_COMM, mesh_n, mesh_f);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_READ_COORD,0,0,0,0);CHKERRQ(ierr);

  if( flag_print & (1<<PRINT_VTKPART)){
    sprintf(vtkfile_n,"%s_part_%d.vtk",myname,rank_mac);
    ierr = spu_vtk_partition( vtkfile_n, &MACRO_COMM );CHKERRQ(ierr);
  }

  ierr = list_init(&physical_list, sizeof(physical_t), NULL);CHKERRQ(ierr);
  ierr = list_init(&function_list, sizeof(physical_t), NULL);CHKERRQ(ierr);
  /* Read materials  */
  ierr = parse_material(MACRO_COMM, input_n);
  if(ierr){
    ierr = PetscPrintf(MACRO_COMM,"Problem parsing materials from input file\n");
    goto end_mac_1;
  }
  /* Read Physical entities */
  ierr = read_physical_entities(MACRO_COMM, mesh_n, mesh_f);
  if(ierr){
    ierr = PetscPrintf(MACRO_COMM,"Problem parsing physical entities from mesh file\n");
    goto end_mac_1;
  }
  /* Read functions */
  ierr = parse_function(MACRO_COMM, input_n);
  if(ierr){
    ierr = PetscPrintf(MACRO_COMM,"Problem parsing functions from input file\n");
    goto end_mac_1;
  }
  /* Read boundaries */
  ierr = macro_parse_boundary(MACRO_COMM, input_n);
  if(ierr){
    ierr = PetscPrintf(MACRO_COMM,"Problem reading boundaries from input file\n");
    goto end_mac_1;
  }
  ierr = set_id_on_material_and_boundary(MACRO_COMM);
  if(ierr){
    ierr = PetscPrintf(MACRO_COMM,"Problem determing ids on materials and boundaries\n");
    goto end_mac_1;
  }
  ierr = CheckPhysicalID();
  if(ierr){
    ierr = PetscPrintf(MACRO_COMM,"Problem checking physical ids\n");
    goto end_mac_1;
  }
  
  ierr = read_boundary(MACRO_COMM, mesh_n, mesh_f);CHKERRQ(ierr);
  ierr = mac_init_boundary(MACRO_COMM, &boundary_list);

  /*
     Allocate matrices & vectors
  */ 
  ierr = PetscLogEventBegin(EVENT_ALLOC_MATVEC,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscPrintf(MACRO_COMM, "allocating matrices & vectors\n");CHKERRQ(ierr);
  ierr = mac_alloc(MACRO_COMM);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_ALLOC_MATVEC,0,0,0,0);CHKERRQ(ierr);

  /*
     Setting solver options 
  */
  ierr = KSPCreate(MACRO_COMM,&ksp); CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPCG); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A); CHKERRQ(ierr);

  /*
     Init Gauss point shapes functions and derivatives
  */
  ierr = PetscLogEventBegin(EVENT_INIT_GAUSS,0,0,0,0);CHKERRQ(ierr);
  ierr = fem_inigau(); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_INIT_GAUSS,0,0,0,0);CHKERRQ(ierr);

  /*
     Begin time dependent loop
  */
  int      i, j, nr_its = -1, kspits = -1;
  int      time_step = 0;
  double   t = t0;
  double   norm = -1.0, kspnorm = -1.0;
  double   strain_mac[6] = {0.1, 0.1, 0.2, 0.0, 0.0, 0.0};
  double   stress_mac[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  // Initial condition <x> = 0
  ierr = VecZeroEntries(x);CHKERRQ(ierr);

  if(flag_testcomm == TESTCOMM_STRAIN)
  {
    /*
       It is test. Sends a calculating strain to micro and obtain the stress
    */
    for(i=0;i<6;i++){
      memset(strain_mac,0.0,nvoi*sizeof(double));
      strain_mac[i] = 0.005;
      ierr = mac_send_signal(WORLD_COMM, MAC2MIC_STRAIN);CHKERRQ(ierr);
      ierr = mac_send_strain(WORLD_COMM, strain_mac);CHKERRQ(ierr);
      ierr = mac_recv_stress(WORLD_COMM, stress_mac);CHKERRQ(ierr);
      ierr = PetscPrintf(MACRO_COMM,"\nstress_ave = ");CHKERRQ(ierr);
      for(j=0;j<6;j++)
	ierr = PetscPrintf(MACRO_COMM,"%e ",stress_mac[j]);CHKERRQ(ierr);
      ierr = PetscPrintf(MACRO_COMM,"\n");CHKERRQ(ierr);
    }
  }
  else{
    /*
       It is real calculation
    */

    if(macmic.type==COUP_1){
      ierr = PetscPrintf(MACRO_COMM,"\ncalculating homo_cij\n", time_step, t);CHKERRQ(ierr);
      ierr = mac_calc_homo_cij(((mac_coup_1_t*)macmic.coup)->homo_cij);
      for(i=0;i<6;i++){
	for(j=0;j<6;j++){
	  ierr = PetscPrintf(MACRO_COMM,"%e ",((mac_coup_1_t*)macmic.coup)->homo_cij[i*6+j]);CHKERRQ(ierr);
	}ierr = PetscPrintf(MACRO_COMM,"\n");CHKERRQ(ierr);
      }ierr = PetscPrintf(MACRO_COMM,"\n");CHKERRQ(ierr);
    }

    while( t < (tf + 1.0e-10))
    {
      ierr = PetscPrintf(MACRO_COMM, "\nTime step %3d %e seg\n", time_step, t);CHKERRQ(ierr);

      /*
	 Setting Displacement on Dirichlet Indeces on <x>
       */
      ierr = PetscLogEventBegin(EVENT_SET_DISP_BOU,0,0,0,0);CHKERRQ(ierr);
      ierr = MacroSetDisplacementOnBoundary( t, &x);
      if( flag_print & (1<<PRINT_PETSC) ){
	ierr = PetscViewerASCIIOpen(MACRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
	ierr = VecView(x,viewer); CHKERRQ(ierr);
      }
      ierr = PetscLogEventEnd(EVENT_SET_DISP_BOU,0,0,0,0);CHKERRQ(ierr);

      /*
	 If the Residual Norm is bigger than <NormTol>
	 we should iterate
       */
      nr_its = 0; norm = 2*nr_norm_tol;
      while( nr_its < nr_max_its && norm > nr_norm_tol )
      {
	/*
	   Assemblying Residual
	 */
	ierr = PetscLogEventBegin(EVENT_ASSEMBLY_RES,0,0,0,0);CHKERRQ(ierr);
	ierr = PetscPrintf(MACRO_COMM, "Assembling Residual ");
	ierr = assembly_residual_sd( &x, &b);CHKERRQ(ierr);
	ierr = MacroSetBoundaryOnResidual( &b ); CHKERRQ(ierr);
	if( flag_print & (1<<PRINT_PETSC) ){
	  ierr = PetscViewerASCIIOpen(MACRO_COMM,"b.dat",&viewer); CHKERRQ(ierr);
	  ierr = VecView(b,viewer); CHKERRQ(ierr);
	}
	ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
	ierr = PetscPrintf(MACRO_COMM,"|b| = %e\n",norm);CHKERRQ(ierr);
	ierr = VecScale(b,-1.0); CHKERRQ(ierr);
	ierr = PetscLogEventEnd(EVENT_ASSEMBLY_RES,0,0,0,0);CHKERRQ(ierr);
	if( norm < nr_norm_tol )break;
	/*
	   Assemblying Jacobian
	 */
	ierr = PetscLogEventBegin(EVENT_ASSEMBLY_JAC,0,0,0,0);CHKERRQ(ierr);
	ierr = PetscPrintf(MACRO_COMM, "Assembling Jacobian\n");
	ierr = assembly_jacobian_sd(&A);
	ierr = MacroSetBoundaryOnJacobian( &A ); CHKERRQ(ierr);
	if( flag_print & (1<<PRINT_PETSC) ){
	  ierr = PetscViewerASCIIOpen(MACRO_COMM,"A.dat",&viewer); CHKERRQ(ierr);
	  ierr = MatView(A,viewer); CHKERRQ(ierr);
	}
	ierr = PetscLogEventEnd(EVENT_ASSEMBLY_JAC,0,0,0,0);CHKERRQ(ierr);
	/*
	   Solving Problem
	 */
	ierr = PetscLogEventBegin(EVENT_SOLVE_SYSTEM,0,0,0,0);CHKERRQ(ierr);
	ierr = PetscPrintf(MACRO_COMM, "Solving Linear System ");
	ierr = KSPSolve(ksp,b,dx);CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(ksp,&kspits);CHKERRQ(ierr);
	ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
	ierr = KSPGetResidualNorm(ksp,&kspnorm);CHKERRQ(ierr);
	ierr = VecAXPY( x, 1.0, dx); CHKERRQ(ierr);
	if( flag_print == PRINT_PETSC ){
	  ierr = PetscViewerASCIIOpen(MACRO_COMM,"dx.dat",&viewer); CHKERRQ(ierr);
	  ierr = VecView(dx,viewer); CHKERRQ(ierr);
	  ierr = PetscViewerASCIIOpen(MACRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
	  ierr = VecView(x,viewer); CHKERRQ(ierr);
	}

	switch(reason){
	  case KSP_CONVERGED_RTOL:
	    strcpy(reason_s, "RTOL");
	    break;
	  case KSP_CONVERGED_ATOL:
	    strcpy(reason_s, "ATOL");
	    break;
	  default :
	    strcpy(reason_s, "I DONT KNOW");
	    break;
	}
	ierr = PetscPrintf(MACRO_COMM,"Iterations %D Norm %e reason %s",kspits, kspnorm, reason_s);
	ierr = PetscLogEventEnd(EVENT_SOLVE_SYSTEM,0,0,0,0);CHKERRQ(ierr);

	nr_its ++;
      }

      if(flag_print & (1<<PRINT_VTK | 1<<PRINT_VTU)){ 
	strain = malloc(nelm*6*sizeof(double));
	stress = malloc(nelm*6*sizeof(double));
	energy = malloc(nelm*sizeof(double));
	ierr = assembly_residual_sd( &x, &b);CHKERRQ(ierr);
	ierr = calc_strain_stress_energy(&x, strain, stress, energy);
	if(flag_print & (1<<PRINT_VTK)){ 
	  sprintf(vtkfile_n,"%s_t_%d_%d.vtk",myname,time_step,rank_mac);
	  ierr = write_vtk(MACRO_COMM, vtkfile_n, &x, strain, stress);
	}
	if(flag_print & (1<<PRINT_VTU)){ 
	  sprintf(vtkfile_n,"%s_t_%d",myname,time_step);
	  ierr = write_vtu(MACRO_COMM, vtkfile_n, &x, &b, strain, stress, energy);
	}
	free(stress); free(strain); free(energy);
      }

      t += dt;
      time_step ++;
    }
  }

end_mac_1:
  /*
     Stop signal to micro if it is coupled
  */
  if(flag_coupling){
    ierr = mac_send_signal(WORLD_COMM, MIC_END); CHKERRQ(ierr);
  }


  /*
     Free Memory and close things
  */
  if(rank_mac==0) fclose(file_out); 

  list_clear(&material_list);
  list_clear(&physical_list);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);  
  ierr = VecDestroy(&b);CHKERRQ(ierr);  
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  PetscPrintf(MACRO_COMM,
      "--------------------------------------------------\n"
      "  MACRO: FINISH COMPLETE\n"
      "--------------------------------------------------\n");

  ierr = PetscFinalize();

end_mac_2:
  ierr = MPI_Finalize();

  return 0;
}
