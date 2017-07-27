/*
*
*   MACRO main function
*
*   Program for solving the displacement field inside a solid 
*   structure representing the macrostructure
*
*   Author: Guido Giuntoli
*
 */


static char help[] = "Solves the displacement field inside a solid structure. \
		      It has the capability of being couple with MICRO, a code for solving and RVE problem.n\n"
		      "-p_vtk [0 (no print vtk) | 1 (print partition) | 2 (print displacement,strain & stress)]\n";


#define   FLAG_VTK_NONE 0
#define   FLAG_VTK_PART 1
#define   FLAG_VTK_DISP 2

#include "macro.h"


int main(int argc, char **argv)
{

  int        i, ierr;
  int        flag_print_vtk = FLAG_VTK_NONE;
  char       *myname = strdup("macro");
  PetscBool  set;

  PetscLogEvent  CHECK_ERROR;    /* event number for error checking */
  PetscViewer    viewer;

#if defined(PETSC_USE_LOG)
  PetscLogStage stages[3];
  PetscLogEvent EVENT_READ_MESH_ELEM,
                EVENT_PART_MESH,
                EVENT_CALC_GHOSTS,
                EVENT_REENUMERATE,
                EVENT_READ_COORD,
		EVENT_INIT_GAUSS,
                EVENT_ALLOC_MATVEC,
                EVENT_SET_DISP_BOU,
                EVENT_ASSEMBLY_JAC,
                EVENT_ASSEMBLY_RES,
		EVENT_SOLVE_SYSTEM;
#endif

  world_comm = MPI_COMM_WORLD;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(world_comm, &nproc_wor);
  ierr = MPI_Comm_rank(world_comm, &rank_wor);

  if(argc>1){
    strcpy(input_n,argv[1]);
  }
  else{
    printf("mac_main.c:no input file has been given\n");
  }

  print_flag = false;
  for(i=2;i<argc;i++){
    if(strcmp(argv[i],"-p")==0){
      print_flag = true;
    }
  }
  //
  // File to print information about Mesh:
  // -> number of elements (before and after partition)
  // -> number of nodes for each <boundary_t>
  //

  spu_parse_scheme( input_n );

  /* 
     Stablish a new local communicator and a set of 
     intercommunicators with micro programs 
  */
  mac_comm_init();
  

  // here we collect time from all processes
  if(rank_mac == 0){
    time_vec = calloc(nproc_mac, sizeof(double));
  }

  spu_parse_mesh(input_n);

  /*
     Set PETSc communicator to MACRO_COMM
  */
  PETSC_COMM_WORLD = MACRO_COMM;
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);
  PetscPrintf(MACRO_COMM,
      "--------------------------------------------------\n"
      "  MACRO: COMPOSITE MATERIAL MULTISCALE CODE\n"
      "--------------------------------------------------\n");

  FileOutputStructures = NULL;
  if(rank_mac==0) FileOutputStructures = fopen("macro_structures.dat","w");

#if defined(PETSC_USE_LOG)
  ierr = PetscLogEventRegister("Read Elems of Mesh"    ,PETSC_VIEWER_CLASSID,&EVENT_READ_MESH_ELEM);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Partition Mesh"        ,PETSC_VIEWER_CLASSID,&EVENT_PART_MESH);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Calculate Ghosts Nodes",PETSC_VIEWER_CLASSID,&EVENT_CALC_GHOSTS);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Reenumerates Nodes"    ,PETSC_VIEWER_CLASSID,&EVENT_REENUMERATE);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Read Coordinates"      ,PETSC_VIEWER_CLASSID,&EVENT_READ_COORD);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Init Gauss Points"     ,PETSC_VIEWER_CLASSID,&EVENT_INIT_GAUSS);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Allocate Mat and Vec"  ,PETSC_VIEWER_CLASSID,&EVENT_ALLOC_MATVEC);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Set Displ on Bou "     ,PETSC_VIEWER_CLASSID,&EVENT_SET_DISP_BOU);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Assembly Jacobian"     ,PETSC_VIEWER_CLASSID,&EVENT_ASSEMBLY_JAC);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Assembly Residual"     ,PETSC_VIEWER_CLASSID,&EVENT_ASSEMBLY_RES);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Solve Linear System"   ,PETSC_VIEWER_CLASSID,&EVENT_SOLVE_SYSTEM);CHKERRQ(ierr);
#endif

  /*
     Get command line arguments
  */
  ierr = PetscOptionsGetInt(NULL, NULL, "-p_vtk", &flag_print_vtk, &set); CHKERRQ(ierr); 

  //
  // Register various stages for profiling
  //
  ierr = PetscLogStageRegister("Read Mesh Elements",&stages[0]);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Linear System 1",&stages[1]);CHKERRQ(ierr);
  ierr = PetscLogStageRegister("Linear System 2",&stages[2]);CHKERRQ(ierr);

  //
  // Register a user-defined event for profiling (error checking).
  //
  CHECK_ERROR = 0;
  ierr        = PetscLogEventRegister("Check Error",KSP_CLASSID,&CHECK_ERROR);CHKERRQ(ierr);

  ierr = PetscLogStagePush(stages[0]);CHKERRQ(ierr);

  //
  // read mesh
  //    
  ierr = PetscLogEventBegin(EVENT_READ_MESH_ELEM,0,0,0,0);CHKERRQ(ierr);
  PetscPrintf(MACRO_COMM,"MACRO: Reading mesh elements\n");
  strcpy(mesh_f,"gmsh");
  read_mesh_elmv(&MACRO_COMM, myname, mesh_n, mesh_f);
  ierr = PetscLogEventEnd(EVENT_READ_MESH_ELEM,0,0,0,0);CHKERRQ(ierr);

  //
  // End curent profiling stage
  //
  ierr = PetscLogStagePop();CHKERRQ(ierr);

  nelm = elmdist[rank_mac+1] - elmdist[rank_mac];
  part = (int*)malloc(nelm * sizeof(int));

  //
  // partition the mesh
  //
  PetscPrintf(MACRO_COMM,"MACRO: Partitioning and distributing mesh\n");
  ierr = PetscLogEventBegin(EVENT_PART_MESH,0,0,0,0);CHKERRQ(ierr);
  part_mesh_PARMETIS(&MACRO_COMM, time_fl, myname, NULL, PARMETIS_MESHKWAY );
  ierr = PetscLogEventEnd(EVENT_PART_MESH,0,0,0,0);CHKERRQ(ierr);
  //
  // We delete repeated nodes and save the <NAllMyNod> values on <AllMyNodOrig> in order
  // this should be inside part_mesh_PARMETIS ?
  //
  clean_vector_qsort(&MACRO_COMM, myname, eptr[nelm], eind, &AllMyNodOrig, &NAllMyNod);
  //
  // Calculate <*ghosts> and <nghosts> 
  //
  PetscPrintf(MACRO_COMM,"MACRO: Calculating Ghost Nodes\n");
  ierr = PetscLogEventBegin(EVENT_CALC_GHOSTS,0,0,0,0);CHKERRQ(ierr);
  calculate_ghosts(&MACRO_COMM, myname);
  ierr = PetscLogEventEnd(EVENT_CALC_GHOSTS,0,0,0,0);CHKERRQ(ierr);
  //
  // Reenumerate Nodes
  //
  PetscPrintf(MACRO_COMM,"MACRO: Reenumering nodes\n");
  ierr = PetscLogEventBegin(EVENT_REENUMERATE,0,0,0,0);CHKERRQ(ierr);
  reenumerate_PETSc(&MACRO_COMM);
  ierr = PetscLogEventEnd(EVENT_REENUMERATE,0,0,0,0);CHKERRQ(ierr);
  //
  // Coordinate Reading
  //
  PetscPrintf(MACRO_COMM,"MACRO: Reading Coordinates\n");
  ierr = PetscLogEventBegin(EVENT_READ_COORD,0,0,0,0);CHKERRQ(ierr);
  read_mesh_coord(&MACRO_COMM, myname, mesh_n, mesh_f);
  ierr = PetscLogEventEnd(EVENT_READ_COORD,0,0,0,0);CHKERRQ(ierr);

  char  vtkfile_n[NBUF];

  sprintf(vtkfile_n,"%s_part_%d.vtk",myname,rank_mac);
  spu_vtk_partition( vtkfile_n, &MACRO_COMM );
  /*
     read materials and physical entities from input and mehs files
  */
  ierr = list_init(&physical_list, sizeof(physical_t), NULL); CHKERRQ(ierr);
  ierr = list_init(&function_list, sizeof(physical_t), NULL); CHKERRQ(ierr);
  ierr = SpuParseMaterials( &MACRO_COMM, input_n ); CHKERRQ(ierr);            
  ierr = SpuParsePhysicalEntities( &MACRO_COMM, mesh_n ); CHKERRQ(ierr);
  ierr = SpuParseFunctions( &MACRO_COMM, input_n ); CHKERRQ(ierr); 
  ierr = SpuParseBoundary(&MACRO_COMM, input_n ); CHKERRQ(ierr); 
  ierr = SetGmshIDOnMaterialsAndBoundaries(MACRO_COMM); CHKERRQ(ierr); 
  ierr = CheckPhysicalID(); CHKERRQ(ierr);
  ierr = SpuReadBoundary(MACRO_COMM, mesh_n, mesh_f, FileOutputStructures );CHKERRQ(ierr);

  /*
     Allocate matrices & vectors
  */ 
  ierr = PetscLogEventBegin(EVENT_ALLOC_MATVEC,0,0,0,0);CHKERRQ(ierr);
  PetscPrintf(MACRO_COMM, "allocating matrices & vectors\n");
  MacroAllocMatrixVector( MACRO_COMM, NMyNod*3, NTotalNod*3);
  ierr = PetscLogEventEnd(EVENT_ALLOC_MATVEC,0,0,0,0);CHKERRQ(ierr);

  /*
     Setting solver options 
  */
  ierr = KSPCreate(MACRO_COMM,&ksp); CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPCG); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A); CHKERRQ(ierr);

  ierr = PetscLogEventBegin(EVENT_INIT_GAUSS,0,0,0,0);CHKERRQ(ierr);
  fem_inigau();
  ierr = PetscLogEventEnd(EVENT_INIT_GAUSS,0,0,0,0);CHKERRQ(ierr);

  /*
     Begin time dependent loop
  */
  int    nr_its = -1, kspits = -1;
  int    time_step = 0;
  double t0 = 0.0, tf = 1.0, dt = 1.0, t = t0;
  double norm = -1.0, NormTol = 1.0e-8, NRMaxIts = 3, kspnorm = -1.0;


  // Initial condition <x> = 0
  ierr = VecZeroEntries(x);CHKERRQ(ierr);

  while( t < (tf + 1.0e-10))
  {
    ierr = PetscPrintf(MACRO_COMM, "\nTime step %3d %e seg\n", time_step, t);CHKERRQ(ierr);

    /*
       Setting Displacement on Dirichlet Indeces on <x>
    */
    ierr = PetscLogEventBegin(EVENT_SET_DISP_BOU,0,0,0,0);CHKERRQ(ierr);
    ierr = SputnikSetDisplacementOnBoundary( t, &x);
    if(print_flag){
      ierr = PetscViewerASCIIOpen(MACRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(x,viewer); CHKERRQ(ierr);
    }
    ierr = PetscLogEventEnd(EVENT_SET_DISP_BOU,0,0,0,0);CHKERRQ(ierr);

    /*
       If the Residual Norm is bigger than <NormTol>
       we should iterate
    */

    nr_its = 0; norm = 2*NormTol;
    while( nr_its < NRMaxIts && norm > NormTol )
    {

      /*
	 Assemblying Residual
       */
      ierr = PetscLogEventBegin(EVENT_ASSEMBLY_RES,0,0,0,0);CHKERRQ(ierr);
      ierr = PetscPrintf(MACRO_COMM, "Assembling Residual ");CHKERRQ(ierr);
      ierr = AssemblyResidualSmallDeformation( &x, &b);CHKERRQ(ierr);
      ierr = SputnikSetBoundaryOnResidual( &b ); CHKERRQ(ierr);
      if(print_flag){
	ierr = PetscViewerASCIIOpen(MACRO_COMM,"b.dat",&viewer); CHKERRQ(ierr);
	ierr = VecView(b,viewer); CHKERRQ(ierr);
      }
      ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
      ierr = PetscPrintf(MACRO_COMM,"|b| = %e\n",norm);CHKERRQ(ierr);
      ierr = VecScale(b,-1.0); CHKERRQ(ierr);
      ierr = PetscLogEventEnd(EVENT_ASSEMBLY_RES,0,0,0,0);CHKERRQ(ierr);
      if( !(norm > NormTol) )break;
      /*
	 Assemblying Jacobian
       */
      ierr = PetscLogEventBegin(EVENT_ASSEMBLY_JAC,0,0,0,0);CHKERRQ(ierr);
      ierr = PetscPrintf(MACRO_COMM, "Assembling Jacobian\n");
      ierr = AssemblyJacobianSmallDeformation(&A);
      ierr = SputnikSetBoundaryOnJacobian( &A ); CHKERRQ(ierr);
      if(print_flag){
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
      if(print_flag){
	ierr = PetscViewerASCIIOpen(MACRO_COMM,"dx.dat",&viewer); CHKERRQ(ierr);
	ierr = VecView(dx,viewer); CHKERRQ(ierr);
	ierr = PetscViewerASCIIOpen(MACRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
	ierr = VecView(x,viewer); CHKERRQ(ierr);
      }
      ierr = PetscPrintf(MACRO_COMM,"Iterations %D Norm %e reason %d\n",kspits, kspnorm, reason);CHKERRQ(ierr);
      ierr = PetscLogEventBegin(EVENT_SOLVE_SYSTEM,0,0,0,0);CHKERRQ(ierr);

      nr_its ++;
    }

    if(flag_print_vtk & FLAG_VTK_DISP){ 
      sprintf(vtkfile_n,"%s_displ_%d.vtk",myname,rank_mac);
      SpuVTKPlot_Displ_Strain_Stress(MACRO_COMM, vtkfile_n, &x, &Strain, &Stress);
    }

    t += dt;
    time_step ++;
  }

  /*
     Free Memory and close things
  */
  if(rank_mac==0) fclose(FileOutputStructures); 

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
  ierr = MPI_Finalize();

  return 0;
}
