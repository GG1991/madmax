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
		      It has the capability of being couple with MICRO, a code for solving and RVE problem.n\n";


#include "macro.h"


int main(int argc, char **argv)
{

  int        i, ierr;
  char       *myname = strdup("macro");

  PetscLogEvent  CHECK_ERROR;    /* event number for error checking */
  PetscViewer    viewer,viewer1;

#if defined(PETSC_USE_LOG)
  PetscLogStage stages[3];
  PetscLogEvent EVENT_READ_MESH_ELEM,
                EVENT_PART_MESH,
                EVENT_CALC_GHOSTS,
                EVENT_REENUMERATE,
                EVENT_READ_COORD,
		EVENT_INIT_GAUSS,
                EVENT_ALLOC_MATVEC,
                EVENT_ASSEMBLY_JAC,
                EVENT_ASSEMBLY_RES;
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

  // 
  // Stablish a new local communicator and a set of 
  // intercommunicators with micro programs 
  //
  mac_comm_init();
  

  // here we collect time from all processes
  if(rank_mac == 0){
    time_vec = calloc(nproc_mac, sizeof(double));
  }

  spu_parse_mesh(input_n);

  //
  // Set PETSc communicator to MACRO_COMM
  //
  PETSC_COMM_WORLD = MACRO_COMM;
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);
  PetscPrintf(MACRO_COMM,"--------------------------------------------------\n"
      "  MACRO: COMPOSITE MATERIAL MULTISCALE CODE\n"
      "--------------------------------------------------\n");

  FileOutputStructures = NULL;
  if(rank_mac==0) FileOutputStructures = fopen("macro_structures.dat","w");

#if defined(PETSC_USE_LOG)
  ierr = PetscLogEventRegister("read_mesh_elmv     Event",PETSC_VIEWER_CLASSID,&EVENT_READ_MESH_ELEM);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("part_mesh_PARMETIS Event",PETSC_VIEWER_CLASSID,&EVENT_PART_MESH);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("calculate_ghosts   Event",PETSC_VIEWER_CLASSID,&EVENT_CALC_GHOSTS);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("reenumerate_PETSc  Event",PETSC_VIEWER_CLASSID,&EVENT_REENUMERATE);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("read_mesh_coord    Event",PETSC_VIEWER_CLASSID,&EVENT_READ_COORD);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("fem_inigau         Event",PETSC_VIEWER_CLASSID,&EVENT_INIT_GAUSS);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("AllocMatrixVector  Event",PETSC_VIEWER_CLASSID,&EVENT_ALLOC_MATVEC);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Assembly Jacobian  Event",PETSC_VIEWER_CLASSID,&EVENT_ASSEMBLY_JAC);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("Assembly Residual  Event",PETSC_VIEWER_CLASSID,&EVENT_ASSEMBLY_RES);CHKERRQ(ierr);
#endif

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
  //
  // read materials and physical entities from input and mehs files
  //
  ierr = list_init(&physical_list, sizeof(physical_t), NULL); CHKERRQ(ierr);
  ierr = list_init(&function_list, sizeof(physical_t), NULL); CHKERRQ(ierr);
  ierr = SpuParseMaterials( &MACRO_COMM, input_n ); CHKERRQ(ierr);            
  ierr = SpuParsePhysicalEntities( &MACRO_COMM, mesh_n ); CHKERRQ(ierr);
  ierr = SpuParseFunctions( &MACRO_COMM, input_n ); CHKERRQ(ierr); 
  ierr = SpuParseBoundary(&MACRO_COMM, input_n ); CHKERRQ(ierr); 
  ierr = SetGmshIDOnMaterialsAndBoundaries(); CHKERRQ(ierr); 
  ierr = CheckPhysicalID(); CHKERRQ(ierr);
  ierr = SpuReadBoundary( &MACRO_COMM, mesh_n, mesh_f, FileOutputStructures );CHKERRQ(ierr);
  //
  // Allocate matrices & vectors
  // 
  ierr = PetscLogEventBegin(EVENT_ALLOC_MATVEC,0,0,0,0);CHKERRQ(ierr);
  PetscPrintf(MACRO_COMM, "allocating matrices & vectors\n");
  AllocMatrixVector( MACRO_COMM, NMyNod*3, NTotalNod*3, &A, &x, &b);
  //
  // Currently, all PETSc parallel matrix formats are partitioned by
  // contiguous chunks of rows across the processors.  Determine which
  // rows of the matrix are locally owned.
  //
  int Istart, Iend;
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  if( Istart != StartIndexRank[rank_mac]*3 ){
    printf("mac_main: error on indeces set for matrix and vector.\n");
    return 1;
  }
  if(rank_mac<nproc_mac-1){
    if( Iend != StartIndexRank[rank_mac+1]*3 ){
      printf("mac_main: error on indeces set for matrix and vector.\n");
      return 1;
    }
  }
  else{
    if( Iend != NTotalNod*3 ){
      printf("mac_main: error on indeces set for matrix and vector.\n");
      return 1;
    }
  }
  ierr = PetscLogEventEnd(EVENT_ALLOC_MATVEC,0,0,0,0);CHKERRQ(ierr);

  //
  // Setting solver options 
  //
  ierr = KSPCreate(MACRO_COMM,&ksp); CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPCG); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A); CHKERRQ(ierr);

  ierr = PetscLogEventBegin(EVENT_INIT_GAUSS,0,0,0,0);CHKERRQ(ierr);
  fem_inigau();
  ierr = PetscLogEventEnd(EVENT_INIT_GAUSS,0,0,0,0);CHKERRQ(ierr);

  double time;

//  ierr = SetBoundaryConditionsOnX( &x, time);

  int KspIterationNum;
  //
  // Assemblying Jacobian
  //
  ierr = PetscLogEventBegin(EVENT_ASSEMBLY_JAC,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscPrintf(MACRO_COMM, "Assembling Jacobian\n");
  ierr = AssemblyJacobianSmallDeformation(&A);
  ierr = PetscLogEventEnd(EVENT_ASSEMBLY_JAC,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(MACRO_COMM,"A.dat",&viewer); CHKERRQ(ierr);
  ierr = MatView(A,viewer); CHKERRQ(ierr);
  //
  // Assemblying Residual
  //
  ierr = PetscLogEventBegin(EVENT_ASSEMBLY_RES,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscPrintf(MACRO_COMM, "Assembling Residual\n");CHKERRQ(ierr);
  ierr = AssemblyResidualSmallDeformation( &x, &b);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_ASSEMBLY_RES,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(MACRO_COMM,"b.dat",&viewer1); CHKERRQ(ierr);
  ierr = VecView(b,viewer1); CHKERRQ(ierr);
  //
  // Solving Problem
  //
  int its;	
  double norm;
  ierr = PetscPrintf(MACRO_COMM, "Solving Linear System\n");
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
  ierr = PetscPrintf(MACRO_COMM,"Norm of error %g Iterations %D reason %d\n",norm,its,reason);CHKERRQ(ierr);
  ierr = PetscPrintf(MACRO_COMM, "OOKK !\n");

  //
  // Free Memory and close things
  //
  if(rank_mac==0) fclose(FileOutputStructures); 

  list_clear(&material_list);
  list_clear(&physical_list);

  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);  
  ierr = VecDestroy(&b);CHKERRQ(ierr);  
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  ierr = PetscFinalize();
  ierr = MPI_Finalize();

  return 0;
}
