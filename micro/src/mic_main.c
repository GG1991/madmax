/*

   MICRO main function

   Program for solving the microscopic problem 
   for multi-scale approach.

   Author: Guido Giuntoli
   Date: 28-07-2017

 */


static char help[] = "Solves an RVE problem.\n\n";



#include "micro.h"



int main(int argc, char **argv)
{

  int        i, ierr;
  char       *myname = strdup("micro");
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

  if(argc>1){
    strcpy(input_n,argv[1]);
  }
  else{
    printf("mic_main.c:no input file has been given\n");
  }

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
     Get command line arguments
  */

  ierr = PetscOptionsGetInt(NULL, NULL, "-p_vtk", &flag_print_vtk, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE) flag_print_vtk = FLAG_VTK_NONE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-coupl", &flag_coupling, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE) flag_coupling  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL, NULL, "-p", &print_flag, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE) print_flag  = PETSC_FALSE;

  /* 
     Stablish a new local communicator
   */
  color = MICRO;
  macmic.type = COUP_NULL;
  ierr = MacMicParseScheme(input_n);
  ierr = MacMicColoring(WORLD_COMM, &color, &macmic, &MICRO_COMM); /* color can change */
  ierr = MPI_Comm_size(MICRO_COMM, &nproc_mic);
  ierr = MPI_Comm_rank(MICRO_COMM, &rank_mic);

  ierr = spu_parse_mesh(input_n);

  ierr = PetscFinalize();CHKERRQ(ierr);

  /*
     Set PETSc communicator to MICRO_COMM
     and start again
   */
  PETSC_COMM_WORLD = MICRO_COMM;
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);
  
  if(flag_coupling){
    PetscPrintf(MICRO_COMM,
	"--------------------------------------------------\n"
	"  MICRO: COUPLING \n"
	"--------------------------------------------------\n");
  }
  else{
    PetscPrintf(MICRO_COMM,
	"--------------------------------------------------\n"
	"  MICRO: STANDALONE \n"
	"--------------------------------------------------\n");
  }

  FileOutputStructures = NULL;
  if(rank_mic==0) FileOutputStructures = fopen("micro_structures.dat","w");

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
  strcpy(mesh_f,"gmsh");
  read_mesh_elmv(&MICRO_COMM, myname, mesh_n, mesh_f);
  ierr = PetscLogEventEnd(EVENT_READ_MESH_ELEM,0,0,0,0);CHKERRQ(ierr);

  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /*
     partition the mesh
  */
  ierr = PetscLogEventBegin(EVENT_PART_MESH,0,0,0,0);CHKERRQ(ierr);
  part_mesh_PARMETIS(&MICRO_COMM, time_fl, myname, NULL, PARMETIS_MESHKWAY );
  ierr = PetscLogEventEnd(EVENT_PART_MESH,0,0,0,0);CHKERRQ(ierr);

  /*
     Calculate <*ghosts> and <nghosts> 
  */
  ierr = PetscLogEventBegin(EVENT_CALC_GHOSTS,0,0,0,0);CHKERRQ(ierr);
  calculate_ghosts(&MICRO_COMM, myname);
  ierr = PetscLogEventEnd(EVENT_CALC_GHOSTS,0,0,0,0);CHKERRQ(ierr);

  /*
     Reenumerate Nodes
  */
  ierr = PetscLogEventBegin(EVENT_REENUMERATE,0,0,0,0);CHKERRQ(ierr);
  reenumerate_PETSc(&MICRO_COMM);
  ierr = PetscLogEventEnd(EVENT_REENUMERATE,0,0,0,0);CHKERRQ(ierr);

  /*
     Coordinate Reading
  */
  ierr = PetscLogEventBegin(EVENT_READ_COORD,0,0,0,0);CHKERRQ(ierr);
  read_mesh_coord(&MICRO_COMM, myname, mesh_n, mesh_f);
  ierr = PetscLogEventEnd(EVENT_READ_COORD,0,0,0,0);CHKERRQ(ierr);

  char  vtkfile_n[NBUF];

  sprintf(vtkfile_n,"%s_part_%d.vtk",myname,rank_mic);
  spu_vtk_partition( vtkfile_n, &MICRO_COMM );

  /*
     read materials and physical entities from input and mehs files
  */
  ierr = list_init(&physical_list, sizeof(physical_t), NULL); CHKERRQ(ierr);
  ierr = list_init(&function_list, sizeof(physical_t), NULL); CHKERRQ(ierr);
  ierr = SpuParseMaterials( &MICRO_COMM, input_n ); CHKERRQ(ierr);            
  ierr = SpuParsePhysicalEntities( &MICRO_COMM, mesh_n ); CHKERRQ(ierr);
  ierr = SpuParseFunctions( &MICRO_COMM, input_n ); CHKERRQ(ierr); 
  ierr = SpuParseBoundary(&MICRO_COMM, input_n ); CHKERRQ(ierr); 
  ierr = SetGmshIDOnMaterialsAndBoundaries(MICRO_COMM); CHKERRQ(ierr); 
  ierr = CheckPhysicalID(); CHKERRQ(ierr);
  ierr = SpuReadBoundary(MICRO_COMM, mesh_n, mesh_f, FileOutputStructures );CHKERRQ(ierr);
  ierr = MacMicInitGaussStructure(eptr, nelm);CHKERRQ(ierr);

  /*
     Allocate matrices & vectors
  */ 
  ierr = PetscLogEventBegin(EVENT_ALLOC_MATVEC,0,0,0,0);CHKERRQ(ierr);
  ierr = MicroAllocMatrixVector( MICRO_COMM, NMyNod*3, NTotalNod*3);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_ALLOC_MATVEC,0,0,0,0);CHKERRQ(ierr);

  /*
     Setting solver options 
  */
  ierr = KSPCreate(MICRO_COMM,&ksp); CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPCG); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A); CHKERRQ(ierr);

  /*
     Init Gauss point shapes functions and derivatives
  */
  ierr = PetscLogEventBegin(EVENT_INIT_GAUSS,0,0,0,0);CHKERRQ(ierr);
  ierr = fem_inigau(); CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_INIT_GAUSS,0,0,0,0);CHKERRQ(ierr);

  if(flag_coupling){
    ierr = MicCommWaitStartSignal( WORLD_COMM );CHKERRQ(ierr);
  }

  ierr = PetscFinalize();
  ierr = MPI_Finalize();

  return 0;
}
