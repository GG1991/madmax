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
"It has the capability of being couple with MICRO, a code for solving and RVE problem.n\n"
"-coupl    [0 (no coupling ) | 1 (coupling with micro)]\n"
"-testcomm [0 (no test) | 1 (sends a strain value and receive a stress calculated from micro)]\n"
"-print [0 (no print) | 1 (print PETSc structures) | 2 (print VTK output)]\n"
"-p_vtk [0 (no print vtk) | 1 (print partition) | 2 (print displacement,strain & stress)]\n";

#include "macro.h"

int main(int argc, char **argv)
{

  int        ierr;
  char       *myname = strdup("macro");
  PetscBool  set;

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
    printf("mac_main.c:no input file has been given\n");
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
  ierr = PetscOptionsGetInt(NULL, NULL, "-testcomm", &flag_testcomm, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE) flag_testcomm  = TESTCOMM_NULL;
  ierr = PetscOptionsGetString(NULL, NULL, "-mesh", mesh_n, 128, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE) SETERRQ(MICRO_COMM,1,"MICRO:mesh file not given on command line.");

  /* 
     Stablish a new local communicator
  */
  color = MACRO;
  macmic.type = COUP_NULL;
  ierr = MacMicParseScheme(input_n);
  ierr = MacMicColoring(WORLD_COMM, &color, &macmic, &MACRO_COMM);
  ierr = MPI_Comm_size(MACRO_COMM, &nproc_mac);
  ierr = MPI_Comm_rank(MACRO_COMM, &rank_mac);
  
  ierr = PetscFinalize();CHKERRQ(ierr);

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

  FileOutputStructures = NULL;
  if(rank_mac==0) FileOutputStructures = fopen("macro_structures.dat","w");

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
  PetscPrintf(MACRO_COMM,"MACRO: Reading mesh elements\n");
  strcpy(mesh_f,"gmsh");
  read_mesh_elmv(&MACRO_COMM, myname, mesh_n, mesh_f);
  ierr = PetscLogEventEnd(EVENT_READ_MESH_ELEM,0,0,0,0);CHKERRQ(ierr);

  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /*
     partition the mesh
  */
  PetscPrintf(MACRO_COMM,"MACRO: Partitioning and distributing mesh\n");
  ierr = PetscLogEventBegin(EVENT_PART_MESH,0,0,0,0);CHKERRQ(ierr);
  part_mesh_PARMETIS(&MACRO_COMM, time_fl, myname, NULL, PARMETIS_MESHKWAY );
  ierr = PetscLogEventEnd(EVENT_PART_MESH,0,0,0,0);CHKERRQ(ierr);

  /*
     Calculate <*ghosts> and <nghosts> 
  */
  PetscPrintf(MACRO_COMM,"MACRO: Calculating Ghost Nodes\n");
  ierr = PetscLogEventBegin(EVENT_CALC_GHOSTS,0,0,0,0);CHKERRQ(ierr);
  calculate_ghosts(&MACRO_COMM, myname);
  ierr = PetscLogEventEnd(EVENT_CALC_GHOSTS,0,0,0,0);CHKERRQ(ierr);

  /*
     Reenumerate Nodes
  */
  PetscPrintf(MACRO_COMM,"MACRO: Reenumering nodes\n");
  ierr = PetscLogEventBegin(EVENT_REENUMERATE,0,0,0,0);CHKERRQ(ierr);
  reenumerate_PETSc(&MACRO_COMM);
  ierr = PetscLogEventEnd(EVENT_REENUMERATE,0,0,0,0);CHKERRQ(ierr);

  /*
     Coordinate Reading
  */
  PetscPrintf(MACRO_COMM,"MACRO: Reading Coordinates\n");
  ierr = PetscLogEventBegin(EVENT_READ_COORD,0,0,0,0);CHKERRQ(ierr);
  read_mesh_coord(&MACRO_COMM, myname, mesh_n, mesh_f);
  ierr = PetscLogEventEnd(EVENT_READ_COORD,0,0,0,0);CHKERRQ(ierr);

  char  vtkfile_n[NBUF];

  sprintf(vtkfile_n,"%s_part_%d.vtk",myname,rank_mac);
  spu_vtk_partition( vtkfile_n, &MACRO_COMM );

  /*
     Read materials, physical entities, boundaries from input and mesh file
  */
  ierr = list_init(&physical_list, sizeof(physical_t), NULL); CHKERRQ(ierr);
  ierr = list_init(&function_list, sizeof(physical_t), NULL); CHKERRQ(ierr);
  ierr = SpuParseMaterials( &MACRO_COMM, input_n ); CHKERRQ(ierr);            
  ierr = SpuParsePhysicalEntities( &MACRO_COMM, mesh_n ); CHKERRQ(ierr);
  ierr = SpuParseFunctions( &MACRO_COMM, input_n ); CHKERRQ(ierr); 
  ierr = MacroParseBoundary(&MACRO_COMM, input_n ); CHKERRQ(ierr); 
  ierr = SetGmshIDOnMaterialsAndBoundaries(MACRO_COMM); CHKERRQ(ierr); 
  ierr = CheckPhysicalID(); CHKERRQ(ierr);
  ierr = SpuReadBoundary(MACRO_COMM, MACRO, mesh_n, mesh_f);CHKERRQ(ierr);
  ierr = MacroFillBoundary(MACRO_COMM);
  ierr = MacMicInitGaussStructure(eptr, nelm);CHKERRQ(ierr);

  /*
     Allocate matrices & vectors
  */ 
  ierr = PetscLogEventBegin(EVENT_ALLOC_MATVEC,0,0,0,0);CHKERRQ(ierr);
  ierr = PetscPrintf(MACRO_COMM, "allocating matrices & vectors\n");CHKERRQ(ierr);
  ierr = MacroAllocMatrixVector( MACRO_COMM, NMyNod*3, NTotalNod*3);CHKERRQ(ierr);
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
    ierr = MacroSetDisplacementOnBoundary( t, &x);
    if( flag_print == PRINT_PETSC ){
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
      ierr = MacroSetBoundaryOnResidual( &b ); CHKERRQ(ierr);
      if( flag_print == PRINT_PETSC ){
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
      ierr = MacroSetBoundaryOnJacobian( &A ); CHKERRQ(ierr);
      if( flag_print == PRINT_PETSC ){
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
      ierr = PetscPrintf(MACRO_COMM,"Iterations %D Norm %e reason %d\n",kspits, kspnorm, reason);CHKERRQ(ierr);
      ierr = PetscLogEventBegin(EVENT_SOLVE_SYSTEM,0,0,0,0);CHKERRQ(ierr);

      nr_its ++;
    }

    if(flag_print_vtk & FLAG_VTK_DISP){ 
      strain = malloc(nelm*6*sizeof(double));
      stress = malloc(nelm*6*sizeof(double));
      ierr = SpuCalcStressOnElement(&x, strain, stress);
      sprintf(vtkfile_n,"%s_displ_%d.vtk",myname,rank_mac);
      ierr = SpuVTKPlot_Displ_Strain_Stress(MACRO_COMM, vtkfile_n, &x, strain, stress);
      free(stress); free(strain);
    }

    t += dt;
    time_step ++;
  }

  /*
     Routine to send a calculating strain to micro and obtain the stress
  */
  if(flag_coupling && (flag_testcomm == TESTCOMM_STRAIN)){
    double strain[6] = {0.1, 0.1, 0.2, 0.0, 0.0, 0.0};
    double stress[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    ierr = MacCommSendSignal( WORLD_COMM, SIGNAL_SEND_STRAIN);CHKERRQ(ierr);
    ierr = MacCommSendStrain( WORLD_COMM, strain);CHKERRQ(ierr);
    ierr = MacCommRecvStress( WORLD_COMM, stress);CHKERRQ(ierr);
  }
  /*
     Stop signal to micro if it is coupled
  */
  if(flag_coupling){
    ierr = MacCommSendSignal( WORLD_COMM, SIGNAL_MICRO_END); CHKERRQ(ierr);
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
