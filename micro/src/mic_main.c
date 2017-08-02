/*

   MICRO main function

   Program for solving the microscopic problem 
   for multi-scale approach.

   Author> Guido Giuntoli
   Date> 28-07-2017

 */

static char help[] = 
"MICRO MULTISCALE CODE\n"
"Solves the RVE problem inside a solid structure. \n"
"It has the capability of being couple with MACRO.\n"
"-coupl    [0 (no coupling ) | 1 (coupling with micro)]\n"
"-testcomm [0 (no test) | 1 (sends a strain value and receive a stress calculated from micro)]\n"
"-print [0 (no print) | 1 (print PETSc structures) | 2 (print VTK output)]\n"
"-p_vtk [0 (no print vtk) | 1 (print partition) | 2 (print displacement,strain & stress)]\n";

#include "micro.h"

int main(int argc, char **argv)
{

  int    ierr;
  char   *myname = strdup("micro");
  char   vtkfile_n[NBUF];
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
  flag_print = PETSC_FALSE;
  ierr = PetscOptionsHasName(NULL,NULL,"-print_petsc",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_PETSC);
  ierr = PetscOptionsHasName(NULL,NULL,"-print_disp",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_VTKDISP);
  ierr = PetscOptionsHasName(NULL,NULL,"-print_part",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_VTKPART);
  ierr = PetscOptionsHasName(NULL,NULL,"-print_all",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_ALL);
  ierr = PetscOptionsGetBool(NULL, NULL, "-coupl", &flag_coupling, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE) flag_coupling  = PETSC_FALSE;
  ierr = PetscOptionsGetString(NULL, NULL, "-mesh", mesh_n, 128, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE) SETERRQ(MICRO_COMM,1,"MICRO:mesh file not given on command line.");
  ierr = PetscOptionsGetString(NULL, NULL, "-input", input_n, 128, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE) SETERRQ(MICRO_COMM,1,"MICRO:input file not given.");

  /* 
     Stablish a new local communicator
   */
  color = MICRO;
  macmic.type = COUP_NULL;
  ierr = MacMicParseScheme(input_n);
  ierr = MacMicColoring(WORLD_COMM, &color, &macmic, &MICRO_COMM); /* color can change */
  ierr = MPI_Comm_size(MICRO_COMM, &nproc_mic);
  ierr = MPI_Comm_rank(MICRO_COMM, &rank_mic);

  ierr = PetscFinalize();CHKERRQ(ierr);

  /*
     Set PETSc communicator to MICRO_COMM
     and start again
   */
  PETSC_COMM_WORLD = MICRO_COMM;
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);
  
  if(!flag_coupling){
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
  if(!flag_coupling)
    PetscPrintf(MICRO_COMM,"MACRO: Reading mesh elements\n");
  strcpy(mesh_f,"gmsh");
  read_mesh_elmv(&MICRO_COMM, myname, mesh_n, mesh_f);
  ierr = PetscLogEventEnd(EVENT_READ_MESH_ELEM,0,0,0,0);CHKERRQ(ierr);

  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /*
     partition the mesh
  */
  if(!flag_coupling)
    PetscPrintf(MICRO_COMM,"MICRO: Partitioning and distributing mesh\n");
  ierr = PetscLogEventBegin(EVENT_PART_MESH,0,0,0,0);CHKERRQ(ierr);
  part_mesh_PARMETIS(&MICRO_COMM, time_fl, myname, NULL, PARMETIS_MESHKWAY );
  ierr = PetscLogEventEnd(EVENT_PART_MESH,0,0,0,0);CHKERRQ(ierr);

  /*
     Calculate <*ghosts> and <nghosts> 
  */
  if(!flag_coupling)
    PetscPrintf(MICRO_COMM,"MICRO: Calculating Ghost Nodes\n");
  ierr = PetscLogEventBegin(EVENT_CALC_GHOSTS,0,0,0,0);CHKERRQ(ierr);
  calculate_ghosts(&MICRO_COMM, myname);
  ierr = PetscLogEventEnd(EVENT_CALC_GHOSTS,0,0,0,0);CHKERRQ(ierr);

  /*
     Reenumerate Nodes
  */
  if(!flag_coupling)
    PetscPrintf(MICRO_COMM,"MICRO: Reenumering nodes\n");
  ierr = PetscLogEventBegin(EVENT_REENUMERATE,0,0,0,0);CHKERRQ(ierr);
  reenumerate_PETSc(&MICRO_COMM);
  ierr = PetscLogEventEnd(EVENT_REENUMERATE,0,0,0,0);CHKERRQ(ierr);

  /*
     Coordinate Reading
  */
  if(!flag_coupling)
    PetscPrintf(MICRO_COMM,"MICRO: Reading Coordinates\n");
  ierr = PetscLogEventBegin(EVENT_READ_COORD,0,0,0,0);CHKERRQ(ierr);
  read_mesh_coord(&MICRO_COMM, myname, mesh_n, mesh_f);
  ierr = PetscLogEventEnd(EVENT_READ_COORD,0,0,0,0);CHKERRQ(ierr);

  if(flag_print | (1<<PRINT_VTKPART)){
    sprintf(vtkfile_n,"%s_part_%d.vtk",myname,rank_mic);
    ierr = spu_vtk_partition( vtkfile_n, &MICRO_COMM );
  }

  /*
     Read materials, physical entities, boundaries from input and mesh file
  */
  ierr = list_init(&physical_list, sizeof(physical_t), NULL); CHKERRQ(ierr);
  ierr = list_init(&function_list, sizeof(physical_t), NULL); CHKERRQ(ierr);
  ierr = SpuParseMaterials( &MICRO_COMM, input_n ); CHKERRQ(ierr);            
  ierr = SpuParsePhysicalEntities( &MICRO_COMM, mesh_n ); CHKERRQ(ierr);
  ierr = SpuParseFunctions( &MICRO_COMM, input_n ); CHKERRQ(ierr); 

  /*
     Check if Physical Entities <P000 P100 P010 X0 X1 Y0 Y1 Z0 Z1> exist in mesh
     Creates the boundary_list with <P000 P100 P010 X0 X1 Y0 Y1 Z0 Z1>
  */
  ierr = MicroCheckPhysicalEntities(&physical_list);CHKERRQ(ierr);
  ierr = MicroCreateBoundary(&boundary_list);CHKERRQ(ierr);
  ierr = SetGmshIDOnMaterialsAndBoundaries(MICRO_COMM); CHKERRQ(ierr); 
  ierr = CheckPhysicalID(); CHKERRQ(ierr);
  ierr = SpuReadBoundary(MICRO_COMM, mesh_n, mesh_f);CHKERRQ(ierr);
  ierr = MicroSetBoundary(MICRO_COMM, &boundary_list );CHKERRQ(ierr);
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

  /*
     micro main coupling loop
   */
  double strain_bc[6], stress_ave[6], ttensor[36];
  LX = LY = LZ = 1.0;

  if(flag_coupling){
    /*
       Coupling
     */
    int signal = SIGNAL_NULL;
    while(signal != SIGNAL_MICRO_END ){
      ierr = MicCommWaitSignal( WORLD_COMM, &signal );CHKERRQ(ierr);
      switch( signal ){
	case SIGNAL_SEND_STRAIN:
	  /*
	     Wait for strain
	  */
	  ierr = MicCommRecvStrain( WORLD_COMM, strain );
	  /*
	     Performs the micro calculation
	  */
	  ierr = MicroSetDisplacementOnBoundary( 0, strain[0], LX, LY, LZ, &x );
	  if( flag_print | (1<<PRINT_PETSC) ){
	    ierr = PetscViewerASCIIOpen(MICRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
	    ierr = VecView(x,viewer); CHKERRQ(ierr);
	  }
	  stress_ave[0] = 1.0; stress_ave[1] = 1.0; stress_ave[2] = 1.0; stress_ave[3] = -1.0; stress_ave[4] = 1.0; stress_ave[5] = 1.0;
	  /*
	     Send Stress
	  */
	  ierr = MicCommSendStress( WORLD_COMM, stress_ave );
	  break;
	case SIGNAL_MICRO_END:
	  break;
	default:
	  SETERRQ(MICRO_COMM,1,"MICRO:can no identify recv signal.");
      }
    }
  }
  else{
    /*
       Standalone
       Performs 6 experiment
     */
    int    nr_its = -1, kspits = -1;
    int    dir = 0;
    double norm = -1.0, NormTol = 1.0e-8, NRMaxIts = 3, kspnorm = -1.0;

    strain_bc[0] = 0.005; strain_bc[1] = 0.005; strain_bc[2] = 0.005; strain_bc[3] = 0.005; strain_bc[4] = 0.005; strain_bc[5] = 0.005;
    ierr = PetscLogEventBegin(EVENT_SET_DISP_BOU,0,0,0,0);CHKERRQ(ierr);
    ierr = MicroSetDisplacementOnBoundary( 0, strain_bc[0], LX, LY, LZ, &x );
    if( flag_print | (1<<PRINT_PETSC) ){
      ierr = PetscViewerASCIIOpen(MICRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
      ierr = VecView(x,viewer); CHKERRQ(ierr);
    }
    ierr = PetscLogEventEnd(EVENT_SET_DISP_BOU,0,0,0,0);CHKERRQ(ierr);

    nr_its = 0; norm = 2*NormTol;
    while( nr_its < NRMaxIts && norm > NormTol )
    {
      /*
	 Assemblying Residual
       */
      ierr = PetscLogEventBegin(EVENT_ASSEMBLY_RES,0,0,0,0);CHKERRQ(ierr);
      ierr = PetscPrintf(MICRO_COMM, "Assembling Residual ");CHKERRQ(ierr);
      ierr = AssemblyResidualSmallDeformation( &x, &b);CHKERRQ(ierr);
      ierr = MicroSetBoundaryOnResidual(dir, &b); CHKERRQ(ierr);
      if( flag_print | (1<<PRINT_PETSC) ){
	ierr = PetscViewerASCIIOpen(MICRO_COMM,"b.dat",&viewer); CHKERRQ(ierr);
	ierr = VecView(b,viewer); CHKERRQ(ierr);
      }
      ierr = VecNorm(b,NORM_2,&norm);CHKERRQ(ierr);
      ierr = PetscPrintf(MICRO_COMM,"|b| = %e\n",norm);CHKERRQ(ierr);
      ierr = VecScale(b,-1.0); CHKERRQ(ierr);
      ierr = PetscLogEventEnd(EVENT_ASSEMBLY_RES,0,0,0,0);CHKERRQ(ierr);
      if( !(norm > NormTol) )break;
      /*
	 Assemblying Jacobian
       */
      ierr = PetscLogEventBegin(EVENT_ASSEMBLY_JAC,0,0,0,0);CHKERRQ(ierr);
      ierr = PetscPrintf(MICRO_COMM, "Assembling Jacobian\n");
      ierr = AssemblyJacobianSmallDeformation(&A);
      ierr = MicroSetBoundaryOnJacobian(dir, &A); CHKERRQ(ierr);
      if( flag_print == (1<<PRINT_PETSC) ){
	ierr = PetscViewerASCIIOpen(MICRO_COMM,"A.dat",&viewer); CHKERRQ(ierr);
	ierr = MatView(A,viewer); CHKERRQ(ierr);
      }
      ierr = PetscLogEventEnd(EVENT_ASSEMBLY_JAC,0,0,0,0);CHKERRQ(ierr);
      /*
	 Solving Problem
       */
      ierr = PetscLogEventBegin(EVENT_SOLVE_SYSTEM,0,0,0,0);CHKERRQ(ierr);
      ierr = PetscPrintf(MICRO_COMM, "Solving Linear System ");
      ierr = KSPSolve(ksp,b,dx);CHKERRQ(ierr);
      ierr = KSPGetIterationNumber(ksp,&kspits);CHKERRQ(ierr);
      ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
      ierr = KSPGetResidualNorm(ksp,&kspnorm);CHKERRQ(ierr);
      ierr = VecAXPY( x, 1.0, dx); CHKERRQ(ierr);
      if( flag_print == (1<<PRINT_PETSC) ){
	ierr = PetscViewerASCIIOpen(MICRO_COMM,"dx.dat",&viewer); CHKERRQ(ierr);
	ierr = VecView(dx,viewer); CHKERRQ(ierr);
	ierr = PetscViewerASCIIOpen(MICRO_COMM,"x.dat",&viewer); CHKERRQ(ierr);
	ierr = VecView(x,viewer); CHKERRQ(ierr);
      }
      ierr = PetscPrintf(MICRO_COMM,"Iterations %D Norm %e reason %d\n",kspits, kspnorm, reason);CHKERRQ(ierr);
      ierr = PetscLogEventBegin(EVENT_SOLVE_SYSTEM,0,0,0,0);CHKERRQ(ierr);

      nr_its ++;
    }

    if(flag_print | (1<<PRINT_VTKDISP)){ 
      strain = malloc(nelm*6*sizeof(double));
      stress = malloc(nelm*6*sizeof(double));
      ierr = SpuCalcStressOnElement(&x, strain, stress);
      sprintf(vtkfile_n,"%s_displ_%d.vtk",myname,rank_mic);
      ierr = SpuVTKPlot_Displ_Strain_Stress(MICRO_COMM, vtkfile_n, &x, strain, stress);
      free(stress); free(strain);
    }

  }

  if(!flag_coupling){
    PetscPrintf(MICRO_COMM,
	"--------------------------------------------------\n"
	"  MICRO: FINISH COMPLETE\n"
	"--------------------------------------------------\n");
  }

  ierr = PetscFinalize();
  ierr = MPI_Finalize();

  return 0;
}
