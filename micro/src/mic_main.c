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
"-print_petsc\n"
"-print_vtk\n"
"-print_part\n"
"-print_vtu\n"
"-print_all\n"
"[-homo_taylor -homo_linear]\n";

#include "micro.h"

int main(int argc, char **argv)
{

  int ierr;
  char *myname = strdup("micro");
  char vtkfile_n[NBUF];
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
  /*
     Mesh and Input Options
  */
  mesh_f = FORMAT_NULL;
  ierr = PetscOptionsHasName(NULL,NULL,"-mesh_gmsh",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) mesh_f = FORMAT_GMSH;
  ierr = PetscOptionsHasName(NULL,NULL,"-mesh_alya",&set);CHKERRQ(ierr);
  if(set == PETSC_TRUE) mesh_f = FORMAT_ALYA;
  if(mesh_f == FORMAT_NULL)SETERRQ(MICRO_COMM,1,"mesh format not given on command line.");
  ierr = PetscOptionsGetString(NULL, NULL, "-mesh", mesh_n, 128, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE) SETERRQ(MICRO_COMM,1,"mesh file not given on command line.");
  ierr = PetscOptionsGetString(NULL, NULL, "-input", input_n, 128, &set); CHKERRQ(ierr); 
  if(set == PETSC_FALSE) SETERRQ(MICRO_COMM,1,"input file not given.");
  /*
     Homogenization Options
  */
  homo_type=0;
  ierr = PetscOptionsHasName(NULL,NULL,"-homo_taylor",&set);CHKERRQ(ierr);
  if(set==PETSC_TRUE) homo_type=HOMO_TAYLOR;
  ierr = PetscOptionsHasName(NULL,NULL,"-homo_linear",&set);CHKERRQ(ierr);
  if(set==PETSC_TRUE) homo_type=HOMO_LINEAR;
  ierr = PetscOptionsHasName(NULL,NULL,"-homo_linear_hexa",&set);CHKERRQ(ierr);
  if(set==PETSC_TRUE) homo_type=HOMO_LINEAR_HEXA;
  if(homo_type==0)SETERRQ(MICRO_COMM,1,"no homogenization option specified");

  /* 
     Stablish a new local communicator
   */
  color = MICRO;
  ierr = macmic_coloring(WORLD_COMM, &color, &macmic, &MICRO_COMM); /* color can change */

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
    PetscPrintf(MICRO_COMM,"Reading mesh elements\n");
  ierr = read_mesh_elmv(MICRO_COMM, myname, mesh_n, mesh_f);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_READ_MESH_ELEM,0,0,0,0);CHKERRQ(ierr);

  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /*
     partition the mesh
  */
  if(!flag_coupling)
    PetscPrintf(MICRO_COMM,"Partitioning and distributing mesh\n");
  ierr = PetscLogEventBegin(EVENT_PART_MESH,0,0,0,0);CHKERRQ(ierr);
  ierr = part_mesh_PARMETIS(&MICRO_COMM, time_fl, myname, NULL, PARMETIS_MESHKWAY );CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_PART_MESH,0,0,0,0);CHKERRQ(ierr);

  /*
     Calculate <*ghosts> and <nghosts> 
  */
  if(!flag_coupling)
    PetscPrintf(MICRO_COMM,"Calculating Ghost Nodes\n");
  ierr = PetscLogEventBegin(EVENT_CALC_GHOSTS,0,0,0,0);CHKERRQ(ierr);
  ierr = calculate_ghosts(&MICRO_COMM, myname);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_CALC_GHOSTS,0,0,0,0);CHKERRQ(ierr);

  /*
     Reenumerate Nodes
  */
  if(!flag_coupling)
    PetscPrintf(MICRO_COMM,"Reenumering nodes\n");
  ierr = PetscLogEventBegin(EVENT_REENUMERATE,0,0,0,0);CHKERRQ(ierr);
  ierr = reenumerate_PETSc(MICRO_COMM);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_REENUMERATE,0,0,0,0);CHKERRQ(ierr);

  /*
     Coordinate Reading
  */
  if(!flag_coupling)
    PetscPrintf(MICRO_COMM,"Reading Coordinates\n");
  ierr = PetscLogEventBegin(EVENT_READ_COORD,0,0,0,0);CHKERRQ(ierr);
  ierr = read_mesh_coord(MICRO_COMM, mesh_n, mesh_f);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(EVENT_READ_COORD,0,0,0,0);CHKERRQ(ierr);

  if(flag_print & (1<<PRINT_VTKPART)){
    sprintf(vtkfile_n,"%s_part_%d.vtk",myname,rank_mic);
    ierr = spu_vtk_partition( vtkfile_n, &MICRO_COMM );
  }

  /*
     Read materials, physical entities, boundaries from input and mesh file
  */
  ierr = list_init(&physical_list, sizeof(physical_t), NULL);CHKERRQ(ierr);
  ierr = list_init(&function_list, sizeof(physical_t), NULL);CHKERRQ(ierr);
  ierr = parse_material(MICRO_COMM, input_n);CHKERRQ(ierr);
  ierr = read_physical_entities(MICRO_COMM, mesh_n, mesh_f);CHKERRQ(ierr);
  ierr = parse_function(MICRO_COMM, input_n);CHKERRQ(ierr); 

  /*
     Check if Physical Entities <P000 P100 P010 X0 X1 Y0 Y1 Z0 Z1> exist in mesh
     Creates the boundary_list with <P000 P100 P010 X0 X1 Y0 Y1 Z0 Z1>
  */
  ierr = micro_check_physical_entities(&physical_list);CHKERRQ(ierr);
  ierr = micro_init_boundary_list(&boundary_list);CHKERRQ(ierr);
  ierr = set_id_on_material_and_boundary(MICRO_COMM); CHKERRQ(ierr); 
  ierr = CheckPhysicalID(); CHKERRQ(ierr);
  ierr = read_boundary(MICRO_COMM, mesh_n, mesh_f);CHKERRQ(ierr);
  ierr = micro_init_boundary(MICRO_COMM, &boundary_list );CHKERRQ(ierr);

  /*
     Allocate matrices & vectors
  */ 
  ierr = PetscLogEventBegin(EVENT_ALLOC_MATVEC,0,0,0,0);CHKERRQ(ierr);
  ierr = MicroAllocMatrixVector( MICRO_COMM, nmynods*3, NTotalNod*3);CHKERRQ(ierr);
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
  double strain_mac[6], strain_ave[6], stress_ave[6], ttensor[36];

  ierr = get_bbox_limit_lengths(MICRO_COMM,coord,nmynods,&LX,&LY,&LZ);CHKERRQ(ierr);
  ierr = PetscPrintf(MICRO_COMM,"LX=%e LY=%e LZ=%e\n",LX,LY,LZ);CHKERRQ(ierr);

  if(flag_coupling){

    /*
       COUPLING EXECUTION

       In this mode <micro> is used to homogenize RVE properties 
       and send them to <macro> program

       1) waits instruction 
       2) execute instruction
       3) finish if instruction = SIGNAL_MICRO_END  

     */
    int signal=-1;

    while(signal!=MIC_END)
    {
      ierr = mic_recv_signal(WORLD_COMM, &signal);CHKERRQ(ierr);

      switch(signal)
      {
	case MAC2MIC_STRAIN:
	  /*
	     Wait for strain
	  */
	  ierr = mic_recv_strain(WORLD_COMM, strain_mac);CHKERRQ(ierr);
	  /*
	     Performs the micro homogenization
	  */
	  ierr = micro_homogenize(MICRO_COMM, strain_mac, strain_ave, stress_ave);CHKERRQ(ierr);
	  /*
	     Send Stress
	  */
	  ierr = mic_send_stress(WORLD_COMM, stress_ave);CHKERRQ(ierr);
	  break;
	case MIC_END:
	  break;
	default:
	  SETERRQ(MICRO_COMM,1,"can no identify recv signal.");
      }
    }
  }
  else{

    /*
       STANDALONE EXECUTION

       We perform 6 homogenization using :
       
       e0 = (    0.005 0 0 0 0 0     )
       e1 = (    0 0.005 0 0 0 0     )
       e2 = (    0 0 0.005 0 0 0     )
       e3 = (    0 0 0 0.005 0 0     )
       e4 = (    0 0 0 0 0 0 0.005   )
       e5 = (    0 0 0 0 0 0 0 0.005 )
     */
    int i, j;
    double strain_mac[6];

    memset(ttensor,0.0,36*sizeof(double));
    for(i=0;i<6;i++){

      memset(strain_mac,0.0,6*sizeof(double));strain_mac[i]=0.005;
      ierr = micro_homogenize(MICRO_COMM, strain_mac, strain_ave, stress_ave);CHKERRQ(ierr);

      ierr = PetscPrintf(MICRO_COMM,"\nstrain_ave = ");CHKERRQ(ierr);
      for(j=0;j<6;j++){
	ierr = PetscPrintf(MICRO_COMM,"%e ",strain_ave[j]);CHKERRQ(ierr);
      }
      ierr = PetscPrintf(MICRO_COMM,"\nstress_ave = ");CHKERRQ(ierr);
      for(j=0;j<6;j++){
	ierr = PetscPrintf(MICRO_COMM,"%e ",stress_ave[j]);CHKERRQ(ierr);
      }
      ierr = PetscPrintf(MICRO_COMM,"\n");CHKERRQ(ierr);
      for(j=0;j<6;j++){
	ttensor[j*6+i] = stress_ave[j] / strain_ave[i];
      }

      if(flag_print & (1<<PRINT_VTK | 1<<PRINT_VTU)){ 
	strain = malloc(nelm*6*sizeof(double));
	stress = malloc(nelm*6*sizeof(double));
	ierr = SpuCalcStressOnElement(&x, strain, stress);
	if(flag_print & (1<<PRINT_VTK)){ 
	  sprintf(vtkfile_n,"%s_displ_exp%d_%d.vtk",myname,i,rank_mic);
	  ierr = SpuVTKPlot_Displ_Strain_Stress(MICRO_COMM, vtkfile_n, &x, strain, stress);
	}
	if(flag_print & (1<<PRINT_VTU)){ 
	  sprintf(vtkfile_n,"%s_displ_exp%d",myname,i);
	  ierr = write_vtu(MICRO_COMM, vtkfile_n, &x, strain, stress);CHKERRQ(ierr);
	}
	free(stress); free(strain);
      }

    }
    ierr = PetscPrintf(MICRO_COMM,"\nConstitutive Average Tensor\n");CHKERRQ(ierr);
    for(i=0;i<6;i++){
      for(j=0;j<6;j++){
	ierr = PetscPrintf(MICRO_COMM,"%e ",(fabs(ttensor[i*6+j])>1.0)?ttensor[i*6+j]:0.0);CHKERRQ(ierr);
      }ierr = PetscPrintf(MICRO_COMM,"\n");CHKERRQ(ierr);
    }ierr = PetscPrintf(MICRO_COMM,"\n");CHKERRQ(ierr);

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
