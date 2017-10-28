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
"-mat_fiber_t0  [E,nu] : material \"FIBER\" type_0 given by command line\n"
"-mat_matrix_t0 [E,nu] : material \"FIBER\" type_0 given by command line\n"
"-homo_taylor_s     : c =  vi ci + vm cm            (serial)\n"
"-homo_taylor_p     : c = (vi ci^-1 + vm cm^-1)^-1  (parallel)\n"
"-homo_us           : homogenization using uniform strains approach\n"
"-fiber_cilin <r,dx,dy,dz>\n"
"-fiber_nx <nx>\n"
"-fiber_ny <ny>\n"
"-struct_mesh [<nx,ny>] if dim = 2\n"
"-print_petsc\n"
"-print_vtk\n"
"-print_part\n"
"-print_vtu\n"
"-print_all\n";

#include "micro.h"

int main(int argc, char **argv)
{

  int        i, ierr, ierr_1=0;
  int        nval;
  char       vtkfile_n[NBUF];
  PetscBool  set;

  myname            = strdup("micro");
  flag_linear_micro = 0;
  first_time_homo   = 1;
  energy_interp     = NULL;
  flag_struct_mesh  = false;

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
  if(set == PETSC_TRUE){
    flag_coupling = PETSC_TRUE;
    macmic.type = COUP_1;
  }

  /* Stablish a new local communicator */
  color = MICRO;
  ierr = macmic_coloring(WORLD_COMM, &color, &macmic, &MICRO_COMM); /* color can change */
  if(ierr){
    ierr_1 = 1;
    PetscPrintf(MICRO_COMM,"micro: problem during coloring\n");
    goto end_mic_0;
  }

  ierr = MPI_Comm_size(MICRO_COMM, &nproc_mic);
  ierr = MPI_Comm_rank(MICRO_COMM, &rank_mic);
  
end_mic_0:
  ierr = PetscFinalize();CHKERRQ(ierr);
  if(ierr_1) goto end_mic_2;

  /*
     Set PETSc communicator to MICRO_COMM
     and start again
   */
  PETSC_COMM_WORLD = MICRO_COMM;
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);


  /* 
     Dimension, mesh and input Options 

     The problem can be execute using the options
     -struct_mesh <nx,ny> (if dim=2)
     -struct_size <lx,ly> (if dim=2)

     or given unstructured meshes as gmsh format
   
   */
  PetscOptionsGetInt(NULL, NULL, "-dim", &dim, &set);
  if(set == PETSC_FALSE){
    PetscPrintf(MPI_COMM_SELF,"dimension (-dim <dim>) not given\n");
    ierr_1 = 1;
    goto end_mic_1;
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
      goto end_mic_1;
  }

  {
    /* struct mesh */

    int    nval_expect;
    int    struct_mesh_n[3];
    double struct_mesh_l[3];

    if( dim == 2 ) nval_expect = nval = 2;
    if( dim == 3 ) nval_expect = nval = 3;
    PetscOptionsGetIntArray(NULL, NULL, "-struct_mesh",struct_mesh_n,&nval,&set);
    flag_struct_mesh  = true;
    if( set == PETSC_TRUE )
    {
      if( nval != nval_expect ){
	PetscPrintf(MPI_COMM_SELF,"-struct_mesh should include %d double arguments\n", nval_expect);
	ierr_1 = 1;
	goto end_mic_0;
      }
      nx = struct_mesh_n[0];
      ny = struct_mesh_n[1];
      if( dim == 3 ) nz = struct_mesh_n[2];
    }

    PetscOptionsGetRealArray(NULL, NULL, "-struct_size",struct_mesh_l,&nval,&set);
    if( set == PETSC_TRUE )
    {
      if( nval != nval_expect ){
	PetscPrintf(MPI_COMM_SELF,"-struct_size should include %d double arguments\n", nval_expect);
	ierr_1 = 1;
	goto end_mic_0;
      }
      nx = struct_mesh_n[0];
      ny = struct_mesh_n[1];
      if( dim == 3 ) nz = struct_mesh_n[0];
    }
    else if( flag_struct_mesh == true ){
	PetscPrintf(MPI_COMM_SELF,"-struct_size should be include because -struct_mesh was specified\n");
	ierr_1 = 1;
	goto end_mic_0;
    }
    lx = struct_mesh_l[0];
    ly = struct_mesh_l[1];
    if( dim == 3 ) lz = struct_mesh_l[2];
  }
  if( flag_struct_mesh == false){

    /* mesh from gmsh */

    mesh_f = FORMAT_NULL;
    PetscOptionsHasName(NULL,NULL,"-mesh_gmsh",&set);CHKERRQ(ierr);
    if(set == PETSC_TRUE) mesh_f = FORMAT_GMSH;
    PetscOptionsHasName(NULL,NULL,"-mesh_alya",&set);CHKERRQ(ierr);
    if(set == PETSC_TRUE) mesh_f = FORMAT_ALYA;
    if(mesh_f == FORMAT_NULL)SETERRQ(MICRO_COMM,1,"mesh format not given on command line.");
    PetscOptionsGetString(NULL, NULL, "-mesh",mesh_n,128,&set);
    if(set == PETSC_FALSE) SETERRQ(MICRO_COMM,1,"mesh file not given on command line.");
    PetscOptionsGetString(NULL, NULL, "-input", input_n,128,&set);
    if(set == PETSC_FALSE) SETERRQ(MICRO_COMM,1,"input file not given.");
  }


  /* Mesh partition algorithms */
  partition_algorithm = PARMETIS_MESHKWAY;
  PetscOptionsHasName(NULL,NULL,"-part_meshkway",&set);
  if(set==PETSC_TRUE) partition_algorithm = PARMETIS_MESHKWAY;
  PetscOptionsHasName(NULL,NULL,"-part_geom",&set);
  if(set==PETSC_TRUE) partition_algorithm = PARMETIS_GEOM;

  /* Homogenization Options */
  homo_type=0;
  PetscOptionsHasName(NULL,NULL,"-homo_taylor_p",&set);
  if(set==PETSC_TRUE) homo_type = TAYLOR_P;
  PetscOptionsHasName(NULL,NULL,"-homo_taylor_s",&set);
  if(set==PETSC_TRUE) homo_type = TAYLOR_S;
  PetscOptionsHasName(NULL,NULL,"-homo_us",&set);
  if(set==PETSC_TRUE) homo_type = UNIF_STRAINS;
  if(homo_type==0){
    PetscPrintf(MPI_COMM_SELF,"no homogenization option specified\n");
    ierr_1 = 1;
    goto end_mic_0;
  }

  /* Solver Options */
  PetscOptionsGetInt(NULL, NULL, "-nr_max_its", &nr_max_its, &set);
  if(set==PETSC_FALSE) nr_max_its=5;
  PetscOptionsGetReal(NULL, NULL, "-nr_norm_tol", &nr_norm_tol, &set);
  if(set==PETSC_FALSE) nr_norm_tol=1.0e-7;

  /* 
     Geometry specifications or material distribution
     inside the rve
   */

  /* Fiber in the middle */
  flag_fiber_cilin = 0;
  nval = 4;
  PetscOptionsGetRealArray(NULL, NULL, "-fiber_cilin", fiber_cilin_vals, &nval,&set);
  if(set==PETSC_TRUE) {
    flag_fiber_cilin=1;
    fiber_cilin_r=fiber_cilin_vals[0];
    if(nval==1){
      for(i=0;i<dim;i++){
	fiber_cilin_center_devi[i]=0.0;
      }
    }
    else if(nval==0){
      PetscPrintf(MPI_COMM_SELF,"-fiber_cilin specified with no argument\n");
      ierr_1 = 1;
      goto end_mic_0;
    }
    else{
      for(i=0;i<nval;i++){
	fiber_cilin_center_devi[i]=fiber_cilin_vals[1+i];
      }
    }
  }
  PetscOptionsGetInt(NULL, NULL, "-fiber_nx", &nx_fibers, &set);
  if(set==PETSC_FALSE) nx_fibers = 1;
  PetscOptionsGetInt(NULL, NULL, "-fiber_ny", &ny_fibers, &set);
  if(set==PETSC_FALSE) ny_fibers = 1;

  {
    /* 
       Materials by command line 

       We expect a fiber and a matrix material only
     */
    material_t mat;
    double E, v;
    list_init(&material_list,sizeof(material_t),NULL);

    nval = 3;
    PetscOptionsGetRealArray(NULL, NULL, "-mat_fiber_t0",mat_fiber_t0,&nval,&set);
    if( set == PETSC_TRUE )
    {
      if( nval != 3 ){
	PetscPrintf(MPI_COMM_SELF,"-mat_fiber_t0 should include 3 double arguments\n");
	ierr_1 = 1;
	goto end_mic_0;
      }

      mat.type_id = TYPE_0;
      mat.name = strdup("FIBER");
      mat.type = malloc(sizeof(type_0));
      ((type_0*)mat.type)->rho         = mat_fiber_t0[0];
      E = ((type_0*)mat.type)->young   = mat_fiber_t0[1];
      v = ((type_0*)mat.type)->poisson = mat_fiber_t0[2];
      ((type_0*)mat.type)->lambda      = (E*v)/((1+v)*(1-2*v));
      ((type_0*)mat.type)->mu          = E/(2*(1+v));

      list_insertlast( &material_list , &mat );
    }

    PetscOptionsGetRealArray(NULL,NULL,"-mat_matrix_t0",mat_matrix_t0,&nval,&set);
    if( set == PETSC_TRUE )
    {
      if( nval != 3 ){
	PetscPrintf(MPI_COMM_SELF,"-mat_matrix_t0 should include 3 double arguments\n");
	ierr_1 = 1;
	goto end_mic_0;
      }

      mat.type_id = TYPE_0;
      mat.name = strdup("MATRIX");
      mat.type = malloc(sizeof(type_0));
      ((type_0*)mat.type)->rho         = mat_matrix_t0[0];
      E = ((type_0*)mat.type)->young   = mat_matrix_t0[1];
      v = ((type_0*)mat.type)->poisson = mat_matrix_t0[2];
      ((type_0*)mat.type)->lambda      = (E*v)/((1+v)*(1-2*v));
      ((type_0*)mat.type)->mu          = E/(2*(1+v));

      list_insertlast( &material_list , &mat );
    }
  }

  /* Printing Options */
  flag_print = 0;
  PetscOptionsHasName(NULL,NULL,"-print_petsc",&set);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_PETSC);
  PetscOptionsHasName(NULL,NULL,"-print_vtk",&set);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_VTK);
  PetscOptionsHasName(NULL,NULL,"-print_part",&set);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_VTKPART);
  PetscOptionsHasName(NULL,NULL,"-print_vtu",&set);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_VTU);
  PetscOptionsHasName(NULL,NULL,"-print_all",&set);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_ALL);
  //**************************************************

  if(!flag_coupling){
    PetscPrintf(MICRO_COMM,
	"--------------------------------------------------\n"
	"  MICRO: STANDALONE \n"
	"--------------------------------------------------\n");
  }

  file_out = NULL;
  if(rank_mic==0) file_out = fopen("micro_structures.dat","w");

#if defined(PETSC_USE_LOG)
  PetscLogEventRegister("ASSEMBLY_JAC",PETSC_VIEWER_CLASSID,&EVENT_ASSEMBLY_JAC);
  PetscLogEventRegister("ASSEMBLY_RES",PETSC_VIEWER_CLASSID,&EVENT_ASSEMBLY_RES);
  PetscLogEventRegister("SOLVE_SYSTEM",PETSC_VIEWER_CLASSID,&EVENT_SOLVE_SYSTEM);
  PetscLogEventRegister("ASSEMBLY_JAC",PETSC_VIEWER_CLASSID,&EVENT_INIT);
#endif

  PetscLogEventBegin(EVENT_INIT,0,0,0,0);
  if( flag_struct_mesh == false)
  {
    /* read mesh */    
    if(!flag_coupling)
      PetscPrintf(MICRO_COMM,"Reading mesh elements\n");
    ierr = read_mesh_elmv(MICRO_COMM, myname, mesh_n, mesh_f);
    if(ierr){
      goto end_mic_1;
    }

    /* Partition the mesh */
    if(!flag_coupling)
      PetscPrintf(MICRO_COMM,"Partitioning and distributing mesh\n");
    ierr = part_mesh_PARMETIS(&MICRO_COMM, time_fl, myname, NULL);

    /* Calculate <*ghosts> and <nghosts> */
    if(!flag_coupling)
      PetscPrintf(MICRO_COMM,"Calculating Ghost Nodes\n");
    ierr = calculate_ghosts(&MICRO_COMM, myname);

    /* Reenumerate Nodes */
    if(!flag_coupling)
      PetscPrintf(MICRO_COMM,"Reenumering nodes\n");
    ierr = reenumerate_PETSc(MICRO_COMM);

    /* Coordinate Reading */
    if(!flag_coupling)
      PetscPrintf(MICRO_COMM,"Reading Coordinates\n");
    ierr = read_mesh_coord(MICRO_COMM, mesh_n, mesh_f);

    if(flag_print & (1<<PRINT_VTKPART)){
      sprintf(vtkfile_n,"%s_part_%d.vtk",myname,rank_mic);
      ierr = spu_vtk_partition( vtkfile_n, &MICRO_COMM );
    }

    list_init(&physical_list, sizeof(physical_t), NULL);
    list_init(&boundary_list, sizeof(boundary_t), NULL);

    /* Read Physical entities */
    ierr = read_physical_entities(MICRO_COMM, mesh_n, mesh_f);
    if(ierr){
      PetscPrintf(MICRO_COMM,"Problem parsing physical entities from mesh file\n");
      goto end_mic_1;
    }
    ierr = mic_parse_boundary(MICRO_COMM, input_n);
    if(ierr){
      PetscPrintf(MICRO_COMM,"Problem parsing physical entities from mesh file\n");
      goto end_mic_1;
    }
    ierr = set_id_on_material_and_boundary(MICRO_COMM);
    ierr = check_elm_id();
    ierr = mic_check_linear_material();

    /* Read boundary */
    ierr = read_boundary(MICRO_COMM, mesh_n, mesh_f);

    /* Allocate matrices & vectors */ 
    ierr = mic_alloc(MICRO_COMM);

    /* Setting solver options */
    ierr = KSPCreate(MICRO_COMM,&ksp); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp,A,A); CHKERRQ(ierr);

    /* Init Gauss point shapes functions and derivatives */
    ierr = fem_inigau();

    /* micro main coupling loop */

    ierr = get_bbox_limit_lengths(MICRO_COMM,coord,nmynods,&LX,&LY,&LZ);
    PetscPrintf(MICRO_COMM,"LX=%e LY=%e LZ=%e\n",LX,LY,LZ);
    ierr = get_domain_center(MICRO_COMM, coord, nmynods, center_domain);
    PetscPrintf(MICRO_COMM,"center = %e %e %e\n",center_domain[0],center_domain[1],center_domain[2]);
    ierr = calc_rho(MICRO_COMM, &rho);
    PetscPrintf(MICRO_COMM,"density = %e\n", rho);
  }

  PetscLogEventEnd(EVENT_INIT,0,0,0,0);

  double strain_mac[6], strain_ave[6], stress_ave[6], c_homo[36];

  if(flag_coupling){

    /*
       COUPLING EXECUTION

       1) waits instruction 
       2) execute instruction
       3) finish if instruction = SIGNAL_MICRO_END  
     */
    int signal=-1;

    while(signal!=MIC_END){

      /* Receive instruction */
      ierr = mic_recv_signal(WORLD_COMM, &signal);

      switch(signal){

	case MAC2MIC_STRAIN:
	  /* Wait for strain */
	  ierr = mic_recv_strain(WORLD_COMM, strain_mac);
	  /* Performs the micro localization + homogenization */
	  ierr = mic_calc_stress_ave(MICRO_COMM, strain_mac, strain_ave, stress_ave);
	  /* Send Stress */
	  ierr = mic_send_stress(WORLD_COMM, stress_ave);
	  break;

	case C_HOMO:
	  /* Wait for strain */
	  ierr = mic_recv_strain(WORLD_COMM, strain_mac);
	  /* Wait for macro_gp number */
	  ierr = mic_recv_macro_gp(WORLD_COMM, &macro_gp);
	  /* Calculates C homogenized */
	  ierr = mic_calc_c_homo(MICRO_COMM, strain_mac, c_homo);
	  /* Send C homogenized */
	  ierr = mic_send_c_homo(WORLD_COMM, c_homo);
	  break;

	case RHO:
	  /* Send rho homogenized */
	  ierr = mic_send_rho(WORLD_COMM, &rho);
	  break;

	case MIC_END:
	  break;

	default:
	  PetscPrintf(MICRO_COMM,"MICRO:signal %d not identified\n",signal);
	  goto end_mic_1;

      }
    }
  }
  else{

    /*
       STANDALONE EXECUTION

       We perform 6 homogenization using :

       dim = 2

       e0 = (  0.005 0 0  )
       e1 = (  0 0.005 0  )
       e2 = (  0 0 0.005  )

       dim = 3
       
       e0 = (  0.005 0 0 0 0 0     )
       e1 = (  0 0.005 0 0 0 0     )
       e2 = (  0 0 0.005 0 0 0     )
       e3 = (  0 0 0 0.005 0 0     )
       e4 = (  0 0 0 0 0 0 0.005   )
       e5 = (  0 0 0 0 0 0 0 0.005 )
     */
    int j;
    double strain_mac[6];

    memset(c_homo,0.0,36*sizeof(double));
    for(i=0;i<nvoi;i++){

      memset(strain_mac,0.0,nvoi*sizeof(double)); strain_mac[i]=0.005;
      ierr = mic_homogenize(MICRO_COMM, strain_mac, strain_ave, stress_ave);
      if(ierr){
	goto end_mic_1;
      }

      ierr = PetscPrintf(MICRO_COMM,"\nstrain_ave = ");
      for(j=0;j<nvoi;j++){
	ierr = PetscPrintf(MICRO_COMM,"%e ",strain_ave[j]);
      }
      ierr = PetscPrintf(MICRO_COMM,"\nstress_ave = ");
      for(j=0;j<nvoi;j++){
	ierr = PetscPrintf(MICRO_COMM,"%e ",stress_ave[j]);
      }
      ierr = PetscPrintf(MICRO_COMM,"\n");
      for(j=0;j<nvoi;j++){
	c_homo[j*nvoi+i] = stress_ave[j] / strain_ave[i];
      }

      if(flag_print & (1<<PRINT_VTK | 1<<PRINT_VTU) && homo_type==UNIF_STRAINS){
	strain = malloc(nelm*nvoi*sizeof(double));
	stress = malloc(nelm*nvoi*sizeof(double));
	energy = malloc(nelm*sizeof(double));
	ierr = assembly_residual_sd(&x, &b);CHKERRQ(ierr);
	ierr = calc_strain_stress_energy(&x, strain, stress, energy);
	if(flag_print & (1<<PRINT_VTK)){ 
	  sprintf(vtkfile_n,"%s_exp%d_%d.vtk",myname,i,rank_mic);
	  ierr = write_vtk(MICRO_COMM, vtkfile_n, &x, strain, stress);
	}
	if(flag_print & (1<<PRINT_VTU)){ 
	  sprintf(vtkfile_n,"%s_exp%d",myname,i);
	  ierr = write_vtu(MICRO_COMM, vtkfile_n, &x, &b, strain, stress, energy);
	  if(ierr){
	    PetscPrintf(MICRO_COMM,"Problem writing vtu file\n");
	    goto end_mic_1;
	  }
	}
	free(stress); free(strain); free(energy);
      }

    }
    PetscPrintf(MICRO_COMM,"\nConstitutive Average Tensor\n");
    for(i=0;i<nvoi;i++){
      for(j=0;j<nvoi;j++){
	PetscPrintf(MICRO_COMM,"%e ",(fabs(c_homo[i*nvoi+j])>1.0)?c_homo[i*nvoi+j]:0.0);
      }
      PetscPrintf(MICRO_COMM,"\n");
    }
    PetscPrintf(MICRO_COMM,"\n");

    /*
       Experiment to test if the homogenization with <strain_mac>
       gives the same <stress_ave> than doing <c_homo>*<strain_mac>
     */

    strain_mac[0] = 0.01; strain_mac[1] = -0.02; strain_mac[2] = +0.03;
    strain_mac[3] = 0.01; strain_mac[4] = -0.02; strain_mac[5] = +0.03;
    ierr = mic_homogenize(MICRO_COMM, strain_mac, strain_ave, stress_ave);
    ierr = PetscPrintf(MICRO_COMM,"\nstrain_ave = ");
    for(j=0;j<nvoi;j++){
      PetscPrintf(MICRO_COMM,"%e ",strain_ave[j]);
    }
    PetscPrintf(MICRO_COMM,"\nstress_ave = ");
    for(j=0;j<nvoi;j++){
      PetscPrintf(MICRO_COMM,"%e ",stress_ave[j]);
    }
    for(i=0;i<nvoi;i++){
      stress_ave[i] = 0.0;
      for(j=0;j<nvoi;j++){
	stress_ave[i] +=  c_homo[i*nvoi+j] * strain_mac[j];
      }
    }
    PetscPrintf(MICRO_COMM,"\nstress_ave = ");
    for(j=0;j<nvoi;j++){
      ierr = PetscPrintf(MICRO_COMM,"%e ",stress_ave[j]);
    }
    PetscPrintf(MICRO_COMM," (c_homo*strain_mac)\n");

  }

end_mic_1:

  if(!flag_coupling){
    PetscPrintf(MICRO_COMM,
	"--------------------------------------------------\n"
	"  MICRO: FINISH COMPLETE\n"
	"--------------------------------------------------\n");
  }

  if(rank_mic==0) fclose(file_out); 

  ierr = PetscFinalize();

end_mic_2:
  ierr = MPI_Finalize();


  return 0;
}
