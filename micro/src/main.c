/*

   MICRO main function

   Program for solving the microscopic problem 
   for multi-scale approach.

   Author> Guido Giuntoli
   Date> 28-07-2017

 */

#include "micro.h"

static char help[] = 
"MICRO MULTISCALE CODE\n"
"Solves the RVE problem inside a solid structure.                                             \n"
"It has the capability of being couple with MACRO.                                            \n"
"-coupl    [0 (no coupling ) | 1 (coupling with micro)]                                       \n"
"-testcomm [0 (no test) | 1 (sends a strain value and receive a stress calculated from micro)]\n"
"-mat_fiber_t0  [E,nu] : material \"FIBER\" type_0 given by command line                      \n"
"-mat_matrix_t0 [E,nu] : material \"FIBER\" type_0 given by command line                      \n"
"-homo_ts     : c =  vi ci + vm cm            (serial)                                        \n"
"-homo_tp     : c = (vi ci^-1 + vm cm^-1)^-1  (parallel)                                      \n"
"-homo_us     : homogenization using uniform strains approach                                 \n"
"-fiber_cilin <r,dx,dy,dz>                                                                    \n"
"-fiber_nx <nx>                                                                               \n"
"-fiber_ny <ny>                                                                               \n"
"-struct_n [<nx,ny>] if dim = 2                                                               \n"
"-struct_l [<lx,ly>] if dim = 2                                                               \n"
"-print_petsc                                                                                 \n"
"-print_vtk                                                                                   \n"
"-print_vtu                                                                                   \n"
"-print_all                                                                                   \n";

int main(int argc, char **argv)
{

  int        i, j, ierr, ierr_1 = 0;
  int        nval;
  PetscBool  set;

  myname            = strdup("micro");
  flag_linear_micro = 0;
  first_time_homo   = 1;
  flag_struct_mesh  = false;
  flag_first_alloc  = true;

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

  MPI_Comm_size(MICRO_COMM, &nproc_mic);
  MPI_Comm_rank(MICRO_COMM, &rank_mic);
  
end_mic_0:
  ierr = PetscFinalize();
  if(ierr_1) goto end_mic_2;

  /* Set PETSc communicator to MICRO_COMM and start */

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
  nvoi = (dim == 2) ? 3 : 6;

  {
    /* struct mesh */

    int    nval_expect;
    int    struct_mesh_n[3];
    double struct_mesh_l[3];
    nx = ny = nz = 1;
    hx = hy = hz = -1;
    lx = ly = lz = -1;

    if( dim == 2 ) nval_expect = nval = 2;
    if( dim == 3 ) nval_expect = nval = 3;
    PetscOptionsGetIntArray(NULL, NULL, "-struct_n",struct_mesh_n,&nval,&set);
    if( set == PETSC_TRUE ){

      if( nval != nval_expect ){
	PetscPrintf(MPI_COMM_SELF,"-struct_n should include %d double arguments\n", nval_expect);
	ierr_1 = 1;
	goto end_mic_0;
      }
      nx   = struct_mesh_n[0];
      ny   = struct_mesh_n[1];
      if( dim == 3 ) nz = struct_mesh_n[2];
      nn   = nx*ny*nz;
      nex  = (nx-1);
      ney  = (ny-1) / nproc_mic + (((ny-1) % nproc_mic > rank_mic) ? 1:0); // local number of elements in y direction
      nez  = (nz-1);
      nelm = ( dim == 2 ) ? nex*ney : nex*ney*nez;                         // number of local elements
      nyl  = ( rank_mic == 0 ) ? ney+1 : ney;                              // local number of nodes in y direction
      nl   = ( dim == 2 ) ? nyl*nx : nyl*nx*nz;                            // local number of nodes

      int *nyl_arr = malloc(nproc_mic * sizeof(int));
      ierr = MPI_Allgather( &nyl, 1, MPI_INT, nyl_arr, 1, MPI_INT, MICRO_COMM); if(ierr) return 1;
      ny_inf = 0;
      for( i = 0 ; i < rank_mic ; i++ ){
	ny_inf += nyl_arr[i];
      }
      free(nyl_arr);

      npe  = ( dim == 2 ) ? 4 : 8;
      ngp  = ( dim == 2 ) ? 4 : 8;
      if( !( ny > nproc_mic ) ){
	PetscPrintf(MPI_COMM_SELF,"ny %d not large enough to be executed with %d processes\n", ny, nproc_mic);
	ierr_1 = 1;
	goto end_mic_0;
      }

    }
    else{
	PetscPrintf(MPI_COMM_SELF,"-struct_n is request\n");
	ierr_1 = 1;
	goto end_mic_0;
    }

    PetscOptionsGetRealArray(NULL, NULL, "-struct_l",struct_mesh_l,&nval,&set);
    if( set == PETSC_TRUE )
    {
      if( nval != nval_expect ){
	PetscPrintf(MPI_COMM_SELF,"-struct_l should include %d double arguments\n", nval_expect);
	ierr_1 = 1;
	goto end_mic_0;
      }
      lx = struct_mesh_l[0];
      ly = struct_mesh_l[1];
      if( dim == 3 ) lz = struct_mesh_l[2];
    }
    else{
	PetscPrintf(MPI_COMM_SELF,"-struct_l is request\n");
	ierr_1 = 1;
	goto end_mic_0;
    }

    /* set the elements' size */
    hx = lx/nex;
    hy = ly/(ny-1);
    if( dim == 3 ) hz = lz/nez;

    /* set volumes and center*/
    center_coor = malloc ( dim * sizeof(double));
    if ( dim == 2 ){
      center_coor[0] = lx / 2;
      center_coor[1] = ly / 2;
      vol_elem = hx*hy;
      vol_tot  = lx * ly;
      vol_loc  = lx * (hy*ney);
    }
    else{
      center_coor[0] = lx / 2;
      center_coor[1] = ly / 2;
      center_coor[2] = lz / 2;
      vol_elem = hx*hy*hz;
      vol_tot  = lx * ly       * lz;
      vol_loc  = lx * (hy*ney) * lz;
    }

  }

  /* Homogenization Options */
  homo_type=0;
  PetscOptionsHasName(NULL,NULL,"-homo_tp",&set);
  if(set==PETSC_TRUE) homo_type = TAYLOR_P;
  PetscOptionsHasName(NULL,NULL,"-homo_ts",&set);
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
  if(set==PETSC_FALSE) nr_max_its = 5;
  PetscOptionsGetReal(NULL, NULL, "-nr_norm_tol", &nr_norm_tol, &set);
  if(set==PETSC_FALSE) nr_norm_tol = 1.0e-5;

  /* 
     Geometry specifications or material distribution
     inside the rve

     we should assign a value to "micro_type"
   */


  {
    /* 
       CIRCULAR_FIBER 
       -fiber_cilin x,y,z,r
       -fiber_nx <n> (optional)
       -fiber_ny <n> (optional)
     */
    int     nval = 4;
    double  fiber_cilin_vals[4];
    PetscOptionsGetRealArray(NULL, NULL, "-fiber_cilin", fiber_cilin_vals, &nval,&set);

    if( set == PETSC_TRUE ){
      micro_type = CIRCULAR_FIBER;
      cilin_fiber.radius = fiber_cilin_vals[0];
      if(nval==1){
	for( i = 0 ; i < dim ; i++ )
	  cilin_fiber.deviation[i] = 0.0;
      }
      else if(nval==0){
	PetscPrintf(MPI_COMM_SELF,"-fiber_cilin specified with no argument\n");
	ierr_1 = 1;
	goto end_mic_0;
      }
      else{
	for( i = 0 ; i < dim ; i++ )
	  cilin_fiber.deviation[i] = fiber_cilin_vals[1+i];
      }
    }
    PetscOptionsGetInt(NULL, NULL, "-fiber_nx", &cilin_fiber.nx, &set);
    if( set == PETSC_FALSE ) cilin_fiber.nx = 1;
    PetscOptionsGetInt(NULL, NULL, "-fiber_ny", &cilin_fiber.ny, &set);
    if( set == PETSC_FALSE ) cilin_fiber.ny = 1;
  }

  /**************************************************/
  {
    /* Materials by command line 
       (expect a fiber and a matrix material only)
     */
    double E, v;
    material_t mat;
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
  /**************************************************/

  /**************************************************/
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
  /**************************************************/

  if(!flag_coupling){
    PetscPrintf(MICRO_COMM,
	"--------------------------------------------------\n"
	"  MICRO: STANDALONE \n"
	"--------------------------------------------------\n");
  }

  /**************************************************/
  /* initialize global variables*/
  A   = NULL;
  b   = NULL;
  x   = NULL;
  dx  = NULL;

  /* alloc variables*/
  loc_elem_index = malloc( dim*npe * sizeof(int));
  glo_elem_index = malloc( dim*npe * sizeof(int));
  elem_disp      = malloc( dim*npe * sizeof(double));
  stress_gp      = malloc( nvoi    * sizeof(double));
  strain_gp      = malloc( nvoi    * sizeof(double));
  if( flag_print & ( 1 << PRINT_VTU ) ){
    elem_strain = malloc( nelm*nvoi * sizeof(double));
    elem_stress = malloc( nelm*nvoi * sizeof(double));
    elem_energy = malloc( nelm      * sizeof(double));
    elem_type   = malloc( nelm      * sizeof(int));
  }
  /* Initilize shape functions, derivatives, jacobian, b_matrix */

  init_shapes( &struct_sh, &struct_dsh, &struct_wp );

  /* alloc the B matrix */
  struct_bmat = malloc( nvoi * sizeof(double**));
  for( i = 0 ; i < nvoi  ; i++ ){
    struct_bmat[i] = malloc( npe*dim * sizeof(double*));
    for( j = 0 ; j < npe*dim ; j++ )
      struct_bmat[i][j] = malloc( ngp * sizeof(double));
  }

  /* calc B matrix */
  int is, gp;
  for( gp = 0; gp < ngp ; gp++ ){
    for( is = 0; is < npe ; is++ ){
      if( dim == 2 ){
	struct_bmat[0][is*dim + 0][gp] = struct_dsh[is][0][gp];
	struct_bmat[0][is*dim + 1][gp] = 0;
	struct_bmat[1][is*dim + 0][gp] = 0;
	struct_bmat[1][is*dim + 1][gp] = struct_dsh[is][1][gp];
	struct_bmat[2][is*dim + 0][gp] = struct_dsh[is][1][gp];
	struct_bmat[2][is*dim + 1][gp] = struct_dsh[is][0][gp];
      }
    }
  }

  init_trace( MICRO_COMM, "micro_trace.dat" );
  /**************************************************/

  /* Setting solver options */
  ierr = KSPCreate(MICRO_COMM,&ksp); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

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
	  ierr = mic_send_c_homo(WORLD_COMM, nvoi, c_homo);
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

    /* STANDALONE EXECUTION */

    double strain_mac[6];

    memset(c_homo,0.0,36*sizeof(double));
    for( i = 0 ; i < nvoi ; i++ )
    {
      memset(strain_mac,0.0,nvoi*sizeof(double)); strain_mac[i]=0.005;
      ierr = mic_homogenize(MICRO_COMM, strain_mac, strain_ave, stress_ave);
      if(ierr) goto end_mic_1;

      PetscPrintf(MICRO_COMM,"\nstrain_ave = ");
      for( j = 0 ; j < nvoi ; j++ )
	PetscPrintf(MICRO_COMM,"%e ",strain_ave[j]);

      PetscPrintf(MICRO_COMM,"\nstress_ave = ");
      for( j = 0 ; j < nvoi ; j++ )
	PetscPrintf(MICRO_COMM,"%e ",stress_ave[j]);
      PetscPrintf(MICRO_COMM,"\n");

      for( j = 0 ; j < nvoi ; j++ )
	c_homo[j*nvoi+i] = stress_ave[j] / strain_ave[i];

      if( flag_print & ( 1 << PRINT_VTU ) && homo_type == UNIF_STRAINS )
      {
        get_elem_properties();
	sprintf(filename,"micro_exp%d",i);
	ierr = micro_pvtu( filename );
	if(ierr){
	  PetscPrintf(MICRO_COMM,"Problem writing vtu file\n");
	  goto end_mic_1;
	}
      }

    }
    PetscPrintf(MICRO_COMM,"\nConstitutive Average Tensor\n");
    for( i = 0 ; i < nvoi ; i++ ){
      for( j = 0 ; j < nvoi ; j++ )
	PetscPrintf(MICRO_COMM,"%e ",(fabs(c_homo[i*nvoi+j])>1.0)?c_homo[i*nvoi+j]:0.0);
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
    PetscPrintf(MICRO_COMM,"\nstrain_ave = ");
    for( j = 0 ; j < nvoi ; j++ )
      PetscPrintf(MICRO_COMM,"%e ",strain_ave[j]);
    PetscPrintf(MICRO_COMM,"\nstress_ave = ");
    for( j = 0 ; j < nvoi ; j++ )
      PetscPrintf(MICRO_COMM,"%e ",stress_ave[j]);
    for( i = 0 ; i < nvoi ; i++ ){
      stress_ave[i] = 0.0;
	for( j = 0 ; j < nvoi ; j++ )
	stress_ave[i] +=  c_homo[i*nvoi+j] * strain_mac[j];
    }
    PetscPrintf(MICRO_COMM,"\nstress_ave = ");
    for( j = 0 ; j < nvoi ; j++ )
      PetscPrintf(MICRO_COMM,"%e ",stress_ave[j]);
    PetscPrintf(MICRO_COMM," (c_homo*strain_mac)\n");

  }

  micro_print_info();

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
  for( i = 0 ; i < nvoi  ; i++ ){
    for( j = 0 ; j < npe*dim ; j++ )
      free(struct_bmat[i][j]);
    free(struct_bmat[i]);
  }
  free(struct_bmat);

  end_trace( MICRO_COMM );
  /**************************************************/

end_mic_1:

  if(!flag_coupling){
    PetscPrintf(MICRO_COMM,
	"--------------------------------------------------\n"
	"  MICRO: FINISH COMPLETE\n"
	"--------------------------------------------------\n");
  }

  ierr = PetscFinalize();

end_mic_2:
  ierr = MPI_Finalize();


  return 0;
}