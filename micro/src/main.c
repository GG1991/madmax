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
    printf_p(&MICRO_COMM,"micro: problem during coloring\n");
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

  /**************************************************/

  nx = ny = nz = -1;
  hx = hy = hz = -1;
  lx = ly = lz = -1;

  /**************************************************/

  /* dimension */

  PetscOptionsGetInt(NULL, NULL, "-dim", &dim, &set);
  if( set == PETSC_FALSE ){
    printf_p( &MICRO_COMM, "dimension (-dim <dim>) not given\n" );
    ierr_1 = 1;
    goto end_mic_1;
  }
  nvoi = (dim == 2) ? 3 : 6;

  /**************************************************/

  /* micro_struct */

  char format[256];

  PetscOptionsGetString( NULL, NULL, "-micro_struct", format, 256, &set );
  if( set == PETSC_FALSE ){
    printf_p( &MICRO_COMM, "micro structure ( -micro_struct <format> ) not given\n" );
    ierr_1 = 1;
    goto end_mic_1;
  }
  micro_struct_init( dim, format, &micro_struct );

  lx = micro_struct.size[0];
  ly = micro_struct.size[1];
  lz = ( dim == 3 ) ? micro_struct.size[2] : -1;

  /**************************************************/

  /* struct mesh */

  int    nval_expect;
  int    struct_mesh_n[3];

  if( dim == 2 ) nval_expect = nval = 2;
  if( dim == 3 ) nval_expect = nval = 3;

  PetscOptionsGetIntArray( NULL, NULL, "-struct_n", struct_mesh_n, &nval, &set );
  if( set == PETSC_TRUE ){

    if( nval != nval_expect ){
      printf_p(&MICRO_COMM,"-struct_n should include %d arguments\n", nval_expect);
      ierr_1 = 1;
      goto end_mic_0;
    }
    nx   = struct_mesh_n[0];
    ny   = struct_mesh_n[1];
    nz   = ( dim == 3 ) ? struct_mesh_n[2] : 1;

    nn   = nx*ny*nz;
    nex  = (nx-1);

    /* local number of elements in y direction */
    ney  = (ny-1)/nproc_mic + (((ny-1) % nproc_mic > rank_mic) ? 1:0); 
    nez  = (nz-1);
    nelm = ( dim == 2 ) ? nex*ney : nex*ney*nez;       // number of local elements
    nyl  = ( rank_mic == 0 ) ? ney+1 : ney;            // local number of nodes in y direction
    nl   = ( dim == 2 ) ? nyl*nx : nyl*nx*nz;          // local number of nodes

    /* set the elements' size */
    hx   = lx/nex;
    hy   = ly/(ny-1);
    hz   = ( dim == 3 ) ? (lz/nez) : -1;

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
      printf_p( &MICRO_COMM, "ny %d not large enough to be executed with %d processes\n", ny, nproc_mic);
      ierr_1 = 1;
      goto end_mic_0;
    }

  }
  else{
    printf_p(&MICRO_COMM,"-struct_n is request\n");
    ierr_1 = 1;
    goto end_mic_0;
  }


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

  /**************************************************/

  {
    /* Materials by command line */
    int    nval = 4;
    char   *string[nval];
    char   *data;
    material_t mat;
    list_init( &material_list, sizeof(material_t), NULL );

    PetscOptionsGetStringArray( NULL, NULL, "-material", string, &nval, &set );
    if( set == PETSC_TRUE )
    {
      for( i = 0 ; i < nval ; i++ )
      {
	data = strtok( string[i] , " \n" );
	mat.name = strdup( data );
	data = strtok( NULL , " \n" );
	if( strcmp( data, "TYPE_0" ) == 0 )
	{
	  double E, v;
	  mat.type_id = TYPE_0;
	  mat.type    = malloc(sizeof(type_0));
	  data = strtok( NULL , " \n" );
	  ((type_0*)mat.type)->rho         = atof(data);
	  data = strtok( NULL , " \n" );
	  E = ((type_0*)mat.type)->young   = atof(data);
	  data = strtok( NULL , " \n" );
	  v = ((type_0*)mat.type)->poisson = atof(data);
	  ((type_0*)mat.type)->lambda      = (E*v)/((1+v)*(1-2*v));
	  ((type_0*)mat.type)->mu          = E/(2*(1+v));
	}
	else if ( strcmp( data, "TYPE_1" ) == 0 )
	{
	  printf_p( MACRO_COMM, "TYPE_1 not allowed in micro code.\n" );
	  goto end_mic_1;
	}
	else
	{
	  printf_p( MACRO_COMM, "type %s not known.\n", data );
	  goto end_mic_1;
	}

	list_insertlast( &material_list , &mat );
      }
    }

  }

  /**************************************************/

  /* homogenization options */

  homo_type=0;
  PetscOptionsHasName(NULL,NULL,"-homo_tp",&set);
  if(set==PETSC_TRUE) homo_type = TAYLOR_P;
  PetscOptionsHasName(NULL,NULL,"-homo_ts",&set);
  if(set==PETSC_TRUE) homo_type = TAYLOR_S;
  PetscOptionsHasName(NULL,NULL,"-homo_us",&set);
  if(set==PETSC_TRUE) homo_type = UNIF_STRAINS;
  if(homo_type==0){
    printf_p(&MICRO_COMM,"no homogenization option specified\n");
    ierr_1 = 1;
    goto end_mic_0;
  }

  /**************************************************/

  /* solver options */

  PetscOptionsGetInt(NULL, NULL, "-nr_max_its", &nr_max_its, &set);
  if(set==PETSC_FALSE) nr_max_its = 5;
  PetscOptionsGetReal(NULL, NULL, "-nr_norm_tol", &nr_norm_tol, &set);
  if(set==PETSC_FALSE) nr_norm_tol = 1.0e-5;

  /**************************************************/

  /* printing options */

  flag_print = 0;
  PetscOptionsHasName(NULL,NULL,"-print_petsc",&set);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_PETSC);
  PetscOptionsHasName(NULL,NULL,"-print_vtu",&set);
  if(set == PETSC_TRUE) flag_print = flag_print | (1<<PRINT_VTU);

  /**************************************************/

  if(!flag_coupling){
    printf_p(&MICRO_COMM,
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
  }
  elem_type   = malloc( nelm      * sizeof(int));
  micro_struct_init_elem_type( &micro_struct, dim, nelm, &get_elem_centroid, elem_type );

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
	  printf_p(&MICRO_COMM,"MICRO:signal %d not identified\n",signal);
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

      printf_p(&MICRO_COMM,"\nstrain_ave = ");
      for( j = 0 ; j < nvoi ; j++ )
	printf_p(&MICRO_COMM,"%e ",strain_ave[j]);

      printf_p(&MICRO_COMM,"\nstress_ave = ");
      for( j = 0 ; j < nvoi ; j++ )
	printf_p(&MICRO_COMM,"%e ",stress_ave[j]);
      printf_p(&MICRO_COMM,"\n");

      for( j = 0 ; j < nvoi ; j++ )
	c_homo[j*nvoi+i] = stress_ave[j] / strain_ave[i];

      if( flag_print & ( 1 << PRINT_VTU ) && homo_type == UNIF_STRAINS )
      {
        get_elem_properties();
	sprintf(filename,"micro_exp%d",i);
	ierr = micro_pvtu( filename );
	if(ierr){
	  printf_p(&MICRO_COMM,"Problem writing vtu file\n");
	  goto end_mic_1;
	}
      }

    }
    printf_p(&MICRO_COMM,"\nConstitutive Average Tensor\n");
    for( i = 0 ; i < nvoi ; i++ ){
      for( j = 0 ; j < nvoi ; j++ )
	printf_p(&MICRO_COMM,"%e ",(fabs(c_homo[i*nvoi+j])>1.0)?c_homo[i*nvoi+j]:0.0);
      printf_p(&MICRO_COMM,"\n");
    }
    printf_p(&MICRO_COMM,"\n");

    /*
       Experiment to test if the homogenization with <strain_mac>
       gives the same <stress_ave> than doing <c_homo>*<strain_mac>
     */

    strain_mac[0] = 0.01; strain_mac[1] = -0.02; strain_mac[2] = +0.03;
    strain_mac[3] = 0.01; strain_mac[4] = -0.02; strain_mac[5] = +0.03;
    ierr = mic_homogenize(MICRO_COMM, strain_mac, strain_ave, stress_ave);
    printf_p(&MICRO_COMM,"\nstrain_ave = ");
    for( j = 0 ; j < nvoi ; j++ )
      printf_p(&MICRO_COMM,"%e ",strain_ave[j]);
    printf_p(&MICRO_COMM,"\nstress_ave = ");
    for( j = 0 ; j < nvoi ; j++ )
      printf_p(&MICRO_COMM,"%e ",stress_ave[j]);
    for( i = 0 ; i < nvoi ; i++ ){
      stress_ave[i] = 0.0;
	for( j = 0 ; j < nvoi ; j++ )
	stress_ave[i] +=  c_homo[i*nvoi+j] * strain_mac[j];
    }
    printf_p(&MICRO_COMM,"\nstress_ave = ");
    for( j = 0 ; j < nvoi ; j++ )
      printf_p(&MICRO_COMM,"%e ",stress_ave[j]);
    printf_p(&MICRO_COMM," (c_homo*strain_mac)\n");

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
    printf_p(&MICRO_COMM,
	"--------------------------------------------------\n"
	"  MICRO: FINISH COMPLETE\n"
	"--------------------------------------------------\n");
  }

  ierr = PetscFinalize();

end_mic_2:
  ierr = MPI_Finalize();


  return 0;
}

/****************************************************************************************************/

int micro_print_info( void )
{

  FILE *fm = fopen("mic_info.dat","w");

  int  *i_data;
  int   i , ierr;

  if( rank_mic == 0 ){
    fprintf(fm,"-----------\n");
    fprintf(fm,"nproc %d\n", nproc_mic);
    fprintf(fm,"-----------\n");
    fprintf(fm,"nx    ny    nz\n");
    fprintf(fm,"%2d    %2d    %2d\n", nx, ny, nz);
    fprintf(fm,"lx    ly    lz\n");
    fprintf(fm,"%lf    %lf    %lf\n", lx, ly, lz);
    fprintf(fm,"hx    hy    hz\n");
    fprintf(fm,"%lf    %lf    %lf\n", hx, hy, hz);
    fprintf(fm,"nex   ney   nez\n");
    fprintf(fm,"%2d    %2d    %2d\n", nex, ny-1, nez);
    fprintf(fm,"-----------\n");
  }

  if( rank_mic == 0 ){
    i_data = malloc( nproc_mic * sizeof(double) );
    fprintf(fm,"%-20s","rank ");
    for( i = 0 ; i < nproc_mic ; i++  )
      fprintf(fm,"%d ", i);
    fprintf(fm,"\n");
    fprintf(fm,"%-20s","");
    for( i = 0 ; i < nproc_mic ; i++ )
      fprintf(fm,"--");
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&nelm, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return 1;

  if( rank_mic == 0 ){
    fprintf(fm,"%-20s","nelm");
    for( i = 0 ; i < nproc_mic ; i++ )
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&ney, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return 1;

  if( rank_mic == 0 ){
    fprintf(fm,"%-20s","eyl");
    for( i = 0 ; i < nproc_mic ; i++ )
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&nyl, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return 1;

  if( rank_mic == 0 ){
    fprintf(fm,"%-20s","nyl");
    for( i = 0 ; i < nproc_mic ; i++ )
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&ny_inf, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return 1;

  if( rank_mic == 0 ){
    fprintf(fm,"%-20s","ny_inf");
    for( i = 0 ; i < nproc_mic ; i++ )
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&ngho, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return 1;

  if( rank_mic == 0 ){
    fprintf(fm,"%-20s","ngho");
    for( i = 0 ; i < nproc_mic ; i++ )
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  double norm;
 
  if( b != NULL ){
    VecNorm( b, NORM_2, &norm );
    if( rank_mic == 0 )
      fprintf(fm,"|b| = %lf \n", norm);
  }

  if( x != NULL ){
    VecNorm( x, NORM_2, &norm );
    if( rank_mic == 0 )
      fprintf(fm,"|x| = %lf \n", norm);
  }

  if( A != NULL ){
    MatNorm( A , NORM_FROBENIUS , &norm );
    if( rank_mic == 0 )
      fprintf(fm,"|A| = %lf \n",norm);
  }

  fclose(fm);
  return 0;
}

/****************************************************************************************************/

int micro_pvtu( char *name )
{

  FILE    *fm;
  char     file_name[NBUF];
  double  *xvalues;
  Vec      xlocal;

  if( rank_mic == 0 ){

    /* rank 0 writes the .pvtu file first */
    strcpy(file_name,name);
    strcat(file_name,".pvtu");
    fm = fopen(file_name,"w");

    fprintf(fm, "<?xml version=\"1.0\"?>\n"
	"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
	"<PUnstructuredGrid GhostLevel=\"0\">\n"
	"<PPoints>\n"
	"<PDataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\"/>\n"
	"</PPoints>\n"
	"<PCells>\n"
	"<PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>\n"
	"</PCells>\n" 

	"<PPointData Vectors=\"displ\">\n" 
	"<PDataArray type=\"Float64\" Name=\"displ\"    NumberOfComponents=\"3\" />\n"
	"<PDataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" />\n"
	"</PPointData>\n"

	"<PCellData>\n"
	"<PDataArray type=\"Int32\"   Name=\"part\"   NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\"/>\n"
	"<PDataArray type=\"Int32\"   Name=\"elem_type\" NumberOfComponents=\"1\"/>\n"
	"</PCellData>\n" , nvoi , nvoi);

    int i;
    for( i = 0 ; i < nproc_mic ; i++ ){
      sprintf(file_name,"%s_%d",name,i);
      fprintf(fm,	"<Piece Source=\"%s.vtu\"/>\n",file_name);
    }
    fprintf(fm,	"</PUnstructuredGrid>\n" 
      "</VTKFile>\n" );

    fclose(fm);

  } // rank = 0

  sprintf( file_name, "%s_%d.vtu", name, rank_mic);
  fm = fopen(file_name,"w"); 
  if(!fm){
    printf_p(PETSC_COMM_WORLD,"Problem trying to opening file %s for writing\n", file_name);
    return 1;
  }

  fprintf(fm, 
      "<?xml version=\"1.0\"?>\n"
      "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      "<UnstructuredGrid>\n");
  fprintf(fm,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nl + ngho, nelm);
  fprintf(fm,"<Points>\n");
  fprintf(fm,"<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");

  double *coord = malloc( dim * sizeof(double));
  int    n , d;

  for( n = 0 ; n < nl ; n++ ){
    get_node_local_coor( n , coord );
    for( d = 0 ; d < dim ; d++ )
      fprintf(fm,"%e ",  coord[d] );
    for( d = dim ; d < 3 ; d++ )
      fprintf(fm,"%e ",0.0);
    fprintf(fm,"\n");
  }
  for( n = 0 ; n < ngho ; n++ ){
    get_node_ghost_coor( n , coord );
    for( d = 0 ; d < dim ; d++ )
      fprintf(fm,"%e ",  coord[d] );
    for( d = dim ; d < 3 ; d++ )
      fprintf(fm,"%e ",0.0);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");
  fprintf(fm,"</Points>\n");
  fprintf(fm,"<Cells>\n");

  int e;

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for ( e = 0 ; e < nelm ; e++ ){
    get_local_elem_node( e , loc_elem_index );
    for ( n = 0 ; n < npe ; n++ )
      fprintf(fm,"%d ", loc_elem_index[n]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  int ce = npe;
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for ( e = 0 ; e < nelm ; e++ ){
    fprintf(fm,"%d ", ce);
    ce += npe;
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for ( e = 0 ; e < nelm ; e++ )
    fprintf(fm, "%d ",vtkcode( dim , npe ) );  
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"</Cells>\n");
  
  fprintf(fm,"<PointData Vectors=\"displ\">\n"); // Vectors inside is a filter we should not use this here

  /* <displ> */
  VecGhostUpdateBegin( x , INSERT_VALUES , SCATTER_FORWARD );
  VecGhostUpdateEnd(   x , INSERT_VALUES , SCATTER_FORWARD );
  VecGhostGetLocalForm(x , &xlocal );

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
  VecGetArray( xlocal , &xvalues );
  for( n = 0 ; n < (nl + ngho) ; n++ ){
    for( d = 0 ; d < dim ; d++ )
      fprintf(fm, "%lf ", xvalues[ n * dim + d ]);
    for( d = dim ; d < 3 ; d++ )
      fprintf(fm,"%lf ",0.0);
    fprintf(fm,"\n");
  }
  VecRestoreArray( xlocal , &xvalues );
  fprintf(fm,"</DataArray>\n");

  /* <residual> */
  VecGhostUpdateBegin( b , INSERT_VALUES, SCATTER_FORWARD );
  VecGhostUpdateEnd  ( b , INSERT_VALUES, SCATTER_FORWARD );
  VecGhostGetLocalForm(b , &xlocal);

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
  VecGetArray(xlocal, &xvalues);
  for( n = 0 ; n < (nl + ngho) ; n++ ){
    for( d = 0 ; d < dim ; d++ )
      fprintf(fm, "%lf ", xvalues[ n * dim + d ]);
    for( d = dim ; d < 3 ; d++ )
      fprintf(fm, "%lf ", 0.0);
    fprintf(fm,"\n");
  }
  VecRestoreArray( xlocal , &xvalues );
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"</PointData>\n");
  fprintf(fm,"<CellData>\n");

  /* <part> */
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"part\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for( e = 0; e < nelm ; e++ )
    fprintf( fm, "%d ", rank_mic );  
  fprintf( fm, "\n");
  fprintf( fm, "</DataArray>\n");

  int v;

  /* <strain> */
  fprintf(fm,"<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for( e = 0; e < nelm ; e++ ){
    for( v = 0 ; v < nvoi ; v++ )
      fprintf(fm, "%lf ", elem_strain[ e*nvoi + v ]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  /* <stress> */
  fprintf(fm,"<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for( e = 0; e < nelm ; e++ ){
    for( v = 0 ; v < nvoi ; v++ )
      fprintf(fm, "%lf ", elem_stress[ e*nvoi + v ]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  /* <elem_type> */
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"elem_type\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for( e = 0; e < nelm ; e++ )
    fprintf( fm, "%d ", elem_type[e] );
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  /* <energy> */
//  fprintf(fm,"<DataArray type=\"Float64\" Name=\"energy\" NumberOfComponents=\"1\" format=\"ascii\">\n");
//  for (i=0;i<nelm;i++){
//    fprintf(fm, "%lf ", energy[i]);
//  }
//  fprintf(fm,"\n");
//  fprintf(fm,"</DataArray>\n");

  fprintf(fm,
      "</CellData>\n"
      "</Piece>\n"
      "</UnstructuredGrid>\n"
      "</VTKFile>\n");

  fclose(fm);
  return 0;
}

/****************************************************************************************************/
