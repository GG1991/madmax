#include "micro.h"

static char help[] = 
"micro multiscale code\n"
"-coupl    [0 (no coupling ) | 1 (coupling with micro)]                                       \n"
"-homo_ts     : c =  vi ci + vm cm            (serial)                                        \n"
"-homo_tp     : c = (vi ci^-1 + vm cm^-1)^-1  (parallel)                                      \n"
"-homo_us     : homogenization using uniform strains approach                                 \n"
"-struct_n [<nx,ny>] if dim = 2                                                               \n"
"-print_petsc                                                                                 \n"
"-print_vtu                                                                                   \n";

params_t params;

#define CHECK_AND_GOTO(error){if(error){myio_printf(&MICRO_COMM, "error line %d at %s\n", __LINE__, __FILE__); goto end;}}
#define CHECK_INST_ELSE_GOTO(cond, instr){if(cond){instr}else{myio_printf(&MICRO_COMM, "error line %d at %s\n", __LINE__, __FILE__); goto end;}}
#define PRINT_ARRAY(name_str, array, length){ myio_printf(&MICRO_COMM,"%s ",name_str); for(int i = 0 ; i < length ; i++) myio_printf(&MICRO_COMM, "%e ", array[i]); myio_printf(&MICRO_COMM, "\n");}

int main(int argc, char **argv){

  int i, j, ierr;
  int nval;
  bool found;

  myio_comm_line_init(argc, argv, &command_line);

  WORLD_COMM = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(WORLD_COMM, &nproc_wor);
  MPI_Comm_rank(WORLD_COMM, &rank_wor);

  init_variables(&params, &message);

  myio_comm_line_search_option(&command_line, "-coupl", &found);
  if(found) params.flag_coupling = true;

  macmic.type = COUP_1;
  color = COLOR_MICRO; /* color can change */
  ierr = macmic_coloring(WORLD_COMM, &color, &macmic, &MICRO_COMM, params.flag_coupling); CHECK_AND_GOTO(ierr);

  MPI_Comm_size(MICRO_COMM, &nproc_mic);
  MPI_Comm_rank(MICRO_COMM, &rank_mic);
  
  PETSC_COMM_WORLD = MICRO_COMM;
  ierr = PetscInitialize(&argc,&argv,(char*)0,help);

  nx = ny = nz = -1;
  hx = hy = hz = -1;
  lx = ly = lz = -1;

  myio_comm_line_get_int(&command_line, "-dim");
  CHECK_INST_ELSE_GOTO(command_line.found, dim = command_line.int_val;)

  nvoi = (dim == 2) ? 3 : 6;

  myio_comm_line_get_string(&command_line, "-micro_struct");
  CHECK_INST_ELSE_GOTO(command_line.found, micro_struct_init(dim, command_line.str, &micro_struct);)
  
  lx = micro_struct.size[0];
  ly = micro_struct.size[1];
  lz = (dim == 3) ? micro_struct.size[2] : -1;

  if(dim == 2) nval = 2;
  if(dim == 3) nval = 3;

  
  myio_comm_line_get_int_array(&command_line, nval, "-struct_n");
  if(command_line.found){

    if(nval != command_line.n_int_found){
      myio_printf(&MICRO_COMM,"-struct_n should include %d arguments\n", nval);
      goto end;
    }
    nx   = command_line.int_arr[0];
    ny   = command_line.int_arr[1];
    nz   = ( dim == 3 ) ? command_line.int_arr[2] : 1;

    nn   = nx*ny*nz;
    nex  = (nx-1);

    ney  = (ny-1)/nproc_mic + (((ny-1) % nproc_mic > rank_mic) ? 1:0); 
    nez  = (nz-1);
    nelm = ( dim == 2 ) ? nex*ney : nex*ney*nez;
    nyl  = ( rank_mic == 0 ) ? ney+1 : ney;
    nl   = ( dim == 2 ) ? nyl*nx : nyl*nx*nz;

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
      myio_printf( &MICRO_COMM, "ny %d not large enough to be executed with %d processes\n", ny, nproc_mic);
      goto end;
    }

  }
  else{
    myio_printf(&MICRO_COMM,"-struct_n is request\n");
    goto end;
  }

  center_coor = malloc ( dim * sizeof(double));
  if(dim == 2){
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

  ierr = material_fill_list_from_command_line(&command_line, &material_list); CHECK_AND_GOTO(ierr)

  myio_comm_line_search_option(&command_line, "-homo_us", &found);
  if(found) params.homog_method = HOMOG_METHOD_UNIF_STRAINS;

  myio_comm_line_search_option(&command_line, "-homo_tp", &found);
  if(found) params.homog_method = HOMOG_METHOD_TAYLOR_PARALLEL;

  myio_comm_line_search_option(&command_line, "-homo_ts", &found);
  if(found) params.homog_method = HOMOG_METHOD_TAYLOR_SERIAL;

  myio_comm_line_get_int(&command_line, "-nl_max_its");
  if(command_line.found) params.non_linear_max_its = command_line.int_val;

  myio_comm_line_get_int(&command_line, "-nl_min_norm_tol");
  if(command_line.found) params.non_linear_min_norm_tol = command_line.double_val;

  flag_print = 0;
  myio_comm_line_search_option(&command_line, "-print_petsc", &found);
  if(found) flag_print = flag_print | (1<<PRINT_PETSC);

  myio_comm_line_search_option(&command_line, "-print_vtu", &found);
  if(found) flag_print = flag_print | (1<<PRINT_VTU);

  if(params.flag_coupling == false){
    myio_printf(&MICRO_COMM,
	"--------------------------------------------------\n"
	"  MICRO: STANDALONE \n"
	"--------------------------------------------------\n");
  }

  A   = NULL;
  b   = NULL;
  x   = NULL;
  dx  = NULL;

  loc_elem_index = malloc(dim*npe*sizeof(int));
  glo_elem_index = malloc(dim*npe*sizeof(int));
  elem_disp      = malloc(dim*npe*sizeof(double));
  stress_gp      = malloc(nvoi*sizeof(double));
  strain_gp      = malloc(nvoi*sizeof(double));
  if( flag_print & ( 1 << PRINT_VTU ) ){
    elem_strain = malloc(nelm*nvoi*sizeof(double));
    elem_stress = malloc(nelm*nvoi*sizeof(double));
    elem_energy = malloc(nelm*sizeof(double));
  }
  elem_type   = malloc( nelm      * sizeof(int));

  ierr = micro_struct_init_elem_type(&micro_struct, dim, nelm, &get_elem_centroid, elem_type); CHECK_AND_GOTO(ierr)

  ierr = micro_check_material_and_elem_type(&material_list, elem_type, nelm);
  if(ierr){
    myio_printf(&MICRO_COMM, "error checking elem_type and material_list\n");
    goto end;
  }


  init_shapes( &struct_sh, &struct_dsh, &struct_wp );

  struct_bmat = malloc( nvoi * sizeof(double**));
  for( i = 0 ; i < nvoi ; i++ ){
    struct_bmat[i] = malloc( npe*dim * sizeof(double*));
    for( j = 0 ; j < npe*dim ; j++ )
      struct_bmat[i][j] = malloc( ngp * sizeof(double));
  }

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

  init_trace(MICRO_COMM, "micro_trace.dat");

  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  homogenize_init();

  double strain_mac[6], strain_ave[6], stress_ave[6];

  double *c_tangent_ave = malloc(36*sizeof(double));

  if(params.flag_coupling == true){

    while(message.action != ACTION_MICRO_END){

      ierr = comm_micro_recv(&message);

      switch(message.action){

	case ACTION_MICRO_CALC_STRESS:

	  ARRAY_COPY(strain_mac, message.strain_mac, nvoi)

	  ierr = homogenize_get_strain_stress(strain_mac, strain_ave, stress_ave);

	  ARRAY_COPY(message.stress_ave, stress_ave, nvoi)

	  break;

	case ACTION_MICRO_CALC_C_TANGENT:

	  ARRAY_COPY(strain_mac, message.strain_mac, nvoi)

	  ierr = homogenize_get_c_tangent(strain_mac, &c_tangent_ave);

	  ARRAY_COPY(message.c_tangent_ave, c_tangent_ave, nvoi*nvoi)

	  break;

	case ACTION_MICRO_END:
	  break;

	default:
	  myio_printf(&MICRO_COMM, "MICRO:signal %d not identified\n", signal);
	  goto end;

      }

      ierr = comm_micro_send(&message);
    }
  }
  else{

    double strain_mac[6], c_homo[36];

    ARRAY_SET_TO_ZERO(c_homo, 36)

    for(int i = 0 ; i < nvoi ; i++){

      ARRAY_SET_TO_ZERO(strain_mac, nvoi); strain_mac[i]=0.005;

      ierr = homogenize_get_strain_stress_non_linear(strain_mac, strain_ave, stress_ave);

      PRINT_ARRAY("strain_ave = ", strain_ave, nvoi)
      PRINT_ARRAY("stress_ave = ", stress_ave, nvoi)

      for(int j = 0 ; j < nvoi ; j++ )
	c_homo[j*nvoi+i] = stress_ave[j] / strain_ave[i];

      if(flag_print & (1 << PRINT_VTU) && params.homog_method == HOMOG_METHOD_UNIF_STRAINS){
        get_elem_properties();
	sprintf(filename,"micro_exp%d",i);
	ierr = micro_pvtu( filename );
	if(ierr){
	  myio_printf(&MICRO_COMM,"Problem writing vtu file\n");
	  goto end;
	}
      }

    }

    myio_printf(&MICRO_COMM,"\nConstitutive Average Tensor\n");
    for(int i = 0 ; i < nvoi ; i++){
      for(int j = 0 ; j < nvoi ; j++)
	myio_printf(&MICRO_COMM, "%e ", (fabs(c_homo[i*nvoi+j])>1.0) ? c_homo[i*nvoi+j] : 0.0);
      myio_printf(&MICRO_COMM, "\n");
    }
    myio_printf(&MICRO_COMM, "\n");

    strain_mac[0] = 0.01; strain_mac[1] = -0.02; strain_mac[2] = +0.03;
    strain_mac[3] = 0.01; strain_mac[4] = -0.02; strain_mac[5] = +0.03;

    ierr = homogenize_get_strain_stress_non_linear(strain_mac, strain_ave, stress_ave);

    PRINT_ARRAY("strain_ave = ", strain_ave, nvoi)
    PRINT_ARRAY("stress_ave = ", stress_ave, nvoi)

    for( i = 0 ; i < nvoi ; i++ ){
      stress_ave[i] = 0.0;
      for( j = 0 ; j < nvoi ; j++ )
	stress_ave[i] +=  c_homo[i*nvoi+j] * strain_mac[j];
    }

    PRINT_ARRAY("stress_ave (c_homo*strain_mac) = ", stress_ave, nvoi)

  }

  micro_print_info();

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

  for( i = 0 ; i < nvoi  ; i++ ){
    for( j = 0 ; j < npe*dim ; j++ )
      free(struct_bmat[i][j]);
    free(struct_bmat[i]);
  }
  free(struct_bmat);

  end_trace( MICRO_COMM );

  PRINTF1("--------------------------------------------------\n"
      "  MICRO: FINISH COMPLETE\n"
      "--------------------------------------------------\n");

  ierr = PetscFinalize();

  finalize();

end:
  ierr = MPI_Finalize();


  return 0;
}


int micro_print_info( void ){

  FILE *fm = fopen("mic_info.dat","w");

  int  *i_data;
  int   i , ierr;

  if(rank_mic == 0){
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

  if(rank_mic == 0){
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

  sprintf(file_name,"%s_%d.vtu",name,rank_mic);
  fm = fopen(file_name,"w"); 
  if(!fm){
    myio_printf(&MICRO_COMM,"Problem trying to opening file %s for writing\n", file_name);
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

  for(int n = 0 ; n < nl ; n++){
    get_node_local_coor( n , coord);
    for(int d = 0 ; d < dim ; d++)
      fprintf(fm,"%e ",  coord[d]);
    for(int d = dim ; d < 3 ; d++)
      fprintf(fm,"%e ",0.0);
    fprintf(fm,"\n");
  }
  for(int n = 0 ; n < ngho ; n++ ){
    get_node_ghost_coor( n , coord );
    for(int d = 0 ; d < dim ; d++ )
      fprintf(fm,"%e ",  coord[d] );
    for(int d = dim ; d < 3 ; d++ )
      fprintf(fm,"%e ",0.0);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");
  fprintf(fm,"</Points>\n");
  fprintf(fm,"<Cells>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int e = 0 ; e < nelm ; e++){
    get_local_elem_node( e , loc_elem_index );
    for(int n = 0 ; n < npe ; n++)
      fprintf(fm,"%d ", loc_elem_index[n]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  int ce = npe;
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int e = 0 ; e < nelm ; e++){
    fprintf(fm,"%d ", ce);
    ce += npe;
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int e = 0 ; e < nelm ; e++)
    fprintf(fm, "%d ",vtkcode(dim, npe));  
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"</Cells>\n");
  
  fprintf(fm,"<PointData Vectors=\"displ\">\n"); // Vectors inside is a filter we should not use this here

  VecGhostUpdateBegin(x , INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(x , INSERT_VALUES, SCATTER_FORWARD);
  VecGhostGetLocalForm(x , &xlocal);

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
  VecGetArray( xlocal , &xvalues );
  for(int n = 0 ; n < (nl + ngho) ; n++){
    for(int d = 0 ; d < dim ; d++)
      fprintf(fm, "%lf ", xvalues[ n * dim + d ]);
    for(int d = dim ; d < 3 ; d++)
      fprintf(fm,"%lf ",0.0);
    fprintf(fm,"\n");
  }
  VecRestoreArray( xlocal , &xvalues );
  fprintf(fm,"</DataArray>\n");

  VecGhostUpdateBegin( b , INSERT_VALUES, SCATTER_FORWARD );
  VecGhostUpdateEnd  ( b , INSERT_VALUES, SCATTER_FORWARD );
  VecGhostGetLocalForm(b , &xlocal);

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
  VecGetArray(xlocal, &xvalues);
  for(int n = 0 ; n < (nl + ngho) ; n++){
    for(int d = 0 ; d < dim ; d++)
      fprintf(fm, "%lf ", xvalues[ n * dim + d ]);
    for(int d = dim ; d < 3 ; d++)
      fprintf(fm, "%lf ", 0.0);
    fprintf(fm,"\n");
  }
  VecRestoreArray( xlocal , &xvalues );
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"</PointData>\n");
  fprintf(fm,"<CellData>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"part\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int e = 0; e < nelm ; e++ )
    fprintf( fm, "%d ", rank_mic );  
  fprintf( fm, "\n");
  fprintf( fm, "</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for(int e = 0; e < nelm ; e++){
    for(int v = 0 ; v < nvoi ; v++)
      fprintf(fm, "%lf ", elem_strain[ e*nvoi + v ]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for(int e = 0; e < nelm ; e++){
    for(int v = 0 ; v < nvoi ; v++)
      fprintf(fm, "%lf ", elem_stress[ e*nvoi + v ]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"elem_type\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int e = 0; e < nelm ; e++)
    fprintf(fm, "%d ", elem_type[e]);
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,
      "</CellData>\n""</Piece>\n"
      "</UnstructuredGrid>\n"
      "</VTKFile>\n");

  fclose(fm);

  return 0;
}


int micro_check_material_and_elem_type(list_t *material_list, int *elem_type, int nelm){

  char *word_to_search;

  for(int e ; e < nelm ; e++ ){

    switch(elem_type[e]){

      case ID_FIBER:
	word_to_search = strdup("FIBER");
	break;

      case ID_MATRIX:
	word_to_search = strdup("MATRIX");
	break;

      default:
	return 1;
    }

    material_t  *mat_p;
    node_list_t *pm = material_list->head;

    while(pm != NULL){
      mat_p = (material_t *)pm->data;
      if(strcmp(mat_p->name, word_to_search) == 0) break;
      pm = pm->next;
    }

    if(pm == NULL) return 1;
  }

  return 0;
}
