#include "macro.h"

static char help[] = 
"macro multiscale code                                                                             \n"
"-coupl       : coupled with \"micro\" code for solving multiscale problem                         \n"
"-normal      : normal execution, solves a time dependent boundary condition problem               \n"
"-testcomm    : communication testing with the \"micro\" code                                      \n"
"-eigen       : calculates the eigensystem Mx = -(1/omega)Kx                                       \n"
"-print_petsc : prints petsc structures on files such as Mat and Vec objects                       \n"
"-print_vtu   : prints solutions on .vtu and .pvtu files                                           \n";

params_t params;
flags_t flags;
solver_t solver;

#define CHECK_AND_GOTO(error){if(error){myio_printf(&MACRO_COMM, "error line %d at %s\n", __LINE__, __FILE__);goto end;}}
#define CHECK_INST_ELSE_GOTO(cond,instr){if(cond){instr}else{myio_printf(&MACRO_COMM, "error line %d at %s\n", __LINE__, __FILE__); goto end;}}
#define CHECK_FOUND_GOTO(message){if(found == false){myio_printf(&MACRO_COMM, "%s\n", message); goto end;}}
#define CHECK_ERROR_GOTO(message){if(ierr != 0){myio_printf(&MACRO_COMM, "%s\n", message); goto end;}}

int main(int argc, char **argv){

  int ierr;
  bool found;

  myname = strdup("macro");

  myio_comm_line_init(argc, argv, &command_line);

  WORLD_COMM = MPI_COMM_WORLD;
  MPI_Init( &argc, &argv );
  MPI_Comm_size( WORLD_COMM, &nproc_wor );
  MPI_Comm_rank( WORLD_COMM, &rank_wor );

  init_variables();

  myio_comm_line_search_option(&command_line, "-coupl", &found);
  if(found == true) flags.coupled = true;

  macmic.type = COUP_1;
  color = COLOR_MACRO;
  ierr = macmic_coloring(WORLD_COMM, &color, &macmic, &MACRO_COMM, flags.coupled); CHECK_AND_GOTO(ierr);

  MPI_Comm_size(MACRO_COMM, &nproc_mac);
  MPI_Comm_rank(MACRO_COMM, &rank_mac);

  myio_printf(&MACRO_COMM,
      "--------------------------------------------------\n"
      "  MACRO: COMPOSITE MATERIAL MULTISCALE CODE\n"
      "--------------------------------------------------\n");

  myio_comm_line_search_option(&command_line, "-help", &found);
  if(found == true){
    myio_printf(&MACRO_COMM, "%s", help);
    goto end;
  }

  myio_comm_line_search_option(&command_line, "-normal", &found);
  if(found == true){
    params.calc_mode = CALC_MODE_NORMAL;

    myio_comm_line_get_double(&command_line, "-tf", &params.tf, &found);
    CHECK_FOUND_GOTO("-tf option should be given in -normal mode.\n");

    myio_comm_line_get_double(&command_line, "-dt", &params.dt, &found);
    CHECK_FOUND_GOTO("-dt option should be given in -normal mode.\n");
  }

  myio_comm_line_search_option(&command_line, "-testcomm", &found);
  if(found == true){
    params.calc_mode = CALC_MODE_TEST;
  }

  myio_comm_line_search_option(&command_line, "-eigen", &found);
  if(found == true){
    params.calc_mode = CALC_MODE_EIGEN;
    myio_comm_line_get_double(&command_line, "-energy_stored", &params.energy_stored, &found);
  }

  mesh_f = FORMAT_GMSH;

  myio_comm_line_get_string(&command_line, "-mesh", mesh_n, &found);
  if(found == false){
    myio_printf(&MACRO_COMM,"mesh file not given on command line.\n");
    goto end;
  }

  FILE *fm = fopen(mesh_n, "r");
  if(fm == NULL){
    myio_printf(&MACRO_COMM,"mesh file not found.\n");
    goto end;
  }
  myio_comm_line_get_int(&command_line, "-dim", &dim, &found);
  if(found == false){
    myio_printf(&MACRO_COMM,"-dim not given on command line.\n");
    goto end;
  }

  nvoi    = (dim == 2) ? 3 : 6;
  npe_max = (dim == 2) ? 4 : 8;
  ngp_max = npe_max;

  myio_comm_line_search_option(&command_line, "-print_matrices", &found);
  if(found == true)
    flags.print_matrices = true;

  myio_comm_line_search_option(&command_line, "-print_vectors", &found);
  if(found == true)
    flags.print_vectors = true;

  myio_comm_line_search_option(&command_line, "-print_pvtu", &found);
  if(found == true)
    flags.print_pvtu = true;

  myio_comm_line_get_int(&command_line, "-nl_max_its", &params.non_linear_max_its, &found);

  myio_comm_line_get_double(&command_line, "-nl_min_norm_tol", &params.non_linear_min_norm_tol, &found);

  ierr = function_fill_list_from_command_line(&command_line, &function_list); CHECK_AND_GOTO(ierr);

  ierr = mesh_fill_boundary_list_from_command_line(&command_line, &boundary_list); CHECK_AND_GOTO(ierr);

  ierr = material_fill_list_from_command_line(&command_line, &material_list); CHECK_AND_GOTO(ierr);

  partition_algorithm = PARMETIS_GEOM;
  myio_comm_line_search_option(&command_line, "-part_kway", &found);
  if(found) partition_algorithm = PARMETIS_MESHKWAY;

  myio_comm_line_search_option(&command_line, "-part_geom", &found);
  if(found) partition_algorithm = PARMETIS_GEOM;

  myio_printf(&MACRO_COMM, "reading mesh elements\n" );
  ierr = read_mesh_elmv(MACRO_COMM, myname, mesh_n, mesh_f);
  CHECK_AND_GOTO(ierr);

  myio_printf(&MACRO_COMM, "partitioning and distributing mesh\n");
  ierr = part_mesh(MACRO_COMM, myname, NULL); CHECK_AND_GOTO(ierr);

  myio_printf(&MACRO_COMM, "calculating ghost nodes\n");
  ierr = calc_local_and_ghost(MACRO_COMM, nallnods, allnods, &ntotnod, &nmynods, &mynods, &nghost, &ghost ); CHECK_AND_GOTO(ierr);

  myio_printf(&MACRO_COMM, "reenumering nodes\n");
  ierr = reenumerate_PETSc( MACRO_COMM ); CHECK_AND_GOTO(ierr);

  myio_printf(&MACRO_COMM, "reading Coordinates\n");
  ierr = read_coord(mesh_n, nmynods, mynods, nghost , ghost, &coord ); CHECK_AND_GOTO(ierr);

  ierr = read_bc(); CHECK_AND_GOTO(ierr);

  list_init(&physical_list, sizeof(physical_t), NULL );
  gmsh_get_physical_list(mesh_n, &physical_list);

  myio_printf(&MACRO_COMM, "allocating ");
  ierr = alloc_memory();

  ierr = fem_inigau();

  if(params.calc_mode == CALC_MODE_EIGEN){

    EPS eps;

    VecZeroEntries(x);

    VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

    ierr = assembly_AM_petsc();
    CHECK_ERROR_GOTO("problem during matrix assembly\n")

    int nconv;
    double error;

    EPSCreate(MACRO_COMM, &eps);
    EPSSetOperators(eps, M, A);
    EPSSetProblemType(eps, EPS_GHEP);
    EPSSetFromOptions(eps);
    EPSGetDimensions(eps, &params.num_eigen_vals, NULL, NULL);
    params.eigen_vals = malloc( params.num_eigen_vals*sizeof(double));
    myio_printf(&MACRO_COMM,"Number of requested eigenvalues: %d\n", params.num_eigen_vals);

    EPSSolve(eps);
    EPSGetConverged(eps, &nconv);
    myio_printf(&MACRO_COMM,"Number of converged eigenpairs: %d\n",nconv);

    for(int i = 0 ; i < params.num_eigen_vals ; i++){

      EPSGetEigenpair( eps, i, &params.eigen_vals[i], NULL, x, NULL );
      EPSComputeError( eps, i, EPS_ERROR_RELATIVE, &error );
      myio_printf(&MACRO_COMM, "omega %d = %e   error = %e\n", i, params.eigen_vals[i], error);

      if(flags.print_pvtu == true){
	get_elem_properties();
	sprintf( filename, "macro_eigen_%d", i);
	macro_pvtu( filename );
      }

    }

    EPSDestroy(&eps);

  }
  else if(params.calc_mode == CALC_MODE_NORMAL){

    KSP ksp;

    KSPCreate(MACRO_COMM, &ksp);
    KSPSetFromOptions( ksp );

    params.t = 0.0;
    params.ts = 0;

    VecZeroEntries(x);
    VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

    while(params.t < (params.tf + 1.0e-10)){

      myio_printf(&MACRO_COMM,"\ntime step %-3d %-e seg\n", params.ts, params.t);

      update_boundary(params.t, &function_list, &boundary_list);

      Vec x_loc;
      double *x_arr;

      VecGhostGetLocalForm(x, &x_loc);
      VecGetArray(x_loc, &x_arr);

      node_list_t * pn = boundary_list.head;
      while(pn){
	mesh_boundary_t *bou = (mesh_boundary_t *)pn->data;
	for(int i = 0 ; i < bou->ndirix ; i++)
	  x_arr[bou->dir_loc_ixs[i]] = bou->dir_val[i];
	pn = pn->next;
      }

      VecRestoreArray(x_loc, &x_arr);
      VecGhostRestoreLocalForm(x, &x_loc);

      VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
      VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

      params.non_linear_its = 0; params.residual_norm = 2*params.non_linear_min_norm_tol;

      while(params.non_linear_its < params.non_linear_max_its && params.residual_norm > params.non_linear_min_norm_tol){

	myio_printf(&MACRO_COMM, "MACRO: assembling residual\n" );
	assembly_b_petsc();

	VecNorm(b, NORM_2, &params.residual_norm);

	myio_printf(&MACRO_COMM,"MACRO: |b| = %e\n", params.residual_norm );
	if(params.residual_norm < params.non_linear_min_norm_tol) break;

	myio_printf(&MACRO_COMM, "MACRO: assembling jacobian\n");
	assembly_A_petsc();

	myio_printf(&MACRO_COMM, "MACRO: solving system\n" );
	KSPSetOperators(ksp, A, A);
	KSPSolve(ksp, b, dx);
	print_ksp_info( MACRO_COMM, ksp);
	myio_printf(&MACRO_COMM, "\n");

	VecAXPY(x, 1.0, dx);
	VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
	VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

	params.non_linear_its ++;
      }

      if(flags.print_pvtu == true){
	get_elem_properties();
	sprintf(filename, "macro_t_%d", params.ts);
	macro_pvtu(filename);
      }

      params.t += params.dt;
      params.ts ++;
    }
    KSPDestroy(&ksp);

  }else if(params.calc_mode == CALC_MODE_TEST){


  }

end:

  myio_printf(&MACRO_COMM,
      "--------------------------------------------------\n"
      "  MACRO: FINISH COMPLETE\n"
      "--------------------------------------------------\n");

  finalize();

  return 0;
}


int read_bc(){

  mesh_boundary_t    *bou;

  node_list_t *pn = boundary_list.head;
  while(pn){
    int *ix, n;
    bou = (mesh_boundary_t *)pn->data;
    int ierr = gmsh_get_node_index(mesh_n, bou->name, nmynods, mynods, dim, &n, &ix);
    if(ierr){
      myio_printf(&MACRO_COMM, "problem finding nodes of boundary %s on msh file\n", bou->name );
      return 1;
    }
    bou->ndir        = n;
    bou->ndirix      = bou->ndir * bou->ndirpn;
    bou->dir_val     = malloc( bou->ndirix * sizeof(double));
    bou->dir_loc_ixs = malloc( bou->ndirix * sizeof(int));
    bou->dir_glo_ixs = malloc( bou->ndirix * sizeof(int));
    for(int i = 0 ; i < n ; i++ ){
      int da = 0;
      int * p = bsearch( &ix[i], mynods, nmynods, sizeof(int), cmpfunc );
      for(int d = 0 ; d < dim ; d++ )
	if( bou->kind & (1<<d) ) {
	  bou->dir_loc_ixs[i* (bou->ndirpn) + da] = (p - mynods) * dim + d;
	  bou->dir_glo_ixs[i* (bou->ndirpn) + da] = loc2petsc[(p - mynods)] * dim + d;
	  da++;
	}
    }
    free(ix);
    pn = pn->next;
  }

  return 0;
}


int read_coord( char *mesh_n, int nmynods, int *mynods, int nghost , int *ghost, double **coord ){

  (*coord) = malloc( ( nmynods + nghost )*dim * sizeof(double));

  int ierr = gmsh_read_coord_parall( mesh_n, dim, nmynods, mynods, nghost , ghost, *coord );

  return ierr;
}


int get_strain(int e , int gp, int *loc_elem_index, double ***dsh_gp,  double ***bmat, double *strain_gp){

  double  *x_arr; 
  Vec      x_loc; 
  VecGhostGetLocalForm(x, &x_loc);
  VecGetArray(x_loc, &x_arr);

  int  npe = eptr[e+1] - eptr[e];
  for(int i = 0 ; i < npe*dim ; i++)
    elem_disp[i] = x_arr[loc_elem_index[i]];

  VecRestoreArray(x_loc , &x_arr);
  VecGhostRestoreLocalForm(x, &x_loc);

  for(int v = 0; v < nvoi ; v++ ){
    strain_gp[v] = 0.0;
    for(int i = 0 ; i < npe*dim ; i++ )
      strain_gp[v] += bmat[v][i][gp] * elem_disp[i];
    strain_gp[v] = ( fabs(strain_gp[v]) < 1.0e-6 ) ? 0.0 : strain_gp[v];
  }

  return 0;
}


int get_stress(int e, int gp, double *strain_gp, double *stress_gp){

  char        name_s[64];
  material_t  *mat_p;
  get_mat_name(elm_id[e], name_s);

  node_list_t *pn = material_list.head;
  while(pn != NULL){
    mat_p = (material_t *)pn->data;
    if(strcmp(name_s, mat_p->name) == 0) break;
    pn = pn->next;
  }
  if(pn == NULL){
    myio_printf(&MACRO_COMM, "Material %s corresponding to element %d not found on material list\n", name_s, e);
    return 1;
  }

  if(mat_p->type_id == MAT_MICRO){

    message.action = ACTION_MICRO_CALC_STRESS;

    ARRAY_COPY(message.strain_mac, strain_gp, nvoi);

    comm_macro_send(&message);
    comm_macro_recv(&message);

    ARRAY_COPY(stress_gp, message.stress_ave, nvoi);

  }
  else
    material_get_stress(mat_p, dim, strain_gp, stress_gp);

  return 0;
}


int get_c_tan(const char *name, int e, int gp, double *strain_gp, double *c_tan){

  char        name_s[64];
  material_t  *mat_p;
  get_mat_name(elm_id[e], name_s);

  node_list_t *pn = material_list.head;
  while( pn ){
    mat_p = (material_t *)pn->data;
    if(strcmp(name_s, mat_p->name) == 0) break;
    pn = pn->next;
  }
  if(pn == NULL){
    myio_printf(&MACRO_COMM, "Material %s corresponding to element %d not found on material list\n", name_s, e );
    return 1;
  }

  if(mat_p->type_id == MAT_MICRO){

    message.action = ACTION_MICRO_CALC_C_TANGENT;

    ARRAY_COPY(message.strain_mac, strain_gp, nvoi);

    comm_macro_send(&message);
    comm_macro_recv(&message);

    ARRAY_COPY(c_tan, message.c_tangent_ave, nvoi*nvoi);

  }
  else
    material_get_c_tang(mat_p, dim, strain_gp, c_tan);

  return 0;
}


int get_rho(const char *name, int e, double *rho){

  material_t  *mat_p;
  char         name_s[64];
  int          ierr;

  get_mat_name( elm_id[e], name_s );

  node_list_t *pn = material_list.head;
  while( pn ){
    mat_p = ( material_t * )pn->data;
    if( strcmp( name_s, mat_p->name ) == 0 ) break;
    pn = pn->next;
  }
  if( pn == NULL ){
    myio_printf(&MACRO_COMM, "Material %s corresponding to element %d not found on material list\n", name_s, e );
    return 1;
  }

  if(mat_p->type_id == MAT_MICRO){
    ierr = mac_send_signal(WORLD_COMM, RHO); if(ierr) return 1;
    ierr = mac_recv_rho(WORLD_COMM, rho);
  }
  else
    material_get_rho( mat_p, dim, rho );

  return 0;
}


int get_mat_name( int id , char * name_s )
{

  node_list_t *pn;
  physical_t  *phy_p;
  
  pn = physical_list.head;
  while( pn ){
    phy_p = ( physical_t * )pn->data;
    if( id == phy_p->id ) break;
    pn = pn->next;
  }
  if( pn == NULL ) return 1;

  strcpy( name_s, phy_p->name );

  return 0;
}


int get_global_elem_index( int e, int * glo_elem_index )
{

  int  n, d;
  int  npe = eptr[e+1] - eptr[e];

  for( n = 0 ; n < npe ; n++ ){
    for( d = 0 ; d < dim ; d++ )
      glo_elem_index[ n * dim + d ] = loc2petsc[ eind[ eptr[e] + n ] ] * dim + d;
  }
  return 0;
}


int get_local_elem_index( int e, int * loc_elem_index )
{

  int  n, d;
  int  npe = eptr[e+1] - eptr[e];

  for( n = 0 ; n < npe ; n++ ){
    for( d = 0 ; d < dim ; d++ )
      loc_elem_index[ n * dim + d ] = eind[ eptr[e] + n ] * dim + d;
  }
  return 0;
}


int get_dsh( int e, int *loc_elem_index, double ***dsh, double *detj )
{

  double ***dsh_master;
  int       i, gp;
  int       npe = eptr[e+1] - eptr[e];
  int       ngp = npe;

  for( i = 0 ; i < npe*dim ; i++ )
    elem_coor[i] = coord[loc_elem_index[i]];

  for( gp = 0; gp < ngp ; gp++ ){

    fem_get_dsh_master( npe, dim, &dsh_master );

    fem_calc_jac( dim, npe, gp, elem_coor, dsh_master, jac );
    fem_invjac( dim, jac, jac_inv, &detj[gp] );
    fem_trans_dsh( dim, npe, gp, jac_inv, dsh_master, dsh );

  }

  return 0;
}


int get_bmat( int e, double ***dsh, double ***bmat )
{


  int       i, gp;
  int       npe = eptr[e+1] - eptr[e];
  int       ngp = npe;

  if( dim == 2 ){
    for( i = 0 ; i < npe ; i++ ){
      for( gp = 0; gp < ngp ; gp++ ){
	bmat[0][i*dim + 0][gp] = dsh[i][0][gp];
	bmat[0][i*dim + 1][gp] = 0             ;
	bmat[1][i*dim + 0][gp] = 0             ;
	bmat[1][i*dim + 1][gp] = dsh[i][1][gp];
	bmat[2][i*dim + 0][gp] = dsh[i][1][gp];
	bmat[2][i*dim + 1][gp] = dsh[i][0][gp];
      }
    }
  }

  return 0;
}


int get_sh( int dim, int npe, double ***sh )
{

  fem_get_sh( npe, dim, sh );

  return 0;
}


int get_wp( int dim, int npe, double **wp )
{

  fem_get_wp( npe, dim, wp );

  return 0;
}


int get_elem_properties( void )
{

  int      e, v, gp;
  double  *strain_aux = malloc( nvoi * sizeof(double) );
  double  *stress_aux = malloc( nvoi * sizeof(double) );
  double  *wp;

  for ( e = 0 ; e < nelm ; e++ ){

    int     npe = eptr[e+1] - eptr[e];
    int     ngp = npe;
    double  vol_elem = 0.0;

    for ( v = 0 ; v < nvoi ; v++ ) 
      strain_aux[v] = stress_aux[v] = 0.0;

    get_local_elem_index (e, loc_elem_index);

    get_dsh( e, loc_elem_index, dsh, detj);
    get_bmat( e, dsh, bmat );
    get_wp( dim, npe, &wp );

    for ( gp = 0 ; gp < ngp ; gp++ ){

      detj[gp] = fabs( detj[gp] );

      get_strain( e , gp, loc_elem_index, dsh, bmat, strain_gp );
      get_stress( e , gp, strain_gp, stress_gp );
      for ( v = 0 ; v < nvoi ; v++ ){
	strain_aux[v] += strain_gp[v] * detj[gp] * wp[gp];
	stress_aux[v] += stress_gp[v] * detj[gp] * wp[gp];
      }
      vol_elem += detj[gp] * wp[gp];
    }
    for ( v = 0 ; v < nvoi ; v++ ){
      elem_strain[ e*nvoi + v ] = strain_aux[v] / vol_elem;
      elem_stress[ e*nvoi + v ] = stress_aux[v] / vol_elem;
    }

    physical_t * phy;
    node_list_t * pn = physical_list.head;
    while ( pn )
    {
      phy = pn->data;
      if( phy->id == elm_id[e] ) break;
      pn = pn->next;
    }
    if( !pn ) return 1;

    int type = 0;
    pn = material_list.head;
    while ( pn )
    {
      material_t *mat = pn->data;
      if( strcmp( phy->name , mat->name) == 0 ) break;
      pn = pn->next;
      type ++;
    }
    if( !pn ) return 1;

    elem_type[e] = type;
  }

  return 0;
}


int update_boundary( double t , list_t * function_list, list_t * boundary_list )
{

  node_list_t * pn = boundary_list->head;
  while( pn )
  {
    mesh_boundary_t * bou = ( mesh_boundary_t * ) pn->data;
    function_t   * function = NULL;
    int i, d;
    for( d = 0 ; d < dim ; d++ ){
      function_get_from_list( bou->fnum[d] , function_list , &function );
      double val;
      function_eval( t , function , &val );
      for( i = 0 ; i < bou->ndir ; i++ )
	bou->dir_val[ i* (bou->ndirpn) + d ] = val;
    }
    pn = pn->next;
  }

  return 0;
}


int macro_pvtu( char *name )
{

  FILE    *fm;
  char    file_name[NBUF];
  double  *xvalues;
  Vec     xlocal;

  if( rank_mac == 0 ){

    strcpy(file_name,name);
    strcat(file_name,".pvtu");
    fm = fopen(file_name,"w");

    fprintf(fm, "<?xml version=\"1.0\"?>\n"
	"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
	"<PUnstructuredGrid GhostLevel=\"0\">\n"
	"<PPoints>\n"
	"<PDataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\"/>\n"
	"</PPoints>\n"
	"<PCells>\n"
	"<PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>\n"
	"</PCells>\n");

    fprintf(fm, "<PPointData Vectors=\"displ\">\n");
    if( x != NULL )
    fprintf(fm, "<PDataArray type=\"Float64\" Name=\"displ\"    NumberOfComponents=\"3\" />\n");
    if( b != NULL )
      fprintf(fm, "<PDataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" />\n");

    fprintf(fm, "</PPointData>\n"
	"<PCellData>\n"
	"<PDataArray type=\"Int32\"   Name=\"part\"   NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\"/>\n"
	"<PDataArray type=\"Int32\"   Name=\"elem_type\" NumberOfComponents=\"1\"/>\n"
	"</PCellData>\n" , nvoi , nvoi);

    int i;
    for( i = 0 ; i < nproc_mac ; i++ ){
      sprintf(file_name,"%s_%d", name, i );
      fprintf(fm, "<Piece Source=\"%s.vtu\"/>\n", file_name );
    }
    fprintf(fm,	"</PUnstructuredGrid>\n</VTKFile>\n" );

    fclose(fm);

  } // rank = 0

  sprintf( file_name, "%s_%d.vtu", name, rank_mac);
  fm = fopen(file_name,"w"); 
  if(!fm){
    myio_printf(&PETSC_COMM_WORLD,"Problem trying to opening file %s for writing\n", file_name);
    return 1;
  }

  fprintf(fm, 
      "<?xml version=\"1.0\"?>\n"
      "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      "<UnstructuredGrid>\n");
  fprintf(fm,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nallnods, nelm);
  fprintf(fm,"<Points>\n");
  fprintf(fm,"<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");

  int    n , d;

  for( n = 0 ; n < nallnods ; n++ ){
    for( d = 0 ; d < dim ; d++ )
      fprintf( fm,"% 01.6e ",  coord[n*dim + d] );
    for( d = dim ; d < 3 ; d++ )
      fprintf( fm, "% 01.6e ", 0.0 );
    fprintf( fm, "\n" );
  }
  fprintf(fm,"</DataArray>\n");
  fprintf(fm,"</Points>\n");
  fprintf(fm,"<Cells>\n");

  int npe, e;

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for ( e = 0 ; e < nelm ; e++ ){
    npe = eptr[e+1] - eptr[e];
    for ( n = 0 ; n < npe ; n++ )
      fprintf(fm,"%-6d ", eind[eptr[e]+n]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  int ce = 0.0;
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for ( e = 0 ; e < nelm ; e++ ){
    npe = eptr[e+1] - eptr[e];
    ce += npe;
    fprintf(fm,"%d ", ce);
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for ( e = 0 ; e < nelm ; e++ )
    fprintf(fm, "%-3d ", vtkcode( dim , npe ) );
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"</Cells>\n");

  fprintf(fm,"<PointData Vectors=\"displ\">\n"); // Vectors inside is a filter we should not use this here

  if( x != NULL ){
    VecGhostUpdateBegin( x , INSERT_VALUES , SCATTER_FORWARD);
    VecGhostUpdateEnd(   x , INSERT_VALUES , SCATTER_FORWARD);
    VecGhostGetLocalForm(x , &xlocal );

    fprintf(fm,"<DataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
    VecGetArray( xlocal , &xvalues );
    for( n = 0 ; n < nallnods ; n++ ){
      for( d = 0 ; d < dim ; d++ )
	fprintf(fm, "% 01.6e ", xvalues[ n * dim + d ]);
      for( d = dim ; d < 3 ; d++ )
	fprintf(fm,"% 01.6e ",0.0);
      fprintf(fm,"\n");
    }
    VecRestoreArray( xlocal , &xvalues );
    fprintf(fm,"</DataArray>\n");
  }

  if( b != NULL ){
    VecGhostUpdateBegin( b , INSERT_VALUES,SCATTER_FORWARD);
    VecGhostUpdateEnd  ( b , INSERT_VALUES,SCATTER_FORWARD);
    VecGhostGetLocalForm(b , &xlocal);

    fprintf(fm,"<DataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
    VecGetArray(xlocal, &xvalues);
    for( n = 0 ; n < nallnods ; n++ ){
      for( d = 0 ; d < dim ; d++ )
	fprintf(fm, "% 01.6e ", xvalues[ n * dim + d ]);
      for( d = dim ; d < 3 ; d++ )
	fprintf(fm, "% 01.6e ", 0.0);
      fprintf(fm,"\n");
    }
    VecRestoreArray( xlocal , &xvalues );
    fprintf(fm,"</DataArray>\n");

  }
  fprintf(fm,"</PointData>\n");
  fprintf(fm,"<CellData>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"part\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for( e = 0; e < nelm ; e++ )
    fprintf( fm, "%d ", rank_mac );  
  fprintf( fm, "\n");
  fprintf( fm, "</DataArray>\n");

  int v;

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for( e = 0 ; e < nelm ; e++ ){
    for( v = 0 ; v < nvoi ; v++ )
      fprintf( fm, "% 01.6e ", elem_strain[ e*nvoi + v ]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for( e = 0; e < nelm ; e++ ){
    for( v = 0 ; v < nvoi ; v++ )
      fprintf(fm, "% 01.6e ", elem_stress[ e*nvoi + v ]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"elem_type\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for( e = 0; e < nelm ; e++ )
    fprintf( fm, "%d ", elem_type[e] );
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,
      "</CellData>\n"
      "</Piece>\n"
      "</UnstructuredGrid>\n"
      "</VTKFile>\n");

  fclose(fm);
  return 0;
}
