#include "macro.h"

static char help[] =
"macro multiscale code \n"
"-coupl       : coupled with \"micro\" code for solving multiscale problem \n"
"-normal      : normal execution, solves a time dependent boundary condition problem \n"
"-testcomm    : communication testing with the \"micro\" code \n"
"-eigen       : calculates the eigensystem Mx = -(1/omega)Kx \n"
"-print_matrices \n"
"-print_vectors \n"
"-print_pvtu \n";

params_t params;
flags_t flags;
solver_t solver;
gmsh_mesh_t gmsh_mesh;
mesh_t mesh;
comm_t comm;
list_t physical_list;
list_t boundary_list;

#define CHECK_FOUND_GOTO(message){\
  if(found == false){\
    myio_printf(MACRO_COMM, "%s\n", message);\
    goto end;}}

#define CHECK_ERROR_GOTO(ierr, message){\
  if(ierr != 0){\
    myio_printf(MACRO_COMM, "%s\n", message);\
    goto end;}}

int main(int argc, char **argv){

  int ierr;
  bool found;

  myio_comm_line_init(argc, argv, &command_line);

  init_variables();

  myio_comm_line_search_option(&command_line, "-coupl", &found);
  if(found == true) flags.coupled = true;

  if(flags.coupled == true){
    comm.color = COLOR_MACRO;
    ierr = comm_coloring(WORLD_COMM, &comm, &MACRO_COMM);
    if(ierr != 0){
      flags.coupled = false;
      myio_printf(MACRO_COMM, RED "error in coloring" NORMAL "\n");
      goto end_no_message;
    }
  }

  MPI_Comm_size(MACRO_COMM, &nproc_mac);
  MPI_Comm_rank(MACRO_COMM, &rank_mac);

  myio_printf(MACRO_COMM, GREEN
      "--------------------------------------------------\n"
      "  MACRO: START\n"
      "--------------------------------------------------" NORMAL "\n\n");

  myio_comm_line_search_option(&command_line, "-help", &found);
  if(found == true){
    myio_printf(MACRO_COMM, "%s", help);
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

  myio_comm_line_get_string(&command_line, "-mesh", mesh_n, &found);
  if(found == false){
    myio_printf(MACRO_COMM,"mesh file not given on command line.\n");
    goto end;
  }

  FILE *fm = fopen(mesh_n, "r");
  if(fm == NULL){
    myio_printf(MACRO_COMM,"mesh file not found.\n");
    goto end;
  }

  myio_comm_line_get_int(&command_line, "-dim", &dim, &found);
  if(found == false){
    myio_printf(MACRO_COMM,"-dim not given on command line.\n");
    goto end;
  }
  mesh.dim = dim;
  gmsh_mesh.dim = dim;

  nvoi = (dim == 2) ? 3 : 6;
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

  ierr = function_fill_list_from_command_line(&command_line, &function_list);
  CHECK_ERROR_GOTO(ierr, RED "error parsing functions from command line" NORMAL "\n");

  ierr = mesh_fill_boundary_list_from_command_line(&command_line, &boundary_list, &mesh);
  CHECK_ERROR_GOTO(ierr, RED "error parsing boundaries from command line" NORMAL "\n");

  ierr = material_fill_list_from_command_line(&command_line, &material_list);
  CHECK_ERROR_GOTO(ierr, RED "error parsing materials from command line" NORMAL "\n");

  mesh.partition = PARMETIS_GEOM;
  myio_comm_line_search_option(&command_line, "-part_kway", &found);
  if(found == true) mesh.partition = PARMETIS_MESHKWAY;

  myio_comm_line_search_option(&command_line, "-part_geom", &found);
  if(found == true) mesh.partition = PARMETIS_GEOM;

  ierr = gmsh_read_mesh(MACRO_COMM, mesh_n, &gmsh_mesh);
  CHECK_ERROR_GOTO(ierr, RED "error reading gmsh mesh" NORMAL "\n")
  copy_gmsh_to_mesh(&gmsh_mesh, &mesh);

  if(nproc_mac > 1){
    ierr = mesh_do_partition(MACRO_COMM, &mesh);
    CHECK_ERROR_GOTO(ierr, RED "error partitioning mesh" NORMAL "\n")
  }

  ierr = mesh_calc_local_and_ghost(MACRO_COMM, &mesh);
  CHECK_ERROR_GOTO(ierr, RED "error calculating ghost nodes" NORMAL "\n")

  ierr = mesh_reenumerate(MACRO_COMM, &mesh);
  CHECK_ERROR_GOTO(ierr, RED "error reenumbering nodes" NORMAL "\n")

  ierr = read_bc();
  CHECK_ERROR_GOTO(ierr, RED "error reading boundaries from mesh" NORMAL "\n")

  list_init(&physical_list, sizeof(physical_t), NULL );
  gmsh_get_physical_list(mesh_n, &physical_list);

  myio_printf(MACRO_COMM, "allocating ");
  ierr = alloc_memory();

  ierr = fem_init();

  if(params.calc_mode == CALC_MODE_EIGEN){

    EPS eps;
    VecZeroEntries(x);
    VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

    ierr = assembly_AM_petsc();
    CHECK_ERROR_GOTO(ierr, "problem during matrix assembly\n")

    int nconv;
    double error;

    EPSCreate(MACRO_COMM, &eps);
    EPSSetOperators(eps, M, A);
    EPSSetProblemType(eps, EPS_GHEP);
    EPSSetFromOptions(eps);
    EPSGetDimensions(eps, &params.num_eigen_vals, NULL, NULL);
    params.eigen_vals = malloc( params.num_eigen_vals*sizeof(double));
    myio_printf(MACRO_COMM,"Number of requested eigenvalues: %d\n", params.num_eigen_vals);

    EPSSolve(eps);
    EPSGetConverged(eps, &nconv);
    myio_printf(MACRO_COMM,"Number of converged eigenpairs: %d\n",nconv);

    for(int i = 0 ; i < params.num_eigen_vals ; i++){

      EPSGetEigenpair( eps, i, &params.eigen_vals[i], NULL, x, NULL );
      EPSComputeError( eps, i, EPS_ERROR_RELATIVE, &error );
      myio_printf(MACRO_COMM, "omega %d = %e   error = %e\n", i, params.eigen_vals[i], error);

      if(flags.print_pvtu == true){
	get_elem_properties();
	char filename[64];
	sprintf(filename, "macro_eigen_%d", i);
	macro_pvtu(filename);
      }

    }

    EPSDestroy(&eps);

  }
  else if(params.calc_mode == CALC_MODE_NORMAL){

    KSP ksp;
    KSPCreate(MACRO_COMM, &ksp);
    KSPSetFromOptions( ksp );

    VecZeroEntries(x);
    VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
    myio_printf(MACRO_COMM, "\n");

    while(params.t < (params.tf + 1.0e-10)){

      myio_printf(MACRO_COMM,"\ntime step %-3d %-e seg\n", params.ts, params.t);

      update_boundary(params.t, &function_list, &boundary_list);

      Vec x_loc;
      double *x_arr;
      VecGhostGetLocalForm(x, &x_loc);
      VecGetArray(x_loc, &x_arr);

      node_list_t * pn = boundary_list.head;
      while(pn != NULL){
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

	assembly_b_petsc();
	VecNorm(b, NORM_2, &params.residual_norm);
	myio_printf(MACRO_COMM, GREEN "|b| = %e" NORMAL " ", params.residual_norm);

	if(params.residual_norm < params.non_linear_min_norm_tol)
	  break;

	assembly_A_petsc();

	KSPSetOperators(ksp, A, A);
	KSPSolve(ksp, b, dx);
	print_petsc_ksp_info( MACRO_COMM, ksp);
	myio_printf(MACRO_COMM, "\n");

	VecAXPY(x, 1.0, dx);
	VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
	VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

	params.non_linear_its ++;
      }
      myio_printf(MACRO_COMM, "\n");

      if(flags.print_pvtu == true){
	get_elem_properties();
	char filename[64];
	sprintf(filename, "macro_t_%d", params.ts);
	macro_pvtu(filename);
      }

      params.t += params.dt;
      params.ts ++;
    }
    KSPDestroy(&ksp);

  }else if(params.calc_mode == CALC_MODE_TEST){


  }
  myio_printf(MACRO_COMM, "\n");

end:

  myio_printf(MACRO_COMM, GREEN
      "--------------------------------------------------\n"
      "  MACRO: FINISH\n"
      "--------------------------------------------------" NORMAL "\n");

end_no_message:

  finalize();

  return 0;
}


int copy_gmsh_to_mesh(gmsh_mesh_t *gmsh_mesh, mesh_t *mesh){

  mesh->nelm_local = gmsh_mesh->nelm_local;
  mesh->nelm_total = gmsh_mesh->nelm_total;

  mesh->eptr = malloc((mesh->nelm_local+1)*sizeof(int));
  ARRAY_COPY(mesh->eptr, gmsh_mesh->eptr, mesh->nelm_local + 1);

  mesh->elm_id = malloc(mesh->nelm_local*sizeof(int));
  ARRAY_COPY(mesh->elm_id, gmsh_mesh->elm_id, mesh->nelm_local);

  mesh->npe = malloc(mesh->nelm_local*sizeof(int));
  for(int i = 0 ; i < mesh->nelm_local ; i++)
    mesh->npe[i] = mesh->eptr[i+1] - mesh->eptr[i];

  mesh->eind = malloc(mesh->eptr[mesh->nelm_local]*sizeof(int));
  for(int i = 0 ; i < mesh->eptr[mesh->nelm_local] ; i++)
    mesh->eind[i] = gmsh_mesh->eind[i] - 1;

  mesh->nelm_dist = malloc(nproc_mac*sizeof(int));
  ARRAY_COPY(mesh->nelm_dist, gmsh_mesh->nelm_dist, nproc_mac);

  mesh->dim = gmsh_mesh->dim;
  mesh->nnods_total = gmsh_mesh->nnods;

  mesh->coord = malloc(mesh->nnods_total*dim*sizeof(double));
  ARRAY_COPY(mesh->coord, gmsh_mesh->coord, mesh->nnods_total*dim)

  return 0;
}


int read_bc(void){

  int ierr;
  mesh_boundary_t *bou;
  node_list_t *pn = boundary_list.head;

  while(pn != NULL){

    int *ix, n;
    bou = (mesh_boundary_t *)pn->data;
    ierr = gmsh_get_node_index(mesh_n, bou->name, mesh.nnods_local, mesh.local_nods, dim, &n, &ix);

    bou->ndir = n;
    bou->ndirix = bou->ndir * bou->ndirpn;
    bou->dir_val = malloc(bou->ndirix * sizeof(double));
    bou->dir_loc_ixs = malloc(bou->ndirix * sizeof(int));
    bou->dir_glo_ixs = malloc(bou->ndirix * sizeof(int));

    for(int i = 0 ; i < n ; i++){

      int da = 0;
      int *p = bsearch(&ix[i], mesh.local_nods, mesh.nnods_local, sizeof(int), mesh_cmpfunc);

      for(int d = 0 ; d < dim ; d++){
	if(bou->kind & (1<<d)){
	  bou->dir_loc_ixs[i*(bou->ndirpn) + da] = (p - mesh.local_nods) * dim + d;
	  bou->dir_glo_ixs[i*(bou->ndirpn) + da] = mesh.local_to_global[(p - mesh.local_nods)] * dim + d;
	  da++;
	}
      }
    }

    free(ix);
    pn = pn->next;
  }

  return ierr;
}


int get_strain(int e , int gp, int *loc_elem_index, double ***dsh_gp,  double ***bmat, double *strain_gp){

  double  *x_arr;
  Vec      x_loc;
  VecGhostGetLocalForm(x, &x_loc);
  VecGetArray(x_loc, &x_arr);

  int  npe = mesh.eptr[e+1] - mesh.eptr[e];
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
  get_mat_name(mesh.elm_id[e], name_s);

  node_list_t *pn = material_list.head;
  while(pn != NULL){
    mat_p = (material_t *)pn->data;
    if(strcmp(name_s, mat_p->name) == 0) break;
    pn = pn->next;
  }
  if(pn == NULL){
    myio_printf(MACRO_COMM, "Material %s corresponding to element %d not found on material list\n", name_s, e);
    return 1;
  }

  if(mat_p->type_id == MAT_MICRO){

    message.action = ACTION_MICRO_CALC_STRESS;

    ARRAY_COPY(message.strain_mac, strain_gp, nvoi);

    comm_macro_send(&message, &comm);
    comm_macro_recv(&message, &comm);

    ARRAY_COPY(stress_gp, message.stress_ave, nvoi);

  }
  else
    material_get_stress(mat_p, dim, strain_gp, stress_gp);

  return 0;
}


int get_c_tan(const char *name, int e, int gp, double *strain_gp, double *c_tan){

  char name_s[64];
  material_t *mat_p;
  get_mat_name(mesh.elm_id[e], name_s);

  node_list_t *pn = material_list.head;
  while(pn != NULL){
    mat_p = (material_t *)pn->data;
    if(strcmp(name_s, mat_p->name) == 0) break;
    pn = pn->next;
  }
  if(pn == NULL){
    myio_printf(MACRO_COMM, "Material %s corresponding to element %d not found on material list\n", name_s, e );
    return 1;
  }

  if(mat_p->type_id == MAT_MICRO){

    message.action = ACTION_MICRO_CALC_C_TANGENT;

    ARRAY_COPY(message.strain_mac, strain_gp, nvoi);

    comm_macro_send(&message, &comm);
    comm_macro_recv(&message, &comm);

    ARRAY_COPY(c_tan, message.c_tangent_ave, nvoi*nvoi);

  }
  else
    material_get_c_tang(mat_p, dim, strain_gp, c_tan);

  return 0;
}


int get_rho(const char *name, int e, double *rho){

  char name_s[64];
  int ierr;

  get_mat_name(mesh.elm_id[e], name_s);

  material_t *mat_p;
  node_list_t *pn = material_list.head;
  while(pn != NULL){
    mat_p = (material_t *)pn->data;
    if(strcmp(name_s, mat_p->name) == 0) break;
    pn = pn->next;
  }
  if(pn == NULL){
    myio_printf(MACRO_COMM, "Material %s corresponding to element %d not found on material list\n", name_s, e );
    return 1;
  }

  if(mat_p->type_id == MAT_MICRO){

    message.action = ACTION_MICRO_CALC_RHO;

    ierr = comm_macro_send(&message, &comm);
    ierr = comm_macro_recv(&message, &comm);

    *rho = message.rho;
  }
  else
    material_get_rho( mat_p, dim, rho );

  return ierr;
}


int get_mat_name(int id, char *name_s){

  node_list_t *pn;
  physical_t  *phy_p;

  pn = physical_list.head;
  while(pn != NULL){
    phy_p = ( physical_t * )pn->data;
    if( id == phy_p->id ) break;
    pn = pn->next;
  }
  if(pn == NULL) return 1;

  strcpy(name_s, phy_p->name);

  return 0;
}


int get_global_elem_index(int e, int *glo_elem_index){

  int  npe = mesh.eptr[e+1] - mesh.eptr[e];

  for(int n = 0 ; n < npe ; n++){
    for(int d = 0 ; d < dim ; d++)
      glo_elem_index[n*dim + d] = mesh.local_to_global[mesh.eind[mesh.eptr[e] + n]]*dim + d;
  }
  return 0;
}


int get_local_elem_index(int e, int *loc_elem_index){

  int  npe = mesh.eptr[e+1] - mesh.eptr[e];

  for(int n = 0 ; n < npe ; n++){
    for(int d = 0 ; d < dim ; d++)
      loc_elem_index[n*dim + d] = mesh.eind[mesh.eptr[e] + n]*dim + d;
  }
  return 0;
}


int get_dsh(int e, int *loc_elem_index, double ***dsh, double *detj){

  double ***dsh_master;
  int npe = mesh.eptr[e+1] - mesh.eptr[e];
  int ngp = npe;

  for(int i = 0 ; i < npe*dim ; i++)
    elem_coor[i] = mesh.coord_local[loc_elem_index[i]];

  fem_get_dsh_master(npe, dim, &dsh_master);

  for(int gp = 0; gp < ngp ; gp++){
    fem_calc_jac(dim, npe, gp, elem_coor, dsh_master, jac);
    fem_invjac(dim, jac, jac_inv, &detj[gp]);
    fem_trans_dsh(dim, npe, gp, jac_inv, dsh_master, dsh);
  }

  return 0;
}


int get_bmat(int e, double ***dsh, double ***bmat){

  int npe = mesh.eptr[e+1] - mesh.eptr[e];
  int ngp = npe;

  if(dim == 2){
    for(int i = 0 ; i < npe ; i++){
      for(int gp = 0; gp < ngp ; gp++){
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


int get_sh(int dim, int npe, double ***sh){

  return fem_get_sh(npe, dim, sh);
}


int get_wp(int dim, int npe, double **wp){

  return fem_get_wp(npe, dim, wp);
}


int get_elem_properties(void){

  double *strain_aux = malloc(nvoi*sizeof(double));
  double *stress_aux = malloc(nvoi*sizeof(double));
  double *wp;

  for(int e = 0 ; e < mesh.nelm_local ; e++){

    int npe = mesh.eptr[e+1] - mesh.eptr[e];
    int ngp = npe;
    double vol_elem = 0.0;

    for(int v = 0 ; v < nvoi ; v++)
      strain_aux[v] = stress_aux[v] = 0.0;

    get_local_elem_index(e, loc_elem_index);

    get_dsh(e, loc_elem_index, dsh, detj);
    get_bmat(e, dsh, bmat);
    get_wp(dim, npe, &wp);

    for(int gp = 0 ; gp < ngp ; gp++){

      detj[gp] = fabs(detj[gp]);

      get_strain(e, gp, loc_elem_index, dsh, bmat, strain_gp);
      get_stress(e, gp, strain_gp, stress_gp);
      for(int v = 0 ; v < nvoi ; v++){
	strain_aux[v] += strain_gp[v] * detj[gp] * wp[gp];
	stress_aux[v] += stress_gp[v] * detj[gp] * wp[gp];
      }
      vol_elem += detj[gp] * wp[gp];
    }
    for(int v = 0 ; v < nvoi ; v++){
      elem_strain[e*nvoi + v] = strain_aux[v] / vol_elem;
      elem_stress[e*nvoi + v] = stress_aux[v] / vol_elem;
    }

    physical_t * phy;
    node_list_t * pn = physical_list.head;
    while(pn != NULL){
      phy = pn->data;
      if(phy->id == mesh.elm_id[e]) break;
      pn = pn->next;
    }
    if(pn == NULL) return 1;

    int type = 0;
    pn = material_list.head;
    while(pn != NULL){
      material_t *mat = pn->data;
      if(strcmp(phy->name, mat->name) == 0) break;
      pn = pn->next;
      type ++;
    }
    if(pn == NULL) return 1;

    elem_type[e] = type;
  }

  return 0;
}


int update_boundary(double t, list_t *function_list, list_t *boundary_list){

  node_list_t * pn = boundary_list->head;
  while(pn != NULL){

    mesh_boundary_t * bou = ( mesh_boundary_t * ) pn->data;
    function_t   * function = NULL;
    for(int d = 0 ; d < dim ; d++){
      function_get_from_list(bou->fnum[d], function_list, &function);
      double val;
      function_eval(t, function, &val);
      for(int i = 0 ; i < bou->ndir ; i++)
	bou->dir_val[i*(bou->ndirpn) + d] = val;
    }
    pn = pn->next;
  }

  return 0;
}
