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

#define CHECK_FOUND_GOTO(message) {\
  if (found == false) {\
    myio_printf(MACRO_COMM, "%s\n", message);\
    goto end;}}

#define CHECK_ERROR_GOTO(ierr, message) {\
  if (ierr != 0) {\
    myio_printf(MACRO_COMM, "%s\n", message);\
    goto end;}}

int main(int argc, char **argv) {

  int ierr;
  bool found;

  myio_comm_line_init(argc, argv, &command_line);

  init_variables();

  myio_comm_line_search_option(&command_line, "-coupl", &found);
  if (found == true) flags.coupled = true;

  if (flags.coupled == true) {
    comm.color = COLOR_MACRO;
    ierr = comm_coloring(WORLD_COMM, &comm, &MACRO_COMM);
    if (ierr != 0) {
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
  if (found == true) {
    myio_printf(MACRO_COMM, "%s", help);
    goto end;
  }

  myio_comm_line_search_option(&command_line, "-normal", &found);
  if (found == true) {
    params.calc_mode = CALC_MODE_NORMAL;

    myio_comm_line_get_double(&command_line, "-tf", &params.tf, &found);
    CHECK_FOUND_GOTO("-tf option should be given in -normal mode.\n");

    myio_comm_line_get_double(&command_line, "-dt", &params.dt, &found);
    CHECK_FOUND_GOTO("-dt option should be given in -normal mode.\n");
  }

  myio_comm_line_search_option(&command_line, "-testcomm", &found);
  if (found == true) {
    params.calc_mode = CALC_MODE_TEST;
  }

  myio_comm_line_search_option(&command_line, "-eigen", &found);
  if (found == true) {
    params.calc_mode = CALC_MODE_EIGEN;
    myio_comm_line_get_double(&command_line, "-energy_stored", &params.energy_stored, &found);
  }

  myio_comm_line_get_string(&command_line, "-mesh", mesh_n, &found);
  if (found == false) {
    myio_printf(MACRO_COMM,"mesh file not given on command line.\n");
    goto end;
  }

  FILE *fm = fopen(mesh_n, "r");
  if (fm == NULL) {
    myio_printf(MACRO_COMM,"mesh file not found.\n");
    goto end;
  }

  myio_comm_line_get_int(&command_line, "-dim", &dim, &found);
  if (found == false) {
    myio_printf(MACRO_COMM,"-dim not given on command line.\n");
    goto end;
  }
  mesh.dim = dim;
  gmsh_mesh.dim = dim;

  nvoi = (dim == 2) ? 3 : 6;
  npe_max = (dim == 2) ? 4 : 8;
  ngp_max = npe_max;

  myio_comm_line_search_option(&command_line, "-print_matrices", &found);
  if (found == true) flags.print_matrices = true;

  myio_comm_line_search_option(&command_line, "-print_vectors", &found);
  if (found == true) flags.print_vectors = true;

  myio_comm_line_search_option(&command_line, "-print_pvtu", &found);
  if (found == true) flags.print_pvtu = true;

  myio_comm_line_get_int(&command_line, "-nl_max_its", &params.non_linear_max_its, &found);

  myio_comm_line_get_double(&command_line, "-nl_min_norm_tol", &params.non_linear_min_norm_tol, &found);

  ierr = function_fill_list_from_command_line(&command_line, &function_list);
  CHECK_ERROR_GOTO(ierr, RED "error parsing functions from command line" NORMAL "\n");

  ierr = mesh_fill_boundary_list_from_command_line(&command_line, &boundary_list, &mesh);
  CHECK_ERROR_GOTO(ierr, RED "error parsing boundaries from command line" NORMAL "\n");

  ierr = material_fill_list_from_command_line(&command_line, &material_list);
  CHECK_ERROR_GOTO(ierr, RED "error parsing materials from command line" NORMAL "\n");

  myio_comm_line_search_option(&command_line, "-part_kway", &found);
  if (found == true) mesh.partition = PARMETIS_MESHKWAY;

  myio_comm_line_search_option(&command_line, "-part_geom", &found);
  if (found == true) mesh.partition = PARMETIS_GEOM;

  ierr = gmsh_read_mesh(MACRO_COMM, mesh_n, &gmsh_mesh);
  CHECK_ERROR_GOTO(ierr, RED "error reading gmsh mesh" NORMAL "\n")
  copy_gmsh_to_mesh(&gmsh_mesh, &mesh);

  if (nproc_mac > 1) {
    ierr = mesh_do_partition(MACRO_COMM, &mesh);
    CHECK_ERROR_GOTO(ierr, RED "error partitioning mesh" NORMAL "\n")
  }

  ierr = mesh_calc_local_and_ghost(MACRO_COMM, &mesh);
  CHECK_ERROR_GOTO(ierr, RED "error calculating ghost nodes" NORMAL "\n")

  ierr = mesh_reenumerate(MACRO_COMM, &mesh);
  CHECK_ERROR_GOTO(ierr, RED "error reenumbering nodes" NORMAL "\n")

  ierr = boundary_read();
  CHECK_ERROR_GOTO(ierr, RED "error reading boundaries from mesh" NORMAL "\n")

  list_init(&physical_list, sizeof(physical_t), NULL );
  gmsh_get_physical_list(mesh_n, &physical_list);

  myio_printf(MACRO_COMM, "allocating ");
  ierr = alloc_memory();

  ierr = fem_init();

  if (params.calc_mode == CALC_MODE_EIGEN) {

    EPS eps;
    VecZeroEntries(x);
    VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

    ierr = assembly_AM();
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
    myio_printf(MACRO_COMM, "Number of converged eigenpairs: %d\n", nconv);

    for (int i = 0 ; i < params.num_eigen_vals ; i++) {

      EPSGetEigenpair(eps, i, &params.eigen_vals[i], NULL, x, NULL);
      EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &error);
      myio_printf(MACRO_COMM, "omega %d = %e   error = %e\n", i, params.eigen_vals[i], error);

      if (flags.print_pvtu == true) {
	assembly_get_elem_properties();
	char filename[64];
	sprintf(filename, "macro_eigen_%d", i);
	macro_pvtu(filename);
      }

    }

    EPSDestroy(&eps);

  }else if (params.calc_mode == CALC_MODE_NORMAL) {

    KSP ksp;
    KSPCreate(MACRO_COMM, &ksp);
    KSPSetFromOptions( ksp );

    VecZeroEntries(x);
    VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
    myio_printf(MACRO_COMM, "\n");

    while (params.t < (params.tf + 1.0e-10)) {

      myio_printf(MACRO_COMM,"\ntime step %-3d %-e seg\n", params.ts, params.t);

      boundary_update(params.t);

      boundary_setx();

      params.non_linear_its = 0; params.residual_norm = 2*params.non_linear_min_norm_tol;

      while (params.non_linear_its < params.non_linear_max_its && params.residual_norm > params.non_linear_min_norm_tol) {

	assembly_b(&params.residual_norm);
	myio_printf(MACRO_COMM, GREEN "|b| = %e" NORMAL " ", params.residual_norm);

	if (params.residual_norm < params.non_linear_min_norm_tol)
	  break;

	assembly_A();

	KSPSetOperators(ksp, A, A);
	KSPSolve(ksp, b, dx);
	solvers_print_petsc_ksp_info( MACRO_COMM, ksp);
	myio_printf(MACRO_COMM, "\n");

	VecAXPY(x, 1.0, dx);
	VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
	VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

	params.non_linear_its ++;
      }
      myio_printf(MACRO_COMM, "\n");

      if (flags.print_pvtu == true) {
	assembly_get_elem_properties();
	char filename[64];
	sprintf(filename, "macro_t_%d", params.ts);
	macro_pvtu(filename);
      }

      params.t += params.dt;
      params.ts ++;
    }
    KSPDestroy(&ksp);

  }else if (params.calc_mode == CALC_MODE_TEST) {


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


int copy_gmsh_to_mesh(gmsh_mesh_t *gmsh_mesh, mesh_t *mesh) {

  mesh->nelm_local = gmsh_mesh->nelm_local;
  mesh->nelm_total = gmsh_mesh->nelm_total;

  mesh->eptr = malloc((mesh->nelm_local+1)*sizeof(int));
  ARRAY_COPY(mesh->eptr, gmsh_mesh->eptr, mesh->nelm_local + 1);

  mesh->elm_id = malloc(mesh->nelm_local*sizeof(int));
  ARRAY_COPY(mesh->elm_id, gmsh_mesh->elm_id, mesh->nelm_local);

  mesh->npe = malloc(mesh->nelm_local*sizeof(int));
  for (int i = 0 ; i < mesh->nelm_local ; i++)
    mesh->npe[i] = mesh->eptr[i+1] - mesh->eptr[i];

  mesh->eind = malloc(mesh->eptr[mesh->nelm_local]*sizeof(int));
  for (int i = 0 ; i < mesh->eptr[mesh->nelm_local] ; i++)
    mesh->eind[i] = gmsh_mesh->eind[i] - 1;

  mesh->nelm_dist = malloc(nproc_mac*sizeof(int));
  ARRAY_COPY(mesh->nelm_dist, gmsh_mesh->nelm_dist, nproc_mac);

  mesh->dim = gmsh_mesh->dim;
  mesh->nnods_total = gmsh_mesh->nnods;

  mesh->coord = malloc(mesh->nnods_total*dim*sizeof(double));
  ARRAY_COPY(mesh->coord, gmsh_mesh->coord, mesh->nnods_total*dim)

  return 0;
}
