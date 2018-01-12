#include "micro.h"

static char help[] =
"micro multiscale code \n"
"-coupl    [0 (no coupling ) | 1 (coupling with micro)] \n"
"-homo_ts     : c =  vi ci + vm cm            (serial) \n"
"-homo_tp     : c = (vi ci^-1 + vm cm^-1)^-1  (parallel) \n"
"-homo_us     : homogenization using uniform strains approach \n"
"-struct_n [<nx,ny>] if dim = 2 \n"
"-print_matrices \n"
"-print_vectors \n"
"-print_pvtu \n";

params_t params;
flags_t flags;
solver_t solver;
comm_t comm;
mesh_struct_t mesh_struct;

#define CHECK_AND_GOTO(error) {\
if (error) {\
  myio_printf(MICRO_COMM, "error line %d at %s\n", __LINE__, __FILE__);\
  goto end;}}

#define CHECK_ERROR_GOTO(message) {\
  if (ierr != 0) {\
    myio_printf(MICRO_COMM, "%s\n", message);\
    goto end;}}

int main(int argc, char **argv)
{
  int ierr;
  bool found;

  myio_comm_line_init(argc, argv, &command_line);

  init_variables();

  myio_comm_line_search_option(&command_line, "-help", &found);
  if (found == true) {
    myio_printf(MACRO_COMM, "%s", help);
    goto end;
  }

  ierr = comm_line_set_flags();

  comm.color = COLOR_MICRO; /* color changes */
  if (flags.coupled == true) {
    ierr = comm_coloring(WORLD_COMM, &comm, &MICRO_COMM); CHECK_AND_GOTO(ierr);
  }

  MPI_Comm_size(MICRO_COMM, &nproc_mic);
  MPI_Comm_rank(MICRO_COMM, &rank_mic);

  if(nproc_mic != 1){
    myio_printf(MICRO_COMM, RED "error during coloring, more than one process detected in micro\n" NORMAL);
    goto end_no_message;
  }

  ierr = myio_comm_line_get_int(&command_line, "-dim", &dim, &found);
  if (found == false) {
    myio_printf(MICRO_COMM,"-dim not given on command line.\n");
    goto end;
  }

  nvoi = (dim == 2) ? 3 : 6;

  char string_buf[NBUF];
  ierr = myio_comm_line_get_string(&command_line, "-micro_struct", string_buf, &found);
  ierr = micro_struct_init(dim, string_buf, &micro_struct);

  int nval_found, nval_expect = (dim == 2) ? 2 : 3;
  int values_int[3];
  myio_comm_line_get_int_array(&command_line, "-struct_n", values_int, nval_expect, &nval_found, &found);

  if (found == true) {
    if (nval_found != nval_expect) {
      myio_printf(MICRO_COMM, RED "-struct_n should include %d arguments\n" NORMAL, nval_expect);
      return 1;
    }
  }else{
    myio_printf(MICRO_COMM, RED "-struct_n is request\n" NORMAL, nval_expect);
    return 1;
  }

  mesh_struct_init(dim, values_int, micro_struct.size, &mesh_struct);

  ngp = (dim == 2) ? 4 : 8;

  ierr = material_fill_list_from_command_line(&command_line, &material_list); CHECK_AND_GOTO(ierr)

  PRINTF1(GREEN
      "--------------------------------------------------\n"
      "  MICRO: START\n"
      "--------------------------------------------------" NORMAL "\n\n");

  ierr = alloc_memory();

  ierr = micro_struct_init_elem_type(&micro_struct, dim, mesh_struct.nelm, &get_elem_centroid, elem_type); CHECK_AND_GOTO(ierr)

  ierr = micro_check_material_and_elem_type(&material_list, elem_type, mesh_struct.nelm);
  CHECK_ERROR_GOTO("error checking elem_type and material_list");

  double h[3]; h[0] = mesh_struct.hx; h[1] = mesh_struct.hy; h[2] = mesh_struct.hz;
  fem_init();
  fem_init_struct(&struct_sh, &struct_dsh, &struct_wp, h, mesh_struct.dim);

  for (int gp = 0; gp < ngp ; gp++) {
    for (int is = 0; is < mesh_struct.npe ; is++) {
      if (dim == 2) {
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

  homogenize_init();

  double strain_mac[6], strain_ave[6], stress_ave[6];

  double *c_tangent_ave = malloc(36*sizeof(double));

  if (flags.coupled == true) {

    while (message.action != ACTION_MICRO_END) {

      ierr = comm_micro_recv(&message, &comm);

      switch(message.action) {

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

	case ACTION_MICRO_CALC_RHO:

	  message.rho = params.rho;
	  break;

	case ACTION_MICRO_END:
	  break;

	default:
	  myio_printf(MICRO_COMM, "MICRO:signal %d not identified\n", signal);
	  goto end;

      }

      if (message.action != ACTION_MICRO_END)
	ierr = comm_micro_send(&message, &comm);
    }
  }
  else{

    myio_printf(MICRO_COMM,"\nConstitutive Average Tensor\n");
    for (int i = 0 ; i < nvoi ; i++) {
      for (int j = 0 ; j < nvoi ; j++)
	myio_printf(MICRO_COMM, "%e ", (fabs(params.c_tangent_linear[i*nvoi+j])>1.0) ? params.c_tangent_linear[i*nvoi+j] : 0.0);
      myio_printf(MICRO_COMM, "\n");
    }
    myio_printf(MICRO_COMM, "\n");

  }

  micro_print_info();

  end_trace(MICRO_COMM);

end:

  PRINTF1(GREEN
      "--------------------------------------------------\n"
      "  MICRO: FINISH\n"
      "--------------------------------------------------" NORMAL "\n");

end_no_message:

  ierr = finalize();

  return ierr;
}


int micro_check_material_and_elem_type(list_t *material_list, int *elem_type, int nelm) {

  char *word_to_search;

  for (int e ; e < nelm ; e++ ) {

    switch(elem_type[e]) {

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

    while (pm != NULL) {
      mat_p = (material_t *)pm->data;
      if (strcmp(mat_p->name, word_to_search) == 0) break;
      pm = pm->next;
    }

    if (pm == NULL) return 1;
  }

  return 0;
}
