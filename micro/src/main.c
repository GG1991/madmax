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

#define CHECK_AND_GOTO(error){if(error){myio_printf(MICRO_COMM, "error line %d at %s\n", __LINE__, __FILE__); goto end;}}
#define CHECK_INST_ELSE_GOTO(cond, instr){if(cond){instr}else{myio_printf(MICRO_COMM, "error line %d at %s\n", __LINE__, __FILE__); goto end;}}
#define CHECK_ERROR_GOTO(message){if(ierr != 0){myio_printf(MICRO_COMM, "%s\n", message); goto end;}}

int main(int argc, char **argv){

  int ierr;
  int values_i[10];
  char string_buf[NBUF];
  int nval_expect, nval_found;
  bool found;

  myio_comm_line_init(argc, argv, &command_line);

  init_variables();

  MPI_Comm_size(WORLD_COMM, &nproc_wor);
  MPI_Comm_rank(WORLD_COMM, &rank_wor);

  myio_comm_line_search_option(&command_line, "-help", &found);
  if(found == true){
    myio_printf(MACRO_COMM, "%s", help);
    goto end;
  }

  myio_comm_line_search_option(&command_line, "-coupl", &found);
  if(found == true) flags.coupled = true;

  macmic.type = COUP_1;
  color = COLOR_MICRO; /* color can change */
  ierr = macmic_coloring(WORLD_COMM, &color, &macmic, &MICRO_COMM, flags.coupled); CHECK_AND_GOTO(ierr);

  MPI_Comm_size(MICRO_COMM, &nproc_mic);
  MPI_Comm_rank(MICRO_COMM, &rank_mic);
  
  nx = ny = nz = -1;
  hx = hy = hz = -1;
  lx = ly = lz = -1;

  ierr = myio_comm_line_get_int(&command_line, "-dim", &dim, &found); CHECK_AND_GOTO(ierr)

  nvoi = (dim == 2) ? 3 : 6;

  ierr = myio_comm_line_get_string(&command_line, "-micro_struct", string_buf, &found); CHECK_AND_GOTO(ierr)
  micro_struct_init(dim, string_buf, &micro_struct);
  
  lx = micro_struct.size[0];
  ly = micro_struct.size[1];
  lz = (dim == 3) ? micro_struct.size[2] : -1;

  if(dim == 2) nval_expect = 2;
  if(dim == 3) nval_expect = 3;

  myio_comm_line_get_int_array(&command_line, "-struct_n", values_i, nval_expect, &nval_found, &found);
  if(found){

    if(nval_found != nval_expect){
      myio_printf(MICRO_COMM,"-struct_n should include %d arguments\n", nval_expect);
      goto end;
    }
    nx   = values_i[0];
    ny   = values_i[1];
    nz   = (dim == 3) ? values_i[2] : 1;

    nn   = nx*ny*nz;
    nex  = (nx-1);

    ney  = (ny-1)/nproc_mic + (((ny-1) % nproc_mic > rank_mic) ? 1:0); 
    nez  = (nz-1);
    nelm = (dim == 2) ? nex*ney : nex*ney*nez;
    nyl  = (rank_mic == 0) ? ney+1 : ney;
    nl   = (dim == 2) ? nyl*nx : nyl*nx*nz;

    hx   = lx/nex;
    hy   = ly/(ny-1);
    hz   = (dim == 3) ? (lz/nez) : -1;

    int *nyl_arr = malloc(nproc_mic * sizeof(int));
    ierr = MPI_Allgather(&nyl, 1, MPI_INT, nyl_arr, 1, MPI_INT, MICRO_COMM); if(ierr) return 1;
    ny_inf = 0;
    for(int i = 0 ; i < rank_mic ; i++)
      ny_inf += nyl_arr[i];
    free(nyl_arr);

    npe  = (dim == 2) ? 4 : 8;
    ngp  = (dim == 2) ? 4 : 8;
    if(ny < nproc_mic){
      myio_printf(MICRO_COMM, "ny %d not large enough to be executed with %d processes\n", ny, nproc_mic);
      goto end;
    }

  }
  else{
    myio_printf(MICRO_COMM,"-struct_n is request\n");
    goto end;
  }

  center_coor = malloc(dim*sizeof(double));
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

  myio_comm_line_get_int(&command_line, "-nl_max_its", &params.non_linear_max_its, &found);
  myio_comm_line_get_int(&command_line, "-nl_min_norm_tol", &params.non_linear_max_its, &found);

  myio_comm_line_search_option(&command_line, "-print_matrices", &found);
  if(found == true)
    flags.print_matrices = true;

  myio_comm_line_search_option(&command_line, "-print_vectors", &found);
  if(found == true)
    flags.print_vectors = true;

  myio_comm_line_search_option(&command_line, "-print_pvtu", &found);
  if(found == true &&
      params.homog_method != HOMOG_METHOD_TAYLOR_PARALLEL &&
      params.homog_method != HOMOG_METHOD_TAYLOR_SERIAL)
    flags.print_pvtu = true;

  PRINTF1(
      "--------------------------------------------------\n"
      "  MICRO: START\n"
      "--------------------------------------------------\n\n");

  PRINTF1("allocating...\n")
  ierr = alloc_memory();

  ierr = micro_struct_init_elem_type(&micro_struct, dim, nelm, &get_elem_centroid, elem_type); CHECK_AND_GOTO(ierr)

  ierr = micro_check_material_and_elem_type(&material_list, elem_type, nelm);
  CHECK_ERROR_GOTO("error checking elem_type and material_list");

  init_shapes(&struct_sh, &struct_dsh, &struct_wp);

  for(int gp = 0; gp < ngp ; gp++){
    for(int is = 0; is < npe ; is++){
      if(dim == 2){
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

  if(flags.coupled == true){

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

	case ACTION_MICRO_CALC_RHO:

	  message.rho = params.rho;
	  break;

	case ACTION_MICRO_END:
	  break;

	default:
	  myio_printf(MICRO_COMM, "MICRO:signal %d not identified\n", signal);
	  goto end;

      }

      if(message.action != ACTION_MICRO_END)
	ierr = comm_micro_send(&message);
    }
  }
  else{

    myio_printf(MICRO_COMM,"\nConstitutive Average Tensor\n");
    for(int i = 0 ; i < nvoi ; i++){
      for(int j = 0 ; j < nvoi ; j++)
	myio_printf(MICRO_COMM, "%e ", (fabs(params.c_tangent_linear[i*nvoi+j])>1.0) ? params.c_tangent_linear[i*nvoi+j] : 0.0);
      myio_printf(MICRO_COMM, "\n");
    }
    myio_printf(MICRO_COMM, "\n");

  }

  micro_print_info();

  end_trace(MICRO_COMM);

end:

  PRINTF1(
      "--------------------------------------------------\n"
      "  MICRO: FINISH COMPLETE\n"
      "--------------------------------------------------\n");

  ierr = finalize();

  return ierr;
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
