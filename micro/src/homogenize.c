#include "micro.h"

int homog_get_strain_stress(double *strain_mac, double *strain_ave, double *stress_ave)
{
  int ierr;

  if (flags.linear_materials == true && flags.c_linear_calculated == true) {

    for (int i = 0 ; i < nvoi ; i++) {
      strain_ave[i] = strain_mac[i];
      stress_ave[i] = 0.0;
      for (int j = 0 ; j < nvoi ; j++)
	stress_ave[i] += params.c_tangent_linear[i*nvoi + j] * strain_mac[j];
    }
  }else{

    ierr = homog_get_strain_stress_non_linear(strain_mac, strain_ave, stress_ave);
    if (ierr) return 1;
  }
  return 0;
}

int homog_get_strain_stress_non_linear(double *strain_mac, double *strain_ave, double *stress_ave)
{
  int ierr = 0;

  if (params.multis_method == MULTIS_MIXP || params.multis_method == MULTIS_MIXS)
    ierr = homog_taylor(strain_mac, strain_ave, stress_ave);
  else if (params.multis_method == MULTIS_FE2)
    ierr = homog_fe2(strain_mac, strain_ave, stress_ave);

  return ierr;
}

int homog_get_c_tangent(double *strain_mac, double **c_tangent)
{
  int ierr = 0;

  if (flags.linear_materials == true && flags.linear_materials == true)
    (*c_tangent) = params.c_tangent_linear;

  else if (flags.linear_materials == false) {

    ierr = homog_calculate_c_tangent(strain_mac, params.c_tangent);
    (*c_tangent) = params.c_tangent;
  }

  return ierr;
}

int homog_calculate_c_tangent_around_zero(double *c_tangent)
{
  double strain_zero[MAX_NVOIGT];
  ARRAY_SET_TO_ZERO(strain_zero, nvoi);

  return homog_calculate_c_tangent(strain_zero, c_tangent);
}

int homog_calculate_c_tangent(double *strain_mac, double *c_tangent)
{
  double strain_1[MAX_NVOIGT], strain_2[MAX_NVOIGT];
  double stress_1[MAX_NVOIGT], stress_2[MAX_NVOIGT];
  double strain_aux[MAX_NVOIGT];

  PRINTF1("calc stress in ") PRINT_ARRAY("strain", strain_1, nvoi)
  ARRAY_COPY(strain_1, strain_mac, nvoi);

  int ierr = homog_get_strain_stress(strain_1, strain_aux, stress_1);

  for (int i = 0 ; i < nvoi ; i++) {

    PRINTF2("exp %d\n", i)
    ARRAY_COPY(strain_2, strain_mac, nvoi);

    strain_2[i] = strain_2[i] + HOMOGENIZE_DELTA_STRAIN;

    ierr = homog_get_strain_stress(strain_2, strain_aux, stress_2);

    for (int j = 0 ; j < nvoi ; j++)
      c_tangent[j*nvoi + i] = (stress_2[j] - stress_1[j]) / (strain_2[i] - strain_1[i]);

    if (flags.print_pvtu == true) {
      get_elem_properties();
      char filename[64];
      sprintf(filename,"micro_exp%d",i);
      ierr = micro_pvtu(filename);
      if (ierr != 0) {
	myio_printf(MICRO_COMM,"Problem writing vtu file\n");
	return ierr;
      }
    }
  }
  return ierr;
}

int homog_taylor(double *strain_mac, double *strain_ave, double *stress_ave)
{
  double *c_i = malloc(nvoi*nvoi * sizeof(double));
  double *c_m = malloc(nvoi*nvoi * sizeof(double));
  double *c = malloc(nvoi*nvoi * sizeof(double));
  double vol_i = 0.0, vol_m = 0.0;

  int ne_i = 0, ne_m = 0;
  for (int e = 0; e < mesh_struct.nelm; e++) {
    if (elem_type[e] == ID_FIBER) {
      vol_i += mesh_struct.vol_elm;
      ne_i++;
    } else {
      vol_m += mesh_struct.vol_elm;
      ne_m++;
    }
  }

  vi = vol_i / mesh_struct.vol;
  vm = vol_m / mesh_struct.vol;
  PRINTF2("vi = %lf \n", vi);
  PRINTF2("vm = %lf \n", vm);

  get_c_tan("FIBER" , -1, -1, NULL, c_i);
  get_c_tan("MATRIX", -1, -1, NULL, c_m);

  if (params.multis_method == MULTIS_MIXP) {

    for (int i = 0; i < nvoi; i++) {
      for (int j = 0; j < nvoi; j++)
	c[i*nvoi + j] = vi * c_i[i*nvoi + j] + vm * c_m[i*nvoi + j];
    }

  } else if (params.multis_method == MULTIS_MIXS) {

    double *c_mi = malloc(nvoi*nvoi*sizeof(double));
    double *c_ii = malloc(nvoi*nvoi*sizeof(double));
    gsl_matrix_view  gsl_c_m , gsl_c_i ;
    gsl_matrix_view  gsl_c_mi, gsl_c_ii;

    gsl_permutation *p = gsl_permutation_alloc(nvoi);

    gsl_c_i = gsl_matrix_view_array(c_i , nvoi, nvoi);
    gsl_c_m = gsl_matrix_view_array(c_m , nvoi, nvoi);
    gsl_c_ii = gsl_matrix_view_array(c_ii, nvoi, nvoi);
    gsl_c_mi = gsl_matrix_view_array(c_mi, nvoi, nvoi);

    int s;
    gsl_linalg_LU_decomp(&gsl_c_m.matrix, p, &s);
    gsl_linalg_LU_invert(&gsl_c_m.matrix, p, &gsl_c_mi.matrix);
    gsl_linalg_LU_decomp(&gsl_c_i.matrix, p, &s);
    gsl_linalg_LU_invert(&gsl_c_i.matrix, p, &gsl_c_ii.matrix);

    for (int i = 0; i < nvoi*nvoi; i++)
      c_i[i] = vi*c_ii[i] + vm*c_mi[i];

    gsl_linalg_LU_decomp(&gsl_c_i.matrix, p, &s);
    gsl_linalg_LU_invert(&gsl_c_i.matrix, p, &gsl_c_mi.matrix);

    for (int i = 0; i < nvoi; i++) {
      for (int j = 0; j < nvoi; j++)
	c[i*nvoi + j] = gsl_matrix_get( &gsl_c_mi.matrix, i, j );
    }

    gsl_permutation_free(p);
    free(c_mi);
    free(c_ii);
  }

  for (int i = 0; i < nvoi; i++) {
    strain_ave[i] = strain_mac[i];
    stress_ave[i] = 0.0;
    for (int j = 0; j < nvoi; j++)
      stress_ave[i] += c[i*nvoi + j] * strain_mac[j];
  }

  return 0;
}

int homog_fe2(double *strain_mac, double *strain_ave, double *stress_ave)
{
  clock_t start, end;
  double time_ass_b = 0.0, time_ass_A = 0.0, time_sol = 0.0;

  VecZeroEntries(x);

  params.non_linear_its = 0;
  params.residual_norm = params.non_linear_min_norm_tol * 10;
  while (params.non_linear_its < params.non_linear_max_its && params.residual_norm > params.non_linear_min_norm_tol) {

    save_event(MICRO_COMM, "ass_0");

    start = clock();
    assembly_b_petsc(&params.residual_norm, strain_mac);
    end = clock();
    time_ass_b = ((double) (end - start)) / CLOCKS_PER_SEC;
    PRINTF2(GREEN "|b| = %lf" NORMAL "\n", params.residual_norm);

    if (params.residual_norm < params.non_linear_min_norm_tol) break;

    VecScale(b, -1.0);

    start = clock();
    assembly_A_petsc();
    end = clock();
    time_ass_A = ((double) (end - start)) / CLOCKS_PER_SEC;

    save_event(MICRO_COMM, "ass_1");

    save_event(MICRO_COMM, "sol_0");
    start = clock();
    KSPSetOperators(ksp, A, A);
    KSPSolve(ksp, b, dx);
    end = clock();
    time_sol = ((double) (end - start)) / CLOCKS_PER_SEC;
    save_event(MICRO_COMM, "sol_1");

    VecAXPY(x, 1.0, dx);

    params.non_linear_its ++;
  }
  save_event(MICRO_COMM, "ass_1");

  PRINTF2(BLUE "time ass_b = %lf" NORMAL "\n", time_ass_b);
  PRINTF2(BLUE "time ass_A = %lf" NORMAL "\n", time_ass_A);
  PRINTF2(BLUE "time sol   = %lf" NORMAL "\n", time_sol);

  get_averages(strain_ave, stress_ave);

  return 0;
}

int strain_x_coord(double *strain, double *coord, double *u)
{
  if (dim == 2) {
    u[0] = strain[0]   * coord[0] + strain[2]/2 * coord[1];
    u[1] = strain[2]/2 * coord[0] + strain[1]   * coord[1];
  }
  return 0;
}

int assembly_b_petsc(double *norm, double *strain_mac)
{
  double *b_arr;
  VecZeroEntries(b);
  VecGetArray(b, &b_arr);

  int npe = mesh_struct.npe;
  int dim = mesh_struct.dim;
  double *res_elem = malloc(dim*npe*sizeof(double));

  for (int e = 0 ; e < mesh_struct.nelm ; e++) {

    ARRAY_SET_TO_ZERO(res_elem, npe*dim);
    mesh_struct_get_elem_indeces(&mesh_struct, e, elem_index);

    for (int gp = 0; gp < ngp; gp++) {

      get_strain(e, gp, strain_gp);
      get_stress(e, gp, strain_gp, stress_gp);

      for (int i = 0; i < npe*dim; i++) {
	for (int j = 0; j < nvoi; j++)
	  res_elem[i] += struct_bmat[j][i][gp] * stress_gp[j] * struct_wp[gp];
      }

    }

    for (int i = 0; i < npe*dim; i++ )
      b_arr[elem_index[i]] += res_elem[i];
  }

  double *x_arr;
  int nn = mesh_struct.nn;
  VecGetArray(x, &x_arr);
  if (dim == 2) {
    if (params.fe2_bc == BC_USTRAIN) {
      double displ[2];
      for (int n = 0; n < mesh_struct.nnods_boundary ; n++ ) {
	strain_x_coord(strain_mac, &mesh_struct.boundary_coord[n*dim], displ);
	for (int d = 0; d < dim ; d++)
	  b_arr[mesh_struct.boundary_indeces[n*dim + d]] = x_arr[mesh_struct.boundary_indeces[n*dim + d]] - displ[d];
      }
    }
    else if (params.fe2_bc == BC_PERIODIC) {

      /*       | fi              |
       *  r =  | fe - µ          |
       *       | u+ - u- - e x dl|
       */
      
      for (int n = 0; n < mesh_struct.ny - 2 ; n++) { // fe - µ en X0 y X1
	for (int d = 0; d < dim ; d++) {
	  b_arr[mesh_struct.nods_x0[n]*dim + d] -= x_arr[nn*dim + n*dim + d];
	  b_arr[mesh_struct.nods_x1[n]*dim + d] += x_arr[nn*dim + n*dim + d];
	}
      }
      for (int n = 0; n < mesh_struct.nx - 2 ; n++) { // fe - µ en Y0 y Y1
	for (int d = 0; d < dim ; d++) {
	  b_arr[mesh_struct.nods_y0[n]*dim + d] -= x_arr[nn*dim + n*dim + d];
	  b_arr[mesh_struct.nods_y1[n]*dim + d] += x_arr[(nn + mesh_struct.ny - 2)*dim + n*dim + d];
	}
      }

      double displ[2], coor[2];

      // we set displacements u - e.x in corners
      coor[0] = 0.0; coor[1] = 0.0; // x0y0
      strain_x_coord(strain_mac, coor, displ);
      for (int d = 0; d < dim ; d++)
	b_arr[mesh_struct.nod_x0y0*dim + d] = x_arr[mesh_struct.nod_x0y0*dim + d] - displ[d];

      coor[0] = mesh_struct.lx; coor[1] = 0.0; // x1y0
      strain_x_coord(strain_mac, coor, displ);
      for (int d = 0; d < dim ; d++)
	b_arr[mesh_struct.nod_x1y0*dim + d] = x_arr[mesh_struct.nod_x1y0*dim + d] - displ[d];

      coor[0] = mesh_struct.lx; coor[1] = mesh_struct.ly; // x1y1
      strain_x_coord(strain_mac, coor, displ);
      for (int d = 0; d < dim ; d++)
	b_arr[mesh_struct.nod_x1y1*dim + d] = x_arr[mesh_struct.nod_x1y1*dim + d] - displ[d];

      coor[0] = 0.0; coor[1] = mesh_struct.ly; // x0y1
      strain_x_coord(strain_mac, coor, displ);
      for (int d = 0; d < dim ; d++)
	b_arr[mesh_struct.nod_x0y1*dim + d] = x_arr[mesh_struct.nod_x0y1*dim + d] - displ[d];

      double delta_length[2];

      delta_length[0] = mesh_struct.lx; delta_length[1] = 0.0; // set u+ - u- - e*dx en x0/x1
      strain_x_coord(strain_mac, delta_length, displ);
      for (int n = 0; n < mesh_struct.ny - 2 ; n++) {
	for (int d = 0; d < dim ; d++)
	  b_arr[nn*dim + n*dim + d] = x_arr[mesh_struct.nods_x1[n]*dim + d] - x_arr[mesh_struct.nods_x0[n]*dim + d] - displ[d];
      }

      delta_length[0] = 0.0; delta_length[1] = mesh_struct.ly; // set u+ - u- - e*dx en y0/y1
      strain_x_coord(strain_mac, delta_length, displ);
      for (int n = 0; n < mesh_struct.nx - 2 ; n++) {
	for (int d = 0; d < dim ; d++)
	  b_arr[(nn + mesh_struct.ny - 2)*dim + n*dim + d] = x_arr[mesh_struct.nods_y1[n]*dim + d] - x_arr[mesh_struct.nods_y0[n]*dim + d] - displ[d];
      }
    }
  }
  VecRestoreArray(x, &x_arr);

  VecRestoreArray(b, &b_arr);
  VecNorm(b, NORM_2, norm);

  if (flags.print_vectors == true) {
    PetscViewer viewer;
    PetscViewerASCIIOpen(MICRO_COMM,"b.dat" ,&viewer); VecView(b ,viewer);
    PetscViewerASCIIOpen(MICRO_COMM,"x.dat" ,&viewer); VecView(x ,viewer);
    PetscViewerASCIIOpen(MICRO_COMM,"dx.dat",&viewer); VecView(dx,viewer);
  }

  return 0;
}

int assembly_A_petsc(void)
{
  MatZeroEntries(A);

  int npe = mesh_struct.npe;
  int dim = mesh_struct.dim;
  double *k_elem = malloc(dim*npe*dim*npe*sizeof(double));
  double *c = malloc(nvoi*nvoi*sizeof(double));

  for (int e = 0; e < mesh_struct.nelm; e++) {

    ARRAY_SET_TO_ZERO(k_elem, npe*dim*npe*dim)
    mesh_struct_get_elem_indeces(&mesh_struct, e, elem_index);

    for (int gp = 0; gp < ngp; gp++) {

      get_strain(e, gp, strain_gp);
      get_c_tan(NULL, e, gp, strain_gp, c);

      for (int i = 0 ; i < npe*dim; i++)
	for (int j = 0 ; j < npe*dim; j++)
	  for (int k = 0; k < nvoi ; k++)
	    for (int h = 0; h < nvoi; h++)
	      k_elem[ i*npe*dim + j] += struct_bmat[h][i][gp]*c[h*nvoi + k]*struct_bmat[k][j][gp]*struct_wp[gp];

    }
    MatSetValues(A, npe*dim, elem_index, npe*dim, elem_index, k_elem, ADD_VALUES);
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  if (params.fe2_bc == BC_USTRAIN) {
    MatZeroRowsColumns(A, mesh_struct.nnods_boundary * mesh_struct.dim, mesh_struct.boundary_indeces, 1.0, NULL, NULL);
  }
  else if (params.fe2_bc == BC_PERIODIC) {

    /*       | K   1 |
     * jac = |       |
     *       | 1   0 |
     */
    int index[3], row_index;
    double values[3];
    if (dim == 2) {
      index[0] = mesh_struct.nod_x0y0 * dim + 0;
      index[1] = mesh_struct.nod_x0y0 * dim + 1;
      MatZeroRows(A, 2, index, 1.0, NULL, NULL);
      index[0] = mesh_struct.nod_x1y0 * dim + 0;
      index[1] = mesh_struct.nod_x1y0 * dim + 1;
      MatZeroRows(A, 2, index, 1.0, NULL, NULL);
      index[0] = mesh_struct.nod_x1y1 * dim + 0;
      index[1] = mesh_struct.nod_x1y1 * dim + 1;
      MatZeroRows(A, 2, index, 1.0, NULL, NULL);
      index[0] = mesh_struct.nod_x0y1 * dim + 0;
      index[1] = mesh_struct.nod_x0y1 * dim + 1;
      MatZeroRows(A, 2, index, 1.0, NULL, NULL);

      values[0] = -1.0;
      values[1] = +1.0;

      int col_index;

      for (int n = 0; n < mesh_struct.ny - 2 ; n++) { // f + µ in x0
	for (int d = 0; d < dim ; d++){
	  row_index = mesh_struct.nods_x1[n]*dim + d;
	  col_index = mesh_struct.nn*dim + n*dim + d;
	  MatSetValues(A, 1, &row_index, 1, &col_index, &values[1], INSERT_VALUES);
	}
      }
      for (int n = 0; n < mesh_struct.ny - 2 ; n++) { // f - µ in x0
	for (int d = 0; d < dim ; d++){
	  row_index = mesh_struct.nods_x0[n]*dim + d;
	  col_index = mesh_struct.nn*dim + n*dim + d;
	  MatSetValues(A, 1, &row_index, 1, &col_index, &values[0], INSERT_VALUES);
	}
      }
      for (int n = 0; n < mesh_struct.nx - 2 ; n++) { // f + µ in y1
	for (int d = 0; d < dim ; d++){
	  row_index = mesh_struct.nods_y1[n]*dim + d;
	  col_index = (mesh_struct.nn + mesh_struct.ny - 2)*dim + n*dim + d;
	  MatSetValues(A, 1, &row_index, 1, &col_index, &values[1], INSERT_VALUES);
	}
      }
      for (int n = 0; n < mesh_struct.nx - 2 ; n++) { // f - µ in y0
	for (int d = 0; d < dim ; d++){
	  row_index = mesh_struct.nods_y0[n]*dim + d;
	  col_index = (mesh_struct.nn + mesh_struct.ny - 2)*dim + n*dim + d;
	  MatSetValues(A, 1, &row_index, 1, &col_index, &values[0], INSERT_VALUES);
	}
      }

      for (int n = 0; n < mesh_struct.ny - 2 ; n++) { // u+ - u- - e x dl en x0 x1
	for (int d = 0; d < dim ; d++){
	  row_index = mesh_struct.nn*dim + n*dim + d;
	  index[0] = mesh_struct.nods_x0[n]*dim + d;
	  index[1] = mesh_struct.nods_x1[n]*dim + d;
	  MatSetValues(A, 1, &row_index, 2, index, values, INSERT_VALUES);
	}
      }
      for (int n = 0; n < mesh_struct.nx - 2 ; n++) {
	for (int d = 0; d < dim ; d++){
	  row_index = (mesh_struct.nn + mesh_struct.ny - 2)*dim + n*dim + d;
	  index[0] = mesh_struct.nods_y0[n]*dim + d;
	  index[1] = mesh_struct.nods_y1[n]*dim + d;
	  MatSetValues(A, 1, &row_index, 2, index, values, INSERT_VALUES);
	}
      }
    }
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  if (flags.print_matrices == true) {
    PetscViewer viewer;
    PetscViewerASCIIOpen(MICRO_COMM,"A.dat" ,&viewer); MatView(A ,viewer);
  }
  return 0;
}

int get_averages(double *strain_ave, double *stress_ave)
{
  ARRAY_SET_TO_ZERO(strain_ave, nvoi);
  ARRAY_SET_TO_ZERO(stress_ave, nvoi);

  for (int e = 0 ; e < mesh_struct.nelm ; e++) {

    for (int gp = 0; gp < ngp; gp++) {

      get_strain(e, gp, strain_gp);
      get_stress(e, gp, strain_gp, stress_gp);

      for (int i = 0; i < nvoi; i++) {
	stress_ave[i] += stress_gp[i] * struct_wp[gp];
	strain_ave[i] += strain_gp[i] * struct_wp[gp];
      }
    }
  }

  for (int i = 0; i < nvoi; i++) {
    stress_ave[i] /= mesh_struct.vol;
    strain_ave[i] /= mesh_struct.vol;
  }
  return 0;
}

int get_elem_properties(void)
{
  for (int e = 0; e < mesh_struct.nelm; e++) {

    double strain_aux[MAX_NVOIGT];
    double stress_aux[MAX_NVOIGT];
    ARRAY_SET_TO_ZERO(strain_aux, nvoi);
    ARRAY_SET_TO_ZERO(stress_aux, nvoi);

    for (int gp = 0 ; gp < ngp ; gp++) {

      get_strain( e , gp, strain_gp );
      get_stress( e , gp, strain_gp, stress_gp );
      for (int v = 0 ; v < nvoi ; v++ ) {
	strain_aux[v] += strain_gp[v] * struct_wp[gp];
	stress_aux[v] += stress_gp[v] * struct_wp[gp];
      }

    }
    for (int v = 0 ; v < nvoi ; v++) {
      elem_strain[ e*nvoi + v ] = strain_aux[v] / mesh_struct.vol_elm;
      elem_stress[ e*nvoi + v ] = stress_aux[v] / mesh_struct.vol_elm;
    }
  }
  return 0;
}

int get_strain(int e, int gp, double *strain_gp)
{
  double *x_arr;
  VecGetArray(x, &x_arr);

  mesh_struct_get_elem_indeces(&mesh_struct, e, elem_index);

  for (int i = 0; i < mesh_struct.npe * mesh_struct.dim; i++)
    elem_disp[i] = x_arr[elem_index[i]];

  for (int v = 0; v < nvoi; v++) {
    strain_gp[v] = 0.0;
    for (int i = 0; i < mesh_struct.npe*dim; i++)
      strain_gp[v] += struct_bmat[v][i][gp] * elem_disp[i];
  }

  return 0;
}


int get_stress(int e, int gp, double *strain_gp, double *stress_gp)
{
  char *word_to_search;

  switch (elem_type[e]) {

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
  node_list_t *pm = material_list.head;

  while (pm != NULL) {
    mat_p = (material_t *)pm->data;
    if (strcmp(mat_p->name, word_to_search) == 0) break;
    pm = pm->next;
  }

  if (pm == NULL) return 1;

  return material_get_stress(mat_p, dim, strain_gp, stress_gp);
}

int get_c_tan(const char *name, int e, int gp, double *strain_gp, double *c_tan)
{
  char *word_to_search;

  if (name != NULL)
    word_to_search = strdup(name);
  else if ( elem_type[e] == ID_FIBER )
    word_to_search = strdup("FIBER");
  else if ( elem_type[e] == ID_MATRIX )
    word_to_search = strdup("MATRIX");

  material_t  *mat_p;
  node_list_t *pm = material_list.head;

  while (pm != NULL){
    mat_p = (material_t *)pm->data;
    if (strcmp(mat_p->name, word_to_search) == 0) break;
    pm = pm->next;
  }
  if (pm == NULL) return 1;

  return material_get_c_tang(mat_p, dim, strain_gp, c_tan);
}

int get_elem_centroid(int e, int dim, double *centroid)
{
  if (dim == 2) {
    centroid[0] = (e % mesh_struct.nex + 0.5) * mesh_struct.hx;
    centroid[1] = (e / mesh_struct.nex + 0.5) * mesh_struct.hy;
  }
  return 0;
}
