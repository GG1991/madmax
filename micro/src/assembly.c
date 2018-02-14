#include "micro.h"

int assembly_jac(void)
{
  switch (params.solver) {
    case SOL_PETSC:
      assembly_jac_petsc();
      break;
    case SOL_ELL:
      assembly_jac_ell();
      break;
    default:
      return 1;
  }
  return 0;
}

int assembly_res(double *norm, double *strain_mac)
{
  switch (params.solver) {
    case SOL_PETSC:
      assembly_res_petsc(norm, strain_mac);
      break;
    case SOL_ELL:
      assembly_res_ell(norm, strain_mac);
      break;
    default:
      return 1;
  }
  return 0;
}

int assembly_res_ell(double *norm, double *strain_mac)
{
  int npe = mesh_struct.npe;
  int dim = mesh_struct.dim;
  int nn = mesh_struct.nn;
  int nelm = mesh_struct.nelm;
  double *res_e = malloc((dim*npe) * sizeof(double));

  for (int i = 0 ; i < (nn*dim) ; i++) res_ell[i] = 0.0;

  for (int e = 0 ; e < nelm ; e++) {
    for (int i = 0 ; i < (npe*dim) ; i++) res_e[i] = 0.0;
    mesh_struct_get_elem_indeces(&mesh_struct, e, elem_index);
    
    for (int gp = 0; gp < ngp; gp++) { // integration
      get_strain(e, gp, strain_gp);
      get_stress(e, gp, strain_gp, stress_gp);
      for (int i = 0; i < npe*dim; i++) {
	for (int j = 0; j < nvoi; j++)
	  res_e[i] += struct_bmat[j][i][gp] * stress_gp[j] * struct_wp[gp];
      }
    }
    
    for (int i = 0; i < (npe*dim); i++ ) res_ell[elem_index[i]] += res_e[i]; // assembly
  }

  if (params.fe2_bc == BC_USTRAIN) {
    for (int d = 0; d < dim ; d++) {
      for (int n = 0; n < mesh_struct.ny - 2 ; n++) {
	res_ell[mesh_struct.nods_x0[n]*dim + d] = 0.0;
	res_ell[mesh_struct.nods_x1[n]*dim + d] = 0.0;
      }
      for (int n = 0; n < mesh_struct.nx - 2 ; n++) {
	res_ell[mesh_struct.nods_y0[n]*dim + d] = 0.0;
	res_ell[mesh_struct.nods_y1[n]*dim + d] = 0.0;
      }
      res_ell[mesh_struct.nod_x0y0*dim + d] = 0.0;
      res_ell[mesh_struct.nod_x1y0*dim + d] = 0.0;
      res_ell[mesh_struct.nod_x1y1*dim + d] = 0.0;
      res_ell[mesh_struct.nod_x0y1*dim + d] = 0.0;
    }
  }

  *norm = 0;
  for (int i = 0 ; i < (nn*dim) ; i++)
    *norm += res_ell[i] * res_ell[i];
  *norm = sqrt(*norm);

  return 0;
}

int assembly_jac_ell(void)
{
  int npe = mesh_struct.npe;
  int dim = mesh_struct.dim;
  int nelm = mesh_struct.nelm;
  double *jac_e = malloc((npe*dim*npe*dim) * sizeof(double));
  double *ctang = malloc(nvoi*nvoi*sizeof(double));

  ell_set_zero_mat(&jac_ell);
  for (int e = 0 ; e < nelm ; e++) {
    for (int i = 0 ; i < (npe*dim*npe*dim) ; i++) jac_e[i] = 0.0;
    mesh_struct_get_elem_indeces(&mesh_struct, e, elem_index);

    for (int gp = 0; gp < ngp; gp++) { // integration

      get_strain(e, gp, strain_gp);
      get_c_tan(NULL, e, gp, strain_gp, ctang);

      for (int i = 0 ; i < npe*dim; i++)
	for (int j = 0 ; j < npe*dim; j++)
	  for (int k = 0; k < nvoi ; k++)
	    for (int h = 0; h < nvoi; h++)
	      jac_e[ i*(npe*dim) + j] += struct_bmat[h][i][gp] * ctang[h*nvoi + k] * struct_bmat[k][j][gp] * struct_wp[gp];

    }
    ell_add_vals(&jac_ell, elem_index, npe*dim, elem_index, npe*dim, jac_e); // assembly
  }

  if (params.fe2_bc == BC_USTRAIN) {
    for (int d = 0; d < dim ; d++) {
      for (int n = 0; n < mesh_struct.ny - 2 ; n++) {
	ell_set_zero_row (&jac_ell, mesh_struct.nods_x0[n]*dim + d, 1.0);
	ell_set_zero_col (&jac_ell, mesh_struct.nods_x0[n]*dim + d, 1.0);
	ell_set_zero_row (&jac_ell, mesh_struct.nods_x1[n]*dim + d, 1.0);
	ell_set_zero_col (&jac_ell, mesh_struct.nods_x1[n]*dim + d, 1.0);
      }
      for (int n = 0; n < mesh_struct.nx - 2 ; n++) {
	ell_set_zero_row (&jac_ell, mesh_struct.nods_y0[n]*dim + d, 1.0);
	ell_set_zero_col (&jac_ell, mesh_struct.nods_y0[n]*dim + d, 1.0);
	ell_set_zero_row (&jac_ell, mesh_struct.nods_y1[n]*dim + d, 1.0);
	ell_set_zero_col (&jac_ell, mesh_struct.nods_y1[n]*dim + d, 1.0);
      }
      ell_set_zero_row (&jac_ell, mesh_struct.nod_x0y0*dim + d, 1.0);
      ell_set_zero_col (&jac_ell, mesh_struct.nod_x0y0*dim + d, 1.0);
      ell_set_zero_row (&jac_ell, mesh_struct.nod_x1y0*dim + d, 1.0);
      ell_set_zero_col (&jac_ell, mesh_struct.nod_x1y0*dim + d, 1.0);
      ell_set_zero_row (&jac_ell, mesh_struct.nod_x1y1*dim + d, 1.0);
      ell_set_zero_col (&jac_ell, mesh_struct.nod_x1y1*dim + d, 1.0);
      ell_set_zero_row (&jac_ell, mesh_struct.nod_x0y1*dim + d, 1.0);
      ell_set_zero_col (&jac_ell, mesh_struct.nod_x0y1*dim + d, 1.0);
    }
  }
  return 0;
}

int assembly_res_petsc(double *norm, double *strain_mac)
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
    for (int i = 0; i < npe*dim; i++) b_arr[elem_index[i]] += res_elem[i];
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
    else if (params.fe2_bc == BC_PER_LM) {

      /*       | fi              |      |u      |
       *  r =  | fe - µ          |  x = |µ en y0| 
       *       | u+ - u- - e x dl|      |µ en x0|
       */
      
      for (int n = 0; n < mesh_struct.nx - 2 ; n++) { // fe - µ en Y0 y Y1
	for (int d = 0; d < dim ; d++) {
	  b_arr[mesh_struct.nods_y0[n]*dim + d] -= x_arr[nn*dim + n*dim + d];
	  b_arr[mesh_struct.nods_y1[n]*dim + d] += x_arr[nn*dim + n*dim + d];
	}
      }
      for (int n = 0; n < mesh_struct.ny - 2 ; n++) { // fe - µ en X0 y X1
	for (int d = 0; d < dim ; d++) {
	  b_arr[mesh_struct.nods_x0[n]*dim + d] -= x_arr[(nn + mesh_struct.nx - 2)*dim + n*dim + d];
	  b_arr[mesh_struct.nods_x1[n]*dim + d] += x_arr[(nn + mesh_struct.nx - 2)*dim + n*dim + d];
	}
      }

      // we set displacements u = e.x in corners
      for (int d = 0; d < dim ; d++) // x0y0
	b_arr[mesh_struct.nod_x0y0*dim + d] = 0.0;
      for (int d = 0; d < dim ; d++) // x1y0
	b_arr[mesh_struct.nod_x1y0*dim + d] = 0.0;
      for (int d = 0; d < dim ; d++) // x1y1
	b_arr[mesh_struct.nod_x1y1*dim + d] = 0.0;
      for (int d = 0; d < dim ; d++) // x0y1
	b_arr[mesh_struct.nod_x0y1*dim + d] = 0.0;

      double delta_length[2], displ[2];

      delta_length[0] = 0.0; delta_length[1] = mesh_struct.ly; // set u+ - u- - e*dx en y0/y1
      strain_x_coord(strain_mac, delta_length, displ);
      for (int n = 0; n < mesh_struct.nx - 2 ; n++) {
	for (int d = 0; d < dim ; d++)
	  b_arr[nn*dim + n*dim + d] = x_arr[mesh_struct.nods_y1[n]*dim + d] - x_arr[mesh_struct.nods_y0[n]*dim + d] - displ[d];
      }

      delta_length[0] = mesh_struct.lx; delta_length[1] = 0.0; // set u+ - u- - e*dx en x0/x1
      strain_x_coord(strain_mac, delta_length, displ);
      for (int n = 0; n < mesh_struct.ny - 2 ; n++) {
	for (int d = 0; d < dim ; d++)
	  b_arr[(nn + mesh_struct.nx - 2)*dim + n*dim + d] = x_arr[mesh_struct.nods_x1[n]*dim + d] - x_arr[mesh_struct.nods_x0[n]*dim + d] - displ[d];
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


int assembly_jac_petsc(void)
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
    MatZeroRowsColumns(A, mesh_struct.nnods_boundary*mesh_struct.dim, mesh_struct.boundary_indeces, 1.0, NULL, NULL);
  }
  else if (params.fe2_bc == BC_PER_LM) {

    /*       | K   1 |
     * jac = |       |
     *       | 1   0 |
     */
    int index[3], row_index;
    double values[3];
    if (dim == 2) {

      // corners
      index[0] = mesh_struct.nod_x0y0 * dim + 0;
      index[1] = mesh_struct.nod_x0y0 * dim + 1;
//      MatZeroRows(A, 2, index, 1.0, NULL, NULL);
      MatZeroRowsColumns(A, 2, index, 1.0, NULL, NULL);
      index[0] = mesh_struct.nod_x1y0 * dim + 0;
      index[1] = mesh_struct.nod_x1y0 * dim + 1;
//      MatZeroRows(A, 2, index, 1.0, NULL, NULL);
      MatZeroRowsColumns(A, 2, index, 1.0, NULL, NULL);
      index[0] = mesh_struct.nod_x1y1 * dim + 0;
      index[1] = mesh_struct.nod_x1y1 * dim + 1;
//      MatZeroRows(A, 2, index, 1.0, NULL, NULL);
      MatZeroRowsColumns(A, 2, index, 1.0, NULL, NULL);
      index[0] = mesh_struct.nod_x0y1 * dim + 0;
      index[1] = mesh_struct.nod_x0y1 * dim + 1;
//      MatZeroRows(A, 2, index, 1.0, NULL, NULL);
      MatZeroRowsColumns(A, 2, index, 1.0, NULL, NULL);

      values[0] = -1.0;
      values[1] = +1.0;

      int col_index;

      for (int n = 0; n < mesh_struct.nx - 2 ; n++) { // f + µ in y1
	for (int d = 0; d < dim ; d++){
	  row_index = mesh_struct.nods_y1[n]*dim + d;
	  col_index = mesh_struct.nn*dim + n*dim + d;
	  MatSetValues(A, 1, &row_index, 1, &col_index, &values[1], INSERT_VALUES);
	}
      }
      for (int n = 0; n < mesh_struct.nx - 2 ; n++) { // f - µ in y0
	for (int d = 0; d < dim ; d++){
	  row_index = mesh_struct.nods_y0[n]*dim + d;
	  col_index = mesh_struct.nn*dim + n*dim + d;
	  MatSetValues(A, 1, &row_index, 1, &col_index, &values[0], INSERT_VALUES);
	}
      }
      for (int n = 0; n < mesh_struct.ny - 2 ; n++) { // f + µ in x0
	for (int d = 0; d < dim ; d++){
	  row_index = mesh_struct.nods_x1[n]*dim + d;
	  col_index = (mesh_struct.nn + mesh_struct.nx - 2)*dim + n*dim + d;
	  MatSetValues(A, 1, &row_index, 1, &col_index, &values[1], INSERT_VALUES);
	}
      }
      for (int n = 0; n < mesh_struct.ny - 2 ; n++) { // f - µ in x0
	for (int d = 0; d < dim ; d++){
	  row_index = mesh_struct.nods_x0[n]*dim + d;
	  col_index = (mesh_struct.nn + mesh_struct.nx - 2)*dim + n*dim + d;
	  MatSetValues(A, 1, &row_index, 1, &col_index, &values[0], INSERT_VALUES);
	}
      }

      double zero = 0.0;
      for (int n = 0; n < mesh_struct.nx - 2 ; n++) { // u+ - u- - e x dl en y0 y1
	for (int d = 0; d < dim ; d++){
	  row_index = mesh_struct.nn*dim + n*dim + d;
	  index[0] = mesh_struct.nods_y0[n]*dim + d;
	  index[1] = mesh_struct.nods_y1[n]*dim + d;
	  MatSetValues(A, 1, &row_index, 2, index, values, INSERT_VALUES);
	  MatSetValues(A, 1, &row_index, 1, &row_index, &zero, INSERT_VALUES);
	}
      }
      for (int n = 0; n < mesh_struct.ny - 2 ; n++) { // u+ - u- - e x dl en x0 x1
	for (int d = 0; d < dim ; d++){
	  row_index = (mesh_struct.nn + mesh_struct.nx - 2)*dim + n*dim + d;
	  index[0] = mesh_struct.nods_x0[n]*dim + d;
	  index[1] = mesh_struct.nods_x1[n]*dim + d;
	  MatSetValues(A, 1, &row_index, 2, index, values, INSERT_VALUES);
	  MatSetValues(A, 1, &row_index, 1, &row_index, &zero, INSERT_VALUES);
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

