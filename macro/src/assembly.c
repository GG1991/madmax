#include "macro.h"

#ifdef PETSC

int assembly_b(double *norm) {

  double *wp;
  double *b_arr;
  Vec b_loc;
  int ierr;

  VecGhostGetLocalForm(b, &b_loc);
  VecGetArray(b_loc, &b_arr);

  ARRAY_SET_TO_ZERO(b_arr, mesh.nnods_local_ghost*dim)

  for (int e = 0 ; e < mesh.nelm_local ; e++) {

    assembly_get_local_elem_index(e, loc_elem_index);
    int npe = mesh.eptr[e+1] - mesh.eptr[e];
    int ngp = npe;

    ARRAY_SET_TO_ZERO(res_elem, npe*dim)

    ierr = assembly_get_dsh(e, loc_elem_index, dsh, detj);
    ierr = assembly_get_bmat(e, dsh, bmat);
    ierr = assembly_get_wp(dim, npe, &wp);

    for (int gp = 0; gp < ngp ; gp++) {

      if (detj[gp] < 0.0) flag_neg_detj = 1;
      detj[gp] = fabs(detj[gp]);

      assembly_get_strain(e, gp, loc_elem_index, dsh, bmat, strain_gp);
      assembly_get_stress(e, gp, strain_gp, stress_gp);

      for (int i = 0 ; i < npe*dim ; i++) {
	for (int j = 0; j < nvoi ; j++)
	  res_elem[i] += bmat[j][i][gp]*stress_gp[j]*wp[gp]*detj[gp];
      }
    }

    for (int i = 0 ; i < npe*dim ; i++)
      b_arr[loc_elem_index[i]] += res_elem[i];

  }

  VecRestoreArray(b_loc, &b_arr);
  VecGhostRestoreLocalForm(b, &b_loc);

  VecGhostUpdateBegin(b, ADD_VALUES, SCATTER_REVERSE);
  VecGhostUpdateEnd(b, ADD_VALUES, SCATTER_REVERSE);
  VecGhostUpdateBegin(b, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(b, INSERT_VALUES, SCATTER_FORWARD);

  if (flag_neg_detj == 1)
    myio_printf(MACRO_COMM, "MACRO: warning negative jacobian detected\n");

  VecGhostGetLocalForm(b, &b_loc);
  VecGetArray(b_loc, &b_arr);

  node_list_t *pn = boundary_list.head;
  while (pn) {
    mesh_boundary_t *bou = (mesh_boundary_t *)pn->data;
    for (int i = 0 ; i < bou->ndirix ; i++)
      b_arr[bou->dir_loc_ixs[i]] = 0.0;
    pn = pn->next;
  }

  VecRestoreArray(b_loc, &b_arr);
  VecGhostRestoreLocalForm(b, &b_loc);
  VecGhostUpdateBegin(b, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(b, INSERT_VALUES, SCATTER_FORWARD);
  VecScale(b, -1.0);

  if (flags.print_vectors == true) {
    PetscViewer viewer;
    PetscViewerASCIIOpen(MACRO_COMM, "b.dat", &viewer); VecView(b ,viewer);
  }

  VecNorm(b, NORM_2, norm);

  return ierr;
}


int assembly_AM(void) {

  MatZeroEntries(A);
  MatZeroEntries(M);

  int ierr;
  double rho_gp;
  double **sh;
  double *wp;

  for (int e = 0 ; e < mesh.nelm_local ; e++ ) {

    int npe = mesh.eptr[e+1] - mesh.eptr[e];
    int ngp = npe;

    ARRAY_SET_TO_ZERO(k_elem, npe*dim*npe*dim);
    ARRAY_SET_TO_ZERO(m_elem, npe*dim*npe*dim);

    assembly_get_local_elem_index(e, loc_elem_index);
    assembly_get_global_elem_index(e, glo_elem_index);

    assembly_get_sh(dim, npe, &sh);
    assembly_get_dsh(e, loc_elem_index, dsh, detj);
    assembly_get_bmat(e, dsh, bmat);
    assembly_get_wp(dim, npe, &wp);

    for (int gp = 0; gp < ngp ; gp++) {

      detj[gp] = fabs(detj[gp]);

      ierr = assembly_get_strain(e , gp, loc_elem_index, dsh, bmat, strain_gp);
      ierr = assembly_get_c_tan(NULL, e, gp, strain_gp, c); if (ierr != 0) return 1;
      ierr = assembly_get_rho(NULL, e, &rho_gp); if (ierr != 0) return 1;

      for (int i = 0 ; i < npe*dim ; i++)
	for (int j = 0 ; j < npe*dim ; j++)
	  for (int k = 0; k < nvoi ; k++)
	    for (int h = 0; h < nvoi ; h++)
	      k_elem[i*npe*dim + j] += bmat[h][i][gp]*c[h*nvoi + k]*bmat[k][j][gp]*wp[gp]*detj[gp];

      for (int d = 0 ; d < dim ; d++)
	for (int i = 0 ; i < npe; i++)
	  for (int j = 0 ; j < npe; j++)
	    m_elem[(i*dim)*(npe*dim) + j*dim + (d*dim*npe + d)] += rho_gp*sh[i][gp]*sh[j][gp]*wp[gp]*detj[gp];

    }
    MatSetValues(A, npe*dim, glo_elem_index, npe*dim, glo_elem_index, k_elem, ADD_VALUES);
    MatSetValues(M, npe*dim, glo_elem_index, npe*dim, glo_elem_index, m_elem, ADD_VALUES);

  }

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);

  node_list_t *pn = boundary_list.head;
  while (pn != NULL) {
    mesh_boundary_t * bou = (mesh_boundary_t * )pn->data;
    MatZeroRowsColumns(M, bou->ndirix, bou->dir_glo_ixs, 1.0, NULL, NULL);
    MatZeroRowsColumns(A, bou->ndirix, bou->dir_glo_ixs, 1.0, NULL, NULL);
    pn = pn->next;
  }
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  if (flags.print_vectors == true) {
    PetscViewer viewer;
    PetscViewerASCIIOpen(MACRO_COMM, "M.dat", &viewer); MatView(M, viewer);
    PetscViewerASCIIOpen(MACRO_COMM, "A.dat", &viewer); MatView(A, viewer);
  }

  return ierr;
}


int assembly_A(void) {

  MatZeroEntries(A);

  int ierr;
  double *wp;

  for (int e = 0 ; e < mesh.nelm_local ; e++ ) {

    int npe = mesh.eptr[e+1] - mesh.eptr[e];
    int ngp = npe;

    ARRAY_SET_TO_ZERO(k_elem, npe*dim*npe*dim);

    assembly_get_local_elem_index (e, loc_elem_index);
    assembly_get_global_elem_index(e, glo_elem_index);

    assembly_get_dsh(e, loc_elem_index, dsh, detj);
    assembly_get_bmat(e, dsh, bmat);
    assembly_get_wp(dim, npe, &wp);

    for (int gp = 0; gp < ngp ; gp++) {

      detj[gp] = fabs(detj[gp]);

      assembly_get_strain(e, gp, loc_elem_index, dsh, bmat, strain_gp);

      ierr = assembly_get_c_tan(NULL , e , gp , strain_gp , c); if (ierr != 0) return ierr;

      for (int i = 0 ; i < npe*dim ; i++)
	for (int j = 0 ; j < npe*dim ; j++)
	  for (int k = 0; k < nvoi ; k++)
	    for (int h = 0; h < nvoi ; h++)
	      k_elem[ i*npe*dim + j] += bmat[h][i][gp]*c[h*nvoi + k]*bmat[k][j][gp]*wp[gp]*detj[gp];

    }
    MatSetValues(A, npe*dim, glo_elem_index, npe*dim, glo_elem_index, k_elem, ADD_VALUES);
  }

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  node_list_t *pn = boundary_list.head;
  while (pn != NULL) {
    mesh_boundary_t *bou = (mesh_boundary_t *)pn->data;
    MatZeroRowsColumns(A, bou->ndirix, bou->dir_glo_ixs, 1.0, NULL, NULL);
    pn = pn->next;
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  if (flags.print_vectors == true) {
    PetscViewer viewer;
    PetscViewerASCIIOpen(MACRO_COMM, "A.dat", &viewer); MatView(A ,viewer);
  }

  return ierr;
}

#endif


int assembly_get_dsh(int e, int *loc_elem_index, double ***dsh, double *detj)
{
  double ***dsh_master;
  int npe = mesh.eptr[e+1] - mesh.eptr[e];
  int ngp = npe;

  for (int i = 0 ; i < npe*dim ; i++)
    elem_coor[i] = mesh.coord_local[loc_elem_index[i]];

  fem_get_dsh_master(npe, dim, &dsh_master);

  for (int gp = 0; gp < ngp ; gp++) {
    fem_calc_jac(dim, npe, gp, elem_coor, dsh_master, jac);
    fem_invjac(dim, jac, jac_inv, &detj[gp]);
    fem_trans_dsh(dim, npe, gp, jac_inv, dsh_master, dsh);
  }

  return 0;
}


int assembly_get_bmat(int e, double ***dsh, double ***bmat)
{
  int npe = mesh.eptr[e+1] - mesh.eptr[e];
  int ngp = npe;

  if (dim == 2) {
    for (int i = 0 ; i < npe ; i++) {
      for (int gp = 0; gp < ngp ; gp++) {
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


int assembly_get_sh(int dim, int npe, double ***sh)
{
  return fem_get_sh(npe, dim, sh);
}


int assembly_get_wp(int dim, int npe, double **wp)
{
  return fem_get_wp(npe, dim, wp);
}


int assembly_get_global_elem_index(int e, int *glo_elem_index)
{
  int  npe = mesh.eptr[e+1] - mesh.eptr[e];
  for (int n = 0 ; n < npe ; n++)
    for (int d = 0 ; d < dim ; d++)
      glo_elem_index[n*dim + d] = mesh.local_to_global[mesh.eind[mesh.eptr[e] + n]]*dim + d;
  return 0;
}


int assembly_get_local_elem_index(int e, int *loc_elem_index)
{
  int  npe = mesh.eptr[e+1] - mesh.eptr[e];
  for (int n = 0 ; n < npe ; n++)
    for (int d = 0 ; d < dim ; d++)
      loc_elem_index[n*dim + d] = mesh.eind[mesh.eptr[e] + n]*dim + d;
  return 0;
}


int assembly_get_strain(int e , int gp, int *loc_elem_index, double ***dsh_gp,  double ***bmat, double *strain_gp)
{
  double *x_arr; Vec x_loc;
  VecGhostGetLocalForm(x, &x_loc);
  VecGetArray(x_loc, &x_arr);

  int  npe = mesh.eptr[e+1] - mesh.eptr[e];
  for (int i = 0 ; i < npe*dim ; i++)
    elem_disp[i] = x_arr[loc_elem_index[i]];

  VecRestoreArray(x_loc , &x_arr);
  VecGhostRestoreLocalForm(x, &x_loc);

  for (int v = 0; v < nvoi ; v++ ) {
    strain_gp[v] = 0.0;
    for (int i = 0 ; i < npe*dim ; i++ )
      strain_gp[v] += bmat[v][i][gp] * elem_disp[i];
    strain_gp[v] = ( fabs(strain_gp[v]) < 1.0e-6 ) ? 0.0 : strain_gp[v];
  }

  return 0;
}


int assembly_get_stress(int e, int gp, double *strain_gp, double *stress_gp)
{
  char name_s[64];
  material_t *mat_p;
  assembly_get_mat_name(mesh.elm_id[e], name_s);

  node_list_t *pn = material_list.head;
  while (pn != NULL) {
    mat_p = (material_t *)pn->data;
    if (strcmp(name_s, mat_p->name) == 0) break;
    pn = pn->next;
  }
  if (pn == NULL) {
    myio_printf(MACRO_COMM, "Material %s corresponding to element %d not found on material list\n", name_s, e);
    return 1;
  }

  if (mat_p->type_id == MAT_MICRO) {

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


int assembly_get_c_tan(const char *name, int e, int gp, double *strain_gp, double *c_tan)
{
  char name_s[64];
  material_t *mat_p;
  assembly_get_mat_name(mesh.elm_id[e], name_s);

  node_list_t *pn = material_list.head;
  while (pn != NULL) {
    mat_p = (material_t *)pn->data;
    if (strcmp(name_s, mat_p->name) == 0) break;
    pn = pn->next;
  }
  if (pn == NULL) {
    myio_printf(MACRO_COMM, "Material %s corresponding to element %d not found on material list\n", name_s, e );
    return 1;
  }

  if (mat_p->type_id == MAT_MICRO) {

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


int assembly_get_rho(const char *name, int e, double *rho)
{
  char name_s[64];
  int ierr;

  assembly_get_mat_name(mesh.elm_id[e], name_s);

  material_t *mat_p;
  node_list_t *pn = material_list.head;
  while (pn != NULL) {
    mat_p = (material_t *)pn->data;
    if (strcmp(name_s, mat_p->name) == 0) break;
    pn = pn->next;
  }
  if (pn == NULL) {
    myio_printf(MACRO_COMM, "Material %s corresponding to element %d not found on material list\n", name_s, e );
    return 1;
  }

  if (mat_p->type_id == MAT_MICRO) {

    message.action = ACTION_MICRO_CALC_RHO;
    ierr = comm_macro_send(&message, &comm);

    ierr = comm_macro_recv(&message, &comm);
    *rho = message.rho;
  }
  else
    material_get_rho( mat_p, dim, rho );

  return ierr;
}


int assembly_get_mat_name(int id, char *name_s)
{
  physical_t *phy_p;
  node_list_t *pn = physical_list.head;
  while (pn != NULL) {
    phy_p = ( physical_t * )pn->data;
    if ( id == phy_p->id ) break;
    pn = pn->next;
  }
  if (pn == NULL) return 1;

  strcpy(name_s, phy_p->name);

  return 0;
}


int assembly_get_elem_properties(void)
{
  double *strain_aux = malloc(nvoi*sizeof(double));
  double *stress_aux = malloc(nvoi*sizeof(double));
  double *wp;

  for (int e = 0 ; e < mesh.nelm_local ; e++) {

    int npe = mesh.eptr[e+1] - mesh.eptr[e];
    int ngp = npe;
    double vol_elem = 0.0;

    for (int v = 0 ; v < nvoi ; v++)
      strain_aux[v] = stress_aux[v] = 0.0;

    assembly_get_local_elem_index(e, loc_elem_index);

    assembly_get_dsh(e, loc_elem_index, dsh, detj);
    assembly_get_bmat(e, dsh, bmat);
    assembly_get_wp(dim, npe, &wp);

    for (int gp = 0 ; gp < ngp ; gp++) {

      detj[gp] = fabs(detj[gp]);

     assembly_get_strain(e, gp, loc_elem_index, dsh, bmat, strain_gp);
     assembly_get_stress(e, gp, strain_gp, stress_gp);
      for (int v = 0 ; v < nvoi ; v++) {
	strain_aux[v] += strain_gp[v] * detj[gp] * wp[gp];
	stress_aux[v] += stress_gp[v] * detj[gp] * wp[gp];
      }
      vol_elem += detj[gp] * wp[gp];
    }
    for (int v = 0 ; v < nvoi ; v++) {
      elem_strain[e*nvoi + v] = strain_aux[v] / vol_elem;
      elem_stress[e*nvoi + v] = stress_aux[v] / vol_elem;
    }

    physical_t * phy;
    node_list_t * pn = physical_list.head;
    while (pn != NULL) {
      phy = pn->data;
      if (phy->id == mesh.elm_id[e]) break;
      pn = pn->next;
    }
    if (pn == NULL) return 1;

    int type = 0;
    pn = material_list.head;
    while (pn != NULL) {
      material_t *mat = pn->data;
      if (strcmp(phy->name, mat->name) == 0) break;
      pn = pn->next;
      type ++;
    }
    if (pn == NULL) return 1;

    elem_type[e] = type;
  }

  return 0;
}
