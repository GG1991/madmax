#include "macro.h"


int boundary_read(void)
{
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

    for(int i = 0; i < n; i++){

      int da = 0;
      int *p = bsearch(&ix[i], mesh.local_nods, mesh.nnods_local, sizeof(int), mesh_cmpfunc);

      for(int d = 0; d < dim; d++){
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


int boundary_update(double time)
{
  double value;
  node_list_t *pn = boundary_list.head;

  while (pn != NULL) {

    mesh_boundary_t *bou = (mesh_boundary_t *)pn->data;
    function_t *function = NULL;

    for (int d = 0; d < dim; d++) {
      function_get_from_list(bou->fnum[d], &function_list, &function);
      function_eval(time, function, &value);
      for (int i = 0 ; i < bou->ndir ; i++) bou->dir_val[i*(bou->ndirpn) + d] = value;
    }
    pn = pn->next;
  }

  return 0;
}


int boundary_setx(void)
{ 
  Vec x_loc; double *x_arr;

  VecGhostGetLocalForm(x, &x_loc);
  VecGetArray(x_loc, &x_arr);

  node_list_t * pn = boundary_list.head;
  while (pn != NULL) {
    mesh_boundary_t *bou = (mesh_boundary_t *)pn->data;
    for (int i = 0; i < bou->ndirix; i++) x_arr[bou->dir_loc_ixs[i]] = bou->dir_val[i];
    pn = pn->next;
  }

  VecRestoreArray(x_loc, &x_arr);
  VecGhostRestoreLocalForm(x, &x_loc);

  VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);

  return 0;
}
