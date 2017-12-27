#include "macro.h"


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
