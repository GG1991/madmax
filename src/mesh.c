#include "mesh.h"


int mesh_fill_boundary_list_from_command_line(command_line_t *command_line, list_t *boundary_list, mesh_t *mesh) {

  bool found;
  int num_string_found;
  char string_arr[MAX_NUM_OF_BOUNDARIES][128];

  myio_comm_line_get_string_array(command_line, "-boundary", string_arr, MAX_NUM_OF_BOUNDARIES, &num_string_found, &found);

  if (found == false || num_string_found == 0)
    return 1;

  list_init(boundary_list, sizeof(mesh_boundary_t), NULL);

  for (int i = 0 ; i < num_string_found ; i++) {

    char *str_token = strtok(string_arr[i]," \n");

    mesh_boundary_t bou;
    bou.name = strdup(str_token);
    str_token = strtok(NULL," \n");
    bou.kind = strbin2dec(str_token);

    if (mesh->dim == 2) {
      if (bou.kind == 1 || bou.kind == 2) bou.ndirpn = 1;
      if (bou.kind == 3) bou.ndirpn =  2;
      if (bou.kind == 0) bou.ndirpn =  0;
    }

    bou.nneupn   = mesh->dim - bou.ndirpn;
    bou.fnum     = malloc(mesh->dim * sizeof(int));

    for (int d = 0 ; d < mesh->dim ; d++) {
      str_token = strtok(NULL," \n");
      bou.fnum[d] = atoi(str_token);
    }
    bou.dir_loc_ixs = NULL;
    bou.dir_val = NULL;
    list_insertlast(boundary_list, &bou);
  }

  return 0;
}


int mesh_do_partition(MPI_Comm COMM, mesh_t *mesh) {

  int rank, nproc, ierr;
  MPI_Comm_size(COMM, &nproc);
  MPI_Comm_rank(COMM, &rank);

  if (nproc == 1) return 0;

  int *part = malloc(mesh->nelm_local * sizeof(int));

  if (mesh->partition == PARMETIS_GEOM || mesh->partition == PARMETIS_MESHKWAY) {

#ifdef PARTMETIS

    idx_t *elmwgt;
    idx_t wgtflag;
    idx_t numflag;
    idx_t ncon;
    idx_t ncommonnodes;
    idx_t nparts;
    real_t *tpwgts;
    real_t *ubvec;
    idx_t options[3];
    idx_t edgecut;

    elmwgt = NULL;
    wgtflag = 0;
    numflag = 0;
    nparts = nproc;
    ncon = 1;
    tpwgts = malloc(ncon*nparts*sizeof(real_t));

    for (int i = 0 ; i < ncon*nparts ; i++)
      tpwgts[i] = 1.0 / nparts;

    ncommonnodes = 3;

    options[0] = 0;
    options[1] = 0;
    options[2] = 0;

    ubvec = malloc( ncon * sizeof(real_t) );
    for (int i = 0 ; i < ncon ; i++)
      ubvec[i] = 1.05;

    if (mesh->partition == PARMETIS_GEOM) {

      ierr = ParMETIS_V3_PartGeom( elmdist, &dim, (real_t*)elmv_centroid, part, &COMM );

    }else if (mesh->partition == PARMETIS_MESHKWAY) {

      ierr = ParMETIS_V3_PartMeshKway(mesh->elm_dist, mesh->eptr, mesh->eind, elmwgt, &wgtflag, &numflag,
	  &ncon, &ncommonnodes, &nparts, tpwgts, ubvec, options, &edgecut, part, &COMM);
    }
    free(ubvec);
    free(tpwgts);

#else
    return 1;
#endif
  }

  int *eind_swi = malloc(mesh->eptr[mesh->nelm_local]*sizeof(int));
  int *npe_swi = malloc(mesh->nelm_local*sizeof(int));
  int *elm_id_swi = malloc(mesh->nelm_local*sizeof(int));
  int *eind_swi_size = malloc(nproc*sizeof(int));
  int *npe_swi_size = malloc(nproc*sizeof(int));
  int *eind_size_new = malloc(nproc*sizeof(int));
  int *npe_size_new = malloc(nproc*sizeof(int));

  ierr = swap_vectors_SCR(part, nproc, mesh->nelm_local,
      mesh->npe, mesh->eptr, mesh->eind, mesh->elm_id,
      npe_swi, eind_swi, elm_id_swi,
      npe_swi_size, eind_swi_size );

  ierr = MPI_Alltoall(npe_swi_size, 1, MPI_INT, npe_size_new, 1, MPI_INT, COMM);
  if (ierr != 0)
    return 1;

  ierr = MPI_Alltoall(eind_swi_size, 1, MPI_INT, eind_size_new, 1, MPI_INT, COMM);
  if (ierr != 0)
    return 1;

  int npe_size_new_tot = 0, eind_size_new_tot = 0;

  for (int i = 0 ; i < nproc ; i++) {
    npe_size_new_tot += npe_size_new[i];
    eind_size_new_tot += eind_size_new[i];
  }
  mesh->nelm_local = npe_size_new_tot;

  free(mesh->npe);
  free(mesh->eptr);
  free(mesh->eind);
  free(mesh->elm_id);
  free(part);

  mesh->npe = malloc(mesh->nelm_local*sizeof(int));
  mesh->eptr = malloc((mesh->nelm_local + 1)*sizeof(int));
  mesh->eind = malloc(eind_size_new_tot*sizeof(int));
  mesh->elm_id = malloc(mesh->nelm_local*sizeof(int));

  int *sdispls = malloc(nproc*sizeof(int));
  int *rdispls = malloc(nproc*sizeof(int));

  for (int i = 0 ; i < nproc ; i++) {
    sdispls[i] = 0;
    for (int j = 0 ; j < i ; j++)
      sdispls[i] += npe_swi_size[j];
  }
  for (int i = 0 ; i < nproc ; i++) {
    rdispls[i] = 0;
    for (int j = 0 ; j < i ; j++)
      rdispls[i] += npe_size_new[j];
  }

  ierr = MPI_Alltoallv(npe_swi, npe_swi_size, sdispls, MPI_INT,
      mesh->npe, npe_size_new, rdispls, MPI_INT, COMM);
  if (ierr != 0)
    return ierr;

  ierr = MPI_Alltoallv(elm_id_swi, npe_swi_size, sdispls, MPI_INT,
      mesh->elm_id, npe_size_new, rdispls, MPI_INT, COMM);
  if (ierr != 0)
    return ierr;

  mesh->eptr[0] = 0;
  for (int i = 0 ; i < mesh->nelm_local ; i++)
    mesh->eptr[i+1] = mesh->eptr[i] + mesh->npe[i];

  for (int i = 0 ; i < nproc ; i++) {
    sdispls[i] = 0;
    for (int j = 0 ; j < i ; j++)
      sdispls[i] += eind_swi_size[j];
  }
  for (int i = 0 ; i < nproc ; i++) {
    rdispls[i] = 0;
    for (int j = 0 ; j < i ; j++)
      rdispls[i] += eind_size_new[j];
  }

  ierr = MPI_Alltoallv(eind_swi, eind_swi_size, sdispls, MPI_INT,
      mesh->eind, eind_size_new, rdispls, MPI_INT, COMM);

  free(eind_swi);
  free(npe_swi);
  free(npe_swi_size);
  free(npe_size_new);
  free(eind_swi_size);
  free(eind_size_new);
  free(sdispls);
  free(rdispls);
  free(elm_id_swi);

  return 0;
}


int swap_vector(int *swap, int n, int *vector, int *new_vector, int *cuts) {

  /*
     swap       = [ 0 1 0 0 1 2 2 ]
     vector     = [ 0 1 2 3 4 5 6 ]
     new_vector = [ 0 2 3 1 4 5 6 ]
     cut        = [ 3 2 2 ]
   */

  int *aux_vector;

  if (n == 0) return 0;

  if (vector == NULL || cuts == NULL) return 1;

  if (new_vector == NULL)
    aux_vector = vector;
  else
    aux_vector = new_vector;

  int j = 0;
  for (int p = 0 ; p < n ; p++) {
    cuts[p] = 0;
    for (int i = 0 ; i < n ; i++) {
      if (swap[i] == p) {
	int aux = vector[i];
	aux_vector[i] = vector[j];
	aux_vector[j] = aux;
	j++;
	cuts[p] ++;
      }
    }
  }
  return 0;
}


int swap_vectors_SCR(int *swap, int nproc, int n,  int *npe,
    int *eptr, int *eind, int *elm_id,
    int *npe_swi, int *eind_swi, int *elm_id_swi,
    int *npe_size, int *eind_size) {

  /*
     swap        = [ 0 2 1 0 ] (swap will be generally the "part" array)
     npe         = [ 3 2 3 1 ]
     elm_id  = [ 0 0 1 2 ]
     eind        = [ 3 2 0 | 1 2 | 1 0 1 |3 ]

     npe_swi     = [ 3 1 3 2 ]
     elm_id  = [ 0 2 1 0 ]
     eind_swi    = [ 3 2 0 | 3 | 1 0 1 | 1 2 ]
   */

  int lp, pi, c;

  if (n == 0) return 0;

  if (!npe || !eind || !elm_id ||
      !eind_swi || !npe_swi || !elm_id_swi ||
      !npe_size || !eind_size) {
    return 1;
  }

  int j = pi = lp = 0;
  for (int p = 0 ; p < nproc ; p++) {

    npe_size[p] = 0;
    for (int e = 0 ; e < n ; e++) {

      if (swap[e] == p) {

	// swap npe
	npe_swi[j] = npe[e];
	elm_id_swi[j] = elm_id[e];
	j ++;

	// swap eind
	pi = eptr[e];

	for (int i = 0 ; i < npe[e] ; i++) {
	  eind_swi[lp] = eind[ pi + i ];
	  lp ++;
	}
	npe_size[p] ++;
      }
    }
  }

  c = 0;
  for (int i = 0 ; i < nproc ; i++) {
    eind_size[i] = 0;
    for (j = 0 ; j < npe_size[i] ; j++) {
      eind_size[i] += npe_swi[c];
      c++;
    }
  }

  return 0;
}


int mesh_calc_local_and_ghost(MPI_Comm COMM, mesh_t *mesh) {

  int ierr;

  int *peer_sizes, mysize, *peer_nod_glo;
  int **rep_matrix, *nrep;

  MPI_Request *request;

  int  rank, nproc;
  MPI_Comm_rank(COMM, &rank);
  MPI_Comm_size(COMM, &nproc);

  util_clean_and_sort_vector(mesh->eind, mesh->eptr[mesh->nelm_local], &mesh->local_ghost_nods, &mesh->nnods_local_ghost);

  mysize = mesh->nnods_local_ghost;
  peer_sizes = NULL;

  peer_sizes = malloc(nproc*sizeof(int));
  request = malloc(nproc*sizeof(MPI_Request));
  rep_matrix = calloc(nproc,sizeof(int*) );
  nrep = calloc(nproc,sizeof(int) );

  ierr = MPI_Allgather(&mysize, 1, MPI_INT, peer_sizes, 1, MPI_INT, COMM);

  for (int i = 0 ; i < nproc ; i++) {
    if (i != rank) {
      ierr = MPI_Isend(mesh->local_ghost_nods, mesh->nnods_local_ghost, MPI_INT, i, 0, COMM, &request[i]);
      if (ierr != 0) return 1;
    }
  }
  for (int i = 0 ; i < nproc ; i++) {
    if (i != rank) {
      peer_nod_glo = malloc(peer_sizes[i]*sizeof(int));
      ierr = MPI_Recv(peer_nod_glo, peer_sizes[i], MPI_INT, i, 0, COMM, MPI_STATUS_IGNORE );
      util_sort_vector_intersec(mesh->local_ghost_nods, mesh->nnods_local_ghost, peer_nod_glo, peer_sizes[i], &rep_matrix[i], &nrep[i]);
      free(peer_nod_glo);
    }
  }

  int nreptot = 0;
  for (int i = 0 ; i < nproc ; i++) {
    if (i != rank)
      nreptot += nrep[i];
  }
  int *rep_array = malloc(nreptot*sizeof(int));

  int c = 0;
  for (int i = 0 ; i < nproc ; i++) {
    if (i != rank) {
      for (int j = 0 ; j < nrep[i] ; j++)
	rep_array[c++] = rep_matrix[i][j];
    }
  }

  int *rep_array_clean = NULL, nreptot_clean = 0;
  ierr = util_clean_and_sort_vector(rep_array, nreptot, &rep_array_clean, &nreptot_clean);
  free(rep_array);

  int rep_count = 0, remote_rank;

  if (nreptot_clean != 0) {

    mesh->nnods_local = 0;
    mesh->nnods_ghost = 0;
    for (int i = 0 ; i < mesh->nnods_local_ghost ; i++) {

      if (rep_count < nreptot_clean) {

	if (mesh->local_ghost_nods[i] == rep_array_clean[rep_count]) {

	  int ismine = mesh_ownership_selection_rule(COMM, rep_matrix, nrep, mesh->local_ghost_nods[i], &remote_rank);
	  if (ismine != 1)
	    mesh->nnods_local++;
	  else
	    mesh->nnods_ghost++;

	  rep_count++;

	}else
	  mesh->nnods_local++;
      }
      else
	mesh->nnods_local++;
    }
  }else{
    mesh->nnods_local = mesh->nnods_local_ghost;
    mesh->nnods_ghost = 0;
  }

  mesh->local_nods = malloc(mesh->nnods_local*sizeof(int));
  mesh->ghost_nods = malloc(mesh->nnods_ghost*sizeof(int));

  int local_count = 0; int ghost_count = 0; rep_count = 0;
  if (nreptot_clean != 0) {

    for (int i = 0 ; i < mesh->nnods_local_ghost ; i++) {

      if (rep_count < nreptot_clean) {

	if (mesh->local_ghost_nods[i] == rep_array_clean[rep_count]) {

	  int ismine = mesh_ownership_selection_rule(COMM, rep_matrix, nrep, mesh->local_ghost_nods[i], &remote_rank);

	  if (ismine)
	    mesh->local_nods[local_count++] = mesh->local_ghost_nods[i];
	  else
	    mesh->ghost_nods[ghost_count++] = mesh->local_ghost_nods[i];

	  rep_count++;

	}else
	  mesh->local_nods[local_count++] = mesh->local_ghost_nods[i];

      }else
	mesh->local_nods[local_count++] = mesh->local_ghost_nods[i];
    }

  }else{
    for (int i = 0 ;i < mesh->nnods_local_ghost ; i++)
      mesh->local_nods[i] = mesh->local_ghost_nods[i];
  }


  for (int i = 0 ; i < nproc ; i++) if (i != rank) free(rep_matrix[i]);
  free(rep_matrix);
  free(nrep);
  free(request);
  free(peer_sizes);

  mesh->coord_local = malloc(mesh->nnods_local_ghost*mesh->dim*sizeof(double));
  for (int i = 0 ; i < mesh->nnods_local ; i++) {
    for (int d = 0 ; d < mesh->dim ; d++)
      mesh->coord_local[i*mesh->dim + d] = mesh->coord[mesh->local_nods[i]*mesh->dim + d];
  }
  for (int i = 0 ; i < mesh->nnods_ghost ; i++) {
    for (int d = 0 ; d < mesh->dim ; d++)
      mesh->coord_local[(i + mesh->nnods_local)*mesh->dim + d] = mesh->coord[mesh->ghost_nods[i]*mesh->dim + d];
  }

  return 0;
}


int mesh_reenumerate(MPI_Comm COMM, mesh_t *mesh) {

  int rank, nproc;
  MPI_Comm_rank(COMM, &rank);
  MPI_Comm_size(COMM, &nproc);

  int *rem_nnod = malloc(nproc*sizeof(int));
  int ierr = MPI_Allgather(&mesh->nnods_local, 1, MPI_INT, rem_nnod, 1, MPI_INT, COMM); if (ierr != 0) return 1;

  int *disp_nods = malloc(nproc*sizeof(int));
  disp_nods[0] = 0;
  for (int i = 1 ; i < nproc ; i++)
    disp_nods[i] = disp_nods[i-1] + rem_nnod[i-1];

  for (int i = 0 ; i < mesh->eptr[mesh->nelm_local] ; i++) {
    int *p = bsearch(&mesh->eind[i], mesh->local_nods, mesh->nnods_local, sizeof(int), mesh_cmpfunc);
    if (p != NULL)
      mesh->eind[i] = p - mesh->local_nods;
    else{
      p = bsearch(&mesh->eind[i], mesh->ghost_nods, mesh->nnods_ghost, sizeof(int), mesh_cmpfunc);
      if (p != NULL)
	mesh->eind[i] = mesh->nnods_local + p - mesh->ghost_nods;
      else
	return 1;
    }
  }

  mesh->local_to_global = malloc(mesh->nnods_local_ghost*sizeof(int));
  for (int i = 0 ; i < mesh->nnods_local ; i++)
    mesh->local_to_global[i] = disp_nods[rank] + i;

  MPI_Request *request = malloc(nproc*sizeof(MPI_Request));
  int *ghost_global_ix = malloc(mesh->nnods_ghost*sizeof(int));
  for (int i = 0 ; i < mesh->nnods_ghost ; i++) ghost_global_ix[i] = -1;

  for (int i = 0 ; i < nproc ; i++) {
    if (i != rank) {
      ierr = MPI_Isend(mesh->local_nods, mesh->nnods_local, MPI_INT, i, 0, COMM, &request[i]); if (ierr != 0) return 1;
    }
  }

  for (int i = 0 ; i < nproc ; i++) {

    if (i != rank) {

      int *rem_nods = malloc(rem_nnod[i]*sizeof(int));
      ierr = MPI_Recv(rem_nods, rem_nnod[i], MPI_INT, i, 0, COMM, MPI_STATUS_IGNORE ); if (ierr != 0) return 1;

      for (int j = 0 ; j < mesh->nnods_ghost ; j++) {
	int *p = bsearch(&mesh->ghost_nods[j], rem_nods, rem_nnod[i], sizeof(int), mesh_cmpfunc);
	if (p != NULL)
	  ghost_global_ix[j] = disp_nods[i] + p - rem_nods;
      }
      free(rem_nods);
    }
  }
  free(rem_nnod);
  free(request);

  for (int i = 0 ; i < mesh->nnods_ghost ; i++) if (ghost_global_ix[i] == -1) return 1;

  for (int i = 0 ; i < mesh->nnods_ghost ; i++)
    mesh->local_to_global[mesh->nnods_local + i] =  ghost_global_ix[i];

  free(ghost_global_ix);
  return 0;
}


int mesh_ownership_selection_rule(MPI_Comm COMM, int **rep_matrix, int *nrep, int node_guess, int *owner_rank) {

  int nproc, rank;
  MPI_Comm_rank(COMM, &rank);
  MPI_Comm_size(COMM, &nproc);

  int rank_guess = node_guess % nproc;

  for (int i = 0 ; i < nproc ; i++) {

    if (rank_guess == rank) {

      *owner_rank = rank_guess; return 1;

    }else{

      if (util_is_in_vector(node_guess, &rep_matrix[rank_guess][0], nrep[rank_guess]) == 1) {

	*owner_rank = rank_guess; return 0;

      }else{

	rank_guess ++; if (rank_guess == nproc) rank_guess = 0;

      }
    }
  }

  return -1;	
}


int mesh_get_bounding_box(mesh_t *mesh, double *x, double *y, double *z) {

  if (mesh->nnods_total == 0) return 0;

  x[0] = mesh->coord_local[0*mesh->dim+0]; x[1] = mesh->coord_local[0*mesh->dim+0];
  y[0] = mesh->coord_local[0*mesh->dim+1]; y[1] = mesh->coord_local[0*mesh->dim+1];
  if (mesh->dim == 2) {
    z[0] = 0.0; z[1] = 0.0;
  }
  else if (mesh->dim == 3) {
    z[0] = mesh->coord_local[0*mesh->dim+2]; z[1] = mesh->coord_local[0*mesh->dim+2];
  }

  for (int i = 1 ; i < mesh->nnods_local_ghost ; i++) {
    if (mesh->coord_local[i*mesh->dim+0] < x[0] ) x[0] = mesh->coord_local[i*mesh->dim+0];
    if (mesh->coord_local[i*mesh->dim+0] > x[1] ) x[1] = mesh->coord_local[i*mesh->dim+0];
    if (mesh->coord_local[i*mesh->dim+1] < y[0] ) y[0] = mesh->coord_local[i*mesh->dim+1];
    if (mesh->coord_local[i*mesh->dim+1] > y[1] ) y[1] = mesh->coord_local[i*mesh->dim+1];
    if (mesh->dim == 3) {
      if (mesh->coord_local[i*mesh->dim+2] < z[0]) z[0] = mesh->coord_local[i*mesh->dim+2];
      if (mesh->coord_local[i*mesh->dim+2] > z[1]) z[1] = mesh->coord_local[i*mesh->dim+2];
    }
  }

  return 0;
}


int mesh_get_domain_center(MPI_Comm COMM, mesh_t *mesh, double center[3]) {

  double x[2] ,y[2] ,z[2] ,x_abs[2] ,y_abs[2] ,z_abs[2];

  int rank, nproc;
  MPI_Comm_size(COMM, &nproc);
  MPI_Comm_rank(COMM, &rank);

  int *x_all = malloc(nproc*2*sizeof(double));
  int *y_all = malloc(nproc*2*sizeof(double));
  int *z_all = malloc(nproc*2*sizeof(double));

  int ierr = mesh_get_bounding_box(mesh, x, y, z);

  ierr = MPI_Allgather(x, 2, MPI_DOUBLE, x_all, 2, MPI_DOUBLE, COMM);
  ierr = MPI_Allgather(y, 2, MPI_DOUBLE, y_all, 2, MPI_DOUBLE, COMM);
  ierr = MPI_Allgather(z, 2, MPI_DOUBLE, z_all, 2, MPI_DOUBLE, COMM);

  x_abs[0] = x_all[0]; x_abs[1] = x_all[1];
  y_abs[0] = y_all[0]; y_abs[1] = y_all[1];
  z_abs[0] = z_all[0]; z_abs[1] = z_all[1];

  for (int i = 1 ; i < nproc ; i++) {
    if (x_all[2*i+0] < x_abs[0]) x_abs[0] = x_all[2*i+0];
    if (x_all[2*i+1] > x_abs[1]) x_abs[1] = x_all[2*i+1];
    if (y_all[2*i+0] < y_abs[0]) y_abs[0] = y_all[2*i+0];
    if (y_all[2*i+1] > y_abs[1]) y_abs[1] = y_all[2*i+1];
    if (z_all[2*i+0] < z_abs[0]) z_abs[0] = z_all[2*i+0];
    if (z_all[2*i+1] > z_abs[1]) z_abs[1] = z_all[2*i+1];
  }

  center[0] = (x_abs[1] + x_abs[0])/2;
  center[1] = (y_abs[1] + y_abs[0])/2;
  center[2] = (z_abs[1] + z_abs[0])/2;

  return ierr;
}


int get_bbox_limit_lengths(MPI_Comm COMM, mesh_t *mesh, double *lx, double *ly, double *lz) {

  int rank, nproc;
  double x[2] ,y[2] ,z[2] ,x_abs[2] ,y_abs[2] ,z_abs[2];

  MPI_Comm_size(COMM, &nproc);
  MPI_Comm_rank(COMM, &rank);

  int *x_all = malloc(nproc*2*sizeof(double));
  int *y_all = malloc(nproc*2*sizeof(double));
  int *z_all = malloc(nproc*2*sizeof(double));

  int ierr = mesh_get_bounding_box(mesh, x, y, z);

  ierr = MPI_Allgather(x, 2, MPI_DOUBLE, x_all, 2, MPI_DOUBLE, COMM);
  ierr = MPI_Allgather(y, 2, MPI_DOUBLE, y_all, 2, MPI_DOUBLE, COMM);
  ierr = MPI_Allgather(z, 2, MPI_DOUBLE, z_all, 2, MPI_DOUBLE, COMM);

  x_abs[0]=x_all[0]; x_abs[1]=x_all[1];
  y_abs[0]=y_all[0]; y_abs[1]=y_all[1];
  z_abs[0]=z_all[0]; z_abs[1]=z_all[1];

  for (int i = 1 ; i < nproc ; i++) {
    if ( x_all[2*i+0] < x_abs[0] ) x_abs[0] = x_all[2*i+0];
    if ( x_all[2*i+1] > x_abs[1] ) x_abs[1] = x_all[2*i+1];
    if ( y_all[2*i+0] < y_abs[0] ) y_abs[0] = y_all[2*i+0];
    if ( y_all[2*i+1] > y_abs[1] ) y_abs[1] = y_all[2*i+1];
    if ( z_all[2*i+0] < z_abs[0] ) z_abs[0] = z_all[2*i+0];
    if ( z_all[2*i+1] > z_abs[1] ) z_abs[1] = z_all[2*i+1];
  }

  *lx = x_abs[1] - x_abs[0];
  *ly = y_abs[1] - y_abs[0];
  *lz = z_abs[1] - z_abs[0];

  free(x_all);
  free(y_all);
  free(z_all);

  return ierr;
}


int build_structured_2d(int **eind, int **eptr, double **coor, double limit[4], int nx, int ny) {

  double   x0 = limit[0];
  double   x1 = limit[1];
  double   y0 = limit[2];
  double   y1 = limit[3];
  int      nnod = nx*ny;
  int      nex = (nx-1);
  int      ney = (ny-1);
  double   dx = (x1-x0)/nex;
  double   dy = (y1-y0)/ney;
  int      nelm = nex*ney;
  int      npe = 4;
  int      i, j, e, n;

  *eind = malloc(nelm*npe*sizeof(int));
  *eptr = malloc((nelm+1)*sizeof(int));
  *coor = malloc(nnod*2*sizeof(double));

  for (i=0;i<nex;i++) {
    for (j=0;j<ney;j++) {
      e = i*nex + j;
      (*eptr)[e+0]     = ( e + 0 )*npe;
      (*eptr)[e+1]     = ( e + 1 )*npe;
      (*eind)[e*npe+0] = j   + i*nx;
      (*eind)[e*npe+1] = j+1 + i*nx;
      (*eind)[e*npe+2] = j+1 + (i+1)*nx;
      (*eind)[e*npe+3] = j   + (i+1)*nx;
    }
  }
  for (i=0;i<nx;i++) {
    for (j=0;j<ny;j++) {
      n = i*nx + j;
      (*coor)[n*2+0] = x0 + dx*j;
      (*coor)[n*2+1] = y0 + dy*i;
    }
  }

  return 0;
}


int get_element_structured_2d(double centroid[2], double limit[4], int nx, int ny, int *es) {

  double   x0 = limit[0];
  double   x1 = limit[1];
  double   y0 = limit[2];
  double   y1 = limit[3];
  int      nex = (nx-1);
  int      ney = (ny-1);
  int      i, j;
  double   dx = (x1-x0)/nex;
  double   dy = (y1-y0)/ney;
  double   x_min;
  double   x_max;
  double   y_min;
  double   y_max;

  j=0;
  while (j<nex) {
    x_min = x0 + dx*j;
    x_max = x0 + dx*(j+1);
    if (x_min < centroid[0] && centroid[0] < x_max) {
      break;
    }
    j++;
  }

  i=0;
  while (i<ney) {
    y_min = y0 + dy*i;
    y_max = y0 + dy*(i+1);
    if (y_min < centroid[1] && centroid[1] < y_max) {
      break;
    }
    i++;
  }

  *es = i*nex + j;
  return 0;
}


int mesh_cmpfunc(const void * a, const void * b) {

  return *(int*)a - *(int*)b;
}


int cmpfunc_for_list(void * a, void * b) {

  return *(int*)a - *(int*)b;
}
