#ifndef MESH_H
#define MESH_H

#include "list.h"
#include "util.h"
#include "myio.h"
#include "geometry.h"

#ifdef PARMETIS
#include "parmetis.h"
#endif

#define PARMETIS_GEOMKWAY   1
#define PARMETIS_GEOM       2
#define PARMETIS_KWAY       3
#define PARMETIS_MESHKWAY   4

#define NBUF 256

#define MAX_ADJ_NODES 30
#define MAX_NUM_OF_BOUNDARIES 4
#define MAX_NPE 8
#define MAX_DIM 3

typedef struct{

  char *name;
  int kind;
  int *fnum;
  int ndirpn;
  int nneupn;
  int ndirix;
  int ndir;
  int *dir_loc_ixs;
  int *dir_glo_ixs;
  double *dir_val;

}mesh_boundary_t;

typedef struct{

  int dim;
  int nnods_total;
  int nnods_local;
  int nnods_ghost;
  int nnods_local_ghost;
  int *local_nods;
  int *ghost_nods;
  int *local_ghost_nods;
  int *local_to_global;
  int nelm_total;
  int nelm_local;
  int *eptr;
  int *eind;
  int **elements;
  int *npe;
  int *nelm_dist;
  int *elm_id;
  double *elm_centroid;
  int partition;
  double *coord;
  double *coord_local;

}mesh_t;

typedef struct{

  int dim;
  int nx, ny, nz;
  int nex, ney, nez;
  int nelm;
  int nn;
  int npe;
  int nnods_boundary;
  int *boundary_nods;
  int *boundary_indeces;
  double *boundary_coord;
  double lx, ly, lz;
  double hx, hy, hz;
  double vol;
  double vol_elm;

}mesh_struct_t;

extern mesh_t mesh;
extern mesh_struct_t mesh_struct;
extern list_t boundary_list;


int mesh_do_partition(MPI_Comm COMM, mesh_t *mesh);
int mesh_reenumerate(MPI_Comm COMM, mesh_t *mesh);
int mesh_calc_local_and_ghost(MPI_Comm COMM, mesh_t *mesh);
int mesh_ownership_selection_rule(MPI_Comm COMM, int **rep_matrix, int *nrep, int node_guess, int *owner_rank);
int mesh_cmpfunc(const void * a, const void * b);
int mesh_check_orientation(mesh_t *mesh);
int mesh_check_orientation_tria(int *nods, double *coords);
int mesh_check_orientation_quad(int *nods, double *coords);

int mesh_struct_init(int dim, int *sizes, double *length, mesh_struct_t *mesh_struct);
int mesh_struct_get_node_coord(mesh_struct_t *mesh_struct, int node, double *coord);
int mesh_struct_get_elem_nods(mesh_struct_t *mesh_struct, int e, int *elem_nods);
int mesh_struct_get_elem_indeces(mesh_struct_t *mesh_struct, int e, int *elem_indeces);

int swap_vectors_SCR(int *swap, int nproc, int n,  int *npe,
    int *eptr, int *eind, int *elm_id,
    int *npe_swi, int *eind_swi, int *elm_id_swi,
    int *cuts_npe, int *cuts_eind);

int mesh_fill_boundary_list_from_command_line(command_line_t *command_line, list_t *boundary_list, mesh_t *mesh);


#endif
