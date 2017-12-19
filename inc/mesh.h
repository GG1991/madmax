#ifndef MESH_H
#define MESH_H

#include "list.h"
#include "util.h"

#define PARMETIS_GEOMKWAY   1
#define PARMETIS_GEOM       2
#define PARMETIS_KWAY       3
#define PARMETIS_MESHKWAY   4

#define MAX_ADJ_NODES 30
#define MAX_NUM_OF_BOUNDARIES 4

int *elmdist;
int nelm;
int *eptr;
int *eind;

int *allnods;
int nallnods;
int *mynods;
int nmynods;
int *ghost;
int nghost;
int ntotnod;

double *coord;

int *loc2petsc;

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
  int nelem;
  int nelm_local;
  int *eptr;
  int *eind;
  int **elements;
  int *local_nodes;
  int *ghost_nodes;
  int *local_ghost_nodes;
  int *local_to_global;
  double *coord;
  list_t boundary_list;

}mesh_t;

extern mesh_t mesh;


typedef struct{

  int dim;
  int nx, ny, nz;
  int nex, ney, nez;
  int nelm;
  int nnods;
  int *eptr;
  int *eind;
  int **elements;
  double lx, ly, lz;
  double hx, hy, hz;
  double *coord;

}mesh_struct_t;

extern mesh_struct_t mesh_struct;


int part_mesh(MPI_Comm COMM, char *myname, double *centroid);
int reenumerate_PETSc(MPI_Comm COMM);

int calc_local_and_ghost( MPI_Comm COMM, int nallnods, int *allnods,
    int *ntotnod, int *nmynods, int **mynods, int *nghost , int **ghost );

int ownership_selec_rule( MPI_Comm COMM, int **repeated, int *nrep, int node, int *remoterank );
int swap_vectors_SCR( int *swap, int nproc, int n,  int *npe,
    int *eptr, int *eind, int *elm_id,
    int *npe_swi, int *eind_swi, int *elm_id_swi,
    int *cuts_npe, int *cuts_eind );

int mesh_fill_boundary_list_from_command_line(command_line_t *command_line, list_t *boundary_list);


#endif
