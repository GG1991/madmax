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

int      *elmdist;     // number of elements in each procesor
int       nelm;        // # of local elements
int      *eptr;        // list of indeces of nodes inside eind
int      *eind;        // list of nodes for elem "i" [eind[eptr[i]] , eind[eptr[i+1]]) (not including)

int      *allnods;     // all nodes including locals and ghosts
int       nallnods;    // <nmynods> + <nghost>
int      *mynods;      // Original (gmsh) numbers of my nodes
int       nmynods;     // number of local nodes
int      *ghost;       // ghosts nodes (global numbering)
int       nghost;      // number of ghost nodes
int       ntotnod;     // number of total nodes in the mesh

double   *coord;       // nodes' coordinates

int      *loc2petsc;   // array of size <nmynods> + <nghost>
                       // returns the position in PETSc matrix & vectors

typedef struct{

  char     *name;
  int       kind;                // boundary kind ( xxx -> to decimal )
  int      *fnum;                // funtion numbers to evaluate
  int       ndirpn;              // dirichlet values per node
  int       nneupn;              // neumann values per node
  int       ndirix;              // number of dir indeces
  int       ndir;                // number of dir nodes
  int      *dir_loc_ixs;         // dirichlet indeces (local)
  int      *dir_glo_ixs;         // dirichlet indeces (global)
  double   *dir_val;             // dirichlet values

}mesh_boundary_t;

typedef struct{

  int dim;
  int n_elem_local;
  int n_elem_total;
  int n_node_local;
  int n_node_total;
  int *local_nodes;
  int *ghost_nodes;
  int *local_ghost_nodes;
  double *coord;
  list_t boundary_list;

}mesh_t;

mesh_t mesh;


int part_mesh( MPI_Comm COMM, char *myname, double *centroid );
int reenumerate_PETSc( MPI_Comm COMM );

int calc_local_and_ghost( MPI_Comm COMM, int nallnods, int *allnods,
    int *ntotnod, int *nmynods, int **mynods, int *nghost , int **ghost );

int ownership_selec_rule( MPI_Comm COMM, int **repeated, int *nrep, int node, int *remoterank );
int swap_vectors_SCR( int *swap, int nproc, int n,  int *npe, 
    int *eptr, int *eind, int *elm_id,
    int *npe_swi, int *eind_swi, int *elm_id_swi,
    int *cuts_npe, int *cuts_eind );
;
int mesh_fill_boundary_list_from_command_line(command_line_t *command_line, list_t *boundary_list);

#endif
