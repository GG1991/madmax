#ifndef MESH_H
#define MESH_H

#define PARMETIS_GEOMKWAY   1
#define PARMETIS_GEOM       2
#define PARMETIS_KWAY       3
#define PARMETIS_MESHKWAY   4

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


int part_mesh( MPI_Comm COMM, char *myname, double *centroid );
int reenumerate_PETSc( MPI_Comm COMM );
int calculate_ghosts( MPI_Comm COMM, char *myname );
int ownership_selec_rule( MPI_Comm COMM, int **repeated, int *nrep, int node, int *remoterank );
int swap_vectors_SCR( int *swap, int nproc, int n,  int *npe, 
    int *eptr, int *eind, int *elm_id,
    int *npe_swi, int *eind_swi, int *elm_id_swi,
    int *cuts_npe, int *cuts_eind );

#endif
