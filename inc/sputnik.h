#ifndef SPUTNIK_H
#define SPUTNIK_H

#include "petscksp.h"
#include "slepceps.h"
#include "list.h"
#include "stdbool.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fem.h"
#include "gmsh.h"
#include "material.h"
#include <gsl/gsl_linalg.h>

#define FORMAT_NULL        0
#define FORMAT_GMSH        1
#define FORMAT_ALYA        2

#define PRINT_PETSC        0
#define PRINT_VTK          1
#define PRINT_VTU          2
#define PRINT_VTKPART      4
#define PRINT_ALL          8

int         partition_algorithm;

char        *myname;

int         rank_wor;               //  rank on world comm
int         nproc_wor;              //  # of processes (WORLD_COMM)
int         dim;                    //  problem dimensions 

int         flag_print;
bool        flag_reactions;
bool        flag_first_alloc;

/* Cilindrical Fiber in quad Matrix */
double      LX, LY, LZ;
int         flag_fiber_cilin;
int         nx_fibers;
int         ny_fibers;
double      fiber_cilin_r;
double      fiber_cilin_center_devi[3];
double      fiber_cilin_vals[4];
double      center_domain[3];

/* Interpolation over structured grids */
int         nx_interp;
int         ny_interp;
int         nz_interp;

// Structures to save de mesh on CSR format 

int          *elm_id;                 // element property number
double       *elmv_centroid;

int          *StartIndexRank;
int          *allnods;                // all nodes including locals and ghosts
int          nallnods;                // <nmynods> + <nghost>
int          *mynods;                 // Original (gmsh) numbers of my nodes
int          nmynods;                 // number of local nodes
int          *ghost;                  // ghosts nodes (global numbering)
int          nghost;                  // number of ghost nodes
int          ntotnod;                 // number of total nodes in the mesh

double       *coord;                  // nodes' coordinates

int          *loc2petsc;              // array of size <nmynods> + <nghost>
                                      // returns the position in PETSc matrix & vectors

int read_mesh_elmv( MPI_Comm COMM, char *myname, char *mesh_n, int mesh_f);
int read_mesh_elmv_CSR_GMSH( MPI_Comm COMM, char *myname, char *mesh_n);
int read_mesh_elmv_CSR_ALYA( MPI_Comm COMM, char *myname, char *mesh_n);

int read_mesh_coord( MPI_Comm COMM, char *mesh_n, int mesh_f);
int read_mesh_coord_GMSH( MPI_Comm COMM, char *mesh_n);
int read_mesh_coord_ALYA( MPI_Comm COMM, char *mesh_n);

int CSR_give_pointer( int e, int *npe, int *eind, int *p);
int clean_vector_qsort(int n, int *input, int **output, int *not_rep);
int give_repvector_qsort( MPI_Comm * comm, char *myname, int n, int *input, int **output, int *nrep);
int give_inter_sort( MPI_Comm comm, char *myname,
    int *array1, int n1,
    int *array2, int n2,
    int **reps, int *nreps);
int is_in_vector(int val, int *vector, int size);
int cmpfunc(const void * a, const void * b);
int gmsh_npe(int code);
int gmsh_is_vol_elm(int code);
int get_bbox_limit_lengths(MPI_Comm PROBLEM_COMM, double *coord, int n, double *lx, double *ly, double *lz);
int get_bbox_local_limits(double *coord, int n, double *x, double *y, double *z);
int get_domain_center(MPI_Comm PROBLEM_COMM, double *coord, int n, double center[3]);
int interpolate_structured_2d(double limit[2], int nx, int ny, double *field, double *var_interp);
int get_element_structured_2d(double centroid[2], double limit[4], int nx, int ny, int *es);
int build_structured_2d(int **eind, int **eptr, double **coor, double limit[4], int nx, int ny);

#endif
