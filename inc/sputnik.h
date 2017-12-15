#ifndef SPUTNIK_H
#define SPUTNIK_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include "petscksp.h"
#include "slepceps.h"
#include "list.h"
#include "stdbool.h"
#include "fem.h"
#include "gmsh.h"
#include "material.h"

#define FORMAT_NULL        0
#define FORMAT_GMSH        1
#define FORMAT_ALYA        2

int         partition_algorithm;

char        *myname;

int         rank_wor;
int         nproc_wor;
int         dim;

int          *elm_id;
double       *elmv_centroid;

int read_mesh_elmv( MPI_Comm COMM, char *myname, char *mesh_n, int mesh_f);
int read_mesh_elmv_CSR_GMSH( MPI_Comm COMM, char *myname, char *mesh_n);

int clean_vector_qsort(int n, int *input, int **output, int *not_rep);
int give_repvector_qsort(MPI_Comm *COMM, char *myname, int n, int *input, int **output, int *nrep);

int give_inter_sort( MPI_Comm COMM, int *array1, int n1, int *array2, int n2, int **reps, int *nreps);

int is_in_vector(int val, int *vector, int size);
int cmpfunc(const void * a, const void * b);
int gmsh_npe(int code);
int get_bbox_limit_lengths(MPI_Comm PROBLEM_COMM, double *coord, int n, double *lx, double *ly, double *lz);
int get_bbox_local_limits(double *coord, int n, double *x, double *y, double *z);
int get_domain_center(MPI_Comm PROBLEM_COMM, double *coord, int n, double center[3]);
int interpolate_structured_2d(double limit[2], int nx, int ny, double *field, double *var_interp);
int get_element_structured_2d(double centroid[2], double limit[4], int nx, int ny, int *es);
int build_structured_2d(int **eind, int **eptr, double **coor, double limit[4], int nx, int ny);

#endif
