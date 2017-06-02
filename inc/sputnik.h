/*****************************************************************************************************
   SPUTNIK external lybraries
*****************************************************************************************************/

#include "petscksp.h"
#include "list.h"
#include "mpi.h"
#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#define MACRO         1
#define MICRO         2
#define MESH_N_LENGTH 64
#define BUF_N_LENGTH  128

// spu_parser.c
int parse_mpi(const char mpi_file[], int nproc[2], int * nkind, int ** nproc_k);

// spu_mesh.c
int read_mesh(char *mesh_n, char *mesh_f, int rank, int nproc, int ** elmdist, int ** eptr, int ** eind);
int read_mesh_CSR_GMSH(char *mesh_n, int rank, int nproc, int ** elmdist, int ** eptr, int ** eind);
