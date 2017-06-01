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
#include "parmetis.h"

#define MACRO 1
#define MICRO 2

// spu_parser.c
int parse_mpi(const char mpi_file[], int nproc[2], int * nkind, int ** nproc_k);

// spu_mesh.c
int read_mesh(char *mesh_n, char *mesh_f);
int read_mesh_GMSH_CSR(char *mesh_n, int *eptr, int *eind);
