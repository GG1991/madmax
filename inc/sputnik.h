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

#define MACRO              1     // MACRO IDs and colors
#define MICRO              2     // MICRO IDs and colors
#define NBUF               256   // Buffers length for read from a file


// Execution schemes
#define SCHEME_MACRO_MICRO_COU     1   // Coupled execution
#define SCHEME_MACRO_PART_MESH     2   // Standalone execution
#define SCHEME_MICRO_PART_MESH     3
#define SCHEME_MACRO_PDE_ALONE     4
#define SCHEME_MICRO_PDE_ALONE     5

#define PARMETIS_GEOMKWAY  1
#define PARMETIS_GEOM      2
#define PARMETIS_KWAY      3
#define PARMETIS_MESHKWAY  4


/*****************************************************************************************************
   SPUTNIK global variables
*****************************************************************************************************/

int *id_vec;    // ID vector size = #proc (info of which ID has each rank
int nproc_mac;  // number of macro processes 
int nproc_mic;  // number of micro processes
int nstru_mic;  // number of micro structures
int scheme;     // communication approach


/*****************************************************************************************************
   SPUTNIK function definitions
*****************************************************************************************************/

// spu_parser.c (common rutines for parser inputs files from MACRO and MICRO)
int spu_parse_comm( char *input_fila);

// spu_mesh.c
int read_mesh(MPI_Comm * comm, char *mesh_n, char *mesh_f, int ** elmdist, int ** eptr, int ** eind);
int read_mesh_CSR_GMSH(MPI_Comm * comm, char *mesh_n, int ** elmdist, int ** eptr, int ** eind);
int part_mesh_PARMETIS(MPI_Comm * comm, int * elmdist, int * eptr, int * eind, double * centroid, int algorithm);
