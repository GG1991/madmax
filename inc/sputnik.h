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
#define MACRO_MICRO         1   // Coupled execution
#define MACRO_ALONE         2   // Standalone execution of MACRO
#define MICRO_ALONE         3   // Standalone execution of MICRO

#define PARMETIS_GEOMKWAY   1
#define PARMETIS_GEOM       2
#define PARMETIS_KWAY       3
#define PARMETIS_MESHKWAY   4


/*****************************************************************************************************
   SPUTNIK global variables
*****************************************************************************************************/

char         input_n[NBUF];            // Input file name

MPI_Comm     world_comm;

int          rank_wor;                 //  rank on world comm
int          nproc_wor;                //  # of processes (world_comm)

int          *id_vec;                  // ID vector size = #proc (info of which ID has each rank
int          nproc_mac_tot;            // number of macro processes total (inside world_comm)
int          nproc_mic_tot;            // number of micro processes total (inside world_comm)
int          nstruc_mic;               // number of micro structures
int          *nproc_per_mic;           // number of processes per micro structure ( size = nstruc_mic )
int          nproc_mic_group;          // number of micro process in a group = sum_i nproc_per_mic[i]
int          nmic_worlds;              // number of micro worlds nproc_mic / nproc_mic_group
int          scheme;                   // communication approach

// Time measurement

double       t0,t1;

// Structures to save de mesh on CSR format 

char         mesh_n[NBUF];             // Mesh file name
char         mesh_f[4];                // Mesh format name

int          *part;
int          *elmdist;                 // number of elements inside each procesor
int           nelm;                    // # of local elements
int          *eptr;                    // list of indeces of nodes inside eind
int          *eind;                    // list of nodes for elem "i" is between 
                                       // eind[eptr[i]] eind[eptr[i+1]] (not including)


/*****************************************************************************************************
   SPUTNIK function definitions
*****************************************************************************************************/

// spu_parser.c (common rutines for parser inputs files from MACRO and MICRO)
int spu_parse_scheme( char *input );
int spu_parse_mesh( char * input );

// spu_mesh.c
int read_mesh(MPI_Comm * comm, char *myname, char *mesh_n, char *mesh_f, int ** elmdist, int ** eptr, int ** eind);
int read_mesh_CSR_GMSH(MPI_Comm *comm, char *myname, char *mesh_n, int **elmdist, int **eptr, int **eind);
int part_mesh_PARMETIS(MPI_Comm *comm, char *myname, int *elmdist, int *eptr, int *eind, int *part, double *centroid, int algorithm);
int swap_vectors_SCR( int *swap, int nproc, int n,  int *npe, int *eptr, int *eind, int *npe_new, int *eind_new, int
*cuts_npe, int *cuts_eind);
int CSR_give_pointer( int e, int *npe, int *eind, int *p);
