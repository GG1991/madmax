/*****************************************************************************************************
   SPUTNIK external lybraries
*****************************************************************************************************/

#include "petscksp.h"
#include "petscsys.h"
#include "list.h"
#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fem.h"
#include "material.h"
#include "boundary.h"

#ifndef SPUTNIK_H
#define SPUTNIK_H

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

#define CHECK_SPU_ERROR(data)                                                          \
     {if(!data){                                                                       \
	 printf("SPUTNIK ERROR on %s line %d\n",__FILE__,__LINE__);                    \
	 return -1;                                                                    \
     }}

#define CHECK_NEG_ERROR(data_int)                                                      \
     {if(data_int<0){                                                                  \
	 printf("SPUTNIK ERROR on %s line %d\n",__FILE__,__LINE__);                    \
	 return -1;                                                                    \
     }}

/*****************************************************************************************************
   SPUTNIK user defined structures
*****************************************************************************************************/

typedef struct _physical_t{
    
    int    dim;
    int    GmshID;
    char   *name; 
    int    FlagFound;

}physical_t;


/*****************************************************************************************************
   SPUTNIK global variables
*****************************************************************************************************/

char         input_n[NBUF];            // Input file name

MPI_Comm     world_comm;
MPI_Status   status;

bool         print_flag;

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

FILE         *time_fl;
FILE         *FileOutputStructures;
double       t0,t1;
double       *time_vec;

// Structures to save de mesh on CSR format 

char          mesh_n[NBUF];            // Mesh file name
char          mesh_f[4];               // Mesh format name

int           *part;
int           *elmdist;                // number of elements inside each procesor
int            nelm;                   // # of local elements
int           *eptr;                   // list of indeces of nodes inside eind
int           *eind;                   // list of nodes for elem "i" is between 
                                       // eind[eptr[i]] eind[eptr[i+1]] (not including)
int           *PhysicalID;             // element property number

int           *StartIndexRank;
int           *AllMyNodOrig;           // Original (gmsh) numbers of my nodes + my ghosts
int           NAllMyNod;               // <NMyNod> + <NMyGhost>
int           *MyNodOrig;              // Original (gmsh) numbers of my nodes
int           NMyNod;                  // Number of my nodes 
int           *MyGhostOrig;            // Original (gmsh) numbers of my ghosts nodes
int           NMyGhost;                // Number of my ghost nodes
int           NTotalNod;               // Number of total nodes in the mesh

double        *coord;                  // nodes' coordinates

int           *loc2petsc;              // array of size <NMyNod>+<NMyGhost>
                                       // returns the position in PETSc matrix & vectors

// List of different utilities
list_t function_list;
list_t physical_list;

/*****************************************************************************************************
   SPUTNIK function definitions
*****************************************************************************************************/

// spu_parser.c (common rutines for parser inputs files from MACRO and MICRO)
int spu_parse_scheme( char *input );
int spu_parse_mesh( char * input );
int SpuParseMaterials(MPI_Comm *PROBLEM_COMM, char * input );
int SpuParsePhysicalEntities( MPI_Comm *PROBLEM_COMM, char *mesh_n );
int SetGmshIDOnMaterialsAndBoundaries(MPI_Comm PROBLEM_COMM);
int CheckPhysicalID(void);
int SpuParseBoundary(MPI_Comm *PROBLEM_COMM, char *input);
int SpuParseFunctions(MPI_Comm *PROBLEM_COMM, char *input );
int StrBin2Dec(char *str);

// spu_mesh.c
int read_mesh_elmv(MPI_Comm * comm, char *myname, char *mesh_n, char *mesh_f);
int read_mesh_elmv_CSR_GMSH(MPI_Comm *comm, char *myname, char *mesh_n);
int read_mesh_coord(MPI_Comm * comm, char *myname, char *mesh_n, char *mesh_f);
int read_mesh_coord_GMSH(MPI_Comm * comm, char *myname, char *mesh_n);
int part_mesh_PARMETIS(MPI_Comm *comm, FILE *time_fl, char *myname, double *centroid, int algorithm);
int swap_vectors_SCR( int *swap, int nproc, int n,  int *npe, 
    int *eptr, int *eind, int *PhysicalID,
    int *npe_swi, int *eind_swi, int *PhysicalID_swi,
    int *cuts_npe, int *cuts_eind );
int CSR_give_pointer( int e, int *npe, int *eind, int *p);
int clean_vector_qsort(MPI_Comm * comm, char *myname, int n, int *input, int **output, int *not_rep);
int give_repvector_qsort(MPI_Comm * comm, char *myname, int n, int *input, int **output, int *nrep);
int give_inter_sort(MPI_Comm *comm, char *myname, int *array1, int n1, int *array2, int n2, int **reps, int *nreps);
int calculate_ghosts(MPI_Comm * comm, char *myname);
int ownership_selec_rule( MPI_Comm *comm, int **repeated, int *nrep, int node, int *remoterank );
int is_in_vector(int val, int *vector, int size);
int reenumerate_PETSc(MPI_Comm *comm);
int search_position_linear(int *array, int size, int val, int *pos);
int search_position_logn(int *array, int size, int val, int *pos);
int SpuReadBoundary(MPI_Comm PROBLEM_COMM, char *mesh_n, char *mesh_f, FILE *outfile );
int SpuReadBoundaryGmsh(MPI_Comm PROBLEM_COMM, char *mesh_n, FILE *outfile);
int cmpfunc (const void * a, const void * b);
int cmpfunc_for_list (void * a, void * b);
int cmpfuncBou (void * a, void * b);
int GmshNodesPerElement(int code);
int GmshIsAsurfaceElement(int code);

// spu_boundary.c
int SputnikSetDisplacementOnBoundary( double time, Vec *x );

// spu_vtk.c
int spu_vtk_partition( char *vtkfile_n, MPI_Comm *comm );
int vtkcode(int dim,int npe);

// spu_time.c
int save_time(MPI_Comm *comm, const char *string, FILE *file, double dt);

// spu_assembly.c
int GetPETScIndeces(int *LocalNod, int n, int *local2PETSc, int *PETScIndex);
int GetElemCoord(int *LocalNod, int n, double ElemCoord[8][3]);
int AssemblyJacobianSmallDeformation(Mat *J);
int AssemblyResidualSmallDeformation(Vec *Displacement_old, Vec *Residue);
int GetShapeDerivs(int gp, int npe, double coor[8][3], double ShapeDerivs[8][3], double *DetJac);
int GetB( int npe, double ShapeDerivs[8][3], double B[6][3*8] );
int GetWeight(int npe, double **wp);
int GetDsDe( int npe, double *ElemDisp, double DsDe[6][6] );
material_t * GetMaterial(int GmshIDToSearch);
int GetElemenDispls( int e, Vec *Displacement, double *ElemDisps );

#endif
