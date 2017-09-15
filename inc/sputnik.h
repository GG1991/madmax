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

#define NBUF               256   // Buffers length for read from a file

// Execution schemes
#define MACRO_MICRO         1   // Coupled execution
#define MACRO_ALONE         2   // Standalone execution of MACRO
#define MICRO_ALONE         3   // Standalone execution of MICRO

#define PARMETIS_GEOMKWAY   1
#define PARMETIS_GEOM       2
#define PARMETIS_KWAY       3
#define PARMETIS_MESHKWAY   4

int partition_algorithm;

#if defined(PETSC_USE_LOG)
  PetscLogStage stages[3];
  PetscLogEvent EVENT_READ_MESH_ELEM,
		EVENT_PART_MESH,
		EVENT_CALC_GHOSTS,
		EVENT_REENUMERATE,
		EVENT_READ_COORD,
		EVENT_INIT_GAUSS,
		EVENT_ALLOC_MATVEC,
		EVENT_SET_DISP_BOU,
		EVENT_ASSEMBLY_JAC,
		EVENT_ASSEMBLY_RES,
		EVENT_SOLVE_SYSTEM;
#endif

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

MPI_Comm     WORLD_COMM;
MPI_Status   status;
PetscViewer  viewer;
int          reason;
char         reason_s[64];

int rank_wor;  //  rank on world comm
int nproc_wor; //  # of processes (WORLD_COMM)
int dim;       //  problem dimensions 
int nvoi;      //  number of voigt components (3 if dim=2, 6 if dim=3)

#define PRINT_PETSC        0
#define PRINT_VTK          1
#define PRINT_VTU          2
#define PRINT_VTKPART      4
#define PRINT_ALL          8

int      flag_print;
bool     flag_reactions;

/*
   Cilindrical Fiber in quad Matrix
*/
double   LX, LY, LZ;
int      flag_fiber_cilin;
int      nx_fibers;
int      ny_fibers;
double   fiber_cilin_r;
double   fiber_cilin_center_devi[3];
double   fiber_cilin_vals[4];
double   center_domain[3];

#define FORMAT_NULL        0
#define FORMAT_GMSH        1
#define FORMAT_ALYA        2

int       mesh_f;                   // Mesh format number

// Time measurement

FILE      *time_fl;
FILE      *file_out;
double    t0,t1;
double    *time_vec;

// Solver options

int       nr_max_its;
double    nr_norm_tol;

// Structures to save de mesh on CSR format 

char      mesh_n[NBUF];           // Mesh file name

int      *part;
int      *elmdist;                // number of elements inside each procesor
int      nelm;                    // # of local elements
int      *eptr;                   // list of indeces of nodes inside eind
int      *eind;                   // list of nodes for elem "i" is between 
                                  // eind[eptr[i]] eind[eptr[i+1]] (not including)
int      *PhysicalID;             // element property number
double   *elmv_centroid;

int *StartIndexRank;
int *allnods;                // all nodes including mynods and ghost
int nallnods;                // <nmynods> + <nghost>
int *mynods;                 // Original (gmsh) numbers of my nodes
int nmynods;                 // Number of <mynods> 
int *ghost;                  // Original numbers of my ghosts nodes
int nghost;                  // Number of my ghost nodes
int NTotalNod;               // Number of total nodes in the mesh

double *coord;               // nodes' coordinates

int *loc2petsc;              // array of size <nmynods>+<nghost>
                             // returns the position in PETSc matrix & vectors

// List of different utilities
list_t function_list;
list_t physical_list;

/*****************************************************************************************************
   SPUTNIK function definitions
*****************************************************************************************************/

// spu_parser.c (common rutines for parser inputs files from MACRO and MICRO)
int spu_parse_scheme( char *input );
int parse_material(MPI_Comm PROBLEM_COMM, char *input);
int CheckPhysicalID(void);
int parse_function(MPI_Comm PROBLEM_COMM, char *input);
int StrBin2Dec(char *str);

// spu_mesh.c
int read_mesh_elmv(MPI_Comm PROBLEM_COMM, char *myname, char *mesh_n, int mesh_f);
int read_mesh_elmv_CSR_GMSH(MPI_Comm PROBLEM_COMM, char *myname, char *mesh_n);
int read_mesh_elmv_CSR_ALYA(MPI_Comm PROBLEM_COMM, char *myname, char *mesh_n);

int read_mesh_coord(MPI_Comm PROBLEM_COMM, char *mesh_n, int mesh_f);
int read_mesh_coord_GMSH(MPI_Comm PROBLEM_COMM, char *mesh_n);
int read_mesh_coord_ALYA(MPI_Comm PROBLEM_COMM, char *mesh_n);

int read_physical_entities(MPI_Comm PROBLEM_COMM, char *mesh_n, int mesh_f);
int read_physical_entities_GMSH(MPI_Comm PROBLEM_COMM, char *mesh_n);
int read_physical_entities_ALYA(MPI_Comm PROBLEM_COMM, char *mesh_n);

int read_boundary(MPI_Comm PROBLEM_COMM, char *mesh_n,int mesh_f);
int read_boundary_GMSH(MPI_Comm PROBLEM_COMM, char *mesh_n);
int read_boundary_ALYA(MPI_Comm PROBLEM_COMM, char *mesh_n);

int part_mesh_PARMETIS(MPI_Comm *comm, FILE *time_fl, char *myname, double *centroid);
int swap_vectors_SCR( int *swap, int nproc, int n,  int *npe, 
    int *eptr, int *eind, int *PhysicalID,
    int *npe_swi, int *eind_swi, int *PhysicalID_swi,
    int *cuts_npe, int *cuts_eind );
int CSR_give_pointer( int e, int *npe, int *eind, int *p);
int clean_vector_qsort(int n, int *input, int **output, int *not_rep);
int give_repvector_qsort(MPI_Comm * comm, char *myname, int n, int *input, int **output, int *nrep);
int give_inter_sort(MPI_Comm *comm, char *myname, int *array1, int n1, int *array2, int n2, int **reps, int *nreps);
int calculate_ghosts(MPI_Comm * comm, char *myname);
int ownership_selec_rule( MPI_Comm *comm, int **repeated, int *nrep, int node, int *remoterank );
int is_in_vector(int val, int *vector, int size);
int reenumerate_PETSc(MPI_Comm PROBLEM_COMM);
int search_position_linear(int *array, int size, int val, int *pos);
int search_position_logn(int *array, int size, int val, int *pos);
int cmpfunc (const void * a, const void * b);
int cmpfunc_for_list (void * a, void * b);
int gmsh_npe(int code);
int gmsh_is_surf_elm(int code);
int gmsh_is_vol_elm(int code);
int set_id_on_material_and_boundary(MPI_Comm PROBLEM_COMM);
int get_bbox_limit_lengths(MPI_Comm PROBLEM_COMM, double *coord, int n, double *lx, double *ly, double *lz);
int get_bbox_local_limits(double *coord, int n, double *x, double *y, double *z);
int get_domain_center(MPI_Comm PROBLEM_COMM, double *coord, int n, double center[3]);

// spu_out.c
int spu_vtk_partition( char *vtkfile_n, MPI_Comm *comm );
int write_vtk(MPI_Comm PROBLEM_COMM, char *vtkfile_n, Vec *Displa, double *Strain, double *Stress);
int write_pvtu(MPI_Comm PROBLEM_COMM, char *name);
int write_vtu(MPI_Comm PROBLEM_COMM, char *name, Vec *x, Vec *b, double *strain, double *stress, double *energy);
int vtkcode(int dim,int npe);

// spu_time.c
int save_time(MPI_Comm *comm, const char *string, FILE *file, double dt);

// spu_assembly.c
int GetPETScIndeces(int *LocalNod, int n, int *local2PETSc, int *PETScIndex);
int GetElemCoord(int *LocalNod, int n, double elem_coor[8][3]);
int assembly_jacobian_sd(Mat *J);
int assembly_residual_sd(Vec *x, Vec *b);
int get_dsh(int gp, int npe, double coor[8][3], double ShapeDerivs[8][3], double *DetJac);
int GetB(int npe, double ShapeDerivs[8][3], double B[6][3*8]);
int GetWeight(int npe, double **wp);
int GetDsDe(int npe, double *ElemDisp, double DsDe[6][6]);
material_t * GetMaterial(int e);
int GetElemenDispls( int e, double *x, double *elem_disp );
int calc_strain_stress_energy(Vec *x, double *strain, double *stress, double *energy);
int calc_ave_strain_stress(MPI_Comm PROBLEM_COMM, Vec *x, double strain_ave[6], double stress_ave[6]);
int get_elem_coor(int e, double elem_coor[8][3]);
int is_inside_fiber_cilin(int e);
int get_centroid(int e, double centroid[3]);

// spu_util.c
int get_nods_bc(int **nods, int *nnods);
int get_nods_index(int *nods_bc, int nnods_bc, int *ix_loc, int *ix_glo);


#endif
