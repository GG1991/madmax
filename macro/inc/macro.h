/*****************************************************************************************************
   MACRO external lybraries
*****************************************************************************************************/

#include "sputnik.h"
#include "comm.h"
#include "util.h"
#include "fun.h"

#define NORMAL      1
#define EIGENSYSTEM 2
#define TEST_COMM   3

#define NBUF        256            // buffer length

int         macro_mode;

int         nr_max_its;          // newton raphson maximum number of iterations
double      nr_norm_tol;         // newton raphson minimum tolerance

double      **struct_sh;          // Shape functions
double      ***struct_dsh;        // Derivative shapes functions
double      *struct_wp;           // Weights
double      **bmat;               // B matrix (Bu = epsilon)
int         npe_max;              // maximum number of nodes per element
double      *elem_disp;
int         *loc_elem_index;      // local elemental index vector for assembly and reading
int         *glo_elem_index;      // global elemental index vector for assembly and reading
double      *strain_gp;           // strain at Gauss point
double      *stress_gp;           // stress at Gauss point
int         *elem_type;           // type in each element
double      *elem_strain;         // strain at each element
double      *elem_stress;         // stress at each element
double      *elem_energy;         // energy at each element

f1d_t       func_bc;              // boundary function

int         mymicro_rank_worker;

int         rank_mac;             //  rank on macro comm
int         nproc_mac;            //  # of macro processes (WORLD_COMM)  

char        filename[NBUF];       //  string for different purposes
char        mesh_n[NBUF];         //  mesh path
int         mesh_f;               //  mesh format

typedef struct bound_t_
{
  char    *name;
  int     kind;
  int     *fnum;
  int     nix;
  int     *disp_ixs;
  double  *disp_val;

}bound_t;

list_t boundary_list;

// mac_main.c 
int main(int argc, char **args);

// mac_comm.c
int mac_comm_init(void);

// mac_alloc.c
int mac_alloc(MPI_Comm PROBLEM_COMM);

// mac_boundary.c
int mac_init_boundary(MPI_Comm PROBLEM_COMM, list_t *boundary_list);
int MacroSetDisplacementOnBoundary( double time, Vec *x );
int MacroSetBoundaryOnJacobian( Mat *J );
int MacroSetBoundaryOnResidual( Vec *b );
int macro_parse_boundary(MPI_Comm PROBLEM_COMM, char *input);
int cmpfunc_mac_bou (void * a, void * b);
int SetGmshIDOnMaterialsAndBoundaries(MPI_Comm PROBLEM_COMM);
