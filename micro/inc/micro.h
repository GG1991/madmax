/*
   MICRO external lybraries
*/
#include "petscksp.h"
#include "petscsys.h"
#include "slepceps.h"
#include "list.h"
#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "material.h"
#include "comm.h"
#include "trace.h"
#include "vtk.h"
#include <gsl/gsl_linalg.h>

#define TAYLOR_S     1
#define TAYLOR_P     2
#define UNIF_STRAINS 3

#define NBUF         256

char        filename[NBUF];    //  string for different purposes

/*
   Internal Variables on Gauss Points are going to be saved on 
   this vector
*/
double      *gauss_param_d;

/* Variables Send by <macro> */

int         rank_mic;          //  rank on macro comm
int         nproc_mic;         //  # of micro processes (MICRO_COMM)
int         dim;               //  problem dimensions
int         nvoi;              //  number of voigt components (3 if dim=2, 6 if dim=3)

char        *myname;
int         nr_max_its;
double      nr_norm_tol;

int         micro_type;
FILE        *fm_info;

/* structured mesh */
bool        flag_struct_mesh;
int         nx, ny, nz, nn;       // number of nodes 
int         nl;                   // number of local nodes
int         nelm;
int         nex , ney , nez;      // total number of elements per direction 
int         nyl;                  // local number of nodes in y direction
int         ny_inf;               // inferior numeration in y direction
double      lx, ly, lz;           // dimain lenght
double      hx, hy, hz;           // element sizes 
double      vol_tot;              // total volume
double      vol_loc;              // local volume
double      vol_elem;             // elemental volume
int         npe;                  // Nodes per element
int         ngp;                  // Number of gauss points per element
int         ngho;                 // Number of ghosts nodes
double      **struct_sh;          // Shape functions
double      ***struct_dsh;        // Derivative shapes functions
double      *struct_wp;           // Weights
double      ***struct_bmat;       // B matrix (Bu = epsilon)
double      *elem_disp;
int         *loc_elem_index;      // local elemental index vector for assembly and reading
int         *glo_elem_index;      // global elemental index vector for assembly and reading
int         istart, iend;         // starting and ending index of matrices
int         nstart, nend;         // starting and ending index of nodes
int         *dir_ix_loc;          // dirichlet indeces (local numeration)
int         *dir_ix_glo;          // dirichlet indeces (global numeration)
int         ndir_ix;              // number of dirichlet indeces
double      *coor_dir;            // coordinates of dirichlet nodes
double      *strain_gp;           // strain at Gauss point
double      *stress_gp;           // stress at Gauss point
int         *elem_type;           // type in each element (depends on the micro structure "micro_type")
double      *elem_strain;         // strain at each element
double      *elem_stress;         // stress at each element
double      *elem_energy;         // energy at each element
double      *center_coor;         // coordinates of the center

/* Taylor homogenization variables */

double       vi;                  // volumetric fraction of inclusions
double       vm;                  // volumetric fraction of matrix

int          homo_type;
int          flag_linear_micro;   // 1 if all materials in micro are linear 0 otherwise
int          macro_gp;
int          first_time_homo;
double       c_homo_lineal[36];
double       rho;

double       mat_fiber_t0[3];     // array of properties for type_0 fiber
double       mat_matrix_t0[3];    // array of properties for type_0 matrix

Mat          A, J;                // petsc matrices
Vec          x, b;                // petsc vectors
Vec          dx;
KSP          ksp;                 // linear solver context


#define PRINT_PETSC        0
#define PRINT_VTU          1

int         flag_print;
bool        flag_first_alloc;

/* Micro figure type */

#define CIRCULAR_FIBER 1

double      center_domain[3];

typedef struct _cilin_fiber_t
{
  double radius;
  double deviation[3];
  int    nx;
  int    ny;

}cilin_fiber_t;

cilin_fiber_t cilin_fiber;

int mic_homogenize(MPI_Comm COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6]);
int mic_homogenize_taylor(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6]);
int mic_homog_us(MPI_Comm COMM, double strain_bc[6], double strain_ave[6], double stress_ave[6]);
int voigt2mat(double voigt[6], double matrix[3][3]);

int micro_print_info( void );
int micro_pvtu( char *name );

int mic_calc_stress_ave(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6]);
int mic_calc_c_homo(MPI_Comm MICRO_COMM, double strain_mac[6], double c_homo[36]);
int mic_calc_c_homo_lineal(MPI_Comm MICRO_COMM, double c_homo_lineal[36]);

int get_local_elem_index( int e, int *loc_index );
int get_global_elem_index( int e, int *glo_elem_index );
int assembly_b( void );
int assembly_A( void );
int get_stress( int e , int gp, double *strain_gp , double *stress_gp );
int get_strain( int e , int gp, double *strain_gp );
int get_c_tan( const char * name , int e , int gp, double *strain_gp , double *c_tan );
int get_centroid_struct( int e, double *centroid );
int is_in_fiber( int e );
int strain_x_coord( double * strain , double * coord , double * u );
int get_node_local_coor( int n , double * coord );
int get_node_ghost_coor( int n , double * coord );
int get_local_elem_node( int e , int * n_loc );
int local_to_global_index( int local );
int get_averages( double * strain_ave, double * stress_ave );
int get_elem_type( int e , int *type );
int get_elem_properties( void );
int init_shapes( double ***sh, double ****dsh, double **wp );
