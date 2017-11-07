/*****************************************************************************************************
   MICRO external lybraries
*****************************************************************************************************/
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
#include <gsl/gsl_linalg.h>

/*****************************************************************************************************
   MICRO global variables 
*****************************************************************************************************/

/*
   Internal Variables on Gauss Points are going to be saved on 
   this vector
*/
double       *gauss_param_d;

/* Variables Send by <macro> */

int          rank_mic;          //  rank on macro comm
int          nproc_mic;         //  # of micro processes (MICRO_COMM)
int          dim;               //  problem dimensions
int          nvoi;              //  number of voigt components (3 if dim=2, 6 if dim=3)

#define NBUF         256        // Buffers length for read from a file
#define TAYLOR_S     1
#define TAYLOR_P     2
#define UNIF_STRAINS 3

char         *myname;
int          nr_max_its;
double       nr_norm_tol;


int          micro_type;
FILE         *fm_info;

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
#define PRINT_VTK          1
#define PRINT_VTU          2
#define PRINT_VTKPART      4
#define PRINT_ALL          8

int         flag_print;
bool        flag_first_alloc;

/* Micro figure type */

#define CIRCULAR_FIBER 1

int         nx_fibers;
int         ny_fibers;
double      fiber_cilin_r;
double      fiber_cilin_center_devi[3];
double      center_domain[3];


/*****************************************************************************************************
   MICRO function definitions
*****************************************************************************************************/

// mic_main.c 
int main(int argc, char **args);

// mac_comm.c
int mic_comm_init(void);

// mic_alloc.c
int mic_alloc(MPI_Comm comm);

// mic_boundary.c
int mic_parse_boundary(MPI_Comm PROBLEM_COMM, char *input);
int mic_init_boundary_list(list_t *boundary_list);
int micro_check_physical_entities( list_t *physical_list );

int mic_homogenize(MPI_Comm COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6]);
int mic_homogenize_taylor(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6]);
int mic_homogenize_unif_strain(MPI_Comm COMM, double strain_bc[6], double strain_ave[6], double stress_ave[6]);
int mic_homog_us(MPI_Comm COMM, double strain_bc[6], double strain_ave[6], double stress_ave[6]);
int mic_homogenize_ld_lagran(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6]);
int voigt2mat(double voigt[6], double matrix[3][3]);

int mic_calc_stress_ave(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6]);
int mic_calc_c_homo(MPI_Comm MICRO_COMM, double strain_mac[6], double c_homo[36]);
int mic_calc_c_homo_lineal(MPI_Comm MICRO_COMM, double c_homo_lineal[36]);
int mic_check_l_us(void);

int micro_pvtu( char *name );

void micro_print_info( void );
int  get_elem_properties( void );
int  init_shapes( double ***sh, double ****dsh, double **wp );

