/*****************************************************************************************************
   MICRO external lybraries
*****************************************************************************************************/

#include "sputnik.h"

/*****************************************************************************************************
   MICRO global variables 
*****************************************************************************************************/

/*
   Internal Variables on Gauss Points are going to be saved on 
   this vector
*/
double       *gauss_param_d;

/* Variables Send by <macro> */

MPI_Comm     MICRO_COMM;
int          rank_mic;          //  rank on macro comm
int          nproc_mic;         //  # of micro processes (MICRO_COMM)

#define TAYLOR_S     1
#define TAYLOR_P     2
#define UNIF_STRAINS 3


/* Micro figure type */
#define CIRCULAR_FIBER 1

int          micro_type;
FILE         *fm_info;

/* structured mesh */
bool        flag_struct_mesh;
int         nx, ny, nz, nn;       // number of nodes 
int         nl;                   // number of local nodes
int         nex , ney , nez;      // total number of elements per direction 
int         nyl;                  // local number of nodes in y direction
int         ny_inf;               // inferior numeration in y direction
double      lx, ly, lz;           // dimain lenght
double      hx, hy, hz;           // element sizes 
double      vol_tot;              // total volume
double      vol_loc;              // local volume
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

Mat          J;                   // extended matrix for lagrange multipliers boundary setting
Vec          xe, re;              // extended distributed vectors for lagrange multipliers boundary setting
Vec          b1;

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

int micro_pvtu( char *name, double *strain, double *stress, double *energy);

void micro_print_info( void );
int  get_elem_properties( double * stress, double * strain, double * energy );
int  init_shapes( double ***sh, double ****dsh, double **wp );
