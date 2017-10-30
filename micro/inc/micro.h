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


/* structured mesh */
bool        flag_struct_mesh;
int         nx, ny, nz, nn;       // number of nodes 
int         nex, ney, nez;        // number of elements per direction 
int         nyl;
double      lx, ly, lz;           // dimain lenght
double      hx, hy, hz;           // element sizes 
int         npe;                  // Nodes per element
int         ngp;                  // Number of gauss points per element
double      **struct_sh;          // Shape functions
double      ***struct_dsh;        // Derivative shapes functions
double      *struct_wp;           // Weights
double      **struct_bmat;        // B matrix (Bu = epsilon)
int         *loc_index;           // local index vector for reading and writting
int         istart, iend;         // starting and ending index of matrices
int         nstart, nend;         // starting and ending index of nodes

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
