/*****************************************************************************************************
   MICRO external lybraries
*****************************************************************************************************/

#include "sputnik.h"
#include "macmic.h"  

/*****************************************************************************************************
   MICRO global variables 
*****************************************************************************************************/

/*
   Internal Variables on Gauss Points are going to be saved on 
   this vector
*/
double       *gauss_param_d;

/*
   Variables Send by <macro>
*/

MPI_Comm     MICRO_COMM;
int          rank_mic;          //  rank on macro comm
int          nproc_mic;         //  # of micro processes (MICRO_COMM)

#define TAYLOR       1
#define UNIF_STRAINS 2

int          homo_type;
int          flag_linear_micro; // 1 if all materials in micro are linear 0 otherwise
int          macro_gp;
int          first_time_c_homo_lineal_ask;
double       c_homo_lineal[36];

Mat          J;                // extended matrix for lagrange multipliers boundary setting
Vec          xe, re;           // extended distributed vectors for lagrange multipliers boundary setting
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
int mic_homogenize_taylor(MPI_Comm COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6]);
int mic_homogenize_unif_strain(MPI_Comm COMM, double strain_bc[6], double strain_ave[6], double stress_ave[6]);
int mic_homogenize_ld_lagran(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6]);
int voigt2mat(double voigt[6], double matrix[3][3]);

int mic_calc_c_homo(MPI_Comm MICRO_COMM, double strain_mac[6], double c_homo[36]);
int mic_calc_c_homo_lineal(MPI_Comm MICRO_COMM, double strain_mac[6], double c_homo_lineal[36]);
