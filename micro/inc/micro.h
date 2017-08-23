/*****************************************************************************************************
   MICRO external lybraries
*****************************************************************************************************/

#include "sputnik.h"
#include "macmic.h"  /* Routines inside sputinik that are common for <macro> and <micro> only */

/*****************************************************************************************************
   MICRO global variables 
*****************************************************************************************************/

/*
   Internal Variables on Gauss Points are going to be saved on 
   this vector
*/
double          *gauss_param_d;

/*
   Variables Send by <macro>
*/

MPI_Comm     MICRO_COMM;
int          rank_mic;          //  rank on macro comm
int          nproc_mic;         //  # of micro processes (MICRO_COMM)

/*
   indeces for specifying boundary conditions
*/
double LX, LY, LZ;
int *index_x0_ux, *index_x0_uy, *index_x0_uz, nnods_x0; 
int *index_y0_ux, *index_y0_uy, *index_y0_uz, nnods_y0; 
int *index_z0_ux, *index_z0_uy, *index_z0_uz, nnods_z0; 
int *index_x1_ux, *index_x1_uy, *index_x1_uz, nnods_x1; 
int *index_y1_ux, *index_y1_uy, *index_y1_uz, nnods_y1; 
int *index_z1_ux, *index_z1_uy, *index_z1_uz, nnods_z1; 
double *value_x0_ux, *value_x0_uy, *value_x0_uz;   
double *value_y0_ux, *value_y0_uy, *value_y0_uz; 
double *value_z0_ux, *value_z0_uy, *value_z0_uz; 
double *value_x1_ux, *value_x1_uy, *value_x1_uz; 
double *value_y1_ux, *value_y1_uy, *value_y1_uz; 
double *value_z1_ux, *value_z1_uy, *value_z1_uz; 
int P000[3], P000_ismine, P100[3], P100_ismine, P010[3], P010_ismine;
double PVAL[3];

/*
   Homogenization Variables
*/

#define HOMO_TAYLOR  1
#define HOMO_LINEAR  2
#define HOMO_PERIOD  3
#define HOMO_EXP     9

int homo_type;


/*****************************************************************************************************
   MICRO function definitions
*****************************************************************************************************/

// mic_main.c 
int main(int argc, char **args);

// mac_comm.c
int mic_comm_init(void);

// mic_alloc.c
int MicroAllocMatrixVector(MPI_Comm comm, int nlocal, int ntotal);

// mic_boundary.c
int micro_init_boundary_list(list_t *boundary_list);
int micro_init_boundary(MPI_Comm PROBLEM_COMM, list_t *boundary_list);
int micro_check_physical_entities( list_t *physical_list );

#define DISPLACE    0
#define JACOBIAN    1
#define RESIDUAL    2
#define SET_DISPLACE    1
#define SET_JACOBIAN    2
#define SET_RESIDUAL    4
#define SET_JACRES      6

int micro_apply_bc(int dir, double strain[6], Vec *x, Mat *J, Vec *b, int flag);
int micro_apply_bc_linear(double strain_mac[6], Vec *x, Mat *J, Vec *b, int flag);
int micro_homogenize(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6]);
int micro_homogenize_taylor(MPI_Comm MICRO_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6]);
int micro_homogenize_linear(MPI_Comm MICRO_COMM, int i, double strain_bc[6], double strain_ave[6], double stress_ave[6]);
int voigt2mat(double voigt[6], double matrix[3][3]);
