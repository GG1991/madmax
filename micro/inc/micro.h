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

int             rank_mic;          //  rank on macro comm
int             nproc_mic;         //  # of micro processes (MICRO_COMM)

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
int *value_x0_ux, *value_x0_uy, *value_x0_uz;   
int *value_y0_ux, *value_y0_uy, *value_y0_uz; 
int *value_z0_ux, *value_z0_uy, *value_z0_uz; 
int *value_x1_ux, *value_x1_uy, *value_x1_uz; 
int *value_y1_ux, *value_y1_uy, *value_y1_uz; 
int *value_z1_ux, *value_z1_uy, *value_z1_uz; 


/*****************************************************************************************************
   MICRO function definitions
*****************************************************************************************************/

// mic_main.c 
int main(int argc, char **args);

// mac_comm.c
int mic_comm_init(void);

// mic_alloc.c
int MicroAllocMatrixVector(MPI_Comm comm, int nlocal, int ntotal);
