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

/*****************************************************************************************************
   MICRO function definitions
*****************************************************************************************************/

// mic_main.c 
int main(int argc, char **args);

// mac_comm.c
int mic_comm_init(void);

// mic_alloc.c
int MicroAllocMatrixVector(MPI_Comm comm, int nlocal, int ntotal);
