/*****************************************************************************************************
   MICRO external lybraries
*****************************************************************************************************/

#include "sputnik.h"
#include "macmic.h"  /* Routines inside sputinik that are common for <macro> and <micro> only */

#define  DIM   3

/*****************************************************************************************************
   MICRO global variables 
*****************************************************************************************************/

MPI_Comm         MICRO_COMM;

int              nproc_mac;         //  # of macro processes (world_comm)  
int              nproc_mic;         //  # of micro processes (micro_comm)
int              rank_mic;          //  rank on macro comm
int             *remote_ranks;      //  remote ranks if micro processes

int              nev;


/*****************************************************************************************************
   MICRO function definitions
*****************************************************************************************************/

// mic_main.c 
int main(int argc, char **args);

// mac_comm.c
int mic_comm_init(void);

// mic_alloc.c
int MicroAllocMatrixVector(MPI_Comm comm, int nlocal, int ntotal);
