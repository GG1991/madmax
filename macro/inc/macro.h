/*****************************************************************************************************
   MACRO external lybraries
*****************************************************************************************************/

#include "sputnik.h"
#include "macmic.h"  /* Routines inside sputinik that are common for <macro> and <micro> only */

/*****************************************************************************************************
   MACRO global variables 
*****************************************************************************************************/

int              mymicro_rank_worker;

int              nev;

int              rank_mac;          //  rank on macro comm
int              nproc_mac;         //  # of macro processes (WORLD_COMM)  


/*****************************************************************************************************
   MACRO function definitions
*****************************************************************************************************/

// mac_main.c 
int main(int argc, char **args);

// mac_comm.c
int mac_comm_init(void);

// mac_alloc.c
int MacroAllocMatrixVector(MPI_Comm comm, int nlocal, int ntotal);

