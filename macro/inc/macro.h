/*****************************************************************************************************
   MACRO external lybraries
*****************************************************************************************************/

#include "sputnik.h"

#define  DIM   3

/*****************************************************************************************************
   MACRO global variables 
*****************************************************************************************************/

MPI_Comm         MACRO_COMM;

int              nproc_mac;         //  # of processes (macro_comm)
int              rank_mac;          //  rank on macro comm
int             *remote_ranks;      //  remote ranks if micro processes
                                    //  to build the macmic_comm

int              nev;


/*****************************************************************************************************
   MACRO function definitions
*****************************************************************************************************/

// mac_main.c 
int main(int argc, char **args);

// mac_comm.c
int mac_comm_init(void);

// mac_boundary.c
int MacroSetBoundaryConditions( double time, Vec *x );
