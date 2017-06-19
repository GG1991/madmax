/*****************************************************************************************************
   MACRO external lybraries
*****************************************************************************************************/

#include "sputnik.h"

#define  DIM   3

/*****************************************************************************************************
   MACRO global variables 
*****************************************************************************************************/

//PetscErrorCode   ierr;
//PetscBool        couple_fl; //  couple flag = 0|1 (o:not coupled, 1:coupled)

MPI_Comm         macro_comm;
MPI_Comm        *macmic_comm;

int              nproc_wor;         //  # of processes (world_comm)
int              nproc_mac;         //  # of processes (macro_comm)
int              rank_wor;          //  rank on world comm
int              rank_mac;          //  rank on macro comm
int             *remote_ranks;      //  remote ranks if micro processes
                                    //  to build the macmic_comm

int              nev;

char             input_n[NBUF];     // Input file name
char             mesh_n[NBUF];      // Mesh file name
char             mesh_f[4];         // Mesh format name

/*****************************************************************************************************
   MACRO function definitions
*****************************************************************************************************/

// mac_main.c 
int main(int argc, char **args);

// mac_comm.c
int mac_comm_init(void);
