/*****************************************************************************************************
   MICRO external lybraries
*****************************************************************************************************/

#include "petscksp.h"
#include "sputnik.h"
#include <stdbool.h>

#define  DIM   3

/*****************************************************************************************************
   MICRO global variables 
*****************************************************************************************************/

MPI_Comm         micro_comm;

int              nproc_mac;         //  # of macro processes (world_comm)  
int              nproc_mic;         //  # of micro processes (micro_comm)
int              rank_mic;          //  rank on macro comm
int             *remote_ranks;      //  remote ranks if micro processes

int              nev;

char             input_n[NBUF];     // Input file name
char             mesh_n[NBUF];      // Mesh file name
char             mesh_f[4];         // Mesh format name

/*****************************************************************************************************
   MICRO function definitions
*****************************************************************************************************/

// mic_main.c 
int main(int argc, char **args);

// mac_comm.c
int mic_comm_init(void);
