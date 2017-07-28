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

/*
   Internal Variables on Gauss Points are going to be saved on 
   this vector
*/
double          *gauss_param_d;

/*
   Variables Send by <macro>
*/

int             MyMacroRankLeader;

int             mac_gp;
double          mac_eps[6];
double          mic_stress[6];
double          mic_ttensor[9];
double          *mic_stress_ttensor;


/*****************************************************************************************************
   MICRO function definitions
*****************************************************************************************************/

// mic_main.c 
int main(int argc, char **args);

// mac_comm.c
int mic_comm_init(void);
int MicCommWaitStartSignal( MPI_Comm WORLD_COMM );
int MicCommRecvStrain( MPI_Comm WORLD_COMM );
int MicCommRecvGPnum( MPI_Comm WORLD_COMM );
int MicCommSendAveStressAndTanTensor( MPI_Comm WORLD_COMM );
int MicCommSendAveStress( MPI_Comm WORLD_COMM );
int MicCommSendAveTTensor( MPI_Comm WORLD_COMM );

// mic_alloc.c
int MicroAllocMatrixVector(MPI_Comm comm, int nlocal, int ntotal);
