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
int              color;

int              nev;

char             input_n[NBUF];     // Input file name
char             mesh_n[NBUF];      // Mesh file name
char             mesh_f[4];         // Mesh format name

/*****************************************************************************************************
 Structures to performe de mesh partition with ParMetis
 We are going to adopt the same one to store the mesh
*****************************************************************************************************/

int            * elmdist;   // number of elements inside each procesor
int            * eptr;      // list of indeces of nodes inside eind
int            * eind;      // list of nodes for elem "i" is between 
                            // eind[eptr[i]] eind[eptr[i+1]] (not including)

/*****************************************************************************************************
   MACRO function definitions
*****************************************************************************************************/

// mac_main.c 
int main(int argc, char **args);

// mac_comm.c
int mac_comm_init(void);
