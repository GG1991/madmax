/*****************************************************************************************************
   MICRO external lybraries
*****************************************************************************************************/

#include "petscksp.h"
#include "list.h"
#include <stdbool.h>

#define  DIM   3

/*****************************************************************************************************
   MICRO global variables 
*****************************************************************************************************/

PetscErrorCode   ierr;
PetscBool        couple_fl; //  couple flag = 0|1 (o:not coupled, 1:coupled)

int              nproc;     //  # of processors
int              rank;     

int              nev;

/*****************************************************************************************************
Structures to performe de mesh partition with ParMetis
*****************************************************************************************************/

int            * vtxdist;   // number of elements inside each procesor
int            * eptr;      // list of indeces of nodes inside eind
int            * eind;      /* list of nodes for elem "i" is between 
                               eind[eptr[i]] eind[eptr[i+1]] (not including) */

/*****************************************************************************************************
   MICRO function definitions
*****************************************************************************************************/

// mic_main.c 
int main(int argc, char **args);

// mic_mesh.c
int mic_rmsh(char *mesh_n, char *mesh_f);
int mic_readev_gmsh(char * mesh_n);

// mic_util.c
int cmp_int(void *a, void *b);

// mic_parser.c
