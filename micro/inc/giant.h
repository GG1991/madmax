/*****************************************************************************************************
   giant external lybraries
*****************************************************************************************************/

#include "petscksp.h"
#include "list.h"
#include <stdbool.h>

#define  DIM   3

/*****************************************************************************************************
   giant global variables 
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
   giant function definitions
*****************************************************************************************************/

// gia_main.c 
int main(int argc, char **args);

// gia_mesh.c
int gia_rmsh(char *mesh_n, char *mesh_f);
int gia_readev_gmsh(char * mesh_n);

// gia_util.c
int cmp_int(void *a, void *b);

// gia_parser.c
int gia_parse_mpi(const char mpi_file[], int * nproc_mac, int * nsubs_mac, int * nkind_mic, int ** mpinfo_mic);
