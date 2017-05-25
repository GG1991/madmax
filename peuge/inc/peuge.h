/*
   peuge external lybraries
*/

//#include "petscksp.h"
#include "list.h"

#define  DIM   3

/*
   peunt global variables 
*/

//PetscErrorCode   ierr;
//PetscBool        couple_fl; //  couple flag = 0|1 (o:not coupled, 1:coupled)

int              nproc;     //  # of processors
int              grank;     //  rank number global
int              lrank;     //  rank number local 
int              color;

int              nev;

//
// Structures to performe de mesh partition with ParMetis
// We are going to adopt the same one to store the mesh
//
int            * vtxdist;   // number of elements inside each procesor
int            * eptr;      // list of indeces of nodes inside eind
int            * eind;      // list of nodes for elem "i" is between eind[eptr[i]] eind[eptr[i+1]] (not including)

/*
   peuge function definitions
*/

// peu_main.c 
int main(int argc, char **args);

// peu_mesh.c
int peu_rmsh(char *mesh_n, char *mesh_f);
int peu_readev_gmsh(char * mesh_n);

// peu_util.c
int cmp_int(void *a, void *b);
