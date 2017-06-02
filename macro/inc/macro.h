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

int              nproc;     //  # of processors
int              rank;      //  rank number global
int              lrank;     //  rank number local 
int              color;

int              nev;

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

// mac_util.c
int cmp_int(void *a, void *b);

// mac_parser.c
