/*****************************************************************************************************
   MACRO external lybraries
*****************************************************************************************************/

#include "sputnik.h"
#include "macmic.h"  /* Routines inside sputinik that are common for <macro> and <micro> only */

#define TESTCOMM_NULL   -1
#define TESTCOMM_STRAIN  1

int    flag_testcomm;

/*****************************************************************************************************
   MACRO global variables 
*****************************************************************************************************/

int              mymicro_rank_worker;

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
int MacroAllocMatrixVector(MPI_Comm PROBLEM_COMM, int nlocal, int ntotal);

// mac_boundary.c
int MacroFillBoundary(MPI_Comm PROBLEM_COMM, list_t *boundary_list);
int MacroSetDisplacementOnBoundary( double time, Vec *x );
int MacroSetBoundaryOnJacobian( Mat *J );
int MacroSetBoundaryOnResidual( Vec *b );
int macro_parse_boundary(MPI_Comm PROBLEM_COMM, char *input);
int cmpfunc_mac_bou (void * a, void * b);
int SetGmshIDOnMaterialsAndBoundaries(MPI_Comm PROBLEM_COMM);
