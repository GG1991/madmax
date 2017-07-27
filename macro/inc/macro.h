/*****************************************************************************************************
   MACRO external lybraries
*****************************************************************************************************/

#include "sputnik.h"

#define  FLAG_VTK_NONE 0
#define  FLAG_VTK_PART 1
#define  FLAG_VTK_DISP 2

/*****************************************************************************************************
   MACRO global variables 
*****************************************************************************************************/

MPI_Comm         MACRO_COMM;

int              nproc_mac;         //  # of processes (macro_comm)
int              rank_mac;          //  rank on macro comm
int             *remote_ranks;      //  remote ranks if micro processes
                                    //  to build the macmic_comm

int              nev;

// Matrices and vectors

Mat           A;                    /* Jacobian Matrix          */
Vec           x, dx, b;             /* Vectors unknowns and RHS */
KSP           ksp;                  /* linear solver context    */
KSPConvergedReason  reason;

double        *stress, *strain;     // Averange strain and stress on each element

int           flag_print_vtk;


/*****************************************************************************************************
   MACRO function definitions
*****************************************************************************************************/

// mac_main.c 
int main(int argc, char **args);

// mac_comm.c
int mac_comm_init(void);

// mac_alloc.c
int MacroAllocMatrixVector(MPI_Comm comm, int nlocal, int ntotal);

