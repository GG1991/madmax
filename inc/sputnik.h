/*****************************************************************************************************
   SPUTNIK external lybraries
*****************************************************************************************************/

#include "petscksp.h"
#include "list.h"
#include "mpi.h"
#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#define MACRO           1
#define MICRO           2
#define MESH_N_LENGTH   64
#define BUF_N_LENGTH    128

#define APPROACH_MACRO     1
#define APPROACH_MICRO     2
#define APPROACH_MACMIC_1  3


/*****************************************************************************************************
   SPUTNIK structures
*****************************************************************************************************/

typedef struct spu_comm_t_{

   int    approach_type;
   void * approach;

}spu_comm_t;

typedef struct approach_1_t_{


}approach_1_t;

typedef struct approach_2_t_{


}approach_2_t;

typedef struct approach_3_t_{

    int nproc_mac;
    int nproc_mic;

}approach_3_t;


/*****************************************************************************************************
   SPUTNIK function definitions
*****************************************************************************************************/

// spu_parser.c
int parse_mpi( const char mpi_file[], spu_comm_t * spu_comm );

// spu_mesh.c
int read_mesh(MPI_Comm comm, char *mesh_n, char *mesh_f, int ** elmdist, int ** eptr, int ** eind);
int read_mesh_CSR_GMSH(MPI_Comm comm, char *mesh_n, int ** elmdist, int ** eptr, int ** eind);
