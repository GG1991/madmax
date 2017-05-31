/* 
   SPUTNIK common definitions 
*/

#define MACRO 1
#define MICRO 2

// spu_parse.c
int parse_mpi(const char mpi_file[], int nproc[2], int * nkind, int ** nproc_k);
