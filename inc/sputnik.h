/* 
   Montseny common definitions 
*/

#define MACRO 1
#define MICRO 2

// mon_parse.c
int parse_mpi(const char mpi_file[], int * nproc_mac, int * nsubs_mac, int * nkind_mic, int ** mpinfo_mic);
