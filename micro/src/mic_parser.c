
/*

   This routine parse a file that specifies the connectivities between processes

   Recommend to call this file "mpi.dat"

   The format is:

   # N_M <#MACRO> <#SUBMAC_2>
   1  1
   (we save this on mpinfo_mic array and check dimensions)

   # <#KINDS> 
   2

   # kind_1 <#MICRO_1> <#SUBMIC_1> kind_2 <#MICRO_2> <#SUBMIC_2> ...
   1 1  1 1
   (we save this on mpinfo_mic array and check dimensions)

*/

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#define SIZEBUF 128

