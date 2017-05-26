
/*

   This routine parse a file that specifies the conectivities between processes

   Recommend to call this file "mpi.dat"

   The format is:

   # N_M <#MACRO> <#SUBMAC_2>
   1  1
   (we save this on mpinfo_mic array and check dimensions)

   # <#KINDS> 
   1

   # kind_1 <#MICRO_1> <#SUBMIC_1> -> kind_2 <#MICRO_2> <#SUBMIC_2> ...
   1 1  1 1
   (we save this on mpinfo_mic array and check dimensions)

*/

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#define SIZEBUF 128

int peu_parse_mpi(const char mpi_file[], int * nproc_mac, int * nsubs_mac, int * nkind_mic, int ** mpinfo_mic)
{
    FILE * file = fopen(mpi_file,"r");
    char   buf[SIZEBUF];
    char * a;
    int    flg = 0, i;

    if(!file){
	return 1;
    }
    
    while(fgets(buf,SIZEBUF,file) != NULL)
    {
	a = strtok(buf," \n");
	if(a[0] != '#' && a != NULL){

	    if(flg == 0){

                /* read number of process that should be executing macro structure */
                *nproc_mac = atoi(a); 
		
		a = strtok(NULL," \n");
		if(a == NULL){
		    return 1;
		}
                *nsubs_mac = atoi(a);

		flg = 1;
		break;

	    }
	    else if(flg == 1){

                /* read number of microscopic kinds that are going to be read and allocate mpinfo_mic */
                *nkind_mic  = atoi(a); 
		*mpinfo_mic = calloc( *nkind_mic * 2,sizeof(int));
		flg = 2;
		break;

	    }
	    else if(flg == 2){

                /* read microscopic kinds */
		i = 0;
		while(a != NULL){
		    *mpinfo_mic[i] = atoi(a);
		    a = strtok(NULL," \n");
		    i++; 
		}
		if(i % 3 != 0){
		    return 1;
		}
		if(i / 3 != *nkind_mic){
		    return 1;
		}
		break;

	    }
	}
    }
    return 0;
}

