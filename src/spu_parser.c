/*

   This routine parse a file that specifies the conectivities between processes

   Recommend to call this file "mpi.dat"

   The format is:

   # MACRO INFORMATION
   # <# PROCESSES> 
      2  
   
   # MICRO INFORMATION
   # <# PROCESSES> 
      6  
   
   # <#KINDS> 
      3
   
   # <# PROCESSES KIND 0> <# PROCESSES KIND 1> <# PROCESSES KIND 2>
      1                    1                    1

   (we save this on nproc_k array and check dimensions)

*/

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#define SIZEBUF 128

int parse_mpi(const char mpi_file[], int nproc[2], int * nkind, int ** nproc_k)
{

    /*
       nproc[0]   : # of processes macro code
       nproc[1]   : # of processes micro code
       nkind      : # of microscopic kinds
       nproc_k[i] : # of process kind "i"
     */

    FILE * file = fopen(mpi_file,"r");
    char   buf[SIZEBUF];
    char * a;
    int    flg = 0, i, sum;

    if(!file){
	return 1;
    }
    
    while(fgets(buf,SIZEBUF,file) != NULL)
    {
	a = strtok(buf," \n");
	if(a != NULL){
	    if(a[0] != '#'){

		if(flg == 0){

		    /* read number of process that should be executing the macro structure */
		    nproc[0] = atoi(a); 
		    flg = 1;

		}
		else if(flg == 1){

		    /* read number of process that should be executing the micro structure */
		    nproc[1] = atoi(a); 
		    flg = 2;

		}
		else if(flg == 2){

		    /* read number of microscopic kinds that are going to be read and allocate nproc_k */
		    *nkind   = atoi(a); 
		    *nproc_k = malloc( *nkind * sizeof(int));
		    flg = 3;

		}
		else if(flg == 3){

		    /* read microscopic kinds */
		    i   = 0;
		    sum = 0;
		    while((a != NULL) && (i < *nkind)){
			(*nproc_k)[i] = atoi(a);
			sum += (*nproc_k)[i];
			a = strtok(NULL," \n");
			i++; 
		    }

                    // is all 0K ?
		    if(i != *nkind){
			return 1;
		    }

		    // checks if # of microscopic process fits for the process per kind defined 
		    if(nproc[1] != sum * nproc[0]){
		      return 1;
		    }
		    flg = 4;
		}
	    }
	}
    }

    // check if all the filds have been read
    if(flg != 4){
      return 1;
    }

    return 0;
}

