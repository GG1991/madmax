/*

   This routine parse a file that specifies the conectivities between processes.
   The "mpi_file" is read only if the option -c <mpi_file> is given.

   Recommendation: call this file "spu_comm.dat"

   The formats are classified with the key word APPROACH <APPROACH>:

   ************************************************** 
   APPROACH MACMIC_1
   NPROC_MAC <NPROC_MAC> 
   NPROC_MIC <NPROC_MIC> 

   Note: This approach is to perform a coupling between 
         the "macro" code and the "micro" code with two 
	 communicators as:

   
   SPUTNIK ( MPI_COMM_WORLD )
                |
	  _____/ \____
   MACRO               MICRO  
   |                   |    
   |_rank 0            |_rank 0
   |_rank 1            |_rank 1
   |_rank 2            |_rank 2
   |_ ...              |_ ...
   |_rank NPROC_MAC    |_rank NPROC_MIC


   ************************************************** 
   
*/

#include "sputnik.h"

int parse_spu_comm( const char spu_comm_file[], spu_comm_t * spu_comm )
{

    /*

       Parse the communication file in order to 
       stablish communications between processes

     */

    FILE * file = fopen(spu_comm_file,"r");
    char   buf[BUF_N_LENGTH];
    char * data;
    int    ln;

    if(!file){
	return 1;
    }
    
    ln = 0;
    while(fgets(buf,BUF_N_LENGTH,file) != NULL)
    {
	ln ++;
	data = strtok(buf," \n");
	if(data != NULL){
	    if(data[0] != '#'){

		if(strcmp(data,"MACRO")){
		    approach_1_t approach_1;
		    spu_comm->approach_type = APPROACH_MACRO;
		    // parseamos y llenamos approach_1 (en este caso no tiene nada)
		    spu_comm->approach = (void*)malloc(sizeof(approach_1_t));
		    memcpy((void*)spu_comm->approach,(void*)&approach_1,sizeof(approach_1));
		}
		else if(strcmp(data,"MICRO")){
		    approach_2_t approach_2;
		    spu_comm->approach_type = APPROACH_MICRO;
		    // parseamos y llenamos approach_2 (en este caso no tiene nada)
		    spu_comm->approach = (void*)malloc(sizeof(approach_2_t));
		    memcpy((void*)&spu_comm->approach,(void*)&approach_2,sizeof(approach_2));
		}
		else{
		    PetscPrintf(PETSC_COMM_WORLD,"parse_mpi.c: non valid option at line %d\n",ln);
		}

	    }
	}
    }

    return 0;
}

