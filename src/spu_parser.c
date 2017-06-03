/*

   This routine parse a file that specifies the conectivities between processes

   Recommend to call this file "mpi.dat"

   The formats are classified with the key word APPROACH <APPROACH>:


   ************************************************** 
   APPROACH MACRO

   Note: This approach is for test the "macro" code
         It does not have any other argument

   ************************************************** 
   APPROACH MICRO

   Note: This approach is for test the "micro" code
         It does not have any other argument

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

int parse_mpi( const char mpi_file[], mpi_comm_t * mpi_comm )
{

    /*

       Parse the communication file in order to 
       stablish communications between processes

     */

    FILE * file = fopen(mpi_file,"r");
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
		    mpi_comm->approach_type = APPROACH_MACRO;
		    // parseamos y llenamos approach_1 (en este caso no tiene nada)
		    mpi_comm->approach = (void*)malloc(sizeof(approach_1_t));
		    memcpy(mpi_comm->approach,approach_1);
		}
		else if(strcmp(data,"MICRO")){
		    approach_2_t approach_2;
		    mpi_comm->approach_type = APPROACH_MICRO;
		    // parseamos y llenamos approach_2 (en este caso no tiene nada)
		    mpi_comm->approach = (void*)malloc(sizeof(approach_2_t));
		    memcpy(mpi_comm->approach,approach_2);
		}
		else{
		    PetscPrintf(PETSC_COMM_WORLD,"parse_mpi.c: non valid option at line %d\n",ln);
		}

	    }
	}
    }

    return 0;
}

