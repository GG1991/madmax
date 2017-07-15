/*
   Routines related to the communicators creation

   This is going to be done taking the mpi communication 
   strategy taken from the input.

   ************************************************** 

   MACRO_ALONE

   Note: This approach is for test the "macro" code
         It does not have any other argument

   ************************************************** 

   MICRO_ALONE

   Note: This approach is for test the "micro" code
         It does not have any other argument

   ************************************************** 

   MACRO_MICRO

   Note:

   SPUTNIK ( MPI_COMM_WORLD )
          |
	 / \
   MACRO     MICRO  (colors with code on sputnik.h)
   |         |    
   |_color 1 |_color 2
   |_color 1 |_color 2
   |_color 1 |_color 2
   |_...     |_...
   |_color 1 |_color 2

   color 1 (( MACRO_R0 ))   (( MICRO_K1_R0 -> MICRO_K1_R1 -> MICRO_K1_R2 ))  color 2  
	                   
                            (( MICRO_K2_R0 -> MICRO_K2_R1 ))                 color 3 


   color 1 (( MACRO_R1 ))   (( MICRO_K1_R0 -> MICRO_K1_R1 -> MICRO_K1_R2 ))  color 4  
	                    
                            (( MICRO_K2_R0 -> MICRO_K2_R1 ))                 color 5 


   color 1 (( MACRO_R2 ))   (( MICRO_K1_R0 -> MICRO_K1_R1 -> MICRO_K1_R2 ))  color 6
	                    
                            (( MICRO_K2_R0 -> MICRO_K2_R1 ))                 color 7


   Author: Guido Giuntoli
   
*/

#include "sputnik.h"    
#include "macro.h"    

int mac_comm_init(void)
{

  /* 

     Performs the creation of the new communicator "MACRO_COMM" with its 
     inter-communicators with micro_comms if the scheme is MACRO_MICRO

     Are defined:

     id_vec        : vector of size nproc_tot that 
                     id_vec[rank_wor] = MACRO|MICRO

     macro_world   : communicator this communicator holds 
                     all those process that are going to solve the 
                     macro structure in a distributed way
     
   */

  int  i, ierr, c, m;
  int  color;

  color = MACRO;

  if(scheme == MACRO_MICRO){

    // fills the id_vec array (collective with micro code)
    id_vec = malloc(nproc_wor * sizeof(int));
    ierr = MPI_Allgather(&color,1,MPI_INT,id_vec,1,MPI_INT,world_comm);
    if(ierr){
      return 1;
    }

    // MACRO_COMM creation
    MPI_Comm_split(world_comm, color, 0, &MACRO_COMM);
    ierr = MPI_Comm_size(MACRO_COMM, &nproc_mac);
    ierr = MPI_Comm_rank(MACRO_COMM, &rank_mac);

    // we count the number of processes 
    // that are with micro and macro from id_vec
    nproc_mic_tot = nproc_mac_tot = 0;
    for(i=0;i<nproc_wor;i++){
        if(id_vec[i] == MICRO){
	  nproc_mic_tot ++;
	}
	else{
	  nproc_mac_tot++;
	}
    }
    if(!rank_mac){
      printf("MACRO: # of macro process (id_vec) = %d\n",nproc_mac_tot);
      printf("MACRO: # of micro process (id_vec) = %d\n",nproc_mic_tot);
    }
    if(nproc_mic_tot % nproc_mic_group != 0){
      printf("MACRO: mod(nproc_mic,nproc_mic_group) = %d\n",nproc_mic_tot % nproc_mic_group);
      return 1;
    }
    nmic_worlds = nproc_mic_tot / nproc_mic_group;

    // remote_rank array is filled 
    // these ranks correspond to the 
    // micro processes' leaders
    remote_ranks = malloc(nmic_worlds * sizeof(int));
    m = c = 0;
    for(i=0;i<nproc_wor;i++){
        if(id_vec[i] == MICRO){
	  if(c == 0){
	    remote_ranks[m] = i;
	    c = nproc_per_mic[m % nstruc_mic];
	    m ++;
	  }
	  else{
	    c--;
	  }
	}
    }

  }

  else{
      return 1;
  }

  return 0;
}
