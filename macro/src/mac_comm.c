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
   |_color 0 |_color 3
   |_color 1 |_color 4
   |_color 2 |_color 5
             |_color 6 
	     |_color 7
             |_color 8

   Strategy the first "num_macro" process will correspond to MACRO, example : 

                             d=0           d=1           d=2
   color 0 MACRO_1 --> (( MICRO_K1_1 -> MICRO_K1_2 -> MICRO_K1_3 ))  color 3  b=0
	           |   
                    -> (( MICRO_K2_1 -> MICRO_K2_2 ))                color 4  b=0              
		             d=3           d=4


                             d=0           d=1           d=2
   color 1 MACRO_2 --> (( MICRO_K1_1 -> MICRO_K1_2 -> MICRO_K1_3 ))  color 5  b=1
	           |   
                    -> (( MICRO_K2_1 -> MICRO_K2_2 ))                color 6  b=1            
		             d=3           d=4


                             d=0           d=1           d=2
   color 2 MACRO_3 --> (( MICRO_K1_1 -> MICRO_K1_2 -> MICRO_K1_3 ))  color 7  b=2
	           |   
                    -> (( MICRO_K2_1 -> MICRO_K2_2 ))                color 8  b=2             
		             d=3           d=4


   Author: Guido Giuntoli
   
*/

#include "sputnik.h"    
#include "macro.h"    

int mac_comm_init(void)
{

  /* 

     Performs the creation of the new communicator "macro_comm" with its 
     inter-communicators with micro_comms if the scheme is MACRO_MICRO

     Are defined:

     id_vec        : vector of size nproc_tot that 
                     id_vec[rank_wor] = MACRO|MICRO

     macro_world   : communicator this communicator holds 
                     all those process that are going to solve the 
                     macro structure in a distributed way
     
     comm_macmic   : array of inter-comunicators for message interchange 
                     between the "macro_comm" and "micro_comm" communicators 
                     (all the micro worlds)

   */

  int  i, ierr, c, m;
  int  color;

  color = MACRO;

  if(scheme == MACRO_MICRO){

    MPI_Comm_split(world_comm, color, 0, &macro_comm);

    ierr = MPI_Comm_size(macro_comm, &nproc_mac);
    ierr = MPI_Comm_size(macro_comm, &rank_mac);

    // create the intercommunicators

    // fills the id_vec array (collective with micro code)
    id_vec = malloc(nproc_wor * sizeof(int));
    ierr = MPI_Allgather(&color,1,MPI_INT,id_vec,1,MPI_INT,world_comm);
    if(ierr){
      return 1;
    }

    // remote_rank array is filled 
    remote_ranks = malloc(nmic_worlds * sizeof(int));
    nproc_mic = nproc_mac = m = c = 0;
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
	  nproc_mic ++;
	}
	else{
	  nproc_mac++;
	}
    }
    if(!rank_mac){
      printf("number of macro process = %d",nproc_mac);
      printf("number of micro process = %d",nproc_mic);
    }

    // array of intercommunicator with micro-world
    macmic_comm = (MPI_Comm*)malloc(nmic_worlds * sizeof(MPI_Comm));
    for(m=0;m<nmic_worlds;m++){
      // args  =         local comm, local rank, global comm, remote rank, TAG, inter comm
      MPI_Intercomm_create(macro_comm, 0, world_comm, remote_ranks[m], 0, &macmic_comm[m]);
    }

  }
  else if(scheme == MACRO_MICRO){
     // TODO
  }
  else if(scheme == MICRO_ALONE){
     // TODO
  }
  else{
      return 1;
  }

  return 0;
}
