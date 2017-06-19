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

int mac_init_comm(void)
{

  /* 

     Performs the creation of the new communicator "macro_comm" with its 
     inter-communicators with micro_comms if the scheme is MACRO_MICRO

     Are defined:

     macro_world   : communicator this communicator holds 
     all those process that are going to solve the 
     macro structure in a distributed way
     
     comm_macmic   : array of inter-comunicators for message interchange 
     between the "macro_comm" and "micro_comm" communicators 
     (all the micro worlds)

   */

  int  i;

  // we are going to work with a "local" communicator saved in "world"
  // we create a newone called localworld
  MPI_Comm_size(world, &nproc);
  MPI_Comm_rank(world, &rank);

  if(scheme == MACRO_MICRO){

    MPI_Comm_split(world_comm, MACRO, 0, macro_comm);

    ierr = MPI_Comm_size(macro_comm, &nproc_mac);
    ierr = MPI_Comm_size(macro_comm, &rank_mac);

    // create the intercommunicators
    macmic_comm = (MPI_Comm*)malloc(nmic_worlds * sizeof(MPI_Comm));
//    MPI_Intercomm_create(*macro, 0, world, grank_remote, 0, &macmic_comm[i]);

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
