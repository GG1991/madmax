/*
   Routine to performe coloring on macro processes

   This is going to being done taking the mpi communication 
   strategy taken from the input.

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
   MACRO               MICRO  -> colors to assign
   |                   |    
   |_rank 0            |_rank 0
   |_rank 1            |_rank 1
   |_rank 2            |_rank 2
   |_ ...              |_ ...
   |_rank NPROC_MAC    |_rank NPROC_MIC


   ************************************************** 

   APPROACH MACMIC_2

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

   
*/

#include "sputnik.h"    // SPUTNIK common feature such as  variables to performe coloring

int mac_color(MPI_Comm world, spu_comm_t  spu_comm, MPI_Comm  * macro, MPI_Comm ** macmic_comm)
{
  
  /* 

     Author: Guido Giuntoli

     Performs the creation of the new communicator "macro" with its 
     inter-communicators (if exist according to "spu_comm") with 
     "micro"'s communicators


     Input: 

     world        : World communicator (macro & micro live here)

     spu_comm     : Communication approaches that stablish the rules 
                    for creating the communications between those 
		    communicators.


     Output:      

     macro        : "macro" communicator this communicator holds 
                    all those process that are going to solve the 
	            macro structure in a distributed way

     macmic_comm  : array of inter-comunicators for message interchange between 
                    the "macro" & and "micro" communicators


   */

    int  i, nproc, rank;
    int  sum_micro_proc;     // partial sum of microprocess per macroscopic process (vary with "i")
    int  sum_micro_proc_tot; // total sum of microprocess per macroscopic process
    int  grank_remote;
    int  color;

    // we are going to work with a "local" communicator saved in "world"
    // we create a newone called localworld
    MPI_Comm_size(world, &nproc);
    MPI_Comm_rank(world, &rank);

    if(spu_comm.approach_type == APPROACH_MACRO){
      //  In this approach it should exist only one communicator for all
      //  the macro processes, no micro process should exist here
      //  so all ranks have the same color = MACRO
      color = MACRO;
      MPI_Comm_split(world, color, 0, macro);

    }
    else{

      return 1;
    }


//    // we need an array of intercommunicator for each microscopic kind
//    macmic_comm = (MPI_Comm*)malloc(nkind_mic*sizeof(MPI_Comm));
//
//    // we performe the splitting (it should be made on all the processes)
//    // all the process that enters here have a different color {0 1 2 nproc_mac-1} 
//    *color = rank;
//    MPI_Comm_split(world, *color, 0, macro);
//
//    // compute total quantity of microscopic processes per macroscopic processes
//    sum_micro_proc_tot = 0;
//    for(i=0;i<nkind_mic;i++){
//	sum_micro_proc_tot += mpinfo_mic[i*2];
//    }
//
//    // verify that total number of process is valid for the communicators structure file
//    if(nproc != sum_micro_proc_tot * nproc_mac + nproc_mac){
//	return 1;
//    }
//
//    // we have to create the intercommunicators
//    sum_micro_proc = 0;
//    for( i=0 ; i < nkind_mic ; i++ ){
//        // the global rank formula aim to the first microscopic process that should be in 
//	// that remote communicator to perform the intercommunication
//	grank_remote = nproc_mac + rank * sum_micro_proc_tot + sum_micro_proc;
//	sum_micro_proc += mpinfo_mic[i*2]; 
//	MPI_Intercomm_create(*macro, 0, world, grank_remote, 0, &macmic_comm[i]);
//    }
//
//    printf("giant: rank = %d - nproc = %d - color = %d - rem_rank(to macro) %d\n",rank,nproc,*color,grank_remote);

    return 0;
}
