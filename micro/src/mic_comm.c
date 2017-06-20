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
#include "micro.h"    

int mic_comm_init(void)
{

  /* 

     Performs the creation of the new communicator "micro_comm" with its 
     inter-communicators with macro_comms if the scheme is MACRO_MICRO

     Are defined:

     id_vec        : vector of size nproc_tot that 
                     id_vec[rank_wor] = MACRO|MICRO

     <--------nproc_wor-------------->

     [ 1 1 2 1 2 2 2 1 2 1 2 2 ... 2 ]

     micro_world   : this communicator holds all those processes
                     that belong to the same micro-world are going
		     to solve the micro structure in a distributed way
     
     micmac_inter_comm   : array of inter-comunicators for message interchange 
                     between the "micro_comm" and "macro_comm" communicators 
                     (all the macro communicators)

   */

  int  i, ierr, c, m;
  int  color;
  int  size_inter, rank_inter;

  color = MICRO;

  if(scheme == MACRO_MICRO){

    // fills the id_vec array (collective with macro code)
    id_vec = malloc(nproc_wor * sizeof(int));
    ierr = MPI_Allgather(&color,1,MPI_INT,id_vec,1,MPI_INT,world_comm);
    if(ierr){
      return 1;
    }

    // determines nproc_mac_tot and nproc_mic_tot
    nproc_mic_tot = nproc_mac_tot = 0;
    for(i=0;i<nproc_wor;i++){
      if(id_vec[i] == MACRO){
	nproc_mac_tot++;
      }
      else if(id_vec[i] == MICRO){
	nproc_mic_tot++;
      }
      else{
	return 1;
      }
    }
    if(nproc_mic_tot % nproc_mic_group != 0){
      printf("mod(nproc_mic,nproc_mic_group) = %d\n",nproc_mic_tot % nproc_mic_group);
      return 1;
    }
    nmic_worlds = nproc_mic_tot / nproc_mic_group;
    
    // color for MICRO
    m = c = 0;
    for(i=0;i<rank_wor;i++){
      if(id_vec[i] == MICRO){
	if( m == nproc_per_mic[c % nstruc_mic] ){
	  c ++;
	  m = 0;
	}
	else{
	  m ++;
	}
      }
    }
    color += c;

    // micro_comm creation
    MPI_Comm_split(world_comm, color, 0, &micro_comm);
    ierr = MPI_Comm_size(micro_comm, &nproc_mic);
    ierr = MPI_Comm_rank(micro_comm, &rank_mic);

    // remote ranks
    // these remote ranks correspond to
    // the remote micro leaders
    remote_ranks = malloc(nproc_mac * sizeof(int));
    m = 0;
    for(i=0;i<nproc_wor;i++){
      if(id_vec[i] == MACRO){
	remote_ranks[m] = i;
	m++;
      }
    }

  }
  else if(scheme == MACRO_ALONE){
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
