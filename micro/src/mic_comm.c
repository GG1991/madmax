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

     Performs the creation of the new communicator "MICRO_COMM" with its 
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
                     between the "MICRO_COMM" and "macro_comm" communicators 
                     (all the macro communicators)

   */

  int  i, ierr, c, m;
  int  color;

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

    // MICRO_COMM creation
    MPI_Comm_split(world_comm, color, 0, &MICRO_COMM);
    ierr = MPI_Comm_size(MICRO_COMM, &nproc_mic);
    ierr = MPI_Comm_rank(MICRO_COMM, &rank_mic);

    /*
       these remote ranks correspond to
       the remote micro leaders
     */
    remote_ranks = malloc(nproc_mac_tot * sizeof(int));
    m = 0;
    for(i=0;i<nproc_wor;i++){
      if(id_vec[i] == MACRO){
	remote_ranks[m] = i;
	m++;
      }
    }
    printf("nproc_mac = %d m = %d\n",nproc_mac_tot,m);

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

/****************************************************************************************************/

int MicCommWaitStartSignal( MPI_Comm WORLD_COMM )
{

  /*
     The processes will wait here until they receive a signal
  */

  int ierr, signal;
  MPI_Status status;

  ierr = MPI_Recv(&signal, 1, MPI_INT, MyMacroRankLeader, 0, WORLD_COMM, &status); CHKERRQ(ierr);

  return(signal == MPI_MICRO_START) ? 0 : 1;
}

/****************************************************************************************************/

int MicCommRecvStrain( MPI_Comm WORLD_COMM )
{

  /*
     The processes will wait here until they receive a signal
  */

  int ierr;
  MPI_Status status;

  ierr = MPI_Recv(mac_eps, 6, MPI_DOUBLE, MyMacroRankLeader, 0, WORLD_COMM, &status); CHKERRQ(ierr);

  return 0;
}

/****************************************************************************************************/

int MicCommRecvGPnum( MPI_Comm WORLD_COMM )
{

  /*
     The processes will wait here until they receive a signal
  */

  int ierr;
  MPI_Status status;

  ierr = MPI_Recv(&mac_gp, 1, MPI_INT, MyMacroRankLeader, 0, WORLD_COMM, &status); CHKERRQ(ierr);

  return 0;
}

/****************************************************************************************************/

int MicCommSendAveStressAndTanTensor( MPI_Comm WORLD_COMM )
{

  /*
     Sends to macro leader the contiguos vector that has the averange 
     stress and tangent strain tensor calculated here
  */

  int ierr;
  MPI_Status status;

  if(rank_mic==0){
    ierr = MPI_Ssend(mic_stress_ttensor, 6+81, MPI_DOUBLE, MyMacroRankLeader, 0, WORLD_COMM); CHKERRQ(ierr);
  }

  return 0;
}

/****************************************************************************************************/

int MicCommSendAveStress( MPI_Comm WORLD_COMM )
{

  /*
     Sends to macro leader the averange Stress tensor calculated here
  */

  int ierr;
  MPI_Status status;

  if(rank_mic==0){
    ierr = MPI_Ssend(mic_stress, 6, MPI_DOUBLE, MyMacroRankLeader, 0, WORLD_COMM); CHKERRQ(ierr);
  }

  return 0;
}

/****************************************************************************************************/

int MicCommSendAveTTensor( MPI_Comm WORLD_COMM )
{

  /*
     Sends to macro leader the averange Tangent Tensor calculated here
  */

  int ierr;
  MPI_Status status;

  if(rank_mic==0){
    ierr = MPI_Ssend(mic_ttensor, 9, MPI_DOUBLE, MyMacroRankLeader, 0, WORLD_COMM); CHKERRQ(ierr);
  }

  return 0;
}
