/*
   Common functions for <macro> & <micro> programs
 */

#include "macmic.h"

#define NBUF 64

int MacMicInitGaussStructure(int *eptr, int nelm)
{
  /*
     Allocates memory for Gauss point global structure <gauss>
     of type <gauss_t>
   */
  int i, ngp, ngauss = 0;
  for(i=1;i<(nelm+1);i++){
    // we assume that each element has the same number of gauss points as vertices 
    ngp = eptr[i] - eptr[i-1]; 
    ngauss += ngp; 
  }
  gauss = malloc(ngauss * sizeof(gauss_t)); if(!gauss) return 1;

  for(i=0;i<ngauss;i++){
    gauss[i].param_d = NULL;
  }

  return 0;
}

/****************************************************************************************************/

int MacMicParseScheme( char *input )
{

  /*
     Parse the communication scheme from input file

     Examples>

     $scheme
     scheme  COUP_1
     $end_scheme
   */

  FILE * file = fopen(input,"r");
  char   buf[NBUF];
  char * data;
  int    ln, i;


  if(!file) return 1;

  ln = 0;
  while(fgets(buf,NBUF,file) != NULL)
  {

    ln ++;
    data = strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$scheme")){

	fgets(buf,NBUF,file);
	ln ++;
	data = strtok(buf," \n"); 
	if(!data) return 1;
	if(strcmp(data,"scheme")) return 1;

	data = strtok(NULL," \n"); 
	if(!data) return 1;
	if(!strcmp(data,"COUP_1")){
	  macmic.type = COUP_1;
	  return 0;
	}
	else{
	  return -1;
	}

	return -1; //no encontro $end_scheme

      } // inside $scheme
    }

  }
  return 1;
}

/****************************************************************************************************/

int MacMicColoring(int id)
{

  /* 

     Creates the new communicators "MACRO_COMM" & "MICRO_COMM" 

     Are defined:

     id_vec        > vector of size nproc_tot that 
     id_vec[rank_wor] = MACRO|MICRO

     <--------nproc_wor-------------->

     [ 1 1 2 1 2 2 2 1 2 1 2 2 ... 2 ]

     macro_world (=1) > this communicator holds all those processes
     that belong to the same macro program are going
     to solve the micro structure in a distributed way

     micro_world (=2) > this communicator holds all those processes
     that belong to the same micro program and are going
     to solve the micro structure in a distributed way. They are joined
     in groups that are going to be coubled with macro program.

   */

  int  i, ierr, c, m;
  int  color;

  if(id == MICRO){

    color = MICRO;

    if(macmic.type == COUP_1){

      // fills the id_vec array (collective with macro code)
      id_vec = malloc(nproc_wor * sizeof(int));
      ierr = MPI_Allgather(&color,1,MPI_INT,id_vec,1,MPI_INT,WORLD_COMM);
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
      MPI_Comm_split(WORLD_COMM, color, 0, &MICRO_COMM);
      ierr = MPI_Comm_size(MICRO_COMM, &nproc_mic);
      ierr = MPI_Comm_rank(MICRO_COMM, &rank_mic);

      /*
	 these remote ranks correspond to
	 the remote MACRO leaders
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

      /* 
	 Initialize structure <coupling_t macmic>
       */
      macmic.type = COUP_1;
      macmic.coup = malloc(sizeof(coupMic_1_t));
      ((coupMic_1_t*)macmic.coup)->myMacLeader = remote_ranks[rank_mic];
      ((coupMic_1_t*)macmic.coup)->imMicLeader = 0;

    }
    else
      return 1;

  }else if(id==MACRO){

    color = MACRO;

    if(macmic.type == COUP_1){

      // fills the id_vec array (collective with MICRO code)
      id_vec = malloc(nproc_wor * sizeof(int));
      ierr = MPI_Allgather(&color,1,MPI_INT,id_vec,1,MPI_INT,WORLD_COMM);CHKERRQ(ierr);

      // MACRO_COMM creation
      MPI_Comm_split(WORLD_COMM, color, 0, &MACRO_COMM);
      ierr = MPI_Comm_size(MACRO_COMM, &nproc_mac);
      ierr = MPI_Comm_rank(MACRO_COMM, &rank_mac);

      /* 
	 we count the number of processes 
	 that are with micro and macro from id_vec
       */
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

      /* 
	 remote_rank array is filled 
	 these ranks correspond to the 
	 micro processes' leaders
       */ 
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

      /* 
	 Initialize structure <coupling_t macmic>
       */
      macmic.type = COUP_1;
      macmic.coup = malloc(sizeof(coupMac_1_t));
      ((coupMac_1_t*)macmic.coup)->myMicWorker = remote_ranks[rank_mac];

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

  if(rank_mic==0){
    ierr = MPI_Ssend(mic_ttensor, 9, MPI_DOUBLE, MyMacroRankLeader, 0, WORLD_COMM); CHKERRQ(ierr);
  }

  return 0;
}
