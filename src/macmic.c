/*
   Common functions for <macro> & <micro> programs

   Author: Guido Giuntoli
   Date : 31-07-2017

 */

#include "macmic.h"
#include "micro.h"
#include "macro.h"

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

	while(fgets(buf,NBUF,file) != NULL){
	  ln ++;
	  data = strtok(buf," \n"); 
	  if(data){ 
	    if(strcmp(data,"scheme")) return 1;

	    data = strtok(NULL," \n"); if(!data) return 1;

	    if(!strcmp(data,"COUP_1")){
	      macmic.type = COUP_1;
	      return 0;
	    }
	    else{
	      return -1;
	    }
	  }
	}

	return -1; //no encontro $end_scheme

      } // inside $scheme
    }

  }
  return 1;
}
/****************************************************************************************************/
int MacMicColoring(MPI_Comm WORLD_COMM, int *color, coupling_t *macmic, MPI_Comm *LOCAL_COMM)
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

  int  i, ierr, c;
  int  nproc_wor, rank_wor;
  int  nmic_worlds;
  int  nproc_mac_tot = 0, nproc_mic_tot = 0, mic_nproc_group;
  int  nproc_local, rank_local;

  ierr = MPI_Comm_size(WORLD_COMM, &nproc_wor);
  ierr = MPI_Comm_rank(WORLD_COMM, &rank_wor);

  int  *id_vec = malloc(nproc_wor * sizeof(int));

  // Allgather of the id_vec array 
  ierr = MPI_Allgather(color,1,MPI_INT,id_vec,1,MPI_INT,WORLD_COMM);CHKERRQ(ierr);

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

  if(flag_coupling == true && (nproc_mic_tot==0 || nproc_mac_tot==0))
      SETERRQ(PETSC_COMM_WORLD,1,"SPUTNIK: Want to coupling executing ONLY ONE code.");


  if(flag_coupling == false){

    // LOCAL_COMM creation
    ierr = MPI_Comm_split(WORLD_COMM, *color, 0, LOCAL_COMM);CHKERRQ(ierr);

  }
  else if(macmic->type == COUP_1){

    if(nproc_mic_tot % nproc_mac_tot != 0){
      printf("mod(nproc_mic_tot,nproc_mac_tot) = %d\n",nproc_mic_tot % nproc_mac_tot);
      return 1;
    }
    mic_nproc_group = nproc_mic_tot / nproc_mac_tot;

    int im_leader;
    int mic_pos = -1;

    if(*color == MICRO){

      // determine MICRO color 
      c = -1;
      for(i=0;i<=rank_wor;i++){
	if(id_vec[i] == MICRO){
	  mic_pos++;
	  if( mic_pos % mic_nproc_group == 0){
	    c ++; 
	  }
	}
      }
      *color += c;

      im_leader = (mic_pos % mic_nproc_group == 0) ? 1 : 0;

      // determine MACRO leaders
      int mac_rank;
      i = 0;
      c = -1;
      while( i<nproc_wor ){
	if(id_vec[i] == MACRO){
	  c++;
	}
	if(c == mic_pos/mic_nproc_group){
	  mac_rank = i; 
	  break;
	}
	i++;
      }

      macmic->coup = malloc(sizeof(coupMic_1_t));
      ((coupMic_1_t*)macmic->coup)->mac_rank = mac_rank;
      ((coupMic_1_t*)macmic->coup)->im_leader = im_leader;

    } // in MICRO
    else{

      /*
	 The color is only one here (MACRO)
      */

      // determine MICRO leaders
      int mac_pos = 0;
      i = 0;
      while( i<rank_wor ){
	if(id_vec[i] == MACRO){
	  mac_pos ++;
	}
	i++;
      }

      int mic_rank;
      i = 0; c = 0; mic_pos = 0;
      while( i<nproc_wor ){
	if(id_vec[i] == MICRO){
	  if(c == mac_pos){
	    mic_rank = i;
	    break;
	  }
	  if(mic_pos % mic_nproc_group == 0) c++;
	  mic_pos ++;
	}
	i++;
      }

      macmic->coup = malloc(sizeof(coupMac_1_t));
      ((coupMac_1_t*)macmic->coup)->mic_rank = mic_rank;

    }

    // LOCAL_COMM creation
    ierr = MPI_Comm_split(WORLD_COMM, *color, 0, LOCAL_COMM);CHKERRQ(ierr);

  }
  else{
    return 1;
  }

  free(id_vec);

  return 0;
}
/****************************************************************************************************/
int MicCommWaitSignal( MPI_Comm WORLD_COMM, int *signal )
{

  /*
     The processes will wait here until they receive the signal
  */
  int ierr, remote_rank;
  MPI_Status status;
  *signal = -1;
  if(macmic.type == COUP_1){
    remote_rank = ((coupMic_1_t*)macmic.coup)->mac_rank;
    ierr = MPI_Recv(signal, 1, MPI_INT, remote_rank, 0, WORLD_COMM, &status);CHKERRQ(ierr);
  }
  else{
    return 1;
  }
  return 0;
}
/****************************************************************************************************/
int MacCommSendSignal( MPI_Comm WORLD_COMM, int signal )
{
  /*
     The processes will wait here until they receive the signal
  */
  int ierr, remote_rank;
  if(macmic.type == COUP_1){
    remote_rank = ((coupMac_1_t*)macmic.coup)->mic_rank;
    ierr = MPI_Ssend(&signal, 1, MPI_INT, remote_rank, 0, WORLD_COMM);CHKERRQ(ierr);
  }
  else{
    return 1;
  }
  return 0;
}
/****************************************************************************************************/
int MicCommRecvStrain( MPI_Comm WORLD_COMM, double strain[6] )
{
  /*
     The processes will wait here until they receive a signal
  */
  int ierr, remote_rank;
  MPI_Status status;
  if(macmic.type == COUP_1){
    remote_rank = ((coupMic_1_t*)macmic.coup)->mac_rank;
    ierr = MPI_Recv(strain, 6, MPI_DOUBLE, remote_rank, 0, WORLD_COMM, &status); CHKERRQ(ierr);
  }
  else{
    return 1;
  }
  return 0;
}
/****************************************************************************************************/
int MacCommSendStrain( MPI_Comm WORLD_COMM, double strain[6] )
{
  /*
     The processes will wait here until they receive the signal
  */
  int ierr, remote_rank;
  if(macmic.type == COUP_1){
    remote_rank = ((coupMac_1_t*)macmic.coup)->mic_rank;
    ierr = MPI_Ssend(strain, 6, MPI_DOUBLE, remote_rank, 0, WORLD_COMM);CHKERRQ(ierr);
  }
  else{
    return 1;
  }
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

//  ierr = MPI_Recv(&mac_gp, 1, MPI_INT, MyMacroRankLeader, 0, WORLD_COMM, &status); CHKERRQ(ierr);

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
//    ierr = MPI_Ssend(mic_stress_ttensor, 6+81, MPI_DOUBLE, MyMacroRankLeader, 0, WORLD_COMM); CHKERRQ(ierr);
  }

  return 0;
}
/****************************************************************************************************/
int MicCommSendStress( MPI_Comm WORLD_COMM, double stress[6] )
{
  /*
     Sends to macro leader the averange Stress tensor calculated here
  */
  int ierr, remote_rank;
  if(macmic.type == COUP_1){
    if(((coupMic_1_t*)macmic.coup)->im_leader){
      // only the micro leader sends the stress
      remote_rank = ((coupMic_1_t*)macmic.coup)->mac_rank;
      ierr = MPI_Ssend(stress, 6, MPI_DOUBLE, remote_rank, 0, WORLD_COMM);CHKERRQ(ierr);
    }
  }
  else{
    return 1;
  }
  return 0;
}
/****************************************************************************************************/
int MacCommRecvStress( MPI_Comm WORLD_COMM, double stress[6] )
{
  /*
     The processes will wait here until they receive the stress
  */
  int ierr, remote_rank;
  MPI_Status status;
  if(macmic.type == COUP_1){
    remote_rank = ((coupMac_1_t*)macmic.coup)->mic_rank;
    ierr = MPI_Recv(stress, 6, MPI_DOUBLE, remote_rank, 0, WORLD_COMM, &status); CHKERRQ(ierr);
  }
  else{
    return 1;
  }
  return 0;
}
/****************************************************************************************************/
int MicCommSendTTensor( MPI_Comm WORLD_COMM )
{

  /*
     Sends to macro leader the averange Tangent Tensor calculated here
  */

  int ierr;

  if(rank_mic==0){
//    ierr = MPI_Ssend(mic_ttensor, 9, MPI_DOUBLE, MyMacroRankLeader, 0, WORLD_COMM); CHKERRQ(ierr);
  }

  return 0;
}
