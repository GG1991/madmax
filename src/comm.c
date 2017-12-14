#include "comm.h"


int comm_init_message(message_t *message){

  message->action = ACTION_NULL;

  const int nitems = 5;
  int block_lengths[5] = {1, 1, MAX_VOIGT, MAX_VOIGT, MAX_VOIGT*MAX_VOIGT};
  MPI_Datatype types[5] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint offsets[5];

  offsets[0] = offsetof(message_t, action);
  offsets[1] = offsetof(message_t, num_voigt);
  offsets[2] = offsetof(message_t, strain_mac);
  offsets[3] = offsetof(message_t, stress_ave);
  offsets[4] = offsetof(message_t, c_tangent_ave);


  int ierr = 0;
  ierr = MPI_Type_create_struct(nitems, block_lengths, offsets, types, &mpi_message_t);
  ierr = MPI_Type_commit(&mpi_message_t);

  return ierr;
}


int comm_finalize_message(void){

  int ierr = MPI_Type_free(&mpi_message_t);

  return ierr;
}


int comm_macro_send(message_t *message){

  int remote_rank = ((mac_coup_1_t*)macmic.coup)->mic_rank;
  int ierr = MPI_Ssend(message, 1, mpi_message_t, remote_rank, 0, WORLD_COMM);

  return ierr;
}


int comm_macro_recv(message_t *message){

  int remote_rank = ((mac_coup_1_t*)macmic.coup)->mic_rank;
  int ierr = MPI_Recv(message, 1, mpi_message_t, remote_rank, 0, WORLD_COMM, MPI_STATUS_IGNORE);

  return ierr;
}


int comm_micro_send(message_t *message){

  int ierr = 0;

  if(((mic_coup_1_t*)macmic.coup)->im_leader){
    int remote_rank = ((mic_coup_1_t*)macmic.coup)->mac_rank;
    ierr = MPI_Ssend(message, 1, mpi_message_t, remote_rank, 0, WORLD_COMM);
    if(ierr == 1) return 1;
  }

  return ierr;
}


int comm_micro_recv(message_t *message){

  int ierr = 0;

  if(((mic_coup_1_t*)macmic.coup)->im_leader){
    int remote_rank = ((mic_coup_1_t*)macmic.coup)->mac_rank;
    ierr = MPI_Recv(message, 1, mpi_message_t, remote_rank, 0, WORLD_COMM, MPI_STATUS_IGNORE);
    if(ierr == 1) return 1;
  }

  ierr = MPI_Bcast(message, 1, mpi_message_t, 0, MICRO_COMM);

  return ierr;
}


int macmic_coloring(MPI_Comm WORLD_COMM, int *color, coupling_t *macmic, MPI_Comm *LOCAL_COMM, bool flag_coupling){

  int  i, ierr, c;
  int  nproc_wor, rank_wor;
  int  nproc_mac_tot = 0, nproc_mic_tot = 0, mic_nproc_group;

  ierr = MPI_Comm_size(WORLD_COMM, &nproc_wor);
  ierr = MPI_Comm_rank(WORLD_COMM, &rank_wor);

  int  *id_vec = malloc(nproc_wor * sizeof(int));

  ierr = MPI_Allgather(color,1,MPI_INT,id_vec,1,MPI_INT,WORLD_COMM);
  if(ierr)
    return 1;

  nproc_mic_tot = nproc_mac_tot = 0;
  for(i=0;i<nproc_wor;i++){
    if(id_vec[i] == COLOR_MACRO){
      nproc_mac_tot++;
    }
    else if(id_vec[i] == COLOR_MICRO){
      nproc_mic_tot++;
    }
    else{
      return 1;
    }
  }

  if(flag_coupling == true && (nproc_mic_tot==0 || nproc_mac_tot==0)){
    return 1;
  }


  if(flag_coupling == false){

    ierr = MPI_Comm_split(WORLD_COMM, *color, 0, LOCAL_COMM); if(ierr) return 1;

  }
  else if(macmic->type == COUP_1){

    if(nproc_mic_tot % nproc_mac_tot != 0){
      return 1;
    }
    mic_nproc_group = nproc_mic_tot / nproc_mac_tot;

    int im_leader;
    int mic_pos = -1;

    if(*color == COLOR_MICRO){

      // determine MICRO color 
      c = -1;
      for(i=0;i<=rank_wor;i++){
	if(id_vec[i] == COLOR_MICRO){
	  mic_pos++;
	  if( mic_pos % mic_nproc_group == 0){
	    c ++; 
	  }
	}
      }
      *color += c;

      im_leader = (mic_pos % mic_nproc_group == 0) ? 1 : 0;

      // determine MACRO leaders
      int mac_rank = -1;
      i = 0;
      c = -1;
      while( i<nproc_wor ){
	if(id_vec[i] == COLOR_MACRO){
	  c++;
	}
	if(c == mic_pos/mic_nproc_group){
	  mac_rank = i; 
	  break;
	}
	i++;
      }
      if(mac_rank < 0) return 1;

      macmic->coup = malloc(sizeof(mic_coup_1_t));
      ((mic_coup_1_t*)macmic->coup)->mac_rank = mac_rank;
      ((mic_coup_1_t*)macmic->coup)->im_leader = im_leader;

    } // in MICRO
    else{

      int mac_pos = 0;
      i = 0;
      while( i<rank_wor ){
	if(id_vec[i] == COLOR_MACRO){
	  mac_pos ++;
	}
	i++;
      }

      int mic_rank = -1;
      i = 0; c = 0; mic_pos = 0;
      while( i<nproc_wor ){
	if(id_vec[i] == COLOR_MICRO){
	  if(c == mac_pos){
	    mic_rank = i;
	    break;
	  }
	  if(mic_pos % mic_nproc_group == 0) c++;
	  mic_pos ++;
	}
	i++;
      }
      if(mic_rank < 0) return 1;

      macmic->coup = malloc(sizeof(mac_coup_1_t));
      ((mac_coup_1_t*)macmic->coup)->mic_rank = mic_rank;

    }

    ierr = MPI_Comm_split(WORLD_COMM, *color, 0, LOCAL_COMM); if(ierr) return 1;

  }
  else{
    return 1;
  }

  free(id_vec);

  return 0;
}
