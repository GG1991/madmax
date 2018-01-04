#include "comm.h"

int comm_init_message(message_t *message)
{
  message->action = ACTION_NULL;

  const int nitems = 6;
  int block_lengths[6] = {1, 1, MAX_VOIGT, MAX_VOIGT, MAX_VOIGT*MAX_VOIGT, 1};
  MPI_Datatype types[6] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint offsets[6];

  offsets[0] = offsetof(message_t, action);
  offsets[1] = offsetof(message_t, num_voigt);
  offsets[2] = offsetof(message_t, strain_mac);
  offsets[3] = offsetof(message_t, stress_ave);
  offsets[4] = offsetof(message_t, c_tangent_ave);
  offsets[5] = offsetof(message_t, rho);


  int ierr = 0;
  ierr = MPI_Type_create_struct(nitems, block_lengths, offsets, types, &mpi_message_t);
  ierr = MPI_Type_commit(&mpi_message_t);

  return ierr;
}

int comm_finalize_message(void)
{
  return MPI_Type_free(&mpi_message_t);
}

int comm_macro_send(message_t *message, comm_t *comm)
{
  return MPI_Ssend(message, 1, mpi_message_t, comm->micro_slave, 0, WORLD_COMM);
}

int comm_macro_recv(message_t *message, comm_t *comm)
{
  return MPI_Recv(message, 1, mpi_message_t, comm->micro_slave, 0, WORLD_COMM, MPI_STATUS_IGNORE);
}

int comm_micro_send(message_t *message, comm_t *comm)
{
  return MPI_Ssend(message, 1, mpi_message_t, comm->macro_leader, 0, WORLD_COMM);
}

int comm_micro_recv(message_t *message, comm_t *comm)
{
  return MPI_Recv(message, 1, mpi_message_t, comm->macro_leader, 0, WORLD_COMM, MPI_STATUS_IGNORE);
}

int comm_coloring(MPI_Comm WORLD_COMM, comm_t *comm, MPI_Comm *LOCAL_COMM)
{
  int  nproc_world, rank_world;
  MPI_Comm_size(WORLD_COMM, &nproc_world);
  MPI_Comm_rank(WORLD_COMM, &rank_world);

  int *id_vec = malloc(nproc_world*sizeof(int));
  int ierr = MPI_Allgather(&comm->color, 1, MPI_INT, id_vec, 1, MPI_INT, WORLD_COMM);

  int nproc_macro = 0;  int nproc_micro = 0;
  for (int i = 0 ; i < nproc_world ; i++) {
    if (id_vec[i] == COLOR_MACRO)
      nproc_macro++;
    else if (id_vec[i] == COLOR_MICRO)
      nproc_micro++;
    else
      return 1;
  }

  if (nproc_micro != nproc_macro) return 1;

  if (comm->color == COLOR_MICRO) {

    int micro_position = 0;
    for (int i = 0 ; i < rank_world ; i++)
      if (id_vec[i] == COLOR_MICRO) micro_position++;
    comm->color += micro_position;

    int macro_count = 0; int macro_leader = 0;
    while (macro_leader < nproc_world) {

      if (id_vec[macro_leader] == COLOR_MACRO) {
	if (macro_count == micro_position) break;
	macro_count++;
      }
      macro_leader++;
    }
    comm->macro_leader = macro_leader;

  }else{

    int macro_position = 0;
    for (int i = 0 ; i < rank_world ; i++)
      if (id_vec[i] == COLOR_MACRO) macro_position++;

    int micro_count = 0; int micro_slave = 0;
    while (micro_slave < nproc_world) {

      if (id_vec[micro_slave] == COLOR_MICRO) {
	if (micro_count == macro_position) break;
	micro_count++;
      }
      micro_slave++;
    }
    comm->micro_slave = micro_slave;
  }

  ierr = MPI_Comm_split(WORLD_COMM, comm->color, 0, LOCAL_COMM);
  free(id_vec);

  return ierr;
}
