/*
   Routine to performe coloring on peuge processes
   
   Strategy the first "nmacro" process will correspond to this peuge

   
*/

#include "mpi.h"
#include "monts.h"    // montseny common feature such as  variables to performe coloring
#include "peuge.h"

int peu_colors(MPI_Comm world, MPI_Comm intercomm)
{

    // we are going to work with a "local" communicator saved in "world"
    // we create a newone called localworld
    MPI_Comm_size(world, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &grank);
    MPI_Comm_rank(world, &lrank);

    if(lrank == 0)
	printf("peuge:grank = %d - nproc = %d",grank,nproc);

    // we performe the splitting (it should be made on all the processes)
    color = MACRO;
    MPI_Comm_split(world, color, 0, &localworld);

    return 0;
}
