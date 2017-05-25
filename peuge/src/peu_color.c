/*
   Routine to performe coloring on peuge processes

   MPI_COMM_WORLD
         |
	/ \
   MACRO   MICRO  (colors with code on monts.h)

   Strategy the first "nmacro" process will correspond to this peuge, example : 

   MACRO_1 --> MICRO_K1_1 -> MICRO_K1_2 -> MICRO_K1_3    (i=0 - sum_micro_proc = 0 - lrank = 0)
	   \   
            -> MICRO_K2_1 -> MICRO_K2_2                  (i=1 - sum_micro_proc = 3 - lrank = 0)

   MACRO_2 --> MICRO_K1_1 -> MICRO_K1_2 -> MICRO_K1_3    (i=0 - sum_micro_proc = 0 - lrank = 1)
	   \   
            -> MICRO_K2_1 -> MICRO_K2_2                  (i=1 - sum_micro_proc = 3 - lrank = 1)
   
*/

#include "mpi.h"
#include "monts.h"    // montseny common feature such as  variables to performe coloring
#include "stdlib.h"
#include "stdio.h"

int peu_colors(MPI_Comm world, MPI_Comm * macro, MPI_Comm * macmic_comm, 
	int   nproc_mac, 
	int   nkind_mic, 
	int * color,
	int * lrank,
	int * grank,
	int * mpinfo_mic)
{

    int  i, nproc;
    int  sum_micro_proc;     // partial sum of microprocess per macroscopic process (vary with "i")
    int  sum_micro_proc_tot; // total sum of microprocess per macroscopic process
    int  grank_remote;

    // we are going to work with a "local" communicator saved in "world"
    // we create a newone called localworld
    MPI_Comm_size(world, &nproc);
    if(nproc != nproc_mac)
	return 1;

    MPI_Comm_rank(MPI_COMM_WORLD, grank);
    MPI_Comm_rank(world, lrank);

    // we need one intercomm per kind_mic
    macmic_comm = (MPI_Comm*)malloc(nkind_mic*sizeof(MPI_Comm));

    if(lrank == 0)
	printf("peuge:lrank = %d - grank = %d - nproc = %d",*lrank,*grank,nproc);

    // we performe the splitting (it should be made on all the processes)
    *color = MACRO;
    MPI_Comm_split(world, *color, 0, macro);

    // compute total quantity of process per macroscopic processes
    sum_micro_proc_tot = 0;
    for(i=0;i<nkind_mic;i++){
	sum_micro_proc_tot += mpinfo_mic[i*3+1];
    }

    // we have to create the intercommunicators
    sum_micro_proc = 0;
    for( i=0 ; i < nkind_mic ; i++ ){
        // the global rank formula aim to the first microscopic process that should be in 
	// that remote communicator for perform the interconection
	sum_micro_proc += mpinfo_mic[i*3+1]; 
	grank_remote = nproc_mac + (*lrank) * sum_micro_proc_tot + sum_micro_proc;
	MPI_Intercomm_create(*macro, 0, world, grank_remote, 0, &macmic_comm[i]);
    }

    return 0;
}
