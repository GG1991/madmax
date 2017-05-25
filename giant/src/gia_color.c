/*
   Routine to performe coloring on peuge processes

   MPI_COMM_WORLD
          |
	 / \
   MACRO     MICRO  (colors with code on monts.h)
   |         |    
   |_color 0 |_color 3
   |_color 1 |_color 4
   |_color 2 |_color 5
             |_color 6 
	     |_color 7
             |_color 8

   Strategy the first "nmacro" process will correspond to this peuge, example : 

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

#include "mpi.h"
#include "monts.h"    // montseny common feature such as  variables to performe coloring
#include "stdlib.h"
#include "stdio.h"

int gia_color(MPI_Comm world, MPI_Comm * micro, MPI_Comm * micmac_comm, 
	int   nproc_mac, 
	int   nkind_mic, 
	int * color,
	int * lrank,
	int * grank,
	int * mpinfo_mic)
{

    int  i, nproc; 
    int  a, b, c, d, e, f, g; // variables for calculating the color in an "easy" way
    int  sum_micro_proc;      // partial sum of microprocess per macroscopic process (vary with "i")
    int  sum_micro_proc_tot;  // total sum of microprocess per macroscopic process
    int  grank_remote;

    // we are going to work with a "local" communicator saved in "world"
    // we create a newone called localworld
    MPI_Comm_size(world, &nproc);
    MPI_Comm_rank(world, rank);

    // we calculate the color correspoding to this process
    sum_micro_proc_tot = 0;
    for( i=0 ; i < nkind_mic ; i++ ){
	sum_micro_proc_tot += mpinfo_mic[i*2]; 
    }
    a = rank - nproc_mac;
    c = sum_micro_proc_tot;
    b = a / c;
    d = a % c;
    f = 0;

    sum_micro_proc = 0;
    for( i=0 ; i < nkind_mic ; i++ ){
	if( d < sum_micro_proc ){
	    break;
	}
	sum_micro_proc += mpinfo_mic[i*2]; 
	f++;
    }
    e = nproc_mac + b + f ;

    *color = e;

    // we performe the splitting (it should be made on all the processes)
    MPI_Comm_split(world, *color, 0, micro);

    // verify that total number of process is valid for the communicators structure file
    if(nproc != sum_micro_proc_tot * nproc_mac + nproc_mac){
	return 1;
    }

    // create the intercommunicators
    // that remote communicator for perform the interconection
    grank_remote = (e - nproc_mac) % nproc_mac;
    MPI_Intercomm_create(*micro, 0, world, grank_remote, 0, &micmac_comm);

    printf("giant: rank = %d - nproc = %d - color = %d - rem_rank(to macro) %d\n",rank,nproc,e,grank_remote);

    return 0;
}
