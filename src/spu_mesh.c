/*

   SPUTNIK mesh treatment functions

*/

#include "sputnik.h"
#include "parmetis.h"

int part_mesh_PARMETIS(MPI_Comm *comm, FILE *time_fl, char *myname, int *elmdist, int *eptr, int *eind, int *part, double *centroid, int algorithm)
{

    /*
       Performes the mesh partition mesh saved on the mesh structures.

       a) First it builds the dual graph (nodes are elements)  of the 
          original (nodes are nodes)

       b) Do the partition 

       c) distribute the graph to processes

     */

    int        rank, nproc, i,j, ierr;
    int        nelm;             // number of elements of this process
    idx_t    * elmwgt;           // (inp) Element weights
    idx_t      wgtflag;          // (inp) Element weight flag (0 desactivated)
    idx_t      numflag;          // (inp) Numeration ( 0 in C, 1 in Fortran)
    idx_t      ncon;             // (inp) number of constrains of the graph ?
    idx_t      ncommonnodes;     // (inp) degree of connectivities among vertices in dual graph
    idx_t      nparts;           // (inp) number of partitions
    real_t   * tpwgts;           // (inp) array of size "ncon" x "npart" fraction of vertex for each subdomain
    real_t   * ubvec;            // (inp) array of size "ncon"
    idx_t      options[3];       // (inp) option parameters
    idx_t      edgecut;          // (out) number of edges cut of the partition
//    idx_t          * part;             // (out) Array for storing the solution
    double     t0_loc,t1_loc;

    MPI_Comm_size(*comm, &nproc);
    MPI_Comm_rank(*comm, &rank);

    nelm = elmdist[rank+1] - elmdist[rank];

    //**************************************************
    //
    // Set up some options
    //    
    elmwgt  = NULL; // no weights per elements
    wgtflag = 0;    // no weights per elements
    numflag = 0;    // C numeration
    
    nparts = nproc; // number of partitions 

    ncon = 1;
    tpwgts = (real_t*)malloc(ncon * nparts * sizeof(real_t));
    for(i=0; i < ncon * nparts ;i++){
      // uniform distribution of vertex in all processes
      tpwgts[i] = 1.0 / nparts;
    }
    
    ncommonnodes = 8;
    
    options[0] = 0; // options (1,2) : 0 default, 1 considered
    options[1] = 0; // level of information returned
    options[2] = 0; // random seed

    ubvec = (real_t*)malloc(ncon * sizeof(real_t));
    for(i=0;i<ncon;i++){
      ubvec[i] = 1.05;
    }

    //**************************************************

    if(algorithm == PARMETIS_GEOMKWAY){

    }
    else if(algorithm == PARMETIS_GEOM){

    }
    else if(algorithm == PARMETIS_KWAY){

    }
    else if(algorithm == PARMETIS_MESHKWAY){

      // Performe the partition with no weights
      t0_loc = MPI_Wtime();

      ParMETIS_V3_PartMeshKway (
	  elmdist, eptr, eind, elmwgt, &wgtflag, &numflag,
	  &ncon, &ncommonnodes, &nparts, tpwgts, ubvec,
	  options, &edgecut, part, comm );

      t1_loc = MPI_Wtime() - t0_loc;
      ierr = MPI_Gather(&t1_loc, 1, MPI_DOUBLE, time_vec, 1, MPI_DOUBLE, 0, *comm);
      if(ierr){
	return 1;
      }
      if(rank == 0){
	fprintf(time_fl,"%-20s", "mesh partition");
	for(i=0;i<nproc;i++){
	  fprintf(time_fl," %e",time_vec[i]);
	}
	fprintf(time_fl,"\n");
      }

      /* 
       * Graph distribution
       *
       * First we create an array "npe" that follows
       * the same function as eptr but is a global 
       * reference 
       *
       * eptr = [ 0 3 5 8 9 ]
       * npe  = [ 3 2 3 1 ]    (npe[i] = eptr[i+1] - eptr[i])
       *
       * Then vectors are switched acording to "part"
       * we create npe_swi, eind_swi, npe_swi_size
       *
       * We do MPI_Alltoall of : "npe_swi_size"  ->  "npe_size_new"
       *                         "eind_swi_size" ->  "eind_size_new"
       *
       * we free and reallocate memory for : "npe" using "npe_size" 
       *                                     "eind" using "eind_size" 
       *
       */

      t0_loc = MPI_Wtime();

      int *eind_swi, *eind_swi_size, *eind_size_new;
      int *npe_swi, *npe_swi_size, *npe_size_new;         
      int *npe;

      npe = malloc(nelm*sizeof(int));
      for(i=0;i<nelm;i++){
	npe[i] = eptr[i+1] - eptr[i];
      }

      eind_swi      = malloc(eptr[nelm]*sizeof(int)); 
      npe_swi       = malloc(nelm*sizeof(int)); 
      eind_swi_size = malloc(nproc*sizeof(int)); 
      npe_swi_size  = malloc(nproc*sizeof(int)); 
      eind_size_new = malloc(nproc*sizeof(int)); 
      npe_size_new  = malloc(nproc*sizeof(int)); 
      
      // swap npe and eind
      swap_vectors_SCR( part, nproc, nelm, npe, eptr, eind, npe_swi, eind_swi, npe_swi_size, eind_swi_size );

      printf("%-6s r%2d %-14s :", myname, rank, "npe_swi_size");
      for(i=0;i<nproc;i++){
	printf("%8d ",npe_swi_size[i]);
      }
      printf("\n");

      printf("%-6s r%2d %-14s :", myname, rank, "eind_swi_size");
      for(i=0;i<nproc;i++){
	printf("%8d ",eind_swi_size[i]);
      }
      printf("\n");

      // free & reallocate memory for "npe" & "eind"
      int npe_size_new_tot;
      int eind_size_new_tot;

      ierr = MPI_Alltoall(npe_swi_size, 1, MPI_INT, npe_size_new, 1, MPI_INT, *comm);
      if(ierr){
	return 1;
      }

      ierr = MPI_Alltoall(eind_swi_size, 1, MPI_INT, eind_size_new, 1, MPI_INT, *comm);
      if(ierr){
	return 1;
      }

      printf("%-6s r%2d %-14s :", myname, rank, "npe_size_new");
      for(i=0;i<nproc;i++){
	printf("%8d ",npe_size_new[i]);
      }
      printf("\n");

      printf("%-6s r%2d %-14s :", myname, rank, "eind_size_new");
      for(i=0;i<nproc;i++){
	printf("%8d ",eind_size_new[i]);
      }
      printf("\n");

      npe_size_new_tot = 0;
      for(i=0;i<nproc;i++){
	npe_size_new_tot += npe_size_new[i];
      }

      eind_size_new_tot = 0;
      for(i=0;i<nproc;i++){
	eind_size_new_tot += eind_size_new[i];
      }
      
      printf("%-6s r%2d %-14s : %8d\n", myname, rank, "npe_size_new_tot", npe_size_new_tot);
      printf("%-6s r%2d %-14s : %8d\n", myname, rank, "eind_size_new_tot", eind_size_new_tot);

      free(npe);
      free(eind);
      npe = malloc(npe_size_new_tot*sizeof(int));
      eind = malloc(eind_size_new_tot * sizeof(int));

      /* performe the MPI_Alltoall operation for calculating "npe" & "eind"
       *
       * for "npe"
       * sdispls = npe_swi_size
       * rdispls = npe_size_new
       *
       * for "eind"
       * sdispls = eind_swi_size
       * rdispls = eind_size_new
       */

      int *sdispls, *rdispls;

      sdispls = malloc(nproc*sizeof(int)); 
      rdispls = malloc(nproc*sizeof(int)); 

      for(i=0;i<nproc;i++){
	sdispls[i] = 0;
	for(j=0;j<i;j++){
	  sdispls[i] += npe_swi_size[j];
	}
      }
      for(i=0;i<nproc;i++){
	rdispls[i] = 0;
	for(j=0;j<i;j++){
	  rdispls[i] += npe_size_new[j];
	}
      }

      ierr = MPI_Alltoallv(npe_swi, npe_swi_size, sdispls, MPI_INT, 
	  npe, npe_size_new, rdispls, MPI_INT, *comm);

      for(i=0;i<nproc;i++){
	sdispls[i] = 0;
	for(j=0;j<i;j++){
	  sdispls[i] += eind_swi_size[j];
	}
      }
      for(i=0;i<nproc;i++){
	rdispls[i] = 0;
	for(j=0;j<i;j++){
	  rdispls[i] += eind_size_new[j];
	}
      }

      ierr = MPI_Alltoallv(eind_swi, eind_swi_size, sdispls, MPI_INT, 
	  eind, eind_size_new, rdispls, MPI_INT, *comm);

      free(eind_swi);
      free(npe_swi);

      t1_loc = MPI_Wtime() - t0_loc;
      ierr = MPI_Gather(&t1_loc, 1, MPI_DOUBLE, time_vec, 1, MPI_DOUBLE, 0, *comm);
      if(ierr){
	return 1;
      }
      if(rank == 0){
	fprintf(time_fl,"%-20s", "mesh distribution");
	for(i=0;i<nproc;i++){
	  fprintf(time_fl," %e",time_vec[i]);
	}
	fprintf(time_fl,"\n");
      }

    }
    else{

      return 1;
    }
    free(tpwgts);

    return 0;
}

int swap_vector( int *swap, int n, int *vector, int *new_vector, int *cuts )
{

  /* 
     swaps a vector

     example:

     swap       = [ 0 1 0 0 1 2 2 ]
     vector     = [ 0 1 2 3 4 5 6 ]
     n = 7

     new_vector = [ 0 2 3 1 4 5 6 ] 
     cut        = [ 3 2 2 ]

     Notes:

     -> if new_vector = NULL the result is saved on vector
     -> swap should have values in [0,n)
  */

  int *aux_vector;
  int i,p,aux,j;

  if(n==0){
    return 0;
  }
  else if(vector == NULL || cuts == NULL){
    return 1;
  }
 
  if(new_vector == NULL){
    aux_vector = vector;
  }
  else{
    aux_vector = new_vector;
  }

  j = 0;
  for(p=0;p<n;p++){
    cuts[p] = 0;
    for(i=0;i<n;i++){
      if(swap[i] == p){
        aux=vector[i];
	aux_vector[i] = vector[j];
	aux_vector[j] = aux;
	j ++;
	cuts[p] ++;
      }
    }
  }

  return 0;
}

int swap_vectors_SCR( int *swap, int nproc, int n,  int *npe, 
    int *eptr, int *eind, int *npe_swi, int *eind_swi, 
    int *cuts_npe, int *cuts_eind )
{

  /* 
     swaps a vectors in SCR format 

     example:

     swap = [ 0 2 1 0 ] (swap will be generally the "part" array)
     npe  = [ 3 2 3 1 ]             
     eind = [ 3 2 0 | 1 2 | 1 0 1 |3 ]   
     | 
     (swap operation with swap_vectors_CSR)
     | 
     npe  = [ 3 1 3 2 ]
     eind = [ 3 2 0 | 3 | 1 0 1 | 1 2 ]   

     Notes:

     -> "n" is the length of "npe"
     -> results are saved on "eind_swi" and "npe_swi" (memory is duplicated)
     -> swap should have values in [0,n)
  */

  int e, p, i, j, lp, pi;

  if(n==0){
    return 0;
  }
  if(!npe || !eind || !eind_swi || !npe_swi || !cuts_npe || !cuts_eind){
    return 1;
  }
  
  j = pi = lp = 0;
  for(p=0;p<nproc;p++){

    cuts_npe[p] = 0;
    for(e=0;e<n;e++){

      if(swap[e] == p){

	// swap npe
        npe_swi[j] = npe[e];
	j ++;

	// swap eind
	// CSR_give_pointer( e, npe, eind, &pi);
        pi = eptr[e];

	for(i=0;i<npe[e];i++){
	  eind_swi[lp] = eind[ pi + i ];
	  lp ++;
	}
	cuts_npe[p] ++;
      }
    }
  }

  for(i=0;i<nproc;i++){
    cuts_eind[i] = 0;
    for(j=0;j<cuts_npe[i];j++){
      if(i==0){
	cuts_eind[i] += npe_swi[j];
      }
      else{
	cuts_eind[i] += npe_swi[cuts_npe[i-1] + j];
      }
    }
  }

  return 0;
}

int CSR_give_pointer( int e, int *npe, int *eind, int *p)
{
  /*
   * returns the position "p" inside "eind" of element "e"
   *
   */

  int m;
  *p = 0;
  for(m=0;m<e;m++){
    *p += npe[m];
  }

  return 0;
}


int read_mesh(MPI_Comm * comm, char *myname, char *mesh_n, char *mesh_f, int ** elmdist, int ** eptr, int ** eind)
{

    /*

       Reads the mesh according to the format specified
       and performs the partition if it is required

       returns: 0 success
                1 not found
	       -1 failed

     */

    if(strcmp(mesh_f,"gmsh") == 0){
	return read_mesh_CSR_GMSH(comm, myname, mesh_n, elmdist, eptr, eind);
    }else{
      return 1;
    }
}

/****************************************************************************************************/

int read_mesh_CSR_GMSH(MPI_Comm * comm, char *myname, char *mesh_n, int ** elmdist, int ** eptr, int ** eind)
{

    /* 

       Info:   Reads the elements with the nodes conectivities and saves on 
               "elmdist[]", "eptr[]" and "eind[]" in CSR format (same names
	       that parmetis)

       Input: 
       char   * mesh_n   : file name with path
       MPI_Comm comm     : the communicator of these processes
       
       Output:
       int  ** elmdist  : number of elements for each process            (MAH)
       int  ** eptr     : array of indeces for "eind" (CSR format)       (MAH)
       int  ** eind     : element conectivities with nodes (CSR format)	 (MAH)


       1) first counts the total number of volumetric element on the mesh nelm_tot

       2) calculates nelm = nelm_tot/nproc (elements assigned to this process)
          calculates the vector elmdist in order to know how many elems will be for each process 

       3) read the mesh again, each process reads its own group of elements and see
          element types determines "npe" and fills "eptr[nelm+1]"
          finally alloc memory for "eind[eptr[nelm]]"

       4) reads the mesh again and fill "eind[]"

       Notes:

       a) rank and nproc are going to be respect to the communicator "comm"

       b) all processes do fopen and fread up to their corresponding position
          in the file

       Author: Guido Giuntoli

    */

    FILE               * fm;
    unsigned long int    offset;

    int                  nelm, nelm_tot;
    int                  npe;
    int                  total;
    int                  resto;
    int                  i, d, n; 
    int                  len;               // strlen(buf) for adding to offset
    int                  ln;                // line counter
    int                  ierr;
    int                  ntag;              // ntag to read gmsh element conectivities
    int                  rank;
    int                  nproc;
    
    char                 buf[NBUF];   
    char               * data;

    MPI_Comm_size(*comm, &nproc);
    MPI_Comm_rank(*comm, &rank);

    fm = fopen(mesh_n,"r");
    if(!fm){
	ierr = PetscPrintf(*comm,"file not found : %s\n",mesh_n);CHKERRQ(ierr);
	return 1;
    }

    /**************************************************/
    //  count the total number of volumetric elements 
    //  nelm_tot
    //
    offset   = 0;
    nelm_tot = 0;
    while(fgets(buf,NBUF,fm)!=NULL){
        offset += strlen(buf); 
	data=strtok(buf," \n");
	//
	// leemos hasta encontrar $Elements
	//
	if(strcmp(data,"$Elements")==0){
	    //
	    // leemos el numero total pero no lo usamos 
	    // (incluye elementos de superficie y de volumen)
	    //
	    fgets(buf,NBUF,fm);
	    offset += strlen(buf); 
	    data  = strtok(buf," \n");
	    total = atoi(data);

	    //
	    // leemos hasta $EndElements
	    // y contamos el numero total de los elementos volumen
	    //
	    for(i=0; i<total; i++){
	        fgets(buf,NBUF,fm); 
		len = strlen(buf);
		data=strtok(buf," \n");
		data=strtok(NULL," \n");
		if(atoi(data) == 4 || atoi(data) == 5 || atoi(data) == 6){
		  nelm_tot ++;
		}
		else{
		  // is a surface element so be continue counting
		  ln ++;
		  offset += len; 
		}
	    }
	    ierr = PetscPrintf(*comm,"%-6s %-8s   :  %d\n", myname, "nelm_tot", nelm_tot);
	    break;
	}
	ln ++;
    }
    //
    /**************************************************/

    /**************************************************/
    //  armamos el vector elmdist. 
    //  example: P0 tiene sus elementos entre 
    //  elmdist[0] a elemdist[1] (not included)
    //  los elementos que sobran a la division 
    //  nelm_tot/nproc los repartimos uno por 
    //  uno entre los primeros procesos
    //
    ierr = PetscPrintf(*comm,"%-6s %-8s   : ", myname, "elmdist");
    *elmdist = (int*)calloc( nproc + 1 ,sizeof(int));
    resto = nelm_tot % nproc;
    (*elmdist)[0] = 0;
    ierr = PetscPrintf(*comm,"%d ",(*elmdist)[0]);CHKERRQ(ierr);
    for(i=1; i < nproc + 1; i++){
	(*elmdist)[i] = i * nelm_tot / nproc;
        if(resto>0){
	  (*elmdist)[i] += 1;
	  resto --;
	}
	ierr = PetscPrintf(*comm,"%d ",(*elmdist)[i]);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(*comm,"\n");CHKERRQ(ierr);

    // ya podemos allocar el vector "eptr" su dimension es :
    // número de elementos locales + 1 = nelm + 1
    nelm = (*elmdist)[rank+1] - (*elmdist)[rank];
    *eptr = (int*)calloc( nelm + 1 ,sizeof(int));
    //
    /**************************************************/

    /**************************************************/
    //   
    // we read the file again (from offset) 
    // to count number of nodes 
    // per element and fill "eptr"
    // with this vector we can alloc memory for "eind"
    //    
    fseek( fm, offset, SEEK_SET);         // we go up to the first volumetric element
    for(i=0; i<(*elmdist)[rank]; i++){    // we go to the first element we have to store
      fgets(buf,NBUF,fm); 
      offset += strlen(buf); 
    }
    (*eptr)[0] = 0;
    for(i=1; i<nelm+1; i++){
      fgets(buf,NBUF,fm); 
      data=strtok(buf," \n");
      data=strtok(NULL," \n");
      switch(atoi(data)){
	case 4:
	  npe = 4;
	  break;
	case 5:
	  npe = 8;
	  break;
	case 6:
	  npe = 6;
	  break;
	default:
	  break;
      }
      (*eptr)[i] = (*eptr)[i-1] + npe; 
    }
    *eind = (int*)calloc( (*eptr)[nelm] ,sizeof(int));
    //
    /**************************************************/


    /**************************************************/
    //
    // repetimos el proceso pero esta vez leemos los 
    // nodos y completamos el vector "eind[eptr[nelm]]"
    // empezamos a leer desde "offset"
    //
    fseek( fm, offset, SEEK_SET);         // we go up to the first volumetric element
    n = 0;
    for(i=0; i < nelm ; i++){
      fgets(buf,NBUF,fm); 
      data=strtok(buf," \n");
      data=strtok(NULL," \n");
      switch(atoi(data)){
	case 4:
	  npe = 4;
	  break;
	case 5:
	  npe = 8;
	  break;
	case 6:
	  npe = 6;
	  break;
	default:
	  break;
      }
      data=strtok(NULL," \n");
      ntag = atoi(data);
      d = 0;
      while(d<ntag){
	data = strtok(NULL," \n");
	d++;
      }
      d = 0;
      while(d<npe){
	data = strtok(NULL," \n");
	(*eind)[n+d] = atoi(data); 
	d++;
      }
      n += npe;
    }
    //
    /**************************************************/

    return 0;   
}

/****************************************************************************************************/
