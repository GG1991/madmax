/*

   SPUTNIK mesh treatment functions

*/

#include "sputnik.h"
#include "parmetis.h"

int part_mesh_PARMETIS(MPI_Comm * comm, int * elmdist, int * eptr, int * eind, double * centroid, int algorithm)
{

    /*
       Performes the mesh partition mesh saved on the mesh structures.

       a) First it builds the dual graph (nodes are elements)  of the 
          original (nodes are nodes)

       b) Do the partition 

     */
    int                  rank, nproc, i;
    int                  nelm;             // number of elements of this process
    idx_t              * elmwgt;           // (inp) Element weights
    idx_t                wgtflag;          // (inp) Element weight flag (0 desactivated)
    idx_t                numflag;          // (inp) Numeration ( 0 in C, 1 in Fortran)
    idx_t                ncon;             // (inp) number of constrains of the graph ?
    idx_t                ncommonnodes;     // (inp) degree of connectivities among vertices in dual graph
    idx_t                nparts;           // (inp) number of partitions
    real_t             * tpwgts;           // (inp) array of size "ncon" x "npart" fraction of vertex for each subdomain
    real_t             * ubvec;            // (inp) array of size "ncon"
    idx_t                options[3];       // (inp) option parameters
    idx_t                edgecut;          // (out) number of edges cut of the partition
    idx_t              * part;             // (out) Array for storing the solution

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
    
    ncommonnodes = 2;
    
    options[0] = 0; // options (1,2) : 0 default, 1 considered
    options[1] = 0; // level of information returned
    options[2] = 0; // random seed

    part = (idx_t*)malloc(nelm * sizeof(idx_t));

    ubvec = (real_t*)malloc(ncon * sizeof(real_t));
    for(i=0;i<ncon;i++){
      ubvec[i] = 1.05;
    }
    //
    //**************************************************

    if(algorithm == PARMETIS_GEOMKWAY){
    //  ParMETIS_V32_Mesh2Dual (
    //      idx t *elmdist, idx t *eptr, idx t *eind, idx t *numflag, idx t *ncommonnodes,
    //      idx t **xadj, idx t **adjncy, MPI Comm *comm
    //      )

    }
    else if(algorithm == PARMETIS_GEOM){

    }
    else if(algorithm == PARMETIS_KWAY){

    }
    else if(algorithm == PARMETIS_MESHKWAY){

      // Performe the partition with no weights
      ParMETIS_V3_PartMeshKway (
	  elmdist, eptr, eind, elmwgt, &wgtflag, &numflag,
	  &ncon, &ncommonnodes, &nparts, tpwgts, ubvec,
	  options, &edgecut, part, comm );

    }
    else{

      return 1;
    }

    return 0;
}

int read_mesh(MPI_Comm * comm, char *mesh_n, char *mesh_f, int ** elmdist, int ** eptr, int ** eind)
{

    /*

       Reads the mesh according to the format specified
       and performs the partition if it is required

     */

    if(strcmp(mesh_f,"gmsh") == 0){
	if(read_mesh_CSR_GMSH(comm, mesh_n, elmdist, eptr, eind))
	    return 1;
    }

    return 0;
}

/****************************************************************************************************/

int read_mesh_CSR_GMSH(MPI_Comm * comm, char *mesh_n, int ** elmdist, int ** eptr, int ** eind)
{

    /* 

       Author: Guido Giuntoli
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
    
    char                 buf[BUF_N_LENGTH];   
    char               * data;

    MPI_Comm_size(*comm, &nproc);
    MPI_Comm_rank(*comm, &rank);

    fm = fopen(mesh_n,"r");
    if(!fm){
	ierr = PetscPrintf(PETSC_COMM_WORLD,"file not found : %s\n",mesh_n);CHKERRQ(ierr);
	return 1;
    }

    /**************************************************/
    //  count the total number of volumetric elements 
    //  nelm_tot
    //
    offset   = 0;
    nelm_tot = 0;
    while(fgets(buf,BUF_N_LENGTH,fm)!=NULL){
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
	    fgets(buf,BUF_N_LENGTH,fm);
	    offset += strlen(buf); 
	    data  = strtok(buf," \n");
	    total = atoi(data);

	    //
	    // leemos hasta $EndElements
	    // y contamos el numero total de los elementos volumen
	    //
	    for(i=0; i<total; i++){
	        fgets(buf,BUF_N_LENGTH,fm); 
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
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"nelm_tot    : %d\n",nelm_tot);CHKERRQ(ierr);
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
    ierr = PetscPrintf(PETSC_COMM_WORLD,"elmdist     : ");CHKERRQ(ierr);
    *elmdist = (int*)calloc( nproc + 1 ,sizeof(int));
    resto = nelm_tot % nproc;
    (*elmdist)[0] = 0;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"%d ",(*elmdist)[0]);CHKERRQ(ierr);
    for(i=1; i < nproc + 1; i++){
	(*elmdist)[i] = i * nelm_tot / nproc;
        if(resto>0){
	  (*elmdist)[i] += 1;
	  resto --;
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"%d ",(*elmdist)[i]);CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);

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
      fgets(buf,BUF_N_LENGTH,fm); 
      offset += strlen(buf); 
    }
    (*eptr)[0] = 0;
    for(i=1; i<nelm+1; i++){
      fgets(buf,BUF_N_LENGTH,fm); 
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
      fgets(buf,BUF_N_LENGTH,fm); 
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
