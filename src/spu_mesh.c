/*

   SPUTNIK mesh treatment functions

*/

#include "sputnik.h"

int part_mesh(void)
{

    /*
       Performes the mesh partition mesh saved on the mesh structures.
     */

    //  ParMETIS_V32_Mesh2Dual (
    //      idx t *elmdist, idx t *eptr, idx t *eind, idx t *numflag, idx t *ncommonnodes,
    //      idx t **xadj, idx t **adjncy, MPI Comm *comm
    //      )

    return 0;
}

int read_mesh(char *mesh_n, char *mesh_f)
{
    /*
       Reads the mesh according to the format specified
       and performs the partition if it is required
     */

    if(strcmp(mesh_f,"gmsh") == 0){
	if(readev_gmsh(mesh_n))
	    return 1;
    }

    return 0;
}

/****************************************************************************************************/

int read_mesh_CSR_GMSH(char *mesh_n)
{

    /* 
       Reads the elements with the nodes conectivities and saves on "eptr[]" and "eind[]" in CSR format

       1) first counts the total number of volumetric element on the mesh nelm_tot

       2) calculates nelm = nelm_tot/nproc
          calculates the vector elmdist in order to know how many elems will be for each process 

       3) read the mesh again and seeing element types, so "npe", is it possible to fill "eptr[nelm+1]"
          finally alloc memory for "eind[eptr[nelm]]"

       4) reads the mesh again and fill "eind[]"
    */

    FILE   * fm;

    int      nelm, nelm_tot;
    int      nnod_tot;
    int      npe;
    int      total;
    int      resto;
    int      i, d, ei; 
    int      ln;               // line counter
    int      ntve;             // number of lines up to the first volume element
    int      ierr;
    int      nproc;
    
    char     buf[128];   
    char   * data;

    fm = fopen(mesh_n,"r");
    if(!fm){
	ierr = PetscPrintf(PETSC_COMM_WORLD,"file not found : %s\n",mesh_n);CHKERRQ(ierr);
	return 1;
    }


    /**************************************************/
    //  count the total number of volumetric elements 
    //  nelm_tot
    //
    ln = 0;
    nelm_tot = 0;
    while(fgets(buf,128,fm)!=NULL){
	data=strtok(buf," \n");
	//
	// leemos hasta encontrar $Elements
	//
	if(strcmp(data,"$Elements")==0){
	    //
	    // leemos el numero total pero no lo usamos 
	    // (incluye elementos de superficie y de volumen)
	    //
	    fgets(buf,128,fm);
	    data  = strtok(buf," \n");
	    total = atoi(data);

	    //
	    // leemos hasta $EndElements
	    // y contamos el numero total de los elementos volumen
	    //
	    for(i=0; i<total; i++){
	        fgets(buf,128,fm); 
		data=strtok(buf," \n");
		data=strtok(NULL," \n");
		if(atoi(data) == 4 || atoi(data) == 5 || atoi(data) == 6){
		  nelm_tot ++;
		}
		else{
		  // is a surface element so be continue counting
		  ln ++;
		}
	    }
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"nelm_tot  : %d\n",nelm_tot);CHKERRQ(ierr);
	    break;
	}
	ln ++;
    }
    rewind(fm);
    ntev = ln;
    //
    /**************************************************/

    /**************************************************/
    //  armamos el vector elmdist. 
    //  example: P0 tiene sus elementos entre vtxdist[0] 
    //  a vtxdist[1] (not included)
    //  los elementos que sobran a la division 
    //  nelm_tot/nproc los repartimos uno por 
    //  uno entre los primeros procesos
    //
    ierr = PetscPrintf(PETSC_COMM_WORLD,"vtxdist     : ");CHKERRQ(ierr);
    elmdist = (int*)calloc( nproc + 1 ,sizeof(int));
    resto = nelm_tot % nproc;
    d = 1;
    for(i=0; i < nproc + 1; i++){
	elmdist[i] = i * nelm_tot / nproc + d - 1;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"%d ",elmdist[i]);CHKERRQ(ierr);
	if(i == resto - 1) d = 0;
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);

    // ya podemos allocar el vector "eptr" su dimension es :
    // número de elementos locales + 1 = nelm + 1
    nelm = elmdist[rank+1] - elmdist[rank];
    eptr = (int*)calloc( nelm + 1 ,sizeof(int));
    //
    /**************************************************/

    /**************************************************/
    //   
    // we read the file again to count number of nodes 
    // per element and fill "eptr"
    // with this vector we can alloc memory for "eind"
    //    
    ln = 0;
    while(fgets(buf,128,fm)!=NULL){
        
	//
	// leemos hasta "ntve" 
	//
	if(ln < ntve){

	    ei = 0;
	    for(i=0; i<nelmtot; i++){
	        fgets(buf,128,fm); 
		data=strtok(buf," \n");
		data=strtok(NULL," \n");
		isvol = 1;
		switch(atoi(data)){
		  case 4:
		    npe = 4;
		    break;
		  case 5:
		    npe = 6;
		    break;
		  case 6:
		    npe = 8;
		    break;
		  default:
		    isvol = 0;
		    break;
		}
		if( isvol == 1 ){
		  if(i >= elmdist[rank] && i < elmdist[rank+1]){
		    if( i != elmdist[rank] ){
		      eptr[i] = eptr[i-1] + npe; 
		    }else{
		      //  i = 0
		      eptr[i] = 0;
		    }
		    ei ++;
		  }
		}
	    }
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"nelm_tot  : %d\n",nelm_tot);CHKERRQ(ierr);
	    break;
	}
	ln ++;
    }
    rewind(fm);
    eind = (int*)calloc( eptr[nelm] ,sizeof(int));


    //
    /**************************************************/

    

    /**************************************************/
    //
    // repetimos el proceso pero esta vez leemos los 
    // nodos y completamos los vectores 
    //
    while(fgets(buf,128,fm)!=NULL){
	data=strtok(buf," \n");
	//
	// leemos hasta encontrar $Elements
	//
	if(strcmp(data,"$Elements")==0){
	    //
	    // leemos el numero total pero no lo usamos (incluye elementos de superficie
	    // y de volumen
	    //
	    fgets(buf,128,fm);
	    data  = strtok(buf," \n");
	    total = atoi(data);

	    //
	    // leemos hasta $EndElements
	    // y contamos el número total de los de volumen
	    //
	    for(i=0; i<total; i++){
	        fgets(buf,128,fm); 
		data=strtok(buf," \n");
		data=strtok(NULL," \n");
		if(atoi(data) == 4 || atoi(data) == 5 || atoi(data) == 6) nelm_tot ++;
	    }
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"nelm_tot  : %d\n",nelm_tot);CHKERRQ(ierr);
	    break;
	}
    }
    rewind(fm);
    //
    /**************************************************/

    return 0;   
}

/****************************************************************************************************/
