
/*
   Giant mesh partition routines
 */


#include "parmetis.h"
#include "giant.h"


int gia_part(void)
{

    /*
       Performes the partition of the mesh saved on the mesh structures.
     */

    //  ParMETIS_V32_Mesh2Dual (
    //      idx t *elmdist, idx t *eptr, idx t *eind, idx t *numflag, idx t *ncommonnodes,
    //      idx t **xadj, idx t **adjncy, MPI Comm *comm
    //      )

    return 0;
}

int gia_rmsh(char *mesh_n, char *mesh_f)
{
    /*
       Reads the mesh according to the format specified
       and performs the partition if it is required
     */

    

    if(strcmp(mesh_f,"gmsh") == 0){
	if(gia_readev_gmsh(mesh_n))
	    return 1;
    }

    return 0;
}

int gia_readev_gmsh(char *mesh_n)
{

    /* 
       Reads the elements with the nodes conectivities and saves on list_elemv
       1) first counts the total number of volumetric element on the mesh n_elem_tot
       2) calculates n_elem = n_elem_tot / nproc
       3) calculates the vector vtxdist in order to know how many elem will be for each process 
       4) reads its own group of elements
    */

    FILE   * fm;
    int      n_elem_tot;
    int      total;
    int      i;
    
    char     buf[128];   
    char   * data;

    fm = fopen(mesh_n,"r");
    if(!fm){
	ierr = PetscPrintf(PETSC_COMM_WORLD,"file not found : %s\n",mesh_n);CHKERRQ(ierr);
	return 1;
    }

    n_elem_tot = 0;

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
	    // y contamos el numero total de los de volumen
	    //
	    for(i=0; i<total; i++){
	        fgets(buf,128,fm); 
		data=strtok(buf," \n");
		data=strtok(NULL," \n");
		if(atoi(data) == 4 || atoi(data) == 5 || atoi(data) == 6) n_elem_tot ++;
	    }
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"n_elem_tot  : %d\n",n_elem_tot);CHKERRQ(ierr);
	}
    }

    return 0;   
}
