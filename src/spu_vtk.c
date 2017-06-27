/*
 * VTK format file routines
 *
 *
 */

#include "sputnik.h"

int spu_vtk_partition( char *myname, char *mesh_n, MPI_Comm *comm )
{

  /* This function is for plotting a partition
   * rank = 0 is the responsible of writing the vtk
   * file
   *
   */

  int     rank, nproc, i, j, ierr;

  MPI_Comm_size(*comm, &nproc);
  MPI_Comm_rank(*comm, &rank);

  if(rank==0){ 

    FILE    *vtkfl, *meshfl;
    int     nnod;
    double  vald;
    char    buf[NBUF], filevtk[32], ending[5], *data;

    strcpy(filevtk,myname);
    filevtk[strlen(myname)]='\0';
    sprintf(ending,".vtk");
    strcat(filevtk,ending);

    vtkfl = fopen(filevtk,"w");
    if(!vtkfl){
      return 1;
    }
    meshfl = fopen(mesh_n,"r");
    if(!meshfl){
      return 1;
    }

    /************************************************************/
    /*Mesh geometry data                                        */
    /************************************************************/    
    fprintf(vtkfl, "# vtk DataFile Version 2.0\n");
    fprintf(vtkfl, "SPUTNIK\n");
    fprintf(vtkfl, "ASCII\n");
    fprintf(vtkfl, "DATASET UNSTRUCTURED_GRID\n");

    while(fgets(buf,NBUF,meshfl)!=NULL){
      data=strtok(buf," \n");
      //
      // leemos hasta encontrar $Nodes
      //
      if(strcmp(data,"$Nodes")==0){
	//
	// leemos el numero total de nodos 
	//
	fgets(buf,NBUF,meshfl);
	data  = strtok(buf," \n");
	nnod = atoi(data);
	//
	// leemos hasta $EndNodes
	//
	for(i=0; i<nnod; i++){
	  fgets(buf,NBUF,meshfl); 
	  data=strtok(buf," \n");
	  for(j=0;j<3;j++){
	    data=strtok(NULL," \n");
	    fprintf(vtkfl," %e",atof(data)); 
	  }
	  fprintf(vtkfl,"\n");
	}
	break;
      }
    }
    fclose(meshfl);
    fprintf(vtkfl, "POINTS %d double\n", nnod);

    fclose(vtkfl);
  }
  else{

  }

  return 0;
}   
