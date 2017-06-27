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
  int     length;

  MPI_Comm_size(*comm, &nproc);
  MPI_Comm_rank(*comm, &rank);

  if(rank==0){ 

    FILE    *vtkfl, *meshfl;
    int     nnod, *length_vec;
    int     *eptr_a, *eind_a;
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

    length_vec = malloc(nproc*sizeof(int));

    length = nelm + 1;
    ierr = MPI_Gather(&length, 1, MPI_INT, length_vec, 1, MPI_INT, 0, *comm);
    if(ierr){
      return 1;
    }

    // from each process we receive a their mesh and we write on the file
    for(i=1;i<nproc;i++){ 
      eptr_a = malloc((length_vec[i])*sizeof(int));
      ierr = MPI_Recv(eptr_a, length_vec[i], MPI_INT, i, 0, *comm, &status);
      if(ierr){
	return 1;
      }


      free(eptr_a);
    }

    length = eptr[nelm];
    ierr = MPI_Gather(&length, 1, MPI_INT, length_vec, 1, MPI_INT, 0, *comm);
    if(ierr){
      return 1;
    }

    // from each process we receive a their mesh and we write on the file
    for(i=1;i<nproc;i++){ 
      eind_a = malloc((length_vec[i])*sizeof(int));
      ierr = MPI_Recv(eind_a, length_vec[i], MPI_INT, i, 0, *comm, &status);
      if(ierr){
	return 1;
      }


      free(eind_a);
    }

    fclose(vtkfl);
  }
  else{

    // I'm not rank 0 so I should send my mesh to
    // him and wait

     
    length = nelm + 1;
    ierr = MPI_Gather(&length, 1, MPI_INT, NULL, 1, MPI_INT, 0, *comm);
    MPI_Ssend(eptr,length,MPI_INT,0,0,*comm);

    length = eptr[nelm];
    ierr = MPI_Gather(&length, 1, MPI_INT, NULL, 1, MPI_INT, 0, *comm);
    MPI_Ssend(eind,length,MPI_INT,0,0,*comm);

  }

  return 0;
}   
