/*
 * VTK routines
 *
 */

#include "sputnik.h"

int spu_vtk_partition( char *vtkfile_n, MPI_Comm *comm )
{

  /* 
   * Function for plotting the partition
   */

  FILE    *vtkfl;

  int     rank, nproc, i, j, ierr;
  int     length;
  int     nnod, *length_vec;
  int     *eptr_a, *eind_a;
  char    buf[NBUF], *data;

  MPI_Comm_size(*comm, &nproc);
  MPI_Comm_rank(*comm, &rank);


//  strcpy(filevtk,myname);
//  filevtk[strlen(myname)]='\0';
//  sprintf(ending,".vtk");
//  strcat(filevtk,ending);

  vtkfl = fopen(vtkfile_n,"w");
  if(!vtkfl){
    return 1;
  }

  /************************************************************/
  /*Mesh geometry data                                        */
  /************************************************************/    
  fprintf(vtkfl, "# vtk DataFile Version 2.0\n");
  fprintf(vtkfl, "SPUTNIK\n");
  fprintf(vtkfl, "ASCII\n");
  fprintf(vtkfl, "DATASET UNSTRUCTURED_GRID\n");

  fclose(vtkfl);

  return 0;
}   
