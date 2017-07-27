/*
 * VTK routines
 *
 */

#define   VTK_POINT         1
#define   VTK_LINE          3
#define   VTK_TRIANGLE      5
#define   VTK_QUADRANGLE    9
#define   VTK_TETRAHEDRON   10
#define   VTK_HEXAHEDRON    12
#define   VTK_6N_PRISM      13

#include "sputnik.h"

int spu_vtk_partition( char *vtkfile_n, MPI_Comm *comm )
{

  /* 
   * Function for plotting the partition
   */

  FILE    *vtkfl;

  int     rank, nproc; 
  int     count, d, n, e;

  MPI_Comm_size(*comm, &nproc);
  MPI_Comm_rank(*comm, &rank);

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
  fprintf(vtkfl, "POINTS %d double\n", NAllMyNod);

  for (n=0;n<NAllMyNod;n++){
    for(d=0;d<3;d++){
      fprintf(vtkfl, "%lf ", coord[n*3 + d]);
    }
    fprintf(vtkfl, "\n");
  }

  count=0;
  for(e=0;e<nelm;e++){
    count += eptr[e+1] - eptr[e] + 1;
  }
  fprintf(vtkfl, "CELLS %d %d\n", nelm, count);
  for (e=0;e<nelm;e++){
    fprintf(vtkfl, "%d ", eptr[e+1] - eptr[e]);
    for (n=0;n<(eptr[e+1] - eptr[e]);n++){
      fprintf(vtkfl, "%d ", eind[eptr[e] + n]);
    }
    fprintf(vtkfl, "\n");
  }

  fprintf(vtkfl, "CELL_TYPES %i\n", nelm);
  for (e=0;e<nelm;e++){
    fprintf(vtkfl, "%d\n",vtkcode(3,eptr[e+1] - eptr[e]));  
  }

  fprintf(vtkfl, "CELL_DATA %i\n",nelm);

  fprintf(vtkfl, "SCALARS part FLOAT\n");
  fprintf(vtkfl, "LOOKUP_TABLE default\n");
  for (e=0;e<nelm;e++){
    fprintf(vtkfl, "%lf\n",rank*1.0);
  }

  fclose(vtkfl);

  return 0;
}   

/****************************************************************************************************/

int SpuVTKPlot_Displ_Strain_Stress(MPI_Comm PROBLEM_COMM, char *vtkfile_n, Vec *Displa, Vec *Strain, Vec *Stress)
{

  /* 
     Plots in ASCII VTK > 

     Displa (On nodes),  
     Strain (On elements)
     Stress (On elements)

   */

  FILE    *vtkfl;

  int     rank, nproc; 
  int     count, d, n, e;
  double  index, value, *xvalues;
  int     ierr;

  Vec xlocal;

  MPI_Comm_size(PROBLEM_COMM, &nproc);
  MPI_Comm_rank(PROBLEM_COMM, &rank);

//  xvalues = malloc(NAllMyNod*3*sizeof(double));

  ierr = VecGhostUpdateBegin(*Displa,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecGhostUpdateEnd(*Displa,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecGhostGetLocalForm(*Displa,&xlocal);

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
  fprintf(vtkfl, "POINTS %d double\n", NAllMyNod);

  for (n=0;n<NAllMyNod;n++){
    for(d=0;d<3;d++){
      fprintf(vtkfl, "%lf ", coord[n*3 + d]);
    }
    fprintf(vtkfl, "\n");
  }

  count=0;
  for(e=0;e<nelm;e++){
    count += eptr[e+1] - eptr[e] + 1;
  }
  fprintf(vtkfl, "CELLS %d %d\n", nelm, count);
  for (e=0;e<nelm;e++){
    fprintf(vtkfl, "%d ", eptr[e+1] - eptr[e]);
    for (n=0;n<(eptr[e+1] - eptr[e]);n++){
      fprintf(vtkfl, "%d ", eind[eptr[e] + n]);
    }
    fprintf(vtkfl, "\n");
  }

  fprintf(vtkfl, "CELL_TYPES %i\n", nelm);
  for (e=0;e<nelm;e++){
    fprintf(vtkfl, "%d\n",vtkcode(3,eptr[e+1] - eptr[e]));  
  }

  fprintf(vtkfl, "POINT_DATA %i\n",NAllMyNod);

  fprintf(vtkfl, "VECTORS Displa FLOAT\n");
  fprintf(vtkfl, "LOOKUP_TABLE default\n");
  ierr = VecGetArray(xlocal, &xvalues); CHKERRQ(ierr);
  for (n=0;n<NAllMyNod;n++){
    for (d=0;d<3;d++){
//      index = n*3+d;
//      ierr = VecGetValues(xlocal, 1, &index, &value); CHKERRQ(ierr);
      fprintf(vtkfl, "%lf ", xvalues[n*3+d]);
    }
    fprintf(vtkfl,"\n");
  }
  VecRestoreArray(xlocal,&xvalues); CHKERRQ(ierr);

//  fprintf(vtkfl, "CELL_DATA %i\n",nelm);
//
//  fprintf(vtkfl, "SCALARS Strain FLOAT\n");
//  fprintf(vtkfl, "LOOKUP_TABLE default\n");
//  for (e=0;e<nelm;e++){
//    fprintf(vtkfl, "%lf\n",rank*1.0);
//  }

  fclose(vtkfl);

  return 0;
}   

/****************************************************************************************************/

int vtkcode(int dim,int npe)
{

  switch(dim){
    case 1:
      switch(npe){
        case 2 :
          return VTK_LINE;
        default:
          return -1;
      }
    case 2:
      switch(npe){
        case 3 :
          return VTK_TRIANGLE;
        case 4 :
          return VTK_QUADRANGLE;
        default:
          return -1;
      }
    case 3:
      switch(npe){
        case 4 :
          return VTK_TETRAHEDRON;
        case 6 :
          return VTK_6N_PRISM;  
        case 8 :
          return VTK_HEXAHEDRON;  
        default:
          return -1;
      }
    default:
      return -1;
  }
}
