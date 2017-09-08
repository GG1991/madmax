/*
   Output routines
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
     Function for plotting the partition
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
  fprintf(vtkfl, "POINTS %d double\n", nallnods);

  for (n=0;n<nallnods;n++){
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
    fprintf(vtkfl, "%d\n",vtkcode(dim,eptr[e+1] - eptr[e]));  
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
int SpuVTKPlot_Displ_Strain_Stress(MPI_Comm PROBLEM_COMM, char *vtkfile_n, Vec *Displa, double *Strain, double *Stress)
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
  double  *xvalues;
  int     ierr;

  Vec xlocal;

  MPI_Comm_size(PROBLEM_COMM, &nproc);
  MPI_Comm_rank(PROBLEM_COMM, &rank);

//  xvalues = malloc(nallnods*3*sizeof(double));

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
  fprintf(vtkfl, "POINTS %d double\n", nallnods);

  for (n=0;n<nallnods;n++){
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
    fprintf(vtkfl, "%d\n",vtkcode(dim,eptr[e+1] - eptr[e]));  
  }

  fprintf(vtkfl, "POINT_DATA %i\n",nallnods);

  fprintf(vtkfl, "VECTORS Displa FLOAT\n");
//  fprintf(vtkfl, "LOOKUP_TABLE default\n");
  ierr = VecGetArray(xlocal, &xvalues); CHKERRQ(ierr);
  for (n=0;n<nallnods;n++){
    for (d=0;d<3;d++){
      fprintf(vtkfl, "%lf ", xvalues[n*3+d]);
    }
    fprintf(vtkfl,"\n");
  }
  VecRestoreArray(xlocal,&xvalues); CHKERRQ(ierr);

  fprintf(vtkfl, "CELL_DATA %i\n",nelm);

  fprintf(vtkfl, "TENSORS Strain FLOAT\n");
//  fprintf(vtkfl, "LOOKUP_TABLE default\n");
  for (e=0;e<nelm;e++){
    fprintf(vtkfl, "%lf %lf %lf\n", Strain[e*6+0],Strain[e*6+3],Strain[e*6+5]);
    fprintf(vtkfl, "%lf %lf %lf\n", Strain[e*6+3],Strain[e*6+1],Strain[e*6+4]);
    fprintf(vtkfl, "%lf %lf %lf\n", Strain[e*6+5],Strain[e*6+4],Strain[e*6+2]);
    fprintf(vtkfl, "\n"); 
  }

  fprintf(vtkfl, "TENSORS Stress FLOAT\n");
//  fprintf(vtkfl, "LOOKUP_TABLE default\n");
  for (e=0;e<nelm;e++){
    fprintf(vtkfl, "%lf %lf %lf\n", Stress[e*6+0],Stress[e*6+3],Stress[e*6+5]);
    fprintf(vtkfl, "%lf %lf %lf\n", Stress[e*6+3],Stress[e*6+1],Stress[e*6+4]);
    fprintf(vtkfl, "%lf %lf %lf\n", Stress[e*6+5],Stress[e*6+4],Stress[e*6+2]);
    fprintf(vtkfl, "\n"); 
  }

  fclose(vtkfl);

  return 0;
}   
/****************************************************************************************************/
int write_pvtu(MPI_Comm PROBLEM_COMM, char *name)
{
  /*
     Rank 0 writes this file
  */
  int  rank, nproc; 

  MPI_Comm_size(PROBLEM_COMM, &nproc);
  MPI_Comm_rank(PROBLEM_COMM, &rank);

  if(!rank){

    FILE *fm;
    char file_name[NBUF];
    int  i;

    strcpy(file_name,name);
    strcat(file_name,".pvtu");
    fm = fopen(file_name,"w"); if(!fm)SETERRQ1(PROBLEM_COMM,1,"file %s could not be opened",file_name);

    fprintf(fm, "<?xml version=\"1.0\"?>\n"
	"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
	"<PUnstructuredGrid GhostLevel=\"0\">\n"
	"<PPoints>\n"
	"<PDataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\"/>\n"
	"</PPoints>\n"
	"<PointData> Vectors=\"displ\"\n"
	"<PDataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" />"
	"</PointData>\n"
	"<PCells>\n"
	"<PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>\n"
	"</PCells>\n" 
	"<PCellData Scalars=\"Stress\">"
	"<PDataArray type=\"Int32\" Name=\"Material\" NumberOfComponents=\"1\"/>\n"
	"</PCellData>\n"); 
    for(i=0;i<nproc;i++){
      sprintf(file_name,"%s_%d",name,i);
      fprintf(fm,	"<Piece Source=\"%s.vtu\"/>\n",file_name);
    }
    fprintf(fm,	"</PUnstructuredGrid>\n" 
      "</VTKFile>\n"
      );

    fclose(fm);
  }
  return 0;
}
/****************************************************************************************************/
int write_vtu(MPI_Comm PROBLEM_COMM, char *name, Vec *x, Vec *b, double *strain, double *stress, double *energy)
{
  int  rank, nproc, ierr; 

  MPI_Comm_size(PROBLEM_COMM, &nproc);
  MPI_Comm_rank(PROBLEM_COMM, &rank);

  FILE *fm;
  char file_name[NBUF];
  int  i, d;
  Vec  xlocal;
  double *xvalues;
 
  /* 
     rank 0 writes the .pvtu file first
   */
  if(!rank){

    strcpy(file_name,name);
    strcat(file_name,".pvtu");
    fm = fopen(file_name,"w"); if(!fm)SETERRQ1(PROBLEM_COMM,1,"file %s could not be opened",file_name);

    fprintf(fm, "<?xml version=\"1.0\"?>\n"
	"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
	"<PUnstructuredGrid GhostLevel=\"0\">\n"
	"<PPoints>\n"
	"<PDataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\"/>\n"
	"</PPoints>\n"
	"<PCells>\n"
	"<PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>\n"
	"</PCells>\n" 

	"<PPointData Vectors=\"displ\">\n" // Vectors inside is a filter we should not use this here
	"<PDataArray type=\"Float64\" Name=\"displ\"    NumberOfComponents=\"3\" />\n"
	"<PDataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" />\n"
	"<PDataArray type=\"Int32\"   Name=\"bc_kinds\" NumberOfComponents=\"1\" />\n"
	"</PPointData>\n"

	"<PCellData>\n"
	"<PDataArray type=\"Int32\"   Name=\"part\" NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"9\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"9\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"energy\" NumberOfComponents=\"1\"/>\n"
	"</PCellData>\n"); 
    for(i=0;i<nproc;i++){
      sprintf(file_name,"%s_%d",name,i);
      fprintf(fm,	"<Piece Source=\"%s.vtu\"/>\n",file_name);
    }
    fprintf(fm,	"</PUnstructuredGrid>\n" 
      "</VTKFile>\n"
      );

    fclose(fm);
  }

  sprintf(file_name,"%s_%d.vtu",name,rank);
  fm = fopen(file_name,"w"); if(!fm)SETERRQ1(PROBLEM_COMM,1,"file %s could not be opened",file_name);

  fprintf(fm, 
      "<?xml version=\"1.0\"?>\n"
      "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      "<UnstructuredGrid>\n");
  fprintf(fm,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nallnods, nelm);
  fprintf(fm,"<Points>\n");
  fprintf(fm,"<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(i=0;i<nallnods;i++){
    for(d=0;d<dim;d++){
      fprintf(fm,"%e ",coord[i*3+d]);
    }
    for(d=dim;d<3;d++){
      fprintf(fm,"%e ",0.0);
    }
    fprintf(fm,"\n ");
  }
  fprintf(fm,"</DataArray>\n");
  fprintf(fm,"</Points>\n");
  fprintf(fm,"<Cells>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (i=0;i<nelm;i++){
    for (d=0;d<(eptr[i+1] - eptr[i]);d++){
      fprintf(fm,"%d ",eind[eptr[i] + d]);
    }
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (i=1;i<nelm+1;i++){
    fprintf(fm,"%d ",eptr[i]);
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (i=0;i<nelm;i++){
    fprintf(fm, "%d ",vtkcode(dim,eptr[i+1] - eptr[i]));  
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"</Cells>\n");
  
  fprintf(fm,"<PointData Vectors=\"displ\">\n"); // Vectors inside is a filter we should not use this here

  /*
     <displ>
   */
  ierr = VecGhostUpdateBegin(*x,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecGhostUpdateEnd(*x,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecGhostGetLocalForm(*x,&xlocal);
  fprintf(fm,"<DataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
  ierr = VecGetArray(xlocal, &xvalues); CHKERRQ(ierr);
  for(i=0;i<nallnods;i++){
    for(d=0;d<3;d++){
      fprintf(fm, "%lf ", xvalues[i*3+d]);
    }
    fprintf(fm,"\n");
  }
  VecRestoreArray(xlocal,&xvalues); CHKERRQ(ierr);
  fprintf(fm,"</DataArray>\n");

  /*
     <residual>
   */
  ierr = VecGhostUpdateBegin(*b,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecGhostUpdateEnd(*b,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecGhostGetLocalForm(*b,&xlocal);
  fprintf(fm,"<DataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
  ierr = VecGetArray(xlocal, &xvalues); CHKERRQ(ierr);
  for(i=0;i<nallnods;i++){
    for(d=0;d<dim;d++){
      fprintf(fm, "%lf ", xvalues[i*3+d]);
    }
    for(d=dim;d<3;d++){
      fprintf(fm, "%lf ", 0.0);
    }
    fprintf(fm,"\n");
  }
  VecRestoreArray(xlocal,&xvalues); CHKERRQ(ierr);
  fprintf(fm,"</DataArray>\n");

  /*
     <bc_kinds>
   */
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"bc_kinds\" NumberOfComponents=\"1\" format=\"ascii\" >\n");
  for(i=0;i<nmynods;i++){
    fprintf(fm, "%d ", bc_kinds[i]);
  }
  for(i=0;i<nghost;i++){
    fprintf(fm, "%d ", -1);
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");
  fprintf(fm,"</PointData>\n");

  fprintf(fm,"<CellData>\n");

  /*
     <part>
   */
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"part\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (i=0;i<nelm;i++){
    fprintf(fm,"%d ",rank);  
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  /*
     <strain>
   */
  if(dim==2){
    fprintf(fm,"<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"4\" format=\"ascii\">\n");
    for (i=0;i<nelm;i++){
      fprintf(fm, "%lf %lf ", strain[i*3+0],strain[i*3+2]);
      fprintf(fm, "%lf %lf ", strain[i*3+2],strain[i*3+1]);
    }
    fprintf(fm,"\n");
  }
  else if(dim==3){
    fprintf(fm,"<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"9\" format=\"ascii\">\n");
    for (i=0;i<nelm;i++){
      fprintf(fm, "%lf %lf %lf ", strain[i*6+0],strain[i*6+3],strain[i*6+5]);
      fprintf(fm, "%lf %lf %lf ", strain[i*6+3],strain[i*6+1],strain[i*6+4]);
      fprintf(fm, "%lf %lf %lf ", strain[i*6+5],strain[i*6+4],strain[i*6+2]);
    }
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  /*
     <stress>
   */
  if(dim==2){
    fprintf(fm,"<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"4\" format=\"ascii\">\n");
    for (i=0;i<nelm;i++){
      fprintf(fm, "%lf %lf ", stress[i*3+0],stress[i*3+2]);
      fprintf(fm, "%lf %lf ", stress[i*3+2],stress[i*3+1]);
    }
    fprintf(fm,"\n");
  }
  else if(dim==3){
    fprintf(fm,"<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"9\" format=\"ascii\">\n");
    for (i=0;i<nelm;i++){
      fprintf(fm, "%lf %lf %lf ", stress[i*6+0],stress[i*6+3],stress[i*6+5]);
      fprintf(fm, "%lf %lf %lf ", stress[i*6+3],stress[i*6+1],stress[i*6+4]);
      fprintf(fm, "%lf %lf %lf ", stress[i*6+5],stress[i*6+4],stress[i*6+2]);
    }
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  /*
     <stress>
   */
  fprintf(fm,"<DataArray type=\"Float64\" Name=\"energy\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (i=0;i<nelm;i++){
    fprintf(fm, "%lf ", energy[i]);
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,
      "</CellData>\n"
      "</Piece>\n"
      "</UnstructuredGrid>\n"
      "</VTKFile>\n");

  fclose(fm);
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
/****************************************************************************************************/
//  if(outfile!=NULL){
//    // queremos escribir algo
//    if(rank==0){
//      fprintf(outfile,"boundary_list\n");
//      pBound = boundary_list.head;
//      while(pBound){
//	fprintf(outfile,"name: %-8s NNod: %6d\n",
//	    ((boundary_t*)pBound->data)->name,((boundary_t*)pBound->data)->NNods);
//	for(i=0;i<((boundary_t*)pBound->data)->NNods;i++){
//	  fprintf(outfile,"%6d ", ((boundary_t*)pBound->data)->Nods[i]);
//	}
//	fprintf(outfile,"\n");
//	pBound=pBound->next;
//      }
//    }
//  }
