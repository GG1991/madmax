#include "micro.h"


int micro_print_info(void){

  FILE *fm = fopen("micro_info.dat","w");
  int *i_data;

  if(rank_mic == 0){
    fprintf(fm,"-----------\n");
    fprintf(fm,"nproc %d\n", nproc_mic);
    fprintf(fm,"-----------\n");
    fprintf(fm,"nx    ny    nz\n");
    fprintf(fm,"%2d    %2d    %2d\n", nx, ny, nz);
    fprintf(fm,"lx    ly    lz\n");
    fprintf(fm,"%lf    %lf    %lf\n", lx, ly, lz);
    fprintf(fm,"hx    hy    hz\n");
    fprintf(fm,"%lf    %lf    %lf\n", hx, hy, hz);
    fprintf(fm,"nex   ney   nez\n");
    fprintf(fm,"%2d    %2d    %2d\n", nex, ny-1, nez);
    fprintf(fm,"-----------\n");
  }

  if(rank_mic == 0){
    i_data = malloc( nproc_mic * sizeof(double));
    fprintf(fm,"%-20s","rank ");
    for(int i = 0 ; i < nproc_mic ; i++)
      fprintf(fm,"%d ", i);
    fprintf(fm,"\n");
    fprintf(fm,"%-20s","");
    for(int i = 0 ; i < nproc_mic ; i++)
      fprintf(fm,"--");
    fprintf(fm,"\n");
  }

  int ierr = MPI_Gather(&nelm, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return 1;

  if(rank_mic == 0){
    fprintf(fm,"%-20s","nelm");
    for(int i = 0 ; i < nproc_mic ; i++)
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&ney, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return 1;

  if(rank_mic == 0){
    fprintf(fm,"%-20s","eyl");
    for(int i = 0 ; i < nproc_mic ; i++)
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&nyl, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr != 0) return 1;

  if(rank_mic == 0){
    fprintf(fm,"%-20s","nyl");
    for(int i = 0 ; i < nproc_mic ; i++)
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&ny_inf, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return 1;

  if(rank_mic == 0){
    fprintf(fm,"%-20s","ny_inf");
    for(int i = 0 ; i < nproc_mic ; i++)
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  ierr = MPI_Gather(&ngho, 1, MPI_INT, (!rank_mic)?i_data:NULL, 1, MPI_INT, 0, MICRO_COMM);
  if(ierr) return 1;

  if(rank_mic == 0){
    fprintf(fm,"%-20s","ngho");
    for(int i = 0 ; i < nproc_mic ; i++)
      fprintf(fm,"%d ", i_data[i]);
    fprintf(fm,"\n");
  }

  double norm;
 
  if(b != NULL){
    VecNorm( b, NORM_2, &norm );
    if( rank_mic == 0 )
      fprintf(fm,"|b| = %lf \n", norm);
  }

  if(x != NULL){
    VecNorm(x, NORM_2, &norm);
    if(rank_mic == 0)
      fprintf(fm,"|x| = %lf \n", norm);
  }

  if(A != NULL){
    MatNorm(A, NORM_FROBENIUS, &norm);
    if(rank_mic == 0)
      fprintf(fm,"|A| = %lf \n",norm);
  }

  fclose(fm);
  return 0;
}

int micro_pvtu(char *name){

  FILE    *fm;
  char     file_name[NBUF];
  double  *xvalues;
  Vec      xlocal;

  if(rank_mic == 0){

    strcpy(file_name,name);
    strcat(file_name,".pvtu");
    fm = fopen(file_name,"w");

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

	"<PPointData Vectors=\"displ\">\n" 
	"<PDataArray type=\"Float64\" Name=\"displ\"    NumberOfComponents=\"3\" />\n"
	"<PDataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" />\n"
	"</PPointData>\n"

	"<PCellData>\n"
	"<PDataArray type=\"Int32\"   Name=\"part\"   NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\"/>\n"
	"<PDataArray type=\"Int32\"   Name=\"elem_type\" NumberOfComponents=\"1\"/>\n"
	"</PCellData>\n" , nvoi , nvoi);

    for(int i = 0 ; i < nproc_mic ; i++ ){
      sprintf(file_name,"%s_%d",name,i);
      fprintf(fm,	"<Piece Source=\"%s.vtu\"/>\n",file_name);
    }
    fprintf(fm,	"</PUnstructuredGrid>\n" 
      "</VTKFile>\n" );

    fclose(fm);

  }

  sprintf(file_name,"%s_%d.vtu",name,rank_mic);
  fm = fopen(file_name,"w"); 
  if(!fm){
    myio_printf(MICRO_COMM,"Problem trying to opening file %s for writing\n", file_name);
    return 1;
  }

  fprintf(fm, 
      "<?xml version=\"1.0\"?>\n"
      "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      "<UnstructuredGrid>\n");
  fprintf(fm,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", nl + ngho, nelm);
  fprintf(fm,"<Points>\n");
  fprintf(fm,"<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");

  double *coord = malloc( dim * sizeof(double));

  for(int n = 0 ; n < nl ; n++){
    get_node_local_coor( n , coord);
    for(int d = 0 ; d < dim ; d++)
      fprintf(fm,"%e ",  coord[d]);
    for(int d = dim ; d < 3 ; d++)
      fprintf(fm,"%e ",0.0);
    fprintf(fm,"\n");
  }
  for(int n = 0 ; n < ngho ; n++ ){
    get_node_ghost_coor( n , coord );
    for(int d = 0 ; d < dim ; d++ )
      fprintf(fm,"%e ",  coord[d] );
    for(int d = dim ; d < 3 ; d++ )
      fprintf(fm,"%e ",0.0);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");
  fprintf(fm,"</Points>\n");
  fprintf(fm,"<Cells>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int e = 0 ; e < nelm ; e++){
    get_local_elem_node( e , loc_elem_index );
    for(int n = 0 ; n < npe ; n++)
      fprintf(fm,"%d ", loc_elem_index[n]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  int ce = npe;
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int e = 0 ; e < nelm ; e++){
    fprintf(fm,"%d ", ce);
    ce += npe;
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int e = 0 ; e < nelm ; e++)
    fprintf(fm, "%d ",vtkcode(dim, npe));  
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"</Cells>\n");
  
  fprintf(fm,"<PointData Vectors=\"displ\">\n"); // Vectors inside is a filter we should not use this here

  VecGhostUpdateBegin(x , INSERT_VALUES, SCATTER_FORWARD);
  VecGhostUpdateEnd(x , INSERT_VALUES, SCATTER_FORWARD);
  VecGhostGetLocalForm(x , &xlocal);

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
  VecGetArray( xlocal , &xvalues );
  for(int n = 0 ; n < (nl + ngho) ; n++){
    for(int d = 0 ; d < dim ; d++)
      fprintf(fm, "%lf ", xvalues[ n * dim + d ]);
    for(int d = dim ; d < 3 ; d++)
      fprintf(fm,"%lf ",0.0);
    fprintf(fm,"\n");
  }
  VecRestoreArray( xlocal , &xvalues );
  fprintf(fm,"</DataArray>\n");

  VecGhostUpdateBegin( b , INSERT_VALUES, SCATTER_FORWARD );
  VecGhostUpdateEnd  ( b , INSERT_VALUES, SCATTER_FORWARD );
  VecGhostGetLocalForm(b , &xlocal);

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
  VecGetArray(xlocal, &xvalues);
  for(int n = 0 ; n < (nl + ngho) ; n++){
    for(int d = 0 ; d < dim ; d++)
      fprintf(fm, "%lf ", xvalues[ n * dim + d ]);
    for(int d = dim ; d < 3 ; d++)
      fprintf(fm, "%lf ", 0.0);
    fprintf(fm,"\n");
  }
  VecRestoreArray( xlocal , &xvalues );
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"</PointData>\n");
  fprintf(fm,"<CellData>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"part\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int e = 0; e < nelm ; e++ )
    fprintf( fm, "%d ", rank_mic );  
  fprintf( fm, "\n");
  fprintf( fm, "</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for(int e = 0; e < nelm ; e++){
    for(int v = 0 ; v < nvoi ; v++)
      fprintf(fm, "%lf ", elem_strain[ e*nvoi + v ]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for(int e = 0; e < nelm ; e++){
    for(int v = 0 ; v < nvoi ; v++)
      fprintf(fm, "%lf ", elem_stress[ e*nvoi + v ]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"elem_type\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int e = 0; e < nelm ; e++)
    fprintf(fm, "%d ", elem_type[e]);
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,
      "</CellData>\n""</Piece>\n"
      "</UnstructuredGrid>\n"
      "</VTKFile>\n");

  fclose(fm);
  return 0;
}


