#include "macro.h"


int macro_pvtu(char *name) {

  FILE    *fm;
  char    file_name[NBUF];
  double  *xvalues;
  Vec     xlocal;

  if (rank_mac == 0) {

    strcpy(file_name,name);
    strcat(file_name,".pvtu");
    fm = fopen(file_name,"w");

    fprintf(fm, "<?xml version=\"1.0\"?>\n"
	"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
	"<PUnstructuredGrid GhostLevel=\"0\">\n"
	"<PPoints>\n"
	"<PDataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\"/>\n"
	"</PPoints>\n"
	"<PCells>\n"
	"<PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>\n"
	"</PCells>\n");

    fprintf(fm, "<PPointData Vectors=\"displ\">\n");
    if (x != NULL)
    fprintf(fm, "<PDataArray type=\"Float64\" Name=\"displ\"    NumberOfComponents=\"3\" />\n");
    if (b != NULL)
      fprintf(fm, "<PDataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" />\n");

    fprintf(fm, "</PPointData>\n"
	"<PCellData>\n"
	"<PDataArray type=\"Int32\"   Name=\"part\"   NumberOfComponents=\"1\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\"/>\n"
	"<PDataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\"/>\n"
	"<PDataArray type=\"Int32\"   Name=\"elem_type\" NumberOfComponents=\"1\"/>\n"
	"</PCellData>\n" , nvoi , nvoi);

    for (int i = 0 ; i < nproc_mac ; i++) {
      sprintf(file_name,"%s_%d", name, i);
      fprintf(fm, "<Piece Source=\"%s.vtu\"/>\n", file_name );
    }
    fprintf(fm,	"</PUnstructuredGrid>\n</VTKFile>\n" );

    fclose(fm);
  }

  sprintf( file_name, "%s_%d.vtu", name, rank_mac);
  fm = fopen(file_name,"w");
  if (!fm) {
    myio_printf(PETSC_COMM_WORLD,"Problem trying to opening file %s for writing\n", file_name);
    return 1;
  }

  fprintf(fm,
      "<?xml version=\"1.0\"?>\n"
      "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      "<UnstructuredGrid>\n");
  fprintf(fm,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", mesh.nnods_local_ghost, mesh.nelm_local);
  fprintf(fm,"<Points>\n");
  fprintf(fm,"<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");

  for (int n = 0 ; n < mesh.nnods_local_ghost ; n++) {
    for (int d = 0 ; d < dim ; d++)
      fprintf(fm,"% 01.6e ",  mesh.coord_local[n*dim + d]);
    for (int d = dim ; d < 3 ; d++)
      fprintf(fm, "% 01.6e ", 0.0 );
    fprintf(fm, "\n" );
  }
  fprintf(fm,"</DataArray>\n");
  fprintf(fm,"</Points>\n");
  fprintf(fm,"<Cells>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (int e = 0 ; e < mesh.nelm_local ; e++) {
    int npe = mesh.eptr[e+1] - mesh.eptr[e];
    for (int n = 0 ; n < npe ; n++)
      fprintf(fm,"%-6d ", mesh.eind[mesh.eptr[e] + n]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  int ce = 0.0;
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (int e = 0 ; e < mesh.nelm_local ; e++) {
    int npe = mesh.eptr[e+1] - mesh.eptr[e];
    ce += npe;
    fprintf(fm, "%d ", ce);
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (int e = 0 ; e < mesh.nelm_local ; e++) {
    int npe = mesh.eptr[e+1] - mesh.eptr[e];
    fprintf(fm, "%-3d ", vtkcode(dim , npe));
  }
  fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"</Cells>\n");

  fprintf(fm,"<PointData Vectors=\"displ\">\n"); // Vectors inside is a filter we should not use this here

  if (x != NULL) {
    VecGhostUpdateBegin(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostGetLocalForm(x, &xlocal);

    fprintf(fm,"<DataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
    VecGetArray( xlocal , &xvalues );
    for (int n = 0 ; n < mesh.nnods_local_ghost ; n++) {
      for (int d = 0 ; d < dim ; d++)
	fprintf(fm, "% 01.6e ", xvalues[n*dim + d]);
      for (int d = dim ; d < 3 ; d++)
	fprintf(fm,"% 01.6e ",0.0);
      fprintf(fm,"\n");
    }
    VecRestoreArray(xlocal , &xvalues);
    fprintf(fm,"</DataArray>\n");
  }

  if (b != NULL) {
    VecGhostUpdateBegin(b , INSERT_VALUES,SCATTER_FORWARD);
    VecGhostUpdateEnd(b , INSERT_VALUES,SCATTER_FORWARD);
    VecGhostGetLocalForm(b , &xlocal);

    fprintf(fm,"<DataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
    VecGetArray(xlocal, &xvalues);
    for (int n = 0 ; n < mesh.nnods_local_ghost ; n++) {
      for (int d = 0 ; d < 3 ; d++) fprintf(fm, "% 01.6e ", (d < dim) ? xvalues[n*dim + d] : 0);
      fprintf(fm,"\n");
    }
    VecRestoreArray(xlocal, &xvalues);
    fprintf(fm,"</DataArray>\n");
  }
  fprintf(fm,"</PointData>\n<CellData>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"part\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (int e = 0; e < mesh.nelm_local ; e++) fprintf(fm, "%d ", rank_mac ); fprintf(fm, "\n");
  fprintf(fm, "</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for (int e = 0 ; e < mesh.nelm_local ; e++) {
    for (int v = 0 ; v < nvoi ; v++) fprintf(fm, "% 01.6e ", elem_strain[e*nvoi + v]); fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for (int e = 0; e < mesh.nelm_local ; e++) {
    for (int v = 0 ; v < nvoi ; v++) fprintf(fm, "% 01.6e ", elem_stress[e*nvoi + v]); fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"elem_type\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (int e = 0; e < mesh.nelm_local ; e++) fprintf(fm, "%d ", elem_type[e]); fprintf(fm,"\n");
  fprintf(fm,"</DataArray>\n");

  fprintf(fm, "</CellData>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");

  fclose(fm);
  return 0;
}
