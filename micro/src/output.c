#include "micro.h"


int micro_print_info(void)
{
  FILE *fm = fopen("micro_info.dat","w");

  fprintf(fm,"-----------\n");
  fprintf(fm,"nproc %d\n", nproc_mic);
  fprintf(fm,"-----------\n");
  fprintf(fm,"nx    ny    nz\n");
  fprintf(fm,"%2d    %2d    %2d\n",mesh_struct.nx, mesh_struct.ny, mesh_struct.nz);
  fprintf(fm,"lx    ly    lz\n");
  fprintf(fm,"%lf    %lf    %lf\n",mesh_struct.lx, mesh_struct.ly, mesh_struct.lz);
  fprintf(fm,"hx    hy    hz\n");
  fprintf(fm,"%lf    %lf    %lf\n",mesh_struct.hx, mesh_struct.hy, mesh_struct.hz);
  fprintf(fm,"nex   ney   nez\n");
  fprintf(fm,"%2d    %2d    %2d\n",mesh_struct.nex, mesh_struct.ney, mesh_struct.nez);
  fprintf(fm,"-----------\n");


  double norm;

  if (b != NULL) {
    VecNorm(b, NORM_2, &norm);
    fprintf(fm,"|b| = %lf \n", norm);
  }

  if (x != NULL) {
    VecNorm(x, NORM_2, &norm);
    fprintf(fm,"|x| = %lf \n", norm);
  }

  if (A != NULL) {
    MatNorm(A, NORM_FROBENIUS, &norm);
    fprintf(fm,"|A| = %lf \n",norm);
  }

  fclose(fm);
  return 0;
}

int micro_pvtu(char *name)
{
  FILE *fm;
  char file_name[NBUF];
  double *xvalues;

  strcpy(file_name, name);
  strcat(file_name, ".pvtu");
  fm = fopen(file_name, "w");

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

  sprintf(file_name,"%s_0",name);
  fprintf(fm, "<Piece Source=\"%s.vtu\"/>\n",file_name);
  fprintf(fm, "</PUnstructuredGrid>\n</VTKFile>\n");

  fclose(fm);

  sprintf(file_name, "%s_%d.vtu", name, rank_mic);
  fm = fopen(file_name, "w");
  if (!fm) {
    myio_printf(MICRO_COMM,"Problem trying to opening file %s for writing\n", file_name);
    return 1;
  }

  fprintf(fm,
      "<?xml version=\"1.0\"?>\n"
      "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      "<UnstructuredGrid>\n");
  fprintf(fm,"<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", mesh_struct.nn, mesh_struct.nelm);
  fprintf(fm,"<Points>\n");
  fprintf(fm,"<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n");

  double *coord = malloc(dim*sizeof(double));
  int npe = mesh_struct.npe;
  int dim = mesh_struct.npe;

  for (int n = 0; n < mesh_struct.nn; n++) {
    mesh_struct_get_node_coord(&mesh_struct, n, coord);
    for (int d = 0 ; d < 3 ; d++)
      fprintf(fm,"%e ",(d < mesh_struct.dim) ? coord[d] : 0.0);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n</Points>\n<Cells>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (int e = 0; e < mesh_struct.nelm; e++) {
    mesh_struct_get_elem_nods(&mesh_struct, e, elem_nods);
    for (int n = 0; n < npe; n++)
      fprintf(fm,"%d ", elem_nods[n]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  int ce = npe;
  fprintf(fm,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (int e = 0 ; e < mesh_struct.nelm ; e++) {
    fprintf(fm,"%d ", ce);
    ce += npe;
  }
  fprintf(fm,"\n</DataArray>\n");

  fprintf(fm,"<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (int e = 0 ; e < mesh_struct.nelm ; e++)
    fprintf(fm, "%d ", vtkcode(dim, npe));
  fprintf(fm,"\n</DataArray>\n");

  fprintf(fm,"</Cells>\n");

  fprintf(fm,"<PointData Vectors=\"displ\">\n"); // Vectors inside is a filter we should not use this here

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"displ\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
  VecGetArray(x, &xvalues);
  for (int n = 0 ; n < mesh_struct.nn ; n++) {
    for (int d = 0 ; d < 3 ; d++)
      fprintf(fm,"%e ",(d < dim) ? xvalues[n*dim + d] : 0.0);
    fprintf(fm,"\n");
  }
  VecRestoreArray(x, &xvalues);
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"residual\" NumberOfComponents=\"3\" format=\"ascii\" >\n");
  VecGetArray(b, &xvalues);
  for (int n = 0; n < mesh_struct.nn; n++) {
    for (int d = 0; d < 3; d++)
      fprintf(fm,"%e ",(d < dim) ? xvalues[n*dim + d] : 0.0);
    fprintf(fm,"\n");
  }
  VecRestoreArray(x, &xvalues);
  fprintf(fm,"</DataArray>\n</PointData>\n<CellData>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"part\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (int e = 0; e < mesh_struct.nelm ; e++)
    fprintf(fm, "%d ", rank_mic);
  fprintf( fm, "\n</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"strain\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for (int e = 0; e < mesh_struct.nelm; e++) {
    for (int v = 0; v < nvoi; v++)
      fprintf(fm, "%lf ", elem_strain[e*nvoi + v]);
    fprintf(fm, "\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Float64\" Name=\"stress\" NumberOfComponents=\"%d\" format=\"ascii\">\n",nvoi);
  for (int e = 0; e < mesh_struct.nelm; e++) {
    for (int v = 0; v < nvoi; v++)
      fprintf(fm, "%lf ", elem_stress[e*nvoi + v]);
    fprintf(fm,"\n");
  }
  fprintf(fm,"</DataArray>\n");

  fprintf(fm,"<DataArray type=\"Int32\" Name=\"elem_type\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for (int e = 0; e < mesh_struct.nelm ; e++)
    fprintf(fm, "%d ", elem_type[e]);
  fprintf(fm,"\n</DataArray>\n");

  fprintf(fm, "</CellData>\n""</Piece>\n</UnstructuredGrid>\n</VTKFile>\n");

  fclose(fm);
  return 0;
}
