#ifndef VTK_H_
#define VTK_H_

#define   VTK_POINT         1
#define   VTK_LINE          3
#define   VTK_TRIANGLE      5
#define   VTK_QUADRANGLE    9
#define   VTK_TETRAHEDRON   10
#define   VTK_HEXAHEDRON    12
#define   VTK_6N_PRISM      13

int vtkcode( int dim,int npe );

#endif
