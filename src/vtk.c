#include "vtk.h"

/****************************************************************************************************/

int vtkcode( int dim, int npe )
{

  switch(dim) {
    case 1:
      switch(npe) {
        case 2 :
          return VTK_LINE;
        default:
          return -1;
      }
    case 2:
      switch(npe) {
        case 3 :
          return VTK_TRIANGLE;
        case 4 :
          return VTK_QUADRANGLE;
        default:
          return -1;
      }
    case 3:
      switch(npe) {
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
