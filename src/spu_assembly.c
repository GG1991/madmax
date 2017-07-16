/*
 * Routines for matrix assembly
 *
 */

#include "sputnik.h"

int AssemblyJac(Mat *J)
{

  /*    Assembly the Jacobian for Small Deformation
   *    approach.
   *
   */

  int e, npe;
  int PETScIndex[8*3];

  double coord[8*3];

  for(e=0;e<nelm;e++){
    npe = eptr[e+1]-eptr[e];
    GetPETScIndeces( &eind[eptr[e]], npe, loc2petsc, PETScIndex);




  }
  return 0;

}

int GetPETScIndeces(int *LocalNod, int n, int *local2PETSc, int *PETScIndex)
{
   
    /* 
     * Gives the indeces to Gather fields and assembly on
     * PETSc matrices and vectors
     *
     * Input:
     *   LocalNod    :  array with local node numeration
     *   n           : number of nodes of the element
     *   local2PETSc : renumbering vetor to transform from local to PETSc
     *
     * Output:
     *    PETScIndex : Array with PETSc numeration
     *
     */

    int i, d;

    for(i=0;i<n;i++){
      for(d=0;d<3;d++){
	PETScIndex[i*3+d] = local2PETSc[LocalNod[i]]*3+d;
      }
    }
    return 0;
}
