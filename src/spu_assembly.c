/*
 * Routines for matrix assembly
 *
 */

#include "sputnik.h"

/****************************************************************************************************/

int AssemblyJac(Mat *J)
{

  /*    Assembly the Jacobian for Small Deformation
   *    approach.
   *
   */

  int e, gp, ngp, npe;
  int PETScIndex[8*3];

  double ElemCoord[8][3];
  double ShapeDerivs[8][3];
  double DetJac;

  for(e=0;e<nelm;e++){
    npe = eptr[e+1]-eptr[e];
    ngp = npe;
    GetPETScIndeces( &eind[eptr[e]], npe, loc2petsc, PETScIndex);
    GetElemCoord(&eind[eptr[e]], npe, ElemCoord);

    for(gp=0;gp<ngp;gp++){

      GetShapeDerivs(gp, npe, ElemCoord, ShapeDerivs, &DetJac);

    }




  }
  return 0;

}

/****************************************************************************************************/
int GetShapeDerivs(int gp, int npe, double coor[8][3], double ShapeDerivs[8][3], double *DetJac)
{

  double ***ShapeDerivsMaster;

  double jac[3][3];
  double ijac[3][3];

  ShapeDerivsMaster = FemGetPointer2ShapeDerivsMaster(npe, 3);
  if(ShapeDerivsMaster == NULL) return 1;

  FemCalculateJac3D( coor, ShapeDerivsMaster, npe, gp, jac);
  FemInvertJac3D( jac, ijac, DetJac);
  FemGiveShapeDerivs( ijac, npe, gp, ShapeDerivsMaster, ShapeDerivs);



  return 0;
}

/****************************************************************************************************/

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

/****************************************************************************************************/

int GetElemCoord(int *LocalNod, int n, double ElemCoord[8][3])
{
   
  /*
   *  Returns the coordinates of the vertex of that element
   */

  int i, d;

  for(i=0;i<n;i++){
    for(d=0;d<3;d++){
      ElemCoord[i][d] = coord[ LocalNod[i]*3 + d ];
    }
  }

  return 0;
}
