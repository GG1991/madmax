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

  int    i, j, k, e, gp, ngp, npe;
  int    PETScIdx[8*3];
  int    ierr;

  double ElemCoord[8][3];
  double ShapeDerivs[8][3];
  double DetJac;
  double ElemMatrix[8*3 * 8*3];
  double B[6][3*8], Baux[6][3*8];
  double DsDe[6][6];
  double *wp = NULL;
  double ElemDispls[8*3];

  ierr = MatZeroEntries(*J);CHKERRQ(ierr);

  for(e=0;e<nelm;e++){
    npe = eptr[e+1]-eptr[e];
    ngp = npe;
    GetPETScIndeces( &eind[eptr[e]], npe, loc2petsc, PETScIdx);
    GetElemCoord(&eind[eptr[e]], npe, ElemCoord);
    GetWeight(npe, wp);

    // calculate <elemMat> by numerical integration

    for(gp=0;gp<ngp;gp++){

      GetShapeDerivs(gp, npe, ElemCoord, ShapeDerivs, &DetJac);
      GetB( npe, ShapeDerivs, B );
      GetDsDe( npe, ElemDispls, DsDe );

      for(i=0;i<6;i++){
	for(j=0;j<npe*3;j++){
	  Baux[i][j]=0.0;
	  for(k=0;k<6;k++){
	    Baux[i][j] += DsDe[i][k]*B[k][j];
	  }
	}
      }

      for(i=0;i<npe*3;i++){
	for(j=0;j<npe*3;j++){
	  for(k=0;k<3;k++){
	    ElemMatrix[i*npe*3+j] += B[k][i]*Baux[k][j] * DetJac * wp[gp];
	  }
	}
      }


    }
    ierr = MatSetValues(*J, npe*3, PETScIdx, npe*3, PETScIdx, ElemMatrix, ADD_VALUES);CHKERRQ(ierr);

  }
  ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  return 0;

}

/****************************************************************************************************/

int GetDsDe( int e, double *ElemDisp, double DsDe[6][6] )
{

  /*  Calculates constitutive tensor
   *  according to the element type
   */

  int npe, type;

  npe = eptr[e+1]-eptr[e];

  //  type = GetMaterialType(PhysicalID[e]);

  switch(type){

    case 0:
      /* 
       * Elástico lineal 
       */

      break;

    default:
      break;

  }

  return 0;
}

/****************************************************************************************************/

int GetWeight(int npe, double *wp)
{

  wp = FemGetPointer2Weight(npe, 3);
  return 0;
}

/****************************************************************************************************/

int GetB( int npe, double ShapeDerivs[8][3], double B[6][3*8] )
{

  int i;

  for(i=0;i<npe;i++){

    B[0][i*3+0] = ShapeDerivs[i][0]; 
    B[0][i*3+1] = 0.0         ;
    B[0][i*3+2] = 0.0         ; 

    B[1][i*3+0] = 0.0         ;
    B[1][i*3+1] = ShapeDerivs[i][1];
    B[1][i*3+2] = 0.0         ; 

    B[2][i*3+0] = 0.0         ;
    B[2][i*3+1] = 0.0         ;
    B[2][i*3+2] = ShapeDerivs[i][2]; 

    B[3][i*3+0] = ShapeDerivs[i][1];
    B[3][i*3+1] = ShapeDerivs[i][0];
    B[3][i*3+2] = 0.0         ; 

    B[4][i*3+0] = 0.0         ;
    B[4][i*3+1] = ShapeDerivs[i][2];
    B[4][i*3+2] = ShapeDerivs[i][1];

    B[5][i*3+0] = ShapeDerivs[i][2];
    B[5][i*3+1] = 0.0         ;
    B[5][i*3+2] = ShapeDerivs[i][0]; 

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
