/*
   Routines for performing homogenization on RVE

   Author > Guido Giuntoli
   Date   > 18-08-2017

 */

#include "micro.h"

int micro_homogenize_taylor(MPI_Comm PROBLEM_COMM, double strain_mac[6], double strain_ave[6], double stress_ave[6])
{

  int    i, k, e, gp, ngp, npe, ierr;

  double ElemCoord[8][3];
  double ShapeDerivs[8][3];
  double DetJac;
  double stress_gp[6];
  double DsDe[6][6];
  double *wp = NULL;
  double ElemDispls[8*3];
  double vol = -1.0, vol_tot = -1.0;
  double stress_aux[6], strain_aux[6];
  register double wp_eff;

  for(i=0;i<6;i++) strain_aux[i] = stress_aux[i] = 0.0;
 
  vol = 0.0;

  for(e=0;e<nelm;e++){

    npe = eptr[e+1]-eptr[e];
    ngp = npe;

    GetElemCoord(&eind[eptr[e]], npe, ElemCoord);
    GetWeight(npe, &wp);

    // calculate <ElemResidue> by numerical integration

    for(gp=0;gp<ngp;gp++){

      memset(stress_gp,0.0,6*sizeof(double));

      GetShapeDerivs(gp, npe, ElemCoord, ShapeDerivs, &DetJac);
      GetDsDe( e, ElemDispls, DsDe );
      wp_eff = DetJac*wp[gp];

      for(i=0;i<6;i++){
	for(k=0;k<6;k++){
	  stress_gp[i] += DsDe[i][k]*strain_mac[k];
	}
      }

      for(i=0;i<6;i++){
	stress_aux[i] += stress_gp[i] * wp_eff;
	strain_aux[i] += strain_mac[i] * wp_eff;
      }

      vol += wp_eff;
    }

  }

  ierr = MPI_Allreduce(stress_aux, stress_ave, 6, MPI_DOUBLE, MPI_SUM, PROBLEM_COMM);CHKERRQ(ierr);
  ierr = MPI_Allreduce(strain_aux, strain_ave, 6, MPI_DOUBLE, MPI_SUM, PROBLEM_COMM);CHKERRQ(ierr);
  ierr = MPI_Allreduce(&vol, &vol_tot, 1, MPI_DOUBLE, MPI_SUM, PROBLEM_COMM);CHKERRQ(ierr);

  for(i=0;i<6;i++){
    stress_ave[i] /= vol_tot;
    strain_ave[i] /= vol_tot;
  }

  return 0;
}
/****************************************************************************************************/
