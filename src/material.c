#include "material.h"

int mat_get_stress( material_t *mat_p, int dim , double * strain_gp, double * stress_gp )
{

  if( mat_p->type_id == TYPE_0 ){

    /* is a linear material stress = C * strain */

    double  young   = ((type_0*)mat_p->type)->young;
    double  poisson = ((type_0*)mat_p->type)->poisson;
    int     i , j;

    if( dim == 2 ){

      /* Plane strain ( long fibers case ) */
      int     nvoi = 3;
      double  c[3][3];
      c[0][0]=1.0-poisson; c[0][1]=poisson    ; c[0][2]=0.0            ;
      c[1][0]=poisson    ; c[1][1]=1.0-poisson; c[1][2]=0.0            ;
      c[2][0]=0.0        ; c[2][1]=0.0        ; c[2][2]=(1-2*poisson)/2;

      for( i = 0; i < nvoi ; i++ ){
	for( j = 0 ; j < nvoi ; j++ )
	  c[i][j] *= young/((1+poisson)*(1-2*poisson));
      }
      for( i = 0; i < nvoi ; i++ ){
	stress_gp[i] = 0.0;
	for( j = 0 ; j < nvoi ; j++ )
	  stress_gp[i] += c[i][j] * strain_gp[j];
      }

    }
  }
  return 0;
}

/****************************************************************************************************/

int mat_get_c_tang( material_t *mat_p, int dim , double * strain_gp, double * c_tan )
{

  if( mat_p->type_id == TYPE_0 ){

    /* is a linear material stress = C * strain */
    double  young   = ((type_0*)mat_p->type)->young;
    double  poisson = ((type_0*)mat_p->type)->poisson;
    int     i , j;
    int     nvoi = 3;

    if(dim==2){

      /* Plane strain ( long fibers case ) */
      c_tan[0*nvoi+0]=1.0-poisson; c_tan[0*nvoi+1]=poisson    ; c_tan[0*nvoi+2]=0.0            ;
      c_tan[1*nvoi+0]=poisson    ; c_tan[1*nvoi+1]=1.0-poisson; c_tan[1*nvoi+2]=0.0            ;
      c_tan[2*nvoi+0]=0.0        ; c_tan[2*nvoi+1]=0.0        ; c_tan[2*nvoi+2]=(1-2*poisson)/2;

      for( i = 0; i < nvoi ; i++ ){
	for( j = 0 ; j < nvoi ; j++ )
	  c_tan[i*nvoi + j] *= young/((1+poisson)*(1-2*poisson));
      }

    }
  }
  return 0;
}

/****************************************************************************************************/

int mat_get_rho( material_t *mat_p, int dim , double * rho )
{

  if( mat_p->type_id == TYPE_0 )
    *rho   = ((type_0*)mat_p->type)->rho;

  return 0;
}

/****************************************************************************************************/
