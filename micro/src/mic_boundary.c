/*
   Routines to impose Boundary condition on RVE cells

   Author > Guido Giuntoli
   Date > 31-07-2017
*/

#include "micro.h"

int MicroSetDisplacementOnBoundary( int dir, double strain_dir, double LX, double LY, double LZ, Vec *x )
{
  /*  
      Dirichlet Boundary condition set is set on <x> 
      usamos VecSetValuesLocal aqui ya que vamos a modificar valores locales unicamente
   */
  int ierr;
  switch(dir){
    case 0:
      /* CARA X0-UX y X1-UX e11 */
      memset(value_x0_ux, 0.0, nnods_x0 );
      ierr = VecSetValues( *x, nnods_x0, index_x0_ux, value_x0_ux, INSERT_VALUES); CHKERRQ(ierr);
      memset(value_x1_ux, strain_dir*LX, nnods_x1 );
      ierr = VecSetValues( *x, nnods_x1, index_x1_ux, value_x1_ux, INSERT_VALUES); CHKERRQ(ierr);
      break;
    case 1:
      /* CARA Y0-UY y Y1-UY e22 */
      memset(value_y0_uy, 0.0, nnods_y0 );
      ierr = VecSetValues( *x, nnods_y0, index_y0_uy, value_y0_uy, INSERT_VALUES); CHKERRQ(ierr);
      memset(value_y1_uy, strain_dir*LY, nnods_y1 );
      ierr = VecSetValues( *x, nnods_y1, index_y1_uy, value_y1_uy, INSERT_VALUES); CHKERRQ(ierr);
      break;
    case 2:
      /* CARA Z0-UZ y Z1-UZ e33 */
      memset(value_z0_uy, 0.0, nnods_z0 );
      ierr = VecSetValues( *x, nnods_z0, index_z0_uz, value_z0_uz, INSERT_VALUES); CHKERRQ(ierr);
      memset(value_z1_uz, strain_dir*LZ, nnods_z1 );
      ierr = VecSetValues( *x, nnods_z1, index_z1_uz, value_z1_uz, INSERT_VALUES); CHKERRQ(ierr);
      break;
    case 3:
      /* CARA X0-UX y Z1-UX e12 */
      memset(value_z0_uy, 0.0, nnods_z0 );
      ierr = VecSetValues( *x, nnods_z0, index_z0_uz, value_z0_uz, INSERT_VALUES); CHKERRQ(ierr);
      memset(value_z1_uz, strain_dir*LZ, nnods_z1 );
      ierr = VecSetValues( *x, nnods_z1, index_z1_uz, value_z1_uz, INSERT_VALUES); CHKERRQ(ierr);
      break;
    case 4:
      /* CARA X0-UX y Z1-UX e12 */
      memset(value_z0_uy, 0.0, nnods_z0 );
      ierr = VecSetValues( *x, nnods_z0, index_z0_uz, value_z0_uz, INSERT_VALUES); CHKERRQ(ierr);
      memset(value_z1_uz, strain_dir*LZ, nnods_z1 );
      ierr = VecSetValues( *x, nnods_z1, index_z1_uz, value_z1_uz, INSERT_VALUES); CHKERRQ(ierr);
      break;
    case 5:
      /* CARA X0-UX y Z1-UX e12 */
      memset(value_z0_uy, 0.0, nnods_z0 );
      ierr = VecSetValues( *x, nnods_z0, index_z0_uz, value_z0_uz, INSERT_VALUES); CHKERRQ(ierr);
      memset(value_z1_uz, strain_dir*LZ, nnods_z1 );
      ierr = VecSetValues( *x, nnods_z1, index_z1_uz, value_z1_uz, INSERT_VALUES); CHKERRQ(ierr);
      break;
    default:
      break;
  }
  /* communication between processes */
  ierr = VecAssemblyBegin(*x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*x);CHKERRQ(ierr);
  return 0;
}
