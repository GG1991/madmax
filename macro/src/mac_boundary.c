#include "sputnik.h"

int MacroSetBoundaryDisplacement( double time, Vec *x )
{

  /* 
     Go over all <boundary_list> elements and set on 
     all nodes of that that <boundary_t> elements the
     value of the boundary condition at that time 
     (<time>) evaluating the <f1d_t> functions.

     1) primero rellenamos el array boundary.values con el valor
     de la funcion evaluada para cada elemento de la lista

     2) luego procedemos con VecSetValues ya tenemos los indices guardados

     Nota: No deberiamos hacer mallocs aqui dentro

   */

  int    i, d, numnodes, kind;
  int    *pToIndeces;
  double *pToValues;
  double ValueToSet;
  int    ierr;

  node_list_t *pBound;
  f1d_t  *f1d_aux;

  pBound = boundary_list.head;
  while(pBound)
  {
    numnodes   = ((boundary_t*)pBound->data)->NNods;
    kind       = ((boundary_t*)pBound->data)->kind;
    pToValues  = ((boundary_t*)pBound->data)->values;
    pToIndeces = ((boundary_t*)pBound->data)->indeces;

    for(d=0;d<3;d++){
      /* Barremos primero las dirección x -> y -> z */
      if( (kind & (1<<d)) == (1<<d) ){
	/* es Dirichlet */
	switch(d){
	  case 0:
	    f1d_aux = ((boundary_t*)pBound->data)->fx;
	    break;
	  case 1:
	    f1d_aux = ((boundary_t*)pBound->data)->fy;
	    break;
	  case 2:
	    f1d_aux = ((boundary_t*)pBound->data)->fz;
	    break;
	  default:
	    return 1;
	}
	f1d_eval( time, f1d_aux, &ValueToSet );
	/* pToValues = [ valx valy valz valx valy valz ... valx valy valz ] */
	for(i=0;i<numnodes;i++) pToValues[i*3+d] = ValueToSet;
      }
    }
    //  usamos VecSetValuesLocal aqui ya que vamos a modificar valores locales unicamente
    ierr = VecSetValuesLocal( *x, numnodes, pToIndeces, pToValues, INSERT_VALUES); CHKERRQ(ierr);
    pBound = pBound->next;
  }
  /* communication between processes */
  ierr = VecAssemblyBegin(*x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*x);CHKERRQ(ierr);
  return 0;
}
