#include "sputnik.h"

int SputnikSetDisplacementOnBoundary( double time, Vec *x )
{

  /* 
     Go over all <boundary_list> elements and set on 
     all nodes of that that <boundary_t> elements the
     value of the boundary condition at that time 
     (<time>) evaluating the <f1d_t> functions.

     1) primero rellenamos el array boundary.values con el valor
     de la funcion evaluada para cada elemento de la lista

     2) luego procedemos con VecSetValues ya tenemos los indices guardados

     Nota > No deberiamos hacer mallocs aqui dentro

   */

  int    i, d, numnodes, kind;
  int    *pToDirIndeces, *pToNeuIndeces;
  int    NDirIndeces, NNeuIndeces;
  int    *pToIndeces;
  double *pToDirValues, *pToNeuValues;
  double ValueToSet;
  int    ierr;
  int    ofs_neu, ofs_dir;
  int    NDirPerNode;
  int    NNeuPerNode;

  node_list_t *pBound;
  f1d_t  *f1d_aux;

  pBound = boundary_list.head;
  while(pBound)
  {
    numnodes      = ((boundary_t*)pBound->data)->NNods;
    NDirPerNode   = ((boundary_t*)pBound->data)->NDirPerNode;
    NNeuPerNode   = ((boundary_t*)pBound->data)->NNeuPerNode;
    pToDirIndeces = ((boundary_t*)pBound->data)->DirichletIndeces;
    NDirIndeces   = ((boundary_t*)pBound->data)->NDirIndeces;
    pToNeuIndeces = ((boundary_t*)pBound->data)->NeumannIndeces;
    NNeuIndeces   = ((boundary_t*)pBound->data)->NNeuIndeces;
    numnodes      = ((boundary_t*)pBound->data)->NNods;
    kind          = ((boundary_t*)pBound->data)->kind;
    pToDirValues  = ((boundary_t*)pBound->data)->DirichletValues;
    pToNeuValues  = ((boundary_t*)pBound->data)->NeumannValues;
    ofs_neu = ofs_dir = 0;

    for(d=0;d<3;d++)
    {
      /* Barremos primero las dirección x -> y -> z */
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

      if( (kind & (1<<d)) == (1<<d) ){
	/* es Dirichlet */
	for(i=0;i<numnodes;i++){ pToDirValues[i*NDirPerNode + ofs_dir] = ValueToSet;}
	ofs_dir++;
      }
      else{
	/* es Neumann */
	for(i=0;i<numnodes;i++){ pToNeuValues[i*NNeuPerNode + ofs_neu] = ValueToSet;}
	ofs_neu++;
      }
      /* pToValues = [ valx valy valz valx valy valz ... valx valy valz ] */
    }

    /*  
	Dirichlet Boundary condition set is set on <x> 
	usamos VecSetValuesLocal aqui ya que vamos a modificar valores locales unicamente
     */
    ierr = VecSetValues( *x, NDirIndeces, pToDirIndeces, pToDirValues, INSERT_VALUES); CHKERRQ(ierr);
    pBound = pBound->next;
  }
  /* communication between processes */
  ierr = VecAssemblyBegin(*x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*x);CHKERRQ(ierr);
  return 0;
}

/****************************************************************************************************/

int SputnikSetBoundaryOnJacobian( Mat *J )
{
  /* 
     Sets 1's on the diagonal corresponding to Dirichlet indeces and 0's
     on the rest of the row and column 
  */

  int    i;
  int    *pToDirIndeces;
  int    NDirIndeces;
  int    ierr;

  node_list_t *pBound;

  pBound = boundary_list.head;
  while(pBound)
  {
    pToDirIndeces = ((boundary_t*)pBound->data)->DirichletIndeces;
    NDirIndeces   = ((boundary_t*)pBound->data)->NDirIndeces;
    ierr = MatZeroRowsColumns(*J, NDirIndeces, pToDirIndeces, 1.0, NULL, NULL); CHKERRQ(ierr);
    pBound = pBound->next;
  }

  /* communication between processes */
  ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  return 0;
}

/****************************************************************************************************/

int SputnikSetBoundaryOnResidual( Vec *b )
{
  /* 
     Sets 0's on the Dirichlet indeces over the Residual <b>
  */

  int    *pToDirIndeces, NDirIndeces;
  double *pToDirValues;
  int    ierr;

  node_list_t *pBound;

  pBound = boundary_list.head;
  while(pBound)
  {
    pToDirIndeces = ((boundary_t*)pBound->data)->DirichletIndeces;
    pToDirValues  = ((boundary_t*)pBound->data)->DirichletValues;
    NDirIndeces   = ((boundary_t*)pBound->data)->NDirIndeces;

    memset(pToDirValues, 0.0, NDirIndeces*sizeof(double));
    ierr = VecSetValues( *b, NDirIndeces, pToDirIndeces, pToDirValues, INSERT_VALUES); CHKERRQ(ierr);
    pBound = pBound->next;
  }
  /* communication between processes */
  ierr = VecAssemblyBegin(*b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*b);CHKERRQ(ierr);
  return 0;
}
