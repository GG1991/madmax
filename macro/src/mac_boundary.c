/*
   MACRO boundary management routines

   Author> Guido Giuntoli
   Date> 01-08-2017
 */

#include "macro.h"

int MacroFillBoundary(MPI_Comm PROBLEM_COMM, list_t *boundary_list_aux, list_t *boundary_list)
{
  /* 
     Here we fill the boundary_list structure acording to our problem 
     macro where it contains boundary_t elements
   */
  // ahora metemos los nodos de las listas en la estructura posta <boundary_list>
  int NodeOrig, NodeLocal, NodeGlobal, numnodes, kind, NDirPerNode= -1, NNeuPerNode = -1;
  int DirCount, NeuCount, *p, d, n;
  int nproc, rank;

  MPI_Comm_size(PROBLEM_COMM, &nproc);
  MPI_Comm_rank(PROBLEM_COMM, &rank);

  node_list_t *pBound = boundary_list->head;
  node_list_t *pAuxBound = boundary_list_aux->head;
  while(pBound)
  {
    /* 
       asignamos el NNods, allocamos memory y guardamos los nods 
       (ya van a estar ordenados y sin repetir)
     */
    numnodes = ((AuxBoundary_t *)pAuxBound->data)->Nods.sizelist;
    kind     = ((boundary_t *)pBound->data)->kind;

    if(kind==0)                  { NDirPerNode = 0; NNeuPerNode = 3;}
    if(kind==1||kind==2||kind==4){ NDirPerNode = 1; NNeuPerNode = 2;}
    if(kind==3||kind==5||kind==6){ NDirPerNode = 2; NNeuPerNode = 1;}
    if(kind==7)                  { NDirPerNode = 3; NNeuPerNode = 0;}

    ((boundary_t*)pBound->data)->NNods             = numnodes;
    ((boundary_t*)pBound->data)->Nods              = malloc( numnodes * sizeof(int) );
    ((boundary_t*)pBound->data)->indeces           = malloc( numnodes * 3 * sizeof(int) );
    ((boundary_t*)pBound->data)->values            = malloc( numnodes * 3 * sizeof(double) );
    ((boundary_t*)pBound->data)->NDirPerNode       = NDirPerNode;
    ((boundary_t*)pBound->data)->NNeuPerNode       = NNeuPerNode;
    ((boundary_t*)pBound->data)->NDirIndeces       = numnodes * NDirPerNode;
    ((boundary_t*)pBound->data)->DirichletIndeces  = malloc( numnodes * NDirPerNode * sizeof(int) );
    ((boundary_t*)pBound->data)->DirichletValues   = malloc( numnodes * NDirPerNode * sizeof(double) );
    ((boundary_t*)pBound->data)->NNeuIndeces       = numnodes * NNeuPerNode;
    ((boundary_t*)pBound->data)->NeumannIndeces    = malloc( numnodes * NNeuPerNode * sizeof(int) );
    ((boundary_t*)pBound->data)->NeumannValues     = malloc( numnodes * NNeuPerNode * sizeof(double) );

    n=0; DirCount = 0; NeuCount = 0;
    while(n<numnodes)
    {
      /* we set the Original node numeration first */
      NodeOrig = *(int*)(((AuxBoundary_t *)pAuxBound->data)->Nods.head->data); 
      ((boundary_t*)pBound->data)->Nods[n] = NodeOrig;

      p = bsearch(&NodeOrig, MyNodOrig, NMyNod, sizeof(int), cmpfunc); 
      if(!p){SETERRQ2(PROBLEM_COMM,1,
	  "A boundary node (%d) seems now to not belong to this process (rank:%d)",NodeOrig,rank);}

      NodeLocal  = p - MyNodOrig;        // Local numeration
      NodeGlobal = loc2petsc[NodeLocal]; // PETSc numeration

      for(d=0;d<3;d++){
	/* Dirichlet */
	if( (kind & (1<<d)) == (1<<d) ){
	  ((boundary_t*)pBound->data)->DirichletIndeces[DirCount] = NodeGlobal*3 + d; DirCount++;
	}
	/* Neumann */
	else{
	  ((boundary_t*)pBound->data)->NeumannIndeces[NeuCount]   = NodeGlobal*3 + d; NeuCount++;
	}
      }
      n++;
      list_delfirst( &(((AuxBoundary_t *)pAuxBound->data)->Nods) ) ;
    }

    if(((AuxBoundary_t *)pAuxBound->data)->Nods.head)
      SETERRQ(PROBLEM_COMM,1,"It's seems that there some more nodes in the list.");

    list_clear(&(((AuxBoundary_t *)pAuxBound->data)->Nods));
    pBound = pBound->next;	pAuxBound = pAuxBound->next;
  }
  list_clear(boundary_list_aux);
  return 0;
}

