#include "util.h" 

//int get_nods_bc(int **nods, int *nnods)
//{
//
//  int          nnods_bc=-1, nnods_bc_a=-1, c, n_orig; 
//  int          i, *aux_nod;
//  node_list_t  *pb,*pn;
//
//  /* count the number of total node with repetitions on boundary list */
//  nnods_bc=0;
//  pb = boundary_list.head;
//  while(pb){
//    nnods_bc += ((boundary_t*)pb->data)->Nods.sizelist;
//    pb=pb->next;
//  }
//  aux_nod = malloc(nnods_bc*sizeof(int));
//
//  i=0;
//  pb = boundary_list.head;
//  while(pb){
//    pn = ((boundary_t*)pb->data)->Nods.head;
//    while(pn){
//      n_orig     = *(int*)pn->data;
//      aux_nod[i] = n_orig;
//      i++;
//      pn=pn->next;
//    }
//    pb=pb->next;
//  }
//
//  /*
//     we order the vector of boundary nodes in order to delete those 
//     nodes which are repeated
//   */
//  qsort(aux_nod, nnods_bc, sizeof(int), cmpfunc);
//  for(i=0;i<nnods_bc;i++){
//    if(i!=0){
//      if(aux_nod[i]!=aux_nod[i-1]){
//	nnods_bc_a++;
//      }
//    }
//    else if(i==0){
//      nnods_bc_a=1;
//    }
//  }
//
//  *nnods = nnods_bc_a;
//  *nods  = malloc(nnods_bc_a*sizeof(int));
//
//  c = 0;
//  for(i=0;i<nnods_bc;i++){
//    if(i!=0){
//      if(aux_nod[i]!=aux_nod[i-1]){
//	(*nods)[c] = aux_nod[i];
//	c++;
//      }
//    }
//    else if(i==0){
//      (*nods)[c] = aux_nod[i];
//      c++;
//    }
//  }
//
//  return 0;
//}
///****************************************************************************************************/
//int get_nods_index(int *nods_bc, int nnods_bc, int *ix_loc, int *ix_glo)
//{
//  /* 
//     ix_loc -> size =  nnods_bc x dim
//     ix_glo -> size =  nnods_bc x dim
//   */
//  int i, d, *p;
//
//  if(!ix_loc || !ix_glo) return 1;
//
//  for(i=0;i<nnods_bc;i++){
//    p = bsearch(&nods_bc[i], mynods, nmynods, sizeof(int), cmpfunc); 
//    if(!p){
//      PetscPrintf(PETSC_COMM_SELF,"bc node %d not belong <mynods>",nods_bc[i]);
//      return 1;
//    }
//    for(d=0;d<dim;d++){
//      ix_loc[i*dim+d] = (p - mynods)*dim + d;
//      ix_glo[i*dim+d] = loc2petsc[p - mynods]*dim + d;
//    }
//  }
//  return 0;
//}

/****************************************************************************************************/

#ifdef PETSC
int print_ksp_info(MPI_Comm COMM, KSP ksp)
{
 
  int    kspits, reason;
  double kspnorm;
  char   *reason_s;
  
  KSPGetIterationNumber( ksp, &kspits );
  KSPGetConvergedReason( ksp, &reason );
  KSPGetResidualNorm   ( ksp, &kspnorm);
  switch(reason)
  {
    case KSP_CONVERGED_RTOL:
      reason_s = strdup( "RTOL" );
      break;
    case KSP_CONVERGED_ATOL:
      reason_s = strdup( "ATOL" );
      break;
    default :
      reason_s = strdup( "UNKNOW" );
      break;
  }
  PetscPrintf( COMM,"kspits %D kspnorm %e kspreason %s\n", kspits, kspnorm, reason_s );
  return 0;
}
#endif

/****************************************************************************************************/
