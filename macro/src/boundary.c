/*
   MACRO boundary management routines

   Author> Guido Giuntoli
   Date> 01-08-2017
 */

#include "macro.h"

int mac_init_boundary(MPI_Comm PROBLEM_COMM, list_t *boundary_list)
{
  /* 
     Here we fill the <boundary_list> structure acording to our problem 
     macro where it contains <mac_boundary_t> elements
   */

  int nproc, rank;

  MPI_Comm_size(PROBLEM_COMM, &nproc);
  MPI_Comm_rank(PROBLEM_COMM, &rank);

  /* 
     Terminamos de leer el archivo ahora quitamos 
     los repetidos en order decendente 
   */
  int i, j, h, k;
  node_list_t  *nbb, *nba, *nnb, *nna;

  for(i=boundary_list->sizelist;i>0;i--){
    nbb=boundary_list->head;
    for(j=0;j<i-1;j++){
      nbb=nbb->next;
    }
    nba=boundary_list->head;
    for(j=0;j<i-1;j++){
      nna=((boundary_t*)nba->data)->Nods.head;
      for(h=0;h<((boundary_t*)nba->data)->Nods.sizelist;h++){
	nnb=((boundary_t*)nbb->data)->Nods.head;
	for(k=0;k<((boundary_t*)nbb->data)->Nods.sizelist;k++){
	  if(*(int*)nnb->data==*(int*)nna->data){
	    if(list_del(&((boundary_t*)nbb->data)->Nods,nnb)){
	      return 1;
	    }
	    break;
	  }
	  nnb=nnb->next;
	}
	nna=nna->next;
      }
      nba=nba->next;
    }
  }

  // ahora metemos los nodos de las listas en la estructura posta <boundary_list>
  int NodeOrig, NodeLocal, NodeGlobal, numnodes, kind, ndir_pn=-1, nneu_pn=-1;
  int DirCount, NeuCount, *p, d, n;

  mac_boundary_t *mac_boundary;

  bc_kinds=malloc(nmynods*sizeof(int));
  for(i=0;i<nmynods;i++)bc_kinds[i]=-1;

  node_list_t *pBound = boundary_list->head;
  while(pBound)
  {
    /* 
       asignamos el nnod, allocamos memory y guardamos los nods 
       (ya van a estar ordenados y sin repetir)
     */
    numnodes = ((boundary_t *)pBound->data)->Nods.sizelist;

    mac_boundary = ((boundary_t *)pBound->data)->bvoid;
    kind = mac_boundary->kind;
    if(dim==2){
      if(kind==0)                  { ndir_pn = 0; nneu_pn = 2;}
      if(kind==1||kind==2)         { ndir_pn = 1; nneu_pn = 1;}
      if(kind==3)                  { ndir_pn = 2; nneu_pn = 0;}
    }
    else if(dim==3){
      if(kind==0)                  { ndir_pn = 0; nneu_pn = 3;}
      if(kind==1||kind==2||kind==4){ ndir_pn = 1; nneu_pn = 2;}
      if(kind==3||kind==5||kind==6){ ndir_pn = 2; nneu_pn = 1;}
      if(kind==7)                  { ndir_pn = 3; nneu_pn = 0;}
    }

    mac_boundary->nnod     = numnodes;
    mac_boundary->nods     = malloc(numnodes * sizeof(int));
    mac_boundary->ndir_pn  = ndir_pn;
    mac_boundary->nneu_pn  = nneu_pn;
    mac_boundary->ndir     = numnodes*ndir_pn;
    mac_boundary->nneu     = numnodes*nneu_pn;
    mac_boundary->dir_idx  = malloc(numnodes * ndir_pn * sizeof(int));
    mac_boundary->dir_val  = malloc(numnodes * ndir_pn * sizeof(double));
    mac_boundary->neu_idx  = malloc(numnodes * nneu_pn * sizeof(int));
    mac_boundary->neu_val  = malloc(numnodes * nneu_pn * sizeof(double));

    node_list_t  *pn = ((boundary_t*)pBound->data)->Nods.head;
    n=0; DirCount = 0; NeuCount = 0;
    while(n<numnodes)
    {
      /* we set the Original node numeration first */
      NodeOrig = *(int*)(pn->data); 
      mac_boundary->nods[n] = NodeOrig;

      p = bsearch(&NodeOrig, mynods, nmynods, sizeof(int), cmpfunc); 
      if(!p){
	return 1;
      }

      NodeLocal  = p - mynods;           // Local numeration
      NodeGlobal = loc2petsc[NodeLocal]; // PETSc numeration
      bc_kinds[NodeLocal] = kind;

      for(d=0;d<dim;d++){
	if( (kind & (1<<d)) == (1<<d) ){ /* Dirichlet */
	  mac_boundary->dir_idx[DirCount] = NodeGlobal*dim + d; DirCount++;
	}
	else{ /* Neumann */
	  mac_boundary->neu_idx[NeuCount] = NodeGlobal*dim + d; NeuCount++;
	}
      }
      n++;
      pn = pn->next;
    }

    pBound = pBound->next;
  }
  return 0;
}
/****************************************************************************************************/
int macro_parse_boundary(MPI_Comm PROBLEM_COMM, char *input)
{
  /*
     Parse the boundary of the problem

     Searchs for keywords>

     $Boundary
     <name1> <order> <kind> <fnumx> <fnumy> <fnumz>
     <name2> <order> <kind> <fnumx> <fnumy> <fnumz>
     ...
     $EndBoundary
   */

  FILE   *file = fopen(input,"r"); if(!file) SETERRQ(PETSC_COMM_SELF,1,"File not found");
  char   buf[NBUF], *data;
  int    ln=0, flag_start_boundary=0;

  mac_boundary_t mac_boundary;

  mac_boundary.kind    = -1;
  mac_boundary.order   = -1;
  mac_boundary.nfx     = mac_boundary.nfy = mac_boundary.nfz = -1;
  mac_boundary.fx      = mac_boundary.fy  = mac_boundary.fz  = NULL;
  mac_boundary.nnod    = -1;
  mac_boundary.nods    = NULL;
  mac_boundary.ndir    = mac_boundary.ndir_pn = -1;
  mac_boundary.dir_idx = NULL;
  mac_boundary.dir_val = NULL;
  mac_boundary.nneu    = mac_boundary.nneu_pn = -1;
  mac_boundary.neu_idx = NULL;
  mac_boundary.neu_val = NULL;

  boundary_t boundary;
  list_init(&boundary_list, sizeof(boundary_t), cmpfunc_mac_bou);
  list_init(&boundary.Nods, sizeof(int), cmpfunc_for_list);

  while(fgets(buf,NBUF,file) != NULL)
  {

    ln ++;
    data = strtok(buf," \n");
    if(data){ 
      if(!strcmp(data,"$Boundary")){

	flag_start_boundary=1;
	while(fgets(buf,NBUF,file) != NULL)
	{
	  ln ++;

	  // <name>
	  data = strtok(buf," \n"); 
	  if(!strcmp(data,"$EndBoundary")) break;
	  boundary.name = strdup(data);

	  // <order> 
	  data = strtok(NULL," \n"); 
	  mac_boundary.order = atoi(data);

	  // <kind> 
	  data = strtok(NULL," \n");
	  mac_boundary.kind = StrBin2Dec(data);
	  if(mac_boundary.kind<0 || mac_boundary.kind>7){
	    SETERRQ(PETSC_COMM_SELF,1,"Bad <kind> code on boundary element.");
	  }

	  // <nfz> 
	  data = strtok(NULL," \n");
	  mac_boundary.nfz = atoi(data);

	  // <nfy> 
	  data = strtok(NULL," \n");
	  mac_boundary.nfy = atoi(data);

	  // <nfx> 
	  data = strtok(NULL," \n");
	  mac_boundary.nfx = atoi(data);

	  mac_boundary.fx = GetFunctionPointer(&function_list,mac_boundary.nfx);
	  mac_boundary.fy = GetFunctionPointer(&function_list,mac_boundary.nfy);
	  mac_boundary.fz = GetFunctionPointer(&function_list,mac_boundary.nfz);

	  // si llegamos hasta acá esta todo 0K lo insertamos en la lista 
	  boundary.bvoid = malloc(sizeof(mac_boundary_t));
	  memcpy(boundary.bvoid,&mac_boundary,sizeof(mac_boundary_t));
	  list_insert_se(&boundary_list, &boundary);
	}
      } // inside $Boundary

      if(!strcmp(data,"$EndBoundary")){
	CHKERRQ(!flag_start_boundary);
	return 0;
      }
    }
  }
  SETERRQ(PETSC_COMM_SELF,1,"any boundary condition found");
  return 1;
}
/****************************************************************************************************/
int MacroSetDisplacementOnBoundary( double time, Vec *x )
{

  /* 
     Go over all <boundary_list> elements and set on 
     all nodes of that that <mac_boundary_t> elements the
     value of the boundary condition at that time 
     (<time>) evaluating the <f1d_t> functions.

     1) primero rellenamos el array boundary.values con el valor
     de la funcion evaluada para cada elemento de la lista

     2) luego procedemos con VecSetValues ya tenemos los indices guardados

     Nota > No deberiamos hacer mallocs aqui dentro

   */

  int    i, d, nnod, kind;
  int    *dir_idx;
  int    ndir;
  double *dir_val, *neu_val;
  double ValueToSet;
  int    ierr;
  int    ofs_neu, ofs_dir;
  int    ndir_pn;
  int    nneu_pn;

  node_list_t *pBound;
  f1d_t  *f1d_aux;
  mac_boundary_t *mac_boundary;

  pBound = boundary_list.head;
  while(pBound)
  {
    mac_boundary  = ((boundary_t*)pBound->data)->bvoid;
    nnod          = mac_boundary->nnod;
    ndir_pn       = mac_boundary->ndir_pn;
    nneu_pn       = mac_boundary->nneu_pn;
    dir_idx       = mac_boundary->dir_idx;
    ndir          = mac_boundary->ndir;
    kind          = mac_boundary->kind;
    dir_val       = mac_boundary->dir_val;
    neu_val       = mac_boundary->neu_val;
    ofs_neu       = ofs_dir = 0;

    for(d=0;d<dim;d++){
      /* Barremos primero las dirección x -> y -> z */
      switch(d){
	case 0:
	  f1d_aux = mac_boundary->fx;
	  break;
	case 1:
	  f1d_aux = mac_boundary->fy;
	  break;
	case 2:
	  f1d_aux = mac_boundary->fz;
	  break;
	default:
	  return 1;
      }
      f1d_eval( time, f1d_aux, &ValueToSet );

      if( (kind & (1<<d)) == (1<<d) ){
	/* es Dirichlet */
	for(i=0;i<nnod;i++){ dir_val[i*ndir_pn + ofs_dir] = ValueToSet;}
	ofs_dir++;
      }
      else{
	/* es Neumann */
	for(i=0;i<nnod;i++){ neu_val[i*nneu_pn + ofs_neu] = ValueToSet;}
	ofs_neu++;
      }
      /* pToValues = [ valx valy valz valx valy valz ... valx valy valz ] */
    }

    /*  
	Dirichlet Boundary condition set is set on <x> 
	usamos VecSetValuesLocal aqui ya que vamos a modificar valores locales unicamente
     */
    ierr = VecSetValues( *x, ndir, dir_idx, dir_val, INSERT_VALUES); CHKERRQ(ierr);
    pBound = pBound->next;
  }
  /* communication between processes */
  ierr = VecAssemblyBegin(*x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*x);CHKERRQ(ierr);
  return 0;
}
/****************************************************************************************************/
int MacroSetBoundaryOnJacobian( Mat *J )
{
  /* 
     Sets 1's on the diagonal corresponding to Dirichlet indeces and 0's
     on the rest of the row and column 
  */

  int    *dir_idx;
  int    ndir;
  int    ierr;

  node_list_t *pBound;
  mac_boundary_t *mac_boundary;

  pBound = boundary_list.head;
  while(pBound){
    mac_boundary  = ((boundary_t*)pBound->data)->bvoid;
    dir_idx = mac_boundary->dir_idx;
    ndir = mac_boundary->ndir;
    ierr = MatZeroRowsColumns(*J, ndir, dir_idx, 1.0, NULL, NULL); CHKERRQ(ierr);
    pBound = pBound->next;
  }

  /* communication between processes */
  ierr = MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  return 0;
}
/****************************************************************************************************/
int MacroSetBoundaryOnResidual( Vec *b )
{
  /* 
     Sets 0's on the Dirichlet indeces over the Residual <b>
  */
  int    ierr;
  node_list_t *pBound;
  mac_boundary_t *mac_boundary;

  pBound = boundary_list.head;
  while(pBound){
    mac_boundary  = ((boundary_t*)pBound->data)->bvoid;
    memset(mac_boundary->dir_val, 0.0, mac_boundary->ndir*sizeof(double));
    ierr = VecSetValues(*b,mac_boundary->ndir,mac_boundary->dir_idx,mac_boundary->dir_val,INSERT_VALUES);CHKERRQ(ierr);
    pBound = pBound->next;
  }
  /* communication between processes */
  ierr = VecAssemblyBegin(*b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(*b);CHKERRQ(ierr);
  return 0;
}
/****************************************************************************************************/
int cmpfunc_mac_bou (void * a, void * b)
{
  return ( ((mac_boundary_t*)(((boundary_t*)a)->bvoid))->order - ((mac_boundary_t*)(((boundary_t*)b)->bvoid))->order );
}
/****************************************************************************************************/
