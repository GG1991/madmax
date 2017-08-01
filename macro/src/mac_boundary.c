/*
   MACRO boundary management routines

   Author> Guido Giuntoli
   Date> 01-08-2017
 */

#include "macro.h"

int MacroFillBoundary(MPI_Comm PROBLEM_COMM, list_t *boundary_list)
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
  int NodeOrig, NodeLocal, NodeGlobal, numnodes, kind, NDirPerNode= -1, NNeuPerNode = -1;
  int DirCount, NeuCount, *p, d, n;

  mac_boundary_t *mac_boundary;

  node_list_t *pBound = boundary_list->head;
  while(pBound)
  {
    /* 
       asignamos el NNods, allocamos memory y guardamos los nods 
       (ya van a estar ordenados y sin repetir)
     */
    numnodes = ((boundary_t *)pBound->data)->Nods.sizelist;

    mac_boundary = ((boundary_t *)pBound->data)->bvoid;
    kind = mac_boundary->kind;

    if(kind==0)                  { NDirPerNode = 0; NNeuPerNode = 3;}
    if(kind==1||kind==2||kind==4){ NDirPerNode = 1; NNeuPerNode = 2;}
    if(kind==3||kind==5||kind==6){ NDirPerNode = 2; NNeuPerNode = 1;}
    if(kind==7)                  { NDirPerNode = 3; NNeuPerNode = 0;}

    mac_boundary->NNods             = numnodes;
    mac_boundary->Nods              = malloc( numnodes * sizeof(int) );
    mac_boundary->indeces           = malloc( numnodes * 3 * sizeof(int) );
    mac_boundary->values            = malloc( numnodes * 3 * sizeof(double) );
    mac_boundary->NDirPerNode       = NDirPerNode;
    mac_boundary->NNeuPerNode       = NNeuPerNode;
    mac_boundary->NDirIndeces       = numnodes * NDirPerNode;
    mac_boundary->DirichletIndeces  = malloc( numnodes * NDirPerNode * sizeof(int) );
    mac_boundary->DirichletValues   = malloc( numnodes * NDirPerNode * sizeof(double) );
    mac_boundary->NNeuIndeces       = numnodes * NNeuPerNode;
    mac_boundary->NeumannIndeces    = malloc( numnodes * NNeuPerNode * sizeof(int) );
    mac_boundary->NeumannValues     = malloc( numnodes * NNeuPerNode * sizeof(double) );

    n=0; DirCount = 0; NeuCount = 0;
    while(n<numnodes)
    {
      /* we set the Original node numeration first */
      NodeOrig = *(int*)(((boundary_t*)pBound->data)->Nods.head->data); 
      mac_boundary->Nods[n] = NodeOrig;

      p = bsearch(&NodeOrig, MyNodOrig, NMyNod, sizeof(int), cmpfunc); 
      if(!p){SETERRQ2(PROBLEM_COMM,1,
	  "A boundary node (%d) seems now to not belong to this process (rank:%d)",NodeOrig,rank);}

      NodeLocal  = p - MyNodOrig;        // Local numeration
      NodeGlobal = loc2petsc[NodeLocal]; // PETSc numeration

      for(d=0;d<3;d++){
	if( (kind & (1<<d)) == (1<<d) ){ /* Dirichlet */
	  mac_boundary->DirichletIndeces[DirCount] = NodeGlobal*3 + d; DirCount++;
	}
	else{ /* Neumann */
	  mac_boundary->NeumannIndeces[NeuCount]   = NodeGlobal*3 + d; NeuCount++;
	}
      }
      n++;
      list_delfirst( &(((boundary_t *)pBound->data)->Nods) ) ; 
    }

    if(((boundary_t*)pBound->data)->Nods.head) // simple check 
      SETERRQ(PROBLEM_COMM,1,"It's seems that there some more nodes in the list.");

    list_clear(&(((boundary_t *)pBound->data)->Nods)); 

    /* completamos el <bvoid> */
    //    ((boundary_t *)pBound->data)->bvoid = malloc(sizeof(mac_boundary_t));
    //    memcpy(((boundary_t *)pBound->data)->bvoid, &mac_boundary, sizeof(mac_boundary_t)); 

    pBound = pBound->next;
  }
  return 0;
}
/****************************************************************************************************/
int MacroParseBoundary(MPI_Comm *PROBLEM_COMM, char *input )
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
  char   buf[NBUF];
  char   *data;
  int    ln = 0;
  int    flag_start_boundary = 0;

  mac_boundary_t mac_boundary;
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
	  if(mac_boundary.kind<0 || mac_boundary.kind>7)SETERRQ(PETSC_COMM_SELF,1,"Bad <kind> code on boundary element.");

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
	//      PetscPrintf(*PROBLEM_COMM, "# of boundaries found in %s : %d\n", input, boundary_list.sizelist);
	return 0;
      }
    }
  }
  // any boundary condition found
  // printf("SpuParseBoundary: Any boundary found on input file\n");
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

  int    i, d, numnodes, kind;
  int    *pToDirIndeces;// *pToNeuIndeces;
  int    NDirIndeces;// NNeuIndeces;
//  int    *pToIndeces;
  double *pToDirValues, *pToNeuValues;
  double ValueToSet;
  int    ierr;
  int    ofs_neu, ofs_dir;
  int    NDirPerNode;
  int    NNeuPerNode;

  node_list_t *pBound;
  f1d_t  *f1d_aux;
  mac_boundary_t *mac_boundary;

  pBound = boundary_list.head;
  while(pBound)
  {
    mac_boundary  = (mac_boundary_t*)((boundary_t*)pBound->data)->bvoid;
    numnodes      = mac_boundary->NNods;
    NDirPerNode   = mac_boundary->NDirPerNode;
    NNeuPerNode   = mac_boundary->NNeuPerNode;
    pToDirIndeces = mac_boundary->DirichletIndeces;
    NDirIndeces   = mac_boundary->NDirIndeces;
//    pToNeuIndeces = mac_boundary->NeumannIndeces;
//    NNeuIndeces   = mac_boundary->NNeuIndeces;
    kind          = mac_boundary->kind;
    pToDirValues  = mac_boundary->DirichletValues;
    pToNeuValues  = mac_boundary->NeumannValues;
    ofs_neu = ofs_dir = 0;

    for(d=0;d<3;d++)
    {
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
int MacroSetBoundaryOnJacobian( Mat *J )
{
  /* 
     Sets 1's on the diagonal corresponding to Dirichlet indeces and 0's
     on the rest of the row and column 
  */

  int    *pToDirIndeces;
  int    NDirIndeces;
  int    ierr;

  node_list_t *pBound;
  mac_boundary_t *mac_boundary;

  pBound = boundary_list.head;
  while(pBound){
    mac_boundary  = (mac_boundary_t*)((boundary_t*)pBound->data)->bvoid;
    pToDirIndeces = mac_boundary->DirichletIndeces;
    NDirIndeces   = mac_boundary->NDirIndeces;
    ierr = MatZeroRowsColumns(*J, NDirIndeces, pToDirIndeces, 1.0, NULL, NULL); CHKERRQ(ierr);
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

  int    *pToDirIndeces, NDirIndeces;
  double *pToDirValues;
  int    ierr;

  node_list_t *pBound;

  pBound = boundary_list.head;
  while(pBound)
  {
    pToDirIndeces = ((mac_boundary_t*)pBound->data)->DirichletIndeces;
    pToDirValues  = ((mac_boundary_t*)pBound->data)->DirichletValues;
    NDirIndeces   = ((mac_boundary_t*)pBound->data)->NDirIndeces;

    memset(pToDirValues, 0.0, NDirIndeces*sizeof(double));
    ierr = VecSetValues( *b, NDirIndeces, pToDirIndeces, pToDirValues, INSERT_VALUES); CHKERRQ(ierr);
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
