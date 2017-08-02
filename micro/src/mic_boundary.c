/*
   Routines to impose Boundary condition on RVE cells

   Author > Guido Giuntoli
   Date > 31-07-2017
*/

#include "micro.h"

/****************************************************************************************************/
int MicroCreateBoundary(list_t *boundary_list)
{
  /*
     Creates the boundary list with names
     P000 P100 P010 X0 X1 Y0 Y1 Z0 Z1 
   */
  int i = 0;
  boundary_t boundary;

  list_init(boundary_list, sizeof(boundary_t), NULL);
  list_init(&boundary.Nods,sizeof(int), cmpfunc_for_list);
  while(i<9)
  {
    switch(i){
      case 0:
	boundary.name = strdup("P000");break;
      case 1:
	boundary.name = strdup("P100");break;
      case 2:
	boundary.name = strdup("P010");break;
      case 3:
	boundary.name = strdup("X0")  ;break;
      case 4:
	boundary.name = strdup("X1")  ;break;
      case 5:
	boundary.name = strdup("Y0")  ;break;
      case 6:
	boundary.name = strdup("Y1")  ;break;
      case 7:
	boundary.name = strdup("Z0")  ;break;
      case 8:
	boundary.name = strdup("Z1")  ;break;
      default:
	break;
    }
    list_insertlast(boundary_list, &boundary);
    i++;
  }
  return 0;
}
/****************************************************************************************************/
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
/****************************************************************************************************/
int MicroCheckPhysicalEntities( list_t *physical_list )
{
  /*
     Checks if the physical entities defined on mesh file are
     P000 P100 P010 X0 X1 Y0 Y1 Z0 Z1 (at least)
   */
  int i = 0, flag = 0, flag_pn = 0;
  char *name;
  while(i<9)
  {
    node_list_t * pn = physical_list->head;
    flag_pn=0;
    while(pn && !flag_pn)
    {
      name = ((physical_t*)pn->data)->name;
      switch(i){
	case 0:
	  if(!strcmp(name,"P000")){flag=flag|(1<<0);flag_pn=1;}break;
	case 1:
	  if(!strcmp(name,"P100")){flag=flag|(1<<1);flag_pn=1;}break;
	case 2:
	  if(!strcmp(name,"P010")){flag=flag|(1<<2);flag_pn=1;}break;
	case 3:
	  if(!strcmp(name,"X0")  ){flag=flag|(1<<3);flag_pn=1;}break;
	case 4:
	  if(!strcmp(name,"X1")  ){flag=flag|(1<<4);flag_pn=1;}break;
	case 5:
	  if(!strcmp(name,"Y0")  ){flag=flag|(1<<5);flag_pn=1;}break;
	case 6:
	  if(!strcmp(name,"Y1")  ){flag=flag|(1<<6);flag_pn=1;}break;
	case 7:
	  if(!strcmp(name,"Z0")  ){flag=flag|(1<<7);flag_pn=1;}break;
	case 8:
	  if(!strcmp(name,"Z1")  ){flag=flag|(1<<8);flag_pn=1;}break;
	default:
	  break;
      }
      pn=pn->next;
    }
    i++;
  }
  if(flag != 511)SETERRQ(MICRO_COMM,1, "MICRO:physical entity not found (P000 P100 P010 X0 X1 Y0 Y1 Z0 Z1)");
  return 0;
}
/****************************************************************************************************/
int MicroSetBoundary(MPI_Comm PROBLEM_COMM, list_t *boundary_list)
{
  /*
     Checks if "P000" "P100" "P010" "X0" "X1" "Y0" "Y1" "Z0" "Z1"
     are define in the boundary_list
  */
  int i = 0, n = -1, *p, flag = 0, flag_pn = 0, nnods = 0;
  int node_petsc, node_orig, *px, *py, *pz;
  char *name;
  int nproc, rank;

  MPI_Comm_size(PROBLEM_COMM, &nproc);
  MPI_Comm_rank(PROBLEM_COMM, &rank);

  node_list_t *pb, *pn;
  while(i<9)
  {
    pb = boundary_list->head;
    flag_pn=0;
    while(pb && !flag_pn)
    {
      nnods = ((boundary_t*)pb->data)->Nods.sizelist;
      name  = ((boundary_t*)pb->data)->name;
      switch(i){
	case 0:
	  if(!strcmp(name,"P000")){flag=flag|(1<<0);flag_pn=1;}break;
	case 1:
	  if(!strcmp(name,"P100")){flag=flag|(1<<1);flag_pn=1;}break;
	case 2:
	  if(!strcmp(name,"P010")){flag=flag|(1<<2);flag_pn=1;}break;
	case 3:
	  if(!strcmp(name,"X0")){
	    flag=flag|(1<<3);
	    nnods_x0 = nnods;
	    index_x0_ux = malloc( nnods_x0 * sizeof(int)); value_x0_ux = malloc( nnods_x0 * sizeof(double));
	    index_x0_uy = malloc( nnods_x0 * sizeof(int)); value_x0_uy = malloc( nnods_x0 * sizeof(double));
	    index_x0_uz = malloc( nnods_x0 * sizeof(int)); value_x0_uz = malloc( nnods_x0 * sizeof(double));
	    px = index_x0_ux; py = index_x0_uy; pz = index_x0_uz;
	    flag_pn=1;
	  }
	  break;
	case 4:
	  if(!strcmp(name,"X1")){
	    flag=flag|(1<<4);
	    nnods_x1 = nnods;
	    index_x1_ux = malloc( nnods_x1 * sizeof(int)); value_x1_ux = malloc( nnods_x1 * sizeof(double));
	    index_x1_uy = malloc( nnods_x1 * sizeof(int)); value_x1_uy = malloc( nnods_x1 * sizeof(double));
	    index_x1_uz = malloc( nnods_x1 * sizeof(int)); value_x1_uz = malloc( nnods_x1 * sizeof(double));
	    px = index_x1_ux; py = index_x1_uy; pz = index_x1_uz;
	    flag_pn=1;
	  }
	  break;
	case 5:
	  if(!strcmp(name,"Y0")){
	    flag=flag|(1<<5);
	    nnods_y0 = nnods;
	    index_y0_ux = malloc( nnods_y0 * sizeof(int)); value_y0_ux = malloc( nnods_y0 * sizeof(double));
	    index_y0_uy = malloc( nnods_y0 * sizeof(int)); value_y0_uy = malloc( nnods_y0 * sizeof(double));
	    index_y0_uz = malloc( nnods_y0 * sizeof(int)); value_y0_uz = malloc( nnods_y0 * sizeof(double));
	    px = index_y0_ux; py = index_y0_uy; pz = index_y0_uz;
	    flag_pn=1;
	  }
	  break;
	case 6:
	  if(!strcmp(name,"Y1")){
	    flag=flag|(1<<6);
	    nnods_y1 = nnods;
	    index_y1_ux = malloc( nnods_y1 * sizeof(int)); value_y1_ux = malloc( nnods_y1 * sizeof(double));
	    index_y1_uy = malloc( nnods_y1 * sizeof(int)); value_y1_uy = malloc( nnods_y1 * sizeof(double));
	    index_y1_uz = malloc( nnods_y1 * sizeof(int)); value_y1_uz = malloc( nnods_y1 * sizeof(double));
	    px = index_y1_ux; py = index_y1_uy; pz = index_y1_uz;
	    flag_pn=1;
	  }
	  break;
	case 7:
	  if(!strcmp(name,"Z0")){
	    flag=flag|(1<<7);
	    nnods_z0 = nnods;
	    index_z0_ux = malloc( nnods_z0 * sizeof(int)); value_z0_ux = malloc( nnods_z0 * sizeof(double));
	    index_z0_uy = malloc( nnods_z0 * sizeof(int)); value_z0_uy = malloc( nnods_z0 * sizeof(double));
	    index_z0_uz = malloc( nnods_z0 * sizeof(int)); value_z0_uz = malloc( nnods_z0 * sizeof(double));
	    px = index_z0_ux; py = index_z0_uy; pz = index_z0_uz;
	    flag_pn=1;
	  }
	  break;
	case 8:
	  if(!strcmp(name,"Z1")){
	    flag=flag|(1<<8);
	    nnods_z1 = nnods;
	    index_z1_ux = malloc( nnods_z1 * sizeof(int)); value_z1_ux = malloc( nnods_z1 * sizeof(double));
	    index_z1_uy = malloc( nnods_z1 * sizeof(int)); value_z1_uy = malloc( nnods_z1 * sizeof(double));
	    index_z1_uz = malloc( nnods_z1 * sizeof(int)); value_z1_uz = malloc( nnods_z1 * sizeof(double));
	    px = index_z1_ux; py = index_z1_uy; pz = index_z1_uz;
	    flag_pn=1;
	  }
	  break;
	default:
	  break;
      }

      /*
	 Fill index_xx_xx arrays
      */
      pn = ((boundary_t*)pn->data)->Nods.head; n = 0;
      while(pn)
      {
	node_orig = *(int*)(pn->data); 
	p = bsearch(&node_orig, MyNodOrig, NMyNod, sizeof(int), cmpfunc); 
	if(!p){SETERRQ2(PROBLEM_COMM,1,
	    "A boundary node (%d) seems now to not belong to this process (rank:%d)",node_orig,rank);}
	node_petsc  = loc2petsc[p - MyNodOrig];  // PETSc numeration

	px[n] = node_petsc*3 + 0;
	py[n] = node_petsc*3 + 1;
	pz[n] = node_petsc*3 + 2;
	pn=pn->next; n ++;
      }

      pb=pb->next;
    }
    i++;
  }
  if(flag != 511)SETERRQ(MICRO_COMM,1, "MICRO:One entity not found P000 P100 P010 X0 X1 Y0 Y1 Z0 Z1.");
  return 0;
}
/****************************************************************************************************/
