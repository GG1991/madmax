/*

   Routines to create boundary structures and to impose boundary condition on RVE cells

   Author > Guido Giuntoli
   Date > 31-07-2017

*/

#include "micro.h"

int mic_parse_boundary(MPI_Comm PROBLEM_COMM, char *input)
{
  /*
     Parse boundary names

     $BoundaryMic
     <name1>
     <name2>
     ...
     $EndBoundaryMic
   */

  FILE        *file;
  char        buf[NBUF], *data;
  int         ln=0, flag_start_boundary=0;
  boundary_t  boundary;

  file = fopen(input,"r");
  if(!file){
    PetscPrintf(PETSC_COMM_WORLD,"input file %s not found.\n", input);
    return 1;
  }

  list_init(&boundary.Nods, sizeof(int), cmpfunc_for_list);
  while(fgets(buf,NBUF,file) != NULL)
  {

    ln ++;
    data = strtok(buf," \n");
    if(data){ 
      if(!strcmp(data,"$BoundaryMic")){

	flag_start_boundary=1;
	while(fgets(buf,NBUF,file) != NULL)
	{
	  ln ++;

	  // <name>
	  data = strtok(buf," \n"); 
	  if(!strcmp(data,"$EndBoundaryMic")) break;
	  boundary.name = strdup(data);

	  // si llegamos hasta acá esta todo 0K lo insertamos en la lista 
	  boundary.bvoid = malloc(sizeof(mic_boundary_t));
	  list_insertlast(&boundary_list, &boundary);
	}
      } 

      if(!strcmp(data,"$EndBoundaryMic")){
	if(!flag_start_boundary){
	  SETERRQ(PETSC_COMM_SELF,1,"$EndBoundaryMic found but not $Boundary");
	}
	return 0;
      }
    }
  }
  SETERRQ(PETSC_COMM_SELF,1,"any boundary condition found");
  return 1;
}
/****************************************************************************************************/
int voigt2mat(double voigt[6], double matrix[3][3])
{
  if(!voigt || !matrix) return 1;
  matrix[0][0] = voigt[0]; matrix[0][1] = voigt[3]; matrix[0][2] = voigt[5];
  matrix[1][0] = voigt[3]; matrix[1][1] = voigt[1]; matrix[1][2] = voigt[4];
  matrix[2][0] = voigt[5]; matrix[2][1] = voigt[4]; matrix[2][2] = voigt[2];
  return 0;
}
/****************************************************************************************************/
int mic_init_boundary_list(list_t *boundary_list)
{
  /*
     Creates the boundary list with names
     (¿ NO ?) P000 P100 P010 X0 X1 Y0 Y1 Z0 Z1 (at least)
     X0 X1 Y0 Y1 Z0 Z1 (at least)
   */
  int i = 0;
  boundary_t boundary;

  list_init(boundary_list, sizeof(boundary_t), NULL);
  list_init(&boundary.Nods,sizeof(int), cmpfunc_for_list);
  while(i<6)
  {
    switch(i){
      case 0:
	boundary.name = strdup("X0");break;
      case 1:
	boundary.name = strdup("X1");break;
      case 2:
	boundary.name = strdup("Y0");break;
      case 3:
	boundary.name = strdup("Y1");break;
      case 4:
	boundary.name = strdup("Z0");break;
      case 5:
	boundary.name = strdup("Z1");break;
      default:
	break;
    }
    list_insertlast(boundary_list, &boundary);
    i++;
  }
  return 0;
}
/****************************************************************************************************/
int micro_check_physical_entities( list_t *physical_list )
{
  /*
     Checks if the physical entities defined on mesh file are
     (¿ NO ?) P000 P100 P010 X0 X1 Y0 Y1 Z0 Z1 (at least)
     X0 X1 Y0 Y1 Z0 Z1 (at least)
   */
  int i = 0, flag = 0, flag_pn = 0;
  char *name;
  while(i<6)
  {
    node_list_t * pn = physical_list->head;
    flag_pn=0;
    while(pn && !flag_pn)
    {
      name = ((physical_t*)pn->data)->name;
      switch(i){
	case 0:
	  if(!strcmp(name,"X0")  ){flag=flag|(1<<0);flag_pn=1;}break;
	case 1:                                   
	  if(!strcmp(name,"X1")  ){flag=flag|(1<<1);flag_pn=1;}break;
	case 2:                                   
	  if(!strcmp(name,"Y0")  ){flag=flag|(1<<2);flag_pn=1;}break;
	case 3:                                   
	  if(!strcmp(name,"Y1")  ){flag=flag|(1<<3);flag_pn=1;}break;
	case 4:                                   
	  if(!strcmp(name,"Z0")  ){flag=flag|(1<<4);flag_pn=1;}break;
	case 5:                                   
	  if(!strcmp(name,"Z1")  ){flag=flag|(1<<5);flag_pn=1;}break;
	default:
	  break;
      }
      pn=pn->next;
    }
    i++;
  }
  if(flag!=63)SETERRQ(MICRO_COMM,1,"MICRO:physical entity not found (X0 X1 Y0 Y1 Z0 Z1)");
  //if(flag != 511)SETERRQ(MICRO_COMM,1, "MICRO:physical entity not found (P000 P100 P010 X0 X1 Y0 Y1 Z0 Z1)");
  return 0;
}
/****************************************************************************************************/
