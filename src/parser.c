/*

   Routines common for MICRO & MACRO for parsing the input 
   file searching for the keywords.

 */

#include "sputnik.h"

int parse_material(MPI_Comm PROBLEM_COMM, char * input)
{
  /*
     Parse the materials of the problem

     Example>

     $materials
     <PhysicalName> <TYPEXX> <options>
     IRON  TYPE_0 E=1.0e6 v=1.0e6
     MICRO TYPE01 E=1.0e6 v=1.0e6
     $end_materials
   */

  FILE         *file; 
  char         buf[NBUF], *data;
  int          ln=0;
  int          flag_start_material=0;
  material_t   material;

  file = fopen(input,"r");
  if(!file){
    PetscPrintf(PETSC_COMM_WORLD,"input file %s not found.\n", input);
    return 1;
  }

  list_init(&material_list, sizeof(material_t), NULL); 

  while(fgets(buf,NBUF,file) != NULL)
  {

    ln ++;
    data = strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$Materials")){

	flag_start_material=1;
	while(fgets(buf,NBUF,file) != NULL)
	{
	  ln ++;

	  // <name>
	  data = strtok(buf," \n"); 
	  if(!data) SETERRQ(PROBLEM_COMM,1,"<name> expected.");

	  if(data[0]!='#'){

	    if(!strcmp(data,"$EndMaterials")) break;
	    //	  strcpy(material.name,data);
	    material.name = strdup(data);

	    // <type> & <options>
	    data = strtok(NULL," \n");
	    if(!data) SETERRQ(PROBLEM_COMM,1,"<type> expected.");

	    if(!strcmp(data,"TYPE_0")){

	      material.type_id = TYPE_0;
	      material.GmshID = -1;
	      material.type = malloc(sizeof(type_0));

	      // densidad rho
	      data = strtok(NULL," \n");
	      if(!data)return 1;
	      if(strncmp(data,"rho=",4)){
		PetscPrintf(PETSC_COMM_WORLD,"<rho=<value>> expected\n", input);
		return 1;
	      }
	      if(data[4]=='\0') return 1;
	      ((type_0*)material.type)->rho = atof(&data[4]);

	      // módulo de young
	      data = strtok(NULL," \n");
	      if(!data)return 1;
	      if(strncmp(data,"E=",2)){
		PetscPrintf(PETSC_COMM_WORLD,"<E=<value>> expected\n", input);
		return 1;
	      }
	      if(data[2]=='\0') return 1;
	      ((type_0*)material.type)->young = atof(&data[2]);

	      // módulo de poisson
	      data = strtok(NULL," \n");
	      if(!data) return 1;
	      if(strncmp(data,"v=",2)){
		PetscPrintf(PETSC_COMM_WORLD,"<v=<value>> expected\n", input);
		return 1;
	      }
	      if(data[2]=='\0') return 1;
	      ((type_0*)material.type)->poisson = atof(&data[2]);

	      // calculamos parametros derivados
	      double E, v;
	      E = ((type_0*)material.type)->young;
	      v = ((type_0*)material.type)->poisson;
	      ((type_0*)material.type)->lambda = (E*v)/((1+v)*(1-2*v));
	      ((type_0*)material.type)->mu = E/(2*(1+v));
	    }
	    else if(!strcmp(data,"MICRO")){
	      material.type_id = MICRO;
	      material.GmshID = -1;
	      material.type = NULL;
	    }
	    else{
	      SETERRQ1(PROBLEM_COMM,1,"material type %s not valid.",data);
	    }
	    list_insertlast(&material_list, &material);
	  }
	}
      } 

      if(!strcmp(data,"$EndMaterials")){
	if(!flag_start_material){
	  return 1;
	}
	fclose(file);
	return 0;
      }
    } 
  }
  return 1;
}
/****************************************************************************************************/
int is_linear_material(int material)
{
  switch(material){
    case TYPE_0:
      return 1;
  }
  return 0;
}
/****************************************************************************************************/
int check_elm_id(void)
{
  /* 
     Checks if all the elements of <elm_id>
     have their correspondent material in <material_list>
   */
  int e;

  node_list_t *pm;

  for(e=0;e<nelm;e++){
    pm = material_list.head;
    while(pm){
      if( ((material_t*)pm->data)->GmshID == elm_id[e] ) break;
      pm = pm->next;
    }
    if(!pm){
      return 1;
    }
  }

  return 0;
}
/****************************************************************************************************/
int parse_function(MPI_Comm PROBLEM_COMM, char *input )
{
  /*
     Parse the functions of the problem
     
     $Function
     <fnum> <inter> <n>
     x1 y1
     x2 y2
     $EndFunction
   */

  FILE   *file;
  char   buf[NBUF], *data;
  int    ln=0, n;
  int    flag_start_function=0;
  f1d_t  f1d;

  file = fopen(input,"r");
  if(!file){
    PetscPrintf(PETSC_COMM_WORLD,"input file %s not found.\n", input);
    return 1;
  }

  list_init(&function_list, sizeof(f1d_t), NULL);

  while(fgets(buf,NBUF,file) != NULL)
  {

    ln ++;
    data = strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$Function")){
        
	if(flag_start_function)return 1;
	flag_start_function=1;
	fgets(buf,NBUF,file); ln ++;

	// <fnum>
	data = strtok(buf," \n");
	if(!data)return 1;
	f1d.fnum = atoi(data);

	// <inter> 
	data = strtok(NULL," \n");
	if(!data)return 1;

	if(!strcmp(data,"INTER1")){
	    f1d.inter = INTER1 ;
	}
	else{
	  return 1;
	}

	// <n>
	data = strtok(NULL," \n");
	if(!data)SETERRQ1(PROBLEM_COMM,1,"format error on %s",input);
	f1d.n = atoi(data);
	f1d.x = malloc(f1d.n*sizeof(double));
	f1d.y = malloc(f1d.n*sizeof(double));
        
	n = 0;
	while( fgets(buf,NBUF,file) != NULL ){
	  if(n >= f1d.n) break;

	  data = strtok(buf," \n");
	  if(!data)return 1;
	  f1d.x[n] = atof(data);

	  data = strtok(NULL," \n");
	  if(!data)return 1;
	  f1d.y[n] = atof(data);
	  n++;
	}
	data = strtok(buf," \n");
	if(!data)return 1;

	if(strcmp(data,"$EndFunction"))return 1;
	list_insertlast(&function_list, &f1d);
      } 

      if(!strcmp(data,"$EndFunction")){
	if(!flag_start_function)SETERRQ1(PROBLEM_COMM,1,"format error on %s",input);
	flag_start_function = 0;
      }
    }
  }
  return 0;
}
/****************************************************************************************************/
int StrBin2Dec(char *str)
{
  /* Converts the string in "str" that suppose to have a continuue
     sequence of "0" and "1" to its decimal representation.
   */

  int i;
  int dec = 0; 

  for(i=strlen(str)-1;i>=0;i--){
    if(str[i]=='0' || str[i]=='1'){
      if(str[i]=='1')
	dec+=(int)pow(2,strlen(str)-1-i);
    }else{
      return -1;
    }
  }
  return dec;
}
/****************************************************************************************************/
int isfloat(char *s)
{
  int i = 0;

  if(!s) return 1;

  if( s[i] == '\0' ) return 1;

  while( s[i] != '\0' && i < strlen(s) ) {

    if (s[i] != 'e' ||
	s[i] != 'E' ||
	s[i] != '.' ||
	s[i] != '0' ||
	s[i] != '1' ||
	s[i] != '2' ||
	s[i] != '3' ||
	s[i] != '4' ||
	s[i] != '5' ||
	s[i] != '6' ||
	s[i] != '7' ||
	s[i] != '8' ||
	s[i] != '9' )
      return 1;

    i++;
  }

  return 0;
}
