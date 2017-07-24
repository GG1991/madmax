/*

   Routines common for MICRO & MACRO for parsing the input 
   file searching for the keywords.

 */

#include "sputnik.h"

#define CHECK_INPUT_ERROR(data_char)                                                   \
     {if(!(data_char)){                                                                \
	 printf("INPUT ERROR on %s line %d\n",__FILE__,__LINE__);                      \
	 return -1;                                                                    \
     }}

int spu_parse_scheme( char * input )
{

  /*
   * 
   *  Parse the communication scheme from input file
   * 
   *  returns: 0 success
   *  -1 failed
   *  1 not found 
   *  
   *  Searchs for keywords:
   *  
   *  $scheme
   *  ...
   *  $end_scheme
   *  
   *  Examples:
   *  
   *  scheme          MACRO_MICRO
   *  nstruc_mic      <ns>
   *  nproc_per_mic   <np1> <np1> <np1> ... <npns> 
   * 
   */

  FILE * file = fopen(input,"r");
  char   buf[NBUF];
  char * data;
  int    ln, i;


  nproc_per_mic = NULL;	
  nstruc_mic    = -1;
  scheme        = -1;

  if(!file){
    return 1;
  }

  ln = 0;
  while(fgets(buf,NBUF,file) != NULL)
  {

    ln ++;
    data = strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$scheme")){

	if(fgets(buf,NBUF,file) == NULL){
	  printf("scheme section incomplete at line %d\n",ln);
	  printf("scheme keyword expected at line %d\n",ln);
	  return -1;
	}
	ln ++;
	data = strtok(buf," \n"); CHECK_INPUT_ERROR(data);
	if(strcmp(data,"scheme")){
	  printf("scheme keyword expected at line %d\n",ln);
	  return -1;
	}
	data = strtok(NULL," \n"); CHECK_INPUT_ERROR(data);
	if(!strcmp(data,"MACRO_MICRO")){
	  scheme = MACRO_MICRO;
	}
	else if(!strcmp(data,"MACRO_ALONE")){
	  scheme = MACRO_ALONE;
	}
	else if(!strcmp(data,"MACRO_ALONE")){
	  scheme = MACRO_ALONE;
	}
	else{
	  printf("scheme type %s not valid at line %d\n",data,ln);
	  return -1;
	}

	if(scheme==MACRO_MICRO){ // MACRO_MICRO scheme

	  if(fgets(buf,NBUF,file) == NULL){ // num_micro_structures
	    printf("scheme section incomplete at line %d\n",ln);
	    printf("num_micro_struct keyword expected at line %d\n",ln);
	    return -1;
	  }
	  ln ++;
	  data = strtok(buf," \n");
	  if(strcmp(data,"num_micro_struct")){
	    printf("nproc_per_mic keyword expected at line %d\n",ln);
	    return -1;
	  }
	  data = strtok(NULL," \n");
	  nstruc_mic = atoi(data);

	  if(fgets(buf,NBUF,file) == NULL){ // nproc_micro_structures
	    printf("scheme section incomplete at line %d\n",ln);
	    printf("nproc_per_mic keyword expected at line %d\n",ln);
	    return -1;
	  }
	  ln ++;
	  data = strtok(buf," \n");
	  if(strcmp(data,"nproc_per_mic")){
	    printf("nproc_per_mic keyword expected at line %d\n",ln);
	    return -1;
	  }

	  nproc_per_mic = malloc(nstruc_mic * sizeof(int));
	  data = strtok(NULL," \n");
	  nproc_mic_group = i = 0;
	  while(data){
	    nproc_per_mic[i] = atoi(data);
	    nproc_mic_group += nproc_per_mic[i];
	    data = strtok(NULL," \n");
	    i++;
	  }
	  if(i!=nstruc_mic){
	    printf("nproc_per_mic given are more then the specified previouly: %d at line %d\n", nstruc_mic, ln);
	    return -1;
	  }

	  if(fgets(buf,NBUF,file) == NULL){ // nproc_micro_structures
	    printf("scheme section incomplete at line %d\n",ln);
	    printf("$end_scheme keyword expected at line %d\n",ln);
	    return -1;
	  }
	  data = strtok(buf," \n");
	  if(!strcmp(data,"$end_scheme")){
	    return 0;
	  }

	}

	return -1; //no encontro $end_scheme

      } // inside $scheme
    }

  }
  return 1;
}

/****************************************************************************************************/

int spu_parse_mesh( char * input )
{

  /*
   *   Parse the mesh name from input file
   * 
   *   returns: 0 success
   *   -1 failed
   *   1 not found 
   *   
   *   Searchs for keywords:
   *   
   *   $mesh
   *   mesh <mesh_file>
   *   $end_mesh
   * 
   */

  FILE * file = fopen(input,"r");
  char   buf[NBUF];
  char * data;
  int    ln;

  if(!file){
    return 1;
  }

  ln = 0;
  while(fgets(buf,NBUF,file) != NULL)
  {

    ln ++;
    data = strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$mesh")){

	if(fgets(buf,NBUF,file) == NULL){
	  printf("mesh section incomplete at line %d\n",ln);
	  printf("mesh keyword expected at line %d\n",ln);
	  return -1;
	}
	ln ++;
	data = strtok(buf," \n");
	if(strcmp(data,"mesh")){
	  printf("mesh keyword expected at line %d\n",ln);
	  return -1;
	}
	data = strtok(NULL," \n");
	strcpy(mesh_n,data);

	if(fgets(buf,NBUF,file) == NULL){ 
	  printf("mesh section incomplete at line %d\n",ln);
	  printf("$end_mesh keyword expected at line %d\n",ln);
	  return -1;
	}
	data = strtok(buf," \n");
	if(!strcmp(data,"$end_mesh")){
	  return 0;
	}
	return -1; //no se encontro $end_mesh

      } // inside $mesh
    }
  }
  return 1;
}

/****************************************************************************************************/

int SpuParseMaterials(MPI_Comm *PROBLEM_COMM, char * input )
{

  /*
   * Parse the materials of the problem
   * 
   *    returns: 0 success
   *            -1 failed
   *             1 not found 
   *    
   *    Searchs for keywords:
   *    
   *    $materials
   *    <PhysicalName> <TYPEXX> <options>
   *    IRON TYPE00 E=1.0e6 v=1.0e6
   *    $end_materials
   * 
   */

  FILE   *file = fopen(input,"r");
  char   buf[NBUF];
  char   *data;
  int    ln = 0;
  int    flag_start_material = 0;

  material_t material;

  if(!file){
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
	  if(!data) SETERRQ(PETSC_COMM_SELF,1,"SpuParseMaterials: <name> expected.");

	  if(data[0]!='#'){

	    if(!strcmp(data,"$EndMaterials")) break;
	    //	  strcpy(material.name,data);
	    material.name = strdup(data);

	    // <type> & <options>
	    data = strtok(NULL," \n");
	    if(!data) SETERRQ(PETSC_COMM_SELF,1,"SpuParseMaterials: <type> expected.");

	    if(!strcmp(data,"TYPE00")){

	      material.typeID = TYPE00;
	      material.GmshID = -1;
	      material.type = malloc(sizeof(type_00));

	      // módulo de young
	      data = strtok(NULL," \n");CHECK_INPUT_ERROR(data);
	      if(strncmp(data,"E=",2)){
		printf("SpuParseMaterials: <E=<value>> expected\n");
		return 1;
	      }
	      ((type_00*)material.type)->young = atof(&data[2]);

	      // módulo de poisson
	      data = strtok(NULL," \n");CHECK_INPUT_ERROR(data);
	      if(strncmp(data,"v=",2)){
		printf("SpuParseMaterials: <v=<value>> expected\n");
		return 1;
	      }
	      ((type_00*)material.type)->poisson = atof(&data[2]);

	      // calculamos parametros derivados
	      double E, v;
	      E = ((type_00*)material.type)->young;
	      v = ((type_00*)material.type)->poisson;
	      ((type_00*)material.type)->lambda = (E*v)/((1+v)*(1-2*v));
	      ((type_00*)material.type)->mu = E/(2*(1+v));

	      // lo insertamos en la lista 
	      list_insertlast(&material_list, &material);
	    }
	    else{
	      printf("SpuParseMaterials: %s unknown.\n", data);
	      return 1;
	    }

	  }
	}
      } // inside $Materials

      if(!strcmp(data,"$EndMaterials")){
	if(!flag_start_material)SETERRQ(PETSC_COMM_SELF,1,"$EndMaterials detected without $Materials above.");
	PetscPrintf(*PROBLEM_COMM, "# of materials found in %s : %d\n", input, material_list.sizelist);
	fclose(file);
	return 0;
      }
    } // data != NULL
  }
  // any $Material found
  SETERRQ(PETSC_COMM_SELF,1,"$Materials section not found on input file."); 
}

/****************************************************************************************************/

int SpuParsePhysicalEntities( MPI_Comm *PROBLEM_COMM, char *mesh_n )
{

  /* 
   * Info:   Reads the physical entities from the GMSH file
   *         and save them on <physical_list>
   *
   * Input: 
   * char   * mesh_n   : file name with path
   * 
   * Output:
   * list_t   physical_list : physical entities list
   *
   */

  FILE                 *fm;

  int                  i, ntot; 
  int                  ln = 0;  // line counter and offset for moving faster in the file
  int                  flag_start_physical;

  char                 buf[NBUF];   
  char                 *data;

  physical_t physical;

  fm = fopen(mesh_n,"r");
  if(!fm){
    printf("file not found : %s\n",mesh_n);
    return 1;
  }


  /**************************************************/
  //  go to "$PhysicalNames" and then read them
  //  filling physical_list
  //
  while(fgets(buf,NBUF,fm)!=NULL)
  {
    ln++;
    data=strtok(buf," \n");
    //
    // leemos hasta encontrar $PhysicalNames
    //
    if(strcmp(data,"$PhysicalNames")==0){
      //
      // leemos la cantidad de physical entities
      // solo para hacer verificaciones
      //
      flag_start_physical=1;
      fgets(buf,NBUF,fm);
      data  = strtok(buf," \n");
      ntot = atoi(data);

      while(fgets(buf,NBUF,fm)!=NULL)
      {
	ln++;

	// dimension
	data=strtok(buf," \n");
	if(!data){
	  printf("format error at line %d on %s\n", ln, mesh_n);
	  return -1;
	}
	if(!strcmp(data,"$EndPhysicalNames")) break;
	physical.dim = atoi(data);

	// GmshID
	data=strtok(NULL," \n");CHECK_INPUT_ERROR(data);
	physical.GmshID = atoi(data);

	// name (le sacamos los \" de las puntas)
	data=strtok(NULL," \n");CHECK_INPUT_ERROR(data);
	physical.name = malloc((strlen(data)-1)*sizeof(char)); 
	for(i=1;i<strlen(data)-1;i++){
	  physical.name[i-1] = data[i];
	}
	physical.name[i-1] = '\0';

	list_insertlast(&physical_list,&physical);

      }
      if(!strcmp(data,"$EndPhysicalNames")){
	CHECK_INPUT_ERROR(flag_start_physical);
	CHECK_INPUT_ERROR(physical_list.sizelist == ntot);
	PetscPrintf(*PROBLEM_COMM, "# of physical found in %s : %d\n", mesh_n, physical_list.sizelist);
	return 0;
      }
    }
    ln ++;
  }
  return 0;
}

/****************************************************************************************************/

int SetGmshIDOnMaterialsAndBoundaries(void)
{

  /* For each material on <material_list> 
   * Searchs for the <GmshID> in the 
   * <physical_list>
   */

   node_list_t *pm, *pp;

   pm = material_list.head;
   while(pm){
     pp = physical_list.head;
     while(pp){
       if( !strcmp( ((physical_t*)pp->data)->name, ((material_t*)pm->data)->name ) ){
	 ((material_t*)pm->data)->GmshID = ((physical_t*)pp->data)->GmshID;
	 break;
       }
       pp = pp->next;
     }
     if(!pp){ 
       SETERRQ1(PETSC_COMM_SELF,1,"Material %s not found in Gmsh File.",((material_t*)pm->data)->name);
     }
     pm = pm->next;
   }

   pm = boundary_list.head;
   while(pm){
     pp = physical_list.head;
     while(pp){
       if( !strcmp( ((physical_t*)pp->data)->name, ((boundary_t*)pm->data)->name ) ){
	 ((boundary_t*)pm->data)->GmshID = ((physical_t*)pp->data)->GmshID;
	 break;
       }
       pp = pp->next;
     }
     if(!pp){ 
       SETERRQ1(PETSC_COMM_SELF,1,"Boundary %s not found in Gmsh File.",((boundary_t*)pm->data)->name);
     }
     pm = pm->next;
   }

  return 0;
}

/****************************************************************************************************/

int CheckPhysicalID(void)
{

  /* Checks if all the elements of <PhysicalID>
   * have their correspondent material in <material_list>
   */

  int e;

  node_list_t *pm;

  for(e=0;e<nelm;e++){
    pm = material_list.head;
    while(pm){
      if( ((material_t*)pm->data)->GmshID == PhysicalID[e] ) break;
      pm = pm->next;
    }
    if(!pm) return 1;
  }

  return 0;
}

/****************************************************************************************************/

int SpuParseBoundary(MPI_Comm *PROBLEM_COMM, char *input )
{

  /*
   * Parse the boundary of the problem
   * 
   *    returns: 0 success
   *             1 failed
   *    
   *    Searchs for keywords:
   *    
   *    $Boundary
   *    <name1> <order> <kind> <fnumx> <fnumy> <fnumz>
   *    <name2> <order> <kind> <fnumx> <fnumy> <fnumz>
   *    ...
   *    $EndBoundary
   * 
   */

  FILE   *file = fopen(input,"r");
  char   buf[NBUF];
  char   *data;
  int    ln = 0;
  int    flag_start_boundary = 0;

  if(!file){
    return 1;
  }

  boundary_t boundary;
  list_init(&boundary_list, sizeof(boundary_t), cmpfuncBou);

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
	  data = strtok(buf," \n"); CHECK_INPUT_ERROR(data);
	  if(!strcmp(data,"$EndBoundary")) break;
	  boundary.name = strdup(data);

	  // <order> 
	  data = strtok(NULL," \n"); CHECK_INPUT_ERROR(data);
	  boundary.order = atoi(data);

	  // <kind> 
	  data = strtok(NULL," \n"); CHECK_INPUT_ERROR(data);
	  boundary.kind = atoi(data);

	  // <nfz> 
	  data = strtok(NULL," \n"); CHECK_INPUT_ERROR(data);
	  boundary.nfz = atoi(data);

	  // <nfy> 
	  data = strtok(NULL," \n"); CHECK_INPUT_ERROR(data);
	  boundary.nfy = atoi(data);

	  // <nfx> 
	  data = strtok(NULL," \n"); CHECK_INPUT_ERROR(data);
	  boundary.nfx = atoi(data);

	  boundary.fx = GetFunctionPointer(&function_list,boundary.nfx);CHECK_INPUT_ERROR(boundary.fx);
	  boundary.fy = GetFunctionPointer(&function_list,boundary.nfy);CHECK_INPUT_ERROR(boundary.fy);
	  boundary.fz = GetFunctionPointer(&function_list,boundary.nfz);CHECK_INPUT_ERROR(boundary.fz);

	  // si llegamos hasta acá esta todo 0K lo insertamos en la lista 
	  list_insert_se(&boundary_list, &boundary);
	}
      } // inside $Boundary

      if(!strcmp(data,"$EndBoundary")){
	CHECK_INPUT_ERROR(flag_start_boundary);
	PetscPrintf(*PROBLEM_COMM, "# of boundaries found in %s : %d\n", input, boundary_list.sizelist);
	return 0;
      }
    }
  }
  // any boundary condition found
  printf("SpuParseBoundary: Any boundary found on input file\n");
  return 1;
}

/****************************************************************************************************/

int SpuParseFunctions(MPI_Comm *PROBLEM_COMM, char *input )
{

  /*
   * Parse the functions of the problem
   * 
   *    returns: 0 success
   *             1 failed
   *    
   *    Searchs for keywords:
   *    
   *    $Function
   *    <fnum> <inter> <n>
   *    x1 y1
   *    x2 y2
   *    ...
   *    xn yn
   *    $EndFunction
   * 
   */

  FILE   *file = fopen(input,"r");
  char   buf[NBUF];
  char   *data;
  int    ln = 0, n;
  int    flag_start_function = 0;

  if(!file){
    return 1;
  }

  f1d_t f1d;
  list_init(&function_list, sizeof(f1d_t), NULL);

  while(fgets(buf,NBUF,file) != NULL)
  {

    ln ++;
    data = strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$Function")){
        
	if(flag_start_function) CHECK_INPUT_ERROR(NULL);
	flag_start_function=1;
	CHECK_INPUT_ERROR(fgets(buf,NBUF,file)); ln ++;

	// <fnum>
	data = strtok(buf," \n"); CHECK_INPUT_ERROR(data);
	f1d.fnum = atoi(data);

	// <inter> 
	data = strtok(NULL," \n"); CHECK_INPUT_ERROR(data);
	if(!strcmp(data,"INTER1")){
	    f1d.inter = INTER1 ;
	}
	else{
	  return 1;
	}

	// <n>
	data = strtok(NULL," \n"); CHECK_INPUT_ERROR(data);
	f1d.n = atoi(data);
	f1d.x = malloc(f1d.n*sizeof(double));
	f1d.y = malloc(f1d.n*sizeof(double));
        
	n = 0;
	while( fgets(buf,NBUF,file) != NULL ){
	  if(n >= f1d.n) break;
	  data = strtok(buf," \n"); CHECK_INPUT_ERROR(data); 
	  f1d.x[n] = atof(data);
	  data = strtok(NULL," \n"); CHECK_INPUT_ERROR(data); 
	  f1d.y[n] = atof(data);
	  n++;
	}
	data = strtok(buf," \n"); CHECK_INPUT_ERROR(data); 
	if(strcmp(data,"$EndFunction"))return 1;
	// si llegamos hasta acá esta todo 0K lo insertamos en la lista 
	list_insertlast(&function_list, &f1d);
      } // inside $Function

      if(!strcmp(data,"$EndFunction")){
	CHECK_INPUT_ERROR(flag_start_function);
	flag_start_function = 0;
      }
    }
  }
  PetscPrintf(*PROBLEM_COMM, "# of functions found in %s : %d\n", input, function_list.sizelist);
  return 0;
}

/****************************************************************************************************/


int cmpfuncBou (void * a, void * b)
{
  return ( ((boundary_t *)a)->order - ((boundary_t *)b)->order );
}
