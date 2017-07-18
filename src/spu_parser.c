/*

   Routines common for MICRO & MACRO for parsing the input 
   file searching for the keywords.

 */

#include "sputnik.h"

#define CHECK_INPUT_ERROR(data_char)                                                   \
     {if(!(data_char)){                                                                  \
	 printf("INPUT ERROR on %s line %d\n",__FILE__,__LINE__);                      \
	 return -1;                                                                    \
     }}

int spu_parse_scheme( char * input )
{

  /*

     Parse the communication scheme from input file

returns: 0 success
-1 failed
1 not found 

Searchs for keywords:

$scheme
...
$end_scheme

Examples:

scheme          MACRO_MICRO
nstruc_mic      <ns>
nproc_per_mic   <np1> <np1> <np1> ... <npns> 

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

    if(!strcmp(data,"$scheme")){

      if(fgets(buf,NBUF,file) == NULL){
	printf("scheme section incomplete at line %d\n",ln);
	printf("scheme keyword expected at line %d\n",ln);
	return -1;
      }
      ln ++;
      data = strtok(buf," \n");
      if(strcmp(data,"scheme")){
	printf("scheme keyword expected at line %d\n",ln);
	return -1;
      }
      data = strtok(NULL," \n");
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
  return 1;
}

/****************************************************************************************************/

int spu_parse_mesh( char * input )
{

  /*

     Parse the mesh name from input file

returns: 0 success
-1 failed
1 not found 

Searchs for keywords:

$mesh
mesh <mesh_file>
$end_mesh

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
  return 1;
}

/****************************************************************************************************/

int SpuParseMaterials(MPI_Comm *PROBLEM_COMM, char * input )
{

  /*

     Parse the materials of the problem

returns: 0 success
-1 failed
1 not found 

Searchs for keywords:

$materials
<PhysicalName> <TYPEXX> <options>
IRON TYPE00 E=1.0e6 v=1.0e6
$end_materials

   */

  FILE   *file = fopen(input,"r");
  char   buf[NBUF];
  char   *data;
  int    ln = 0;
  //    int    flag_end_material   = 0;
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

    if(!strcmp(data,"$materials")){

      flag_start_material=1;
      while(fgets(buf,NBUF,file) != NULL)
      {
	ln ++;

	// <name>
	data = strtok(buf," \n");
	if(!data){
	  printf("SpuParseMaterials: <name> expected\n");
	  return 1;
	}
	if(!strcmp(data,"$end_materials")) break;
	//	  strcpy(material.name,data);
	material.name = strdup(data);

	// <type> & <options>
	data = strtok(NULL," \n");
	if(!data){
	  printf("SpuParseMaterials: <name> expected\n");
	  return 1;
	}
	if(!strcmp(data,"TYPE00")){

	  material.typeID = TYPE00;
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
	  ((type_00*)material.type)->young = atof(&data[2]);

	  // lo insertamos en la lista 
	  list_insertlast(&material_list, &material);
	}
	else{
	  printf("SpuParseMaterials: %s unknown.\n", data);
	  return 1;
	}

      }
    } // inside $mesh

    if(!strcmp(data,"$end_materials")){
      CHECK_INPUT_ERROR(flag_start_material);
      PetscPrintf(*PROBLEM_COMM, "# of materials found in %s : %d\n", input, material_list.sizelist);
      return 0;
    }
  }
  // any material found 
  printf("SpuParseMaterials: Any material found on input file\n");
  return 0;
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
	PetscPrintf(*PROBLEM_COMM, "# of materials found in %s : %d\n", mesh_n, physical_list.sizelist);
	return 0;
      }
    }
    ln ++;
  }
  return 0;
}
