/*

   Routines common for MICRO & MACRO for parsing the input 
   file searching for the keywords.

*/

#include "sputnik.h"

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

int SpuParseMaterials( char * input )
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
	    data = strtok(NULL," \n");
	    if(!data){
	      printf("SpuParseMaterials: <E=<value>> expected\n");
	      return 1;
	    }
	    if(strncmp(data,"E=",2)){
	      printf("SpuParseMaterials: <E=<value>> expected\n");
	      return 1;
	    }
	    ((type_00*)material.type)->young = atof(&data[2]);

            // módulo de poisson
	    data = strtok(NULL," \n");
	    if(!data){
	      printf("SpuParseMaterials: <v=<value>> expected\n");
	      return 1;
	    }
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
	if(flag_start_material==0){
	  printf("format error at line %d\n",ln);
	  return -1;
	}
	return 0;
      }
    }
    // any material found 
    printf("SpuParseMaterials: Any material found on input file\n");
    return 0;
}
