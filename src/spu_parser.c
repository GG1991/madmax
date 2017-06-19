/*

   Routines common for MICRO & MACRO for parsing the input 
   file searching for the keywords.

*/

#include "sputnik.h"

int spu_parse_scheme( char * input )
{

    /*

       Parse the communication scheme from input file

       Searchs for keywords:

       $scheme
         ...
       $end_scheme

       defines the global variables:

       int   scheme;
       int   nstruc_mic;
       int   nproc_per_mic[nstruc_mic];

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

      if(strcmp(data,"$scheme")){

	if(fgets(buf,NBUF,file) == NULL){
	  printf("scheme section incomplete at line %d\n",ln);
	  printf("scheme keyword expected at line %d\n",ln);
	  return 1;
	}
	ln ++;
	data = strtok(buf," \n");
	if(strcmp(data,"scheme")){
	  printf("scheme keyword expected at line %d\n",ln);
	  return 1;
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
	  return 1;
	}

	if(scheme==MACRO_MICRO){ // MACRO_MICRO scheme

	  if(fgets(buf,NBUF,file) == NULL){ // num_micro_structures
	    printf("scheme section incomplete at line %d\n",ln);
	    printf("num_micro_struct keyword expected at line %d\n",ln);
	    return 1;
	  }
	  ln ++;
	  data = strtok(buf," \n");
	  if(strcmp(data,"num_micro_struct")){
	    printf("nproc_per_mic keyword expected at line %d\n",ln);
	    return 1;
	  }
	  data = strtok(NULL," \n");
	  nstruc_mic = atoi(data);

	  if(fgets(buf,NBUF,file) == NULL){ // nproc_micro_structures
	    printf("scheme section incomplete at line %d\n",ln);
	    printf("nproc_per_mic keyword expected at line %d\n",ln);
	    return 1;
	  }
	  ln ++;

	  nproc_per_mic = malloc(nstruc_mic * sizeof(int));
	  data = strtok(buf," \n");
	  i = 0;
	  while(data){
	    nproc_per_mic[i] = atoi(data);
	    data = strtok(NULL," \n");
	    i++;
	  }
	  if(i!=nstruc_mic){
	    printf("nproc_per_mic given are more then the specified previouly: %d at line %d\n", nstruc_mic, ln);
	    return 1;
	  }

	}
	else if(scheme==MACRO_ALONE){ // MACRO_ALONE scheme
	  //TODO
	}
	else if(scheme==MICRO_ALONE){ // MICRO_ALONE scheme
	  //TODO
	}

	if(!strcmp(data,"$end_scheme")){
	  return 0;
	}
	return 1; //no encontro $end_scheme

      } // inside $scheme

    }

    return 0;
}

