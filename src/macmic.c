/*
   Common functions for <macro> & <micro> programs
 */

#include "macmic.h"

#define NBUF 64

int MacMicInitGaussStructure(int *eptr, int nelm)
{
  /*
     Allocates memory for Gauss point global structure <gauss>
     of type <gauss_t>
   */
  int i, ngp, ngauss = 0;
  for(i=1;i<(nelm+1);i++){
    // we assume that each element has the same number of gauss points as vertices 
    ngp = eptr[i] - eptr[i-1]; 
    ngauss += ngp; 
  }
  gauss = malloc(ngauss * sizeof(gauss_t)); if(!gauss) return 1;

  for(i=0;i<ngauss;i++){
    gauss[i].param_d = NULL;
  }

  return 0;
}

/****************************************************************************************************/

int MacMicParseScheme( char *input )
{

  /*
     Parse the communication scheme from input file

     Examples>

     $scheme
     scheme  COUP_1
     $end_scheme
   */

  FILE * file = fopen(input,"r");
  char   buf[NBUF];
  char * data;
  int    ln, i;


  if(!file) return 1;

  ln = 0;
  while(fgets(buf,NBUF,file) != NULL)
  {

    ln ++;
    data = strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$scheme")){

	fgets(buf,NBUF,file);
	ln ++;
	data = strtok(buf," \n"); 
	if(!data) return 1;
	if(strcmp(data,"scheme")) return 1;

	data = strtok(NULL," \n"); 
	if(!data) return 1;
	if(!strcmp(data,"COUP_1")){
	  macmic.type = COUP_1;
	  return 0;
	}
	else{
	  return -1;
	}

	return -1; //no encontro $end_scheme

      } // inside $scheme
    }

  }
  return 1;
}

