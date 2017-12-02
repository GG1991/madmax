#include "function.h"

/****************************************************************************************************/

int function_init(double *x, double *y, int n, int inter, function_t *function){

  if( n == 0 || x == NULL || y == NULL || function == NULL )
    return 1;

  function->n     = n;
  function->inter = inter;
  function->fnum  = -1;
  function->x     = calloc(n,sizeof(double));
  function->y     = calloc(n,sizeof(double));
  if( !function->x || !function->y )
    return 1;
  int i;
  for( i = 0 ; i < n ; i++ ){
    function->x[i] = x[i];
    function->y[i] = y[i];
  }
  return 0;
}

/****************************************************************************************************/

int function_eval(double x, function_t *function, double *y)
{


  if( function == NULL ) return 1;

  if( function->n < 1 ) return 1;

  if( function->n == 1 ){
    *y=function->y[0];
    return 0;
  }
  if( x < function->x[0] ){
    *y = function->y[0];
    return 0;
  }
  if( x > function->x[1] ){
    *y = function->y[1];
    return 0;
  }

  int i = 1;
  while( i < function->n ){
    if( function->x[i-1] <= x && x < function->x[i] )
      break;
    i++;
  }
  if( i == function->n ){
    *y=function->y[i-1];
    return 0;
  }

  *y = ( function->y[i] - function->y[i-1] )*( x - function->x[i-1] )/( function->x[i] - function->x[i-1] ) + function->y[i-1];

  return 0;
}

/****************************************************************************************************/

int function_comp(void *a, void *b)
{
  if ( ((function_t *)a)->fnum > ((function_t *)b)->fnum )
    return 1;
  else if( ((function_t*)a)->fnum == ((function_t*)b)->fnum )
    return 0;
  else
    return -1;
  return 1;
}

/****************************************************************************************************/

int function_get_from_list(int fn, list_t *function_list, function_t **function)
{

  /* returns 0 if was found and 1 if not or error */

   if( function_list == NULL ){
     function = NULL;
     return 1;
   }
   if( function_list->sizelist == 0 ){
     function = NULL;
     return 1;
   }

   node_list_t * pn = function_list->head;
   function_t * function_a;
   while(pn)
   {
     function_a = (function_t *)pn->data;
     if( function_a->fnum == fn ) break;
     pn = pn->next;
   }
   if( pn == NULL ){
     function = NULL;
     return 1;
   }
   *function = function_a;

   return 0;
}

/****************************************************************************************************/

int function_fill_list_from_command_line(command_line_t *command_line, list_t *function_list){

  myio_get_string_array_command_line(command_line, MAX_NUM_OF_FUNCTIONS, "-function");

  if(!command_line->found || !command_line->n_str_found)
    return 1;

  list_init(function_list, sizeof(function_t), NULL);
  char *data;

  int i, j;
  function_t fun;
  for( i=0 ; i<command_line->n_str_found ; i++ ){
    data     = strtok(command_line->str_arr[i]," \n");
    fun.fnum = atoi(data);
    data     = strtok(NULL, " \n");
    fun.n    = atoi(data);
    fun.x    = malloc(fun.n*sizeof(double));
    fun.y    = malloc(fun.n*sizeof(double));
    for( j = 0 ; j < fun.n ; j++ ){
      data = strtok(NULL," \n"); fun.x[j] = atof(data);
      data = strtok(NULL," \n"); fun.y[j] = atof(data);
    }
    list_insertlast(function_list, &fun);
  }

  return 0;
}

/****************************************************************************************************/
