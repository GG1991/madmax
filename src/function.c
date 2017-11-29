#include "function.h"

/****************************************************************************************************/

int function_init(double *x, double *y, int n, int inter, f1d_t *f1d){

  if( n == 0 || x == NULL || y == NULL || f1d == NULL )
    return 1;

  f1d->n     = n;
  f1d->inter = inter;
  f1d->fnum  = -1;
  f1d->x     = calloc(n,sizeof(double));
  f1d->y     = calloc(n,sizeof(double));
  if( !f1d->x || !f1d->y )
    return 1;
  int i;
  for( i = 0 ; i < n ; i++ ){
    f1d->x[i] = x[i];
    f1d->y[i] = y[i];
  }
  return 0;
}

/****************************************************************************************************/

int function_eval(double x, f1d_t *f1d, double *y)
{


  if( f1d == NULL ) return 1;

  if( f1d->n < 1 ) return 1;

  if( f1d->n == 1 ){
    *y=f1d->y[0];
    return 0;
  }
  if( x < f1d->x[0] ){
    *y = f1d->y[0];
    return 0;
  }
  if( x > f1d->x[1] ){
    *y = f1d->y[1];
    return 0;
  }

  int i = 1;
  while( i < f1d->n ){
    if( f1d->x[i-1] <= x && x < f1d->x[i] )
      break;
    i++;
  }
  if( i == f1d->n ){
    *y=f1d->y[i-1];
    return 0;
  }

  *y = ( f1d->y[i] - f1d->y[i-1] )*( x - f1d->x[i-1] )/( f1d->x[i] - f1d->x[i-1] ) + f1d->y[i-1];

  return 0;
}

/****************************************************************************************************/

int function_comp(void *a, void *b)
{
  if ( ((f1d_t *)a)->fnum > ((f1d_t *)b)->fnum )
    return 1;
  else if( ((f1d_t*)a)->fnum == ((f1d_t*)b)->fnum )
    return 0;
  else
    return -1;
  return 1;
}

/****************************************************************************************************/

int function_get_from_list(int fn, list_t *function_list, f1d_t **f1d)
{

  /* returns 0 if was found and 1 if not or error */

   if( function_list == NULL ){
     f1d = NULL;
     return 1;
   }
   if( function_list->sizelist == 0 ){
     f1d = NULL;
     return 1;
   }

   node_list_t * pn = function_list->head;
   f1d_t * f1d_a;
   while(pn)
   {
     f1d_a = (f1d_t *)pn->data;
     if( f1d_a->fnum == fn ) break;
     pn = pn->next;
   }
   if( pn == NULL ){
     f1d = NULL;
     return 1;
   }
   *f1d = f1d_a;

   return 0;
}

/****************************************************************************************************/

int function_fill_list_from_command_line(int argc, const char **argv, list_t *function_list)
{

  list_init(function_list, sizeof(f1d_t), NULL);

  char **string_array, *data;
  int    found, n_str_found;
  f1d_t  fun;
  found = myio_get_string_array_command_line(argc, argv, "-function", MAX_NUM_OF_FUNCTIONS, &string_array, &n_str_found);

  if(found || !n_str_found)
    return 1;

  int i, j;
  for( i = 0 ; i < n_str_found ; i++ )
  {
    data     = strtok(string_array[i]," \n");
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
