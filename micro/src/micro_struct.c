/*
 * micro_struct.c - functions to represent the geometrical micro-structure
 * pattern of a solid. The description of what each micro-structure 
 * represents can be found on micro_struc.h.
 * 
 *
 * author : Guido Giuntoli
 * date   : 26 - 11 - 2017
 *
 */

#include "micro_struct.h"

/**********************************************************************/

int micro_struct_init( int dim, const char *string, micro_struct_t *micro_struct )
{

  /*
   * Initializates the micro_struct parsing "string" the first word
   * indicates the kind of "micro_struct" type. Then a convention 
   * is implemented to filled all the structure.
   *
   *
   * string = { "type" , "parameters" }
   *
   * return 0 if success
   * return 1 if non-specific error
   * return 2 if format error
   *
   */

  char   *stra = strdup( string );
  char   *data = strtok( stra, " \n" );
  double *size = malloc( dim * sizeof(double) );
  int     d;

  if( !strcmp(data, "fiber_cilin") )
  {

    /* Reads the format
     * "size"[dim] "nx_fib" "ny_fib" "radio" "desv"[2]
     */

    fiber_cilin_t *fiber_cilin = malloc(sizeof(fiber_cilin_t));

    /* read and write */
    fiber_cilin->desv = malloc(dim*sizeof(double));

    for( d=0 ; d<dim ; d++ )
    {
      data = strtok( NULL, " \n" ); 
      if(!data) return 2;
      size[d] = atof( data );
    }

    data = strtok( NULL, " \n" );
    if(!data) return 2;
    fiber_cilin->nx_fib = atoi( data );

    data = strtok( NULL, " \n" );
    if(!data) return 2;
    fiber_cilin->ny_fib = atoi( data );

    data = strtok( NULL, " \n" );
    if(!data) return 2;
    fiber_cilin->radio = atof( data );

    for( d=0 ; d<2 ; d++ )
    {
      data = strtok( NULL, " \n" );
      if(!data) return 2;
      fiber_cilin->desv[d] = atof( data );
    }

    /* assign to micro_struct */
    micro_struct->type = FIBER_CILIN;
    micro_struct->data = fiber_cilin;

  }
  else if( !strcmp(data, "fiber_line") )
  {

    /* Reads the format
     * "size"[dim] "ntype" "nfib[ntype]" "tetha[ntype]"  "seps[ntype]" "width[ntype]" "desv[ntype]"
     */

    int ntype;

    fiber_line_t *fiber_line = malloc(sizeof(fiber_line_t));

    /* read and write */
    data = strtok( NULL, " \n" );
    if(!data) return 2;
    fiber_line->ntype = ntype = atoi( data );

    fiber_line->theta = malloc(ntype*sizeof(double));
    fiber_line->seps  = malloc(ntype*sizeof(double));
    fiber_line->width = malloc(ntype*sizeof(double));
    fiber_line->desv  = malloc(ntype*sizeof(double));
    fiber_line->nfib  = malloc(ntype*sizeof(double));

    for( d=0 ; d<dim ; d++ )
    {
      data = strtok( NULL, " \n" );
      if(!data) return 2;
      size[d] = atof( data );
    }

    for( d=0 ; d<ntype ; d++ )
    {
      data = strtok( NULL, " \n" );
      if(!data) return 2;
      fiber_line->theta[d] = atof( data );
    }

    for( d=0 ; d<ntype ; d++ )
    {
      data = strtok( NULL, " \n" );
      if(!data) return 2;
      fiber_line->seps[d] = atoi( data );
    }

    for( d=0 ; d<ntype ; d++ )
    {
      data = strtok( NULL, " \n" );
      if(!data) return 2;
      fiber_line->width[d] = atof( data );
    }

    for( d=0 ; d<ntype ; d++ )
    {
      data = strtok( NULL, " \n" );
      if(!data) return 2;
      fiber_line->desv[d] = atof( data );
    }

    /* assign to micro_struct */
    micro_struct->type = FIBER_LINE;
    micro_struct->data = fiber_line;

  }
  else
  {
    return 1;
  }

  micro_struct->size = size;

  return 0;
}

/**********************************************************************/

int micro_struct_get_elem_id( int dim, micro_struct_t *micro_struct, double *elem_centroid, int *elem_id )
{

  /* returns elem_id as a function of the centroid coordinate */

  int    i, j, d;

  if( micro_struct->type == FIBER_CILIN )
  {

      fiber_cilin_t * fiber_cilin = (fiber_cilin_t *)micro_struct->data;

      double deviation[2];
      double center[2];
      double lx = micro_struct->size[0];
      double ly = micro_struct->size[1];
      center[0] = lx / 2;
      center[1] = ly / 2;

      /* as default is in the matrix */
      *elem_id = ID_MATRIX;

      /* check if it is inside one of the fibers */
      for( i = 0 ; i < fiber_cilin->nx_fib ; i++ )
      {
	for( j = 0 ; j < fiber_cilin->ny_fib ; j++ )
	{
	  deviation[0] = fiber_cilin->desv[0] - lx/2 + (lx/fiber_cilin->nx_fib)/2 + i * (lx/fiber_cilin->nx_fib);
	  deviation[1] = fiber_cilin->desv[1] - ly/2 + (ly/fiber_cilin->ny_fib)/2 + j * (ly/fiber_cilin->ny_fib);
	  double l = 0.0;
	  for( d = 0 ; d < 2 ; d++ ){
	    l = l + pow( elem_centroid[d] - (center[d] + deviation[d]), 2 );
	  }
	  l = sqrt(l);
	  if( l <= fiber_cilin->radio ){
	    *elem_id = ID_FIBER ;
	    return 0;
	  }

	}
      }

  }
  else if( micro_struct->type == FIBER_LINE )
  {

    fiber_line_t * fiber_line = (fiber_line_t *)micro_struct->data;

    /* as default is in the matrix */
    *elem_id = ID_MATRIX;

    /* check if it is inside one of the fibers */
    for( i=0 ; i<fiber_line->ntype ; i++ ){
      for( j=0 ; j<fiber_line->nfib[i] ; j++ )
      {

      }
    }

  }

  return 0;
}

/**********************************************************************/

int micro_struct_init_elem_type(
    micro_struct_t *micro_struct,
    int dim,
    int nelm,
    int (*get_centroid)( int e, int dim, double *elem_centroid ),
    int *elem_type )
{

  /*
     fills the elem_type vector according to the micro structure
     > get_centroid( int e, int dim, double *centroid ) is a function
     that returns the centroid of an element
   */

  int     e;
  double *elem_centroid = malloc(dim*sizeof(double));

  for( e = 0 ; e < nelm ; e++ ){

    get_centroid( e, dim, elem_centroid );
    micro_struct_get_elem_id( dim, micro_struct, elem_centroid,	&elem_type[e] );

  }

  return 0;
}

/****************************************************************************************************/
