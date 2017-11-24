#include "micro_struct.h"

/**********************************************************************/

int init_micro_structure( int dim, micro_struct_t * micro_struct, char * format )
{

  char * data = strtok( format, " \n" );
  int    d;

  if( ! strcmp( data, "fiber_cilin" ) )
  {
    /* -micro_struct "fiber_cilin <size[0]> <size[1]> <nx> <ny> <radio> <desv[0]> <desv[1]>" */

    fiber_cilin_t * fiber_cilin;
    micro_struct->type = FIBER_CILIN;
    micro_struct->data = malloc( sizeof(fiber_cilin_t) );
    fiber_cilin = ( fiber_cilin_t * )micro_struct->data;
    fiber_cilin->size = malloc( dim * sizeof(double) );
    fiber_cilin->desv = malloc( dim * sizeof(double) );

    for( d = 0 ; d < dim ; d++ ){
      data = strtok( NULL, " \n" );
      fiber_cilin->size[d] = atof( data );
    }

    data = strtok( NULL, " \n" );
    fiber_cilin->nx = atoi( data );
    data = strtok( NULL, " \n" );
    fiber_cilin->ny = atoi( data );

    data = strtok( NULL, " \n" );
    fiber_cilin->radio = atof( data );

    for( d = 0 ; d < 2 ; d++ ){
      data = strtok( NULL, " \n" );
      fiber_cilin->desv[d] = atof( data );
    }

  }
  else if( ! strcmp( data, "fiber_planar" ) )
  {

   /*
   -micro_struct "fiber_planar <size[0]> <size[1]> <ntype> <angle[0]> <seps[0]> <width[0]> <desv[0]> <numb[0]>
                                                           <angle[1]> <seps[1]> <width[1]> <desv[1]> <numb[1]>
      						     ..."
    */

    micro_struct->type = FIBER_PLANAR;

  }
  else{
    return 1;
  }

  return 0;
}

/**********************************************************************/

int get_elem_id_centroid( int dim, micro_struct_t * micro_struct, double * elem_centroid, int * elem_id )
{

  /* returns elem_id as a function of the centroid coordinate */

  if( micro_struct->type == FIBER_CILIN )
  {

      fiber_cilin_t * fiber_cilin = ( fiber_cilin_t * )micro_struct->data;

      int    i, j, d;
      double deviation[2];
      double center[2];
      double lx = fiber_cilin->size[0];
      double ly = fiber_cilin->size[1];
      center[0] = lx / 2;
      center[1] = ly / 2;
      for( i = 0 ; i < fiber_cilin->nx ; i++ )
      {
	for( j = 0 ; j < fiber_cilin->ny ; j++ )
	{
	  deviation[0] = fiber_cilin->desv[0] - lx/2 + ( lx / fiber_cilin->nx )/2 + i * ( lx / fiber_cilin->nx );
	  deviation[1] = fiber_cilin->desv[1] - ly/2 + ( ly / fiber_cilin->ny )/2 + j * ( ly / fiber_cilin->ny );
	  double l = 0.0;
	  for( d = 0 ; d < 2 ; d++ )
	    l = l + pow( elem_centroid[d] - (center[d] + deviation[d]), 2 );
	  l = sqrt(l);
	  *elem_id = ( l <= fiber_cilin->radio ) ? FIBER_ID : MATRIX_ID;
	}
      }

  }

  return 0;
}

/**********************************************************************/
