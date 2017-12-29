#include "micro_struct.h"

int micro_struct_init(int dim, char *string, micro_struct_t *micro_struct) {

  char *stra = strdup(string );
  char *data = strtok(stra, " \n" );
  double *size = malloc(dim * sizeof(double) );

  if (strcmp(data, "fiber_cilin") == 0) {

    /* Reads "size"[dim] "nx_fib" "ny_fib" "radio" "desv"[2]*/

    fiber_cilin_t *fiber_cilin = malloc(sizeof(fiber_cilin_t));

    fiber_cilin->desv = malloc(dim*sizeof(double));

    for (int d = 0 ; d < dim ; d++) {
      data = strtok( NULL, " \n" );
      if (!data) return 2;
      size[d] = atof( data );
    }

    data = strtok( NULL, " \n" );
    if (!data) return 2;
    fiber_cilin->nx_fib = atoi( data );

    data = strtok( NULL, " \n" );
    if (!data) return 2;
    fiber_cilin->ny_fib = atoi( data );

    data = strtok( NULL, " \n" );
    if (!data) return 2;
    fiber_cilin->radio = atof( data );

    for (int d = 0 ; d < dim ; d++) {
      data = strtok( NULL, " \n" );
      if (!data) return 2;
      fiber_cilin->desv[d] = atof( data );
    }

    micro_struct->type = FIBER_CILIN;
    micro_struct->data = fiber_cilin;

  }else if (strcmp(data, "fiber_line") == 0) {

    /* Reads "size"[dim] "ntype" "nfib[ntype]" "tetha[ntype]"  "seps[ntype]" "width[ntype]" "desv[ntype]"*/

    int ntype;

    fiber_line_t *fiber_line = malloc(sizeof(fiber_line_t));

    for (int d = 0 ; d < dim ; d++) {
      data = strtok( NULL, " \n" );
      if (!data) return 2;
      size[d] = atof( data );
    }

    data = strtok( NULL, " \n" );
    if (!data) return 2;
    fiber_line->ntype = ntype = atoi( data );

    fiber_line->theta = malloc(ntype*sizeof(double));
    fiber_line->sep   = malloc(ntype*sizeof(double));
    fiber_line->width = malloc(ntype*sizeof(double));
    fiber_line->desv  = malloc(ntype*sizeof(double));
    fiber_line->nfib  = malloc(ntype*sizeof(double));

    for (int d = 0 ; d < dim ; d++) {
      data = strtok( NULL, " \n" );
      if (!data) return 2;
      fiber_line->nfib[d] = atoi( data );
    }

    for (int d = 0 ; d < dim ; d++) {
      data = strtok( NULL, " \n" );
      if (!data) return 2;
      fiber_line->theta[d] = atof( data );
    }

    for (int d = 0 ; d < dim ; d++) {
      data = strtok( NULL, " \n" );
      if (!data) return 2;
      fiber_line->sep[d] = atof( data );
    }

    for (int d = 0 ; d < dim ; d++) {
      data = strtok( NULL, " \n" );
      if (!data) return 2;
      fiber_line->width[d] = atof( data );
    }

    for (int d = 0 ; d < dim ; d++) {
      data = strtok( NULL, " \n" );
      if (!data) return 2;
      fiber_line->desv[d] = atof( data );
    }

    micro_struct->type = FIBER_LINE;
    micro_struct->data = fiber_line;

  }else
    return 1;

  micro_struct->size = size;

  return 0;
}


int micro_struct_get_elem_id(int dim, micro_struct_t *micro_struct, double *elem_centroid, int *elem_id) {

  double lx = micro_struct->size[0];
  double ly = micro_struct->size[1];

  if (micro_struct->type == FIBER_CILIN) {

      fiber_cilin_t * fiber_cilin = (fiber_cilin_t *)micro_struct->data;

      double deviation[2];
      double center[2];
      center[0] = lx / 2;
      center[1] = ly / 2;

      *elem_id = ID_MATRIX;

      for (int i = 0 ; i < fiber_cilin->nx_fib ; i++) {
	for (int j = 0 ; j < fiber_cilin->ny_fib ; j++) {
	  deviation[0] = fiber_cilin->desv[0] - lx/2 + (lx/fiber_cilin->nx_fib)/2 + i * (lx/fiber_cilin->nx_fib);
	  deviation[1] = fiber_cilin->desv[1] - ly/2 + (ly/fiber_cilin->ny_fib)/2 + j * (ly/fiber_cilin->ny_fib);
	  double l = 0.0;
	  for (int d = 0 ; d < 2 ; d++) {
	    l = l + pow( elem_centroid[d] - (center[d] + deviation[d]), 2 );
	  }
	  l = sqrt(l);
	  if (l <= fiber_cilin->radio) {
	    *elem_id = ID_FIBER ;
	    return 0;
	  }

	}
      }

  }else if (micro_struct->type == FIBER_LINE) {

    fiber_line_t * fiber_line = (fiber_line_t *)micro_struct->data;

    *elem_id = ID_MATRIX;

    double  p_line[2];
    double  n_line[2];
    int     side_1, side_2;
    double *point = elem_centroid;


    for (int i = 0 ; i < fiber_line->ntype ; i++) {
      for (int j = 0 ; j < fiber_line->nfib[i] ; j++) {

	/* normal vector of the line */
	n_line[0] = sin(fiber_line->theta[i]);
	n_line[1] = cos(fiber_line->theta[i]);

	/* line 1 (below) p_line = ( 0 , ly/2 + fiber_line->sep[i] + desv[i] - width[i]/2 ) */
	p_line[0] = lx/2;
	p_line[1] = ly/2 + fiber_line->sep[i]*(j-fiber_line->nfib[i]/2) \
		    + fiber_line->desv[i] + fiber_line->width[i]/2;
	side_1 = geom_2d_line_side(n_line, p_line, point);

	/* line 1 (upper) p_line = ( 0 , ly/2 + fiber_line->sep[i] + desv[i] + width[i]/2 ) */
	p_line[0] = lx/2;
	p_line[1] = ly/2 + fiber_line->sep[i]*(j-fiber_line->nfib[i]/2) \
		    + fiber_line->desv[i] - fiber_line->width[i]/2;
	side_2 = geom_2d_line_side(n_line, p_line, point);

	if (side_1 == 0 || side_2 == 0 || (side_1*side_2) == -1) {
	    *elem_id = ID_FIBER;
	    return 0;
	}

      }
    }

  }

  return 0;
}


int micro_struct_init_elem_type(micro_struct_t *micro_struct, int dim, int nelm, int (*get_centroid)(int e, int dim, double *elem_centroid), int *elem_type)
{
  double *elem_centroid = malloc(dim*sizeof(double));

  for (int e = 0 ; e < nelm ; e++) {
    get_centroid(e, dim, elem_centroid);
    micro_struct_get_elem_id(dim, micro_struct, elem_centroid, &elem_type[e]);
  }
  return 0;
}
