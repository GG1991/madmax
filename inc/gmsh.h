/* 
   gmsh functions for reading thing grom file
 */

#ifndef GMSH_H
#define GMSH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "list.h"

#define NBUF_GMSH 256

typedef struct _physical_t{
    
    int    dim;
    int    id;
    char   *name;

}physical_t;

int gmsh_get_node_index( const char * mesh_n, const char * phy_name, int nmynods, int *mynods, int dim, int * n, int ** ix );
int gmsh_get_physical_list( char *mesh_n, list_t *physical_list );
int gmsh_which_id( const char * mesh_n, const char * name );
int gmsh_is_surf( int code, int dim );
int gmsh_npe( int code );
int gmsh_funcmp_int_a(void *a, void *b);
int gmsh_funcmp_int_b(const void *a, const void *b);
int gmsh_read_coord_parall( char *mesh_n, int dim, int nmynods, int *mynods, int nghost , int *ghost, double *coord );

#endif
