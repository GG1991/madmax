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

int gmsh_get_node_index( const char * mesh_n, const char * phy_name, int nmynods, int *mynods, int dim, int * n, int ** ix );
int gmsh_which_id( const char * mesh_n, const char * name );
int gmsh_is_surf( int code, int dim );
int gmsh_npe( int code );
int gmsh_funcmp_int_a(void *a, void *b);
int gmsh_funcmp_int_b(const void *a, const void *b);

#endif
