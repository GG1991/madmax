/* Gmsh declarations
 * 
 * 
 */

#ifndef GMSH_H
#define GMSH_H

typedef struct _physical_t{
    
    int    dim;
    int    GmshID;
    char   *name; 

}physical_t;

list_t physical_list;

#endif
