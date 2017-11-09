/*
   functions for reading thing from gmsh file

   Author: Guido Giuntoli
   Date  : 09 - 11 - 17
*/
#include "gmsh.h"

int gmsh_get_node_index( const char * mesh_n, const char * bou_name, int nmynods, int *mynods, int dim, int * n, int ** ix )
{

  FILE         *fm = fopen(mesh_n,"r"); if( fm == NULL ) return 1;
  int          i, d, npe;
  int          flag = 0; 
  int          ntag;
  int          id_s, id; 
  int          ns, *nf; 
  char         buf[NBUF_GMSH], *data;   
  list_t       nod_list;
  node_list_t  *pn;

  list_init( &nod_list, sizeof(int), &gmsh_funcmp_int_a );

  /* searchs in the physical entities the "id" */
  id_s = gmsh_which_id( mesh_n, bou_name ); if( id_s < 0 ) return 1;

  while( fgets( buf, NBUF_GMSH, fm ) != NULL ){

    data = strtok(buf," \n");
    if( flag == 0 ){
      if( strcmp(data,"$Elements") == 0 )
      {
	fgets( buf, NBUF_GMSH, fm );
	flag  = 1;
      }

    }
    else{

      data = strtok( NULL, " \n" );

      if( gmsh_is_surf( atoi(data), dim ) ){

	npe  = gmsh_npe(atoi(data));
	data = strtok(NULL," \n");
	ntag = atoi(data);
	data = strtok(NULL," \n");
	id   = atoi(data);

	if( id == id_s ){

	  // salteamos los tags y nos vamos derecho para los nodos
	  d = 1;
	  while( d < ntag ){
	    data = strtok(NULL," \n");
	    d++;
	  }
	  i = 0;
	  while( i < npe ){
	    data = strtok(NULL," \n");
	    ns   = atoi(data);
	    nf   = bsearch( &ns, mynods, nmynods, sizeof(int), &gmsh_funcmp_int_b );
	    if( nf )
	      list_insert_se( &nod_list, nf );
	    i++;
	  }
	}

      }
      else{
	*n  = nod_list.sizelist;
	*ix = malloc( *n * sizeof(int));
	pn  = nod_list.head;
	i = 0;
	while( pn ){
	  (*ix)[i++] = *( int * ) pn->data;
	  pn = pn->next;
	}
	list_clear( &nod_list );
	fclose(fm);
	return 0;   
      }
    } 
  }

  return 1;
}

/****************************************************************************************************/

int gmsh_which_id( const char * mesh_n, const char * name )
{

  /* returns the "id" of a physical entity with name "name" */

  FILE         *fm = fopen(mesh_n,"r"); if( fm == NULL ) return 1;
  int          id;
  int          i = 0, total;
  int          flag = 0; 
  char         buf[NBUF_GMSH], *data;   

  while( fgets( buf, NBUF_GMSH, fm ) != NULL ){

    data = strtok(buf," \n");

    if( flag == 0 ){
      if( !strcmp(data,"$PhysicalNames"))
      {
	fgets( buf, NBUF_GMSH, fm );
	data  = strtok(buf," \n");
	total = atoi(data);
	flag = 1;
      }
    }
    else{

      data = strtok( NULL, " \n" );
      id   = atoi(data);
      data = strtok( NULL, " \"\n" );
      if( strcmp( name , data ) == 0 )
	return id;
      i++;

    }
    if( i == total ) break;

  }

  return -1;
}

/****************************************************************************************************/

int  gmsh_funcmp_int_a(void *a, void *b){
     if( *(int*)a > *(int*)b ) return  1;
     if( *(int*)a < *(int*)b ) return -1;
     return 0;
}

/****************************************************************************************************/

int  gmsh_funcmp_int_b(const void *a, const void *b){
     if( *(int*)a > *(int*)b ) return  1;
     if( *(int*)a < *(int*)b ) return -1;
     return 0;
}

/****************************************************************************************************/

int gmsh_is_surf( int code , int dim )
{
  // returns 1 if the code corresponds to a surface element, 0 othewhise
  if(dim == 2)
    return (code == 1 || code == 15) ? 1 : 0;
  else if(dim == 3)
    return (code == 1 || code == 2 || code == 3 || code == 15) ? 1 : 0;
  return 1;
}

/****************************************************************************************************/

int gmsh_npe(int code)
{
    switch(code){
	case 1:  
	    return 2; 
	case 2:    
	    return 3;
	case 3:    
	    return 4;
	case 4:    
	    return 4;
	case 5:    
	    return 8;
	case 6:   
	    return 6;
	case 15:   
	    return 1;
	default:    
	    return -1;
    }
}

/****************************************************************************************************/
