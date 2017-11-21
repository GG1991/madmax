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

int gmsh_get_physical_list( char *mesh_n, list_t *physical_list )
{

  FILE  *fm;
  int   i = 0, total;
  int   flag = 0;
  char  buf[NBUF_GMSH], *data;

  physical_t physical;

  fm = fopen(mesh_n,"r"); if( fm == NULL ) return 1;

  while(fgets(buf,NBUF_GMSH,fm)!=NULL)
  {
    data=strtok(buf," \n");
    if( strcmp( data,"$PhysicalNames") == 0 ){
      fgets( buf, NBUF_GMSH, fm );
      data  = strtok(buf," \n");
      total = atoi(data);
      flag  = 1;
    }
    else if( flag == 1 ){

      physical.dim = atoi(data);

      data=strtok(NULL," \n");
      physical.id = atoi(data);

      data=strtok(NULL," \"\n");
      physical.name = strdup( data );

      list_insertlast( physical_list , &physical );
      i ++;

      if( i == total ) break;
    }
  }
  return 0;
}

/****************************************************************************************************/

int gmsh_read_coord_parall( char *mesh_n, int dim, int nmynods, int *mynods, int nghost , int *ghost, double *coord )
{

  int   i, c, d;
  int   ln, offset;  // line counter and offset for moving faster in the file
  int   ntotnod;

  char  buf[NBUF_GMSH];
  char  *data;

  FILE * fm = fopen(mesh_n,"r");
  if( !fm ){
    printf( "file %s not found\n" , mesh_n );
    return 1;
  }

  offset   = 0;
  while(fgets(buf,NBUF_GMSH,fm)!=NULL){
    offset += strlen(buf);
    data=strtok(buf," \n");
    if(strcmp(data,"$Nodes")==0){
      fgets(buf,NBUF_GMSH,fm);
      offset += strlen(buf);
      data  = strtok(buf," \n");
      ntotnod = atoi(data);

      /* start with local nodes */
      i = c = 0;
      while( c < nmynods ){
	while( i< mynods[c] ){
	  fgets(buf,NBUF_GMSH,fm);
	  i++;
	}
	data=strtok(buf," \n");
	for( d=0;d<dim;d++){
	  data=strtok(NULL," \n");
	  coord[c*dim + d] = atof(data);
	}
	c++;
      }
      if( c > ntotnod ){
	printf("gmsh_read_coord : more nodes (%d) in %s than calculated (%d)\n", ntotnod, mesh_n, c);
	return 1;
      }
      break;
    }
    ln ++;
  }

  /* continue with ghost nodes */
  fseek( fm, offset, SEEK_SET);
  i = c = 0;
  while( c < nghost ){
    while( i<ghost[c] ){
      fgets(buf,NBUF_GMSH,fm);
      i++;
    }
    data=strtok(buf," \n");
    for( d=0;d<dim;d++){
      data=strtok(NULL," \n");
      coord[ (nmynods + c)*dim + d] = atof(data);
    }
    c++;
  }
  if( c > ntotnod ){
    printf("gmsh_read_coord : more nodes (%d) in %s than calculated (%d)\n", ntotnod, mesh_n, c);
    return 1;
  }
  fclose(fm);
  return 0;
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
