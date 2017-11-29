#include "material.h"

int material_get_stress( material_t *mat_p, int dim , double * strain_gp, double * stress_gp )
{

  if( mat_p->type_id == MAT_ELASTIC ){

    /* is a linear material stress = C * strain */

    double  young   = ((type_0*)mat_p->type)->young;
    double  poisson = ((type_0*)mat_p->type)->poisson;
    int     i , j;

    if( dim == 2 ){

      /* Plane strain ( long fibers case ) */
      int     nvoi = 3;
      double  c[3][3];
      c[0][0]=1.0-poisson; c[0][1]=poisson    ; c[0][2]=0.0            ;
      c[1][0]=poisson    ; c[1][1]=1.0-poisson; c[1][2]=0.0            ;
      c[2][0]=0.0        ; c[2][1]=0.0        ; c[2][2]=(1-2*poisson)/2;

      for( i = 0; i < nvoi ; i++ ){
	for( j = 0 ; j < nvoi ; j++ )
	  c[i][j] *= young/((1+poisson)*(1-2*poisson));
      }
      for( i = 0; i < nvoi ; i++ ){
	stress_gp[i] = 0.0;
	for( j = 0 ; j < nvoi ; j++ )
	  stress_gp[i] += c[i][j] * strain_gp[j];
      }

    }
  }
  return 0;
}

/****************************************************************************************************/

int material_get_c_tang( material_t *mat_p, int dim , double * strain_gp, double * c_tan )
{

  if( mat_p->type_id == MAT_ELASTIC ){

    /* is a linear material stress = C * strain */
    double  young   = ((type_0*)mat_p->type)->young;
    double  poisson = ((type_0*)mat_p->type)->poisson;
    int     i , j;
    int     nvoi = 3;

    if(dim==2){

      /* Plane strain ( long fibers case ) */
      c_tan[0*nvoi+0]=1.0-poisson; c_tan[0*nvoi+1]=poisson    ; c_tan[0*nvoi+2]=0.0            ;
      c_tan[1*nvoi+0]=poisson    ; c_tan[1*nvoi+1]=1.0-poisson; c_tan[1*nvoi+2]=0.0            ;
      c_tan[2*nvoi+0]=0.0        ; c_tan[2*nvoi+1]=0.0        ; c_tan[2*nvoi+2]=(1-2*poisson)/2;

      for( i = 0; i < nvoi ; i++ ){
	for( j = 0 ; j < nvoi ; j++ )
	  c_tan[i*nvoi + j] *= young/((1+poisson)*(1-2*poisson));
      }

    }
  }
  return 0;
}

/****************************************************************************************************/

int material_get_rho( material_t *mat_p, int dim , double * rho )
{

  if( mat_p->type_id == MAT_ELASTIC )
    *rho   = ((type_0*)mat_p->type)->rho;

  return 0;
}

/****************************************************************************************************/

int material_fill_list_from_command_line(int argc, const char **argv, list_t *material_list)
{

  list_init(material_list, sizeof(material_t), NULL);

  char **string_array;
  int    found, n_str_found;
  found = myio_get_string_array_command_line(argc, argv, "-material", MAX_NUM_OF_MATERIALS, &string_array, &n_str_found);

  if(found || !n_str_found)
    return 1;

  int i;
  for( i = 0 ; i < n_str_found ; i++ )
  {
    char       *data;
    material_t  mat;
    data = strtok(string_array[i]," \n");
    mat.name = strdup( data );
    data = strtok(NULL," \n");
    if(!strcmp(data,"MAT_ELASTIC"))
    {
      double E, v;
      mat.type_id = MAT_ELASTIC;
      mat.type    = malloc(sizeof(type_0));
      data = strtok( NULL , " \n" );
      ((type_0*)mat.type)->rho         = atof(data);
      data = strtok( NULL , " \n" );
      E = ((type_0*)mat.type)->young   = atof(data);
      data = strtok( NULL , " \n" );
      v = ((type_0*)mat.type)->poisson = atof(data);
      ((type_0*)mat.type)->lambda      = (E*v)/((1+v)*(1-2*v));
      ((type_0*)mat.type)->mu          = E/(2*(1+v));
    }
    else if (!strcmp(data,"MAT_MICRO")){
      mat.type_id = MAT_MICRO;
    }
    else
    {
      myio_printf(NULL,"type %s not known.\n",data);
      return 1;
    }

    list_insertlast(material_list, &mat);
  }

  return 0;
}

/****************************************************************************************************/
