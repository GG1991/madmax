#include "material.h"


int material_get_stress( material_t *mat_p, int dim , double * strain_gp, double * stress_gp ) {

  if (mat_p->type_id == MAT_ELASTIC) {

    double  young   = ((type_0*)mat_p->type)->young;
    double  poisson = ((type_0*)mat_p->type)->poisson;

    if (dim == 2) {

      /* plane strain ( long fibers case ) */
      int     nvoi = 3;
      double  c[3][3];
      c[0][0]=1.0-poisson; c[0][1]=poisson    ; c[0][2]=0.0            ;
      c[1][0]=poisson    ; c[1][1]=1.0-poisson; c[1][2]=0.0            ;
      c[2][0]=0.0        ; c[2][1]=0.0        ; c[2][2]=(1-2*poisson)/2;

      for (int i = 0; i < nvoi ; i++) {
	for (int j = 0 ; j < nvoi ; j++)
	  c[i][j] *= young/((1+poisson)*(1-2*poisson));
      }
      for (int i = 0; i < nvoi ; i++) {
	stress_gp[i] = 0.0;
	for (int j = 0 ; j < nvoi ; j++)
	  stress_gp[i] += c[i][j] * strain_gp[j];
      }

    }
  }
  return 0;
}


int material_get_c_tang(material_t *mat_p, int dim, double *strain_gp, double *c_tan) {

  if (mat_p->type_id == MAT_ELASTIC) {

    double  young   = ((type_0*)mat_p->type)->young;
    double  poisson = ((type_0*)mat_p->type)->poisson;
    int     nvoi = 3;

    if (dim == 2) {

      /* plane strain ( long fibers case ) */
      c_tan[0*nvoi+0]=1.0-poisson; c_tan[0*nvoi+1]=poisson    ; c_tan[0*nvoi+2]=0.0            ;
      c_tan[1*nvoi+0]=poisson    ; c_tan[1*nvoi+1]=1.0-poisson; c_tan[1*nvoi+2]=0.0            ;
      c_tan[2*nvoi+0]=0.0        ; c_tan[2*nvoi+1]=0.0        ; c_tan[2*nvoi+2]=(1-2*poisson)/2;

      for (int i = 0; i < nvoi ; i++) {
	for (int j = 0 ; j < nvoi ; j++)
	  c_tan[i*nvoi + j] *= young/((1+poisson)*(1-2*poisson));
      }

    }
  }
  return 0;
}


int material_get_rho(material_t *mat_p, int dim, double *rho) {

  if (mat_p->type_id == MAT_ELASTIC)
    *rho = ((type_0*)mat_p->type)->rho;

  return 0;
}


bool material_are_all_linear(list_t *material_list) {

  node_list_t *node_ptr = material_list->head;
  while (node_ptr != NULL) {
    material_t *mat = node_ptr->data;
    if (mat->is_linear == false) break;
    node_ptr = node_ptr->next;
  }
  return (node_ptr == NULL) ? true : false;

}


int material_fill_list_from_command_line(command_line_t *command_line, list_t *material_list) {

  bool found;
  int num_string_found;
  char string_arr[MAX_NUM_OF_MATERIALS][128];

  myio_comm_line_get_string_array(command_line, "-material", string_arr, MAX_NUM_OF_MATERIALS, &num_string_found, &found);

  if (found == false || num_string_found == 0)
    return 1;

  list_init(material_list, sizeof(material_t), NULL);


  for (int i = 0 ; i < num_string_found ; i++) {

    char *data = strtok(string_arr[i], " \n");

    material_t  mat;
    mat.name = strdup( data );

    data = strtok(NULL, " \n");
    if (strcmp(data, "MAT_ELASTIC") == 0) {

      mat.type_id = MAT_ELASTIC;
      mat.type    = malloc(sizeof(type_0));
      mat.is_linear = true;

      data = strtok( NULL , " \n");
      ((type_0*)mat.type)->rho         = atof(data);

      data = strtok( NULL , " \n");
      double E = ((type_0*)mat.type)->young   = atof(data);

      data = strtok( NULL , " \n");
      double v = ((type_0*)mat.type)->poisson = atof(data);

      ((type_0*)mat.type)->lambda      = (E*v)/((1+v)*(1-2*v));
      ((type_0*)mat.type)->mu          = E/(2*(1+v));

    }
    else if (strcmp(data,"MAT_MICRO") == 0) {
      mat.type_id = MAT_MICRO;
    }
    else
      return 1;

    list_insertlast(material_list, &mat);
  }

  return 0;
}


bool material_check_in_list(list_t *material_list, char *name_to_check) {

  node_list_t *pm = material_list->head;

  while (pm != NULL) {
    material_t *mat_p = (material_t *)pm->data;
    if (strcmp(mat_p->name, name_to_check) == 0) break;
    pm = pm->next;
  }

  return (pm == NULL) ? false : true;
}
