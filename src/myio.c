#include "myio.h"


int myio_printf(void* COMM, const char format[], ...){

  va_list arg;
  int     done;

#ifdef WITH_MPI
  int     rank;
  done = MPI_Comm_rank( *(MPI_Comm*)COMM , &rank );
  if(!rank){
    va_start(arg, format);
    done = vfprintf( stdout, format, arg );
    va_end(arg);
  }
#else
  va_start(arg, format);
  done = vfprintf( stdout, format, arg );
  va_end(arg);
#endif

  return done;
}


int myio_comm_line_init(int argc, char **argv, command_line_t *command_line){

  if(command_line == NULL) return 1;

  command_line->argc = argc;
  command_line->argv = malloc(argc*sizeof(char*));
  for(int i = 0 ; i < argc ; i++)
    command_line->argv[i] = strdup(argv[i]);

  return 0;
}


int myio_comm_line_search_option(command_line_t *command_line, const char *option_name, bool *found){

  *found = false;

  if(command_line->argv == 0 || option_name == 0) return 1;
  if(command_line->argc == 0) return 0;

  int i = 0; SEARCH_ARGV_INDEX(i, option_name)
  if(i < command_line->argc) *found = true;

  return 0;
}


int myio_comm_line_get_int(command_line_t *command_line, const char *option_name, int *value, bool *found){

  *found = false;

  if(command_line->argv == NULL || option_name == NULL) return 1;
  if(command_line->argc == 0) return 0;

  int i = 0; SEARCH_ARGV_INDEX(i, option_name)
  if(i >= command_line->argc - 1) return 1;

  *value = atoi(command_line->argv[i+1]);
  *found = true;

  return 0;
}


int myio_comm_line_get_int_array(command_line_t *command_line, const char *option_name, int *values, int nval_expect, int *nval_found, bool *found){

  *found = false;

  if(command_line->argv == NULL || option_name == NULL) return 1;
  if(command_line->argc == 0) return 0;

  int i = 0; SEARCH_ARGV_INDEX(i, option_name)
  if(i >= command_line->argc - 1) return 1;

  char *argv_dup = strdup(command_line->argv[i+1]);
  char *str_token = strtok(argv_dup, ",\n");

  int n = 0;
  while(str_token){
    values[n] = atoi(str_token);
    str_token = strtok(NULL, ",\n");
    n++;
  }
  free(argv_dup);

  *nval_found = n;
  *found = true;

  return 0;
}


int myio_comm_line_get_double(command_line_t *command_line, const char *option_name, double *value, bool *found){

  *found = false;

  if(command_line->argv == NULL || option_name == NULL || value == NULL) return 1;
  if(command_line->argc == 0) return 0;

  int i = 0; SEARCH_ARGV_INDEX(i, option_name)
  if(i >= command_line->argc - 1) return 1;

  *value = atof(command_line->argv[i+1]);
  *found = true;

  return 0;
}


int myio_comm_line_get_string(command_line_t *command_line, const char *option_name, char *string, bool *found){

  *found = false;

  if(command_line->argv == NULL || option_name == NULL || string == NULL) return 1;
  if(command_line->argc == 0) return 0;

  int i = 0; SEARCH_ARGV_INDEX(i, option_name)
  if(i == command_line->argc - 1) return 1;

  strcpy(string, command_line->argv[i+1]);

  *found = true;

  return 0;
}


int myio_comm_line_get_string_array(command_line_t *command_line, const char *option_name, char string_arr[][STRING_LENGTH], int nval_expect, int *nval_found, bool *found){

  *found = false;

  if(command_line->argv == NULL || option_name == NULL || string_arr == NULL) return 1;

  if(command_line->argc == 0) return 0;

  int i = 0; SEARCH_ARGV_INDEX(i, option_name)
  if(i == command_line->argc - 1) return 1;

  char *str_token, *argv_dup;
  argv_dup = strdup(command_line->argv[i+1]);
  str_token = strtok(argv_dup, ",\n");

  int n = 0;
  while(str_token){
    strcpy(string_arr[n], strdup(str_token));
    str_token = strtok(NULL, ",\n");
    n++;
  }
  free(argv_dup);

  *nval_found = n;
  *found = true;
  
  return 0;
}
