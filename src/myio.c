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

  command_line->str = NULL;
  command_line->str_arr = NULL;
  command_line->int_arr = NULL;

  return 0;
}


int myio_comm_line_search_option(command_line_t *command_line, const char *option_name){

  command_line->found = false;

  if(! command_line->argv || ! option_name)
    return 1;

  if(! command_line->argc)
    return 0;

  int i = 0;

  while(i < command_line->argc){
    if(! strcmp(command_line->argv[i], option_name)){
      command_line->found = true;
      break;
    }
    i++;
  }

  return 0;
}


int myio_comm_line_get_int(command_line_t *command_line, const char *option_name){

  command_line->found = false;

  if(! command_line->argv || ! option_name)
    return 1;

  if(! command_line->argc)
    return 0;

  int i = 0;

  while(i < command_line->argc){
    if(! strcmp(command_line->argv[i], option_name)) break;
    i++;
  }

  if(i < command_line->argc - 1){
    i = i + 1;
    command_line->int_val = atoi(command_line->argv[i]);
    command_line->found = true;
  }
  else{
    command_line->str = NULL;
  }

  return 0;
}


int myio_comm_line_get_int_array(command_line_t *command_line, int n_int_expec, const char *option_name){

  command_line->found = false;

  if(command_line->argv == NULL || option_name == NULL)
    return 1;

  if(command_line->argc == 0)
    return 0;

  if(command_line->int_arr != NULL) free(command_line->int_arr);
  command_line->int_arr = malloc(n_int_expec*sizeof(int));

  int i = 0;

  while(i < command_line->argc){
    if(strcmp(command_line->argv[i], option_name) == 0) break;
    i++;
  }

  if(i < command_line->argc - 1){

    i = i + 1;

    char *argv_dup = strdup(command_line->argv[i]);
    char *str_token = strtok(argv_dup, ",\n");
    command_line->n_str_found = 0;

    while(str_token){
      command_line->int_arr[command_line->n_str_found] = atoi(str_token);
      str_token = strtok(NULL, ",\n");
      command_line->n_str_found++;
    }
    free(argv_dup);

    command_line->found = true;
  }
  else
    command_line->int_arr = NULL;

  return 0;
}


int myio_comm_line_get_double(command_line_t *command_line, const char *option_name){

  command_line->found = false;

  if(! command_line->argv || ! option_name)
    return 1;

  if(! command_line->argc)
    return 0;

  int i = 0;

  while(i < command_line->argc){
    if(! strcmp(command_line->argv[i], option_name)) break;
    i++;
  }

  if(i < command_line->argc - 1){
    i = i + 1;
    command_line->double_val = atof(command_line->argv[i]);
    command_line->found = true;
  }
  else{
    command_line->str = NULL;
  }

  return 0;
}


int myio_comm_line_get_string_array(command_line_t *command_line, int n_str_expec, const char *option_name){

  command_line->found = false;

  if(! command_line->argv || ! option_name)
    return 1;

  if(! command_line->argc)
    return 0;

  myio_free_string_array(command_line->str_arr, command_line->n_str_found);

  command_line->n_str_expec = n_str_expec;
  command_line->str_arr = malloc(command_line->n_str_expec*sizeof(char**));

  int i = 0;

  while(i < command_line->argc){
    if(! strcmp(command_line->argv[i], option_name))
      break;
    i++ ;
  }

  command_line->n_str_found = 0;

  if(i < command_line->argc - 1){
    i = i + 1;
    char *str_token, *argv_dup;
    argv_dup = strdup(command_line->argv[i]);
    str_token = strtok(argv_dup, ",\n");
    while(str_token){
      command_line->str_arr[command_line->n_str_found] = strdup(str_token);
      str_token = strtok(NULL, ",\n");
      command_line->n_str_found++;
    }
    command_line->found = true;
  }
  else{
    command_line->str_arr = NULL;
  }

  return 0;
}


int myio_comm_line_get_string(command_line_t *command_line, const char *option_name){

  command_line->found = false;

  if(! command_line->argv || ! option_name)
    return 1;

  if(! command_line->argc)
    return 0;

  myio_free_string(command_line->str);

  int i = 0;

  while(i < command_line->argc){
    if(! strcmp(command_line->argv[i], option_name)) break;
    i++;
  }

  if(i < command_line->argc - 1){
    i = i + 1;
    char *str_token, *argv_dup;
    argv_dup = strdup(command_line->argv[i]);
    str_token = strtok(argv_dup, "\n");
    command_line->str = strdup(str_token);
    command_line->found = true;
  }
  else{
    command_line->str = NULL;
  }

  return 0;
}


int myio_free_string_array(char **str_arr, int n_str_found){

  if(str_arr == NULL) return 0;

  for(int i = 0 ; i < n_str_found ; i++)
    if(str_arr[i]) free(str_arr[i]);

  free(str_arr);

  return 0;
}


int myio_free_string(char *str){

  if(str)
    free(str);

  return 0;
}
