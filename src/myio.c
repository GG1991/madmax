#include "myio.h"

int myio_printf(void* COMM, const char format[], ...)
{
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

int myio_init_command_line(int argc, char **argv, command_line_t *command_line){

  if(!command_line){
    return 1;
  }
  else{
    command_line->argc = argc;
    command_line->argv = malloc(argc*sizeof(char*));
    int i;
    for( i=0 ; i<argc ; i++ ){
      command_line->argv[i] = strdup(argv[i]);
    }
  }

  return 0;
}

int myio_search_option_in_command_line(command_line_t *command_line, const char *option_name)
{

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


int myio_get_string_array_command_line(command_line_t *command_line, int n_str_expec, const char *option_name)
{

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
    command_line->n_str_found = 0;
    while(str_token){
      command_line->str_arr[command_line->n_str_found] = strdup(str_token);
      str_token = strtok(NULL, ",\n");
      command_line->n_str_found++;
    }
    command_line->found = true;
  }

  return 0;
}

int myio_get_string_command_line(command_line_t *command_line, const char *option_name)
{

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

  int i;
  if(! str_arr) return 0;

  for( i=0 ; i<n_str_found ; i++ ){
    if(str_arr[i]) free(str_arr[i]);
  }
  free(str_arr);

  return 0;
}


int myio_free_string(char *str){

  if(str)
    free(str);

  return 0;
}
