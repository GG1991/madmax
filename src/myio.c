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

int myio_get_string_array_command_line(int argc, char **argv, const char *option_name, int n_str_expect, char ***strings, bool *flag_found, int *n_str_found)
{

  *flag_found = false;

  if(!argc || !n_str_expect){
    strings = NULL;
    return 0;
  }

  if(!argv || !option_name)
    return 1;

  *strings = malloc(n_str_expect*sizeof(char**));

  char *str_token, *argv_dup;
  int i = 0;

  while(i < argc){
    argv_dup = strdup(argv[i]);
    str_token = strtok(argv_dup," \"\n");
    if(!strcmp(str_token, option_name)){
      *flag_found = true;
      break;
    }
    i++ ;
  }

  int n_str_found_aux = 0;

  if(i < (argc-1)){
    i = i + 1;
    argv_dup = strdup(argv[i]);
    str_token = strtok(argv_dup,",\"\n");
    n_str_found_aux = 0;
    while(str_token){
      (*strings)[n_str_found_aux++] = strdup(str_token);
      str_token = strtok(NULL,",\"\n");
    }
  }
  else
    return 1;

  *n_str_found = n_str_found_aux;

  return 0;
}

int myio_search_option_in_command_line(int argc, char **argv, const char *option_name, bool *flag_found)
{

  *flag_found = false;

  if(!argv || !option_name)
    return 1;

  if(!argc)
    return 0;

  char *str_token, *argv_dup;
  int i = 0;

  while(i < argc){
    argv_dup = strdup(argv[i]);
    str_token = strtok(argv_dup," \n");
    if(!strcmp(str_token,option_name)){
      *flag_found = true;
      break;
    }
    i++;
  }

  return 0;
}

int myio_get_string_command_line(int argc, char **argv, const char *option_name, char **string, bool *flag_found)
{

  *flag_found = false;

  if(!argv || !option_name)
    return 1;

  if(!argc)
    return 0;

  char *str_token, *argv_dup;
  int i = 0;

  while(i < argc){
    argv_dup = strdup(argv[i]);
    str_token = strtok(argv_dup," \n");
    if(!strcmp(str_token,option_name)){
      *flag_found = true;
      break;
    }
    i++;
  }

  if(i < (argc-1)){
    i = i + 1;
    *string = strdup(argv[i]);
  }

  return 0;
}
