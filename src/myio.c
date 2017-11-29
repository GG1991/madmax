/* parallel implementation of some io functions */

#include "myio.h"

/************************************************************/

/* 
   we defined myio_printf (printf for MPI case) 
   only rank 0 print on stdout, if MPI is not 
   defined it behaves as printf and can be called
   as 

   myio_printf( NULL, format[], ... );

 */


int myio_printf( void* COMM, const char format[], ... )
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

/************************************************************/

/* Get string array from the command line, they should be after
 * "option_name" and each one should be enclose with \"
 * symbol * and separated by commas (,).
 *
 * int    argc : number of command line arguments
 * char **argv : command line arguments
 * const char *option_name : string to search previous to the string array
 *
 */

int myio_get_string_array_command_line(int argc, const char **argv, const char *option_name, int n_str_expect, char ***strings, int *n_str_found)
{

  char *str_token;

  if(!argc || !n_str_expect){
    strings = NULL;
    return 0;
  }

  if(!argv || !option_name)
    return 1;

  *strings = malloc(n_str_expect*sizeof(char**));

  char  *argv_elem_dup;
  int    i=0;
  while(i<argc){
    argv_elem_dup  = strdup(argv[i]);
    str_token = strtok(argv_elem_dup," \"\n");
    if(!strcmp(str_token, option_name)){
      free(argv_elem_dup);
      break;
    }
    free(argv_elem_dup);
    i++;
  }

  int n_str_found_aux = 0;
  if(i<(argc-1)){
    i = i + 1;
    argv_elem_dup  = strdup(argv[i]);
    str_token = strtok(argv_elem_dup,",\"\n");
    n_str_found_aux = 0;
    while(str_token){
      (*strings)[n_str_found_aux++] = strdup(str_token);
      str_token = strtok(NULL,",\"\n");
    }
    free(argv_elem_dup);
  }
  else
    return 1;

  *n_str_found = n_str_found_aux;

  return 0;
}

/************************************************************/

/* Duplicate argv but change type from (char*) to (const char*) */

int myio_duplicate_argv_char_to_const_char(int argc, char **argv, const char ***argv_dup)
{

  if(!argc){
    *argv_dup = NULL;
    return 0;
  }

  if(!argv)
    return 1;

  *argv_dup = malloc(argc*sizeof(char**));
  int i=0;
  for( i=0 ; i<argc ; i++ )
    (*argv_dup)[i] = strdup(argv[i]);

  return 0;
}
