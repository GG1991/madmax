#ifndef MYIO_H
#define MYIO_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdbool.h>

#ifdef WITH_MPI

#include <mpi.h>

#endif

typedef struct command_line_t_
{
  int argc;
  char **argv;

  int n_str_found;
  char **string_array;

  char *string;

  bool found_flag;

}command_line_t;

command_line_t command_line;


int myio_init_command_line(int argc, char **argv, command_line_t *command_line);

int myio_printf( void *COMM, const char format[], ... );

int myio_get_string_command_line(int argc, char **argv, const char *option_name, char **string, bool *flag_found);

int myio_get_string_array_command_line(int argc, char **argv, const char *option_name, int n_str_expect, char ***strings, bool *flag_found, int *n_str_found);

int myio_search_option_in_command_line(int argc, char **argv, const char *option_name, bool *flag_found);


#endif
