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

typedef struct command_line_t_{

  int argc;
  char **argv;

  char **str_arr;
  int n_str_found;
  int n_str_expec;

  char *str;

  bool found;

}command_line_t;

command_line_t command_line;


int myio_printf(void *COMM, const char format[], ...);

int myio_free_string_array(char **string_arr, int n_str);

int myio_free_string(char *str);

int myio_init_command_line(int argc, char **argv, command_line_t *command_line);

int myio_get_string_command_line(command_line_t *command_line, const char *option_name);

int myio_get_string_array_command_line(command_line_t *command_line, int n_str_expec, const char *option_name);

int myio_search_option_in_command_line(command_line_t *command_line, const char *option_name);


#endif
