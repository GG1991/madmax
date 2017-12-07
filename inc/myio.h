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

  char *str;

  int int_val;
  int n_int_found;
  int n_int_expect;
  int *int_arr;

  double double_val;

  bool found;

}command_line_t;

command_line_t command_line;


int myio_printf(void *COMM, const char format[], ...);

int myio_free_string_array(char **string_arr, int n_str);

int myio_free_string(char *str);

int myio_comm_line_init(int argc, char **argv, command_line_t *command_line);

int myio_comm_line_search_option(command_line_t *command_line, const char *option_name);

int myio_comm_line_get_int(command_line_t *command_line, const char *option_name);
int myio_comm_line_get_int_array(command_line_t *command_line, int n_int_expec, const char *option_name);

int myio_comm_line_get_double(command_line_t *command_line, const char *option_name);

int myio_comm_line_get_string(command_line_t *command_line, const char *option_name);
int myio_comm_line_get_string_array(command_line_t *command_line, const char *option_name, char ***string_arr, int max_num_string_expected, int *num_string_found, bool *found);


#endif
