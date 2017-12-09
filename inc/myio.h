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

#define SEARCH_ARGV_INDEX(i, option_name) { \
  while(i < command_line->argc){ \
    if(strcmp(command_line->argv[i], option_name) == 0) break; \
    i++; \
  } \
}

typedef struct{

  int argc;
  char **argv;

}command_line_t;

command_line_t command_line;

int myio_printf(void *COMM, const char format[], ...);

int myio_comm_line_init(int argc, char **argv, command_line_t *command_line);

int myio_comm_line_search_option(command_line_t *command_line, const char *option_name, bool *found);

int myio_comm_line_get_int(command_line_t *command_line, const char *option_name, int *value, bool *found);
int myio_comm_line_get_int_array(command_line_t *command_line, const char *option_name, int *values, int nval_expect, int *nval_found, bool *found);

int myio_comm_line_get_double(command_line_t *command_line, const char *option_name, double *value, bool *found);

int myio_comm_line_get_string(command_line_t *command_line, const char *option_name, char *string, bool *found);
int myio_comm_line_get_string_array(command_line_t *command_line, const char *option_name, char **string_arr, int nval_expect, int *nval_found, bool *found);


#endif
