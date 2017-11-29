#ifndef MYIO_H
#define MYIO_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

int myio_printf( void *COMM, const char format[], ... );

int myio_get_string_array_command_line(int argc, const char **argv, const char *option_name, int n_str_expect, char ***strings, int *n_str_found);

int myio_duplicate_argv_char_to_const_char(int argc, char **argv, const char ***argv_dup);

#endif
