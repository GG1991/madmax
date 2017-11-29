#ifndef MYIO_H
#define MYIO_H

#include <stdio.h>
#include <stdarg.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

int myio_printf( void *COMM, const char format[], ... );

int myio_get_strings_command_line(
    int          argc,
    const int  **argv,
    const char  *option_name,
    char       **strings,
    int          n_str_expect,
    int          n_str_found
    );

#endif
