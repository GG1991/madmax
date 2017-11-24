/* parallel implementation of some io functions */

#include "myio.h"

/* 
   we defined printf_p (printf parallel) 
   only rank 0 print on stdout, if MPI is not 
   defined it behaves as printf and can be called
   as 

   printf_p( NULL, format[], ... );

 */


int printf_p( void* COMM, const char format[], ... )
{
  va_list arg;
  int     done;

#ifdef WITH_MPI

  int     rank;
  done = MPI_Comm_rank( *(MPI_Comm*)COMM , &rank );

  if( !rank )
  {
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

