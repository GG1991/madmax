/*
 *  Function declarations and prototypes
 *
 */
 

#ifndef FUN_H
#define FUN_H

#include "stdlib.h"
#include "list.h"

enum {INTER1,INTER2}; 

typedef struct _f1d_t{
    
    int n;
    int inter;
    int fnum;    
    
    double *x;
    double *y;
    
}f1d_t;

int f1d_init( double * x , double * y , int n , int inter , f1d_t * f1d );
int f1d_eval( double   x , f1d_t * f1d , double * y );
int cmp_f1d ( void *a, void *b);
int get_f1d ( int fn , list_t * function_list , f1d_t ** f1d );

#endif
