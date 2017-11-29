/* Function declarations and prototypes
 *  
 */
 
#ifndef FUN_H
#define FUN_H

#include "stdlib.h"
#include "list.h"
#include "myio.h"

#define MAX_NUM_OF_FUNCTIONS 4

enum {INTER1,INTER2}; 

typedef struct _f1d_t{
    
    int n;
    int inter;
    int fnum;    
    
    double *x;
    double *y;
    
}f1d_t;

int function_init(double *x, double *y, int n, int inter, f1d_t * f1d);
int function_eval(double x, f1d_t *f1d, double * y);
int function_comp(void *a, void *b);
int function_get_from_list(int fn , list_t *function_list, f1d_t **f1d);
int function_fill_list_from_command_line(int argc, const char **argv, list_t *function_list);

#endif
