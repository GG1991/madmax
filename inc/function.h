#ifndef FUN_H
#define FUN_H

#include "stdlib.h"
#include "list.h"
#include "myio.h"

#define MAX_NUM_OF_FUNCTIONS 4

typedef struct _function_t{

    int n;
    int inter;
    int fnum;

    double *x;
    double *y;

}function_t;

int function_init(double *x, double *y, int n, int inter, function_t * f1d);
int function_eval(double x, function_t *f1d, double * y);
int function_comp(void *a, void *b);
int function_get_from_list(int fn , list_t *function_list, function_t **f1d);
int function_fill_list_from_command_line(command_line_t *command_line, list_t *function_list);

list_t function_list;

#endif
