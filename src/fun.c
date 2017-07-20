/*
 * Function utilities routines
 * 
 */
#include "fun.h"

int f1d_init(double *x, double *y, int n, int inter, f1d_t *f1d){

    int i;
    if(n==0 || !x || !y || !f1d)
        return 1;
    f1d->n=n;
    f1d->inter=inter;
    f1d->fnum=-1;
    f1d->x=(double*)calloc(n,sizeof(double));
    f1d->y=(double*)calloc(n,sizeof(double));
    if(!f1d->x || !f1d->y)
        return 1;
    for(i=0;i<n;i++){
        f1d->x[i]=x[i];
        f1d->y[i]=y[i];
    }
    return 0;
}

int f1d_eval(double x, f1d_t *f1d, double *y){

    int i;

    if(!f1d) 
        return 1;
    if(f1d->n<1)
        return 1;
    if(f1d->n==1){
        *y=f1d->y[0];
        return 0;
    }
    if(x<f1d->x[0]){
        *y=f1d->y[0];
        return 0;
    }
    i=1;
    while(i<f1d->n){
        if(f1d->x[i-1]<=x && x<f1d->x[i])
            break;
        i++;
    }
    if(i==f1d->n){
        *y=f1d->y[i-1];
        return 0;
    }
    *y=(f1d->y[i]-f1d->y[i-1])*(x - f1d->x[i-1])/(f1d->x[i]-f1d->x[i-1]) + f1d->y[i-1];
    return 0;
}

int cmp_f1d(void *a, void *b){
    if ( ((f1d_t *)a)->fnum > ((f1d_t *)b)->fnum ){
	return 1;
    }else if( ((f1d_t*)a)->fnum == ((f1d_t*)b)->fnum ){
	return 0;
    }else{
	return -1;
    }
}

f1d_t * GetFunctionPointer(list_t *function_list, int fnumToSearch)
{
  /*
   * Returns a pointer to <f1d_t> element in function list
   */

   if(!function_list) return NULL;
   if(!function_list->sizelist) return NULL;
   node_list_t *pFun;
   pFun = function_list->head;
   while(pFun){
     if( ((f1d_t*)pFun->data)->fnum == fnumToSearch ) break;
     pFun = pFun->next;
   }
   return (f1d_t*)pFun->data;

}
