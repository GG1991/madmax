#include "util.h"


#ifdef PETSC
int print_petsc_ksp_info(MPI_Comm COMM, KSP ksp) {

  int kspits, reason;
  double kspnorm;
  char *reason_s;

  KSPGetIterationNumber( ksp, &kspits );
  KSPGetConvergedReason( ksp, &reason );
  KSPGetResidualNorm   ( ksp, &kspnorm);
  switch(reason) {

    case KSP_CONVERGED_RTOL:
      reason_s = strdup( "RTOL" );
      break;
    case KSP_CONVERGED_ATOL:
      reason_s = strdup( "ATOL" );
      break;
    default :
      reason_s = strdup( "UNKNOW" );
      break;
  }
  PetscPrintf(COMM,"kspits %D kspnorm %e kspreason %s", kspits, kspnorm, reason_s );

  return 0;
}
#endif


int strbin2dec(char *str) {

  int dec = 0;
  for (int i = strlen(str)-1 ; i >= 0 ; i-- ) {
    if (str[i]=='0' || str[i]=='1') {
      if (str[i] == '1')
	dec += (int) pow(2,strlen(str)-1-i);
    }else
      return -1;
  }

  return dec;
}


int util_is_in_vector(int val, int *vector, int size) {

  int j = 0;
  while (j < size) {
    if (vector[j] == val) break;
    j++;
  }
  return (j == size) ? 0 : 1;

  return -1;
}


int util_clean_and_sort_vector(int *in_vector, int n_in, int **out_vector, int *n_out) {

  int swi, val_o;

  if (n_in == 0) return 0;

  int *aux = malloc(n_in*sizeof(int));
  for (int i = 0 ; i < n_in ; i++) aux[i] = in_vector[i];

  qsort(aux, n_in, sizeof(int), util_cmpfunc);

  val_o = aux[0];
  int c = 1;
  for (int i = 1 ; i < n_in ; i++) {
    swi = 1;
    if (aux[i] == val_o) {
      swi = 0;
    }
    else{
      val_o = aux[i];
      swi = 1;
    }
    if (swi==1) {
      c++;
    }
  }
  (*out_vector) = malloc(c*sizeof(int));

  val_o = aux[0];
  (*out_vector)[0] = aux[0];
  c = 1;
  for (int i = 1 ; i < n_in ; i++) {
    swi = 1;
    if (aux[i] == val_o)
      swi = 0;
    else{
      val_o = aux[i];
      swi = 1;
    }
    if (swi == 1) (*out_vector)[c++] = aux[i];
  }

  free(aux);
  *n_out = c;

  return 0;
}


int util_sort_vector_intersec(int *array_1, int n1, int *array_2, int n2, int **inter_arr, int *n_inter) {

  int index_1 = 0, index_2 = 0, n_inter_aux = 0;

  while (index_1 < n2 && index_2 < n1) {
    if (array_1[index_2] < array_2[index_1])
      index_2 ++;
    else if (array_1[index_2] > array_2[index_1])
      index_1 ++;
    else if (array_1[index_2] == array_2[index_1]) {
      index_2 ++;
      index_1 ++;
      n_inter_aux ++;
    }
  }

  *n_inter = n_inter_aux;
  *inter_arr = malloc(n_inter_aux*sizeof(int));

  index_1 = index_2 = n_inter_aux = 0;
  while (index_1 < n2 && index_2 < n1) {
    if (array_1[index_2] < array_2[index_1])
      index_2 ++;
    else if (array_1[index_2] > array_2[index_1])
      index_1 ++;
    else if (array_1[index_2] == array_2[index_1]) {
      (*inter_arr)[n_inter_aux] = array_2[index_1];
      index_2 ++;
      index_1 ++;
      n_inter_aux ++;
    }
  }

  return 0;
}


int util_cmpfunc(const void * a, const void * b) {

  return (*(int*)a - *(int*)b);
}
