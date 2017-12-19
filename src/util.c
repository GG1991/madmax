#include "util.h"


#ifdef PETSC
int print_petsc_ksp_info(MPI_Comm COMM, KSP ksp){

  int kspits, reason;
  double kspnorm;
  char *reason_s;

  KSPGetIterationNumber( ksp, &kspits );
  KSPGetConvergedReason( ksp, &reason );
  KSPGetResidualNorm   ( ksp, &kspnorm);
  switch(reason){

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


int strbin2dec(char *str){

  int dec = 0;
  for(int i = strlen(str)-1 ; i >= 0 ; i-- ){
    if(str[i]=='0' || str[i]=='1'){
      if(str[i] == '1')
	dec += (int) pow(2,strlen(str)-1-i);
    }else
      return -1;
  }

  return dec;
}
