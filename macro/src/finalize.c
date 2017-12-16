#ifndef FINALIZE_H
#define FINALIZE_H

#include "macro.h"

int finalize(void){

  int ierr;

  if(flags.coupled == true){

    message.action = ACTION_MICRO_END;
    ierr = comm_macro_send(&message);
    if(ierr != 0){
      myio_printf(PETSC_COMM_WORLD, "macro: problem sending MIC_END to micro\n");
      return 1;
    }
  }

  free(loc_elem_index);
  free(glo_elem_index);
  free(elem_disp);
  free(stress_gp);
  free(strain_gp);
  free(elem_strain);
  free(elem_stress);
  free(elem_energy);
  free(elem_type);

  for(int i = 0 ; i < nvoi ; i++){
    for(int j = 0 ; j < npe_max*dim ; j++)
      free(bmat[i][j]);
    free(bmat[i]);
  }
  free(bmat);

  for(int i = 0 ; i < npe_max ; i++){
    for(int j = 0 ; j < dim ; j++)
      free(dsh[i][j]);
    free(dsh[i]);
  }
  free(dsh);

  list_clear(&material_list);
  list_clear(&physical_list);
  list_clear(&function_list);

  MatDestroy(&A);
  VecDestroy(&x);
  VecDestroy(&b);

  comm_finalize_message();

#ifdef SLEPC
  SlepcFinalize();
#elif  PETSC
  PetscFinalize();
#endif

  ierr = MPI_Finalize();

  return ierr;
}

#endif
