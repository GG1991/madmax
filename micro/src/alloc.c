#include "micro.h"

int alloc_memory(void){

  int ierr = 0;

  if(solver.type == SOLVER_PETSC){

#ifdef PETSC
    PETSC_COMM_WORLD = MICRO_COMM;
    PetscInitialize(&command_line.argc, &command_line.argv, (char*)0, NULL);
#else
    return 1;
#endif
  }
  
  loc_elem_index = malloc(dim*npe*sizeof(int));
  glo_elem_index = malloc(dim*npe*sizeof(int));
  elem_disp = malloc(dim*npe*sizeof(double));
  stress_gp = malloc(nvoi*sizeof(double));
  strain_gp = malloc(nvoi*sizeof(double));
  elem_strain = malloc(nelm*nvoi*sizeof(double));
  elem_stress = malloc(nelm*nvoi*sizeof(double));
  elem_energy = malloc(nelm*sizeof(double));
  elem_type = malloc(nelm*sizeof(int));

  struct_bmat = malloc(nvoi*sizeof(double**));
  for(int i = 0 ; i < nvoi ; i++){
    struct_bmat[i] = malloc(npe*dim*sizeof(double*));
    for(int j = 0 ; j < npe*dim ; j++)
      struct_bmat[i][j] = malloc(ngp*sizeof(double));
  }

  return ierr;
}
