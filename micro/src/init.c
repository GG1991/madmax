#include "micro.h"

void init_variables(params_t *params){

  params->homog_method = HOMOG_METHOD_NULL;
  params->non_linear_max_its = 2;
  params->non_linear_min_norm_tol = 1.0e-4;

  params->flag_coupling = false;

  return;
}
