#include "micro.h"

void init_variables(params_t *params, message_t *message){

  params->homog_method = HOMOG_METHOD_NULL;
  params->non_linear_max_its = 2;
  params->non_linear_min_norm_tol = 1.0e-4;

  flags.coupled = false;
  flags.linear_materials = false;
  flags.allocated = false;
  flags.c_linear_calculated = false;
  
  comm_init_message(message);

  return;
}
