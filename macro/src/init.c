#include "macro.h"

void init_variables(params_t *params){

  params->calc_mode = CALC_MODE_NULL;
  params->non_linear_max_its = 1;
  params->non_linear_min_norm_tol = 1.0;
  params->num_eigen_vals = 1;
  params->final_time = 0.0;
  params->delta_time = 0.0;
  params->time = 0.0;
  params->time_step = 0.0;
  params->energy_stored = 1.0;

  return;
}
