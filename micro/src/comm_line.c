#include "micro.h"

int comm_line_set_flags(void)
{
  bool found;
  int ierr;

  ierr = myio_comm_line_search_option(&command_line, "-coupl", &found);
  if (ierr != 0) return 1;
  else if (found == true) flags.coupled = true;

  ierr = myio_comm_line_search_option(&command_line, "-print_matrices", &found);
  if (ierr != 0) return 1;
  else if (found == true) flags.print_matrices = true;

  ierr = myio_comm_line_search_option(&command_line, "-print_vectors", &found);
  if (ierr != 0) return 1;
  else if (found == true) flags.print_vectors = true;

  ierr = myio_comm_line_search_option(&command_line, "-print_pvtu", &found);
  if (ierr != 0) return 1;
  else if (found == true && params.homog_method != HOMOG_METHOD_TAYLOR_PARALLEL && params.homog_method != HOMOG_METHOD_TAYLOR_SERIAL)
    flags.print_pvtu = true;

  myio_comm_line_search_option(&command_line, "-homo_us", &found);
  if (found == true) params.homog_method = HOMOG_METHOD_UNIF_STRAINS;

  myio_comm_line_search_option(&command_line, "-homo_periodic", &found);
  if (found == true) params.homog_method = HOMOG_METHOD_PERIODIC;

  myio_comm_line_search_option(&command_line, "-homo_tp", &found);
  if (found == true) params.homog_method = HOMOG_METHOD_TAYLOR_PARALLEL;

  myio_comm_line_search_option(&command_line, "-homo_ts", &found);
  if (found == true) params.homog_method = HOMOG_METHOD_TAYLOR_SERIAL;

  myio_comm_line_get_int(&command_line, "-nl_max_its", &params.non_linear_max_its, &found);

  myio_comm_line_get_int(&command_line, "-nl_min_norm_tol", &params.non_linear_max_its, &found);

  return 0;
}
