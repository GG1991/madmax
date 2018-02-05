#include "micro.h"

int comm_line_set_flags(void)
{
  bool found;
  int ierr;

  ierr = myio_comm_line_search_option(&command_line, "-coupl", &found);
  if (ierr != 0) return 1;
  else if (found == true) flags.coupled = true;

  ierr = myio_comm_line_search_option(&command_line, "-print_matrices", &found);
  if (ierr != 0) return 2;
  else if (found == true) flags.print_matrices = true;

  ierr = myio_comm_line_search_option(&command_line, "-print_vectors", &found);
  if (ierr != 0) return 3;
  else if (found == true) flags.print_vectors = true;

  ierr = myio_comm_line_search_option(&command_line, "-print_pvtu", &found);
  if (ierr != 0) return 4;
  else if (found == true && params.multis_method == MULTIS_FE2)
    flags.print_pvtu = true;

  myio_comm_line_search_option(&command_line, "-fe2", &found);
  if (found == true) {
    params.multis_method = MULTIS_FE2;

    myio_comm_line_search_option(&command_line, "-bc_ustrain", &found);
    if (found == true) params.fe2_bc = BC_USTRAIN;

    myio_comm_line_search_option(&command_line, "-bc_periodic", &found);
    if (found == true) params.fe2_bc = BC_PERIODIC;

    myio_comm_line_search_option(&command_line, "-bc_ustress", &found);
    if (found == true) params.fe2_bc = BC_USTRESS;

    if (params.fe2_bc == BC_NULL) {
      MIC_PRINTF_0("no bc was given for the fe2 multis method (-bc_ustrain -bc_ustress -bc_ustress)\n");
      return 5;
    }
  }

  myio_comm_line_search_option(&command_line, "-mixp", &found);
  if (found == true) params.multis_method = MULTIS_MIXP;

  myio_comm_line_search_option(&command_line, "-mixs", &found);
  if (found == true) params.multis_method = MULTIS_MIXS;

  myio_comm_line_get_int(&command_line, "-nl_max_its", &params.non_linear_max_its, &found);

  myio_comm_line_get_int(&command_line, "-nl_min_norm_tol", &params.non_linear_max_its, &found);

  return 0;
}
