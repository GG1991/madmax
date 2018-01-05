#include "macro.h"

int set_comm_line_flags(void)
{
  bool found;
  int ierr;

  ierr = myio_comm_line_search_option(&command_line, "-coupl", &found);
  if (ierr != 0)
    return 1;
  else if (found == true)
    flags.coupled = true;

  ierr = myio_comm_line_search_option(&command_line, "-print_matrices", &found);
  if (ierr != 0)
    return 1;
  else if (found == true)
    flags.print_matrices = true;

  ierr = myio_comm_line_search_option(&command_line, "-print_vectors", &found);
  if (ierr != 0)
    return 1;
  else if (found == true)
    flags.print_vectors = true;

  ierr = myio_comm_line_search_option(&command_line, "-print_pvtu", &found);
  if (ierr != 0)
    return 1;
  else if (found == true)
    flags.print_pvtu = true;

  return 0;
}
