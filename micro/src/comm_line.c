#include "micro.h"

int comm_line_set_flags(void)
{
  bool found;

  myio_comm_line_search_option(&command_line, "-coupl", &found); if (found == true) flags.coupled = true;

  myio_comm_line_search_option(&command_line, "-solver_petsc", &found); if (found == true) params.solver = SOL_PETSC;
  myio_comm_line_search_option(&command_line, "-solver_ell", &found); if (found == true) params.solver = SOL_ELL;

  myio_comm_line_search_option(&command_line, "-fe2", &found);
  if (found == true) {
    params.multis_method = MULTIS_FE2;
    myio_comm_line_search_option(&command_line, "-bc_ustrain", &found); if (found == true) params.fe2_bc = BC_USTRAIN;
    myio_comm_line_search_option(&command_line, "-bc_ustress", &found); if (found == true) params.fe2_bc = BC_USTRESS;
    myio_comm_line_search_option(&command_line, "-bc_per_lm", &found); if (found == true) params.fe2_bc = BC_PER_LM;
    myio_comm_line_search_option(&command_line, "-bc_per_ms", &found); if (found == true) params.fe2_bc = BC_PER_MS;
    if (params.fe2_bc == BC_NULL) {
      MIC_PRINTF_0("no bc was given for the fe2 multis method (-bc_ustrain -bc_ustress -bc_ustress)\n");
      return 5;
    }
  }

  myio_comm_line_search_option(&command_line, "-print_matrices", &found); if (found == true) flags.print_matrices = true;
  myio_comm_line_search_option(&command_line, "-print_vectors", &found); if (found == true) flags.print_vectors = true;
  myio_comm_line_search_option(&command_line, "-print_pvtu", &found); if (found == true) flags.print_pvtu = true;

  myio_comm_line_search_option(&command_line, "-mixp", &found); if (found == true) params.multis_method = MULTIS_MIXP;
  myio_comm_line_search_option(&command_line, "-mixs", &found); if (found == true) params.multis_method = MULTIS_MIXS;

  myio_comm_line_get_int(&command_line, "-nl_max_its", &params.nl_max_its, &found);
  myio_comm_line_get_double(&command_line, "-nl_min_norm", &params.nl_min_norm, &found);
  return 0;
}
