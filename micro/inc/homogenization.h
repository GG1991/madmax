/*
   Homogenization Variables
*/

#define HOMO_TAYLOR           1
#define HOMO_LINEAR           2
#define HOMO_LINEAR_HEXA      3
#define LD_LAGRAN_SEQ         4
#define LD_LAGRAN_PAR         5

/*
   Linear deformation using Lagrange multipliers
*/

typedef struct _homog_ld_lagran_t{

  int nnods_bc;
  int *index;
  double *ub_val;

}homog_ld_lagran_t;

typedef struct _homog_t{

  int type;
  int dim;
  void *st;

}homog_t;

homog_t homo;

int mic_init_homo();
