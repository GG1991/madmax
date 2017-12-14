#ifndef MACRO_H
#define MACRO_H

#include "sputnik.h"
#include "comm.h"
#include "util.h"
#include "gmsh.h"
#include "function.h"
#include "vtk.h"
#include "myio.h"
#include "mesh.h"
#include "solvers.h"

#define CALC_MODE_NULL   0
#define CALC_MODE_NORMAL 1
#define CALC_MODE_EIGEN  2
#define CALC_MODE_TEST   3

#define NBUF 256

double    **struct_sh;           // Shape functions
double   ***struct_dsh;          // Derivative shapes functions
double     *struct_wp;           // Weights
double   ***bmat;                // B matrix (Bu = epsilon)
int         npe_max;             // maximum number of nodes per element
int         ngp_max;             // maximum number of gauss points
double     *elem_disp;           // elemental displacements
double     *elem_coor;           // coordinates of the element's verteces
int        *loc_elem_index;      // local elemental index vector for assembly and reading
int        *glo_elem_index;      // global elemental index vector for assembly and reading
double     *strain_gp;           // strain at Gauss point
double     *stress_gp;           // stress at Gauss point
int        *elem_type;           // type in each element
double     *elem_strain;         // strain at each element
double     *elem_stress;         // stress at each element
double     *elem_energy;         // energy at each element
double     *res_elem;            // elemental residue
double     *k_elem;              // elemental steffiness matrix
double     *m_elem;              // elemental mass matrix
double     *c;                   // constitutive tensor
double   ***dsh;                 // shape functions derivatives at gauss points
double     *detj;                // jacobian determinants for each gauss point
double    **jac;                 // elemental jacobian
double    **jac_inv;             // inverse of elemental jacobian

int         nvoi;                // voigt number

int         mymicro_rank_worker;
int         flag_neg_detj;       // negative jacobian flag

int         rank_mac;
int         nproc_mac;

char        filename[NBUF];
char        mesh_n[128];
int         mesh_f;

Mat         A;
Mat         M;
Vec         x, dx, b;

typedef struct{

  int calc_mode;

  int non_linear_max_its;
  double non_linear_min_norm_tol;

  int num_eigen_vals;
  double *eigen_vals;

  double tf;
  double dt;
  double t;
  int ts;

  double energy_stored;
  int non_linear_its;
  double residual_norm;

}params_t;

typedef struct{

  bool coupled;
  bool allocated;
  bool print_pvtu;
  bool print_vectors;
  bool print_matrices;

}flags_t;

extern params_t params;
extern flags_t flags;

list_t boundary_list;
list_t physical_list;

int read_bc(void);
int assembly_A_petsc(void);
int assembly_b_petsc(void);
int assembly_AM_petsc(void);
int update_bound(double t);
int get_global_elem_index(int e, int *glo_elem_index);
int get_local_elem_index(int e, int *loc_elem_index);
int get_elem_properties(void);
int get_strain(int e , int gp, int *loc_elem_index, double ***dsh, double ***bmat, double *strain_gp);
int get_stress(int e , int gp, double *strain_gp , double *stress_gp);
int get_c_tan(const char * name, int e, int gp, double *strain_gp, double *c_tan);
int get_rho(const char * name, int e, double *rho);
int get_sh(int dim, int npe, double ***sh);
int get_dsh(int e, int * loc_elem_index, double ***dsh, double *detj);
int get_wp(int dim, int npe, double **wp);
int get_bmat(int e, double ***dsh, double ***bmat);
int get_mat_name(int id, char * name_s);
int macro_pvtu(char *name);
int update_boundary(double t, list_t *function_list, list_t *boundary_list);
int read_coord(char *mesh_n, int nmynods, int *mynods, int nghost, int *ghost, double **coord);

void init_variables(void);
int finalize(void);

int alloc_memory(void);

#endif
