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
#include <string.h>

#define CALC_MODE_NORMAL 1
#define CALC_MODE_EIGEN  2
#define CALC_MODE_TEST   3

#define NBUF 256

int         macro_mode;

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

char        filename[NBUF];      // string for different purposes
char       *mesh_n;              // mesh path;
int         mesh_f;              // mesh format

Mat         A;
Mat         M;
Vec         x, dx, b;

typedef struct{

  int calc_mode;

  int non_linear_max_its;
  double non_linear_min_norm_tol;

  int num_eigen_vals;
  double *eigen_vals;

  double final_time;
  double delta_time;
  double time;

  int time_step;

  double energy_stored;

}params_t;

extern params_t params;

list_t boundary_list;
list_t physical_list;

#endif
