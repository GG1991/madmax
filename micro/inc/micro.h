#ifndef MICRO_H
#define MICRO_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#include "petscksp.h"
#include "petscsys.h"
#include "slepceps.h"
#include "list.h"
#include "material.h"
#include "comm.h"
#include "trace.h"
#include "myio.h"
#include "vtk.h"
#include "micro_struct.h"
#include "solvers.h"
#include "mesh.h"
#include "fem.h"
#include "ell.h"

enum {MULTIS_NULL, MULTIS_MIXS, MULTIS_MIXP, MULTIS_FE2};
enum {BC_NULL, BC_USTRAIN, BC_USTRESS, BC_PERIODIC};
enum {SOL_PETSC, SOL_ELL};

#define HOMOGENIZE_DELTA_STRAIN      0.005

#define GREEN  "\x1B[32m"
#define RED    "\x1B[31m"
#define BLUE   "\x1B[34m"
#define NORMAL "\x1B[0m"

#define PRINT_ALWAYS 0
#define MIC_PRINTF_0(message){if(flags.coupled == false || PRINT_ALWAYS) printf(message);}
#define MIC_PRINTF_1(message, arg_1){if(flags.coupled == false || PRINT_ALWAYS) printf(message, arg_1);}

#define PRINTF1(message){if(flags.coupled == false || PRINT_ALWAYS) myio_printf(MICRO_COMM, message);}
#define PRINTF2(message, arg_1){if(flags.coupled == false || PRINT_ALWAYS) myio_printf(MICRO_COMM, message,arg_1);}
#define PRINT_ARRAY(name_str, array, length){ PRINTF2("%s ",name_str); for(int i = 0 ; i < length ; i++) PRINTF2("%lf ", array[i]); PRINTF1("\n");}

#define MAX_NVOIGT   6

#define NBUF         256

int rank_mic;
int nproc_mic;
int dim;
int nvoi;

int ngp;
double **struct_sh;
double ***struct_dsh;
double *struct_wp;
double ***struct_bmat;
double *elem_disp;
int *elem_index;
int *elem_nods;
double *strain_gp;
double *stress_gp;
int *elem_type;
double *elem_strain;
double *elem_stress;
double *elem_energy;

double vi;
double vm;

int homo_type;
int macro_gp;

ell_matrix jac_ell;
double *x_ell, *res_ell, *dx_ell;

Mat A, J;
Vec x, b;
Vec dx;
KSP ksp;

typedef struct{

  int multis_method;
  int fe2_bc;
  int nl_max_its;
  int solver;

  double nl_min_norm;
  double c_tangent_linear[MAX_NVOIGT*MAX_NVOIGT];
  double c_tangent[MAX_NVOIGT*MAX_NVOIGT];
  double rho;

}params_t;

typedef struct{

  bool coupled;
  bool allocated;
  bool linear_materials;
  bool c_linear_calculated;
  bool print_pvtu;
  bool print_vectors;
  bool print_matrices;

}flags_t;

extern params_t params;
extern flags_t flags;

double center_domain[3];

int voigt2mat(double voigt[6], double matrix[3][3]);

int micro_print_info( void );
int micro_pvtu( char *name );

int micro_check_material_and_elem_type(list_t *material_list, int *elem_type, int nelm);

int get_local_elem_index(int e, int *loc_index);
int get_global_elem_index(int e, int *glo_elem_index);
int get_stress(int e , int gp, double *strain_gp, double *stress_gp);
int get_strain(int e , int gp, double *strain_gp);
int get_c_tan(const char * name, int e, int gp, double *strain_gp, double *c_tan);
int get_elem_centroid( int e, int dim, double *centroid);
int strain_x_coord( double * strain , double * coord , double *u);
int get_node_local_coor(int n, double *coord);
int get_node_ghost_coor(int n, double *coord);
int get_local_elem_node(int e, int *n_loc);
int local_to_global_index(int local);
int get_averages(double * strain_ave, double *stress_ave);
int get_elem_type(int e, int *type);
int get_elem_properties(void);

int assembly_res(double *norm, double *strain_mac);
int assembly_res_petsc(double *norm, double *strain_mac);
int assembly_jac_petsc(void);
int assembly_jac(void);
int assembly_jac_ell(void);
int assembly_res_ell(double *norm, double *strain_mac);

int init_variables_1(void);
int init_variables_2(void);
int finalize(void);

int alloc_memory(void);

int comm_line_set_flags(void);

int homog_get_c_tangent(double *strain_mac, double **c_tangent);
int homog_get_rho(double *rho);
int homog_calculate_c_tangent(double *strain_mac, double *c_tangent);
int homog_calculate_c_tangent_around_zero(double *c_tangent);
int homog_get_strain_stress(double *strain_mac, double *strain_ave, double *stress_ave);
int homog_get_strain_stress_non_linear(double *strain_mac, double *strain_ave, double *stress_ave);
int homog_taylor(double *strain_mac, double *strain_ave, double *stress_ave);
int homog_fe2(double *strain_mac, double *strain_ave, double *stress_ave);
int negative_res(void);
int add_x_dx(void);
int solve(void);

#endif
