#ifndef MICRO_H
#define MICRO_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include "homogenize.h"
#include "micro_struct.h"
#include "solvers.h"
#include "mesh.h"

#define GREEN "\x1B[32m"
#define RED "\x1B[31m"
#define NORMAL "\x1B[0m"

#define PRINT_ALWAYS 0
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

Mat A, J;
Vec x, b;
Vec dx;
KSP ksp;

typedef struct{

  int homog_method;
  int non_linear_max_its;
  double non_linear_min_norm_tol;
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
int assembly_b_petsc(void);
int assembly_A_petsc(void);
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
int init_shapes(double ***sh, double ****dsh, double **wp);

void homogenize_check_linear_material(void);

void init_variables(void);
int finalize(void);

int alloc_memory(void);

#endif
