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

#define PRINT_ALWAYS 0
#define PRINTF1(message){if(flags.coupled == false || PRINT_ALWAYS) myio_printf(&MICRO_COMM, message);}
#define PRINTF2(message, arg_1){if(flags.coupled == false || PRINT_ALWAYS) myio_printf(&MICRO_COMM, message,arg_1);}

#define MAX_NVOIGT   6

#define NBUF         256

char        filename[NBUF];

int         rank_mic;          //  rank on macro comm
int         nproc_mic;         //  # of micro processes (MICRO_COMM)
int         dim;               //  problem dimensions
int         nvoi;              //  number of voigt components (3 if dim=2, 6 if dim=3)

char        *myname;

int         micro_type;
FILE        *fm_info;

bool        flag_struct_mesh;
int         nx, ny, nz, nn;       // number of nodes 
int         nl;                   // number of local nodes
int         nelm;
int         nex , ney , nez;      // total number of elements per direction 
int         nyl;                  // local number of nodes in y direction
int         ny_inf;               // inferior numeration in y direction
double      lx, ly, lz;           // dimain lenght
double      hx, hy, hz;           // element sizes 
double      vol_tot;              // total volume
double      vol_loc;              // local volume
double      vol_elem;             // elemental volume
int         npe;                  // Nodes per element
int         ngp;                  // Number of gauss points per element
int         ngho;                 // Number of ghosts nodes
double      **struct_sh;          // Shape functions
double      ***struct_dsh;        // Derivative shapes functions
double      *struct_wp;           // Weights
double      ***struct_bmat;       // B matrix (Bu = epsilon)
double      *elem_disp;           // elemental displacements
int         *loc_elem_index;      // local elemental index vector for assembly and reading
int         *glo_elem_index;      // global elemental index vector for assembly and reading
int         istart, iend;         // starting and ending index of matrices
int         nstart, nend;         // starting and ending index of nodes
int         *dir_ix_loc;          // dirichlet indeces (local numeration)
int         *dir_ix_glo;          // dirichlet indeces (global numeration)
int         ndir_ix;              // number of dirichlet indeces
double      *coor_dir;            // coordinates of dirichlet nodes
double      *strain_gp;           // strain at Gauss point
double      *stress_gp;           // stress at Gauss point
int         *elem_type;           // type in each element (depends on the micro_struc)
double      *elem_strain;         // strain at each element
double      *elem_stress;         // stress at each element
double      *elem_energy;         // energy at each element
double      *center_coor;         // coordinates of the center

double       vi;                  // volumetric fraction of inclusions
double       vm;                  // volumetric fraction of matrix

int          homo_type;
int          macro_gp;
double       c_homo_lineal[36];
double       rho;

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

}params_t;

typedef struct{

  bool coupled;
  bool allocated;
  bool linear_materials;
  bool c_linear_calculated;

}flags_t;

extern params_t params;
extern flags_t flags;

#define PRINT_PETSC        0
#define PRINT_VTU          1

int         flag_print;

double      center_domain[3];

int voigt2mat(double voigt[6], double matrix[3][3]);

int micro_print_info( void );
int micro_pvtu( char *name );

int micro_check_material_and_elem_type(list_t *material_list, int *elem_type, int nelm);

int get_local_elem_index(int e, int *loc_index);
int get_global_elem_index(int e, int *glo_elem_index);
int assembly_b(void);
int assembly_A(void);
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
void finalize(void);

int alloc_memory(void);

#endif
