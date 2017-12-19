#ifndef FEM_H
#define FEM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int fem_init(void);
int fem_init_struct(double ***sh, double ****dsh, double **wp, double *h, int dim);
int fem_invjac(int dim, double ** jac, double ** ijac, double *det);
int fem_calc_jac(int dim, int npe, int gp, double * coor, double *** dsh, double ** jac);

int fem_trans_dsh(int dim, int nsh, int gp, double **ijac, double ***dsh_master, double ***dsh);
int fem_get_dsh_master(int npe, int dim , double ****dsh);
int fem_get_sh(int npe, int dim, double ***sh);
int fem_get_dsh( int npe, int dim, int gp, double **dsh, double *detj);
int fem_get_wp(int npe, int dim, double **wp);

double **xp_segm_2;
double *wp_segm_2;
double **sh_segm_2;
double ***ds_segm_2;

double **xp_tria_3;
double *wp_tria_3;
double **sh_tria_3;
double ***ds_tria_3;

double **xp_quad_4;
double *wp_quad_4;
double **sh_quad_4;
double ***ds_quad_4;

double **xp_tetra_4;
double *wp_tetra_4;
double **sh_tetra_4;
double ***ds_tetra_4;

double **xp_prism_6;
double *wp_prism_6;
double **sh_prism_6;
double ***ds_prism_6;

double **xp_hexa_8;
double *wp_hexa_8;
double **sh_hexa_8;
double ***ds_hexa_8;


#endif
