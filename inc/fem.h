/* Finite element procedures
 * 
 * 
 */

#ifndef FEM2_H
#define FEM2_H

#include <stdio.h>
#include <stdlib.h>
#include "math.h"

int fem_inigau(void);
int fem_caljac3(double coor[8][3],double ***ds, int npe,int gp,double jac[3][3]);
int fem_invjac3(double jac[3][3],double ijac[3][3],double *det);
int fem_calode(int npe, int dim, double ****oder);
int fem_calwei(int npe, int dim, double **wp);
int fem_calder3(double ijac[3][3],int nsh,int gp,double ***oder,double der[3][3]);
int fem_calare(double **pts, int npe, int dim, double *area);
int fem_vecmod(double *vec, int n, double *mod);
int fem_dotdsh(int i, int j, double **derivs, int dim, double *p);
int fem_vcross(double *v1, double *v2, double *vr);

// Segment 2 nodes
double **xp_segm_2;
double *wp_segm_2;
double **sh_segm_2;
double ***ds_segm_2;

// Triangle 3 nodes
double **xp_tria_3;
double *wp_tria_3;
double **sh_tria_3;
double ***ds_tria_3;

// Quadrangle 4 nodes
double **xp_quad_4;
double *wp_quad_4;
double **sh_quad_4;
double ***ds_quad_4;

// Tetrahedron 4 nodes
double **xp_tetra_4;
double *wp_tetra_4;
double **sh_tetra_4;
double ***ds_tetra_4;

// Prism 6 nodes
double **xp_prism_6;
double *wp_prism_6;
double **sh_prism_6;
double ***ds_prism_6;
          
// Hexahedron 8 nodes
double **xp_hexa_8;
double *wp_hexa_8;
double **sh_hexa_8;
double ***ds_hexa_8;

#endif
