#include "fem.h"


int fem_init_struct(double ***sh, double ****dsh, double **wp, double *h, int dim) {

  int nsh = ( dim == 2 ) ? 4 : 8;
  int ngp = ( dim == 2 ) ? 4 : 8;
  double *xp = malloc( ngp*dim * sizeof(double));

  if (dim == 2) {
    xp[0] = -0.577350269189626;   xp[1]= -0.577350269189626;
    xp[2] = +0.577350269189626;   xp[3]= -0.577350269189626;
    xp[4] = +0.577350269189626;   xp[5]= +0.577350269189626;
    xp[6] = -0.577350269189626;   xp[7]= +0.577350269189626;
  }

  *sh = malloc( nsh * sizeof(double*));
  for (int is = 0 ; is < nsh ; is++) {
    (*sh)[is] = malloc(ngp*sizeof(double));
  }

  *dsh  = malloc( nsh * sizeof(double**));
  for (int is = 0 ; is < nsh ; is++) {
    (*dsh)[is] = malloc(dim * sizeof(double*));
    for (int d = 0 ; d < dim ; d++) {
      (*dsh)[is][d] = malloc(ngp*sizeof(double));
    }
  }

  *wp = malloc(ngp*sizeof(double));

  if (dim == 2) {

    for (int gp = 0 ; gp < ngp ; gp++) {
      (*sh)[0][gp] = (1 - xp[2*gp]) * (1 - xp[2*gp+1])/4;
      (*sh)[1][gp] = (1 + xp[2*gp]) * (1 - xp[2*gp+1])/4;
      (*sh)[2][gp] = (1 + xp[2*gp]) * (1 + xp[2*gp+1])/4;
      (*sh)[3][gp] = (1 - xp[2*gp]) * (1 + xp[2*gp+1])/4;
    }

    double hx = h[0], hy = h[1];
    for (int gp = 0 ; gp < ngp ; gp++) {
      (*dsh)[0][0][gp] = -1 * (1 - xp[2*gp+1]) /4 * 2/hx; // d phi / d x
      (*dsh)[1][0][gp] = +1 * (1 - xp[2*gp+1]) /4 * 2/hx;
      (*dsh)[2][0][gp] = +1 * (1 + xp[2*gp+1]) /4 * 2/hx;
      (*dsh)[3][0][gp] = -1 * (1 + xp[2*gp+1]) /4 * 2/hx;
      (*dsh)[0][1][gp] = -1 * (1 - xp[2*gp+0]) /4 * 2/hy; // d phi / d y
      (*dsh)[1][1][gp] = -1 * (1 + xp[2*gp+0]) /4 * 2/hy;
      (*dsh)[2][1][gp] = +1 * (1 + xp[2*gp+0]) /4 * 2/hy;
      (*dsh)[3][1][gp] = +1 * (1 - xp[2*gp+0]) /4 * 2/hy;
    }

    double vol = hx*hy;
    for (int gp = 0 ; gp < ngp ; gp++)
      (*wp)[gp] = vol / ngp;

  }

  free(xp);

  return 0;
}


int fem_init(void) {

    xp_segm_2 = calloc(2,sizeof(double*));
    wp_segm_2 = calloc(2,sizeof(double));
    sh_segm_2 = calloc(2,sizeof(double*));
    ds_segm_2 = calloc(2,sizeof(double**));

    for (int i = 0 ; i < 2 ; i++)
      ds_segm_2[i] = calloc(1,sizeof(double*));

    for (int i = 0 ; i < 2 ; i++) {
      xp_segm_2[i] = calloc(1,sizeof(double));
      sh_segm_2[i] = calloc(2,sizeof(double));
      for (int j = 0 ; j < 1 ; j++)
	ds_segm_2[i][j] = (double *)calloc(2,sizeof(double));
    }

    xp_segm_2[0][0] = -0.577350269189626; wp_segm_2[0]= +1.0;
    xp_segm_2[1][0] = +0.577350269189626; wp_segm_2[1]= +1.0;

    for (int gp = 0 ; gp < 2 ; gp++) {
      sh_segm_2[0][gp]= +0.5*(1.0-xp_segm_2[gp][0]);
      ds_segm_2[0][0][gp]= -0.5;

      sh_segm_2[1][gp]= +0.5*(1.0+xp_segm_2[gp][0]);
      ds_segm_2[1][0][gp]= +0.5;
    }

    xp_tria_3 = calloc(3,sizeof(double*));
    wp_tria_3 = calloc(3,sizeof(double));
    sh_tria_3 = calloc(3,sizeof(double*));
    ds_tria_3 = calloc(3,sizeof(double**));

    for (int i = 0 ; i < 3 ; i++)
      ds_tria_3[i]=(double **)calloc(2,sizeof(double*));

    for (int i = 0 ; i < 3 ; i++) {
      xp_tria_3[i]=(double *)calloc(2,sizeof(double));
      sh_tria_3[i]=(double *)calloc(3,sizeof(double));
      for (int j = 0 ; j < 2 ; j++)
	ds_tria_3[i][j]=(double *)calloc(3,sizeof(double));
    }

    xp_tria_3[0][0] = +0.166666666666666;
    xp_tria_3[0][1] = +0.166666666666666;
    wp_tria_3[0] = 1/6.0; // detj[gp] = 2*Area => wp[gp] = 1/(2*3)

    xp_tria_3[1][0] = +0.666666666666666;
    xp_tria_3[1][1] = +0.166666666666666;
    wp_tria_3[1] = 1/6.0;

    xp_tria_3[2][0] = +0.166666666666666;
    xp_tria_3[2][1] = +0.666666666666666;
    wp_tria_3[2] = 1/6.0;


    for (int gp = 0 ; gp < 3 ; gp++) {

      sh_tria_3[0][gp] = +1.0 - xp_tria_3[gp][0] - xp_tria_3[gp][1];
      ds_tria_3[0][0][gp] = -1.0;
      ds_tria_3[0][1][gp] = -1.0;

      sh_tria_3[1][gp] = +xp_tria_3[gp][0];
      ds_tria_3[1][0][gp] = +1.0;
      ds_tria_3[1][1][gp] = +0.0;

      sh_tria_3[2][gp] = +xp_tria_3[gp][1];
      ds_tria_3[2][0][gp] = +0.0;
      ds_tria_3[2][1][gp] = +1.0;
    }


    xp_quad_4 = calloc(4,sizeof(double*));
    wp_quad_4 = calloc(4,sizeof(double));
    sh_quad_4 = calloc(4,sizeof(double*));
    ds_quad_4 = calloc(4,sizeof(double**));

    for (int i = 0 ; i < 4 ; i++)
      ds_quad_4[i] = calloc(2,sizeof(double*));

    for (int i = 0 ; i < 4 ; i++) {
      xp_quad_4[i] = calloc(2,sizeof(double));
      sh_quad_4[i] = calloc(4,sizeof(double));
      for (int j = 0 ; j < 2 ; j++)
	ds_quad_4[i][j] = calloc(4,sizeof(double));
    }

    xp_quad_4[0][0] = -0.577350269189626;  xp_quad_4[0][1] = -0.577350269189626;  wp_quad_4[0] = 1.0;
    xp_quad_4[1][0] = +0.577350269189626;  xp_quad_4[1][1] = -0.577350269189626;  wp_quad_4[1] = 1.0;
    xp_quad_4[2][0] = +0.577350269189626;  xp_quad_4[2][1] = +0.577350269189626;  wp_quad_4[2] = 1.0;
    xp_quad_4[3][0] = -0.577350269189626;  xp_quad_4[3][1] = +0.577350269189626;  wp_quad_4[3] = 1.0; // detj[gp] = A/4 => wp[gp] = 1

    for (int gp = 0 ; gp < 4 ; gp++) {

      sh_quad_4[0][gp] = (1.0-xp_quad_4[gp][0])*(1.0-xp_quad_4[gp][1])*0.25;
      ds_quad_4[0][0][gp] = -1.0*(1.0-xp_quad_4[gp][1])*0.25;
      ds_quad_4[0][1][gp] = -1.0*(1.0-xp_quad_4[gp][0])*0.25;

      sh_quad_4[1][gp] = (1.0+xp_quad_4[gp][0])*(1.0-xp_quad_4[gp][1])*0.25;
      ds_quad_4[1][0][gp] = +1.0*(1.0-xp_quad_4[gp][1])*0.25;
      ds_quad_4[1][1][gp] = -1.0*(1.0+xp_quad_4[gp][0])*0.25;

      sh_quad_4[2][gp] = (1.0+xp_quad_4[gp][0])*(1.0+xp_quad_4[gp][1])*0.25;
      ds_quad_4[2][0][gp] = +1.0*(1.0+xp_quad_4[gp][1])*0.25;
      ds_quad_4[2][1][gp] = +1.0*(1.0+xp_quad_4[gp][0])*0.25;

      sh_quad_4[3][gp] = (1.0-xp_quad_4[gp][0])*(1.0+xp_quad_4[gp][1])*0.25;
      ds_quad_4[3][0][gp] = -1.0*(1.0+xp_quad_4[gp][1])*0.25;
      ds_quad_4[3][1][gp] = +1.0*(1.0-xp_quad_4[gp][0])*0.25;
    }


    xp_tetra_4 = calloc(4,sizeof(double*));
    wp_tetra_4 = calloc(4,sizeof(double));
    sh_tetra_4 = calloc(4,sizeof(double*));
    ds_tetra_4 = calloc(4,sizeof(double**));

    for (int i = 0 ; i < 4 ; i++)
      ds_tetra_4[i] = calloc(3,sizeof(double*));

    for (int i = 0 ; i < 4 ; i++) {
      xp_tetra_4[i] = calloc(3,sizeof(double));
      sh_tetra_4[i] = calloc(4,sizeof(double));
      for (int j = 0 ; j < 3 ; j++)
	ds_tetra_4[i][j] = calloc(4,sizeof(double));
    }

    xp_tetra_4[0][0] = +0.138196601125011;
    xp_tetra_4[0][1] = +0.138196601125011;
    xp_tetra_4[0][2] = +0.138196601125011;
    wp_tetra_4[0] = +0.041666666666666;

    xp_tetra_4[1][0] = +0.585410196624968;
    xp_tetra_4[1][1] = +0.138196601125011;
    xp_tetra_4[1][2] = +0.138196601125011;
    wp_tetra_4[1] = +0.041666666666666;

    xp_tetra_4[2][0] = +0.138196601125011;
    xp_tetra_4[2][1] = +0.585410196624968;
    xp_tetra_4[2][2] = +0.138196601125011;
    wp_tetra_4[2] = +0.041666666666666;

    xp_tetra_4[3][0] = +0.138196601125011;
    xp_tetra_4[3][1] = +0.138196601125011;
    xp_tetra_4[3][2] = +0.585410196624968;
    wp_tetra_4[3] = +0.041666666666666;

    for (int gp = 0 ; gp < 4 ; gp++) {

      sh_tetra_4[0][gp] = 1.0 - xp_tetra_4[gp][0] - xp_tetra_4[gp][1] - xp_tetra_4[gp][2];
      ds_tetra_4[0][0][gp] = -1.0;
      ds_tetra_4[0][1][gp] = -1.0;
      ds_tetra_4[0][2][gp] = -1.0;

      sh_tetra_4[1][gp] = xp_tetra_4[gp][0];
      ds_tetra_4[1][0][gp] = +1.0;
      ds_tetra_4[1][1][gp] = +0.0;
      ds_tetra_4[1][2][gp] = +0.0;

      sh_tetra_4[2][gp] = xp_tetra_4[gp][1];
      ds_tetra_4[2][0][gp] = +0.0;
      ds_tetra_4[2][1][gp] = +1.0;
      ds_tetra_4[2][2][gp] = +0.0;

      sh_tetra_4[3][gp]= xp_tetra_4[gp][2];
      ds_tetra_4[3][0][gp]= +0.0;
      ds_tetra_4[3][1][gp]= +0.0;
      ds_tetra_4[3][2][gp]= +1.0;
    }

    xp_prism_6 = calloc(6,sizeof(double*));
    wp_prism_6 = calloc(6,sizeof(double));
    sh_prism_6 = calloc(6,sizeof(double*));
    ds_prism_6 = calloc(6,sizeof(double**));

    for (int i = 0 ; i < 6 ; i++)
      ds_prism_6[i] = calloc(3,sizeof(double*));

    for (int i = 0 ; i < 6 ; i++) {
      xp_prism_6[i] = calloc(3,sizeof(double));
      sh_prism_6[i] = calloc(6,sizeof(double));
      for (int j = 0 ; j < 3 ; j++)
	ds_prism_6[i][j] = calloc(6,sizeof(double));
    }

    xp_prism_6[0][0] = +0.166666666666666;
    xp_prism_6[0][1] = +0.166666666666666;
    xp_prism_6[0][2] = -0.577350269189626;//+0.166666666666666;
    wp_prism_6[0] = +0.166666666666666;

    xp_prism_6[1][0] = +0.666666666666666;
    xp_prism_6[1][1] = +0.166666666666666;
    xp_prism_6[1][2] = -0.577350269189626;//+0.166666666666666;
    wp_prism_6[1] = +0.166666666666666;

    xp_prism_6[2][0] = +0.166666666666666;
    xp_prism_6[2][1] = +0.666666666666666;
    xp_prism_6[2][2] = -0.577350269189626;//+0.166666666666666;
    wp_prism_6[2] = +0.166666666666666;

    xp_prism_6[3][0] = +0.166666666666666;
    xp_prism_6[3][1] = +0.166666666666666;
    xp_prism_6[3][2] = +0.577350269189626;//+0.666666666666666;
    wp_prism_6[3] = +0.166666666666666;

    xp_prism_6[4][0] = +0.666666666666666;
    xp_prism_6[4][1] = +0.166666666666666;
    xp_prism_6[4][2] = +0.577350269189626;//+0.666666666666666;
    wp_prism_6[4] = +0.166666666666666;

    xp_prism_6[5][0] = +0.166666666666666;
    xp_prism_6[5][1] = +0.666666666666666;
    xp_prism_6[5][2] = +0.577350269189626;//+0.666666666666666;
    wp_prism_6[5] = +0.166666666666666;

    for (int gp = 0 ; gp < 6 ; gp++) {

      sh_prism_6[0][gp] = (1.0-xp_prism_6[gp][0]-xp_prism_6[gp][1])*(1.0-xp_prism_6[gp][2])*0.5;
      ds_prism_6[0][0][gp] = -1.0*(1.0-xp_prism_6[gp][2])*0.5;
      ds_prism_6[0][1][gp] = -1.0*(1.0-xp_prism_6[gp][2])*0.5;
      ds_prism_6[0][2][gp] = -1.0*(1.0-xp_prism_6[gp][0]-xp_prism_6[gp][1])*0.5;

      sh_prism_6[1][gp] =  xp_prism_6[gp][0]*(1.0-xp_prism_6[gp][2])*0.5;
      ds_prism_6[1][0][gp] = +1.0*(1.0-xp_prism_6[gp][2])*0.5;
      ds_prism_6[1][1][gp] = +0.0;
      ds_prism_6[1][2][gp] = -1.0*xp_prism_6[gp][0]*0.5;

      sh_prism_6[2][gp] = xp_prism_6[gp][1]*(1.0-xp_prism_6[gp][2])*0.5;
      ds_prism_6[2][0][gp] = +0.0;
      ds_prism_6[2][1][gp] = +1.0*(1.0-xp_prism_6[gp][2])*0.5;
      ds_prism_6[2][2][gp] = -1.0*xp_prism_6[gp][1]*0.5;

      sh_prism_6[3][gp] = (1.0-xp_prism_6[gp][0]-xp_prism_6[gp][1])*(1.0+xp_prism_6[gp][2])*0.5;
      ds_prism_6[3][0][gp] = -1.0*(1.0+xp_prism_6[gp][2])*0.5;
      ds_prism_6[3][1][gp] = -1.0*(1.0+xp_prism_6[gp][2])*0.5;
      ds_prism_6[3][2][gp] = +1.0*(1.0-xp_prism_6[gp][0]-xp_prism_6[gp][1])*0.5;

      sh_prism_6[4][gp] = xp_prism_6[gp][0]*(1.0 + xp_prism_6[gp][2])*0.5;
      ds_prism_6[4][0][gp] = +1.0*(1.0+xp_prism_6[gp][2])*0.5;
      ds_prism_6[4][1][gp] = +0.0;
      ds_prism_6[4][2][gp] = +1.0*xp_prism_6[gp][0]*0.5;

      sh_prism_6[5][gp] = xp_prism_6[gp][1]*(1.0+xp_prism_6[gp][2])*0.5;
      ds_prism_6[5][0][gp] = +0.0;
      ds_prism_6[5][1][gp] = +1.0*(1.0+xp_prism_6[gp][2])*0.5;
      ds_prism_6[5][2][gp] = +1.0*xp_prism_6[gp][1]*0.5;
    }

    xp_hexa_8 = calloc(8,sizeof(double*));
    wp_hexa_8 = calloc(8,sizeof(double));
    sh_hexa_8 = calloc(8,sizeof(double*));
    ds_hexa_8 = calloc(8,sizeof(double**));

    for (int i = 0 ; i < 8 ; i++)
      ds_hexa_8[i] = calloc(3,sizeof(double*));

    for (int i = 0 ; i < 8 ; i++) {
      xp_hexa_8[i] = calloc(3,sizeof(double));
      sh_hexa_8[i] = calloc(8,sizeof(double));
      for (int j = 0 ; j < 3 ; j++)
	ds_hexa_8[i][j] = calloc(8,sizeof(double));
    }

    xp_hexa_8[0][0] = -0.577350269189626;
    xp_hexa_8[0][1] = -0.577350269189626;
    xp_hexa_8[0][2] = -0.577350269189626;
    wp_hexa_8[0] = +1.0;

    xp_hexa_8[1][0]= +0.577350269189626;
    xp_hexa_8[1][1]= -0.577350269189626;
    xp_hexa_8[1][2]= -0.577350269189626;
    wp_hexa_8[1] = +1.0;

    xp_hexa_8[2][0]= +0.577350269189626;
    xp_hexa_8[2][1]= +0.577350269189626;
    xp_hexa_8[2][2]= -0.577350269189626;
    wp_hexa_8[2] = +1.0;

    xp_hexa_8[3][0]= -0.577350269189626;
    xp_hexa_8[3][1]= +0.577350269189626;
    xp_hexa_8[3][2]= -0.577350269189626;
    wp_hexa_8[3] = +1.0;

    xp_hexa_8[4][0]= -0.577350269189626;
    xp_hexa_8[4][1]= -0.577350269189626;
    xp_hexa_8[4][2]= +0.577350269189626;
    wp_hexa_8[4] = +1.0;

    xp_hexa_8[5][0]= +0.577350269189626;
    xp_hexa_8[5][1]= -0.577350269189626;
    xp_hexa_8[5][2]= +0.577350269189626;
    wp_hexa_8[5] = +1.0;

    xp_hexa_8[6][0]= +0.577350269189626;
    xp_hexa_8[6][1]= +0.577350269189626;
    xp_hexa_8[6][2]= +0.577350269189626;
    wp_hexa_8[6] = +1.0;

    xp_hexa_8[7][0]= -0.577350269189626;
    xp_hexa_8[7][1]= +0.577350269189626;
    xp_hexa_8[7][2]= +0.577350269189626;
    wp_hexa_8[7] = +1.0;

    for (int gp = 0 ; gp < 8 ; gp++) {

      sh_hexa_8[0][gp] = ( 1.0 - xp_hexa_8[gp][0] )*(1.0 - xp_hexa_8[gp][1] )*(1.0 - xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[0][0][gp] = -1.0*(1.0 - xp_hexa_8[gp][1])*(1.0 - xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[0][1][gp] = -1.0*(1.0 - xp_hexa_8[gp][0])*(1.0 - xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[0][2][gp] = -1.0*(1.0 - xp_hexa_8[gp][0])*(1.0 - xp_hexa_8[gp][1])*0.125;

      sh_hexa_8[1][gp] = ( 1.0 + xp_hexa_8[gp][0] )*(1.0 - xp_hexa_8[gp][1] )*(1.0 - xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[1][0][gp] = +1.0*(1.0 - xp_hexa_8[gp][1])*(1.0 - xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[1][1][gp] = -1.0*(1.0 + xp_hexa_8[gp][0])*(1.0 - xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[1][2][gp] = -1.0*(1.0 + xp_hexa_8[gp][0])*(1.0 - xp_hexa_8[gp][1])*0.125;

      sh_hexa_8[2][gp] = ( 1.0 + xp_hexa_8[gp][0] )*(1.0 + xp_hexa_8[gp][1] )*(1.0 - xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[2][0][gp] = +1.0*(1.0 + xp_hexa_8[gp][1])*(1.0 - xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[2][1][gp] = +1.0*(1.0 + xp_hexa_8[gp][0])*(1.0 - xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[2][2][gp] = -1.0*(1.0 + xp_hexa_8[gp][0])*(1.0 + xp_hexa_8[gp][1])*0.125;

      sh_hexa_8[3][gp] = ( 1.0 - xp_hexa_8[gp][0] )*(1.0 + xp_hexa_8[gp][1] )*(1.0 - xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[3][0][gp] = -1.0*(1.0 + xp_hexa_8[gp][1])*(1.0 - xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[3][1][gp] = +1.0*(1.0 - xp_hexa_8[gp][0])*(1.0 - xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[3][2][gp] = -1.0*(1.0 - xp_hexa_8[gp][0])*(1.0 + xp_hexa_8[gp][1])*0.125;

      sh_hexa_8[4][gp] = ( 1.0 - xp_hexa_8[gp][0] )*(1.0 - xp_hexa_8[gp][1] )*(1.0 + xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[4][0][gp] = -1.0*(1.0 - xp_hexa_8[gp][1])*(1.0 + xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[4][1][gp] = -1.0*(1.0 - xp_hexa_8[gp][0])*(1.0 + xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[4][2][gp] = +1.0*(1.0 - xp_hexa_8[gp][0])*(1.0 - xp_hexa_8[gp][1])*0.125;

      sh_hexa_8[5][gp] = ( 1.0 + xp_hexa_8[gp][0] )*(1.0 - xp_hexa_8[gp][1] )*(1.0 + xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[5][0][gp] = +1.0*(1.0 - xp_hexa_8[gp][1])*(1.0 + xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[5][1][gp] = -1.0*(1.0 + xp_hexa_8[gp][0])*(1.0 + xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[5][2][gp] = +1.0*(1.0 + xp_hexa_8[gp][0])*(1.0 - xp_hexa_8[gp][1])*0.125;

      sh_hexa_8[6][gp] = ( 1.0 + xp_hexa_8[gp][0] )*(1.0 + xp_hexa_8[gp][1] )*(1.0 + xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[6][0][gp] = +1.0*(1.0 + xp_hexa_8[gp][1])*(1.0 + xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[6][1][gp] = +1.0*(1.0 + xp_hexa_8[gp][0])*(1.0 + xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[6][2][gp] = +1.0*(1.0 + xp_hexa_8[gp][0])*(1.0 + xp_hexa_8[gp][1])*0.125;

      sh_hexa_8[7][gp] = ( 1.0 - xp_hexa_8[gp][0] )*(1.0 + xp_hexa_8[gp][1] )*(1.0 + xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[7][0][gp] = -1.0*(1.0 + xp_hexa_8[gp][1])*(1.0 + xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[7][1][gp] = +1.0*(1.0 - xp_hexa_8[gp][0])*(1.0 + xp_hexa_8[gp][2])*0.125;
      ds_hexa_8[7][2][gp] = +1.0*(1.0 - xp_hexa_8[gp][0])*(1.0 + xp_hexa_8[gp][1])*0.125;
    }

    return 0;
}


int fem_calc_jac(int dim, int npe, int gp, double * coor, double *** dsh, double ** jac) {

  if (jac == NULL || coor == NULL || dsh == NULL)
    return 1;

  for (int i = 0 ; i < dim ; i++) {
    for (int j = 0 ; j < dim ; j++) {
      jac[i][j] = 0.0;
      for (int n = 0 ; n < npe ; n++)
	jac[i][j] += dsh[n][i][gp]*coor[n*dim + j];
    }
  }
  return 0;
}


int fem_invjac(int dim, double **jac, double **ijac, double *det) {

  double c00, c01, c02, c10, c11, c12, c20, c21, c22;

  if (jac == NULL || ijac == NULL)
    return 1;

  if (dim == 2) {
    c00 = +jac[1][1];
    c01 = -jac[1][0];
    c10 = -jac[0][1];
    c11 = +jac[0][0];

    (*det) = jac[0][0]*c00 + jac[0][1]*c01;
    ijac[0][0] = c00/(*det);
    ijac[0][1] = c10/(*det);
    ijac[1][0] = c01/(*det);
    ijac[1][1] = c11/(*det);

  }else if (dim == 3) {
    c00 = +jac[1][1]*jac[2][2]-jac[2][1]*jac[1][2];
    c01 = -jac[1][0]*jac[2][2]+jac[2][0]*jac[1][2];
    c02 = +jac[1][0]*jac[2][1]-jac[2][0]*jac[1][1];
    c10 = -jac[0][1]*jac[2][2]+jac[2][1]*jac[0][2];
    c11 = +jac[0][0]*jac[2][2]-jac[2][0]*jac[0][2];
    c12 = -jac[0][0]*jac[2][1]+jac[2][0]*jac[0][1];
    c20 = +jac[0][1]*jac[1][2]-jac[1][1]*jac[0][2];
    c21 = -jac[0][0]*jac[1][2]+jac[1][0]*jac[0][2];
    c22 = +jac[0][0]*jac[1][1]-jac[1][0]*jac[0][1];

    (*det) = jac[0][0]*c00 + jac[0][1]*c01 + jac[0][2]*c02;
    ijac[0][0] = c00/(*det);
    ijac[0][1] = c10/(*det);
    ijac[0][2] = c20/(*det);
    ijac[1][0] = c01/(*det);
    ijac[1][1] = c11/(*det);
    ijac[1][2] = c21/(*det);
    ijac[2][0] = c02/(*det);
    ijac[2][1] = c12/(*det);
    ijac[2][2] = c22/(*det);
  }

  return 0;
}


int fem_trans_dsh(int dim, int nsh, int gp, double **ijac, double ***dsh_master, double ***dsh) {

  if (ijac == NULL || dsh_master == NULL || dsh == NULL) return 1;

  for (int sh = 0 ; sh < nsh ; sh++) {
    for (int i = 0 ; i < dim ; i++) {
      dsh[sh][i][gp] = 0.0;
      for (int j = 0 ; j < dim ; j++)
        dsh[sh][i][gp] += ijac[i][j] * dsh_master[sh][j][gp];
    }
  }
  return 0;
}


int fem_get_dsh_master(int npe, int dim, double ****dsh) {

  switch(dim) {

    case 2:

      switch(npe) {

	case 3:
	  *dsh = ds_tria_3;
	  break;

	case 4:
	  *dsh = ds_quad_4;
	  break;

	default:
	  return 1;
      }
      break;

    case 3:

      switch(npe) {

	case 4:
	  *dsh = ds_tetra_4;
	  break;

	case 6:
	  *dsh = ds_prism_6;
	  break;

	case 8:
	  *dsh = ds_hexa_8;
	  break;

	default:
	  return 1;
      }
      break;

    default:
      return 1;
  }

  return 0;
}


int fem_get_sh(int npe, int dim, double ***sh) {

  switch(dim) {

    case 2:

      switch(npe) {

        case 3:
          *sh = sh_tria_3;
	  break;

        case 4:
          *sh = sh_quad_4;
	  break;

        default:
          return 1;
      }
      break;

    case 3:

      switch(npe) {

        case 4:
          *sh = sh_tetra_4;
	  break;

        case 6:
          *sh = sh_prism_6;
	  break;

        case 8:
          *sh = sh_hexa_8;
	  break;

        default:
          return 1;
      }
      break;

    default:
      return 1;
  }

  return 0;
}


int fem_get_wp(int npe, int dim, double **wp) {

  switch(dim) {

    case 2:

      switch(npe) {

        case 3:
          *wp = wp_tria_3;
	  break;

        case 4:
          *wp = wp_quad_4;
	  break;

        default:
	  return 1;
      }
      break;

    case 3:

      switch(npe) {

        case 4:
          *wp = wp_tetra_4;
	  break;

        case 6:
          *wp = wp_prism_6;
	  break;

        case 8:
          *wp = wp_hexa_8;
	  break;

        default:
	  return 1;
      }
      break;

    default:
      return 1;
  }

  return 0;
}
