#ifndef HOMOGENIZE_H
#define HOMOGENIZE_H

#define HOMOG_METHOD_NULL         0
#define HOMOG_METHOD_TAYLOR_S     1
#define HOMOG_METHOD_TAYLOR_P     2
#define HOMOG_METHOD_UNIF_STRAINS 3

int homogenize_init(void);
int homogenize_get_average_c_tangent(double *strain_mac, double **c_tangent);
int homogenize_get_average_c_tangent_non_linear(double *strain_mac, double *c_tangent);
int homogenize_get_average_strain_stress(double *strain_mac, double *strain_ave, double *stress_ave);
int homogenize_get_average_strain_stress_non_linear(double *strain_mac, double *strain_ave, double *stress_ave);
int mic_homogenize_taylor(double *strain_mac, double *strain_ave, double *stress_ave);
int mic_homog_us(double *strain_mac, double *strain_ave, double *stress_ave);

#endif
