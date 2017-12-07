#ifndef HOMOGENIZE_H
#define HOMOGENIZE_H

#define HOMOG_METHOD_NULL         0
#define HOMOG_METHOD_TAYLOR_S     1
#define HOMOG_METHOD_TAYLOR_P     2
#define HOMOG_METHOD_UNIF_STRAINS 3

#define HOMOGENIZE_DELTA_STRAIN 0.005

int homogenize_init(void);

int homogenize_get_c_tangent(double *strain_mac, double **c_tangent);
int homogenize_calculate_c_tangent(double *strain_mac, double *c_tangent);
int homogenize_calculate_c_tangent_around_zero(double *c_tangent);

int homogenize_get_strain_stress(double *strain_mac, double *strain_ave, double *stress_ave);
int homogenize_get_strain_stress_non_linear(double *strain_mac, double *strain_ave, double *stress_ave);

int mic_homogenize_taylor(double *strain_mac, double *strain_ave, double *stress_ave);
int mic_homog_us(double *strain_mac, double *strain_ave, double *stress_ave);

#endif
