#ifndef HYDRO1D_GEOM
#define HYDRO1D_GEOM

double (*CM)(double x0, double x1);
double (*DA)(double x);
double (*DV)(double x0, double x1);

double CM_cart(double x0, double x1);
double DA_cart(double x);
double DV_cart(double x0, double x1);

double CM_cyl(double x0, double x1);
double DA_cyl(double x);
double DV_cyl(double x0, double x1);

double CM_sph(double x0, double x1);
double DA_sph(double x);
double DV_sph(double x0, double x1);

#endif
