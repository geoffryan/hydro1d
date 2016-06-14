#include <math.h>
#include "geom.h"

//Cartesian Geometry
double CM_cart(double x0, double x1)
{
    return 0.5*(x0+x1);
}

double DA_cart(double x)
{
    return 1.0;
}

double DV_cart(double x0, double x1)
{
    return x1-x0;
}

//Cylindrical Geometry
double CM_cyl(double x0, double x1)
{
    return 0.5*(x0+x1);
}

double DA_cyl(double x)
{
    return 2*M_PI*x;
}

double DV_cyl(double x0, double x1)
{
    return M_PI*(x0+x1)*(x1-x0);
}

//Spherical Geometry
double CM_sph(double x0, double x1)
{
    return 2*(x1*x1+x1*x0+x0*x0)/(3*(x1+x0));
}

double DA_sph(double x)
{
    return 4*M_PI*x*x;
}

double DV_sph(double x0, double x1)
{
    return 4.0*M_PI*(x1*x1+x0*x1+x0*x0)*(x1-x0)/3.0;
}
