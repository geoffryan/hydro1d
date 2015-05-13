#ifndef HYDRO1D_GEOM
#define HYDRO1D_GEOM

double CM(double x0, double x1);
double DA(double x);
double DV(double x0, double x1);

#if GEOM == 0
//Cartesian Geometry
double CM(double x0, double x1)
{
    return 0.5*(x0+x1);
}

double DA(double x)
{
    return 1.0;
}

double DV(double x0, double x1)
{
    return x1-x0;
}

#elif GEOM == 1
//Cylindrical Geometry
double CM(double x0, double x1)
{
    return 0.5*(x0+x1);
}

double DA(double x)
{
    return 2*M_PI*x;
}

double DV(double x0, double x1)
{
    return M_PI*(x0+x1)*(x1-x0);
}

#elif GEOM == 2
//Spherical Geometry
double CM(double x0, double x1)
{
    return 2*(x1*x1+x1*x0+x0*x0)/(3*(x1+x0));
}

double DA(double x)
{
    return 4*M_PI*x*x;
}

double DV(double x0, double x1)
{
    return 4.0*M_PI*(x1*x1+x0*x1+x0*x0)*(x1-x0)/3.0;
}
#endif


#endif
