#ifndef HYDRO1D_HYDRO
#define HYDRO1D_HYDRO

void prim2cons(double *prim, double *cons, double x, double dV,
                struct parList *pars);
void cons2prim(double *cons, double *prim, double x, double dV,
                struct parList *pars);
void flux(double *prim, double *F, double x, struct parList *pars);
void add_source(double *prim, double *cons, double x, double dVdt, 
                struct parList *pars);
void wave_speeds(double *prim1, double *prim2, double *sL, double *sR,
                    double *sC, double x, struct parList *pars);
double mindt(double *prim, double x, double dx, struct parList *pars);

#endif
