#ifndef HYDRO1D_RIEMANN
#define HYDRO1D_RIEMANN

void riemann_solve(double *primL, double *consL, double xL,
                    double *primR, double *consR, double xR, 
                    double x, double dt, struct parList *pars);

#endif
