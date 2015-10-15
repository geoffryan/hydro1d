#ifndef HYDRO1D_RIEMANN
#define HYDRO1D_RIEMANN

void (*riemann_flux)(double primL[], double primR[], double F[], int nq,
                    double x, double w, struct parList *pars);

int set_riemann_solver(struct parList *pars);

void lax_friedrichs_flux(double primL[], double primR[], double F[], int nq,
                            double x, double w, struct parList *pars);
void hll_flux(double primL[], double primR[], double F[], int nq,
                 double x, double w, struct parList *pars);

#endif
