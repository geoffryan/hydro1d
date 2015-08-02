#ifndef HYDRO1D_INITIAL
#define HYDRO1D_INITIAL

void (*initial_value)(double *, double, struct parList *);

int set_initial(struct parList *pars);
void initialize_grid(struct grid *g, struct parList *pars);
void initial_uniform(double *prim, double x, struct parList *pars);
void initial_shocktube(double *prim, double x, struct parList *pars);
void initial_shocktube_transverse(double *prim, double x, 
                                struct parList *pars);
void initial_isentrope(double *prim, double x, struct parList *pars);
void initial_disc(double *prim, double x, struct parList *pars);

#endif
