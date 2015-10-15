#ifndef HYDRO1D_GRID
#define HYDRO1D_GRID

struct parList;

struct grid
{
    int nx;
    int nx_int;
    int ng;
    int nc;
    int nq;
    double xmin;
    double xmax;
    double t;
    double *x;
    double *x_rk;
    double *w;
    double *prim;
    double *cons;
    double *cons_rk;
    double *grad;
    double PLM;
};

const static struct grid GRID_DEFAULT = {
    .nx = 0,
    .nx_int = 0,
    .ng = 0,
    .nc = 0,
    .nq = 0,
    .xmin = 0.0,
    .xmax = 1.0,
    .t = 0.0,
    .x = NULL,
    .x_rk = NULL,
    .w = NULL,
    .prim = NULL,
    .cons = NULL,
    .cons_rk = NULL,
    .grad = NULL,
    .PLM = 1.5
};

void (*reconstruction)(struct grid *);

int set_reconstruction(struct parList *pars);
void make_grid(struct grid *g, struct parList *pars);
void free_grid(struct grid *g);

void interpolate_constant(struct grid *g);
void interpolate_plm(struct grid *g);
void copy_to_rk(struct grid *g);
void update_cons(struct grid *g, double fac1, double fac2);
void update_cons_rk(struct grid *g, double fac1, double fac2);
void update_x(struct grid *g, double fac1, double fac2);
void update_x_rk(struct grid *g, double fac1, double fac2);

#endif
