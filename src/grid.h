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
    double *x;
    double *prim;
    double *cons;
    double *cons_rk;
    double *grad;
};

const static struct grid GRID_DEFAULT = {
    .nx = 0,
    .nx_int = 0,
    .ng = 0,
    .nc = 0,
    .nq = 0,
    .xmin = 0.0,
    .xmax = 1.0,
    .x = NULL,
    .prim = NULL,
    .cons = NULL,
    .cons_rk = NULL,
    .grad = NULL
};

void make_grid(struct grid *g, struct parList *pars);
void print_grid(struct grid *g, char *filename);
void free_grid(struct grid *g);

#endif
