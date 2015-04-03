#ifndef GRIDH
#define GRIDH

struct grid
{
    int nx;
    int ng;
    double xmin;
    double xmax;
    double *x;
    double *prim;
    double *cons;
}

void make_grid(struct parList *pars);

#endif
